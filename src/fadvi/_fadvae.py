from __future__ import annotations

from typing import TYPE_CHECKING

import torch
import torch.nn.functional as F
from scvi import REGISTRY_KEYS
from scvi.module._classifier import Classifier
from scvi.module._vae import VAE
from scvi.module.base import LossOutput, SupervisedModuleClass
from scvi.nn import Encoder
from torch.distributions import Normal
from torch.distributions import kl_divergence as kl

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Literal

    from torch.distributions import Distribution


class GradientReversalFunction(torch.autograd.Function):
    """Gradient reversal layer for adversarial training."""

    @staticmethod
    def forward(ctx, x, alpha):
        ctx.alpha = alpha
        return x.view_as(x)

    @staticmethod
    def backward(ctx, grad_output):
        return grad_output.neg() * ctx.alpha, None


def gradient_reversal(x, alpha=1.0):
    """Apply gradient reversal to input tensor."""
    return GradientReversalFunction.apply(x, alpha)


class FADVAE(SupervisedModuleClass, VAE):
    """Factor Disentanglement Variational Autoencoder.

    This model disentangles batch-related variation (z_b), label-related variation (z_l),
    and residual variation (z_r) using adversarial training and cross-correlation penalties.

    Parameters
    ----------
    n_input
        Number of input genes
    n_batch
        Number of batches
    n_labels
        Number of labels
    n_hidden
        Number of nodes per hidden layer
    n_latent_b
        Dimensionality of batch latent space
    n_latent_l
        Dimensionality of label latent space
    n_latent_r
        Dimensionality of residual latent space
    n_layers
        Number of hidden layers used for encoder and decoder NNs
    n_continuous_cov
        Number of continuous covariates
    n_cats_per_cov
        Number of categories for each extra categorical covariate
    dropout_rate
        Dropout rate for neural networks
    dispersion
        Dispersion parameter option
    log_variational
        Log(data+1) prior to encoding for numerical stability
    gene_likelihood
        One of 'nb', 'zinb', 'poisson'
    use_observed_lib_size
        If True, use observed library size
    beta
        KL divergence weight
    lambda_b
        Weight for batch classification loss
    lambda_l
        Weight for label classification loss
    alpha_bl
        Weight for adversarial loss (label prediction from batch latents)
    alpha_lb
        Weight for adversarial loss (batch prediction from label latents)
    alpha_rb
        Weight for adversarial loss (batch prediction from residual latents)
    alpha_rl
        Weight for adversarial loss (label prediction from residual latents)
    gamma
        Weight for cross-correlation penalty
    use_batch_norm
        Whether to use batch norm in layers
    use_layer_norm
        Whether to use layer norm in layers
    **vae_kwargs
        Keyword args for VAE
    """

    def __init__(
        self,
        n_input: int,
        n_batch: int = 0,
        n_labels: int = 0,
        n_hidden: int = 128,
        n_latent_b: int = 30,
        n_latent_l: int = 30,
        n_latent_r: int = 10,
        n_layers: int = 1,
        n_continuous_cov: int = 0,
        n_cats_per_cov: Iterable[int] | None = None,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        log_variational: bool = True,
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb",
        use_observed_lib_size: bool = True,
        beta: float = 1.0,
        lambda_b: float = 50,
        lambda_l: float = 50,
        alpha_bl: float = 1.0,
        alpha_lb: float = 1.0,
        alpha_rb: float = 1.0,
        alpha_rl: float = 1.0,
        gamma: float = 1.0,
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        unlabeled_category_id: int | None = None,
        **vae_kwargs,
    ):
        # Calculate total latent dimension for VAE parent
        n_latent_total = n_latent_b + n_latent_l + n_latent_r

        super().__init__(
            n_input,
            n_hidden=n_hidden,
            n_latent=n_latent_total,
            n_layers=n_layers,
            n_continuous_cov=n_continuous_cov,
            n_cats_per_cov=n_cats_per_cov,
            dropout_rate=dropout_rate,
            n_batch=n_batch,
            dispersion=dispersion,
            log_variational=log_variational,
            gene_likelihood=gene_likelihood,
            use_observed_lib_size=use_observed_lib_size,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            **vae_kwargs,
        )

        # Store dimensions
        self.n_latent_b = n_latent_b
        self.n_latent_l = n_latent_l
        self.n_latent_r = n_latent_r
        self.n_labels = n_labels
        self.n_cats_per_cov = n_cats_per_cov
        self._unlabeled_category_id = unlabeled_category_id

        # Store loss weights
        self.beta = beta
        self.lambda_b = lambda_b
        self.lambda_l = lambda_l
        self.alpha_bl = alpha_bl
        self.alpha_lb = alpha_lb
        self.alpha_rb = alpha_rb
        self.alpha_rl = alpha_rl
        self.gamma = gamma

        use_batch_norm_encoder = use_batch_norm == "encoder" or use_batch_norm == "both"
        use_layer_norm_encoder = use_layer_norm == "encoder" or use_layer_norm == "both"

        # Construct categorical list for encoders (following VAE pattern)
        cat_list = [n_batch] + list([] if n_cats_per_cov is None else n_cats_per_cov)

        # Separate encoders for each latent factor
        self.encoder_b = Encoder(
            n_input,
            n_latent_b,
            n_cat_list=cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            return_dist=True,
        )

        self.encoder_l = Encoder(
            n_input,
            n_latent_l,
            n_cat_list=cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            return_dist=True,
        )

        self.encoder_r = Encoder(
            n_input,
            n_latent_r,
            n_cat_list=cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            return_dist=True,
        )

        # Supervised classification heads
        self.head_batch = Classifier(
            n_latent_b,
            n_labels=n_batch,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            logits=True,
        )

        self.head_label = Classifier(
            n_latent_l,
            n_labels=n_labels,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            logits=True,
        )

        # Adversarial classification heads
        self.adv_head_batch_from_l = Classifier(
            n_latent_l,
            n_labels=n_batch,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            logits=True,
        )

        self.adv_head_label_from_b = Classifier(
            n_latent_b,
            n_labels=n_labels,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            logits=True,
        )

        self.adv_head_batch_from_r = Classifier(
            n_latent_r,
            n_labels=n_batch,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            logits=True,
        )

        self.adv_head_label_from_r = Classifier(
            n_latent_r,
            n_labels=n_labels,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            logits=True,
        )

    def _get_inference_input(self, tensors):
        """Get input for inference."""
        x = tensors[REGISTRY_KEYS.X_KEY]

        input_dict = {"x": x}

        # Add batch index if present
        if REGISTRY_KEYS.BATCH_KEY in tensors.keys():
            input_dict.update({"batch_index": tensors[REGISTRY_KEYS.BATCH_KEY]})

        # Add categorical covariates if present
        if REGISTRY_KEYS.CAT_COVS_KEY in tensors.keys():
            input_dict.update({"cat_covs": tensors[REGISTRY_KEYS.CAT_COVS_KEY]})

        # Add continuous covariates if present
        if REGISTRY_KEYS.CONT_COVS_KEY in tensors.keys():
            input_dict.update({"cont_covs": tensors[REGISTRY_KEYS.CONT_COVS_KEY]})

        return input_dict

    def _get_generative_input(self, tensors, inference_outputs):
        """Get input for generative model."""
        z_b = inference_outputs["z_b"]
        z_l = inference_outputs["z_l"]
        z_r = inference_outputs["z_r"]
        library = inference_outputs["library"]

        # Concatenate all latent factors
        z = torch.cat([z_b, z_l, z_r], dim=-1)

        input_dict = {"z": z, "library": library}

        # Add batch information if present
        if REGISTRY_KEYS.BATCH_KEY in tensors.keys():
            input_dict.update({"batch_index": tensors[REGISTRY_KEYS.BATCH_KEY]})

        # Add categorical covariates if present
        if REGISTRY_KEYS.CAT_COVS_KEY in tensors.keys():
            input_dict.update({"cat_covs": tensors[REGISTRY_KEYS.CAT_COVS_KEY]})

        # Add continuous covariates if present
        if REGISTRY_KEYS.CONT_COVS_KEY in tensors.keys():
            input_dict.update({"cont_covs": tensors[REGISTRY_KEYS.CONT_COVS_KEY]})

        return input_dict

    def inference(
        self, x, batch_index=None, cat_covs=None, cont_covs=None, n_samples=1
    ):
        """Run inference to get latent representations."""
        # Prepare input
        x_ = x
        if self.log_variational:
            x_ = torch.log1p(x_)

        # Ensure all tensors are on the same device as x
        device = x.device

        # Move batch_index to correct device if provided
        if batch_index is not None:
            batch_index = batch_index.to(device)

        # Move categorical covariates to correct device if provided
        if cat_covs is not None:
            cat_covs = [
                cov.to(device) if isinstance(cov, torch.Tensor) else cov
                for cov in cat_covs
            ]
            categorical_input = tuple(cat_covs)
        else:
            categorical_input = ()

        # Move continuous covariates to correct device if provided
        if cont_covs is not None:
            cont_covs = [
                cov.to(device) if isinstance(cov, torch.Tensor) else cov
                for cov in cont_covs
            ]

        # Encode to separate latent factors
        qb, z_b = self.encoder_b(x_, batch_index, *categorical_input)
        ql, z_l = self.encoder_l(x_, batch_index, *categorical_input)
        qr, z_r = self.encoder_r(x_, batch_index, *categorical_input)

        # Sample from distributions if needed
        if n_samples > 1:
            z_b = qb.rsample((n_samples,))
            z_l = ql.rsample((n_samples,))
            z_r = qr.rsample((n_samples,))

        # Library size encoding (from parent VAE)
        if not self.use_observed_lib_size:
            ql_lib, library_encoded = self.l_encoder(
                x_, batch_index, *categorical_input
            )
            if n_samples > 1:
                library = ql_lib.rsample((n_samples,))
            else:
                library = library_encoded
        else:
            library = torch.log(x.sum(1)).unsqueeze(1)
            ql_lib = None

        outputs = {
            "z_b": z_b,
            "z_l": z_l,
            "z_r": z_r,
            "qb": qb,
            "ql": ql,
            "qr": qr,
            "library": library,
        }

        if not self.use_observed_lib_size:
            outputs["ql_lib"] = ql_lib

        return outputs

    def generative(self, z, library, batch_index=None, cat_covs=None, cont_covs=None):
        """Run generative model to reconstruct data."""
        # Ensure all tensors are on the same device as z
        device = z.device

        # Move library to correct device
        if library is not None:
            library = library.to(device)

        # Move batch_index to correct device if provided
        if batch_index is not None:
            batch_index = batch_index.to(device)

        # Move categorical covariates to correct device if provided
        if cat_covs is not None:
            cat_covs = [
                cov.to(device) if isinstance(cov, torch.Tensor) else cov
                for cov in cat_covs
            ]

        # Move continuous covariates to correct device if provided
        if cont_covs is not None:
            cont_covs = [
                cov.to(device) if isinstance(cov, torch.Tensor) else cov
                for cov in cont_covs
            ]

        # Use the parent VAE's decoder with concatenated latents
        return super().generative(z, library, batch_index, cat_covs, cont_covs)

    def cross_covariance_penalty(self, z_b, z_l, z_r):
        """Compute cross-covariance penalty to enforce independence."""
        # Center the latent representations
        z_b_centered = z_b - z_b.mean(dim=0, keepdim=True)
        z_l_centered = z_l - z_l.mean(dim=0, keepdim=True)
        z_r_centered = z_r - z_r.mean(dim=0, keepdim=True)

        n = z_b.size(0)

        # Compute cross-covariances
        cov_bl = torch.mm(z_b_centered.t(), z_l_centered) / (n - 1)
        cov_br = torch.mm(z_b_centered.t(), z_r_centered) / (n - 1)
        cov_lr = torch.mm(z_l_centered.t(), z_r_centered) / (n - 1)

        # Return sum of squared cross-covariances (Frobenius norm)
        penalty = (cov_bl**2).sum() + (cov_br**2).sum() + (cov_lr**2).sum()

        return penalty

    def loss(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor | Distribution | None],
        generative_outputs: dict[str, Distribution | None],
        kl_weight: float = 1.0,
        labelled_tensors: dict[str, torch.Tensor] | None = None,
        **kwargs,  # Accept additional keyword arguments for compatibility
    ):
        """Compute the total loss including ELBO, supervised, adversarial, and decorrelation terms."""
        px: Distribution = generative_outputs["px"]
        x: torch.Tensor = tensors[REGISTRY_KEYS.X_KEY]
        batch_index: torch.Tensor = tensors[REGISTRY_KEYS.BATCH_KEY]
        labels: torch.Tensor = tensors[REGISTRY_KEYS.LABELS_KEY]

        # Get latent representations and distributions
        z_b = inference_outputs["z_b"]
        z_l = inference_outputs["z_l"]
        z_r = inference_outputs["z_r"]
        qb = inference_outputs["qb"]
        ql = inference_outputs["ql"]
        qr = inference_outputs["qr"]

        # Reconstruction loss - handle different gene likelihoods properly
        if self.gene_likelihood == "zinb" or self.gene_likelihood == "nb":
            # For negative binomial distributions, ensure x is non-negative integers
            x_safe = torch.clamp(torch.round(x), min=0.0)
        elif self.gene_likelihood == "poisson":
            # For Poisson, ensure x is non-negative
            x_safe = torch.clamp(x, min=0.0)
        else:
            # For other distributions, use general clamping
            x_safe = torch.clamp(x, min=1e-8)

        reconst_loss = -px.log_prob(x_safe).sum(-1)

        # KL divergences for each latent factor
        kl_divergence_b = kl(
            qb, Normal(torch.zeros_like(qb.loc), torch.ones_like(qb.scale))
        ).sum(dim=-1)
        kl_divergence_l = kl(
            ql, Normal(torch.zeros_like(ql.loc), torch.ones_like(ql.scale))
        ).sum(dim=-1)
        kl_divergence_r = kl(
            qr, Normal(torch.zeros_like(qr.loc), torch.ones_like(qr.scale))
        ).sum(dim=-1)

        kl_divergence = kl_divergence_b + kl_divergence_l + kl_divergence_r

        # Library size KL if not using observed library size
        if not self.use_observed_lib_size:
            ql_lib = inference_outputs["ql_lib"]
            (
                local_library_log_means,
                local_library_log_vars,
            ) = self._compute_local_library_params(batch_index)

            kl_divergence_lib = kl(
                ql_lib,
                Normal(local_library_log_means, torch.sqrt(local_library_log_vars)),
            ).sum(dim=1)
            kl_divergence += kl_divergence_lib

        # ELBO loss
        elbo_loss = torch.mean(reconst_loss + self.beta * kl_divergence * kl_weight)

        # Supervised classification losses
        pred_b = self.head_batch(z_b)
        pred_l = self.head_label(z_l)

        sup_loss_b = F.cross_entropy(pred_b, batch_index.squeeze())

        # For label prediction, handle potential unlabeled category
        # If the classifier has reduced outputs (n_labels < original),
        # the labels might need remapping to exclude unlabeled category
        labels_squeezed = labels.squeeze()
        if (
            hasattr(self, "_unlabeled_category_id")
            and self._unlabeled_category_id is not None
        ):
            # Filter out unlabeled data for classification loss
            labeled_mask = labels_squeezed != self._unlabeled_category_id
            if labeled_mask.sum() > 0:  # Only compute if there are labeled samples
                pred_l_labeled = pred_l[labeled_mask]
                labels_labeled = labels_squeezed[labeled_mask]

                # If labels need remapping (classifier outputs < total categories)
                if pred_l_labeled.size(1) < (self._unlabeled_category_id + 1):
                    # Remap labels: all indices >= unlabeled_category_id should be shifted down
                    labels_remapped = labels_labeled.clone()
                    labels_remapped[labels_labeled > self._unlabeled_category_id] -= 1
                    sup_loss_l = F.cross_entropy(pred_l_labeled, labels_remapped)
                else:
                    sup_loss_l = F.cross_entropy(pred_l_labeled, labels_labeled)
            else:
                sup_loss_l = torch.tensor(0.0, device=pred_l.device, requires_grad=True)
        else:
            sup_loss_l = F.cross_entropy(pred_l, labels_squeezed)

        sup_loss = self.lambda_b * sup_loss_b + self.lambda_l * sup_loss_l

        # Adversarial losses with gradient reversal
        adv_pred_bl = self.adv_head_label_from_b(gradient_reversal(z_b, self.alpha_bl))
        adv_pred_lb = self.adv_head_batch_from_l(gradient_reversal(z_l, self.alpha_lb))
        adv_pred_rb = self.adv_head_batch_from_r(gradient_reversal(z_r, self.alpha_rb))
        adv_pred_rl = self.adv_head_label_from_r(gradient_reversal(z_r, self.alpha_rl))

        # Apply same logic to adversarial losses
        if (
            hasattr(self, "_unlabeled_category_id")
            and self._unlabeled_category_id is not None
        ):
            labeled_mask = labels_squeezed != self._unlabeled_category_id
            if labeled_mask.sum() > 0:
                adv_pred_bl_labeled = adv_pred_bl[labeled_mask]
                adv_pred_rl_labeled = adv_pred_rl[labeled_mask]
                labels_labeled = labels_squeezed[labeled_mask]

                if adv_pred_bl_labeled.size(1) < (self._unlabeled_category_id + 1):
                    labels_remapped = labels_labeled.clone()
                    labels_remapped[labels_labeled > self._unlabeled_category_id] -= 1
                    adv_loss_bl = F.cross_entropy(adv_pred_bl_labeled, labels_remapped)
                    adv_loss_rl = F.cross_entropy(adv_pred_rl_labeled, labels_remapped)
                else:
                    adv_loss_bl = F.cross_entropy(adv_pred_bl_labeled, labels_labeled)
                    adv_loss_rl = F.cross_entropy(adv_pred_rl_labeled, labels_labeled)
            else:
                adv_loss_bl = torch.tensor(
                    0.0, device=adv_pred_bl.device, requires_grad=True
                )
                adv_loss_rl = torch.tensor(
                    0.0, device=adv_pred_rl.device, requires_grad=True
                )
        else:
            adv_loss_bl = F.cross_entropy(adv_pred_bl, labels_squeezed)
            adv_loss_rl = F.cross_entropy(adv_pred_rl, labels_squeezed)

        adv_loss_lb = F.cross_entropy(adv_pred_lb, batch_index.squeeze())
        adv_loss_rb = F.cross_entropy(adv_pred_rb, batch_index.squeeze())

        adv_loss = adv_loss_bl + adv_loss_lb + adv_loss_rb + adv_loss_rl

        # Cross-covariance penalty
        xcov_penalty = self.cross_covariance_penalty(z_b, z_l, z_r)
        xcov_loss = self.gamma * xcov_penalty

        # Total loss
        total_loss = elbo_loss + sup_loss + adv_loss + xcov_loss

        # For training plan metrics, filter out unlabeled data to match classifier outputs
        if (
            hasattr(self, "_unlabeled_category_id")
            and self._unlabeled_category_id is not None
        ):
            labels_squeezed = labels.squeeze()
            labeled_mask = labels_squeezed != self._unlabeled_category_id
            if labeled_mask.sum() > 0:
                # Filter logits and labels for metrics computation
                filtered_logits = pred_l[labeled_mask]
                filtered_labels = labels_squeezed[labeled_mask]

                # Remap labels if classifier has reduced outputs
                if filtered_logits.size(1) < (self._unlabeled_category_id + 1):
                    filtered_labels_remapped = filtered_labels.clone()
                    filtered_labels_remapped[
                        filtered_labels > self._unlabeled_category_id
                    ] -= 1
                    metrics_labels = filtered_labels_remapped
                else:
                    metrics_labels = filtered_labels

                metrics_logits = filtered_logits
                metrics_classification_loss = sup_loss
            else:
                # No labeled data in this batch - set classification_loss to None to skip metrics
                # metrics_logits = pred_l[:0]  # Empty tensor with correct dimensions
                # metrics_labels = labels.squeeze()[:0]  # Empty tensor
                # metrics_classification_loss = None
                return LossOutput(
                    loss=total_loss,
                    reconstruction_loss=reconst_loss,
                    kl_local=kl_divergence,
                    extra_metrics={
                        "elbo_loss": elbo_loss,
                        "sup_loss": sup_loss,
                        "sup_loss_b": sup_loss_b,
                        "sup_loss_l": sup_loss_l,
                        "adv_loss": adv_loss,
                        "adv_loss_bl": adv_loss_bl,
                        "adv_loss_lb": adv_loss_lb,
                        "adv_loss_rb": adv_loss_rb,
                        "adv_loss_rl": adv_loss_rl,
                        "xcov_loss": xcov_loss,
                        "z_b": z_b,
                        "z_l": z_l,
                        "z_r": z_r,
                    },
                )
        else:
            # No unlabeled category handling needed
            metrics_logits = pred_l
            metrics_labels = labels.squeeze()
            metrics_classification_loss = sup_loss

        return LossOutput(
            loss=total_loss,
            reconstruction_loss=reconst_loss,
            kl_local=kl_divergence,
            classification_loss=metrics_classification_loss,
            true_labels=metrics_labels,
            logits=metrics_logits,
            extra_metrics={
                "elbo_loss": elbo_loss,
                "sup_loss": sup_loss,
                "sup_loss_b": sup_loss_b,
                "sup_loss_l": sup_loss_l,
                "adv_loss": adv_loss,
                "adv_loss_bl": adv_loss_bl,
                "adv_loss_lb": adv_loss_lb,
                "adv_loss_rb": adv_loss_rb,
                "adv_loss_rl": adv_loss_rl,
                "xcov_loss": xcov_loss,
                "z_b": z_b,
                "z_l": z_l,
                "z_r": z_r,
            },
        )
