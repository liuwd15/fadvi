from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import torch
import torchmetrics.functional as tmf
from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LabelsWithUnlabeledObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
)
from scvi.model._utils import _init_library_size
from scvi.model.base import (
    ArchesMixin,
    BaseMinifiedModeModelClass,
    RNASeqMixin,
    SemisupervisedTrainingMixin,
    VAEMixin,
)
from scvi.module.base import BaseModuleClass, LossOutput
from scvi.train import SemiSupervisedTrainingPlan, TrainingPlan
from scvi.train._constants import METRIC_KEYS
from scvi.train._metrics import ElboMetric

from ._fadvae import FADVAE

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData

logger = logging.getLogger(__name__)


class SemiSupervisedTrainingPlanFixed(SemiSupervisedTrainingPlan):
    def __init__(
        self,
        module: BaseModuleClass,
        n_classes: int,
        *,
        lr: float = 1e-3,
        weight_decay: float = 1e-6,
        n_steps_kl_warmup: int | None = None,
        n_epochs_kl_warmup: int | None = 400,
        reduce_lr_on_plateau: bool = False,
        lr_factor: float = 0.6,
        lr_patience: int = 30,
        lr_threshold: float = 0.0,
        lr_scheduler_metric: Literal[
            "elbo_validation", "reconstruction_loss_validation", "kl_local_validation"
        ] = "elbo_validation",
        compile: bool = False,
        compile_kwargs: dict | None = None,
        **loss_kwargs,
    ):
        super().__init__(
            module=module,
            n_classes=n_classes,
            lr=lr,
            weight_decay=weight_decay,
            n_steps_kl_warmup=n_steps_kl_warmup,
            n_epochs_kl_warmup=n_epochs_kl_warmup,
            reduce_lr_on_plateau=reduce_lr_on_plateau,
            lr_factor=lr_factor,
            lr_patience=lr_patience,
            lr_threshold=lr_threshold,
            lr_scheduler_metric=lr_scheduler_metric,
            compile=compile,
            compile_kwargs=compile_kwargs,
            **loss_kwargs,
        )

    def compute_and_log_metrics(
        self, loss_output: LossOutput, metrics: dict[str, ElboMetric], mode: str
    ):
        """Computes and logs metrics."""
        TrainingPlan.compute_and_log_metrics(self, loss_output, metrics, mode)
        # no labeled observations in minibatch
        if loss_output.classification_loss is None:
            return

        classification_loss = loss_output.classification_loss
        true_labels = loss_output.true_labels
        logits = loss_output.logits
        predicted_labels = torch.argmax(logits, dim=-1)

        accuracy = tmf.classification.multiclass_accuracy(
            predicted_labels,
            true_labels,
            self.n_classes,
            average="micro",
        )

        f1 = tmf.classification.multiclass_f1_score(
            predicted_labels,
            true_labels,
            self.n_classes,
            average="micro",
        )
        ce = tmf.classification.multiclass_calibration_error(
            logits,
            true_labels,
            self.n_classes,
        )

        self.log_with_mode(
            METRIC_KEYS.CLASSIFICATION_LOSS_KEY,
            classification_loss,
            mode,
            on_step=False,
            on_epoch=True,
            batch_size=loss_output.n_obs_minibatch,
        )
        self.log_with_mode(
            METRIC_KEYS.ACCURACY_KEY,
            accuracy,
            mode,
            on_step=False,
            on_epoch=True,
            batch_size=loss_output.n_obs_minibatch,
        )
        self.log_with_mode(
            METRIC_KEYS.F1_SCORE_KEY,
            f1,
            mode,
            on_step=False,
            on_epoch=True,
            batch_size=loss_output.n_obs_minibatch,
        )
        self.log_with_mode(
            METRIC_KEYS.CALIBRATION_ERROR_KEY,
            ce,
            mode,
            on_step=False,
            on_epoch=True,
            batch_size=loss_output.n_obs_minibatch,
        )


class FADVI(
    RNASeqMixin,
    SemisupervisedTrainingMixin,
    VAEMixin,
    ArchesMixin,
    BaseMinifiedModeModelClass,
):
    """Factor Disentanglement Variational Inference model.

    This model disentangles batch-related variation, label-related variation,
    and residual variation using adversarial training and cross-correlation penalties.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.FADVI.setup_anndata`.
    n_hidden
        Number of nodes per hidden layer.
    n_latent_b
        Dimensionality of the batch latent space.
    n_latent_l
        Dimensionality of the label latent space.
    n_latent_r
        Dimensionality of the residual latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    dropout_rate
        Dropout rate for neural networks.
    dispersion
        One of the following:

        * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - dispersion can differ between different batches
        * ``'gene-label'`` - dispersion can differ between different labels
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell
    gene_likelihood
        One of:

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
        * ``'poisson'`` - Poisson distribution
    use_observed_lib_size
        If ``True``, use the observed library size for RNA as the scaling factor in the mean of the
        conditional distribution.
    beta
        Weight for KL divergence in ELBO.
    lambda_b
        Weight for batch classification loss.
    lambda_l
        Weight for label classification loss.
    alpha_bl
        Weight for adversarial loss (label prediction from batch latents).
    alpha_lb
        Weight for adversarial loss (batch prediction from label latents).
    alpha_rb
        Weight for adversarial loss (batch prediction from residual latents).
    alpha_rl
        Weight for adversarial loss (label prediction from residual latents).
    gamma
        Weight for cross-correlation penalty.
    **model_kwargs
        Keyword args for :class:`~fadvi.FADVAE`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> fadvi.FADVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    >>> model = fadvi.FADVI(adata)
    >>> model.train()
    >>> adata.obsm["X_fadvi_b"] = model.get_latent_representation(representation="b")
    >>> adata.obsm["X_fadvi_l"] = model.get_latent_representation(representation="l")
    >>> adata.obsm["X_fadvi_r"] = model.get_latent_representation(representation="r")
    """

    _training_plan_cls = SemiSupervisedTrainingPlanFixed

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent_b: int = 30,
        n_latent_l: int = 30,
        n_latent_r: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
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
        **model_kwargs,
    ):
        super().__init__(adata)

        self._set_indices_and_labels()
        self._set_batch_mapping()

        n_batch = self.summary_stats.n_batch

        # Following SCANVI pattern: reduce n_labels for module if unlabeled category exists
        n_labels = self.summary_stats.n_labels
        unlabeled_category_id = None
        if (
            self.unlabeled_category_ is not None
            and self.unlabeled_category_ in self.labels_
        ):
            n_labels = self.summary_stats.n_labels - 1
            # Find the ID of the unlabeled category in the label mapping
            unlabeled_category_id = np.where(
                self._label_mapping == self.unlabeled_category_
            )[0]
            if len(unlabeled_category_id) > 0:
                unlabeled_category_id = int(unlabeled_category_id[0])
            else:
                unlabeled_category_id = None

        # Store both original and reduced n_labels
        self.n_labels_original = (
            self.summary_stats.n_labels
        )  # For training plan compatibility
        self.n_labels = n_labels  # Reduced for module (actual prediction classes)
        n_continuous_cov = self.summary_stats.get("n_extra_continuous_covs", 0)
        n_cats_per_cov = (
            self.adata_manager.get_state_registry(
                REGISTRY_KEYS.CAT_COVS_KEY
            ).n_cats_per_key
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )

        # Complain if unlabeled
        use_size_factor_key = (
            REGISTRY_KEYS.SIZE_FACTOR_KEY in self.adata_manager.data_registry
        )
        library_log_means, library_log_vars = None, None
        if not use_size_factor_key:
            library_log_means, library_log_vars = _init_library_size(
                self.adata_manager, n_batch
            )

        self.module = FADVAE(
            n_input=self.summary_stats.n_vars,
            n_batch=n_batch,
            n_labels=n_labels,
            n_continuous_cov=n_continuous_cov,
            n_cats_per_cov=n_cats_per_cov,
            n_hidden=n_hidden,
            n_latent_b=n_latent_b,
            n_latent_l=n_latent_l,
            n_latent_r=n_latent_r,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            use_observed_lib_size=use_observed_lib_size,
            library_log_means=library_log_means,
            library_log_vars=library_log_vars,
            beta=beta,
            lambda_b=lambda_b,
            lambda_l=lambda_l,
            alpha_bl=alpha_bl,
            alpha_lb=alpha_lb,
            alpha_rb=alpha_rb,
            alpha_rl=alpha_rl,
            gamma=gamma,
            unlabeled_category_id=unlabeled_category_id,
            **model_kwargs,
        )
        self._model_summary_string = (
            f"FADVI Model with the following params: \nn_hidden: {n_hidden}, "
            f"n_latent_b: {n_latent_b}, n_latent_l: {n_latent_l}, n_latent_r: {n_latent_r}, "
            f"n_layers: {n_layers}, dropout_rate: {dropout_rate}, "
            f"dispersion: {dispersion}, gene_likelihood: {gene_likelihood}, "
            f"beta: {beta}, lambda_b: {lambda_b}, lambda_l: {lambda_l}, "
            f"alpha_bl: {alpha_bl}, alpha_lb: {alpha_lb}, alpha_rb: {alpha_rb}, alpha_rl: {alpha_rl}, "
            f"gamma: {gamma}"
        )
        self.init_params_ = self._get_init_params(locals())

    def _set_batch_mapping(self):
        """Set up batch mapping for converting codes to original batch labels."""
        if REGISTRY_KEYS.BATCH_KEY in self.adata_manager.data_registry:
            batch_state_registry = self.adata_manager.get_state_registry(
                REGISTRY_KEYS.BATCH_KEY
            )
            self._batch_mapping = batch_state_registry.categorical_mapping
            self._code_to_batch = dict(enumerate(self._batch_mapping))
        else:
            self._batch_mapping = None
            self._code_to_batch = None

    def _get_label_mapping_for_predictions(self):
        """Get label mapping excluding unlabeled category for predictions."""
        if self._label_mapping is None:
            return None

        # Filter out unlabeled category if it exists
        if (
            self.unlabeled_category_ is not None
            and self.unlabeled_category_ in self._label_mapping
        ):
            # Return all categories except the unlabeled one
            filtered_mapping = [
                cat for cat in self._label_mapping if cat != self.unlabeled_category_
            ]
            return np.array(filtered_mapping)
        else:
            return self._label_mapping

    def _get_code_to_label_for_predictions(self):
        """Get code-to-label mapping excluding unlabeled category for predictions."""
        filtered_mapping = self._get_label_mapping_for_predictions()
        if filtered_mapping is None:
            return None
        return dict(enumerate(filtered_mapping))

    @classmethod
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        batch_key: str | None = None,
        labels_key: str | None = None,
        unlabeled_category: str = "Unknown",
        size_factor_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        **kwargs,
    ) -> AnnDataManager | None:
        """Set up AnnData object for FADVI model.

        A mapping will be created between data fields used by FADVI and AnnData objects.
        None of the data in adata are modified. Only adds fields to uns.

        Parameters
        ----------
        AnnData object.
            Rows represent cells, columns represent features.
        layer
            If not `None`, uses this as the key in `adata.layers` for raw count data.
        batch_key
            key in `adata.obs` for batch information. Categories will automatically be converted into
            integer categories and saved to `adata.obs['_scvi_batch']`. If `None`, assigns the same batch
            to all the data.
        labels_key
            key in `adata.obs` for label information. Categories will automatically be converted into
            integer categories and saved to `adata.obs['_scvi_labels']`. If `None`, assigns the same label
            to all the data.
        unlabeled_category
            value in `adata.obs[labels_key]` that indicates unlabeled observations.
        size_factor_key
            key in `adata.obs` for size factor information. Instead of using library size as a size factor,
            the provided size factor column will be used as offset in the mean of the likelihood. Assumed
            to be on linear scale.
        categorical_covariate_keys
            keys in `adata.obs` that correspond to categorical data.
            These covariates can be added in addition to the batch covariate and are also treated as
            nuisance factors (i.e., the model tries to minimize their effects on the latent space). Thus,
            these should not be used for biologically-relevant factors that you do _not_ want to correct
            for.
        continuous_covariate_keys
            keys in `adata.obs` that correspond to continuous data.
            These covariates can be added in addition to the batch covariate and are also treated as
            nuisance factors (i.e., the model tries to minimize their effects on the latent space). Thus,
            these should not be used for biologically-relevant factors that you do _not_ want to correct
            for.

        Returns
        -------
        None. Adds the following fields:

            .uns['_scvi']
                `scvi` setup dictionary
            .obs['_scvi_labels']
                labels encoded as integers
            .obs['_scvi_batch']
                batch encoded as integers
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            LabelsWithUnlabeledObsField(
                REGISTRY_KEYS.LABELS_KEY, labels_key, unlabeled_category
            ),
            NumericalObsField(
                REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False
            ),
            CategoricalJointObsField(
                REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys
            ),
            NumericalJointObsField(
                REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys
            ),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    def get_latent_representation(
        self,
        adata: AnnData | None = None,
        indices: list[int] | None = None,
        give_mean: bool = True,
        mc_samples: int = 5000,
        batch_size: int | None = None,
        return_dist: bool = False,
        representation: Literal[
            "full", "b", "batch", "l", "label", "r", "residual", "lr", "label_residual"
        ] = "label",
    ) -> dict[str, torch.Tensor] | torch.Tensor:
        """Return the latent representation for each cell.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        give_mean
            Give mean of distribution or sample from it.
        mc_samples
            For distributions that have no closed-form mean (e.g. `LogNormal`), how many Monte Carlo
            samples to take for computing mean.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_dist
            If `True`, a mapping will be returned with key "dist" that contains distributional
            samples of the latent variables.
        representation
            Which latent representation to return:
            - "full": concatenated representation [z_b, z_l, z_r]
            - "b" or "batch": batch representation only
            - "l" or "label": label representation only
            - "r" or "residual": residual representation only
            - "lr" or "label_residual": concatenated label and residual representation

        Returns
        -------
        latent_representation
            Low-dimensional representation for each cell or dict of tensors if `return_dist` is True.
        """
        if self.is_trained_ is False:
            raise RuntimeError("Please train the model first.")

        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )

        latent = []
        for tensors in scdl:
            # Move tensors to model device
            device = next(self.module.parameters()).device
            tensors = {
                k: v.to(device) if isinstance(v, torch.Tensor) else v
                for k, v in tensors.items()
            }

            inference_inputs = self.module._get_inference_input(tensors)
            outputs = self.module.inference(**inference_inputs)

            if representation in ["b", "batch"]:
                z = outputs["z_b"]
                qz = outputs["qb"]
            elif representation in ["l", "label"]:
                z = outputs["z_l"]
                qz = outputs["ql"]
            elif representation in ["r", "residual"]:
                z = outputs["z_r"]
                qz = outputs["qr"]
            elif representation in ["lr", "label_residual"]:
                z_l = outputs["z_l"]
                z_r = outputs["z_r"]
                z = torch.cat([z_l, z_r], dim=-1)
                qz = None
            else:  # "full"
                z_b = outputs["z_b"]
                z_l = outputs["z_l"]
                z_r = outputs["z_r"]
                z = torch.cat([z_b, z_l, z_r], dim=-1)
                # For full representation, we don't return qz for now
                qz = None

            if not give_mean and qz is not None:
                z = qz.sample()

            latent += [z.cpu()]

        if return_dist and qz is not None:
            return {"mean": torch.cat(latent).detach().numpy(), "dist": qz}

        return torch.cat(latent).detach().numpy()

    def get_normalized_expression(
        self,
        adata: AnnData | None = None,
        indices: list[int] | None = None,
        transform_batch: str | int | None = None,
        gene_list: list[str] | None = None,
        library_size: float | Literal["latent"] = 1,
        n_samples: int = 1,
        n_samples_overall: int | None = None,
        batch_size: int | None = None,
        return_mean: bool = True,
        return_numpy: bool | None = None,
    ):
        """Return the normalized (decoded) gene expression.

        This is denoted as :math:`\\rho_n` in the scVI paper.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        transform_batch
            Batch to condition on.
            If transform_batch is:

            - None, then real observed batch is used.
            - int, then batch transform_batch is used.
        gene_list
            Return frequencies of expression for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
        library_size
            Scale the expression frequencies to a common library size.
            This allows gene expression levels to be interpreted on a common scale of relevant
            magnitude. If set to `"latent"`, use the latent library size.
        n_samples
            Number of posterior samples to use for estimation.
        n_samples_overall
            Number of posterior samples to use for estimation. Overrides `n_samples`.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a :class:`~numpy.ndarray` instead of a :class:`~pandas.DataFrame`. DataFrame includes
            gene names as columns. If either `n_samples=1` or `return_mean=True`, defaults to `False`.
            Otherwise, it defaults to `True`.

        Returns
        -------
        normalized_expression
            If `n_samples` > 1 and `return_mean` is False, then the shape is `(samples, cells, genes)`.
            Otherwise, shape is `(cells, genes)`. In this case, return type is
            :class:`~pandas.DataFrame` unless `return_numpy` is True.
        """
        return super().get_normalized_expression(
            adata=adata,
            indices=indices,
            transform_batch=transform_batch,
            gene_list=gene_list,
            library_size=library_size,
            n_samples=n_samples,
            n_samples_overall=n_samples_overall,
            batch_size=batch_size,
            return_mean=return_mean,
            return_numpy=return_numpy,
        )

    @torch.inference_mode()
    def predict(
        self,
        adata: AnnData | None = None,
        indices: list[int] | None = None,
        prediction_mode: Literal["b", "batch", "l", "label"] = "label",
        soft: bool = False,
        batch_size: int | None = None,
        return_numpy: bool = True,
    ):
        """Predict batch or label categories using the supervised classification heads.

        This method uses the trained encoders and classification heads to predict
        either batch categories (using batch latent b) or label categories (using
        label latent l).

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        prediction_mode
            What to predict:
            - "b" or "batch": Predict batch categories using the batch latent (b) and batch classifier
            - "l" or "label": Predict label categories using the label latent (l) and label classifier
        soft
            If `True`, return class probabilities. If `False`, return class predictions.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_numpy
            Return a :class:`~numpy.ndarray` instead of a :class:`~torch.Tensor`.

        Returns
        -------
        predictions
            If `soft=True`, returns class probabilities with shape `(n_cells, n_classes)`.
            If `soft=False`, returns class predictions with shape `(n_cells,)`.
            If `return_numpy=True`, returns numpy array, otherwise torch tensor.
        """
        if self.module.training:
            self.module.eval()

        adata = self._validate_anndata(adata)
        if indices is None:
            indices = np.arange(adata.n_obs)

        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )

        predictions = []
        for tensors in scdl:
            # Move tensors to model device
            device = next(self.module.parameters()).device
            tensors = {
                k: v.to(device) if isinstance(v, torch.Tensor) else v
                for k, v in tensors.items()
            }

            x = tensors[REGISTRY_KEYS.X_KEY]
            batch_index = tensors.get(REGISTRY_KEYS.BATCH_KEY, None)
            cat_covs = tensors.get(REGISTRY_KEYS.CAT_COVS_KEY, None)
            cont_covs = tensors.get(REGISTRY_KEYS.CONT_COVS_KEY, None)

            # Get the appropriate latent representation and predictions
            if prediction_mode in ["b", "batch"]:
                # Use inference method to get batch latent (device-aware)
                inference_outputs = self.module.inference(
                    x, batch_index=batch_index, cat_covs=cat_covs, cont_covs=cont_covs
                )
                z_b = inference_outputs["z_b"]

                # Get batch predictions from batch classifier
                batch_logits = self.module.head_batch(z_b)

                if soft:
                    batch_pred = torch.softmax(batch_logits, dim=1)
                else:
                    batch_pred = torch.argmax(batch_logits, dim=1)

                predictions.append(batch_pred.cpu())

            elif prediction_mode in ["l", "label"]:
                # Use inference method to get label latent (device-aware)
                inference_outputs = self.module.inference(
                    x, batch_index=batch_index, cat_covs=cat_covs, cont_covs=cont_covs
                )
                z_l = inference_outputs["z_l"]

                # Get label predictions from label classifier
                label_logits = self.module.head_label(z_l)

                if soft:
                    label_pred = torch.softmax(label_logits, dim=1)
                    # Slice to exclude unlabeled category probabilities if present
                    filtered_mapping = self._get_label_mapping_for_predictions()
                    if (
                        filtered_mapping is not None
                        and len(filtered_mapping) < label_pred.shape[1]
                    ):
                        # Remove the unlabeled category probability (assumed to be last)
                        label_pred = label_pred[:, :-1]
                else:
                    label_pred = torch.argmax(label_logits, dim=1)

                predictions.append(label_pred.cpu())
            else:
                raise ValueError(
                    f"prediction_mode must be 'b', 'batch', 'l', or 'label', got {prediction_mode}"
                )

        predictions = torch.cat(predictions, dim=0)

        # Convert to numpy if tensor
        if isinstance(predictions, torch.Tensor):
            predictions = predictions.numpy()

        # Map numerical predictions back to original categorical labels
        if not soft:
            # Hard predictions - map indices to original labels
            if prediction_mode in ["b", "batch"] and self._code_to_batch is not None:
                predictions = np.array(
                    [self._code_to_batch[int(p)] for p in predictions]
                )
            elif prediction_mode in ["l", "label"]:
                # Use filtered mapping that excludes unlabeled category
                code_to_label_filtered = self._get_code_to_label_for_predictions()
                if code_to_label_filtered is not None:
                    # Handle case where classifier might predict unlabeled category index
                    # Find the index of unlabeled category in the original mapping
                    unlabeled_idx = None
                    if self.unlabeled_category_ is not None:
                        unlabeled_indices = np.where(
                            self._label_mapping == self.unlabeled_category_
                        )[0]
                        if len(unlabeled_indices) > 0:
                            unlabeled_idx = unlabeled_indices[0]

                    mapped_predictions = []
                    for p in predictions:
                        p_int = int(p)
                        if p_int in code_to_label_filtered:
                            # Valid prediction, use filtered mapping
                            mapped_predictions.append(code_to_label_filtered[p_int])
                        elif unlabeled_idx is not None and p_int == unlabeled_idx:
                            # Classifier predicted unlabeled category - this shouldn't happen in a well-trained model
                            # For now, map to the most frequent class or handle gracefully
                            # We'll map to index 0 (first valid category) as a fallback
                            mapped_predictions.append(code_to_label_filtered[0])
                        else:
                            # Unexpected index, fallback to first valid category
                            mapped_predictions.append(code_to_label_filtered[0])

                    predictions = np.array(mapped_predictions)
        else:
            # Soft predictions - create DataFrame with original label names as columns
            if prediction_mode in ["b", "batch"] and self._batch_mapping is not None:
                n_categories = len(predictions[0])
                predictions = pd.DataFrame(
                    predictions,
                    columns=self._batch_mapping[:n_categories],
                    index=adata.obs_names[indices],
                )
            elif prediction_mode in ["l", "label"]:
                # Use filtered mapping that excludes unlabeled category
                label_mapping_filtered = self._get_label_mapping_for_predictions()
                if label_mapping_filtered is not None:
                    n_categories = len(predictions[0])
                    predictions = pd.DataFrame(
                        predictions,
                        columns=label_mapping_filtered[:n_categories],
                        index=adata.obs_names[indices],
                    )

        # Return numpy array or pandas DataFrame as requested
        if not return_numpy and not soft:
            return predictions  # Keep as array for hard predictions
        elif soft:
            return predictions  # DataFrame for soft predictions
        else:
            return predictions  # numpy array for hard predictions
