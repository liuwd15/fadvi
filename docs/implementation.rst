Implementation Details
======================

This page documents how FADVI is trained, why it is stable, and how to choose
its hyperparameters. It is intended both as practical guidance and as a
reproducibility reference.

.. contents::
   :local:
   :depth: 2


Model structure
---------------

FADVI encodes each cell into three latent subspaces:

* ``z_b`` — batch-associated variation,
* ``z_l`` — biological (label-associated) variation,
* ``z_r`` — residual variation not explained by either.

The decoder reconstructs expression from all three. Separation is enforced by
four mechanisms working together, all optimised jointly in a single objective:

.. code-block:: text

   total loss = ELBO
              + lambda_b * CE(batch  | z_b)     # supervised: batch head
              + lambda_l * CE(label  | z_l)     # supervised: label head
              + adversarial terms via gradient reversal
              + gamma * cross-covariance penalty

**Supervised heads** pull the intended signal into its subspace.
**Adversarial heads** push the *unintended* signal out: they attempt to predict
the label from ``z_b``, the batch from ``z_l``, and both from ``z_r``, with a
gradient-reversal layer in between so the encoder learns to make those
predictions *fail*. The **cross-covariance penalty** is a non-adversarial term
that directly decorrelates the subspaces, providing a second, convex route to
the same goal.


Gradient reversal and training stability
----------------------------------------

Adversarial objectives in VAEs have a reputation for instability. That
reputation comes from *alternating* minimax optimization, in which a generator
and a discriminator are updated in separate steps — the setting that produces
oscillation and mode collapse in GAN training.

**FADVI does not use alternating optimization.** It uses a gradient-reversal
layer (GRL; Ganin & Lempitsky, 2015):

.. code-block:: python

   class GradientReversalFunction(torch.autograd.Function):
       @staticmethod
       def forward(ctx, x, alpha):
           ctx.alpha = alpha
           return x.view_as(x)          # identity on the forward pass

       @staticmethod
       def backward(ctx, grad_output):
           return grad_output.neg() * ctx.alpha, None   # sign-flipped gradient

The forward pass is the identity; the backward pass negates and rescales the
gradient. The encoder and the auxiliary classifiers are therefore trained
**jointly in a single backward pass under ordinary SGD** — the "min" and the
"max" happen simultaneously rather than as a two-player game.

Three further design choices reinforce stability:

1. **The adversarial weight is small and fixed.** ``alpha_* = 1`` by default and
   is never scheduled, so the adversary cannot dominate an objective anchored by
   reconstruction and supervision (``lambda_* = 50``).
2. **Disentanglement is redundantly enforced.** The cross-covariance penalty
   (``gamma``) targets the same goal without any adversarial dynamics, so
   separation does not hinge on the adversary alone.
3. **Empirically verified.** Sweeping the gradient-reversal strength from 0 to
   50× the default on a 330k-cell atlas, every run converged smoothly, the
   adversarial loss stayed bounded, no run produced a NaN or Inf, and the
   largest single-epoch increase in training loss was at most 0.13% of its
   initial value.

Monitoring training
~~~~~~~~~~~~~~~~~~~

Every loss component is logged separately and is available in
``model.history``:

.. code-block:: python

   model.train(max_epochs=30)
   model.history.keys()
   # 'train_loss', 'elbo_loss_train', 'sup_loss_train', 'adv_loss_train',
   # 'xcov_loss_train', 'adv_loss_bl_train', 'adv_loss_lb_train',
   # 'adv_loss_rb_train', 'adv_loss_rl_train', ...

Healthy training shows the total loss decreasing smoothly and the adversarial
loss remaining **bounded** — it should *not* collapse toward zero (the adversary
has been defeated outright, so no pressure remains) nor grow without limit.

.. note::

   scvi-tools only runs the validation loop if you ask it to. Pass
   ``check_val_every_n_epoch=1`` (or enable early stopping) if you want
   ``*_validation`` metrics recorded; otherwise no validation curve is logged.


Hyperparameters
---------------

.. list-table::
   :header-rows: 1
   :widths: 16 12 72

   * - Parameter
     - Default
     - Role
   * - ``beta``
     - 1.0
     - Weight on the KL term of the ELBO. Higher values regularise the latent
       space more strongly, trading biological conservation for batch mixing.
   * - ``lambda_b``
     - 50
     - Batch-classification weight; drives batch signal into ``z_b``.
   * - ``lambda_l``
     - 50
     - Label-classification weight; drives biological signal into ``z_l``.
       Set to 0 for fully unsupervised training.
   * - ``alpha_bl``, ``alpha_lb``, ``alpha_rb``, ``alpha_rl``
     - 1.0
     - Gradient-reversal strengths for the four adversarial heads
       (label-from-``z_b``, batch-from-``z_l``, batch-from-``z_r``,
       label-from-``z_r``).
   * - ``gamma``
     - 1.0
     - Cross-covariance (decorrelation) penalty between subspaces.
   * - ``n_latent_b`` / ``n_latent_l`` / ``n_latent_r``
     - 30 / 30 / 10
     - Dimensionality of each subspace. ``z_r`` is deliberately smaller, since
       it should capture only leftover variation.

Sensitivity
~~~~~~~~~~~

A one-at-a-time sweep of all five weights across six single-cell atlases found
the total integration score varied by at most **0.015–0.041** (≤ ~5% of the
score) over wide ranges — ``beta`` ∈ [0.1, 5], ``lambda_*`` ∈ [10, 200],
``alpha``/``gamma`` ∈ [0, 10]. **The defaults are near-optimal and no
per-dataset tuning is required.**

The trends that do exist are smooth and interpretable:

* increasing ``beta`` shifts the balance from biological conservation toward
  batch correction;
* increasing ``lambda_l`` improves biological conservation;
* ``alpha`` and ``gamma`` have the smallest effect of all — turning the
  adversary off entirely, or up tenfold, changes little.

If you do tune, adjust ``beta`` first, then ``lambda_l``.


Training procedure and reproducibility
--------------------------------------

.. code-block:: python

   import scanpy as sc
   from fadvi import FADVI
   import scvi

   scvi.settings.seed = 0                  # reproducible initialisation

   FADVI.setup_anndata(
       adata,
       batch_key="batch",
       labels_key="cell_type",
       unlabeled_category="Unknown",       # semi-supervised: mask unknown cells
       layer="counts",                     # RAW counts, not normalised values
   )
   model = FADVI(adata, n_latent_b=30, n_latent_l=30, n_latent_r=10, n_layers=2)
   model.train(max_epochs=30, batch_size=256)

   adata.obsm["X_fadvi"] = model.get_latent_representation(representation="l")

Practical notes
~~~~~~~~~~~~~~~

* **Always train on raw counts.** The likelihood (ZINB by default) models
  counts; passing normalised values will degrade results.
* **Feature selection.** Around 2,000 highly variable genes is a good default.
* **Epochs.** 30 epochs is adequate for most datasets. On a 330k-cell atlas the
  validation ELBO plateaus at roughly 110 epochs if you train to full
  convergence.
* **Batch size.** 256 works well; larger batches speed up training on big
  datasets at some cost in gradient noise.
* **A single-sample final minibatch will crash batch normalization.** If your
  cell count leaves a remainder of one, pass
  ``datasplitter_kwargs={"drop_last": True}``.
* **Which representation to use.** Use ``representation="l"`` for integration,
  clustering and UMAP. Use ``"lr"`` in fully unsupervised mode.

Expected runtime and memory
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Measured on a single NVIDIA L40S, 2,000 genes, batch size 256:

.. list-table::
   :header-rows: 1
   :widths: 20 26 26 28

   * - Cells
     - Seconds / epoch
     - Minutes to convergence
     - Peak GPU / host memory
   * - 10,000
     - 1.1
     - 1.8
     - 0.10 GB / 3.8 GB
   * - 50,000
     - 5.4
     - 9.4
     - 0.10 GB / 4.7 GB
   * - 100,000
     - 10.8
     - 19.0
     - 0.10 GB / 6.7 GB
   * - 300,000
     - 32.5
     - 57.2
     - 0.10 GB / 12.7 GB

Runtime scales **linearly** with cell number, and GPU memory is independent of
dataset size because training is minibatched.


See :doc:`tutorials/spatial_and_single_cell` for recommended practice with
high-resolution spatial data.
