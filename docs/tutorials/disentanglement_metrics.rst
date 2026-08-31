Evaluating Disentanglement
========================================

FADVI factorizes the latent space into three subspaces. This tutorial shows how to
verify it quantitatively: whether batch information really is confined to
``z_b``, whether biological identity really is confined to ``z_l``, and whether
the residual subspace ``z_r`` carries its own signal rather than absorbing the
other two.

All of the functions below are importable directly from ``fadvi``, so they can be
applied to latent representations from any source -- not only to a fitted FADVI
model.

Why measure it
----------------------------------------

A good integration score is not evidence of disentanglement. A model can score
well while still leaving batch structure in the biological subspace, and feature
attribution on a subspace that has absorbed the wrong signal will produce
misleading gene lists. The metrics here test the factorization directly.

Three questions, three kinds of metric:

* **Predictability** -- can each factor be recovered from its own subspace, and
  *not* from the others?
* **Independence** -- how much linear structure do the subspaces share?
* **Residual interpretability** -- does ``z_r`` behave as a signal sink?

Quick evaluation
----------------------------------------

The convenience methods on a fitted model extract the three subspaces for you:

.. code-block:: python

   from fadvi import FADVI

   FADVI.setup_anndata(
       adata, batch_key="batch", labels_key="cell_type", layer="counts"
   )
   model = FADVI(adata)
   model.train(max_epochs=30)

   results = model.evaluate_disentanglement(adata)

   print(results["batch_predictability"]["accuracy"])       # z_b -> batch
   print(results["biological_predictability"]["accuracy"])  # z_l -> cell type
   print(results["cross_correlations"])                     # between subspaces
   print(results["overall_quality"]["overall_score"])

Focused analyses are available as separate methods:

.. code-block:: python

   model.analyze_subspace_specificity(adata)      # per-subspace predictability
   model.compute_orthogonality_scores(adata)      # CCA-based independence
   model.analyze_residual_interpretability(adata) # is z_r a signal sink?

Using the metrics on any representation
----------------------------------------

The same functions work on plain arrays, which is what you need to compare FADVI
against another method or to score embeddings loaded from disk:

.. code-block:: python

   import numpy as np
   from fadvi import (
       compute_cross_correlation,
       compute_orthogonality_score,
       evaluate_subspace_predictability,
       evaluate_disentanglement_quality,
   )

   z_b = model.get_latent_representation(adata, representation="b")
   z_l = model.get_latent_representation(adata, representation="l")
   z_r = model.get_latent_representation(adata, representation="r")

   # 1. Can each factor be recovered from its own subspace?
   batch_from_b = evaluate_subspace_predictability(z_b, adata.obs["batch"])
   type_from_l  = evaluate_subspace_predictability(z_l, adata.obs["cell_type"])
   print(batch_from_b["accuracy"], type_from_l["accuracy"])

   # 2. ...and NOT from the others? (this is the informative direction)
   batch_from_l = evaluate_subspace_predictability(z_l, adata.obs["batch"])
   type_from_b  = evaluate_subspace_predictability(z_b, adata.obs["cell_type"])

   # 3. How much structure do the subspaces share?
   mean_corr, corr_matrix = compute_cross_correlation(z_b, z_l)
   orth = compute_orthogonality_score(z_b, z_l)      # method="cca" by default

Interpreting the numbers
----------------------------------------

Read accuracy against the majority class, not against chance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the single most common way to misread these metrics. Class frequencies
in real atlases are far from uniform, so the reference point is the
**majority-class rate** -- the accuracy of always predicting the commonest class
-- and *not* ``1 / n_classes``.

Concretely: in a dataset where 92% of cells come from one batch, a classifier
that recovers batch from ``z_r`` at 0.92 accuracy has learned **nothing**; it is
reproducing the class prior. Against a ``1 / n_classes`` reference of 0.04 the
same number would look like severe leakage.

.. code-block:: python

   import numpy as np

   def majority_rate(labels):
       _, counts = np.unique(np.asarray(labels), return_counts=True)
       return counts.max() / counts.sum()

   leak = evaluate_subspace_predictability(z_r, adata.obs["batch"])
   baseline = majority_rate(adata.obs["batch"])
   print(f"z_r -> batch: {leak['accuracy']:.3f}  (baseline {baseline:.3f}, "
         f"excess {leak['accuracy'] - baseline:+.3f})")

Only the **excess over the baseline** indicates real information. Predictability
is scored with cross-validation, so it is not inflated by overfitting.

Cross-correlation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``compute_cross_correlation`` returns the mean absolute correlation over *all*
pairs of dimensions, plus the full matrix. Lower is better. Note that because
the mean is taken over every pair, even two identical inputs score well below
1.0 -- only the matching diagonal pairs correlate strongly. Compare values
across subspace pairs rather than against an absolute cut-off, and inspect
``corr_matrix`` when you want to know *which* dimensions are entangled.

Orthogonality
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``compute_orthogonality_score`` returns a canonical-correlation based score in
``[0, 1]``, where higher means the two subspaces share less linear structure. It
is a continuous quantity and there is no community-standard threshold that
separates "orthogonal" from "entangled", so report the value rather than a
pass/fail verdict. In practice this score varies more across datasets than
cross-correlation does, so the two are worth reporting together.

Checking the residual subspace
----------------------------------------

The residual subspace is only weakly constrained, so it is the most likely place
for leaked signal to accumulate. That matters if you intend to run feature
attribution, since attributing to a subspace that has absorbed the wrong signal
gives a meaningless gene list.

.. code-block:: python

   from fadvi import analyze_residual_interpretability

   res = analyze_residual_interpretability(
       z_r, z_l, z_b, adata,
       batch_key="batch",        # explicit keys are safer than the fallbacks
       labels_key="cell_type",
   )

   print(res["batch_leakage_accuracy"], res["label_leakage_accuracy"])
   print(res["residual_variance_fraction"])   # z_r's share of latent variance
   print(res["is_signal_sink"])               # overall verdict

.. note::

   Pass ``batch_key`` and ``labels_key`` explicitly whenever you know them. If
   omitted, the function falls back to searching for ``'batch'`` then
   ``'method'``, and ``'cell_type'`` then ``'annotation'`` then ``'labels'``.
   The fallback is convenient but silent: if none of the candidates is present,
   the corresponding leakage test is skipped and its key is simply absent from
   the returned dictionary. The keys ``batch_key_used`` and ``labels_key_used``
   record which columns were actually used, so you can confirm.

Comparing several methods
----------------------------------------

``compare_disentanglement_methods`` scores multiple representations in one call
and returns a tidy DataFrame. Any method that yields separate batch and
biological embeddings can be included:

.. code-block:: python

   from fadvi import compare_disentanglement_methods

   comparison = compare_disentanglement_methods(
       {
           "FADVI": {"z_b": z_b, "z_l": z_l, "z_r": z_r},
           "other_method": {"batch": other_z_b, "label": other_z_l},
       },
       adata,
       batch_key="batch",
       label_key="cell_type",
   )
   print(comparison)

Both naming conventions (``z_b``/``z_l``/``z_r`` and
``batch``/``label``/``residual``) are accepted. If no residual embedding is
supplied, a zero array is substituted so that the remaining metrics can still be
computed; residual-specific results are then not meaningful for that method.

Practical notes
----------------------------------------

* **Cost.** Predictability fits a cross-validated classifier per call, and
  ``analyze_residual_interpretability`` correlates every residual dimension
  against every gene. On a large atlas, subsample cells (a few tens of thousands
  is usually ample) before evaluating.
* **Which representation.** These metrics describe the subspaces themselves. Use
  ``representation="l"`` for integration, clustering and visualisation.
* **Report values, not verdicts.** The boolean fields
  (``is_signal_sink``, ``meets_quality_criteria``) apply fixed internal
  thresholds and exist for convenience. The underlying numbers, read against the
  majority-class baseline, are what should go into a figure or a manuscript.

See also
----------------------------------------

* :doc:`../implementation` -- loss terms, hyperparameters and their sensitivity
* :doc:`advanced_usage` -- feature attribution on the latent subspaces
* :doc:`../api/metrics` -- full API reference for these functions
