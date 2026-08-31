Disentanglement Metrics
=======================

Quantify how cleanly the three latent subspaces are separated: whether batch
information is confined to ``z_b``, biological identity to ``z_l``, and whether
the unsupervised residual subspace ``z_r`` carries its own distinct signal
rather than absorbing the other two.

The most convenient entry points are the methods on
:class:`~fadvi.FADVI` itself:

.. code-block:: python

   results = model.evaluate_disentanglement(adata)
   print(results["batch_predictability"]["accuracy"])       # z_b -> batch
   print(results["biological_predictability"]["accuracy"])  # z_l -> cell type
   print(results["cross_correlations"])                     # between subspaces

   # Focused analyses
   model.analyze_subspace_specificity(adata)
   model.compute_orthogonality_scores(adata)
   model.analyze_residual_interpretability(adata)

Interpreting the numbers:

* **Predictability** — accuracy of recovering a factor from a subspace.
  A well-disentangled model predicts batch from ``z_b`` and cell type from
  ``z_l`` with high accuracy, and each factor poorly from the other subspace.
  Accuracy is computed with cross-validation, so it is not inflated by
  overfitting.
* **Cross-correlation** — mean absolute correlation between two subspaces.
  Lower is better; values below ~0.3 indicate good separation. Note the score
  averages over *all* dimension pairs, so even identical inputs score well
  below 1.
* **Orthogonality** — canonical-correlation based score in ``[0, 1]``, where
  higher means the two subspaces share less linear structure.

.. automodule:: fadvi._disentanglement_metrics
   :members:
   :undoc-members:
   :show-inheritance:
