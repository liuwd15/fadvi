"""FADVI: Factor Disentanglement Variational Inference.

A VAE framework that factorises single-cell and spatial omics data into
batch-associated (``z_b``), biological/label-associated (``z_l``) and residual
(``z_r``) latent subspaces.

The model classes :class:`FADVI` and :class:`FADVAE` are the main entry points.
The disentanglement metrics are also exported as standalone functions, so they
can be applied to latent representations from any source -- including other
integration methods -- and not only to a fitted FADVI model. The equivalent
methods on :class:`FADVI` (``evaluate_disentanglement``,
``analyze_subspace_specificity``, ``compute_orthogonality_scores`` and
``analyze_residual_interpretability``) are thin wrappers that extract the three
subspaces for you.
"""

from ._fadvae import FADVAE
from ._fadvi import FADVI
from ._disentanglement_metrics import (
    compute_cross_correlation,
    compute_mutual_information,
    compute_orthogonality_score,
    evaluate_subspace_predictability,
    analyze_residual_interpretability,
    evaluate_disentanglement_quality,
    compare_disentanglement_methods,
)

__all__ = [
    # models
    "FADVI",
    "FADVAE",
    # disentanglement metrics
    "compute_cross_correlation",
    "compute_mutual_information",
    "compute_orthogonality_score",
    "evaluate_subspace_predictability",
    "analyze_residual_interpretability",
    "evaluate_disentanglement_quality",
    "compare_disentanglement_methods",
]
