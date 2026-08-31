"""
Disentanglement validation metrics for FADVI latent subspaces.

This module provides quantitative methods to validate that:
1. Biological signals concentrate in z_l (label latent space)
2. Batch effects concentrate in z_b (batch latent space)
3. Residual subspace z_r captures meaningful orthogonal variation (not a signal sink)
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple, Union
import warnings

import numpy as np
import pandas as pd
import torch
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import LabelEncoder
from sklearn.decomposition import PCA
from sklearn.feature_selection import mutual_info_regression, mutual_info_classif
from sklearn.cross_decomposition import CCA

import scanpy as sc
import anndata

warnings.filterwarnings('ignore')


def compute_cross_correlation(
    z1: np.ndarray,
    z2: np.ndarray,
    method: str = "pearson"
) -> Tuple[float, np.ndarray]:
    """
    Compute cross-correlation between two latent subspaces.

    Lower correlation indicates better disentanglement (target: < 0.3).

    Parameters
    ----------
    z1 : np.ndarray
        First latent representation (n_cells, n_latent1)
    z2 : np.ndarray
        Second latent representation (n_cells, n_latent2)
    method : str
        Correlation method ('pearson' or 'spearman')

    Returns
    -------
    mean_corr : float
        Mean absolute correlation across all factor pairs
    corr_matrix : np.ndarray
        Full correlation matrix (n_latent1, n_latent2)
    """
    n_cells1, n_latent1 = z1.shape
    n_cells2, n_latent2 = z2.shape

    assert n_cells1 == n_cells2, "Latent representations must have same number of cells"

    corr_func = pearsonr if method == "pearson" else spearmanr
    corr_matrix = np.zeros((n_latent1, n_latent2))

    for i in range(n_latent1):
        for j in range(n_latent2):
            corr, _ = corr_func(z1[:, i], z2[:, j])
            corr_matrix[i, j] = corr if not np.isnan(corr) else 0.0

    mean_corr = np.mean(np.abs(corr_matrix))
    return mean_corr, corr_matrix


def compute_mutual_information(
    z_continuous: np.ndarray,
    z_discrete: np.ndarray,
    discrete_threshold: float = 0.1
) -> float:
    """
    Compute mutual information between continuous and discrete latent factors.

    Parameters
    ----------
    z_continuous : np.ndarray
        Continuous latent factors (n_cells, n_factors)
    z_discrete : np.ndarray
        Discrete latent factors or labels (n_cells,)
    discrete_threshold : float
        Threshold for discretizing continuous variables if needed

    Returns
    -------
    float
        Mean mutual information across all factor combinations
    """
    if z_continuous.ndim == 1:
        z_continuous = z_continuous.reshape(-1, 1)

    n_factors = z_continuous.shape[1]
    mi_scores = []

    for i in range(n_factors):
        # Use regression MI for continuous vs discrete
        mi = mutual_info_regression(
            z_continuous[:, i].reshape(-1, 1),
            z_discrete,
            random_state=42
        )[0]
        mi_scores.append(mi)

    return np.mean(mi_scores)


def evaluate_subspace_predictability(
    latent_repr: np.ndarray,
    target_labels: np.ndarray,
    task_type: str = "classification",
    cv_folds: int = 5
) -> Dict[str, float]:
    """
    Evaluate how well a latent subspace can predict target labels.

    High predictability indicates the subspace captures information
    relevant to the target (batch effects, cell types, etc.).

    Uses cross-validation for primary accuracy metric (unbiased) and
    train/test split for additional metrics to avoid overfitting bias.

    Parameters
    ----------
    latent_repr : np.ndarray
        Latent representation (n_cells, n_latent)
    target_labels : np.ndarray
        Target labels to predict (n_cells,)
    task_type : str
        'classification' or 'regression'
    cv_folds : int
        Number of cross-validation folds

    Returns
    -------
    Dict[str, float]
        Dictionary with predictability metrics:
        - accuracy: Cross-validated accuracy/R2 (primary unbiased metric)
        - f1_score: Test set F1 score (classification only)
        - auc: Test set AUC score (classification only)
        - r2_score: Test set R2 score (regression only)
        - n_classes: Number of unique classes/targets
    """
    # Handle string labels
    if target_labels.dtype == object or isinstance(target_labels[0], str):
        le = LabelEncoder()
        target_encoded = le.fit_transform(target_labels)
    else:
        target_encoded = target_labels

    # Skip if too few unique labels
    n_unique = len(np.unique(target_encoded))
    if n_unique < 2:
        return {"accuracy": 0.0, "f1_score": 0.0, "auc": 0.0}

    if task_type == "classification":
        from sklearn.model_selection import train_test_split

        if n_unique > 10:  # Many classes, use simpler model
            model = LogisticRegression(random_state=42, max_iter=1000)
        else:
            model = RandomForestClassifier(n_estimators=50, random_state=42, n_jobs=2)

        # Cross-validation scores (unbiased primary metric)
        cv_scores = cross_val_score(model, latent_repr, target_encoded,
                                  cv=min(cv_folds, n_unique), scoring='accuracy')
        accuracy = cv_scores.mean()

        # Train/test split for additional metrics (avoid overfitting bias)
        try:
            X_train, X_test, y_train, y_test = train_test_split(
                latent_repr, target_encoded,
                test_size=0.3,
                random_state=42,
                stratify=target_encoded if n_unique > 1 else None
            )

            # Fit on training data only
            model.fit(X_train, y_train)
            pred_proba = model.predict_proba(X_test)
            predictions = model.predict(X_test)

            f1 = f1_score(y_test, predictions, average='weighted')

            # AUC for binary or use macro for multiclass
            if n_unique == 2:
                auc = roc_auc_score(y_test, pred_proba[:, 1])
            else:
                auc = roc_auc_score(y_test, pred_proba,
                                   multi_class='ovr', average='macro')
        except (ValueError, Exception):
            # Fallback to CV-only metrics if train/test split fails
            f1 = accuracy  # Use CV accuracy as proxy
            auc = accuracy

        return {
            "accuracy": accuracy,  # Primary metric (CV-based, unbiased)
            "f1_score": f1,       # Test set F1 (unbiased)
            "auc": auc,           # Test set AUC (unbiased)
            "n_classes": n_unique
        }

    else:  # regression
        from sklearn.ensemble import RandomForestRegressor
        from sklearn.metrics import r2_score
        from sklearn.model_selection import train_test_split

        model = RandomForestRegressor(n_estimators=50, random_state=42, n_jobs=2)

        # Cross-validation scores (unbiased primary metric)
        cv_scores = cross_val_score(model, latent_repr, target_encoded,
                                  cv=cv_folds, scoring='r2')
        r2_cv = cv_scores.mean()

        # Train/test split for additional validation
        try:
            X_train, X_test, y_train, y_test = train_test_split(
                latent_repr, target_encoded,
                test_size=0.3,
                random_state=42
            )

            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            r2_test = r2_score(y_test, y_pred)
        except Exception:
            r2_test = r2_cv  # Fallback to CV score

        return {"r2_score": r2_test, "accuracy": r2_cv}  # Use CV r2 as primary accuracy


def compute_orthogonality_score(
    z1: np.ndarray,
    z2: np.ndarray,
    method: str = "cca"
) -> Dict[str, float]:
    """
    Compute orthogonality between two latent subspaces.

    Parameters
    ----------
    z1 : np.ndarray
        First latent subspace (n_cells, n_latent1)
    z2 : np.ndarray
        Second latent subspace (n_cells, n_latent2)
    method : str
        Method to use ('cca' for canonical correlation, 'cosine' for cosine similarity)

    Returns
    -------
    Dict[str, float]
        Orthogonality metrics
    """
    if method == "cca":
        # Canonical Correlation Analysis
        n_components = min(z1.shape[1], z2.shape[1], z1.shape[0] // 10)
        n_components = max(1, n_components)

        try:
            cca = CCA(n_components=n_components)
            cca.fit(z1, z2)

            # Transform and compute correlations
            z1_c, z2_c = cca.transform(z1, z2)

            correlations = []
            for i in range(n_components):
                corr, _ = pearsonr(z1_c[:, i], z2_c[:, i])
                if not np.isnan(corr):
                    correlations.append(abs(corr))

            if correlations:
                max_canonical_corr = max(correlations)
                mean_canonical_corr = np.mean(correlations)
            else:
                max_canonical_corr = mean_canonical_corr = 0.0

            # Orthogonality score (1 - correlation, higher is more orthogonal)
            orthogonality_score = 1 - max_canonical_corr

        except Exception as e:
            print(f"CCA failed: {e}")
            orthogonality_score = max_canonical_corr = mean_canonical_corr = 0.0

        return {
            "orthogonality_score": orthogonality_score,
            "max_canonical_correlation": max_canonical_corr,
            "mean_canonical_correlation": mean_canonical_corr,
            "n_components": n_components
        }

    elif method == "cosine":
        # Cosine similarity between subspace centroids
        centroid1 = np.mean(z1, axis=0)
        centroid2 = np.mean(z2, axis=0)

        # Normalize
        centroid1 = centroid1 / (np.linalg.norm(centroid1) + 1e-8)
        centroid2 = centroid2 / (np.linalg.norm(centroid2) + 1e-8)

        # Pad to same dimension
        max_dim = max(len(centroid1), len(centroid2))
        centroid1 = np.pad(centroid1, (0, max_dim - len(centroid1)))
        centroid2 = np.pad(centroid2, (0, max_dim - len(centroid2)))

        cosine_sim = np.dot(centroid1, centroid2)
        orthogonality_score = 1 - abs(cosine_sim)

        return {
            "orthogonality_score": orthogonality_score,
            "cosine_similarity": cosine_sim
        }


def analyze_residual_interpretability(
    z_residual: np.ndarray,
    z_biological: np.ndarray,
    z_batch: np.ndarray,
    adata: anndata.AnnData,
    n_top_genes: int = 100,
    batch_key: Optional[str] = None,
    labels_key: Optional[str] = None,
) -> Dict[str, Union[float, List[str]]]:
    """
    Analyze whether the residual subspace captures meaningful variation
    or acts as a signal sink.

    Parameters
    ----------
    z_residual : np.ndarray
        Residual latent factors (n_cells, n_latent_r)
    z_biological : np.ndarray
        Biological latent factors (n_cells, n_latent_l)
    z_batch : np.ndarray
        Batch latent factors (n_cells, n_latent_b)
    adata : anndata.AnnData
        Expression data for gene correlation analysis
    n_top_genes : int
        Number of top correlated genes to identify per residual factor
    batch_key : str, optional
        Column in ``adata.obs`` holding the batch variable. If omitted, falls
        back to ``'batch'`` then ``'method'``; the latter is what spatial
        datasets use for the profiling platform.
    labels_key : str, optional
        Column in ``adata.obs`` holding the cell-type labels. If omitted, falls
        back to ``'cell_type'``, ``'annotation'``, then ``'labels'``.

    Returns
    -------
    Dict[str, Union[float, List[str]]]
        Analysis results including signal leakage and interpretability metrics
    """
    results = {}

    # 1. Test for signal leakage (residual shouldn't predict batch/labels well)
    #
    # Resolve the obs columns rather than assuming their names. Passing
    # `batch_key` explicitly is preferred; otherwise fall back to a candidate
    # list. The fallback matters for spatial datasets, where the batch variable
    # is `method` (the profiling platform) and not `batch` -- assuming 'batch'
    # silently skipped the batch-leakage test on those datasets.
    batch_candidates = [batch_key] if batch_key else ['batch', 'method']
    batch_col = next((c for c in batch_candidates if c and c in adata.obs), None)
    if batch_col is not None:
        batch_pred = evaluate_subspace_predictability(z_residual, adata.obs[batch_col])
        results["batch_leakage_accuracy"] = batch_pred["accuracy"]
        results["batch_leakage_acceptable"] = batch_pred["accuracy"] < 0.7  # Should be low
        results["batch_key_used"] = batch_col

    label_candidates = ([labels_key] if labels_key
                        else ['cell_type', 'annotation', 'labels'])
    label_col = next((c for c in label_candidates if c and c in adata.obs), None)
    if label_col is not None:
        label_pred = evaluate_subspace_predictability(z_residual, adata.obs[label_col])
        results["label_leakage_accuracy"] = label_pred["accuracy"]
        results["label_leakage_acceptable"] = label_pred["accuracy"] < 0.7
        results["labels_key_used"] = label_col

    # 2. Orthogonality with biological and batch subspaces
    bio_orthogonality = compute_orthogonality_score(z_residual, z_biological)
    batch_orthogonality = compute_orthogonality_score(z_residual, z_batch)

    results["biological_orthogonality"] = bio_orthogonality["orthogonality_score"]
    results["batch_orthogonality"] = batch_orthogonality["orthogonality_score"]
    results["orthogonality_acceptable"] = (
        bio_orthogonality["orthogonality_score"] > 0.7 and
        batch_orthogonality["orthogonality_score"] > 0.7
    )

    # 3. Gene correlation analysis - what does residual capture?
    X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X

    residual_gene_correlations = []
    top_genes_per_factor = []

    for i in range(z_residual.shape[1]):
        gene_corrs = []
        for j in range(X.shape[1]):
            corr, p_val = pearsonr(z_residual[:, i], X[:, j])
            if not np.isnan(corr):
                gene_corrs.append((abs(corr), adata.var.index[j], corr, p_val))

        # Sort by absolute correlation
        gene_corrs.sort(reverse=True)
        residual_gene_correlations.extend([gc[0] for gc in gene_corrs[:n_top_genes]])

        # Top genes for this factor
        top_genes = [gc[1] for gc in gene_corrs[:20]]
        top_genes_per_factor.append(top_genes)

    results["mean_gene_correlation"] = np.mean(residual_gene_correlations) if residual_gene_correlations else 0.0
    results["top_correlated_genes"] = top_genes_per_factor

    # 4. Variance explained comparison
    residual_var = np.var(z_residual, axis=0).sum()
    bio_var = np.var(z_biological, axis=0).sum()
    batch_var = np.var(z_batch, axis=0).sum()
    total_var = residual_var + bio_var + batch_var

    results["residual_variance_fraction"] = residual_var / total_var if total_var > 0 else 0.0
    results["variance_reasonable"] = 0.1 <= results["residual_variance_fraction"] <= 0.4  # Should be meaningful but not dominant

    # 5. Overall interpretability assessment
    interpretability_checks = [
        results.get("batch_leakage_acceptable", True),
        results.get("label_leakage_acceptable", True),
        results["orthogonality_acceptable"],
        results["variance_reasonable"],
        results["mean_gene_correlation"] > 0.05  # Should correlate with some genes
    ]

    results["interpretability_score"] = sum(interpretability_checks) / len(interpretability_checks)
    results["is_signal_sink"] = results["interpretability_score"] < 0.6

    return results


def evaluate_disentanglement_quality(
    z_batch: np.ndarray,
    z_biological: np.ndarray,
    z_residual: np.ndarray,
    adata: anndata.AnnData,
    batch_key: str = "batch",
    label_key: str = "cell_type"
) -> Dict[str, Union[float, Dict, bool]]:
    """
    Comprehensive disentanglement quality evaluation.

    Validates that FADVI achieves proper factor separation:
    - z_batch should predict batch effects well
    - z_biological should predict cell types well
    - z_residual should be orthogonal and interpretable

    Parameters
    ----------
    z_batch : np.ndarray
        Batch latent factors
    z_biological : np.ndarray
        Biological latent factors
    z_residual : np.ndarray
        Residual latent factors
    adata : anndata.AnnData
        Data with batch and label annotations
    batch_key : str
        Column name for batch information
    label_key : str
        Column name for biological labels

    Returns
    -------
    Dict[str, Union[float, Dict, bool]]
        Comprehensive disentanglement evaluation results
    """
    results = {
        "batch_predictability": {},
        "biological_predictability": {},
        "cross_correlations": {},
        "orthogonality_scores": {},
        "residual_analysis": {},
        "overall_quality": {}
    }

    # 1. Batch predictability (z_batch should predict batch well)
    if batch_key in adata.obs:
        batch_pred = evaluate_subspace_predictability(z_batch, adata.obs[batch_key])
        results["batch_predictability"] = batch_pred
        batch_quality = batch_pred["accuracy"] > 0.8
    else:
        batch_quality = False
        print(f"Warning: {batch_key} not found in adata.obs")

    # 2. Biological predictability (z_biological should predict cell types well)
    label_col = label_key
    if label_key not in adata.obs:
        # Try common alternatives
        for alt in ['annotation', 'cell_type', 'labels', 'celltype']:
            if alt in adata.obs:
                label_col = alt
                break

    if label_col in adata.obs:
        bio_pred = evaluate_subspace_predictability(z_biological, adata.obs[label_col])
        results["biological_predictability"] = bio_pred
        bio_quality = bio_pred["accuracy"] > 0.8
    else:
        bio_quality = False
        print(f"Warning: No cell type labels found in adata.obs")

    # 3. Cross-correlations (should be low between subspaces)
    batch_bio_corr, _ = compute_cross_correlation(z_batch, z_biological)
    batch_res_corr, _ = compute_cross_correlation(z_batch, z_residual)
    bio_res_corr, _ = compute_cross_correlation(z_biological, z_residual)

    results["cross_correlations"] = {
        "batch_biological": batch_bio_corr,
        "batch_residual": batch_res_corr,
        "biological_residual": bio_res_corr
    }

    corr_quality = all(corr < 0.3 for corr in [batch_bio_corr, batch_res_corr, bio_res_corr])

    # 4. Orthogonality scores
    batch_bio_orth = compute_orthogonality_score(z_batch, z_biological)
    batch_res_orth = compute_orthogonality_score(z_batch, z_residual)
    bio_res_orth = compute_orthogonality_score(z_biological, z_residual)

    results["orthogonality_scores"] = {
        "batch_biological": batch_bio_orth,
        "batch_residual": batch_res_orth,
        "biological_residual": bio_res_orth
    }

    orth_quality = all(
        orth["orthogonality_score"] > 0.7
        for orth in [batch_bio_orth, batch_res_orth, bio_res_orth]
    )

    # 5. Residual analysis
    residual_analysis = analyze_residual_interpretability(
        z_residual, z_biological, z_batch, adata
    )
    results["residual_analysis"] = residual_analysis
    residual_quality = not residual_analysis["is_signal_sink"]

    # 6. Overall quality assessment
    quality_checks = {
        "batch_specificity": batch_quality,
        "biological_specificity": bio_quality,
        "low_cross_correlation": corr_quality,
        "high_orthogonality": orth_quality,
        "residual_interpretable": residual_quality
    }

    overall_score = sum(quality_checks.values()) / len(quality_checks)

    results["overall_quality"] = {
        "individual_checks": quality_checks,
        "overall_score": overall_score,
        "disentanglement_successful": overall_score >= 0.8,
        "meets_quality_criteria": (
            batch_quality and bio_quality and
            corr_quality and residual_quality
        )
    }

    return results


def compare_disentanglement_methods(
    embeddings_dict: Dict[str, Dict[str, np.ndarray]],
    adata: anndata.AnnData,
    batch_key: str = "batch",
    label_key: str = "cell_type"
) -> pd.DataFrame:
    """
    Compare disentanglement quality across multiple methods.

    Parameters
    ----------
    embeddings_dict : Dict[str, Dict[str, np.ndarray]]
        Nested dict with structure: {method_name: {subspace_name: embedding}}
        e.g., {"FADVI": {"z_b": batch_embed, "z_l": bio_embed, "z_r": res_embed}}
    adata : anndata.AnnData
        Data for evaluation
    batch_key : str
        Batch annotation column
    label_key : str
        Cell type annotation column

    Returns
    -------
    pd.DataFrame
        Comparison results across methods
    """
    comparison_results = []

    for method_name, subspaces in embeddings_dict.items():
        print(f"Evaluating {method_name}...")

        # Extract subspaces (handle different naming conventions).
        # NOTE: use an explicit None check rather than `a or b` -- numpy arrays
        # raise "truth value of an array is ambiguous" under boolean `or`.
        def _pick(*names):
            for nm in names:
                arr = subspaces.get(nm)
                if arr is not None:
                    return arr
            return None

        z_batch = _pick("z_b", "batch")
        z_biological = _pick("z_l", "biological", "label")
        z_residual = _pick("z_r", "residual")

        if z_batch is None or z_biological is None:
            print(f"  Warning: Missing required subspaces for {method_name}")
            continue

        # Use zeros if no residual provided
        if z_residual is None:
            z_residual = np.zeros((z_batch.shape[0], 1))

        try:
            results = evaluate_disentanglement_quality(
                z_batch, z_biological, z_residual, adata, batch_key, label_key
            )

            # Flatten results for comparison table
            row = {
                "method": method_name,
                "batch_accuracy": results["batch_predictability"].get("accuracy", 0),
                "biological_accuracy": results["biological_predictability"].get("accuracy", 0),
                "batch_bio_correlation": results["cross_correlations"]["batch_biological"],
                "batch_residual_correlation": results["cross_correlations"]["batch_residual"],
                "bio_residual_correlation": results["cross_correlations"]["biological_residual"],
                "batch_bio_orthogonality": results["orthogonality_scores"]["batch_biological"]["orthogonality_score"],
                "residual_interpretable": not results["residual_analysis"]["is_signal_sink"],
                "overall_score": results["overall_quality"]["overall_score"],
                "meets_criteria": results["overall_quality"]["meets_quality_criteria"]
            }

            comparison_results.append(row)

        except Exception as e:
            print(f"  Error evaluating {method_name}: {e}")
            continue

    if comparison_results:
        return pd.DataFrame(comparison_results)
    else:
        print("No methods could be evaluated successfully")
        return pd.DataFrame()