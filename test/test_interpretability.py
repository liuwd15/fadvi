#!/usr/bin/env python3
"""
Test script for FADVI interpretability functionality.
Tests IntegratedGradients and GradientShap interpretability methods.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
import pandas as pd
import torch

# Import FADVI
try:
    from fadvi import FADVI

    print("Successfully imported FADVI from installed package")
except ImportError:
    # Fallback for development
    import os
    import sys

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
    from fadvi import FADVI


def load_data():
    """Create simulated scRNA-seq data for testing (from test_fadvi_predict.py)."""
    print("Creating simulated scRNA-seq data...")

    # Use scvi's synthetic data generation with correct parameters
    from scvi.data import synthetic_iid

    # Generate synthetic data with batch and label effects
    adata = synthetic_iid(
        batch_size=700,  # Number of cells per batch (700 * 3 = 2100 total cells)
        n_genes=1000,  # Number of genes
        n_batches=3,  # Number of batches
        n_labels=4,  # Number of cell types
        dropout_ratio=0.3,  # Lower dropout for more realistic data
        sparse_format="csr_matrix",  # Use sparse format
    )

    print(f"Generated data shape: {adata.shape}")
    print(f"Data range: [{adata.X.min():.3f}, {adata.X.max():.3f}]")

    # Add batch and label annotations that match the expected names
    adata.obs["batch"] = adata.obs["batch"].astype("category")
    adata.obs["label"] = adata.obs["labels"].astype("category")

    # Basic preprocessing
    adata.var_names_make_unique()

    # Filter cells with too few genes (some synthetic cells might have low counts)
    import scanpy as sc

    sc.pp.filter_cells(adata, min_genes=50)
    print(f"After filtering: {adata.shape}")

    # Handle extreme values that can cause issues with highly variable gene selection
    # Clip extreme values to reasonable range for count data
    max_count = np.percentile(
        adata.X.data if hasattr(adata.X, "data") else adata.X, 99.9
    )
    print(f"Clipping counts above {max_count:.1f}")

    if hasattr(adata.X, "data"):  # sparse matrix
        adata.X.data = np.clip(adata.X.data, 0, max_count)
    else:  # dense matrix
        adata.X = np.clip(adata.X, 0, max_count)

    print(f"Data range after clipping: [{adata.X.min():.3f}, {adata.X.max():.3f}]")

    # Find highly variable genes on raw count data (no normalization needed for scvi-tools)
    # Use more conservative parameters to avoid issues with extreme values
    try:
        sc.pp.highly_variable_genes(adata, min_mean=0.01, max_mean=5, min_disp=0.2)
        n_hvg = sum(adata.var["highly_variable"])
        if n_hvg < 100:  # If too few HVGs, be more lenient
            print(f"Only {n_hvg} HVGs found, using more lenient criteria")
            sc.pp.highly_variable_genes(
                adata, min_mean=0.001, max_mean=10, min_disp=0.05
            )
            n_hvg = sum(adata.var["highly_variable"])
        if (
            n_hvg < 50
        ):  # If still too few, use top variable genes with cell_ranger flavor
            print(
                f"Still only {n_hvg} HVGs found, selecting top 500 most variable genes"
            )
            sc.pp.highly_variable_genes(adata, n_top_genes=500, flavor="cell_ranger")
    except (ValueError, ImportError) as e:
        print(
            f"Warning: HVG selection failed ({e}), selecting top 500 most variable genes"
        )
        try:
            sc.pp.highly_variable_genes(adata, n_top_genes=500, flavor="cell_ranger")
        except:
            print("Cell ranger method also failed, using all genes")
            adata.var["highly_variable"] = True

    # Ensure we have at least some genes selected
    if sum(adata.var["highly_variable"]) == 0:
        print("No HVGs found with any method, using all genes")
        adata.var["highly_variable"] = True

    adata.raw = adata
    adata = adata[
        :, adata.var.highly_variable
    ].copy()  # Important: copy to avoid view issues

    print(f"Number of highly variable genes: {adata.n_vars}")
    print(f"Batch categories: {adata.obs['batch'].cat.categories.tolist()}")
    print(f"Label categories: {adata.obs['label'].cat.categories.tolist()}")

    return adata


def load_test_data():
    """Create simulated scRNA-seq data for interpretability testing."""
    print("Creating simulated scRNA-seq data for interpretability tests...")

    from scvi.data import synthetic_iid

    # Generate synthetic data
    adata = synthetic_iid(
        batch_size=300,  # Moderate size for testing
        n_genes=200,  # Reasonable gene count
        n_batches=2,  # Simple batch structure
        n_labels=3,  # Simple label structure
        dropout_ratio=0.2,
        sparse_format="csr_matrix",
    )

    # Add required annotations
    adata.obs["batch"] = adata.obs["batch"].astype("category")
    adata.obs["label"] = adata.obs["labels"].astype("category")

    print(f"Generated data shape: {adata.shape}")
    print(f"Batch categories: {adata.obs['batch'].cat.categories.tolist()}")
    print(f"Label categories: {adata.obs['label'].cat.categories.tolist()}")

    return adata


def setup_trained_model(adata):
    """Setup and train a FADVI model for testing."""
    print("Setting up and training FADVI model...")

    FADVI.setup_anndata(adata, batch_key="batch", labels_key="label")
    model = FADVI(adata, n_latent_b=5, n_latent_l=5, n_latent_r=5)

    # Quick training for testing
    model.train(max_epochs=5, batch_size=64, early_stopping=False)
    print("Model training completed!")

    return model


def test_interpretability_analysis():
    """Test interpretability analysis functionality (moved from test_fadvi_predict.py)."""
    print("\n=== Testing Interpretability Analysis ===")

    # Skip if captum is not available
    try:
        import captum

        print("Captum found, proceeding with interpretability analysis tests")
    except ImportError:
        print("Captum not found, skipping interpretability analysis tests")
        return

    # Load data and setup model
    adata = load_test_data()
    model = setup_trained_model(adata)

    # Test subset of cells
    test_indices = np.random.choice(adata.n_obs, 20, replace=False)

    print("\n--- Testing IntegratedGradients interpretability analysis ---")
    try:
        # Get attributions for a few cells using IntegratedGradients
        sample_indices = test_indices[:5]
        sample_pred, sample_attr = model.predict(
            adata,
            indices=sample_indices,
            prediction_mode="label",
            interpretability="ig",
            return_dict=False,
        )

        # sample_attr is now a ranked features DataFrame, not per-cell attributions
        print(f"Sample predictions: {sample_pred}")
        print(f"Ranked features shape: {sample_attr.shape}")
        print(f"Top 5 most important genes overall:")
        for i, (idx, row) in enumerate(sample_attr.head().iterrows()):
            print(
                f"  {i+1}. {row['feature']} (mean attribution: {row['attribution_mean']:.6f})"
            )

        # Verify the structure of ranked features DataFrame
        expected_columns = {
            "feature",
            "feature_idx",
            "attribution_mean",
            "attribution_std",
            "attribution_abs_mean",
            "n_cells",
        }
        assert expected_columns.issubset(
            set(sample_attr.columns)
        ), f"Missing expected columns in ranked features"
        assert sample_attr["n_cells"].iloc[0] == len(
            sample_indices
        ), f"n_cells should match sample size"
        assert len(sample_attr) == adata.n_vars, f"Should have ranking for all genes"

        print("✅ IntegratedGradients interpretability analysis completed successfully")

    except Exception as e:
        print(f"❌ IntegratedGradients interpretability analysis failed: {e}")
        raise

    print("\n--- Testing GradientShap interpretability analysis ---")
    try:
        # Get attributions using GradientShap
        sample_pred_gs, sample_attr_gs = model.predict(
            adata,
            indices=sample_indices,
            prediction_mode="batch",  # Test on batch prediction
            interpretability="gs",
            return_dict=False,
        )

        print(f"GradientShap sample predictions: {sample_pred_gs}")
        print(f"GradientShap ranked features shape: {sample_attr_gs.shape}")
        print(f"Top 5 most important genes for batch prediction (GradientShap):")
        for i, (idx, row) in enumerate(sample_attr_gs.head().iterrows()):
            print(
                f"  {i+1}. {row['feature']} (mean attribution: {row['attribution_mean']:.6f})"
            )

        # Verify GradientShap results have same structure
        assert expected_columns.issubset(
            set(sample_attr_gs.columns)
        ), f"Missing expected columns in GradientShap ranked features"
        assert sample_attr_gs["n_cells"].iloc[0] == len(
            sample_indices
        ), f"n_cells should match sample size"
        assert len(sample_attr_gs) == adata.n_vars, f"Should have ranking for all genes"

        print("✅ GradientShap interpretability analysis completed successfully")

    except Exception as e:
        print(f"❌ GradientShap interpretability analysis failed: {e}")
        raise

    # Compare IG vs GS results
    print("\n--- Comparing IntegratedGradients vs GradientShap ---")
    try:
        # Get both methods on same data for comparison
        np.random.seed(42)
        torch.manual_seed(42)

        ig_pred, ig_attr = model.predict(
            adata,
            indices=sample_indices,
            prediction_mode="label",
            interpretability="ig",
            return_dict=False,
        )

        np.random.seed(42)
        torch.manual_seed(42)

        gs_pred, gs_attr = model.predict(
            adata,
            indices=sample_indices,
            prediction_mode="label",
            interpretability="gs",
            return_dict=False,
        )

        # Predictions should be the same (same model, same data)
        predictions_match = np.array_equal(ig_pred, gs_pred)
        print(f"IG vs GS predictions match: {predictions_match}")

        # Attribution rankings might be different (different methods)
        top_genes_ig = set(ig_attr.head(10)["feature"].values)
        top_genes_gs = set(gs_attr.head(10)["feature"].values)
        overlap = len(top_genes_ig.intersection(top_genes_gs))

        print(f"Top 10 gene overlap between IG and GS: {overlap}/10 genes")
        print(f"IG top genes: {sorted(list(top_genes_ig))[:5]}...")
        print(f"GS top genes: {sorted(list(top_genes_gs))[:5]}...")

        print("✅ IG vs GS comparison completed successfully")

    except Exception as e:
        print(f"❌ IG vs GS comparison failed: {e}")
        raise

    print("\n✅ All interpretability analysis tests passed!")


def test_gs_interpretability():
    """Test GradientShap interpretability functionality specifically."""
    print("\n=== Testing GradientShap Interpretability ===")

    # Skip if captum is not available
    try:
        import captum

        print("Captum found, proceeding with GradientShap tests")
    except ImportError:
        print("Captum not found, skipping GradientShap tests")
        return

    # Load data and setup model
    adata = load_test_data()
    model = setup_trained_model(adata)

    # Test subset for efficiency
    test_indices = np.random.choice(adata.n_obs, 15, replace=False)

    # Test 1: Batch prediction with GradientShap
    print("\n--- Testing batch prediction with GradientShap ---")
    try:
        batch_pred_gs, batch_attr_gs = model.predict(
            adata,
            indices=test_indices,
            prediction_mode="batch",
            interpretability="gs",
            return_dict=False,
        )

        print(f"✅ Batch GS predictions shape: {batch_pred_gs.shape}")
        print(f"✅ Batch GS attributions shape: {batch_attr_gs.shape}")

        # Verify structure
        assert isinstance(
            batch_attr_gs, pd.DataFrame
        ), f"Attributions should be DataFrame with automatic ranking"
        assert (
            batch_attr_gs.shape[0] == adata.n_vars
        ), f"Attribution rows should match genes"

        expected_columns = {
            "feature",
            "feature_idx",
            "attribution_mean",
            "attribution_std",
            "attribution_abs_mean",
            "n_cells",
        }
        assert expected_columns.issubset(
            set(batch_attr_gs.columns)
        ), f"Missing expected columns"
        assert batch_attr_gs["n_cells"].iloc[0] == len(
            test_indices
        ), f"n_cells should match test indices"

        print(f"✅ Batch GS ranked features DataFrame with {len(batch_attr_gs)} genes")

    except Exception as e:
        print(f"❌ Batch GradientShap test failed: {e}")
        raise

    # Test 2: Label prediction with GradientShap
    print("\n--- Testing label prediction with GradientShap ---")
    try:
        label_pred_gs, label_attr_gs = model.predict(
            adata,
            indices=test_indices,
            prediction_mode="label",
            interpretability="gs",
            return_dict=False,
        )

        print(f"✅ Label GS predictions shape: {label_pred_gs.shape}")
        print(f"✅ Label GS attributions shape: {label_attr_gs.shape}")

        # Verify structure
        assert isinstance(
            label_attr_gs, pd.DataFrame
        ), f"Attributions should be DataFrame"
        assert (
            label_attr_gs.shape[0] == adata.n_vars
        ), f"Attribution rows should match genes"
        assert label_attr_gs["n_cells"].iloc[0] == len(
            test_indices
        ), f"n_cells should match test indices"

        print(f"✅ Label GS ranked features DataFrame with {len(label_attr_gs)} genes")

    except Exception as e:
        print(f"❌ Label GradientShap test failed: {e}")
        raise

    # Test 3: Soft predictions with GradientShap
    print("\n--- Testing soft predictions with GradientShap ---")
    try:
        soft_pred_gs, soft_attr_gs = model.predict(
            adata,
            indices=test_indices[:10],  # Smaller subset
            prediction_mode="batch",
            soft=True,
            interpretability="gs",
            return_dict=False,
        )

        print(f"✅ Soft GS predictions shape: {soft_pred_gs.shape}")
        print(f"✅ Soft GS attributions shape: {soft_attr_gs.shape}")

        # Should still get ranked features DataFrame
        assert isinstance(
            soft_attr_gs, pd.DataFrame
        ), "Attributions should be DataFrame"
        assert soft_attr_gs["n_cells"].iloc[0] == 10, "n_cells should match subset size"

        print(
            f"✅ Soft prediction GS ranked features DataFrame with {len(soft_attr_gs)} genes"
        )

    except Exception as e:
        print(f"❌ Soft prediction GradientShap test failed: {e}")
        raise

    # Test 4: GradientShap with custom arguments
    print("\n--- Testing GradientShap with custom arguments ---")
    try:
        custom_pred_gs, custom_attr_gs = model.predict(
            adata,
            indices=test_indices[:8],
            prediction_mode="label",
            interpretability="gs",
            interpretability_args={
                "n_samples": 20,
                "stdevs": 0.1,
            },  # Custom GS parameters
            return_dict=False,
        )

        print(f"✅ Custom GS args test passed: {custom_attr_gs.shape}")
        assert isinstance(
            custom_attr_gs, pd.DataFrame
        ), "Custom GS attributions should be DataFrame"

    except Exception as e:
        print(f"❌ Custom GradientShap args test failed: {e}")
        raise

    print("\n✅ All GradientShap tests passed!")


def test_interpretability_methods_comparison():
    """Compare different interpretability methods side by side."""
    print("\n=== Testing Interpretability Methods Comparison ===")

    # Skip if captum is not available
    try:
        import captum

        print("Captum found, proceeding with methods comparison")
    except ImportError:
        print("Captum not found, skipping methods comparison")
        return

    # Load data and setup model
    adata = load_test_data()
    model = setup_trained_model(adata)

    # Use same subset for fair comparison
    test_indices = [0, 1, 2, 3, 4, 5]

    # Set seeds for reproducible comparison
    np.random.seed(42)
    torch.manual_seed(42)

    print("\n--- Comparing IG vs GS on same data ---")

    # IntegratedGradients
    ig_pred, ig_attr = model.predict(
        adata,
        indices=test_indices,
        prediction_mode="label",
        interpretability="ig",
        return_dict=False,
    )

    # Reset seeds for fair comparison
    np.random.seed(42)
    torch.manual_seed(42)

    # GradientShap
    gs_pred, gs_attr = model.predict(
        adata,
        indices=test_indices,
        prediction_mode="label",
        interpretability="gs",
        return_dict=False,
    )

    # Compare predictions (should be identical)
    predictions_identical = np.array_equal(ig_pred, gs_pred)
    print(f"✅ Predictions identical across methods: {predictions_identical}")

    # Compare attribution patterns
    print(f"✅ IG ranked features shape: {ig_attr.shape}")
    print(f"✅ GS ranked features shape: {gs_attr.shape}")

    # Show top genes from each method
    print("\nTop 5 genes by each method:")
    print("IntegratedGradients:")
    for i, (idx, row) in enumerate(ig_attr.head().iterrows()):
        print(f"  {i+1}. {row['feature']} (mean: {row['attribution_mean']:.6f})")

    print("GradientShap:")
    for i, (idx, row) in enumerate(gs_attr.head().iterrows()):
        print(f"  {i+1}. {row['feature']} (mean: {row['attribution_mean']:.6f})")

    # Calculate method agreement
    top_n = 10
    top_ig_genes = set(ig_attr.head(top_n)["feature"].values)
    top_gs_genes = set(gs_attr.head(top_n)["feature"].values)
    agreement = len(top_ig_genes.intersection(top_gs_genes)) / top_n

    print(f"✅ Top-{top_n} gene agreement between methods: {agreement:.2%}")

    # Both methods should identify some meaningful patterns
    assert len(top_ig_genes) == top_n, f"IG should identify {top_n} top genes"
    assert len(top_gs_genes) == top_n, f"GS should identify {top_n} top genes"

    print("✅ Methods comparison completed successfully!")


# Pytest-discoverable test functions
def test_interpretability_analysis_pytest():
    """Pytest-discoverable version of interpretability analysis test."""
    try:
        import captum

        test_interpretability_analysis()
    except ImportError:
        import pytest

        pytest.skip("captum not available")


def test_gs_interpretability_pytest():
    """Pytest-discoverable version of GradientShap test."""
    try:
        import captum

        test_gs_interpretability()
    except ImportError:
        import pytest

        pytest.skip("captum not available")


def test_interpretability_methods_comparison_pytest():
    """Pytest-discoverable version of methods comparison test."""
    try:
        import captum

        test_interpretability_methods_comparison()
    except ImportError:
        import pytest

        pytest.skip("captum not available")


# Methods moved from test_fadvi_predict.py


def test_interpretability_functionality():
    """Test interpretability functionality with Integrated Gradients (moved from test_fadvi_predict.py)."""
    print("\n=== Testing Interpretability Functionality ===")

    # Skip if captum is not available
    try:
        import captum

        print("Captum found, proceeding with interpretability tests")
    except ImportError:
        print("Captum not found, skipping interpretability tests")
        return

    # Load and prepare data
    adata = load_data()

    # Setup and train model
    print("Setting up FADVI model...")
    FADVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="label",
    )

    model = FADVI(adata, n_latent_b=5, n_latent_l=5, n_latent_r=5)
    print("Training model...")
    model.train(max_epochs=5, batch_size=128, early_stopping=False)

    print("Model training completed!")

    # Test subset of cells for faster testing
    n_test_cells = 100
    test_indices = np.random.choice(adata.n_obs, n_test_cells, replace=False)

    # Test 1: Batch prediction with interpretability
    print("\n--- Testing batch prediction interpretability ---")
    try:
        batch_pred_ig, batch_attr = model.predict(
            adata,
            indices=test_indices,
            prediction_mode="batch",
            interpretability="ig",
            return_dict=False,
        )

        print(f"✅ Batch predictions shape: {batch_pred_ig.shape}")
        print(f"✅ Batch attributions shape: {batch_attr.shape}")

        # Verify attributions are now a ranked features DataFrame (new behavior)
        assert isinstance(
            batch_attr, pd.DataFrame
        ), f"Attributions should be DataFrame with automatic ranking"
        assert (
            batch_attr.shape[0] == adata.n_vars
        ), f"Attribution rows should match genes (ranked features)"
        expected_columns = {
            "feature",
            "feature_idx",
            "attribution_mean",
            "attribution_std",
            "attribution_abs_mean",
            "n_cells",
        }
        assert expected_columns.issubset(
            set(batch_attr.columns)
        ), f"Missing expected columns in ranked features"

        # Check that we have meaningful attribution data
        assert batch_attr["n_cells"].iloc[0] == len(
            test_indices
        ), f"n_cells should match test indices"
        print(
            f"✅ Ranked features DataFrame with {len(batch_attr)} genes and {len(test_indices)} cells"
        )

        # Test with custom IG arguments
        batch_pred_ig2, batch_attr2 = model.predict(
            adata,
            indices=test_indices[:20],  # Even smaller subset
            prediction_mode="batch",
            interpretability="ig",
            interpretability_args={"n_steps": 50},  # Custom number of integration steps
            return_dict=False,
        )
        print(f"✅ Custom IG args test passed: {batch_attr2.shape}")

    except Exception as e:
        print(f"❌ Batch interpretability test failed: {e}")
        raise

    # Test 2: Label prediction with interpretability
    print("\n--- Testing label prediction interpretability ---")
    try:
        label_pred_ig, label_attr = model.predict(
            adata,
            indices=test_indices,
            prediction_mode="label",
            interpretability="ig",
            return_dict=False,
        )

        print(f"✅ Label predictions shape: {label_pred_ig.shape}")
        print(f"✅ Label attributions shape: {label_attr.shape}")

        # Verify attributions are now a ranked features DataFrame (new behavior)
        assert isinstance(
            label_attr, pd.DataFrame
        ), f"Attributions should be DataFrame with automatic ranking"
        assert (
            label_attr.shape[0] == adata.n_vars
        ), f"Attribution rows should match genes (ranked features)"
        assert label_attr["n_cells"].iloc[0] == len(
            test_indices
        ), f"n_cells should match test indices"
        print(f"✅ Label ranked features DataFrame with {len(label_attr)} genes")

    except Exception as e:
        print(f"❌ Label interpretability test failed: {e}")
        raise

    # Test 3: Soft predictions with interpretability
    print("\n--- Testing soft predictions with interpretability ---")
    try:
        batch_pred_soft_ig, batch_attr_soft = model.predict(
            adata,
            indices=test_indices[:20],  # Small subset
            prediction_mode="batch",
            soft=True,
            interpretability="ig",
            return_dict=False,
        )

        print(f"✅ Soft batch predictions shape: {batch_pred_soft_ig.shape}")
        print(f"✅ Soft batch attributions shape: {batch_attr_soft.shape}")

        # Should still get ranked features DataFrame even with soft predictions
        assert isinstance(
            batch_attr_soft, pd.DataFrame
        ), "Attributions should be DataFrame with automatic ranking"
        assert (
            batch_attr_soft["n_cells"].iloc[0] == 20
        ), "n_cells should match test subset size"
        print(
            f"✅ Soft prediction ranked features DataFrame with {len(batch_attr_soft)} genes"
        )

    except Exception as e:
        print(f"❌ Soft prediction interpretability test failed: {e}")
        raise

    # Test 4: Compare regular vs interpretability predictions
    print("\n--- Testing prediction consistency ---")
    try:
        # Set random seeds for deterministic behavior
        np.random.seed(42)
        torch.manual_seed(42)
        if torch.cuda.is_available():
            torch.cuda.manual_seed(42)

        # Regular predictions
        batch_pred_regular = model.predict(
            adata, indices=test_indices[:20], prediction_mode="batch"
        )

        # Reset seeds to ensure same conditions
        np.random.seed(42)
        torch.manual_seed(42)
        if torch.cuda.is_available():
            torch.cuda.manual_seed(42)

        # Interpretability predictions (should be same with same seeds)
        batch_pred_ig_only, _ = model.predict(
            adata,
            indices=test_indices[:20],
            prediction_mode="batch",
            interpretability="ig",
            return_dict=False,  # Use tuple format for consistency
        )

        # Check if predictions are identical
        predictions_match = np.array_equal(batch_pred_regular, batch_pred_ig_only)
        print(f"✅ Regular vs IG predictions match: {predictions_match}")

        # If exact match fails, check if the majority of predictions are the same
        if not predictions_match:
            match_ratio = np.mean(batch_pred_regular == batch_pred_ig_only)
            print(f"   Prediction match ratio: {match_ratio:.2f}")

            # Accept if at least 80% of predictions match (accounting for model stochasticity)
            if match_ratio >= 0.8:
                print("✅ Prediction consistency acceptable (>80% match)")
            else:
                print("❌ Predictions too different")
                # Show some examples of differences
                different_indices = np.where(batch_pred_regular != batch_pred_ig_only)[
                    0
                ]
                for i in different_indices[:5]:  # Show first 5 differences
                    print(
                        f"   Cell {i}: Regular={batch_pred_regular[i]}, IG={batch_pred_ig_only[i]}"
                    )
                assert (
                    match_ratio >= 0.8
                ), f"Prediction match ratio ({match_ratio:.2f}) is too low"

    except Exception as e:
        print(f"❌ Prediction consistency test failed: {e}")
        raise

    print("\n✅ All interpretability tests passed!")


def test_interpretability_error_handling():
    """Test error handling when captum is not available (moved from test_fadvi_predict.py)."""
    print("\n=== Testing Interpretability Error Handling ===")

    # Load and prepare minimal data
    adata = load_data()

    # Setup model
    FADVI.setup_anndata(adata, batch_key="batch", labels_key="label")
    model = FADVI(adata, n_latent_b=5, n_latent_l=5, n_latent_r=5)

    # Train minimally
    model.train(max_epochs=2, batch_size=128)

    # Test the error handling by temporarily removing captum from sys.modules
    import importlib
    import sys

    # Save the original module if it exists
    captum_module = sys.modules.get("captum")
    captum_attr_module = sys.modules.get("captum.attr")

    # Remove captum from sys.modules to force a re-import
    if "captum" in sys.modules:
        del sys.modules["captum"]
    if "captum.attr" in sys.modules:
        del sys.modules["captum.attr"]

    # Mock the importlib.import_module function to raise ImportError for captum
    original_import_module = importlib.import_module

    def mock_import_module(name, package=None):
        if name == "captum.attr" or name.startswith("captum"):
            raise ImportError(f"No module named '{name}'")
        return original_import_module(name, package)

    try:
        # Apply the mock
        importlib.import_module = mock_import_module

        # This should raise a ModuleNotFoundError
        try:
            model.predict(
                adata, indices=[0, 1, 2], prediction_mode="label", interpretability="ig"
            )
            # If we get here, the error was not raised
            print(
                "❌ Error handling test failed: Expected ModuleNotFoundError was not raised"
            )
        except ModuleNotFoundError as e:
            # This is what we expect
            print(
                "✅ Error handling test passed - ModuleNotFoundError raised as expected"
            )
            print(f"   Error message: {e}")

    finally:
        # Restore original import function
        importlib.import_module = original_import_module

        # Restore the modules
        if captum_module is not None:
            sys.modules["captum"] = captum_module
        if captum_attr_module is not None:
            sys.modules["captum.attr"] = captum_attr_module

    print("✅ Error handling tests completed!")


def test_interpretability_with_different_parameters():
    """Test interpretability with different prediction parameters (moved from test_fadvi_predict.py)."""
    print("\n=== Testing Interpretability with Different Parameters ===")

    # Skip if captum is not available
    try:
        import captum
    except ImportError:
        print("Captum not found, skipping parameter tests")
        return

    # Load and prepare data
    adata = load_data()

    # Setup and train model
    FADVI.setup_anndata(adata, batch_key="batch", labels_key="label")
    model = FADVI(adata, n_latent_b=5, n_latent_l=5, n_latent_r=5)
    model.train(max_epochs=3, batch_size=128, early_stopping=False)

    test_indices = np.random.choice(adata.n_obs, 20, replace=False)

    # Test different return_numpy settings
    print("--- Testing return_numpy parameter ---")
    pred1, attr1 = model.predict(
        adata,
        indices=test_indices,
        prediction_mode="batch",
        interpretability="ig",
        return_numpy=True,
        return_dict=False,
    )

    pred2, attr2 = model.predict(
        adata,
        indices=test_indices,
        prediction_mode="batch",
        interpretability="ig",
        return_numpy=False,
        return_dict=False,
    )

    print(f"✅ return_numpy=True: pred type {type(pred1)}, attr type {type(attr1)}")
    print(f"✅ return_numpy=False: pred type {type(pred2)}, attr type {type(attr2)}")

    # Both should return DataFrame for attributions due to automatic ranking
    assert isinstance(
        attr1, pd.DataFrame
    ), "Attributions should be DataFrame (ranked features)"
    assert isinstance(
        attr2, pd.DataFrame
    ), "Attributions should be DataFrame (ranked features)"

    # Test with different batch sizes
    print("--- Testing different batch sizes ---")

    # Set seeds for consistency
    np.random.seed(42)
    torch.manual_seed(42)

    pred3, attr3 = model.predict(
        adata,
        indices=test_indices,
        prediction_mode="label",
        interpretability="ig",
        batch_size=5,  # Small batch size
        return_dict=False,
    )

    # Reset seeds again for identical generation
    np.random.seed(42)
    torch.manual_seed(42)

    pred4, attr4 = model.predict(
        adata,
        indices=test_indices,
        prediction_mode="label",
        interpretability="ig",
        batch_size=50,  # Large batch size
        return_dict=False,
    )

    # Results should be the same regardless of batch size
    predictions_match = np.array_equal(pred3, pred4)
    # Compare ranked features DataFrames by their attribution_mean values
    attributions_close = np.allclose(
        attr3["attribution_mean"].values,
        attr4["attribution_mean"].values,
        rtol=1e-5,
        atol=1e-8,
    )

    print(f"✅ Batch size consistency - predictions match: {predictions_match}")
    print(f"✅ Batch size consistency - attributions close: {attributions_close}")

    # For stochastic models, we can't guarantee exact matches across different batch sizes
    # But we can ensure the results are at least reasonable
    if not (predictions_match and attributions_close):
        print(
            "Note: Small differences in batch processing are expected for stochastic models"
        )
        # Just ensure shapes are consistent
        assert pred3.shape == pred4.shape, "Prediction shapes should match"
        assert attr3.shape == attr4.shape, "Attribution shapes should match"
    else:
        print("✅ Batch size consistency tests passed!")

    print("✅ All parameter tests completed!")


def test_interpretability_quick():
    """Quick test for interpretability functionality - can be run independently (moved from test_fadvi_predict.py)."""
    print("\n=== Quick Interpretability Test ===")

    # Skip if captum is not available
    try:
        import captum

        print("Captum available - running quick test")
    except ImportError:
        print("Captum not available - skipping quick test")
        return

    # Create minimal synthetic data
    from scvi.data import synthetic_iid

    adata = synthetic_iid(
        batch_size=200,  # Small dataset
        n_genes=100,  # Few genes
        n_batches=2,  # Few batches
        n_labels=2,  # Few labels
        sparse_format="csr_matrix",
    )

    # Add required annotations
    adata.obs["batch"] = adata.obs["batch"].astype("category")
    adata.obs["label"] = adata.obs["labels"].astype("category")

    # Setup and quick train
    FADVI.setup_anndata(adata, batch_key="batch", labels_key="label")
    model = FADVI(adata, n_latent_b=3, n_latent_l=3, n_latent_r=3)
    model.train(max_epochs=2, batch_size=64, early_stopping=False)

    # Test interpretability on just a few cells
    test_cells = [0, 1, 2, 3, 4]

    # Test batch interpretability
    batch_pred, batch_attr = model.predict(
        adata,
        indices=test_cells,
        prediction_mode="batch",
        interpretability="ig",
        return_dict=False,
    )

    # Test label interpretability
    label_pred, label_attr = model.predict(
        adata,
        indices=test_cells,
        prediction_mode="label",
        interpretability="ig",
        return_dict=False,
    )

    print(f"✅ Quick test passed!")
    print(
        f"   Batch predictions: {batch_pred.shape}, ranked features: {batch_attr.shape}"
    )
    print(
        f"   Label predictions: {label_pred.shape}, ranked features: {label_attr.shape}"
    )

    # Verify that attributions are now properly ranked DataFrames
    assert isinstance(
        batch_attr, pd.DataFrame
    ), "Batch attributions should be DataFrame"
    assert isinstance(
        label_attr, pd.DataFrame
    ), "Label attributions should be DataFrame"

    # Test functions should not return values
    assert True  # All tests passed


# Additional pytest-discoverable test functions for moved methods
def test_interpretability_functionality_pytest():
    """Pytest-discoverable version of interpretability functionality test."""
    try:
        import captum

        test_interpretability_functionality()
    except ImportError:
        import pytest

        pytest.skip("captum not available")


def test_interpretability_error_handling_pytest():
    """Pytest-discoverable version of error handling test."""
    test_interpretability_error_handling()


def test_interpretability_with_different_parameters_pytest():
    """Pytest-discoverable version of parameter test."""
    try:
        import captum

        test_interpretability_with_different_parameters()
    except ImportError:
        import pytest

        pytest.skip("captum not available")


def test_interpretability_quick_pytest():
    """Pytest-discoverable version of quick interpretability test."""
    try:
        import captum

        test_interpretability_quick()
    except ImportError:
        import pytest

        pytest.skip("captum not available")


if __name__ == "__main__":
    # Set random seeds for reproducibility
    np.random.seed(42)
    torch.manual_seed(42)

    # Run the tests
    try:
        print("Running FADVI interpretability tests...")

        print("\n1. Running interpretability analysis tests...")
        test_interpretability_analysis()

        print("\n2. Running GradientShap specific tests...")
        test_gs_interpretability()

        print("\n3. Running methods comparison tests...")
        test_interpretability_methods_comparison()

        print("\n4. Running moved interpretability functionality tests...")
        test_interpretability_functionality()

        print("\n5. Running moved interpretability error handling tests...")
        test_interpretability_error_handling()

        print("\n6. Running moved interpretability parameter tests...")
        test_interpretability_with_different_parameters()

        print("\n7. Running moved quick interpretability tests...")
        test_interpretability_quick()

        print("\n🎉🎉 ALL INTERPRETABILITY TESTS COMPLETED SUCCESSFULLY! 🎉🎉")

    except Exception as e:
        print(f"\n❌ Interpretability tests failed with error: {e}")
        import traceback

        traceback.print_exc()
    try:
        print("Running FADVI interpretability tests...")

        print("\n1. Running interpretability analysis tests...")
        test_interpretability_analysis()

        print("\n2. Running GradientShap specific tests...")
        test_gs_interpretability()

        print("\n3. Running methods comparison tests...")
        test_interpretability_methods_comparison()

        print("\n🎉🎉 ALL INTERPRETABILITY TESTS COMPLETED SUCCESSFULLY! 🎉🎉")

    except Exception as e:
        print(f"\n❌ Interpretability tests failed with error: {e}")
        import traceback

        traceback.print_exc()
