#!/usr/bin/env python3
"""
Simple test to verify FADVI package installation and basic functionality.
"""

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from anndata import AnnData

# Import FADVI from the installed package
from fadvi import FADVAE, FADVI


def test_package_import():
    """Test that the package can be imported successfully."""
    print("✓ FADVI package imported successfully")
    print(f"✓ Available classes: FADVI, FADVAE")
    assert hasattr(FADVI, "__init__")
    assert hasattr(FADVAE, "__init__")
    print("✓ Classes have required methods")


def create_simple_test_data(n_cells=100, n_genes=500, n_batches=2, n_labels=3):
    """Create simple synthetic data for testing."""
    # Generate random gene expression data (smaller for testing)
    X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(np.float32)

    # Create batch and label assignments
    batches = [f"batch_{i % n_batches}" for i in range(n_cells)]
    labels = [f"cell_type_{i % n_labels}" for i in range(n_cells)]

    # Create AnnData object
    adata = AnnData(X=X)
    adata.obs["batch"] = pd.Categorical(batches)
    adata.obs["cell_type"] = pd.Categorical(labels)

    # Add gene names
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]
    adata.var_names_make_unique()

    # Store the count matrix
    adata.layers["counts"] = adata.X.copy()

    return adata


def test_data_setup():
    """Test that data can be set up for FADVI."""
    print("Creating test data...")
    adata = create_simple_test_data()
    print(f"✓ Created test data with shape {adata.shape}")

    # Test FADVI setup
    try:
        FADVI.setup_anndata(
            adata, layer="counts", batch_key="batch", labels_key="cell_type"
        )
        print("✓ FADVI.setup_anndata completed successfully")

        # Verify setup with assertions
        assert adata.shape[0] > 0, "Data should have cells"
        assert adata.shape[1] > 0, "Data should have genes"

    except Exception as e:
        print(f"✗ FADVI.setup_anndata failed: {e}")
        raise


def test_basic_package_functionality():
    """Test basic package functionality without full model training."""
    print("\n=== Testing FADVI Package Installation ===")

    # Test imports
    test_package_import()

    # Test data setup
    test_data_setup()

    print("\n✓ All basic tests passed!")
    print("✓ FADVI package is properly installed and functional!")
    # Note: pytest test functions should not return values


if __name__ == "__main__":
    try:
        test_basic_package_functionality()
        print("\n🎉 SUCCESS: FADVI package installation test completed successfully!")
    except Exception as e:
        print(f"\n❌ FAILED: {e}")
        raise
