"""
Example usage of the FADVI model for disentangling batch and label effects.
"""

import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc

# Add the src directory to the path for local development
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import scvi
from anndata import AnnData

# Import FADVI directly from the installed package
try:
    from fadvi import FADVI

    print("Successfully imported FADVI from installed package")
except ImportError:
    # Fallback for development
    import os
    import sys

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
    from fadvi import FADVI


# Simulate some simple data for testing
def create_test_data(n_cells=1000, n_genes=2000, n_batches=3, n_labels=4):
    """Create synthetic data for testing FADVI."""
    # Generate random gene expression data
    X = np.random.negative_binomial(20, 0.3, size=(n_cells, n_genes))

    # Create batch and label assignments
    batches = np.random.choice(n_batches, size=n_cells)
    labels = np.random.choice(n_labels, size=n_cells)

    # Create AnnData object
    adata = AnnData(X=X.astype(np.float32))
    adata.obs["batch"] = pd.Categorical(batches.astype(str))
    adata.obs["labels"] = pd.Categorical(labels.astype(str))

    # Add some gene names
    adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
    adata.obs_names = [f"Cell_{i}" for i in range(n_cells)]

    return adata


def test_fadvi():
    """Test the FADVI model implementation."""
    print("Creating synthetic data...")
    adata = create_test_data()

    print("Setting up FADVI model...")
    FADVI.setup_anndata(adata, batch_key="batch", labels_key="labels")

    print("Initializing FADVI model...")
    model = FADVI(adata)

    print("Model summary:")
    print(model)

    # Debug: Check what's in the registry
    print("Registry keys:", list(model.adata_manager.registry.keys()))
    if "labels" in model.adata_manager.registry:
        print("Labels mapping:", model.adata_manager.registry["labels"])

    print("Training model (small number of epochs for testing)...")
    model.train(max_epochs=5, batch_size=128)

    print("Getting latent representations...")
    # Get different latent representations
    latent_full = model.get_latent_representation(representation="full")
    latent_b = model.get_latent_representation(representation="b")
    latent_l = model.get_latent_representation(representation="l")
    latent_r = model.get_latent_representation(representation="r")

    print(f"Full latent shape: {latent_full.shape}")
    print(f"Batch latent (b) shape: {latent_b.shape}")
    print(f"Label latent (l) shape: {latent_l.shape}")
    print(f"Residual latent (r) shape: {latent_r.shape}")

    # Add to AnnData object
    adata.obsm["X_fadvi_full"] = latent_full
    adata.obsm["X_fadvi_b"] = latent_b
    adata.obsm["X_fadvi_l"] = latent_l
    adata.obsm["X_fadvi_r"] = latent_r

    # Verify the test completed successfully with assertions
    assert (
        latent_full.shape[0] == adata.n_obs
    ), "Full latent representation has wrong number of cells"
    assert (
        latent_b.shape[0] == adata.n_obs
    ), "Batch latent representation has wrong number of cells"
    assert (
        latent_l.shape[0] == adata.n_obs
    ), "Label latent representation has wrong number of cells"
    assert (
        latent_r.shape[0] == adata.n_obs
    ), "Residual latent representation has wrong number of cells"

    print("FADVI model test completed successfully!")
    # Note: pytest test functions should not return values


if __name__ == "__main__":
    test_fadvi()
