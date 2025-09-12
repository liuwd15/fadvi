# FADVI: Factor Disentanglement Variational Inference

FADVI is a deep learning method for single-cell RNA sequencing analysis that disentangles batch-related variation, label-related variation, and residual variation using adversarial training and cross-correlation penalties.

Read the [documentation](https://fadvi.readthedocs.io/en/latest/) for usage and demo.

## Features

- **Factor Disentanglement**: Separates batch effects, cell type effects, and residual variation in single-cell data
- **Integration with scvi-tools**: Built on top of the scvi-tools framework for scalable analysis
- **Batch Correction**: Removes unwanted batch effects while preserving biological signal
- **Cell Type Classification**: Performs supervised learning for cell type prediction

## Installation

### Install from PyPI

```bash
pip install fadvi
```


## Quick Start

```python

import scanpy as sc
import scvi
from fadvi import FADVI

# Load your single-cell data
adata = sc.read_h5ad("your_data.h5ad")

# Setup the model
FADVI.setup_anndata(
    adata,
    batch_key="batch",
    labels_key="cell_type",
    unlabeled_category="Unknown",
    layer="counts"
)

# Create and train the model
model = FADVI(adata)
model.train(max_epochs=30)

# Get latent representations
latent_l = model.get_latent_representation(representation="label")
latent_b = model.get_latent_representation(representation="batch")

# Get label predictions
prediction_label = model.predict(prediction_mode="label")

```

## Model Architecture

FADVI uses a variational autoencoder architecture with three latent subspaces:

- **z_b**: Batch-related latent factors
- **z_l**: Label-related latent factors  
- **z_r**: Residual latent factors

The model uses adversarial training to ensure proper disentanglement between these factor subspaces.


## Citation

Comming soon...

If you use FADVI in your research, please cite:

```bibtex
@article{fadvi2025,
  title={FADVI: disentangled representation learning for robust integration of single-cell and spatial omics data},
  author={Wendao Liu},
  journal={bioRxiv},
  year={2025}
}
```
