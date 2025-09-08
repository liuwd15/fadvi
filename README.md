# FADVI: Factor Disentanglement Variational Inference

FADVI is a deep learning method for single-cell RNA sequencing analysis that disentangles batch-related variation, label-related variation, and residual variation using adversarial training and cross-correlation penalties.

## Features

- **Factor Disentanglement**: Separates batch effects, cell type effects, and residual variation in single-cell data
- **Integration with scvi-tools**: Built on top of the scvi-tools framework for scalable analysis
- **Batch Correction**: Removes unwanted batch effects while preserving biological signal
- **Cell Type Classification**: Performs supervised learning for cell type prediction

## Installation

### Install from source (development)

```bash
# Clone the repository
git clone https://github.com/liuwd15/fadvi.git
cd fadvi

# Install in development mode
pip install -e .
```


## Quick Start

```python
import fadvi
import scanpy as sc
import scvi

# Load your single-cell data
adata = sc.read_h5ad("your_data.h5ad")

# Setup the model
fadvi.FADVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="batch",
    labels_key="cell_type"
)

# Create and train the model
model = fadvi.FADVI(adata)
model.train(max_epochs=40)

# Get latent representations
latent_l = model.get_latent_representation(representation="label")
latent_b = model.get_latent_representation(representation="batch")

# Get label predictions
prediction_label = model.predict(prediction_mode="label")

# Perform batch correction
corrected = model.get_normalized_expression()
```

## Model Architecture

FADVI uses a variational autoencoder architecture with three latent spaces:

- **z_b**: Batch-related latent factors
- **z_l**: Label-related latent factors  
- **z_r**: Residual latent factors

The model uses adversarial training to ensure proper disentanglement between these factor spaces.

## Key Components

### FADVAE
The core variational autoencoder module that implements the disentanglement architecture.

### FADVI
The high-level model class that provides the user interface and integrates with scvi-tools.


## Citation

If you use FADVI in your research, please cite:

```bibtex
@article{fadvi2025,
  title={FADVI: Factor Disentanglement Variational Inference for Single-Cell Analysis},
  author={Wendao Liu},
  journal={bioRxiv},
  year={2025}
}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.
