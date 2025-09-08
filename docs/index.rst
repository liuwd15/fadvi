FADVI: Factor Disentanglement Variational Inference
==================================================

FADVI (Factor Disentanglement Variational Inference) is a deep learning framework for disentangling batch effects and biological factors in single-cell RNA sequencing data. Built on top of scvi-tools, FADVI provides a robust solution for analyzing complex single-cell datasets.

Features
--------

* **Batch Effect Correction**: Effectively removes technical batch effects while preserving biological signal
* **Factor Disentanglement**: Separates biological factors from technical confounders
* **scvi-tools Integration**: Built on the proven scvi-tools framework for reliability and performance
* **Easy to Use**: Simple API for both beginners and advanced users

Quick Start
-----------

Install FADVI using pip:

.. code-block:: bash

   pip install fadvi

Basic usage:

.. code-block:: python

   import fadvi
   import scanpy as sc
   
   # Load your data
   adata = sc.read_h5ad("your_data.h5ad")
   
   # Setup data registration
   fadvi.FADVI.setup_anndata(adata, batch_key="batch", labels_key="cell_type", unlabeled_category="Unknown")
   
   # Initialize and train FADVI model
   model = fadvi.FADVI(adata)
   model.train(max_epochs=30)

   # Get different latent representations
   adata.obsm["X_fadvi_batch"] = model.get_latent_representation(representation="b")
   adata.obsm["X_fadvi_label"] = model.get_latent_representation(representation="l") 
   adata.obsm["X_fadvi_residual"] = model.get_latent_representation(representation="r")

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   installation
   tutorials/index
   examples/index
   api/index
   contributing
   changelog

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

