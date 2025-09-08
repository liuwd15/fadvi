Basic Usage
===========

This tutorial covers the fundamental concepts and workflow of FADVI. FADVI is built on scvi-tools framework and follows similar design principles.

Understanding FADVI
-------------------

FADVI (Factor Disentanglement Variational Inference) is designed to:

* **Separate batch effects** from biological signal
* **Disentangle factors** that contribute to gene expression variation
* **Provide interpretable latent representations** for downstream analysis
* **Semi-supervised learning** with same input as scANVI

Key Concepts
~~~~~~~~~~~~

* **Batch Key**: Column in `adata.obs` that identifies technical batches
* **Labels Key**: Column in `adata.obs` that identifies biological conditions/cell types
* **Latent Space**: Low-dimensional representation learned by the model
* **Disentanglement**: Separation of batch effects from biological factors

Data Preparation
----------------

* The input anndata can either be scRNA-seq data or other modalities like scATAC-seq data (peaks, windows, gene activities).

Preparing AnnData Object
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import scanpy as sc
   import pandas as pd
   
   # Load your data
   
   adata = sc.read_h5ad("your_data.h5ad")
   
   # Ensure required metadata is present
   print("Batch column:", 'batch' in adata.obs.columns)
   print("Labels column:", 'labels' in adata.obs.columns)
   
   # Basic preprocessing
   sc.pp.filter_cells(adata, min_genes=200)
   sc.pp.filter_genes(adata, min_cells=3)
   
   # Highly variable genes (recommended)
   # If using scATAC-seq data, consider using different feature selection methods or select more peaks/windows (e.g. ~50k)
   sc.pp.highly_variable_genes(adata, n_top_genes=2000)

Model Initialization
-------------------

Basic Setup
~~~~~~~~~~~

.. code-block:: python

   import fadvi
   
   # Initialize model with default parameters
   model = fadvi.FADVI(
       adata,
       batch_key='batch',
       labels_key='cell_type'  # or whatever your label column is called
       unlabeled_category='Unknown',  # Optional: specify unlabeled category
       layer='counts'  # The input must be raw counts
   )

Advanced Configuration
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Initialize with custom parameters
   model = fadvi.FADVI(
       adata,
       batch_key='batch',
       labels_key='cell_type',
       unlabeled_category='Unknown',
       n_latent_b=30,         # Batch latent dimensions (default 30)
       n_latent_l=30,         # Label latent dimensions (default 30)
       n_latent_r=10,         # Residual latent dimensions (default 10)
       n_hidden=256,          # Hidden layer size (default 128)
       n_layers=2,            # Number of hidden layers (default 1)
       dropout_rate=0.1,      # Dropout rate (default 0.1)
       gene_likelihood='zinb' # Gene likelihood (zinb/nb/poisson)
   )

Training the Model
-----------------

Basic Training
~~~~~~~~~~~~~~

.. code-block:: python

   # Train with default settings
   model.train(max_epochs=30) # 30 epoches should be good for most datasets

Custom Training
~~~~~~~~~~~~~~~

.. code-block:: python

   # Train with custom parameters
   model.train(
       max_epochs=30,
       lr=1e-3,
       batch_size=256,
       check_val_every_n_epoch=10,
       early_stopping=True,
       early_stopping_patience=20
   )

Getting Results
--------------

Latent Representation
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Get different latent representations
   latent_b = model.get_latent_representation(representation="b")  # Batch latents
   latent_l = model.get_latent_representation(representation="l")  # Label latents  
   latent_r = model.get_latent_representation(representation="r")  # Residual latents
   
   # Get combined latent representation (default)
   latent_combined = model.get_latent_representation()  # All latents concatenated
   
   print(f"Batch latent shape: {latent_b.shape}")      # (n_cells, n_latent_b)
   print(f"Label latent shape: {latent_l.shape}")      # (n_cells, n_latent_l)
   print(f"Residual latent shape: {latent_r.shape}")   # (n_cells, n_latent_r)
   print(f"Combined latent shape: {latent_combined.shape}")  # (n_cells, n_latent_b+n_latent_l+n_latent_r)
   
   # Add to original AnnData
   adata.obsm['X_fadvi_b'] = latent_b
   adata.obsm['X_fadvi_l'] = latent_l
   adata.obsm['X_fadvi_r'] = latent_r
   adata.obsm['X_fadvi'] = latent_combined


Batch and Label Predictions
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Predict batch effects
   batch_pred = model.predict(prediction_mode="b")
   
   # Predict biological labels
   label_pred = model.predict(prediction_mode="l")
   
   # Add predictions to metadata
   adata.obs['batch_pred'] = batch_pred
   adata.obs['label_pred'] = label_pred

Downstream Analysis
------------------

Visualization
~~~~~~~~~~~~

.. code-block:: python

   import scanpy as sc
   
   # UMAP on different FADVI latent representations
   
   # Option 1: Label latents only (biological variation)
   sc.pp.neighbors(adata, use_rep='X_fadvi_l')
   sc.tl.umap(adata, key_added='X_umap_label')  

   # Option 2: Batch latents only
   sc.pp.neighbors(adata, use_rep='X_fadvi_b')
   sc.tl.umap(adata, key_added='X_umap_batch')  

   # Option 3: Residual latents (batch-corrected)
   sc.pp.neighbors(adata, use_rep='X_fadvi_r')
   sc.tl.umap(adata, key_added='X_umap_residual')

   # Option 4: Combined latent representation (all factors)
   sc.pp.neighbors(adata, use_rep='X_fadvi')
   sc.tl.umap(adata, key_added='X_umap_combined')

   # Plot results
   sc.pl.umap(adata, color=['batch', 'cell_type', 'batch_pred', 'label_pred'], basis='X_umap_combined')
   sc.pl.umap(adata, color=['batch', 'cell_type'], basis='X_umap_label', title='Label Latents')
   sc.pl.umap(adata, color=['batch', 'cell_type'], basis='X_umap_residual', title='Residual Latents')

Quality Assessment
~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Calculate batch mixing metrics
   from sklearn.metrics import adjusted_rand_score
   
   # Batch correction quality (lower is better)
   batch_ari = adjusted_rand_score(adata.obs['batch'], adata.obs['batch_pred'])
   print(f"Batch ARI: {batch_ari:.3f}")
   
   # Biological preservation (higher is better)
   label_ari = adjusted_rand_score(adata.obs['cell_type'], adata.obs['label_pred'])
   print(f"Label ARI: {label_ari:.3f}")

Saving and Loading Models
------------------------

Save Model
~~~~~~~~~~

.. code-block:: python

   # Save trained model
   model.save("fadvi_save", overwrite=True, save_anndata=True) 

Load Model
~~~~~~~~~~

.. code-block:: python

   # Load pre-trained model
   loaded_model = fadvi.FADVI.load("fadvi_save", adata)

Next Steps
----------

* Explore :doc:`advanced_usage` for more sophisticated use cases
* Check the :doc:`../api/index` for detailed parameter descriptions
* See example notebooks in our GitHub repository
