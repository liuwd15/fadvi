Quick Start Guide
==================================================

This tutorial will get you started with FADVI in just a few minutes.

Installation
--------------------------------------------------

First, install FADVI:

.. code-block:: bash

   pip install fadvi

Basic Usage
---------------------------------------------------

1. **Import libraries**

.. code-block:: python

   import fadvi
   import scanpy as sc

2. **Load data**

.. code-block:: python

   # Load your data
   adata = sc.read_h5ad("your_data.h5ad")

3. **Initialize and train the model**

.. code-block:: python

    fadvi.FADVI.setup_anndata(adata,
       batch_key="batch",
       labels_key="cell_type",
       unlabeled_category="Unknown",
       layer="counts"
   )

   # Create FADVI model
   model = fadvi.FADVI(adata)
   
   # Train the model
   model.train(max_epochs=30)

4. **Get results**

.. code-block:: python

   # Get latent representation
   latent = model.get_latent_representation()
   adata.obsm["X_fadvi_l"] = latent


Next Steps
-----------------------------------------------

* Learn more about the :doc:`basic_usage` workflow
* Explore :doc:`advanced_usage` features
* Check out the :doc:`../api/index` for detailed parameter descriptions
