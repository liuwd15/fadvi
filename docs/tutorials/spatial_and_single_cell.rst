Integrating Spatial Transcriptomics and Single-Cell Data
===============================================================

This tutorial covers the best practices for integrating spatial transcriptomics (ST) and single-cell RNA-seq (scRNA-seq) data using FADVI.

What input ST data should be used?
---------------------------------------------------------------

FADVI is **NOT** designed to deconvolve cell types from large-spot-based ST data. 
Instead, it leverages the high-resolution, single-cell level ST data with matched scRNA-seq data.
FADVI will output joint representations of ST and scRNA-seq data, which can be used for downstream analysis like clustering and visualization.

* **Supported technologies**: 10X Genomics Visium HD, 10X Genomics Xenium, CosMx, Stereo-seq and other high-resolution ST platforms
* **Unsupported technologies**: 10X Genomics Visium, and other large-spot-based ST data

Can I use unlabeled ST data?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Yes**: FADVI can integrate unlabeled ST data with labeled scRNA-seq data.
* **However**, labeled ST data is recommended for optimal integration performance.

(Recommended) Integrating labeled ST data with scRNA-seq data
---------------------------------------------------------------

The first step is to annotate cell type labels in the ST data using matched scRNA-seq data. 
Any methods for label transfer can be employed. `TACCO <https://github.com/simonwm/tacco>`_ is demonstrated in this tutorial.

Annotating ST data with TACCO
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import scanpy as sc
   import anndata as ad
   import tacco as tc
   
   # Load your data
   
   adata_sc = sc.read_h5ad("scRNAseq_data.h5ad")
   adata_st = sc.read_h5ad("st_data.h5ad")
   
   # TACCO annotation, the annotated label will be in adata_st.obs["cell_type"]
   tc.tl.annotate(adata_st, adata_sc, annotation_key="cell_type", result_key="cell_type")

   # concatenate ST and scRNA-seq data
   adata = ad.concat({"scRNA-seq": adata_sc, "spatial": adata_st}, label="tech")


Then ST and scRNA-seq data can be integrated using FADVI to obtain a joint representation.

FADVI integration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import fadvi

   # Set up AnnData
   fadvi.FADVI.setup_anndata(adata,
       batch_key="tech",
       labels_key="cell_type",
       layer="counts"
   )

   # Initialize model with default parameters
   model = fadvi.FADVI(adata)

   # Train with default settings
   model.train(max_epochs=30) # 30 epoches should be good for most datasets

   # Get latent representation
   latent = model.get_latent_representation()
   adata.obsm["X_fadvi_l"] = latent

(Alternative) Integrating unlabeled ST data with scRNA-seq data
---------------------------------------------------------------

It is possible to directly integrate unlabeled ST data with labeled scRNA-seq data, but sometimes the results may be suboptimal. This is due to the large discrepancy in the transcription profiles of same cell types in different technologies.


.. code-block:: python

   import scanpy as sc
   import anndata as ad
   import tacco as tc
   import fadvi

   # Load your data
   
   adata_sc = sc.read_h5ad("scRNAseq_data.h5ad")
   adata_st = sc.read_h5ad("st_data.h5ad")

   # Assign "Unknown" label for all ST data
   adata_st.obs["cell_type"] = "Unknown"

   # concatenate ST and scRNA-seq data
   adata = ad.concat({"scRNA-seq": adata_sc, "spatial": adata_st}, label="tech")

   # Set up AnnData
   fadvi.FADVI.setup_anndata(adata,
       batch_key="tech",
       labels_key="cell_type",
       unlabeled_category="Unknown",
       layer="counts"
   )

   # Initialize model with default parameters
   model = fadvi.FADVI(adata)

   # Train with default settings
   model.train(max_epochs=30) # 30 epoches should be good for most datasets

   # Get latent representation
   latent = model.get_latent_representation()
   adata.obsm["X_fadvi_l"] = latent


Then ST and scRNA-seq data can be integrated using FADVI to obtain a joint representation.

Next Steps
------------------------------------------

* Explore :doc:`advanced_usage` for more sophisticated use cases
* Check the :doc:`../api/index` for detailed parameter descriptions
* See example notebooks
