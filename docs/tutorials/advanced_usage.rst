Advanced Usage
================================

This tutorial covers advanced features and customization options for power users.

Custom Loss Functions
--------------------------------

You can customize the loss function behavior:

.. code-block:: python

   # Initialize model with custom loss weights
   model = fadvi.FADVI(
       adata,
       batch_key="batch",
       labels_key="cell_type",
       beta=1.0,           # KL divergence weight
       lambda_b=25.0,      # Batch classification weight (default 50)
       lambda_l=75.0,      # Label classification weight (default 50)
       alpha_bl=0.5,       # Adversarial: label from batch (default 1.0)
       alpha_lb=0.5,       # Adversarial: batch from label (default 1.0)
       alpha_rb=1.5,       # Adversarial: batch from residual (default 1.0)
       alpha_rl=1.5,       # Adversarial: label from residual (default 1.0)
       gamma=0.3           # Cross-correlation penalty (default 1.0)
   )

Starting with reference data and query data
--------------------------------------------

A common use case is to start with a labeled reference dataset and an unlabeled query dataset. 
You should concatenate both datasets into a single AnnData object, and specify the appropriate keys.

.. code-block:: python

   import scanpy as sc

   data_ref = sc.read_h5ad("reference_data.h5ad")
   data_query = sc.read_h5ad("query_data.h5ad")

   print(data_ref.obs["cell_type"].unique()) # Make sure you have the required labels in the reference data
   print(data_ref.obs["batch"].unique()) # Make sure you have the required batches in the reference data

   data_query.obs["cell_type"] = "Unknown" # Assign a placeholder label for the unlabeled query data

   # If data_query has batch information, use it and make sure it"s name is same as data_ref; otherwise, assign a default
   if "batch" not in data_query.obs.columns:
       data_query.obs["batch"] = "query_batch" # Assign a batch label for the query data if not present

   # Concatenate reference and query data
   adata = data_ref.concatenate(data_query)

   # Initialize model with concatenated data
   model = fadvi.FADVI(
       adata,
       batch_key="batch",
       labels_key="cell_type",
       unlabeled_category="Unknown"  # Specify the category for unlabeled query data
   )

Integration with scvi-tools Ecosystem
---------------------------------------

FADVI is built on scvi-tools and can be used with other scvi-tools modules:

.. code-block:: python

   import scvi
   
   # Use with scVI data loaders
   scvi.data.setup_anndata(adata, batch_key="batch", labels_key="cell_type")
   
   model.train(plan_kwargs={"lr": 1e-4, "weight_decay": 1e-4})


Custom Data Splitting
------------------------------------

Control train/validation splits:

.. code-block:: python
   
   # Train with custom indices
   model.train(
       max_epochs=100,
       train_size=0.8,  # Or use indices directly if supported
       validation_size=0.2
   )



Export for Other Tools
------------------------------------

Export results for use with other analysis tools:

.. code-block:: python

   # Export anndata
   adata.write_h5ad("fadvi_results.h5ad")
   
   # Export different latent representations as CSV
   import pandas as pd
   
   # Export batch latents
   latent_b_df = pd.DataFrame(
       model.get_latent_representation(representation="b"),
       index=adata.obs.index,
       columns=[f"FADVI_batch_{i}" for i in range(model.module.n_latent_b)]
   )
   latent_b_df.to_csv("fadvi_latent_batch.csv")
   
   # Export label latents
   latent_l_df = pd.DataFrame(
       model.get_latent_representation(representation="l"),
       index=adata.obs.index,
       columns=[f"FADVI_label_{i}" for i in range(model.module.n_latent_l)]
   )
   latent_l_df.to_csv("fadvi_latent_label.csv")
   
   # Export residual latents
   latent_r_df = pd.DataFrame(
       model.get_latent_representation(representation="r"),
       index=adata.obs.index,
       columns=[f"FADVI_residual_{i}" for i in range(model.module.n_latent_r)]
   )
   latent_r_df.to_csv("fadvi_latent_residual.csv")
   
   # Export normalized expression
   normalized = model.get_normalized_expression()
   normalized_df = pd.DataFrame(
       normalized,
       index=adata.obs.index,
       columns=adata.var.index
   )
   normalized_df.to_csv("fadvi_normalized.csv") # Be cautious with large datasets


Next Steps
-----------------------

* Explore :doc:`spatial_and_single_cell` for integrating spatial transcriptomics data with single-cell data
* Contribute to the project on `GitHub <https://github.com/your-username/fadvi>`_
* Report issues or request features
* Check out the :doc:`../api/index` for complete API documentation
