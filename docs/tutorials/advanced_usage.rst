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
   
   # Export normalized expression, be cautious with large datasets
   normalized = model.get_normalized_expression()
   normalized_df = pd.DataFrame(
       normalized,
       index=adata.obs.index,
       columns=adata.var.index
   )
   normalized_df.to_csv("fadvi_normalized.csv") 


Interpretability Analysis
---------------------------

FADVI provides interpretability analysis to understand which genes contribute most to batch and label predictions.
This feature uses attribution methods from the `Captum <https://captum.ai/>`_ library to identify important features.

Prerequisites
~~~~~~~~~~~~~

Install the interpretability dependencies:

.. code-block:: bash

   pip install captum

Basic Interpretability Usage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Get feature attributions for batch predictions:

.. code-block:: python

   # Predict with interpretability for specific cells
   cell_indices = [0, 1, 2, 3, 4]  # Indices of cells to analyze
   
   # Get batch predictions with feature attributions
   predictions, attributions = model.predict(
       adata,
       indices=cell_indices,
       prediction_mode="batch",
       interpretability="ig",  # Use Integrated Gradients
       return_dict=False  # Return as tuple (predictions, attributions)
   )
   
   print(f"Batch predictions: {predictions}")
   print(f"Top 5 most important genes:")
   print(attributions.head())

Label Prediction Interpretability
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Analyze which genes are important for cell type predictions:

.. code-block:: python

   # Get label predictions with feature attributions
   label_predictions, label_attributions = model.predict(
       adata,
       indices=cell_indices,
       prediction_mode="label",
       interpretability="ig",
       return_dict=False
   )
   
   # The attributions DataFrame contains ranked features
   print("Top 10 genes for label prediction:")
   for i, (idx, row) in enumerate(label_attributions.head(10).iterrows()):
       print(f"{i+1:2d}. {row['feature']:15s} "
             f"(mean attribution: {row['attribution_mean']:8.6f}, "
             f"std: {row['attribution_std']:8.6f})")

Available Attribution Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FADVI supports two interpretability methods:

**Integrated Gradients:**

.. code-block:: python

   predictions, attributions = model.predict(
       adata,
       indices=cell_indices,
       prediction_mode="batch",
       interpretability="ig",
       return_dict=False
   )

**GradientShap:**

.. code-block:: python

   predictions, attributions = model.predict(
       adata,
       indices=cell_indices,
       prediction_mode="label", 
       interpretability="gs",
       return_dict=False
   )

Working with Attribution Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The attribution results are returned as a ranked pandas DataFrame with the following columns:

* **feature**: Gene name/identifier
* **feature_idx**: Index of the gene in the original data
* **attribution_mean**: Mean attribution across analyzed cells
* **attribution_std**: Standard deviation of attributions
* **attribution_abs_mean**: Mean absolute attribution (magnitude of importance)
* **n_cells**: Number of cells included in the analysis

.. code-block:: python

   # Analyze attribution results
   print(f"Analyzed {attributions['n_cells'].iloc[0]} cells")
   print(f"Total genes ranked: {len(attributions)}")
   
   # Get top positive and negative contributors
   top_positive = attributions.nlargest(10, 'attribution_mean')
   top_negative = attributions.nsmallest(10, 'attribution_mean')
   
   print("Top positive contributors:")
   print(top_positive[['feature', 'attribution_mean', 'attribution_std']])
   
   print("Top negative contributors:")
   print(top_negative[['feature', 'attribution_mean', 'attribution_std']])
   
   # Export results for further analysis
   attributions.to_csv("fadvi_attributions.csv")

Soft Predictions with Interpretability
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can also get probability distributions along with attributions:

.. code-block:: python

   # Get soft predictions (probabilities) with attributions
   soft_predictions, soft_attributions = model.predict(
       adata,
       indices=cell_indices,
       prediction_mode="batch",
       soft=True,  # Return probabilities instead of hard predictions
       interpretability="ig",
       return_dict=False
   )
   
   print("Prediction probabilities:")
   print(soft_predictions)  # Shape: (n_cells, n_batches)
   
   # Attributions are still computed and ranked across all cells
   print("Feature attributions (same format as hard predictions):")
   print(soft_attributions.head())

Batch Processing for Large Datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For large datasets, process cells in batches to manage memory usage:

.. code-block:: python

   # Process cells in smaller batches
   batch_size = 50
   all_attributions = []
   
   for i in range(0, len(cell_indices), batch_size):
       batch_indices = cell_indices[i:i+batch_size]
       
       _, batch_attrs = model.predict(
           adata,
           indices=batch_indices,
           prediction_mode="label",
           interpretability="ig",
           batch_size=128,  # Internal batch size for attribution computation
           return_dict=False
       )
       
       # Store results (each batch gives same genes, just different stats)
       batch_attrs['batch_id'] = i // batch_size
       all_attributions.append(batch_attrs)
   
   # Combine results if needed for analysis
   print(f"Processed {len(all_attributions)} batches")

Visualization and Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~

Combine interpretability results with visualization:

.. code-block:: python

   import matplotlib.pyplot as plt
   import seaborn as sns
   
   # Plot top contributing genes
   top_genes = attributions.head(20)
   
   plt.figure(figsize=(10, 8))
   sns.barplot(data=top_genes, y='feature', x='attribution_mean', 
               xerr=top_genes['attribution_std'])
   plt.title('Top 20 Genes Contributing to Batch Predictions')
   plt.xlabel('Mean Attribution Score')
   plt.ylabel('Gene')
   plt.tight_layout()
   plt.savefig('fadvi_top_genes.png', dpi=300, bbox_inches='tight')
   plt.show()
   
   # Correlation with highly variable genes
   if 'highly_variable' in adata.var.columns:
       hvg_genes = set(adata.var[adata.var['highly_variable']].index)
       top_attributed = set(attributions.head(100)['feature'])
       
       overlap = len(hvg_genes.intersection(top_attributed))
       print(f"Overlap between top 100 attributed genes and HVGs: {overlap}/100")


Next Steps
-----------------------

* Explore :doc:`spatial_and_single_cell` for integrating spatial transcriptomics data with single-cell data
* Contribute to the project on `GitHub <https://github.com/your-username/fadvi>`_
* Report issues or request features
* Check out the :doc:`../api/index` for complete API documentation
