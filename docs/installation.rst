Installation
=============================================

Requirements
---------------------------------------------

FADVI requires Python 3.11 or higher and the following dependencies:

* scvi-tools >= 1.3.0
* PyTorch >= 1.8.0
* NumPy
* Pandas
* scanpy
* anndata

Installation Methods
---------------------------------------------

From PyPI (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pip install fadvi

From Source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To install the latest development version:

.. code-block:: bash

   git clone https://github.com/liuwd15/fadvi.git
   cd fadvi
   pip install -e .

Verify Installation
---------------------------------------------

To verify that FADVI is installed correctly:

.. code-block:: python

   import fadvi
   print(fadvi.__version__)

If you encounter any issues, please open an issue on our `GitHub repository <https://github.com/liuwd15/fadvi>`_.

GPU Support
---------------------------------------------

FADVI automatically detects and uses GPU acceleration when available. To check if GPU is being used:

.. code-block:: python

   import torch
   print(f"CUDA available: {torch.cuda.is_available()}")
   print(f"Number of GPUs: {torch.cuda.device_count()}")

For optimal performance with large datasets, we recommend using a GPU with at least 8GB of memory.
