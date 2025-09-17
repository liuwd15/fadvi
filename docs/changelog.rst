Changelog
=====================================

All notable changes to FADVI will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

[Released]
-------------------------------------

[0.2.0] - 2025-09-15
-------------------------------------

Added
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Interpretability Analysis**: Feature attribution support using Captum library
  
  * Integrated Gradients (IG) method for feature importance analysis
  * GradientShap (GS) method as alternative attribution approach
  * Automatic feature ranking with statistical summaries
  * Support for both batch and label prediction interpretability
  * Integration with ``predict()`` method via ``interpretability`` parameter
  * Customizable attribution method parameters

* **Enhanced Prediction Interface**:
  
  * ``get_ranked_features()`` method for automatic feature ranking
  * Automatic integration of interpretability with prediction workflow
  * Flexible return formats (dict/tuple) via ``return_dict`` parameter
  * Comprehensive attribution result DataFrames with statistical summaries

* **Documentation Improvements**:
  
  * Advanced usage tutorial with interpretability examples
  * Complete code examples for attribution methods
  * Performance optimization guidelines
  * Visualization and analysis workflows
  * Error handling best practices

* **Test Suite Enhancements**:
  
  * Dedicated interpretability test module (``test_interpretability.py``)
  * Comprehensive testing for both IG and GS methods  
  * Test coverage for error handling and edge cases
  * Parameter validation and consistency testing
  * Integration tests for prediction + interpretability workflows

Improved
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Code Organization**: 
  
  * Reorganized test suite with clear separation of concerns
  * ``test_fadvi_predict.py`` now focuses solely on basic prediction testing
  * All interpretability functionality consolidated in ``test_interpretability.py``
  * Improved test discoverability with pytest-compatible naming

* **API Consistency**:
  
  * Standardized return formats across prediction methods
  * Enhanced parameter validation and error messages
  * Backward-compatible API with new optional parameters

* **Performance**:
  
  * Memory-efficient batch processing for interpretability analysis
  * Optimized attribution computation for large datasets
  * Configurable batch sizes for different hardware configurations

Technical
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Added Captum as optional dependency for interpretability features
* Enhanced type hints for interpretability-related methods
* Improved error handling with informative messages when dependencies missing
* Thread-safe attribution computation

Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* captum (optional, required for interpretability features)
* All existing dependencies maintained


[0.1.0.post1] - 2025-09-12
-------------------------------------

Added
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Release on PyPI
* Documentation hosting on ReadTheDocs

Technical
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Python 3.10+ compatibility

[Unreleased]
-------------------------------------

[0.1.0] - 2025-09-09
-------------------------------------

Added
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Initial release of FADVI
* Core FADVI model implementation
* Factor disentanglement for batch effects and biological labels
* Integration with scvi-tools ecosystem
* Comprehensive API documentation
* Tutorial guides and examples
* Test suite with >90% coverage
* Support for synthetic data generation
* GPU acceleration support

Features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **FADVI Model**: Main interface for factor disentanglement
* **FADVAE**: Underlying VAE implementation with disentanglement
* **Batch Effect Correction**: Remove technical batch effects
* **Label Preservation**: Maintain biological signal during correction
* **Flexible Architecture**: Customizable network architecture
* **Multiple Likelihoods**: Support for ZINB, NB, and Poisson likelihoods
* **Training Utilities**: Built-in training loops with early stopping
* **Evaluation Metrics**: Batch mixing and label preservation metrics

Technical
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Python 3.11+ compatibility
* scvi-tools >=1.3.0 integration
* PyTorch backend
* Comprehensive type hints
* Modular design for extensibility

Documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Complete API reference
* Quick start guide
* Basic and advanced tutorials
* Installation instructions
* Contributing guidelines
* Code examples and best practices

Testing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Unit tests for all major components
* Integration tests for full workflows
* Synthetic data generation for testing
* Continuous integration setup
* >90% test coverage

Known Issues
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* None at release

Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* scvi-tools >=1.3.0
* torch >=1.8.0
* numpy
* pandas  
* scanpy
* anndata

[0.0.1] - 2025-09-04
-------------------------------------

Added
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Initial project setup
* Basic package structure
* Core model skeleton
