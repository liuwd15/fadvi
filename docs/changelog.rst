Changelog
=====================================

All notable changes to FADVI will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

[Released]
-------------------------------------

[0.3.0] - 2026-08-30
-------------------------------------

Added
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Disentanglement metrics** (``fadvi._disentanglement_metrics``), exposed as
  methods on :class:`~fadvi.FADVI`:

  * ``evaluate_disentanglement()`` — overall quantification of subspace separation
  * ``analyze_subspace_specificity()`` — how well each subspace predicts its factor
  * ``compute_orthogonality_scores()`` — canonical-correlation based independence
  * ``analyze_residual_interpretability()`` — checks whether ``z_r`` acts as a
    signal sink. Accepts explicit ``batch_key`` and ``labels_key`` arguments and
    reports the columns actually used (``batch_key_used``, ``labels_key_used``);
    if omitted, the obs columns are resolved from a candidate list.
  * Predictability is scored with cross-validation, so it is not inflated by
    overfitting
  * Orthogonality and cross-correlation are reported as continuous values. No
    pass/fail threshold is applied to them, since neither has a
    community-standard cut-off; interpret predictability against the
    majority-class rate rather than ``1 / n_classes``.

* **Implementation documentation**: a new :doc:`implementation` page covering the
  loss terms, gradient-reversal layer and training-stability evidence,
  hyperparameter meanings and sensitivity, reproducibility guidance, and measured
  runtime/memory at 10k–300k cells.

* **Tests**: ``test_metrics.py`` covering the disentanglement metrics against
  synthetic data with known structure.

* **Top-level exports**: the disentanglement metrics are now importable directly
  from ``fadvi``, so they can be applied to latent representations from any
  source rather than only to a fitted model:
  ``compute_cross_correlation``, ``compute_mutual_information``,
  ``compute_orthogonality_score``, ``evaluate_subspace_predictability``,
  ``analyze_residual_interpretability``, ``evaluate_disentanglement_quality``
  and ``compare_disentanglement_methods``.

* **Tutorial**: a new :doc:`tutorials/disentanglement_metrics` page covering how
  to evaluate the factorisation, including why predictability must be read
  against the majority-class rate rather than ``1 / n_classes``.

Fixed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Malformed reStructuredText in the ``predict()`` docstring that produced Sphinx
  warnings and rendered incorrectly.
* Duplicate autodoc entries in the API reference.

[0.2.0] - 2025-11-04
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
