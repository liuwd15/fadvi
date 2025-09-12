Changelog
=====================================

All notable changes to FADVI will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

[Released]
-------------------------------------

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
