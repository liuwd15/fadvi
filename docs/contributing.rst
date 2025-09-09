Contributing to FADVI
======================================================

We welcome contributions to FADVI! This guide will help you get started.

Getting Started
------------------------------------------------------

Development Setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Fork the repository on GitHub
2. Clone your fork locally:

.. code-block:: bash

   git clone https://github.com/liuwd15/fadvi.git
   cd fadvi

3. Create a development environment:

.. code-block:: bash

   conda create -n fadvi-dev python=3.11
   conda activate fadvi-dev
   pip install -e ".[dev]"

4. Install pre-commit hooks:

.. code-block:: bash

   pre-commit install

Types of Contributions
------------------------------------------------------

Bug Reports
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When reporting bugs, please include:

* Your operating system and Python version
* Steps to reproduce the bug
* Expected vs actual behavior
* Error messages or stack traces
* Minimal code example

Feature Requests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For new features, please:

* Check if it already exists in the issues
* Explain the use case and motivation
* Provide examples of how it would be used
* Consider implementation complexity

Code Contributions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Bug fixes
* New features
* Performance improvements
* Documentation improvements
* Test coverage improvements

Development Workflow
------------------------------------------------------

1. **Create a branch** for your feature/fix:

.. code-block:: bash

   git checkout -b feature/new-awesome-feature

2. **Make changes** and write tests
3. **Run tests** to ensure everything works:

.. code-block:: bash

   pytest test/ -v

4. **Check code quality**:

.. code-block:: bash

   black src/
   flake8 src/

5. **Commit changes** with descriptive messages:

.. code-block:: bash

   git commit -m "Add new awesome feature"

6. **Push and create a pull request**

Code Standards
------------------------------------------------------

Style Guide
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Follow PEP 8 style guidelines
* Use Black for code formatting
* Use type hints where appropriate
* Write docstrings for all public functions

Documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Use Google-style docstrings
* Include examples in docstrings
* Update documentation for new features
* Keep README.md up to date

Testing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Write tests for new features
* Maintain test coverage above 80%
* Use pytest for testing
* Include integration tests for major features

Pull Request Process
------------------------------------------------------

1. **Update documentation** if needed
2. **Add tests** for new functionality
3. **Ensure all tests pass** locally
4. **Update changelog** with your changes
5. **Submit pull request** with clear description

Pull Request Template
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please include in your PR description:

* **What**: Brief description of changes
* **Why**: Motivation for the changes
* **How**: Implementation approach
* **Testing**: How you tested the changes
* **Breaking changes**: Any breaking changes
* **Related issues**: Link to related issues

Code Review Process
------------------------------------------------------

* All pull requests require review
* Reviewers will check:
  
  * Code quality and style
  * Test coverage
  * Documentation
  * Performance impact
  * Breaking changes

* Address reviewer feedback promptly
* Be respectful and constructive

Release Process
------------------------------------------------------

FADVI follows semantic versioning (SemVer):

* **Major** (X.0.0): Breaking changes
* **Minor** (0.X.0): New features, backwards compatible
* **Patch** (0.0.X): Bug fixes, backwards compatible

Community Guidelines
------------------------------------------------------

* Be respectful and inclusive
* Help newcomers get started
* Follow the code of conduct
* Share knowledge and best practices

Getting Help
------------------------------------------------------

* **Documentation**: Check the docs first
* **GitHub Issues**: Search existing issues
* **Discussions**: Use GitHub Discussions for questions
* **Email**: Contact maintainers directly for sensitive issues

Recognition
------------------------------------------------------

All contributors will be:

* Listed in the CONTRIBUTORS.md file
* Mentioned in release notes
* Added to the documentation credits

Thank you for contributing to FADVI! 🎉
