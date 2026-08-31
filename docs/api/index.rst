API Reference
=============

This section contains the complete API reference for FADVI.

Main Classes
------------

:class:`~fadvi.FADVI`
    The user-facing model: data setup, training, latent representations,
    prediction and feature attribution. See :doc:`fadvi`.

:class:`~fadvi.FADVAE`
    The underlying variational module implementing the three subspaces, the
    auxiliary heads and the gradient-reversal layer. See :doc:`model`.

Latent subspaces
----------------

A trained :class:`~fadvi.FADVI` model exposes three latent subspaces through
``get_latent_representation(representation=...)``:

.. list-table::
   :header-rows: 1
   :widths: 12 22 66

   * - Key
     - Subspace
     - Contents
   * - ``"b"``
     - ``z_b`` (batch)
     - Technical / batch-associated variation, isolated by the batch classifier
       and adversarial terms.
   * - ``"l"``
     - ``z_l`` (biological)
     - Cell-type / label-associated biological variation. **This is the
       representation used for integration and for UMAP visualisation.**
   * - ``"r"``
     - ``z_r`` (residual)
     - Remaining variation not attributed to batch or label. Unsupervised, but
       constrained to be decorrelated from the other two.
   * - ``"lr"``
     - ``z_l`` + ``z_r``
     - Concatenation, used in fully unsupervised mode where label and residual
       cannot be separated.

Detailed pages
--------------

.. toctree::
   :maxdepth: 2

   fadvi
   model
   metrics
