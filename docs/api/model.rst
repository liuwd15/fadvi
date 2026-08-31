Module Implementation
=====================

The underlying :class:`~fadvi.FADVAE` module implements the three-subspace
variational autoencoder, the auxiliary classifier heads, the gradient-reversal
layer and the cross-covariance penalty. Most users interact with
:class:`~fadvi.FADVI` instead; this page is a reference for the internals.

See :doc:`../implementation` for a narrative description of the loss terms and
training procedure.

.. autoclass:: fadvi.FADVAE
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

Gradient reversal
-----------------

.. autofunction:: fadvi._fadvae.gradient_reversal

.. autoclass:: fadvi._fadvae.GradientReversalFunction
   :members:
