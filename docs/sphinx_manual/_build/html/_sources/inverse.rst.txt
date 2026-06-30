Toy Inverse Problem
===================

The script ``scripts/inverse_toy/run_toy_inverse.m`` demonstrates the
least-squares parameter estimation method used to fit the model to
EXPORTS-NA observations. The toy problem estimates two parameters ---
the surface value :math:`P_0` and the attenuation coefficient :math:`b`
--- of a simple exponential depth profile :math:`P(z) = P_0\,e^{-bz}`
from six noisy synthetic observations. The cost function is

.. math::
   J(P_0, b) = \sum_{k=1}^{6}
   \left(\frac{P_\mathrm{obs}^k - P_\mathrm{mod}^k}{\sigma_d^k}\right)^2
   + \left(\frac{P_0 - P_0^\mathrm{prior}}{\sigma_{P_0}}\right)^2
   + \left(\frac{b - b^\mathrm{prior}}{\sigma_b}\right)^2,

where :math:`\sigma_d^k` is the observation error and :math:`\sigma_{P_0}`,
:math:`\sigma_b` are the prior uncertainties (Tikhonov regularization).

The script uses ``fminsearch`` to minimize the cost function and computes an
approximate posterior covariance from the finite-difference Hessian. The
true parameter values are :math:`P_0 = 5.0` and :math:`b = 0.008`
m\ :sup:`-1`; a typical run recovers :math:`P_0 \approx 5.2` and
:math:`b \approx 0.0082` m\ :sup:`-1`, well within one posterior standard
deviation of the true values.
