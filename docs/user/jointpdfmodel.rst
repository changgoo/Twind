.. _model:

Joint PDF Model
===============

.. note::
    See `Kim et al. (2020b) <link>`_ for details.

Cool outflow model
------------------

The cool outflow (:math:`T<2\times10^4{\rm K}`) in the TIGRESS suite is well
described by a model combining `log-normal <https://en.wikipedia.org/wiki/Log-normal_distribution>`_ and `generalized gamma distribution <https://en.wikipedia.org/wiki/Generalized_gamma_distribution>`_:

.. math:: \tilde{f}_{M}^{\rm cool}(u,w) = A \left(\frac{v_{\rm out}}{v_{\rm out,0}}\right)^2
    \exp\left[-\left(\frac{v_{\rm out}}{v_{\rm out,0}}\right)\right]
    \exp\left[-\frac{1}{2}\left(\frac{\ln(c_{\rm s}/c_{\rm s,0})}{\sigma}\right)^2\right]

where :math:`A=(\ln 10)^2/(2\pi\sigma^2)^{1/2}=2.12/\sigma`.

.. math:: \frac{v_{\rm out,0}}{{\rm km/s}} = 25
    \left(\frac{\Sigma_{\rm SFR}}{M_{\odot}{\rm kpc^{-2}yr^{-1}}}\right)^{0.23}+3

Hot outflow model
-----------------

The hot outflow (:math:`T>5\times10^5{\rm K}`) in the TIGRESS suite is well
described by a model combining `two generalized gamma distributions
<https://en.wikipedia.org/wiki/Generalized_gamma_distribution>`_:

.. math:: \tilde{f}_{M}^{\rm hot}(u,w) = B \left(\frac{v_{\mathcal{B},z}}{v_{\mathcal{B},0}}\right)^2
    \exp\left[-\left(\frac{v_{\mathcal{B},z}}{v_{\mathcal{B},0}}\right)^4\right]
    \left(\frac{\mathcal{M}}{\mathcal{M}_0}\right)^3
    \exp\left[-\left(\frac{\mathcal{M}}{\mathcal{M}_0}\right)\right]

where :math:`B\equiv 2(\ln 10)^2/\Gamma(1/2)=5.98`,
:math:`v_{\mathcal{B},z}\equiv(v_{\rm out}^2+c_s^2)^{1/2}`, and
:math:`\mathcal{M}=v_{\rm out}/c_s`.

.. math:: \frac{v_{\mathcal{B},0}}{10^3{\rm km/s}} = 2.4
    \left(\frac{\Sigma_{\rm SFR,1}^{1/2}}{2+\Sigma_{\rm SFR,1}^{1/2}}\right)+0.8

where :math:`\Sigma_{\rm SFR,1}\equiv \Sigma_{\rm SFR}/(M_{\odot}{\rm kpc^{-2}yr^{-1}})`
