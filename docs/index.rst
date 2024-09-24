Twind
=====

**Twind** is a Python prototype of
`TIGRESS Multiphase Galactic Wind Launching Model <https://github.com/changgoo/twind>`_. The model is built on `the TIGRESS
simulation suite <https://ui.adsabs.harvard.edu/abs/2020ApJ...900...61K/>`_
developed as part of `the SMAUG
project <https://www.simonsfoundation.org/flatiron/center-for-computational-astrophysics/galaxy-formation/smaug2/papersplash1/>`_.

Basic Usage
-----------

If you want to obtain scaling relations between wind (mass, momentum, energy,
and metal) loading factors and star formation rate surface density at different
escape velocity cuts, you would construct PDFs with **Twind**:

.. code-block:: python

    import twind

    tw = twind.TigressWindModel(z0='H',verbose=True)
    tw.set_axes(verbose=True)
    pdf = tw.build_model(renormalize=True,energy_bias=True)

This will return 3D PDFs in the :math:`(v_{\rm out}, c_s, \Sigma_{\rm SFR})` space.
A more complete example is available in the :ref:`quickstart` tutorial.

``pdf`` stores all information using `xarray <http://xarray.pydata.org/en/stable/>`_.
Then, additional manipulations are easy.
If you want to apply velocity cuts to get loading factors with a selected condition,
you would do something like:

.. code-block:: python

    dbinsq = pdf.attrs['dlogcs']*pdf.attrs['dlogvout']
    cdf_over_vB100 = pdf['Mpdf'].where(pdf['vBz']>100).sum(dim=['logcs','logvout'])*dbinsq
    etaM_over_vB100 = pdf['etaM']*cdf_over_vB100

You can get a quick and dirty plot for the mass loading factor:

.. code-block:: python

    etaM_over_vB100.plot()

A more complete example is available in the :ref:`loading_sfr` tutorial.


.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user/install
   user/jointpdfmodel

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   tutorials/quickstart
   tutorials/loading_sfr
   tutorials/sampling

.. toctree::
   :maxdepth: 1
   :caption: TIGRESS Simulation Suite

   tutorials/model_table
   tutorials/simulation_pdfs

.. toctree::
   :maxdepth: 1
   :caption: Examples

   tutorials/paper_figures


.. toctree::
    :maxdepth: 4
    :caption: API Reference

    api/twind

License & Attribution
---------------------
**Twind** is free software made available under the MIT License.

If you make use of **Twind** in your work, please cite our papers:

-   ``Kim et al. 2020a, ApJ, 900, 61`` `First results from SMAUG: Characterization of
    Multiphase Galactic Outflows from a Suite of Local Star-Forming Galactic Disk
    Simulations` [`arXiv <https://arxiv.org/abs/1202.3665>`__,
    `ADS <https://ui.adsabs.harvard.edu/abs/2020ApJ...900...61K/abstract>`__,
    `BibTeX <https://ui.adsabs.harvard.edu/abs/2020ApJ...900...61K/exportcitation>`__]
-   ``Kim et al. 2020b, ApJL, 903, 34`` `A Framework for Multiphase Galactic Wind
    Launching using TIGRESS` [`arXiv <https://arxiv.org/abs/2010.09090>`__,
    `ADS <https://ui.adsabs.harvard.edu/abs/2020ApJ...903L..34K/abstract>`__,
    `BibTex <https://ui.adsabs.harvard.edu/abs/2020ApJ...903L..34K/exportcitation>`__]
