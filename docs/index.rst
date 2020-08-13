Twind
=====

**Twind** is a Python prototyping of `TIGRESS Multiphase Wind Launching
Model <https://github.com/changgoo/twind>`_. The model is built on `the TIGRESS
simulation suite <https://arxiv.org/abs/2006.16315>`_ developed as part
of `the SMAUG
project <https://www.simonsfoundation.org/flatiron/center-for-computational-astrophysics/galaxy-formation/smaug2/papersplash1/>`_.

Basic Usage
-----------

If you wanted to obtain scaling relations between wind (mass, momentum, energy,
and metal) loading factors and star formation rate surface density at different
escape velocity cuts, you would do something like:

.. code-block:: python

    import twind

    tw=twind.TigressWindModel(z0='H',verbose=True)
    tw.set_axes(verbose=True)
    pdf=tw.build_model(renormalize=True,energy_bias=True)

This will return 3D PDFs of :math:`(v_{\rm out}, c_s, \Sigma_{\rm SFR})` in
the `xarray.Dataset <http://xarray.pydata.org/en/stable/>`_ format.
A more complete example is available in the :ref:`quickstart` tutorial.


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
    :maxdepth: 4
    :caption: API Reference

    api/twind

License & Attribution
---------------------
**Twind** is free software made available under the MIT License.

If you make use of **Twind** in your work, please cite our papers
Kim et al. (2020a) [`arXiv <https://arxiv.org/abs/1202.3665>`_,
`ADS <https://ui.adsabs.harvard.edu/abs/2020arXiv200616315K>`_,
`BibTeX <https://ui.adsabs.harvard.edu/abs/2020arXiv200616315K/exportcitation>`_]
and Kim et al. (2020b) [links].
