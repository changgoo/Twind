**Note:** This tutorial was generated from an IPython notebook that can be downloaded
`here <https://github.com/changgoo/Twind/tree/master/docs/_static/notebooks/paper_figures.ipynb>`_.

.. _paper_figures:



Figures in Paper II
===================

This tutorial shows how ``Twind`` package can be used to produce figures
in Paper II. See the source code
`twind/plotter.py <https://github.com/changgoo/Twind/tree/master/twind/plotter.py>`__
to know what is happening under the scean.

.. code:: python

    import twind
    from twind.plotter import *

Figure 1: display simulation PDFs
---------------------------------

.. code:: python

    # read in simulated PDF
    sim = twind.TigressSimLoader('R4','H')
    sim.load(download=True)
    sim.set_axes(sim.simpdf)

.. code:: python

    with plt.style.context([{'axes.grid':False}]):
        fig=plot_flux_pdfs_yZ(sim.simpdf,grid=True)



.. image:: paper_figures_files/paper_figures_6_0.png


This figure shows the joint PDFs of
:math:`u\equiv \log_{10} v_{\rm out}` and :math:`w\equiv \log_{10} c_s`
for Model R4 at \|z\| = H. (a) Mass loading PDF, (b) energy loading PDF,
and (c) enrichment factor. The red and blue dashed lines denote
temperature cuts to separate cool (:math:`T < 2 \times 10^4` K),
intermediate (:math:`2\times10^4` K :math:`< T < 5 \times 10^5` K), and
hot (:math:`5 \times 10^5` K :math:`< T`) phases. The dotted gray lines
denote loci of constant Bernoulli velocity (labeled in (a))

.. math:: v_{B,z}\equiv (v_{\rm out}^2 + 5c_s^2)^{1/2}

and Mach number (labeled in (b))

.. math:: \mathcal{M} \equiv v_{\rm out}/c_s

**Notice that both mass and energy loading PDFs are distributed widely
in large range of u and w, and there is clear correlation between the
enrichment factor** :math:`\zeta(u,w)` **and Bernoulli velocity**
:math:`v_{B,z}`.

Figure 2: reconstruct PDFs from mass loading PDF
------------------------------------------------

As the joint PDF we are after is a function of the outflow velocity
:math:`v_{\rm out}` and sound speed :math:`c_s`, the momentum, energy,
and metal loading PDFs can be reconstructed from the mass loading PDF.
We present a reconstruction procedure in Section 3 of the paper. How
good is this reconstruction procedure? See below.

.. code:: python

    sims=twind.TigressSimContainer(z0='H')

.. code:: python

    fig=flux_reconstruction(sims.sims)



.. image:: paper_figures_files/paper_figures_10_0.png


**(a)** Examples of PDFs for loading factors projected onto
:math:`\log v_{B,z}` for model R4 at :math:`|z|=H`. Thick lines show
direct measurements of all PDFs, while thin lines with the same color
(overlying the thick lines almost everywhere) show reconstructions from
the mass PDF of momentum, energy, and metal PDFs.

**(b)** The ratios of reconstructed PDFs to the original PDFs for all
models at :math:`|z|=H`. The mean ratio at a given :math:`\log v_{B,z}`
is obtained by a massflux weighted average. (Thick lines correspond to
model R4, shown to left.

Figure 3: comparison between model and simulation PDFs
------------------------------------------------------

:ref:`model` gives a mass loading PDF model for cool and hot outflows
separately. The model PDF is entirely determined with two parameters:
:math:`\Sigma_{\rm SFR}` and :math:`Z_{\rm ISM}`. How good are the
simple ``Twind`` model PDFs in comparison with the complex simulation
PDFs? Here, we show multi-axes projection to give a quantitative view.

.. code:: python

    figM=comparison_pdfs(sims.sims['R4'],q='M')
    figE=comparison_pdfs(sims.sims['R4'],q='E')



.. image:: paper_figures_files/paper_figures_13_0.png



.. image:: paper_figures_files/paper_figures_13_1.png


Comparison between simulated and model PDFs for R4: **(a)** mass loading
and **(b)** energy loading. In each row, the first column shows full
joint PDFs in logarithmic color scale
(:math:`\log\, f_{M,E}\,[\rm dex^{-2}]`) from the simulation (color) and
model (contour). The remaining four panels are histograms showing
projections onto (from left to right) **outflow velocity**
:math:`v_{\rm out}`, **sound speed** :math:`c_s`, **Bernoulli velocity**
:math:`v_{B,z}`, **and Mach number :math:`\mathcal{M}`** axes. Model
PDFs are separated into cool (blue) and hot (orange) components. The sum
of the two (yellow) matches simulated PDFs (black lines) well
(especially for dominating components).

Figure 4: Loading factor scalings
---------------------------------

What is the model prediction for outflow properties as a function of
:math:`\Sigma_{\rm SFR}`? How much mass, momentum, energy, and metals
can travel far from the launching position in a given galactic halo?
Beyond the velocity-integrated loading factor scalings presented in
`Paper
I <https://ui.adsabs.harvard.edu/abs/2020ApJ...900...61K/abstract>`__,
it is important to ask how much of outflows have specific energy
:math:`v_{B}^2/2` large enough to climb up the potential well:

.. math::

   \newcommand\vout{v_{\rm out}}
   \newcommand\vesc{v_{\rm esc}}
   \newcommand\vBz{v_{\mathcal{B},z}}
       \eta_q(\vBz>\vesc)\equiv \tilde{\eta}_q
       \int_{\vBz=\vesc}^\infty \tilde{f}_{q}(u,w)dudw,

Depending on specific questions, one can use
:math:`v_{\rm esc} \equiv \sqrt{2\Delta\Phi}` for gravitational
potential difference between any distance, e.g.,
:math:`\Delta \Phi= \Phi(R_{\rm vir}) - \Phi(H)`.

.. code:: python

    tw=twind.TigressWindModel(z0=sims.z0,verbose=False)
    tw.set_axes(verbose=False)
    modelpdf=tw.build_model(renormalize=True,energy_bias=True,verbose=False)

.. code:: python

    fig = show_loading(modelpdf,sims=sims.sims)



.. image:: paper_figures_files/paper_figures_17_0.png


Loading factors for outflows with :math:`v_{B,z}>v_{\rm esc}`. Filled
circles are directly calculated from the simulation PDFs, while solid
lines are from the model PDFs. Solid and dashed lines in (d) denote the
model loading factors for :math:`Z_{\rm ISM}=0.02` and 0, respectively.

**Overall, the model tracks the general behavior of the simulation
results.**

Beyond the result at :math:`|z|=H` presented in Paper II, ``Twind``
includes model and simulation PDFs at :math:`|z|=2H`, 500 pc, and 1 kpc.
Given its simplicity, the agreement with the simulation PDFs is
stunning!

Figure 4 at z=2H
~~~~~~~~~~~~~~~~

.. code:: python

    sims=twind.TigressSimContainer(z0='2H')
    tw=twind.TigressWindModel(z0=sims.z0,verbose=False)
    tw.set_axes(verbose=False)
    modelpdf=tw.build_model(renormalize=True,energy_bias=True,verbose=False)
    fig = show_loading(modelpdf,sims=sims.sims)
    plt.setp(fig.axes,'ylim',(1.e-3,1.e2))




.. parsed-literal::

    [0.001, 100.0, 0.001, 100.0, 0.001, 100.0, 0.001, 100.0]




.. image:: paper_figures_files/paper_figures_20_1.png


Figure 4 at z=500pc
~~~~~~~~~~~~~~~~~~~

.. code:: python

    sims=twind.TigressSimContainer(z0='500')
    tw=twind.TigressWindModel(z0=sims.z0,verbose=False)
    tw.set_axes(verbose=False)
    modelpdf=tw.build_model(renormalize=True,energy_bias=True,verbose=False)
    fig = show_loading(modelpdf,sims=sims.sims)
    plt.setp(fig.axes,'ylim',(1.e-3,1.e2))




.. parsed-literal::

    [0.001, 100.0, 0.001, 100.0, 0.001, 100.0, 0.001, 100.0]




.. image:: paper_figures_files/paper_figures_22_1.png


Figure 4 at z=1kpc
~~~~~~~~~~~~~~~~~~

.. code:: python

    sims=twind.TigressSimContainer(z0='1000')
    tw=twind.TigressWindModel(z0=sims.z0,verbose=False)
    tw.set_axes(verbose=False)
    modelpdf=tw.build_model(renormalize=True,energy_bias=True,verbose=False)
    fig = show_loading(modelpdf,sims=sims.sims)
    plt.setp(fig.axes,'ylim',(1.e-3,1.e2))




.. parsed-literal::

    [0.001, 100.0, 0.001, 100.0, 0.001, 100.0, 0.001, 100.0]




.. image:: paper_figures_files/paper_figures_24_1.png


Figure 5: sampling example
--------------------------

Can we use ``Twind`` in cosmological simulations? Yes! We are working
hard on developing subgrid model based on ``Twind`` in `the SMAUG
collaboration <https://www.simonsfoundation.org/flatiron/center-for-computational-astrophysics/galaxy-formation/smaug/>`__.
Here’s a quick demonstration of wind particle sampling. (See Appendix B
of Paper II or the `source
code <https://github.com/changgoo/Twind/tree/master/twind/sampler.py#L166>`__
for the procedure in detail.)

.. code:: python

    # read in time series
    sim = twind.TigressSimLoader('R8','H')
    sim.load(download=True,time_series=True)
    sim.set_axes(sim.simpdf)

.. code:: python

    fig = sampling_from_simulation_sfr(sim)



.. image:: paper_figures_files/paper_figures_28_0.png


Model sampling demonstration for **(a)** mass outflow rate of cool gas
and **(b)** energy outflow rate of hot gas. The simulation result (black
solid) is compared to the model for three different particle mass
choices (colored lines; see keys). The input to the model is
:math:`\Sigma_{\rm SFR}(t)` from TIGRESS simulation R8, where
SFR\ :math:`=\Sigma_{\rm SFR} L_x L_y` is shown as the grey shaded
region in (a) and the corresponding SN energy injection rate is shown as
the grey region in (b). **(c)** Distributions of cool (circles) and hot
(squares) outflow particles sampled over :math:`t =220` -
:math:`440`\ Myr from the different mass sampling cases (number of
particles drawn is shown in the legend). The simulation PDF over the
same time interval is shown in the background.
