**Note:** This tutorial was generated from an IPython notebook that can be downloaded
`here <https://github.com/changgoo/Twind/tree/master/docs/_static/notebooks/model_table.ipynb>`_.

.. _model_table:



Simulation Model Tables
=======================

We construct model PDFs based on the results from the TIGRESS simulation
suite presented in `Paper
I <https://ui.adsabs.harvard.edu/abs/2020ApJ...900...61K/abstract>`__.
We summarize model parameters and some integrated outflow propertes here
using table files made available at
`zenodo <https://doi.org/10.5281/zenodo.3872048>`__ or
`github <https://github.com/changgoo/tigress-wind-figureset/tree/v1.0.1>`__.
Mainly, the results are phase separated (three large bins in temperature
or :math:`c_s`) but outflow velocity integrated (:math:`v_{\rm out}>0`).

You can download the original
`notebook <https://github.com/changgoo/tigress-wind-figureset/blob/paper1/tables/Example_scripts.ipynb>`__
to reproduce tables and figures in `Paper
I <https://ui.adsabs.harvard.edu/abs/2020ApJ...900...61K/abstract>`__.

Download and Prepare Tables
---------------------------

.. code:: python

    # Download Tables
    import urllib.request
    import os
    
    repo_url='https://changgoo.github.io/tigress-wind-figureset'
    tbl_files=['table-mean.ecsv','table-mean-err.ecsv']
    if not os.path.isdir('tables/'): os.mkdir('tables/')
    for f in tbl_files:
        if not os.path.isfile(f):
            urllib.request.urlretrieve('{}/tables/{}'.format(repo_url,f),'tables/'+f)

.. code:: python

    # Read Tables with astropy:
    
    from astropy.table import QTable,Table
    tmean=Table.read('tables/table-mean.ecsv')
    terr=Table.read('tables/table-mean-err.ecsv')

.. code:: python

    # add additional time scales for Table 2 in Paper I
    import astropy.constants as ac
    import astropy.units as au
    tmean['torb']=(2*np.pi/tmean['Omega_0'].quantity).to('Myr')
    tmean['tosca']=(2*np.pi/np.sqrt(4*np.pi*ac.G*tmean['rho_tot'].quantity)).to('Myr')
    tmean['toscn']=(2.0*np.pi*tmean['H'].quantity/tmean['sigma_eff'].quantity).to('Myr')
    
    
    # set format for more compact display
    for k in tmean.keys():
        if tmean[k].info.dtype == 'float64':
            tmean[k].info.format = '15.2g'
            if k in terr: terr[k].info.format = '15.2g'

Table 1: Model Parameters
-------------------------

.. code:: python

    table1_varlist=['model','Sigma_gas0','Sigma_star','rho_dm','Omega_0','z_star','R_0']
    for k in table1_varlist:
        if tmean[k].info.dtype == 'float64':
            tmean[k].info.format = '15.3g'
    
    tbl1=tmean[(tmean['z']=='H') & (tmean['phase']=='whole') ][table1_varlist]
    
    tbl1.pprint_all()


.. parsed-literal::

    model    Sigma_gas0      Sigma_star        rho_dm         Omega_0          z_star           R_0      
           solMass / pc2   solMass / pc2   solMass / pc3    km / (kpc s)         pc             kpc      
    ----- --------------- --------------- --------------- --------------- --------------- ---------------
       R2             150             450            0.08             100             245               2
       R4              50             208           0.024            53.7             245               4
       R8              12              42          0.0064              28             245               8
      R16            2.49            1.71         0.00143            11.9             245              16
     LGR2             150             110           0.015              50             500               2
     LGR4              60              50           0.005              30             500               4
     LGR8              12              10          0.0016              15             500               8


-  **Sigma_gas0**: initial gas surface density,
   :math:`\Sigma_\text{gas,0}`
-  **Sigma_star**: stellar surface density, :math:`\Sigma_{*}`
-  **rho_dm**: midplane dark matter density, :math:`\rho_\text{dm}`
-  **Omega**: angular velocity of galactic rotation, :math:`\Omega`
-  **R_0**: galactocentric radius, :math:`R_0`
-  **z_star**: scale height of stellar disk, :math:`z_*`

Table 2: Time Scales
--------------------

.. code:: python

    table2_varlist=['model','torb','toscn','tosca','tdep40','surf','sfr40']
    tbl2=tmean[(tmean['z']=='H') & (tmean['phase']=='whole') ][table2_varlist]
    
    tbl2.pprint_all()


.. parsed-literal::

    model       torb           toscn           tosca           tdep40           surf             sfr40       
                Myr             Myr             Myr             Myr        solMass / pc2  solMass / (kpc2 yr)
    ----- --------------- --------------- --------------- --------------- --------------- -------------------
       R2              61              32              23              66              74                 1.1
       R4         1.1e+02              51              38         2.4e+02              29                0.12
       R8         2.2e+02         1.2e+02              75         2.1e+03              11              0.0051
      R16         5.2e+02         4.6e+02         3.1e+02         3.1e+04             2.5               8e-05
     LGR2         1.2e+02              52              48         1.5e+02              74                0.49
     LGR4           2e+02              87              80         4.2e+02              38                0.09
     LGR8         4.1e+02         2.2e+02         1.7e+02         3.3e+03              10              0.0032


-  **torb**: orbit time, :math:`t_\text{orb}=2\pi/\Omega`
-  **toscn**: vertical oscillation time derived from numerical measures,
   :math:`t_\text{osc,n}=2\pi H/\sigma_{\rm z,eff}`
-  **tosca**: vertical oscillation time derived from input parameters,
   :math:`t_\text{osc,a}=2\pi/(4\pi G\rho_{\rm tot})^{1/2}`
-  **tdep40**: gas depletion time with SFR surface density in 40 Myr,
   :math:`t_\text{dep,40}=\Sigma_\text{gas}/\Sigma_\text{SFR,40}`
-  **surf**: mean gas surface density, :math:`\Sigma_\text{gas}`
-  **sfr40**: mean SFR surface density from star particles young than 40
   Myr, :math:`\Sigma_\text{SFR,40}`

*mean and error are determined from bootstrap resampling with a sample
size of 10 for time series over* :math:`0.5<t/t_{\rm orb}<1.5`

Table 3-1: Fluxes
-----------------

.. code:: python

    z0='H' # height can be ('H','2H','500','1000')
    table3_varlist1=['model','phase','mass','mom','energy','metal','metal_sn']
    tbl3=tmean[tmean['z']==z0][table3_varlist1]
    
    tbl3.pprint_all()


.. parsed-literal::

    model phase         mass                  mom                 energy            metal              metal_sn     
                solMass / (kpc2 yr) km solMass / (kpc2 s yr) erg / (kpc2 yr) solMass / (kpc2 yr) solMass / (kpc2 yr)
    ----- ----- ------------------- ------------------------ --------------- ------------------- -------------------
       R2  cool                0.74                       50         7.2e+46               0.029              0.0032
       R2   int               0.063                       10         2.8e+46              0.0026             0.00056
       R2   hot                0.13                  1.4e+02         2.8e+48              0.0096              0.0062
       R2 whole                0.94                    2e+02         2.9e+48               0.041                0.01
       R4  cool                0.26                       12           1e+46              0.0081             0.00042
       R4   int               0.014                      1.8         4.1e+45             0.00047             7.1e-05
       R4   hot               0.026                       18         2.2e+47              0.0013             0.00058
       R4 whole                 0.3                       32         2.3e+47              0.0098               0.001
       R8  cool               0.032                     0.78         4.4e+44             0.00071             2.1e-05
       R8   int              0.0012                     0.12         2.3e+44             2.9e-05             2.9e-06
       R8   hot              0.0013                     0.67         5.5e+45             4.1e-05             1.5e-05
       R8 whole               0.035                      1.6         6.2e+45             0.00078             3.8e-05
      R16  cool              0.0055                    0.085         2.3e+43             0.00011             2.5e-09
      R16   int             3.6e-05                   0.0028         3.7e+42             7.7e-07             5.2e-08
      R16   hot             1.4e-05                   0.0093         6.1e+43             4.4e-07             1.8e-07
      R16 whole              0.0055                    0.097         8.8e+43             0.00011             8.4e-08
     LGR2  cool                0.55                       26         2.8e+46               0.018              0.0015
     LGR2   int               0.026                      3.6         8.8e+45             0.00097             0.00019
     LGR2   hot               0.055                       48         6.8e+47              0.0033              0.0018
     LGR2 whole                0.63                       78         7.1e+47               0.023              0.0034
     LGR4  cool                0.45                       14         8.3e+45               0.012             0.00021
     LGR4   int                0.01                      1.2         2.5e+45              0.0003             3.7e-05
     LGR4   hot               0.015                       10         1.1e+47             0.00065             0.00028
     LGR4 whole                0.47                       25         1.2e+47               0.013             0.00048
     LGR8  cool                0.04                     0.86         3.6e+44             0.00087             7.9e-06
     LGR8   int             0.00074                    0.073         1.3e+44             1.7e-05             1.5e-06
     LGR8   hot             0.00089                     0.44         3.3e+45             2.7e-05             8.6e-06
     LGR8 whole               0.042                      1.4         3.8e+45             0.00092             1.8e-05


-  **mass**: mass flux, :math:`\overline{\mathcal{F}}_M`
-  **mom**: momentum flux, :math:`\overline{\mathcal{F}}_p`
-  **energy**: energy flux, :math:`\overline{\mathcal{F}}_E`
-  **metal**: metal flux, :math:`\overline{\mathcal{F}}_Z`
-  **metal_sn**: SN-origin metal flux,
   :math:`\overline{\mathcal{F}}_Z^{SN}`

*mean and error are determined from bootstrap resampling with a sample
size of 10 for time series over* :math:`0.5<t/t_{\rm orb}<1.5`

Table 3-2: Loading Factors
--------------------------

.. code:: python

    z0='H' # height can be ('H','2H','500','1000')
    table3_varlist2=['model','phase','mass_loading','mom_loading',
                     'energy_loading','metal_loading','metal_sn_loading',]
    tbl3=tmean[tmean['z']==z0][table3_varlist2]
    
    tbl3.pprint_all()


.. parsed-literal::

    model phase   mass_loading    mom_loading    energy_loading  metal_loading  metal_sn_loading
                                                                                                
    ----- ----- --------------- --------------- --------------- --------------- ----------------
       R2  cool            0.68           0.035          0.0064             1.3             0.14
       R2   int           0.058          0.0071          0.0025            0.11            0.025
       R2   hot            0.12             0.1            0.24            0.42             0.27
       R2 whole            0.86            0.14            0.25             1.8             0.44
       R4  cool             2.2           0.075           0.008             3.2             0.17
       R4   int            0.12           0.012          0.0032            0.19            0.028
       R4   hot            0.22            0.12            0.17             0.5             0.23
       R4 whole             2.5             0.2            0.18             3.9              0.4
       R8  cool             6.3            0.12          0.0081             6.6             0.19
       R8   int            0.24           0.018          0.0043            0.27            0.027
       R8   hot            0.25           0.099             0.1            0.38             0.14
       R8 whole             6.8            0.23            0.11             7.3             0.36
      R16  cool              56            0.66           0.022              54           0.0012
      R16   int            0.37           0.022          0.0036            0.38            0.025
      R16   hot            0.14           0.072            0.06            0.22            0.087
      R16 whole              56            0.75           0.086              55            0.041
     LGR2  cool             1.2           0.042          0.0056             1.9             0.15
     LGR2   int           0.054          0.0058          0.0018           0.098             0.02
     LGR2   hot            0.12           0.077            0.14            0.33             0.18
     LGR2 whole             1.3            0.13            0.14             2.3             0.35
     LGR4  cool               5            0.12          0.0088             6.2             0.11
     LGR4   int            0.11            0.01          0.0027            0.16             0.02
     LGR4   hot            0.17           0.085            0.11            0.35             0.15
     LGR4 whole             5.3            0.21            0.12             6.8             0.26
     LGR8  cool              12             0.2           0.011              13             0.12
     LGR8   int            0.23           0.017           0.004            0.26            0.022
     LGR8   hot            0.28             0.1           0.099             0.4             0.13
     LGR8 whole              13            0.32            0.11              14             0.27


-  **mass_loading**: mass loading factor, :math:`\eta_M`
-  **mom_loading**: mom loading factor, :math:`\eta_p`
-  **energy_loading**: energy loading factor, :math:`\eta_E`
-  **metal_loading**: mass loading factor, :math:`\eta_Z`
-  **metal_sn_loading**: SN-origin metal loading factor,
   :math:`\eta_Z^{SN}`

*mean and error are determined from bootstrap resampling with a sample
size of 10 for time series over* :math:`0.5<t/t_{\rm orb}<1.5`

Table 4: Velocities and Metals
------------------------------

.. code:: python

    z0='H' # height can be ('H','2H','500','1000')
    
    table4_varlist=['model','phase','vout_flux','vB','Z','enrichment','fmass_sn','fmetal_sn']
    tbl4=tmean[tmean['z']==z0][table4_varlist]
    
    tbl4.pprint_all()


.. parsed-literal::

    model phase    vout_flux           vB              Z           enrichment       fmass_sn       fmetal_sn   
                     km / s          km / s                                                                    
    ----- ----- --------------- --------------- --------------- --------------- --------------- ---------------
       R2  cool              69           1e+02           0.039             1.1           0.026            0.14
       R2   int         1.4e+02         2.1e+02           0.042             1.2           0.044            0.21
       R2   hot         5.8e+02         1.4e+03           0.072             2.1            0.23            0.63
       R2 whole         1.6e+02         5.6e+02           0.044             1.3           0.059            0.27
       R4  cool              47              67           0.032             1.1           0.011           0.068
       R4   int         1.1e+02         1.6e+02           0.034             1.1           0.023            0.13
       R4   hot         3.8e+02         8.2e+02           0.046             1.6           0.095             0.4
       R4 whole           1e+02         3.2e+02           0.034             1.1           0.024            0.14
       R8  cool              20              37           0.022               1          0.0035           0.032
       R8   int              69         1.3e+02           0.024             1.1           0.012             0.1
       R8   hot         2.4e+02           6e+02           0.031             1.4           0.054            0.34
       R8 whole              34         1.4e+02           0.023             1.1          0.0066           0.057
      R16  cool             7.9              20            0.02               1         7.3e-06         7.9e-05
      R16   int              36              95           0.022             1.1          0.0063           0.071
      R16   hot         1.3e+02         5.5e+02           0.032             1.6           0.051            0.37
      R16 whole             8.4              32            0.02               1         6.3e-05         0.00068
     LGR2  cool              44              68           0.035             1.1           0.015           0.084
     LGR2   int         1.1e+02         1.8e+02           0.039             1.2           0.036            0.19
     LGR2   hot         4.2e+02           1e+03           0.057             1.8            0.15            0.51
     LGR2 whole              92         3.4e+02           0.038             1.2           0.031            0.16
     LGR4  cool              30              45           0.028               1          0.0046           0.032
     LGR4   int              92         1.5e+02            0.03             1.1           0.018            0.12
     LGR4   hot         3.1e+02         7.4e+02           0.041             1.5            0.08            0.38
     LGR4 whole              47         1.7e+02           0.028             1.1          0.0085           0.058
     LGR8  cool              13              26           0.022               1          0.0014           0.013
     LGR8   int              50         1.2e+02           0.024             1.1           0.014            0.11
     LGR8   hot         1.6e+02         4.6e+02           0.029             1.4           0.039            0.29
     LGR8 whole              17              72           0.022               1          0.0032           0.027


-  **vout_flux**: characteristic outflow velocity,
   :math:`\overline{v}_\text{out}`
-  **vB**: Bernoulli velocity, :math:`\overline{v}_{\mathcal{B}}`
-  **Z**: outflow metallicity, :math:`\overline{Z}`
-  **enrichment**: metal enrichment factor, :math:`\zeta`
-  **fmass_sn**: fraction of SN-origin mass flux, :math:`f_M^{SN}`
-  **fmetal_sn**: fraction of SN-origin metal flux, :math:`f_Z^{SN}`

*mean and error are determined from bootstrap resampling with a sample
size of 10 for time series over* :math:`0.5<t/t_{\rm orb}<1.5`
