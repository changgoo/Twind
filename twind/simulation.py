import numpy as np
import matplotlib.pyplot as plt
import astropy.units as au
import astropy.constants as ac

import xarray as xr
from scipy.interpolate import interp1d
from scipy.stats import poisson

class TigressSimContainer(object):
    def __init__(self,z0='H',verbose=True):
        tigmodels=['R2','R4','R8','R16','LGR2','LGR4','LGR8']
        self.z0=z0
        simpdfs=dict()
        modelpdfs=dict()
        wm=TigressWindModel(z0=z0)
        for k in tigmodels:
            pdffile='../data/pdfs/{}-{}.nc'.format(k,z0)
            if verbose:
                print('================================================')
                print('Reading from {}'.format(pdffile))
            with xr.open_dataset(pdffile) as dset:
                simpdfs[k]=dset.load()
            wm.set_axes(simpdfs[k],verbose=verbose)
            modelpdfs[k]=wm.build_model(verbose=verbose)
        self.tigmodels = tigmodels
        self.simpdfs = simpdfs
        self.modelpdfs = modelpdfs
        self.windmodel = wm
