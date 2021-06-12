import numpy as np
import matplotlib.pyplot as plt
import astropy.units as au
import astropy.constants as ac
import os

import xarray as xr
from scipy.interpolate import interp1d
from scipy.stats import poisson

from .models import TigressWindModel

__all__ = [ "TigressSimContainer", "TigressSimLoader" ]

class TigressSimContainer(object):
    """Simulation PDF container for the TIGRESS simulation suite

    Load all models at a given height.

    Parameters
    ----------
    z0 : ['H','2H','500','1000']
    modelnames : ['R2','R4','R8','R16','LGR2','LGR4','LGR8']
        list of model names to load

    Examples
    --------
    >>> from twind import *
    >>> sim = TigressSimContainer(z0='H')

    """
    def __init__(self,z0='H',modelnames=['R2','R4','R8','R16','LGR2','LGR4','LGR8'],
                 basedir='./'):
        self._set_z0(z0)

        sims=dict()
        for name in modelnames:
            sim = TigressSimLoader(name,z0,basedir=basedir)
            sim.load(download=True)
            sim.set_axes(sim.simpdf)
            sim.pdf_reconstruction()
            sims[name] = sim

        self.sims = sims

    def _set_z0(self,z0):
        """Initialize z0"""
        if z0 in ['H','2H','500','1000']:
            self.z0 = z0
        else:
            raise ValueError('z0 must be one of {}'.format(self.z0list))


class TigressSimLoader(TigressWindModel):
    """Simulated PDF loader for the TIGRESS simulation suite

    Parameters
    ----------
    name : ['R2','R4','R8','R16','LGR2','LGR4','LGR8']
    z0 : ['H','2H','500','1000']

    Examples
    --------
    >>> from twind import *
    >>> sim = TigressSimContainer(model='R4',z0='H')

    """
    def __init__(self,name='R4',z0='H',basedir='./'):
        TigressWindModel.__init__(self, z0, False)

        self.modelnames=['R2','R4','R8','R16','LGR2','LGR4','LGR8']
        self._set_name(name)
        self.basedir=os.path.join(basedir,'data')

        fname = '{}-{}.nc'.format(self.name,self.z0)
        self.pdffile=os.path.join(self.basedir,'pdfs',fname)
        self.tsfile=os.path.join(self.basedir,'time_series',fname)

    def _set_name(self,name):
        """Initialize name"""
        if name in self.modelnames:
            self.name = name
        else:
            raise ValueError('name must be one of {}'.format(self.modelnames))

    def load(self,download=False,time_series=False):
        """Load simulation PDF

        Parameters
        ----------
        download : bool
            automatically call download() method if file doesn't exist

        """
        if not os.path.isfile(self.pdffile) and not download:
            print("Please download simulation data:")
            print("  call download() first or try load(download=True)")
            return
        if download:
            if not os.path.isfile(self.pdffile):
                self.download()
            if not os.path.isfile(self.tsfile) and time_series:
                self.download(time_series=True)

        with xr.open_dataset(self.pdffile) as dset:
            simpdf = dset.load()

        # rename to be consistent with model pdf
        simpdf = simpdf.rename(massflux='Mpdf',metalflux='Zpdf',
                energyflux='Epdf',momflux='ppdf').drop_vars(['poyntingflux'])

        # calculate loading factor
        ref = self.ref_params
        par = self.params
        attrs = simpdf.attrs

        # convert flux to loading
        attrs['mass_ref'] = (ref['Mref']/par['mstar']).cgs.value
        attrs['mom_ref'] = (ref['pref']/par['mstar']).to('km/s').value
        attrs['energy_ref'] = (ref['Eref']/par['mstar']).to('erg/Msun').value
        attrs['metal_ref'] = (ref['Zref']/par['mstar']).cgs.value

        for fl,loading in zip(['mass','mom','energy','metal'],
                ['etaM','etap','etaE','etaZ']):
            conv_fact = attrs[fl+'flux_unit']/attrs['NxNyNt']/attrs['sfr']/attrs[fl+'_ref']
            attrs[loading] = attrs[fl+'flux']*conv_fact

        # add metal related fields
        simpdf['Z'] = simpdf['Zpdf']/simpdf['Mpdf']*attrs['metalflux']/attrs['massflux']
        simpdf['yZ'] = simpdf['Z']/attrs['ZISM']

        # add Mach number for convenience
        simpdf['Mach'] = 10.**(-simpdf.logcs + simpdf.logvout)

        self.simpdf = simpdf

        # load time_series if asked
        if time_series:
            with xr.open_dataset(self.tsfile) as dset:
                ts = dset.load()
            self.time_series = ts

    def download(self,source='tigressdata',time_series=False):
        """Download simulation pdf data

        Parameters
        ----------
        source : ['tigressdata','dataverse','cca']

        Note
        ----
        'cca' server is not yet available as a source

        """
        if source == 'dataverse':
            file_id = dict(LGR2_1000=4063040,LGR2_2H=4063041,
                           LGR2_500=4063046,LGR2_H=4063029,
                           LGR4_1000=4063034,LGR4_2H=4063028,
                           LGR4_500=4063048,LGR4_H=4063035,
                           LGR8_1000=4063037,LGR8_2H=4063051,
                           LGR8_500=4063026,LGR8_H=4063039,
                           R16_1000=4063053,R16_2H=4063044,
                           R16_500=4063032,R16_H=4063052,
                           R2_1000=4063045,R2_2H=4063043,R2_500=4063031,R2_H=4063027,
                           R4_1000=4063042,R4_2H=4063047,R4_500=4063036,R4_H=4063038,
                           R8_1000=4063033,R8_2H=4063049,R8_500=4063050,R8_H=4063030)
            url = 'https://dataverse.harvard.edu/api/access/datafile/'
            url += '{}'.format(file_id['_'.join([self.name,self.z0])])
        elif source == 'tigressdata':
            url ='https://tigress-web.princeton.edu/~changgoo/'
            url += 'TIGRESS_example_data/wind-paper/pdfs/'
            url += '-'.join([self.name,self.z0])+'.nc'
        else:
            raise ValueError("source = ['dataverse', 'tigressdata']")

        directories = [self.basedir,
                       os.path.join(self.basedir,'pdfs'),
                       os.path.join(self.basedir,'time_series')]
        for d in directories:
            if not os.path.isdir(d):
                print('creating folder {}'.format(d))
                os.mkdir(d)

        #if exists(url):
        try:
            import urllib.request
            if time_series:
                msg = 'downloading time_series/{}-{}.nc from {}'
                msg = msg.format(self.name,self.z0,source)
                print(msg,end='... ')
                tsurl = url.replace('pdfs/','time_series/')
                urllib.request.urlretrieve(tsurl, self.tsfile)
                print('complete!')
            else:
                msg = 'downloading pdfs/{}-{}.nc from {}'
                msg = msg.format(self.name,self.z0,source)
                print(msg,end='... ')
                urllib.request.urlretrieve(url, self.pdffile)
                print('complete!')
        except:
            print('failed download from {}'.format(url))

    def pdf_reconstruction(self):
        """PDF reconstruction from mass loading PDF

        Parameters
        ----------

        """
        if not hasattr(self,'simpdf'):
            print("Please call load() function first")
            return

        pdf = self.simpdf
        vout = 10.**pdf.logvout
        cs = 10.**pdf.logcs
        attrs = pdf.attrs

        # Eq. XX
        flratio = attrs['massflux']/attrs['momflux']
        pdf['ppdf_r'] = pdf['Mpdf']*(vout+cs**2/vout)*flratio

        # Eq. XX and Eq. XX for b
        flratio = (attrs['massflux']/attrs['energyflux'])
        energy_bias = self._energy_bias(pdf['vBz'])
        pdf['Epdf_r'] = pdf['Mpdf']*0.5*pdf['vBz']**2*flratio/energy_bias

        # Eq. XX and Eq. XX for tildeZ
        flratio = (attrs['massflux']/attrs['metalflux'])
        Zmodel = self._Zmodel(pdf['vBz'],attrs['sfr'],attrs['ZISM'])
        pdf['Zpdf_r'] = pdf['Mpdf']*Zmodel*flratio
