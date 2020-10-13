import numpy as np
import matplotlib.pyplot as plt
import astropy.units as au
import astropy.constants as ac

import xarray as xr
from scipy.interpolate import interp1d
from scipy.stats import poisson

import logging

__all__ = [ "TigressWindModel" ]

class TigressWindModel(object):
    """TIGRESS Wind Launching Model class

    Parameters
    ----------
    z0 : ['H','2H','500','1000']

    Examples
    --------
    >>> from twind import *
    >>> tw = TigressWindModel(z0='H')
    >>> tw.set_axes()
    >>> pdf = tw.build_Mpdf()
    """

    def __init__(self, z0='H', verbose=False):
        self.ncomp = 2 # fixed
        self.scaling_field = 'sfr' # fixed
        self.z0list=['H','2H','500','1000']
        self._set_z0(z0)
        self._set_model_parameters(verbose=verbose)

    def _set_z0(self,z0):
        """Initialize z0"""
        if z0 in self.z0list:
            self.z0 = z0
        else:
            raise ValueError('z0 must be one of {}'.format(self.z0list))

    def _set_model_parameters(self, verbose=False):
        """Initialize PDF model parameters.

        Parameters
        ----------
        verbose : bool
            call `show_parameters()` to print parameters
        """
        from scipy.special import gamma

        z0 = self.z0

        # set parameters that are constants
        p_v, d_v, cs0, sigma, vout0 = (1, 2, 6.7, 0.1, 25.0)
        p_vB, d_vB, Mach0, p_M, d_M = (4, 2, 0.5, 1, 3)

        # calculate amplitudes that make the pdf integrate to 1
        A_v = np.log(10)*p_v/gamma(d_v/p_v)
        A_cs = np.log(10)/np.sqrt(2*np.pi)/sigma
        A_vB = np.log(10)*p_vB/gamma(d_vB/p_vB)
        A_M = np.log(10)*p_M/gamma(d_M/p_M)

        # store them in dictionaries
        self.cool_params = dict(A_v=A_v, p_v=p_v, d_v=d_v,
                                A_cs=A_cs, cs0=cs0, sigma=sigma, vout0=vout0)
        self.hot_params = dict(A_vB=A_vB, p_vB=p_vB, d_vB=d_vB,
                   A_M=A_M, Mach0=Mach0,p_M=p_M,d_M=d_M)
        # SN related parameters that set the reference values for loading factors
        self.params = dict(Esn=1.e51*au.erg, mstar=95.5*au.M_sun, vcool=200*au.km/au.s,
                           Mej=10.*au.M_sun, ZSN=0.2, ZISM0=0.02)
        self.params['vej'] = np.sqrt(2.0*self.params['Esn']/self.params['Mej']).to('km/s')
        self.ref_params = dict(Mref=self.params['mstar'],
            pref=self.params['Esn']/(2*self.params['vcool']),
            Eref=self.params['Esn'],
            Zref=self.params['Mej']*self.params['ZSN'])

        # coefficients used in conversion from mass to other PDFs
        self.vp = (self.ref_params['pref']/self.params['mstar']).to('km/s').value
        self.vE = np.sqrt(self.ref_params['Eref']/self.params['mstar']).to('km/s').value
        self.Ze = (self.ref_params['Zref']/self.params['mstar']).cgs.value

        # parameters for scaling relations from Paper~I
        a = np.array(fit_alpha[z0])
        b = np.array(fit_beta[z0])

        self.scaling_params = dict(a=a, b=b)
        if z0 == '2H':
            self.cool_params['vout0'] = 45
            self.cool_params['cs0'] = 7.5
        elif z0 == '500':
            self.cool_params['vout0'] = 45
            self.cool_params['cs0'] = 8.5
        elif z0 == '1000':
            self.cool_params['vout0'] = 60
            self.cool_params['cs0'] = 10.0
        self.scaling_params['A'] = np.round(10.**(np.array(self.scaling_params['a'])),2)
        self.scaling_params['p'] = 1.+np.array(self.scaling_params['b'])
        self.enum=dict(M_cool=0, M_int=1, M_hot=2, M_total=3,
                            p_cool=4, p_int=5, p_hot=6, p_total=7,
                            E_cool=8, E_int=9, E_hot=10, E_total=11,
                            Z_cool=12, Z_int=13, Z_hot=14, Z_total=15)

        # print parameters
        if verbose:
            self.show_parameters()

    def reset_parameters(self,z0):
        """Reset parameters for different z0

        This sets `z0` attribute and calls `_set_parameters()` method.

        Parameters
        ----------
        z0 : ['H','2H','500','1000']
        """
        self._set_z0(z0)
        self._set_parameters()

    def show_parameters(self):
        """Print all parameters in readable forms
        """
        with np.printoptions(precision=3, suppress=True):
            print('number of wind phase = {}'.format(self.ncomp))
            print('galactic parameter = {}'.format(self.scaling_field))
            print('reference height = {}'.format(self.z0))
            for p in ['cool_params','hot_params','params','ref_params','scaling_params']:
                params = getattr(self,p)
                print(p)
                for k,v in params.items():
                    print('    {} = {}'.format(k,v))

    def set_axes(self,pdf=None,sfr=(-6,2,100),vout=(0,4,500),cs=(0,4,500),verbose=False):
        """Set axes (`vout`, `cs`, `sfr`) using `xarray` for convenient manipulations

        If a simulated pdf is an input, model axes are set to be identical to
        those of the input pdf otherwise, axes are set for a given (log min,
        log max, N bins)

        Key attributes, `u` = `logvout`, `w` = `logcs`, `logsfr`, `vBz`, and
        `Mach`, will be set

        Parameters
        ----------
        pdf : xarray.Dataset or xarray.DataArray
            a joint pdf from simulation
        sfr : float, tuple, list
            a single value of SFR surface density, or
            log values of min, max, and No. of bins for SFR axis
        vout : tuple, list
            log values of min, max, and No. of bins for vout axis
        cs : tuple, list
            log values of min, max, and No. of bins for cs axis
        verbose : bool
            print input ranges
        """
        if pdf is not None:
            if verbose:
                print('Setting up from simulation PDF...')
                attrs = pdf.attrs
                u = pdf.logvout.data
                x1,x2,dbin = u.min(), u.max(), attrs['dbin']
                print('  u in ({:.1f},{:.1f}) with du = {:.2f}'.format(x1,x2,dbin))
                w = pdf.logcs.data
                x1,x2,dbin = w.min(), w.max(), attrs['dbin']
                print('  w in ({:.1f},{:.1f}) with dw = {:.2f}'.format(x1,x2,dbin))
                print('  Sigma_SFR = {:.3g},'.format(attrs['sfr']), end=' ')
                print('ZISM = {:.3g}'.format(attrs['ZISM']))
                for fl in ['Mpdf','ppdf','Epdf','Zpdf']:
                    # log value of 2.e4 and 5.5e5
                    T1,T2=(1.1854266752455402,1.9121660193614398)
                    c=pdf[fl].sel(logcs=slice(0,T1)).sum().data*dbin**2
                    i=pdf[fl].sel(logcs=slice(T1,T2)).sum().data*dbin**2
                    h=pdf[fl].sel(logcs=slice(T2,4)).sum().data*dbin**2
                    t=pdf[fl].sel().sum().data*dbin**2
                    msg = '  {:5s}:'.format(fl)
                    for ph, fph in zip(['cool','int','hot','total'],[c,i,h,t]):
                        msg += ' {}={:.3f}'.format(ph,fph)
                    print(msg)
            self.logvout = pdf.logvout
            self.logcs = pdf.logcs
            self.dlogvout = pdf.attrs['dbin']
            self.dlogcs = pdf.attrs['dbin']
            self.sfr = pdf.attrs['sfr']
            self.logsfr = np.log10(self.sfr)
            self.vout = 10.**self.logvout
            self.cs = 10.**self.logcs
            self.params['ZISM0']=pdf.attrs['ZISM']
        else:
            ranges=dict(cs=cs,vout=vout)
            if hasattr(sfr, '__len__'):
                if len(sfr) == 3:
                    ranges['sfr']=sfr
                else:
                    raise ValueError('sfr should either be an array/list/tuple of'+
                                     'three elements (log min, log max, N), '+
                                     'but len(sfr)={}'.format(len(sfr)))
            else: # scalar
                self.sfr=sfr
                if sfr>0:
                    self.logsfr=np.log10(sfr)
                else:
                    raise ValueError('sfr must be positive, but sfr={}'.format(sfr))
                if verbose: print('sfr={}'.format(sfr))

            for f in ranges:
                if len(ranges[f]) != 3:
                    raise ValueError('{} should either be array-like with '.format(f)+
                                     'three elements (log min, log max, N), '+
                                     'but len({})={}'.format(f,len(ranges[f])))

                x1,x2,N = ranges[f]
                if verbose: print('{}: min={}, max={}, N={}'.format(f,x1,x2,N))
                x = np.linspace(x1,x2,N)
                x_da = xr.DataArray(x,coords=[x],dims=['log'+f])
                setattr(self,'dlog'+f,x[1]-x[0])
                setattr(self,'log'+f,getattr(x_da,'log'+f))
                setattr(self,f,10.**getattr(self,'log'+f))

        self.u = self.logvout
        self.w = self.logcs
        self.vBz = np.sqrt(5.0*self.cs**2+self.vout**2)
        self.Mach = 1/self.cs*self.vout

    def _vout0(self,x):
        """Scaling relation for the characteristic outflow velocity of cool gas
        vs SFR surface density

        This is model specific (i.e., v0 varies with height `z0`)
        from Equation XX

        The mean value of the adopted distribution is 2*vout0

        Parameters
        ----------
        x : float, array
            SFR surface density

        Returns
        -------
        vout0 : float, array
            charateristic outflow velocity of cool gas
        """
        v0 = self.cool_params['vout0']
        return v0*x**0.23+3

    def _vB0(self,x):
        """Scaling relation for the characteristic Bernoulli velocity (outgoing
        component) of hot gas vs SFR surface density

        vBz = (vout^2+5*cs^2) with gamma = 5/3 is assumed
        This can be model specific (now it is fixed)
        from Equation XX

        The mean value of the adopted distribution is 0.69*vB0

        Parameters
        ----------
        x : float, array
            SFR surface density

        Returns
        -------
        vB0 : float, array
            characteristic Bernoulli velocity (outgoing component) of hot gas
        """
        return 2.4e3*x**0.5/(2+x**0.5)+8.e2

    def _eta_sfr_scaling(self,x,q):
        """Loading factor scaling relation

        Parameters
        ----------
        x : float, array
            SFR surface density
        q : str
            `quantity` and `phase`, e.g., `M_cool`

        Returns
        -------
        eta : float, array
            loading factor of the `quantity` and `phase`
        """
        i = self.enum[q]
        A = self.scaling_params['A'][i]
        b = self.scaling_params['b'][i]
        return A*x**b

    def _flux_sfr_scaling(self,x,q):
        """Flux scaling relation

        Parameters
        ----------
        x : float, array
            SFR surface density
        q : str
            `quantity` and `phase`, e.g., `M_cool`

        Returns
        -------
        flux : float, array
            outgoing flux of the `quantity` and `phase`
        """
        i = self.enum[q]
        A = self.scaling_params['A'][i]
        p = self.scaling_params['p'][i]
        return A*x**p

    def _etaM_cool(self,x):
        """Cool mass loading factor
        """
        return self._eta_sfr_scaling(x,'M_cool')

    def _etaM_hot(self,x):
        """Hot mass loading factor
        """
        return self._eta_sfr_scaling(x,'M_hot')

    def _etaM(self,x):
        """Total mass loading factor
        """
        return self._etaM_cool(x) + self._etaM_hot(x)

    def _etap(self,x):
        """Total momentum loading factor
        """
        return self._eta_sfr_scaling(x,'p_cool') + self._eta_sfr_scaling(x,'p_hot')

    def _etaE_cool(self,x):
        """Cool energy loading factor
        """
        return self._eta_sfr_scaling(x,'E_cool')


    def _etaE_hot(self,x):
        """Hot energy loading factor
        """
        return self._eta_sfr_scaling(x,'E_hot')


    def _etaE(self,x):
        """Total energy loading factor
        """
        return self._etaE_cool(x) + self._etaE_hot(x)

    def _etaZ(self,x):
        """Total metal loading factor
        """
        return self._eta_sfr_scaling(x,'Z_cool')+self._eta_sfr_scaling(x,'Z_hot')

    def _energy_bias(self,vBz):
        """Model of the ratio of specific energy in outgoing component to total

        Notes
        -----
        .. math::

            b(v_{\mathcal{B}, z}) \equiv \frac{v_{\mathcal{B}, z}^2}{v_{\mathcal{B}}^2}
            = \frac{v_{\rm out}^2 + 5c_s^2)}{v^2 + 5c_s^2}
            \approx = 0.1 \log (v_{\mathcal{B}, z}) + 0.6

        vBz must be smaller than 10^4 km/s
        Equation (XX) in the paper
        """
        return np.clip(0.1*np.log10(vBz)+0.6,0,1)

    def _Zmodel(self,vBz,sfr,ZISM):
        """Model of the outflow metallicity

        Notes
        -----
        .. math::

            Z(v_{\mathcal{B}, z}) \equiv \tilde{\zeta}(v_{\mathcal{B}, z})Z_{\rm ISM}
            = \left(\frac{v_{\mathcal{B}, z}}{3.2\times10^3{\rm km/s}}\right)^1.7
              (0.2 - Z_{\rm ISM}) + Z_{\rm ISM}

        Equation (XX) in the paper
        """
        yZcool = 1.0
        expo = 1.7

        vmax = self.params['vej'].to('km/s').value
        ZSN = self.params['ZSN']
        Mej = self.params['Mej']

        Zmodel = np.clip((vBz/vmax)**expo*(ZSN-ZISM)+ZISM*yZcool,None,ZSN)
        return Zmodel

    def CoolMassFluxPDF(self,u,w,sfr=1.0,params=None):
        """Model of mass loading/flux PDF for cool gas

        This utilizes generalized gamma (vout) and log-normal (cs) distributions.

        Parameters
        ----------
        u : array_like (xarray.DataArray)
            log vout
        w : array_like (xarray.DataArray)
            log cs
        sfr : float, array_like
            SFR surface density
        params : array_like or None
            (p_v, d_v, cs0, sigma)
            if None, `cool_params` attribute will be used

        Returns
        -------
        pdf : array_like (xarray.DataArray)

        Notes
        -----
        see :ref:`model`
        """

        vout = 10.**u
        cs = 10.**w

        if params is None:
            A_v = self.cool_params['A_v']
            p_v = self.cool_params['p_v']
            d_v = self.cool_params['d_v']
            A_cs = self.cool_params['A_cs']
            cs0 = self.cool_params['cs0']
            sigma = self.cool_params['sigma']
        else:
            from scipy.special import gamma
            p_v,d_v,cs0,sigma = params
            A_v = np.log(10)*p_v/gamma(d_v/p_v)
            A_cs = np.log(10)/np.sqrt(2*np.pi)/sigma

        vout0 = self._vout0(sfr)
        v=(vout/vout0)
        PDF_v = A_v*v**d_v*np.exp(-v**p_v)
        PDF_cs = A_cs*np.exp(-0.5*(np.log(cs/cs0)/sigma)**2)

        return PDF_cs*PDF_v

    def HotMassFluxPDF(self,u,w,sfr=1.0,params=None):
        """Model of mass loading/flux PDF for hot gas

        This utilizes generalized gamma distributions in vBz and Mach,
        where vBz = sqrt(vout^2 + 5*cs^2) and Mach = vout/cs

        Parameters
        ----------
        u : array_like (xarray.DataArray)
            log vout
        w : array_like (xarray.DataArray)
            log cs
        sfr : float, array_like
            SFR surface density
        params : array_like or None
            (p_v, d_v, cs0, sigma)
            if None, `cool_params` attribute will be used

        Returns
        -------
        pdf : array_like (xarray.DataArray)

        Notes
        -----
        see :ref:`model`
        """

        vout = 10.**u
        cs = 10.**w

        vB = np.sqrt(5.0*cs**2+vout**2)
        Mach = vout/cs

        if params is None:
            A_vB = self.hot_params['A_vB']
            p_vB = self.hot_params['p_vB']
            d_vB = self.hot_params['d_vB']
            A_M = self.hot_params['A_M']
            p_M = self.hot_params['p_M']
            d_M = self.hot_params['d_M']
            Mach0 = self.hot_params['Mach0']
        else:
            from scipy.special import gamma
            p_vB,d_vB,Mach0,p_M,d_M = params
            A_vB = np.log(10)*p_vB/gamma(d_vB/p_vB)
            A_M = np.log(10)*p_M/gamma(d_M/p_M)

        vB0 = self._vB0(sfr)
        v=(vB/vB0)
        PDF_v = A_vB*v**d_vB*np.exp(-v**p_vB)
        m=(Mach/Mach0)
        PDF_M = A_M*m**d_M*np.exp(-m**p_M)

        return PDF_v*PDF_M

    def build_model(self,ZISM=None,renormalize=True,energy_bias=True,verbose=False):
        """Build full PDFs for mass, momentum, energy, and metal PDFs

        This will use axes attributes (`logvout`, `logcs`, `sfr`) set by `set_axes` method

        Parameters
        ----------
        ZISM : float
            set ZISM for metal PDF (if None, ZISM=0.02)
        renormalize : bool
            if True, momentum, energy, and metal PDFs are renormalized
        energy_bias : bool
            if True, apply energy bias factor in building the energy PDF
        verbose : bool
            print integrations of both cool and hot PDFs

        Returns
        -------
        pdfs : xarray.Dataset
        """

        if (not hasattr(self,'u')) or (not hasattr(self,'w')) or (not hasattr(self,'sfr')):
            raise AttributeError("axes are not set. Call set_axes() first")

        if ZISM is None:
            ZISM = self.params['ZISM0']

        pdf_dset = self.build_Mpdf(verbose)
        self._build_ppdf(pdf_dset,renormalize)
        self._build_Epdf(pdf_dset,energy_bias,renormalize)
        self._build_Zpdf(pdf_dset,ZISM,renormalize)

        return pdf_dset

    def build_Mpdf(self,verbose=False):
        """Build mass loading/flux PDF

        This will use axes attributes (`logvout`, `logcs`, `sfr`) set by `set_axes` method

        Parameters
        ----------
        verbose : bool
            print integrations of both cool and hot PDFs

        Returns
        -------
        pdfs : xarray.Dataset
        """

        if (not hasattr(self,'u')) or (not hasattr(self,'w')) or (not hasattr(self,'sfr')):
            raise AttributeError("axes are not set. Call set_axes() first")

        # place holder for xarray dataset
        pdf_dset = xr.Dataset()

        # copy information to dataset

        pdf_dset['vBz']=self.vBz
        pdf_dset['Mach']=self.Mach
        pdf_dset.attrs['dlogcs']=self.dlogcs
        pdf_dset.attrs['dlogvout']=self.dlogvout

        dbinsq = self.dlogcs*self.dlogvout

        # Mass flux PDF
        mpdfc = self.CoolMassFluxPDF(self.u,self.w,sfr=self.sfr)
        mpdfh = self.HotMassFluxPDF(self.u,self.w,sfr=self.sfr)
        if verbose: # sanity check
            msg = 'Mass PDFs are integrated to:'
            fcool = mpdfc.sum(dim=['logcs','logvout']).min().data*dbinsq
            fhot = mpdfh.sum(dim=['logcs','logvout']).min().data*dbinsq
            msg += ' cool={:.3g}'.format(fcool)
            msg += ' hot={:.3g}'.format(fhot)
            print(msg)

        # renormalization
        etac = self._etaM_cool(self.sfr)
        etah = self._etaM_hot(self.sfr)
        etaM = self._etaM(self.sfr)
        pdf_dset['etaM-hot']=etah
        pdf_dset['etaM-cool']=etac
        pdf_dset['etaM']=etaM

        mpdfc *= etac/etaM
        mpdfh *= etah/etaM
        mpdf = mpdfc + mpdfh

        pdf_dset['Mpdf-cool']=mpdfc
        pdf_dset['Mpdf-hot']=mpdfh
        pdf_dset['Mpdf'] = mpdf

        return pdf_dset

    def _build_ppdf(self,pdf_dset,renormalize):
        """Build momentum loading/flux PDF
        """

        if (not hasattr(self,'u')) or (not hasattr(self,'w')) or (not hasattr(self,'sfr')):
            raise AttributeError("axes are not set. Call set_axes() first")

        dbinsq = self.dlogcs*self.dlogvout

        # Momentum flux PDF
        etaM = pdf_dset['etaM'] # in Msun/kpc^2/yr
        etap = self._etap(self.sfr) # in (Msun*km/s)/kpc^2/yr
        pdf_dset['etap'] = etap

        pfact = (self.vout**2+self.cs**2)/(self.vp*self.vout)
        ppdfc = etaM/etap*pdf_dset['Mpdf-cool']*pfact
        ppdfh = etaM/etap*pdf_dset['Mpdf-hot']*pfact
        ppdf = ppdfc + ppdfh

        if renormalize:
            renorm = ppdf.sum(dim=['logcs','logvout'])*dbinsq
            ppdfc = ppdfc/renorm
            ppdfh = ppdfh/renorm
            ppdf = ppdf/renorm
            pdf_dset['p_renorm'] = renorm

        pdf_dset['ppdf-cool'] = ppdfc
        pdf_dset['ppdf-hot'] = ppdfh
        pdf_dset['etap-cool'] = pdf_dset['etap']*ppdfc.sum(dim=['logcs','logvout'])*dbinsq
        pdf_dset['etap-hot'] = pdf_dset['etap']*ppdfh.sum(dim=['logcs','logvout'])*dbinsq
        pdf_dset['ppdf'] = ppdf

    def _build_Epdf(self,pdf_dset,energy_bias,renormalize):
        """Build energy loading/flux PDF
        """

        if (not hasattr(self,'u')) or (not hasattr(self,'w')) or (not hasattr(self,'sfr')):
            raise AttributeError("axes are not set. Call set_axes() first")

        dbinsq = self.dlogcs*self.dlogvout

        # Energy flux PDF
        if energy_bias:
            bias = self._energy_bias(self.vBz)
            pdf_dset.attrs['energy_bias_correction']=True
        else:
            bias = 1.0
            pdf_dset.attrs['energy_bias_correction']=False

        etaM = pdf_dset['etaM']
        etaE = self._etaE(self.sfr)
        pdf_dset['etaE'] = etaE

        Efact = 0.5*self.vBz**2/self.vE**2/bias

        epdfc = etaM/etaE*pdf_dset['Mpdf-cool']*Efact
        epdfh = etaM/etaE*pdf_dset['Mpdf-hot']*Efact
        epdf = epdfc + epdfh

        if renormalize:
            renorm = epdf.sum(dim=['logcs','logvout'])*dbinsq
            epdfc = epdfc/renorm
            epdfh = epdfh/renorm
            epdf = epdf/renorm
            pdf_dset['E_renorm'] = renorm

        pdf_dset['Epdf-cool'] = epdfc
        pdf_dset['Epdf-hot'] = epdfh
        pdf_dset['etaE-cool'] = pdf_dset['etaE']*epdfc.sum(dim=['logcs','logvout'])*dbinsq
        pdf_dset['etaE-hot'] = pdf_dset['etaE']*epdfh.sum(dim=['logcs','logvout'])*dbinsq
        pdf_dset['Epdf'] = epdf

    def _build_Zpdf(self,pdf_dset,ZISM,renormalize):
        """Build metal loading/flux PDF
        """

        if (not hasattr(self,'u')) or (not hasattr(self,'w')) or (not hasattr(self,'sfr')):
            raise AttributeError("axes are not set. Call set_axes() first")

        pdf_dset.attrs['ZISM']=ZISM

        dbinsq = self.dlogcs*self.dlogvout

        # Metal flux PDF
        etaM = pdf_dset['etaM']
        etaZ = self._etaZ(self.sfr)
        pdf_dset['etaZ']=etaZ

        Zfact = self._Zmodel(self.vBz,self.sfr,self.params['ZISM0'])/self.Ze
        Zpdfc = etaM/etaZ*pdf_dset['Mpdf-cool']*Zfact
        Zpdfh = etaM/etaZ*pdf_dset['Mpdf-hot']*Zfact
        Zpdf = Zpdfc + Zpdfh

        if renormalize:
            renorm = Zpdf.sum(dim=['logcs','logvout'])*dbinsq
            Zpdfc = Zpdfc/renorm
            Zpdfh = Zpdfh/renorm
            Zpdf = Zpdf/renorm
            pdf_dset['Z_renorm'] = renorm

        if ZISM != self.params['ZISM0']:
            Z = self._Zmodel(self.vBz,self.sfr,ZISM)
            Zplist = [Zpdfc, Zpdfh, Zpdf]
            for Zp in Zplist:
                Zp *= Z*self.Ze/Zfact

        pdf_dset['Zpdf-cool'] = Zpdfc
        pdf_dset['Zpdf-hot'] = Zpdfh
        pdf_dset['etaZ-cool'] = pdf_dset['etaZ']*Zpdfc.sum(dim=['logcs','logvout'])*dbinsq
        pdf_dset['etaZ-hot'] = pdf_dset['etaZ']*Zpdfh.sum(dim=['logcs','logvout'])*dbinsq
        pdf_dset['Zpdf'] = Zpdf

# hard coded fitting results from Kim et al. 2020
# https://ui.adsabs.harvard.edu/abs/2020arXiv200616315K/abstract
# linear fits to log Sigma_SFR,40 and log eta_q
# intercept:
#   column -- (cool, intermediate, hot, whole) phases
#   row -- eta_M, eta_p, eta_E, eta_Z
fit_alpha = {'H':[-0.0671392 , -1.21580154, -0.85672715,  0.0131208 ,
                  -1.42578049, -2.15089666, -1.01273815, -0.87942345,
                  -2.22758553, -2.62695091, -0.69832133, -0.6946117 ,
                   0.17085313, -0.94227817, -0.37524992, 0.28085058 ],
             '2H':[-0.48178938, -1.39292468, -1.06803198, -0.35677985,
                   -1.6584143 , -2.28362473, -1.32440513, -1.14772651,
                   -2.36198454, -2.76012799, -1.10862539, -1.07917999,
                   -0.23027555, -1.10167822, -0.63811413, -0.0635826 ],
             '500':[-0.3758737 , -1.3527475 , -0.99574129, -0.26085479,
                    -1.60938683, -2.26403553, -1.19451471, -1.05262572,
                    -2.35185742, -2.73542988, -0.92237569, -0.90241354,
                    -0.13504726, -1.06994167, -0.53412128, 0.03099716 ],
             '1000':[-0.73548645, -1.56370501, -1.21368564, -0.59947458,
                     -1.83809179, -2.41726079, -1.52892506, -1.35397334,
                     -2.52731165, -2.88037334, -1.3744673 , -1.33644891,
                     -0.47739876, -1.26198951, -0.79526589, -0.30588739]}

# slope:
#   column -- (cool, intermediate, hot, whole) phases
#   row -- eta_M, eta_p, eta_E, eta_Z
fit_beta = {'H':[-0.44065478, -0.22666131, -0.06933898, -0.41838774,
                 -0.28653679, -0.14879847,  0.01576121, -0.15650136,
                 -0.11676794, -0.07649367,  0.13593299,  0.10784259,
                 -0.36418699, -0.14614239,  0.03989157, -0.33538503],
            '2H':[-0.40946307, -0.23313865, -0.05772923, -0.37688503,
                  -0.23729182, -0.15618482,  0.0139299 , -0.13221861,
                  -0.06530682, -0.08637603,  0.13041081,  0.10240284,
                  -0.33124967, -0.14576976,  0.04706089, -0.28998012],
            '500':[-0.50891164, -0.30985307, -0.14488489, -0.4789282 ,
                   -0.32633961, -0.22712803, -0.07896586, -0.21185669,
                   -0.14646706, -0.15378422,  0.02831253,  0.01220397,
                   -0.43370472, -0.2242491 , -0.03971429, -0.38875551],
            '1000':[-0.51368956, -0.29764828, -0.12523562, -0.4842106 ,
                    -0.31296749, -0.20059778, -0.07786858, -0.22549327,
                    -0.12941253, -0.11355882,  0.02750263,  0.00155861,
                    -0.4365193 , -0.20517516, -0.02569534, -0.39478495]}
