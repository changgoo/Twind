import numpy as np
import matplotlib.pyplot as plt
import astropy.units as au
import astropy.constants as ac

import xarray as xr
from scipy.interpolate import interp1d
from scipy.stats import poisson

from .models import TigressWindModel

__all__ = [ "TigressWindSampler", "to_time_series"]

@np.vectorize
def GGD(x,d=2,p=1):
    """Two parameter generalized gamma distribution (GGD)

    Parameters
    ----------
    x : array_like (positive)
    d : float (positive)
    p : float (positive)

    Returns
    -------
    pdf : array_like

    Notes
    -----
    .. math::

        G(x;d,p) = \frac{p}{\Gamma(d/p)}x^{d-1}\exp{-x^p}

    where Gamma() is the gamma function
    """

    from scipy.special import gamma
    return p/gamma(d/p)*x**(d-1)*np.exp(-x**p)

def GGD_CDF(d=2,p=1,log=False):
    """Tabulate cumulative distribution function (CDF) of a GGD

    Parameters
    ----------
    d : float (positive)
    p : float (positive)
    log : bool
        if True, CDF is tabulated with uniform inteval in log y

    Returns
    -------
    y : array_like
        range over which the CDF is calculated
    cdf : array_like

    Notes
    -----
    .. math::

        CDF(y) = \int_0^y G(x;d,p) dx

    where G(x;d,p) is a GGD
    """

    if log:
        dlogy = 0.01
        logy = np.arange(-4,2,dlogy)
        y0 = 10.**logy
        pdf = np.log(10)*y0*GGD(y0,d=d,p=p)
        cdf = pdf.cumsum()*dlogy
    else:
        dy = 0.01
        y0 = np.arange(0,20,dy)
        pdf = GGD(y0,d=d,p=p)
        cdf = pdf.cumsum()*dy

    return y0,cdf

class TigressWindSampler(TigressWindModel):
    """Particle sampler for the TIGRESS Wind Model

    Parameters
    ----------
    z0 : ['H','2H','500','1000']

    Examples
    --------
    >>> from twind import *
    >>> sampler = TigressWindSampler()
    >>> cool,hot=sampler.draw_mass(sfr0,mcool,mhot,area=area,dt=dt)

    """
    def __init__(self, z0='H', verbose=False):
        TigressWindModel.__init__(self, z0, verbose)

        # a conversion factor between (erg/Msun) and (km/s)^2
        self.vEsq=(1.0*au.erg/au.M_sun).to('km^2/s^2').value

        # for vout cool
        p = self.cool_params
        y0,vout_cdf = GGD_CDF(d=p['d_v'],p=p['p_v'])
        self.vout_cdf = vout_cdf

        # for cs cool
        self.cs0 = p['cs0']
        self.sigma = p['sigma']

        # for vB hot
        p = self.hot_params
        y0,vB_cdf = GGD_CDF(d=p['d_vB'],p=p['p_vB'])
        self.vB_cdf = vB_cdf

        # for Mach hot
        y0,Mach_cdf = GGD_CDF(d=p['d_M'],p=p['p_M'])
        self.Mach0 = p['Mach0']
        self.Mach_cdf = Mach_cdf
        self.y0 = y0

        # set some constants for convenience
        self.ZISM = self.params['ZISM0']
        self.mstar = self.params['mstar'].to('Msun').value
        self.Eref = self.ref_params['Eref'].to('erg').value


    def get_refs(self,sfr):
        """Obtain reference rates and loading factors
        for a given SFR surface density using scaling relations

        Parameters
        ----------
        sfr : array_like
            SFR surface density

        Returns
        -------
        refs : array_like
            reference mass, momemtum, energy, metal outflow rates
        eta : array_like
            mass, momemtum, energy, metal loading factors for total gas
        eta_cool : array_like
            mass, momemtum, energy, metal loading factors for cool gas
        eta_hot : array_like
            mass, momemtum, energy, metal loading factors for hot gas
        """
        snr=sfr/self.params['mstar'].to('Msun').value
        mr=snr*self.ref_params['Mref'].to('Msun').value
        pr=snr*self.ref_params['pref'].to('Msun*km/s').value
        er=snr*self.ref_params['Eref'].to('erg').value
        Zr=snr*self.ref_params['Zref'].to('Msun').value
        refs=[mr,pr,er,Zr]

        etas=[]
        for name in ['M_total','p_total','E_total','Z_total']:
            etas.append(self._eta_sfr_scaling(sfr,name))

        etas_cool=[]
        for name in ['M_cool','p_cool','E_cool','Z_cool']:
            etas_cool.append(self._eta_sfr_scaling(sfr,name))

        etas_hot=[]
        for name in ['M_hot','p_hot','E_hot','Z_hot']:
            etas_hot.append(self._eta_sfr_scaling(sfr,name))

        return refs, etas, etas_cool, etas_hot

    def draw_mass(self,sfr,mcool,mhot,area=1.0,dt=1.e3):
        """Draw particles with fixed particle mass quanta

        Parameters
        ----------
        sfr : float, array_like
            SFR surface density in Msun/yr/kpc^2
        mcool : float
            Mass of cool gas in Msun
        mhot : float
            Mass of hot gas in Msun
        area : float
            area in kpc^2
        dt : float, array_like
            time interval in yr over which particle is sampled

        Returns
        -------
        cool, hot : dicts
            dicts containg particle mass, 3 component velocity, sound speed, metallicity,
            and index of each particle in the corresponding input SFR surface density array,
            which will be used for reconstuction of time series
        """

        # Step 0: preparation
        sfr_ = np.atleast_1d(sfr)
        dt_ = np.atleast_1d(dt)
        mstar_ = sfr_*dt_*area

        # Step 1: obatin the mass of the wind in each gas phase
        mcool_out = self._etaM_cool(sfr)*mstar_
        mcool_out[mstar_ == 0] = 0.
        mhot_out = self._etaM_hot(sfr)*mstar_
        mhot_out[mstar_ == 0] = 0.

        # Step 2: draw an integer random variate for number of particles
        # expected number of particles
        ncool_ = mcool_out/mcool
        nhot_ = mhot_out/mhot

        # Step 3-6:
        cool, hot = self._sample_particles(ncool_,nhot_,sfr_)

        # Make mass as array
        cool['mass'] = mcool*np.ones_like(cool['vz'])
        hot['mass'] = mhot*np.ones_like(hot['vz'])

        return cool,hot


    def draw_energy(self,sfr,ecool,ehot,area=1.0,dt=1.e3):
        """Draw particles with fixed particle energy quanta

        Parameters
        ----------
        sfr : float, array_like
            SFR surface density in Msun/yr/kpc^2
        ecool : float
            energy of cool gas in 10^51 erg
        ehot : float
            energy of hot gas in 10^51 erg
        area : float
            area in kpc^2
        dt : float, array_like
            time interval over which particle is sampled

        Returns
        -------
        cool, hot : dicts
            dicts containg particle mass, 3 component velocity, sound speed, metallicity,
            and index of each particle in the corresponding input SFR surface density array,
            which will be used for reconstuction of time series
        """

        # Step 0: preparation
        sfr_ = np.atleast_1d(sfr)
        dt_ = np.atleast_1d(dt)
        nsn_ = sfr_*dt_*area/self.mstar
        Einj_ = nsn_*self.Eref/1.e51

        # Step 1: obatin the energy of the wind in each gas phase
        ecool_out = self._etaE_cool(sfr)*Einj_
        ehot_out = self._etaE_hot(sfr)*Einj_

        # Step 2: draw an integer random variate for number of particles
        # expected number of particles
        ncool_ = ecool_out/ecool
        nhot_ = ehot_out/ehot

        # Step 3-6:
        cool, hot = self._sample_particles(ncool_,nhot_,sfr_)

        # get mass from energy
        vsqc = 0.5*((cool['vx']**2+cool['vy']**2+cool['vz']**2) + 5*cool['cs']**2)
        vsqh = 0.5*((hot['vx']**2+hot['vy']**2+hot['vz']**2) + 5*hot['cs']**2)
        mcool = ecool*1.e51/vsqc*self.vEsq
        mhot = ehot*1.e51/vsqh*self.vEsq
        cool['mass'] = mcool
        hot['mass'] = mhot

        return cool,hot


    def _sample_particles(self,ncool_,nhot_,sfr_):
        """Sampling particles for a given number of particles and SFR surface density
        """

        # get integer number of particles using poisson sampler
        ncool = np.atleast_1d(poisson.rvs(ncool_))
        nhot = np.atleast_1d(poisson.rvs(nhot_))

        Nc = ncool.sum()
        Nh = nhot.sum()

        # Step 3.0: prepare to draw particle's velocity and sound speed
        # this step is required to avoid for loops in actual sampling step
        # maybe there will be a better pythonic way for this, but
        # at least this is working and not too slow...

        # store indices of SFR that has non-zero number of particles
        coolidx=[]
        hotidx=[]
        for i,nc,nh in zip(range(len(sfr_)),ncool,nhot):
            for j in range(nc):
                coolidx.append(i)
            for k in range(nh):
                hotidx.append(i)

        # SFR surface density information for particles
        sfrcool = sfr_[coolidx]
        sfrhot = sfr_[hotidx]

        # Steps 3 and 4: Obtain particle velocity and sound speed
        vzc, cc = self._draw_cool(Nc,sfrcool)
        vzh, ch = self._draw_hot(Nh,sfrhot)

        # calculate vBz
        vBzc = np.sqrt(vzc**2 + 5*cc**2)
        vBzh = np.sqrt(vzh**2 + 5*ch**2)

        # Step 5: Assign metallicity
        Zc = self._Zmodel(vBzc,sfrcool,self.ZISM)
        Zh = self._Zmodel(vBzh,sfrhot,self.ZISM)

        # Step 6: Assign transverse velocity
        # calculate the magnitude of transverse velocity from the energy bias model
        bc = self._energy_bias(vBzc)
        bh = self._energy_bias(vBzh)

        vperpc = np.sqrt((1-bc)/bc)*vBzc
        vperph = np.sqrt((1-bh)/bh)*vBzh

        # draw uniform random number to assign vx and vy
        theta = np.random.rand(Nc)*2*np.pi
        vxc = vperpc*np.cos(theta)
        vyc = vperpc*np.sin(theta)

        theta = np.random.rand(Nh)*2*np.pi
        vxh = vperph*np.cos(theta)
        vyh = vperph*np.sin(theta)

        cool = dict(vx=vxc, vy=vyc, vz=vzc, cs=cc, Z=Zc, idx=coolidx)
        hot = dict(vx=vxh, vy=vyh, vz=vzh, cs=ch, Z=Zh, idx=hotidx)

        return cool,hot


    def _draw_hot(self,N,sfr,log=False):
        """Sample outflow velocity and sound speed of hot gas
        """
        xi = np.random.rand(N)
        vB = np.interp(xi,self.vB_cdf,self.y0)*self._vB0(sfr)
        xi = np.random.rand(N)
        Mach = np.interp(xi,self.Mach_cdf,self.y0)*self.Mach0

        cs = np.sqrt(vB**2/(Mach**2+5))
        vout = Mach*cs
        if log:
            u = np.log10(vout)
            w = np.log10(cs)
            return u,w

        return vout,cs

    def _draw_cool(self,N,sfr,log=False):
        """Sample outflow velocity and sound speed of cool gas
        """
        xi = np.random.rand(N)
        vout = np.interp(xi,self.vout_cdf,self.y0)*self._vout0(sfr)
        lncs = np.random.randn(N)*self.sigma + np.log(self.cs0)
        cs = np.exp(lncs)

        if log:
            u = np.log10(vout)
            w = np.log10(cs)
            return u,w
        return vout,cs

def to_time_series(p,time):
    """Function to convert the particle data into time series

    Parameters
    ----------
    p : dict
        paticle data as returned by `TigressWindSampler.draw` method
    time : array_like
        time array corresponding to SFR time series used to sample particles

    Returns
    -------
    out : (m, p, E, mZ)
        time series of mass, momemtum, energy, and metals carried by sampled particles
    """
    msum = np.zeros_like(time)
    psum = np.zeros_like(time)
    Esum = np.zeros_like(time)
    mZsum = np.zeros_like(time)

    to_erg = (ac.M_sun*(au.km/au.s)**2).to('erg').value
    value,index,counts = np.unique(p['idx'],return_index=True,return_counts=True)
    m = p['mass']*np.ones_like(p['vz'])
    mom = p['mass']*p['vz']+p['mass']*p['cs']**2/p['vz']
    e = 0.5*p['mass']*(p['vx']**2+p['vy']**2+p['vz']**2+5*p['cs']**2)*to_erg
    mZ = p['mass']*p['Z']

    #print(m.mean(),mom.mean(),e.mean(),mZ.mean())
    for i,i0,n in zip(value,index,counts):
        msum[i] = m[i0:i0+n].sum()
        psum[i] = mom[i0:i0+n].sum()
        Esum[i] = e[i0:i0+n].sum()
        mZsum[i] = mZ[i0:i0+n].sum()

    return msum,psum,Esum,mZsum
