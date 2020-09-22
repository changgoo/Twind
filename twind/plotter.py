from matplotlib.colors import ListedColormap,LogNorm
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import cmasher as cma
import pandas as pd

from matplotlib.ticker import LogLocator,AutoLocator,AutoMinorLocator,MaxNLocator

from .sampler import to_time_series, TigressWindSampler

discrete_cmap = ListedColormap(sns.color_palette('tab20c',n_colors=20,desat=0.5).as_hex())
pdf_cmap = cma.fall_r
pdfmin = -2.5
pdfmax = 1

def pdf_projection(pdf,wpdf=None,dvB=0.02,dMach=0.02):
    """Obtain 1D PDFs prjected onto u, w, log vB, log Mach

    Parameters
    ----------
    pdf : xarray.Dataset
        2D pdf of u and w (from TigressWindModel)
    wpdf : xarray.Dataset, optional
        weight field
    dvB : float
        log vB bin (the default is 0.02)
    dMach : float
        log Mach bin (the default is 0.02)

    Returns
    -------
    bins, pdf : list of tuples (bin, pdf) for all four variables
        [(u, pdf_u), (w, pdf_w), (logvB, pdf_logvB) (log Mach, pdf_logMach)]
    """
    dlogcs = np.diff(pdf.logcs)[0]
    dlogvout = np.diff(pdf.logvout)[0]
    dbinsq = dlogcs*dlogvout

    if wpdf is None:
        pdf_u = pdf.sum(dim=['logcs'])*dlogcs
        pdf_w = pdf.sum(dim=['logvout'])*dlogvout
    else:
        pdf_u = (pdf*wpdf).sum(dim=['logcs'])/wpdf.sum(dim=['logcs'])
        pdf_w = (pdf*wpdf).sum(dim=['logvout'])/wpdf.sum(dim=['logvout'])

    vBz_bins = np.arange(0,4,dvB)
    Mach_bins = np.arange(-2,2,dMach)

    cs = 10.**pdf.logcs
    vout = 10.**pdf.logvout
    vBz = np.sqrt(5.0*cs**2+vout**2)
    Mach = 1/cs*vout

    vBz.name = 'vBz'
    Mach.name = 'Mach'

    if wpdf is None:
        pdf_vBz = pdf.groupby_bins(np.log10(vBz),vBz_bins)
        pdf_vBz = pdf_vBz.sum(dim=['stacked_logcs_logvout'])*dbinsq/dvB
        pdf_Mach = pdf.groupby_bins(np.log10(Mach),Mach_bins)
        pdf_Mach = pdf_Mach.sum(dim=['stacked_logcs_logvout'])*dbinsq/dMach
    else:
        wpdf_vBz = wpdf.groupby_bins(np.log10(vBz),vBz_bins)
        wpdf_vBz = wpdf_vBz.sum(dim=['stacked_logcs_logvout'])
        wpdf_Mach = wpdf.groupby_bins(np.log10(Mach),Mach_bins)
        wpdf_Mach = wpdf_Mach.sum(dim=['stacked_logcs_logvout'])
        pdf_vBz = (pdf*wpdf).groupby_bins(np.log10(vBz),vBz_bins)
        pdf_vBz = pdf_vBz.sum(dim=['stacked_logcs_logvout'])/wpdf_vBz
        pdf_Mach = (pdf*wpdf).groupby_bins(np.log10(Mach),Mach_bins)
        pdf_Mach = pdf_Mach.sum(dim=['stacked_logcs_logvout'])/wpdf_Mach

    pdf_list=[(pdf.logvout,pdf_u),(pdf.logcs,pdf_w),
              (vBz_bins[:-1]+0.5*dvB,pdf_vBz.T),
              (Mach_bins[:-1]+0.5*dMach,pdf_Mach.T)]

    return pdf_list

def scifmt(value,fmt=':9.1e'):
    """Formatting float in scientific format"""
    maxdigits=int(fmt.split('.')[1][0])+1
    sp='{{{}}}'.format(fmt).format(value).split('e')
    digits=int(sp[1])
    if abs(digits) < maxdigits:
        return'{{{}}}'.format(':9.{}f'.format(maxdigits-digits-1)).format(value)
    else:
        if eval(sp[0]) == 1.0:
            return '10^{{{}}}'.format(int(sp[1]))
        else:
            return sp[0]+'\\cdot 10^{{{}}}'.format(int(sp[1]))

def add_vBM_grid(ax,pdf,vB_labels=True,M_labels=True):
    """Add contours of log vB and log Mach

    Parameters
    ----------
    ax : matplotlib.axes
    pdf : xarray.Dataset
        Dataset created by TigressWindModel or TigressSimLoader
    vB_labels, M_labels : bool
        add contour labels
    """
    levels=[1,1.5,2,2.5,3,3.5]
    u = pdf.logvout
    w = pdf.logcs
    vB = pdf['vBz']
    ct=ax.contour(u,w,np.log10(vB),levels=levels,colors='gray',alpha=1.0,linestyles=':')
    if vB_labels:
        ax.clabel(ct,[3],inline=1, inline_spacing=40,manual=[(1.0,3)],
                  fmt=r'$\log_{10} v_{\mathcal{B},z}=%1.1f$', fontsize='x-small')
        ax.clabel(ct,[1,1.5,2,2.5,3.5],inline=1, fmt=r'$%1.1f$', fontsize='x-small')

    levels=[-2,-1,0,1,2]
    Mach = pdf['Mach']
    ct=ax.contour(u,w,np.log10(Mach),levels=levels,colors='gray',alpha=1.0,linestyles=':')

    if M_labels:
        ax.clabel(ct,[-1],inline=1, inline_spacing=40,manual=[(0.5,2.5)],
                  fmt=r'$\log_{10} \mathcal{M}=%1.0f$', fontsize='x-small')
        ax.clabel(ct,[-2,0,1,2],inline=1,manual=[(0.5,2.5),(0.7,0.5),(1.5,0.4),(2.5,0.5)],
                 fmt=r'$%1.0f$', fontsize='x-small')

def add_phase_lines(ax,labels=True):
    """Add horizontal lines for phase separation
    """
    T1=2.e4
    T2=5.e5
    from .tigress_tools import T2cs
    ax.axhline(T2cs(T1),color='C0',ls='--',lw=1)
    ax.axhline(T2cs(T2),color='C1',ls='--',lw=1)
    if labels:
        ax.annotate('cool',(3.4,T2cs(T1)-0.02),color='C0',
                    ha='right',va='top',fontsize='small')
        ax.annotate('hot',(3.4,T2cs(T2)+0.02),color='C1',
                    ha='right',va='bottom',fontsize='small')
        ax.annotate('int.',(3.4,0.5*(T2cs(T1)+T2cs(T2))),color='C2',
                    ha='right',va='center',fontsize='small')

def toggle_xticks(axes,visible=False):
    plt.setp([ax.get_xticklabels() for ax in axes],visible=visible)
def toggle_yticks(axes,visible=False):
    plt.setp([ax.get_yticklabels() for ax in axes],visible=visible)

def show2d(pdf,xlabel=True,ylabel=True,label=None,vBM=[False,False],**kwargs):
    """Display 2D joint PDF as image
    """
    dbin = np.diff(pdf.logvout)[0]
    extent=[pdf.logvout.min(),pdf.logvout.max()+dbin,pdf.logcs.min(),pdf.logcs.max()+dbin]
    im=plt.imshow(np.log10(pdf),vmin=pdfmin,vmax=pdfmax,origin='lower',
                  extent=extent,cmap=pdf_cmap,**kwargs)
    ax=plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim(0,3.5)
    ax.set_ylim(0,3.5)
    if ylabel: ax.set_ylabel(r'$w\equiv\log_{10}\,c_s\,[\rm km/s]$')
    if xlabel: ax.set_xlabel(r'$u\equiv\log_{10}\,v_{\rm out}\,[\rm km/s]$')

    ax.set_yticks([0,1,2,3])
    ax.set_xticks([0,1,2,3])
    if label is not None:
        ax.annotate(label,(0.02,0.98),xycoords='axes fraction',
                    ha='left',va='top',fontsize='small')

    if any(vBM):
        add_vBM_grid(ax,pdf,vB_labels=vBM[0],M_labels=vBM[1])

    return im

def show2d_ct(pdf,xlabel=True,ylabel=True,label=None,**kwargs):
    """Display 2D joint PDF as contours
    """
    dbin = np.diff(pdf.logvout)[0]
    extent=[pdf.logvout.min(),pdf.logvout.max()+dbin,pdf.logcs.min(),pdf.logcs.max()+dbin]
    ct=plt.contour(pdf.logvout,pdf.logcs,np.log10(pdf),
                   vmin=pdfmin,vmax=pdfmax,
                   cmap=pdf_cmap,levels=np.arange(pdfmin,pdfmax+0.5,0.5),
                   linewidths=3,**kwargs)
    ax=plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim(0,3.5)
    ax.set_ylim(0,3.5)
    if ylabel: ax.set_ylabel(r'$w\equiv\log_{10}\,c_s\,[\rm km/s]$')
    if xlabel: ax.set_xlabel(r'$u\equiv\log_{10}\,v_{\rm out}\,[\rm km/s]$')

    ax.set_yticks([0,1,2,3])
    ax.set_xticks([0,1,2,3])
    if label is not None:
        ax.annotate(label,(0.02,0.98),xycoords='axes fraction',
                    ha='left',va='top',fontsize='small')

    return ct

def plot_flux_pdfs_yZ(pdf,grid=False):
    """Three-panel figure for Mpdf, Epdf, and yZ

    Figure 1 of the wind launching paper

    Parameters
    ----------
    pdf : xarray.Dataset
        a joint pdf from simulation (TigressSimLoader)
    grid : bool
        toggle grid of vBz and Mach
    """
    fig,axes=plt.subplots(3,1,figsize=(5,12),sharex=True,
                          gridspec_kw=dict(left=0.06,right=0.90,
                                           bottom=0.15,top=0.98,
                                           hspace=0.0))
    fields = ['Mpdf','Epdf','yZ']
    labels = [[r'(a) mass loading',r'$\log_{10}\, f_{M}(u,w)\,[{\rm dex^{-2}}]$'],
              [r'(b) energy loading',r'$\log_{10}\, f_{E}(u,w)\,[{\rm dex^{-2}}]$'],
              [r'(c) enrichment factor',r'$\zeta(u,w)$']]

    for i,ax,wf,lab in zip(range(len(fields)),axes[:],fields,labels):
        attrs=pdf.attrs
        pdfdata=pdf[wf]
        dbin=attrs['dbin']
        u = pdf.logvout
        w = pdf.logcs
        extent=[u.min(),u.max()+dbin,w.min(),w.max()+dbin]

        if wf.startswith('Z'):
            im=ax.imshow(pdfdata,origin='lower',extent=extent,
                         vmin=0.02,vmax=0.12,cmap=discrete_cmap)
        elif wf.startswith('yZ'):
            im=ax.imshow(pdfdata,origin='lower',extent=extent,
                         vmin=1,vmax=3.5,cmap=discrete_cmap)
        else:
            im=ax.imshow(np.log10(pdfdata),origin='lower',extent=extent,
                         vmin=pdfmin,vmax=pdfmax,cmap=pdf_cmap)

        ax.set_aspect('equal')
        ax.set_xlim(0,3.5)
        ax.set_ylim(0,3.5)
        ax.set_xticks([0,1,2,3])
        ax.set_yticks([0,1,2,3])
        if i == 2: ax.set_xlabel(r'$u\equiv\log_{10}\,v_{\rm out}\,[\rm km/s]$')

        add_phase_lines(ax,labels=(i==0))
        if grid:
            add_vBM_grid(ax,pdf,vB_labels=(i==0),M_labels=(i==1))
        else:
            add_vB_lines(ax,pdf,labels=(i==1))

        ax.annotate(lab[0],(0.02,0.98),xycoords='axes fraction',ha='left',va='top',fontsize='small')
        cbar=plt.colorbar(im,ax=ax,pad=0)
        cbar.set_label(lab[1])
        if i < 2:
            cbar.set_ticks([-3,-2,-1,0,1])
            cbar.ax.set_yticklabels(['',-2,-1,0,1])

        ax.set_ylabel(r'$w\equiv\log_{10}\,c_s\,[\rm km/s]$')

    return fig

def plot_vB_projection(sim,reconstructed=True,**kwargs):
    """Plot four 1D pdfs along log vB from simulation

    Figure 2(a) of the wind launching paper.

    Parameters
    ----------
    sim : TigressSimLoader

    """

    fields=['Mpdf','ppdf','Epdf','Zpdf']
    labels=[r'$f_M$',
            r'$f_p$',
            r'$f_E$',
            r'$f_Z$',]
    colors=['C3','C0','C1','C2']

    pdf = sim.simpdf
    for k,label,c in zip(fields,labels,colors):
        pdf1d = pdf_projection(pdf[k])
        x, y = pdf1d[2]
        l, = plt.step(x,y,label=label,lw=3,color=c,alpha=0.5)
        if reconstructed:
            if k+'_r' in pdf:
                pdf1d = pdf_projection(pdf[k+'_r'])
                x, y = pdf1d[2]
                plt.step(x,y,lw=1,color=l.get_color())

    plt.xlim(1,3.5)
    plt.ylabel(r'$f_q(\log_{10}\,v_{\mathcal{B},z})\,[{\rm dex}^{-1}]$')
    plt.xlabel(r'$\log_{10}\, v_{\mathcal{B},z}\,{\rm [km/s]}$')

def plot_vB_projection_ratio(sim,**kwargs):
    """Plot ratios of recunstructed and original 1D pdfs along log vB

    Figure 2(b) of the wind launching paper.

    Parameters
    ----------
    sim : TigressSimLoader

    """

    fields=['ppdf','Epdf','Zpdf']
    labels=[r'$f_{p}^r/f_{p}$',
            r'$f_{E}^r/f_{E}$',
            r'$f_{Z}^r/f_{Z}$',]

    colors=['C0','C1','C2','C3','C4','C5','C6','C7']

    for k,label,c in zip(fields,labels,colors):
        ratio = sim.simpdf[k+'_r']/sim.simpdf[k]
        pdf1d = pdf_projection(ratio,wpdf=sim.simpdf['Mpdf'])
        x, y = pdf1d[2]
        l,=plt.plot(x,y,label=label,color=c,**kwargs)

    plt.xlim(1,3.5)
    plt.ylabel('Ratio')
    plt.xlabel(r'$\log_{10}\, v_{\mathcal{B},z}\,{\rm [km/s]}$')

def flux_reconstruction(sims,stdmodel='R4'):
    """Two-panel figure for comparison of orignal and reconstructed PDFs

    Figure 2 of the wind launching paper

    Parameters
    ----------
    sims : TigressSimContainer
    stdmodel : ['R2','R4','R8','R16','LGR2','LGR4','LGR8']

    """
    fig,axes = plt.subplots(1,2,figsize=(10,4))

    plt.sca(axes[0])
    plot_vB_projection(sims[stdmodel])
    plt.legend(loc=9,fontsize='x-small')
    plt.xlim(1,3.5)
    plt.xticks([1,2,3,])
    plt.yscale('log')
    plt.ylim(1.e-2,10)
    plt.annotate('(a)',(0.03,0.97),xycoords='axes fraction',ha='left',va='top')

    plt.sca(axes[1])
    plot_vB_projection_ratio(sims[stdmodel],lw=3)
    plt.legend(loc=4,fontsize='small')
    plt.annotate('(b)',(0.03,0.97),xycoords='axes fraction',ha='left',va='top')

    for k,sim in sims.items():
        plot_vB_projection_ratio(sim,lw=1,alpha=0.5)

    plt.ylim(0,1.3)
    plt.xlim(1,3.5)
    plt.xticks([1,2,3,])
    plt.axhline(1,ls=':')
    plt.axhline(0,ls=':')

    plt.tight_layout()

    return fig


def plot_proj(pdf1d,axes,**kwargs):
    """Plot all four 1D histograms

    Figure 4 right panels.

    Parameters
    ----------
    pdf1d : list
        bin and 1d pdf prjected on different axes
        created by pdf_projection
    axes : list
        axes to plot

    """

    xlabels=[r'$v_{\rm out} [{\rm km/s}]$',
             r'$c_s [{\rm km/s}]$',
             r'$v_{\mathcal{B},z} [{\rm km/s}]$',
             r'$\mathcal{M}$']
    for ax,(x,y),xlab in zip(axes,pdf1d,xlabels):
        plt.sca(ax)
        plt.step(x,y,where='mid',**kwargs)
        plt.yscale('log')
        plt.ylim(1.e-3,10)


def comparison_pdfs(sim,q='M'):
    """Five-panel figure of 2D and 1D PDFs

    Figure 3 of the wind launching paper

    Parameters
    ----------
    sim : TigressSimLoader
    q : ['M','p','E','Z']
        varialble to plot

    """

    modelpdf=sim.build_model()

    def set_fivepanel_axes(nrows=2):
        '''
            construct custom axes for one 2D PDFs and four 1D histograms
        '''

        fig=plt.figure(figsize=(20,4*nrows),constrained_layout=False)
        outergrid=fig.add_gridspec(nrows,2,width_ratios=[1.2,4],wspace=0.2)
        axes=[]
        for i in range(nrows):
            axes_row=[]
            innergrid=outergrid[i*2].subgridspec(1,2,width_ratios=[1,0.05],wspace=0)
            for ig in innergrid:
                axes_row.append(fig.add_subplot(ig))

            innergrid=outergrid[i*2+1].subgridspec(1,4,width_ratios=[1,1,1,1],wspace=0)
            for ig in innergrid:
                axes_row.append(fig.add_subplot(ig))
            axes.append(axes_row)
        if nrows == 1:
            return fig,axes[0]
        else:
            return fig,axes

    fig,_axes = set_fivepanel_axes(nrows=1)

    # draw simulation PDF
    plt.sca(_axes[0])
    im = show2d(sim.simpdf[q+'pdf'],alpha=0.7)
    plt.colorbar(im,cax=_axes[1])

    # draw 1D projections
    pdf1d = pdf_projection(sim.simpdf[q+'pdf'])
    plot_proj(pdf1d,_axes[2:],color='k',lw=1,label='sim.')

    # draw model PDF contour
    plt.sca(_axes[0])
    ct = show2d_ct(modelpdf[q+'pdf'])
    cbar = plt.colorbar(ct,cax=_axes[1])
    cbar.set_ticks([-2,-1,0,1])

    # draw 1D projections
    pdf1d = pdf_projection(modelpdf[q+'pdf-cool'])
    plot_proj(pdf1d,_axes[2:],lw=2,label='cool',zorder=1)
    pdf1d = pdf_projection(modelpdf[q+'pdf-hot'])
    plot_proj(pdf1d,_axes[2:],lw=2,label='hot',zorder=1)
    pdf1d = pdf_projection(modelpdf[q+'pdf'])
    plot_proj(pdf1d,_axes[2:],lw=3,label='total',zorder=0)

    toggle_yticks(_axes[3:])

    xlabels=[r'$\log_{10}\,v_{\rm out}$',r'$\log_{10}\,c_s$',
            r'$\log_{10}\,v_{\mathcal{B},z}$',r'$\log_{10}\,\mathcal{M}$']
    xlims=[(0,3.5),(0,3.5),(0.5,4),(-1.5,2)]
    for ax,xlab,xlim in zip(_axes[2:],xlabels,xlims):
        ax.set_xlabel(xlab+r'$\,[\rm km/s]$')
        ax.set_xlim(xlim)
    _axes[2].set_ylabel(r'$f_{}\,[{{\rm dex}}^{{-1}}]$'.format(q))

    _axes[-1].legend(fontsize='xx-small')
    for ax in _axes[2:]:
        ax.yaxis.set_major_locator(LogLocator(numticks=15))
        ax.yaxis.set_minor_locator(LogLocator(subs=np.arange(2, 10) * .1, numticks=15))
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    return fig

def show_loading(model,vlist=[0,30,1000,300,100],sims=None):
    """Loading factor scaling as a function of Sigma_SFR

    Figure 4 of the wind launching paper

    Parameters
    ----------
    model : TigressWindModel
    vlist : list
        list of escape velocities
    sims : TigressSimContainer (optional)
        if passed, overplot simulation data point

    """
    fig,axes = plt.subplots(1,4,sharey='row',sharex='col',figsize=(18,5),
                           gridspec_kw=dict(hspace=0.1,wspace=0.1))

    sfr=10.**model.logsfr
    vBz=model.vBz

    for ax,q in zip(axes,'MpEZ'):
        for i,vesc in enumerate(vlist):
            dbinsq=model.attrs['dlogvout']*model.attrs['dlogcs']
            eta = model['eta'+q]
            pdf = model[q+'pdf']
            cdf = pdf.where(vBz>vesc).sum(dim=['logcs','logvout'])*dbinsq
            loading = cdf*eta

            plt.sca(ax)
            if sims is not None:
                for k,sim in sims.items():
                    simpdf = sim.simpdf
                    simeta = sim.simpdf.attrs['eta'+q]
                    dbin = sim.simpdf.attrs['dbin']
                    simsfr = sim.simpdf.attrs['sfr']
                    simcdf = simpdf[q+'pdf'].where(sim.simpdf['vBz']>vesc).sum()*dbin**2
                    simloading = simcdf*simeta
                    plt.plot(simsfr,simloading,'o',color='C{}'.format(i))

            plt.plot(sfr,loading,color='C{}'.format(i))

        plt.xscale('log')

    axes[0].set_ylabel(r'$\eta\,(v_{\mathcal{B},z}>v_{\rm esc})$')
    plt.setp(axes,'yscale','log')
    plt.setp(axes,'xscale','log')
    plt.setp(axes,'xlim',(2.e-5,2))
    plt.setp(axes,'ylim',(5.e-3,1.e2))

    ax=axes[0]
    ax.annotate(r'$v_{\rm esc}=0$',(5.e-4,50),ha='left',va='top',color='C0',fontsize='small')
    ax.annotate('$30$',(1.e-4,5),ha='left',va='bottom',color='C1',fontsize='small')
    ax.annotate('$100$',(1.5e-4,0.35),ha='left',va='bottom',color='C4',fontsize='small')
    ax.annotate('$300$',(1.5e-4,0.25),ha='left',va='top',color='C3',fontsize='small')
    ax.annotate(r'$10^3$',(2.e-4,0.015),ha='left',va='bottom',color='C2',fontsize='small')
    plt.setp(axes,'xlabel',r'$\Sigma_{\rm SFR}\,[{\rm M_\odot\,kpc^{-2}\,yr^{-1}}]$')

    return fig

def sampling_from_simulation_sfr(sim,tdelay=10):
    """Three-panel figure for
    (a) mass outflow rate of cool gas
    (b) energy outflow rate of hot gas
    (c) distribution of sampled particles

    Parameters
    ----------
    sim : TigressSimLoader
    tdelay : float
        time delay to be applied to sampled data in Myr
    """

    sampler=TigressWindSampler(z0=sim.z0)

    ts = sim.time_series
    sfr = ts['sfr10']
    time = ts['tMyr']
    dt = 1.e6
    area = 1.024
    mc = (1.e4,1.e5,1.e6)
    mh = (1.e2,1.e3,1.e4)

    # determine indices that match the desired time range
    torb = ts.torb
    idx = np.arange(len(torb))
    idx = idx[(torb>1) & (torb<2)]
    imin, imax = idx.min(), idx.max()

    fig=plt.figure(figsize=(18,8),constrained_layout=False)
    og=fig.add_gridspec(2,2,width_ratios=[1.5,1],wspace=0.2,hspace=0.1)

    axes=[]

    axes.append(fig.add_subplot(og[0,0]))
    axes.append(fig.add_subplot(og[1,0]))
    axes.append(fig.add_subplot(og[:,1]))

    refs,etas,etac,etah = sampler.get_refs(sfr)

    axes[0].fill_between(time,0,refs[0]*area,color='k',alpha=0.2,lw=0)
    axes[1].fill_between(time,0,refs[2]*area,color='k',alpha=0.2,lw=0)

    for m1, m2, alpha in zip(mc,mh,[0.5,0.7,0.9]):
        cool,hot = sampler.draw_mass(sfr,m1,m2,area=area,dt=dt)

        ts_cool = to_time_series(cool,time)
        ts_hot = to_time_series(hot,time)

        # plot mass outflow rate of cool gas
        plt.sca(axes[0])

        sr = pd.Series(ts_cool[0]/dt,time)
        plt.plot(time+tdelay,sr.rolling(10).mean(),
                 lw=2,label=r'${}M_\odot$'.format(scifmt(m1)))
        plt.yscale('log')
        plt.ylim(1.e-3,1)
        plt.legend(loc=1,ncol=3,title=r'$m^{{\rm cool}}$',
                   fontsize='x-small',framealpha=1.0)
        plt.ylabel(r'$\dot{M}_{\rm out}\,[M_\odot\,{\rm yr^{-1}}]$')
        plt.annotate('(a)',(0.02,0.95),xycoords='axes fraction',ha='left',va='top')
        toggle_xticks([axes[0]])

        # plot energy outflow rate of hot gas
        plt.sca(axes[1])
        sr = pd.Series(ts_hot[2]/dt,time)
        plt.plot(time+tdelay,sr.rolling(10).mean(),
                 lw=2,label=r'${}M_\odot$'.format(scifmt(m2))+
                            r'$('+scifmt(ts_hot[2].mean(),fmt=':9.1e')+
                            r'{{\rm erg}})$')
        plt.ylabel(r'$\dot{E}_{\rm out}\,[{\rm erg\,yr^{-1}}]$')
        plt.xlabel(r'${\rm time [Myr]}$')
        plt.legend(loc=1,fontsize='x-small',framealpha=1.0,ncol=3,
                   title=r'$m^{\rm hot}( \overline{e^{{\rm hot}}})$')
        plt.annotate('(b)',(0.02,0.95),xycoords='axes fraction',ha='left',va='top')


        plt.setp(axes[:2],'xlim',(0,600))

        plt.yscale('log')

        # plot sampled particle distribution on top of simulation distribution
        # only particles in 1<t/torb<2

        plt.sca(axes[2])

        show2d(sim.simpdf['Mpdf'],alpha=0.5)

        cool_idx = (cool['idx']>imin) & (cool['idx']<imax)
        hot_idx = (hot['idx']>imin) & (hot['idx']<imax)

        uc = np.log10(cool['vz'][cool_idx])
        wc = np.log10(cool['cs'][cool_idx])
        uh = np.log10(hot['vz'][hot_idx])
        wh = np.log10(hot['cs'][hot_idx])

        l,=plt.plot(uc,wc,'o',alpha=alpha,markeredgewidth=0,
                    label=r'$N^{{\rm cool}}={}$'.format(len(uc)))
        plt.plot(uh,wh,'s',alpha=alpha,color=l.get_color(),
                 markeredgewidth=0,label=r'$N^{{\rm hot}}={}$'.format(len(uh)))

        plt.annotate('(c)',(0.05,0.95),xycoords='axes fraction',ha='left',va='top')
        plt.legend(fontsize='small',loc=4)

    axes[0].set_ylim(1.e-3,1)
    axes[0].plot(ts['tMyr'],ts['mass_whole'],color='k')
    axes[1].set_ylim(1.e44,2.e48)
    axes[1].plot(ts['tMyr'],ts['energy_whole'],color='k')

    return fig
