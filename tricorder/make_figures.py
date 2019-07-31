#### Make figures for paper (in eps/pdf format)
from __future__ import division
import numpy as np

import matplotlib
matplotlib.use('agg')

import matplotlib.pyplot as plt
import chainconsumer
from scipy.optimize import minimize
import plottools
reload(plottools)
import palettable
from matplotlib.colors import LogNorm

from astropy.io import fits
import paths
import skymapper as skm

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rc('xtick', labelsize=12) 
matplotlib.rc('ytick', labelsize=12) 
matplotlib.rc('font', family='serif')

def make_easy_figures():
    ## 
    ##Footprint of simulated gals
    ##
    HD = fits.getdata(paths.rm_y1[0])
    HL = fits.getdata(paths.rm_y1_HL[0])
    HHL = fits.getdata(paths.rm_y1_HHL[0])

    slice1 = HD[((HD['ZREDMAGIC'] > .15) & (HD['ZREDMAGIC'] < .3))]
    slice2 = HD[((HD['ZREDMAGIC'] > .3) & (HD['ZREDMAGIC'] < .45))]
    slice3 = HD[((HD['ZREDMAGIC'] > .45) & (HD['ZREDMAGIC'] < .6))]
    slice4 = HL[((HL['ZREDMAGIC'] > .6) & (HL['ZREDMAGIC'] < .75))]
    slice5 = HHL[((HHL['ZREDMAGIC'] > .75) & (HHL['ZREDMAGIC'] < .9))]

    sample = np.concatenate([slice1,slice2,slice3,slice4,slice5])
    skm.plotDensity(sample['RA'],sample['DEC'],sep=10)

    plt.savefig('./figures/skymap.pdf',dpi=300,bbox_inches='tight')
    plt.savefig('./figures/skymap.png',dpi=300,bbox_inches='tight')

    plt.figure()
    
    ##
    ## Redmagic n(z) histogram
    ##
    plt.hist(slice1['ZSPEC'],histtype='step',bins=100,range=(0,1))
    plt.hist(slice2['ZSPEC'],histtype='step',bins=100,range=(0,1))
    plt.hist(slice3['ZSPEC'],histtype='step',bins=100,range=(0,1))
    plt.hist(slice4['ZSPEC'],histtype='step',bins=100,range=(0,1))
    plt.hist(slice5['ZSPEC'],histtype='step',bins=100,range=(0,1));
    plt.axvspan(.15,.3,alpha=.2,color='b')
    plt.axvspan(.3,.45,alpha=.2,color='orange')
    plt.axvspan(.45,.6,alpha=.2,color='g')
    plt.axvspan(.6,.75,alpha=.2,color='r')
    plt.axvspan(.75,.9,alpha=.2,color='violet')
    plt.xlabel('True Redshift')
    plt.ylabel('N(z)')

    plt.savefig('./figures/tomobins.pdf',dpi=300,bbox_inches='tight')
    plt.savefig('./figures/tomobins.png',dpi=300,bbox_inches='tight')

    plt.figure()

    ##
    ## Redmagic errors vs redshift
    ##

    plt.hist2d(sample['ZREDMAGIC'],sample['ZREDMAGIC_E'],bins=100,norm=LogNorm());
    plt.vlines([.3,.45,.6,.75],0.007,0.1,linestyle='--',color='r')
    plt.ylim(0.007,0.1);
    plt.xlabel(r'$z_{RM}$')
    plt.ylabel(r'$\sigma(z_{RM})$')

    plt.savefig('./figures/redmagic_errors.pdf',dpi=300,bbox_inches='tight')
    plt.savefig('./figures/redmagic_errors.png',dpi=300,bbox_inches='tight')

def make_triangle_diagrams():
    colors = palettable.colorbrewer.qualitative.Paired_5.hex_colors
    def drawTriangleScatter(min_r,max_r,min_u,max_u,min_v,max_v,color=0,num=20):
        test = np.array(np.meshgrid(np.linspace(min_r,max_r,num=num),
                np.linspace(min_u,max_u,num=num),
                np.linspace(min_v,max_v,num=10*num),))
        d2 = test[0,:,:,:]
        d3 = test[0,:,:,:]*test[1,:,:,:]
        d1 = np.abs(test[2,:,:,:])*d3+d2 
        cosine = (d2**2 + d3**2 - d1**2)/(2*d2*d3+1e-9)
        sine = np.sqrt(1-cosine**2)
        points = d3*cosine,d3*sine
        plt.scatter(np.array(points).reshape(2,-1)[0,:],np.array(points).reshape(2,-1)[1,:],color=colors[color],alpha=.1,marker='.')
        plt.plot([0,20],[0,0],color='k')
        plt.plot([0,-4,20],[0,12,0],color='k',linestyle='--')
        plt.plot(20,0,color='r',marker='s')
        plt.xlim(-30,25)
        plt.ylim(-5,25)

    drawTriangleScatter(19.9,20.1,.3,1,0,.2, color=0,num=50)
    drawTriangleScatter(19.9,20.1,.3,1,.2,.4,color=1,num=50)
    drawTriangleScatter(19.9,20.1,.3,1,.4,.6,color=2,num=50)
    drawTriangleScatter(19.9,20.1,.3,1,.6,.8,color=3,num=50)
    drawTriangleScatter(19.9,20.1,.3,1,.8,1, color=4,num=100)
    plt.xlabel('RA (arcmin)')
    plt.ylabel('DEC (arcmin)')   

    plt.savefig('./figures/perfect_binning.pdf',dpi=300,bbox_inches='tight')
    plt.savefig('./figures/perfect_binning.png',dpi=300,bbox_inches='tight')

    plt.figure() 

    drawTriangleScatter(16,24,.3,1,0,.2, color=0,num=50)
    drawTriangleScatter(16,24,.3,1,.2,.4,color=1,num=50)
    drawTriangleScatter(16,24,.3,1,.4,.6,color=2,num=50)
    drawTriangleScatter(16,24,.3,1,.6,.8,color=3,num=50)
    drawTriangleScatter(16,24,.3,1,.8,1, color=4,num=120)    

    plt.savefig('./figures/tolerance_20.pdf',dpi=300,bbox_inches='tight')
    plt.savefig('./figures/tolerance_20.png',dpi=300,bbox_inches='tight')


def make_covariance_fig():

    path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/photoz_syst/'

    zmin = .45
    zmax = zmin+.15
    data = plottools.load_res_indep(path, 'dm','newpaper13.1','ZSPEC',zmin,zmax,'12x20',sigma=0)
    data2 = plottools.load_res_indep(path,'dm','newpaper13.1','ZREDMAGIC',zmin,zmax,'12x20',sigma=0)
    galdata = plottools.load_res_indep(path, 'newbuzzardrm2','newpaper13.1','ZSPEC',zmin,zmax,'20',sigma=0)
    galdata2 = plottools.load_res_indep(path,'newbuzzardrm2','newpaper13.1','ZREDMAGIC',zmin,zmax,'20',sigma=0)
    plt.figure()
    fig,axarr = plt.subplots(2,2,figsize=(6,6))
    
    vmin= 0
    vmax= 1
    
    axarr[0,0].imshow(np.corrcoef(plottools.compress_dv(galdata['Q'].values.reshape(-1,10)).T),vmin=vmin,vmax=vmax,origin='lower',extent=(0,1,0,1))
    axarr[0,0].set_title(r'Galaxies: True n(z)')
    axarr[0,0].set_xticks([.1,.3,.5,.7,.9])
    axarr[0,0].set_yticks([.1,.3,.5,.7,.9])
    #axarr[0,0].set_xlabel('v')
    axarr[0,0].set_ylabel('v')
    axarr[0,1].imshow(np.corrcoef(plottools.compress_dv(galdata2['Q'].values.reshape(-1,10)).T),vmin=vmin,vmax=vmax,origin='lower',extent=(0,1,0,1))
    axarr[0,1].set_title(r'Galaxies: RedMaGiC n(z)')
    #axarr[0,1].set_xlabel('v')
    #axarr[0,1].set_ylabel('v')
    axarr[0,1].set_xticks([.1,.3,.5,.7,.9])
    axarr[0,1].set_yticks([.1,.3,.5,.7,.9])
    axarr[1,0].imshow(np.corrcoef(plottools.compress_dv(data['Q'].values.reshape(-1,10)).T),vmin=vmin,vmax=vmax,origin='lower',extent=(0,1,0,1))
    axarr[1,0].set_title(r'DM: True n(z)')
    axarr[1,0].set_xticks([.1,.3,.5,.7,.9])
    axarr[1,0].set_yticks([.1,.3,.5,.7,.9])
    axarr[1,0].set_xlabel('v')
    axarr[1,0].set_ylabel('v')
    im = axarr[1,1].imshow(np.corrcoef(plottools.compress_dv(data2['Q'].values.reshape(-1,10)).T),vmin=vmin,vmax=vmax,origin='lower',extent=(0,1,0,1))
    axarr[1,1].set_title(r'DM: RedMaGiC n(z)')
    axarr[1,1].set_xlabel('v')
    #axarr[1,1].set_ylabel('v')
    axarr[1,1].set_xticks([.1,.3,.5,.7,.9])
    axarr[1,1].set_yticks([.1,.3,.5,.7,.9])
    
    #fig.colorbar(im, ax=axarr.ravel().tolist(),shrink=.4,aspect=10)
    #plt.suptitle(str(zmin)+r'$<z<$'+str(zmax))
    plt.tight_layout()
    cbar_ax = fig.add_axes([.97, 0.2, .04, 0.6])
    fig.colorbar(im, cax=cbar_ax, orientation="vertical")

    plt.savefig('./figures/corrmat_examples.pdf',dpi=300,bbox_inches='tight')
    plt.savefig('./figures/corrmat_examples.png',dpi=300,bbox_inches='tight') 

def make_data_vector_plots():

    path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/photoz_syst/'

    for i,zmin in enumerate([.15, .3, .45, .6]):
        zmax = zmin+.15
        data =  plottools.load_res_indep(path,'dm','newpaper13.1','ZSPEC',zmin,zmax,'12x20',sigma=0)
        galdata = plottools.load_res_indep(path,'newbuzzardrm2','newpaper13.1','ZSPEC',zmin,zmax,'20',sigma=0)
        plt.figure()
        plottools.plot_dv(data,'Q', offset=0.,indiv_runs=True,color='Blue',label=r'Dark Matter',compressed=True)
        plottools.plot_dv(galdata,'Q', offset=0.,indiv_runs=True,color='Red',label='Galaxies',compressed=True)
        plt.xlabel('v')
        plt.ylabel('Q')
        plt.title(str(zmin)+r'$<z<$'+str(zmax))
        plt.legend() 

        plt.savefig('./figures/spec_bin'+str(i+1)+'_dv.pdf',dpi=300,bbox_inches='tight')
        plt.savefig('./figures/spec_bin'+str(i+1)+'_dv.png',dpi=300,bbox_inches='tight')
        plt.figure()


def get_zspec(is11k=False):
    summaries = []
    for i,zmin in enumerate([.15,.3,.45,.6]):
        if zmin != .45:
            path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/tolerance_syst/'
            config_fname = 'newpaper14.1'
        else:
            path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/photoz_syst/'
            config_fname = 'newpaper13.1'
        config = plottools.load_config(config_fname)
        zmax = zmin+.15
        cc = chainconsumer.ChainConsumer()
        if ((zmin == .15) | (zmin == .45)):
            data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'12x20',sigma=0)
            data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'20',sigma=0)
        else:
            data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0)
            data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'True redshifts')
        
        summary = cc.analysis.get_summary()
        summaries.append(summary)
       
        if ((zmin == .15) | (zmin == .45)):
            xi1 = plottools.load_res_xi_indep(path,'dm',config_fname,zmin,zmax,'12x20')
            xi2 = plottools.load_res_xi_indep(path,'newbuzzardrm2',config_fname,zmin,zmax,'20')
        else:
            xi1 = plottools.load_res_xi_indep(path,'dm',config_fname,zmin,zmax,'10x10')
            xi2 = plottools.load_res_xi_indep(path,'newbuzzardrm2',config_fname,zmin,zmax,'10')
        #b1 = np.sqrt(np.mean(xi2,axis=0)/np.mean(xi1,axis=0))
        b1 = np.sqrt(xi2/np.mean(xi1,axis=0))
        r = np.logspace(config['2PCF']['min_sep'],np.log10(config['2PCF']['max_sep']),num=config['2PCF']['nbins'])
        b1_mean = np.mean(b1[:,(r > config['3PCF']['min_sep']*config['3PCF']['min_u']) & (r < config['3PCF']['max_sep']*config['3PCF']['max_u'])])
        b1_std = np.std(b1[:,(r > config['3PCF']['min_sep']*config['3PCF']['min_u']) & (r < config['3PCF']['max_sep']*config['3PCF']['max_u'])])
        
        figure = cc.plotter.plot(figsize='column',extents=[(0,3),(-3,3)]);
        figure.axes[2].axvspan(b1_mean-b1_std, b1_mean+b1_std, alpha=0.3, color='b')
        figure.axes[2].axvspan(b1_mean-2*b1_std, b1_mean+2*b1_std, alpha=0.3, color='b')

        figure.axes[2].plot(np.linspace(0,3,100),plottools.get_lazeyras_schmidt_2015(np.linspace(0,3,100)),label='Lazeyras+Schmidt 2015')
        figure.axes[2].plot(np.linspace(0,3,100),plottools.get_hoffman_2015(np.linspace(0,3,100)),label='Hoffman 2015')

        if is11k:
            str11k = '_11k'
        else:
            str11k = ''
 
        axarr = figure.get_axes()
        #plt.suptitle(str(zmin)+r'$ < z < $'+str(zmax))
        plt.legend()
        plt.savefig('./figures/spec_bin'+str(i+1)+str11k+'.pdf',dpi=300,bbox_inches='tight')
        plt.savefig('./figures/spec_bin'+str(i+1)+str11k+'.png',dpi=300,bbox_inches='tight')

def get_zspec_both():
    summaries = []
    for i,zmin in enumerate([.15,.3,.45,.6]):
        if zmin != .45:
            path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/tolerance_syst/'
            config_fname = 'newpaper14.1'
        else:
            path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/photoz_syst/'
            config_fname = 'newpaper13.1'
        config = plottools.load_config(config_fname)
        zmax = zmin+.15
        cc = chainconsumer.ChainConsumer()
        if ((zmin == .15) | (zmin == .45)):
            data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'12x20',sigma=0)
            data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'20',sigma=0)
        else:
            data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0)
            data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=False)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'DES Y1-like',color='b')
        samples2 = plottools.make_inference(red_qdm,red_qrm,is11k=True)
        cc.add_chain(samples2.flatchain,parameters=['b1','b2'],name=r'12,000 sq. deg. survey', color='k')
        
        summary = cc.analysis.get_summary()
        summaries.append(summary)
       
        if ((zmin == .15) | (zmin == .45)):
            xi1 = plottools.load_res_xi_indep(path,'dm',config_fname,zmin,zmax,'12x20')
            xi2 = plottools.load_res_xi_indep(path,'newbuzzardrm2',config_fname,zmin,zmax,'20')
        else:
            xi1 = plottools.load_res_xi_indep(path,'dm',config_fname,zmin,zmax,'10x10')
            xi2 = plottools.load_res_xi_indep(path,'newbuzzardrm2',config_fname,zmin,zmax,'10')
        #b1 = np.sqrt(np.mean(xi2,axis=0)/np.mean(xi1,axis=0))
        b1 = np.sqrt(xi2/np.mean(xi1,axis=0))
        r = np.logspace(config['2PCF']['min_sep'],np.log10(config['2PCF']['max_sep']),num=config['2PCF']['nbins'])
        b1_mean = np.mean(b1[:,(r > config['3PCF']['min_sep']*config['3PCF']['min_u']) & (r < config['3PCF']['max_sep']*config['3PCF']['max_u'])])
        b1_std = np.std(b1[:,(r > config['3PCF']['min_sep']*config['3PCF']['min_u']) & (r < config['3PCF']['max_sep']*config['3PCF']['max_u'])])
        
        figure = cc.plotter.plot(figsize='column',extents=[(0,3),(-3,3)]);
        figure.axes[2].axvspan(b1_mean-b1_std, b1_mean+b1_std, alpha=0.3, color='b')
        figure.axes[2].axvspan(b1_mean-2*b1_std, b1_mean+2*b1_std, alpha=0.3, color='b')

        figure.axes[2].plot(np.linspace(0,3,100),plottools.get_lazeyras_schmidt_2015(np.linspace(0,3,100)),label='Lazeyras+Schmidt 2015')
        figure.axes[2].plot(np.linspace(0,3,100),plottools.get_hoffman_2015(np.linspace(0,3,100)),label='Hoffman 2015')

        if is11k:
            str11k = '_11k'
        else:
            str11k = ''
 
        axarr = figure.get_axes()
        plt.suptitle(str(zmin)+r'$ < z < $'+str(zmax))
        plt.legend()
        plt.savefig('./figures/both_bin'+str(i+1)+'.pdf',dpi=300,bbox_inches='tight')
        plt.savefig('./figures/both_bin'+str(i+1)+'.png',dpi=300,bbox_inches='tight')

def get_zspec_zrm(is11k=False):
    path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/photoz_syst/'
    summaries = []
    for i,zmin in enumerate([.15,.3,.45,.6]):
        if zmin != .45:
            path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/tolerance_syst/'
            config_fname = 'newpaper14.1'
        else:
            path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/photoz_syst/'
            config_fname = 'newpaper13.1'
        config = plottools.load_config(config_fname)
        zmax = zmin+.15
        cc = chainconsumer.ChainConsumer()
        if ((zmin == .15) | (zmin == .45)):
            data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'12x20',sigma=0)
            data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'20',sigma=0)
        else:
            data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0)
            data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'True redshifts')
        
        if zmin != .45:
            data1 = plottools.load_res_indep(path,'dm',config_fname,'ZREDMAGIC',zmin,zmax,'10x10',sigma=0)
            data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZREDMAGIC',zmin,zmax,'10',sigma=0)
        else:
            data1 = plottools.load_res_indep(path,'dm',config_fname,'ZREDMAGIC',zmin,zmax,'12x20',sigma=0)
            data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZREDMAGIC',zmin,zmax,'20',sigma=0)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'redMaGiC redshifts')
        
        summary = cc.analysis.get_summary()
        summaries.append(summary)
       
        if is11k:
            str11k = '_11k'
        else:
            str11k = ''
 
        figure = cc.plotter.plot(figsize='column',extents=[(0,3),(-3,3)]);
        axarr = figure.get_axes()
        plt.suptitle(str(zmin)+r'$ < z < $'+str(zmax))
        plt.legend()
        plt.savefig('./figures/rm_bin'+str(i+1)+str11k+'.pdf',dpi=300,bbox_inches='tight')
        plt.savefig('./figures/rm_bin'+str(i+1)+str11k+'.png',dpi=300,bbox_inches='tight')

def get_tolerance_figs(is11k=False):
    oldpath = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/photoz_syst/'
    path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/tolerance_syst/'
    for i,zmin in enumerate([.15,.3,.45,.6]):
        #config_fname = 'newpaper13.1'
        #config = plottools.load_config(config_fname)
        zmax = zmin+.15
        cc = chainconsumer.ChainConsumer()
        
        data1 = plottools.load_res_indep(oldpath,'dm','newpaper13.1','ZSPEC',zmin,zmax,'12x20',sigma=0)
        data2 = plottools.load_res_indep(oldpath,'newbuzzardrm2','newpaper13.1','ZSPEC',zmin,zmax,'20',sigma=0)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'Tolerance $= 20\%$')
        
        data1 = plottools.load_res_indep(path,'dm','newpaper13.1_up','ZSPEC',zmin,zmax,'12x20',sigma=0)
        data2 = plottools.load_res_indep(path,'newbuzzardrm2','newpaper13.1_up','ZSPEC',zmin,zmax,'20',sigma=0)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'Tolerance $= 25\%$')
        
        data1 = plottools.load_res_indep(path,'dm','newpaper13.1_down','ZSPEC',zmin,zmax,'12x20',sigma=0)
        data2 = plottools.load_res_indep(path,'newbuzzardrm2','newpaper13.1_down','ZSPEC',zmin,zmax,'20',sigma=0)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'Tolerance $= 15\%$')
        
        data1 = plottools.load_res_indep(path,'dm','newpaper13.1_down2','ZSPEC',zmin,zmax,'12x20',sigma=0)
        data2 = plottools.load_res_indep(path,'newbuzzardrm2','newpaper13.1_down2','ZSPEC',zmin,zmax,'20',sigma=0)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'Tolerance $= 10\%$')
        
        data1 = plottools.load_res_indep(path,'dm','newpaper13.1_up2','ZSPEC',zmin,zmax,'12x20',sigma=0)
        data2 = plottools.load_res_indep(path,'newbuzzardrm2','newpaper13.1_up2','ZSPEC',zmin,zmax,'20',sigma=0)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'Tolerance $= 30\%$')
        
        figure = cc.plotter.plot(figsize='column',extents=[(0,3),(-3,3)]);
        axarr = figure.get_axes()
        plt.suptitle(str(zmin)+r'$ < z < $'+str(zmax))
        plt.legend()

        if is11k:
            str11k = '_11k'
        else:
            str11k = ''

        plt.savefig('./figures/tol_bin'+str(i+1)+str11k+'.pdf',dpi=300,bbox_inches='tight')
        plt.savefig('./figures/tol_bin'+str(i+1)+str11k+'.png',dpi=300,bbox_inches='tight')

def get_u_figs(is11k=False): 
    oldpath = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/photoz_syst/'
    path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/tolerance_syst/'   
    #testing choice of u
    for i,zmin in enumerate([.15,.3,.45,.6]):
        config_fname = 'newpaper13.1'
        zmax = zmin+.15
        cc = chainconsumer.ChainConsumer()
        
        data1 = plottools.load_res_indep(oldpath,'dm','newpaper13.1','ZSPEC',zmin,zmax,'12x20',sigma=0)
        data2 = plottools.load_res_indep(oldpath,'newbuzzardrm2','newpaper13.1','ZSPEC',zmin,zmax,'20',sigma=0)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'$.3 < u < 1$')
        
        data1 = plottools.load_res_indep(path,'dm','newpaper13.1_newu','ZSPEC',zmin,zmax,'12x20',sigma=0)
        data2 = plottools.load_res_indep(path,'newbuzzardrm2','newpaper13.1_newu','ZSPEC',zmin,zmax,'20',sigma=0)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'$.6 < u < 1$')
        
        data1 = plottools.load_res_indep(path,'dm','newpaper13.1_newu2','ZSPEC',zmin,zmax,'12x20',sigma=0)
        data2 = plottools.load_res_indep(path,'newbuzzardrm2','newpaper13.1_newu2','ZSPEC',zmin,zmax,'20',sigma=0)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'$.9 < u < 1$')
        
        figure = cc.plotter.plot(figsize='column',extents=[(0,3),(-3,3)]);
        axarr = figure.get_axes()
        plt.suptitle(str(zmin)+r'$ < z < $'+str(zmax))
        plt.legend()
        
        if is11k:
            str11k = '_11k'
        else:
            str11k = ''

        plt.savefig('./figures/u_bin'+str(i+1)+str11k+'.pdf',dpi=300,bbox_inches='tight')
        plt.savefig('./figures/u_bin'+str(i+1)+str11k+'.png',dpi=300,bbox_inches='tight')


def get_gaussian_photoz(is11k=False):
    path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/photoz_syst/'
    summaries = []
    for i,zmin in enumerate([.15,.3,.45,.6]):
        config_fname = 'newpaper13.1'
        config = plottools.load_config(config_fname)
        zmax = zmin+.15
        cc = chainconsumer.ChainConsumer()
        data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'12x20',sigma=0)
        data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'20',sigma=0)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'$\sigma_z = 0$')
        
        data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'12x20',sigma=0.005)
        data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'20',sigma=0.005)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'$\sigma_z = 0.005$')
        
        data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'12x20',sigma=0.01)
        data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'20',sigma=0.01)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'$\sigma_z = 0.01$')
        
        data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'12x20',sigma=0.02)
        data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'20',sigma=0.02)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'$\sigma_z = 0.02$')
        
        summary = cc.analysis.get_summary()
        summaries.append(summary)
        
        figure = cc.plotter.plot(figsize='column',extents=[(0,3),(-3,3)]);
        axarr = figure.get_axes()
        plt.suptitle(str(zmin)+r'$ < z < $'+str(zmax))
        plt.legend()

        if is11k:
            str11k = '_11k'
        else:
            str11k = ''

        plt.savefig('./figures/gauss_bin'+str(i+1)+str11k+'.pdf',dpi=300,bbox_inches='tight')
        plt.savefig('./figures/gauss_bin'+str(i+1)+str11k+'.png',dpi=300,bbox_inches='tight')
        

if __name__ == '__main__':
    #make_easy_figures()
    make_data_vector_plots()
    #make_triangle_diagrams()
    make_covariance_fig()
    get_zspec_both()
    get_zspec()
    get_zspec(is11k=True)
    get_zspec_zrm()
    get_zspec_zrm(is11k=True)
    get_gaussian_photoz()
    get_gaussian_photoz(is11k=True)
    get_tolerance_figs()
    get_tolerance_figs(is11k=True)
    get_u_figs()
    get_u_figs(is11k=True)
