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
        b1 = np.sqrt(np.mean(xi2,axis=0)/np.mean(xi1,axis=0))
        r = np.logspace(config['2PCF']['min_sep'],np.log10(config['2PCF']['max_sep']),num=config['2PCF']['nbins'])
        b1_mean = np.mean(b1[(r > config['3PCF']['min_sep']*config['3PCF']['min_u']) & (r < config['3PCF']['max_sep']*config['3PCF']['max_u'])])
        b1_std = np.std(b1[(r > config['3PCF']['min_sep']*config['3PCF']['min_u']) & (r < config['3PCF']['max_sep']*config['3PCF']['max_u'])])
        
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
        plt.savefig('./figures/spec_bin'+str(i+1)+str11k+'.pdf',dpi=300,bbox_inches='tight')

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
        

if __name__ == '__main__':
    get_zspec()
    get_zspec(is11k=True)
    get_zspec_zrm()
    get_zspec_zrm(is11k=True)
    get_gaussian_photoz()
    get_gaussian_photoz(is11k=True)
    #get_tolerance_figs()
    #get_tolerance_figs(is11k=True)
    #get_u_figs()
    #get_u_figs(is11k=True)
