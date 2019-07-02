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

def get_zspec_zrm(is11k=False):
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
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'ZSPEC')
        
        if zmin != .45:
            data1 = plottools.load_res_indep(path,'dm',config_fname,'ZREDMAGIC',zmin,zmax,'10x10',sigma=0)
            data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZREDMAGIC',zmin,zmax,'10',sigma=0)
        else:
            data1 = plottools.load_res_indep(path,'dm',config_fname,'ZREDMAGIC',zmin,zmax,'12x20',sigma=0)
            data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZREDMAGIC',zmin,zmax,'20',sigma=0)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'ZREDMAGIC')
        
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
        
if __name__ == '__main__':
    get_zspec_zrm()
    get_zspec_zrm(is11k=True)
