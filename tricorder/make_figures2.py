#### Make figures for paper (in eps/pdf format)
from __future__ import division
import matplotlib
matplotlib.use('agg')

import numpy as np
import chainconsumer
from scipy.optimize import minimize
import plottools
reload(plottools)
import palettable
path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/paper_validation/'
path2 = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/paper_validation2/'

import matplotlib.pyplot as plt
import chainconsumer
from scipy.optimize import minimize
from matplotlib.colors import LogNorm

from astropy.io import fits
import paths
import skymapper as skm

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rc('xtick', labelsize=12) 
matplotlib.rc('ytick', labelsize=12) 
matplotlib.rc('font', family='serif')

# Covmats 3 cols, 4 rows: cols are ZSPEC, Gaussian z, and ZRM
# 	Or maybe 4x4 since we might need spec gals, spec dm, rm gals, rmdm
# Contours: zspec+zrm for DES Y1 are bad
# But 11k good: Contours: zspec, zrm, 2pt, lazeyras
# Scale
# Tolerance
# U min

b1_2pt = [1.71,1.78,1.9,2.22]
b1_err_2pt = [.064,.064,.071,.074]

b1_spec = [1.51,1.71,1.71,2.1]
b1_up_spec = [.21,.17,.18,.21]
b1_down_spec = [.22,.15,.19,.19]

b2_spec = [.21,.49,.35,.7]
b2_up_spec = [.37,.29,.31,.5]
b2_down_spec = [.40,.34,.35,.43]

b1_rm = [1.41,1.87,1.69,2.33]
b1_up_rm = [.38,.37,.28,.98]
b1_down_rm = [.28,.3,.2,.54]

b2_rm = [-.18,.56,.37,1]
b2_up_rm = [.52,.83,.61,3.5]
b2_down_rm = [.42,.58,.49,1.1]

sigmas_list = [0.02,0.03,0.03,0.03]

# ## UNCOMMENT STARTING HERE
# # data vectors ZSPEC

# plt.figure()
# config_fname = 'paper_3dval4'
# for i,zmin in enumerate([.15,.3,.45,.6]):
#     sigma = 0
#     zmax = zmin+.15
#     data =  plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=sigma,use_alt_randoms=True)
#     galdata = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=sigma,use_alt_randoms=True)
#     print len(data), len(galdata)
#     plottools.plot_dv(data,'Q', offset=0.,indiv_runs=True,color='Blue',label=r'Dark Matter',compressed=True)
#     plottools.plot_dv(galdata,'Q', offset=0.,indiv_runs=True,color='Red',label=r'Galaxies',compressed=True)
#     plt.xlabel('v')
#     plt.ylabel('Q')
#     plt.title(str(zmin)+r'$<z<$'+str(zmax))
#     plt.ylim(-.75,3)
#     plt.legend()

#     plt.savefig('./figures/spec_bin'+str(i+1)+'_dv.pdf',dpi=300,bbox_inches='tight')
#     plt.savefig('./figures/spec_bin'+str(i+1)+'_dv.png',dpi=300,bbox_inches='tight')
#     plt.figure()

# # data vectors ZREDMAGIC

# config_fname = 'paper_3dval4'
# for i,zmin in enumerate([.15,.3,.45,.6]):
#     sigma = 0
#     zmax = zmin+.15
#     data =  plottools.load_res_indep(path2,'dm',config_fname,'ZREDMAGIC',zmin,zmax,'10x10',sigma=sigma,use_alt_randoms=False)
#     galdata = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZREDMAGIC',zmin,zmax,'10',sigma=sigma,use_alt_randoms=False)
#     print len(data), len(galdata)
#     plottools.plot_dv(data,'Q', offset=0.,indiv_runs=True,color='Blue',label=r'Dark Matter',compressed=True)
#     plottools.plot_dv(galdata,'Q', offset=0.,indiv_runs=True,color='Red',label=r'Galaxies',compressed=True)
#     plt.xlabel('v')
#     plt.ylabel('Q')
#     plt.title(str(zmin)+r'$<z<$'+str(zmax))
#     plt.ylim(-.75,3)
#     plt.legend()

#     plt.savefig('./figures/rm_bin'+str(i+1)+'_dv.pdf',dpi=300,bbox_inches='tight')
#     plt.savefig('./figures/rm_bin'+str(i+1)+'_dv.png',dpi=300,bbox_inches='tight')
#     plt.figure()

# #covmats

# def make_corr_fig(zmin):

#     zmax = zmin+.15
#     if zmin == .15:
#         placeholder_sigma = 0.01
#     else: 
#         placeholder_sigma = 0.02
#     data = plottools.load_res_indep(path, 'dm','paper_3dval4','ZSPEC',zmin,zmax,'10x10',sigma=0)
#     data2 = plottools.load_res_indep(path,'dm','paper_3dval4','ZSPEC',zmin,zmax,'10x10',sigma=placeholder_sigma)
#     data3 = plottools.load_res_indep(path2,'dm','paper_3dval4','ZREDMAGIC',zmin,zmax,'10x10',sigma=0)
    
#     galdata = plottools.load_res_indep(path, 'newbuzzardrm2','paper_3dval4','ZSPEC',zmin,zmax,'10',sigma=0)
#     galdata2 = plottools.load_res_indep(path,'newbuzzardrm2','paper_3dval4','ZSPEC',zmin,zmax,'10',sigma=placeholder_sigma)
#     galdata3 = plottools.load_res_indep(path,'newbuzzardrm2','paper_3dval4','ZREDMAGIC',zmin,zmax,'10',sigma=0)
#     plt.figure()
    
#     print len(data),len(data2),len(data3)
#     print len(galdata), len(galdata2), len(galdata3)
    
#     fig,axarr = plt.subplots(2,3,figsize=(8,7))
    
#     vmin= -1
#     vmax= 1
#     cmap_choice = 'viridis'
    
#     axarr[0,0].imshow(np.corrcoef(plottools.compress_dv(galdata['Q'].values.reshape(-1,10)).T),vmin=vmin,vmax=vmax,origin='lower',extent=(0,1,0,1),cmap=cmap_choice)
#     axarr[0,0].set_title(r'True redshifts')
#     axarr[0,0].set_xticks([.1,.3,.5,.7,.9])
#     axarr[0,0].set_yticks([.1,.3,.5,.7,.9])
#     #axarr[0,0].set_xlabel('v')
#     axarr[0,0].set_ylabel('v')
#     axarr[0,1].imshow(np.corrcoef(plottools.compress_dv(galdata2['Q'].values.reshape(-1,10)).T),vmin=vmin,vmax=vmax,origin='lower',extent=(0,1,0,1),cmap=cmap_choice)
#     axarr[0,1].set_title(r'Galaxies \\ \\ $\sigma_z = $'+str(placeholder_sigma))
#     #axarr[0,1].set_xlabel('v')
#     #axarr[0,1].set_ylabel('v')
#     axarr[0,1].set_xticks([.1,.3,.5,.7,.9])
#     axarr[0,1].set_yticks([.1,.3,.5,.7,.9])
    
#     axarr[0,2].imshow(np.corrcoef(plottools.compress_dv(galdata3['Q'].values.reshape(-1,10)).T),vmin=vmin,vmax=vmax,origin='lower',extent=(0,1,0,1),cmap=cmap_choice)
#     axarr[0,2].set_title(r'RedMaGiC redshifts')
#     #axarr[0,1].set_xlabel('v')
#     #axarr[0,1].set_ylabel('v')
#     axarr[0,2].set_xticks([.1,.3,.5,.7,.9])
#     axarr[0,2].set_yticks([.1,.3,.5,.7,.9])
    
#     axarr[1,0].imshow(np.corrcoef(plottools.compress_dv(data['Q'].values.reshape(-1,10)).T),vmin=vmin,vmax=vmax,origin='lower',extent=(0,1,0,1),cmap=cmap_choice)
#     axarr[1,0].set_title(r'True redshifts')
#     axarr[1,0].set_xticks([.1,.3,.5,.7,.9])
#     axarr[1,0].set_yticks([.1,.3,.5,.7,.9])
#     axarr[1,0].set_xlabel('v')
#     axarr[1,0].set_ylabel('v')
#     im = axarr[1,1].imshow(np.corrcoef(plottools.compress_dv(data2['Q'].values.reshape(-1,10)).T),vmin=vmin,vmax=vmax,origin='lower',extent=(0,1,0,1),cmap=cmap_choice)
#     axarr[1,1].set_title(r'Dark Matter \\ \\ $\sigma_z = $'+str(placeholder_sigma))
#     axarr[1,1].set_xlabel('v')
#     #axarr[1,1].set_ylabel('v')
#     axarr[1,1].set_xticks([.1,.3,.5,.7,.9])
#     axarr[1,1].set_yticks([.1,.3,.5,.7,.9])
    
#     axarr[1,2].imshow(np.corrcoef(plottools.compress_dv(data3['Q'].values.reshape(-1,10)).T),vmin=vmin,vmax=vmax,origin='lower',extent=(0,1,0,1),cmap=cmap_choice)
#     axarr[1,2].set_title(r'RedMaGiC redshifts')
#     axarr[1,2].set_xlabel('v')
#     #axarr[1,1].set_ylabel('v')
#     axarr[1,2].set_xticks([.1,.3,.5,.7,.9])
#     axarr[1,2].set_yticks([.1,.3,.5,.7,.9])
    
#     #fig.colorbar(im, ax=axarr.ravel().tolist(),shrink=.4,aspect=10)
#     #plt.suptitle(str(zmin)+r'$<z<$'+str(zmax))
#     #plt.tight_layout()
#     cbar_ax = fig.add_axes([.97, 0.2, .04, 0.6])
#     fig.colorbar(im, cax=cbar_ax, orientation="vertical")
#     plt.suptitle(str(zmin)+r'$< z <$'+str(zmax))

# for i,zmin in enumerate([.15,.3,.45,.6]):
#     make_corr_fig(zmin)
#     plt.savefig('./figures/cov_bin'+str(i+1)+'.pdf',dpi=300,bbox_inches='tight')
#     plt.savefig('./figures/cov_bin'+str(i+1)+'.png',dpi=300,bbox_inches='tight')
#     plt.figure()

# # DES Y1 contours (spec+photoz), no 2pt or lazeyras comparison

# is11k = False
# for i,zmin in enumerate([.15,.3,.45,.6]):
#     zmax = zmin+.15
#     cc = chainconsumer.ChainConsumer()
    
#     for sigma in [0]:
#         config_fname = 'paper_3dval4'
#         data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=sigma,use_alt_randoms=False)
#         data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=sigma,use_alt_randoms=False)
#         red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
#         red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
#         samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
#         cc.add_chain(samples.flatchain,parameters=['b1','b2'],name='Spectroscopic Redshifts')
        
#     data1 = plottools.load_res_indep(path2,'dm',config_fname,'ZREDMAGIC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
#     data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
#     red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
#     red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
#     samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
#     cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'RedMaGiC Photometric Redshifts')
    
#     cc.configure(legend_kwargs={"loc": "lower right"},label_font_size=14,tick_font_size=14)
#     figure = cc.plotter.plot(figsize='column',extents=[(0,10),(-10,10)]);
#     axarr = figure.get_axes()
#     figure.axes[2].axvspan(b1_2pt[i]-b1_err_2pt[i], b1_2pt[i]+b1_err_2pt[i], alpha=0.3, color='b')
#     figure.axes[2].axvspan(b1_2pt[i]-2*b1_err_2pt[i], b1_2pt[i]+2*b1_err_2pt[i], alpha=0.3, color='b')
#     plt.suptitle(str(zmin)+r'$ < z < $'+str(zmax))
    
#     plt.savefig('./figures/bothy1_bin'+str(i+1)+'_inf.pdf',dpi=300,bbox_inches='tight')
#     plt.savefig('./figures/bothy1_bin'+str(i+1)+'_inf.png',dpi=300,bbox_inches='tight')
#     plt.figure()
    
# # 11k degree contours, 2pt and lazeyras comparison.

# is11k = True
# for i,zmin in enumerate([.15,.3,.45,.6]):
#     zmax = zmin+.15
#     cc = chainconsumer.ChainConsumer()
    
#     for sigma in [0]:
#         config_fname = 'paper_3dval4'
#         data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=sigma,use_alt_randoms=False)
#         data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=sigma,use_alt_randoms=False)
#         red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
#         red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
#         samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
#         cc.add_chain(samples.flatchain,parameters=['b1','b2'],name='Spectroscopic Redshifts')

#     #replace later once RM is fixed:   
#     #data1 = plottools.load_res_indep(path2,'dm',config_fname,'ZREDMAGIC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
#     if zmin != .15:
#         data1 = plottools.load_res_indep(path2,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=sigmas_list[i],use_alt_randoms=False)
#     else:
#         data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=sigmas_list[i],use_alt_randoms=False)
        
#     data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
    
#     red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
#     red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
#     samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
#     cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'RedMaGiC Photometric Redshifts')

#     # add 2pt
#     cc.configure(legend_kwargs={"loc": "lower right"},label_font_size=14,tick_font_size=14)
#     figure = cc.plotter.plot(figsize='column',extents=[(0,4),(-4,4)]);
#     axarr = figure.get_axes()

#     figure.axes[2].axvspan(b1_2pt[i]-b1_err_2pt[i], b1_2pt[i]+b1_err_2pt[i], alpha=0.3, color='b')
#     figure.axes[2].axvspan(b1_2pt[i]-2*b1_err_2pt[i], b1_2pt[i]+2*b1_err_2pt[i], alpha=0.3, color='b')

#     figure.axes[2].plot(np.linspace(0,4,100),plottools.get_lazeyras_schmidt_2015(np.linspace(0,4,100)),label='Lazeyras+Schmidt 2015')
#     figure.axes[2].plot(np.linspace(0,4,100),plottools.get_hoffman_2015(np.linspace(0,4,100)),label='Hoffman 2015')

#     plt.suptitle(str(zmin)+r'$ < z < $'+str(zmax))
    
#     plt.savefig('./figures/both11k_bin'+str(i+1)+'_inf.pdf',dpi=300,bbox_inches='tight')
#     plt.savefig('./figures/both11k_bin'+str(i+1)+'_inf.png',dpi=300,bbox_inches='tight')
#     plt.figure()
    
# # tol dv redmagic
# colors = palettable.colorbrewer.sequential.Reds_3.hex_colors
# colors2 = palettable.colorbrewer.sequential.Blues_3.hex_colors
# for i,zmin in enumerate([.15,.3,.45,.6]):
#     zmax = zmin+.15
#     galdata = plottools.load_res_indep(path2,'newbuzzardrm2','fiducial3d_toldown','ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
#     galdata1 = plottools.load_res_indep(path,'newbuzzardrm2','paper_3dval4','ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
#     galdata2 = plottools.load_res_indep(path2,'newbuzzardrm2','fiducial3d_tolup','ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
    
#     data = plottools.load_res_indep(path2,'dm','fiducial3d_toldown','ZREDMAGIC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
#     data1 = plottools.load_res_indep(path2,'dm','paper_3dval4','ZREDMAGIC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
#     data2 = plottools.load_res_indep(path2,'dm','fiducial3d_tolup','ZREDMAGIC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
    
#     plottools.plot_dv(data,'Q', offset=0.,indiv_runs=False,color=colors2[0]   ,compressed=True)
#     plottools.plot_dv(data1,'Q', offset=0.01,indiv_runs=False,color=colors2[1],compressed=True)
#     plottools.plot_dv(data2,'Q', offset=0.02,indiv_runs=False,color=colors2[2],compressed=True)
    
    
#     plottools.plot_dv(galdata,'Q', offset=0.,indiv_runs=False,color=colors[0]   ,label=r'Tolerance = $\pm 2.5 Mpc$' ,compressed=True)
#     plottools.plot_dv(galdata1,'Q', offset=0.01,indiv_runs=False,color=colors[1],label=r'Tolerance = $\pm 5 Mpc$' ,compressed=True)
#     plottools.plot_dv(galdata2,'Q', offset=0.02,indiv_runs=False,color=colors[2],label=r'Tolerance = $\pm 7.5 Mpc$' ,compressed=True)
    
#     plt.xlabel('v')
#     plt.ylabel('Q')
#     plt.title(str(zmin)+r'$<z<$'+str(zmax))
#     plt.legend()
#     plt.savefig('./figures/tol_bin'+str(i+1)+'_dvrm.pdf',dpi=300,bbox_inches='tight')
#     plt.savefig('./figures/tol_bin'+str(i+1)+'_dvrm.png',dpi=300,bbox_inches='tight')
#     plt.figure()

    
# # tol inference spec

# is11k = True
# tol_labels = ['2.5','5','7.5']
# for i,zmin in enumerate([.15,.3,.45,.6]):
#     # delete when jobs are finished
#     if zmin == .45: continue
#     zmax = zmin+.15
#     cc = chainconsumer.ChainConsumer()
    
#     for j,config_fname in enumerate(['fiducial3d_toldown','paper_3dval4','fiducial3d_tolup']):
#         data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
#         data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
#         red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
#         red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
#         samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
#         cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'Tolerance = $\pm$ '+tol_labels[j]+' Mpc')
    
#     cc.configure(legend_kwargs={"loc": "lower right"},label_font_size=14,tick_font_size=14)
#     figure = cc.plotter.plot(figsize='column',extents=[(0,4),(-4,4)]);
#     axarr = figure.get_axes()
#     plt.suptitle(str(zmin)+r'$ < z < $'+str(zmax))
    
#     plt.savefig('./figures/tol_bin'+str(i+1)+'_inf.pdf',dpi=300,bbox_inches='tight')
#     plt.savefig('./figures/tol_bin'+str(i+1)+'_inf.png',dpi=300,bbox_inches='tight')
#     plt.figure()
    

# # tol inference ZREDMAGIC

# is11k = True
# for i,zmin in enumerate([.15,.3,.45,.6]):
#     zmax = zmin + 0.15
#     cc = chainconsumer.ChainConsumer()
    
#     data1 = plottools.load_res_indep(path2,'dm','fiducial3d_toldown','ZREDMAGIC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
#     data2 = plottools.load_res_indep(path2,'newbuzzardrm2','fiducial3d_toldown','ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
    
#     red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
#     red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
#     samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
#     cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'Tolerance = $\pm 2.5 $ Mpc')
    
#     data1 = plottools.load_res_indep(path2,'dm','paper_3dval4','ZREDMAGIC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
#     data2 = plottools.load_res_indep(path,'newbuzzardrm2','paper_3dval4','ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
    
#     red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
#     red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
#     samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
#     cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'Tolerance = $\pm$ 5 Mpc')
    
#     data1 = plottools.load_res_indep(path2,'dm','fiducial3d_tolup','ZREDMAGIC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
#     data2 = plottools.load_res_indep(path2,'newbuzzardrm2','fiducial3d_tolup','ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
    
#     red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
#     red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
#     samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
#     cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'Tolerance = $\pm$ 7.5 Mpc')
    
#     cc.configure(legend_kwargs={"loc": "lower right"},label_font_size=14,tick_font_size=14)
#     figure = cc.plotter.plot(figsize='column',extents=[(0,4),(-4,4)]);
#     axarr = figure.get_axes()
#     plt.suptitle(str(zmin)+r'$ < z < $'+str(zmax))

#     plt.savefig('./figures/tol_bin'+str(i+1)+'_infrm.pdf',dpi=300,bbox_inches='tight')
#     plt.savefig('./figures/tol_bin'+str(i+1)+'_infrm.png',dpi=300,bbox_inches='tight')
#     plt.figure()

# # u data vectors

# colors = palettable.colorbrewer.sequential.Reds_3.hex_colors
# colors2 = palettable.colorbrewer.sequential.Blues_3.hex_colors
# for i,zmin in enumerate([.15,.3,.45,.6]):
#     plt.figure()
#     zmax = zmin+.15
#     galdata = plottools.load_res_indep(path,'newbuzzardrm2','fiducial3d_halfu','ZSPEC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
#     galdata1 = plottools.load_res_indep(path,'newbuzzardrm2','fiducial3d_75u','ZSPEC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
#     galdata2 = plottools.load_res_indep(path,'newbuzzardrm2','paper_3dval4','ZSPEC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
    
#     data = plottools.load_res_indep(path,'dm','fiducial3d_halfu','ZSPEC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
#     data1 = plottools.load_res_indep(path,'dm','fiducial3d_75u','ZSPEC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
#     data2 = plottools.load_res_indep(path,'dm','paper_3dval4','ZSPEC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
    
#     plottools.plot_dv(data,'Q', offset=0.,indiv_runs=False,color=colors2[0]    ,compressed=True)
#     plottools.plot_dv(data1,'Q', offset=0.01,indiv_runs=False,color=colors2[1],compressed=True)
#     plottools.plot_dv(data2,'Q', offset=0.02,indiv_runs=False,color=colors2[2],compressed=True)
    
    
#     plottools.plot_dv(galdata,'Q', offset=0.,indiv_runs=False,color=colors[0]   ,label=r'$u > .5$' ,compressed=True)
#     plottools.plot_dv(galdata1,'Q', offset=0.01,indiv_runs=False,color=colors[1],label=r'$u > .75$' ,compressed=True)
#     plottools.plot_dv(galdata2,'Q', offset=0.02,indiv_runs=False,color=colors[2],label=r'$u > .9$' ,compressed=True)
    
#     plt.xlabel('v')
#     plt.ylabel('Q')
#     plt.title(str(zmin)+r'$<z<$'+str(zmax))
#     plt.legend()
#     plt.savefig('./figures/u_bin'+str(i+1)+'_dv.pdf',dpi=300,bbox_inches='tight')
#     plt.savefig('./figures/u_bin'+str(i+1)+'_dv.png',dpi=300,bbox_inches='tight')
#     plt.figure()

# # u data vectors ZREDMAGIC

# colors = palettable.colorbrewer.sequential.Reds_3.hex_colors
# colors2 = palettable.colorbrewer.sequential.Blues_3.hex_colors
# for i,zmin in enumerate([.15,.3,.45,.6]):
#     plt.figure()
#     zmax = zmin+.15
#     galdata = plottools.load_res_indep(path2,'newbuzzardrm2','fiducial3d_halfu','ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
#     galdata1 = plottools.load_res_indep(path2,'newbuzzardrm2','fiducial3d_75u','ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
#     galdata2 = plottools.load_res_indep(path,'newbuzzardrm2','paper_3dval4','ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
    
#     data = plottools.load_res_indep(path2,'dm','fiducial3d_halfu','ZREDMAGIC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
#     data1 = plottools.load_res_indep(path2,'dm','fiducial3d_75u', 'ZREDMAGIC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
#     data2 = plottools.load_res_indep(path2,'dm','paper_3dval4',   'ZREDMAGIC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
    
#     plottools.plot_dv(data,'Q', offset=0.,indiv_runs=False,color=colors2[0]    ,compressed=True)
#     plottools.plot_dv(data1,'Q', offset=0.01,indiv_runs=False,color=colors2[1],compressed=True)
#     plottools.plot_dv(data2,'Q', offset=0.02,indiv_runs=False,color=colors2[2],compressed=True)
    
    
#     plottools.plot_dv(galdata,'Q', offset=0.,indiv_runs=False,color=colors[0]   ,label=r'$u > .5$' ,compressed=True)
#     plottools.plot_dv(galdata1,'Q', offset=0.01,indiv_runs=False,color=colors[1],label=r'$u > .75$' ,compressed=True)
#     plottools.plot_dv(galdata2,'Q', offset=0.02,indiv_runs=False,color=colors[2],label=r'$u > .9$' ,compressed=True)
    
#     plt.xlabel('v')
#     plt.ylabel('Q')
#     plt.title(str(zmin)+r'$<z<$'+str(zmax))
#     plt.legend()
#     plt.savefig('./figures/u_bin'+str(i+1)+'_dvrm.pdf',dpi=300,bbox_inches='tight')
#     plt.savefig('./figures/u_bin'+str(i+1)+'_dvrm.png',dpi=300,bbox_inches='tight')
#     plt.figure()

# # u inferences SPEC

# is11k = True
# tol_labels = ['.5','.75','.9']
# for i,zmin in enumerate([.15,.3,.45,.6]):
#     zmax = zmin+.15
#     cc = chainconsumer.ChainConsumer()
    
#     for j,config_fname in enumerate(['fiducial3d_halfu','fiducial3d_75u','paper_3dval4']):
#         data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
#         data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
#         red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
#         red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
#         samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
#         cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'$u >$ '+tol_labels[j])
    
#     cc.configure(legend_kwargs={"loc": "lower right"},label_font_size=14,tick_font_size=14)
#     figure = cc.plotter.plot(figsize='column',extents=[(0,4),(-4,4)]);
#     axarr = figure.get_axes()
#     plt.suptitle(str(zmin)+r'$ < z < $'+str(zmax))
    
#     plt.savefig('./figures/u_bin'+str(i+1)+'_inf.pdf',dpi=300,bbox_inches='tight')
#     plt.savefig('./figures/u_bin'+str(i+1)+'_inf.png',dpi=300,bbox_inches='tight')
#     plt.figure()
    
# # u inferences ZRM

# is11k = True
# tol_labels = ['.5','.75','.9']
# for i,zmin in enumerate([.15,.3,.45,.6]):
#     zmax = zmin+.15
#     cc = chainconsumer.ChainConsumer()
    
#     for j,config_fname in enumerate(['fiducial3d_halfu','fiducial3d_75u','paper_3dval4']):
#         data1 = plottools.load_res_indep(path2,'dm',config_fname,'ZREDMAGIC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
#         if config_fname == 'paper_3dval4':
#             data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
#         else:
#             data2 = plottools.load_res_indep(path2,'newbuzzardrm2',config_fname,'ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
#         red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
#         red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
#         samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
#         cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'$u >$ '+tol_labels[j])
    
#     cc.configure(legend_kwargs={"loc": "lower right"},label_font_size=14,tick_font_size=14)
#     figure = cc.plotter.plot(figsize='column',extents=[(0,4),(-4,4)]);
#     axarr = figure.get_axes()
#     plt.suptitle(str(zmin)+r'$ < z < $'+str(zmax))
#     plt.savefig('./figures/u_bin'+str(i+1)+'_infrm.pdf',dpi=300,bbox_inches='tight')
#     plt.savefig('./figures/u_bin'+str(i+1)+'_infrm.png',dpi=300,bbox_inches='tight')
#     plt.figure()

# # gaussian z dv

# colors = palettable.colorbrewer.sequential.Reds_7.hex_colors
# config_fname = 'paper_3dval4'
# zmin = .15
# plt.figure()
# zmax = zmin+.15
# galdata = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
# galdata1 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.001,use_alt_randoms=False)
# galdata2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.003,use_alt_randoms=False)
# galdata3 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.005,use_alt_randoms=False)
# galdata4 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.01,use_alt_randoms=False)
# galdata8 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)

# plottools.plot_dv(galdata,'Q', offset=0.,indiv_runs=False,color=colors[0]   ,label=r'$\sigma(z) = 0$' ,compressed=True)
# plottools.plot_dv(galdata1,'Q', offset=0.01,indiv_runs=False,color=colors[1],label=r'$\sigma(z) = 0.001$' ,compressed=True)
# plottools.plot_dv(galdata2,'Q', offset=0.02,indiv_runs=False,color=colors[2],label=r'$\sigma(z) = 0.003$' ,compressed=True)
# plottools.plot_dv(galdata3,'Q', offset=0.03,indiv_runs=False,color=colors[3],label=r'$\sigma(z) = 0.005$' ,compressed=True)
# plottools.plot_dv(galdata4,'Q', offset=0.04,indiv_runs=False,color=colors[4],label=r'$\sigma(z) = 0.01$' ,compressed=True)
# plottools.plot_dv(galdata8,'Q', offset=0.05,indiv_runs=False,color='k',label='ZREDMAGIC',compressed=True)

# plt.xlabel('v')
# plt.ylabel('Q')
# plt.title(str(zmin)+r'$<z<$'+str(zmax))
# plt.legend()
# plt.savefig('./figures/gaussrm_bin1_galdv.pdf',dpi=300,bbox_inches='tight')
# plt.savefig('./figures/gaussrm_bin1_galdv.png',dpi=300,bbox_inches='tight')
# plt.figure()

# zmin = .3
# plt.figure()
# zmax = zmin+.15
# galdata = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
# galdata1 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.001,use_alt_randoms=False)
# galdata2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.003,use_alt_randoms=False)
# galdata3 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.005,use_alt_randoms=False)
# galdata4 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.01,use_alt_randoms=False)
# galdata5 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.02,use_alt_randoms=False)
# galdata8 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)

# plottools.plot_dv(galdata,'Q', offset=0.,indiv_runs=False,color=colors[0]   ,label=r'$\sigma(z) = 0$' ,compressed=True)
# plottools.plot_dv(galdata1,'Q', offset=0.01,indiv_runs=False,color=colors[1],label=r'$\sigma(z) = 0.001$' ,compressed=True)
# plottools.plot_dv(galdata2,'Q', offset=0.02,indiv_runs=False,color=colors[2],label=r'$\sigma(z) = 0.003$' ,compressed=True)
# plottools.plot_dv(galdata3,'Q', offset=0.03,indiv_runs=False,color=colors[3],label=r'$\sigma(z) = 0.005$' ,compressed=True)
# plottools.plot_dv(galdata4,'Q', offset=0.04,indiv_runs=False,color=colors[4],label=r'$\sigma(z) = 0.01$' ,compressed=True)
# plottools.plot_dv(galdata5,'Q', offset=0.05,indiv_runs=False,color=colors[5],label=r'$\sigma(z) = 0.02$' ,compressed=True)
# plottools.plot_dv(galdata8,'Q', offset=0.06,indiv_runs=False,color='k',label='ZREDMAGIC',compressed=True)

# plt.xlabel('v')
# plt.ylabel('Q')
# plt.title(str(zmin)+r'$<z<$'+str(zmax))
# plt.legend()
# plt.savefig('./figures/gaussrm_bin2_galdv.pdf',dpi=300,bbox_inches='tight')
# plt.savefig('./figures/gaussrm_bin2_galdv.png',dpi=300,bbox_inches='tight')
# plt.figure()

# zmin = .45
# plt.figure()
# zmax = zmin+.15
# galdata = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
# galdata1 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.001,use_alt_randoms=False)
# galdata2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.003,use_alt_randoms=False)
# galdata3 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.005,use_alt_randoms=False)
# galdata4 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.01,use_alt_randoms=False)
# galdata5 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.02,use_alt_randoms=False)
# galdata8 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)

# plottools.plot_dv(galdata,'Q', offset=0.,indiv_runs=False,color=colors[0]   ,label=r'$\sigma(z) = 0$' ,compressed=True)
# plottools.plot_dv(galdata1,'Q', offset=0.01,indiv_runs=False,color=colors[1],label=r'$\sigma(z) = 0.001$' ,compressed=True)
# plottools.plot_dv(galdata2,'Q', offset=0.02,indiv_runs=False,color=colors[2],label=r'$\sigma(z) = 0.003$' ,compressed=True)
# plottools.plot_dv(galdata3,'Q', offset=0.03,indiv_runs=False,color=colors[3],label=r'$\sigma(z) = 0.005$' ,compressed=True)
# plottools.plot_dv(galdata4,'Q', offset=0.04,indiv_runs=False,color=colors[4],label=r'$\sigma(z) = 0.01$' ,compressed=True)
# plottools.plot_dv(galdata5,'Q', offset=0.05,indiv_runs=False,color=colors[5],label=r'$\sigma(z) = 0.02$' ,compressed=True)
# plottools.plot_dv(galdata8,'Q', offset=0.06,indiv_runs=False,color='k',label='ZREDMAGIC',compressed=True)

# plt.xlabel('v')
# plt.ylabel('Q')
# plt.title(str(zmin)+r'$<z<$'+str(zmax))
# plt.legend()
# plt.savefig('./figures/gaussrm_bin3_galdv.pdf',dpi=300,bbox_inches='tight')
# plt.savefig('./figures/gaussrm_bin3_galdv.png',dpi=300,bbox_inches='tight')

# plt.figure()

# zmin = .6
# plt.figure()
# zmax = zmin+.15
# galdata = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
# galdata1 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.001,use_alt_randoms=False)
# galdata2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.003,use_alt_randoms=False)
# galdata3 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.005,use_alt_randoms=False)
# galdata4 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.01,use_alt_randoms=False)
# galdata5 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.02,use_alt_randoms=False)
# galdata6 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.03,use_alt_randoms=False)
# galdata8 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)

# plottools.plot_dv(galdata,'Q', offset=0.,indiv_runs=False,color=colors[0]   ,label=r'$\sigma(z) = 0$' ,compressed=True)
# plottools.plot_dv(galdata1,'Q', offset=0.01,indiv_runs=False,color=colors[1],label=r'$\sigma(z) = 0.001$' ,compressed=True)
# plottools.plot_dv(galdata2,'Q', offset=0.02,indiv_runs=False,color=colors[2],label=r'$\sigma(z) = 0.003$' ,compressed=True)
# plottools.plot_dv(galdata3,'Q', offset=0.03,indiv_runs=False,color=colors[3],label=r'$\sigma(z) = 0.005$' ,compressed=True)
# plottools.plot_dv(galdata4,'Q', offset=0.04,indiv_runs=False,color=colors[4],label=r'$\sigma(z) = 0.01$' ,compressed=True)
# plottools.plot_dv(galdata5,'Q', offset=0.05,indiv_runs=False,color=colors[5],label=r'$\sigma(z) = 0.02$' ,compressed=True)
# plottools.plot_dv(galdata6,'Q', offset=0.06,indiv_runs=False,color=colors[6],label=r'$\sigma(z) = 0.03$' ,compressed=True)
# plottools.plot_dv(galdata8,'Q', offset=0.07,indiv_runs=False,color='k',label='ZREDMAGIC',compressed=True)

# plt.xlabel('v')
# plt.ylabel('Q')
# plt.title(str(zmin)+r'$<z<$'+str(zmax))
# plt.legend()

# plt.savefig('./figures/gaussrm_bin4_galdv.pdf',dpi=300,bbox_inches='tight')
# plt.savefig('./figures/gaussrm_bin4_galdv.png',dpi=300,bbox_inches='tight')
# plt.figure()

# # dv data vector figures

# colors = palettable.colorbrewer.sequential.Blues_7.hex_colors
# config_fname = 'paper_3dval4'
# zmin = .15
# plt.figure()
# zmax = zmin+.15
# galdata = plottools.load_res_indep(path, 'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
# galdata1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.001,use_alt_randoms=False)
# galdata2 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.003,use_alt_randoms=False)
# galdata3 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.005,use_alt_randoms=False)
# galdata4 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.01,use_alt_randoms=False)
# galdata8 = plottools.load_res_indep(path2,'dm',config_fname,'ZREDMAGIC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)

# plottools.plot_dv(galdata,'Q', offset=0.,indiv_runs=False,color=colors[0]   , label=r'$\sigma(z) = 0$',compressed=True)
# plottools.plot_dv(galdata1,'Q', offset=0.01,indiv_runs=False,color=colors[1], label=r'$\sigma(z) = 0.001$',compressed=True)
# plottools.plot_dv(galdata2,'Q', offset=0.02,indiv_runs=False,color=colors[2], label=r'$\sigma(z) = 0.003$',compressed=True)
# plottools.plot_dv(galdata3,'Q', offset=0.03,indiv_runs=False,color=colors[3], label=r'$\sigma(z) = 0.005$',compressed=True)
# plottools.plot_dv(galdata4,'Q', offset=0.04,indiv_runs=False,color=colors[4], label=r'$\sigma(z) = 0.01$',compressed=True)
# plottools.plot_dv(galdata8,'Q', offset=0.05,indiv_runs=False,color='k', label=r'ZREDMAGIC',compressed=True)

# plt.xlabel('v')
# plt.ylabel('Q')
# plt.title(str(zmin)+r'$<z<$'+str(zmax))
# plt.legend()
# plt.savefig('./figures/gaussrm_bin1_dmdv.pdf',dpi=300,bbox_inches='tight')
# plt.savefig('./figures/gaussrm_bin1_dmdv.png',dpi=300,bbox_inches='tight')
# plt.figure()

# zmin = .3
# plt.figure()
# zmax = zmin+.15
# galdata = plottools.load_res_indep(path, 'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
# galdata1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.001,use_alt_randoms=False)
# galdata2 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.003,use_alt_randoms=False)
# galdata3 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.005,use_alt_randoms=False)
# galdata4 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.01,use_alt_randoms=False)
# galdata5 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.02,use_alt_randoms=False)
# galdata8 = plottools.load_res_indep(path2,'dm',config_fname,'ZREDMAGIC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)

# plottools.plot_dv(galdata,'Q', offset=0.,indiv_runs=False,color=colors[0]   , label=r'$\sigma(z) = 0$',compressed=True)
# plottools.plot_dv(galdata1,'Q', offset=0.01,indiv_runs=False,color=colors[1], label=r'$\sigma(z) = 0.001$',compressed=True)
# plottools.plot_dv(galdata2,'Q', offset=0.02,indiv_runs=False,color=colors[2], label=r'$\sigma(z) = 0.003$',compressed=True)
# plottools.plot_dv(galdata3,'Q', offset=0.03,indiv_runs=False,color=colors[3], label=r'$\sigma(z) = 0.005$',compressed=True)
# plottools.plot_dv(galdata4,'Q', offset=0.04,indiv_runs=False,color=colors[4], label=r'$\sigma(z) = 0.01$',compressed=True)
# plottools.plot_dv(galdata5,'Q', offset=0.05,indiv_runs=False,color=colors[5], label=r'$\sigma(z) = 0.02$',compressed=True)
# plottools.plot_dv(galdata8,'Q', offset=0.06,indiv_runs=False,color='k', label=r'ZREDMAGIC',compressed=True)

# plt.xlabel('v')
# plt.ylabel('Q')
# plt.title(str(zmin)+r'$<z<$'+str(zmax))
# plt.legend()
# plt.savefig('./figures/gaussrm_bin2_dmdv.pdf',dpi=300,bbox_inches='tight')
# plt.savefig('./figures/gaussrm_bin2_dmdv.png',dpi=300,bbox_inches='tight')
# plt.figure()

# zmin = .45
# plt.figure()
# zmax = zmin+.15
# galdata = plottools.load_res_indep(path, 'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
# galdata1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.001,use_alt_randoms=False)
# galdata2 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.003,use_alt_randoms=False)
# galdata3 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.005,use_alt_randoms=False)
# galdata4 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.01,use_alt_randoms=False)
# galdata5 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.02,use_alt_randoms=False)
# galdata8 = plottools.load_res_indep(path2,'dm',config_fname,'ZREDMAGIC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)

# plottools.plot_dv(galdata,'Q', offset=0.,indiv_runs=False,color=colors[0]   , label=r'$\sigma(z) = 0$',compressed=True)
# plottools.plot_dv(galdata1,'Q', offset=0.01,indiv_runs=False,color=colors[1], label=r'$\sigma(z) = 0.001$',compressed=True)
# plottools.plot_dv(galdata2,'Q', offset=0.02,indiv_runs=False,color=colors[2], label=r'$\sigma(z) = 0.003$',compressed=True)
# plottools.plot_dv(galdata3,'Q', offset=0.03,indiv_runs=False,color=colors[3], label=r'$\sigma(z) = 0.005$',compressed=True)
# plottools.plot_dv(galdata4,'Q', offset=0.04,indiv_runs=False,color=colors[4], label=r'$\sigma(z) = 0.01$',compressed=True)
# plottools.plot_dv(galdata5,'Q', offset=0.05,indiv_runs=False,color=colors[5], label=r'$\sigma(z) = 0.02$',compressed=True)
# plottools.plot_dv(galdata8,'Q', offset=0.06,indiv_runs=False,color='k', label=r'ZREDMAGIC',compressed=True)

# plt.xlabel('v')
# plt.ylabel('Q')
# plt.title(str(zmin)+r'$<z<$'+str(zmax))
# plt.legend()
# plt.savefig('./figures/gaussrm_bin3_dmdv.pdf',dpi=300,bbox_inches='tight')
# plt.savefig('./figures/gaussrm_bin3_dmdv.png',dpi=300,bbox_inches='tight')
# plt.figure()

# zmin = .6
# plt.figure()
# zmax = zmin+.15
# galdata = plottools.load_res_indep(path, 'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
# galdata1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.001,use_alt_randoms=False)
# galdata2 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.003,use_alt_randoms=False)
# galdata3 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.005,use_alt_randoms=False)
# galdata4 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.01,use_alt_randoms=False)
# galdata5 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.02,use_alt_randoms=False)
# galdata6 = plottools.load_res_indep(path2,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.03,use_alt_randoms=False)
# galdata8 = plottools.load_res_indep(path2,'dm',config_fname,'ZREDMAGIC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)

# plottools.plot_dv(galdata,'Q', offset=0.,indiv_runs=False,color=colors[0]   , label=r'$\sigma(z) = 0$',compressed=True)
# plottools.plot_dv(galdata1,'Q', offset=0.01,indiv_runs=False,color=colors[1], label=r'$\sigma(z) = 0.001$',compressed=True)
# plottools.plot_dv(galdata2,'Q', offset=0.02,indiv_runs=False,color=colors[2], label=r'$\sigma(z) = 0.003$',compressed=True)
# plottools.plot_dv(galdata3,'Q', offset=0.03,indiv_runs=False,color=colors[3], label=r'$\sigma(z) = 0.005$',compressed=True)
# plottools.plot_dv(galdata4,'Q', offset=0.04,indiv_runs=False,color=colors[4], label=r'$\sigma(z) = 0.01$',compressed=True)
# plottools.plot_dv(galdata5,'Q', offset=0.05,indiv_runs=False,color=colors[5], label=r'$\sigma(z) = 0.02$',compressed=True)
# plottools.plot_dv(galdata6,'Q', offset=0.06,indiv_runs=False,color=colors[6], label=r'$\sigma(z) = 0.03$',compressed=True)
# plottools.plot_dv(galdata8,'Q', offset=0.07,indiv_runs=False,color='k', label=r'ZREDMAGIC',compressed=True)

# plt.xlabel('v')
# plt.ylabel('Q')
# plt.title(str(zmin)+r'$<z<$'+str(zmax))
# plt.legend()
# plt.savefig('./figures/gaussrm_bin4_dmdv.pdf',dpi=300,bbox_inches='tight')
# plt.savefig('./figures/gaussrm_bin4_dmdv.png',dpi=300,bbox_inches='tight')
# plt.figure()

#Gaussian photoz + rm inference

is11k = True
for i,zmin in enumerate([.15,.3,.45,.6]):
    zmax = zmin+.15
    cc = chainconsumer.ChainConsumer()
    
    config_fname = 'paper_3dval4'
    data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
    data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
    red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
    red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
    samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
    cc.add_chain(samples.flatchain,parameters=['b1','b2'],name='Spectroscopic Redshifts',color='k',linewidth=3)
    
    data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.005,use_alt_randoms=False)
    data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.005,use_alt_randoms=False)
    red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
    red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
    samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
    cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'$\sigma(z) = 0.005$')
    
    data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.01,use_alt_randoms=False)
    data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.01,use_alt_randoms=False)
    red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
    red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
    samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
    cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'$\sigma(z) = 0.01$')
    
    if zmin > 0.15:
        data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.02,use_alt_randoms=False)
        data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.02,use_alt_randoms=False)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'$\sigma(z) = 0.02$')
    
    if zmin > .45:
        data1 = plottools.load_res_indep(path2,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0.03,use_alt_randoms=False)
        data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0.03,use_alt_randoms=False)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'$\sigma(z) = 0.03$')
        
    if zmin != .15:
        data1 = plottools.load_res_indep(path2,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=sigmas_list[i],use_alt_randoms=False)
    else:
        data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=sigmas_list[i],use_alt_randoms=False)
    
    #fix to use RM DM later...
    #data1 = plottools.load_res_indep(path2,'dm',config_fname,'ZREDMAGIC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
    data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZREDMAGIC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
    red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
    red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
    samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
    cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'RedMaGiC Photometric Redshifts',linewidth=3)
    
    cc.configure(legend_kwargs={"loc": "lower right"},label_font_size=14,tick_font_size=14)
    figure = cc.plotter.plot(figsize='column',extents=[(0,4),(-4,4)]);
    axarr = figure.get_axes()
    plt.suptitle(str(zmin)+r'$ < z < $'+str(zmax))
    
    plt.savefig('./figures/gaussrm_bin'+str(i+1)+'_inf.pdf',dpi=300,bbox_inches='tight')
    plt.savefig('./figures/gaussrm_bin'+str(i+1)+'_inf.png',dpi=300,bbox_inches='tight')

    plt.figure()

# scale inferences

is11k = True
for i,zmin in enumerate([.15,.3,.45,.6]):
    zmax = zmin+.15
    cc = chainconsumer.ChainConsumer()
    
    config_fname = 'paper_3dval1'
    data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
    data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
    red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
    red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
    samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
    cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'4 Mpc')
    
    config_fname = 'paper_3dval2'
    data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
    data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
    red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
    red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
    samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
    cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'10 Mpc')
    
    config_fname = 'paper_3dval3'
    data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
    data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
    red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
    red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
    samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
    cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'20 Mpc',color='k',linestyle='--')
    
    if zmin == .6:
        config_fname = 'fiducial3d_25Mpc'
        data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
        data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'25 Mpc')

    config_fname = 'paper_3dval4'
    data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
    data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
    print len(data1),len(data2)
    red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
    red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
    samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
    cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'30 Mpc',color='k')
    
    if zmin != .45:
        config_fname = 'fiducial3d_35Mpc'
        data1 = plottools.load_res_indep(path,'dm',config_fname,'ZSPEC',zmin,zmax,'10x10',sigma=0,use_alt_randoms=False)
        data2 = plottools.load_res_indep(path,'newbuzzardrm2',config_fname,'ZSPEC',zmin,zmax,'10',sigma=0,use_alt_randoms=False)
        print len(data1),len(data2)
        red_qdm = plottools.compress_dv(data1['Q'].values.reshape(-1,10))
        red_qrm = plottools.compress_dv(data2['Q'].values.reshape(-1,10))
        samples = plottools.make_inference(red_qdm,red_qrm,is11k=is11k)
        cc.add_chain(samples.flatchain,parameters=['b1','b2'],name=r'35 Mpc')
    
    cc.configure(legend_kwargs={"loc": "lower right"},label_font_size=14,tick_font_size=14)
    figure = cc.plotter.plot(figsize='column',extents=[(0,10),(-10,10)]);
    axarr = figure.get_axes()
    plt.suptitle(str(zmin)+r'$ < z < $'+str(zmax))
    
    plt.savefig('./figures/scale_bin'+str(i+1)+'_inf.pdf',dpi=300,bbox_inches='tight')
    plt.savefig('./figures/scale_bin'+str(i+1)+'_inf.png',dpi=300,bbox_inches='tight')
    plt.figure()

z = np.array([.15,.3,.45,.6])+0.075
line1 = plt.errorbar(z-0.01,b1_spec,yerr=np.array([b1_down_spec,b1_up_spec]),linestyle='None',marker='o',color='b',label='Spectroscopic Redshifts')
line2 = plt.errorbar(z+0.01,b1_rm,yerr=np.array([b1_down_rm,b1_up_rm]),linestyle='None',marker='o',color='r',label='RedMaGiC Redshifts')

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Patch

errorboxes = []

for i in [0,1,2,3]:
    rect = Rectangle((z[i] - 0.02, b1_2pt[i] - b1_err_2pt[i]), 0.04, 2*b1_err_2pt[i])
    errorboxes.append(rect)

pc = PatchCollection(errorboxes, facecolor='k', alpha=0.5,
                         edgecolor='None')
ax = plt.gca()
ax.add_collection(pc)

red_patch = Patch(color='k', alpha=0.5, label='Linear bias from 2PCFs')
plt.legend(handles=[line1,line2,red_patch],loc=2)

plt.xlabel('z')
plt.ylabel('Linear bias')
plt.xlim(.15,.75)

plt.savefig('./figures/b1_evolution.pdf',dpi=300,bbox_inches='tight')
plt.savefig('./figures/b1_evolution.png',dpi=300,bbox_inches='tight')
plt.figure()

plt.errorbar(z-0.01,b2_spec,yerr=np.array([b2_down_spec,b2_up_spec]),linestyle='None',marker='o',color='b',label='Spectroscopic Redshifts')
plt.errorbar(z+0.01,b2_rm,yerr=np.array([b2_down_rm,b2_up_rm]),linestyle='None',marker='o',color='r',label='RedMaGiC Redshifts')
plt.hlines(0,0.2,0.7,linestyles=['--'],colors=['k'])

plt.xlabel('z')
plt.ylabel('Quadratic bias')
plt.legend(loc=2)
plt.xlim(.15,.75)
plt.savefig('./figures/b2_evolution.pdf',dpi=300,bbox_inches='tight')
plt.savefig('./figures/b2_evolution.png',dpi=300,bbox_inches='tight')