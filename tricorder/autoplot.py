from __future__ import print_function
from __future__ import division
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.cosmology import Planck15 as cosmo
from scipy.stats import binned_statistic
import treecorr
import sys
import numpy as np
import json
import os

outdir = '/nfs/slac/g/ki/ki19/des/mbaumer/3pt_runs/'

def plot_run(runname):
    with open(outdir+runname+'.config') as f:
        config = json.loads(f.read())
        
    if 'random_set_id' not in config.keys():
        config['random_set_id'] = 0
    if 'jk_id' not in config.keys():
        config['jk_id'] = -1
    if 'data_z_var' not in config.keys():
        config['data_z_var'] = 'ZSPEC'
    
    data = fits.getdata(config['datapath'])
    randoms = fits.getdata(config['randompath'])
    
    data = data[((data[config['data_z_var']] > config['min_z']) & (data[config['data_z_var']] < config['max_z']))]
    data = data[((data['RA'] > config['min_ra']) & (data['RA'] < config['max_ra']))]
    data = data[((data['DEC'] > config['min_dec']) & (data['DEC'] < config['max_dec']))]

    randoms = randoms[((randoms['Z'] > config['min_z']) & (randoms['Z'] < config['max_z']))]
    randoms = randoms[((randoms['RA'] > config['min_ra']) & (randoms['RA'] < config['max_ra']))]
    randoms = randoms[((randoms['DEC'] > config['min_dec']) & (randoms['DEC'] < config['max_dec']))]
    
    template = treecorr.NNNCorrelation(config=config)

    #make angular plots
    for scale in [3,6,9,12]:
        for ratio in [.5,.333,.25]:
            for tolerance in [.1,.2,.3]:
                for nbins in [8,16,24,32]:
                    mode = 'angle'
                    bins, binned = analyze_single_run(mode,scale=scale,ratio=ratio,tolerance=tolerance,nbins=nbins)
                    for name,var in binned.iteritems():
                        plt.figure()
                        plt.plot(bins,var)
                        plt.xlabel('Angle (degrees)')
                        plt.ylabel(name)
                        plt.title(str(config['min_z'])+'<'config['data_z_var']+str(config['max_z'])+'  '+str(scale*ratio)+':'+str(scale)+' +/- '+str(100*tolerance)+'%')
                        plt.figtext(.7,.7,)
                        #plt.savefig(plotdir+name+'.png')


def computeXvsAngle(ddd,var,stat='mean',scale=6,ratio=.5,tolerance=.1,nbins=15):
    transition_angle = np.arccos(.25)/np.pi*180 #angle at which elongate becomes collapsed
    N_low_bins = np.floor(transition_angle/180*nbins)
    coll_bins = np.linspace(0,transition_angle,num=N_low_bins) 
    elong_bins = np.linspace(transition_angle,180,num=nbins-N_low_bins) 

    collapsed_angles = computeAngularBins(np.exp(ddd.logr),ddd.u,ddd.v,collapsed=True)
    elongated_angles = computeAngularBins(np.exp(ddd.logr),ddd.u,ddd.v,collapsed=False)
    isRightSize = (np.exp(ddd.logr)*ddd.u > scale*ratio-scale*ratio*tolerance) & (np.exp(ddd.logr)*ddd.u < scale*ratio+scale*ratio*tolerance)
    isCollapsed = (((ddd.u*np.abs(ddd.v))*np.exp(ddd.logr)+np.exp(ddd.logr) > scale-scale*tolerance) & ((ddd.u*np.abs(ddd.v))*np.exp(ddd.logr)+np.exp(ddd.logr) < scale*scale*tolerance))
    isElongated = ((np.exp(ddd.logr) > scale-scale*tolerance) & (np.exp(ddd.logr) < scale+scale*tolerance))
    out1,bins1,_ = binned_statistic(elongated_angles[np.where(isRightSize & isElongated)],var[np.where(isRightSize & isElongated)],bins=elong_bins,statistic=stat)
    out2,bins2,_ = binned_statistic(collapsed_angles[np.where(isRightSize & isCollapsed)],var[np.where(isRightSize & isCollapsed)],bins=coll_bins,statistic=stat)
    full_var = np.concatenate((out2,out1))
    bins1 += (bins1[1]-bins1[0])/2 #make edges centers
    bins2 += (bins2[1]-bins2[0])/2
    full_bins = np.concatenate((bins2[:-1],bins1[:-1]))
    return full_var, full_bins

def compute_x_vs_side_length(ddd,var,stat='mean',nbins=15,tolerance=.1):
    isEquilateral = (ddd.u > 1-tolerance) & (np.abs(ddd.v) < tolerance)
    res, b, _ = binned_statistic(ddd.logr[isEquilateral],var[isEquilateral],bins=nbins,statistic=stat)
    b += (b[1]-b[0])/2
    b = b[:-1]
    return res, b

def analyze_single_run(mode,**kwargs):
    ddd,drr,rdr,rrd,ddr,drd,rdd,rrr = load_data_for_run(config[jk_id])
    template = treecorr.NNNCorrelation(config=config)
    
    if mode == 'angle':
        get_binned_stat = computeXvsAngle
    if mode == 'equi':
        get_binned_stat = compute_x_vs_side_length
    
    binned = {}
    binned['ddd'], bins = get_binned_stat(template,ddd,stat='sum')
    binned['ddr'], bins = get_binned_stat(template,ddr,stat='sum')
    binned['drd'], bins = get_binned_stat(template,drd,stat='sum')
    binned['rdd'], bins = get_binned_stat(template,rdd,stat='sum')
    binned['rrd'], bins = get_binned_stat(template,rrd,stat='sum')
    binned['drr'], bins = get_binned_stat(template,drr,stat='sum')
    binned['rdr'], bins = get_binned_stat(template,rdr,stat='sum')
    binned['rrr'], bins = get_binned_stat(template,rrr,stat='sum')

    binned['d1'], bins = get_binned_stat(template,template.u*np.abs(template.v)*np.exp(template.logr)+np.exp(template.logr))
    binned['d2'], bins = get_binned_stat(template,np.exp(template.logr))
    binned['d3'], bins = get_binned_stat(template,template.u*np.exp(template.logr))

    datatot = len(data)
    randtot = len(randoms)
    dddtot = float(datatot)**3/6
    drrtot = float(datatot)*float(randtot)**2/6
    rdrtot = float(datatot)*float(randtot)**2/6
    rrdtot = float(datatot)*float(randtot)**2/6
    ddrtot = float(datatot)**2*float(randtot)/6
    drdtot = float(datatot)**2*float(randtot)/6
    rddtot = float(datatot)**2*float(randtot)/6
    rrrtot = float(randtot)**3/6

    binned['zeta'] = (binned['ddd']+dddtot*(-binned['ddr']/ddrtot-binned['drd']/drdtot-binned['rdd']/rddtot+binned['rrd']/rrdtot+binned['rdr']/rdrtot+binned['drr']/drrtot-binned['rrr']/rrrtot))/(binned['rrr']*dddtot/rrrtot)  
    binned['denom'] = get_two_point_expectation(binned['d1'],binned['d2'],binned['d3'],metric=config['metric'])
    binned['q'] = binned['zeta']/binned['denom']
    return bins, binned


def get_two_point_expectation(d1bins,d2bins,d3bins):
    if config['metric'] == 'Rperp':
        cat = treecorr.Catalog(ra=data['RA'], dec=data['DEC'], ra_units='degrees', dec_units='degrees',r=cosmo.comoving_distance(data['ZSPEC'])/cosmo.h)
        random_cat = treecorr.Catalog(ra=randoms['RA'], dec=randoms['DEC'], ra_units='degrees', dec_units='degrees',r=cosmo.comoving_distance(randoms['Z'])/cosmo.h)
        dd = treecorr.NNCorrelation(min_sep=.1,max_sep=25,nbins=20,bin_slop=0.1,metric='Rperp')
        dr = treecorr.NNCorrelation(min_sep=.1,max_sep=25,nbins=20,bin_slop=0.1,metric='Rperp')
        rr = treecorr.NNCorrelation(min_sep=.1,max_sep=25,nbins=20,bin_slop=0.1,metric='Rperp')
    elif config['metric'] == 'Euclidean':
        cat = treecorr.Catalog(ra=data['RA'], dec=data['DEC'], ra_units='degrees', dec_units='degrees')
        random_cat = treecorr.Catalog(ra=randoms['RA'], dec=randoms['DEC'], ra_units='degrees', dec_units='degrees')
        dd = treecorr.NNCorrelation(min_sep=.1,max_sep=25,nbins=20,bin_slop=0.1,sep_units='arcmin',metric='Euclidean')
        dr = treecorr.NNCorrelation(min_sep=.1,max_sep=25,nbins=20,bin_slop=0.1,sep_units='arcmin',metric='Euclidean')
        rr = treecorr.NNCorrelation(min_sep=.1,max_sep=25,nbins=20,bin_slop=0.1,sep_units='arcmin',metric='Euclidean')
    else:
        raise ValueError('invalid metric specified')
    dd.process(cat)
    dr.process(cat,random_cat)
    rr.process(random_cat)
    xi, varxi = dd.calculateXi(rr=rr,dr=dr)
    
    coeffs = np.polyfit(dd.logr,np.log(xi),deg=1)
    poly = np.poly1d(coeffs)
    yfit = lambda x: np.exp(poly(np.log(x)))
    
    xi1 = yfit(d1bins)
    xi2 = yfit(d2bins)
    xi3 = yfit(d3bins)
    denom_bins = (xi1*xi2+xi2*xi3+xi3*xi1)
    return denom_bins

def load_data_for_run(jk_id):
    ddd = np.load(outdir+runname+'_'+str(random_set_id)+'_'+str(jk_id)+'_'+'ddd.npy')
    ddr = np.load(outdir+runname+'_'+str(random_set_id)+'_'+str(jk_id)+'_'+'ddr.npy')
    drd = np.load(outdir+runname+'_'+str(random_set_id)+'_'+str(jk_id)+'_'+'drd.npy')
    rdd = np.load(outdir+runname+'_'+str(random_set_id)+'_'+str(jk_id)+'_'+'rdd.npy')
    rrd = np.load(outdir+runname+'_'+str(random_set_id)+'_'+str(jk_id)+'_'+'rrd.npy')
    drr = np.load(outdir+runname+'_'+str(random_set_id)+'_'+str(jk_id)+'_'+'drr.npy')
    rdr = np.load(outdir+runname+'_'+str(random_set_id)+'_'+str(jk_id)+'_'+'rdr.npy')
    rrr = np.load(outdir+runname+'_'+str(random_set_id)+'_'+str(jk_id)+'_'+'rrr.npy')
    return ddd,drr,rdr,rrd,ddr,drd,rdd,rrr

def computeAngularBins(r,u,v,collapsed=False):
    #if v < 0: collapsed = not collapsed
    v = np.abs(v)
    d2 = r
    d3 = u*r
    d1 = v*d3+d2
    #law of cosines
    if not collapsed:
        cosine = (d2**2 + d3**2 - d1**2)/(2*d2*d3+1e-9)
    else: 
        cosine = (d1**2 + d3**2 - d2**2)/(2*d1*d3+1e-9)
    bins = np.arccos(cosine)/np.pi*180
    return bins

