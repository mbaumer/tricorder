from __future__ import print_function
from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
buzzard_cosmo = FlatLambdaCDM(68.81,.295)
from scipy.stats import binned_statistic
import subprocess
import pandas as pd
import treecorr
import sys
import numpy as np
import yaml
import os

outdir = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new_triplet_counts/'
plotdir = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/plots/'

def computeXvsAngle(ddd,var,stat='mean',scale=6,ratio=.5,tolerance=.1,nbins=15,**kwargs):
    transition_angle = np.arccos(.25)/np.pi*180 #angle at which elongate becomes collapsed
    N_low_bins = np.floor(transition_angle/180*nbins)
    coll_bins = np.linspace(0,transition_angle,num=N_low_bins) 
    elong_bins = np.linspace(transition_angle,180,num=nbins-N_low_bins) 

    collapsed_angles = computeAngularBins(np.exp(ddd.logr),ddd.u,ddd.v,collapsed=True)
    elongated_angles = computeAngularBins(np.exp(ddd.logr),ddd.u,ddd.v,collapsed=False)
    isRightSize = (np.exp(ddd.logr)*ddd.u > scale*ratio-scale*ratio*tolerance) & (np.exp(ddd.logr)*ddd.u < scale*ratio+scale*ratio*tolerance)
    isCollapsed = (((ddd.u*np.abs(ddd.v))*np.exp(ddd.logr)+np.exp(ddd.logr) > scale-scale*tolerance) & ((ddd.u*np.abs(ddd.v))*np.exp(ddd.logr)+np.exp(ddd.logr) < scale+scale*tolerance))
    isElongated = ((np.exp(ddd.logr) > scale-scale*tolerance) & (np.exp(ddd.logr) < scale+scale*tolerance))
    out1,bins1,_ = binned_statistic(elongated_angles[np.where(isRightSize & isElongated)],var[np.where(isRightSize & isElongated)],bins=elong_bins,statistic=stat)
    out2,bins2,_ = binned_statistic(collapsed_angles[np.where(isRightSize & isCollapsed)],var[np.where(isRightSize & isCollapsed)],bins=coll_bins,statistic=stat)
    full_var = np.concatenate((out2,out1))
    bins1 += (bins1[1]-bins1[0])/2 #make edges centers
    bins2 += (bins2[1]-bins2[0])/2
    full_bins = np.concatenate((bins2[:-1],bins1[:-1]))
    return full_var, full_bins

def compute_x_vs_side_length(ddd,var,stat='mean',nbins=15,tolerance=.1,**kwargs):
    isEquilateral = (ddd.u > 1-tolerance) & (np.abs(ddd.v) < tolerance)
    res, b, _ = binned_statistic(ddd.logr[isEquilateral],var[isEquilateral],bins=nbins,statistic=stat)
    b += (b[1]-b[0])/2
    b = b[:-1]
    return res, b

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

class NNNPlotter (object):

    def __init__(self,zvar,min_z,delta_z,metric):
        self.zvar = zvar
        self.min_z = min_z
        self.delta_z = delta_z
        self.max_z = self.min_z + self.delta_z
        self.metric = metric
        self.runname = self.zvar+str(self.min_z)+'_deltaz'+str(self.delta_z)+'_'+self.metric
        with open(outdir+self.runname+'.yaml') as f:
            self.config = yaml.load(f.read())
        self.data = np.load(self.config['data_path'])
        self.randoms = np.load(self.config['randoms_path'])
        assert self.runname == self.config['runname']
    
    def load_data_for_run(self):
        if self.zvar == 'DISTANCE':
            self.data = self.data[((self.data[self.zvar] > (buzzard_cosmo.h)*buzzard_cosmo.comoving_distance(self.min_z).value) & (self.data[self.zvar] < (buzzard_cosmo.h)*buzzard_cosmo.comoving_distance(self.max_z).value))]
            self.randoms = self.randoms[((self.randoms[self.zvar] > (buzzard_cosmo.h)*buzzard_cosmo.comoving_distance(self.min_z).value) & (self.randoms[self.zvar] < (buzzard_cosmo.h)*buzzard_cosmo.comoving_distance(self.max_z).value))]
        else:
            self.data = self.data[((self.data[self.zvar] > self.min_z) & (self.data[self.zvar] < self.max_z))]
            self.randoms = self.randoms[((self.randoms[self.zvar] > self.min_z) & (self.randoms[self.zvar] < self.max_z))]

        self.ddd = np.load(outdir+self.runname+'_'+'ddd.npy')
        self.ddr = np.load(outdir+self.runname+'_'+'ddr.npy')
        self.drd = np.load(outdir+self.runname+'_'+'drd.npy')
        self.rdd = np.load(outdir+self.runname+'_'+'rdd.npy')
        self.rrd = np.load(outdir+self.runname+'_'+'rrd.npy')
        self.drr = np.load(outdir+self.runname+'_'+'drr.npy')
        self.rdr = np.load(outdir+self.runname+'_'+'rdr.npy')
        self.rrr = np.load(outdir+self.runname+'_'+'rrr.npy')

    def analyze_single_run(self,mode,**kwargs):

        template = treecorr.NNNCorrelation(config=self.config)
        
        if mode == 'angle':
            get_binned_stat = computeXvsAngle
        if mode == 'equi':
            get_binned_stat = compute_x_vs_side_length
        
        binned = {}
        binned['ddd'], bins = get_binned_stat(template,self.ddd,stat='sum',**kwargs)
        binned['ddr'], bins = get_binned_stat(template,self.ddr,stat='sum',**kwargs)
        binned['drd'], bins = get_binned_stat(template,self.drd,stat='sum',**kwargs)
        binned['rdd'], bins = get_binned_stat(template,self.rdd,stat='sum',**kwargs)
        binned['rrd'], bins = get_binned_stat(template,self.rrd,stat='sum',**kwargs)
        binned['drr'], bins = get_binned_stat(template,self.drr,stat='sum',**kwargs)
        binned['rdr'], bins = get_binned_stat(template,self.rdr,stat='sum',**kwargs)
        binned['rrr'], bins = get_binned_stat(template,self.rrr,stat='sum',**kwargs)

        binned['d1'], bins = get_binned_stat(template,template.u*np.abs(template.v)*np.exp(template.logr)+np.exp(template.logr),**kwargs)
        binned['d2'], bins = get_binned_stat(template,np.exp(template.logr),**kwargs)
        binned['d3'], bins = get_binned_stat(template,template.u*np.exp(template.logr),**kwargs)

        datatot = len(self.data)
        randtot = len(self.randoms)
        dddtot = float(datatot)**3/6
        drrtot = float(datatot)*float(randtot)**2/6
        rdrtot = float(datatot)*float(randtot)**2/6
        rrdtot = float(datatot)*float(randtot)**2/6
        ddrtot = float(datatot)**2*float(randtot)/6
        drdtot = float(datatot)**2*float(randtot)/6
        rddtot = float(datatot)**2*float(randtot)/6
        rrrtot = float(randtot)**3/6

        binned['zeta'] = (binned['ddd']+dddtot*(-binned['ddr']/ddrtot-binned['drd']/drdtot-binned['rdd']/rddtot+binned['rrd']/rrdtot+binned['rdr']/rdrtot+binned['drr']/drrtot-binned['rrr']/rrrtot))/(binned['rrr']*dddtot/rrrtot)  
        binned['denom'] = self.get_two_point_expectation(binned['d1'],binned['d2'],binned['d3'])
        binned['q'] = binned['zeta']/binned['denom']
        return bins, binned

    def get_two_point_expectation(self,d1bins,d2bins,d3bins):
        if self.metric == 'Euclidean':
            cat = treecorr.Catalog(ra=self.data['RA'], dec=self.data['DEC'], ra_units='degrees', dec_units='degrees')
            random_cat = treecorr.Catalog(ra=self.randoms['RA'], dec=self.randoms['DEC'], ra_units='degrees', dec_units='degrees')
            dd = treecorr.NNCorrelation(min_sep=1,max_sep=30,nbins=30,bin_slop=0.1,sep_units='arcmin',metric=self.metric)
            dr = treecorr.NNCorrelation(min_sep=1,max_sep=30,nbins=30,bin_slop=0.1,sep_units='arcmin',metric=self.metric)
            rr = treecorr.NNCorrelation(min_sep=1,max_sep=30,nbins=30,bin_slop=0.1,sep_units='arcmin',metric=self.metric)
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

    def plot_run(self):
        
        results = pd.DataFrame()
        self.load_data_for_run()
        
        #make angular plots
        for scale in [10,15,20,25,30]:
            for ratio in [.5]:
                for tolerance in [.1,.2,.3]:
                    for nbins in [8,16,100]:
                        print (scale,ratio,tolerance,nbins)
                        sys.stdout.flush()
                        if ratio == 1:
                            mode = 'equi'
                        else:
                            mode = 'angle'
                        bins, binned = self.analyze_single_run(mode,scale=scale,ratio=ratio,tolerance=tolerance,nbins=nbins)
                        
                        this_res = pd.DataFrame.from_dict(binned)
                        this_res['bins'] = bins
                        this_res['scale'] = scale
                        this_res['ratio'] = ratio
                        this_res['tolerance'] = tolerance
                        this_res['nbins']= nbins
                        
                        results = results.append(this_res)
                        
                        if False:
                            for name,var in binned.iteritems():
                                fig = plt.figure()
                                plt.plot(bins,var)
                                if mode == 'angle':
                                    plt.xlabel('Angle (degrees)')
                                else:
                                    plt.xlabel('Scale (arcmin)')
                                plt.ylabel(name)
                                plt.title(str(self.min_z)+'<'+self.zvar+'<'+str(self.max_z)+' '+str(scale*ratio)+':'+str(scale)+' +/- '+str(100*tolerance)+'%')
                                fig.savefig(plotdir+name+'_'+mode+'_'+str(scale)+'_'+str(ratio)+'_'+str(tolerance)+'_'+str(nbins)+'.png')
                                
        results.to_csv(outdir+self.runname+'.csv')

def runall(min_z, max_z, delta_z, zvar, metric, do3D):
    for lower_z_lim in np.arange(min_z,max_z,delta_z):
        print ("bsub", "-W", "08:00", "python", "-c" ,"import autoplot; plotter = autoplot.NNNPlotter('"+zvar+"',"+str(lower_z_lim)+","+str(delta_z)+",'Euclidean'); plotter.plot_run()")
        subprocess.call(["bsub", "-W", "08:00", "python", "-c" ,"import autoplot; plotter = autoplot.NNNPlotter('"+zvar+"',"+str(lower_z_lim)+","+str(delta_z)+",'Euclidean'); plotter.plot_run()"])