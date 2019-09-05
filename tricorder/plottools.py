from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import paths
from astropy.io import fits

import yaml
import treecorr
import emcee
import chainconsumer 
from sklearn.decomposition import PCA

from glob import glob

def load_config(config_fname):
    config_path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/'+config_fname+'.config'
    with open(config_path) as f:
        return yaml.load(f.read())
    
def get_q(config_name,fname):
    config = load_config(config_name)
    corr = treecorr.NNNCorrelation(config=config['3PCF'])
    dd = treecorr.NNCorrelation(config=config['2PCF'])
    zeta = np.load(fname)
    xi = np.load(fname[:-8]+'xi.npy')
    
    from scipy.interpolate import UnivariateSpline
    yfit = UnivariateSpline(dd.logr[xi != 0],np.log(xi)[xi != 0],k=4)
    xi1 = np.exp(yfit(corr.logr*(1+corr.u*np.abs(corr.v))))
    xi2 = np.exp(yfit(corr.logr))
    xi3 = np.exp(yfit(corr.logr*corr.u))
    denom = (xi1*xi2+xi2*xi3+xi3*xi1)
    return zeta/denom

def plot_qarr(qvec,**kwargs):
    q_mean = np.mean(qvec,axis=0)
    q_std = np.std(qvec,axis=0)
    v = np.linspace(-1,1,qvec.shape[2])[:,np.newaxis]*np.ones_like(q_mean.T)
    for row in np.arange(v.shape[1]):
        plt.errorbar(v[:,row]+.03*row,q_mean.T[:,row],yerr=q_std.T[:,row],**kwargs)
    plt.xlabel('v')
    plt.ylabel('q')

def infer_bias(q_dm_infer,q_gal_infer,icov,use_covmat=True):
    def lnprior(theta):
        b1, b2 = theta
        if 0 < b1 < 10 and -10 < b2 < 10:
            return 0.0
        return -np.inf

    def lnprob_noerror(bias):
        ln_prior = lnprior(bias)
        b1 = bias[0]
        b2 = bias[1]
        return -.5*np.sum((((q_gal_infer-q_dm_infer/b1-b2/(b1**2))**2)/icov)) + ln_prior

    def lnprob(bias):
        ln_prior = lnprior(bias)
        b1 = bias[0]
        b2 = bias[1]
        resid = (q_gal_infer-q_dm_infer/b1-b2/(b1**2))
        return -.5*np.matmul(resid,np.matmul(icov,resid)) + ln_prior
        #return -.5*np.sum(np.dot(np.dot(resid,icov),resid)) + ln_prior

    ndim, nwalkers = 2, 50
    p0 = [np.array([10*np.random.rand(),20*np.random.rand()-10]) for i in range(nwalkers)]

    if use_covmat:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
    else:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_noerror)
    sampler.run_mcmc(p0, 10000)
    
    cc = chainconsumer.ChainConsumer()
    cc.add_chain(sampler.flatchain,parameters=['b1','b2'],name='3pt ZSPEC')
    summary = cc.analysis.get_summary()
    b1=summary['b1'][1]
    b2=summary['b2'][1]
    resid = (q_gal_infer-q_dm_infer/b1-b2/(b1**2))
    print 'chisq is: ', str(.5*np.sum(np.dot(np.dot(resid,icov),resid)))
    print 'alt: ', .5*np.matmul(np.matmul(resid,icov),resid)
    
    #cc = chainconsumer.ChainConsumer()
    #cc.add_chain(sampler.flatchain,parameters=['b1','b2'])
    
    #cc.plotter.plot(figsize='column',extents=[(0,3),(-1,1)]);
    return sampler

def plot_data_vectors(qdm,qrm,v=None,b1=1,b2=0,b1MAP=1,b2MAP=0,rm_color='g',dm_color='b'):
    qdm_mean = np.mean(qdm,axis=0)
    qdm_std = np.std(qdm,axis=0)
    
    qrm_mean = np.mean(qrm,axis=0)
    qrm_std = np.std(qrm,axis=0)
    
    if v is None:
        v = np.arange(qrm.shape[1])
        plt.xlabel('Data Vector Index')
    else:    
        plt.xlabel('v')
        
    plt.errorbar(v,qdm_mean/b1+b2/(b1**2),yerr=qdm_std/np.sqrt(qdm.shape[0]),label='DM',color=dm_color)
    plt.errorbar(v,qrm_mean,yerr=qrm_std,label='Galaxies',color=rm_color)
    if b1 != 1:
        plt.title('Best-fit dark matter')
    else:
        plt.title('Data Vectors')
    plt.ylabel('q')
    plt.legend()
    
def make_inference(red_qdm,red_qrm,covmat_src=None,max_pca_comps=None,is11k=False):
    
    if max_pca_comps is None:
        rmcov = np.cov(red_qrm.T)
        dmcov = np.cov(red_qdm.T)
    else: 
        pca = PCA(n_components = max_pca_comps)
        pca.fit(red_qrm)
        print pca.explained_variance_ratio_
        rmcov = pca.get_covariance()
    
    if covmat_src is None:
        #icov = np.linalg.inv(dmcov/red_qdm.shape[0]+rmcov)
        if is11k:
            #icov = np.linalg.inv(rmcov/len(red_qrm)+dmcov/len(red_qdm))
            # corrected for https://arxiv.org/pdf/astro-ph/0608064.pdf 
            icov = np.linalg.inv(rmcov/len(red_qrm)+dmcov/len(red_qdm))*(len(red_qrm)-red_qrm.shape[1]-1)/(len(red_qrm)-1)
            #Let's see how it goes taking out the DM errors:
            #icov = np.linalg.inv(rmcov/len(red_qrm))*(len(red_qrm)-red_qrm.shape[1]-1)/(len(red_qrm)-1)
        else:
            #icov = np.linalg.inv(rmcov+dmcov/len(red_qdm))
            # corrected for https://arxiv.org/pdf/astro-ph/0608064.pdf 
            icov = np.linalg.inv(rmcov+dmcov/len(red_qdm))*(len(red_qrm)-red_qrm.shape[1]-1)/(len(red_qrm)-1)
            #Let's see how it goes taking out the DM errors:
            #icov = np.linalg.inv(rmcov)*(len(red_qrm)-red_qrm.shape[1]-1)/(len(red_qrm)-1)
        qrm_mean = np.mean(red_qrm,axis=0)
        qdm_mean = np.mean(red_qdm,axis=0)
    else:
        qrm_mean = red_qrm
        qdm_mean = red_qdm
        icov = np.linalg.inv(covmat_src)
    samples = infer_bias(qdm_mean,qrm_mean,icov,use_covmat=True)
    
    #cc = chainconsumer.ChainConsumer()
    #cc.add_chain(samples.flatchain,parameters=['b1','b2'],name='3pt ZSPEC')
    #cc.plotter.plot(figsize='column',extents=[(0,3),(-1,1)]);
    
    return samples

from glob import glob
import os.path
import matplotlib.pyplot as plt
import treecorr
import yaml
import paths
import numpy as np
import pandas as pd
def load_config(config_fname):
    config_path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/'+config_fname+'.config'
    with open(config_path) as f:
        return yaml.load(f.read())
    
def compress_dv(qdm1):
    return ((qdm1[:,5:] + np.flip(qdm1[:,:5],axis=1))/2.0).reshape(qdm1.shape[0],5)[:,-5:]




def load_res(path,dset_name,dset_id,config_fname,zvar,min_z,max_z,rsamp_str,sigma=0,norm=True):
    first_done = False
    for jk_id in range(15):
        try:
            this = load_files(path,dset_name,dset_id,config_fname,zvar,min_z,max_z,rsamp_str,sigma,return_all_norm=norm,jk_id=jk_id,config=config_fname)
        except IOError:
            continue
        this['JK'] = jk_id
        if not first_done:
            res1 = this
            first_done = True
        else:
            res1 = pd.concat([res1,this],ignore_index=True)
    return res1

def load_res_indep(path,dset_name,config_fname,zvar,min_z,max_z,rsamp_str,norm=True,sigma=0,use_alt_randoms=True):
    first_done = False
    for dset_id in range(24):
        try:
            this = load_files(path,dset_name,dset_id,config_fname,zvar,min_z,max_z,rsamp_str,sigma,return_all_norm=norm,config=config_fname,use_alt_randoms=use_alt_randoms)
        except IOError:
            continue
        this['DSET'] = dset_id
        if not first_done:
            res1 = this
            first_done = True
        else:
            res1 = pd.concat([res1,this],ignore_index=True)
    return res1

def load_res_indep2(path,dset_name,config_fname,zvar,min_z,max_z,rsamp_str,norm=False):
    first_done = False
    for dset_id in range(24):
        if norm:
            try:
                this = load_files(path+config_fname+'_'+dset_name+'dset0_jk'+str(dset_id)+'_sigma0_'+zvar+'_'+str(min_z)+'_'+str(max_z)+'_rsamp'+rsamp_str,return_all_norm=True)
            except IOError:
                continue
        else:
            try:
                this = load_files(path+config_fname+'_'+dset_name+'dset0_jk'+str(dset_id)+'_sigma0_'+zvar+'_'+str(min_z)+'_'+str(max_z)+'_rsamp'+rsamp_str,return_all=True)
            except IOError:
                continue
        this['DSET'] = dset_id
        if not first_done:
            res1 = this
            first_done = True
        else:
            res1 = pd.concat([res1,this],ignore_index=True)
    return res1

def load_res_xi(path,dset_name,dset_id,config_fname,min_z,max_z,rsamp_str):
    xilist = []
    for jk_id in range(15):
        try:
            this = np.load(path+config_fname+'_'+dset_name+'dset'+str(dset_id)+'_jk'+str(jk_id)+'_sigma0_ZSPEC_'+str(min_z)+'_'+str(max_z)+'_rsamp'+rsamp_str+'.xi.npy')
        except IOError:
            continue
        xilist.append(this.flatten())
    return xilist

def load_res_xi_indep(path,dset_name,config_fname,min_z,max_z,rsamp_str,zvar='ZSPEC',sigma=0):
    xilist = []
    for dset_id in range(24):
        try:
            this = np.load(path+config_fname+'_'+dset_name+'dset'+str(dset_id)+'_jk-1_sigma'+str(sigma)+'_'+zvar+'_'+str(min_z)+'_'+str(max_z)+'_rsamp'+rsamp_str+'.xi.npy')
        except IOError:
            continue
        xilist.append(this.flatten())
    return xilist

def load_files(path,dset_name,dset_id,config_fname,zvar,min_z,max_z,rsamp_str,sigma=0,get_q=True,jk_id=-1,config='newpaper13.1',return_all=False,return_all_norm=False,use_alt_randoms=True):
    runname = config_fname+'_'+dset_name+'dset'+str(dset_id)+'_jk'+str(jk_id)+'_sigma'+str(sigma)+'_'+zvar+'_'+str(min_z)+'_'+str(max_z)+'_rsamp'+rsamp_str
    ddd = treecorr.NNNCorrelation(config=load_config(config)['3PCF'])
    ddd.read(path+runname+'.ddd')
    ddr = treecorr.NNNCorrelation(config=load_config(config)['3PCF'])
    ddr.read(path+runname+'.ddr')
    drd = treecorr.NNNCorrelation(config=load_config(config)['3PCF'])
    drd.read(path+runname+'.drd')
    rdd = treecorr.NNNCorrelation(config=load_config(config)['3PCF'])
    rdd.read(path+runname+'.rdd')
    rrd = treecorr.NNNCorrelation(config=load_config(config)['3PCF'])
    rrd.read(path+runname+'.rrd')
    drr = treecorr.NNNCorrelation(config=load_config(config)['3PCF'])
    drr.read(path+runname+'.drr')
    rdr = treecorr.NNNCorrelation(config=load_config(config)['3PCF'])
    rdr.read(path+runname+'.rdr')
    rrr = treecorr.NNNCorrelation(config=load_config(config)['3PCF'])
    try:
        rrr.read(path+runname+'.rrr')
    except IOError:
        print 'missing randoms'
        if use_alt_randoms:
            print 'using new randoms'
            newrunnamelist = glob(path+config_fname+'_'+dset_name+'dset*_jk'+str(jk_id)+'_sigma'+str(sigma)+'_'+zvar+'_'+str(min_z)+'_'+str(max_z)+'_rsamp'+rsamp_str+'.rrr')
            if newrunnamelist != []:
                newrunname = newrunnamelist[-1]
                print newrunname
                rrr.read(newrunname)
                print rrr.ntri
            else:
                print 'couldnt find alternate rrr file'
                raise IOError
        else:
            raise IOError
        
    corr = treecorr.NNNCorrelation(config = load_config(config)['3PCF'])
    dd = treecorr.NNCorrelation(config = load_config(config)['2PCF'])
    zeta = (ddd.ntri/ddd.tot - ddr.ntri/ddr.tot - drd.ntri/drd.tot - rdd.ntri/rdd.tot + rrd.ntri/rrd.tot + drr.ntri/drr.tot + rdr.ntri/rdr.tot - rrr.ntri/rrr.tot)/(rrr.ntri/rrr.tot)
    dddrrr = ddd.ntri/ddd.tot/(rrr.ntri/rrr.tot)
    oldest = (ddd.ntri/ddd.tot - ddr.ntri/ddr.tot - drd.ntri/drd.tot - rdd.ntri/rdd.tot)/(rrr.ntri/rrr.tot) + 2

    xi = np.load(path+runname+'.xi.npy')
    
    from scipy.interpolate import UnivariateSpline
    
    yfit = UnivariateSpline(dd.logr[xi != 0],np.log(xi)[xi != 0],k=4)
    xi1 = np.exp(yfit(corr.logr*(1+corr.u*np.abs(corr.v))))
    xi2 = np.exp(yfit(corr.logr))
    xi3 = np.exp(yfit(corr.logr*corr.u))
    denom = (xi1*xi2+xi2*xi3+xi3*xi1)
    
    #insert elisabeth stuff here
    
    yfit = UnivariateSpline(dd.logr[xi != 0],np.log(xi)[xi != 0],k=4)
    xi1 = np.exp(yfit(np.log(np.exp(corr.logr)*(1+corr.u*np.abs(corr.v)))))
    xi2 = np.exp(yfit(np.log(np.exp(corr.logr))))
    xi3 = np.exp(yfit(np.log(np.exp(corr.logr)*corr.u)))
    correct_denom = (xi1*xi2+xi2*xi3+xi3*xi1)
    
    if return_all:
        return pd.DataFrame(np.array([ddd.ntri.flatten().T,ddr.ntri.flatten().T,drd.ntri.flatten().T,
                             rdd.ntri.flatten().T,drr.ntri.flatten().T,rdr.ntri.flatten().T, 
                             rrd.ntri.flatten().T,rrr.ntri.flatten().T,zeta.flatten().T,(zeta/denom).flatten().T, dddrrr.flatten().T, oldest.flatten().T]).T, columns=['DDD','DDR','DRD','RDD','DRR','RDR','RRD','RRR','ZETA','Q','DDD/RRR','OLDEST'])
    
    if return_all_norm:
        return pd.DataFrame(np.array([ddd.ntri.flatten().T/ddd.tot,ddr.ntri.flatten().T/ddr.tot,drd.ntri.flatten().T/drd.tot,
                             rdd.ntri.flatten().T/rdd.tot,drr.ntri.flatten().T/drr.tot,rdr.ntri.flatten().T/rdr.tot, 
                             rrd.ntri.flatten().T/rrd.tot,rrr.ntri.flatten().T/rrr.tot,zeta.flatten().T,denom.flatten().T, correct_denom.flatten().T, (zeta/correct_denom).flatten().T, dddrrr.flatten().T, oldest.flatten().T]).T, columns=['DDD','DDR','DRD','RDD','DRR','RDR','RRD','RRR','ZETA','DENOM','CORRECTDENOM','Q','DDD/RRR','OLDEST'])
    
    if get_q:
        return zeta/denom
    else:
        return zeta

def old_load_files(runname,rsamp,get_q=True,config='paper13.1',return_all=False):
    path = ''
    print len(glob(path+runname+'*'))
    ddd = treecorr.NNNCorrelation(config=load_config(config)['3PCF'])
    ddd.read(path+runname+'.ddd')
    ddr = treecorr.NNNCorrelation(config=load_config(config)['3PCF'])
    ddr.read(path+runname+'.ddr')
    drd = treecorr.NNNCorrelation(config=load_config(config)['3PCF'])
    drd.read(path+runname+'.drd')
    rdd = treecorr.NNNCorrelation(config=load_config(config)['3PCF'])
    rdd.read(path+runname+'.rdd')
    rrd = treecorr.NNNCorrelation(config=load_config(config)['3PCF'])
    rrd.read(path+runname+'.rrd')
    drr = treecorr.NNNCorrelation(config=load_config(config)['3PCF'])
    drr.read(path+runname+'.drr')
    rdr = treecorr.NNNCorrelation(config=load_config(config)['3PCF'])
    rdr.read(path+runname+'.rdr')
    rrr = treecorr.NNNCorrelation(config=load_config(config)['3PCF'])
    rrr.read(path+runname+'.rrr')

    corr = treecorr.NNNCorrelation(config = load_config(config)['3PCF'])
    dd = treecorr.NNCorrelation(config = load_config(config)['2PCF'])
    zeta = (ddd.ntri - ddr.ntri/rsamp - drd.ntri/rsamp - rdd.ntri/rsamp + rrd.ntri/rsamp**2 + drr.ntri/rsamp**2 + rdr.ntri/rsamp**2 - rrr.ntri/rsamp**3)/(rrr.ntri/rsamp**3)

    xi = np.load(path+runname+'.xi.npy')
    
    from scipy.interpolate import UnivariateSpline
    yfit = UnivariateSpline(dd.logr[xi != 0],np.log(xi)[xi != 0],k=4)
    xi1 = np.exp(yfit(corr.logr*(1+corr.u*np.abs(corr.v))))
    xi2 = np.exp(yfit(corr.logr))
    xi3 = np.exp(yfit(corr.logr*corr.u))
    denom = (xi1*xi2+xi2*xi3+xi3*xi1)
    
    if return_all:
        print zeta.flatten()
        return pd.DataFrame(np.array([ddd.ntri.flatten().T,ddr.ntri.flatten().T,drd.ntri.flatten().T,
                             rdd.ntri.flatten().T,drr.ntri.flatten().T,rdr.ntri.flatten().T, 
                             rrd.ntri.flatten().T,rrr.ntri.flatten().T,zeta.flatten().T,(zeta/denom).flatten().T]).T, columns=['DDD','DDR','DRD','RDD','DRR','RDR','RRD','RRR','ZETA','Q'])
    
    if get_q:
        return zeta/denom
    else:
        return zeta

def load_files2(runname,rsamp):
    ddd = treecorr.NNNCorrelation(config=load_config('paper13.1')['3PCF'])
    ddd.read(path+runname+'.ddd')
    rrr = treecorr.NNNCorrelation(config=load_config('paper13.1')['3PCF'])
    rrr.read(path+runname+'.rrr')
    return (ddd.ntri - rrr.ntri/rsamp**3)/(rrr.ntri/rsamp**3)

def load_files3(runname,rsamp):
    ddd = treecorr.NNNCorrelation(config=load_config('paper13.1')['3PCF'])
    ddd.read(path+runname+'.ddd')
    ddr = treecorr.NNNCorrelation(config=load_config('paper13.1')['3PCF'])
    ddr.read(path+runname+'.ddr')
    drd = treecorr.NNNCorrelation(config=load_config('paper13.1')['3PCF'])
    drd.read(path+runname+'.drd')
    rdd = treecorr.NNNCorrelation(config=load_config('paper13.1')['3PCF'])
    rdd.read(path+runname+'.rdd')
    rrr = treecorr.NNNCorrelation(config=load_config('paper13.1')['3PCF'])
    rrr.read(path+runname+'.rrr')
    return (ddd.ntri - ddr.ntri/(rsamp) - drd.ntri/(rsamp) - rdd.ntri/(rsamp))/(rrr.ntri/rsamp**3) + 2


def plot_data_vectors(qdm,qrm,v=None,b1=1,b2=0,b1MAP=1,b2MAP=0,rm_color='g',dm_color='b'):
    qdm_mean = np.mean(qdm,axis=0)
    qdm_std = np.std(qdm,axis=0)
    
    qrm_mean = np.mean(qrm,axis=0)
    qrm_std = np.std(qrm,axis=0)
    
    if v is None:
        v = np.arange(qrm.shape[1])
        plt.xlabel('Data Vector Index')
    else:    
        plt.xlabel('v')
        
    plt.errorbar(v,qdm_mean/b1+b2/(b1**2),yerr=qdm_std/np.sqrt(qdm.shape[0]),label='DM',color=dm_color)
    plt.errorbar(v,qrm_mean,yerr=qrm_std,label='Galaxies',color=rm_color)
    if b1 != 1:
        plt.title('Best-fit dark matter')
    else:
        plt.title('Data Vectors')
    plt.ylabel('q')
    plt.legend()

def plot_dv(res,var,indiv_runs=True,offset=0,compressed=False,**kwargs):
    if compressed:
        if indiv_runs:
            plt.plot(np.linspace(.1,.9,5),compress_dv(res[var].values.reshape(-1,10)).T,alpha=.3,color=kwargs['color'])
        plt.errorbar(np.linspace(.1,.9,5)+offset,np.mean(compress_dv(res[var].values.reshape(-1,10)),axis=0),
                     yerr=np.std(compress_dv(res[var].values.reshape(-1,10)),axis=0),**kwargs)
    else:
        if indiv_runs:
            plt.plot(np.linspace(-.9,.9,10),res[var].values.reshape(-1,10).T,alpha=.3,color=kwargs['color'])
        plt.errorbar(np.linspace(-.9,.9,10)+offset,np.mean(res[var].values.reshape(-1,10),axis=0),
                     yerr=np.std(res[var].values.reshape(-1,10),axis=0),**kwargs)
    
def get_max_like(red_qdm,red_qrm):
    b1 = np.linspace(0, 3, 300)
    b2 = np.linspace(-3, 3, 300)
    z = np.zeros((len(b2),len(b1)))
    icov = np.linalg.inv(np.cov(red_qrm.T)+np.cov(red_qdm.T)/len(red_qdm))
    for i,tb1 in enumerate(b1):
        for j,tb2 in enumerate(b2):
            resid = (np.mean(red_qrm,axis=0)-np.mean(red_qdm,axis=0)/tb1-tb2/(tb1**2))
            z[j,i] = .5*np.matmul(np.matmul(resid,icov),resid)
    print 'manual chisq: ', min(z[~np.isnan(z)])
    return b1, b2, z

def get_max_like_diag(red_qdm,red_qrm):
    b1 = np.arange(0, 3, 0.01)
    b2 = np.arange(-1, 1, 0.01)
    z2 = np.zeros((len(b2),len(b1)))
    icov = np.linalg.inv(np.diag(np.cov(red_qrm.T)+np.cov(red_qdm.T)/len(red_qdm))*np.eye(5))
    for i,tb1 in enumerate(b1):
        for j,tb2 in enumerate(b2):
            resid = (np.mean(red_qrm,axis=0)-np.mean(red_qdm,axis=0)/tb1-tb2/(tb1**2))
            z2[j,i] = .5*np.matmul(resid,np.matmul(resid,icov))
    print 'diag chisq: ', min(z2[~np.isnan(z2)])
    return b1, b2, z2

def get_lazeyras_schmidt_2015(b1):
    return .412-2.143*b1+.929*b1**2+.008*b1**3

def get_hoffman_2015(b1):
    return .51-2.21*b1+1*b1**2
