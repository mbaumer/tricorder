#!/home/mbaumer/anaconda2/bin/python2.7
from __future__ import division
import sys, os.path
sys.path.append('/home/mbaumer/anaconda2/bin/')
from astropy.io import fits
from astropy.cosmology import Planck15 as cosmo
import treecorr
import sys
import numpy as np
import time
import json
import platform

y1_main = {}
y1_main['min_ra'] = 0
y1_main['max_ra'] = 360
y1_main['min_dec'] = -70
y1_main['max_dec'] = -35

##User settings!
footprint = y1_main
datapath = '/nfs/slac/g/ki/ki19/des/mbaumer/3pt_data/redmagic_'
outdir = '/nfs/slac/g/ki/ki19/des/mbaumer/3pt_runs/'
metric = 'Rperp'

#if 'sh-' in platform.node():
#    datapath = '/scratch/PI/kipac/mbaumer/des/data/redmagic_'
#    outdir = '/scratch/PI/kipac/mbaumer/des/3pt_results/'

class NNNProcessor (object):

    def __init__(self,runname, random_set_id):
        
        configdict = {}
        #treecorr ignores irrelevant keys

        self.runname = runname
        self.random_set_id = int(random_set_id)

        configdict['datapath'] = datapath+'data.fits'
        configdict['randompath'] = datapath+'randoms_'+str(self.random_set_id)+'.fits'

        for key,value in footprint.iteritems():
            configdict[key] = value
    
        configdict['metric'] = metric

        configdict['min_z'] = .5
        configdict['max_z'] = .7

        configdict['runname'] = runname
        configdict['min_sep'] = 1
        configdict['max_sep'] = 25
        configdict['nbins'] = 200
        
        configdict['min_u'] = .4
        configdict['max_u'] = 1
        configdict['nubins'] = 60
        
        configdict['min_v'] = 0
        configdict['max_v'] = 1
        configdict['nvbins'] = 100

        configdict['bin_slop'] = 1
        configdict['sep_units'] = 'arcmin'

        self.config = configdict

        #write it out so we remember what we did
        config_fname = outdir+self.runname+'.config'
        #if (!os.path.exists(config_fname)): #not atomic; hard code for now
        if (self.random_set_id == 0): #just write out for first one
            f = open(config_fname,'w')
            f.write(json.dumps(self.config))
            f.close()

    def applyCuts(self,cat):
        if 'ZREDMAGIC' in cat.dtype.names:
            cat = cat[((cat['ZREDMAGIC'] > self.config['min_z']) & (cat['ZREDMAGIC'] < self.config['max_z']))]
        elif 'Z' in cat.dtype.names:
            cat = cat[((cat['Z'] > self.config['min_z']) & (cat['Z'] < self.config['max_z']))]
        else:
            print 'one or more input cats have no redshift data'
        cat = cat[((cat['RA'] > footprint['min_ra']) & (cat['RA'] < footprint['max_ra']))]
        cat = cat[((cat['DEC'] > footprint['min_dec']) & (cat['DEC'] < footprint['max_dec']))]
        #later, color, magnitude cuts...
        return cat

    def prepareCatalogs(self):
        data = self.applyCuts(fits.getdata(self.config['datapath']))
        randoms = self.applyCuts(fits.getdata(self.config['randompath']))

        joint_ra_table = np.hstack([data['RA'],randoms['RA']])
        joint_dec_table = np.hstack([data['DEC'],randoms['DEC']])  
        joint_z_table = np.hstack([data['ZREDMAGIC'],randoms['Z']])

        print 'joint cat has ', joint_ra_table.shape
        print 'randoms have ', len(randoms['RA'])
        
        wt_factor = len(data['RA'])/len(randoms['RA'])
        weights = np.hstack([np.ones_like(data['RA']),-(wt_factor)*np.ones_like(randoms['RA'])])
        print 'sum of weights (should be close to zero): ', np.sum(weights)

        joint_cat = treecorr.Catalog(ra=joint_ra_table, dec=joint_dec_table, 
            ra_units='degrees', dec_units='degrees',  
            r=cosmo.comoving_distance(joint_z_table)/cosmo.h, w=weights)
        random_cat = treecorr.Catalog(ra=randoms['RA'], dec=randoms['DEC'], 
            ra_units='degrees', dec_units='degrees',
            r=cosmo.comoving_distance(randoms['Z'])/cosmo.h)

        return joint_cat, random_cat

    def run(self):
        joint_cat, random_cat = self.prepareCatalogs()

        nnn = treecorr.NNNCorrelation(config=self.config)

        print 'starting numerator!'
        toc = time.time()
        nnn.process(joint_cat,metric=metric)
        tic = time.time()
        print 'numerator took', tic-toc

        fname = outdir+self.config['runname']+str(self.random_set_id)+'.out'
        nnn.write(fname,file_type='FITS')

        if self.random_set_id == 0:
            print 'also doing denominator'
            rrr = treecorr.NNNCorrelation(config=self.config)
            toc = time.time()
            rrr.process(random_cat,metric=metric)
            tic = time.time()
            print 'denominator took', tic-toc
            fname = outdir+self.config['runname']+'RRR'+str(self.random_set_id)+'.out'
            rrr.write(fname,file_type='FITS')


if __name__ == '__main__':
    handler = NNNProcessor(sys.argv[1],sys.argv[2])
    handler.run()
