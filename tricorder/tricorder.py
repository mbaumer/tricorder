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
metric = 'Euclidean'
data_z_var = 'Z_SPEC'
random_z_var = 'Z'
do3D = False

#if 'sh-' in platform.node():
#    datapath = '/scratch/PI/kipac/mbaumer/des/data/redmagic_'
#    outdir = '/scratch/PI/kipac/mbaumer/des/3pt_results/'

class NNNProcessor (object):

    def __init__(self,runname, random_set_id):
        
        configdict = {}
        #treecorr ignores irrelevant keys

        self.runname = runname
        self.random_set_id = int(random_set_id)

        self.data_z_var = data_z_var
        self.random_z_var = random_z_var
        self.do3D = do3D

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
        
        configdict['min_u'] = 0
        configdict['max_u'] = 1
        configdict['nubins'] = 100
        
        configdict['min_v'] = -1
        configdict['max_v'] = 1
        configdict['nvbins'] = 200

        configdict['bin_slop'] = 1
        if not self.do3D:
            configdict['sep_units'] = 'arcmin'
        else: 
            try:
                del configdict['sep_units']
            except KeyError:
                pass

        self.config = configdict

        #write it out so we remember what we did
        config_fname = outdir+self.runname+'.config'
        #if (!os.path.exists(config_fname)): #not atomic; hard code for now
        if (self.random_set_id == 0): #just write out for first one
            f = open(config_fname,'w')
            f.write(json.dumps(self.config))
            f.close()

    def applyCuts(self,cat,zvar):
        try:
            cat = cat[((cat[zvar] > self.config['min_z']) & (cat[zvar] < self.config['max_z']))]
        except KeyError:
            print 'specified zvar: ', zvar, 'not found in ', cat.columns
            raise
        cat = cat[((cat['RA'] > footprint['min_ra']) & (cat['RA'] < footprint['max_ra']))]
        cat = cat[((cat['DEC'] > footprint['min_dec']) & (cat['DEC'] < footprint['max_dec']))]
        #later, color, magnitude cuts...
        return cat

    def prepare_data_cat(self):
        data = self.applyCuts(fits.getdata(self.config['datapath']),self.data_z_var)

        if self.do3D: 
            data_cat = treecorr.Catalog(ra=data['RA'], dec=data['DEC'], 
                ra_units='degrees', dec_units='degrees',
                r=cosmo.comoving_distance(data[self.data_z_var])/cosmo.h)
        else: 
            data_cat = treecorr.Catalog(ra=data['RA'], dec=data['DEC'], 
                ra_units='degrees', dec_units='degrees')
            
        return data_cat

    def prepare_random_cat(self):
        randoms = self.applyCuts(fits.getdata(self.config['randompath']),self.random_z_var)
        if self.do3D:
            random_cat = treecorr.Catalog(ra=randoms['RA'], dec=randoms['DEC'], 
                ra_units='degrees', dec_units='degrees',
                r=cosmo.comoving_distance(randoms[self.random_z_var])/cosmo.h)
        else:
            random_cat = treecorr.Catalog(ra=randoms['RA'], dec=randoms['DEC'], 
                ra_units='degrees', dec_units='degrees')
        return random_cat

    def prepare_joint_cat(self):
        data = self.applyCuts(fits.getdata(self.config['datapath']),self.data_z_var)
        randoms = self.applyCuts(fits.getdata(self.config['randompath']),self.random_z_var)

        joint_ra_table = np.hstack([data['RA'],randoms['RA']])
        joint_dec_table = np.hstack([data['DEC'],randoms['DEC']])  
        joint_z_table = np.hstack([data[self.data_z_var],randoms[self.random_z_var]])

        print 'joint cat shape: ', joint_ra_table.shape
        print 'randoms len: ', len(randoms['RA'])
        
        wt_factor = float(len(data['RA']))/float(len(randoms['RA']))
        weights = np.hstack([np.ones_like(data['RA']),-(wt_factor)*np.ones_like(randoms['RA'])])
        print 'sum of weights (should be close to zero): ', np.sum(weights)

        if self.do3D:
            joint_cat = treecorr.Catalog(ra=joint_ra_table, dec=joint_dec_table, 
                ra_units='degrees', dec_units='degrees',  
                r=cosmo.comoving_distance(joint_z_table)/cosmo.h, w=weights)
        else:
            joint_cat = treecorr.Catalog(ra=joint_ra_table, dec=joint_dec_table, 
                ra_units='degrees', dec_units='degrees', w=weights)

        return joint_cat

    def run(self,set1,set2,set3):

        cat = self.prepare_data_cat()
        random_cat = self.prepare_random_cat()
        joint_cat = self.prepare_joint_cat()

        setdict = {'d':cat,'r':random_cat,'n':joint_cat}
        nnn = treecorr.NNNCorrelation(config=self.config)

        print 'starting processing!'
        toc = time.time()

        if ((set1 == set2) and (set2 == set3)):
            #this is faster even though it shouldn't in theory...
            nnn.process(setdict[set1],metric=metric)
        else:
            nnn.process(setdict[set1],setdict[set2],setdict[set3],metric=metric)

        tic = time.time()
        print 'that took', tic-toc

        fname = outdir+self.config['runname']+'_'+str(self.random_set_id)+'_'+set1+set2+set3+'.out'
        nnn.write(fname,file_type='FITS')

def run_3pt_ana(runname, random_set_id, set1, set2, set3):
        handler = NNNProcessor(runname,random_set_id)
        handler.run(set1,set2,set3)
