#!/home/mbaumer/anaconda/bin/python2.7
from __future__ import division
import sys
sys.path.append('/home/mbaumer/anaconda2/bin/')
from astropy.io import fits
#from astropy.cosmology import Planck15 as cosmo
import treecorr
import numpy as np
import time
import yaml
from os.path import expandvars

from astropy.cosmology import FlatLambdaCDM
buzzard_cosmo = FlatLambdaCDM(68.81,.295)

class NNNProcessor (object):

    def __init__(self,config_fname):
        try: 
            with open(config_fname) as f:
                configdict = yaml.load(f.read())
        except IOError:
            print 'config file '+config_fname+' not found.'
            raise

        self.config = configdict
        self.data = np.load(self.config['data_path'])
        self.randoms = np.load(self.config['randoms_path'])

    def make_catalog(self,cat):

        #redshift cut
        if self.config['zvar'] == 'DISTANCE':
            cat = cat[((cat[self.config['zvar']] > (buzzard_cosmo.h)*buzzard_cosmo.comoving_distance(self.config['min_z']).value) & (cat[self.config['zvar']] < (buzzard_cosmo.h)*buzzard_cosmo.comoving_distance(self.config['max_z']).value))]
        else:
            cat = cat[((cat[self.config['zvar']] > self.config['min_z']) & (cat[self.config['zvar']] < self.config['max_z']))]

        print 'Catalog length after redshift cut:', len(cat)
        sys.stdout.flush()

        if self.config['param_3D'] != 0: 
            if self.config['zvar'] == 'DISTANCE':
                out_cat = treecorr.Catalog(ra=cat['RA'], dec=cat['DEC'], 
                ra_units='degrees', dec_units='degrees',
                r=cat[self.config['zvar']]/cosmo.h)
            else:
                out_cat = treecorr.Catalog(ra=cat['RA'], dec=cat['DEC'], 
                ra_units='degrees', dec_units='degrees',
                r=cosmo.comoving_distance(cat[self.config['zvar']])/cosmo.h)
        else: 
            out_cat = treecorr.Catalog(ra=cat['RA'], dec=cat['DEC'], 
                ra_units='degrees', dec_units='degrees')
            
        return out_cat

    def run(self,set1,set2,set3):

        cat = self.make_catalog(self.data)
        random_cat = self.make_catalog(self.randoms)

        setdict = {'d':cat,'r':random_cat}
        nnn = treecorr.NNNCorrelation(config=self.config)

        print 'starting processing!'
        sys.stdout.flush()

        toc = time.time()
        nnn.process(setdict[set1],setdict[set2],setdict[set3])
        tic = time.time()

        print 'that took', tic-toc
        sys.stdout.flush()

        fname = self.config['outdir']+self.config['runname']+'_'+set1+set2+set3+'.npy'
        np.save(fname,nnn.ntri)

def run_3pt_ana(config_fname, set1, set2, set3):
    print 'constructing analysis'
    handler = NNNProcessor(config_fname)
    print 'handler constructed'
    handler.run(set1,set2,set3)
