#!/home/mbaumer/anaconda/bin/python2.7
from __future__ import division
import sys
sys.path.append('/home/mbaumer/anaconda2/bin/')
from astropy.io import fits
from astropy.cosmology import Planck15 as cosmo
import treecorr
import numpy as np
import time
import yaml
from os.path import expandvars

#needs to stay same as in submit.py
outdir = expandvars('$DES_DATA')+'/new_3pt_runs/'

class NNNProcessor (object):

    def __init__(self,config_fname):
        config_path = outdir+config_fname
        try: 
            with open(config_path) as f:
                configdict = json.loads(f.read())
        except IOError:
            print 'config file '+config_fname+' not found in '+outdir
            raise

        self.config = configdict
        self.data = np.load(self.config['data_path'])
        self.randoms = np.load(self.config['randoms_path'])

    def make_catalog(self):

        #redshift cut
        if zvar == 'DISTANCE':
            cat = cat[((cat[self.config['zvar']] > cosmo.comoving_distance(self.config['min_z']).value) & (cat[self.config['zvar']] < cosmo.comoving_distance(self.config['max_z']).value))]
        else:
            cat = cat[((cat[self.config['zvar']] > self.config['min_z']) & (cat[self.config['zvar']] < self.config['max_z']))]

        print 'Catalog length after redshift cut:', len(cat)
        sys.stdout.flush()

        if self.config['do3D']: 
            if self.config['zvar'] == 'DISTANCE':
                data_cat = treecorr.Catalog(ra=data['RA'], dec=data['DEC'], 
                ra_units='degrees', dec_units='degrees',
                r=data[self.config['zvar']]/cosmo.h)
            else:
                data_cat = treecorr.Catalog(ra=data['RA'], dec=data['DEC'], 
                ra_units='degrees', dec_units='degrees',
                r=cosmo.comoving_distance(data[self.config['zvar']])/cosmo.h)
        else: 
            data_cat = treecorr.Catalog(ra=data['RA'], dec=data['DEC'], 
                ra_units='degrees', dec_units='degrees')
            
        return data_cat

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

        fname = outdir+self.config['runname']+'_'+set1+set2+set3+'.npy'
        np.save(fname,nnn.ntri)

def run_3pt_ana(config_fname, set1, set2, set3):
	print 'constructing analysis'
    handler = NNNProcessor(config_fname)
	print 'handler constructed'
    handler.run(set1,set2,set3)
