#!/home/mbaumer/anaconda2/bin/python2.7
import sys, os.path
sys.path.append('/home/mbaumer/anaconda2/bin/')
from astropy.io import fits
import treecorr
import sys
import numpy as np
import time
import json
import platform

if 'ki-ls' in platform.node():
    datapath = '/nfs/slac/g/ki/ki19/des/mbaumer/3pt_data/redmagic_'
    outdir = '/nfs/slac/g/ki/ki19/des/mbaumer/3pt_results/'

if 'sh-' in platform.node():
    datapath = '/scratch/PI/kipac/mbaumer/des/data/redmagic_'
    outdir = '/scratch/PI/kipac/mbaumer/des/3pt_results/'

runType = 'Y1_sims'

y1_main_footprint = {}
y1_main_footprint['min_ra'] = 0
y1_main_footprint['max_ra'] = 360
y1_main_footprint['min_dec'] = -70
y1_main_footprint['max_dec'] = -35

class NNNProcessor (object):

    def __init__(self,runname, random_set_id):
        
        configdict = {}
        #treecorr ignores irrelevant keys

        self.runname = runname
        self.random_set_id = random_set_id
        configdict['runType'] = runType

        if configdict['runType'] == 'Y1_sims':
            configdict['datapath'] = datapath+'data.fits'
            configdict['randompath'] = datapath+'randoms_'+str(self.random_set_id)+'.fits'
            self.footprint = y1_main_footprint

        else:
            raise IOError('invalid runType')
    
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

        configdict['bin_slop'] = 0.1
        configdict['sep_units'] = 'arcmin'

        self.config = configdict

        self.cat = None
        self.random_cat = None

        #write it out so we remember what we did
        config_fname = outdir+self.runname+'.config'
        #if (!os.path.exists(config_fname)): #not atomic; hard code for now
        if (self.random_set_id == 0): #just write out for first one
            f = open(config_fname,'w')
            f.write(json.dumps(self.config))
            f.close()

    def applyCuts(self,cat):
        if 'ZREDMAGIC' in cat.columns:
            cat = cat[((cat['ZREDMAGIC'] > self.config['min_z']) & (cat['ZREDMAGIC'] < self.config['max_z']))]
        elif 'Z' in cat.columns:
            cat = cat[((cat['Z'] > self.config['min_z']) & (cat['Z'] < self.config['max_z']))]
        else:
            print 'one or more input cats have no redshift data'
        cat = cat[((cat['RA'] > self.footprint['min_ra']) & (cat['RA'] < self.footprint['max_ra']))]
        cat = cat[((cat['DEC'] > self.footprint['min_dec']) & (cat['DEC'] < self.footprint['max_dec']))]
        #later, color, magnitude cuts...
        return cat

    def prepareCatalogs(self):
        data = self.applyCuts(fits.getdata(self.config['datapath']))
        randoms = self.applyCuts(fits.getdata(self.config['randompath']))

        joint_ra_table = np.hstack([data['RA'],randoms['RA']])
        joint_dec_table = np.hstack([data['DEC'],randoms['DEC']])  
        
        wt_factor = len(data['RA'])/len(randoms['RA'])
        weights = np.hstack([np.ones_like(data['RA']),-(wt_factor)*np.ones_like(randoms['RA'])])

        joint_cat = treecorr.Catalog(ra=joint_ra_table, dec=joint_dec_table, ra_units='degrees', dec_units='degrees', w=weights)
        random_cat = treecorr.Catalog(ra=randoms['RA'], dec=randoms['DEC'], ra_units='degrees', dec_units='degrees')

        return joint_cat, random_cat

    def run(self):
        joint_cat, random_cat = self.prepareCatalogs()

        nnn = treecorr.NNNCorrelation(config=self.config)
        rrr = treecorr.NNNCorrelation(config=self.config)

        print 'starting numerator!'
        toc = time.time()
        nnn.process(joint_cat)
        tic = time.time()
        print 'numerator took', tic-toc

        print 'starting denominator!'
        toc = time.time()
        rrr.process(random_cat)
        tic = time.time()
        print 'denominator took', tic-toc

        fname = outdir+self.config['runname']+str(self.random_set_id)+'.out'
        nnn.write(fname,nnn,rrr=rrr,file_type='FITS')

if __name__ == '__main__':
    handler = NNNProcessor(sys.argv[1],sys.argv[2])
    handler.run()
