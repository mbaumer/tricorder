#!/home/mbaumer/anaconda2/bin/python2.7
import sys, os.path
sys.path.append('/home/mbaumer/anaconda2/bin/')
from astropy.io import fits
import treecorr
import sys
import numpy as np
import time
import json
from sklearn.cluster import KMeans
from sklearn.neighbors import KNeighborsRegressor

outdir = '/scratch/PI/kipac/mbaumer/des/3pt_results/'

class NNNProcessor (object):

    def __init__(self,runname):
        config_fname = outdir+runname+'.config'
        
        configdict = {}
        #treecorr ignores irrelevant keys
        configdict['runname'] = runname
        configdict['datapath'] = '/scratch/PI/kipac/mbaumer/des/data/redmagic_sv_data.fits'
        configdict['randompath'] = '/scratch/PI/kipac/mbaumer/des/data/redmagic_sv_5x_randoms.fits'
	#configdict['datapath'] = '/scratch/PI/kipac/mbaumer/des/data/redmagic_Y1_sims_data.fits'
        #configdict['randompath'] = '/scratch/PI/kipac/mbaumer/des/data/redmagic_Y1_sims_5x_randoms.fits'

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

        #write it out so we remember what we did
        f = open(config_fname,'w')
        f.write(json.dumps(configdict))
        f.close()

        self.config = configdict

    def run(self,set1,set2,set3):
        cat, random_cat = self.prepareCatalog()
        setdict = {'d':cat,'r':random_cat}
        nnn = treecorr.NNNCorrelation(config=self.config)

        print 'starting processing!'
        toc = time.time()

        if ((set1 == set2) and (set2 == set3)):
            #this is faster even though it shouldn't in theory...
            nnn.process(setdict[set1])
        else:
            nnn.process(setdict[set1],setdict[set2],setdict[set3])

        tic = time.time()
        print 'that took', tic-toc

        fname = outdir+self.config['runname']+set1+set2+set3+'.out'
        nnn.write(fname,nnn)

if __name__ == '__main__':
    handler = NNNProcessor(sys.argv[4])
    handler.run(sys.argv[1], sys.argv[2], sys.argv[3])
