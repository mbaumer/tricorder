#!/home/mbaumer/anaconda2/bin/python2.7
import sys, os.path
sys.path.append('/home/mbaumer/anaconda2/bin/')
from astropy.io import fits
import treecorr
import sys
import numpy as np
import time
import json

class NNNHandler (object):

    def __init__(self,runname):
        self.runname = runname
        self.datapath = '/scratch/PI/kipac/mbaumer/des/data/redmagic_sv_data.fits'
        self.randompath = '/scratch/PI/kipac/mbaumer/des/data/redmagic_sv_randoms.fits'
        self.outdir = '/scratch/PI/kipac/mbaumer/des/3pt_results/'
        #write config file
        config_fname = outdir+runname+'.config'

        #if runname is already done, read in results
        if os.path.exists(config_fname):
            with open(config_fname) as f:
                configdict = json.loads(f.read())
        #otherwise, make a new one
        else:
            configdict = {}
            configdict['min_sep'] = 1
            configdict['max_sep'] = 50
            configdict['nbins'] = 11
            
            configdict['min_u'] = 0.9
            configdict['max_u'] = 1.0
            configdict['nubins'] = 1
            
            configdict['min_v'] = -0.1
            configdict['max_v'] = 0.1
            configdict['nvbins'] = 1

            configdict['bin_slop'] = 0.1
            configdict['sep_units'] = 'arcmin'

            #rest of this is for other parts of analysis; treecorr doesn't care
            configdict['min_z'] = .5
            configdict['max_z'] = .6

            #write it out so we remember what we did
            f = open(config_fname,'w')
            f.write(json.dumps(configdict))
            f.close()

        self.config = configdict

    def load(self):
        [ddd,ddr,drd,rdd,drr,rdr,rrd,rrr] = 8*[treecorr.NNNCorrelation(config=self.config)]
        for nnn in [ddd,ddr,drd,rdd,drr,rdr,rrd,rrr]:
            #do something

    def run(self,set1,set2,set3):
        data = fits.getdata(self.datapath)
        randoms = fits.getdata(self.randompath)

        cat = treecorr.Catalog(ra=data['RA'], dec=data['DEC'], ra_units='degrees', dec_units='degrees')
        random_cat = treecorr.Catalog(ra=randoms['RA'], dec=randoms['DEC'], ra_units='degrees', dec_units='degrees')
        setdict = {'d':cat,'r':random_cat}
        nnn = treecorr.NNNCorrelation(config=self.config)
        
        #TODO:
        #write out config dict used for processing runname

        print 'starting processing!'
        toc = time.time()
        if ((set1 == set2) and (set2 == set3)):
            nnn.process(setdict[set1])
            #this saves time even though it shouldn't in theory...
        else:
            nnn.process(setdict[set1],setdict[set2],setdict[set3])
        tic = time.time()
        print 'that took', tic-toc
        fname = self.outdir+self.runname+set1+set2+set3+'.out'
        nnn.write(fname,nnn)

if __name__ == '__main__':
    handler = NNNHandler(sys.argv[4])
    handler.run(sys.argv[1], sys.argv[2], sys.argv[3])
