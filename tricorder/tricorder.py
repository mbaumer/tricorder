#!/home/mbaumer/anaconda2/bin/python2.7
import sys
sys.path.append('/home/mbaumer/anaconda2/bin/')
from astropy.io import fits
import treecorr
import sys
import numpy as np
import time

class NNNHandler (object):

    def __init__(self,runname):
        self.runname = runname
        self.datapath = '/scratch/PI/kipac/mbaumer/des/data/redmagic_sv_data.fits'
        self.randompath = '/scratch/PI/kipac/mbaumer/des/data/redmagic_sv_randoms.fits'
        self.outdir = '/scratch/PI/kipac/mbaumer/des/3pt_results/'
        #write config file

        configdict = {}
        configdict['min_r'] = 1
        configdict['max_r'] = 50
        configdict['n_rbins'] = 11
        
        configdict['min_u'] = 0.9
        configdict['max_u'] = 1.0
        configdict['n_ubins'] = 1
        
        configdict['min_v'] = -0.1
        configdict['max_v'] = 0.1
        configdict['n_vbins'] = 1

        self.config = configdict

    def load(self,runname):
        self.config = None #TODO load config dict!
        [ddd,ddr,drd,rdd,drr,rdr,rrd,rrr] = 8*[treecorr.NNNCorrelation(min_sep=self.config['min_r'], max_sep=self.config['max_r'], nbins=self.config['n_rbins'], min_u=self.config['min_u'], 
          max_u=self.config['max_u'], nubins=self.config['n_ubins'], min_v=self.config['min_v'], max_v=self.config['max_v'], nvbins=self.config['n_vbins'], bin_slop=1, sep_units='arcmin')]
        for nnn in [ddd,ddr,drd,rdd,drr,rdr,rrd,rrr]:
            #do something

    def run(self,set1,set2,set3):
        data = fits.getdata(self.datapath)
        randoms = fits.getdata(self.randompath)

        cat = treecorr.Catalog(ra=data['RA'], dec=data['DEC'], ra_units='degrees', dec_units='degrees')
        random_cat = treecorr.Catalog(ra=randoms['RA'], dec=randoms['DEC'], ra_units='degrees', dec_units='degrees')
        setdict = {'d':cat,'r':random_cat}
        nnn = treecorr.NNNCorrelation(min_sep=self.config['min_r'], max_sep=self.config['max_r'], nbins=self.config['n_rbins'], min_u=self.config['min_u'], 
          max_u=self.config['max_u'], nubins=self.config['n_ubins'], min_v=self.config['min_v'], max_v=self.config['max_v'], nvbins=self.config['n_vbins'], bin_slop=1, sep_units='arcmin')
        
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
