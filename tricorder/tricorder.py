#!/home/mbaumer/anaconda2/bin/python2.7
import sys, os.path
sys.path.append('/home/mbaumer/anaconda2/bin/')
from astropy.io import fits
import treecorr
import sys
import numpy as np
import time
import json

outdir = '/scratch/PI/kipac/mbaumer/des/3pt_results/'

class NNNPlotter (object):
    def __init__(self,runname):
        config_fname = outdir+runname+'.config'
        try: 
            with open(config_fname) as f:
                configdict = json.loads(f.read())
        except IOError:
            print 'results for '+runname+' not found in '+outdir

    def load(self):
        [ddd,ddr,drd,rdd,drr,rdr,rrd,rrr] = 8*[treecorr.NNNCorrelation(config=self.config)]
        for nnn in [ddd,ddr,drd,rdd,drr,rdr,rrd,rrr]:
            pass

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

        configdict['min_z'] = .5
        configdict['max_z'] = .7
        #SV SPT-E footprint
        configdict['min_ra'] = 60
        configdict['max_ra'] = 92
        configdict['min_dec'] = -61
        configdict['max_dec'] = -40
	#Y1 main footprint
	#configdict['min_ra'] = 0
        #configdict['max_ra'] = 360
        #configdict['min_dec'] = -70
        #configdict['max_dec'] = -35


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

    def prepareCatalog(self):
        data = fits.getdata(self.config['datapath'])
        randoms = fits.getdata(self.config['randompath'])

        data = data[((data['ZREDMAGIC'] > self.config['min_z']) & (data['ZREDMAGIC'] < self.config['max_z']))]
        data = data[((data['RA'] > self.config['min_ra']) & (data['RA'] < self.config['max_ra']))]
        data = data[((data['DEC'] > self.config['min_dec']) & (data['DEC'] < self.config['max_dec']))]

        randoms = randoms[((randoms['Z'] > self.config['min_z']) & (randoms['Z'] < self.config['max_z']))]
        randoms = randoms[((randoms['RA'] > self.config['min_ra']) & (randoms['RA'] < self.config['max_ra']))]
        randoms = randoms[((randoms['DEC'] > self.config['min_dec']) & (randoms['DEC'] < self.config['max_dec']))]

        cat = treecorr.Catalog(ra=data['RA'], dec=data['DEC'], ra_units='degrees', dec_units='degrees')
        random_cat = treecorr.Catalog(ra=randoms['RA'], dec=randoms['DEC'], ra_units='degrees', dec_units='degrees')
        return cat, random_cat

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
