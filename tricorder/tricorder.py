#!/home/mbaumer/anaconda2/bin/python2.7
import sys
sys.path.append('/home/mbaumer/anaconda2/bin/')
from astropy.io import fits
import treecorr
import sys
import numpy as np
import time

def run3pt(set1='d',set2='d',set3='d'):
    data = fits.getdata('/scratch/PI/kipac/mbaumer/des/data/redmagic_sv_data.fits')
    randoms = fits.getdata('/scratch/PI/kipac/mbaumer/des/data/redmagic_sv_randoms.fits')
    #data = fits.getdata('/scratch/PI/kipac/mbaumer/des/data/redmagic_Y1_sims_data.fits')
    #randoms = fits.getdata('/scratch/PI/kipac/mbaumer/des/data/redmagic_Y1_sims_randoms.fits')
    #mask = [np.random.uniform(size=len(randoms)) < 0.05]
    #randoms = randoms[mask]
    cat = treecorr.Catalog(ra=data['RA'], dec=data['DEC'], ra_units='degrees', dec_units='degrees')
    random_cat = treecorr.Catalog(ra=randoms['RA'], dec=randoms['DEC'], ra_units='degrees', dec_units='degrees')
    setdict = {'d':cat,'r':random_cat}
    nnn = treecorr.NNNCorrelation(min_sep=1, max_sep=50, nbins=11, min_u=0.9, 
      max_u=1.0, nubins=1, min_v=-0.1, max_v=0.1, nvbins=1, bin_slop=0.01, sep_units='arcmin')
    print 'starting processing!'
    toc = time.time()
    if ((set1 == set2) and (set2 == set3)):
        nnn.process(setdict[set1])
        #this saves time even though it shouldn't in theory...
    else:
        nnn.process(setdict[set1],setdict[set2],setdict[set3])
    tic = time.time()
    print 'that took', tic-toc
    fname = '/scratch/PI/kipac/mbaumer/des/3pt_results/'+set1+set2+set3+'sv_test.out'
    nnn.write(fname,nnn)

if __name__ == '__main__':
    run3pt(sys.argv[1], sys.argv[2], sys.argv[3])
