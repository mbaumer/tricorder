#!/u/ki/mbaumer/anaconda/bin/python
from astropy.io import fits
import treecorr
import sys
import numpy as np

def run3pt(set1='d',set2='d',set3='d',Rfrac='.95'):
    data =
fits.getdata('/nfs/slac/g/ki/ki18/des/cpd/cross_correlation_taskforce_buzzard_datasets/redmagic.fits')
    randoms =
fits.getdata('/nfs/slac/g/ki/ki18/des/cpd/cross_correlation_taskforce_buzzard_datasets/redmagic_randoms.fits')
    mask = [np.random.uniform(size=len(randoms)) < float(Rfrac)]
    randoms = randoms[mask]
    cat = treecorr.Catalog(ra=data['RA'], dec=data['DEC'], ra_units='degrees',
dec_units='degrees')
    random_cat = treecorr.Catalog(ra=randoms['RA'], dec=randoms['DEC'],
ra_units='degrees', dec_units='degrees')
    setdict = {'d':cat,'r':random_cat}
    nnn = treecorr.NNNCorrelation(min_sep=1., max_sep=30., nbins=8, min_u=0.9,
max_u=1.1, nubins=1, min_v=-0.1, max_v=0.1, nvbins=1, sep_units='arcmin')
    print 'starting processing!'
    nnn.process(setdict[set1],setdict[set2],setdict[set3])
    fname =
'/nfs/slac/g/ki/ki19/des/mbaumer/3pt_runs/redmagic_sim_'+set1+set2+set3+'_Rfrac'+Rfrac+'.out'
    nnn.write(fname,nnn)

if __name__ == '__main__':
    run3pt(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
