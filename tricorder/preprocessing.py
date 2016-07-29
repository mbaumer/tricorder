import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from sklearn.cluster import KMeans
from sklearn.neighbors import KNeighborsRegressor

outdir = '/scratch/PI/kipac/mbaumer/des/3pt_results/'

#define footprint dicts
sv_spt_footprint['min_ra'] = 60
sv_spt_footprint['max_ra'] = 92
sv_spt_footprint['min_dec'] = -61
sv_spt_footprint['max_dec'] = -40

y1_main_footprint['min_ra'] = 0
y1_main_footprint['max_ra'] = 360
y1_main_footprint['min_dec'] = -70
y1_main_footprint['max_dec'] = -35

class NNNMunger(object):
    
    def __init__(self):
        self.footprint = y1_main_footprint
        self.config['min_z'] = .5
        self.config['max_z'] = .7


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

def makeJKRegionLabels(cat,random_cat,n_jackknife=15):
    data = zip(cat['RA'],cat['DEC'])
    randoms = zip(random_cat['RA'],random_cat['DEC'])
    
    finder = KMeans(n_clusters=n_jackknife)
    data_indices = finder.fit_predict(data)
    
    nbrs = KNeighborsRegressor(n_neighbors=1)
    nbrs.fit(data,data_indices)
    random_indices = nbrs.predict(randoms)
    
    return data_indices, random_indices

def plotRegion(cat, indices, i):
    locs = np.where(indices == i)
    these_pts = cat[locs]
    plt.scatter(these_pts['RA'],these_pts['DEC'],color=np.random.rand(3,1),linewidths=0)
    
def plotFootprint(cat,indices):
    for i in np.unique(indices):
        plotRegion(cat, indices, i)