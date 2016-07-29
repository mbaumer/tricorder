import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from sklearn.cluster import KMeans
from sklearn.neighbors import KNeighborsRegressor

outdir = '/scratch/PI/kipac/mbaumer/des/3pt_results/'
runType = 'Y1_sims'
n_jackknife = 15

#define footprint dicts
sv_spt_footprint['min_ra'] = 60
sv_spt_footprint['max_ra'] = 92
sv_spt_footprint['min_dec'] = -61
sv_spt_footprint['max_dec'] = -40

y1_main_footprint['min_ra'] = 0
y1_main_footprint['max_ra'] = 360
y1_main_footprint['min_dec'] = -70
y1_main_footprint['max_dec'] = -35

def plotRegion(cat, indices, i):
    locs = np.where(indices == i)
    these_pts = cat[locs]
    plt.scatter(these_pts['RA'],these_pts['DEC'],color=np.random.rand(3,1),linewidths=0)
    
def plotFootprint(cat,indices):
    for i in np.unique(indices):
        plotRegion(cat, indices, i)

class NNNMunger(object):
    
    def __init__(self,runname):
        config_fname = outdir+runname+'.config'
        self.config = {}
        self.runname = runname
        self.config['runType'] = runType
        self.config['n_jackknife'] = n_jackknife

        if self.config['runType'] == 'SV':
            self.config['datapath'] = '/scratch/PI/kipac/mbaumer/des/data/redmagic_sv_data.fits'
            self.config['randompath'] = '/scratch/PI/kipac/mbaumer/des/data/redmagic_sv_5x_randoms.fits'
            self.footprint = sv_spt_footprint

        elif self.config['runType'] == 'Y1_sims':
            self.config['datapath'] = '/scratch/PI/kipac/mbaumer/des/data/redmagic_Y1_sims_data.fits'
            self.config['randompath'] = '/scratch/PI/kipac/mbaumer/des/data/redmagic_Y1_sims_5x_randoms.fits'
            self.footprint = y1_main_footprint

        else:
            raise IOError('invalid runType')
    
        self.config['min_z'] = .5
        self.config['max_z'] = .7

        #write it out so we remember what we did
        f = open(config_fname,'w')
        f.write(json.dumps(self.config))
        f.close()

    def applyCuts(self,cat):
        cat = cat[((cat['ZREDMAGIC'] > self.config['min_z']) & (cat['ZREDMAGIC'] < self.config['max_z']))]
        cat = cat[((cat['RA'] > self.footprint['min_ra']) & (cat['RA'] < self.footprint['max_ra']))]
        cat = cat[((cat['DEC'] > self.footprint['min_dec']) & (cat['DEC'] < self.footprint['max_dec']))]
        #later, color, magnitude cuts...
        return cat

    def prepareCatalog(self):
        data = fits.getdata(self.config['datapath'])
        randoms = fits.getdata(self.config['randompath'])

        data = applyCuts(data)
        randoms = applyCuts(randoms)

        self.cat = cat
        self.random_cat = random_cat

    def addJKRegionLabels(self):
        data = zip(self.cat['RA'],self.cat['DEC'])
        randoms = zip(self.random_cat['RA'],self.random_cat['DEC'])
        
        finder = KMeans(n_clusters=self.config['n_jackknife'])
        self.data_jk_indices = finder.fit_predict(data)
        
        nbrs = KNeighborsRegressor(n_neighbors=1)
        nbrs.fit(data,self.data_jk_indices)
        self.random_jk_indices = nbrs.predict(randoms)

    def saveInputsWithJKLabels(self):
        #god i hate fits tables
        new_data_col = fits.Column(name='jk_region', format='I', array=self.data_jk_indices)
        appended_data = fits.BinTableHDU.from_columns(self.cat.columns+new_data_col)
        appended_data.saveto(outdir+self.runname+'input_data.fits')

        new_random_col = fits.Column(name='jk_region', format='I', array=self.random_jk_indices)
        appended_randoms = fits.BinTableHDU.from_columns(self.random_cat.columns+new_random_col)
        appended_randoms.saveto(outdir+self.runname+'input_randoms.fits')

# ------------------------------------------------------

if __name__ == '__main__':
    munger = NNNProcessor(sys.argv[1])
    munger.prepareCatalog()
    munger.makeJKRegionLabels()
    munger.saveInputsWithJKLabels()
    


