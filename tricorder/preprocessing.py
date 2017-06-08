import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from sklearn.cluster import KMeans
from sklearn.neighbors import KNeighborsRegressor

n_jackknife = 30

def plotRegion(cat, indices, i):
    locs = np.where(indices == i)
    these_pts = cat[locs]
    plt.scatter(these_pts['RA'],these_pts['DEC'],color=np.random.rand(3,1),linewidths=0)
    
def plotFootprint(cat,indices):
    for i in np.unique(indices):
        plotRegion(cat, indices, i)

class NNNMunger(object):
    
    def __init__(self,runname):
        def __init__(self,zvar,min_z,delta_z,metric):
        self.zvar = zvar
        self.min_z = min_z
        self.delta_z = delta_z
        self.max_z = self.min_z + self.delta_z
        self.metric = metric
        self.runname = self.zvar+str(self.min_z)+'_deltaz'+str(self.delta_z)+'_'+self.metric
        with open(outdir+self.runname+'.yaml') as f:
            self.config = yaml.load(f.read())
        self.data = np.load(self.config['data_path'])
        self.randoms = np.load(self.config['randoms_path'])
        assert self.runname == self.config['runname']

    def addJKRegionLabels(self):
        data = zip(self.data['RA'],self.data['DEC'])
        randoms = zip(self.randoms['RA'],self.randoms['DEC'])
        
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
    


