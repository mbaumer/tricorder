#!/u/ki/mbaumer/anaconda/bin/python
import readGadgetSnapshot
from astropy import coordinates
from glob import glob
import numpy as np
from astropy.io import fits

box_size = 2600 #1050, 2600, 4000, or 6000
fname_vec = glob('/nfs/slac/g/ki/ki21/cosmo/BCCSims/Chinchillas-midway/Chinchilla-1/Lb'+str(box_size)+'/lightcone000/snapshot_*')

def make_gadget_catalog(frac):
    firstFile = True 
    for i,fname in enumerate(fname_vec):
        print i, fname
        data = readGadgetSnapshot.readGadgetSnapshot(fname,read_pos=True)[1]
        coords = coordinates.SkyCoord(x=data[:,0],y=data[:,1],z=data[:,2],unit='Mpc', representation='cartesian')
        coords = coords[np.random.rand(len(coords)) < frac]
        if firstFile:
            fullcat = coords
            firstFile = False
        else:
            newcat = coords
            fullcat = coordinates.concatenate([fullcat,newcat])
    return fullcat

def save_dm_catalog(fullcat,frac):
    hdu = fits.open('/nfs/slac/g/ki/ki19/des/mbaumer/3pt_data/redmagic_data.fits')
    c1 = fits.Column(name='RA', format='D', array=fullcat.ra)
    c2 = fits.Column(name='DEC', format='D', array=fullcat.dec)
    c3 = fits.Column(name='DISTANCE', format='D', array=fullcat.distance)
    coldefs = fits.ColDefs([c1, c2, c3])
    tbhdu = fits.BinTableHDU.from_columns(coldefs)
    hdu[1] = tbhdu
    hdu.writeto('/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/dark_matter/dm_cat_'+str(box_size)+'h-1Mpc_data_10pct_downsampled.fits',clobber=True)
    
def save_random_catalog(cat,oversample_factor):
    random_ra = 90*np.random.rand(oversample_factor*len(cat.ra))
    random_dec = 90/np.pi*np.arccos(np.random.rand(oversample_factor*len(cat.ra))-1)
    hdu = fits.open('/nfs/slac/g/ki/ki19/des/mbaumer/3pt_data/redmagic_data.fits')
    c1 = fits.Column(name='RA', format='D', array=random_ra)
    c2 = fits.Column(name='DEC', format='D', array=random_dec)
    rand_distances = np.random.choice(cat.distance,len(random_ra))
    c3 = fits.Column(name='DISTANCE', format='D', array=rand_distances)
    coldefs = fits.ColDefs([c1, c2, c3])
    tbhdu2 = fits.BinTableHDU.from_columns(coldefs)
    hdu[1] = tbhdu2
    hdu.writeto('dm_cat_'+str(box_size)+'Mpc_randoms_0.fits',clobber=True)
    
def make_DM_catalogs(frac=0.1,oversample_factor=5):
    gadget_cat = make_gadget_catalog(frac)
    save_dm_catalog(gadget_cat,frac)
    #save_random_catalog(gadget_cat,oversample_factor)
    
if __name__ == '__main__':
    make_DM_catalogs()