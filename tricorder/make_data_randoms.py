'''
This code is run from the command line to generate masked data and random sets with a fixed oversampling w.r.t. the input data.

Parameters
==========

data_filename : str
    The relative (or global) path of the data file to use.
    
mask_filename : str
    The relative (or global) path of the combined zlim/fracgood mask.
    
zvar : str
    The redshift variable to use. Typically 'ZSPEC' or 'ZREDMAGIC'
    
oversampling : float
    The factor relative to the amount of data to determine the amount of randoms to be produced. Typically between 1-20.
    
outname : str
    The name for the output datasets, and the diagnostic figures.
    
Notes
=====

This should probably become object-oriented as some point, with data production managed by a DatasetMaker class, e.g.

'''

import sys
import numpy as np
from astropy.io import fits
import healpy as hp
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
buzzard_cosmo = FlatLambdaCDM(68.81,.295)
#redmine: fiducial Y1 cosmology

def load_mask_data(mask_filename):
    mask = fits.getdata(mask_filename)
    map_fracgood_indices, map_fracgood = convert_pixel_to_map(mask,valkey='FRACGOOD',pixkey='HPIX',return_seen=False)
    map_zmax_indices, map_zmax = convert_pixel_to_map(mask,valkey='ZMAX',pixkey='HPIX',return_seen=False)

    #mask pixels with fracgood < .95 (not sure where this cut comes from, but Chris did it.)
    map_fracgood_indices = map_fracgood_indices[map_fracgood > 0.95]
    map_fracgood = map_fracgood[map_fracgood > 0.95]
    return map_fracgood_indices, map_fracgood, map_zmax_indices, map_zmax

def convert_pixel_to_map(pixel, valkey='SIGNAL', pixkey='PIXEL',
                         nested=False, nside=4096, return_seen=True):
    ycoord = pixel[pixkey]
    if nested:
        # convert to ring
        ycoord = hp.nest2ring(nside, ycoord)

    m = np.zeros(hp.nside2npix(nside), dtype=pixel[pixkey].dtype) + hp.UNSEEN
    m[ycoord] = pixel[valkey]

    if return_seen:
        indices, = np.where(m != hp.UNSEEN)
        m = m[indices]
        return indices, m
    else:
        return np.arange(m.size), m

def index_to_radec(index, nside):
    theta,phi=hp.pixelfunc.pix2ang(nside, index)

    # dec, ra
    return -np.degrees(theta-np.pi/2.), np.degrees(phi)

def radec_to_index(dec, ra, nside):
    return hp.pixelfunc.ang2pix(nside, np.radians(-dec+90.), np.radians(ra))

def generate_randoms_radec(minra, maxra, mindec, maxdec, Ngen, raoffset=0):
    r = 1.0
    # this z is not redshift!
    zmin = r * np.sin(np.pi * mindec / 180.)
    zmax = r * np.sin(np.pi * maxdec / 180.)

    # parity transform from usual, but let's not worry about that
    phimin = np.pi / 180. * (minra - 180 + raoffset)
    phimax = np.pi / 180. * (maxra - 180 + raoffset)

    # generate ra and dec
    z_coord = np.random.uniform(zmin, zmax, Ngen)  # not redshift!
    phi = np.random.uniform(phimin, phimax, Ngen)
    dec_rad = np.arcsin(z_coord / r)

    # convert to ra and dec
    ra = phi * 180 / np.pi + 180 - raoffset
    dec = dec_rad * 180 / np.pi
    return ra, dec

def generate_randoms(dat, zvar, raoffset, Ntot,
                     map_fracgood_indices, map_fracgood,
                     map_zmax_indices, map_zmax, nside,
                     Ngen=1000000, Ntries_max=10000):

    Ncurrent = 0
    Ntry = 0

    minra = min(dat['RA'])
    maxra = max(dat['RA'])
    mindec = min(dat['DEC'])
    maxdec = max(dat['DEC'])
    
    zdist = dat[zvar]

    random_ra = []
    random_dec = []
    random_z = []
    while ((Ncurrent < Ntot) & (Ntry < Ntries_max)):
        if Ntry % 100 == 1:
            print(Ntry, Ntries_max, Ncurrent, Ntot, Ngen)
        # generate random ra and dec
        ra_i, dec_i = generate_randoms_radec(minra, maxra, mindec, maxdec, Ngen, raoffset=raoffset)
        indices_i = radec_to_index(dec_i, ra_i, nside)
        # cut ra and dec on fracgood
        good_index = np.in1d(indices_i, map_fracgood_indices)
        ra_i = ra_i[good_index]
        dec_i = dec_i[good_index]
        indices_i = indices_i[good_index]
        # generate zdist
        z_i = np.random.choice(zdist, len(ra_i))
        # get associated zmax_i
        # magic. assumes map_zmax_indices already sorted (it is)
        # simple way
        zmax_indices_i = np.searchsorted(map_zmax_indices, indices_i)
        # harder way if something goes wrong...
        # zmax_indices_i = np.searchsorted(map_zmax_indices_sorted, indices_i)
        # zmax_indices = np.take(argsort_map_zmax_indices, zmax_indices_i, mode='clip')
        # zmax_indices_i = zmax_indices_i[map_zmax_indices[zmax_indices] == indices_i]

        zmax_i = map_zmax[zmax_indices_i]
        # cut on zmax
        if zvar == 'DISTANCE':
            good_index = buzzard_cosmo.comoving_distance(zmax_i).value*buzzard_cosmo.h >= z_i
        else:
            good_index = zmax_i >= z_i
        
        ra_i = ra_i[good_index]
        dec_i = dec_i[good_index]
        z_i = z_i[good_index]

        random_ra += list(ra_i)
        random_dec += list(dec_i)
        random_z += list(z_i)

        Ncurrent = len(random_ra)
        Ntry += 1
    if Ntry >= Ntries_max:
        print('Warning! We gave up after {0} tries, finding {1} objects instead of the desired {2} objects!'.format(Ntry, Ncurrent, Ntot))

    randoms = np.zeros(len(random_ra), dtype=[('RA', '>f4'), ('DEC', '>f4'), (zvar, '>f4')])
    randoms['RA'] = random_ra
    randoms['DEC'] = random_dec
    randoms[zvar] = random_z
    
    if len(randoms) > Ntot:
        inds_to_keep = np.random.choice(np.arange(len(randoms)), size=Ntot,
                                 replace=False)
        randoms = randoms[inds_to_keep]
        
    return randoms

def mask_data(data, zvar, map_fracgood_indices, map_fracgood, map_zmax_indices, map_zmax, nside):
    indices = radec_to_index(data['DEC'], data['RA'], nside=nside)
    good_index = np.in1d(indices, map_fracgood_indices[map_fracgood > 0.95])
    data = data[good_index] #cut data on simulated depth map
    indices = indices[good_index]

    map_zmax_indices = np.searchsorted(map_zmax_indices, indices)
    map_zmax = map_zmax[map_zmax_indices]
    
    if zvar == 'DISTANCE':
        good_index = buzzard_cosmo.comoving_distance(map_zmax).value*buzzard_cosmo.h >= data[zvar]
    else:
        good_index = map_zmax >= data[zvar]
    
    data = data[good_index]
    indices = indices[good_index]
    
    masked_data = np.zeros(len(data), dtype=[('RA', '>f4'), ('DEC', '>f4'), (zvar, '>f4')])
    masked_data['RA'] = data['RA']
    masked_data['DEC'] = data['DEC']
    masked_data[zvar] = data[zvar]
        
    return masked_data
    
def test_data_randoms(data,randoms,zvar,dname,rname,nbins=200):
    fig, axarr = plt.subplots(2,2,figsize=(12,8))
    plt.suptitle('data: '+dname+'\n'+'mask: '+rname)
    axarr[0,0].hist2d(data['RA'],data['DEC'],bins=nbins)
    axarr[0,0].set_title('Data')
    axarr[0,1].hist2d(randoms['RA'],randoms['DEC'],bins=nbins)
    axarr[0,1].set_title('Randoms')
    axarr[1,0].hist(data[zvar])
    axarr[1,0].set_title('Data '+zvar)
    axarr[1,1].hist(randoms[zvar])
    axarr[1,1].set_title('Randoms '+zvar)
    return fig
    
def generate_datasets(data_filename, mask_filename, zvar, oversampling,outname):
    #mask_filename = 'raw/simulation/buzzard-v1.1-y1a1-full_run_redmapper_v6.4.13_redmagic_highdens_0.5-10_vlim_zmask.fit'
    
    nside = fits.getheader(mask_filename,1)['NSIDE']
    assert nside == 4096
    
    #data_filename = 'raw/simulation/buzzard-v1.1-y1a1-full_run_redmapper_v6.4.13_redmagic_highdens_0.5-10.fit'
    data = fits.getdata(data_filename)
    
    if zvar == 'DISTANCE':
        data = data[np.random.rand(len(data)) < .1] #total effective downsample of .1% (1% from initial gadget read)
        data = data[data['DISTANCE'] < 2600] # distance < 2600 Mpc/h (units already correct)
        data['DEC'] = -data['DEC']
        
    data = data[(data['RA'] > 0) & (data['RA'] < 90)]
    data = data[(data['DEC'] > -60) & (data['DEC'] < -40)]
    
    map_fracgood_indices, map_fracgood, map_zmax_indices, map_zmax = load_mask_data(mask_filename)
    masked_data = mask_data(data, zvar, map_fracgood_indices, map_fracgood, map_zmax_indices, map_zmax, nside)
    
    Ntot = oversampling*len(masked_data)
    randoms = generate_randoms(data, zvar, 0, Ntot,
                     map_fracgood_indices, map_fracgood,
                     map_zmax_indices, map_zmax, nside,
                     Ngen=1000000, Ntries_max=10000)
    
    np.save(outname+'_data.npy',masked_data)
    np.save(outname+'_randoms.npy',randoms)
    
    fig = test_data_randoms(masked_data,randoms,zvar,data_filename,mask_filename)
    fig.savefig(outname+'.png',dpi=fig.dpi)
    
if __name__ == '__main__':
    if len(sys.argv)-1 != 5:
        print 'only got', len(sys.argv)-1, 'arguments.'
        raise IndexError('takes 5 args: data, mask, zvar, oversampling multiplier, outname')
    generate_datasets(sys.argv[1],sys.argv[2],sys.argv[3],int(sys.argv[4]),sys.argv[5])
