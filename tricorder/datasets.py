'''datasets.py

This module loads, pre-processes, and plots metrics of datasets used
in the 3PCF analysis.

Later, we will have multiple classes inherit from BaseDataset, as in
    $ class RedmagicDataset (BaseDataset):
'''
from __future__ import division

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from scipy.stats import itemfreq

# upgrade to this at some point:
# import healpix_util as hu
from make_data_randoms import index_to_radec, radec_to_index

buzzard_cosmo = FlatLambdaCDM(68.81, .295)

zvar_labels = {'ZSPEC': r'$z_{true}$',
               'ZREDMAGIC': r'z_{RM}',
               }


class BaseDataset (object):

    def __init__(self, datapath, maskpath, use_spec_z=True):

        self.data = fits.getdata(datapath)
        self.mask = hp.read_map(maskpath, 0, partial=True)
        self.zmask = hp.read_map(maskpath, 1, partial=True)

        self.nside = None
        self.nbar = None
        self.pixelized = None  # takes format (pixel ra, pixel dec, counts)

        if use_spec_z == True:
            self.zvar = 'ZSPEC'
        else:
            self.zvar = 'ZREDMAGIC'

    def pixelize_at_target_nside(self, nside):
        # this also effectively applies the mask to the data
        self.nside = nside
        mask_targetnside = hp.pixelfunc.ud_grade(
            self.mask, pess=True, nside_out=self.nside)
        gal_index_targetnside = radec_to_index(
            self.data['DEC'], self.data['RA'], self.nside)
        counts = itemfreq(gal_index_targetnside)
        full_map_counts = np.zeros(hp.nside2npix(self.nside))
        full_map_counts[counts[:, 0]] = counts[:, 1]

        good_counts = full_map_counts[np.where(mask_targetnside != hp.UNSEEN)]
        good_fracs = mask_targetnside[np.where(mask_targetnside != hp.UNSEEN)]

        self.nbar = np.average(good_counts, weights=good_fracs)
        print 'nbar is', self.nbar, 'galaxies per pixel'
        pixels_to_count = np.where(mask_targetnside != hp.UNSEEN)
        dec, ra = index_to_radec(pixels_to_count, self.nside)
        final_counts = full_map_counts[pixels_to_count]
        self.pixelized = (ra[0], dec[0], final_counts)

    def _apply_footprint_data(self, min_ra, max_ra, min_dec, max_dec):
        self.data = self.data[self.data['RA'] > min_ra]
        self.data = self.data[self.data['RA'] < max_ra]
        self.data = self.data[self.data['DEC'] > min_dec]
        self.data = self.data[self.data['DEC'] < max_dec]

    def _apply_footprint_mask(self, min_ra, max_ra, min_dec, max_dec):

        self.mask[self.mask < .95] = hp.UNSEEN
        self.mask[self.zmask < .5]= hp.UNSEEN

        dec, ra = index_to_radec(np.arange(hp.nside2npix(4096),dtype='int64'), 4096)  # masks have nside 4096
        bad_ra_dec= np.where(~((dec > min_dec) & (dec < max_dec)
                                & (ra > min_ra) & (ra < max_ra)))
        self.mask[bad_ra_dec]= hp.UNSEEN

    def apply_footprint(self, min_ra, max_ra, min_dec, max_dec):
        self._apply_footprint_data(min_ra, max_ra, min_dec, max_dec)
        self._apply_footprint_mask(min_ra, max_ra, min_dec, max_dec)

    def make_sky_map(self,**kwargs):
        plt.hist2d(self.data['RA'], self.data['DEC'], bins=200,**kwargs)
        plt.xlabel('RA')
        plt.ylabel('DEC')

    def plot_n_z(self,**kwargs):
        plt.hist(self.data[self.zvar], bins=50, **kwargs)
        plt.xlabel(zvar_labels[self.zvar])

    def apply_z_cut(self, min_z, max_z):
        self.data= self.data[self.data[self.zvar] > min_z]
        self.data= self.data[self.data[self.zvar] < max_z]
