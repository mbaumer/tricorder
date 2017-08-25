"""Module for handling datasets used in the 3PCF analysis.

This module loads, pre-processes, and plots metrics of datasets used
in the 3PCF analysis.

Later, we will have multiple classes inherit from BaseDataset, as in
    $ class RedmagicDataset (BaseDataset):
"""
from __future__ import division

import cPickle as pickle

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from scipy.stats import itemfreq
from sklearn.cluster import KMeans

# upgrade to this at some point:
# import healpix_util as hu
from make_data_randoms import index_to_radec, radec_to_index

buzzard_cosmo = FlatLambdaCDM(68.81, .295)

zvar_labels = {'ZSPEC': r'$z_{true}$',
               'ZREDMAGIC': r'z_{RM}',
               }

output_path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/redmagicHD/'


class BaseDataset (object):

    def __init__(self, datapath, maskpath, use_spec_z=True):

        self.datapath = datapath
        self.maskpath = maskpath

        self.data = None
        self.mask = None
        self.zmask = None

        self.n_jackknife = 30
        self.jk_labels = None

        self.min_z = None
        self.max_z = None

        self.nside = None
        self.nbar = None
        self.pixelized = None  # tuple: (pixel ra, pixel dec, counts)
        if use_spec_z:
            self.zvar = 'ZSPEC'
        else:
            self.zvar = 'ZREDMAGIC'

    @classmethod
    def fromfilename(cls, filename):
        """Initialize a BaseDataset from a pickle written by self.write."""
        data = pickle.load(open(filename, 'rb'))
        return data

    def load_data(self):
        self.data = fits.getdata(self.datapath)
        self.mask = hp.read_map(self.maskpath, 0, partial=True)
        self.zmask = hp.read_map(self.maskpath, 1, partial=True)

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
        self.mask[self.zmask < .5] = hp.UNSEEN

        dec, ra = index_to_radec(np.arange(hp.nside2npix(
            4096), dtype='int64'), 4096)  # masks have nside 4096
        bad_ra_dec = np.where(~((dec > min_dec) & (dec < max_dec)
                                & (ra > min_ra) & (ra < max_ra)))
        self.mask[bad_ra_dec] = hp.UNSEEN

    def apply_footprint(self, min_ra, max_ra, min_dec, max_dec):
        self._apply_footprint_data(min_ra, max_ra, min_dec, max_dec)
        self._apply_footprint_mask(min_ra, max_ra, min_dec, max_dec)

    def make_sky_map(self, **kwargs):
        """Plot a simple 2D histogram of RA,DEC for a sample.

        Mostly useful as a quick check of the footprint.
        """
        plt.hist2d(self.data['RA'], self.data['DEC'], bins=200, **kwargs)
        plt.xlabel('RA')
        plt.ylabel('DEC')

    def plot_n_z(self, **kwargs):
        """Plot the distribution of redshifts in the sample.

        Currently uses the zvar specified in the __init__
        """
        plt.hist(self.data[self.zvar], bins=50, **kwargs)
        plt.xlabel(zvar_labels[self.zvar])

    def apply_z_cut(self, min_z, max_z):
        self.min_z = min_z
        self.max_z = max_z
        self.data = self.data[self.data[self.zvar] > self.min_z]
        self.data = self.data[self.data[self.zvar] < self.max_z]

    def compute_new_jk_regions(self):
        """Add jackknife labels to the healpix pixels of the data locs.

        I don't think it should matter that these regions will be
        different for each redshift slice, dataset, etc.
        """
        data = zip(self.pixelized[0], self.pixelized[1])
        finder = KMeans(n_clusters=self.n_jackknife)
        self.jk_labels = finder.fit_predict(data)

    def write(self):
        """Save a BaseDataset instance as a pickle.

        These will be read in later by the 3PCF analysis.
        """
        # don't actually pickle out this huge stuff
        del self.data
        del self.mask
        del self.zmask

        name = output_path + str(self.zvar) + str(self.min_z) + '_' + \
            str(self.max_z) + 'nside' + str(self.nside) + 'nJack' \
            + str(self.n_jackknife)

        with open(name, 'wb') as pickle_file:
            pickle.dump(self, pickle_file)
