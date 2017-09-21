"""Module for handling datasets used in the 3PCF analysis.

This module loads, pre-processes, and plots metrics of datasets used
in the 3PCF analysis.

Later, we will have multiple classes inherit from BaseDataset, as in
    $ class RedmagicDataset (BaseDataset):
"""
from __future__ import division

import cPickle as pickle
from copy import copy

# upgrade to this at some point:
import healpix_util as hu
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from scipy.stats import itemfreq
from sklearn.cluster import KMeans

#from make_data_randoms import index_to_radec, radec_to_index

buzzard_cosmo = FlatLambdaCDM(68.81, .295)

zvar_labels = {'ZSPEC': r'$z_{true}$',
               'ZREDMAGIC': r'z_{RM}',
               'DISTANCE': r' Comoving Distance (Mpc/$h$)',
               'REDSHIFT': r'$z_{LSS}$'
               }

mock = 'Buzzard_v1.1'

output_path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/' + mock + '/'


class BaseDataset (object):

    def __init__(self):

        if not hasattr(self, 'sample_type'):
            self.sample_type = 'unk'

        self.output_path = output_path + self.sample_type + '/data/'

        self.data = None
        self.mask = None
        self.zmask = None

        self.n_jackknife = 30
        self.jk_labels = None

        self.min_z = None
        self.max_z = None

        self.randoms = None

        self.nside = None
        self.nbar = None
        self.pixelized = None  # tuple: (pixel ra, pixel dec, counts)

    @classmethod
    def fromfilename(cls, filename):
        """Initialize a BaseDataset from a pickle written by self.write."""
        data = pickle.load(open(filename, 'rb'))
        return data

    def load_data(self):
        raise NotImplementedError(
            'Need to specify sample class to know how to load!')

    def pixelize_at_target_nside(self, nside):
        # this also effectively applies the mask to the data
        self.nside = nside
        hpix = hu.HealPix('ring', self.nside)
        mask_targetnside = hp.pixelfunc.ud_grade(
            self.mask, pess=True, nside_out=self.nside)
        gal_index_targetnside = hpix.eq2pix(self.data['RA'], self.data['DEC'])

        # prune data that's in a bad part of the mask
        self.data = self.data[mask_targetnside[gal_index_targetnside] != hp.UNSEEN]

        counts = itemfreq(gal_index_targetnside)
        full_map_counts = np.zeros(hp.nside2npix(self.nside))
        full_map_counts[counts[:, 0]] = counts[:, 1]

        good_counts = full_map_counts[np.where(mask_targetnside != hp.UNSEEN)]
        good_fracs = mask_targetnside[np.where(mask_targetnside != hp.UNSEEN)]

        self.nbar = np.average(good_counts, weights=good_fracs)
        print 'nbar is', self.nbar, 'galaxies per pixel'
        pixels_to_count = np.where(mask_targetnside != hp.UNSEEN)
        ra, dec = hpix.pix2eq(pixels_to_count)
        final_counts = full_map_counts[pixels_to_count]
        self.pixelized = (ra, dec, final_counts)

    def _apply_footprint_data(self, min_ra, max_ra, min_dec, max_dec):
        self.data = self.data[self.data['RA'] > min_ra]
        self.data = self.data[self.data['RA'] < max_ra]
        self.data = self.data[self.data['DEC'] > min_dec]
        self.data = self.data[self.data['DEC'] < max_dec]

    def _apply_footprint_mask(self, min_ra, max_ra, min_dec, max_dec):

        self.mask[self.mask < .95] = hp.UNSEEN
        self.mask[self.zmask < .6] = hp.UNSEEN

        hpix = hu.HealPix('ring', 4096)
        ra, dec = hpix.pix2eq(np.arange(hp.nside2npix(
            4096), dtype='int64'))  # masks have nside 4096
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

    def generate_randoms(self, oversamp):
        """Generate randoms with oversamp times the density of the data."""
        # uses more memory, but don't want this method to edit real mask.
        #zero_mask = copy(self.mask)
        zero_mask = self.mask
        zero_mask[zero_mask == hp.UNSEEN] = 0
        dmask = hu.DensityMap('ring', self.mask)
        randoms = {}
        N_randoms = len(self.data) * oversamp
        randoms['RA'], randoms['DEC'] = dmask.genrand(N_randoms)
        randoms[self.zvar] = np.random.choice(self.data[self.zvar], N_randoms)
        self.randoms = randoms

    def make(self, nside, oversamp, min_z, max_z, min_ra, max_ra, min_dec, max_dec, doWrite=False):
        self.load_data()
        self.apply_z_cut(min_z, max_z)
        self.apply_footprint(min_ra, max_ra, min_dec, max_dec)
        self.pixelize_at_target_nside(nside)
        self.generate_randoms(oversamp)
        self.compute_new_jk_regions()
        if doWrite:
            self.write()

    def write(self):
        """Save a BaseDataset instance as a pickle.

        These will be read in later by the 3PCF analysis.
        """
        # don't actually pickle out this huge stuff
        # del self.data  # don't delete randoms either
        del self.mask
        del self.zmask

        name = self.output_path + str(self.zvar) + \
            str(self.min_z) + '_' + str(self.max_z) + \
            'nside' + str(self.nside) + 'nJack' \
            + str(self.n_jackknife)

        with open(name, 'wb') as pickle_file:
            pickle.dump(self, pickle_file)


class RedmagicDataset(BaseDataset):
    def __init__(self, use_spec_z=True):
        self.sample_type = 'redmagicHD'
        if use_spec_z:
            self.zvar = 'ZSPEC'
        else:
            self.zvar = 'ZREDMAGIC'
        self.datapath = '/nfs/slac/des/fs1/g/sims/jderose/addgals/catalogs/Buzzard/Catalog_v1.1/y1a1_mock_analysis/mock1/redmagic/buzzard-v1.1-y1a1-spt_mock1_redmapper_v6.4.13_redmagic_higherlum_1.5-01.fit'
        self.maskpath = '/nfs/slac/des/fs1/g/sims/jderose/addgals/catalogs/Buzzard/Catalog_v1.1/y1a1_mock_analysis/mock1/redmagic/buzzard-v1.1-y1a1-spt_mock1_run_redmapper_v6.4.13_redmagic_1.5_vlim_zmask.fit'
        super(RedmagicDataset, self).__init__()

    def load_data(self):
        self.data = fits.getdata(self.datapath)
        self.mask = hp.read_map(self.maskpath, 0, partial=True)
        self.zmask = hp.read_map(self.maskpath, 1, partial=True)


class DMDataset(BaseDataset):
    def __init__(self, datapath, maskpath):
        self.sample_type = 'dark_matter'
        self.zvar = 'DISTANCE'
        self.datapath = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/dark_matter/dm_cat_2600Mpc_data.fits'
        self.maskpath = '/nfs/slac/des/fs1/g/sims/jderose/addgals/catalogs/Buzzard/Catalog_v1.1/y1a1_mock_analysis/mock1/redmagic/buzzard-v1.1-y1a1-spt_mock1_run_redmapper_v6.4.13_redmagic_1.5_vlim_zmask.fit'
        super(DMDataset, self).__init__()

    def apply_z_cut(self, min_z, max_z):
        # compute dist limits in Mpc/h to agree w DM DISTANCE
        self.min_z = min_z
        self.max_z = max_z
        min_dist = buzzard_cosmo.comoving_distance(
            min_z).value * buzzard_cosmo.h
        max_dist = buzzard_cosmo.comoving_distance(
            max_z).value * buzzard_cosmo.h
        self.data = self.data[self.data[self.zvar] > min_dist]
        self.data = self.data[self.data[self.zvar] < max_dist]

    def load_data(self):
        self.data = fits.getdata(self.datapath)
        # to move DM octant into south to overlap with DES mask.
        self.data['DEC'] = -self.data['DEC']
        self.mask = hp.read_map(self.maskpath, 0, partial=True)
        self.zmask = hp.read_map(self.maskpath, 1, partial=True)


class LSSDataset(BaseDataset):
    def __init__(self):
        self.sample_type = 'lss_sample'
        self.zvar = 'REDSHIFT'
        super(LSSDataset, self).__init__()

    def load_data(self):
        mock = 1
        gold = fits.getdata(
            '/nfs/slac/des/fs1/g/sims/jderose/addgals/catalogs/Buzzard/Catalog_v1.1/y1a1_mock_analysis/mock{0}/mergedcats/Buzzard_v1.1_{0}_gold.fits.gz'.format(mock))
        lss = gold[(gold['SAMPLE'] == 1) + (gold['SAMPLE'] == 3)]
        c1 = fits.Column(name='RA', array=lss['ra'], format='E')
        c2 = fits.Column(name='DEC', array=lss['dec'], format='E')
        c3 = fits.Column(name='REDSHIFT', array=lss['redshift'], format='E')
        t = fits.BinTableHDU.from_columns([c1, c2, c3])
        self.data = t.data

        footprint = hp.read_map(
            '/nfs/slac/g/ki/ki23/des/jderose/SkyFactory/chinchilla-herd/Chinchilla-1/sampleselection/y1a1_gold_1.0.2_wide_footprint_4096.fits.gz')
        badmask = hp.read_map(
            '/nfs/slac/g/ki/ki23/des/jderose/SkyFactory/chinchilla-herd/Chinchilla-1/sampleselection/y1a1_gold_1.0.2_wide_badmask_4096.fits.gz')
        nodup = hp.read_map(
            '/nfs/slac/des/fs1/g/sims/jderose/addgals/catalogs/Buzzard/Catalog_v1.1/depthmaps/nodup_ssmask.fits')
        lss_mask = (footprint >= 1) * (badmask == 0) * (nodup != -9999)
        new_mask = hp.UNSEEN * np.ones(hp.nside2npix(4096))
        new_mask[lss_mask] = 1
        self.mask = new_mask

    def _apply_footprint_mask(self, min_ra, max_ra, min_dec, max_dec):

        hpix = hu.HealPix('ring', 4096)
        ra, dec = hpix.pix2eq(np.arange(hp.nside2npix(
            4096), dtype='int64'))  # masks have nside 4096
        bad_ra_dec = np.where(~((dec > min_dec) & (dec < max_dec)
                                & (ra > min_ra) & (ra < max_ra)))
        self.mask[bad_ra_dec] = hp.UNSEEN
