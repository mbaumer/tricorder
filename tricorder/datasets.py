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

from make_data_randoms import (generate_randoms_radec, index_to_radec,
                               radec_to_index)

buzzard_cosmo = FlatLambdaCDM(68.81, .295)

zvar_labels = {'ZSPEC': r'$z_{true}$',
               'ZREDMAGIC': r'z_{RM}',
               'DISTANCE': r' Comoving Distance (Mpc/$h$)',
               'REDSHIFT': r'$z_{true}$',
               }

mock = 'Buzzard_v1.6_Y3_0_a'

output_path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/' + mock + '/'


class BaseDataset (object):

    def __init__(self):

        if not hasattr(self, 'sample_type'):
            self.sample_type = 'unk'

        self.output_path = output_path + self.sample_type + '/'

        self.data = None
        self.randoms = None
        self.mask = None

        self.n_jackknife = 30
        self.jk_labels = None

        self.min_z = None
        self.max_z = None

        self.nside = None
        self.nbar = None
        self.pixelized = None  # tuple: (pixel ra, pixel dec, counts)

    @classmethod
    def fromfilename(cls, filename):
        """Initialize a BaseDataset from a pickle written by self.write."""
        data = pickle.load(open(filename, 'rb'))
        data.data = np.load(filename + '_data.npy')
        data.randoms = np.load(filename + '_randoms.npy')
        return data

    def load_data(self):
        raise NotImplementedError(
            'Need to specify sample class to know how to load!')

    def pixelize_at_target_nside(self, nside):
        # this also effectively applies the mask to the data
        self.nside = nside
        # so averages are computed in downgrade
        self.mask[self.mask == hp.UNSEEN] = 0  # step 1
        mask_targetnside = hp.pixelfunc.ud_grade(
            self.mask, pess=False, nside_out=self.nside)
        gal_index_targetnside = radec_to_index(
            self.data['DEC'], self.data['RA'], self.nside)
        mask_targetnside[mask_targetnside == 0] = hp.UNSEEN
        self.mask[self.mask == 0] = hp.UNSEEN  # undo step 1

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
        dec, ra = index_to_radec(pixels_to_count, self.nside)
        final_counts = full_map_counts[pixels_to_count]
        self.pixelized = (ra[0], dec[0], final_counts)

    def _apply_footprint_data(self, min_ra, max_ra, min_dec, max_dec):
        self.data = self.data[self.data['RA'] > min_ra]
        self.data = self.data[self.data['RA'] < max_ra]
        self.data = self.data[self.data['DEC'] > min_dec]
        self.data = self.data[self.data['DEC'] < max_dec]

    def _apply_footprint_mask(self, min_ra, max_ra, min_dec, max_dec):

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

    def generate_randoms(self, oversamp,
                         Ngen=1000000, Ntries_max=10000):

        Ncurrent = 0
        Ntry = 0
        Ntot = len(self.data) * oversamp

        minra = min(self.data['RA'])
        maxra = max(self.data['RA'])
        mindec = min(self.data['DEC'])
        maxdec = max(self.data['DEC'])

        zdist = self.data[self.zvar]

        random_ra = []
        random_dec = []
        random_z = []
        while ((Ncurrent < Ntot) & (Ntry < Ntries_max)):
            if Ntry % 100 == 1:
                print(Ntry, Ntries_max, Ncurrent, Ntot, Ngen)
            # generate random ra and dec
            ra_i, dec_i = generate_randoms_radec(minra, maxra,
                                                 mindec, maxdec, Ngen)
            indices_i = radec_to_index(dec_i, ra_i, 4096)

            good_index = np.in1d(indices_i, np.where(self.mask != hp.UNSEEN))
            ra_i = ra_i[good_index]
            dec_i = dec_i[good_index]
            indices_i = indices_i[good_index]

            z_i = np.random.choice(zdist, len(ra_i))

            random_ra += list(ra_i)
            random_dec += list(dec_i)
            random_z += list(z_i)

            Ncurrent = len(random_ra)
            Ntry += 1
        if Ntry >= Ntries_max:
            print('Warning! We gave up after {0} tries, finding {1} objects instead of the desired {2} objects!'.format(
                Ntry, Ncurrent, Ntot))

        randoms = np.zeros(len(random_ra), dtype=[
                           ('RA', '>f4'), ('DEC', '>f4'), (self.zvar, '>f4')])
        randoms['RA'] = random_ra
        randoms['DEC'] = random_dec
        randoms[self.zvar] = random_z

        if len(randoms) > Ntot:
            inds_to_keep = np.random.choice(np.arange(len(randoms)), size=Ntot,
                                            replace=False)
            randoms = randoms[inds_to_keep]

        self.randoms = randoms

    def make(self, nside, oversamp, min_z, max_z, min_ra,
             max_ra, min_dec, max_dec, write=True):
        self.load_data()
        self.apply_z_cut(min_z, max_z)
        self.apply_footprint(min_ra, max_ra, min_dec, max_dec)
        self.pixelize_at_target_nside(nside)
        self.generate_randoms(oversamp)
        self.compute_new_jk_regions()
        if write:
            self.write()

    def write(self):
        """Save a BaseDataset instance as a pickle.

        These will be read in later by the 3PCF analysis.
        """
        # don't actually pickle out this huge stuff
        del self.mask

        name = self.output_path + 'data/' + str(self.zvar) + \
            str(self.min_z) + '_' + str(self.max_z) + \
            'nside' + str(self.nside) + 'nJack' \
            + str(self.n_jackknife) + '.dset'

        if self.data is not None:
            np.save(name + '_data.npy', self.data)
        del self.data

        if self.randoms is not None:
            np.save(name + '_randoms.npy', self.randoms)
        del self.randoms

        with open(name, 'wb') as pickle_file:
            pickle.dump(self, pickle_file, protocol=2)


class RedmagicDataset(BaseDataset):
    def __init__(self, use_spec_z=True):
        self.sample_type = 'redmagicHD'
        if use_spec_z:
            self.zvar = 'ZSPEC'
        else:
            self.zvar = 'ZREDMAGIC'
        self.datapath = '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-0/a/buzzard-0_1.6_y3_run_redmapper_v6.4.20_redmagic_highdens_0.5-10.fit'
        self.maskpath = '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-0/a/buzzard-0_1.6_y3_run_redmapper_v6.4.20_redmagic_highdens_0.5_vlim_zmask.fit'

        super(RedmagicDataset, self).__init__()

    def load_data(self):
        self.data = fits.getdata(self.datapath)
        self.mask = hp.read_map(self.maskpath, 0, partial=True)
        zmask = hp.read_map(self.maskpath, 1, partial=True)
        self.mask[self.mask < .95] = hp.UNSEEN
        self.mask[zmask < .6] = hp.UNSEEN

class MICERedmagicDataset(BaseDataset):
    def __init__(self, use_spec_z=True):
        self.sample_type = 'redmagicHD'
        if use_spec_z:
            self.zvar = 'ZSPEC'
        else:
            self.zvar = 'ZREDMAGIC'
        self.datapath = '/nfs/slac/g/ki/ki23/des/jderose/des/MICE/redmagic/mice2_des_run_redmapper_v6.4.16_redmagic_highdens_0.5-10.fit'
        self.maskpath = None

        super(MICERedmagicDataset, self).__init__()

    def load_data(self):
        self.data = fits.getdata(self.datapath)
        self.data['DEC'] = -self.data['DEC']
        self.mask = np.ones(hp.nside2npix(4096))
            
class DMDataset(BaseDataset):
    def __init__(self, use_spec_z=True):
        self.sample_type = 'dark_matter'
        self.zvar = 'REDSHIFT'
        self.datapath = '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-0/a/downsampled_particles.fits.downsample'
        self.maskpath = '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/a/buzzard_1.6-6a_run_redmapper_v6.4.18_redmagic_highdens_0.5_vlim_zmask.fit'
        super(DMDataset, self).__init__()

    def load_data(self):
        data = fits.getdata(self.datapath)
        c1 = fits.Column(name='RA', array=data['azim_ang'], format='E')
        c2 = fits.Column(name='DEC', array=data['polar_ang'], format='E')
        c3 = fits.Column(name='REDSHIFT', array=data['redshift'], format='E')
        t = fits.BinTableHDU.from_columns([c1, c2, c3])
        self.data = t.data

        self.mask = hp.read_map(self.maskpath, 0, partial=True)
        zmask = hp.read_map(self.maskpath, 1, partial=True)
        self.mask[self.mask < .95] = hp.UNSEEN
        self.mask[zmask < .6] = hp.UNSEEN

class MICEDMDataset(BaseDataset):
    def __init__(self, use_spec_z=True):
        self.sample_type = 'dark_matter'
        self.zvar = 'REDSHIFT'
        self.datapath = '/nfs/slac/g/ki/ki23/des/jderose/des/MICE/dm/dm_v0.fits'
        self.maskpath = None
        super(MICEDMDataset, self).__init__()

    def load_data(self):
        data = fits.getdata(self.datapath)
        data['dec'] = -data['dec']
        c1 = fits.Column(name='RA', array=data['ra'], format='E')
        c2 = fits.Column(name='DEC', array=data['dec'], format='E')
        c3 = fits.Column(name='REDSHIFT', array=data['redshift'], format='E')
        t = fits.BinTableHDU.from_columns([c1, c2, c3])
        self.data = t.data
        
        self.mask = np.ones(hp.nside2npix(4096))

class LSSDataset(BaseDataset):
    def __init__(self, use_spec_z=True):
        self.sample_type = 'lss_sample'
        self.zvar = 'REDSHIFT'
        self.datapath = '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y3/buzzard/flock/buzzard-0/a/Buzzard_v1.6_Y3_gold.fits'
        self.maskpath = ['/nfs/slac/g/ki/ki23/des/jderose/SkyFactory/chinchilla-herd/Chinchilla-1/sampleselection/y1a1_gold_1.0.2_wide_footprint_4096.fits.gz',
                         '/nfs/slac/g/ki/ki23/des/jderose/SkyFactory/chinchilla-herd/Chinchilla-1/sampleselection/y1a1_gold_1.0.2_wide_badmask_4096.fits.gz'
                         ]
        super(LSSDataset, self).__init__()

    def load_data(self):
        gold = fits.getdata(self.datapath)
        lss = gold[gold['lss-sample'] == 1]
        c1 = fits.Column(name='RA', array=lss['ra'], format='E')
        c2 = fits.Column(name='DEC', array=lss['dec'], format='E')
        c3 = fits.Column(name='REDSHIFT', array=lss['redshift'], format='E')
        t = fits.BinTableHDU.from_columns([c1, c2, c3])
        self.data = t.data

        footprint = hp.read_map(self.maskpath[0])
        badmask = hp.read_map(self.maskpath[1])
        lss_mask = (footprint >= 1) * (badmask == 0)
        new_mask = hp.UNSEEN * np.ones(hp.nside2npix(4096))
        new_mask[lss_mask] = 1
        self.mask = new_mask
