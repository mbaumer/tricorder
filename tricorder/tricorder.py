"""A module for computing correlation functions on datasets.

I think this will end up handling 2- and 3-point correlations.
"""

from __future__ import division

import subprocess
import time
from sys import stdout

import numpy as np
import treecorr
import yaml

import datasets

output_path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/'


def write_default_config(runname):

    metric = 'Euclidean'
    do3D = True
    scale_angle_w_redshift = False

    config_2pt = {}
    config_3pt = {}
    configdict = {'2PCF': config_2pt, '3PCF': config_3pt, 
                  'do3D': do3D,
                  'scale_angle_w_redshift' : scale_angle_w_redshift}

    config_2pt['min_sep'] = 1
    config_2pt['max_sep'] = 50
    config_2pt['nbins'] = 20
    if not do3D:
        config_2pt['sep_units'] = 'arcmin'
    config_2pt['bin_slop'] = 0
    config_2pt['metric'] = metric

    #config_2pt['min_rpar'] = -20
    #config_2pt['max_rpar'] = 20

    # 3pt params
    config_3pt['min_sep'] = 16
    config_3pt['max_sep'] = 20
    config_3pt['nbins'] = 1
    config_3pt['min_u'] = .25
    config_3pt['max_u'] = .95
    config_3pt['nubins'] = 7
    config_3pt['min_v'] = -1
    config_3pt['max_v'] = 1
    config_3pt['nvbins'] = 20
    config_3pt['bin_slop'] = 0
    if not do3D:
        config_3pt['sep_units'] = 'arcmin'
    config_3pt['metric'] = metric

    #config_3pt['min_rpar'] = -20
    #config_3pt['max_rpar'] = 20

    config_fname = output_path + 'configs/' + runname + '.config'
    f = open(config_fname, 'w')
    f.write(yaml.dump(configdict))
    f.close()


class BaseCorrelation (object):
    """A base class for three-point analyses.

    Mostly overridden by inherited methods.
    """

    def __init__(self):
        raise NotImplementedError(
            'need to specify pixel- or point-based class.')


class PixelCorrelation (BaseCorrelation):

    def __init__(self, dset_fname, config_fname, jk_to_omit=-1):
        try:
            with open(config_fname) as f:
                configdict = yaml.load(f.read())
        except IOError:
            print 'config file ' + config_fname + ' not found.'
            raise

        self.dset_fname = dset_fname
        self.config_fname = config_fname
        self.jk_to_omit = jk_to_omit
        self.dataset = datasets.BaseDataset.fromfilename(dset_fname)
        # drop the abspath and .config
        name = config_fname.split('/')[-1][:-7]
        self.name = name + '_' + str(self.jk_to_omit)
        
        if configdict['scale_angle_w_redshift']:
            avg_redshift = (self.dataset.min_z + self.dataset.max_z)/2
            arcmin_per_mpc = datasets.buzzard_cosmo.arcsec_per_kpc_comoving(avg_redshift).value/60*1000
            
            configdict['2PCF']['min_sep'] = configdict['2PCF']['min_sep']/datasets.buzzard_cosmo.h*arcmin_per_mpc
            configdict['2PCF']['max_sep'] = configdict['2PCF']['max_sep']/datasets.buzzard_cosmo.h*arcmin_per_mpc
            configdict['3PCF']['min_sep'] = configdict['3PCF']['min_sep']/datasets.buzzard_cosmo.h*arcmin_per_mpc
            configdict['3PCF']['max_sep'] = configdict['3PCF']['max_sep']/datasets.buzzard_cosmo.h*arcmin_per_mpc
        
        self.config_2pt = configdict['2PCF']
        self.config_3pt = configdict['3PCF']
        self.cat = None
        self.kk = None
        self.kkk = None

    def reweight_z_dist(self, target_file, nbins=100):
        gal_dset = datasets.RedmagicDataset.fromfilename(target_file)

        dm_z = self.data['REDSHIFT']
        gal_z = gal_dset.data['ZSPEC']

        dm_counts, edges = np.histogram(dm_z, bins=nbins, normed=True, range=(self.min_z, self.max_z));
        gal_counts, _ = np.histogram(gal_z, bins=nbins, normed=True, range=(self.min_z, self.max_z);

        # we need to weight the DM
        dm_idx = np.digitize(dm_z, edges)
        # set weights of overflow and underflow bins to zero
        dm_weights = np.concatenate(([0], gal_counts / dm_counts, [0]))
        return dm_weights[dm_idx]

    def make_treecorr_cat(self):

        weights = None
        str_list = self.dset_fname.split('/')
        if str_list[-3] == 'dark_matter':
            str_list[-3] = 'redmagicHD'
            str_list[-1] = 'ZSPEC'+str(self.min_z)+'_'+str(self.max_z)+'nside1024nJack30.dset'
            target_file = "/".join(str_list)
            weights = reweight_z_dist(target_file)

        if self.jk_to_omit != -1:
            inds_to_keep = np.where(self.dataset.jk_labels != self.jk_to_omit)
        else:
            inds_to_keep = np.arange(len(self.dataset.pixelized[0]))
        ra_to_use = self.dataset.pixelized[0][inds_to_keep]
        dec_to_use = self.dataset.pixelized[1][inds_to_keep]
        counts_to_use = self.dataset.pixelized[2][inds_to_keep]

        # Need to recompute nbar each time
        nbar = np.sum(counts_to_use) / len(counts_to_use)
        print 'nbar for this run is: ', nbar
        kappa_est = (counts_to_use / nbar) - 1
        self.cat = treecorr.Catalog(ra=ra_to_use,
                                    dec=dec_to_use,
                                    ra_units='degrees', dec_units='degrees',
                                    k=kappa_est
                                    w = weights)

    def compute_2pt_pix(self):
        kk = treecorr.KKCorrelation(config=self.config_2pt)
        toc = time.time()
        kk.process(self.cat)
        tic = time.time()
        print '2PCF took', tic - toc
        stdout.flush()
        self.kk = kk

    def compute_3pt_pix(self):
        kkk = treecorr.KKKCorrelation(config=self.config_3pt)
        toc = time.time()
        kkk.process(self.cat)
        tic = time.time()
        print '3PCF took', tic - toc
        stdout.flush()
        self.kkk = kkk

    def write(self):
        dataname = self.dset_fname.split('/')[-1]
        res_dir = 'pix_results/'
        results_prefix = output_path + datasets.mock + \
            '/' + self.dataset.sample_type + '/' + res_dir
        np.save(results_prefix + self.name + '_' +
                dataname + '.zeta', self.kkk.zeta)
        np.save(results_prefix + self.name + '_' +
                dataname + '.weight', self.kkk.weight)
        np.save(results_prefix + self.name +
                '_' + dataname + '.xi', self.kk.xi)

    def run(self):
        self.make_treecorr_cat()
        self.compute_2pt_pix()
        self.compute_3pt_pix()
        self.write()

    def submit(self):
        command_str = "import tricorder; corr = tricorder.PixelCorrelation('" + \
            self.dset_fname + "', '" + self.config_fname + \
            "'," + str(self.jk_to_omit) + "); corr.run()"
        subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=2000]",
                         "python", "-c", command_str])


class PointCorrelation (BaseCorrelation):
    """This class will implement the code I've written before.

    Work in progress re-factoring it.
    """

    def __init__(self, dset_fname, config_fname, set_str, jk_to_omit=-1):
        try:
            with open(config_fname) as f:
                configdict = yaml.load(f.read())
        except IOError:
            print 'config file ' + config_fname + ' not found.'
            raise

        self.dset_fname = dset_fname
        self.config_fname = config_fname
        self.set_str = set_str
        self.jk_to_omit = jk_to_omit
        self.dataset = datasets.BaseDataset.fromfilename(dset_fname)
        # drop the abspath and .config
        name = config_fname.split('/')[-1][:-7]
        self.name = name + '_' + str(self.jk_to_omit)
        
        if configdict['scale_angle_w_redshift']:
            #avg_redshift = (self.dataset.min_z + self.dataset.max_z)/2
            avg_redshift = np.median(self.dataset.data[self.dataset.zvar])
            arcmin_per_mpc = datasets.buzzard_cosmo.arcsec_per_kpc_comoving(avg_redshift).value/60*1000
            
            configdict['2PCF']['min_sep'] = configdict['2PCF']['min_sep']/datasets.buzzard_cosmo.h*arcmin_per_mpc
            configdict['2PCF']['max_sep'] = configdict['2PCF']['max_sep']/datasets.buzzard_cosmo.h*arcmin_per_mpc
            configdict['3PCF']['min_sep'] = configdict['3PCF']['min_sep']/datasets.buzzard_cosmo.h*arcmin_per_mpc
            configdict['3PCF']['max_sep'] = configdict['3PCF']['max_sep']/datasets.buzzard_cosmo.h*arcmin_per_mpc
        
        self.config_2pt = configdict['2PCF']
        self.config_3pt = configdict['3PCF']
        self.do3D = configdict['do3D']
        self.cat = None
        self.random_cat = None
        self.xi = None
        self.nnn = None

    def make_treecorr_cat(self):
        if self.jk_to_omit != -1:
            inds_to_keep = np.where(self.dataset.jk_labels != self.jk_to_omit)
        else:
            inds_to_keep = np.arange(len(self.dataset.data['RA']))
        ra_to_use = self.dataset.data['RA'][inds_to_keep]
        rand_ra_to_use = self.dataset.randoms['RA'][inds_to_keep]
        dec_to_use = self.dataset.data['DEC'][inds_to_keep]
        rand_dec_to_use = self.dataset.randoms['DEC'][inds_to_keep]
        if self.do3D and self.dataset.zvar == 'DISTANCE':
            dist_to_use = self.dataset.data[self.dataset.zvar][inds_to_keep]
            rand_dist_to_use = self.dataset.randoms[self.dataset.zvar][inds_to_keep]
        elif self.do3D:
            dist_to_use = datasets.buzzard_cosmo.h*datasets.buzzard_cosmo.comoving_distance(
                self.dataset.data[self.dataset.zvar]).value[inds_to_keep]
            rand_dist_to_use = datasets.buzzard_cosmo.h*datasets.buzzard_cosmo.comoving_distance(
                self.dataset.randoms[self.dataset.zvar]).value[inds_to_keep]
        else:
            dist_to_use = None
            rand_dist_to_use = None

        self.cat = treecorr.Catalog(ra=ra_to_use,
                                    dec=dec_to_use,
                                    ra_units='degrees', dec_units='degrees',
                                    r=dist_to_use)
        self.random_cat = treecorr.Catalog(ra=rand_ra_to_use,
                                           dec=rand_dec_to_use,
                                           ra_units='degrees', dec_units='degrees',
                                           r=rand_dist_to_use)

    def compute_2pt_raw(self):
        dd = treecorr.NNCorrelation(config=self.config_2pt)
        dr = treecorr.NNCorrelation(config=self.config_2pt)
        rr = treecorr.NNCorrelation(config=self.config_2pt)
        toc = time.time()
        dd.process(self.cat, metric=self.config_2pt['metric'])
        dr.process(self.cat, self.random_cat,
                   metric=self.config_2pt['metric'])
        rr.process(self.random_cat, metric=self.config_2pt['metric'])
        self.xi, varxi = dd.calculateXi(dr=dr, rr=rr)
        tic = time.time()
        print '2PCF took', tic - toc
        stdout.flush()

    def compute_3pt_raw(self):
        nnn = treecorr.NNNCorrelation(config=self.config_3pt)
        toc = time.time()
        setdict = {'d': self.cat, 'r': self.random_cat}
        nnn.process(setdict[self.set_str[0]],
                    setdict[self.set_str[1]], setdict[self.set_str[2]],
                    metric=self.config_3pt['metric'])
        tic = time.time()
        print '3PCF took', tic - toc
        stdout.flush()
        self.nnn = nnn

    def write(self):
        dataname = self.dset_fname.split('/')[-1]
        if self.do3D and self.config_3pt['metric'] == 'Rperp':
            res_dir = 'point_proj_results/'
        elif self.do3D and self.config_3pt['metric'] == 'Euclidean':
            res_dir = 'point_3D_results/'
        else:
            res_dir = 'point_ang_results/'
        results_prefix = output_path + datasets.mock + \
            '/' + self.dataset.sample_type + '/' + res_dir

        if self.set_str == 'ddd':
            self.nnn.write(results_prefix + self.name +
                           '_' + dataname + '.ddd')
            np.save(results_prefix + self.name +
                    '_' + dataname + '.xi', self.xi)
        else:
            np.save(results_prefix + self.name + '_' +
                    dataname + '.' + self.set_str + 'weight', self.nnn.weight)

    def run(self):
        self.make_treecorr_cat()
        if self.set_str == 'ddd':
            self.compute_2pt_raw()
        self.compute_3pt_raw()
        self.write()

    def submit(self):
        command_str = "import tricorder; corr = tricorder.PointCorrelation('" + \
            self.dset_fname + "', '" + self.config_fname + \
            "', '" + self.set_str + "', " + \
            str(self.jk_to_omit) + "); corr.run()"
        subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=2000]",
                         "python", "-c", command_str])
