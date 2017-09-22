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

output_path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/results/'


def write_default_config(runname):

    metric = 'Euclidean'
    do3D = False

    config_2pt = {}
    config_3pt = {}
    configdict = {'2PCF': config_2pt, '3PCF': config_3pt, 'do3D': do3D}

    config_2pt['min_sep'] = 5
    config_2pt['max_sep'] = 135
    config_2pt['nbins'] = 20
    if not do3D:
        config_2pt['sep_units'] = 'arcmin'
    config_2pt['bin_slop'] = 0
    config_2pt['metric'] = metric

    # 3pt params
    config_3pt['min_sep'] = 68
    config_3pt['max_sep'] = 112
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

    config_fname = output_path + runname + '.config'
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
        self.config_2pt = configdict['2PCF']
        self.config_3pt = configdict['3PCF']
        self.cat = None
        self.kk = None
        self.kkk = None

    def make_treecorr_cat(self):
        inds_to_keep = np.where(self.dataset.jk_labels != self.jk_to_omit)
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
                                    k=kappa_est)

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
        np.save(output_path + self.name + '_' +
                dataname + '.zeta', self.kkk.zeta)
        np.save(output_path + self.name + '_' +
                dataname + '.weight', self.kkk.weight)
        np.save(output_path + self.name + '_' + dataname + '.xi', self.kk.xi)

    def run(self):
        self.make_treecorr_cat()
        self.compute_2pt_pix()
        self.compute_3pt_pix()
        self.write()

    def submit(self):
        command_str = "import tricorder; corr = tricorder.PixelCorrelation('" + \
            self.dset_fname + "', '" + self.config_fname + \
            "'," + str(self.jk_to_omit) + "); corr.run()"
        subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=8000]",
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
        self.config_2pt = configdict['2PCF']
        self.config_3pt = configdict['3PCF']
        self.do3D = configdict['do3D']
        self.cat = None
        self.random_cat = None
        self.xi = None
        self.nnn = None

    def make_treecorr_cat(self):
        inds_to_keep = np.where(self.dataset.jk_labels != self.jk_to_omit)
        ra_to_use = self.dataset.data['RA'][inds_to_keep]
        rand_ra_to_use = self.dataset.randoms['RA'][inds_to_keep]
        dec_to_use = self.dataset.data['DEC'][inds_to_keep]
        rand_dec_to_use = self.dataset.randoms['DEC'][inds_to_keep]
        if self.do3D and self.dataset.zvar == 'DISTANCE':
            dist_to_use = self.dataset.data[self.dataset.zvar][inds_to_keep]
            rand_dist_to_use = self.dataset.randoms[self.dataset.zvar][inds_to_keep]
        elif self.do3D:
            dist_to_use = datasets.buzzard_cosmo.comoving_distance(
                self.dataset.data[self.dataset.zvar])[inds_to_keep]
            rand_dist_to_use = datasets.buzzard_cosmo.comoving_distance(
                self.dataset.randoms[self.dataset.zvar])[inds_to_keep]
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
        dd.process(self.cat)
        dr.process(self.cat, self.random_cat)
        rr.process(self.random_cat)
        self.xi, varxi = dd.calculateXi(dr=dr, rr=rr)
        tic = time.time()
        print '2PCF took', tic - toc
        stdout.flush()

    def compute_3pt_raw(self):
        nnn = treecorr.NNNCorrelation(config=self.config_3pt)
        toc = time.time()
        setdict = {'d': self.cat, 'r': self.random_cat}
        nnn.process(setdict[self.set_str[0]],
                    setdict[self.set_str[1]], setdict[self.set_str[2]])
        tic = time.time()
        print '3PCF took', tic - toc
        stdout.flush()
        self.nnn = nnn

    def write(self):
        dataname = self.dset_fname.split('/')[-1]
        np.save(output_path + self.name + '_' +
                dataname + '.' + self.set_str + 'weight', self.nnn.weight)
        if self.set_str == 'ddd':
            self.nnn.write(output_path + self.name + '_' + dataname + '.ddd')
            np.save(output_path + self.name + '_' + dataname + '.xi', self.xi)

    def run(self):
        self.make_treecorr_cat()
        if self.set_str == 'ddd':
            self.compute_2pt_raw()
        self.compute_3pt_raw()
        self.write()

    def submit(self):
        command_str = "import tricorder; corr = tricorder.PointCorrelation('" + \
            self.dset_fname + "', '" + self.config_fname + \
            "'," + str(self.jk_to_omit) + "); corr.run()"
        subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=8000]",
                         "python", "-c", command_str])
