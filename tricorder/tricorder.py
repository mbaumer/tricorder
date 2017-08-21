"""A module for computing correlation functions on datasets.

I think this will end up handling 2- and 3-point correlations.
"""

from __future__ import division

import time
from sys import stdout

import numpy as np
import treecorr
import yaml

import datasets

output_path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/results/'


def write_default_config(runname):
    config_2pt = {}
    config_3pt = {}
    configdict = {'2PCF': config_2pt, '3PCF': config_3pt}

    config_2pt['min_sep'] = 5
    config_2pt['max_sep'] = 120
    config_2pt['sep_units'] = 'arcmin'
    config_3pt['bin_slop'] = 0.1

    # 3pt params
    config_3pt['min_sep'] = 10
    config_3pt['max_sep'] = 60
    config_3pt['nbins'] = 50
    config_3pt['min_u'] = 0
    config_3pt['max_u'] = 1
    config_3pt['nubins'] = 100
    config_3pt['min_v'] = -1
    config_3pt['max_v'] = 1
    config_3pt['nvbins'] = 100
    config_3pt['bin_slop'] = 0.1

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

    def __init__(self, dset_fname, config_fname):
        try:
            with open(config_fname) as f:
                configdict = yaml.load(f.read())
        except IOError:
            print 'config file ' + config_fname + ' not found.'
            raise

        self.dataset = datasets.BaseDataset.fromfilename(dset_fname)
        self.name = config_fname[:-7]  # drop the .config
        self.config_2pt = configdict['2PCF']
        self.config_3pt = configdict['3PCF']
        self.cat = None
        self.kk = None
        self.kkk = None

    def make_treecorr_cat(self):
        kappa_est = (self.dataset.pixelized[2] / self.dataset.nbar) - 1
        self.cat = treecorr.Catalog(ra=self.dataset.pixelized[0],
                                    dec=self.dataset.pixelized[1],
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
        np.save(output_path + self.name + '.zeta', self.kkk.zeta)
        np.save(output_path + self.name + '.weight', self.kkk.weight)
        np.save(output_path + self.name + '.xi', self.kk.xi)

    def run(self):
        self.make_treecorr_cat()
        self.compute_2pt_pix()
        self.compute_3pt_pix()
        self.write()

    def submit(self):
        print (["bsub", "-W", "47:00", "-R", "rusage[mem=8000]",
                "python", "-c", "import tricorder; \
                corr = tricorder.PixelCorrelation('" + self.dset_fname
                + "'," + self.config_fname + "); corr.run()"])
        """
        subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=8000]",
                         "python", "-c", "import tricorder; \
                        corr = tricorder.PixelCorrelation('" + self.dset_fname
                         + "'," + self.config_fname + "); corr.run()"])
        """


class PointCorrelation (BaseCorrelation):
    """This class will implement the code I've written before.

    Work in progress re-factoring it.
    """

    def __init__(self):
        raise NotImplementedError('will put this in soon from old code.')
