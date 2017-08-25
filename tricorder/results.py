"""A module for plotting the results of 2 and 3pt analyses.

Will also include the bias inference code.
"""

import numpy as np
import treecorr
import yaml
from scipy.stats import binned_statistic

output_path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/results/'


class Results(object):
    """Analyze and plot results of 3pt analyses."""

    def __init__(self, dataname, runname, n_jackknife):
        self.dataname = dataname
        self.runname = runname
        config_fname = output_path + self.runname + '.config'
        try:
            with open(config_fname) as f:
                self.configdict = yaml.load(f.read())
        except IOError:
            print 'config file ' + config_fname + ' not found.'
            raise

        self.kk = treecorr.KKCorrelation(config=self.configdict['2PCF'])
        self.kkk = treecorr.KKKCorrelation(config=self.configdict['3PCF'])

        self.n_jackknife = n_jackknife

        self.zetas = []
        self.weights = []
        self.xis = []

    def _load_data_single_run(self, jk_id):
        """Load data from single jk run."""
        base_filename = output_path + self.dataname + \
            self.runname + '_' + str(jk_id)
        zeta_filename = base_filename + '.zeta.npy'
        weight_filename = base_filename + '.weight.npy'
        xi_filename = base_filename + '.xi.npy'
        self.zetas.append(np.load(zeta_filename))
        self.weights.append(np.load(weight_filename))
        self.xis.append(np.load(xi_filename))

    def load_all_data(self):
        for i in xrange(self.n_jackknife):
            self._load_data_single_run(i)

    def analyze(self, mode, **kwargs):
        results = []
        for i in xrange(self.n_jackknife):
            bins, binned = self._analyze_single_run(mode, i, **kwargs)
            results.append(binned)
        return results, bins

    def _analyze_single_run(self, mode, jk_id, **kwargs):

        if mode == 'angle':
            get_binned_stat = self._computeXvsAngle
        if mode == 'equi':
            get_binned_stat = self._compute_x_vs_side_length
        if mode == 'v':
            get_binned_stat = self._computeXvsV

        binned = {}

        binned['d1'], bins = get_binned_stat(
            self.kkk.u * np.abs(self.kkk.v) * np.exp(self.kkk.logr)
            + np.exp(self.kkk.logr), **kwargs)
        binned['d2'], bins = get_binned_stat(np.exp(self.kkk.logr), **kwargs)
        binned['d3'], bins = get_binned_stat(
            self.kkk.u * np.exp(self.kkk.logr), **kwargs)

        unweighted_zetas, bins = get_binned_stat(
            self.zetas[jk_id] * self.weights[jk_id], stat='sum', **kwargs)
        binned['zeta'] = unweighted_zetas / \
            get_binned_stat(self.weights[jk_id], stat='sum', **kwargs)[0]
        binned['weights'], bins = get_binned_stat(
            self.weights[jk_id], stat='sum', **kwargs)
        binned['denom'] = self._get_two_point_expectation(
            binned['d1'], binned['d2'], binned['d3'], jk_id)
        binned['q'] = binned['zeta'] / binned['denom']
        return bins, binned

    def _get_two_point_expectation(self, d1bins, d2bins, d3bins, jk_id):

        xi = self.xis[jk_id]

        coeffs = np.polyfit(self.kk.logr, np.log(xi), deg=1)
        poly = np.poly1d(coeffs)

        def yfit(x): return np.exp(poly(np.log(x)))

        xi1 = yfit(d1bins)
        xi2 = yfit(d2bins)
        xi3 = yfit(d3bins)
        denom_bins = (xi1 * xi2 + xi2 * xi3 + xi3 * xi1)
        return denom_bins

    def _computeAngularBins(self, collapsed=False):
        # if v < 0: collapsed = not collapsed
        self.kkk.v = np.abs(self.kkk.v)
        d2 = np.exp(self.kkk.logr)
        d3 = self.kkk.u * np.exp(self.kkk.logr)
        d1 = self.kkk.v * d3 + d2
        # law of cosines
        if not collapsed:
            cosine = (d2**2 + d3**2 - d1**2) / (2 * d2 * d3 + 1e-9)
        else:
            cosine = (d1**2 + d3**2 - d2**2) / (2 * d1 * d3 + 1e-9)
        bins = np.arccos(cosine) / np.pi * 180
        return bins

    def _compute_x_vs_side_length(self, var, stat='mean', nbins=15,
                                  tolerance=.1, **kwargs):
        isEquilateral = (self.kkk.u > 1 -
                         tolerance) & (np.abs(self.kkk.v) < tolerance)
        res, b, _ = binned_statistic(
            self.kkk.logr[isEquilateral], var[isEquilateral],
            bins=nbins, statistic=stat)
        b += (b[1] - b[0]) / 2
        b = b[:-1]
        return res, b

    def _computeXvsV(self, var, stat='mean', nbins=15, scale=6, ratio=.5,
                     tolerance=.1, **kwargs):
        isCorrectRatio = ((self.kkk.u > (ratio - tolerance * ratio)) &
                          (self.kkk.u < (ratio + tolerance * ratio)))
        res, b, _ = binned_statistic(
            np.abs(self.kkk.v[isCorrectRatio]), var[isCorrectRatio],
            bins=nbins, statistic=stat)
        b += (b[1] - b[0]) / 2
        b = b[:-1]
        return res, b

    def _computeXvsAngle(self, var, stat='mean', scale=6, ratio=.5,
                         tolerance=.1, nbins=15, **kwargs):

        # angle at which elongate becomes collapsed
        # TODO maybe different for ratios other than 0.5?
        transition_angle = np.arccos(.25) / np.pi * 180
        N_low_bins = np.floor(transition_angle / 180 * nbins)
        coll_bins = np.linspace(0, transition_angle, num=N_low_bins)
        elong_bins = np.linspace(transition_angle, 180, num=nbins - N_low_bins)

        collapsed_angles = self._computeAngularBins(collapsed=True)
        elongated_angles = self._computeAngularBins(collapsed=False)
        isRightSize = (np.exp(self.kkk.logr) * self.kkk.u > scale * ratio - scale * ratio * tolerance) & (
            np.exp(self.kkk.logr) * self.kkk.u < scale * ratio + scale * ratio * tolerance)
        isCollapsed = (((self.kkk.u * np.abs(self.kkk.v)) * np.exp(self.kkk.logr) + np.exp(self.kkk.logr) > scale - scale * tolerance)
                       & ((self.kkk.u * np.abs(self.kkk.v)) * np.exp(self.kkk.logr) + np.exp(self.kkk.logr) < scale + scale * tolerance))
        isElongated = ((np.exp(self.kkk.logr) > scale - scale * tolerance)
                       & (np.exp(self.kkk.logr) < scale + scale * tolerance))
        out1, bins1, _ = binned_statistic(elongated_angles[np.where(isRightSize & isElongated)], var[np.where(
            isRightSize & isElongated)], bins=elong_bins, statistic=stat)
        out2, bins2, _ = binned_statistic(collapsed_angles[np.where(isRightSize & isCollapsed)], var[np.where(
            isRightSize & isCollapsed)], bins=coll_bins, statistic=stat)
        full_var = np.concatenate((out2, out1))
        bins1 += (bins1[1] - bins1[0]) / 2  # make edges centers
        bins2 += (bins2[1] - bins2[0]) / 2
        full_bins = np.concatenate((bins2[:-1], bins1[:-1]))
        return full_var, full_bins
