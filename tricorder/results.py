"""A module for plotting the results of 2 and 3pt analyses.

Will also include the bias inference code.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import lines
import numpy as np
import treecorr
import yaml
from astropy.cosmology import FlatLambdaCDM
from scipy.stats import binned_statistic
import cPickle as pickle

import datasets

buzzard_cosmo = FlatLambdaCDM(68.81, .295)


output_path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/results/'
plot_path = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/plots/'

def add_triangles(alpha=1,lw=1,color='k'):
    # new clear axis overlay with 0-1 limits
    ax2 = plt.axes([0,0,1,1], axisbg=(1,1,1,0))
    ax2.set_axis_off()

    #upper right
    x,y = np.array([[.75, 0.65, .55,.75], [.8, 0.8, 0.83,.8]])
    line = lines.Line2D(x, y, lw=lw, color=color, alpha=alpha)
    ax2.add_line(line)
    #upper middle
    x,y = np.array([[.5, 0.4, 0.44,.5], [.8, 0.8, 0.9,.8]])
    line = lines.Line2D(x, y, lw=lw, color=color, alpha=alpha)
    ax2.add_line(line)
    #upper left
    x,y = np.array([[0.1, 0.2,.3,.1], [.8, 0.8, 0.83,.8]])
    line = lines.Line2D(x, y, lw=lw, color=color, alpha=alpha)
    ax2.add_line(line)
    #lower right
    x,y = np.array([[.75, 0.65, .6,.75], [0.15, 0.15, 0.165,.15]])
    line = lines.Line2D(x, y, lw=lw, color=color, alpha=alpha)
    ax2.add_line(line)
    #lower middle
    x,y = np.array([[.5, 0.4, 0.42,.5], [0.15, 0.15, 0.18,0.15]])
    line = lines.Line2D(x, y, lw=lw, color=color, alpha=alpha)
    ax2.add_line(line)
    #lower left
    x,y = np.array([[0.1, 0.2,.25,.1], [0.15, 0.15, 0.165,0.15]])
    line = lines.Line2D(x, y, lw=lw, color=color, alpha=alpha)
    ax2.add_line(line)

class Results(object):
    """Analyze and plot results of 3pt analyses."""

    def __init__(self, dataname, runname, sample_type, n_jackknife):
        self.sample_type = sample_type #redmagicHD, dark_matter
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
        base_filename = output_path + self.runname + '_' + str(jk_id) \
            + '_' + self.dataname
        zeta_filename = base_filename + '.zeta.npy'
        weight_filename = base_filename + '.weight.npy'
        xi_filename = base_filename + '.xi.npy'
        self.zetas.append(np.load(zeta_filename))
        self.weights.append(np.load(weight_filename))
        self.xis.append(np.load(xi_filename))

    def load_all_data(self):
        for i in xrange(self.n_jackknife):
            try:
                self._load_data_single_run(i)
            except IOError:
                print 'file for jk region '+str(i)+'not found'

    def analyze(self, mode, **kwargs):
        results = []
        for i in xrange(len(self.zetas)):
            bins, binned = self._analyze_single_run(mode, i, **kwargs)
            results.append(binned)
        return results, bins

    def analyze_many_scales(self, **kwargs):
        self.load_all_data()
        for scale in [15, 30, 45, 60, 75]:
            for ratio in [.5, 1.0]:
                for tolerance in [.2]:
                    for nbins in [16]:
                        out, bins = self.analyze(
                            'angle', scale=scale, ratio=ratio,
                            tolerance=tolerance, nbins=nbins)
                        self.make_plots('q', bins, out, scale=scale,
                                        ratio=ratio, tolerance=tolerance, nbins=nbins)
                        self.make_plots('zeta', bins, out, scale=scale,
                                        ratio=ratio, tolerance=tolerance, nbins=nbins)
                        self.make_plots('denom', bins, out, scale=scale,
                                        ratio=ratio, tolerance=tolerance, nbins=nbins)
                        self.make_plots('weights', bins, out, scale=scale,
                                        ratio=ratio, tolerance=tolerance, nbins=nbins,
                                        make_covmat=False)
                        self.make_plots('d1', bins, out, scale=scale,
                                        ratio=ratio, tolerance=tolerance, nbins=nbins,
                                        make_covmat=False)
                        self.make_plots('d2', bins, out, scale=scale,
                                        ratio=ratio, tolerance=tolerance, nbins=nbins,
                                        make_covmat=False)
                        self.make_plots('d3', bins, out, scale=scale,
                                        ratio=ratio, tolerance=tolerance, nbins=nbins,
                                        make_covmat=False)
                        base_name = plot_path + self.dataname + '_' + self.runname + '_' + str(scale) + \
                            '_' + str(ratio) + '_' + str(tolerance) + \
                            '_' + str(nbins)
                        out_tup = (out,bins)
                        with open(base_name+'.res','wb') as f:
                            pickle.dump(out_tup,f)

    def make_plots(self, var, bins, out, scale=30, nbins=16, tolerance=.3,
                   ratio=.5, make_covmat=True):
        out_vec = []
        for res1 in out:
            out_vec.append(res1[var])
        out_vec = np.array(out_vec)
        out_vec = out_vec.T
        plt.plot(bins, out_vec, color='k', alpha=.1)
        # use number of jk analyses that finished.
        jk_factor = (out_vec.shape[1] - 1.0) / out_vec.shape[1]
        plt.errorbar(bins, np.mean(out_vec, axis=1), color='g', marker='o', linestyle='.',
                     linewidth=3, yerr=np.sqrt(jk_factor * np.var(out_vec, axis=1)))
        plt.xlabel('Angle')
        plt.ylabel(var)
        dset = datasets.BaseDataset.fromfilename(
            datasets.output_path + self.sample_type + '/' + self.dataname)
        mean_z = np.mean([dset.max_z, dset.min_z])
        lg_dist = np.round(buzzard_cosmo.kpc_comoving_per_arcmin(
            mean_z).value / 1000 * scale * buzzard_cosmo.h, 1)
        sm_dist = np.round(lg_dist * ratio, 1)
        plt.title(str(dset.min_z) + '<z<' + str(dset.max_z) + '; ' + str(scale * ratio) + "\':" + str(scale) +
                  "\' = " + str(sm_dist) + ':' + str(lg_dist) + ' Mpc/h +/- ' + str(100 * tolerance) + '%')
        base_name = plot_path + self.dataname + '_' + self.runname + '_' + str(scale) + \
            '_' + str(ratio) + '_' + str(tolerance) + \
            '_' + str(nbins) + '_' + var
        plt.savefig(base_name + '.png')
        if make_covmat:
            plt.figure()
            covmat = jk_factor * np.cov(out_vec)
            plt.imshow(covmat, interpolation='None', origin='lower',
                       extent=[0, 180, 0, 180], cmap='RdBu_r')
            plt.title(self.dataname)
            cbar = plt.colorbar()
            cbar.set_label(var)
            plt.savefig(base_name + '_covmat.png')

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

        def yfit(x):
            return np.exp(poly(np.log(x)))

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
        isCorrectSize = ((np.exp(self.kkk.logr) > scale-scale*tolerance) &
                         (np.exp(self.kkk.logr) < scale+scale*tolerance))
        res, b, _ = binned_statistic(
            np.abs(self.kkk.v[isCorrectRatio & isCorrectSize]), var[isCorrectRatio & isCorrectSize],
            bins=nbins, statistic=stat)
        b += (b[1] - b[0]) / 2
        b = b[:-1]
        return res, b

    def _computeXvsAngle(self, var, stat='mean', scale=6, ratio=.5,
                         tolerance=.1, nbins=15, **kwargs):

        # angle at which elongate becomes collapsed
        # changed fixed 0.25 to ratio/2 so ratio can vary
        transition_angle = np.arccos(ratio/2) / np.pi * 180
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
