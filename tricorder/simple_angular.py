#!/u/ki/mbaumer/anaconda/bin/python
from __future__ import division
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import treecorr
import yaml
import numpy as np
import datasets
import time
import tricorder
import yaml
import sys
import os
import paths
from glob import glob
import pandas as pd
import pickle

from make_data_randoms import (generate_randoms_radec, index_to_radec,
                               radec_to_index)


def load_config(config_fname):
    config_path = os.path.join(
        paths.config_dir, config_fname+'.config')
    with open(config_path) as f:
        return yaml.load(f.read())


def get_zslice(data, min_z, max_z, zvar):
    return data[(data[zvar] > min_z) & (data[zvar] < max_z)]

def downselect_pz(input_data, target_cts, target_bins, input_zvar, oversamp):
    
    input_labels = np.digitize(input_data[input_zvar], target_bins)
    input_cts, _ = np.histogram(
        input_data[input_zvar], bins=target_bins, range=(0, 1))
    input_weights = target_cts/(input_cts+1e-40)
    return input_data[np.random.rand(len(input_data)) < oversamp*input_weights[input_labels-1]]

def downselect(input_data, target, input_zvar, target_zvar, oversamp):
    target_cts, bins = np.histogram(
        target[target_zvar], range=(0, 1), bins=100)
    input_labels = np.digitize(input_data[input_zvar], bins)
    input_cts, _ = np.histogram(
        input_data[input_zvar], bins=bins, range=(0, 1))
    input_weights = target_cts/(input_cts+1e-40)
    return input_data[np.random.rand(len(input_data)) < oversamp*input_weights[input_labels-1]]


def calc_2pt(data, randoms, config_fname, do3D, ra_var='RA', dec_var='DEC',
             random_ra_var='RA', random_dec_var='DEC', data_zvar=None, random_zvar=None,):
    if do3D:
        assert data_zvar is not None
        assert random_zvar is not None
    config_2pt = load_config(config_fname)['2PCF']
    if not do3D:
        cat = treecorr.Catalog(ra=data[ra_var], dec=data[dec_var],
                               dec_units='degrees', ra_units='degrees',
                               )
        random_cat = treecorr.Catalog(ra=randoms[random_ra_var], dec=randoms[random_dec_var],
                                      dec_units='degrees', ra_units='degrees',
                                      )
    else:
        cat = treecorr.Catalog(ra=data[ra_var], dec=data[dec_var],
                               dec_units='degrees', ra_units='degrees',
                               r=datasets.buzzard_cosmo.comoving_distance(data[data_zvar]).value)
        random_cat = treecorr.Catalog(ra=randoms[random_ra_var], dec=randoms[random_dec_var],
                                      dec_units='degrees', ra_units='degrees',
                                      r=datasets.buzzard_cosmo.comoving_distance(randoms[random_zvar]).value)

    print config_2pt
    dd = treecorr.NNCorrelation(config=config_2pt)
    dr = treecorr.NNCorrelation(config=config_2pt)
    rr = treecorr.NNCorrelation(config=config_2pt)
    dd.process(cat, metric=config_2pt['metric'])
    dr.process(cat, random_cat,
               metric=config_2pt['metric'])
    rr.process(random_cat, metric=config_2pt['metric'])
    xi, varxi = dd.calculateXi(dr=dr, rr=rr)
    return xi


def calc_3pt(data, randoms, config_fname, do3D, ra_var='RA',
             dec_var='DEC', random_ra_var='RA', random_dec_var='DEC',
             data_zvar=None, random_zvar=None, outvar='zeta'):

    if do3D:
        assert data_zvar is not None
        assert random_zvar is not None

    config_3pt = load_config(config_fname)['3PCF']
    if not do3D:
        cat = treecorr.Catalog(ra=data[ra_var], dec=data[dec_var],
                               dec_units='degrees', ra_units='degrees',
                               )
        random_cat = treecorr.Catalog(ra=randoms[random_ra_var], dec=randoms[random_dec_var],
                                      dec_units='degrees', ra_units='degrees',
                                      )
    else:
        cat = treecorr.Catalog(ra=data[ra_var], dec=data[dec_var],
                               dec_units='degrees', ra_units='degrees',
                               r=datasets.buzzard_cosmo.comoving_distance(data[data_zvar]).value)
        random_cat = treecorr.Catalog(ra=randoms[random_ra_var], dec=randoms[random_dec_var],
                                      dec_units='degrees', ra_units='degrees',
                                      r=datasets.buzzard_cosmo.comoving_distance(randoms[random_zvar]).value)

    print config_3pt

    if outvar == 'zeta':
        ddd = treecorr.NNNCorrelation(config=config_3pt)
        ddr = treecorr.NNNCorrelation(config=config_3pt)
        drd = treecorr.NNNCorrelation(config=config_3pt)
        rdd = treecorr.NNNCorrelation(config=config_3pt)
        rdr = treecorr.NNNCorrelation(config=config_3pt)
        rrd = treecorr.NNNCorrelation(config=config_3pt)
        drr = treecorr.NNNCorrelation(config=config_3pt)
        rrr = treecorr.NNNCorrelation(config=config_3pt)
        ddd.process(cat, metric=config_3pt['metric'])
        ddr.process(cat, cat,  random_cat, metric=config_3pt['metric'])
        drd.process(cat, random_cat, cat, metric=config_3pt['metric'])
        rdd.process(random_cat, cat, cat, metric=config_3pt['metric'])
        rdr.process(random_cat, cat,  random_cat, metric=config_3pt['metric'])
        rrd.process(random_cat, random_cat,  cat, metric=config_3pt['metric'])
        drr.process(cat, random_cat,  random_cat, metric=config_3pt['metric'])
        rrr.process(random_cat, random_cat,  random_cat,
                    metric=config_3pt['metric'])
        output, varzeta = ddd.calculateZeta(
            ddr=ddr, drd=drd, rdd=rdd, rrd=rrd, rdr=rdr, drr=drr, rrr=rrr)
    else:
        nnn = treecorr.NNNCorrelation(config=config_3pt)
        toc = time.time()
        setdict = {'d': cat, 'r': random_cat}
        nnn.process(setdict[outvar[0]],
                    setdict[outvar[1]], setdict[outvar[2]],
                    metric=config_3pt['metric'])
        tic = time.time()
        print '3PCF took', tic - toc
        output = nnn.ntri
    return output


def calc_3pt_noisy_photoz_lss(dset_id, config_fname, do3D, min_z, max_z, sigma_z, zvar, random_zvar, random_oversamp):
    randoms = fits.getdata(paths.lss_y1_randoms)
    data = fits.getdata(paths.lss_y1[dset_id])
    data = data[data['lss-sample'] == 1]

    ra_var = 'RA'
    dec_var = 'DEC'

    data = data[data[ra_var] < 100]
    data = data[data[dec_var] < -20]

    data[zvar] += np.random.normal(size=len(data), scale=sigma_z)
    data_slice = get_zslice(data, min_z, max_z, zvar)
    randoms_slice = get_zslice(randoms, min_z, max_z, random_zvar)

    randoms_slice = randoms_slice[np.random.rand(len(randoms_slice)) < (
        len(data_slice)/len(randoms_slice)*random_oversamp)]

    xi = calc_2pt(data_slice, randoms_slice, config_fname, do3D,
                  ra_var=ra_var, dec_var=dec_var,
                  data_zvar=zvar, random_zvar=random_zvar,)
    zeta = calc_3pt(data_slice, randoms_slice, config_fname, do3D,
                    ra_var=ra_var, dec_var=dec_var,
                    data_zvar=zvar, random_zvar=random_zvar,)

    xi_file_name = config_fname + \
        '_lssdset'+str(dset_id)+'_sigma'+str(sigma_z) + \
        '_'+str(zvar)+'_'+str(min_z)+'_'+str(max_z) + \
        '_rsamp'+str(random_oversamp)+'.xi'
    zeta_file_name = config_fname+'_lssdset' + \
        str(dset_id)+'_sigma'+str(sigma_z) + \
        '_'+str(zvar)+'_'+str(min_z)+'_'+str(max_z) + \
        '_rsamp'+str(random_oversamp)+'.zeta'

    if not do3D:
        np.save(os.path.join(paths.ang_out_dir, xi_file_name), xi)
        np.save(os.path.join(paths.ang_out_dir, zeta_file_name), zeta)
    else:
        np.save(os.path.join(paths.corr_out_dir, xi_file_name), xi)
        np.save(os.path.join(paths.corr_out_dir, zeta_file_name), zeta)


def calc_3pt_noisy_photoz_mice(dset_id, jk_id, config_fname, do3D, min_z, max_z, sigma_z, zvar, random_zvar, random_oversamp, outvar='zeta'):

    if min_z == .6:
        data = fits.getdata(paths.rm_mice_y1_HL[dset_id])
        randoms = generate_randoms(data, random_oversamp, zvar)
    elif min_z == .75:
        data = fits.getdata(paths.rm_mice_y1_HHL[dset_id])
        randoms = generate_randoms(data, random_oversamp, zvar)
    else:
        data = fits.getdata(paths.rm_mice_y1[dset_id])
        randoms = generate_randoms(data, random_oversamp, zvar)

    ra_var = 'RA'
    dec_var = 'DEC'

    data[dec_var] = -data[dec_var]
    data = data[data[ra_var] > 0]
    data = data[data[ra_var] < 90]
    data = data[data[dec_var] > -60]
    data = data[data[dec_var] < -40]

    data[zvar] += np.random.normal(size=len(data), scale=sigma_z)
    data_slice = get_zslice(data, min_z, max_z, zvar)
    randoms_slice = get_zslice(randoms, min_z, max_z, random_zvar)

    #remove jk region
    jk_classifier = pickle.load( open( "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/jk_classifiers/rectangle_0_90_-60_-40_jk.pkl", "rb" ) )
    data_inds = jk_classifier.predict(zip(data_slice[ra_var],data_slice[dec_var]))
    random_inds = jk_classifier.predict(zip(randoms_slice[ra_var],randoms_slice[dec_var]))
    data_slice = data_slice[data_inds != jk_id]
    randoms_slice = randoms_slice[random_inds != jk_id]

    if (outvar == 'zeta') | (outvar == 'ddd'):
        xi = calc_2pt(data_slice, randoms_slice, config_fname, do3D,
                      ra_var=ra_var, dec_var=dec_var,
                      data_zvar=zvar, random_zvar=random_zvar,)
    output = calc_3pt(data_slice, randoms_slice, config_fname, do3D,
                      ra_var=ra_var, dec_var=dec_var,
                      data_zvar=zvar, random_zvar=random_zvar, outvar=outvar)

    xi_file_name = config_fname + \
        '_MICEdset'+str(dset_id)+'_jk'+str(jk_id)+'_sigma'+str(sigma_z) + \
        '_'+str(zvar)+'_'+str(min_z)+'_'+str(max_z) + \
        '_rsamp'+str(random_oversamp)+'.xi'
    output_file_name = config_fname+'_MICEdset' + \
        str(dset_id)+'_jk'+str(jk_id)+'_sigma'+str(sigma_z) + \
        '_'+str(zvar)+'_'+str(min_z)+'_'+str(max_z) + \
        '_rsamp'+str(random_oversamp)+'.'+outvar

    if not do3D:
        if (outvar == 'zeta') | (outvar == 'ddd'):
            np.save(os.path.join(paths.ang_out_dir, xi_file_name), xi)
        np.save(os.path.join(paths.ang_out_dir, output_file_name), output)
    else:
        if (outvar == 'zeta') | (outvar == 'ddd'):
            np.save(os.path.join(paths.corr_out_dir, xi_file_name), xi)
        np.save(os.path.join(paths.corr_out_dir, output_file_name), output)

def calc_3pt_noisy_photoz_MICEdm(dset_id, jk_id, config_fname, do3D, min_z, max_z, sigma_z, zvar, random_zvar, random_oversamp, rw_scheme, outvar='zeta'):
    data = fits.getdata(paths.dm_mice_y1[dset_id])
    data = data[data['redshift'] < 1.0]

    if min_z == .6:
        weight_data = fits.getdata(paths.rm_mice_y1_HL[0])
    elif min_z == .75:
        weight_data = fits.getdata(paths.rm_mice_y1_HHL[0])
    else:
        weight_data = fits.getdata(paths.rm_mice_y1[0])

    ra_var = 'RA'
    dec_var = 'DEC'
    
    data[dec_var] = -data[dec_var]
    data = data[data[ra_var] > 0]
    data = data[data[ra_var] < 90]
    data = data[data[dec_var] > -60]
    data = data[data[dec_var] < -40]

    weight_data_slice = get_zslice(weight_data, min_z, max_z, rw_scheme)

    if rw_scheme == 'ZSPEC':
        target_cts, target_bins = np.histogram(
            weight_data_slice['ZSPEC'], range=(0, 1), bins=100)   
    else:
        target_cts, target_bins = np.histogram(
            weight_data_slice['ZREDMAGIC']+np.random.normal(scale=weight_data_slice['ZREDMAGIC_E']), range=(0, 1), bins=100)    

    data_slice = downselect_pz(data, target_cts, target_bins, 'redshift', 2)
    randoms = generate_randoms(data_slice, random_oversamp, 'redshift')
    randoms_slice = downselect_pz(randoms, target_cts, target_bins, 'Z', 5)

    #remove jk region
    jk_classifier = pickle.load( open( "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/jk_classifiers/rectangle_0_90_-60_-40_jk.pkl", "rb" ) )
    data_inds = jk_classifier.predict(zip(data_slice[ra_var],data_slice[dec_var]))
    random_inds = jk_classifier.predict(zip(randoms_slice[ra_var],randoms_slice[dec_var]))
    data_slice = data_slice[data_inds != jk_id]
    randoms_slice = randoms_slice[random_inds != jk_id]

    if (outvar == 'zeta') | (outvar == 'ddd'):
        xi = calc_2pt(data_slice, randoms_slice, config_fname, do3D,
                      ra_var=ra_var, dec_var=dec_var,
                      data_zvar=zvar, random_zvar=random_zvar,)
    output = calc_3pt(data_slice, randoms_slice, config_fname, do3D,
                      ra_var=ra_var, dec_var=dec_var,
                      data_zvar=zvar, random_zvar=random_zvar, outvar=outvar)

    xi_file_name = config_fname + \
        '_MICEdmdset'+str(dset_id)+'_jk'+str(jk_id)+'_sigma'+str(sigma_z) + \
        '_'+rw_scheme+'_'+str(min_z)+'_'+str(max_z)+'.xi'
    output_file_name = config_fname+'_MICEdmdset' + \
        str(dset_id)+'_jk'+str(jk_id)+'_sigma'+str(sigma_z) + \
        '_'+rw_scheme+'_'+str(min_z)+'_'+str(max_z)+'.'+outvar

    if not do3D:
        if (outvar == 'zeta') | (outvar == 'ddd'):
            np.save(os.path.join(paths.ang_out_dir, xi_file_name), xi)
        np.save(os.path.join(paths.ang_out_dir, output_file_name), output)
    else:
        if (outvar == 'zeta') | (outvar == 'ddd'):
            np.save(os.path.join(paths.corr_out_dir, xi_file_name), xi)
        np.save(os.path.join(paths.corr_out_dir, output_file_name), output)

def calc_3pt_noisy_photoz_halos(dset_id, config_fname, do3D, min_z, max_z, sigma_z, zvar, random_zvar, random_oversamp, outvar='zeta'):

    data = pd.read_pickle('/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/buzzard_halos/halos-'+str(dset_id)+'.pkl')

    randoms = generate_randoms(data, random_oversamp, 'Z')

    ra_var = 'RA'
    dec_var = 'DEC'

    data[zvar] += np.random.normal(size=len(data), scale=sigma_z)
    data_slice = get_zslice(data, min_z, max_z, zvar)
    randoms_slice = get_zslice(randoms, min_z, max_z, random_zvar)

    # randoms_slice = randoms_slice[np.random.rand(len(randoms_slice)) < (
    #    len(data_slice)/len(randoms_slice)*random_oversamp)]

    if (outvar == 'zeta') | (outvar == 'ddd'):
        xi = calc_2pt(data_slice, randoms_slice, config_fname, do3D,
                      ra_var=ra_var, dec_var=dec_var,
                      data_zvar=zvar, random_zvar=random_zvar,)
    output = calc_3pt(data_slice, randoms_slice, config_fname, do3D,
                      ra_var=ra_var, dec_var=dec_var,
                      data_zvar=zvar, random_zvar=random_zvar, outvar=outvar)

    xi_file_name = config_fname + \
        '_halosdset'+str(dset_id)+'_sigma'+str(sigma_z) + \
        '_'+str(zvar)+'_'+str(min_z)+'_'+str(max_z) + \
        '_rsamp'+str(random_oversamp)+'.xi'
    output_file_name = config_fname+'_halosdset' + \
        str(dset_id)+'_sigma'+str(sigma_z) + \
        '_'+str(zvar)+'_'+str(min_z)+'_'+str(max_z) + \
        '_rsamp'+str(random_oversamp)+'.'+outvar

    if not do3D:
        if (outvar == 'zeta') | (outvar == 'ddd'):
            np.save(os.path.join(paths.ang_out_dir, xi_file_name), xi)
        np.save(os.path.join(paths.ang_out_dir, output_file_name), output)
    else:
        if (outvar == 'zeta') | (outvar == 'ddd'):
            np.save(os.path.join(paths.corr_out_dir, xi_file_name), xi)
        np.save(os.path.join(paths.corr_out_dir, output_file_name), output)


def calc_3pt_noisy_photoz_dm(dset_id, config_fname, do3D, min_z, max_z, sigma_z, zvar, random_zvar, rw_scheme, outvar='zeta'):
    randoms = fits.getdata(paths.dm_y1_randoms)
    data = fits.getdata(paths.dm_y1[dset_id])
    data[zvar] += np.random.normal(size=len(data), scale=sigma_z)

    # do downsampling
    # nbins = 1000
    # data = data[(data['redshift'] > .15) & (data['redshift'] < .8)]
    # a, _ = np.histogram(data['redshift'], bins=nbins, range=(.15, .8))
    # from datasets import buzzard_cosmo
    # z = np.linspace(.15, .8, nbins)
    # vol = buzzard_cosmo.differential_comoving_volume(z)
    # probs_vs_z = vol.value/a
    # p = probs_vs_z[np.digitize(data['redshift'], np.linspace(.15, .8, nbins))]
    # p /= np.sum(p)
    # new_downsample = np.random.choice(data, size=8000000, p=p)
    # c, _ = np.histogram(
    #     new_downsample['redshift'], bins=nbins, range=(.15, .8))
    # new_probs = (
    #     vol.value/c)[np.digitize(new_downsample['redshift'], np.linspace(.15, .8, nbins))]
    # new_probs /= np.sum(new_probs)
    # data = np.random.choice(new_downsample, size=8000000, p=new_probs)
    # end downsampling

    ra_var = 'azim_ang'
    dec_var = 'polar_ang'

    if min_z == .6:
        weight_data = fits.getdata(paths.rm_y1_HL[0])
    else:
        weight_data = fits.getdata(paths.rm_y1[0])
    
    weight_data_slice = get_zslice(weight_data, min_z, max_z, rw_scheme)

    target_cts, target_bins = np.histogram(
        weight_data_slice['ZREDMAGIC']+np.random.normal(scale=weight_data_slice['ZREDMAGIC_E']), range=(0, 1), bins=100)    

    data_slice = downselect_pz(data, target_cts, target_bins, 'redshift', 2)
    randoms_slice = downselect_pz(randoms, target_cts, target_bins, 'Z', 5)

    if (outvar == 'zeta') | (outvar == 'ddd'):
        xi = calc_2pt(data_slice, randoms_slice, config_fname, do3D,
                      ra_var=ra_var, dec_var=dec_var,
                      data_zvar=zvar, random_zvar=random_zvar,)
    output = calc_3pt(data_slice, randoms_slice, config_fname, do3D,
                      ra_var=ra_var, dec_var=dec_var,
                      data_zvar=zvar, random_zvar=random_zvar, outvar=outvar)

    xi_file_name = config_fname + \
        '_dm2x5RWdset'+str(dset_id)+'_sigma'+str(sigma_z) + \
        '_'+rw_scheme+'_'+str(min_z)+'_'+str(max_z)+'.xi'
    output_file_name = config_fname+'_dm2x5RWdset' + \
        str(dset_id)+'_sigma'+str(sigma_z) + \
        '_'+rw_scheme+'_'+str(min_z)+'_'+str(max_z)+'.'+outvar

    if not do3D:
        if (outvar == 'zeta') | (outvar == 'ddd'):
            np.save(os.path.join(paths.ang_out_dir, xi_file_name), xi)
        np.save(os.path.join(paths.ang_out_dir, output_file_name), output)
    else:
        if (outvar == 'zeta') | (outvar == 'ddd'):
            np.save(os.path.join(paths.corr_out_dir, xi_file_name), xi)
        np.save(os.path.join(paths.corr_out_dir, output_file_name), output)


def calc_3pt_noisy_photoz(dset_id, config_fname, do3D, min_z, max_z, sigma_z, zvar, random_zvar, random_oversamp):

    if min_z == .6:
        randoms = fits.getdata(paths.rm_y1_HL_randoms)
        data = fits.getdata(paths.rm_y1_HL[dset_id])
    elif min_z == .75:
        randoms = fits.getdata(paths.rm_y1_HHL_randoms)
        data = fits.getdata(paths.rm_y1_HHL[dset_id])
    else:
        randoms = fits.getdata(paths.rm_y1_randoms)
        data = fits.getdata(paths.rm_y1[dset_id])

    ra_var = 'RA'
    dec_var = 'DEC'
    print len(randoms)

    data[zvar] += np.random.normal(size=len(data), scale=sigma_z)
    data_slice = get_zslice(data, min_z, max_z, zvar)
    randoms_slice = get_zslice(randoms, min_z, max_z, random_zvar)

    randoms_slice = randoms_slice[np.random.rand(len(randoms_slice)) < (
        len(data_slice)/len(randoms_slice)*random_oversamp)]

    xi = calc_2pt(data_slice, randoms_slice, config_fname, do3D,
                  ra_var=ra_var, dec_var=dec_var,
                  data_zvar=zvar, random_zvar=random_zvar,)
    zeta = calc_3pt(data_slice, randoms_slice, config_fname, do3D,
                    ra_var=ra_var, dec_var=dec_var,
                    data_zvar=zvar, random_zvar=random_zvar,)

    xi_file_name = config_fname + \
        '_newbuzzardrmdset'+str(dset_id)+'_sigma'+str(sigma_z) + \
        '_'+str(zvar)+'_'+str(min_z)+'_'+str(max_z) + \
        '_rsamp'+str(random_oversamp)+'.xi'
    zeta_file_name = config_fname+'_newbuzzardrmdset' + \
        str(dset_id)+'_sigma'+str(sigma_z) + \
        '_'+str(zvar)+'_'+str(min_z)+'_'+str(max_z) + \
        '_rsamp'+str(random_oversamp)+'.zeta'

    if not do3D:
        np.save(os.path.join(paths.ang_out_dir, xi_file_name), xi)
        np.save(os.path.join(paths.ang_out_dir, zeta_file_name), zeta)
    else:
        np.save(os.path.join(paths.corr_out_dir, xi_file_name), xi)
        np.save(os.path.join(paths.corr_out_dir, zeta_file_name), zeta)


def calc_3pt_noisy_photoz_y3(dset_id, config_fname, do3D, min_z, max_z, sigma_z, zvar, random_zvar, random_oversamp):

    if min_z == .6:
        randoms = fits.getdata(paths.rm_y3_HL_randoms)
        data = fits.getdata(paths.rm_y3_HL[dset_id])
    elif min_z == .75:
        randoms = fits.getdata(paths.rm_y3_HHL_randoms)
        data = fits.getdata(paths.rm_y3_HHL[dset_id])
    else:
        randoms = fits.getdata(paths.rm_y3_randoms)
        data = fits.getdata(paths.rm_y3[dset_id])

    ra_var = 'RA'
    dec_var = 'DEC'
    print len(randoms)

    data[zvar] += np.random.normal(size=len(data), scale=sigma_z)
    data_slice = get_zslice(data, min_z, max_z, zvar)
    randoms_slice = get_zslice(randoms, min_z, max_z, random_zvar)

    randoms_slice = randoms_slice[np.random.rand(len(randoms_slice)) < (
        len(data_slice)/len(randoms_slice)*random_oversamp)]

    xi = calc_2pt(data_slice, randoms_slice, config_fname, do3D,
                  ra_var=ra_var, dec_var=dec_var,
                  data_zvar=zvar, random_zvar=random_zvar,)
    zeta = calc_3pt(data_slice, randoms_slice, config_fname, do3D,
                    ra_var=ra_var, dec_var=dec_var,
                    data_zvar=zvar, random_zvar=random_zvar,)

    xi_file_name = config_fname + \
        '_rmy3dset'+str(dset_id)+'_sigma'+str(sigma_z) + \
        '_'+str(zvar)+'_'+str(min_z)+'_'+str(max_z) + \
        '_rsamp'+str(random_oversamp)+'.xi'
    zeta_file_name = config_fname+'_rmy3dset' + \
        str(dset_id)+'_sigma'+str(sigma_z) + \
        '_'+str(zvar)+'_'+str(min_z)+'_'+str(max_z) + \
        '_rsamp'+str(random_oversamp)+'.zeta'

    if not do3D:
        np.save(os.path.join(paths.ang_out_dir, xi_file_name), xi)
        np.save(os.path.join(paths.ang_out_dir, zeta_file_name), zeta)
    else:
        np.save(os.path.join(paths.corr_out_dir, xi_file_name), xi)
        np.save(os.path.join(paths.corr_out_dir, zeta_file_name), zeta)


def generate_randoms(data, oversamp, zvar,
                     Ngen=1000000, Ntries_max=10000):

    Ncurrent = 0
    Ntry = 0
    Ntot = int(len(data) * oversamp)

    minra = 0
    maxra = 90
    mindec = -60
    maxdec = -40

    zdist = data[zvar]

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
        ('RA', '>f4'), ('DEC', '>f4'), ('Z', '>f4')])
    randoms['RA'] = random_ra
    randoms['DEC'] = random_dec
    randoms['Z'] = random_z

    if len(randoms) > Ntot:
        inds_to_keep = np.random.choice(
            np.arange(len(randoms)), size=Ntot, replace=False)
        randoms = randoms[inds_to_keep]

    return randoms
