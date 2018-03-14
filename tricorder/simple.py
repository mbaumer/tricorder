#!/u/ki/mbaumer/anaconda/bin/python
import matplotlib.pyplot as plt
from astropy.io import fits
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
from simple_script import out_path, config_dir


def load_config(config_fname):
    config_path = os.path.join(
        simple_script.config_dir, config_fname+'.config')
    with open(config_path) as f:
        return yaml.load(f.read())


def get_zslice(data, min_z, max_z, zvar):
    return data[(data[zvar] > min_z) & (data[zvar] < max_z)]


def calc_2pt(data, randoms, config_fname, zvar, random_zvar, ra_var='RA', dec_var='DEC', random_ra_var='RA', random_dec_var='DEC'):
    config_2pt = load_config(config_fname)['2PCF']
    cat = treecorr.Catalog(ra=data[ra_var], dec=data[dec_var],
                           dec_units='degrees', ra_units='degrees',
                           r=datasets.buzzard_cosmo.comoving_distance(data[zvar]).value*datasets.buzzard_cosmo.h)
    random_cat = treecorr.Catalog(ra=randoms[random_ra_var], dec=randoms[random_dec_var],
                                  dec_units='degrees', ra_units='degrees',
                                  r=datasets.buzzard_cosmo.comoving_distance(randoms[random_zvar]).value*datasets.buzzard_cosmo.h)
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


def calc_3pt(data, randoms, config_fname, zvar, random_zvar, ra_var='RA', dec_var='DEC', random_ra_var='RA', random_dec_var='DEC'):
    config_3pt = load_config(config_fname)['3PCF']
    cat = treecorr.Catalog(ra=data[ra_var], dec=data[dec_var],
                           dec_units='degrees', ra_units='degrees',
                           r=datasets.buzzard_cosmo.comoving_distance(data[zvar]).value*datasets.buzzard_cosmo.h)
    random_cat = treecorr.Catalog(ra=randoms[random_ra_var], dec=randoms[random_dec_var],
                                  dec_units='degrees', ra_units='degrees',
                                  r=datasets.buzzard_cosmo.comoving_distance(randoms[random_zvar]).value*datasets.buzzard_cosmo.h)
    print config_3pt
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
    zeta, varzeta = ddd.calculateZeta(
        ddr=ddr, drd=drd, rdd=rdd, rrd=rrd, rdr=rdr, drr=drr, rrr=rrr)
    return zeta


def calc_3pt_noisy_photoz_lss(dset_id, config_fname, min_z, max_z, sigma_z, zvar, random_zvar):
    randoms = np.load(
        '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y1_0_a/lss_sample/data/REDSHIFT0.6_1nsideNonenJack30.dset_randoms.npy')
    randoms = randoms[np.random.rand(len(randoms)) < 0.2]
    data = fits.getdata(paths.lss_y1[dset_id])
    data = data[data['lss-sample'] == 1]
    ra_var = 'RA'
    dec_var = 'DEC'
    print len(randoms)

    data[zvar] += np.random.normal(size=len(data), scale=sigma_z)
    data_slice = get_zslice(data, min_z, max_z, zvar)
    randoms_slice = get_zslice(randoms, min_z, max_z, random_zvar)

    xi = calc_2pt(data_slice, randoms_slice, config_fname, zvar,
                  random_zvar, ra_var=ra_var, dec_var=dec_var)
    zeta = calc_3pt(data_slice, randoms_slice, config_fname,
                    zvar, random_zvar, ra_var=ra_var, dec_var=dec_var)

    xi_file_name = config_fname + '_lssdset' + \
        str(dset_id)+'_sigma'+str(sigma_z)+'_'+str(min_z)+'_'+str(max_z)+'.xi'
    zeta_file_name = config_fname+'_lssdset' +
        str(dset_id)+'_sigma'+str(sigma_z) + \
            '_'+str(min_z)+'_'+str(max_z)+'.zeta'
    np.save(os.path.join(simple_script.out_dir, xi_file_name), xi)
    np.save(os.path.join(simple_script.out_dir, zeta_file_name), zeta)


def calc_3pt_noisy_photoz_dm(dset_id, config_fname, min_z, max_z, sigma_z, zvar, random_zvar):
    randoms = fits.getdata(
        '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/a/buzzard_1.6-6a_run_redmapper_v6.4.18_redmagic_highdens_0.5-10_randoms.fit')
    randoms = randoms[np.random.rand(len(randoms)) < 0.02]
    data = fits.getdata(paths.dm_y1[dset_id])

    ra_var = 'azim_ang'
    dec_var = 'polar_ang'
    print len(randoms)

    data[zvar] += np.random.normal(size=len(data), scale=sigma_z)
    data_slice = get_zslice(data, min_z, max_z, zvar)
    randoms_slice = get_zslice(randoms, min_z, max_z, random_zvar)

    xi = calc_2pt(data_slice, randoms_slice, config_fname, zvar,
                  random_zvar, ra_var=ra_var, dec_var=dec_var)
    zeta = calc_3pt(data_slice, randoms_slice, config_fname,
                    zvar, random_zvar, ra_var=ra_var, dec_var=dec_var)

    xi_file_name = config_fname +
        '_dmdset'+str(dset_id)+'_sigma'+str(sigma_z) + \
            '_'+str(min_z)+'_'+str(max_z)+'.xi'
    zeta_file_name = config_fname+'_dmdset' +
        str(dset_id)+'_sigma'+str(sigma_z) + \
            '_'+str(min_z)+'_'+str(max_z)+'.zeta'
    np.save(os.path.join(simple_script.out_dir, xi_file_name), xi)
    np.save(os.path.join(simple_script.out_dir, zeta_file_name), zeta)


def calc_3pt_noisy_photoz_rm(dset_id, config_fname, min_z, max_z, sigma_z, zvar, random_zvar):
    randoms = fits.getdata(
        '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/a/buzzard_1.6-6a_run_redmapper_v6.4.18_redmagic_highdens_0.5-10_randoms.fit')
    randoms = randoms[np.random.rand(len(randoms)) < 0.02]
    data = fits.getdata(paths.rm_y1[dset_id])

    if dset_id == 0:
        ra_var = 'azim_ang'
        dec_var = 'polar_ang'
    else:
        ra_var = 'RA'
        dec_var = 'DEC'
    print len(randoms)

    data[zvar] += np.random.normal(size=len(data), scale=sigma_z)
    data_slice = get_zslice(data, min_z, max_z, zvar)
    randoms_slice = get_zslice(randoms, min_z, max_z, random_zvar)

    xi = calc_2pt(data_slice, randoms_slice, config_fname, zvar,
                  random_zvar, ra_var=ra_var, dec_var=dec_var)
    zeta = calc_3pt(data_slice, randoms_slice, config_fname,
                    zvar, random_zvar, ra_var=ra_var, dec_var=dec_var)

    xi_file_name = config_fname +
        '_dset'+str(dset_id)+'_sigma'+str(sigma_z) + \
            '_'+str(min_z)+'_'+str(max_z)+'.xi'
    zeta_file_name = config_fname+'_dset' +
        str(dset_id)+'_sigma'+str(sigma_z) + \
            '_'+str(min_z)+'_'+str(max_z)+'.zeta'
    np.save(os.path.join(simple_script.out_dir, xi_file_name), xi)
    np.save(os.path.join(simple_script.out_dir, zeta_file_name), zeta)
