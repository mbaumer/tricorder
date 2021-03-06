import subprocess
import sys
from glob import glob
import numpy as np

dset_fname_vec = ['/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y1_0_a/redmagicHD/data/ZREDMAGIC0.2_0.25nside1024nJack30.dset',
                  '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y1_0_a/redmagicHD/data/ZREDMAGIC0.35_0.4nside1024nJack30.dset',
                  '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y1_0_a/redmagicHD/data/ZREDMAGIC0.5_0.55nside1024nJack30.dset',
                  '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y1_0_a/redmagicHD/data/ZSPEC0.2_0.3nside1024nJack30.dset',
                  '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y1_0_a/redmagicHD/data/ZSPEC0.35_0.45nside1024nJack30.dset',
                  '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y1_0_a/redmagicHD/data/ZSPEC0.5_0.6nside1024nJack30.dset',
                  '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y1_0_a/redmagicHD/data/ZSPEC0.15_0.3nside1024nJack30.dset',
                  '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y1_0_a/redmagicHD/data/ZSPEC0.3_0.45nside1024nJack30.dset',
                  '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y1_0_a/redmagicHD/data/ZSPEC0.45_0.6nside1024nJack30.dset',
                  '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y1_0_a/redmagicHD/data/ZREDMAGIC0.15_0.3nside1024nJack30.dset',
                  '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y1_0_a/redmagicHD/data/ZREDMAGIC0.3_0.45nside1024nJack30.dset',
                  '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y1_0_a/redmagicHD/data/ZREDMAGIC0.45_0.6nside1024nJack30.dset']

dset_fname_vec += ['/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y3_0_a/redmagicHD/data/ZREDMAGIC0.15_0.25nside1024nJack30.dset',
                   '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y3_0_a/redmagicHD/data/ZREDMAGIC0.35_0.4nside1024nJack30.dset',
                   '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y3_0_a/redmagicHD/data/ZREDMAGIC0.5_0.55nside1024nJack30.dset',
                   '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y3_0_a/redmagicHD/data/ZSPEC0.2_0.3nside1024nJack30.dset',
                   '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y3_0_a/redmagicHD/data/ZSPEC0.35_0.45nside1024nJack30.dset',
                   '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y3_0_a/redmagicHD/data/ZSPEC0.5_0.6nside1024nJack30.dset',
                   '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y3_0_a/redmagicHD/data/ZSPEC0.15_0.3nside1024nJack30.dset',
                   '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y3_0_a/redmagicHD/data/ZSPEC0.3_0.45nside1024nJack30.dset',
                   '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y3_0_a/redmagicHD/data/ZSPEC0.45_0.6nside1024nJack30.dset',
                   '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y3_0_a/redmagicHD/data/ZREDMAGIC0.15_0.3nside1024nJack30.dset',
                   '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y3_0_a/redmagicHD/data/ZREDMAGIC0.3_0.45nside1024nJack30.dset',
                   '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y3_0_a/redmagicHD/data/ZREDMAGIC0.45_0.6nside1024nJack30.dset']


# dset_fname_vec += glob(
#    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y1_0_a/lss_sample/data/RED*.dset')
dset_fname_vec += glob(
    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y3_0_a/dark_matter/data/RED*.dset')

# dset_fname_vec += [
#    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/MICE/dark_matter/data/REDSHIFT0.15_0.3nside1024nJack30.dset',
#    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/MICE/dark_matter/data/REDSHIFT0.3_0.45nside1024nJack30.dset',
#    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/MICE/dark_matter/data/REDSHIFT0.45_0.6nside1024nJack30.dset']


# dset_fname_vec += glob(
#    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y1_0_a/halos/data/RED*.dset')

# dset_fname_vec += glob(
#    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y3_0_a/redmagicHD/data/ZSPEC*.dset')
# dset_fname_vec += glob(
#    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.6_Y3_0_a/lss_sample/data/RED*.dset')

# dset_fname_vec += ['/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/MICE/redmagicHD/data/ZSPEC0.15_0.3nside1024nJack30.dset',
#                   '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/MICE/redmagicHD/data/ZSPEC0.3_0.45nside1024nJack30.dset',
#                   '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/MICE/redmagicHD/data/ZSPEC0.45_0.6nside1024nJack30.dset', ]
# dset_fname_vec += glob(
#    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/MICE/dark_matter/data/RED*.dset')

print len(dset_fname_vec)

#'/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test2_proj_rpar50.config',
#'/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test2_proj_rpar100.config',
#'/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test2_3d.config',

# config_fname_pointvec = []
config_fname_pointvec = ['/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test3_3d.config',
                         '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test7_3d.config',
                         '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test10_3d.config'
                         ]
#                        ]
#                         '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test7_proj_rpar20.config',
#                         ,]

# config_fname_pixvec = []
config_fname_pixvec = ['/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test7_angular.config',
                       '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test5_angular.config',
                       '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test9_angular.config']
#                       ]

for dset_fname in dset_fname_vec:

    for config_fname in config_fname_pixvec:
        for jk_id in [-1]:
            command_str = "import tricorder; corr = tricorder.PixelCorrelation('" + \
                dset_fname + "', '" + config_fname + \
                "'," + str(jk_id) + "); corr.run()"
            subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                             "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])

    for config_fname in config_fname_pointvec:
        for jk_id in [-1]:
            for set_str in ['ddd', 'ddr', 'drd', 'rdd', 'rrd', 'rdr', 'drr', 'rrr']:
                command_str = "import tricorder; corr = tricorder.PointCorrelation('" + \
                    dset_fname + "', '" + config_fname + \
                    "', '" + set_str + "', " + \
                    str(jk_id) + "); corr.run()"
                subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                                 "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])
