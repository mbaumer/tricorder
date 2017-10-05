import subprocess
from glob import glob

dset_fname_vec = glob(
    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.1/lss_sample/data/RED*.dset')
dset_fname_vec += glob(
    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.1/dark_matter/data/DIST*.dset')
dset_fname_vec += glob(
    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.1/redmagicHD/data/ZSPEC*.dset')

#'/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test2_3d.config',
config_fname_vec = ['/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test2_proj_rpar50.config',
                    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test2_proj_rpar100.config',
                    ]

for dset_fname in dset_fname_vec:
    for config_fname in config_fname_vec:
        for set_str in ['ddd', 'ddr', 'drd', 'rdd', 'rrd', 'rdr', 'drr', 'rrr']:
            command_str = "import tricorder; corr = tricorder.PointCorrelation('" + \
                dset_fname + "', '" + config_fname + \
                "', '" + set_str + "', " + \
                str(-1) + "); corr.run()"
            subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=2000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                             "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])
            
"""
config_fname = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test_angular.config'
    command_str = "import tricorder; corr = tricorder.PixelCorrelation('" + \
        dset_fname + "', '" + config_fname + \
        "'," + str(-1) + "); corr.run()"
    subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=2000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                     "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])
"""