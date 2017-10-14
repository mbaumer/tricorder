import subprocess
from glob import glob

dset_fname_vec = glob(
    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.3a/lss_sample/data/RED*.dset')
dset_fname_vec += glob(
    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.3a/dark_matter/data/RED*.dset')
dset_fname_vec += glob(
    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.3a/redmagicHD/data/ZSPEC*.dset')

#'/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test2_proj_rpar50.config',
#'/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test2_proj_rpar100.config',
#'/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test2_3d.config',

config_fname_pointvec = ['/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test4_angular.config']

config_fname_pixvec = ['/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test4_angular.config',
                       '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test5_angular.config',
                       '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test6_angular.config']

for dset_fname in dset_fname_vec:
    
    for config_fname in config_fname_pixvec:
        command_str = "import tricorder; corr = tricorder.PixelCorrelation('" + \
            dset_fname + "', '" + config_fname + \
            "'," + str(-1) + "); corr.run()"
        subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=2000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                     "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])
    
    for config_fname in config_fname_pointvec:
        for set_str in ['ddd', 'ddr', 'drd', 'rdd', 'rrd', 'rdr', 'drr', 'rrr']:
            command_str = "import tricorder; corr = tricorder.PointCorrelation('" + \
                dset_fname + "', '" + config_fname + \
                "', '" + set_str + "', " + \
                str(-1) + "); corr.run()"
            subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=2000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                             "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])
            