from glob import glob
import subprocess

dset_fname_vec = glob('/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.1/lss_sample/data/RED*')
dset_fname_vec += glob('/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.1/dark_matter/data/DIST*')
dset_fname_vec += glob('/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/Buzzard_v1.1/redmagicHD/data/ZSPEC*')
"""
config_fname_vec = ['/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test_3d.config',
                    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test_proj_rpar50.config',
                    '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test_proj_rpar100.config',
                   ]
"""
config_fname_vec = ['/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/test_angular.config']

for dset_fname in dset_fname_vec:
    for config_fname in config_fname_vec:
        command_str = "import tricorder; corr = tricorder.PixelCorrelation('" + \
            dset_fname + "', '" + config_fname + \
            "'," + str(-1) + "); corr.run()"
        subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=2000]",
                         "python", "-c", command_str])
        
        for set_str in ['ddd','ddr','drd','rdd','rrd','rdr','drr','rrr']:
            command_str = "import tricorder; corr = tricorder.PointCorrelation('" + \
            dset_fname + "', '" + config_fname + \
            "', '" + set_str + "', " + \
            str(-1) + "); corr.run()"
            subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=2000]",
                         "python", "-c", command_str])