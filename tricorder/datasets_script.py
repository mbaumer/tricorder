#!/u/ki/mbaumer/anaconda/bin/python
import subprocess

#need args
cmd_list = ["import datasets; dset = datasets.DMDataset('/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/dark_matter/dm_cat_2600Mpc_data.fits','/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/redmagicHD/buzzard-v1.1-y1a1-full_run_redmapper_v6.4.13_redmagic_highdens_0.5-10_vlim_zmask.fit'); dset.make(1024,.15,.2,0, 90,-60,-40)",
            "import datasets; dset = datasets.DMDataset('/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/dark_matter/dm_cat_2600Mpc_data.fits','/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/redmagicHD/buzzard-v1.1-y1a1-full_run_redmapper_v6.4.13_redmagic_highdens_0.5-10_vlim_zmask.fit'); dset.make(1024,.3,.35,0, 90,-60,-40)",
            "import datasets; dset = datasets.DMDataset('/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/dark_matter/dm_cat_2600Mpc_data.fits','/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/redmagicHD/buzzard-v1.1-y1a1-full_run_redmapper_v6.4.13_redmagic_highdens_0.5-10_vlim_zmask.fit'); dset.make(1024,.45,.5,0, 90,-60,-40)",
            "import datasets; dset = datasets.DMDataset('/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/dark_matter/dm_cat_2600Mpc_data.fits','/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/redmagicHD/buzzard-v1.1-y1a1-full_run_redmapper_v6.4.13_redmagic_highdens_0.5-10_vlim_zmask.fit'); dset.make(1024,.6,.65,0, 90,-60,-40)",
            "import datasets; dset = datasets.DMDataset('/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/dark_matter/dm_cat_2600Mpc_data.fits','/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/redmagicHD/buzzard-v1.1-y1a1-full_run_redmapper_v6.4.13_redmagic_highdens_0.5-10_vlim_zmask.fit'); dset.make(1024,.15,.18,0,90,-60,-40)",
            "import datasets; dset = datasets.DMDataset('/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/dark_matter/dm_cat_2600Mpc_data.fits','/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/redmagicHD/buzzard-v1.1-y1a1-full_run_redmapper_v6.4.13_redmagic_highdens_0.5-10_vlim_zmask.fit'); dset.make(1024,.3,.33,0, 90,-60,-40)",
            "import datasets; dset = datasets.DMDataset('/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/dark_matter/dm_cat_2600Mpc_data.fits','/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/redmagicHD/buzzard-v1.1-y1a1-full_run_redmapper_v6.4.13_redmagic_highdens_0.5-10_vlim_zmask.fit'); dset.make(1024,.45,.48,0,90,-60,-40)",
            "import datasets; dset = datasets.DMDataset('/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/dark_matter/dm_cat_2600Mpc_data.fits','/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/redmagicHD/buzzard-v1.1-y1a1-full_run_redmapper_v6.4.13_redmagic_highdens_0.5-10_vlim_zmask.fit'); dset.make(1024,.6,.63,0, 90,-60,-40)",
            "import datasets; dset = datasets.DMDataset('/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/dark_matter/dm_cat_2600Mpc_data.fits','/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/redmagicHD/buzzard-v1.1-y1a1-full_run_redmapper_v6.4.13_redmagic_highdens_0.5-10_vlim_zmask.fit'); dset.make(1024,.15,.16,0,90,-60,-40)",
            "import datasets; dset = datasets.DMDataset('/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/dark_matter/dm_cat_2600Mpc_data.fits','/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/redmagicHD/buzzard-v1.1-y1a1-full_run_redmapper_v6.4.13_redmagic_highdens_0.5-10_vlim_zmask.fit'); dset.make(1024,.3,.31,0, 90,-60,-40)",
            "import datasets; dset = datasets.DMDataset('/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/dark_matter/dm_cat_2600Mpc_data.fits','/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/redmagicHD/buzzard-v1.1-y1a1-full_run_redmapper_v6.4.13_redmagic_highdens_0.5-10_vlim_zmask.fit'); dset.make(1024,.45,.46,0,90,-60,-40)",
            "import datasets; dset = datasets.DMDataset('/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/dark_matter/dm_cat_2600Mpc_data.fits','/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new/redmagicHD/buzzard-v1.1-y1a1-full_run_redmapper_v6.4.13_redmagic_highdens_0.5-10_vlim_zmask.fit'); dset.make(1024,.6,.61,0, 90,-60,-40)",
            ]

"""
cmd_list = ["import datasets; dset = datasets.LSSDataset(); dset.make(1024,.15,.2,0, 100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.3,.35,0, 100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.45,.5,0, 100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.6,.65,0, 100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.15,.18,0,100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.3,.33,0, 100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.45,.48,0,100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.6,.63,0, 100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.15,.16,0,100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.3,.31,0, 100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.45,.46,0,100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.6,.61,0, 100,-60,-40)",
            ]
"""

for cmd_str in cmd_list:
    subprocess.call(["bsub", "-W", "08:00", "-R", "rusage[mem=4000]",
                         "python", "-c", cmd_str])