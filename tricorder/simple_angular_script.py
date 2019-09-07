import subprocess
import paths
from time import sleep
import os.path

do3Ds = [True, True, True, True, True, True]
outlogpath = "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/logs4/%J.out"
errlogpath = "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/logs4/%J.err"
ncpus = "4"
primary_dset_id = 0
walltime = '47:00'
memlimit = '2000'


def nice_job_submit(do3D, outvar, config_fname, dset_flavor, dset_id, jk_id, sigma_z, rw_scheme, min_z, max_z, random_oversamp, dm_oversamp=None, sleep_time=.1):

    if do3D:
        checkpath = paths.corr_out_dir
    else:
        checkpath = paths.ang_out_dir

    if outvar == 'ddd':
        checkoutvar = 'rdd'  # last one to be written in a ddd job.
    else:
        checkoutvar = outvar

    if dset_flavor == 'dm':
        oversamp_str = str(dm_oversamp) + 'x' + str(random_oversamp)
    elif dset_flavor == 'newbuzzardrm2':
        oversamp_str = str(random_oversamp)
    else:
        raise ValueError('Unknown dset_flavor: '+dset_flavor)

    outfile = checkpath + '/' + config_fname + \
        '_'+dset_flavor+'dset'+str(dset_id)+'_jk'+str(jk_id)+'_sigma'+str(sigma_z) + \
        '_'+rw_scheme+'_'+str(min_z)+'_'+str(max_z) + \
        '_rsamp'+oversamp_str + '.' + checkoutvar

    if os.path.exists(outfile):
        print 'already done; continuing'
        return
    else:
        if dset_flavor == 'dm':
            command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_dm(" + str(
                dset_id) + ", " + str(jk_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + ",'redshift','Z',"+str(dm_oversamp)+","+str(random_oversamp)+", '"+rw_scheme+"', outvar='"+outvar+"')"
        else:
            command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz(" + str(
                dset_id) + ", " + str(jk_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + ",'"+rw_scheme+"','Z',"+str(random_oversamp)+", outvar='"+outvar+"')"

        print command_str
        subprocess.call(["bsub", "-W", walltime, "-n", ncpus, "-C", "1", "-R", "span[hosts=1] rusage[mem="+memlimit+"] select[hname!=deft0001 && hname!=deft0002 && hname!=deft0003 && hname!=deft0004 && hname!=deft0005 && hname!=deft0006 && hname!=deft0007 && hname!=deft0008 && hname!=deft0009 && hname!=deft0010 && hname!=deft0011 && hname!=deft0012 && hname!=deft0013 && hname!=deft0014 && hname!=deft0015 && hname!=deft0016 && hname!=deft0017 && hname!=deft0018 && hname!=deft0019 && hname!=deft0020 && hname!=deft0021 && hname!=deft0022 && hname!=deft0023 && hname!=deft0024 && hname!=deft0025 && hname!=deft0026 && hname!=deft0027 && hname!=deft0028 && hname!=kiso0030 && hname!=kiso0032 && hname!=kiso0033 && hname!=kiso0034 && hname!=kiso0035 && hname!=kiso0036 && hname!=kiso0037 && hname!=kiso0038 && hname!=kiso0039 && hname!=kiso0010 && hname!=kiso0011 && hname!=kiso0012 && hname!=kiso0013 && hname!=kiso0014 && hname!=kiso0015 && hname!=kiso0016 && hname!=kiso0060 && hname!=kiso0017 && hname!=kiso0061 && hname!=kiso0018 && hname!=kiso0062 && hname!=kiso0019 && hname!=kiso0063 && hname!=kiso0064 && hname!=kiso0065 && hname!=kiso0067 && hname!=kiso0068 && hname!=kiso0040 && hname!=kiso0041 && hname!=kiso0042 && hname!=kiso0043 && hname!=kiso0044 && hname!=kiso0045 && hname!=kiso0046 && hname!=kiso0047 && hname!=kiso0048 && hname!=kiso0049 && hname!=kiso0020 && hname!=kiso0021 && hname!=kiso0022 && hname!=kiso0024 && hname!=kiso0025 && hname!=kiso0026 && hname!=kiso0027 && hname!=kiso0028 && hname!=kiso0029 && hname!=kiso0002 && hname!=kiso0003 && hname!=kiso0004 && hname!=kiso0005 && hname!=kiso0006 && hname!=kiso0050 && hname!=kiso0051 && hname!=kiso0008 && hname!=kiso0052 && hname!=kiso0054 && hname!=kiso0055 && hname!=kiso0056 && hname!=kiso0057 && hname!=kiso0058 && hname!=kiso0059 && hname!=bubble0003 && hname!=bubble0004 && hname!=bubble0005 && hname!=bubble0006]",
                         "-o", outlogpath,
                         "-e", errlogpath, "python", "-c", command_str])
        sleep(sleep_time)
        return


if __name__ == '__main__':
    for i, config_fname in enumerate(['fiducial3d_halfu', 'fiducial3d_75u', 'fiducial3d_tolup', 'fiducial3d_toldown', 'fiducial3d_25Mpc', 'fiducial3d_35Mpc']):
        #do3D = False
        for z_width in [0.15]:
            for sigma_z in [0]:
                for min_z in [.15, .45, .3, .6, ]:
                    max_z = min_z + z_width
                    for random_oversamp in [10]:
                        for jk_id in [-1]:
                            for outvar in ['ddd', 'rrr', 'drr', 'rdr', 'rrd']:
                                if sigma_z == 0:
                                    zlist = ['ZREDMAGIC']
                                else:
                                    zlist = ['ZSPEC']

                                for rw_scheme in zlist:

                                    for dm_oversamp in [10]:

                                        if jk_id == -1:
                                            dset_ids = range(len(paths.dm_y1))
                                        else:
                                            dset_ids = [primary_dset_id]
                                        for dset_id in dset_ids:

                                            nice_job_submit(do3Ds[i], outvar, config_fname, 'dm', dset_id,
                                                            jk_id, sigma_z, rw_scheme, min_z, max_z, random_oversamp,
                                                            dm_oversamp=dm_oversamp)

                                    # RMY1
                                    if jk_id == -1:
                                        dset_ids = range(len(paths.rm_y1))
                                    else:
                                        dset_ids = [primary_dset_id]
                                    for dset_id in dset_ids:
                                        nice_job_submit(do3Ds[i], outvar, config_fname, 'newbuzzardrm2', dset_id,
                                                        jk_id, sigma_z, rw_scheme, min_z, max_z, random_oversamp)
