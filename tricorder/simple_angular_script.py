import subprocess
import paths
from time import sleep
import os.path

do3Ds = [False, False]
outlogpath = "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/logs4/%J.out"
errlogpath = "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/logs4/%J.err"
ncpus = "4"
primary_dset_id = 0
walltime = '95:00'
memlimit = '2000'


def nice_job_submit(do3D, outvar, config_fname, dset_flavor, dset_id, jk_id, sigma_z, rw_scheme, min_z, max_z, random_oversamp, dm_oversamp=None, sleep_time=20):

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
        subprocess.call(["bsub", "-W", walltime, "-n", ncpus, "-C", "1", "-R", "span[hosts=1] rusage[mem="+memlimit+"]", "-o", outlogpath,
                         "-e", errlogpath, "python", "-c", command_str])
        sleep(sleep_time)
        return


if __name__ == '__main__':
    for i, config_fname in enumerate(['larger_scales1']):
        #do3D = False
        for sigma_z in [0]:
            for min_z in [.45, .3, .6, .15]:
                max_z = min_z + 0.15
                for random_oversamp in [10]:
                    for jk_id in [-1]:
                        for outvar in ['ddd', 'drr', 'rdr', 'rrd', 'rrr']:
                            if sigma_z == 0:
                                zlist = ['ZSPEC']
                            else:
                                zlist = ['ZSPEC']

                            for rw_scheme in zlist:

                                # for dm_oversamp in [10]:

                                #     if jk_id == -1:
                                #         dset_ids = range(len(paths.dm_y1))
                                #     else:
                                #         dset_ids = [primary_dset_id]
                                #     for dset_id in dset_ids:

                                #         nice_job_submit(do3Ds[i], outvar, config_fname, 'dm', dset_id,
                                #                         jk_id, sigma_z, rw_scheme, min_z, max_z, random_oversamp,
                                #                         dm_oversamp=dm_oversamp)

                                # RMY1
                                if jk_id == -1:
                                    dset_ids = range(len(paths.rm_y1))
                                else:
                                    dset_ids = [primary_dset_id]
                                for dset_id in dset_ids:
                                    nice_job_submit(do3Ds[i], outvar, config_fname, 'newbuzzardrm2', dset_id,
                                                    jk_id, sigma_z, rw_scheme, min_z, max_z, random_oversamp)
