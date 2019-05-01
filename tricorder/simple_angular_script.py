import subprocess
import paths

# do3Ds = [True]
outlogpath = "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/logs2/%J.out"
errlogpath = "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/logs2/%J.err"
ncpus = "2"
primary_dset_id = 0

if __name__ == '__main__':
    for i, config_fname in enumerate(['paper13.1', 'paper14.1']):
        do3D = False
        for sigma_z in [0]:
            for min_z in [.15, .3, .45, .6]:
                max_z = min_z + 0.15
                for random_oversamp in [10, 20]:
                    for jk_id in range(15):
                        for outvar in ['ddd', 'drr', 'rdr', 'rrd', 'rrr']:
                            for rw_scheme in ['ZSPEC']:
                                # MICE
                                # command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_mice(" + str(
                                #    dset_id) + ", " + str(jk_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + ",'"+rw_scheme+"','Z',"+str(random_oversamp)+", outvar='"+outvar+"')"
                                #print command_str
                                # subprocess.call(["bsub", "-W", "47:00", "-n", ncpus, "-R", "span[hosts=1]", "-o", outlogpath,
                                #                 "-e", errlogpath, "python", "-c", command_str])
                                for dm_oversamp in [1]:
                                    # command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_MICEdm(" + str(
                                    #     dset_id) + ", " + str(jk_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + ",'redshift','Z',"+str(dm_oversamp)+","+str(random_oversamp)+", '"+rw_scheme+"', outvar='"+outvar+"')"
                                    # print command_str
                                    # subprocess.call(["bsub", "-W", "47:00", "-n", ncpus, "-R", "span[hosts=1]", "-o", outlogpath,
                                    #                  "-e", errlogpath, "python", "-c", command_str])

                                    if jk_id == -1:
                                        dset_ids = range(len(paths.dm_y1))
                                    else:
                                        dset_ids = [primary_dset_id]
                                    for dset_id in dset_ids:
                                        command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_dm(" + str(
                                            dset_id) + ", " + str(jk_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + ",'redshift','Z',"+str(dm_oversamp)+","+str(random_oversamp)+", '"+rw_scheme+"', outvar='"+outvar+"')"
                                        print command_str
                                        subprocess.call(["bsub", "-W", "08:00", "-n", ncpus, "-R", "span[hosts=1]", "-o", outlogpath,
                                                         "-e", errlogpath, "python", "-c", command_str])

                                if jk_id == -1:
                                    dset_ids = range(len(paths.rm_y1))
                                else:
                                    dset_ids = [primary_dset_id]
                                for dset_id in dset_ids:
                                    command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz(" + str(
                                        dset_id) + ", " + str(jk_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + ",'"+rw_scheme+"','Z',"+str(random_oversamp)+", outvar='"+outvar+"')"
                                    print command_str
                                    subprocess.call(["bsub", "-W", "08:00", "-n", ncpus, "-R", "span[hosts=1]", "-o", outlogpath,
                                                     "-e", errlogpath, "python", "-c", command_str])
                            # if jk_id == -1:
                            #    dset_ids = [0, 1, 2, 3]
                            # else:
                            #    dset_ids = [primary_dset_id]
                            # for dset_id in dset_ids:
                            #    command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_halos(" + str(
                            #        dset_id) + ", " + str(jk_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + ",'Z','Z',"+str(random_oversamp)+", outvar='"+outvar+"')"
                            #    print command_str
                            #    subprocess.call(["bsub", "-W", "47:00", "-n", ncpus, "-R", "span[hosts=1]", "-o", outlogpath,
                            #                     "-e", errlogpath, "python", "-c", command_str])
