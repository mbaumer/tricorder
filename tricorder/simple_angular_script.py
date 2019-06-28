import subprocess
import paths

# do3Ds = [True]
outlogpath = "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/logs3/%J.out"
errlogpath = "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/logs3/%J.err"
ncpus = "2"
primary_dset_id = 0

if __name__ == '__main__':
    for i, config_fname in enumerate(['newpaper13.1', 'newpaper13.1_newu', 'newpaper13.1_newu2', 'newpaper14.1', 'paper15.1']):
        do3D = False
        for sigma_z in [0]:
            for min_z in [.15,.3,.45]:
                max_z = min_z + 0.15
                for random_oversamp in [20]:
                    for jk_id in [-1]:
                        # , 'drr', 'rdr', 'rrd', 'rrr']:
                        for outvar in ['rrr']:
                            if sigma_z == 0:
                                zlist = ['ZSPEC']
                            else:
                                zlist = ['ZSPEC']

                            for rw_scheme in zlist:
                                # MICE
                                # command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_mice(" + str(
                                #    dset_id) + ", " + str(jk_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + ",'"+rw_scheme+"','Z',"+str(random_oversamp)+", outvar='"+outvar+"')"
                                #print command_str
                                # subprocess.call(["bsub", "-W", "47:00", "-n", ncpus, "-C", "1", "-R", "span[hosts=1]", "-o", outlogpath,
                                #                 "-e", errlogpath, "python", "-c", command_str])
                                # for dm_oversamp in [12]:
                                    # command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_MICEdm(" + str(
                                    #     dset_id) + ", " + str(jk_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + ",'redshift','Z',"+str(dm_oversamp)+","+str(random_oversamp)+", '"+rw_scheme+"', outvar='"+outvar+"')"
                                    # print command_str
                                    # subprocess.call(["bsub", "-W", "47:00", "-n", ncpus, "-C", "1", "-R", "span[hosts=1]", "-o", outlogpath,
                                    #                  "-e", errlogpath, "python", "-c", command_str])

                                    # if jk_id == -1:
                                    #     dset_ids = range(len(paths.dm_y1))
                                    # else:
                                    #     dset_ids = [primary_dset_id]
                                    # for dset_id in dset_ids:
                                    #     command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_dm(" + str(
                                    #         dset_id) + ", " + str(jk_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + ",'redshift','Z',"+str(dm_oversamp)+","+str(random_oversamp)+", '"+rw_scheme+"', outvar='"+outvar+"')"
                                    #     print command_str
                                    #     subprocess.call(["bsub", "-W", "18:00", "-n", ncpus, "-C", "1", "-R", "span[hosts=1]", "-o", outlogpath,
                                    #                      "-e", errlogpath, "python", "-c", command_str])

                                #     if min_z != .45:
                                #         continue

                                #     if jk_id == -1:
                                #         dset_ids = range(len(paths.dm_y1))
                                #     else:
                                #         dset_ids = [primary_dset_id]
                                #     for dset_id in dset_ids:
                                #         command_str = "import simple_angular; simple_angular.calc_3pt_randxrand(" + str(
                                #             dset_id) + ", " + str(jk_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + ",'Z','Z',"+str(dm_oversamp)+","+str(random_oversamp)+", '"+rw_scheme+"', outvar='"+outvar+"')"
                                #         print command_str
                                #         subprocess.call(["bsub", "-W", "12:00", "-n", ncpus, "-C", "1", "-R", "span[hosts=1]", "-o", outlogpath,
                                #                          "-e", errlogpath, "python", "-c", command_str])

                                # RMY1
                                # if jk_id == -1:
                                #    dset_ids = range(len(paths.rm_y1))
                                # else:
                                dset_ids = [primary_dset_id]
                                for dset_id in dset_ids:
                                    command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz(" + str(
                                        dset_id) + ", " + str(jk_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + ",'"+rw_scheme+"','Z',"+str(random_oversamp)+", outvar='"+outvar+"')"
                                    print command_str
                                    subprocess.call(["bsub", "-W", "24:00", "-n", ncpus, "-C", "1", "-R", "span[hosts=1]", "-o", outlogpath,
                                                     "-e", errlogpath, "python", "-c", command_str])

                                # Y3
                                # if jk_id == -1:
                                #     dset_ids = range(len(paths.rm_y3))
                                # else:
                                #     dset_ids = [primary_dset_id]
                                # for dset_id in dset_ids:
                                #     command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_y3(" + str(
                                #         dset_id) + ", " + str(jk_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + ",'"+rw_scheme+"','Z',"+str(random_oversamp)+", outvar='"+outvar+"')"
                                #     print command_str
                                #     subprocess.call(["bsub", "-W", "24:00", "-n", ncpus, "-C", "1", "-R", "span[hosts=1]", "-o", outlogpath,
                                #                      "-e", errlogpath, "python", "-c", command_str])
                            # if jk_id == -1:
                            #     dset_ids = [0, 1, 2, 3]
                            # else:
                            #     dset_ids = [primary_dset_id]
                            # for dset_id in dset_ids:
                            #     command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_halos(" + str(
                            #         dset_id) + ", " + str(jk_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + ",'Z','Z',"+str(random_oversamp)+", outvar='"+outvar+"')"
                            #     print command_str
                            #     subprocess.call(["bsub", "-W", "12:00", "-n", ncpus, "-C", "1", "-R", "span[hosts=1]", "-o", outlogpath,
                            #                      "-e", errlogpath, "python", "-c", command_str])