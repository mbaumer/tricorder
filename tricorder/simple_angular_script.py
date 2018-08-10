import subprocess
import paths

#do3Ds = [True]

if __name__ == '__main__':
    for i, config_fname in enumerate(['paper1.1', 'paper2.1', 'paper3.1', 'paper4.1', 'paper5.1', 'paper6.1', 'paper7.1', 'paper8.1', 'paper9.1']):
        do3D = False
        for sigma_z in [0]:
            for min_z in [.15, .3, .45]:
                max_z = min_z + 0.15
                for random_oversamp in [1, 2, 3]:
                    for dset_id in [0]:  # range(len(paths.rm_y1)):
                        # command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz(" + str(
                        #     dset_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + "," + "'ZSPEC','Z',"+str(random_oversamp)+")"
                        # print command_str
                        # subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                        #                  "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])

                        #     command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_mice(" + str(
                        #         dset_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + "," + "'ZSPEC','Z',"+str(random_oversamp)+")"
                        #     print command_str
                        #     subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                        #                      "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])

                        #     command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz(" + str(
                        #         dset_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + "," + "'ZREDMAGIC','Z',"+str(random_oversamp)+")"
                        #     print command_str
                        #     subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                        #                      "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])

                        #     command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_mice(" + str(
                        #         dset_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + "," + "'ZREDMAGIC','Z',"+str(random_oversamp)+")"
                        #     print command_str
                        #     subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                        #                      "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])
                        # for dset_id in range(len(paths.dm_y1)):
                        #    for set_str in ['ddd', 'ddr', 'drd', 'rdd', 'rrd', 'rdr', 'drr', 'rrr']:
                        #        command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_dm(" + str(
                        #            dset_id) + ", '" + config_fname + "', "+str(do3D)+",  " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + "," + "'redshift','Z',outvar='" + set_str + "')"
                        #        print command_str
                        #        subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                        #                         "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])
                        # for dset_id in range(len(paths.lss_y1)):
                        command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_lss(" + str(dset_id) + ", '" + config_fname + "', "+str(
                            do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + "," + "'redshift','REDSHIFT')"
                        print command_str
                        subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                                         "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])
