import subprocess
import paths

if __name__ == '__main__':
    for config_fname in ['test7_angular', 'test8_angular']:
        for sigma_z in [0, 0.01, 0.02, 0.03, 0.04]:
            for min_z in [.15, .3, .45]:
                max_z = min_z + 0.15
                #for dset_id in range(len(paths.rm_y1)):
                #    command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz(" + str(
                #        dset_id) + ", '" + config_fname + "', " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + "," + "'ZSPEC','Z')"
                #    print command_str
                #    subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                #                     "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])
                for dset_id in range(len(paths.dm_y1)):
                    command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_dm(" + str(
                        dset_id) + ", '" + config_fname + "', " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + "," + "'redshift','Z')"
                    print command_str
                    subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                                     "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])
                for dset_id in range(len(paths.lss_y1)):
                    command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_lss(" + str(dset_id) + ", '" + config_fname + "', " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + "," + "'redshift','REDSHIFT')"
                    print command_str
                    subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out","-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])
