import subprocess
import paths

do3Ds = [True]

if __name__ == '__main__':
    for i, config_fname in enumerate(['test10_3d']):
        do3D = do3Ds[i]
        for sigma_z in [0, 0.001, 0.002, 0.01, 0.02, 0.03]:
            for min_z in [.15, .3, .45]:
                max_z = min_z + 0.15
                # for dset_id in range(len(paths.rm_y1)):
                #    command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz(" + str(
                #        dset_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + "," + "'ZSPEC','Z')"
                #    print command_str
                #    subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                #                     "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])
                for dset_id in range(len(paths.dm_y1)):
                    for set_str in ['ddd', 'ddr', 'drd', 'rdd', 'rrd', 'rdr', 'drr', 'rrr']:
                        command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_dm(" + str(
                            dset_id) + ", '" + config_fname + "', "+str(do3D)+",  " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + "," + "'redshift','Z',outvar='" + set_str + "')"
                        print command_str
                        subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                                         "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])
                # for dset_id in range(len(paths.lss_y1)):
                #    command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_lss(" + str(
                #        dset_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + "," + "'redshift','REDSHIFT')"
                #    print command_str
                #    subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                #                     "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])
