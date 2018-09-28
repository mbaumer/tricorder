import subprocess
import paths

# do3Ds = [True]

if __name__ == '__main__':
    for i, config_fname in enumerate(['paper13.2', 'paper14.2', 'paper15.2', 'paper13.3', 'paper14.3', 'paper15.3']):
        do3D = False
        for sigma_z in [0]:
            for min_z in [.15, .3, .45, .6, .75]:
                max_z = min_z + 0.15
                for rw_scheme in ['ZSPEC', 'ZREDMAGIC']:
                    for dset_id in range(len(paths.dm_y1)):
                        command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_dm(" + str(
                            dset_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + "," + "'redshift','Z','"+rw_scheme+"')"
                        print command_str
                        subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs2/%J.out",
                                         "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs2/%J.err", "python", "-c", command_str])
                    for dset_id in range(len(paths.rm_y1)):
                        for random_oversamp in [.5, 1]:
                            command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz(" + str(
                                dset_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + "," + "'ZSPEC','Z',"+str(random_oversamp)+")"
                            print command_str
                            subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                                             "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])
