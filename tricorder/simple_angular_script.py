import subprocess
import paths

# do3Ds = [True]

if __name__ == '__main__':
    for i, config_fname in enumerate(['paper13.1', 'paper13.1a', 'paper14.1', 'paper14.1a', 'paper15.1', 'paper15.1a']):
        do3D = False
        for sigma_z in [0]:
            for min_z in [.15, .3, .45, .6, .75]:
                max_z = min_z + 0.15
                for rw_scheme in ['Z']:
                    for dset_id in [0, 1, 2, 3, 5, 21]:
                        for random_oversamp in [3]:
                            command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_halos(" + str(
                                dset_id) + ", '" + config_fname + "', "+str(do3D)+", " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + ",'"+rw_scheme+"','Z',"+str(random_oversamp)+")"
                            print command_str
                            subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out",
                                             "-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])
