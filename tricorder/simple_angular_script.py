import subprocess

for config_fname in ['test7_angular','test8_angular']:
    for sigma_z in [0,0.01,0.02,0.03,0.04]:
        for min_z in [.15,.3,.45]:
            max_z = min_z + 0.15
            command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz(" + str(0) + ", '" + config_fname + "', " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + "," + "'REDSHIFT','Z')"
            print command_str
            subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out","-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])
            for dset_id in range(1,37):
            #for dset_id in range(36):
                command_str = "import simple; simple_angular.calc_3pt_noisy_photoz(" + str(dset_id) + ", '" + config_fname + "', " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + "," + "'ZSPEC','Z')"
                #command_str = "import simple_angular; simple_angular.calc_3pt_noisy_photoz_lss(" + str(dset_id) + ", '" + config_fname + "', " + str(min_z) + "," + str(max_z) + "," + str(sigma_z) + "," + "'redshift','REDSHIFT')"
                print command_str
                subprocess.call(["bsub", "-W", "47:00", "-R", "rusage[mem=4000]", "-o", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.out","-e", "/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/logs/%J.err", "python", "-c", command_str])