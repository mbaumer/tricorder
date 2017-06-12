#!/home/mbaumer/anaconda/bin/python
import sys, subprocess
from os.path import expandvars
import yaml
import numpy as np

outdir = expandvars('$DES_SIMS')+'new_triplet_counts/'

setlist = ['\'d\',\'d\',\'d\'','\'d\',\'d\',\'r\'','\'d\',\'r\',\'d\'','\'r\',\'d\',\'d\'',
    '\'r\',\'r\',\'d\'','\'d\',\'r\',\'r\'','\'r\',\'d\',\'r\'','\'r\',\'r\',\'r\'']

def make_config(lower_z_lim,delta_z,zvar,metric,param_3D):

    configdict = {}

    configdict['zvar'] = zvar
    configdict['param_3D'] = param_3D
    configdict['metric'] = metric

    if metric == 'Rperp':
        assert configdict['param_3D'] != 0
        configdict['min_rpar'] = -1*configdict['param_3D']
        configdict['max_rpar'] = 1*configdict['param_3D']

    configdict['outdir'] = outdir
    if zvar == 'DISTANCE':
        configdict['data_path'] = expandvars('$DES_SIMS')+'/processed_input/dark_matter/downsampled_dm_spt_data.npy'
        configdict['randoms_path'] = expandvars('$DES_SIMS')+'/processed_input/dark_matter/downsampled_dm_spt_randoms.npy'
    elif zvar == 'ZSPEC':
        configdict['data_path'] = expandvars('$DES_SIMS')+'/processed_input/redmagic_gals/buzzard_redmagic_spt_zspec_data.npy'
        configdict['randoms_path'] = expandvars('$DES_SIMS')+'/processed_input/redmagic_gals/buzzard_redmagic_spt_zspec_randoms.npy'
    else:
        configdict['data_path'] = expandvars('$DES_SIMS')+'/processed_input/redmagic_gals/buzzard_redmagic_spt_zredmagic_data.npy'
        configdict['randoms_path'] = expandvars('$DES_SIMS')+'/processed_input/redmagic_gals/buzzard_redmagic_spt_zredmagic_randoms.npy'

    configdict['min_z'] = lower_z_lim
    configdict['max_z'] = lower_z_lim + delta_z

    #treecorr params
    configdict['min_sep'] = 3
    configdict['max_sep'] = 40
    configdict['nbins'] = 200
    configdict['min_u'] = 0
    configdict['max_u'] = 1
    configdict['nubins'] = 100
    configdict['min_v'] = -1
    configdict['max_v'] = 1
    configdict['nvbins'] = 200
    configdict['bin_slop'] = 0.1
    if param_3D == 0:
        configdict['sep_units'] = 'arcmin'
    else: 
        try:
            del configdict['sep_units']
        except KeyError:
            pass
    
    configdict['runname'] = zvar+str(lower_z_lim)+'_deltaz'+str(delta_z)+'_'+metric

    #write it out so we remember what we did
    config_fname = outdir+configdict['runname']+'.yaml'
    #if (!os.path.exists(config_fname)): #not atomic; hard code for now
    f = open(config_fname,'w')
    f.write(yaml.dump(configdict))
    f.close()

    return config_fname

def runall(min_z, max_z, delta_z, zvar, metric, param_3D):
    for lower_z_lim in np.arange(min_z,max_z,delta_z):
        config_fname = make_config(lower_z_lim,delta_z,zvar,metric,param_3D)
        for this_set in setlist:
            print "bsub", "-W", "47:00", "python", "-c" ,"import tricorder; tricorder.run_3pt_ana('"+config_fname+"',"+this_set+")"
            subprocess.call(["bsub", "-W", "47:00", "python", "-c" ,"import tricorder; tricorder.run_3pt_ana('"+config_fname+"',"+this_set+")"])

#if __name__ == '__main__':
#    runall(sys.argv[1],sys.argv[1],sys.argv[1],sys.argv[1],sys.argv[1],sys.argv[1])
