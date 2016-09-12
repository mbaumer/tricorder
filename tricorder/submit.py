#!/home/mbaumer/anaconda2/bin/python
import os, sys
import platform
logdir = '/home/mbaumer/tricorder/tricorder/logs/'

def runall(runname):
	if 'ki-ls' in platform.node():
		for random_set_id in range(31):
			os.system('bsub -W 16:00 python tricorder.py '+runname+' '+str(random_set_id))
	elif 'sh-' in platform.node():
		for random_set_id in range(31):
			f = open('temp.sbatch','w')
			configstrs = ['#!/bin/bash \n', 
				'#SBATCH --time=2:00:00 \n',
				'#SBATCH --qos=iric \n',
				'#SBATCH --nodes=1 \n',
				'#SBATCH --mem-per-cpu=4000 \n',
				'#SBATCH --ntasks-per-node=1 \n',
				'#SBATCH --cpus-per-task=16 \n',
				'#SBATCH -p iric \n',
				'#SBATCH --job-name='+runname+'_'+random_set_id+' \n',
				'#SBATCH --output='+logdir+runname+'_'+random_set_id+'.out'+' \n',
				'#SBATCH --error='+logdir+runname+'_'+random_set_id+'.err'+' \n',
				'srun /home/mbaumer/anaconda2/bin/python2.7 tricorder.py '+runname+' '+str(random_set_id)+'\n' ]
			f.writelines(configstrs)
			f.close()
			os.system('sbatch temp.sbatch')
			os.system('rm temp.sbatch')

if __name__ == '__main__':
	runall(sys.argv[1])
