#!/home/mbaumer/anaconda2/bin/python
import os, sys

def runall(runname):
	logdir = '/home/mbaumer/tricorder/tricorder/logs/'
	runstrs = ['d d d','d r r','r d r','r r d','d d r','d r d','r d d','r r r']
	for runstr in runstrs:
		runchars = runstr.replace(" ","")
		f = open('temp.sbatch','w')
		configstrs = ['#!/bin/bash \n', 
			'#SBATCH --time=12:00:00 \n',
			'#SBATCH --qos=iric \n',
			'#SBATCH --nodes=1 \n',
			'#SBATCH --mem-per-cpu=4000 \n',
			'#SBATCH --ntasks-per-node=1 \n',
			'#SBATCH --cpus-per-task=16 \n',
			'#SBATCH -p iric \n',
			'#SBATCH --job-name='+runname+'_'+runchars+' \n',
			'#SBATCH --output='+logdir+runname+'_'+runchars+'.out'+' \n',
			'#SBATCH --error='+logdir+runname+'_'+runchars+'.err'+' \n',
			'srun /home/mbaumer/anaconda2/bin/python2.7 tricorder.py '+runstr+' '+runname+' \n' ]
		f.writelines(configstrs)
		f.close()
		os.system('sbatch temp.sbatch')
		os.system('rm temp.sbatch')

if __name__ == '__main__':
	runall(sys.argv[1])
