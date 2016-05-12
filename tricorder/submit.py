#!/home/mbaumer/anaconda2/bin/python
import os, sys

def runall(runname):
	logdir = '/home/mbaumer/tricorder/tricorder/logs/'
	rundict = {'d d d': 1}
	#rundict = {'d d d': 1, 'd r r': 12,'r d r' : 12,'r r d': 12,'d d r': 2,'d r d' : 2,'r d d' : 2 ,'r r r' : 47}
	for runstr in rundict.keys():
		timestr = '#SBATCH --time='+str(rundict[runstr])+':00:00 \n'
		runchars = runstr.replace(" ","")
		f = open('temp.sbatch','w')
		configstrs = ['#!/bin/bash \n', 
			timestr,
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
		break
		#os.system('sbatch temp.sbatch')
		#os.system('rm temp.sbatch')

if __name__ == '__main__':
	runall(sys.argv[1])
