#!/home/mbaumer/anaconda/bin/python
import sys, subprocess

Nrandoms = 1
Njackknife = 1
method = 'DDD' # or NNN
if method == 'DDD':
	setlist = ['\'d\',\'d\',\'d\'','\'d\',\'d\',\'r\'','\'d\',\'r\',\'d\'','\'r\',\'d\',\'d\'',
	'\'r\',\'r\',\'d\'','\'d\',\'r\',\'r\'','\'r\',\'d\',\'r\'','\'r\',\'r\',\'r\'']
elif method == 'NNN':
	setlist = ['\'n\',\'n\',\'n\'','\'r\',\'r\',\'r\'']
else:
	raise ValueError('method must be NNN or DDD')

def runall(runname):
	for random_set_id in range(Nrandoms):
		for jk_id in range(Njackknife):
			if Njackknife == 1: jk_id = -1	
			for this_set in setlist:
				subprocess.call(["bsub", "-W", "47:00", "python", "-c" ,"import tricorder; tricorder.run_3pt_ana('"+runname+"',"+str(random_set_id)+","+str(jk_id)+","+this_set+")"])

if __name__ == '__main__':
	runall(sys.argv[1])
