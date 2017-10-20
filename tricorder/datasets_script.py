#!/u/ki/mbaumer/anaconda/bin/python
import subprocess

cmd_list = []
oversamp = 2

for delta_z in [0.05,0.1,0.15]:
    z = 0.15
    while z + delta_z < .61:
        cmd_str = "import datasets; dset = datasets.LSSDataset(); dset.make(1024,"+str(oversamp)+","+str(z)+","+str(z+delta_z)+",0, 360,-90,90,write=True)"
        cmd_list.append(cmd_str)
        #cmd_str = "import datasets; dset = datasets.DMDataset(); dset.make(1024,"+str(oversamp)+","+str(z)+","+str(z+delta_z)+",0, 100,-60,-40,write=True)"
        #cmd_list.append(cmd_str)
        cmd_str = "import datasets; dset = datasets.RedmagicDataset(); dset.make(1024,"+str(oversamp)+","+str(z)+","+str(z+delta_z)+",0, 360,-90,90,write=True)"
        cmd_list.append(cmd_str)
        z += delta_z 

for cmd_str in cmd_list:
    print cmd_str
    subprocess.call(["bsub", "-W", "08:00", "-R", "rusage[mem=8000]",
                         "python", "-c", cmd_str])