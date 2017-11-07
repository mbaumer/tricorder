#!/u/ki/mbaumer/anaconda/bin/python
import subprocess

cmd_list = []
oversamp = 1

for delta_z in [0.05, 0.1, 0.15]:
    z = 0.15
    while z + delta_z < .61:
        #cmd_str = "import datasets; dset = datasets.LSSDataset('Buzzard_v1.6_Y1_0_a'); dset.make(1024,"+str(oversamp)+","+str(z)+","+str(z+delta_z)+",0, 100,-60,-40,write=True)"
        # cmd_list.append(cmd_str)
        #cmd_str = "import datasets; dset = datasets.DMDataset('Buzzard_v1.6_Y1_0_a'); dset.make(1024,"+str(oversamp)+","+str(z)+","+str(z+delta_z)+",0, 100,-60,-40,write=True)"
        # cmd_list.append(cmd_str)
        # cmd_str = "import datasets; dset = datasets.RedmagicDataset('Buzzard_v1.6_Y1_0_a',use_spec_z=False); dset.make(1024," + str(
        #    oversamp) + "," + str(z) + "," + str(z + delta_z) + ",0, 100,-60,-40,write=True)"
        # cmd_list.append(cmd_str)

        #cmd_str = "import datasets; dset = datasets.LSSDataset('Buzzard_v1.6_Y3_0_a'); dset.make(1024,"+str(oversamp)+","+str(z)+","+str(z+delta_z)+",0, 360,-90,90,write=True)"
        # cmd_list.append(cmd_str)
        cmd_str = "import datasets; dset = datasets.DMDataset('Buzzard_v1.6_Y3_0_a'); dset.make(1024," + str(
            oversamp) + "," + str(z) + "," + str(z + delta_z) + ",0, 360,-90,90,write=True)"
        cmd_list.append(cmd_str)

        #cmd_str = "import datasets; dset = datasets.MICEDMDataset('MICE'); dset.make(1024,"+str(oversamp)+","+str(z)+","+str(z+delta_z)+",0, 90,-60,-40,write=True)"
        # cmd_list.append(cmd_str)
        #cmd_str = "import datasets; dset = datasets.MICERedmagicDataset('MICE'); dset.make(1024,"+str(oversamp)+","+str(z)+","+str(z+delta_z)+",0, 90,-60,-40,write=True)"
        # cmd_list.append(cmd_str)

        # cmd_str = "import datasets; dset = datasets.RedmagicDataset('Buzzard_v1.6_Y1_0_b'); dset.make(1024," + str(
        #     oversamp) + "," + str(z) + "," + str(z + delta_z) + ",180, 360,-90,0,write=True)"
        # cmd_list.append(cmd_str)
        # cmd_str = "import datasets; dset = datasets.RedmagicDataset('Buzzard_v1.6_Y1_0_c'); dset.make(1024," + str(
        #     oversamp) + "," + str(z) + "," + str(z + delta_z) + ",180, 360,-90,0,write=True)"
        # cmd_list.append(cmd_str)
        # cmd_str = "import datasets; dset = datasets.RedmagicDataset('Buzzard_v1.6_Y1_1_a'); dset.make(1024," + str(
        #     oversamp) + "," + str(z) + "," + str(z + delta_z) + ",180, 360,-90,0,write=True)"
        # cmd_list.append(cmd_str)
        # cmd_str = "import datasets; dset = datasets.RedmagicDataset('Buzzard_v1.6_Y1_1_b'); dset.make(1024," + str(
        #     oversamp) + "," + str(z) + "," + str(z + delta_z) + ",180, 360,-90,0,write=True)"
        # cmd_list.append(cmd_str)
        # cmd_str = "import datasets; dset = datasets.RedmagicDataset('Buzzard_v1.6_Y1_1_c'); dset.make(1024," + str(
        #     oversamp) + "," + str(z) + "," + str(z + delta_z) + ",180, 360,-90,0,write=True)"
        # cmd_list.append(cmd_str)
        # cmd_str = "import datasets; dset = datasets.RedmagicDataset('Buzzard_v1.6_Y1_2_a'); dset.make(1024," + str(
        #     oversamp) + "," + str(z) + "," + str(z + delta_z) + ",180, 360,-90,0,write=True)"
        # cmd_list.append(cmd_str)
        # cmd_str = "import datasets; dset = datasets.RedmagicDataset('Buzzard_v1.6_Y1_2_b'); dset.make(1024," + str(
        #     oversamp) + "," + str(z) + "," + str(z + delta_z) + ",180, 360,-90,0,write=True)"
        # cmd_list.append(cmd_str)
        # cmd_str = "import datasets; dset = datasets.RedmagicDataset('Buzzard_v1.6_Y1_2_c'); dset.make(1024," + str(
        #     oversamp) + "," + str(z) + "," + str(z + delta_z) + ",180, 360,-90,0,write=True)"
        # cmd_list.append(cmd_str)
        # cmd_str = "import datasets; dset = datasets.RedmagicDataset('Buzzard_v1.6_Y1_21_a'); dset.make(1024," + str(
        #     oversamp) + "," + str(z) + "," + str(z + delta_z) + ",180, 360,-90,0,write=True)"
        # cmd_list.append(cmd_str)

        z += delta_z

for cmd_str in cmd_list:
    print cmd_str
    subprocess.call(["bsub", "-W", "08:00", "-R", "rusage[mem=8000]",
                     "python", "-c", cmd_str])
