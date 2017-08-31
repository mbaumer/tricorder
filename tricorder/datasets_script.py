#!/u/ki/mbaumer/anaconda/bin/python
import subprocess

cmd_list = ["import datasets; dset = datasets.LSSDataset(); dset.make(1024,.15,.3,0,100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.3,.45,0,100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.45,.6,0,100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.6,.75,0,100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.75,.9,0,100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.1,.2,0,100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.2,.3,0,100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.3,.4,0,100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.4,.5,0,100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.5,.6,0,100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.6,.7,0,100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.7,.8,0,100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.8,.9,0,100,-60,-40)",
            "import datasets; dset = datasets.LSSDataset(); dset.make(1024,.9,1.0,0,100,-60,-40)"]

for cmd_str in cmd_list:
    subprocess.call(["bsub", "-W", "08:00", "-R", "rusage[mem=4000]",
                         "python", "-c", cmd_str])