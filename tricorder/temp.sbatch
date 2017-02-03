#!/bin/bash 
#SBATCH --time=1:00:00 
#SBATCH --qos=iric 
#SBATCH --nodes=1 
#SBATCH --mem-per-cpu=4000 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=16 
#SBATCH -p iric 
#SBATCH --job-name=test_ddd 
#SBATCH --output=/home/mbaumer/tricorder/tricorder/logs/test_ddd.out 
#SBATCH --error=/home/mbaumer/tricorder/tricorder/logs/test_ddd.err 
srun /home/mbaumer/anaconda2/bin/python2.7 tricorder.py d d d test 
