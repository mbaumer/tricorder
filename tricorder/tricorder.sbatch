#!/bin/bash 
#SBATCH --job-name=sv
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --time=01:00:00 
#SBATCH --qos=iric
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH -p iric

srun /home/mbaumer/anaconda2/bin/python2.7 tricorder.py d d d 
#srun /home/mbaumer/anaconda2/bin/python2.7 openmp_test.py
