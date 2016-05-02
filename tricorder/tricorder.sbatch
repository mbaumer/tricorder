#!/bin/bash 
#SBATCH --job-name=3pt
#SBATCH --output=3pt.out
#SBATCH --error=3pt.err
#SBATCH --time=01:00 
#formta for time is HH:MM
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -p iric

srun /home/mbaumer/anaconda2/bin/python2.7 tricorder.py r r r
