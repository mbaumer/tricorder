#!/home/mbaumer/anaconda2/bin/python2.7
from __future__ import print_function

import sys
sys.path.append('/home/mbaumer/anaconda2/bin/')

import treecorr
import six

import os
import numpy
import ctypes
_treecorr = numpy.ctypeslib.load_library('_treecorr','/home/mbaumer/anaconda2/lib/python2.7/site-packages/treecorr')
_treecorr.SetOMPThreads.restype = ctypes.c_int
_treecorr.SetOMPThreads.argtypes = [ ctypes.c_int ]
import multiprocessing
num_threads = multiprocessing.cpu_count()
print( 'i have this many CPUs available: ', num_threads)
used_threads = _treecorr.SetOMPThreads(num_threads)
print( 'treecorr will use: ', used_threads)
