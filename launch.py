#!/usr/bin/env python

import subprocess
import sys

def minibatch():
  
  #infile = "3D_spatial_network_shuffled.txt"
  #n = 434874
  #d = 4
  
  infile = "gassensor_shuffled.txt"
  n = 13910
  d = 129
  
  #infile = "USCensus1990.data_shuffled.txt"
  #n = 2458285
  #d = 68
  
  max_no_improvement = 10
  
  for k in [16, 64, 256, 1024]:
  #for k in [1024, 4096, 16384]:
  	for bs in map(lambda x: x * k, [4, 16, 64]):
  		if (bs < n):
  			lc = \
  			"../kmeans_minibatch {infile} {n} {d} {k} {bs} {mni}"\
  			.format(infile = infile, n = n, d = d, k = k, bs = bs, mni = max_no_improvement)
  			print lc
  			subprocess.call(lc, shell=True)

def standard():

  infile = "3D_spatial_network_shuffled.txt"
  n = 434874
  d = 4

  #infile = "gassensor_shuffled.txt"
  #n = 13910
  #d = 129

  #infile = "USCensus1990.data_shuffled.txt"
  #n = 2458285
  #d = 68

  for k in [16, 64, 256, 1024]:
  #for k in [1024, 4096, 16384]:
      lc = \
      "../kmeans_standard {n} {d} 1 {infile} {k}"\
      .format(n = n, d = d, infile = infile, k = k)
      print lc
      subprocess.call(lc, shell=True)


if __name__ == '__main__':
  if len(sys.argv) > 1:
    if sys.argv[1] == 'mini':
      minibatch()
    elif sys.argv[1] == 'standard':
      standard()
    else:
      print "unknown option"
  else:
    print "Usage: " + sys.argv[0] + ' [option]'
    print "option: mini/standard"



    
