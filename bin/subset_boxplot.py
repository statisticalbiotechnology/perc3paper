#!/usr/bin/python

'''
This script contains helper functions to parse percolator in and output files (tab delimited formats)
'''

import os
import sys
import getopt
import csv
import operator
import math
import matplotlib.pyplot as plt

import numpy as np

def main(argv):
  qThresh = 0.01
  baseFN = "/media/hdd/matthew/mergespec/data/Pandey/percolator_tdc_swissprot_partial_digest/tab_subset_scoring_%s/Pandey.percolator.seed%d.tab.%s"
  
  subsetSizes = ["100K", "500K", "1,0M", "5,0M"]
  peptideFull = 298301.0
  psmFull = 7928454.0
  
  peptideCounts, psmCounts = [], []
  for subsetSize in subsetSizes:
    peptideCounts.append([])
    psmCounts.append([])
    for seed in range(1,11):
      peptideCounts[-1].append(countBelowThreshold(baseFN % (subsetSize, seed, "peptides"), qThresh) / peptideFull)
      psmCounts[-1].append(countBelowThreshold(baseFN % (subsetSize, seed, "psms"), qThresh) / psmFull)
  
  plt.subplot(1,2,2)
  plt.plot([0,5], [1, 1], 'k-')
  plt.boxplot(psmCounts)
  plt.xticks([1, 2, 3, 4], subsetSizes)
  plt.ylabel("Fraction PSMs retained, q < 0.01", fontsize = 18)
  plt.xlabel("Subset size", fontsize = 18)
  #plt.ylim([0.96,1.04])
  plt.ylim([0.99,1.01])  

  plt.subplot(1,2,1)
  plt.plot([0,5], [1, 1], 'k-')
  plt.boxplot(peptideCounts)
  plt.xticks([1, 2, 3, 4], subsetSizes)
  plt.ylabel("Fraction peptides retained, q < 0.01", fontsize = 18)
  plt.xlabel("Subset size", fontsize = 18)
  #plt.ylim([0.96,1.04])
  plt.ylim([0.99,1.01])

  plt.tight_layout()
  
  plt.show()

def countBelowThreshold(inputFile, qThresh):
  reader = csv.reader(open(inputFile), delimiter='\t')
  reader.next()
  
  #return np.random.rand()*50+3e5
  count = 0
  for row in reader:
    if float(row[2]) < qThresh:
      count += 1
    else:
      break
  return count

if __name__ == "__main__":
   main(sys.argv[1:])
