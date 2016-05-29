#!/usr/bin/python

'''
This file generates a subset of the PSMs in a percolator input file
'''

import sys
import random
import csv
import os
import heapq
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import math
import protein_inference as prot
import copy
import distinct_colours as dc

def main(argv):
  #percFolder = "/media/storage/mergespec/data/Pandey/percolator_tdc/"
  percFolder = "/media/hdd/matthew/mergespec/data/Pandey/percolator_tdc_swissprot_partial_digest"
  inFile = os.path.join(percFolder, "pin/Pandey.tab")
  totalPsms = 73139505
  force = False
  
  percTabBase = os.path.join(percFolder, "tab_subset_100k", "Pandey.percolator")
  targetProtGroups, decoyProtGroups = prot.getProteinGroupMaps(percTabBase)
  
  numProteins = dict()
  
  numProts = runProteinInference(percTabBase, force, targetProtGroups, decoyProtGroups)
  for method in numProts:
    if method not in numProteins:
      numProteins[method] = []
    numProteins[method].append([numProts[method]]*3)
  
  factors = range(1,13)
  for factor in factors:
    numPsms = totalPsms/(2**factor)
    
    print "numPsms =", numPsms
    
    outFolder = os.path.join(percFolder, "pin/subsample" + str(factor))
    if not os.path.exists(outFolder):
      os.makedirs(outFolder)
    
    for seed in range(1,4):
      baseFN = "Pandey.ss" + str(factor) + ".seed" + str(seed)
      outFile = os.path.join(outFolder, baseFN + ".tab") 
      print outFile
      if not os.path.isfile(outFile) or force:
        createSubset(numPsms, totalPsms, inFile, outFile, seed)
      
      runPercolator(outFile, percFolder, baseFN, force)
      
      percTabBase = os.path.join(percFolder, "tab", baseFN + ".percolator")
      numProts = runProteinInference(percTabBase, force, targetProtGroups, decoyProtGroups)
      for method in numProts:
        if seed == 1:
          numProteins[method].append([])
        numProteins[method][-1].append(numProts[method])
  
  extFactors = [0] + factors
  colors = dc.get_distinct_grey(len(numProteins))
  markers = ['s','D','o','v']
  methods = ["fisher","twopept","bestpept","multPEP"]
  for i, method in enumerate(methods):
    plt.plot(extFactors, [np.mean(x) for x in numProteins[method][::-1]], label = method, color = colors[i], markeredgecolor = colors[i], marker = markers[i], markersize = 8, linewidth = 3)
  
  labelFontSize = 24
  plt.xlabel("Number of PSMs", fontsize = labelFontSize)
  plt.ylabel("Number of protein groups", fontsize = labelFontSize)
  plt.legend(loc = 'lower right', prop={'size':20})
  setAxisFontSize(18)
  plt.xticks(extFactors, [millify(totalPsms/(2**i)) for i in extFactors[::-1]], rotation=60, ha = 'center')
  plt.margins(0.05)
  plt.tight_layout()
  
  plt.show()

def setAxisFontSize(size):
  for tick in plt.gca().xaxis.get_major_ticks():
    tick.label.set_fontsize(size)
  for tick in plt.gca().yaxis.get_major_ticks():
    tick.label.set_fontsize(size)

def millify(n):
  millnames = ['','K','M','B','T']
  
  n = float(n)
  millidx = max(0,min(len(millnames)-1,
                      int(math.floor(0 if n == 0 else math.log10(abs(n))/3))))

  return '{:.0f}{}'.format(n / 10**(3 * millidx), millnames[millidx]) 
  
def createSubset(numPsms, totalPsms, inFile, outFile, seed):
  random.seed(seed)
  
  data = open(inFile,'rb')
  reader = csv.reader(data, delimiter='\t')
    
  header = reader.next()
  defaultDirection = reader.next()
  
  cutoffRand = 0
  randDict = dict()
  writeQueue = list()
  for i, row in enumerate(reader):
    if i % 5000000 == 0:
      print i, "/", totalPsms
    expMass = float(row[3])
    key = (row[0], expMass)
    if key in randDict:
      rand = randDict[key]
      randDict.pop(key, None)
    else:
      rand = random.random()
      randDict[key] = rand
    if len(writeQueue) < numPsms or rand > writeQueue[0][0]:
      if len(writeQueue) < numPsms:
        heapq.heappush(writeQueue, (rand, key))
      else:
        heapq.heappushpop(writeQueue, (rand, key))
  
  keys = set(pair[1] for pair in writeQueue)
  
  data.seek(0)
  reader.next()
  reader.next()
  
  print "Starting to write PSMs"
  writer = csv.writer(open(outFile,'wb'), delimiter='\t')
  writer.writerow(header) # copy header
  writer.writerow(defaultDirection) # copy default direction
  for i, row in enumerate(reader):
    if i % 5000000 == 0:
      print i, "/", totalPsms
    expMass = float(row[3])
    key = (row[0], expMass)
    if key in keys:
      writer.writerow(row)
  
  #while len(writeQueue) > 0:
  #  pair = heapq.heappop(writeQueue)
  #  writer.writerow(pair[1])

def runPercolator(pinFile, percFolder, baseFN, force):
  binFolder = "/home/matthewt/mergespec/software/"
  percolatorBinary = os.path.join(binFolder, "percolator")
  
  if not os.path.isdir(os.path.join(percFolder, "tab")):
    os.mkdir(os.path.join(percFolder, "tab"))
  if not os.path.isdir(os.path.join(percFolder, "stdout")):
    os.mkdir(os.path.join(percFolder, "stdout"))
  
  peptidesOutput = os.path.join(percFolder, "tab", baseFN + ".percolator.tab.peptides")
  peptidesDecoyOutput = os.path.join(percFolder, "tab", baseFN + ".percolator.decoys.tab.peptides")
  logOutput = os.path.join(percFolder, "stdout", baseFN + ".stdout.txt")
  
  print logOutput
  if not os.path.isfile(peptidesOutput) or force:
    cmd = "%s %s -r %s -B %s -N 100000 > %s 2>&1" % (percolatorBinary, pinFile, peptidesOutput, peptidesDecoyOutput, logOutput)
    rc = subprocess.call(cmd, shell=True)
    if rc == 1:
      print("Error while processing " + cmd)
      return 1
    else:
      return 0

def runProteinInference(percTabBase, force, targetProtGroups, decoyProtGroups):
  optionsArray = []
  
  options = {}
  options["method"] = "fisher" # fisher, twopept, bestpept, multPEP
  options["targetDecoyAnalysis"] = "picked" # picked, classic, pval
  options["removeSharedPeptides"] = True
  options["proteinGroupingThreshold"] = 1.0
  optionsArray.append(options)
  
  options = copy.deepcopy(options)
  options["method"] = "twopept" # fisher, twopept, bestpept, multPEP
  optionsArray.append(options)
  
  options = copy.deepcopy(options)
  options["method"] = "bestpept" # fisher, twopept, bestpept, multPEP
  optionsArray.append(options)
  
  options = copy.deepcopy(options)
  options["method"] = "multPEP" # fisher, twopept, bestpept, multPEP
  optionsArray.append(options)
  
  #plt.suptitle("Mimic hm\_yeast", fontsize = 24, fontweight = 'bold')
  
  colors = dc.get_distinct_grey(len(optionsArray))
  numProteins = dict()
  for i, options in enumerate(optionsArray):
    targetOutFN, decoyOutFN = prot.getOutputFN(percTabBase, options)
    
    if not os.path.isfile(targetOutFN) or force:
      prot.writeProteinFdrs(percTabBase, options, targetProtGroups, decoyProtGroups)
    
    qvals, numpos = getQvalues(targetOutFN)
    numProteins[options["method"]] = sum(1 for q in qvals if q < 0.01)
  
  return numProteins

def getQvalues(fileName):
  csv.field_size_limit(sys.maxsize)
  file = open(fileName, 'rb') # The input is the protein output file (.tab) from Percolator (-l)
  reader = csv.reader(file, delimiter='\t')
  reader.next()
  
  qvals, numpos = [], []

  numProteins = 0
  numProteinGroups = 0
  
  for row in reader:
    reportedQvalue = float(row[2])    
    numProteins += len(row[0].split(","))
    numProteinGroups += 1
    
    numpos.append(numProteinGroups)
    qvals.append(reportedQvalue)
  return qvals, numpos
    
if __name__ == "__main__":
   main(sys.argv[1:])
