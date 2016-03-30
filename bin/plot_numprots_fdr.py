import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np
import protein_inference as prot
import os
import copy
import distinct_colours as dc

# FOMR = false omission rate
# FDR = false discovery rate

def main():
  #percTabBase = "/media/storage/mergespec/data/Pandey/percolator_tdc/tab_subset_500k/Pandey.percolator"
  percTabBase = "/media/storage/mergespec/data/Pandey/percolator_tdc_uniprot/tab_uppmax/Pandey.percolator"
  
  force = False
  
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
  
  for i, options in enumerate(optionsArray):
    targetOutFN, decoyOutFN = prot.getOutputFN(percTabBase, options)
    
    if not os.path.isfile(targetOutFN) or force:
      prot.writeProteinFdrs(percTabBase, options)
    
    qvals, numpos = getQvalues(targetOutFN)
    plotQvalues(qvals, numpos, options, colors[i])
  
  labelFontSize = 30
  if "uniprot" in percTabBase:
    maxY = 7000
  else:
    maxY = 14000
  
  plt.figure(1)
  plt.plot([0.01, 0.01], [0, maxY], 'k', linestyle = 'dotted')
  plt.xlim([0, 0.05])
  plt.ylim([0, maxY])
  plt.xlabel("Decoy FDR", fontsize = labelFontSize)
  plt.ylabel("Number of protein groups", fontsize = labelFontSize)
  plt.legend(loc = 'lower right', prop={'size':24})
  setAxisFontSize(20)
  plt.tight_layout()
  
  plt.show()

def setAxisFontSize(size):
  for tick in plt.gca().xaxis.get_major_ticks():
    tick.label.set_fontsize(size)
  for tick in plt.gca().yaxis.get_major_ticks():
    tick.label.set_fontsize(size)
    
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

def plotQvalues(qvals, numpos, options, color):
  figIdx = 1
  
  plt.figure(figIdx)
  plt.plot(qvals, numpos, label = options["method"], linewidth = 2, color = color)
   
if __name__ == "__main__":
  main()

