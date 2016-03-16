import sys
import csv 
import matplotlib.pyplot as plt
import numpy as np
import protein_inference as prot
import os

plotIdx_ = 0

# FOMR = false omission rate
# FDR = false discovery rate

def main():
  percTabBase = "/home/matthew/mergespec/data/103111-Yeast-2hr/percolator_tdc/tab/103111-Yeast-2hr.percolator"
  
  force = False
  
  options = {}
  options["method"] = "bestpept" # fisher, twopept, bestpept, multPEP
  options["targetDecoyAnalysis"] = "picked" # picked, classic, pval
  options["removeSharedPeptides"] = True
  options["proteinGroupingThreshold"] = 1.0
  
  targetOutFN, decoyOutFN = prot.getOutputFN(percTabBase, options)
  
  totalDbProteins = 6721*10
  
  fileName = targetOutFN
  print fileName
  
  if not os.path.isfile(fileName) or force:
    prot.writeProteinFdrs(percTabBase, options)
  
  qvals, fdrs, fomrs = getQvalues(fileName, totalDbProteins)
  numSignificant = sum(1 if qval < 0.01 else 0 for qval in qvals)
  print "#significant proteins Reported 1% FDR:", numSignificant
  print "#significant proteins Observed 1% FDR:", sum(1 if qval < 0.01 else 0 for qval in fdrs)
  print "Observed FDR at reported 1% FDR:", fdrs[numSignificant-1]
  plotQvalues(qvals, fdrs, fomrs)
  
  plt.show()

def getProteinIds(fp):
  name, seq = None, []
  for line in fp:
    line = line.rstrip()
    if line.startswith(">"):
      if name: yield name
      name, seq = line[1:].split(" ")[0], []
    else: seq.append(line)
  if name: yield name 

def getQvalues(fileName, totalProteins):
  file = open(fileName, 'rb') # The input is the protein output file (.tab) from Percolator (-l)
  data = csv.reader(file, delimiter='\t')
  table = [row for row in data]
  x = np.linspace(0, 1, num=1000)
  fp = 0
  tp = 0

  qvals = []
  fdrs = []
  fomrs = []

  fpSeen = False
  pSeen = False
  nextProteinGroup = -1
  
  numAbsentProteins = totalProteins/10*9

  for i in range (1, len(table)):
    reportedQvalue = float(table[i][2])
    
    proteinName = table[i][0] # Protein name
    curProteinGroup = table[i][1] # Protein group ID
    if i < len(table) - 1:
      nextProteinGroup = table[i+1][1]
    else:
      nextProteinGroup = -1
    
    for proteinId in proteinName.split(","):
      proteinAbsent = proteinId.startswith("mimic")
      if proteinAbsent:
        fpSeen = True
      else: 
        pSeen = True
    
    if curProteinGroup == nextProteinGroup: 
      continue 
    else:
      if pSeen: 
        tp = tp + 1
      elif fpSeen:
        fp = fp + 1 # Actual false positives
        if fp < 10:
          print proteinId, float(fp) / (tp + fp), reportedQvalue
        #print float(fp) / (p + fp), p, fp, table[i][0], table[i][2], table[i][4]
      fpSeen = False
      pSeen = False
    
    tn = numAbsentProteins - fp
    fn = totalProteins/10 - tp
    fdr = float(fp) / (tp + fp)
    fomr = float(fn) / (tn + fn)
    # Q-value * Positives = Expected false positives

    qvals.append(reportedQvalue)
    fdrs.append(fdr)
    fomrs.append(fomr)
  return qvals, fdrs, fomrs

def plotQvalues(qvals, fdrs, fomrs):
  global plotIdx_
  
  xlabel = "Reported q value"
  ylabel = "Observed $FDR_A$"

  plotIdx_ += 1
  
  plt.subplot(1,3, plotIdx_)
  x = np.linspace(0, 1, num=1000)
  plt.plot(qvals, fdrs)
  plt.plot(x, x, 'k--')
  plt.axis([0, 1, 0, 1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)
  
  plotIdx_ += 1
  
  plt.subplot(1,3, plotIdx_)
  numIds = sum(1 for fdr in qvals if fdr < 0.01)
  plt.title("Mimic " + " (" + str(numIds) + " proteins at 1\% FDR)", fontsize = 24, fontweight = 'bold')
  plt.plot(qvals, fdrs)
  plt.plot(x, x, 'k--')
  plt.axis([0, 0.1, 0, 0.1])
  plt.xlabel(xlabel, fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)
  
  plotIdx_ += 1
  
  plt.subplot(1,3, plotIdx_)
  plt.plot(fomrs, fdrs)
  plt.axis([0, 0.2, 0, 0.2])
  plt.xlabel("False Ommission Rate", fontsize = 18)
  plt.ylabel(ylabel, fontsize = 18)
   
if __name__ == "__main__":
  main()

