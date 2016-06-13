import sys
import csv
import os
import matplotlib.pyplot as plt
import numpy as np
import protein_inference as prot

def main():
  percTabBase = "/media/storage/mergespec/data/103111-Yeast-2hr/percolator_tdc/tab_subset_scoring/103111-Yeast-2hr.percolator"
  
  force = False
  
  optionsArray = []
  
  options = {}
  options["removeSharedPeptides"] = False
  options["proteinGroupingThreshold"] = 0.1
  options["targetDecoyAnalysis"] = "classic"
  options["method"] = "fido" # fisher, twopept, bestpept, multPEP
  optionsArray.append(options)
  
  for i, options in enumerate(optionsArray):
    targetInFN, decoyInFN = prot.getOutputFN(percTabBase, options)
    
    options["targetDecoyAnalysis"] = "exact"
    targetOutFN, decoyOutFN = prot.getOutputFN(percTabBase, options)
    
    fileName = targetOutFN
    print fileName
    
    qvals, repQvals = getQvalues(targetInFN, decoyInFN, targetOutFN, decoyOutFN)
    plt.plot(qvals, repQvals)
    
  x = np.linspace(0, 1, num=1000)  
  plt.plot(x,x)
  plt.show()

def getQvalues(targetFileName, decoyFileName, targetOutFN, decoyOutFN):
  file = open(targetFileName, 'rb') # The input is the protein output file (.tab) from Percolator (-l)
  data = csv.reader(file, delimiter='\t')
  header = data.next()
  table = [row for row in data]
  
  file2 = open(decoyFileName, 'rb') # The input is the protein output file (.tab) from Percolator (-l)
  data2 = csv.reader(file2, delimiter='\t')
  table.extend([row for row in data2][1:])
  
  targetWriter = csv.writer(open(targetOutFN, 'wb'), delimiter = '\t')
  decoyWriter = csv.writer(open(decoyOutFN, 'wb'), delimiter = '\t')
  
  targetWriter.writerow(header)
  decoyWriter.writerow(header)
  
  table = sorted(table, key = lambda x : float(x[3]))
  fp = 0
  tp = 0

  qvals, repQvals = [], []

  fpSeen = False
  pSeen = False
  prevProteinGroup = -1
  
  qvalSum = 0
  numProteins = 0
  printed = False
  
  for i in range(len(table)):
    reportedQvalue = float(table[i][2])
    reportedPEP = float(table[i][3])
    
    proteinName = table[i][0] # Protein name
    curProteinGroup = table[i][1] # Protein group ID
    if i > 0:
      prevProteinGroup = table[i-1][1]
    else:
      prevProteinGroup = -1
    
    if curProteinGroup != prevProteinGroup:
      if not proteinName.startswith("decoy_"):
        qvalSum += reportedPEP
      numProteins += 1
    qvals.append(qvalSum / numProteins)
    repQvals.append(reportedQvalue)
    
    row = table[i]
    row[2] = qvalSum / numProteins
    if not proteinName.startswith("decoy_"):
      targetWriter.writerow(row)
    else:
      decoyWriter.writerow(row)
    
  return qvals, repQvals

if __name__ == "__main__":
  main()
