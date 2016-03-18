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
  percTabBase = "/home/matthew/mergespec/data/103111-Yeast-2hr/percolator_tdc/tab/103111-Yeast-2hr.percolator"
  
  force = False
  
  optionsArray = []
  
  options = {}
  options["method"] = "fisher" # fisher, twopept, bestpept, multPEP
  options["targetDecoyAnalysis"] = "picked" # picked, classic, pval
  options["removeSharedPeptides"] = False
  options["proteinGroupingThreshold"] = 0.05
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
  
  plt.suptitle("Mimic hm\_yeast", fontsize = 24, fontweight = 'bold')
  
  colors = dc.get_distinct_grey(len(optionsArray))
  
  for i, options in enumerate(optionsArray):
    targetOutFN, decoyOutFN = prot.getOutputFN(percTabBase, options)
    
    totalDbProteins = 6721*10
    
    fileName = targetOutFN
    print fileName
    
    if not os.path.isfile(fileName) or force:
      prot.writeProteinFdrs(percTabBase, options)
    
    qvals, fdrs, fomrs, tpfp = getQvalues(fileName, totalDbProteins)
    numSignificant = sum(1 if qval < 0.01 else 0 for qval in qvals)
    print "#significant proteins Reported 1% FDR:", numSignificant
    print "#significant proteins Observed 1% FDR:", sum(1 if qval < 0.01 else 0 for qval in fdrs)
    print "Observed FDR at reported 1% FDR:", fdrs[numSignificant-1]
    plotQvalues(qvals, fdrs, fomrs, tpfp, options, colors[i])
  
  plt.figure(2)
  x = np.linspace(0, 1, num=1000)
  plt.plot(x, x, 'k--')
  plt.legend(loc = 'lower right', prop={'size':24})
  setAxisFontSize(20)
  plt.tight_layout()
  
  plt.figure(4)
  plt.plot([0.01, 0.01], [0, 2000], 'k', linestyle = 'dotted')
  plt.xlim([0, 0.05])
  plt.legend(loc = 'lower right', prop={'size':24})
  setAxisFontSize(20)
  plt.tight_layout()
  
  plt.figure(5)
  x = np.linspace(-5, 1, num=1000)
  plt.plot(x, x, 'k--')
  plt.legend(loc = 'lower right', prop={'size':24})
  setAxisFontSize(20)
  plt.tight_layout()
  
  #plt.subplots_adjust(top=0.92)
  plt.show()

def setAxisFontSize(size):
  for tick in plt.gca().xaxis.get_major_ticks():
    tick.label.set_fontsize(size)
  for tick in plt.gca().yaxis.get_major_ticks():
    tick.label.set_fontsize(size)
    
def getQvalues(fileName, totalProteins):
  file = open(fileName, 'rb') # The input is the protein output file (.tab) from Percolator (-l)
  data = csv.reader(file, delimiter='\t')
  table = [row for row in data]
  x = np.linspace(0, 1, num=1000)
  fp = 0
  tp = 0

  qvals, fdrs, fomrs, tpfp = [], [], [], []

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
    
    tpfp.append([tp,fp])
    qvals.append(reportedQvalue)
    fdrs.append(fdr)
    fomrs.append(fomr)
  return qvals, fdrs, fomrs, tpfp

def plotQvalues(qvals, fdrs, fomrs, tpfp, options, color):    
  numRows = 2
  numCols = 2
  
  xlabel = "$FDR_{REV}$"
  ylabel = "Observed $FDR_{TRAP}$"
  labelFontSize = 30
  figIdx = 1
  
  #numIds = sum(1 for fdr in qvals if fdr < 0.01)
  #plt.suptitle("Mimic " + " (" + str(numIds) + " proteins at 1\% reported FDR)", fontsize = 24, fontweight = 'bold')
  
  #plt.subplot(numRows, numCols, figIdx)
  plt.figure(figIdx)
  x = np.linspace(0, 1, num=1000)
  plt.plot(qvals, fdrs, label = options["method"], linewidth = 2, color = color)
  plt.plot(x, x, 'k--')
  plt.axis([0, 1, 0, 1])
  plt.xlabel(xlabel, fontsize = labelFontSize)
  plt.ylabel(ylabel, fontsize = labelFontSize)
  
  figIdx += 1
  
  #plt.subplot(numRows, numCols, figIdx)
  plt.figure(figIdx)
  plt.plot(qvals, fdrs, label = options["method"], linewidth = 2, color = color)
  plt.axis([0, 0.1, 0, 0.1])
  plt.xlabel(xlabel, fontsize = labelFontSize)
  plt.ylabel(ylabel, fontsize = labelFontSize)
  
  figIdx += 1
  
  #plt.subplot(numRows, numCols, figIdx)
  plt.figure(figIdx)
  plt.plot(fomrs, fdrs, label = options["method"], linewidth = 2, color = color)
  plt.axis([0, 0.2, 0, 0.2])
  plt.xlabel("False Ommission Rate", fontsize = labelFontSize)
  plt.ylabel(ylabel, fontsize = labelFontSize)
  
  figIdx += 1
  
  #plt.subplot(numRows, numCols, figIdx)
  plt.figure(figIdx)
  plt.plot([x[1]/float(x[1]+x[0]) for x in tpfp], [x[0] for x in tpfp], label = options["method"], linewidth = 2, color = color)
  plt.ylim([0, 1200])
  plt.xlabel(ylabel, fontsize = labelFontSize)
  plt.ylabel("Sample proteins", fontsize = labelFontSize)
  
  figIdx += 1
  
  #plt.subplot(numRows, numCols, figIdx)
  plt.figure(figIdx)
  plt.plot([np.log10(qval + 1e-20) for qval in qvals], [np.log10(fdr+1e-20) for fdr in fdrs], label = options["method"], linewidth = 2, color = color)
  plt.axis([-3, -1, -3, -1])
  plt.xlabel("$\log_{10}$(" + xlabel + ")", fontsize = labelFontSize)
  plt.ylabel("$\log_{10}$(" + ylabel + ")", fontsize = labelFontSize)
   
if __name__ == "__main__":
  main()

