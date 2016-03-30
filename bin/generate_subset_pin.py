#!/usr/bin/python

'''
This file generates a subset of the PSMs in a percolator input file
'''

import sys
import random
import csv

def main(argv):
  inFile = "/media/storage/mergespec/data/Pandey/percolator_tdc/pin/Pandey.73M.tab"
  factor = 12
  numPsms = 73139505/(2**factor)
  
  print "numPsms =", numPsms
  
  outFolder = "/media/storage/mergespec/data/Pandey/percolator_tdc/pin/subsample" + str(factor)
  
  for seed in range(1,11):
    outFile = os.path.join(outFolder, "Pandey.ss" + str(factor) + ".seed" + seed + "tab")  
    random.seed(seed)
    
    reader = csv.reader(open(inFile,'rb'), delimiter='\t')
    writer = csv.writer(open(outFile,'wb'), delimiter='\t')
    writer.writerow(reader.next()) # copy header
    writer.writerow(reader.next()) # copy default direction
    
    lastScannr = -1
    writeQueue = list()
    for row in reader:
      if row[1] == "1":
        if int(row[2]) != lastScannr:
          writeQueue = list()
        if random.random() < frac:
          writer.writerow(row)
          row[1] = "-1"
          writeQueue.append(row[:4])
        lastScannr = int(row[2])
      elif row[:4] in writeQueue:
        writer.writerow(row)
        writeQueue.remove(row[:4])
  
if __name__ == "__main__":
   main(sys.argv[1:])
