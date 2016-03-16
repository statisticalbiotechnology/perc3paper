import sys
import matplotlib.pyplot as plt
import numpy as np
import csv

#peptideFile = "/home/matthew/mergespec/data/103111-Yeast-2hr/percolator_tdc/tab/103111-Yeast-2hr.percolator.tab.peptides"
#peptideFile = "/media/storage/mergespec/data/Linfeng/percolator_tdc/tab/Linfeng.percolator.tab.peptides"
#peptideFile = "/home/matthew/mergespec/data/103111-Yeast-2hr/percolator_tdc_full_digest_n2/tab/103111-Yeast-2hr.percolator.tab.peptides"
peptideFile = "/home/matthew/mergespec/data/103111-Yeast-2hr/percolator_tdc/tab/103111-Yeast-2hr.percolator.tab.peptides"
reader = csv.reader(open(peptideFile, 'rb'), delimiter = '\t')

reader.next()

numProts = []
for row in reader:
  numProts.append(len(row)-5)

plt.hist(numProts, bins = map(lambda x : x - 0.5, range(0,10)), normed = True)
plt.xlim([0,10])
plt.show()
