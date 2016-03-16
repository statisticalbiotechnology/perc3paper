import xml.etree.ElementTree as ET
import sys
import matplotlib.pyplot as plt
import numpy as np
import collections
import scipy.stats as stats

ns = {'pout': 'http://per-colator.com/percolator_out/15'}
tree = ET.parse(sys.argv[1])
root = tree.getroot()

peptides = root.find('pout:peptides', ns)
proteinPvalues = collections.defaultdict(list)
for peptide in peptides:
  isDecoy = peptide.get('{http://per-colator.com/percolator_out/15}decoy') != "false"
  if isDecoy:
    count = 0
    for protein_element in peptide.iter("{http://per-colator.com/percolator_out/15}protein_id"):
      protein = protein_element.text
      count += 1
    if count == 1:
      pvalue = float(peptide.find('pout:p_value', ns).text)
      proteinPvalues[protein].append(pvalue)

pvalues = list()
numPeptides = list()
for protein in proteinPvalues:
  c = -2* sum(map(np.log, proteinPvalues[protein]))
  pvalue = 1 - stats.chi2.cdf(c, len(proteinPvalues[protein]) * 2)
  numPeptides.append(len(proteinPvalues[protein]))
  pvalues.append(pvalue)
  
plt.hist(pvalues, bins = 100)
plt.xlabel("p-value")
plt.ylabel("count")

plt.figure()
x = stats.probplot(pvalues, sparams = (0,1), dist="uniform", fit = False, plot=plt)

plt.figure()
plt.plot(x[0])

plt.figure()
plt.loglog(x[0], x[1], 'bx', label = 'probability plot')
plt.loglog([1e-5, 1], [1e-5, 1], 'k--', label = 'y = x')
plt.loglog([1e-5, 1], [2e-5, 2], 'k', ls = 'dotted', label = 'y = 2x')
plt.loglog([1e-5, 1], [0.5e-5, 0.5], 'k', ls = 'dotted', label = 'y = 0.5x')
plt.xlim([1e-5, 1])
plt.ylim([1e-5, 1])
plt.xlabel('quantiles')
plt.ylabel('ordered protein p-values')
plt.legend(loc = 'upper left')

plt.figure()
plt.hist(numPeptides, bins = range(50))
plt.xlabel("peptides per protein")
plt.ylabel("count")

plt.show()
      
