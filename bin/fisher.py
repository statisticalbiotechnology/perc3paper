import numpy as np
import scipy.stats as stats
import bisect
import csv

def main():
  #pvalues = np.random.uniform(0,1,1000) #[0.11, 0.30]
  pvalues = [1e-12, 1e-5]
  pvalue = getFisherPvalue(pvalues)
  print "Fisher:", pvalue
  
  ks, p = stats.kstest(pvalues, 'uniform', alternative = 'greater')
  print "Kolm-Smir:", p

def getFisherPvalue(pvalues):
  c = -2* sum(map(np.log, pvalues))
  pvalue = 1 - stats.chi2.cdf(c, len(pvalues) * 2)
  return pvalue

def getFisherPvalueFromScores(scores, decoyScores):
  c = -2* sum(map(lambda x : np.log(getPvalue(x, decoyScores)), scores))
  pvalue = 1 - stats.chi2.cdf(c, len(scores) * 2)
  return pvalue
  
def getPvalue(score, decoyScores):
  return 1.0 - (bisect.bisect_left(decoyScores, float(score)) + 0.5) / (len(decoyScores)+1)

def getPercScores(peptFile):
  reader = csv.reader(open(peptFile, 'r'), delimiter = '\t')
  reader.next()
  scores = list()
  for row in reader:
    scores.append(float(row[1]))
  return list(reversed(scores))

def getPvalues(peptFile, decoyScores):
  reader = csv.reader(open(peptFile, 'r'), delimiter = '\t')
  reader.next()
  pvalues, scores = list(), list()
  for row in reader:
    scores.append(float(row[1]))
    pvalues.append(getPvalue(float(row[1]), decoyScores))
  return list(reversed(scores)), list(reversed(pvalues))
  
def getPvalues(decoyFN, targetFN):
  decoyScores = getPercScores(decoyFN)
  targetScores, targetPvalues = getPvalues(targetFN, decoyScores)
  decoyScores, decoyPvalues = getPvalues(decoyFN, decoyScores)
  return targetPvalues, decoyPvalues
