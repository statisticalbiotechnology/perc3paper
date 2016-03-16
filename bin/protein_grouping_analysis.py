import sys
import operator as op
import argparse
import numpy as np
import networkx as nx
import itertools as it
import time
import matplotlib.pyplot as plt

# method for reading the given fasta database
def readFasta(fp):
  name, seq = None, []
  for line in fp:
    line = line.rstrip()
    if line.startswith(">"):
      if name: yield (name, "".join(seq))
      name, seq = line[1:].split(" ")[0], []
    else: seq.append(line)
  if name: yield (name, "".join(seq))

# method for generating the list of peptides
def getPeptideByTrypticDigest(seq, minLength = 6, maxLength = 50):

  lenS, start = len(seq), 0

  for i in range(lenS):
    if (seq[i] == 'K' or seq[i] == 'R') and seq[min(lenS-1,i+1)] != 'P':
      lenP = i - start + 1 
      if lenP >= minLength and lenP <= maxLength : yield (seq[start : i + 1])
      start = i + 1
  lenP = lenS - start
  if lenP >= minLength and lenP <= maxLength : yield (seq[start : ])

def getDecoyPeptide(peptides):
  for peptide in peptides:
    yield peptide
    
def insertToGraph(g, allProteins, peptideSet, proteins, prefix, peptideLabel, lookUp, getPeptide = getPeptideByTrypticDigest, decoy = None):

  # get peptide lists for each protein and build the graph
  for protein in proteins:
    proteinName = prefix + protein
    allProteins.add(proteinName)  
    g.add_node(proteinName, bipartite = 0, pep = None)
    
    dPeptides = set()
    for peptide in getPeptide(lookUp[protein]):
      if not peptide in g: 
        g.add_node(peptide, bipartite = 1, pep = None, label = peptideLabel) # t for target, d for decoy
      g.add_edge(proteinName, peptide)
      
      peptideSet.add(peptide)
      if decoy is not None: dPeptides.add(peptide[::-1])
      
    if decoy is not None: decoy[protein] = dPeptides 
    
#################################################### main method 
def main(**kwargs): 
  
  start_time = time.clock()  
  argList = kwargs.get('args', None)
  #print argList
  
  docu = '''A script for simulating the outcome of a mass spec experiment.  
  The script takes a sequence database as input and randomly assignes a 
  fraction of them as present and the rest of the proteins as
  absent. The proteins constituent peptides are subsequently simulated
  as being matched against a set of (virtual) spectra. In the process
  some of the peptides are simulated as being correctly and some
  incorrectly matched. peptides from absent proteins are simulated as
  allways incorrectly matched, while correctly matched peptides are all
  simulated as steming from present poteins.'''

  parser = argparse.ArgumentParser(description = docu)
  parser.add_argument('fastaDatabase',
                   help = 'The sequence database which we have matched with our simulated search engine')
  parser.add_argument('--decoy', default = None,  metavar = 'PATH',
                   help = 'Output the manifactured decoy database as a fasta file named by PATH')
  parser.add_argument('--outputPath', default = 'simulatorOutput.txt',  metavar = 'OUTPATH',
                   help = 'File path where the simulator output will be saved')
  
  if argList: args = parser.parse_args(argList) 
  #  print args
  else: args = parser.parse_args() # parse sys.argv, if not otherwise specified

  ####################### skip this step if a dict is received as an optional argument 

  allProteinSeqs = kwargs.get('allProteinSeqs', None)
  if allProteinSeqs is None:
    #print "reading fasta database"
    allProteinSeqs = {}
    with open(args.fastaDatabase) as f:
      for proteinName, proteinSeq in readFasta(f):
        allProteinSeqs[proteinName] = proteinSeq  
    #print "file read",
    #print time.clock() - start_time, "seconds"

  allProteins = set(allProteinSeqs)
  numProteins = len(allProteins)

  #print "present/absent sampled", 
  #print time.clock() - start_time, "seconds"

  #initialize the protein-peptide bipartite graph and other sets  
  graph, allProteinsResult, decoyProteins = nx.Graph(), set(), {}
  allPeptides, decoyPeptides = set(), set()

  # insert the proteins and their respective peptides into the graph
  insertToGraph(graph, allProteinsResult, allPeptides, allProteins, "target_", 't', lookUp = allProteinSeqs, decoy = decoyProteins)
  #insertToGraph(graph, allProteins, decoyPeptides, decoyProteins, "decoy_", 'd', lookUp = decoyProteins, getPeptide = getDecoyPeptide, decoy = None)
  
  print "graph is built", 
  print time.clock() - start_time, "seconds"
  print ""
  
  graphs = list(nx.connected_component_subgraphs(graph))
  duplicateProteins, fragmentProteins = 0, 0
  numProteinsWithUniquePeptides, numProteinGroups, numProteinGroupsWithUniquePeptides, numPeptidesAfterCollapsing = 0, 0, 0, 0
  numProteinsDist, numPeptidesDist, numProtPeptDist, numEdgesDist = list(), list(), list(), list()
  peptsPerProtDistBefore, peptsPerProtDistAfterFiltering, peptsPerProtDistAfterCollapsing, peptsPerProtDistAfterCollapsingAndFiltering = list(), list(), list(), list()
  for cc in graphs:
    numProteins, numPeptides = 0, 0
    peptideSets, uniquePeptides = list(), set()
    for n in cc.nodes_iter(data=False):
      if n.startswith("target_"):
        numProteins += 1
        nbs = set(cc.neighbors(n))
        peptideSets.append((nbs, n))
        peptsPerProtDistBefore.append(len(nbs))
      else:
        if len(cc.neighbors(n)) == 1:
          uniquePeptides.add(n)
        numPeptides += 1
    numProteinsDist.append(numProteins)
    numPeptidesDist.append(numPeptides)
    numProtPeptDist.append(numProteins + numPeptides)
    numEdgesDist.append(cc.number_of_edges())      
    
    removeList = list()
    for i in range(len(peptideSets)):
      numUniquePeptides = len(peptideSets[i][0] & uniquePeptides)
      if numUniquePeptides >= 1:
        numProteinsWithUniquePeptides += 1
        peptsPerProtDistAfterFiltering.append(numUniquePeptides)
      
      isDuplicateProtein, isFragmentProtein = False, False
      for j in range(len(peptideSets)):
        if i != j:
          psi = peptideSets[i][0]
          psj = peptideSets[j][0]
          if psi == psj:
            isDuplicateProtein = True
            removeList.append(min([peptideSets[i][1],peptideSets[j][1]]))
            #print "d", psi, psj
          elif len(psi) < len(psj) and psi <= psj:
            isFragmentProtein = True
            removeList.append(peptideSets[i][1])
            #print "f1", psi, psj
      if isFragmentProtein:
        fragmentProteins += 1
      elif isDuplicateProtein:
        duplicateProteins += 1
    
    # collapse duplicate/fragment proteins
    for n in set(removeList):
      cc.remove_node(n)
    peptideSets, uniquePeptides = list(), set()
    
    for n in cc.nodes_iter(data=False):
      if n.startswith("target_"):
        numProteinGroups += 1
        nbs = set(cc.neighbors(n))
        peptideSets.append(nbs)
        peptsPerProtDistAfterCollapsing.append(len(nbs))
      else:
        if len(cc.neighbors(n)) == 1:
          uniquePeptides.add(n)
        numPeptidesAfterCollapsing += 1
    
    for i in range(len(peptideSets)):
      numUniquePeptides = len(peptideSets[i] & uniquePeptides)
      if numUniquePeptides >= 1:
        numProteinGroupsWithUniquePeptides += 1
        peptsPerProtDistAfterCollapsingAndFiltering.append(numUniquePeptides)
        
  
  print "1a. #unique peptides:", len(allPeptides)
  print "1b. #proteins:", len(allProteins)
  print "2. #connected components:", len(graphs)
  
  plotHist(numProteinsDist, "3a. #proteins in connected component", 1, 50, 1)
  plotHist(numPeptidesDist, "3b. #peptides in connected component", 5, 500, 2)
  plotHist(numProtPeptDist, "3c. #peptides+proteins in connected component", 5, 500, 3)
  plotHist(numEdgesDist, "3d. #edges in connected component", 5, 500, 4)
  
  print "4. #proteins fully contained in another protein:", fragmentProteins
  print "5. #proteins with >= 1 other protein with same peptides:", duplicateProteins
  print "6. #proteins with >= 1 unique peptide:", numProteinsWithUniquePeptides
  
  plt.figure()
  plotHist(peptsPerProtDistBefore, "7a. Peptides per protein before", 5, 500, 1)
  plotHist(peptsPerProtDistAfterFiltering, "7b. Peptides per protein after removing non-unique peptides", 5, 500, 2)
  plotHist(peptsPerProtDistAfterCollapsing, "7c. Peptides per protein after collapsing duplicates/fragments", 5, 500, 3)
  plotHist(peptsPerProtDistAfterCollapsingAndFiltering, "7d. Peptides per protein after collapsing duplicates/fragments\nand removing non-unique peptides", 5, 500, 4)
  
  print "8a. #proteins after collapsing duplicates/fragments:", numProteinGroups
  print "8b. #proteins after collapsing duplicates/fragments with >= 1 unique peptide:", numProteinGroupsWithUniquePeptides
  
  plt.show()
  
def plotHist(dist, title, binSize, maxBin, idx):
  plt.subplot(2,2,idx)
  plt.hist(dist, bins = np.arange(0.5,min([maxBin, max(dist)])+0.5, binSize))
  plt.title(title)
  
#############################################################################################################
if __name__ == "__main__":
    main()
