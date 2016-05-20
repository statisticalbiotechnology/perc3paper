import csv
import sys
import os
import subprocess

# runs the stand alone version of fido, but cannot get any remotely decent performance from it so far
def main(argv):
  fidoPath = "/home/matthew/software/fido/bin/FidoChooseParameters"
  inputFolder = "/media/storage/mergespec/data/103111-Yeast-2hr/percolator_tdc/tab"
  targetPeptFile = os.path.join(inputFolder, "103111-Yeast-2hr.percolator.tab.peptides")
  decoyPeptFile = os.path.join(inputFolder, "103111-Yeast-2hr.percolator.decoys.tab.peptides")
  
  qvalThresh = 0.10
  baseFN = os.path.join(inputFolder, "103111-Yeast-2hr.fido_q" + str(int(qvalThresh*100)))
  outputGraphFile = baseFN + "_graph.txt"
  outputProtListFile = baseFN + "_prot_list.txt"
  print outputGraphFile
  print outputProtListFile
  if not os.path.isfile(outputGraphFile) and not os.path.isfile(outputProtListFile):
    generateFidoInput(targetPeptFile, outputGraphFile, outputProtListFile, qvalThresh)
    generateFidoInput(decoyPeptFile, outputGraphFile, outputProtListFile, qvalThresh)
  
  fidoOutputFile = baseFN + "_output.txt"
  print fidoOutputFile
  if not os.path.isfile(fidoOutputFile):
    executeCmd("%s -c 1 %s %s > %s" % (fidoPath, outputGraphFile, outputProtListFile, fidoOutputFile))
    
  fidoPercOutputFile = baseFN + ".tab.proteins"
  print fidoPercOutputFile
  if not os.path.isfile(fidoPercOutputFile):
    generatePercolatorOutput(fidoOutputFile, fidoPercOutputFile)
  
  fidoGroupedOutputFile = baseFN + "_output_grouped.txt"
  print fidoGroupedOutputFile
  if not os.path.isfile(fidoGroupedOutputFile):
    executeCmd("%s -g -c 1 %s %s > %s" % (fidoPath, outputGraphFile, outputProtListFile, fidoGroupedOutputFile))
  
  fidoGroupedPercOutputFile = baseFN + "_grouped.tab.proteins"
  print fidoGroupedPercOutputFile
  if not os.path.isfile(fidoGroupedPercOutputFile):
    generatePercolatorOutput(fidoGroupedOutputFile, fidoGroupedPercOutputFile)
  
def generateFidoInput(peptideFile, outputGraphFile, outputProtListFile, qvalThresh):
  prots = set()
  reader = csv.reader(open(peptideFile, 'rb'), delimiter = '\t')
  reader.next()
  with open(outputGraphFile, 'a') as writer:
    for row in reader:
      qval = float(row[2])
      if qval > qvalThresh:
        break
      PEP = float(row[3])
      peptide = row[4]
      proteins = row[5:]
      
      writer.write("e " + peptide[2:-2] + '\n')
      for protein in proteins:
        writer.write("r " + protein + '\n')
        prots.add(protein)
      writer.write("p " + str(1.0 - PEP) + '\n')
   
  with open(outputProtListFile, 'a') as writer:
    protString = " , ".join(prots)
    writer.write("{ " + protString + " }\n")

def generatePercolatorOutput(fidoOutput, percolatorProtFile):
  writer = csv.writer(open(percolatorProtFile, 'wb'), delimiter = '\t')
  writer.writerow(["ProteinId", "ProteinGroupId", "q-value", "posterior_error_prob", "peptideIds"])
  proteinGroups = 1
  falseDiscoveries = 0.0
  with open(fidoOutput, 'r') as reader:
    for line in reader:
      fields = line[:-1].split(" ")
      proteinPost = float(fields[0])
      proteinPEP = 1.0 - proteinPost
      falseDiscoveries += proteinPEP
      proteins = list()
      for f in fields[1:]:
        if len(f) > 1:
          proteins.append(f)
      writer.writerow([",".join(proteins), proteinGroups, falseDiscoveries / proteinGroups, proteinPEP])
      proteinGroups += 1
  
def executeCmd(cmd, jobType = "local"):
  print(cmd)
  rc = subprocess.call(cmd, shell=True)
  #rc = 0
  if rc == 1:
    print("Error while processing " + cmd)
    return 1
  else:
    return 0    
    
if __name__ == "__main__":
  main(sys.argv)

