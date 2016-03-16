#/bin/bash

subsetSize=5000000
subsetSizeFmt=$(numfmt --to=si $subsetSize)
percolatorBin=/home/matthewt/mergespec/bin/percolator
dataDir=/media/hdd/matthew/mergespec/data/Pandey/
percDir=percolator_tdc_swissprot_partial_digest
tabDir=tab_subset_scoring_$subsetSizeFmt
stdoutDir=stdout_subset_scoring_$subsetSizeFmt

mkdir -p $dataDir/$percDir/{$tabDir,$stdoutDir}

for i in {1..10}; do
  $percolatorBin $dataDir/$percDir/pin/Pandey.tab \
    -r $dataDir/$percDir/$tabDir/Pandey.percolator.seed$i.tab.peptides \
    -B $dataDir/$percDir/$tabDir/Pandey.percolator.seed$i.decoys.tab.peptides \
    -m $dataDir/$percDir/$tabDir/Pandey.percolator.seed$i.tab.psms \
    -M $dataDir/$percDir/$tabDir/Pandey.percolator.seed$i.decoys.tab.psms \
    -N $subsetSize -S $i > \
      $dataDir/$percDir/$stdoutDir/Pandey.seed$i.stdout.txt 2>&1
done
