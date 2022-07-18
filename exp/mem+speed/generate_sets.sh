#!/usr/bin/bash

main="Pandey.25M.tab"
sizes=("25000" "20000" "15000" "10000" "5000" "3000" "1000" "500" "300" "100" "50" "30" "10" "5" "3" "1")

cat $main | gawk 'BEGIN{OFS="\t"}($2==-1){$NF="DECOY_N"}{print}' > pre.$main

for s in ${sizes[*]}; do
    size=$((s*1000))
    cat pre.$main | head -n $size > data/set.${s}k.tab
    echo "Created set " data/set.${s}k.tab " of size:"
    wc -l  data/set.${s}k.tab
done
