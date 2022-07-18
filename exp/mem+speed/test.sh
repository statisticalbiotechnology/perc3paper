#!/usr/bin/bash

percolators=("36" "33" "30" "28")
iterations=("1" "2" "3" "4" "5")

readarray -d '' entries < <(printf '%s\0' data/set.*.tab | sort -zV)
for psms in "${entries[@]}"; do 
    a=${psms#data/set.}
    size=${a%.tab}
    echo -n "$size "
    for ver in ${percolators[*]}; do
        echo -n "$ver "
        percolator=$ver/usr/bin/percolator
        for iter in ${iterations[*]}; do
            log=res/per.$ver.$size.$iter.log
            tf=res/per.$ver.$size.$iter.time
            /usr/bin/time -o $tf $percolator -P DECOY $psms 1>/dev/null 2> $log
        done
    done
    echo ""
done
