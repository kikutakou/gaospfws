#!/bin/bash


trap 'exit 1' INT

case $(hostname | grep -o -E "[0-9]+" | head -n 1) in
1) N=10 ; GOAL=(0.037037 0.02381 0.023148 0.02951375 0.0505135 0.03191500 0.04545500 0.022321 0.0217235 0.0303030) ;;
2) N=15 ; GOAL=(0.029915 0.03750 0.040000 0.03472390 0.0476190 0.03125985 0.05382235 0.052632 0.0375881 0.0287148) ;;
esac

bash ../clean.sh >/dev/null

date "+start  %m/%d %H:%M"
for S in 0 1 2 3 4 5 6 7 8 9 ; do
for i in $(seq 1 10) ; do
for T in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ; do

#bash $(dirname $0)/run.sh $N 10 $S -DA=300 -DGOAL=${GOAL[$S]}
bash $(dirname $0)/run.sh $N 10 $S --thread=$T -DGENMAX=1000
#wait

done
done
done > out.txt
date "+finish %m/%d %H:%M"

cat out.txt


