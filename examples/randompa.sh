#!/bin/sh
 
N="200"        # number of states in model
MINSUC="2"     # minimum number of successors per distribution
MAXSUC="5 7 9 11 13 15 17 19 21"   # maximum number of successors per distribution

ACT="2"    # list of different numbers of actions
MINA="1"       # minimum number of actions per state
MAXA="5"       # maximum number of actions per state
USIM="0.0"     # similarity of actions of one state: 0.0=random, 1.0=all actions the same

# benchmark options
#OPT="-O 0:Space -O 1:Time -O 2:Quotient -O 18:QuoPMF"
OPT="-O 0:Space -O 2:Quotient -O 18:QuoPMF"
DATA="--data all"
LABELS="1"



X="1 2 3 4 5 6 7 8"

for x in $X; do
MODELS=""
# create models
for n in $N; do
 for na in $ACT; do
   for ma in $MAXA; do
     for maxsuc in $MAXSUC; do
        ../src/r-mdp $n $MINSUC $maxsuc $na $MINA $ma $USIM > rpa_$n.$na.$ma.$maxsuc;
        MODELS="$MODELS rpa_$n.$na.$ma.$maxsuc";
     done;
   done;
 done
done

 
../src/benchmark --type pa --avg 1 --time-unit s --space-unit kB --labels $LABELS $DATA --model-info --latex --name rpa $OPT --precision 5  $MODELS > rp1_$x.txt


done