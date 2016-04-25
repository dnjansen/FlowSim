#!/bin/sh

# This script is the basis of the performance comparison in the Information and
# Computation article by Zhang and Jansen, accepted for publication in 2016.

RMDP="../../src/r-mdp"	# location of the r-mdp binary
FLOWSIM="../../src/flowsim"   # location of the flowsim binary

OPTS="-O 49:Crafa/Ranzato -O 1:partitioning -O 2:space-efficient -O 18:time-efficient"
MODELS="--lext .label random.pa"
OUTPUT="--model-info --precision 5 --time-unit s --space-unit MiB --data memory --data usertime --data partitions"

if [ ! -f $FLOWSIM -a ! -x $FLOWSIM ]; then
  echo "Please edit this script and set the variable \$FLOWSIM to the location of the flowsim binary on your system."
  exit
fi

N=200	# number of states
MS=5	# maximum number of transitions per state

for (( D=5 ; D<30 ; D+=5 )) # in 5 10 15 20 25 30 35 40 45
do
	echo D = $D, MS = $MS
	date
	rm -f killme
	for (( i=0 ; i<25 ; i++ ))
	do
		$RMDP $N 2 $D 2 1 $MS > killme.pa
		$FLOWSIM $OUTPUT $OPTS --type pa --avg 1 killme.pa >> killme 2>&1
	done
	sed -f random2015.sed < killme > `date "+%Y%m%d-%H%M"`-random-D$D.csv
done

D=5	# maximum number of arrows per transition

for (( MS=10 ; MS<45 ; MS+=5 ))
do
	echo D = $D, MS = $MS
	date
	rm -f killme
	for (( i=0 ; i<25 ; i++ ))
	do
		$RMDP $N 2 $D 2 1 $MS > killme.pa
		$FLOWSIM $OUTPUT $OPTS --type pa --avg 1 killme.pa >> killme 2>&1
	done
	sed -f random2015.sed < killme > `date "+%Y%m%d-%H%M"`-random-MS$MS.csv
done
rm -f killme*
date
