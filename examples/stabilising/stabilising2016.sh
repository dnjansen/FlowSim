#!/bin/sh

# This script is the basis of the performance comparison in the Information and
# Computation article by Zhang and Jansen, accepted for publication in 2016.

FLOWSIM="../../src/flowsim"   # location of the flowsim binary

OPTS="-O 49:Crafa/Ranzato -O 1:partitioning -O 2:space-efficient -O 18:time-efficient"
MODELS="--lext .label stabilising1.pa stabilising0.pa"
OUTPUT="--model-info --precision 5 --latex --time-unit s --space-unit MiB --data memory --data usertime --data partitions --data finalsize"

if [ ! -f $FLOWSIM -a ! -x $FLOWSIM ]; then
  echo "Please edit this script and set the variable \$FLOWSIM to the location of the flowsim binary on your system."
  exit
fi

date
$FLOWSIM $OUTPUT --name `date "+%Y%m%d-%H%M"`-stabilising $OPTS --type pa --avg 4 $MODELS

OPTS="-O 1:partitioning -O 2:space-efficient -O 18:time-efficient"
MODELS="--lext .label.orig stabilising2.pa stabilising3.pa"

date
$FLOWSIM $OUTPUT --name `date "+%Y%m%d-%H%M"`-stabilising $OPTS --type pa --avg 4 $MODELS

OPTS="-O 49:Crafa/Ranzato"

date
$FLOWSIM $OUTPUT --name `date "+%Y%m%d-%H%M"`-stabilising $OPTS --type pa --avg 4 $MODELS
date
