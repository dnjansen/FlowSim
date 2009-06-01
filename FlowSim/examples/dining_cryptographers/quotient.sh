#!/bin/sh

# This script compares the performance of the quotient automaton approach
# with and without parametric maximum flow, to the performance of the
# refinement loop approach with and without state partitioning.

FLOWSIM="../../src/flowsim"   # location of the flowsim binary

MODELS="crypt3.pa crypt4.pa crypt5.pa"
OPTS="-O 0:Refinement -O 1:StatePart -O 2:Quotient -O 18:QuotientPMF"
OUTPUT="--name quotient --precision 4 --latex"

if [ ! -f $FLOWSIM -a ! -x $FLOWSIM ]; then
  echo "Please edit this script and set the variable \$FLOWSIM to the location of the flowsim binary on your system."
  exit
fi

$FLOWSIM $OUTPUT $OPTS --type pa --avg 4 --data usertime --data memory --data maxflow -- $MODELS
