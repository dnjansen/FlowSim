#!/bin/sh

FLOWSIM="../../src/flowsim"   # location of the flowsim binary

OPTS="-O 0:0000 -O 1:0001 -O 4:0010 -O 5:0011 -O 16:1000 -O 17:1001 -O 28:1110 -O 29:1111"
OUTPUT="--name fanout --plot r"

if [ ! -f $FLOWSIM -a ! -x $FLOWSIM ]; then
  echo "Please edit this script and set the variable \$FLOWSIM to the location of the flowsim binary on your system."
  exit
fi

$FLOWSIM --type random $OUTPUT $OPTS --avg 4 --fp-approx 1e-10 --data usertime --data maxflow fanout.rnd
