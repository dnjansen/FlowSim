#!/bin/sh

FLOWSIM="../../src/flowsim"   # location of the flowsim binary

OPTS="-O 4:0010 -O 5:0011 -O 28:1110 -O 29:1111"

if [ ! -f $FLOWSIM -a ! -x $FLOWSIM ]; then
  echo "Please edit this script and set the variable \$FLOWSIM to the location of the flowsim binary on your system."
  exit
fi

$FLOWSIM -o labels1 -t random $OPTS --avg 4 --fp-approx 1e-10 -d usertime --plot r labels1.rnd
$FLOWSIM -o labels2 -t random $OPTS --avg 4 --fp-approx 1e-10 -d usertime --plot r labels2.rnd
