#!/bin/sh

FLOWSIM="../../src/flowsim"   # location of the flowsim binary

OPTS_P="-O 1:0001 -O 5:0011 -O 21:1011 -O 25:1101"
OPTS_NP="-O 0:0000 -O 4:0010 -O 20:1010 -O 24:1100"
OUTPUT="--plot r --time-unit ms"

if [ ! -f $FLOWSIM -a ! -x $FLOWSIM ]; then
  echo "Please edit this script and set the variable \$FLOWSIM to the location of the flowsim binary on your system."
  exit
fi

$FLOWSIM --name pbias_p --type random $OUTPUT $OPTS_P --avg 4 --fp-approx 1e-10 --data usertime --data maxflow pbias.rnd
$FLOWSIM --name pbias_np --type random $OUTPUT $OPTS_NP --avg 4 --fp-approx 1e-10 --data usertime --data maxflow pbias.rnd
