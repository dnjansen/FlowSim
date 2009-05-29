#!/bin/sh

FLOWSIM="../../src/flowsim"

OPTS="-O 0:0000 -O 1:0001 -O 4:0010 -O 5:0011 -O 16:1000 -O 17:1001 -O 20:1010 -O 21:1011 -O 24:1100 -O 25:1101 -O 28:1110 -O 29:1111"
OUTPUT=" --name uniform --precision 4 --time-unit ms"

if [ ! -f $FLOWSIM -a ! -x $FLOWSIM ]; then
  echo "Please edit this script and set the variable \$FLOWSIM to the location of the flowsim binary on your system."
  exit
fi

$FLOWSIM $OUTPUT $OPTS --type random --data usertime --avg 2 uniform.rnd;
