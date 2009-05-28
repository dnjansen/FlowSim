#!/bin/sh

FLOWSIM=""   # location of the flowsim binary

OPTS="-O 0:0000 -O 1:0001 -O 4:0010 -O 5:0011 -O 16:1000 -O 17:1001 -O 28:1110 -O 29:1111"
MODELS="mgcl50.ctmc mgcl75.ctmc mgcl100.ctmc mgcl125.ctmc mgcl150.ctmc"
OUTPUT=" --name mgcl --latex --model-info --precision 4 --quiet"

if `test -z "$FLOWSIM"`; then
  echo "Please edit this script and set the variable \$FLOWSIM to the location of the flowsim binary on your system."
  exit
fi

$FLOWSIM $OPTS $OUTPUT --type ctmc --data usertime --data memory -- $MODELS
