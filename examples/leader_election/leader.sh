#!/bin/sh

FLOWSIM=""   # location of the flowsim binary

OPTS="-O 0:0000 -O 1:0001 -O 4:0010 -O 5:0011 -O 16:1000 -O 17:1001 -O 28:1110 -O 29:1111"
MODELS="leader3_6.model leader3_8.model leader3_10.model leader3_12.model"
OUTPUT="--model-info --extra-rmap --precision 3 --name leader --latex --quiet"

if `test -z "$FLOWSIM"`; then
  echo "Please edit this script and set the variable \$FLOWSIM to the location of the flowsim binary on your system."
  exit
fi

$FLOWSIM $OUTPUT $OPTS --type dtmc --data usertime --data memory -- $MODELS
