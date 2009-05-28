#!/bin/sh

FLOWSIM=""   # location of the flowsim binary

OPTS="-O 0:None -O 1:Part -O 4:PInv -O 5:PInvPart"
MODELS="crypt3.pa crypt4.pa crypt5.pa"
OUTPUT="--model-info --extra-rmap --precision 5 --name dcrypt --latex --quiet"

if `test -z "$FLOWSIM"`; then
  echo "Please edit this script and set the variable \$FLOWSIM to the location of the flowsim binary on your system."
  exit
fi

$FLOWSIM $OUTPUT $OPTS --type pa --data usertime --data memory -- $MODELS
