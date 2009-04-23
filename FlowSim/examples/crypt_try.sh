#!/bin/sh

#1: state partition, 2: quotient 4 P-invariant, 8: significant arcs, 16: pmf  

AVG="--avg 1"
#OPT="-O 1:1 -O 2:2"
OPT="-O 0:None -O 1:Part -O 4:PInv -O 5:PInvPart  -O 16:PMF -O 17:PMFPart -O 28:PMFSigPar -O 29:all"

../src/benchmark --model-info --extra-rmap --precision 5 -t pa -o crypt $OPT $AVG --latex -d all crypt3.pa crypt4.pa crypt5.pa
