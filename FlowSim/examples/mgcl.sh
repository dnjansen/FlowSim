#!/bin/sh

OPT="-O 0:0000 -O 1:0001 -O 4:0010 -O 5:0011 -O 16:1000 -O 17:1001 -O 28:1110 -O 29:1111"
AVG="--avg 10"
TYPE="-o mgcl -d all --space-unit kb --precision 5 --time-unit ms"

../benchmark $OPT -t ctmc $TYPE $AVG --plot n --model-info --extra-rmap --latex --quiet mgcl10.ctmc mgcl25.ctmc mgcl50.ctmc mgcl75.ctmc mgcl100.ctmc mgcl125.ctmc mgcl150.ctmc
