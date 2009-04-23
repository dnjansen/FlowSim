#!/bin/sh

OPT="-O 0:0000 -O 1:0001 -O 4:0010 -O 5:0011 -O 16:1000 -O 17:1001 -O 28:1110 -O 29:1111"

../benchmark -o labels1 -t random $OPT --avg 4 --fp-approx 1e-10 -d usertime --plot r labels1.rnd
../benchmark -o labels2 -t random $OPT --avg 4 --fp-approx 1e-10 -d usertime --plot r labels2.rnd
../benchmark -o labels3 -t random $OPT --avg 4 --fp-approx 1e-10 -d usertime --plot r labels3.rnd
../benchmark -o labels4 -t random $OPT --avg 4 --fp-approx 1e-10 -d usertime --plot r labels4.rnd
