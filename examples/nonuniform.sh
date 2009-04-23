#!/bin/sh

../benchmark --plot r+ -t random --opt 0:0000 --opt 4:0010 nonuniform.rnd
#gnuplot benchmark_usertime_maxsuc_pbias.gnuplot
#gv benchmark_usertime_maxsuc_pbias.eps
#../benchmark --name bias --plot r+ -t random --time-units s --data usertime --data memory --opt 1 bias.rnd
#../benchmark --name bias --plot r+ -t random --time-units s --data usertime --data memory --opt 17 bias.rnd

