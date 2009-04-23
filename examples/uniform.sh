#!/bin/sh


#the configurations 8,9,12,13 are not there, so take it out first
#OPT="-O 0:0000 -O 1:0001 -O 4:0010 -O 5:0011 -O 8:0100 -O 9:0101 -O 12:0110 -O 13:0111 -O 16:1000 -O 17:1001 -O 20:1010 -O 21:1011 -O 24:1100 -O 25:1101 -O 28:1110 -O 29:1111"
OPT="-O 0:0000 -O 1:0001 -O 4:0010 -O 5:0011 -O 16:1000 -O 17:1001 -O 20:1010 -O 21:1011 -O 24:1100 -O 25:1101 -O 28:1110 -O 29:1111"
AWK_COL='$2 " " $3 " " $4 " " $5 " " $6 " " $7 " " $8 " " $9 " " $10 " " $11 " " $12 " " $13 " " $14 " " $15 " " $16 " " $17'


N="400"     # number of states
rm -f uniform_usertime.data uniform_memory.data uniform_maxflow.data
let A=1
let B=5
while [ $B -le 90 ]
do
  echo "generate models ...A="$A, "B="$B
  ../r-dtmc $N $A $B > rnd$N,$A,$B.model;
  echo "benchmarking ... B="$B
  ../src/benchmark -o "rnd"$N,$A,$B --plot n -t dtmc -d all $OPT --avg 2 --time-unit m -q --fp-approx 1e-5 rnd$N,$A,$B.model;
  echo "insert plot data..."
  cat rnd$N,$A,$B\_usertime.data | awk "{ print $B \" \" $AWK_COL; }" >> uniform_usertime.data;
  cat rnd$N,$A,$B\_memory.data | awk "{ print $B \" \" $AWK_COL; }" >> uniform_memory.data;
  cat rnd$N,$A,$B\_maxflow.data | awk "{ print $B \" \" $AWK_COL; }" >> uniform_maxflow.data;
  cat rnd$N,$A,$B\_iterations.data | awk "{ print $B \" \" $AWK_COL; }" >> uniform_iterations.data;
  let B=$B+5
done

echo "plotting ..."
gnuplot uniform_usertime_1.gnuplot 
gnuplot uniform_usertime_2.gnuplot 
gnuplot uniform_memory.gnuplot
gnuplot uniform_maxflow.gnuplot 
