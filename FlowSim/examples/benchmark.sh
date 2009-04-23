#!/bin/sh

OPT="-O 0:A -O 16:B -O 1:C -O 3:D -O 19:E -O 28:F -O 31:G -O 7:H"
LEADER="leader3_3.model leader3_4.model leader3_5.model leader3_6.model leader3_8.model"
LEADER_LARGE="leader3_10.model leader3_12.model leader3_14.model"
AWK_COL='$2 " " $3'





# Create figure 2 (data source for plot over number of successors; optimzations state part & state+net part)
N="50"     # number of states
A="1"      # min number of successors per state

rm figure_2.data
for B in 1 2 3 4 5 6 7 8 9 10; do
  ../src/r-dtmc $N $A $B > rnd$N,$A,$B.model;
  ../src/benchmark -o "rnd" --plot n -t dtmc -d usertime -O 1 -O 3 -q --fp-approx 1e-15 rnd$N,$A,$B.model;
  cat rnd_usertime.data | awk "{ print $a \" \" $AWK_COL; }" >> figure_2.data;
  rm rnd$N,$A,$B.model rnd_usertime.data rnd_usertime.gnuplot;
done




# Create figure 3 (data source for plot over number of labels; optimizations state part & state+net part)
N="50"     # number of states
A="1"      # min number of successors per state
B="5"      # max number of successors per state

rm figure_3.data
../src/r-dtmc $N $A $B > rnd$N,$A,$B.model;
for L in 1 2 3 4 5 6 7 8 9 10; do
  ../src/benchmark -o "rnd" --plot n --labels $L -t dtmc -d usertime -O 1 -O 3 -q --fp-approx 1e-15 rnd$N,$A,$B.model;
  cat rnd_usertime.data | awk "{ print $L \" \" $AWK_COL; }" >> figure_3.data;
  rm rnd_usertime.data rnd_usertime.gnuplot;  
done
rm rnd$N,$A,$B.model
