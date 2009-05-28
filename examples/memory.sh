#!/bin/sh

OPT="-O 0:A -O 16:B -O 1:C -O 3:D -O 19:E -O 28:F -O 31:G -O 7:H"
LEADER="leader3_3.model leader3_4.model leader3_5.model leader3_6.model leader3_8.model"
LEADER_LARGE="leader3_10.model leader3_12.model leader3_14.model"
#AWK_COL='$2 " " $3 " " $4 " " $5 " " $6 " " $7 " " $8 " " $9'
AWK_COL='$2 " " $3 " " $4 " " $5 " " $6 " " $7 " " $8 " " $9 " " $10 " " $11 " " $12 " " $13 " " $14 " " $15 " " $16 " " $17 " " $18 " " $19 " " $20 " " $21 " " $22 " " $23 " " $24 " " $25'


#N="80"     # number of states
#A="40"      # min number of successors per state
#rm -f figure_2.data
#for B in 51 52 53 54 55 56 57 58 59 60; do
##  ../src/r-dtmc $N $A $B > rnd$N,$A,$B.model;
#  ../src/benchmark -o "rnd" --plot n -t dtmc -d usertime --opt all -q --fp-approx 1e-15 rnd$N,$A,$B.model;
#  cat rnd_usertime.data | awk "{ print $B \" \" $AWK_COL; }" >> figure_2.data;
##  rm rnd$N,$A,$B.model rnd_usertime.data rnd_usertime.gnuplot;
#done

N="80"     # number of states
A="40"      # min number of successors per state
rm -f figure_2.data
# for B in 51 52 53 54 55 56 57 58 59 60; do
let B=51
while [ $B -le 54 ]
do
#  ../src/r-dtmc $N $A $B > rnd$N,$A,$B.model;
  ../src/benchmark -o "rnd" --plot n -t dtmc -d memory --opt all --avg 2 -q --fp-approx 1e-15 rnd$N,$A,$B.model;
  cat rnd_usertime.data | awk "{ print $B \" \" $AWK_COL; }" >> figure_2.data;
#  rm rnd$N,$A,$B.model rnd_usertime.data rnd_usertime.gnuplot;
  let B=$B+1
done

gnuplot uniform.gnuplot 
