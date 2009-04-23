#!/usr/local/bin/bash

# Random model description
# The string "VAR" is replaced by some iterated value
#      N   A B  C Fb Lb Cb Pb Sb
MODEL="VAR 1 1 1 0  0  0  0  0"

# These values will be substituted for VAR
VALUES="1 2 3 4 5"

# Benchmark settings
AVG="4"         # bench each model this many times
OPT="0 16"      # optimizations to benchmark; must specify as list of decimal integers (see manual/benchmark.htm)
LABELS="3"      # number of labels

# Data output
DESC="1 usertime
2 systemtime
3 realtime
4 partitions
5 iterations
6 intitial size
7 final size
8 num maxflow
9 num p-invariant failed
10 num sig-arc deleted
11 memory: relation map
12 memory: partition map
13 memory: relation (array of pairs)
14 memory: maxflow networks
15 memory: model
18 nets saved by parametric maxflow
19 times saved nets were reused"
PRINT="1 5 7"
WIDTH="8"      # for integers, the width of the entire field. for floats, width of non-fractional part
PRECISION="4"  # number of decimals after the point for float values
OUTDIR="examples/strongweak"

# Perform benchmark
for v in $VALUES; do
  PARAM="`echo "$MODEL" | sed -e "s/VAR/$v/g"`";
  echo -n "[" "$PARAM" "] ";
  ./r-dtmc $PARAM > "$OUTDIR/model_$v";
  for o in $OPT; do
    ./benchone_$o strong dtmc $LABELS -1 $AVG "$OUTDIR/model_$v" > "$OUTDIR/reports_$v""_$o";
    ./benchone_$o weak   dtmc $LABELS -1 $AVG "$OUTDIR/model_$v" > "$OUTDIR/reportw_$v""_$o";
  done;
#  rm "$OUTDIR/model_$v";
  echo "";
done

# Format results
rm -f "$OUTDIR/report"
for d in $PRINT; do
  if `test $d -le 3`; then
    RWIDTH="`expr $WIDTH + $PRECISION + 1`";
    RPRECISION="$PRECISION";
  else
    RWIDTH="$WIDTH";
    RPRECISION="0";
  fi
  
  echo "=================================================" >> "$OUTDIR/report";
  echo "$DESC" | egrep "^$d " | sed -e 's/^[0-9]+ //g' >> "$OUTDIR/report";
  echo "=================================================" >> "$OUTDIR/report";
  
  echo -n " OPT |" >> $OUTDIR/report;
  for o in $OPT; do
    printf " %*d |" "`expr \"(\" $RWIDTH \"*\" 2 \")\" + 1`" "$o" >> "$OUTDIR/report";
  done;
  echo "" >> "$OUTDIR/report";
  
  for v in $VALUES; do
    printf "%4d |" $v >> "$OUTDIR/report";
    for o in $OPT; do
      SR="`cat "$OUTDIR/reports_$v""_$o" | awk '{ print $'$d'; }'`";
      WR="`cat "$OUTDIR/reportw_$v""_$o" | awk '{ print $'$d'; }'`";
      printf " %*.*f/%*.*f |" $RWIDTH $RPRECISION "$SR" $RWIDTH $RPRECISION "$WR" >> "$OUTDIR/report";
    done;
    echo "" >> "$OUTDIR/report";
  done;
done
