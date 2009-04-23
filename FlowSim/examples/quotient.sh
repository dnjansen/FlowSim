#!/bin/sh

#1: state partition, 2: quotient 4 P-invariant, 8: significant arcs, 16: pmf  

AVG="--avg 1"
#OPT="-O 1:1 -O 2:2"
#OPT="-O 0:None -O 1:Part -O 2:Quotient -O 18:QuotientPMF  -O 22:QuotientPInvariant -O 26:QuotientSig"
OPT="-O 0:Original -O 1:Part -O 2:Quotient -O 18:QuotientPMF"

../benchmark --model-info --precision 5 --labels 1 -t pa  $OPT $AVG --latex -d all crypt3.pa crypt4.pa  > quo_1.txt
../benchmark --model-info --precision 5 --labels 2 -t pa  $OPT $AVG --latex -d all crypt3.pa crypt4.pa  > quo_2.txt
../benchmark --model-info --precision 5 --time s --labels 3 -t pa  $OPT $AVG --latex -d all crypt3.pa crypt4.pa > quo_3.txt




# #!/bin/sh

# #1: state partition, 2: quotient 4 P-invariant, 8: significant arcs, 16: pmf  

# AVG="--avg 1"
# #OPT="-O 1:1 -O 2:2"
# #OPT="-O 0:None -O 1:Part -O 2:Quotient -O 18:QuotientPMF  -O 22:QuotientPInvariant -O 26:QuotientSig"
# OPT="-O 0:Original -O 2:Quotient -O 18:QuotientPMF"

# ../benchmark --model-info --precision 5 --labels 1 -t pa  $OPT $AVG --latex -d all crypt3.pa crypt4.pa crypt5.pa  > quo_1.txt
# ../benchmark --model-info --precision 5 --labels 2 -t pa  $OPT $AVG --latex -d all crypt3.pa crypt4.pa crypt5.pa  > quo_2.txt
# ../benchmark --model-info --precision 5 --time s --labels 3 -t pa  $OPT $AVG --latex -d all crypt3.pa crypt4.pa crypt5.pa  > quo_3.txt


# #MODELS="stabilising0.pa stabilising1.pa stabilising2.pa stabilising3.pa stabilising4.pa stabilising5.pa"
# #MODELS="mutual3.model"
# #MODELS="rabin3.pa"
# #MODELS="beauquier5.pa beauquier7.pa beauquier9.pa "

# #for model in $MODELS; do
# #../benchmark --model-info --precision 5 --labels 1 -t pa  $OPT $AVG --latex -d all $model 
# #done


