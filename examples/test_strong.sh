#!/bin/sh
 
../src/dtmctest.exe --bench strong 5 aa 
../src/dtmctest_par.exe --bench strong 5 aa_par
../src/dtmctest_inv.exe --bench strong 5 aa_inv
../src/dtmctest_pmf.exe --bench strong 5 aa_pmf
../src/dtmctest_sig.exe --bench strong 5 aa_sig

cat aa |fgrep 'user time' > aa.csv
cat aa_par |fgrep 'user time' > aa_par.csv
cat aa_inv |fgrep 'user time' > aa_inv.csv
cat aa_pmf |fgrep 'user time' > aa_pmf.csv
cat aa_sig |fgrep 'user time' > aa_sig.csv