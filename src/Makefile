#*****************************************************************************
#!
#    Copyright 2009-2014 Jonathan Bogdoll, Holger Hermanns, Lijun Zhang,
#                        David N. Jansen
#
#    This file is part of FlowSim.

#    FlowSim is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    FlowSim is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with FlowSim.  If not, see <http://www.gnu.org/licenses/>.
#
#*****************************************************************************



CXX=g++
#CXXFLAGS=-O0 -ggdb -Wpointer-arith -Winline -Wall -W -pedantic -DDEBUG
CXXFLAGS=-O3 -Wpointer-arith -Winline -Wall -W -DDEBUG -Wextra -Werror -ansi -pedantic -fno-implicit-templates -isystem ../../boost_1_60_0

.cc.o:
	$(CXX) $(CXXFLAGS) -c $<

all: r-dtmc r-mdp flowsim benchsingle modelstat

.PHONY: all depend clean clean-all benchsingle

depend:
	sed '/^# DO NOT DELETE/q' < Makefile > Makefile.$$
	echo >> Makefile.$$
	for i in benchone.cc compactmaxflow.cc StrongMC.cc StrongPA.cc StrongQ.cc ; do \
	  $(CXX) -MM $(CXXFLAGS) -MF killme.$$ -MQ "$$i" "$$i" || exit 1 ; \
	  sed '1s/ '"$$i"'//' < killme.$$ >> Makefile.$$ ; \
	done
	for i in [^bcS]*.cc bench.cc bench[^o]*.cc ; do \
	  $(CXX) -MM $(CXXFLAGS) -MF killme.$$ "$$i" || exit 1 ; \
	  sed '1s/ '"$$i"'//' < killme.$$ >> Makefile.$$ ; \
	done
	rm killme.$$
	mv Makefile.$$ Makefile

clean:
	rm -f *~ *.o *.core

clean-all: clean
	rm -f r-dtmc r-mdp flowsim benchone_* modelstat

r-dtmc: r-dtmc.cc prmodel.o
	$(CXX) $(CXXFLAGS) -o r-dtmc -D_STANDALONE_ r-dtmc.cc prmodel.o

r-mdp: r-mdp.cc prmodel.o
	$(CXX) $(CXXFLAGS) -o r-mdp -D_STANDALONE_ r-mdp.cc prmodel.o

modelstat: modelstat.cc prmodel.o
	$(CXX) $(CXXFLAGS) -o $@ $^

flowsim: bench.o benchmark.o prmodel.o compactmaxflow.cc r-dtmc.o
	$(CXX) $(CXXFLAGS) -o $@ $^

benchsingle: benchone_0 benchone_1 benchone_2 benchone_4 benchone_5 benchone_6 benchone_16 \
             benchone_17 benchone_18 benchone_20 benchone_21 benchone_22 benchone_24 benchone_25 \
             benchone_26 benchone_28 benchone_29 benchone_30

benchone_0 : bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongMC.cc StrongPA.cc
	$(CXX) $(CXXFLAGS)                                                                                          -o $@ $^

benchone_1 : bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongMC.cc StrongPA.cc
	$(CXX) $(CXXFLAGS) -DOPT_PARTITION                                                                          -o $@ $^

benchone_2 : bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongQ.cc
	$(CXX) $(CXXFLAGS)                                                                           -DOPT_QUOTIENT -o $@ $^

benchone_4 : bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongMC.cc StrongPA.cc
	$(CXX) $(CXXFLAGS)                 -DOPT_P_INVARIANT                                                        -o $@ $^

benchone_5 : bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongMC.cc StrongPA.cc
	$(CXX) $(CXXFLAGS) -DOPT_PARTITION -DOPT_P_INVARIANT                                                        -o $@ $^

benchone_6 : bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongQ.cc
	$(CXX) $(CXXFLAGS)                 -DOPT_P_INVARIANT                                         -DOPT_QUOTIENT -o $@ $^

benchone_16: bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongMC.cc StrongPA.cc
	$(CXX) $(CXXFLAGS)                                                          -DOPT_CACHE_NETS                -o $@ $^

benchone_17: bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongMC.cc StrongPA.cc
	$(CXX) $(CXXFLAGS) -DOPT_PARTITION                                          -DOPT_CACHE_NETS                -o $@ $^

benchone_18: bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongQ.cc
	$(CXX) $(CXXFLAGS)                                                          -DOPT_CACHE_NETS -DOPT_QUOTIENT -o $@ $^

benchone_20: bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongMC.cc StrongPA.cc
	$(CXX) $(CXXFLAGS)                 -DOPT_P_INVARIANT                        -DOPT_CACHE_NETS                -o $@ $^

benchone_21: bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongMC.cc StrongPA.cc
	$(CXX) $(CXXFLAGS) -DOPT_PARTITION -DOPT_P_INVARIANT                        -DOPT_CACHE_NETS                -o $@ $^

benchone_22: bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongQ.cc
	$(CXX) $(CXXFLAGS)                 -DOPT_P_INVARIANT                        -DOPT_CACHE_NETS -DOPT_QUOTIENT -o $@ $^

benchone_24: bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongMC.cc StrongPA.cc
	$(CXX) $(CXXFLAGS)                                   -DOPT_SIGNIFICIANT_ARC -DOPT_CACHE_NETS                -o $@ $^

benchone_25: bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongMC.cc StrongPA.cc
	$(CXX) $(CXXFLAGS) -DOPT_PARTITION                   -DOPT_SIGNIFICIANT_ARC -DOPT_CACHE_NETS                -o $@ $^

benchone_26: bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongQ.cc
	$(CXX) $(CXXFLAGS)                                   -DOPT_SIGNIFICIANT_ARC -DOPT_CACHE_NETS -DOPT_QUOTIENT -o $@ $^

benchone_28: bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongMC.cc StrongPA.cc
	$(CXX) $(CXXFLAGS)                 -DOPT_P_INVARIANT -DOPT_SIGNIFICIANT_ARC -DOPT_CACHE_NETS                -o $@ $^

benchone_29: bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongMC.cc StrongPA.cc
	$(CXX) $(CXXFLAGS) -DOPT_PARTITION -DOPT_P_INVARIANT -DOPT_SIGNIFICIANT_ARC -DOPT_CACHE_NETS                -o $@ $^

benchone_30: bench.o prmodel.o r-dtmc.o benchone.cc compactmaxflow.cc StrongQ.cc
	$(CXX) $(CXXFLAGS)                 -DOPT_P_INVARIANT -DOPT_SIGNIFICIANT_ARC -DOPT_CACHE_NETS -DOPT_QUOTIENT -o $@ $^
# DO NOT DELETE

benchone.cc: benchmark.h Simrel.h prmodel.h compactmaxflow.h \
  relationmap.h Strong.h StrongQ.h
compactmaxflow.cc: compactmaxflow.h relationmap.h
StrongMC.cc: Strong.h Simrel.h prmodel.h compactmaxflow.h \
  relationmap.h
StrongPA.cc: Strong.h Simrel.h prmodel.h compactmaxflow.h \
  relationmap.h
StrongQ.cc: StrongQ.h Strong.h Simrel.h prmodel.h \
  compactmaxflow.h relationmap.h
modelstat.o: prmodel.h
prmodel.o: prmodel.h
r-dtmc.o: modelbuilder.h prmodel.h
r-mdp.o: prmodel.h
bench.o: benchmark.h Simrel.h prmodel.h compactmaxflow.h \
  relationmap.h
benchmark.o: benchmark.h Simrel.h prmodel.h compactmaxflow.h \
  relationmap.h
