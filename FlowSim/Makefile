CXX=g++
CXXFLAGS=-O0 -ggdb -Wpointer-arith -Winline -Wall -W -pedantic -DDEBUG -DQLOG
#CXXFLAGS=-O0 -s -Wpointer-arith -Winline -Wall -W -pedantic

.cc.o:
	$(CXX) $(CXXFLAGS) -c $<

all: r-dtmc r-mdp benchmark benchsingle modelstat verify breakpoints

depend:
	makedepend -Y *.cc *.h

clean:
	rm -f *~ *.o *.core

r-dtmc: r-dtmc.cc prmodel.o
	$(CXX) $(CXXFLAGS) -o r-dtmc -D_STANDALONE_ r-dtmc.cc prmodel.o

r-mdp: r-mdp.cc prmodel.o
	$(CXX) $(CXXFLAGS) -o r-mdp -D_STANDALONE_ r-mdp.cc prmodel.o

dtmctest: dtmctest.cc prmodel.o
	$(CXX) $(CXXFLAGS) -o dtmctest dtmctest.cc prmodel.o
	$(CXX) $(CXXFLAGS) -o dtmctest_pmf -DOPT_CACHE_NETS dtmctest.cc prmodel.o

patest: patest.cc prmodel.o StrongQ.cc
	$(CXX) $(CXXFLAGS) -o patest patest.cc prmodel.o StrongQ.cc

verify: verify.cc
	$(CXX) $(CXXFLAGS) -o verify verify.cc prmodel.o

breakpoints: breakpoints.cc Weak.h
	$(CXX) $(CXXFLAGS) -o breakpoints breakpoints.cc

modelstat: modelstat.cc prmodel.o
	$(CXX) $(CXXFLAGS) -o modelstat modelstat.cc prmodel.o

flow: main.o prmodel.o StrongMF.o StrongMC.o StrongDTMC.o StrongCTMC.o StrongPA.o StrongCPA.o r-dtmc.o r-mdp.o
	$(CXX) $(CXXFLAGS) -o flow main.o prmodel.o StrongMF.o StrongMC.o StrongDTMC.o StrongCTMC.o StrongPA.o StrongCPA.o r-dtmc.o r-mdp.o

benchmark: benchmark.o prmodel.o compactmaxflow.cc compactmaxflow.h r-dtmc.o
	$(CXX) $(CXXFLAGS) -o benchmark benchmark.o prmodel.o r-dtmc.o

benchsingle: benchone_0 benchone_1 benchone_2 benchone_4 benchone_5 benchone_6 benchone_16 \
             benchone_17 benchone_18 benchone_20 benchone_21 benchone_22 benchone_24 benchone_25 \
             benchone_26 benchone_28 benchone_29 benchone_30

benchone_0: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark
	$(CXX) $(CXXFLAGS) -o benchone_0 benchone.cc prmodel.o r-dtmc.o

benchone_1: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark
	$(CXX) $(CXXFLAGS) -DOPT_PARTITION -o benchone_1 benchone.cc prmodel.o r-dtmc.o

benchone_2: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark StrongQ.cc
	$(CXX) $(CXXFLAGS) -DOPT_QUOTIENT -o benchone_2 benchone.cc prmodel.o r-dtmc.o StrongQ.cc

benchone_4: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark
	$(CXX) $(CXXFLAGS) -DOPT_P_INVARIANT -o benchone_4 benchone.cc prmodel.o r-dtmc.o

benchone_5: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark
	$(CXX) $(CXXFLAGS) -DOPT_P_INVARIANT -DOPT_PARTITION -o benchone_5 benchone.cc prmodel.o r-dtmc.o

benchone_6: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark StrongQ.cc
	$(CXX) $(CXXFLAGS) -DOPT_P_INVARIANT -DOPT_QUOTIENT -o benchone_6 benchone.cc prmodel.o r-dtmc.o StrongQ.cc

benchone_16: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -o benchone_16 benchone.cc prmodel.o r-dtmc.o

benchone_17: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_PARTITION -o benchone_17 benchone.cc prmodel.o r-dtmc.o

benchone_18: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark StrongQ.cc
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_QUOTIENT -o benchone_18 benchone.cc prmodel.o r-dtmc.o StrongQ.cc

benchone_20: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_P_INVARIANT -o benchone_20 benchone.cc prmodel.o r-dtmc.o

benchone_21: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_PARTITION -DOPT_P_INVARIANT -o benchone_21 benchone.cc prmodel.o r-dtmc.o

benchone_22: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark StrongQ.cc
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_QUOTIENT -DOPT_P_INVARIANT -o benchone_22 benchone.cc prmodel.o r-dtmc.o StrongQ.cc

benchone_24: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_SIGNIFICIANT_ARC -o benchone_24 benchone.cc prmodel.o r-dtmc.o

benchone_25: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_PARTITION -DOPT_SIGNIFICIANT_ARC -o benchone_25 benchone.cc prmodel.o r-dtmc.o

benchone_26: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark StrongQ.cc
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_QUOTIENT -DOPT_SIGNIFICIANT_ARC -o benchone_26 benchone.cc prmodel.o r-dtmc.o StrongQ.cc

benchone_28: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_P_INVARIANT -DOPT_SIGNIFICIANT_ARC -o benchone_28 benchone.cc prmodel.o r-dtmc.o

benchone_29: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_PARTITION -DOPT_P_INVARIANT -DOPT_SIGNIFICIANT_ARC -o benchone_29 benchone.cc prmodel.o r-dtmc.o

benchone_30: prmodel.o r-dtmc.o bench.cc benchone.cc benchmark StrongQ.cc
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_QUOTIENT -DOPT_SIGNIFICIANT_ARC -DOPT_P_INVARIANT -o benchone_30 benchone.cc prmodel.o r-dtmc.o StrongQ.cc
# DO NOT DELETE

StrongCPA.o: Strong.h prmodel.h compactmaxflow.h relationmap.h
StrongCPA.o: compactmaxflow.cc
StrongCTMC.o: Strong.h prmodel.h compactmaxflow.h relationmap.h
StrongDTMC.o: Strong.h prmodel.h compactmaxflow.h relationmap.h
StrongMC.o: Strong.h prmodel.h compactmaxflow.h relationmap.h
StrongMC.o: compactmaxflow.cc
StrongMF.o: StrongMF.h prmodel.h maxflow.h maxflow.cc
StrongPA.o: Strong.h prmodel.h compactmaxflow.h relationmap.h
StrongPA.o: compactmaxflow.cc
bench.o: benchmark.h prmodel.h Strong.h compactmaxflow.h relationmap.h
benchmark.o: benchmark.h prmodel.h Strong.h compactmaxflow.h relationmap.h
benchmark.o: bench.cc
benchone.cc: benchmark.h prmodel.h Strong.h compactmaxflow.h relationmap.h
benchone.cc: StrongMC.cc compactmaxflow.cc StrongDTMC.cc StrongCTMC.cc
benchone.cc: StrongPA.cc StrongCPA.cc bench.cc
benchone.o: benchmark.h prmodel.h Strong.h compactmaxflow.h relationmap.h
benchone.o: StrongMC.cc compactmaxflow.cc StrongDTMC.cc StrongCTMC.cc
benchone.o: StrongPA.cc StrongCPA.cc bench.cc
verify.cc: prmodel.h Strong.h compactmaxflow.h relationmap.h
verify.cc: StrongMC.cc compactmaxflow.cc StrongDTMC.cc StrongCTMC.cc
verify.cc: StrongPA.cc StrongCPA.cc
compactmaxflow.o: compactmaxflow.h relationmap.h
main.o: prmodel.h StrongMF.h maxflow.h Strong.h compactmaxflow.h
main.o: relationmap.h
maxflow.o: maxflow.h
modelstat.o: prmodel.h
prmodel.o: prmodel.h
r-dtmc.o: prmodel.h
r-mdp.o: prmodel.h
test.o: maxflow.cc maxflow.h
Strong.o: prmodel.h compactmaxflow.h relationmap.h
StrongMF.o: prmodel.h maxflow.h
benchmark.o: prmodel.h Strong.h compactmaxflow.h relationmap.h
compactmaxflow.o: relationmap.h
StrongQ.o: StrongQ.cc StrongQ.h Strong.h compactmaxflow.cc compactmaxflow.h
