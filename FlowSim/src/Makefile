CXX=g++
#CXXFLAGS=-O0 -ggdb -Wpointer-arith -Winline -Wall -W -pedantic -DDEBUG
CXXFLAGS=-O0 -s -Wpointer-arith -Winline -Wall -W -DDEBUG

.cc.o:
	$(CXX) $(CXXFLAGS) -c $<

all: r-dtmc r-mdp flowsim benchsingle modelstat

depend:
	makedepend -Y *.cc *.h

clean:
	rm -f *~ *.o *.core

clean-all: clean
	rm -f r-dtmc r-mdp flowsim benchone_* modelstat

r-dtmc: r-dtmc.cc prmodel.o
	$(CXX) $(CXXFLAGS) -o r-dtmc -D_STANDALONE_ r-dtmc.cc prmodel.o

r-mdp: r-mdp.cc prmodel.o
	$(CXX) $(CXXFLAGS) -o r-mdp -D_STANDALONE_ r-mdp.cc prmodel.o

modelstat: modelstat.cc prmodel.o
	$(CXX) $(CXXFLAGS) -o modelstat modelstat.cc prmodel.o

flowsim: benchmark.o prmodel.o compactmaxflow.cc compactmaxflow.h r-dtmc.o
	$(CXX) $(CXXFLAGS) -o flowsim benchmark.o prmodel.o r-dtmc.o

benchsingle: benchone_0 benchone_1 benchone_2 benchone_4 benchone_5 benchone_6 benchone_16 \
             benchone_17 benchone_18 benchone_20 benchone_21 benchone_22 benchone_24 benchone_25 \
             benchone_26 benchone_28 benchone_29 benchone_30

benchone_0: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim
	$(CXX) $(CXXFLAGS) -o benchone_0 benchone.cc prmodel.o r-dtmc.o StrongMC.cc StrongPA.cc

benchone_1: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim
	$(CXX) $(CXXFLAGS) -DOPT_PARTITION -o benchone_1 benchone.cc prmodel.o r-dtmc.o StrongMC.cc StrongPA.cc

benchone_2: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim StrongQ.cc
	$(CXX) $(CXXFLAGS) -DOPT_QUOTIENT -o benchone_2 benchone.cc prmodel.o r-dtmc.o StrongQ.cc

benchone_4: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim
	$(CXX) $(CXXFLAGS) -DOPT_P_INVARIANT -o benchone_4 benchone.cc prmodel.o r-dtmc.o StrongMC.cc StrongPA.cc

benchone_5: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim
	$(CXX) $(CXXFLAGS) -DOPT_P_INVARIANT -DOPT_PARTITION -o benchone_5 benchone.cc prmodel.o r-dtmc.o StrongMC.cc StrongPA.cc

benchone_6: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim StrongQ.cc
	$(CXX) $(CXXFLAGS) -DOPT_P_INVARIANT -DOPT_QUOTIENT -o benchone_6 benchone.cc prmodel.o r-dtmc.o StrongQ.cc

benchone_16: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -o benchone_16 benchone.cc prmodel.o r-dtmc.o StrongMC.cc StrongPA.cc

benchone_17: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_PARTITION -o benchone_17 benchone.cc prmodel.o r-dtmc.o StrongMC.cc StrongPA.cc

benchone_18: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim StrongQ.cc
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_QUOTIENT -o benchone_18 benchone.cc prmodel.o r-dtmc.o StrongQ.cc

benchone_20: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_P_INVARIANT -o benchone_20 benchone.cc prmodel.o r-dtmc.o StrongMC.cc StrongPA.cc

benchone_21: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_PARTITION -DOPT_P_INVARIANT -o benchone_21 benchone.cc prmodel.o r-dtmc.o StrongMC.cc StrongPA.cc

benchone_22: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim StrongQ.cc
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_QUOTIENT -DOPT_P_INVARIANT -o benchone_22 benchone.cc prmodel.o r-dtmc.o StrongQ.cc

benchone_24: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_SIGNIFICIANT_ARC -o benchone_24 benchone.cc prmodel.o r-dtmc.o StrongMC.cc StrongPA.cc

benchone_25: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_PARTITION -DOPT_SIGNIFICIANT_ARC -o benchone_25 benchone.cc prmodel.o r-dtmc.o StrongMC.cc StrongPA.cc

benchone_26: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim StrongQ.cc
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_QUOTIENT -DOPT_SIGNIFICIANT_ARC -o benchone_26 benchone.cc prmodel.o r-dtmc.o StrongQ.cc

benchone_28: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_P_INVARIANT -DOPT_SIGNIFICIANT_ARC -o benchone_28 benchone.cc prmodel.o r-dtmc.o StrongMC.cc StrongPA.cc

benchone_29: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_PARTITION -DOPT_P_INVARIANT -DOPT_SIGNIFICIANT_ARC -o benchone_29 benchone.cc prmodel.o r-dtmc.o StrongMC.cc StrongPA.cc

benchone_30: prmodel.o r-dtmc.o bench.cc benchone.cc flowsim StrongQ.cc
	$(CXX) $(CXXFLAGS) -DOPT_CACHE_NETS -DOPT_QUOTIENT -DOPT_SIGNIFICIANT_ARC -DOPT_P_INVARIANT -o benchone_30 benchone.cc prmodel.o r-dtmc.o StrongQ.cc
# DO NOT DELETE

StrongMC.o: Strong.h Simrel.h prmodel.h compactmaxflow.h relationmap.h
StrongMC.o: compactmaxflow.cc
StrongPA.o: Strong.h Simrel.h prmodel.h compactmaxflow.h relationmap.h
StrongPA.o: compactmaxflow.cc
StrongQ.o: StrongQ.h Strong.h Simrel.h prmodel.h compactmaxflow.h
StrongQ.o: relationmap.h compactmaxflow.cc
bench.o: benchmark.h prmodel.h Strong.h Simrel.h compactmaxflow.h
bench.o: relationmap.h
benchmark.o: benchmark.h prmodel.h Strong.h Simrel.h compactmaxflow.h
benchmark.o: relationmap.h bench.cc compactmaxflow.cc
benchone.o: benchmark.h prmodel.h Strong.h Simrel.h compactmaxflow.h
benchone.o: relationmap.h StrongQ.h bench.cc
compactmaxflow.o: compactmaxflow.h relationmap.h
modelstat.o: prmodel.h
prmodel.o: prmodel.h
r-dtmc.o: prmodel.h modelbuilder.h
r-mdp.o: prmodel.h
Simrel.o: prmodel.h compactmaxflow.h relationmap.h
Strong.o: Simrel.h prmodel.h compactmaxflow.h relationmap.h
StrongQ.o: Strong.h Simrel.h prmodel.h compactmaxflow.h relationmap.h
benchmark.o: prmodel.h Strong.h Simrel.h compactmaxflow.h relationmap.h
compactmaxflow.o: relationmap.h
modelbuilder.o: prmodel.h
