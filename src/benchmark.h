/*****************************************************************************/
/*!
 *   Copyright 2009 Jonathan Bogdoll, Holger Hermanns, Lijun Zhang
 *
 *   This file is part of FlowSim.

 *   FlowSim is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   FlowSim is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with FlowSim.  If not, see <http://www.gnu.org/licenses/>.
 */
/*****************************************************************************/



#ifndef _BENCHMARK_H_
#define _BENCHMARK_H_

#ifndef DEBUG
#error "#define DEBUG" required to build benchmark
#endif

#include "Simrel.h"

struct RandomModel
{
  int n, a, b, c;
  unsigned int avg, steps, xtarget, ztarget, labels;
  double cb, lb, fb, pb, sb, xstart, xend, zstart, zend;
};

struct RandomModelPlotInfo
{
  RandomModel mdl;
  char datasource[160];
  unsigned int id;
};

// Run a benchmark over a specific model using different optimizations
class Benchmark
{
public:
  Benchmark();
  ~Benchmark();
  
  void Bench(ProbabilisticModel*, SimulationRelation*, unsigned long*, unsigned int, unsigned int, FILE* = 0);
  
  inline double GetUserTime(unsigned int i) { if (!tu) return 0.0; if (i >= rows) return 0.0; return tu[i]; }
  inline double GetSystemTime(unsigned int i) { if (!ts) return 0.0; if (i >= rows) return 0.0; return ts[i]; }
  inline double GetRealTime(unsigned int i) { if (!tr) return 0.0; if (i >= rows) return 0.0; return tr[i]; }
  inline const SimulationStatistics* GetStats(unsigned int i) { if (!ss || i >= rows) return 0; return ss + i; }
  inline int GetStates() { return states; }
  inline double GetTransitions() { return transitions; }
  void SetProgressCallback(void(*cb)(unsigned int,unsigned int)) { cb_progress = cb; }
  
private:
  double *tu, *ts, *tr, transitions;
  SimulationStatistics *ss;
  unsigned int rows;
  int states;
  void (*cb_progress)(unsigned int, unsigned int);
};

#endif//_BENCHMARK_H_
