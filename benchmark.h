#ifndef _BENCHMARK_H_
#define _BENCHMARK_H_

#ifndef DEBUG
#error "#define DEBUG" required to build benchmark
#endif

#include "prmodel.h"
#include "Strong.h"
#include <math.h>

#ifndef NAN
#define NAN -1
#endif

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
