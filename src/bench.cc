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



#include "benchmark.h"
#include <sys/resource.h>
#include <sys/time.h>
#include <math.h>

Benchmark::Benchmark()
{
  tu = 0;
  ts = 0;
  tr = 0;
  ss = 0;
  rows = 0;
  cb_progress = 0;
}

Benchmark::~Benchmark()
{
  if (tu) delete [] tu;
  if (ts) delete [] ts;
  if (tr) delete [] tr;
  if (ss) delete [] ss;
}

// Perform the benchmark on the model pm using the given simulator. flagvector
// is an array of bitsets specifying which optimizations to use; see Strong.h.
// flagvectorsize is the number of items in flagvector. averages is the number
// of times the simulation is to be performed and averaged over.
void Benchmark::Bench(ProbabilisticModel *pm, SimulationRelation *simulator, unsigned long*, unsigned int flagvectorsize, unsigned int averages, FILE *result)
{
  rusage r1, r2, rc1, rc2;
  timeval t1, t2;
  double t_usr = 0.0, t_sys = 0.0, t_real = 0.0;
  std::set<std::pair<int,int> > resultset;
  std::set<std::pair<int,int> >::iterator i;
  
  if (tu) delete [] tu;
  if (ts) delete [] ts;
  if (tr) delete [] tr;
  if (ss) delete [] ss;
  
  tu = new double[flagvectorsize];
  ts = new double[flagvectorsize];
  tr = new double[flagvectorsize];
  ss = new SimulationStatistics[flagvectorsize];
  rows = flagvectorsize;
  
  states = pm->States();
  transitions = pm->Transitions();
  
  for (unsigned int m = 0; m < flagvectorsize; ++m)
  {
    t_usr = 0.0;
    t_sys = 0.0;
    t_real = 0.0;
    for (unsigned int n = 0; n < averages; ++n)
    {
      getrusage(RUSAGE_SELF, &r1);
      getrusage(RUSAGE_CHILDREN, &rc1);
      gettimeofday(&t1, 0);
      simulator->Simulate(pm, (result && n == 0 && m == 0 ? &resultset : 0));
      getrusage(RUSAGE_SELF, &r2);
      getrusage(RUSAGE_CHILDREN, &rc2);
      gettimeofday(&t2, 0);
      
      if (n == 0)
      {
        memcpy(ss + m, &simulator->stats, sizeof(SimulationStatistics));
        if (result && m == 0)
        {
          for (i = resultset.begin(); i != resultset.end(); ++i) fprintf(result, "%d %d\n", i->first, i->second);
        }
      }
      else
      {
        if (simulator->stats.min_complexity < (ss + m)->min_complexity) (ss + m)->min_complexity = simulator->stats.min_complexity;
        if (simulator->stats.max_complexity > (ss + m)->max_complexity) (ss + m)->max_complexity = simulator->stats.max_complexity;
      }
      
      t_usr += double((r2.ru_utime.tv_sec) - (r1.ru_utime.tv_sec))
            +  double((rc2.ru_utime.tv_sec) - (rc1.ru_utime.tv_sec))
            + ((double(r2.ru_utime.tv_usec) * 0.000001) - (double(r1.ru_utime.tv_usec) * 0.000001))
            + ((double(rc2.ru_utime.tv_usec) * 0.000001) - (double(rc1.ru_utime.tv_usec) * 0.000001));
      t_sys += double((r2.ru_stime.tv_sec) - (r1.ru_stime.tv_sec))
            +  double((rc2.ru_stime.tv_sec) - (rc1.ru_stime.tv_sec))
            + ((double(r2.ru_stime.tv_usec) * 0.000001) - (double(r1.ru_stime.tv_usec) * 0.000001))
            + ((double(rc2.ru_stime.tv_usec) * 0.000001) - (double(rc1.ru_stime.tv_usec) * 0.000001));
      t_real += double((t2.tv_sec) - (t1.tv_sec)) + ((double(t2.tv_usec) * 0.000001) - (double(t1.tv_usec) * 0.000001));
      
      if (cb_progress) (*cb_progress)(n + (m * averages), averages * flagvectorsize);
    }
    
    tu[m] = t_usr / double(averages);
    ts[m] = t_sys / double(averages);
    tr[m] = t_real / double(averages);
  }
}
