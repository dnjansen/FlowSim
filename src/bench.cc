/*****************************************************************************/
/*!
 *   Copyright 2009-2015 Jonathan Bogdoll, Holger Hermanns, Lijun Zhang,
 *                       David N. Jansen
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
  tu_stdev = NULL;
  ts_stdev = NULL;
  tr_stdev = NULL;
  ss = 0;
  rows = 0;
  cb_progress = 0;
}

Benchmark::~Benchmark()
{
  if (tu) delete [] tu;
  if (ts) delete [] ts;
  if (tr) delete [] tr;
  if (NULL != tu_stdev) delete [] tu_stdev;
  if (NULL != ts_stdev) delete [] ts_stdev;
  if (NULL != tr_stdev) delete [] tr_stdev;
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
  std::set<std::pair<int,int> > resultset;
  std::set<std::pair<int,int> >::iterator i;
  
  if (tu) delete [] tu;
  if (ts) delete [] ts;
  if (tr) delete [] tr;
  if (NULL != tu_stdev) delete [] tu_stdev;
  if (NULL != ts_stdev) delete [] ts_stdev;
  if (NULL != tr_stdev) delete [] tr_stdev;
  if (ss) delete [] ss;
  
  tu = new double[flagvectorsize];
  ts = new double[flagvectorsize];
  tr = new double[flagvectorsize];
  if (averages > 1) {
    tu_stdev = new double[flagvectorsize];
    ts_stdev = new double[flagvectorsize];
    tr_stdev = new double[flagvectorsize];
  } else {
    tu_stdev = NULL;
    ts_stdev = NULL;
    tr_stdev = NULL;
  }
  ss = new SimulationStatistics[flagvectorsize];
  rows = flagvectorsize;
  
  states = pm->States();
  transitions = pm->Transitions();
  
  for (unsigned int m = 0; m < flagvectorsize; ++m)
  {
    double t_usr = 0.0, t_sys = 0.0, t_real = 0.0;
    double t_usr_estimate = 0.0, t_sys_estimate = 0.0, t_real_estimate = 0.0;
    double t_usr_square = 0.0, t_sys_square = 0.0, t_real_square = 0.0;

    for (unsigned int n = 0; n < averages; ++n)
    {
      fputc('.', stderr);
      getrusage(RUSAGE_SELF, &r1);
      getrusage(RUSAGE_CHILDREN, &rc1);
      gettimeofday(&t1, 0);
      simulator->Simulate(pm, (result && n == 0 && m == 0 ? &resultset : 0));
      getrusage(RUSAGE_SELF, &r2);
      getrusage(RUSAGE_CHILDREN, &rc2);
      gettimeofday(&t2, 0);
      
      /* numerically stable calculation of the standard deviation: The time
      of the first experiment is used as an estimate of the mean. Then, the
      variance is calculated for (time - estimate). This generally leads to
      smaller numbers, so there will be less cancellation than in a naive
      implementation of "expectation of square minus square of expectation".

      This algorithm is copied from
      https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Computing_shifted_data
      which refers to T.F. Chan, G.H. Golub and R.J. LeVeque: Algorithms for
      computing the sample variance: Analysis and recommendations. The American
      Statistician, 37(3)1983. pp. 242--247. (See in particular Section 3 of
      that article.) */

      double tempt_usr = double(r2.ru_utime.tv_sec) - double(r1.ru_utime.tv_sec)
            +  double(rc2.ru_utime.tv_sec) - double(rc1.ru_utime.tv_sec)
            + ((double(r2.ru_utime.tv_usec) * 0.000001) - (double(r1.ru_utime.tv_usec) * 0.000001))
            + ((double(rc2.ru_utime.tv_usec) * 0.000001) - (double(rc1.ru_utime.tv_usec) * 0.000001));
      double tempt_sys = double(r2.ru_stime.tv_sec) - double(r1.ru_stime.tv_sec)
            +  double(rc2.ru_stime.tv_sec) - double(rc1.ru_stime.tv_sec)
            + ((double(r2.ru_stime.tv_usec) * 0.000001) - (double(r1.ru_stime.tv_usec) * 0.000001))
            + ((double(rc2.ru_stime.tv_usec) * 0.000001) - (double(rc1.ru_stime.tv_usec) * 0.000001));
      double tempt_real = double(t2.tv_sec) - double(t1.tv_sec) + ((double(t2.tv_usec) * 0.000001) - (double(t1.tv_usec) * 0.000001));

      if (n == 0)
      {
        memcpy(ss + m, &simulator->stats, sizeof(SimulationStatistics));
        if (result && m == 0)
        {
          for (i = resultset.begin(); i != resultset.end(); ++i) fprintf(result, "%d %d\n", i->first, i->second);
        }
        t_usr_estimate  = tempt_usr;
        t_sys_estimate  = tempt_sys;
        t_real_estimate = tempt_real;
      }
      else
      {
        if (simulator->stats.min_complexity < ss[m].min_complexity)
          ss[m].min_complexity = simulator->stats.min_complexity;
        if (simulator->stats.max_complexity > ss[m].max_complexity)
          ss[m].max_complexity = simulator->stats.max_complexity;

        tempt_usr -= t_usr_estimate;
        t_usr += tempt_usr;
        t_usr_square += tempt_usr * tempt_usr;
        tempt_sys -= t_sys_estimate;
        t_sys += tempt_sys;
        t_sys_square += tempt_sys * tempt_sys;
        tempt_real -= t_real_estimate;
        t_real += tempt_real;
        t_real_square += tempt_real * tempt_real;
      }
      
      if (cb_progress) (*cb_progress)(n + (m * averages), averages * flagvectorsize);
    }
    
    tu[m] = t_usr  / averages + t_usr_estimate ;
    ts[m] = t_sys  / averages + t_sys_estimate ;
    tr[m] = t_real / averages + t_real_estimate;
    if (averages > 1) {
      tu_stdev[m] = sqrt((t_usr_square  - t_usr  * t_usr  / averages) / (averages - 1));
      ts_stdev[m] = sqrt((t_sys_square  - t_sys  * t_sys  / averages) / (averages - 1));
      tr_stdev[m] = sqrt((t_real_square - t_real * t_real / averages) / (averages - 1));
    }
  }
}
