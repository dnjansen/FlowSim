/*****************************************************************************/
/*!
 *   Copyright 2009 Jonathan Bogdoll, Holger Hermanns, Lijun Zhang
 *
 *   This file is part of FLowSim.

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



#include "Strong.h"
#include "compactmaxflow.cc"

unsigned int StrongSimulation_CPA::Simulate(ProbabilisticModel *model, std::set<std::pair<int,int> > *result)
{
  unsigned int r = StrongSimulation_PA::Simulate(model, result);
#ifdef DEBUG
  stats.mem_model += m->na * sizeof(double);
#endif//DEBUG
  delete [] dist_sums;
  return r;
}

// Compute the initial relation
int StrongSimulation_CPA::BuildRelationMap()
{
  int k, n, size = 0;
  
  // Compute sum of each distribution
  dist_sums = new double[m->na];
  for (int j, i = 0; i < m->na; ++i)
  {
    dist_sums[i] = 0.0;
    for (j = row_starts[i]; j < row_starts[i+1]; ++j) dist_sums[i] += non_zeros[j];
    for (j = row_starts[i]; j < row_starts[i+1]; ++j) non_zeros[j] /= dist_sums[i];
  }

  for (k = 0; k < n_states - 1; ++k)
  {
    for (n = k + 1; n < n_states; ++n)
    {
      if (Label(k) == Label(n))
      {
        size += 2;
        rmap.Set(k, n);
        rmap.Set(n, k);
      }
    }
  }
  
  rmap.Commit();
  
  return size;
}

// Compute strong simulation on a particular pair in the relation
bool StrongSimulation_CPA::DecideStrongSimulation(Pair *p)
{
  int a1, a2;
  bool sim_one;
  CompactMaxFlow<double> *sim;
  
  // Partitioning sorts the actions and we can scan in sub-quadratic time
  //if (optflags & OPT_PARTITION)
#ifdef OPT_PARTITION
  {
    for (a1 = state_starts[p->x], a2 = state_starts[p->y]; a1 < state_starts[p->x + 1]; ++a1)
    {
      // Seek to the first matching action
           if (actions[a2] <  actions[a1]) while (actions[a2  ] <  actions[a1] && a2 < state_starts[p->y + 1] - 1) ++a2;
      else if (actions[a2] >= actions[a1]) while (actions[a2-1] >= actions[a1] && a2 > state_starts[p->y    ]    ) --a2;
      
      if (actions[a1] != actions[a2]) return false;
      if (!CompactMaxFlow<double>::_Tleq(dist_sums[a1], dist_sums[a2])) continue;
      
      sim_one = false;
      for (; a2 < state_starts[p->y + 1] && actions[a1] == actions[a2]; ++a2)
      {
        sim = ConstructNetwork(a1, a2);
        if (sim)
        {
          sim_one = sim->IsFlowTotal();
          delete sim;
          if (sim_one) break;
        }
      }
      
      if (!sim_one) return false;
    }
  }
  //else
#else
  {
    for (a1 = state_starts[p->x]; a1 < state_starts[p->x + 1]; ++a1)
    {
      sim_one = false;
      for (a2 = state_starts[p->y]; a2 < state_starts[p->y + 1]; ++a2)
      {
        if (!CompactMaxFlow<double>::_Tleq(dist_sums[a1], dist_sums[a2])) continue;
        if (actions[a1] == actions[a2])
        {
          sim = ConstructNetwork(a1, a2);
          if (sim)
          {
            sim_one = sim->IsFlowTotal();
            delete sim;
            if (sim_one) break;
          }
        }
      }
      
      if (!sim_one) return false;
    }
  }
#endif
  
  return true;
}
