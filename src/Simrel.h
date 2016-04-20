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



#ifndef SIMREL_H
#define SIMREL_H

#include <set>
#include "prmodel.h"
#include "compactmaxflow.h"


// Interface class for simulation relation for various model types
class SimulationRelation
{
friend class CompactMaxFlow<double>;

public:
  SimulationRelation() : rmap()
  {
    label_func = 0;
#ifdef DEBUG
    memset(&stats, 0, sizeof(stats));
#endif//DEBUG
  }
  virtual ~SimulationRelation() {}
  
  virtual unsigned int Simulate(ProbabilisticModel*, std::set<std::pair<int,int> >*) = 0;
  
  static inline void SetFPPrecision(double p) { CompactMaxFlow<double>::precision = (p < 0.0 ? -p : p); }
  
  inline void SetLabelFunction(int(*lf)(void*,int),void *ud) { label_func = lf; label_user_data = ud; }
  inline int Label(int s) { return (*label_func)(label_user_data, s); }
  
#ifdef DEBUG
  SimulationStatistics stats;
#endif//DEBUG

protected:
  RelationMap rmap;

  int (*label_func)(void*,int);
  void *label_user_data;
};

#endif//SIMREL_H
