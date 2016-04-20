/*****************************************************************************/
/*!
 *   Copyright 2009-2014 Jonathan Bogdoll, Holger Hermanns, Lijun Zhang,
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



#ifndef _STRONG_QUOTIENT_H_
#define _STRONG_QUOTIENT_H_

#include <vector>
#include <map>
#include "Strong.h"

class StrongSimulation_Quotient : public StrongSimulation
{
protected:
  // Abstraction from model type
  class Simulator
  {
  public:
    Simulator(StrongSimulation_Quotient *ssq) { assert(ssq); base = ssq; }
    virtual ~Simulator() {}
    virtual bool sRs(int,int,bool) = 0;
    virtual bool qRq(int,int,bool) = 0;
    bool muRmu(int,int,bool);
    void Flush();
  protected:
    StrongSimulation_Quotient *base;
  private:
    Simulator() {}
#ifdef OPT_CACHE_NETS
    std::map<std::pair<int,int>,CompactMaxFlow<double>*> cache;
#endif//OPT_CACHE_NETS
  };
  class Simulator_MC : public Simulator
  {
  public:
    Simulator_MC(StrongSimulation_Quotient *ssq) : Simulator(ssq) {}
    virtual ~Simulator_MC() {}
    bool sRs(int,int,bool);
    bool qRq(int,int,bool);
  };
  class Simulator_PA : public Simulator
  {
  public:
    Simulator_PA(StrongSimulation_Quotient *ssq) : Simulator(ssq) {}
    virtual ~Simulator_PA() {}
    bool sRs(int,int,bool);
    bool qRq(int,int,bool);
  };
  
friend class Simulator;
  
public:
  StrongSimulation_Quotient() {}
  virtual ~StrongSimulation_Quotient() {}
  
  // Main interface
  unsigned int Simulate(ProbabilisticModel*, std::set<std::pair<int,int> >*);

#ifdef WITH_VERIFIER
  bool Verify(ProbabilisticModel*, std::set<std::pair<int,int> >&, std::set<std::pair<int,int> >*, std::set<std::pair<int,int> >*) { return false; };
#endif//WITH_VERIFIER
  
protected:
  // Member variables
  ProbabilisticModel *model;
  int nStates, nDistributions, nBlocks;
  Simulator *sim;
  unsigned char *action_masks;
  int action_mask_pitch;
  int *partition, *new_partition;
  int iteration;
  std::vector<std::set<int> > sigma;
  std::multimap<int,std::pair<int,int> > forall;
  
#ifdef QLOG
  FILE *qlog;
#endif
  
  // Model structures
  int *state_starts;
  int *actions;
  int *row_starts;
  int *cols;
  double *non_zeros;
  int *lifted_row_starts;
  int *lifted_cols;
  double *lifted_non_zeros;
  
  // Initialization
  void InitializeActionMasks();
  void InitializeRelation();
  
  // Algorithm subroutines
  void LiftDistributions();
  void FindCommonDistributions_PA();
  void FindCommonDistributions_MC();
  void PurgePartitionRelation();
  
  // Auxiliaries
  inline bool Act_Equal(int s1, int s2)
  {
    if (!action_masks) return true;
    return (memcmp(action_masks + (s1 * action_mask_pitch), action_masks + (s2 * action_mask_pitch), action_mask_pitch) == 0);
  }
  inline bool Act_Subseteq(int s1, int s2)
  {
    s1 *= action_mask_pitch;
    s2 *= action_mask_pitch;
    if (!action_masks) return true;
    for (int i = 0; i < action_mask_pitch; ++i)
    {
      if (action_masks[s1 + i] & ~action_masks[s2 + i]) return false;
    }
    return true;
  }
  inline bool Dist_Equal(int mu1, int mu2)
  {
    int a, b;
    if (lifted_row_starts[mu1 + 1] - lifted_row_starts[mu1] != lifted_row_starts[mu2 + 1] - lifted_row_starts[mu2]) return false;
    for (a = lifted_row_starts[mu1]; a < lifted_row_starts[mu1 + 1]; ++a)
    {
      for (b = lifted_row_starts[mu2]; b < lifted_row_starts[mu2 + 1]; ++b)
      {
        if (lifted_cols[a] == lifted_cols[b] && CompactMaxFlow<double>::_Teq(lifted_non_zeros[a], lifted_non_zeros[b])) break;
      }
      
      if (b == lifted_row_starts[mu2 + 1]) return false;
    }
    return true;
  }
  inline bool Steps_Subseteq(int s1, int s2)
  {
    if (!action_masks) return Dist_Equal(s1, s2);
    int a, b;
    for (a = state_starts[s1]; a < state_starts[s1 + 1]; ++a)
    {
      for (b = state_starts[s2]; b < state_starts[s2 + 1]; ++b)
      {
        if (actions[a] == actions[b] && Dist_Equal(a, b)) break;
      }
      
      if (b == state_starts[s2 + 1]) return false;
    }
    return true;
  }
  inline bool Steps_Equal(int s1, int s2)
  {
    if (!action_masks) return Dist_Equal(s1, s2);
    if (state_starts[s1 + 1] - state_starts[s1] != state_starts[s2 + 1] - state_starts[s2]) return false;
    return Steps_Subseteq(s1, s2) && Steps_Subseteq(s2, s1);
  }
  
  // Compute one iteration of the algorithm
  int Iterate();
};

#endif//_STRONG_QUOTIENT_H_
