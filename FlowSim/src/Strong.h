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



#ifndef STRONG_SIMULATION_H
#define STRONG_SIMULATION_H

#include <set>
#include <iterator>
#include <algorithm>
#include <map>
#include <vector>
#include <utility>
#include <assert.h>
#include "Simrel.h"

// Interface class for strong simulation for various model types
class StrongSimulation : public SimulationRelation
{
friend class CompactMaxFlow<double>;

public:
  StrongSimulation() {}
  virtual ~StrongSimulation() {}
  
  virtual unsigned int Simulate(ProbabilisticModel*, std::set<std::pair<int,int> >*) = 0;

#ifdef WITH_VERIFIER
  virtual bool Verify(ProbabilisticModel*, std::set<std::pair<int,int> >&, std::set<std::pair<int,int> >*, std::set<std::pair<int,int> >*) = 0;
#endif//WITH_VERIFIER

protected:
  // Dummy order on states; overriden by classes that sort states
  virtual bool state_Less(int s1, int s2) { return s1 < s2; }
  
  // Functor defining an order on the states
  class StateOrder
  {
  public:
    StateOrder(StrongSimulation *_s) { s = _s; }
    ~StateOrder() {}
    
    inline bool operator()(int a, int b) { return s->state_Less(a, b); }
    
  private:
    StateOrder();
    StrongSimulation *s;
  };
};

// Private intermediate class containing functions and structures used
// to simulate both DTMCs and CTMCs
class StrongSimulation_MC : public StrongSimulation
{
public:
  StrongSimulation_MC() {}
  ~StrongSimulation_MC() {}
  
  unsigned int Simulate(ProbabilisticModel*, std::set<std::pair<int,int> >*);

#ifdef WITH_VERIFIER
  bool Verify(ProbabilisticModel*, std::set<std::pair<int,int> >&, std::set<std::pair<int,int> >*, std::set<std::pair<int,int> >*);
#endif//WITH_VERIFIER
  
protected:
  struct Pair
  {
  public:
    int x, y;
#ifdef OPT_CACHE_NETS
    CompactMaxFlow<double> *simulation;
#endif
    Pair *next;
  };
  
  // Functor defining an order on the successors of a state (used with std::sort)
  class SuccessorOrder
  {
  public:
    SuccessorOrder(StrongSimulation_MC *_s) { s = _s; }
    ~SuccessorOrder() {}
    
    bool operator()(int a, int b)
    {
      if (s->non_zeros[a] == s->non_zeros[b])
      {
        return s->cols[a] < s->cols[b];
      }
      else return s->non_zeros[a] < s->non_zeros[b];
    }
    
  private:
    SuccessorOrder();
    StrongSimulation_MC *s;
  };
  
  bool DecideStrongSimulation(Pair*);
  CompactMaxFlow<double> *ConstructNetwork(Pair*);
  
  int BuildRelationMap_DTMC();
  int BuildRelationMap_CTMC();
  int IterateRelation(bool);
  int IterateRelation_FirstPartition();
  
  void MakeFirstPartition();
  
  bool state_Less(int, int);
  bool state_Different(int, int);
  
  void SortSuccessors();

  MarkovChain *m;
  double *non_zeros; // list of probabilities for transitions
  int *cols;         // list of successor states
  int *row_starts;   // list indices for successor list
  int n_states;      // number of states
  
  int size_of_relation; //number of pairs in relation
  
  int partitions;
  int *partition;
  int *order;

  //The array of the relation (R).
  Pair *relation;
};

// Strong simulation for PAs
class StrongSimulation_PA : public StrongSimulation
{
friend class InitialCondition_Strong_PA;
friend class InitialCondition_Strong_CPA;

public:
  StrongSimulation_PA() {}
  virtual ~StrongSimulation_PA() {}
  
  virtual unsigned int Simulate(ProbabilisticModel*, std::set<std::pair<int,int> >*);

#ifdef WITH_VERIFIER
  bool Verify(ProbabilisticModel*, std::set<std::pair<int,int> >&, std::set<std::pair<int,int> >*, std::set<std::pair<int,int> >*);
#endif//WITH_VERIFIER
  
protected:
  // One pair of the relation. Stores the two states being related together with their maxflow problem
  class Pair
  {
  public:
    Pair() { x = 0, y = 0; }
    ~Pair() {}
    int x, y;
    Pair *next;
  };
  
  int BuildRelationMap_PA();
  int BuildRelationMap_CPA();
  
  int IterateRelation(bool);
  int IterateRelation_FirstPartition();
  
  bool DecideStrongSimulation(Pair*);
  CompactMaxFlow<double> *ConstructNetwork(int, int);
  
  void MakeFirstPartition();
  void InitializeActionMasks();
  
  bool state_Less(int, int);
  bool state_Different(int, int);
  
  void SortSuccessors();
  
  ProbabilisticAutomaton *m;
  int n_states;
  int *state_starts;
  int *actions;
  int *row_starts;
  double *non_zeros;
  int *cols;
  
  unsigned char *action_masks;
  int action_mask_pitch;
  
  int partitions, *partition, *order;
  
  Pair *relation;
  int size_of_relation;
  
  // Functor defining an order on the actions of a state (used with std::sort)
  class ActionOrder
  {
  public:
    ActionOrder(StrongSimulation_PA *_s) { s = _s; }
    ~ActionOrder() {}
    
    bool operator()(int a, int b)
    {
      int n;
      if (s->actions[a] < s->actions[b]) return true;
      else if (s->actions[a] > s->actions[b]) return false;
      if (s->row_starts[a+1] - s->row_starts[a] < s->row_starts[b+1] - s->row_starts[b]) return true;
      else if (s->row_starts[a+1] - s->row_starts[a] > s->row_starts[b+1] - s->row_starts[b]) return false;
      for (n = 0; n < s->row_starts[a + 1] - s->row_starts[a]; ++n)
      {
        if (s->non_zeros[s->row_starts[a] + n] < s->non_zeros[s->row_starts[b] + n]) return true;
        else if (s->non_zeros[s->row_starts[a] + n] > s->non_zeros[s->row_starts[b] + n]) return false;
        if (s->Label(s->cols[s->row_starts[a] + n]) < s->Label(s->cols[s->row_starts[b] + n])) return true;
        if (s->Label(s->cols[s->row_starts[a] + n]) > s->Label(s->cols[s->row_starts[b] + n])) return false;
      }
      return false;
    }
    
  private:
    ActionOrder();
    StrongSimulation_PA *s;
  };
  
  // Functor defining an order on the successors of a state (used with std::sort)
  class SuccessorOrder
  {
  public:
    SuccessorOrder(StrongSimulation_PA *_s) { s = _s; }
    ~SuccessorOrder() {}
    
    bool operator()(int a, int b)
    {
      if (s->non_zeros[a] == s->non_zeros[b])
      {
        return (s->Label(s->cols[a]) < s->Label(s->cols[b]));
      }
      else return s->non_zeros[a] < s->non_zeros[b];
    }
    
  private:
    SuccessorOrder();
    StrongSimulation_PA *s;
  };
  
  bool EqualActionSets(int s1, int s2)
  {
    return (memcmp(action_masks + (s1 * action_mask_pitch), action_masks + (s2 * action_mask_pitch), action_mask_pitch) == 0);
  }
};

#endif//STRONG_SIMULATION_H
