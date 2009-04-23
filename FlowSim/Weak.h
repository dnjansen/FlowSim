#ifndef WEAK_SIMULATION_H
#define WEAK_SIMULATION_H

#include <set>
#include <iterator>
#include <algorithm>
#include <map>
#include <vector>
#include <queue>
#include <assert.h>
#include "Simrel.h"

// Weak simulation for DTMCs and CTMCs
class WeakSimulation_MC : public SimulationRelation
{
public:
  WeakSimulation_MC() {}
  ~WeakSimulation_MC() {}
  
  unsigned int Simulate(ProbabilisticModel*, std::set<std::pair<int,int> >*);
  
protected:
  class Pair
  {
  public:
    Pair() { x = 0, y = 0; }
    ~Pair() {}
    int x, y;
    
    void Reset() {}
  };
  
  // Functor defining an order on the successors of a state (used with std::sort)
  class SuccessorOrder
  {
  public:
    SuccessorOrder(WeakSimulation_MC *_s) { s = _s; }
    ~SuccessorOrder() {}
    
    bool operator()(int a, int b)
    {
      if (s->non_zeros[a] == s->non_zeros[b])
      {
        if (s->cols[a] == s->cols[b]) return s->Label(s->cols[a]) < s->Label(s->cols[b]);
        else return s->cols[a] < s->cols[b];
      }
      else return s->non_zeros[a] < s->non_zeros[b];
    }
    
  private:
    SuccessorOrder();
    WeakSimulation_MC *s;
  };

  int BuildRelationMap();
  void InitializeRelation();
  
  bool DecideWeakSimulation(Pair*, bool&);
  int ConstructSuccessorPartition(int s1, int s2, int **par_array, int **par_array_index);
  bool ValidFlow(int, int, std::set<double>, int = 1, int* = 0, int* = 0);
  bool CanReach(int, int, bool);
  
  int IterateRelation(Pair**);
  
  //void SortSuccessors();

  MarkovChain *m;
  double *non_zeros; // list of probabilities for transitions
  int *cols;         // list of successor states
  int *row_starts;   // list indices for successor list
  int n_states;      // number of states
  
  int size_of_relation; //number of pairs in relation

  //The array of the relation (R).
  Pair *relation;

  bool hold_hard_pairs;
  int num_easy_pairs;

#ifdef DEBUG
  unsigned long rel_space;
  FILE *simlog;
#endif//DEBUG
};

#endif//WEAK_SIMULATION_H
