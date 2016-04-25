/*****************************************************************************/
/*!
 *   Copyright 2014 David N. Jansen
 *
 *   This file is part of FlowSim.
 *
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

// The simulation algorithm in this file follows closely the publication:
// Crafa, Silvia; Ranzato, Francesco: Bisimulation and simulation algorithms
// on probabilistic transition systems by abstract interpretation.
// Formal Methods in System Design (2012) 40:356-376.
// DOI 10.1007/s10703-012-0147-3


#ifndef _STRONG_CRAFA_RANZATO_H_
#define _STRONG_CRAFA_RANZATO_H_

#include <map>
#include "Strong.h"

class StrongSimulation_CR : public StrongSimulation
{
protected:
  struct int_list {
    // This could perhaps be implemented as std::forward_list<int>, but that
    // type is only available from C++11 on.
    int_list *next;
    int el;
  };
  class in_t
  {
  public:
    in_t() { pre = NULL; Count = NULL; Remove = NULL; }
    ~in_t() {
      while ( NULL != pre ) {
        int_list *temp = pre->next;
        delete pre;
        RegisterMemFree(sizeof(*pre));
        pre = temp;
      }
      // The size of Count is not known here.
      assert(NULL == Count);
      // Remove should be empty now.
      assert(NULL == Remove);
    }
    int_list *pre; /* states that lead to this action */
    int *Count; /* counts the number of distributions that are Rd-successors
                        and also have this in-action */
    /* Another possibility, using perhaps less memory, would be:
    - for state x an array of size ((ProbabilisticAutomaton*)model)->da with
      pointers. (Or rather one big array of size nStates * ...da).
    - The a'th entry in this array of pointers points to an array of size
      nDistributions with integers. (If there is no Count(x,a,...) then this
      pointer would be NULL.)
        These arrays will be large. So we make them smaller by the following
      trick: First, number all distributions that have a as in-action. Then,
      use these numbers as index in the array. The array size can be reduced to
      what is necessary.
        I did try to replace Count by std::map<int,int>. However, in self-
      stabilising 12, this meant that the data structures would use 374871132
      instead of 381588480 bytes, a saving of about 1.76% only.
      */
    int_list *Remove; /* states in Remove_a(d) */
  };
  class Distribution
  {
  public:
    Distribution() { in = NULL; }
    ~Distribution() {
      // The size of in[] is not known here, we have to delete it elsewhere.
      assert(NULL == in);
      for ( std::map<int,CompactMaxFlow<double>*>::iterator e=Rd_inv.begin() ;
                  e != Rd_inv.end() ; e++ )
      {
        delete e->second;
        e->second = NULL;
      }
      CompactMaxFlow<double>::RegisterMemFree(Rd_inv.size()
                  * sizeof(CompactMaxFlow<double>));
      RegisterMemFree(Rd_inv.size() *
                (sizeof(std::pair<int,CompactMaxFlow<double>*>)+MAP_OVERHEAD));
      Rd_inv.clear();
    }
    int row_starts, row_ends; /* start (inclusive) and end (exclusive) index in the cols and non_zeros arrays */
    in_t *in; /* pointer to an array of in-actions. Actually the
    array contains an element for every action; however, only if
    in_action[a].pre != NULL, the action is a real in-action. All these arrays
    together need memory |Act|*|Distr|, which fits in the memory complexity. */
    std::map<int,CompactMaxFlow<double>*> Rd_inv; /* R-predecessors of this distribution, each with a witnessing flow network */
  };
  Distribution *distributions; /* array of all (distinct) distributions used */
  int *distr_index; /* indicates for each state which entry in distributions[] is the right one */
  RelationMap Rdmap;

  struct Listener_t {
    Listener_t *next;
    int d, e;
  };
  Listener_t **Listener;
  /* I tried what would happen if one changed the data type to
  std::map<unsigned int,Listener_t *>. Using that type reduces the memory use
  (for Listener) from 602358648 to 523096472 bytes for self-stabilising 12, a
  reduction of 13%. However, the total memory consumption is reduced by 1.85%
  only. I decided to keep the simple type. */

public:
  StrongSimulation_CR() {
    distributions = NULL;
    distr_index = NULL;
    Listener = NULL;
    model = NULL;
    Deleted = NULL;
    action_masks = NULL;
    state_starts = NULL;
    actions = NULL;
    row_starts = NULL;
    row_ends = NULL;
    cols = NULL;
    non_zeros = NULL;
  }
  virtual ~StrongSimulation_CR() {
    assert(NULL == distributions);
    assert(NULL == distr_index);
    assert(NULL == Listener);
  }
  
  // Main interface
  unsigned int Simulate(ProbabilisticModel*, std::set<std::pair<int,int> >*);

protected:
  // Member variables
  ProbabilisticModel *model;
  int nStates, nDistributions, nActions;
  bool probStable, preStable;
  struct Deleted_t {
    Deleted_t *next;
    int t, w;
  } *Deleted;
  unsigned int *action_masks;
  int action_mask_pitch; // = number of elements in action_masks[] reserved per state.

  // Model structures
  int *state_starts;
  int *actions;
  int *row_starts;
  int *row_ends;
  int *cols;
  double *non_zeros;

  // Initialization
  void Initialize();
  void InitializeActionMasks();

  // Auxiliaries
  inline bool Act_Subseteq(int s1, int s2)
  {
    if ( NULL == action_masks ) return true;
    s1 *= action_mask_pitch;
    s2 *= action_mask_pitch;
    for (int i = 0; i < action_mask_pitch; ++i)
    {
      if ( action_masks[s1 + i] & ~action_masks[s2 + i] )
        return false;
    }
    return true;
  }

  // The algorithm
  void PBis();

  // Subroutines of the algorithm
  bool preStabilize();
  bool probStabilize();
  CompactMaxFlow<double>* Init_SMF(int d, int e);
  bool SMF(int d, int e, int t, int w);
};

#endif//_STRONG_CRAFA_RANZATO_H_
