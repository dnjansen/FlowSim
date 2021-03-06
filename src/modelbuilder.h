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



#ifndef _MODELBUILDER_H_
#define _MODELBUILDER_H_

#include <map>
#include <set>
#include <vector>

#include "prmodel.h"

using namespace std;

// A helper class for creating probabilistic models in memory
// where the numbers of states/actions/transitions aren't known
// in advance.
//
// Repeatedly call Add() to add transitions (do not mix the two
// prototypes), then BuildMC() or BuildPA() when done to return
// objects of type MarkovChain or ProbabilisticAutomaton.
//
// Note: This class is mainly meant for testing and was not
// written with performance in mind.
class ModelBuilder
{
public:
  ModelBuilder() { transitions = 0, dist = 0; }
  ~ModelBuilder() {}
  
  void Add(int from, int to, double p)
  {
    ++transitions;
    _tmc[from].insert(make_pair(to, p));
  }
  void Add(int from, int action, int to, double p)
  {
    ++transitions;
    if (_tpa.find(from) == _tpa.end()) _tpa.insert(make_pair(from, multimap<int,pair<int,double> >()));
    _tpa[from].insert(make_pair(action,make_pair(to,p)));
    actions.insert(action);
  }
  void Clear() { transitions = 0; _tmc.clear(); _tpa.clear(); }
  MarkovChain *BuildMC(int states, double*)
  {
    int n, m;
    map<int,double>::iterator i;
    MarkovChain *mc = new MarkovChain;
    mc->n = states;
    mc->nnz = transitions;
    mc->row_starts = new int[states + 1];
    mc->cols = new int[transitions];
    mc->non_zeros = new double[transitions];
    for (n = 0, m = 0; n < states; ++n)
    {
      mc->row_starts[n] = m;
      for (i = _tmc[n].begin(); i != _tmc[n].end(); ++i)
      {
        mc->cols[m] = i->first;
        mc->non_zeros[m] = i->second;
        ++m;
      }
    }
    mc->row_starts[states] = m;
    Clear();
    return mc;
  }
  ProbabilisticAutomaton *BuildPA(int states, double*)
  {
    int n, m, l, aid;
    map<int,multimap<int,pair<int,double> > >::iterator i;
    multimap<int,pair<int,double> >::iterator k;
    pair<multimap<int,pair<int,double> >::iterator, multimap<int,pair<int,double> >::iterator> kl;
    vector<int> act, row_starts;
    set<int>::iterator a;
    ProbabilisticAutomaton *pa = new ProbabilisticAutomaton;
    pa->n = states;
    pa->nnz = transitions;
    pa->state_starts = new int[states + 1];
    pa->cols = new int[transitions];
    pa->non_zeros = new double[transitions];
    for (n = 0, m = 0, l = 0; n < states; ++n)
    {
      i = _tpa.find(n);
      pa->state_starts[n] = m;
      if (i == _tpa.end()) continue;
      for (a = actions.begin(), aid = 0; a != actions.end(); ++a, ++aid)
      {
        kl = i->second.equal_range(*a);
        if (kl.first == kl.second) continue;
        act.push_back(aid);
        row_starts.push_back(l);
        for (k = kl.first; k != kl.second; ++k)
        {
          pa->cols[l] = k->second.first;
          pa->non_zeros[l] = k->second.second;
          ++l;
        }
        ++m;
      }
    }
    row_starts.push_back(l);
    pa->state_starts[n] = m;
    pa->atable = new int[actions.size()];
    for (n = 0, a = actions.begin(); a != actions.end(); ++a) pa->atable[n++] = *a;
    pa->da = actions.size();
    pa->na = act.size();
    pa->row_starts = new int[pa->na + 1];
    pa->actions = new int[pa->na];
    for (n = 0; n < pa->na; ++n)
    {
      pa->actions[n] = act[n];
      pa->row_starts[n + 1] = row_starts[n + 1];
    }
    pa->row_starts[0] = 0;
    Clear();
    return pa;
  }
  
private:
  map<int,map<int,double> > _tmc;
  map<int,multimap<int,pair<int,double> > > _tpa;
  set<int> actions;
  int transitions, dist;
};

#endif
