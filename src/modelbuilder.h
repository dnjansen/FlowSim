#ifndef _MODELBUILDER_H_
#define _MODELBUILDER_H_

#include <map>
#include <set>
#include <vector>
#include <utility>

#include "prmodel.h"

using namespace std;

class ModelBuilder
{
public:
  ModelBuilder() { transitions = 0, dist = 0; }
  ~ModelBuilder() {}
  
  void Add(int from, int to, double p) { ++transitions; _tmc.insert(make_pair(from,make_pair(to,p))); }
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
    multimap<int,pair<int,double> >::iterator i;
    pair<multimap<int,pair<int,double> >::iterator, multimap<int,pair<int,double> >::iterator> ij;
    MarkovChain *mc = new MarkovChain;
    mc->n = states;
    mc->nnz = transitions;
    mc->row_starts = new int[states + 1];
    mc->cols = new int[transitions];
    mc->non_zeros = new double[transitions];
    for (n = 0, m = 0; n < states; ++n)
    {
      ij = _tmc.equal_range(n);
      mc->row_starts[n] = m;
      for (i = ij.first; i != ij.second; ++i)
      {
        mc->cols[m] = i->second.first;
        mc->non_zeros[m] = i->second.second;
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
    //map<int,multimap<pair<int,int>,pair<int,double> > >::iterator i;
    multimap<int,pair<int,double> >::iterator k;
    pair<multimap<int,pair<int,double> >::iterator, multimap<int,pair<int,double> >::iterator> kl;
    //pair<multimap<pair<int,int>,pair<int,double> >::iterator, multimap<pair<int,int>,pair<int,double> >::iterator> kl;
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
  multimap<int,pair<int,double> > _tmc;
  map<int,multimap<int,pair<int,double> > > _tpa;
  set<int> actions;
  int transitions, dist;
};

#endif
