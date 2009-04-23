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



#ifndef COMPACTMAXFLOW_CC
#define COMPACTMAXFLOW_CC

using namespace std;

#include "compactmaxflow.h"

#ifdef DEBUG
template <typename _T> unsigned long CompactMaxFlow<_T>::global_space = 0;
template <typename _T> unsigned long CompactMaxFlow<_T>::global_instances = 0;
template <typename _T> unsigned long CompactMaxFlow<_T>::global_space_peak = 0;
template <typename _T> unsigned long CompactMaxFlow<_T>::global_times_invoked = 0;
template <typename _T> unsigned long CompactMaxFlow<_T>::global_p_inv_fails = 0;
template <typename _T> unsigned long CompactMaxFlow<_T>::global_sig_arc_fails = 0;
template <typename _T> unsigned int  CompactMaxFlow<_T>::min_complexity = 0;
template <typename _T> unsigned int  CompactMaxFlow<_T>::max_complexity = 0;
#endif//DEBUG

#ifndef DISABLE_FP_APPROXIMATION
// Unless floating point approximation is disabled, define the static variable containing
// the desired precision. By default, this value is zero, which means that absolute
// equality/inequality are required. If this value is changed to something non-zero, values
// which are less than ::precision apart from each other are considered equal. This is
// necessary to correctly simulate some models which are stored with insufficient precision
// to accurately represent the model.
template <typename _T> _T CompactMaxFlow<_T>::precision = 0.0;
double CompactFeasibleFlow::precision = 0.0;
#endif//DISABLE_FP_APPROXIMATION

template <typename _T> CompactMaxFlow<_T>::CompactMaxFlow()
{
  n1 = 0, n2 = 0, set1 = 0, set2 = 0, arcs = 0, storage = 0, valid = false;
#ifdef DEBUG
  space_usage = 0;
  global_space += sizeof(*this);
  if (global_space > global_space_peak) global_space_peak = global_space;
#endif//DEBUG
}

template <typename _T> CompactMaxFlow<_T>::~CompactMaxFlow()
{
  _FreeInternals();
#ifdef DEBUG
  global_space -= sizeof(*this);
#endif//DEBUG
}

// Private routine: free memory used by internal structures
template <typename _T> void CompactMaxFlow<_T>::_FreeInternals()
{
  if (storage)
  {
    delete[] storage;
    storage = 0;
#ifdef DEBUG
    --global_instances;
    global_space -= space_usage;
  }
  space_usage = 0;
#else//DEBUG
  }
#endif//DEBUG
  valid = false;
  n_arcs = 0;
}

// Create a network from sparse matrix data. If known_result is true after the call, the
// return value represents the result of the simulation, otherwise the function returns true
// if the network was successfully created and false if there was a problem.
template <typename _T> bool CompactMaxFlow<_T>::CreateNetwork(int *successors, _T *probabilities, RelationMap *rmap,
                                        int l0, int lc, int r0, int rc, bool &known_result, unsigned long optflags)
{
  return CreateNetwork(successors + l0, successors + r0, probabilities + l0, probabilities + r0, rmap, lc, rc, known_result, optflags);
}

template <typename _T> bool CompactMaxFlow<_T>::CreateNetwork(int *sl, int *sr, _T *pl, _T *pr, RelationMap *rmap,
                                        int lc, int rc, bool &known_result, unsigned long/*optflags*/)
{
  int n, m, size;
  arc *pa, **ppa;
  
#if defined(OPT_SIGNIFICIANT_ARC) || defined(OPT_P_INVARIANT)
  _T *relsum1, *relsum2, lsum, rsum;
  int l, r;
#endif
  
  _FreeInternals();

  // No successors at all is simulated by anything
  if (lc == 0)
  {
    valid = true;
    n_arcs = 0;
    known_result = true;
    return true;
  }
  // No successors at all cannot simulate anything
  if (rc == 0)
  {
    known_result = true;
    return false;
  }
  
  assert(sl && sr && pl && pr && rmap);
  
  n1 = lc;
  n2 = rc;
  
  // Determine number of arcs
  for (n = 0, n_arcs = 0; n < lc; ++n)
  {
    for (m = 0; m < rc; ++m)
    {
      if (sl[n] == sr[m] || (*rmap)(sl[n], sr[m])) ++n_arcs;
    }
  }
  
  if (n_arcs == 0)
  {
    known_result = true;
    return false;
  }
  
  // Allocate one memory area for storage and store pointers to individual sections of that block
  if (storage) delete [] storage;
  size = (sizeof(node) * (n1 + n2)) + (sizeof(arc) * n_arcs) + (sizeof(arc*) * ((n_arcs * 2) + n1 + n2));
  storage = new unsigned char[size];
  set1 = (node*)storage;
  set2 = set1 + n1;
  arcs = (arc*)(set2 + n2);
#ifdef DEBUG
  space_usage = size;
  global_space += space_usage;
  ++global_instances;
  if (global_space > global_space_peak) global_space_peak = global_space;
#endif//DEBUG
  
  memset(storage, 0, size);
  
//  if (optflags & (OPT_SIGNIFICIANT_ARC | OPT_P_INVARIANT))
#if defined(OPT_SIGNIFICIANT_ARC) || defined(OPT_P_INVARIANT)
  {
    relsum1 = new _T[n1 + n2];
    relsum2 = relsum1 + n1;
    lsum = 0, rsum = 0;
  }
#endif
  
  // Fill in node data (for left set only) and arc data and connect arcs to nodes on the left
  pa = arcs;
  ppa = (arc**)(arcs + n_arcs);
  //if (optflags & (OPT_SIGNIFICIANT_ARC | OPT_P_INVARIANT))
#if defined(OPT_SIGNIFICIANT_ARC) || defined(OPT_P_INVARIANT)
  {
    for (n = 0; n < n1; ++n)
    {
      set1[n].arcs = ppa;
      set1[n].excess = pl[n];
      set1[n].id = sl[n];
      relsum1[n] = 0;
      lsum += pl[n];
      for (m = 0; m < n2; ++m)
      {
        if ((*rmap)(sl[n], sr[m]) || sl[n] == sr[m])
        {
          relsum1[n] += pr[m];
          *ppa = (pa++);
          (*ppa)->tail = set1 + n;
          (*ppa)->head = set2 + m;
          (*ppa)->flow = 0;
          ++set1[n].label;
          ++set2[m].label;
          ++ppa;
        }
      }
      *ppa = 0;
      ++ppa;
      
      // There is excess on the left side that has no outgoing arc; this can't work
      if (set1[n].label == 0 && _Tless(0, set1[n].excess))
      {
        delete [] relsum1;
        _FreeInternals();
        known_result = true;
        return false;
      }
    }
  }
#else
  //else//optflags & (OPT_SIGNIFICIANT_ARC | OPT_P_INVARIANT)
  {
    for (n = 0; n < n1; ++n)
    {
      set1[n].arcs = ppa;
      set1[n].excess = pl[n];
      set1[n].id = sl[n];
      for (m = 0; m < n2; ++m)
      {
        if ((*rmap)(sl[n], sr[m]) || sl[n] == sr[m])
        {
          *ppa = (pa++);
          (*ppa)->tail = set1 + n;
          (*ppa)->head = set2 + m;
          (*ppa)->flow = 0;
          ++set1[n].label;
          ++set2[m].label;
          ++ppa;
        }
      }
      *ppa = 0;
      ++ppa;
      
      // There is excess on the left side that has no outgoing arc; this can't work
      if (set1[n].label == 0 && _Tless(0, set1[n].excess))
      {
        _FreeInternals();
        known_result = true;
        return false;
      }
    }
  }
#endif
  
  // Fill in node data (right set) and connect arcs to nodes on the right
  for (n = 0; n < n2; ++n)
  {
    set2[n].excess = pr[n];
    set2[n].arcs = ppa;
    set2[n].id = sr[n];
    pa = arcs;
    while (set2[n].label && pa != arcs + n_arcs) // Label is used to count the number of arcs connected
    {                                            // to this node at this point in the code
      if (pa->head == set2 + n)
      {
        --set2[n].label;
        (*ppa) = pa;
        ++ppa;
      }
      ++pa;
    }
    
    *ppa = 0;
    ++ppa;
    //if (optflags & (OPT_SIGNIFICIANT_ARC | OPT_P_INVARIANT))
#if defined(OPT_SIGNIFICIANT_ARC) || defined(OPT_P_INVARIANT)
    {
      relsum2[n] = 0;
      rsum += pr[n];
      for (m = 0; m < lc; ++m)
      {
        if ((*rmap)(sl[m], sr[n]) || sl[m] == sr[n]) relsum2[n] += pl[m];
      }
    }
#endif
  }
  
  //if (optflags & (OPT_SIGNIFICIANT_ARC | OPT_P_INVARIANT))
#if defined(OPT_SIGNIFICIANT_ARC) || defined(OPT_P_INVARIANT)
  {
    aux = rsum - lsum;
    if (_Tless(aux, 0))
    {
#ifdef DEBUG
      ++global_p_inv_fails;
#endif//DEBUG
      delete [] relsum1;
      _FreeInternals();
      known_result = true;
      return false;
    }
    
    // Determine significiant arcs
    for (n = 0; n < n_arcs; ++n)
    {
      l = arcs[n].tail - set1;
      r = arcs[n].head - set2;
      //if ((optflags & OPT_P_INVARIANT) && (_Tless(relsum1[l], set1[l].excess) || _Tless(relsum2[r] + aux, set2[r].excess)))
#ifdef OPT_P_INVARIANT
      if (_Tless(relsum1[l], set1[l].excess) || _Tless(relsum2[r] + aux, set2[r].excess))
      {
#ifdef DEBUG
        ++global_p_inv_fails;
#endif//DEBUG
        delete [] relsum1;
        _FreeInternals();
        known_result = true;
        return false;
      }
#endif
#ifdef OPT_SIGNIFICIANT_ARC
      //if (optflags & OPT_SIGNIFICIANT_ARC)
      {
        (arcs + n)->significiant = _Tless(relsum1[l]       - set2[r].excess, set1[l].excess)
                                || _Tless(relsum2[r] + aux - set1[l].excess, set2[r].excess);
      }
#endif
    }
    
    delete [] relsum1;
  }
#endif

  tpl = pl;
  tpr = pr;

  valid = true;
  incomplete_flow = false;
  
  for (n = 0; n < n1; ++n) set1[n].label = 1;
  for (n = 0; n < n2; ++n) set2[n].label = 0;
  
#ifdef DEBUG
  complexity = -1;
#endif//DEBUG
  
  known_result = false;
  return true;
}

// Update the network (i.e. delete arcs which are no longer in the relation) based on
// a relation map generated by a strong simulation problem.
template <typename _T> bool CompactMaxFlow<_T>::UpdateNetwork(RelationMap *rmap, bool sigarc)
{
  int n;
  
  if (!storage || !valid) return false;
  if (incomplete_flow) return true;
  
  for (n = 0; n < n_arcs; ++n)
  {
    if ((arcs + n)->tail && (arcs + n)->head)
    {
      if ((arcs + n)->tail->id == (arcs + n)->head->id) continue;
      if (!(*rmap)((arcs + n)->tail->id, (arcs + n)->head->id)) _RemoveArc(arcs + n, sigarc);
    }
  }
  
  return true;
}

// Remove a particular arc (private routine).
// Updates significiance flags on arcs if this optimization is enabled.
template <typename _T> void CompactMaxFlow<_T>::_RemoveArc(arc *a, bool/*sigarc*/)
{
  //if (sigarc)
#ifdef OPT_SIGNIFICIANT_ARC
  int l, r, n;
  _T relsum1 = 0, relsum2 = aux;
  l = a->tail - set1;
  r = a->head - set2;
  
  if (a->significiant)
  {
    _FreeInternals();
    incomplete_flow = true;
#ifdef DEBUG
    ++global_sig_arc_fails;
#endif//DEBUG
    return;
  }
#endif
  
  a->tail->excess += a->flow;
  a->head->excess += a->flow;
  a->tail = 0;
  a->head = 0;
  a->flow = 0;
    
#ifdef OPT_SIGNIFICIANT_ARC
  for (n = 0; *((set1 + l)->arcs + n); ++n)
  {
    if ((*((set1 + l)->arcs + n))->head) relsum1 += tpr[((*((set1 + l)->arcs + n))->head - set2)];
  }
  for (n = 0; *((set1 + l)->arcs + n); ++n)
  {
    if ((*((set1 + l)->arcs + n))->head && !(*((set1 + l)->arcs + n))->significiant)
    {
      (*((set1 + l)->arcs + n))->significiant = _Tless(relsum1 - tpr[((*((set1 + l)->arcs + n))->head - set2)], tpl[l]);
    }
  }
  for (n = 0; *((set2 + r)->arcs + n); ++n)
  {
    if ((*((set2 + r)->arcs + n))->tail) relsum2 += tpl[((*((set2 + r)->arcs + n))->tail - set1)];
  }
  for (n = 0; *((set2 + r)->arcs + n); ++n)
  {
    if ((*((set2 + r)->arcs + n))->tail && !(*((set2 + r)->arcs + n))->significiant)
    {
      (*((set2 + r)->arcs + n))->significiant = _Tless(relsum2 - tpl[((*((set2 + r)->arcs + n))->tail - set1)], tpr[r]);
    }
  }
#endif
}

// Determine if the flow in this network satisfies the condition for strong
// simulation
template <typename _T> bool CompactMaxFlow<_T>::IsFlowTotal(bool restart)
{
  int n, m, mlast, a, min_label, l;
  bool keep_going = false;;
  node *h;
  
#ifdef DEBUG
  unsigned int cmpl = 0;
#endif//DEBUG
  
  if (incomplete_flow || !valid) return false;
  
  if (!storage) return true; // Valid network with NULL storage encodes "always true"
  
#ifdef DEBUG
  ++global_times_invoked;
#endif//DEBUG

  // Restart computation from scratch
  if (restart)
  {
    for (n = 0; n < n_arcs; ++n) arcs[n].flow = 0;
    for (n = 0; n < n1; ++n) set1[n].label = 0, set1[n].excess = tpl[n];
    for (n = 0; n < n2; ++n) set2[n].label = 0, set2[n].excess = tpr[n];
  }
  
  // Initial push: push from left nodes to right nodes but only to saturation, not beyond
  for (n = 0; n < n1; ++n)
  {
    // Iterate through arcs connected to this node, pushing down excess if the head is not saturated
    for (m = 0, a = 0, mlast = -1; *((set1 + n)->arcs + m); ++m)
    {
      if ((*((set1 + n)->arcs + m))->head != 0)
      {
        ++a; // Count only active (not deleted) arcs
        if (_Tleq((set1 + n)->excess, (*((set1 + n)->arcs + m))->head->excess) && _Tless(0, (set1 + n)->excess))
        {
          h = (*((set1 + n)->arcs + m))->head;
          // Push entire excess
          (*((set1 + n)->arcs + m))->head->excess -= (set1 + n)->excess;
          (*((set1 + n)->arcs + m))->flow += (set1 + n)->excess;
          (set1 + n)->excess = 0;
          
          // Even though we're about to leave the loop, we have to keep counting active (non-deleted) arcs
          while (*((set1 + n)->arcs + (++m))) if ((*((set1 + n)->arcs + m))->head != 0) ++a;
          break;
        }
        else if (_Tless(0, (*((set1 + n)->arcs + m))->head->excess) && _Tless(0, (set1 + n)->excess))
        {
          h = (*((set1 + n)->arcs + m))->head;
          // Push part of the excess to saturate arc head
          (*((set1 + n)->arcs + m))->flow += (*((set1 + n)->arcs + m))->head->excess;
          (set1 + n)->excess -= (*((set1 + n)->arcs + m))->head->excess;
          (*((set1 + n)->arcs + m))->head->excess = 0;
          mlast = m; // Remember this arc
        }
        if (mlast == -1) mlast = m;
      }
    }
    
    // Update label.
    if (a == 0 && _Tless(0, (set1 + n)->excess))
    {
      _FreeInternals();
      incomplete_flow = true;
      return false;
    }
    else if (a == 1)
    {
      // This node has only one arc connected; all excess must be pushed through this arc
      // and the label is set to maximum so that no excess can be pushed back to this node
      (set1 + n)->label = n1 + n2;
      if (_Tless(0, (set1 + n)->excess))
      {
        (*((set1 + n)->arcs + mlast))->head->excess -= (set1 + n)->excess;
        (*((set1 + n)->arcs + mlast))->flow += (set1 + n)->excess;
        (set1 + n)->excess = 0;
      }
    }
    else if (_Tless(0, set1[n].excess)) keep_going = true;
  }

  // Secondary loop: push & relabel. Find nodes on the right side which are over-saturated
  // and try to redistribute that excess.
  while (true)
  {
#ifdef DEBUG
    ++cmpl;
#endif//DEBUG
    //keep_going = false;
    for (n = 0; n < n2; ++n)
    {
      while (_Tless((set2 + n)->excess, 0))
      {
        keep_going = true;
        for (m = 0, a = -1, min_label = n1 + n2; *((set2 + n)->arcs + m); ++m)
        {
          if ((*((set2 + n)->arcs + m))->tail && _Tless(0, (*((set2 + n)->arcs + m))->flow))
          {
            l = ((*((set2 + n)->arcs + m))->tail)->label;
            if (l < min_label) min_label = l, a = m;
          }
        }
        if (a == -1 || min_label == n1 + n2)
        {
          // Oversaturated node but no arcs with a positive label found
          incomplete_flow = true;
          _FreeInternals();
          return false;
        }
        (set2 + n)->label = min_label + 1;
        if (_Tless(-(set2 + n)->excess, (*((set2 + n)->arcs + a))->flow))
        {
          (*((set2 + n)->arcs + a))->flow += (set2 + n)->excess;
          ((*((set2 + n)->arcs + a))->tail)->excess -= (set2 + n)->excess;
          (set2 + n)->excess = 0;
        }
        else
        {
          (set2 + n)->excess += (*((set2 + n)->arcs + a))->flow;
          ((*((set2 + n)->arcs + a))->tail)->excess += (*((set2 + n)->arcs + a))->flow;
          (*((set2 + n)->arcs + a))->flow = 0;
        }
      }
    }
    
    if (!keep_going) break;
    
    for (n = 0; n < n1; ++n)
    {
      while (_Tless(0, (set1 + n)->excess))
      {
        for (m = 0, a = -1, min_label = n1 + n2; *((set1 + n)->arcs + m); ++m)
        {
          if ((*((set1 + n)->arcs + m))->head)
          {
            l = ((*((set1 + n)->arcs + m))->head)->label;
            if (l < min_label) min_label = l, a = m;
          }
        }
        if (a == -1 || min_label == n1 + n2)
        {
          // Oversaturated node but no arcs with a positive label found
          incomplete_flow = true;
          _FreeInternals();
          return false;
        }
        (set1 + n)->label = min_label + 1;
        (*((set1 + n)->arcs + a))->flow += (set1 + n)->excess;
        ((*((set1 + n)->arcs + a))->head)->excess -= (set1 + n)->excess;
        (set1 + n)->excess = 0;
      }
    }
    
    keep_going = false;
  }

#ifdef DEBUG
  // Record complexity
  if (complexity == -1) complexity = cmpl;
  if (cmpl < min_complexity) min_complexity = cmpl;
  if (cmpl > max_complexity) max_complexity = cmpl;
#endif//DEBUG
  
  return true;
}

#ifdef DEBUG
// Output the current state of the network
template <typename _T> void CompactMaxFlow<_T>::Dump(const char *s)
{
  int n, m;
  node *h;
  printf("=== %s ===\n", s);
  if (!valid || !storage || incomplete_flow)
  {
    printf("valid=%s storage=%p incomplete_flow=%s\n", (valid ? "true" : "false"), (void*)storage, (incomplete_flow ? "true" : "false"));
    return;
  }
  for (n = 0; n < n1; ++n)
  {
    printf("%4d[%2d] %#.10f\n", (set1 + n)->id, (set1 + n)->label, (set1 + n)->excess);
    for (m = 0; *((set1 + n)-> arcs + m); ++m)
    {
      h = (*((set1 + n)-> arcs + m))->head;
      if (h)
      {
        printf("                 -- %f (%c) --> %4d[%2d] %#.10f\n", (*((set1 + n)-> arcs + m))->flow,
#ifdef OPT_SIGNIFICIANT_ARC
        ((*((set1 + n)-> arcs + m))->significiant ? '!' : ' '),
#else//OPT_SIGNIFICIANT_ARC
        ' ',
#endif//OPT_SIGNIFICIANT_ARC
        h->id, h->label, h->excess);
      }
    }
  }
  printf("\n");
}
#endif//DEBUG

CompactFeasibleFlow::CompactFeasibleFlow()
{
  //lset = 0, rset = 0, arcs = 0;
}

CompactFeasibleFlow::~CompactFeasibleFlow()
{
  _FreeInternals();
}

void CompactFeasibleFlow::_FreeInternals()
{
  /*if (lset) delete [] lset;
  if (rset) delete [] rset;
  if (arcs) delete [] arcs;
  lset = 0, rset = 0, arcs = 0;*/
}

// Transform feasible flow problem to maximum flow problem
bool CompactFeasibleFlow::CreateNetwork(int *successors, double *probabilities,
  RelationMap *rmap, int l0, int lc, int r0, int rc, set<int> &MU1, set<int> &MU2, unsigned long)
{
  int n, m;
  
  arcs.clear();
  sourcecap = 0.0, sinkcap = 0.0;
  sink = 4 + lc + rc;
  
  /*printf("  s'  ->  1\n");
  printf("  s   ->  2\n");
  printf("  t   ->  3\n");
  for (n = 0; n < lc; ++n) printf("L %-2d%c -> %2d\n", successors[l0 + n], (MU1.find(successors[l0 + n]) == MU1.end() ? ' ' : '*'), 4 + n);
  for (n = 0; n < rc; ++n) printf("R %-2d%c -> %2d\n", successors[r0 + n], (MU2.find(successors[r0 + n]) == MU2.end() ? ' ' : '*'), 4 + lc + n);
  printf("  t'  -> %2d\n", sink);*/
  
  // Old sink to old source arc
  arcs.push_back(make_pair(3, make_pair(-1.0, 2)));
  
  // Create infinite capacity arcs between the two parts of the bipartite network
  for (n = 0; n < lc; ++n)
  {
    for (m = 0; m < rc; ++m)
    {
      if ((*rmap)(successors[l0 + n], successors[r0 + m]) || (successors[l0 + n] == successors[r0 + m]))
        arcs.push_back(make_pair(4 + n, make_pair(-1.0, 4 + lc + m)));
    }
  }
  
  // Add arcs from new source to MU1 and from old source to not(MU1)
  for (n = 0; n < lc; ++n)
  {
    if (MU1.find(successors[l0 + n]) == MU1.end())
    {
      arcs.push_back(make_pair(2, make_pair(probabilities[l0 + n], 4 + n)));
    }
    else
    {
      sourcecap += probabilities[l0 + n];
      arcs.push_back(make_pair(1, make_pair(probabilities[l0 + n], 4 + n)));
    }
  }
  
  // Add arcs from MU2 to new sink and from not(MU2) to old sink
  for (n = 0; n < rc; ++n)
  {
    if (MU2.find(successors[r0 + n]) == MU2.end())
    {
      arcs.push_back(make_pair(4 + lc + n, make_pair(-probabilities[r0 + n], 3)));
    }
    else
    {
      sinkcap += probabilities[r0 + n];
      arcs.push_back(make_pair(4 + lc + n, make_pair(-probabilities[r0 + n], sink)));
    }
  }
  
  // Add arcs to/from old source and sink
  arcs.push_back(make_pair(1, make_pair(-sinkcap, 3)));
  arcs.push_back(make_pair(2, make_pair(sourcecap, sink)));
  
  return true;
  
#if 0
  Node *n;
  Arc *a;
  int i, j;
  
  n1 = MU1.size();
  n2 = MU2.size();
  
  _FreeInternals();
  
  lset = new Node[n1];
  rest = new Node[n2];
  CompactFeasibleFlow<_T>::probabilities = probabilities;
  
  laux = 0.0, raux = 0.0, rneed = 0.0;
  
  for (n = lset, i = l0; i < lc; ++i)
  {
    if (MU1.find(successors[i]) == MU1.end()) laux += probabilities[i];
    else
    {
      n->excess = probabilities[i];
      n->id = i;
      n->label = 0;
      n->req = 0.0;
      n->aux = 0.0;
      for (j = 0; j < rc; ++j)
      {
        if ((successors[i] == successors[j] || (*rmap)(successors[i], successors[j]))
            && MU2.find(successors[j]) == MU2.end()) n->aux += probabilities[j];
      }
      ++n;
    }
  }
  
  for (n = rset, i = r0; i < rc; ++i)
  {
    if (MU2.find(successors[i]) == MU2.end()) raux += probabilities[i];
    else
    {
      rneed += probabilities[i];
      n->excess = probabilities[i];
      n->id = i;
      n->label = 0;
      n->req = 0.0;
      n->aux = 0.0;
      for (j = 0; j < lc; ++j)
      {
        if ((successors[j] == successors[i] || (*rmap)(successors[j], successors[i]))
            && MU1.find(successors[j]) == MU1.end()) n->aux += probabilities[j];
      }
      ++n;
    }
  }
  
  n_arcs = 0;
  for (i = 0; i < n1; ++i)
  {
    for (j = 0; j < n2; ++j)
    {
      if (successors[lset[i].id] == successors[rset[j].id] || (*rmap)(successors[lset[i].id], successors[rset[j].id]))
      {
        ++n_arcs;
        ++lset[i].label;
        ++rset[i].label;
      }
    }
  }
  
  arcs = new Arc[n_arcs];
  
  for (i = 0; i < n1; ++i)
  {
    lset[i].arcs = new Arc*[lset[i].label + 1];
    lset[i].arcs[lset[i].label] = 0;
    lset[i].label = 1;
  }
  
  for (i = 0; i < n2; ++i)
  {
    rset[i].arcs = new Arc*[rset[i].label + 1];
    rset[i].arcs[rset[i].label] = 0;
    rset[i].label = 0;
  }
  
  for (a = arcs, i = 0; i < n1; ++i)
  {
    for (j = 0; j < n2; ++j)
    {
      if (successors[lset[i].id] == successors[rset[j].id] || (*rmap)(successors[lset[i].id], successors[rset[j].id]))
      {
        a->tail = lset[i];
        a->head = rset[j];
        a->flow = 0.0;
        lset[i].arcs[lset[i].label++] = a;
        rset[i].arcs[rset[j].label++] = a;
        ++a;
      }
    }
  }
  
  return true;
#endif
}

#ifndef HIPR_PROGRAM
#define HIPR_PROGRAM       "./hipr/hi_pr"
#endif

// Delegate feasible flow check to external hi_pr program
bool CompactFeasibleFlow::IsFlowFeasible(set<double> &params, set<double>::iterator *feasible_param)
{
  set<double>::iterator i;
  vector<pair<int, pair<double, int> > >::iterator vi;
  int fd, from, to;
  double cap;
  FILE *f;
  char tmp[PATH_MAX], hipr[PATH_MAX];
  
  for (i = params.begin(); i != params.end(); ++i)
  {
    strcpy(&tmp[0], "/tmp/feasibleflow.XXXXXXXXX");
    fd = mkstemp(&tmp[0]);
    f = fdopen(fd, "w");
    
    fprintf(f, "p max %u %u\nn %d s\nn %d t\n", sink, arcs.size(), 1, sink);
    
    for (vi = arcs.begin(); vi != arcs.end(); ++vi)
    {
      from = vi->first;
      to = vi->second.second;
      cap = (vi->second.first < 0.0 ? -vi->second.first * (*i) : vi->second.first);
      fprintf(f, "a %d %d %.20f\n", from, to, cap);
    }
    
    fclose(f);
    
    sprintf(&hipr[0], "%s %.20e < %s", HIPR_PROGRAM, precision, &tmp[0]);
    f = popen(&hipr[0], "r");
    if (!f)
    {
      fprintf(stderr, "Couldn't run hipr program! I tried '%s'\n", HIPR_PROGRAM);
      exit(-1);
    }
    
    while (1)
    {
      assert(fgets(&hipr[0], PATH_MAX, f));
      if (!strncmp(&hipr[0], "c flow:", 7))
      {
        cap = strtod(&hipr[7], 0);
        break;
      }
    }
    
    pclose(f);
    
    if (_Teq(cap, sourcecap + ((*i) * sinkcap)))
    {
      if (feasible_param) *feasible_param = i;
      return true;
    }
  }
  
  return false;
  
#if 0
  int i, j, min_label;
  Arc **a;
  _T minparam, maxparam, delta;
  set<_T>::iterator param_it;
  set<Node*> lact, ract;
  set<Node*>::iterator n;
  
  if (!lset || !rset || !arcs) return params.end();
  if (n1 == 0 && n2 == 0) return params.begin();
  if (n_arcs == 0) return params.end();
  
  // Determine the smallest parameter that can possibly satisfy feasible flow using P-Invariant ideas
  for (minparam = 0.0, i = 0; i < n1; ++i)
  {
    if (_Teq(probabilities[lset[i].id], 0.0)) continue;
    
    lset[i].req = lset[i].aux;
    for (a = lset[i].arcs; *a; ++a)
    {
      if (!(*a)->tail || !(*)->head) continue;
      lset[i].req += probabilities[(*a)->head->id];
    }
    
    if (_Teq(lset[i].req, 0.0)) return params.end();
    else if (_Tless(minparam, lset[i].excess / lset[i].req)) minparam = lset[i].excess / lset[i].req;
  }
  
  // Determine the largest parameter that can possibly satisfy feasible flow using P-Invariant ideas
  for (maxparam = -1.0, i = 0; i < n2; ++i)
  {
    if (_Teq(probabilities[rset[i].id], 0.0)) continue;
    
    rset[i].req = rset[i].aux;
    for (a = rset[i].arcs; *a; ++a)
    {
      if (!(*a)->tail || !(*)->head) continue;
      rset[i].req += probabilities[(*a)->tail->id];
    }
    
    if (_Teq(rset[i].req, 0.0)) return params.end();
    else if (maxparam < 0.0 || _Tless(rset[i].excess / rset[i].req, maxparam)) maxparam = rset[i].excess / rset[i].req;
  }
  
  // Initialize excess
  for (i = 0; i < n1; ++i) lset[i].excess = 0.0;
  for (i = 0; i < n2; ++i) rset[i].excess = 0.0;
  
  for (param_id = params.begin(); param_it != params.end(); ++param_it)
  {
    // Ignore parameters outside of the established bounds
    if (_Tless(*param_it, minparam)) continue;
    if (_Tless(maxparam, *param_it)) break;
    
    lact.clear();
    ract.clear();
    
    // Initialize labels and required outflow on left side nodes
    for (i = 0; i < n1; ++i)
    {
      lset[i].req = probabilities[lset[i].id] - (*param_it * lset[i].aux);
      if (_Tless(lset[i].req, 0.0)) lset[i].req = 0.0;
      if (_Tless(lset[i].excess, lset[i].req)) lact.insert(lset + i);
      lset[i].label = 1;
    }
    
    // Initialize labels and required inflow on right side nodes
    for (i = 0; i < n2; ++i)
    {
      rset[i].req = (*param_it * probabilities[rset[i].id]) - rset[i].aux;
      if (_Tless(rest[i].req, 0.0)) rest[i].req = 0.0;
      rset[i].label = 0;
    }
    
    // Preflow
    while (1)
    {
      for (n = lact.begin(); n != lact.end(); ++n)
      {
        min_label = (*n)->label;
        for (a = (*n)->arcs; *a && _Tless((*n)->excess, (*n)->req); ++a)
        {
          if (!(*a)->tail || !(*a)->head) continue;
          if ((*a)->head->label < (*n)->label)
          {
            //delta = (_Tless((*n)->req) ? : );
          }
        }
      }
    }
  }
  
  return params.end();
#endif
}

#if 0
template <typename _T> bool CompactMaxFlow_CheckValid<_T>::CreateNetwork(int *successors, _T *probabilities, RelationMap *rmap,
                                        int l0, int lc, int r0, int rc, bool &known_result, unsigned long optflags)
{
  return CreateNetwork(successors + l0, successors + r0, probabilities + l0, probabilities + r0, rmap, lc, rc, known_result, optflags);
}

template <typename _T> bool CompactMaxFlow_CheckValid<_T>::CreateNetwork(int *sl, int *sr, _T *pl, _T *pr, RelationMap *rmap,
                                        int lc, int rc, bool &known_result, unsigned long optflags)
{
  int n, m, k, *succ;
  bool res;
  
  if (cap) delete [] cap;
  
  orig = pr;
  if (MU1 && MU2)
  {
    cap = new _T[MU1->size() + MU2->size()];
    succ = new int [MU1->size() + MU2->size()];
    if ((int)MU1->size() == lc)
    {
      memcpy(succ, sl, lc * sizeof(int));
      memcpy(cap, pl, lc * sizeof(_T));
      m = lc;
    }
    else for (n = 0, m = 0; n < lc; ++n)
    {
      if (MU1->find(sl[n]) != MU1->end())
      {
        succ[m] = sl[n];
        cap[m] = pl[n];
        ++m;
      }
    }
    if ((int)MU2->size() == rc)
    {
      memcpy(succ + lc, sr, rc * sizeof(int));
      memcpy(cap + lc, pr, rc * sizeof(_T));
      k = m + rc;
    }
    else for (n = 0, k = m; n < rc; ++n)
    {
      if (MU2->find(sr[n]) != MU2->end())
      {
        succ[k] = sr[n];
        cap[k] = pr[n] * param;
        ++k;
      }
    }
    res = CompactMaxFlow<_T>::CreateNetwork(succ, succ + m, cap, cap + m, rmap, m, k - m, known_result, optflags);
    delete [] succ;
    return res;
  }
  else
  {
    cap = new _T[rc];
    for (n = 0; n < rc; ++n) cap[n] = pr[n] * param;
    return CompactMaxFlow<_T>::CreateNetwork(sl, sr, pl, cap, rmap, lc, rc, known_result, optflags);
  }
}

template <typename _T> void CompactMaxFlow_CheckValid<_T>::SetParam(_T p)
{
  int n;
  param = p;
  if (orig && cap) for (n = 0; n < CompactMaxFlow<_T>::n2; ++n)
  {
    cap[CompactMaxFlow<_T>::n1 + n] = orig[n] * param;
  }
}
#endif

#endif//COMPACTMAXFLOW_CC
