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



#ifndef COMPACTMAXFLOW_CC
#define COMPACTMAXFLOW_CC

#include <cstdio>
#include <cstring>
#include "compactmaxflow.h"

// ask explicitly for instantiations of the template:
template class CompactMaxFlow<double>;
template void CompactMaxFlow<double>::_FreeInternals();
template bool CompactMaxFlow<double>::IsFlowTotal(bool restart);
template bool CompactMaxFlow<double>::UpdateNetwork(RelationMap *rmap,
            bool sigarc);
template CompactMaxFlow<double>::CompactMaxFlow();

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
  n1 = n2 = 0;
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
    n1 = n2 = 0;
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
        relsum1 = relsum2 = NULL;
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
      relsum1 = relsum2 = NULL;
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
        relsum1 = relsum2 = NULL;
        _FreeInternals();
        known_result = true;
        return false;
      }
#endif
#ifdef OPT_SIGNIFICIANT_ARC
      //if (optflags & OPT_SIGNIFICIANT_ARC)
      {
        arcs[n].significiant = _Tless(relsum1[l] - set2[r].excess, set1[l].excess)
                                || _Tless(relsum2[r] + aux - set1[l].excess, set2[r].excess);
      }
#endif
    }
    
    delete [] relsum1;
    relsum1 = relsum2 = NULL;
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
  
#if defined(OPT_SIGNIFICIANT_ARC) || defined(OPT_P_INVARIANT)
  assert(NULL == relsum1);
#endif
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
    if (arcs[n].tail && arcs[n].head)
    {
      if (arcs[n].tail->id == arcs[n].head->id) continue;
      if (!(*rmap)(arcs[n].tail->id, arcs[n].head->id)) _RemoveArc(arcs + n, sigarc);
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
  for (n = 0; NULL != set1[l].arcs[n]; ++n)
  {
    if (NULL != set1[l].arcs[n]->head)
      relsum1 += tpr[set1[l].arcs[n]->head - set2];
  }
  for (n = 0; NULL != set1[l].arcs[n]; ++n)
  {
    if (NULL != set1[l].arcs[n]->head && !set1[l].arcs[n]->significiant)
    {
      set1[l].arcs[n]->significiant =
                _Tless(relsum1 - tpr[set1[l].arcs[n]->head - set2], tpl[l]);
    }
  }
  for (n = 0; NULL != set2[r].arcs[n]; ++n)
  {
    if (NULL != set2[r].arcs[n]->tail)
      relsum2 += tpl[set2[r].arcs[n]->tail - set1];
  }
  for (n = 0; NULL != set2[r].arcs[n]; ++n)
  {
    if (NULL != set2[r].arcs[n]->tail && !set2[r].arcs[n]->significiant)
    {
      set2[r].arcs[n]->significiant =
                _Tless(relsum2 - tpl[set2[r].arcs[n]->tail - set1], tpr[r]);
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
    for (m = 0, a = 0, mlast = -1; NULL != set1[n].arcs[m]; ++m)
    {
      if (NULL != set1[n].arcs[m]->head)
      {
        ++a; // Count only active (not deleted) arcs
        if (_Tleq(set1[n].excess, set1[n].arcs[m]->head->excess) &&
                                                    _Tless(0, set1[n].excess))
        {
          h = set1[n].arcs[m]->head;
          // Push entire excess
          set1[n].arcs[m]->head->excess -= set1[n].excess;
          set1[n].arcs[m]->flow += set1[n].excess;
          set1[n].excess = 0;
          
          // Even though we're about to leave the loop, we have to keep counting active (non-deleted) arcs
          while (NULL != set1[n].arcs[++m])
            if (NULL != set1[n].arcs[m]->head) ++a;
          break;
        }
        else if (_Tless(0, set1[n].arcs[m]->head->excess) &&
                                                    _Tless(0, set1[n].excess))
        {
          h = set1[n].arcs[m]->head;
          // Push part of the excess to saturate arc head
          set1[n].arcs[m]->flow += set1[n].arcs[m]->head->excess;
          set1[n].excess -= set1[n].arcs[m]->head->excess;
          set1[n].arcs[m]->head->excess = 0;
          mlast = m; // Remember this arc
        }
        if (mlast == -1) mlast = m;
      }
    }
    
    // Update label.
    if (a == 0 && _Tless(0, set1[n].excess))
    {
      _FreeInternals();
      incomplete_flow = true;
      return false;
    }
    else if (a == 1)
    {
      // This node has only one arc connected; all excess must be pushed through this arc
      // and the label is set to maximum so that no excess can be pushed back to this node
      set1[n].label = n1 + n2;
      if (_Tless(0, set1[n].excess))
      {
        set1[n].arcs[mlast]->head->excess -= set1[n].excess;
        set1[n].arcs[mlast]->flow += set1[n].excess;
        set1[n].excess = 0;
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
      while (_Tless(set2[n].excess, 0))
      {
        keep_going = true;
        for (m = 0, a = -1, min_label = n1 + n2; NULL != set2[n].arcs[m]; ++m)
        {
          if (NULL != set2[n].arcs[m]->tail && _Tless(0,set2[n].arcs[m]->flow))
          {
            l = set2[n].arcs[m]->tail->label;
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
        set2[n].label = min_label + 1;
        if (_Tless(-set2[n].excess, set2[n].arcs[a]->flow))
        {
          set2[n].arcs[a]->flow += set2[n].excess;
          set2[n].arcs[a]->tail->excess -= set2[n].excess;
          set2[n].excess = 0;
        }
        else
        {
          set2[n].excess += set2[n].arcs[a]->flow;
          set2[n].arcs[a]->tail->excess += set2[n].arcs[a]->flow;
          set2[n].arcs[a]->flow = 0;
        }
      }
    }
    
    if (!keep_going) break;
    
    for (n = 0; n < n1; ++n)
    {
      while (_Tless(0, set1[n].excess))
      {
        for (m = 0, a = -1, min_label = n1 + n2; NULL != set1[n].arcs[m]; ++m)
        {
          if (set1[n].arcs[m]->head)
          {
            l = set1[n].arcs[m]->head->label;
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
        set1[n].label = min_label + 1;
        set1[n].arcs[a]->flow += set1[n].excess;
        set1[n].arcs[a]->head->excess -= set1[n].excess;
        set1[n].excess = 0;
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
    printf("%4d[%2d] %#.10f\n", set1[n].id, set1[n].label, set1[n].excess);
    for (m = 0; NULL != set1[n].arcs[m]; ++m)
    {
      h = set1[n].arcs[m]->head;
      if (h)
      {
        printf("                 -- %f (%c) --> %4d[%2d] %#.10f\n", set1[n].arcs[m]->flow,
#ifdef OPT_SIGNIFICIANT_ARC
        (set1[n].arcs[m]->significiant ? '!' : ' '),
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

#endif//COMPACTMAXFLOW_CC
