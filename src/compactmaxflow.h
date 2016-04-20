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



#ifndef COMPACT_MAXFLOW
#define COMPACT_MAXFLOW


#include "relationmap.h"

template <typename _T> class CompactMaxFlow
{
public:
  CompactMaxFlow();
  virtual ~CompactMaxFlow();
  
  bool CreateNetwork(int*, int *, _T*, _T*, RelationMap*, int, int, bool&, unsigned long);
  bool CreateNetwork(int *successors, _T *probabilities, RelationMap *rmap,
                int l0, int lc, int r0, int rc, bool &known_result,
                unsigned long optflags)
  {
    return CreateNetwork(successors + l0, successors + r0, probabilities + l0,
                probabilities + r0, rmap, lc, rc, known_result, optflags);
  }
  
  bool UpdateNetwork(RelationMap*, bool);
  
  bool IsFlowTotal(bool restart = false);
  
  bool IsCreated() { return valid; };
  void Reset() { _FreeInternals(); incomplete_flow = false; }
  
#ifdef DISABLE_FP_APPROXIMATION
  static inline bool _Teq(const _T &v1, const _T &v2) { return v1 == v2; }
  static inline bool _Tless(const _T &v1, const _T &v2) { return v1 < v2; }
  static inline bool _Tleq(const _T &v1, const _T &v2) { return v1 <= v2; }
#else//DISABLE_FP_APPROXIMATION
  static _T precision;
  static inline bool _Teq(const _T &v1, const _T &v2)
  {
    _T v = 2 * (v1 - v2);
    if (v < 0) v = -v;
    return v <= precision;
  }

  static inline bool _Tless(const _T &v1, const _T &v2) { return v2 - v1 > precision; }
  static inline bool _Tleq(const _T &v1, const _T &v2) { return v1 - v2 < precision; }
#endif//DISABLE_FP_APPROXIMATION

protected:
  void _FreeInternals();
  
  struct arc;
  struct node
  {
    arc **arcs;
    _T excess;
    int label, id;
  };
  struct arc
  {
    node *tail, *head;
    _T flow;
#ifdef OPT_SIGNIFICIANT_ARC
    bool significiant;
#endif
  };
  
  int arc_space;
  
  void _RemoveArc(arc*, bool);
  
  int n1, n2, n_arcs;
  node *set1, *set2;
  arc *arcs;
  arc **arc_lists;
  bool incomplete_flow, valid;
  
  _T *tpl, *tpr;
  
  _T aux;
  
#ifdef DEBUG
  unsigned long space_usage;
  int complexity;
  
public:
  static unsigned long global_space, global_instances, global_space_peak, global_times_invoked, global_p_inv_fails, global_sig_arc_fails;
  static unsigned int min_complexity, max_complexity;
  
  static void ResetStats()
  {
    global_space = 0;
    global_instances = 0;
    global_space_peak = 0;
    global_times_invoked = 0;
    global_p_inv_fails = 0;
    global_sig_arc_fails = 0;
    min_complexity = UINT_MAX;
    max_complexity = 0;
  }
  
  void Dump(const char*);
  int GetComplexity() { return complexity; }
#endif//DEBUG
};

#endif//COMPACT_MAXFLOW
