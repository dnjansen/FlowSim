#ifndef COMPACT_MAXFLOW
#define COMPACT_MAXFLOW

#include <vector>
#include <set>
#include <utility>
#include <map>

#include "relationmap.h"

template <typename _T> class CompactMaxFlow
{
public:
  CompactMaxFlow();
  virtual ~CompactMaxFlow();
  
  virtual bool CreateNetwork(int*, _T*, RelationMap*, int, int, int, int, bool&, unsigned long);
  virtual bool CreateNetwork(int*, int *, _T*, _T*, RelationMap*, int, int, bool&, unsigned long);
  
  bool UpdateNetwork(RelationMap*, bool);
  
  virtual bool IsFlowTotal(bool = false);
  
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
  unsigned char *storage;
  node *set1, *set2;
  arc *arcs;
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
    min_complexity = (unsigned int)-1;
    max_complexity = 0;
  }
  
  void Dump(const char*);
  int GetComplexity() { return complexity; }
#endif//DEBUG
};

#if 0
template <typename _T> class CompactMaxFlow_CheckValid : public CompactMaxFlow<_T>
{
public:
  CompactMaxFlow_CheckValid() { orig = 0, cap = 0, param = 1, MU1 = 0, MU2 = 0; }
  virtual ~CompactMaxFlow_CheckValid() { _FreeInternals(); }
  
  virtual bool CreateNetwork(int*, _T*, RelationMap*, int, int, int, int, bool&, unsigned long);
  virtual bool CreateNetwork(int*, int*, _T*, _T*, RelationMap*, int, int, bool&, unsigned long);
  inline void SetMU(std::set<int> *m1, std::set<int> *m2) { MU1 = m1, MU2 = m2; }
  void SetParam(_T);

  inline bool IsFlowTotal(bool = false) { return CompactMaxFlow<_T>::IsFlowTotal(true); }
  
protected:
  void _FreeInternals()
  {
    CompactMaxFlow<_T>::_FreeInternals();
    if (cap) delete [] cap;
    cap = 0;
  }
  
private:
  _T *orig, *cap, param;
  std::set<int> *MU1, *MU2;
};
#endif

class CompactFeasibleFlow
{
public:
  CompactFeasibleFlow();
  ~CompactFeasibleFlow();
  
  bool CreateNetwork(int*, double*, RelationMap*, int, int, int, int, std::set<int>&, std::set<int>&, unsigned long);

  bool IsFlowFeasible(std::set<double>&, std::set<double>::iterator*);
  inline bool IsFlowFeasible(const double &p) { std::set<double> s; s.insert(p); return IsFlowFeasible(s, 0); }
  
#ifdef DISABLE_FP_APPROXIMATION
  static inline bool _Teq(const double &v1, const double &v2) { return v1 == v2; }
  static inline bool _Tless(const double &v1, const double &v2) { return v1 < v2; }
  static inline bool _Tleq(const double &v1, const double &v2) { return v1 <= v2; }
#else//DISABLE_FP_APPROXIMATION
  static double precision;
  static inline bool _Teq(const double &v1, const double &v2)
  {
    double v = 2 * (v1 - v2);
    if (v < 0) v = -v;
    return v <= precision;
  }

  static inline bool _Tless(const double &v1, const double &v2) { return v2 - v1 > precision; }
  static inline bool _Tleq(const double &v1, const double &v2) { return v1 - v2 < precision; }
#endif//DISABLE_FP_APPROXIMATION
  
private:
  /*struct Node;
  struct Arc
  {
    Node *tail, *head;
    _T flow;
  }
  struct Node
  {
    _T excess, req, aux;
    unsigned int label, id;
    Arc **arcs;
  };
  
  int n1, n2, n_arcs;
  _T laux, raux, rneed, *probabilities;
  
  Node *lset, *rset;
  Arc *arcs;*/
  
  void _FreeInternals();
  
  int sink;
  std::vector<std::pair<int, std::pair<double, int> > > arcs;
  double sourcecap, sinkcap;
};

#endif//COMPACT_MAXFLOW
