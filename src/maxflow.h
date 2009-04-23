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



#ifndef MAXFLOW_H
#define MAXFLOW_H

using namespace std;

template <typename _capType> class ParametricMaxflow
{
public:
  ParametricMaxflow();
  ~ParametricMaxflow();

  bool BeginNetwork(int nodes);
  bool EndNetwork();

  bool AddArc(int from, int to, _capType capacity);
  bool AddParametricArc(int from, int to, _capType c_lambda, _capType c_const);
  bool RemoveArc(int from, int to);
  
  bool SetArcCapacity(int from, int to, _capType capacity);
  bool SetParametricArcCapacity(int from, int to, _capType c_lambda, _capType c_const);
  
  _capType Flow();
  
  inline bool IsSimple();
  inline bool IsParametric();
  inline void SetOptimize(bool);
  
  inline int Source() { return -2; }
  inline int Sink() { return -3; }
  
  unsigned long GetSpaceUsage() { return space_usage; }
  unsigned long GetNodes() { return _n; }
  
private:
  // Type declarations for the internal algorithm
  struct node;
  struct capacity
  {
    _capType c_const, c_lambda;/* parametric capacity */
  };
  
  struct arc
  {
    capacity  cap;             /* full capacity */
    _capType  resCap;          /* residual capacity */
    node      *head;           /* arc head */
    arc       *rev;            /* reverse arc */
    bool      parametric;      /* parametric arc */
    bool      reverse;         /* symbolic arc that denotes the reverse of another arc */
  };
  
  struct arclist
  {
    capacity cap;              /* arc capacity */
    int from, to;              /* start and end nodes */
    arclist *next;             /* next item in linked list */
  };
  
  struct node
  {
    arc       *first;          /* first outgoing arc */
    arc       *current;        /* current outgoing arc */
    _capType  excess;          /* excess at the node change to double if needed */
    long      d;               /* distance label */
    node      *bNext;          /* next node in bucket */
    node      *bPrev;          /* previous node in bucket */
  };
  
  struct bucket
  {
    node      *firstActive;    /* first node with positive excess */
    node      *firstInactive;  /* first node with zero excess */
  };
  
  // Class members
  long _n;                 /* number of nodes */
  long _m;                 /* number of arcs */
  long _nm;                /* n + ALPHA * m */
  long _nMin;              /* smallest node id */
  node *_nodes;            /* array of nodes */
  arc *_arcs;              /* array of arcs */
  bucket *_buckets;        /* array of buckets */
  node *_source;           /* source node pointer */
  node *_sink;             /* sink node pointer */
  long _dMax;              /* maximum label */
  long _aMax;              /* maximum active node label */
  long _aMin;              /* minimum active node label */
  _capType _flow;          /* flow value */
  node *sentinelNode;      /* end of the node list marker */
  arc *stopA;              /* used in forAllArcs */
  long _workSinceUpdate;   /* the number of arc scans since last update */
  float _globUpdtFreq;     /* global update frequency */
  long i_dist;
  node *i_next, *i_prev;
  arclist *newArcList;    // Linked list of arcs being added as the network is created
  _capType param;         // Parameter setting for parametric networks
  unsigned char flags;    // Bitset with flags to conserve space, see maxflow.cc for more info
  unsigned char *storage;
  int *node_table;
  unsigned long space_usage;
  
  // Private functions
  void Resynch(); // Executes maxflow algorithm if synched==false
  void Dealloc(); // Free data structures
  int GetNodeID(int, bool); // Get internal node number from user-space node number
  
  // Legacy routines from the IGSYS hi_pr implementation
  int allocDS();
  void Init();
  void ReInit();
  void GlobalUpdate();
  void stageTwo();
  long Relabel(node*);
  int Gap(bucket*);
  void Preflow();
  void Discharge(node*);
  
  inline void aAdd(bucket*, node*);
  inline void aRemove(bucket*, node*);
  inline void iAdd(bucket*, node*);
  inline void iDelete(bucket*, node*);
  inline void swapArcs(arc*, arc*);
};

#endif//MAXFLOW_H
