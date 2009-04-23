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



/* Maximum flow - highest lavel push-relabel algorithm */
/* COPYRIGHT C 1995, 2000 by IG Systems, Inc., igsys@eclipse.net */
/* MODIFIED -- code below is not the original IG Systems code */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define GLOB_UPDT_FREQ 0.5
#define ALPHA 6
#define BETA 12
#define MAXLONG (long)0x40000000

#define WHITE 0
#define GREY 1
#define BLACK 2

// Defines for ParametricMaxflow<_T>::flags
#define GMFF_SYNCHED      0x0001 // Unset if network was modified
#define GMFF_VIRGIN       0x0002 // Set if maxflow must be computed from scratch
#define GMFF_CREATED      0x0004 // Network has been created
#define GMFF_CREATING     0x0008 // Network is being created
#define GMFF_SIMPLE       0x0010 // longest source/sink path is 3 and only source/sink arcs are parametric
#define GMFF_PARAMETRIC   0x0020 // only source/sink arcs are parametric, false=no parametric arcs or not a valid param. network
#define GMFF_OPTIMIZE     0x0040 // Use parametric maxflow optimizations

#include "maxflow.h"

/* macros */

#define forAllNodes(i) for (i = _source; i <= _sink; ++i)
#define forAllArcs(i,a) for (a = i->first, stopA = (i+1)->first; a < stopA; ++a)

#define nNode(i)   ((i) - _nodes)
#define pNode(n)   (_nodes + (n))
#define nArc(a)    ((a == NULL) ? -1 : (a) - _arcs)

#define min(a, b)  (((a) < (b)) ? a : b)

/*
   bucket functions:
   bucket's active node list is singly-linked
     operations aAdd, aRemove (from the front)
   bucket's inactive list is doubly-linked
     operations iAdd, iDelete (from arbitrary position)
*/

template <typename _capType> void ParametricMaxflow<_capType>::aAdd(bucket *l, node *i)
{
  assert(l >= _buckets && l < _buckets + _n + 2);
  assert(i >= _nodes && i < _nodes + _n + 2);
  i->bNext = l->firstActive;
  l->firstActive = i;
  i_dist = i->d;
  if (i_dist < _aMin) _aMin = i_dist;
  if (i_dist > _aMax) _aMax = i_dist;
  if (_dMax < _aMax) _dMax = _aMax;
}

/* i must be the first element */
template <typename _capType> void ParametricMaxflow<_capType>::aRemove(bucket *l, node *i)
{
  assert(l >= _buckets && l < _buckets + _n + 2);
  assert(i >= _nodes && i < _nodes + _n + 2);
  assert(l->firstActive == i);
  l->firstActive = i->bNext;
}

template <typename _capType> void ParametricMaxflow<_capType>::iAdd(bucket *l, node *i)
{
  assert(l >= _buckets && l < _buckets + _n + 2);
  assert(i >= _nodes && i < _nodes + _n + 2);
  i_next = l->firstInactive;
  i->bNext = i_next;
  i->bPrev = sentinelNode;
  i_next->bPrev = i;
  l->firstInactive = i;
}

template <typename _capType> void ParametricMaxflow<_capType>::iDelete(bucket *l, node *i)
{
  assert(l >= _buckets && l < _buckets + _n + 2);
  assert(i >= _nodes && i < _nodes + _n + 2);
  i_next = i->bNext;
  if (l->firstInactive == i)
  {
    l->firstInactive = i_next;
    i_next->bPrev = sentinelNode;
  }
  else
  {
    i_prev = i->bPrev;
    i_prev->bNext = i_next;
    i_next->bPrev = i_prev;
  }
}

template <typename _capType> ParametricMaxflow<_capType>::ParametricMaxflow()
{
  flags = 0;
  _nodes = 0, _arcs = 0, _buckets = 0;
  _dMax = 0;
  storage = 0;
  node_table = 0;
  space_usage = 0;
}

template <typename _capType> ParametricMaxflow<_capType>::~ParametricMaxflow()
{
  Dealloc();
}

template <typename _capType> void ParametricMaxflow<_capType>::Dealloc()
{
  delete [] storage;
  delete [] node_table;
  storage = 0, _nodes = 0, _arcs = 0, _buckets = 0, node_table = 0;
}

// Initialize the maximum flow problem. If the update parameter is true,
// only labels and buckets will be initialized and flow/excess values
// will be left as they are.
template <typename _capType> void ParametricMaxflow<_capType>::Init()
{
  node  *i;        /* current node */
  bucket *l;
  arc *a;
  _capType delta;

  // Set residual capacities
  for (a = _arcs; a != _arcs + (2 * _m) + 1; ++a) a->resCap = a->cap.c_const;
  
  forAllNodes(i)
  {
    i->excess = 0;
    i->current = i->first;
  }
  
  // Saturate all source arcs
  _source->excess = 0;
  forAllArcs(_source, a)
  {
    if (a->head != _source)
    {
      delta = a->resCap;
      a->resCap = (_capType)0;
      (a->rev)->resCap = delta;
      a->head->excess += delta;
    }
  }

  // clear all buckets
  for (l = _buckets; l <= _buckets + _n+1; l++)
  {
    l->firstActive = sentinelNode;
    l->firstInactive = sentinelNode;
  }

  /*  setup labels and buckets */
  l = _buckets + 1;

  _aMax = 0;
  _aMin = _n;

  forAllNodes(i)
  {
    if (i == _sink)
    {
      i->d = 0;
      iAdd(_buckets, i);
      continue;
    }
    // Source is labelled _n for a complete init. To prepare for an update, source
    //if (i == _source) i->d = (update ? 1 : _n); // is labelled 1 like all other nodes.
    if (i == _source) i->d = _n; // is labelled 1 like all other nodes.
    else i->d = 1;
    if (i->excess > 0) aAdd(l, i); /* put into active list */
    else if (i->d < _n) iAdd(l, i); /* put into inactive list */
  }

  _dMax = 1;
  _flow = 0;
}

template <typename _capType> void ParametricMaxflow<_capType>::ReInit()
{
  node  *i;        /* current node */
  bucket *l;

  // clear all buckets
  for (l = _buckets; l <= _buckets + _n+1; l++)
  {
    l->firstActive = sentinelNode;
    l->firstInactive = sentinelNode;
  }

  /*  setup labels and buckets */
  l = _buckets + 1;

  _aMax = 0;
  _aMin = _n;

  forAllNodes(i)
  {
    if (i == _sink)
    {
      i->d = 0;
      iAdd(_buckets, i);
      continue;
    }
    //i->d = 1;
    //if (i->excess > 0) aAdd(l, i); /* put into active list */
    //else if (i->d < _n) iAdd(l, i); /* put into inactive list */
    if (i->excess > 0) aAdd(_buckets + i->d, i); /* put into active list */
    else if (i->d <= _n) iAdd(_buckets + i->d, i); /* put into inactive list */
  }

  _dMax = 1;
  _flow = 0;
}

/* global update via backward breadth first search from the sink */
template <typename _capType> void ParametricMaxflow<_capType>::GlobalUpdate()
{
  node *i, *j;       /* node pointers */
  arc *a;           /* current arc pointers  */
  bucket *l, *jL;          /* bucket */
  long curDist, jD;
  long state;

  /* initialization */
  forAllNodes(i) i->d = _n;
  _sink->d = 0;

  for (l = _buckets; l <= _buckets + _dMax; ++l)
  {
    l->firstActive   = sentinelNode;
    l->firstInactive  = sentinelNode;
  }

  _dMax = _aMax = 0;
  _aMin = _n;

  /* breadth first search */

  // add sink to bucket zero
  iAdd(_buckets, _sink);
  for (curDist = 0; 1; curDist++)
  {
    state = 0;
    l = _buckets + curDist;
    jD = curDist + 1;
    jL = l + 1;

    if ((l->firstActive == sentinelNode) &&
        (l->firstInactive == sentinelNode)) break;

    while (1)
    {
      switch (state)
      {
      case 0:
        i = l->firstInactive;
        state = 1;
        break;
      case 1:
        i = i->bNext;
        break;
      case 2:
        i = l->firstActive;
        state = 3;
        break;
      case 3:
        i = i->bNext;
        break;
      default:
        assert(0);
        break;
      }

      if (i == sentinelNode)
      {
        if (state == 1)
        {
          state = 2;
          continue;
        }
        else
        {
          assert(state == 3);
          break;
        }
      }

      /* scanning arcs incident to node i */
      forAllArcs(i, a)
      {
        if (a->rev->resCap > 0)
        {
          j = a->head;
          if (j->d == _n)
          {
            j->d = jD;
            j->current = j->first;
            if (jD > _dMax) _dMax = jD;

            if (j->excess > 0) aAdd(jL,j); /* put into active list */
            else iAdd(jL,j); /* put into inactive list */
          }
        }
      } /* node i is scanned */
    }
  }
}

/* second stage -- preflow to flow */
/*
   do dsf in the reverse flow graph from nodes with excess
   cancel cycles if found
   return excess flow in topological order

   i->d is used for dfs labels
   i->bNext is used for topological order list
   buckets[i-nodes]->firstActive is used for DSF tree
*/
template <typename _capType> void ParametricMaxflow<_capType>::stageTwo()
{
  node *i, *j, *tos, *bos, *restart, *r;
  arc *a;
  _capType delta;

  /* deal with self-loops */
  forAllNodes(i)
  {
    forAllArcs(i, a)
    {
      if (a->head == i) a->resCap = a->cap.c_const; //FIXME parametric here?
    }
  }

  /* initialize */
  tos = bos = NULL;
  forAllNodes(i)
  {
    i->d = WHITE;
    _buckets[i-_nodes].firstActive = sentinelNode;
    i->current = i->first;
  }

  /* eliminate flow cycles, topologicaly order vertices */
  forAllNodes(i)
  {
    if ((i->d == WHITE) && (i->excess > 0) && (i != _source) && (i != _sink))
    {
      r = i;
      r->d = GREY;

      while (true)
      {
        for (; i->current != (i + 1)->first; ++i->current)
        {
          a = i->current;
          if ((a->cap.c_const == (_capType)0) && (a->resCap > (_capType)0))
          {
            j = a->head;
            if (j->d == WHITE)
            {
              /* start scanning j */
              j->d = GREY;
              _buckets[j-_nodes].firstActive = i;
              i = j;
              break;
            }
            else if (j->d == GREY)
            {
              /* find minimum flow on the cycle */
              delta = a->resCap;
              while (1)
              {
                delta = min(delta, j->current->resCap);
                if (j == i) break;
                else j = j->current->head;
              }

              /* remove delta flow units */
              j = i;
              while (1)
              {
                a = j->current;
                a->resCap -= delta;
                a->rev->resCap += delta;
                j = a->head;
                if (j == i) break;
              }

              /* backup DFS to the first saturated arc */
              restart = i;
              for (j = i->current->head; j != i; j = a->head)
              {
                a = j->current;
                if (( j->d == WHITE ) || ( a->resCap == 0 ))
                {
                  j->current->head->d = WHITE;
                  if (j->d != WHITE) restart = j;
                }
              }

              if (restart != i)
              {
                i = restart;
                i->current++;
                break;
              }
            }
          }
        }

        if (i->current == (i+1)->first)
        {
          /* scan of i complete */
          i->d = BLACK;
          if (i != _source)
          {
            if (bos == NULL)
            {
              bos = i;
              tos = i;
            }
            else
            {
              i->bNext = tos;
              tos = i;
            }
          }

          if (i != r)
          {
            i = _buckets[i-_nodes].firstActive;
            i->current++;
          }
          else break;
        }
      }
    }
  }

  /* return excesses */
  /* note that sink is not on the stack */
  if (bos != NULL)
  {
    for (i = tos; i != bos; i = i->bNext)
    {
      a = i->first;
      while (i->excess > 0)
      {
        if ((a->cap.c_const == (_capType)0) && (a->resCap > (_capType)0))
        {
          if (a->resCap < i->excess) delta = a->resCap;
          else delta = i->excess;
          a->resCap -= delta;
          a->rev->resCap += delta;
          i->excess -= delta;
          a->head->excess += delta;
        }
        a++;
      }
    }
    /* now do the bottom */
    i = bos;
    a = i->first;
    while (i->excess > 0)
    {
      if ((a->cap.c_const == (_capType)0) && (a->resCap > (_capType)0))
      {
        if (a->resCap < i->excess) delta = a->resCap;
        else delta = i->excess;
        a->resCap -= delta;
        a->rev->resCap += delta;
        i->excess -= delta;
        a->head->excess += delta;
      }
      a++;
    }
  }
}

/* gap relabeling */
template <typename _capType> int ParametricMaxflow<_capType>::Gap(bucket *emptyB)
{
  bucket *l;
  node  *i;
  long  r;           /* index of the bucket before l  */

  r = (emptyB - _buckets) - 1;

  /* set labels of nodes beyond the gap to "infinity" */
  for (l = emptyB + 1; l <= _buckets + _dMax; ++l)
  {
    for (i = l->firstInactive; i != sentinelNode; i = i->bNext) i->d = _n;
    l->firstInactive = sentinelNode;
  }

  _dMax = r;
  _aMax = r;

  return (_aMin > r) ? 1 : 0;
}

/*--- relabelling node i */

template <typename _capType> long ParametricMaxflow<_capType>::Relabel(node *i)
{
  node  *j;
  long  minD;     /* minimum d of a node reachable from i */
  arc   *minA;    /* an arc which leads to the node with minimal d */
  arc   *a;

  assert(i->excess > 0);

  _workSinceUpdate += BETA;

  i->d = minD = _n;
  minA = NULL;

  /* find the minimum */
  forAllArcs(i, a)
  {
    _workSinceUpdate++;
    if (a->resCap > 0)
    {
      j = a->head;
      if (j->d < minD)
      {
        minD = j->d;
        minA = a;
      }
    }
  }

  minD++;

  if (minD < _n)
  {
    i->d = minD;
    i->current = minA;
    if (_dMax < minD) _dMax = minD;
  } /* end of minD < n */

  return minD;
}

/* discharge: push flow out of i until i becomes inactive */
template <typename _capType> void ParametricMaxflow<_capType>::Discharge(node *i)
{
  node *j;                 /* sucsessor of i */
  long jD;                 /* d of the next bucket */
  bucket *lj;               /* j's bucket */
  bucket *l;                /* i's bucket */
  arc *a;                 /* current arc (i,j) */
  _capType delta;
  arc *stopA;

  assert(i->excess > 0);
  assert(i != _sink);
  
  while (true)
  {
    jD = i->d - 1;
    l = _buckets + i->d;
    
    assert(i->first);

    /* scanning arcs outgoing from  i  */
    for (a = i->current, stopA = (i+1)->first; a < stopA; a++)
    {
      if (a->resCap > 0)
      {
        j = a->head;

        if (j->d == jD)
        {
          if (a->resCap < i->excess) delta = a->resCap;
          else delta = i->excess;
          a->resCap -= delta;
          a->rev->resCap += delta;

          if (j != _sink)
          {
            lj = _buckets + jD;

            if (j->excess == 0)
            {
              iDelete(lj,j); /* remove j from inactive list */
              aAdd(lj,j); /* add j to active list */
            }
          }

          j->excess += delta;
          i->excess -= delta;

          if (i->excess == 0) break;

        } /* j belongs to the next bucket */
      } /* a  is not saturated */
    } /* end of scanning arcs from  i */

    if (a == stopA)
    {
      /* i must be relabeled */
      Relabel(i);

      if (i->d == _n) break;
      if ((l->firstActive == sentinelNode) &&
          (l->firstInactive == sentinelNode)) Gap(l);

      if (i->d == _n) break;
    }
    else
    {
      /* i no longer active */
      i->current = a;
      /* put i on inactive list */
      iAdd(l,i);
      break;
    }
  }
}

/* first stage  -- maximum preflow*/
template <typename _capType> void ParametricMaxflow<_capType>::Preflow()
{
  node   *i;
  bucket  *l;             /* current bucket */

  _workSinceUpdate = 0;

  /* main loop */
  while (_aMax >= _aMin)
  {
    l = _buckets + _aMax;
    i = l->firstActive;

    if (i == sentinelNode) _aMax--;
    else
    {
      aRemove(l,i);

      assert(i->excess > 0);
      Discharge(i);

      if (_aMax < _aMin) break;

      /* is it time for global update? */
      if (_workSinceUpdate * GLOB_UPDT_FREQ > _nm)
      {
        GlobalUpdate();
        _workSinceUpdate = 0;
      }
    }
  } /* end of the main loop */

  _flow = _sink->excess;
}

// Begin the construction of a new network. The previous network will be lost. This function
// must be called before AddArc() or AddParametricArc(). Call EndNetwork() to indicate that
// all arcs have been specified. The function fails and returns false only if a network is
// already being created.
template <typename _capType> bool ParametricMaxflow<_capType>::BeginNetwork(int num_nodes)
{
  int i;
  
  // Fail if BeginNetwork() has been called before
  if (flags & GMFF_CREATING) return false;
  
  if (num_nodes < 3) return false;

  Dealloc();

  flags |= GMFF_CREATING;

  _n = num_nodes;

  newArcList = 0;
  
  // Initialize translation table for node numbers
  node_table = new int[num_nodes];
  *node_table = Source();
  *(node_table + num_nodes - 1) = Sink();
  for (i = 1; i < num_nodes - 1; ++i) *(node_table + i) = -1;
  
  return true;
}

// Add an arc to the network; must be called between calls to BeginNetwork() and EndNetwork().
// The function fails and returns false if BeginNetwork() was not called before, if the node
// numbers are invalid (must be between 1 and the number of nodes in the network) or if the capacity
// is negative.
template <typename _capType> bool ParametricMaxflow<_capType>::AddArc(int from, int to, _capType capacity)
{
  return AddParametricArc(from, to, capacity, (_capType)0);
}

// Add an arc with parametric capacity. This function fails under the same conditions as
// AddArc(). Both c_lambda and c_const must be positive. The capacity cap of the arc with a parameter p
// is computed by: cap = (c_lambda * p) + c_const
template <typename _capType> bool ParametricMaxflow<_capType>::AddParametricArc(int from, int to, _capType c_const, _capType c_lambda)
{
  arclist *na;

  // Fail if BeginNetwork() hasn't been called
  if ((flags & GMFF_CREATING) == 0) return false;

  // Verify parameters are not out-of-bounds
  if ((from < 0 && from != Source()) || (to < 0 && to != Sink()) || c_lambda < (_capType)0 || c_const < (_capType)0) return false;
  
  // Verify that this arc isn't in the network already
  na = newArcList;
  while (na)
  {
    if (na->from == from && na->to == to) return false;
    na = na->next;
  }

  na = new arclist;
  assert(na);
  assert(na > (arclist*)0x100);
  na->cap.c_lambda = c_lambda;
  na->cap.c_const = c_const;
  na->next = newArcList;
  na->from = from;
  na->to = to;
  
  GetNodeID(from, true);
  GetNodeID(to, true);

  newArcList = na;

  return true;
}

// Finish adding arcs and construct the network
template <typename _capType> bool ParametricMaxflow<_capType>::EndNetwork()
{
  long head, tail, i, last, arc_num, arc_new_num, pos_current = 0, n = _n, m, source, sink;
  long *arc_first, *arc_tail;

  node *ndp, *nodes;
  arc *arc_current = 0, *arc_new, *arcs;

  arclist *al = newArcList, *altmp;

  // Determine number of arcs to add
  m = 0;
  while (al)
  {
    ++m;
    al = al->next;
  }
  
  // Determine if the requested number of nodes is actually fully used
  for (n = 0; n < _n - 1 && *(node_table + n) != -1; ++n);
  if (n < _n - 1)
  {
    // We have less nodes than requested; adjust sink node and total node count
    *(node_table + n) = *(node_table + _n - 1);
    ++n;
    _n = n;
  }

  source = 0;
  sink = _n - 1;

  storage = new unsigned char[(sizeof(node) * (_n + 2)) + (sizeof(arc) * (2 * m + 1)) + (sizeof(bucket) * (_n + 2))];
  nodes = (node*)storage;
  arcs = (arc*)(storage + (sizeof(node) * (_n + 2)));
  arc_tail = new long[(2 * m) + (_n + 2)];
  arc_first = arc_tail + (2 * m);
  
  arc_current = arcs;
  
  memset(arc_tail, 0, sizeof(long) * ((2 * m) + (_n + 2)));

  // Create arcs and delete the intermediary data structures
  al = newArcList;
  while (al)
  {
    tail = GetNodeID(al->from, false);
    head = GetNodeID(al->to, false);
    
    assert(tail != -1);
    assert(head != -1);

    ++arc_first[tail + 1];
    ++arc_first[head + 1];

    /* storing information about the arc */
    arc_tail[pos_current]           = tail;
    arc_tail[pos_current+1]         = head;
    (arc_current    )->head         = nodes + head;
    (arc_current    )->resCap       = al->cap.c_const;
    (arc_current    )->cap          = al->cap;
    (arc_current    )->rev          = arc_current + 1;
    (arc_current    )->parametric   = (al->cap.c_lambda > (_capType)0);
    (arc_current    )->reverse      = false;
    (arc_current + 1)->head         = nodes + tail;
    (arc_current + 1)->resCap       = 0;
    (arc_current + 1)->cap.c_const  = (_capType)0;
    (arc_current + 1)->cap.c_lambda = (_capType)0;
    (arc_current + 1)->rev          = arc_current;
    (arc_current + 1)->parametric   = false;
    (arc_current + 1)->reverse      = true;

    arc_current += 2;
    pos_current += 2;

    altmp = al;
    al = al->next;
    delete altmp;
  }

  /********** ordering arcs - linear time algorithm ***********/

  /* first arc from the first node */
  nodes->first = arcs;

  /* before below loop arc_first[i+1] is the number of arcs outgoing from i;
    after this loop arc_first[i] is the position of the first
    outgoing from node i arcs after they would be ordered;
    this value is transformed to pointer and written to node.first[i]
    */
  
  for (i = 1; i <= sink + 1; ++i)
  {
    arc_first[i] += arc_first[i - 1];
    (nodes + i)->first = arcs + arc_first[i];
  }

  for (i = 0; i <= sink; ++i) // scanning all the nodes
  {                                      //      exept the last*/
    last = ((nodes + i + 1)->first) - arcs;
                            /* arcs outgoing from i must be cited
                              from position arc_first[i] to the position
                              equal to initial value of arc_first[i+1]-1  */

    for (arc_num = arc_first[i]; arc_num < last; ++arc_num)
    {
      tail = arc_tail[arc_num];

      while (tail != i)
      {    /* the arc no  arc_num  is not in place because arc cited here
              must go out from i;
              we'll put it to its place and continue this process
              until an arc in this position would go out from i */

        arc_new_num = arc_first[tail];
        arc_current = arcs + arc_num;
        arc_new     = arcs + arc_new_num;

        /* arc_current must be cited in the position arc_new
          swapping these arcs:                                 */
        swapArcs(arc_new, arc_current);

        arc_tail[arc_num] = arc_tail[arc_new_num];
        arc_tail[arc_new_num] = tail;

        /* we increase arc_first[tail]  */
        ++arc_first[tail];

        tail = arc_tail[arc_num];
      }
    }
    /* all arcs outgoing from  i  are in place */
  }

  /* -----------------------  arcs are ordered  ------------------------- */

  /*----------- constructing lists ---------------*/

  for (ndp = nodes; ndp < nodes + n + 1; ++ndp)
  {
    ndp->first = 0;
    ndp->current = 0;
    ndp->d = 0;
  }

  for (arc_current = arcs + (2 * m - 1); arc_current >= arcs; --arc_current)
  {
    arc_num = arc_current - arcs;
    tail = arc_tail[arc_num];
    ndp = nodes + tail;
    ndp->first = arc_current;
  }

  _m = m;
  _source = nodes + source;
  _sink   = nodes + sink;
  _nodes = nodes;
  _arcs = arcs;
  _nm = ALPHA * _n + _m;
  _buckets = (bucket*)(storage + (sizeof(node) * (_n + 2)) + (sizeof(arc) * (2 * _m + 1)));
  
  sentinelNode = nodes + _n;
  sentinelNode->first = _arcs + (2*_m);

  /* free intermediary memory */
  delete[] arc_tail;
  
  // Set distance label to 1 on all nodes which equal source/sink or are directly connected to them
  _source->d = 1;
  _sink  ->d = 1;
  forAllArcs(_source, arc_current) arc_current->head->d = 1;
  forAllArcs(_sink  , arc_current) arc_current->head->d = 1;
  
  // If all nodes now have a distance label of one, this network is considered simple (Note: This may be
  // changed if the network is found to be parametric below)
  flags |= GMFF_SIMPLE; // Set simple
  for (ndp = _source; ndp <= _sink; ++ndp)
  {
    if (ndp->d != 1)
    {
      flags &= ~GMFF_SIMPLE; // Unset simple
      break;
    }
  }

  // Look for parametric arcs to determine if this network is parametric. If at least one parametric arc
  // is detected that is not directly connected to source or sink, the network is also not simple.
  flags &= ~GMFF_PARAMETRIC; // Unset parametric
  for (arc_current = arcs, arc_num = 0; arc_num < 2 * m; ++arc_current, ++arc_num)
  {
    if (arc_current->parametric)
    {
      flags |= GMFF_PARAMETRIC; // Set parametric
      if (arc_current->head != _source && arc_current->rev->head != _source
       && arc_current->head != _sink   && arc_current->rev->head != _sink)
      {
        flags &= ~GMFF_SIMPLE; // Unset simple
        break;
      }
    }
  }

  flags &= ~GMFF_CREATING; // Unset creating

  // Fail if either source of sink have no arcs connected at all
  if (!_source->first || !_sink->first)
  {
    Dealloc();
    return false;
  }
  
  space_usage = sizeof(ParametricMaxflow<_capType>) + (sizeof(node) * (_n + 2)) + (sizeof(arc) * (2 * m + 1)) + (sizeof(bucket) * (_n + 2)) + (sizeof(int) * _n);
  
  flags = (flags & (~GMFF_SYNCHED)) | GMFF_CREATED | GMFF_VIRGIN; // Unset synched, set created & virgin

  return true;
}

// Swap the contents of two arc structures
template <typename _capType> void ParametricMaxflow<_capType>::swapArcs(arc *a1, arc *a2)
{
  arc tmp, *r1, *r2;
  r1 = a1->rev;
  r2 = a2->rev;
  memcpy(&tmp, a1, sizeof(arc));
  memcpy(a1, a2, sizeof(arc));
  memcpy(a2, &tmp, sizeof(arc));
  if (r1 == a2) // If the two arcs are reverse to each other,
  {             // revert the .rev members
    a1->rev = r1;
    a2->rev = r2;
  }
  a1->rev->rev = a1; // Make sure .rev pointers are set correctly after the swap
  a2->rev->rev = a2;
}

// Private: Updates the flow if the network has been modified.
template <typename _capType> void ParametricMaxflow<_capType>::Resynch()
{
  // If no network has been created or the network is in synch, there's nothing to do (woohoo)
  if (!(flags & GMFF_CREATED) || (flags & GMFF_SYNCHED)) return;
  
  if (flags & GMFF_OPTIMIZE)
  {
    // If the network is marked for computation from scratch or if can't be optimized, prepare
    // for recomputing entirely. Otherwise prepare for update.
    if ((flags & GMFF_VIRGIN) || !(flags & GMFF_SIMPLE)) Init();
    else ReInit();
  }
  else
  {
    Init(); // Prepare for recomputing entirely
  }
  
  Preflow();

  // Unset virgin, set synched
  flags = (flags & (~GMFF_VIRGIN)) | GMFF_SYNCHED;
}

// Returns the maximum flow of the network. Repeatedly calling this function will NOT cause the flow to
// be recomputed if the network was not modified in between calls.
template <typename _capType> _capType ParametricMaxflow<_capType>::Flow()
{
  Resynch();
  
  //node *p;
  //for (p = _source; p <= _sink; ++p) printf("d(%d) = %d   e(%d) = %d\n", nNode(p), p->d, nNode(p), p->excess);
  return _flow;
}

// Remove an arc. This implementation merely sets its capacity to zero rather than actually removing it, but
// this is not a required/guaranteed behavior. The function fails if the arc is not found in the network.
template <typename _capType> bool ParametricMaxflow<_capType>::RemoveArc(int from, int to)
{
  return SetParametricArcCapacity(from, to, (_capType)0, (_capType)0);
}

template <typename _capType> bool ParametricMaxflow<_capType>::SetArcCapacity(int from, int to, _capType capacity)
{
  return SetParametricArcCapacity(from, to, capacity, (_capType)0);
}

// Set the capacity of a parametric arc (or turn a constant-capacity arc into a parametric arc). The capacity
// of an arc is defined as:
//    c(from,to) = c_const + (c_lambda * x)
// where "x" is the network parameter. Passing 0 for c_lambda causes the arc to have constant capacity, otherwise
// it will be parametric. The function fails if the arc can't be found in the network. The nodes may be specified
// in inverse order, i.e. for an arc (a,b), passing (b,a) to this function will correctly modify (a,b). Depending
// on the nature of the change and whether or not optimizations are enabled, the computation may be delayed until
// Flow() is called.
template <typename _capType> bool ParametricMaxflow<_capType>::SetParametricArcCapacity(int from, int to, _capType c_const, _capType c_lambda)
{
  arc *a, *b;
  _capType delta;
  node *f, *t;
  
  from = GetNodeID(from, false);
  to = GetNodeID(to, false);
  
  if (from == -1 || to == -1) return false;
  
  f = pNode(from); // Get pointers to the referenced nodes
  t = pNode(to);
  
  if (f < _source || f > _sink || t < _source || t > _sink) return false;
  
  // Search for the referenced arc
  forAllArcs(f, a)
  {
    if (a->head == t)
    {
      if (a->reverse) a = a->rev; // Reverse arcs have zero base capacity, work on the real arc instead
      
      if (a->cap.c_const == c_const && a->cap.c_lambda == c_lambda) return true; // No change
      
      if (!(flags & GMFF_OPTIMIZE))
      {
        a->cap.c_const = c_const;
        a->cap.c_lambda = c_lambda;
        flags = (flags & (~GMFF_SYNCHED)) | GMFF_VIRGIN; // Unset synched, set virgin
        return true;
      }
      
      delta = a->rev->resCap; // Get current flow through arc
      a->cap.c_const = c_const; // Set new capacity
      a->cap.c_lambda = c_lambda;
      
      // Changing type of arc (parametric or non-parametric) requires rebuild
      if ((c_lambda == (_capType)0 &&  a->parametric)
       || (c_lambda != (_capType)0 && !a->parametric))
      {
        flags = (flags & (~GMFF_SYNCHED)) | GMFF_VIRGIN; // Unset synched, set virgin
        return true;
      }
      
      //printf("changing arc %d --> %d with flow %d from cap %d to cap %d\n", from, to, delta, a->cap.c_const, c_const);
      //a->cap.c_const = c_const; // Set new capacity
      //a->cap.c_lambda = c_lambda;
      
      // New capacity equals flow through arc, just set residual capacity, no further update required
      if (c_const + (c_lambda * (a->parametric ? param : (_capType)0)) == delta)
      {
        a->resCap = (_capType)0;
        return true;
      }
      
      // Is the wew capacity greater than the current flow?
      if (c_const + (c_lambda * (a->parametric ? param : (_capType)0)) > delta)
      {
        if (a->resCap == (_capType)0)
        {
          // Since residual capacity is zero on this arc, the flow might increase with a greater
          // capacity. We have to add more excess.
          a->resCap = c_const + (c_lambda * (a->parametric ? param : (_capType)0)) - a->rev->resCap;
          _source->excess += c_const + (c_lambda * (a->parametric ? param : (_capType)0)) - delta;
          flags &= ~GMFF_SYNCHED; // Unset synched
          return true;
        }
        
        // If residual capacity is greater than zero since there can't be more flow, so we just
        // need to set the new residual capacity.
        a->resCap = c_const + (c_lambda * (a->parametric ? param : (_capType)0)) - a->rev->resCap;
        return true;
      }
      
      // New capacity is lower than flow. This operation can't be optimized on all networks.
      if (!(flags & GMFF_SIMPLE) || a->rev->head == _source)
      {
        // The network isn't of bipartite structure or the arc is connected to the source.
        flags = (flags & (~GMFF_SYNCHED)) | GMFF_VIRGIN; // Unset synched, set virgin
        return true;
      }
      
      // Pull excess back to the arc's tail
      delta -= c_const; // This is the amount of excess that will be moved back
      a->rev->head->excess += delta;
      a->resCap = (_capType)0;
      a->rev->resCap = c_const + (a->parametric ? c_lambda * param : (_capType)0);
      
      // If the arc is connected to the sink, simply subtract the delta from the sink's excess and mark for update
      if (a->head == _sink)
      {
        a->head->excess -= delta;
        flags &= ~GMFF_SYNCHED; // Unset synched
        return true;
      }
      
      // If there is any excess at the arc head, pull that back as far as possible/necessary
      if (a->head->excess > delta)
      {
        a->head->excess -= delta;
        flags &= ~GMFF_SYNCHED;
        return true;
      }
      else if (a->head->excess > 0)
      {
        delta -= a->head->excess;
        a->head->excess = 0;
      }
      
      // Iterate through the arcs to find the sink arc (there can only be one)
      forAllArcs(a->head, b)
      {
        //printf("[a %d->%d cap=%d rescap=%d flow=%d rev=%c]\n", nNode(b->rev->head), nNode(b->head), b->cap.c_const, b->resCap, b->rev->resCap, (b->reverse ? '1' : '0'));
        if (b->head != _sink) continue;
        
        b->resCap += delta;
        b->rev->resCap -= delta;
        b->head->excess -= delta; // b->head is sink, so the result won't be negative
        flags &= ~GMFF_SYNCHED; // Unset synched
        return true;
      }
      
      // The constraints of a valid preflow don't allow this point to be reached. If
      // it is reached, something has gone very wrong.
      assert(false);
      return false;
    }
  }
  
  // Arc was not found
  return false;
}

// Self-explanatory functions to query network info
template <typename _capType> bool ParametricMaxflow<_capType>::IsSimple() { return (flags & GMFF_SIMPLE); }
template <typename _capType> bool ParametricMaxflow<_capType>::IsParametric() { return (flags & GMFF_PARAMETRIC); }
template <typename _capType> void ParametricMaxflow<_capType>::SetOptimize(bool on) { flags = (on ? (flags | GMFF_OPTIMIZE) : (flags & (~GMFF_OPTIMIZE))); }

template <typename _capType> int ParametricMaxflow<_capType>::GetNodeID(int n, bool ins)
{
  int i, f;
  
  if (!(flags & (GMFF_CREATING | GMFF_CREATED)) || !node_table) return -1;
  
  for (i = 0, f = -1; i < _n; ++i)
  {
    if (*(node_table + i) == n) return i;
    else if (*(node_table + i) == -1 && f == -1) f = i;
  }
  
  if (f != -1 && ins)
  {
    *(node_table + f) = n;
    return f;
  }
  
  return -1;
}
