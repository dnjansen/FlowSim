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



#ifndef _RMAP_H_
#define _RMAP_H_

#include <climits>
#include <cstdlib>
#include <cassert>

#include "stats.h"

// Store a relation map as a table of bitsets. Changes have to be committed before
// they are returned by the () operator.
class RelationMap
{
public:
  RelationMap() { buffer[0] = 0, buffer[1] = 0, report = 0; }
  ~RelationMap() { clear_mem(); }

  void clear_mem()
  {
    if (NULL != buffer[0])
    {
      delete [] buffer[0];
      RegisterMemFree(size * sizeof(**buffer));
      buffer[0] = NULL;
    }
    if (NULL != buffer[1])
    {
      delete [] buffer[1];
      RegisterMemFree(size * sizeof(**buffer));
      buffer[1] = NULL;
    }
  }
  
  // Allocate a map with s^2 cells. Attempting to access cells outside of that
  // range will cause an assertion failure if ub is less than s. If ub is greater
  // than s, accessing cells between s and ub will cause more memory to be allocated
  // on demand.  Accessing cells greater than ub will also cause an assertion failure.
  void Create(unsigned int s, unsigned int ub = 0, bool init = false)
  {
    if (s >= (1U << (sizeof(size) * CHAR_BIT / 2)))
    {
      // Strictly speaking, a slightly larger s or ub might be acceptable
      // if one calculates size and upperbound in a special way.
      // (However, the square of any index should fit into an unsigned int. See
      // _index() below.)
      fprintf(stderr, "Error: RelationMap::Create(%u,%u), size too large\n", s,
			ub);
      exit(EXIT_FAILURE);
    }
    if (ub >= (1U << (sizeof(size) * CHAR_BIT / 2)))
    {
      fprintf(stderr, "Warning: RelationMap::Create(%u,%u), upper bound too "
			"large; changed to %u\n", s, ub,
			(1U << (sizeof(size) * CHAR_BIT / 2)) - 1);
      ub = (1U << (sizeof(size) * CHAR_BIT / 2)) - 1;
    }
    if (NULL != buffer[0])
    {
      delete [] buffer[0];
      RegisterMemFree(size * sizeof(**buffer));
    }
    if (NULL != buffer[1])
    {
      delete [] buffer[1];
      RegisterMemFree(size * sizeof(**buffer));
    }
    size = (CHAR_BIT - 1 + s * s) / CHAR_BIT;
    upperbound = (CHAR_BIT - 1 + ub * ub) / CHAR_BIT;
    if (report)
    {
      buffer[0] = NULL;
      RegisterMemAlloc(size * sizeof(**buffer));
    }
    else
    {
      RegisterMemAlloc(size * (sizeof(**buffer) * 2));
      buffer[0] = new unsigned char[size];
      memset(buffer[0], init ? (char) -1 : '\0', sizeof(unsigned char) * size);
    }
    buffer[1] = new unsigned char[size];
    memset(buffer[1], init ? (char) -1 : '\0', sizeof(unsigned char) * size);
  }
  
  // Set whether operator() reports the working state or the committed state
  void ReportCurrent(bool b)
  {
    if (b)
      report = 1;
    else
    {
      if (NULL == buffer[0] && NULL != buffer[1])
      {
        RegisterMemAlloc(size * sizeof(**buffer));
        buffer[0] = new unsigned char[size];
        memset(buffer[0], 0, size);
      }
      report = 0;
    }
  }
  
  // Set/Clear operations
  void Set(unsigned int x, unsigned int y)
  {
    unsigned int i = _index(x, y);
    _assertbyte(i / CHAR_BIT);
    buffer[1][i / CHAR_BIT] |= 1 << (i % CHAR_BIT);
  }
  void Clear(unsigned int x, unsigned int y)
  {
    unsigned int i = _index(x, y);
    if (i / CHAR_BIT < size)
      buffer[1][i / CHAR_BIT] &= ~(1 << (i % CHAR_BIT));
  }
  
  // Update base state with working state
  void Commit()
  {
    if (NULL == buffer[0])
    {
      /* if (report) return; */
      RegisterMemAlloc(size * sizeof(**buffer));
      buffer[0] = new unsigned char[size];
    }
    memcpy(buffer[0], buffer[1], size);
  }
  
  // Check if changes have been made in the working copy since the last commit
  bool MapChanged() {
    return NULL != buffer[0] &&
	       memcmp(buffer[0], buffer[1], sizeof(unsigned char) * size) != 0;
  }
  
  // Retrieve the value of a cell. If out of range (even if below upper bound), false is returned.
  bool operator()(unsigned int x, unsigned int y)
  {
    unsigned int i = _index(x, y);
    if (i / CHAR_BIT >= size) return false;
    return buffer[report][i / CHAR_BIT] & (1 << (i % CHAR_BIT));
  }
  
#ifdef DEBUG
  // Return memory used by this class instance
  void CollectStats(SimulationStatistics *stats)
  {
    size_t memsiz = sizeof(RelationMap)
                + size * (sizeof(**buffer) << (NULL == buffer[0] ? 0 : 1));
    stats->mem_relation_map += memsiz;
    stats->mem_model -= memsiz;
  }
#endif//DEBUG
  
private:
  unsigned char *buffer[2];
  unsigned int size, upperbound, report;
  
  // Map x,y-coordinates to a one-dimensional index in a way that is independent of the map size
  static inline unsigned int _index(unsigned int x, unsigned int y)
  { unsigned int j = (x >= y ? x : y); return (j * j) + (x == j ? y : j + j - x); }
  
  // Assert that a byte-index (!) is valid and allocate more memory if necessary
  inline void _assertbyte(unsigned int i)
  {
    unsigned char *tmp;
    if (i < size) return;
    assert(size < upperbound && i < upperbound);
    
    if (NULL != buffer[0])
    {
      RegisterMemAlloc((i + 1) * sizeof(**buffer));
      tmp = new unsigned char[i + 1];
      memcpy(tmp, buffer[0], sizeof(unsigned char) * size);
      memset(tmp + size, 0, i + 1 - size);
      delete [] buffer[0];
      RegisterMemFree(size * sizeof(**buffer));
      buffer[0] = tmp;
    }
    
    RegisterMemAlloc((i + 1) * sizeof(**buffer));
    tmp = new unsigned char[i + 1];
    memcpy(tmp, buffer[1], sizeof(unsigned char) * size);
    memset(tmp + size, 0, i + 1 - size);
    delete [] buffer[1];
    RegisterMemFree(size * sizeof(**buffer));
    buffer[1] = tmp;
    
    size = i + 1;
  }
};

#endif//_RMAP_H_
