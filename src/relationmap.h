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



#ifndef _RMAP_H_
#define _RMAP_H_

#include <string.h>

// Store a relation map as a table of bitsets. Changes have to be committed before
// they are returned by the () operator.
class RelationMap
{
public:
  RelationMap() { buffer[0] = 0, buffer[1] = 0, report = 0; }
  ~RelationMap() { if (buffer[0]) delete [] buffer[0]; if (buffer[1]) delete [] buffer[1]; }
  
  // Allocate a map with s^2 cells. Attempting to access cells outside of that
  // range will cause an assertion failure if ub is less than s. If ub is greater
  // than s, accessing cells between s and ub will cause more memory to be allocated
  // on demand.  Accessing cells greater than ub will also cause an assertion failure.
  void Create(unsigned int s, unsigned int ub = 0)
  {
    if (buffer[0]) delete [] buffer[0];
    if (buffer[1]) delete [] buffer[1];
    size = ((s * s) >> 3) + (((s * s) & 0x7) ? 1 : 0);
    upperbound = ((ub * ub) >> 3) + (((ub * ub) & 0x7) ? 1 : 0);
    buffer[0] = new unsigned char[size];
    buffer[1] = new unsigned char[size];
    memset(buffer[0], 0, sizeof(unsigned char) * size);
    memset(buffer[1], 0, sizeof(unsigned char) * size);
    report = 0;
  }
  
  // Set whether operator() reports the working state or the committed state
  void ReportCurrent(bool b) { report = (b ? 1 : 0); }
  
  // Set/Clear operations
  void Set(unsigned int x, unsigned int y)
  { unsigned int i = _index(x, y); _assertbyte(i >> 3); buffer[1][i >> 3] |=  (1 << (i & 0x7)); }
  void Clear(unsigned int x, unsigned int y)
  { unsigned int i = _index(x, y); _assertbyte(i >> 3); buffer[1][i >> 3] &= ~(1 << (i & 0x7)); }
  
  // Update base state with working state
  void Commit() { memcpy(buffer[0], buffer[1], sizeof(unsigned char) * size); }
  
  // Check if changes have been made in the working copy since the last commit
  bool MapChanged() { return (memcmp(buffer[0], buffer[1], sizeof(unsigned char) * size) != 0); }
  
  // Retrieve the value of a cell. If out of range (even if below upper bound), false is returned.
  bool operator()(unsigned int x, unsigned int y)
  {
    unsigned int i = _index(x, y);
    if ((i >> 3) >= size) return false;
    return buffer[report][i >> 3] & (1 << (i & 0x7));
  }
  
  // Return memory used by this class instance
  unsigned long MemoryUsage() { return sizeof(RelationMap) + (size * 2 * sizeof(unsigned char)); }
  
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
    
    tmp = new unsigned char[i + 1];
    memcpy(tmp, buffer[0], sizeof(unsigned char) * size);
    memset(tmp + size, 0, i + 1 - size);
    delete [] buffer[0];
    buffer[0] = tmp;
    
    tmp = new unsigned char[i + 1];
    memcpy(tmp, buffer[1], sizeof(unsigned char) * size);
    memset(tmp + size, 0, i + 1 - size);
    delete [] buffer[1];
    buffer[1] = tmp;
    
    size = i + 1;
  }
};

#endif//_RMAP_H_
