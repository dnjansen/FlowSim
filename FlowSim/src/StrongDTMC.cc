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



#include "Strong.h"

// Build the initial map of the relation as a 2D array of booleans.
// Pairs are calculated based on labels.
int StrongSimulation_DTMC::BuildRelationMap()
{
  int m, n, size = 0;
  
  //if (optflags & OPT_P_INVARIANT)
#ifdef OPT_P_INVARIANT
  {
    double *psums = new double[n_states];
    for (int j, i = 0; i < n_states; ++i)
    {
      psums[i] = 0.0;
      for (j = row_starts[i]; j < row_starts[i+1]; ++j) psums[i] += non_zeros[j];
    }

    for (m = 0; m < n_states - 1; ++m)
    {
      for (n = m + 1; n < n_states; ++n)
      {
        if (Label(m) == Label(n))
        {
          if (CompactMaxFlow<double>::_Tleq(psums[m], psums[n])) rmap.Set(m, n), ++size;
          if (CompactMaxFlow<double>::_Tleq(psums[n], psums[m])) rmap.Set(n, m), ++size;
        }
      }
    }
    
    delete [] psums;
  }
  //else
#else
  {
    for (m = 0; m < n_states - 1; ++m)
    {
      for (n = m + 1; n < n_states; ++n)
      {
        if (Label(m) == Label(n))
        {
          size += 2;
          rmap.Set(m, n);
          rmap.Set(n, m);
        }
      }
    }
  }
#endif
  
  rmap.Commit();
  
  return size;
}
