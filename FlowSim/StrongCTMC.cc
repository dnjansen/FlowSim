#include "Strong.h"

// Build the initial map of the relation as a 2D array of booleans.
// Pairs are calculated based on labels.
int StrongSimulation_CTMC::BuildRelationMap()
{
  int m, n, size = 0;
  double *psums = 0;
  
  // Conmpute R(s, S) for all s and normalize transition rates
  psums = new double[n_states];
  for (int j, i = 0; i < n_states; ++i)
  {
    psums[i] = 0.0;
    for (j = row_starts[i]; j < row_starts[i+1]; ++j) psums[i] += non_zeros[j];
    for (j = row_starts[i]; j < row_starts[i+1]; ++j) non_zeros[j] /= psums[i];
  }

  // Initialize relation map
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
  
  rmap.Commit();
  
  return size;
}
