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



#include <string.h>
#include <math.h>
#include "prmodel.h"

int main(int argc, char *argv[])
{
  int n, m, i, l, type, mina, maxa, mins, maxs, *svector, *avector;
  double avgs, min_trans, max_trans, max_delta, dist, min_dist, max_dist;
  FILE *f = 0;
  char *p = 0;
  bool continuous;
  ProbabilisticModel *pm = 0;
  MarkovChain *mc = 0;
  ProbabilisticAutomaton *pa = 0;
  
  if (argc < 2)
  {
    printf("Usage: modelstat [[type:]model ...]\n\n");
    printf("Prints statistics about a probabilistic model (DTMC, CTMC, PA, CPA).\n");
    printf("Model type is determined from file extension; file name must be\n");
    printf("prefixed by model type if the file extension doesn't contain it.\n");
    printf("Example:   modelstat dtmc:mymodel.pm\n");
    printf("     or:   modelstat mymodel.dtmc\n\n");
    return 0;
  }
  
  for (n = 1; n < argc; ++n)
  {
    if (!strncmp(argv[n], "dtmc:", 5))
    {
      p = argv[n] + 5;
      type = 1;
      pm = new MarkovChain;
      continuous = false;
    }
    else if (!strncmp(argv[n], "ctmc:", 5))
    {
      p = argv[n] + 5;
      type = 2;
      pm = new MarkovChain;
      continuous = true;
    }
    else if (!strncmp(argv[n], "pa:", 3))
    {
      p = argv[n] + 3;
      type = 3;
      pm = new ProbabilisticAutomaton;
      continuous = false;
    }
    else if (!strncmp(argv[n], "cpa:", 4))
    {
      p = argv[n] + 4;
      type = 4;
      pm = new ProbabilisticAutomaton;
      continuous = true;
    }
    else if (!strcmp(argv[n] + strlen(argv[n]) - 5, ".dtmc"))
    {
      p = argv[n];
      type = 1;
      pm = new MarkovChain;
      continuous = false;
    }
    else if (!strcmp(argv[n] + strlen(argv[n]) - 5, ".ctmc"))
    {
      p = argv[n];
      type = 2;
      pm = new MarkovChain;
      continuous = true;
    }
    else if (!strcmp(argv[n] + strlen(argv[n]) - 3, ".pa"))
    {
      p = argv[n];
      type = 3;
      pm = new ProbabilisticAutomaton;
      continuous = false;
    }
    else if (!strcmp(argv[n] + strlen(argv[n]) - 4, ".cpa"))
    {
      p = argv[n];
      type = 4;
      pm = new ProbabilisticAutomaton;
      continuous = true;
    }
    else
    {
      printf("%s: model type can't be determined\n", argv[n]);
      continue;
    }
    
    f = fopen(p, "rb");
    pm->Parse(f);
    fclose(f);
    
    printf("Model:              %s\n", p);
    
    if (type == 1 || type == 2)
    {
      mc = (MarkovChain*)pm;
    
      printf("States:             %d\n", mc->n);
      printf("Transitions:        %ld\n", mc->nnz);
    
      mins = mc->n;
      maxs = 0;
      avgs = 0;
      min_trans = -1.0;
      max_delta = 0.0;
      max_trans = 0.0;
      min_dist = -1.0;
      max_dist = 0.0;
      for (m = 0; m < mc->n; ++m)
      {
        for (i = mc->row_starts[m], dist = 0.0; i < mc->row_starts[m + 1]; ++i)
        {
          if (mc->non_zeros[i] > 0.0 && (mc->non_zeros[i] < min_trans || min_trans < 0.0)) min_trans = mc->non_zeros[i];
          if (mc->non_zeros[i] > max_trans) max_trans = mc->non_zeros[i];
          dist += mc->non_zeros[i];
        }
        if (dist > 0.0 && (dist < min_dist || min_dist < 0.0)) min_dist = dist;
        if (dist > max_dist) max_dist = dist;
        if (!continuous)
        {
          dist = fabs(1.0 - dist);
          if (dist > max_delta) max_delta = dist;
        }
        l = mc->row_starts[m + 1] - mc->row_starts[m];
        if (l < mins) mins = l;
        if (l > maxs) maxs = l;
        avgs += l;
      }
      avgs /= mc->n;
    
      printf("Minimum successors: %d\n", mins);
      printf("Maximum successors: %d\n", maxs);
      printf("Average successors: %.2f\n\n", avgs);
    
      svector = new int[maxs + 1];
      for (m = 0; m <= maxs; ++m) svector[m] = 0;
    
      for (m = 0; m < mc->n; ++m) svector[mc->row_starts[m + 1] - mc->row_starts[m]]++;
    
      for (m = 0; m <= maxs; ++m) if (svector[m] > 0) printf("States with %4d successors: %d\n", m, svector[m]);
    
      delete [] svector;
      
      printf("\n");
      
      if (continuous)
      {
        printf("Lowest non-zero transition rate:                %.4E\n", min_trans);
        printf("Highest transition rate:                        %.4E\n", max_trans);
        printf("Lowest non-zero rate (entire distribution):     %.4E\n", min_dist);
        printf("Highest rate (entire distribution):             %.4E\n", max_dist);
        
        printf("\n");
        printf("NOTE: The floating point approximation threshold cannot\n");
        printf("      be computed automatically for continuous-time models,\n");
        printf("      but may still be required for correct operation. The\n");
        printf("      upper bound for the threshold is: %.4E\n", min_trans / max_dist);
        printf("      The lower bound depends upon the model's precision.\n");
      }
      else
      {
        printf("Lowest non-zero transition probability:         %.4E\n", min_trans);
        printf("Highest deviation from stochastic distribution: ");
        if (max_delta == 0.0)
        {
          printf("none\n");
          printf("Recommended floating point approximation:       none needed\n");
        }
        else
        {
          printf("%.4E\n", max_delta);
          
          if (max_delta > min_trans) printf("WARNING: Either this model is not stochastic, or the floating point\n"
                                            "         precision is too low and too much information has been lost\n"
                                            "         to rounding errors. If the model is not stochastic, you\n"
                                            "         must manually specify a suitable approximation threshold;\n"
                                            "         otherwise, the simulation relation cannot be computed on\n"
                                            "         this model due to large rounding errors.\n");
                                                                    // Harmonic mean favors the smaller value
          else printf("Recommended floating point approximation:       %.4E\n", 10 / ((9 / min_trans) + (1 / max_delta)));
        }
      }
      
      printf("\n");
    }
    else if (type == 3 || type == 4)
    {
      pa = (ProbabilisticAutomaton*)pm;
    
      printf("States:             %d\n", pa->n);
      printf("Actions:            %d\n", pa->na);
      printf("Transitions:        %d\n", pa->nnz);
    
      mina = pa->n;
      maxa = 0;
      avgs = 0;
      for (m = 0; m < pa->n; ++m)
      {
        l = pa->state_starts[m + 1] - pa->state_starts[m];
        if (l < mina) mina = l;
        if (l > maxa) maxa = l;
        avgs += l;
      }
      avgs /= pa->n;
    
      printf("Minimum actions:    %d\n", mina);
      printf("Maximum actions:    %d\n", maxa);
      printf("Average actions:    %.2f\n", avgs);
    
      mins = pa->n;
      maxs = 0;
      avgs = 0;
      min_trans = -1.0;
      max_trans = 0.0;
      min_dist = -1.0;
      max_dist = 0.0;
      max_delta = 0.0;
      for (m = 0; m < pa->n; ++m)
      {
        for (i = pa->state_starts[m]; i < pa->state_starts[m+1]; ++i)
        {
          for (l = pa->row_starts[i], dist = 0.0; l < pa->row_starts[i+1]; ++l)
          {
            if (pa->non_zeros[l] > 0.0 && (pa->non_zeros[l] < min_trans || min_trans < 0.0)) min_trans = pa->non_zeros[l];
            if (pa->non_zeros[l] > max_trans) max_trans = pa->non_zeros[l];
            dist += pa->non_zeros[l];
          }
          if (dist > 0.0 && (dist < min_dist || min_dist < 0.0)) min_dist = dist;
          if (dist > max_dist) max_dist = dist;
          if (!continuous)
          {
            dist = fabs(1.0 - dist);
            if (dist > max_delta) max_delta = dist;
          }
          l = pa->row_starts[i + 1] - pa->row_starts[i];
          if (l < mins) mins = l;
          if (l > maxs) maxs = l;
          avgs += l;
        }
      }
      avgs /= pa->na;
    
      printf("Minimum successors: %d\n", mins);
      printf("Maximum successors: %d\n", maxs);
      printf("Average successors: %.2f\n", avgs);
      
      avector = new int[maxa + 1];
      for (m = 0; m <= maxa; ++m) avector[m] = 0;
      svector = new int[maxs + 1];
      for (m = 0; m <= maxs; ++m) svector[m] = 0;
    
      for (m = 0; m < pa->n; ++m)
      {
        avector[pa->state_starts[m + 1] - pa->state_starts[m]]++;
        for (i = pa->state_starts[m]; i < pa->state_starts[m + 1]; ++i)
        {
          svector[pa->row_starts[i + 1] - pa->row_starts[i]]++;
        }
      }
      
      printf("\n");
      for (m = 0; m <= maxa; ++m) if (avector[m] > 0) printf("States with %4d actions: %d\n", m, avector[m]);
      
      printf("\n");
      for (m = 0; m <= maxs; ++m) if (svector[m] > 0) printf("Actions with %4d successors: %d\n", m, svector[m]);
      
      delete [] avector;
      delete [] svector;
      
      printf("\n");
      
      if (continuous)
      {
        printf("Lowest non-zero transition rate:                %.4E\n", min_trans);
        printf("Highest transition rate:                        %.4E\n", max_trans);
        printf("Lowest non-zero rate (entire distribution):     %.4E\n", min_dist);
        printf("Highest rate (entire distribution):             %.4E\n", max_dist);
        
        printf("\n");
        printf("NOTE: The floating point approximation threshold cannot\n");
        printf("      be computed automatically for continuous-time models,\n");
        printf("      but may still be required for correct operation. The\n");
        printf("      upper bound for the threshold is: %.4E\n", min_trans / max_dist);
        printf("      The lower bound depends upon the model's precision.\n");
      }
      else
      {
        printf("Lowest non-zero transition probability:         %.4E\n", min_trans);
        printf("Highest deviation from stochastic distribution: ");
        if (max_delta == 0.0)
        {
          printf("none\n");
          printf("Recommended floating point approximation:       none needed\n");
        }
        else
        {
          printf("%.4E\n", max_delta);
        
          if (max_delta > min_trans) printf("WARNING: The floating point precision on this model is too low. Data\n"
                                            "         has been lost due to rounding errors. Attempting to compute\n"
                                            "         simulation on this model will provide incorrect results.\n");
                                                                    // Harmonic mean favors the smaller value
          else printf("Recommended floating point approximation:       %.4E\n", 10 / ((9 / min_trans) + (1 / max_delta)));
        }
      }
    
      printf("\n");
    }
    
    delete pm;
    pa = 0;
  }
  
  return 0;
}
