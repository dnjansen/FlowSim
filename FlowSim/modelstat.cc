#include <stdio.h>
#include <string.h>
#include <math.h>
#include "prmodel.h"

int main(int argc, char *argv[])
{
  int n, m, i, l, type, mins, maxs, *svector;
  double avgs, min_trans, max_delta, dist;
  FILE *f = 0;
  char *p = 0;
  ProbabilisticModel *pm = 0;
  MarkovChain *mc = 0;
  ProbabilisticAutomaton *pa = 0;
  
  if (argc < 2)
  {
    printf("Usage: modelstat [[type:]model ...]\n\n");
    printf("Prints statistics about a probabilistic model (DTMC, CTMC, PA, CPA).\n");
    printf("Model type is determined from file extension; file name must be\n");
    printf("prefixed by model type if the file extension doesn't contain it.\n");
    printf("Example:   modelstat dtmc:mymodel.pm\n\n");
    return 0;
  }
  
  for (n = 1; n < argc; ++n)
  {
    if (!strncmp(argv[n], "dtmc:", 5))
    {
      p = argv[n] + 5;
      type = 1;
      pm = new MarkovChain;
    }
    else if (!strncmp(argv[n], "ctmc:", 5))
    {
      p = argv[n] + 5;
      type = 2;
      pm = new MarkovChain;
    }
    else if (!strncmp(argv[n], "pa:", 3))
    {
      p = argv[n] + 3;
      type = 3;
      pm = new ProbabilisticAutomaton;
    }
    else if (!strncmp(argv[n], "cpa:", 4))
    {
      p = argv[n] + 4;
      type = 4;
      pm = new ProbabilisticAutomaton;
    }
    else if (!strcmp(argv[n] + strlen(argv[n]) - 5, ".dtmc"))
    {
      p = argv[n];
      type = 1;
      pm = new MarkovChain;
    }
    else if (!strcmp(argv[n] + strlen(argv[n]) - 5, ".ctmc"))
    {
      p = argv[n];
      type = 2;
      pm = new MarkovChain;
    }
    else if (!strcmp(argv[n] + strlen(argv[n]) - 3, ".pa"))
    {
      p = argv[n];
      type = 3;
      pm = new ProbabilisticAutomaton;
    }
    else if (!strcmp(argv[n] + strlen(argv[n]) - 4, ".cpa"))
    {
      p = argv[n];
      type = 4;
      pm = new ProbabilisticAutomaton;
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
      min_trans = 1.0;
      max_delta = 0.0;
      for (m = 0; m < mc->n; ++m)
      {
        for (i = mc->row_starts[m], dist = 0.0; i < mc->row_starts[m + 1]; ++i)
        {
          if (mc->non_zeros[i] < min_trans) min_trans = mc->non_zeros[i];
          dist += mc->non_zeros[i];
        }
        dist = fabs(1.0 - dist);
        if (dist > max_delta) max_delta = dist;
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
      
      printf("Lowest transition probability:                  %.4E\n", min_trans);
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
      
      printf("\n");
    }
    else if (type == 3 || type == 4)
    {
      pa = (ProbabilisticAutomaton*)pm;
    
      printf("States:             %d\n", pa->n);
      printf("Actions:            %d\n", pa->na);
      printf("Transitions:        %d\n", pa->nnz);
    
      mins = pa->n;
      maxs = 0;
      avgs = 0;
      for (m = 0; m < pa->n; ++m)
      {
        l = pa->state_starts[m + 1] - pa->state_starts[m];
        if (l < mins) mins = l;
        if (l > maxs) maxs = l;
        avgs += l;
      }
      avgs /= pa->n;
    
      printf("Minimum actions:    %d\n", mins);
      printf("Maximum actions:    %d\n", maxs);
      printf("Average actions:    %.2f\n", avgs);
    
      mins = pa->n;
      maxs = 0;
      avgs = 0;
      min_trans = 1.0;
      max_delta = 0.0;
      for (m = 0; m < pa->n; ++m)
      {
        for (i = pa->state_starts[m]; i < pa->state_starts[m+1]; ++i)
        {
          for (l = pa->row_starts[i], dist = 0.0; l < pa->row_starts[i+1]; ++l)
          {
            if (pa->non_zeros[l] < min_trans) min_trans = pa->non_zeros[l];
            dist += pa->non_zeros[l];
          }
          if (dist == 0.0) printf("Zero distribution in state %d, %d\n", m, i);
          dist = fabs(1.0 - dist);
          if (dist > max_delta) max_delta = dist;
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
    
      svector = new int[maxs + 1];
      for (m = 0; m <= maxs; ++m) svector[m] = 0;
    
      for (m = 0; m < pa->n; ++m)
      {
        for (i = pa->state_starts[m]; i < pa->state_starts[m + 1]; ++i)
        {
          svector[pa->row_starts[i + 1] - pa->row_starts[i]]++;
        }
      }
    
      for (m = 0; m <= maxs; ++m) if (svector[m] > 0) printf("Actions with %4d successors: %d\n", m, svector[m]);
    
      delete [] svector;
      
      printf("\n");
      
      printf("Lowest transition probability:                  %.4E\n", min_trans);
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
    
      printf("\n");
    }
    
    delete pm;
    pa = 0;
  }
  
  return 0;
}
