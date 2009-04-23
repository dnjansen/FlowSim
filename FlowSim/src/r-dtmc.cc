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



#include <set>
#include <vector>
#include <list>
#include <assert.h>
#include <math.h>
#include "prmodel.h"
#include "modelbuilder.h"

/*
Generate a random DTMC and return it as a sparse matrix
* n : number of nodes
* a : minimum number of successors
* b : maximum number of successors

Probabilities are uniformly distributed between successor states and
P(s,S)=1 for all states.
*/
MarkovChain *DTMC(int n, int a, int b)
{
  int i, j, r, row_start, transitions;
  std::set<int> successors;
  std::set<int>::iterator s_i;
  std::vector<int> cols;
  std::vector<double> non_zeros;
  MarkovChain *dtmc = new MarkovChain;
  
  srand(rand() ^ time(0));
  
  dtmc->n = n;
  dtmc->row_starts = new int[n + 1];
  
  // Generate successor sets of random successors from the state set (no probabilities yet)
  row_start = 0;
  for (i = 0; i < n; ++i)
  {
    dtmc->row_starts[i] = row_start;
    r = a + (rand() % (b - a + 1));
    successors.clear();
    while ((int)successors.size() < r) successors.insert(rand() % n);
    for (s_i = successors.begin(); s_i != successors.end(); ++s_i)
    {
      cols.push_back(*s_i);
      ++row_start;
    }
  }
  
  dtmc->row_starts[i] = row_start;
  
  // Setup uniform distribution of transition probabilities for all states
  for (i = 0; i < n; ++i)
  {
    transitions = dtmc->row_starts[i+1] - dtmc->row_starts[i];
    for (j = 0; j < transitions; ++j) non_zeros.push_back(1.0 / transitions);
  }
  
  assert(cols.size() == non_zeros.size());
  
  // Transfer data into sparse matrix
  dtmc->nnz = non_zeros.size();
  dtmc->cols = new int[dtmc->nnz];
  dtmc->non_zeros = new double[dtmc->nnz];
  
  for (i = 0; i < dtmc->nnz; ++i)
  {
    dtmc->cols[i] = cols[i];
    dtmc->non_zeros[i] = non_zeros[i];
  }
  
  return dtmc;
}

/*
Generate a random DTMC and return it as a sparse matrix
* n : number of nodes
* p : the probability that there is an edge between node i and j for 0 <= i,j < n

Probabilities are uniformly distributed between successor states and
P(s,S)=1 for all states.
*/
MarkovChain *DTMC(int n, double p)
{
  int i, j, row_start, transitions, pp = int(p * 1000);
  std::vector<int> cols;
  std::vector<double> non_zeros;
  MarkovChain *dtmc = new MarkovChain;
  
  srand(time(0));
  
  dtmc->n = n;
  dtmc->row_starts = new int[n + 1];
  
  // Decide which states will have transitions to which other states
  row_start = 0;
  for (i = 0; i < n; ++i)
  {
    dtmc->row_starts[i] = row_start;
    for (j = 0; j < n; ++j)
    {
      // Uncomment the following line if states should not have transitions back to themselves
      //if (i == j) continue;
      if (1 + (rand() % 1000) <= pp)
      {
        cols.push_back(j);
        ++row_start;
      }
    }
  }
  
  dtmc->row_starts[n] = row_start;
  
  // Setup uniform distribution of transition probabilities for all states
  for (i = 0; i < n; ++i)
  {
    transitions = dtmc->row_starts[i+1] - dtmc->row_starts[i];
    for (j = 0; j < transitions; ++j) non_zeros.push_back(1.0 / transitions);
  }
  
  assert(cols.size() == non_zeros.size());
  
  // Transfer data into sparse matrix
  dtmc->nnz = non_zeros.size();
  dtmc->cols = new int[dtmc->nnz];
  dtmc->non_zeros = new double[dtmc->nnz];
  
  for (i = 0; i < dtmc->nnz; ++i)
  {
    dtmc->cols[i] = cols[i];
    dtmc->non_zeros[i] = non_zeros[i];
  }
  
  return dtmc;
}

/*
 * Generate a random DTMC with certain structural characteristics. If
 * all bias values are set to 0, the function behaves like DTMC(n, a, b).
 * Note that the bias values represent tendencies which may be overriden
 * by constraints, e.g. if a state is to have 10 successors but there are
 * only 8 successors which will keep the model linear, the model will
 * be non-linear even if lbias is 1 in order to satisfy 10 successors.
 *
 * n    : number of states
 * a    : minimum number of successors per state
 * b    : maximum number of successors per state
 * fbias: fanout bias; 0.0=random, -1.0=always minimum # of suc, 1.0=always maximum # of suc
 * lbias: linearity bias; 0.0=random, 1.0=linear (no cycles)
 * c    : number of clusters in model (only effective if cbias != 0.0)
 * cbias: cluster bias; 0.0=no clusters, 1.0='c' isolated clusters of states,
 *      : in between: there may be transitions between clusters
 * pbias: probability bias; 0.0=uniform distributions, 1.0=biased distr.
 * sbias: successor bias; 0.0=random successors chosen, 1.0=biased towards some states
 */
MarkovChain *DTMC(int n, int a, int b, double fbias, double lbias, int c, double cbias, double pbias, double sbias)
{
  double *matrix = new double[n*n];
  bool *stencil = new bool[n*n];
  int i, j, k, l, prevcluster, possiblechoices, fanout, subdist, curdist, dstart, dend;
  double clusterexcess, lambda;
  int *clusters = new int[c], *stateorder = new int[n];
  std::list<int> statelist;
  MarkovChain *pmc = 0;
  ModelBuilder mb;

  srand(rand() ^ time(0));

  if (fbias <-1.0) fbias =-1.0;
  if (lbias < 0.0) fbias = 0.0;
  if (cbias < 0.0) fbias = 0.0;
  if (pbias < 0.0) fbias = 0.0;
  if (sbias < 0.0) fbias = 0.0;

  if (fbias > 1.0) fbias = 1.0;
  if (lbias > 1.0) fbias = 1.0;
  if (cbias > 1.0) fbias = 1.0;
  if (pbias > 1.0) fbias = 1.0;
  if (sbias > 1.0) fbias = 1.0;

  // Initialize matrix and generate cluster sizes randomly for c clusters
  for (i = 0; i < n * n; ++i) matrix[i] = 0.0, stencil[i] = true;
  for (i = 0, prevcluster = 0; i < c - 1; ++i) prevcluster += (clusters[i] = (rand() % (n - prevcluster - (c - i - 1))) + 1);
  clusters[c - 1] = n - prevcluster;
  
  // If cluster bias is non-zero, begin clustering
  if (cbias > 0.0)
  {
    // Simplified operation for absolute clustering, clear all inter-cluster transitions in the stencil
    if (cbias == 1.0)
    {
      prevcluster = 0;
      for (i = 0; i < c; ++i)
      {
        for (j = 0; j < prevcluster; ++j)
        {
          for (k = prevcluster; k < n; ++k)
          {
            stencil[(j * n) + k] = false;
          }
        }
        for (j = prevcluster; j < n; ++j)
        {
          for (k = 0; k < prevcluster; ++k)
          {
            stencil[(j * n) + k] = false;
          }
        }
        prevcluster += clusters[i];
      }
    }
    else // Cluster bias is less than 1, use bias as a guideline on how many transitions to clear in the stencil
    {
      prevcluster = 0;
      for (i = 0; i < c - 1; ++i)
      {
        prevcluster += clusters[i];
        clusterexcess = (n - prevcluster) * clusters[i] * 2 * cbias;
        while (clusterexcess >= 1.0)
        {
          if (rand() >= (RAND_MAX >> 1))
          {
            j = prevcluster + (rand() % (n - prevcluster));
            k = prevcluster + (rand() % clusters[i]) - clusters[i];
          }
          else
          {
            j = prevcluster + (rand() % clusters[i]) - clusters[i];
            k = prevcluster + (rand() % (n - prevcluster));
          }
          if (!stencil[(j * n) + k]) continue;
          stencil[(j * n) + k] = false;
          clusterexcess -= 1.0;
        }
        if (rand() < clusterexcess * RAND_MAX)
        {
          do
          {
            if (rand() >= (RAND_MAX >> 1))
            {
              j = prevcluster + (rand() % (n - prevcluster));
              k = prevcluster + (rand() % clusters[i]) - clusters[i];
            }
            else
            {
              j = prevcluster + (rand() % clusters[i]) - clusters[i];
              k = prevcluster + (rand() % (n - prevcluster));
            }
          }
          while (!stencil[(j * n) + k]);
          stencil[(j * n) + k] = false;
        }
      }
    }
  }
  
  // If there is a linearity bias, fill out transition matrix with linearity constraints
  if (lbias > 0.0)
  {
    for (i = 0; i < n; ++i)
    {
      for (possiblechoices = 0, j = 0; j < n; ++j) if (stencil[(i * n) + j]) ++possiblechoices;
      if (lbias == 1.0 && possiblechoices > n - i - 1) possiblechoices = n - i - 1;

      if (fbias == -1.0) fanout = a;
      else if (fbias == 1.0) fanout = b;
      else if (fbias <= 0.0)
      {
        lambda = (5 * -fbias) + 1;
        fanout = (int)round(((b - a) * pow((double)rand() / RAND_MAX, lambda)) + a);
      }
      else
      {
        lambda = (5 * fbias) + 1;
        fanout = (int)round(b - ((b - a) * pow((double)rand() / RAND_MAX, lambda)));
      }
      
      if (fanout < a) fanout = a;
      if (fanout > b) fanout = b;
      
      while (fanout > 0 && possiblechoices > 0)
      {
        if (rand() > lbias * RAND_MAX) j = (int)round(pow((double)rand() / RAND_MAX, (4.0 * sbias) + 1) * n);
        else j = (int)round(pow((double)rand() / RAND_MAX, (4.0 * sbias) + 1) * (n - i - 1)) + i + 1;
        if (j < 0 || j >= n || !stencil[(i * n) + j] || matrix[(i * n) + j] > 0.0) continue;
        matrix[(i * n) + j] = 1.0;
        --fanout;
        --possiblechoices;
      }
    }
  }
  else // Zero linearity bias
  {
    for (i = 0; i < n; ++i)
    {
      for (possiblechoices = 0, j = 0; j < n; ++j) if (stencil[(i * n) + j]) ++possiblechoices;

      if (fbias == -1.0) fanout = a;
      else if (fbias == 1.0) fanout = b;
      else if (fbias <= 0.0)
      {
        lambda = (5 * -fbias) + 1;
        fanout = (int)round(((b - a) * pow((double)rand() / RAND_MAX, lambda)) + a);
      }
      else
      {
        lambda = (5 * fbias) + 1;
        fanout = (int)round(b - ((b - a) * pow((double)rand() / RAND_MAX, lambda)));
      }
      
      if (fanout < a) fanout = a;
      if (fanout > b) fanout = b;
      
      while (fanout > 0 && possiblechoices > 0)
      {
        j = (int)round(pow((double)rand() / RAND_MAX, (4.0 * sbias) + 1) * n);
        if (j < 0 || j >=n || !stencil[(i * n) + j] || matrix[(i * n) + j] > 0.0) continue;
        matrix[(i * n) + j] = 1.0;
        --fanout;
        --possiblechoices;
      }
    }
  }
  
  // Generate a random permutation of the states
  for (i = 0; i < n; ++i)
  {
    if (rand() >= (RAND_MAX >> 1)) statelist.push_back(i);
    else statelist.push_front(i);
  }
  for (i = 0; i < n; ++i)
  {
    stateorder[i] = statelist.front();
    statelist.pop_front();
  }
  
  // Compute transition probabilities
  for (i = 0; i < n; ++i)
  {
    k = 0;
    for (j = 0; j < n; ++j) if (matrix[(i * n) + j] == 1.0) ++k;
    if (pbias == 0.0) subdist = 1;
    else subdist = (int)ceil(pow(double(rand()) / RAND_MAX, 1.0 / pbias) * k);
    dstart = 0, dend = 0, curdist = 0;
    for (j = 0, l = 0; j < n; ++j)
    {
      if (matrix[(i * n) + j] > 0.0)
      {
        if (l == dend)
        {
          ++curdist;
          dstart = dend;
          if (dend == (int)floor(1.0 / subdist * curdist * k)) dend = (int)ceil(1.0 / subdist * curdist * k);
          else dend = (int)((rand() > (RAND_MAX >> 1)) ? floor(1.0 / subdist * curdist * k) : ceil(1.0 / subdist * curdist * k));
        }
        matrix[(i * n) + j] = 2.0 / (subdist * (subdist + 1)) * curdist / (dend - dstart);
        ++l;
      }
    }
  }

  // Create sparse matrix  
  for (i = 0, k = 0; i < n; ++i)
  {
    for (j = 0; j < n; ++j)
    {
      if (matrix[(i * n) + j] > 0.0) ++k;
    }
  }
  
  for (i = 0; i < n; ++i)
  {
    for (j = 0; j < n; ++j)
    {
      if (matrix[(i * n) + j] > 0.0)
      {
        mb.Add(stateorder[i], stateorder[j], matrix[(i * n) + j]);
      }
    }
  }
  
  pmc = mb.BuildMC(n, 0);
  
  // Free memory and return markov chain
  delete [] stateorder;
  delete [] matrix;
  delete [] clusters;
  delete [] stencil;
  
  return pmc;
}

// Compile main function only if r-dtmc is compiled as a standalone program and not linked elsewhere
#ifdef _STANDALONE_
int main(int argc, char *argv[])
{
  if (argc != 3 && argc != 4 && argc != 10)
  {
    printf("Usage: r-dtmc N P\n");
    printf("Usage: r-dtmc N A B\n");
    printf("Usage: r-dtmc N A B C Fb Lb Cb Pb Sb\n");
    printf("\n");
    printf("   N   number of states, >0\n");
    printf("   P   probability of an edge existing between two states, [0,1]\n");
    printf("   A   minimum number of successors per state, A >= 0\n");
    printf("   B   maximum number of successors per state, B >= A\n");
    printf("   C   number of clusters in model, C >= 1\n");
    printf("   Fb  bias towards A or B successors per state, [-1,1]\n");
    printf("   Lb  bias towards linear model (acyclic), [0,1]\n");
    printf("   Cb  bias towards clustering, [0,1]\n");
    printf("   Pb  probability distribution bias, 0=uniform, [0,1]\n");
    printf("   Sb  bias towards certain states being picked as successors\n");
    printf("       more often than others, [0,1]\n");
    printf("\n");
    return -1;
  }
  
  int n = strtol(argv[1], 0, 10);
  double p = strtod(argv[2], 0);
  int a = strtol(argv[2], 0, 10);
  int b = 0;
  int c = 0;
  double fb = 0, lb = 0, cb = 0, pb = 0, sb = 0;
  if (argc >= 4) b = strtol(argv[3], 0, 10);
  if (argc == 10)
  {
    c = strtol(argv[4], 0, 10);
    fb = strtod(argv[5], 0);
    lb = strtod(argv[6], 0);
    cb = strtod(argv[7], 0);
    pb = strtod(argv[8], 0);
    sb = strtod(argv[9], 0);
  }
  
  if (n <= 0 || (argc == 3 && (p < 0.0 || p > 1.0)) || (argc == 4 && (a < 0 || a > b))
    || (argc == 10 && (a < 0 || a > b || c < 0 || fb < -1.0 || fb > 1.0 || lb < 0.0 || lb > 1.0 || cb < 0.0 || cb > 1.0 || pb < 0.0 || pb > 1.0 || sb < 0.0 || sb > 1.0)))
  {
    printf("Usage: r-dtmc N P\n");
    printf("Usage: r-dtmc N A B\n");
    printf("\n");
    printf("   N   number of states, >0\n");
    printf("   P   probability of an edge existing between two states, [0,1]\n");
    printf("   A   minimum number of successors per state, A >= 0\n");
    printf("   B   maximum number of successors per state, B >= A\n");
    printf("   C   number of clusters in model, C >= 1\n");
    printf("   Fb  bias towards A or B successors per state, [-1,1]\n");
    printf("   Lb  bias towards linear model (acyclic), [0,1]\n");
    printf("   Cb  bias towards clustering, [0,1]\n");
    printf("   Pb  probability distribution bias, 0=uniform, [0,1]\n");
    printf("   Sb  bias towards certain states being picked as successors\n");
    printf("       more often than others, [0,1]\n");
    printf("\n");
    return -1;
  }
  
  MarkovChain *s;
  switch (argc)
  {
  case 3: s = DTMC(n, p); break;
  case 4: s = DTMC(n, a, b); break;
  case 10: s = DTMC(n, a, b, fb, lb, c, cb, pb, sb); break;
  }
  
  s->Write();
  
  delete s;
  
  return 0;
}
#endif
