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



#include "prmodel.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <set>
#include <vector>

using namespace std;

MarkovChain::MarkovChain(bool c) : ProbabilisticModel(c)
{
  row_starts = 0;
  cols = 0;
  non_zeros = 0;
  mt = MC;
}

ProbabilisticAutomaton::ProbabilisticAutomaton(bool c) : ProbabilisticModel(c)
{
  row_starts = 0;
  cols = 0;
  non_zeros = 0;
  actions = 0;
  state_starts = 0;
  mt = PA;
  atable = 0;
}

MarkovChain::~MarkovChain()
{
  delete [] row_starts;
  delete [] cols;
  delete [] non_zeros;
}

ProbabilisticAutomaton::~ProbabilisticAutomaton()
{
  delete [] row_starts;
  delete [] state_starts;
  delete [] actions;
  delete [] cols;
  delete [] non_zeros;
  delete [] atable;
}

void MarkovChain::InitMatrix()
{
  if (row_starts) delete [] row_starts;
  if (cols) delete [] cols;
  if (non_zeros) delete [] non_zeros;
  
  row_starts = new int[n+1];
  for(int i=0; i<n+1; i++)
  {
    row_starts[i]=0;
  }
  cols = new int[nnz];
  non_zeros = new double[nnz];
  for(int i=0; i<nnz; i++)
  {
    cols[i] = 0;
    non_zeros[i] = 0.0;
  }
}

void MarkovChain::Parse(FILE *f, bool con_time)
{
  const int MAXLINE = 128;
  char in_line[MAXLINE], *next;
  
  if (!f) f = stdin;

  fgets(in_line, MAXLINE, f);
  sscanf(in_line, "%d %ld", &n, &nnz);
  InitMatrix();
  
  continuous = con_time;

  int source, target;
  double prob;
  long no_lines = 0;
  
  while (fgets(in_line, MAXLINE, f) != NULL)
  {
    switch (in_line[0])
    {
    case 'c':                  /* skip lines with comments */
    case '\n':                 /* skip empty lines   */
    case '\0':                 /* skip empty lines at the end of file */
      break;
    default:                    /* transition description */
      source = strtol(&in_line[0], &next, 10);
      target = strtol(next, &next, 10);
      prob = strtod(next, 0);
      row_starts[source + 1]++;
      cols[no_lines] = target;
      non_zeros[no_lines] = prob;
      no_lines ++;
      break;
    } /* end of switch */
  }   /* end of input loop */

  if (no_lines != nnz)
  {
    printf("There are more transitions than described in the first line, please check!\n");
    exit(1);
  }
  
  // now construct the row_starts
  for(int i=1; i<n+1; i++)
  {
    row_starts[i] += row_starts[i-1];
  }
}

void ProbabilisticAutomaton::InitMatrix()
{
  if (row_starts) delete [] row_starts;
  if (cols) delete [] cols;
  if (non_zeros) delete [] non_zeros;
  if (state_starts) delete [] state_starts;
  if (actions) delete [] actions;
  if (atable) delete [] atable;
  
  row_starts = new int[na+1];
  for(int i=0; i<na+1; i++)
  {
    row_starts[i]=0;
  }
  cols = new int[nnz];
  non_zeros = new double[nnz];
  state_starts = new int[n+1];
  actions = new int[na];
  for(int i=0; i<n+1; i++)
  {
    state_starts[i]=0;
  }
  for(int i=0; i<nnz; i++)
  {
    cols[i] = 0;
    non_zeros[i] = 0.0;
  }
}

void ProbabilisticAutomaton::Parse(FILE *f, bool con_time)
{
  const int MAXLINE = 128;
  char in_line[MAXLINE], *next;
  set<int> seen_actions;
  vector<int> act_table;
  vector<int>::iterator at;
  bool prism_format = true;
  
  if (!f) f = stdin;
  
  continuous = con_time;

  fgets(in_line, MAXLINE, f);
  sscanf(in_line, "%d %d %d", &n, &na, &nnz);
  InitMatrix();

  int x, source, target, action, last_index = -1, index = -1, branch = 0, last_source = -1, last_action = -1;
  double prob;
  long no_lines = 0;
  
  while (fgets(in_line, MAXLINE, f) != NULL)
  {
    switch (in_line[0])
    {
    case 'c':                  /* skip lines with comments */
    case '\n':                 /* skip empty lines   */
    case '\0':                 /* skip empty lines at the end of file */
      break;
    default:                    /* transition description */
      source = strtol(&in_line[0], &next, 10);
      action = strtol(next, &next, 10);
      if (*next == '/') prism_format = false;
      if (prism_format && (source != last_source || action != last_action)) ++index;
      else if (!prism_format) index = strtol(next + 1, &next, 10);
      last_source = source;
      last_action = action;
      target = strtol(next, &next, 10);
      prob = strtod(next, 0);
      
      if (last_index != -1 && last_index != index) ++branch, last_index = index, state_starts[source + 1]++;
      else if (last_index == -1) last_index = index, state_starts[source + 1]++;
      
      row_starts[branch + 1]++;
      
      if (seen_actions.find(action) == seen_actions.end())
      {
        actions[branch] = act_table.size();
        seen_actions.insert(action);
        act_table.push_back(action);
      }
      else
      {
        for (at = act_table.begin(), x = 0; at != act_table.end(); ++at, ++x)
        {
          if (*at == action)
          {
            actions[branch] = x;
            break;
          }
        }
      }
      
      cols[no_lines] = target;
      non_zeros[no_lines] = prob;
      no_lines++;
      
      break;
    } /* end of switch */
  }   /* end of input loop */

  if (no_lines != nnz)
  {
    printf("There are more transitions than described in the first line, please check!\n");
    exit(1);
  }
  
  // now construct the row_starts
  for(int i=1; i<na+1; i++)
  {
    row_starts[i] += row_starts[i-1];
  }
  for(int i=1; i<n+1; i++)
  {
    state_starts[i] += state_starts[i-1];
  }
  
  // set up the action table
  atable = new int[act_table.size()];
  for (at = act_table.begin(), index = 0; at != act_table.end(); ++at, ++index) atable[index] = *at;
  da = act_table.size();
}

// Determine a good floating-point approximation threshold for this model. Returns -1.0
// if the model is not stochastic.
double MarkovChain::PrecisionThreshold()
{
  double dist, min_trans, max_delta, sum;
  int m, i;

  min_trans = 1.0;
  max_delta = 0.0;

  for (m = 0; m < n; ++m)
  {
    for (i = row_starts[m], sum = 0.0; i < row_starts[m + 1]; ++i) sum += non_zeros[i];
    dist = fabs(1.0 - sum);
    if (!ContinuousTimeModel()) sum = 1.0;
    for (i = row_starts[m]; i < row_starts[m + 1]; ++i)
    {
      if (non_zeros[i] < min_trans * sum) min_trans = non_zeros[i] / sum;
    }
    if (dist > max_delta) max_delta = dist;
  }

  if (ContinuousTimeModel()) return min_trans * 0.1;
  else
  {
    if (max_delta > min_trans) return -1.0;
    if (max_delta == 0.0) return min_trans * 0.1;
    return 10 / ((9 / min_trans) + (1 / max_delta));
  }

  return 0.0;
}

// Determine a good floating-point approximation threshold for this model. Returns -1.0
// if the model is not stochastic.
double ProbabilisticAutomaton::PrecisionThreshold()
{
  double dist, min_trans, max_delta, sum;
  int m, i, l;

  min_trans = 1.0;
  max_delta = 0.0;

  for (m = 0; m < n; ++m)
  {
    for (i = state_starts[m]; i < state_starts[m+1]; ++i)
    {
      for (l = row_starts[i], sum = 0.0; l < row_starts[i+1]; ++l) sum += non_zeros[l];
      dist = fabs(1.0 - sum);
      if (!ContinuousTimeModel()) sum = 1.0;
      for (l = row_starts[i]; l < row_starts[i+1]; ++l)
      {
        if (non_zeros[l] < min_trans * sum) min_trans = non_zeros[l] / sum;
      }
      if (dist > max_delta) max_delta = dist;
    }
  }
  if (ContinuousTimeModel()) return min_trans * 0.1;
  else
  {
    if (max_delta > min_trans) return -1.0;
    if (max_delta == 0.0) return 0.0;
    return 10 / ((9 / min_trans) + (1 / max_delta));
  }

  return 0.0;
}

// Dump the model to a file.
void MarkovChain::Write(FILE *out)
{
  int i, j;
  
  if (!out || !row_starts || !cols || !non_zeros) return;
  
  fprintf(out, "%d %ld\n", n, nnz);
  
  for (i = 0; i < n; ++i)
  {
    for (j = row_starts[i]; j < row_starts[i+1]; ++j)
    {
      fprintf(out, "%d %d %#.20e\n", i, cols[j], non_zeros[j]);
    }
  }
}

// Dump the model to a file.
void ProbabilisticAutomaton::Write(FILE *out)
{
  int i, j, k, a0, a1, k0, k1, idx;
  
  fprintf(out, "%d %d %d\n", n, na, nnz);
  
  idx = 0;
  
  for (i = 0; i < n; ++i)
  {
    a0 = state_starts[i  ];
    a1 = state_starts[i+1];
    
    for (j = a0; j < a1; ++j)
    {
      k0 = row_starts[j  ];
      k1 = row_starts[j+1];
      
      for (k = k0; k < k1; ++k)
      {
        fprintf(out, "%d %d/%d %d %.10f\n", i, actions[j], idx, cols[k], non_zeros[k]);
      }
      
      if (k0 != k1) ++idx;
    }
  }
}
