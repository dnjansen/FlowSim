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



#include "Weak.h"
#include "compactmaxflow.cc"
#include "Weak_auxiliary.cc"

using namespace std;

// Compute the simulation relation
unsigned int WeakSimulation_MC::Simulate(ProbabilisticModel *model, std::set<std::pair<int,int> > *result)
{
  int r, iterations = 0;
  Pair *new_relation;
  FILE *reldump;
  
  if (!model || model->Type() != ProbabilisticModel::MC) return 0;

  if (model->ContinuousTimeModel())
  {
    fprintf(stderr, "Continuous time MCs are not supported by this implementation\n");
    return 0;
  }
  
  // Copy sparse matrix structures
  m = (MarkovChain*)model;
  non_zeros = new double[m->nnz];
  memcpy(non_zeros, m->non_zeros, sizeof(double) * m->nnz);
  cols = new int[m->nnz];
  memcpy(cols, m->cols, sizeof(int) * m->nnz);
  
  row_starts = m->row_starts;
  n_states = m->n;
  
  rmap.Create(n_states);
  
  size_of_relation = 0;
  
#ifdef DEBUG
  rel_space = 0;
  memset(&stats, 0, sizeof(stats));
  //stats.mem_model = (sizeof(double) * m->nnz) + (sizeof(int) * m->nnz);
#endif//DEBUG

  InitializeRelation();

#ifdef DEBUG
  rel_space = sizeof(Pair) * size_of_relation;
  stats.mem_relation = rel_space;
  stats.num_initial_pairs = size_of_relation;
  stats.mem_relation_map = rmap.MemoryUsage();
#endif//DEBUG

  hold_hard_pairs = true;

  reldump = fopen("reldump", "wb");
  simlog = fopen("simlog", "wb");
  fprintf(reldump, "==== Initial relation ====\n");
  for (int n = 0; n < size_of_relation; ++n) fprintf(reldump, "%4d %4d\n", relation[n].x, relation[n].y);
  fprintf(reldump, "\n");
  
  while (1)
  {
    ++iterations;
    fprintf(simlog, "================ Iteration %d ================\n", iterations);
    r = IterateRelation(&new_relation);
    
    if (r == -1) break;
    
    fprintf(reldump, "==== After iteration %d ====\n", iterations);
    for (int n = 0; n < r; ++n) fprintf(reldump, "%4d %4d\n", new_relation[n].x, new_relation[n].y);
    if (r == 0) fprintf(reldump, "(empty)\n");
    fprintf(reldump, "\n");
    fprintf(simlog, "\n");
    
    delete [] relation;
    relation = new_relation;
    size_of_relation = r;
    
    if (r == 0) break;
  }
  
#ifdef DEBUG
  stats.num_iterations = iterations;
  stats.num_final_pairs = size_of_relation;
  stats.mem_maxflow = CompactMaxFlow<double>::global_space_peak;
  stats.num_maxflow = CompactMaxFlow<double>::global_times_invoked;
  stats.num_p_invariant_fails = CompactMaxFlow<double>::global_p_inv_fails;
  stats.num_sig_arc_fails = CompactMaxFlow<double>::global_sig_arc_fails;
  stats.min_complexity = CompactMaxFlow<double>::min_complexity;
  stats.max_complexity = CompactMaxFlow<double>::max_complexity;
  CompactMaxFlow<double>::ResetStats();
#endif//DEBUG

  delete [] non_zeros;
  delete [] cols;

  // Return result in the provided set if requested
  if (result)
  {
    result->clear();
    for (int n = 0; n < size_of_relation; ++n) result->insert(std::pair<int,int>(relation[n].x, relation[n].y));
  }
  fprintf(reldump, "==== Final relation ====\n");
  for (int n = 0; n < size_of_relation; ++n) fprintf(reldump, "%4d %4d\n", relation[n].x, relation[n].y);
  if (size_of_relation == 0) fprintf(reldump, "(empty)\n");
  
  // Free relation
  if (relation)
  {
    for (int n = 0; n < size_of_relation; ++n) relation[n].Reset();
    delete [] relation;
  }
  
  return size_of_relation;
}

// Build the initial map of the relation as a 2D array of booleans.
// Pairs are calculated based on labels.
int WeakSimulation_MC::BuildRelationMap()
{
  int m, n, size = 0;
  
  for (m = 0; m < n_states - 1; ++m)
  {
    for (n = m + 1; n < n_states; ++n)
    {
      if (Label(n) == Label(m))
      {
        size += 2;
        rmap.Set(m, n);
        rmap.Set(n, m);
      }
    }
  }
  
  rmap.Commit();
  
  return size;
}

// Create the initial list of pairs based on the relation map computed by
// BuildRelationMap(). The pairs are inserted in a certain order which allows
// a constant time test to see if the reverse pair has been tested already or not.
void WeakSimulation_MC::InitializeRelation()
{
  int n, m, s;
  
  size_of_relation = BuildRelationMap();
  
  // Allocate array of pairs now that we know how many pairs there are total
  relation = new Pair[size_of_relation];
  
  // Iterate through all (n,m) and add the pair if it is in the map
  for (n = 0, s = 0; n < n_states; ++n)
  {
    for (m = 0; m < n_states; ++m)
    {
      if (n == m) continue;
      if (rmap(n, m))
      {
        (relation + s)->x = n;
        (relation + s)->y = m;
        ++s;
      }
    }
  }
}

// Naively iterate the relation by repeatedly invoking WeakSimulation
// on every pair.
int WeakSimulation_MC::IterateRelation(Pair **new_relation)
{
  std::vector<int> keep; // List of pairs to remain in the relation
  int i, new_size;
  bool hard, keep_going = false;
  
  num_easy_pairs = 0;

  for (i = 0; i < size_of_relation; ++i)
  {
    if (DecideWeakSimulation(relation + i, hard)) keep.push_back(i);
    else
    {
      if (!hard) ++num_easy_pairs;
      rmap.Clear(relation[i].x, relation[i].y);
    }
  }
  
  if (num_easy_pairs == 0 && hold_hard_pairs)
  {
    keep_going = true;
    hold_hard_pairs = false;
  }
  
  // Copy the new relation map back into the main buffer
  rmap.Commit();
  
  // No pairs were deleted so the current relation is the correct result
  if (keep.size() == (unsigned int)size_of_relation && !keep_going)
  {
    *new_relation = 0;
    return -1;
  }
  
  // Recreate the relation array without the pairs that were dropped
  new_size = keep.size();
  *new_relation = new Pair[new_size];
  for (i = 0; i < new_size; ++i) (*new_relation)[i] = relation[keep[i]];
  
  return new_size;
}

void print_partition(int h, int *pa, int *pi, FILE *simlog)
{
  int n, m;
  for (n = 0; n < h; ++n)
  {
    fprintf(simlog, "[");
    for (m = pi[n]; m < pi[n + 1]; ++m) fprintf(simlog, (m == pi[n] ? "%d" : " %d"), pa[m]);
    fprintf(simlog, "] ");
  }
  fprintf(simlog, "\n");
}

// Compute strong simulation on a particular pair in the relation
bool WeakSimulation_MC::DecideWeakSimulation(Pair *p, bool &hard)
{
  int n, m, l, h;
  int *par_array, *par_index;
  double Ps1Ai, Ps2Ai, *gamma;
  bool res;
  
  set<pair<pair<int,int>,double> > arcs;
  set<double> breakpoints;
  set<double>::iterator bp;
  
  hard = false;
  
  if (p->x == p->y) return true;
  if (row_starts[p->x + 1] - row_starts[p->x] == 0) return true;
  if (row_starts[p->y + 1] - row_starts[p->y] == 0) return false;
  
  // WS: Line 1
  //printf("Line 1\n");
  for (n = row_starts[p->x]; n < row_starts[p->x + 1]; ++n)
  {
    if (!rmap(cols[n], p->y) && cols[n] != p->y) break;
  }
  if (n == row_starts[p->x + 1])
  {
    fprintf(simlog, "KEEP (%d,%d) because of WS:1\n", p->x, p->y);
    return true;
  }
  
  // WS: Line 2
  //printf("Line 2\n");
  for (n = row_starts[p->y]; n < row_starts[p->y + 1]; ++n)
  {
    if (!rmap(p->x, cols[n]) && p->x != cols[n]) break;
  }
  if (n == row_starts[p->y + 1])
  {
    // WS: Lines 3-6
    //printf("Lines 3-6\n");
    if (CanReach(p->x, p->y, false))
    {
      fprintf(simlog, "KEEP (%d,%d) because of WS:3\n", p->x, p->y);
      return true;
    }
    for (n = row_starts[p->x]; n < row_starts[p->x + 1]; ++n)
    {
      if (rmap(cols[n], p->y)) continue;
      if (!CanReach(cols[n], p->y, true)) break;
    }
    if (n == row_starts[p->x + 1])
    {
      fprintf(simlog, "KEEP (%d,%d) because of WS:5\n", p->x, p->y);
      return true;
    }
    else
    {
      fprintf(simlog, "DROP (%d,%d) because of WS:6\n", p->x, p->y);
      return false;
    }
  }
  
  // Ws_DTMC: Line 1
  //printf("DTMC Line 1\n");
  h = ConstructSuccessorPartition(p->x, p->y, &par_array, &par_index);
  hard = (h == 1);
  /*for (n = row_starts[p->x]; n < row_starts[p->x + 1]; ++n)
  {
    for (m = row_starts[p->y]; m < row_starts[p->y + 1]; ++m)
    {
      if (rmap(cols[n], cols[m]) || n == m) printf("%d <- %d\n", cols[n], cols[m]);
    }
  }*/
  fprintf(simlog, "PART (%d,%d) %d %s: ", p->x, p->y, h, (h == 1 ? "block" : "blocks"));
  print_partition(h, par_array, par_index, simlog);
  
  // Ws_DTMC: Line 2
  //printf("DTMC Line 2\n");
  if (h == 1)
  {
    delete [] par_array;
    delete [] par_index;
    
    if (hold_hard_pairs)
    {
      fprintf(simlog, "HOLD (%d,%d) because h=1\n", p->x, p->y);
      return true;
    }
    //fprintf(stderr, "FPS\n");
    // Set up maxflow network
    //NOTE: This part is temporary and will be rewritten
    for (n = row_starts[p->x]; n < row_starts[p->x + 1]; ++n)
    {
      arcs.insert(make_pair(make_pair(1, cols[n] + 2), non_zeros[n]));
    }
    for (n = row_starts[p->y]; n < row_starts[p->y + 1]; ++n)
    {
      arcs.insert(make_pair(make_pair(cols[n] + n_states + 2, 2 * n_states + 3), non_zeros[n]));
    }
    for (n = row_starts[p->y]; n < row_starts[p->y + 1]; ++n)
    {
      for (m = row_starts[p->x]; m < row_starts[p->x + 1]; ++m)
      {
        if (rmap(cols[m], cols[n]) || cols[m] == cols[n]) arcs.insert(make_pair(make_pair(cols[m] + 2, cols[n] + n_states + 2), 100));
      }
    }
    
    // Ws_FPS: Line 1
    FindBreakpoints(arcs, 1, 2 * n_states + 3, breakpoints);
    
    fprintf(simlog, "[--- (%d,%d) computing breakpoints for network:\n", p->x, p->y);
    AnalyzeNetwork(arcs, 1, 2 * n_states + 3, 1.0, 1.0, 0, 0, 0, simlog, "----");
    
    fprintf(simlog, "---] Breakpoints:");
    for (bp = breakpoints.begin(); bp != breakpoints.end(); ++bp) fprintf(simlog, " %f", *bp);
    fprintf(simlog, "\n");
    
    // Ws_FPS: Line 2
    if (ValidFlow(p->x, p->y, breakpoints))
    {
      fprintf(simlog, "KEEP (%d,%d) because of WS_FPS:2\n", p->x, p->y);
      return true;
    }
    
    // Ws_FPS: Line 3
    fprintf(simlog, "DROP (%d,%d) because of WS_FPS:3\n", p->x, p->y);
    return false;
  }
  
  gamma = new double[h];
  
  // Ws_DTMC: Line 3
  //printf("DTMC Line 3\n");
  //FIXME: This operation needs to be optimized; theoretical complexity
  //       of this code is horrible.
  for (n = 1; n < h; ++n)
  {
    Ps1Ai = 0.0;
    Ps2Ai = 0.0;
    for (m = par_index[n]; m < par_index[n + 1]; ++m)
    {
      for (l = row_starts[p->x]; l < row_starts[p->x + 1]; ++l)
      {
        if (par_array[m] == cols[l])
        {
          Ps1Ai += non_zeros[l];
          break;
        }
      }
      for (l = row_starts[p->y]; l < row_starts[p->y + 1]; ++l)
      {
        if (par_array[m] == cols[l])
        {
          Ps2Ai += non_zeros[l];
          break;
        }
      }
    }
    
    //printf("Ps1A%d = %f   ", n + 1, Ps1Ai);
    //printf("Ps2A%d = %f   ", n + 1, Ps2Ai);
    //printf("Q = %f\n", Ps1Ai / Ps2Ai);
    
    // Ws_DTMC: Line 4
    assert(Ps1Ai != 0.0 || Ps2Ai != 0.0);
    
    // Ws_DTMC: Line 5
    //printf("DTMC Line 5\n");
    if (Ps1Ai == 0.0 || Ps2Ai == 0.0)
    {
      delete [] par_array;
      delete [] par_index;
      delete [] gamma;
      
      fprintf(simlog, "DROP (%d,%d) because of WS_DTMC:5\n", p->x, p->y);
      return false;
    }
    
    // Ws_DTMC: Line 6
    //printf("DTMC Line 6\n");
    gamma[n] = Ps1Ai / Ps2Ai;
    
    // Ws_DTMC: Line 7
    if (n > 1 && gamma[n - 1] != gamma[n])
    {
      fprintf(simlog, "DROP (%d,%d) because of WS_DTMC:7\n", p->x, p->y);
      return false;
    }
  }
  
  // Ws_DTMC: Line 8
  //printf("DTMC Line 8\n");
  breakpoints.clear();
  breakpoints.insert(gamma[h - 1]);
  res = ValidFlow(p->x, p->y, breakpoints, h, par_array, par_index);
  
  delete [] par_array;
  delete [] par_index;
  delete [] gamma;

  fprintf(simlog, "%s (%d,%d) because of WS_DTMC:8\n", (res ? "KEEP" : "DROP"), p->x, p->y);
  return res;
}

#if 0
// Sort the successors of each state by probability (1st) and label (2nd)
void WeakSimulation_MC::SortSuccessors()
{
  int i, suc;
  int *successor_order = new int[m->nnz + 1];
  int *new_cols = new int[m->nnz];
  double *new_nz = new double[m->nnz];
  SuccessorOrder cmp(this);
  
  for (i = 0; i < m->nnz; ++i) successor_order[i] = i;
  
  // Sort successor set for each state
  for (i = 0; i < n_states; ++i)
  {
    suc = row_starts[i+1] - row_starts[i];
    if (suc <= 1) continue;
    std::sort(successor_order + row_starts[i], successor_order + row_starts[i+1], cmp);
  }
  
  // Create new cols and non_zero arrays with the updated order
  for (i = 0; i < m->nnz; ++i)
  {
    new_cols[i] = cols[successor_order[i]];
    new_nz[i] = non_zeros[successor_order[i]];
  }
  
  delete [] successor_order;
  
  // Delete old arrays and copy new ones over
  delete [] cols;
  delete [] non_zeros;
  
  cols = new_cols;
  non_zeros = new_nz;
}
#endif

// Compute the partitions A_1, ..., A_h of the localized relation H
int WeakSimulation_MC::ConstructSuccessorPartition(int s1, int s2, int **par_array, int **par_array_index)
{
  int i, n, s, last, total_states = 0, number_partitions = 0;
  bool finished = false;
  int *explored = new int[n_states]; //indicates whether the state (in state_s) is already explored.
  int *suc_s1 = new int[n_states], *suc_s2 = new int[n_states];
  
  vector<set<int> > part;
  vector<set<int> >::iterator part_i;
  set<int> states;
  set<int>::iterator states_i;
  
  queue<int> qstates; // used for the BFS algoriths.
  int processing = 0; // count the number of processing states. Increased after qstates.push().
  
  for (i = 0; i < n_states; ++i)
  {//set the values back to 0.
      suc_s1[i] = 0;
      suc_s2[i] = 0;
      explored[i] = 0;
  }
  
  for (i = row_starts[s1]; i < row_starts[s1 + 1]; i++)
  { //insert V1 into the qstates
    suc_s1[cols[i]] = 1;
    if (rmap(cols[i], s2))
    {
      qstates.push(cols[i]);
      explored[cols[i]] = 1;
      ++processing;
    }
  }
  
  for (i = row_starts[s2]; i < row_starts[s2 + 1]; i++)
  { //insert V2, if not already inserted, into the qstates
    suc_s2[cols[i]] = 1;
    if (rmap(s1, cols[i]) && (!explored[cols[i]]))
    {
      qstates.push(cols[i]);
      explored[cols[i]] = 1;
      ++processing;
    }
  }
  
  //if qstates is empty, the set An is empty. The special partition is saved at the beginning.
  //In the paper is the last partition.
  if (qstates.empty())
  {
    states.clear();
    part.push_back(states);
  }

  while (!finished)
  {
    if (!qstates.empty()) number_partitions++;
    while (!qstates.empty())
    {
      last = qstates.front();
      qstates.pop();
      states.insert(last);
      if (suc_s1[last])
      {
        //for (i = relation_forward_index[last]; i < relation_forward_index[last + 1]; ++i)
        for (i = 0; i < n_states; ++i)
        {
          //s = relation_forward[i];
          if (!rmap(last, i)) continue;
          s = i;
          if (suc_s2[s] && (!explored[s]))
          {
            qstates.push(s);
            explored[s] = 1;
            ++processing;
          }
        }
      }
      if (suc_s2[last])
      {
        //for (i = relation_backward_index[last]; i < relation_backward_index[last + 1]; ++i)
        for (i = 0; i < n_states; ++i)
        {
          //s = relation_backward[i];
          if (!rmap(i, last)) continue;
          s = i;
          if (suc_s1[s] && (!explored[s]))
          {
            qstates.push(s);
            explored[s] = 1;
            ++processing;
          }
        }
      }
    }//end of the first loop

    if (states.size() > 0) part.push_back(states);
    total_states += states.size();
    states.clear();

    for (i = row_starts[s1]; i < row_starts[s1 + 1]; ++i)
    { //insert state into the qstates
      if (!(explored[cols[i]]))
      {
        qstates.push(cols[i]);
        explored[cols[i]] = 1;
        ++processing;
        break;
      }
    }
    
    if (qstates.empty())
    {
      for (i = row_starts[s2]; i < row_starts[s2 + 1]; ++i)
      { //insert state into the qstates
        if (!(explored[cols[i]]))
        {
          qstates.push(cols[i]);
          explored[cols[i]] = 1;
          ++processing;
          break;
        }
      }
    }
    
    if (qstates.empty()) finished = true;
  }//end of while

  number_partitions = part.size();
  *par_array = new int[total_states];
  *par_array_index = new int[number_partitions + 1];
  
  (*par_array_index)[0] = 0;
  for (part_i = part.begin(), i = 1, n = 0; part_i != part.end(); ++part_i, ++i)
  {
    (*par_array_index)[i] = (*par_array_index)[i - 1] + part_i->size();
    for (states_i = part_i->begin(); states_i != part_i->end(); ++states_i)
    {
      (*par_array)[n++] = *states_i;
    }
  }
  
  delete [] explored;
  delete [] suc_s1;
  delete [] suc_s2;
  
  return number_partitions;
}

void print_set(set<int> &s, const char *l, FILE *simlog)
{
  set<int>::iterator i;
  fprintf(simlog, "%s:", l);
  for (i = s.begin(); i != s.end(); ++i) printf(" %d", *i);
  printf("\n");
}

void print_setd(set<double> &s, const char *l, FILE *simlog)
{
  set<double>::iterator i;
  fprintf(simlog, "%s:", l);
  for (i = s.begin(); i != s.end(); ++i) printf(" %f", *i);
  printf("\n");
}

// Test if at least one of the parameters in the given set
// is valid for the network given by the pair of states
// and the current relation.
bool WeakSimulation_MC::ValidFlow(int s1, int s2, set<double> param, int h, int *par, int *par_idx)
{
  int n, m, k;
  int s1_suc = row_starts[s1 + 1] - row_starts[s1]; // Get number of successor states
  int s2_suc = row_starts[s2 + 1] - row_starts[s2]; //
  CompactFeasibleFlow sim;
  set<int> MU1, MU2;
  set<int>::iterator mui;
  set<double> paramkeep;
  set<double>::iterator pi;
  
  assert(h > 0);
  
  if (h == 1)
  {
    for (m = row_starts[s1]; m < row_starts[s1 + 1]; ++m)
    {
      if (!rmap(cols[m], s2)) { MU1.insert(cols[m]); }
    }
    for (m = row_starts[s2]; m < row_starts[s2 + 1]; ++m)
    {
      if (!rmap(s1, cols[m])) { MU2.insert(cols[m]); }
    }
    if (MU1.size() != 0 && MU2.size() != 0)
    {
      fprintf(simlog, "PART (%d,%d) MU1 [", s1, s2);
      for (mui = MU1.begin(), m = 0; mui != MU1.end(); ++mui) fprintf(simlog, "%s%d", (m++ == 0 ? "" : " "), *mui);
      fprintf(simlog, "] MU2 [");
      for (mui = MU2.begin(), m = 0; mui != MU2.end(); ++mui) fprintf(simlog, "%s%d", (m++ == 0 ? "" : " "), *mui);
      fprintf(simlog, "]\n");
      sim.CreateNetwork(cols, non_zeros, &rmap, row_starts[s1], s1_suc, row_starts[s2], s2_suc, MU1, MU2, 0);
      if (sim.IsFlowFeasible(param, 0)) return true;
    }
  }
  else
  {
    MU1.clear();
    MU2.clear();
    for (m = row_starts[s1]; m < row_starts[s1 + 1]; ++m)
    {
      for (k = 0; k < par_idx[1]; ++k)
      {
        if (cols[m] == par[k] && !rmap(cols[m], s2)) { MU1.insert(cols[m]); break; }
      }
    }
    for (m = row_starts[s2]; m < row_starts[s2 + 1]; ++m)
    {
      for (k = 0; k < par_idx[1]; ++k)
      {
        if (cols[m] == par[k] && !rmap(s1, cols[m])) { MU2.insert(cols[m]); break; }
      }
    }
    //print_set(MU1, "MU1_1", simlog);
    //print_set(MU2, "MU2_1", simlog);
    if (MU1.size() != 0 && MU2.size() != 0)
    {
      sim.CreateNetwork(cols, non_zeros, &rmap, row_starts[s1], s1_suc, row_starts[s2], s2_suc, MU1, MU2, 0);
      
      for (pi = param.begin(); pi != param.end(); ++pi)
      {
        if (sim.IsFlowFeasible(*pi)) paramkeep.insert(*pi);
      }
    }
    else paramkeep = param;

    for (n = 1; n < h && paramkeep.size() > 0; ++n)
    {
      param = paramkeep;
      paramkeep.clear();
    
      MU1.clear();
      MU2.clear();
      for (m = row_starts[s1]; m < row_starts[s1 + 1]; ++m)
      {
        for (k = par_idx[n]; k < par_idx[n + 1]; ++k)
        {
          if (cols[m] == par[k]) { MU1.insert(cols[m]); break; }
        }
      }
      for (m = row_starts[s2]; m < row_starts[s2 + 1]; ++m)
      {
        for (k = par_idx[n]; k < par_idx[n + 1]; ++k)
        {
          if (cols[m] == par[k]) { MU2.insert(cols[m]); break; }
        }
      }
      //print_set(MU1, "MU1_n", simlog);
      //print_set(MU2, "MU2_n", simlog);
      if (MU1.size() == 0 || MU2.size() == 0) continue;
      sim.CreateNetwork(cols, non_zeros, &rmap, row_starts[s1], s1_suc, row_starts[s2], s2_suc, MU1, MU2, 0);
      
      for (pi = param.begin(); pi != param.end(); ++pi)
      {
        if (sim.IsFlowFeasible(*pi)) paramkeep.insert(*pi);
      }
    }
  }
  
  return paramkeep.size() > 0;
}

// Check if there exists a state s, reachable from s2, which simulates s1 (inrel=true)
// or which does not simulate s1 (inrel=false).
bool WeakSimulation_MC::CanReach(int s1, int s2, bool inrel)
{
  queue<int> reachables;
  bool *visited;
  int s, n;
  
  visited = new bool[n_states];
  for (n = 0; n < n_states; ++n) visited[n] = false;
  
  for (n = row_starts[s2]; n < row_starts[s2 + 1]; ++n)
  {
    reachables.push(cols[n]);
    visited[cols[n]] = true;
  }
  
  while (!reachables.empty())
  {
    s = reachables.front();
    
    if (rmap(s1, s) == inrel) break;
    else for (n = row_starts[s]; n < row_starts[s + 1]; ++n)
    {
      if (!visited[cols[n]])
      {
        reachables.push(cols[n]);
        visited[cols[n]] = true;
      }
    }
    
    reachables.pop();
  }
  
  delete [] visited;
  
  return (!reachables.empty());
}
