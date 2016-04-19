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



#include "Strong.h"

// Compute the simulation relation
unsigned int StrongSimulation_MC::Simulate(ProbabilisticModel *model, std::set<std::pair<int,int> > *result)
{
  int r, iterations = 0;
  Pair *pp;
  
  if (!model || model->Type() != ProbabilisticModel::MC) return 0;
  
  // Copy sparse matrix structures
  m = (MarkovChain*)model;
  non_zeros = new double[m->nnz];
  memcpy(non_zeros, m->non_zeros, sizeof(double) * m->nnz);
  cols = new int[m->nnz];
  memcpy(cols, m->cols, sizeof(int) * m->nnz);
  
  row_starts = m->row_starts;
  n_states = m->n;
  
  rmap.Create(1, n_states);
  
  // Build initial relation
  relation = 0;
  if (m->ContinuousTimeModel()) size_of_relation = BuildRelationMap_CTMC();
  else size_of_relation = BuildRelationMap_DTMC();
  
#ifdef DEBUG
  memset(&stats, 0, sizeof(stats));
  CompactMaxFlow<double>::ResetStats();
  stats.mem_model = (sizeof(double) * m->nnz) + (sizeof(int) * m->nnz);
  stats.num_initial_pairs = size_of_relation;
#endif//DEBUG
  
#if defined(OPT_PARTITION)
  // Preprocessing for state partitioning optimization
  SortSuccessors();
  MakeFirstPartition();
#endif

#ifdef DEBUG
#if defined(OPT_PARTITION)
  {
    stats.mem_partition_map = sizeof(char) * partitions * partitions;
    stats.num_partitions = partitions;
  }
#else
  stats.mem_partition_map = 0;
  stats.num_partitions = 0;
#endif
#endif//DEBUG

  // Main refinement loop
  while (1)
  {
    ++iterations;
    
#if defined(OPT_PARTITION)
    if (iterations == 1) r = IterateRelation_FirstPartition();
    else r = IterateRelation(false);
#else
    r = IterateRelation(iterations == 1);
#endif
    
    if (r == -1) break;
    
    size_of_relation = r;
    
    if (r == 0) break;
  }
  
#ifdef DEBUG
  stats.mem_relation_map = rmap.MemoryUsage();
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

  //if (optflags & (OPT_PARTITION | OPT_PARTITION2))
#if defined(OPT_PARTITION)
  {
    delete [] partition;
  }
#endif

  delete [] non_zeros;
  delete [] cols;

  // Return result in the provided set if requested
  if (result)
  {
    result->clear();
    pp = relation;
    while (pp)
    {
      result->insert(std::make_pair(pp->x, pp->y));
      pp = pp->next;
    }
  }
  
  // Free relation
  while (relation)
  {
    pp = relation->next;
    delete relation;
    relation = pp;
  }
  
  return size_of_relation;
}

#ifdef WITH_VERIFIER
// Verify that a set of pairs is a simulation relation for the given model
bool StrongSimulation_MC::Verify(ProbabilisticModel *model, std::set<std::pair<int,int> > &hypothesis,
        std::set<std::pair<int,int> > *false_positives, std::set<std::pair<int,int> > *false_negatives)
{
  bool res, global_res = true, ctmc;
  std::set<std::pair<int,int> >::iterator hi;
  Pair p;
  int i, j;
  double *dist_sums;
  
  if (!model || model->Type() != ProbabilisticModel::MC) return false;
  
  // Copy sparse matrix structures
  m = (MarkovChain*)model;
  non_zeros = m->non_zeros;
  cols = m->cols;
  row_starts = m->row_starts;
  n_states = m->n;
  
  // If we are dealing with a CTMC, we have to normalize transitions first
  if (ctmc = m->ContinuousTimeModel())
  {
    non_zeros = new double[m->Transitions()];
    dist_sums = new double[m->States()];
    for (i = 0; i < n_states; ++i)
    {
      dist_sums[i] = 0.0;
      for (j = row_starts[i]; j < row_starts[i + 1]; ++j) dist_sums[i] += m->non_zeros[j];
      for (j = row_starts[i]; j < row_starts[i + 1]; ++j) non_zeros[j] = m->non_zeros[j] / dist_sums[i];
    }
  }
  
  // Load the hypothesis into the relation map
  rmap.Create(n_states);
  
  size_of_relation = hypothesis.size();
  
  for (hi = hypothesis.begin(); hi != hypothesis.end(); ++hi)
  {
    rmap.Set(hi->first, hi->second);
  }
  
  rmap.Commit();
  
  // Decide simulation for all non-identity pairs, based on the relation
  // being tested, and match the respective results against the hypothesis.
  p.next = 0;
  for (p.x = 0; p.x < n_states; ++p.x)
  {
    for (p.y = 0; p.y < n_states; ++p.y)
    {
      if (p.x == p.y) continue;
#ifdef OPT_CACHE_NETS
      p.simulation = 0;
#endif//OPT_CACHE_NETS
      res = (Label(p.x) == Label(p.y)
             && (ctmc ? CompactMaxFlow<double>::_Tleq(dist_sums[p.x], dist_sums[p.y]) : true)
             && DecideStrongSimulation(&p));
      if (res != rmap(p.x, p.y))
      {
        global_res = false;
        if (res && !rmap(p.x, p.y) && false_negatives) false_negatives->insert(std::make_pair(p.x, p.y));
        else if (!res && rmap(p.x, p.y) && false_positives) false_positives->insert(std::make_pair(p.x, p.y));
      }
#ifdef OPT_CACHE_NETS
      if (res) {
        delete p.simulation;
        p.simulation = NULL;
      }
#endif//OPT_CACHE_NETS
    }
  }
  
  if (ctmc)
  {
    delete [] non_zeros;
    delete [] dist_sums;
  }
  
  return global_res;
}
#endif//WITH_VERIFIER

// Build relation map for discrete time MCs
int StrongSimulation_MC::BuildRelationMap_DTMC()
{
  int m, n, size = 0;
  
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
  
  rmap.Commit();
  
  return size;
}

// Build relation map for continuous time MCs
int StrongSimulation_MC::BuildRelationMap_CTMC()
{
  int m, n, size = 0;
  double *psums = 0;
  
  // Compute R(s, S) for all s and normalize transition rates
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

// Naively iterate the relation by repeatedly invoking StrongSimulation
// on every pair.
int StrongSimulation_MC::IterateRelation(bool first)
{
  int s1, s2, new_size = size_of_relation;
  Pair p, *pp, **anchor = &relation;
  
  // If we are in the first iteration, add those pairs to the relation
  // that simulate and satisfy the initial condition. Otherwise, delete
  // the pairs that don't simulate.
  if (first)
  {
    for (s1 = 0; s1 < n_states; ++s1)
    {
      for (s2 = 0; s2 < n_states; ++s2)
      {
        if (s1 == s2 || !rmap(s1, s2)) continue;
        p.x = s1;
        p.y = s2;
#ifdef OPT_CACHE_NETS
        p.simulation = 0;
#endif//OPT_CACHE_NETS
        if (DecideStrongSimulation(&p))
        {
          pp = new Pair;
          pp->x = s1;
          pp->y = s2;
#ifdef OPT_CACHE_NETS
          pp->simulation = p.simulation;
          p.simulation = NULL;
#endif//OPT_CACHE_NETS
          pp->next = relation;
          relation = pp;
        }
        else
        {
          --new_size;
          rmap.Clear(s1, s2);
        }
      }
    }
#ifdef DEBUG
    stats.mem_relation = sizeof(Pair) * new_size;
#endif//DEBUG

    rmap.Commit();
    return new_size;
  }
  else
  {
    pp = relation;
    while (pp)
    {
      if (!DecideStrongSimulation(pp))
      {
        rmap.Clear(pp->x, pp->y);
        *anchor = pp->next;
#       ifdef OPT_CACHE_NETS
          assert(NULL == pp->simulation);
#       endif
        delete pp;
        pp = *anchor;
        --new_size;
      }
      else
      {
        anchor = &pp->next;
        pp = pp->next;
      }
    }
  }
  
  // Copy the new relation map back into the main buffer
  rmap.Commit();
  
  // No pairs were deleted so the current relation is the correct result
  if (new_size == size_of_relation) return -1;
  
  return new_size;
}

// Compute the first iteration exploiting state partitioning. This must be used
// only for the first iteration and will be incorrect otherwise. MakeFirstPartition()
// must be invoked before this function.
int StrongSimulation_MC::IterateRelation_FirstPartition()
{
  int s1, s2, new_size = size_of_relation;
  char *partition_map, c; // 2D array mapping pairs of partitions to values of unknown, failure or success
  Pair pair, *pp;

  // Allocate and initialize the partition map
  partition_map = new char[partitions * partitions];
  memset(partition_map, 0, sizeof(char) * partitions * partitions);
  
  for (s1 = 0; s1 < n_states; ++s1)
  {
    for (s2 = 0; s2 < n_states; ++s2)
    {
      if (s1 == s2 || !rmap(s1, s2)) continue;
      
      // Identical partitions are known positive
      if (partition[s1] == partition[s2])
      {
        pp = new Pair;
        pp->x = s1;
        pp->y = s2;
#ifdef OPT_CACHE_NETS
        pp->simulation = 0;
#endif//OPT_CACHE_NETS
        pp->next = relation;
        relation = pp;
        continue;
      }
      
      // Look up the current pairs of states in the partition map (by their respective partitions)
      switch (partition_map[(partitions * partition[s1]) + partition[s2]])
      {
      case '-':
        rmap.Clear(s1, s2);
        --new_size;
        break; // Known failure; skip pair
      case '+': // Known success; add pair
        pp = new Pair;
        pp->x = s1;
        pp->y = s2;
#ifdef OPT_CACHE_NETS
        pp->simulation = 0;
#endif//OPT_CACHE_NETS
        pp->next = relation;
        relation = pp;
        break;
      default: // Unknown result; execute simulation and write result to the partition map
        pair.x = s1;
        pair.y = s2;
#ifdef OPT_CACHE_NETS
        pair.simulation = 0;
#endif//OPT_CACHE_NETS
        if (DecideStrongSimulation(&pair))
        {
          pp = new Pair;
          pp->x = s1;
          pp->y = s2;
#ifdef OPT_CACHE_NETS
          pp->simulation = pair.simulation;
          pair.simulation = NULL;
#endif//OPT_CACHE_NETS
          pp->next = relation;
          relation = pp;
          c = '+';
        }
        else
        {
#         ifdef OPT_CACHE_NETS
            assert(NULL == pair.simulation);
#         endif
          c = '-';
          --new_size;
          rmap.Clear(s1, s2);
        }
        partition_map[(partitions * partition[s1]) + partition[s2]] = c;
      }
    }
  }
  
#ifdef DEBUG
  stats.mem_relation = sizeof(Pair) * new_size;
#endif//DEBUG
  
  // Drop partition map; this is only valid for the duration of the first partition
  delete [] partition_map;
  delete [] partition;
  partition = NULL;
  
  // Copy the new relation map back into the main buffer
  rmap.Commit();
  
  return new_size;
}

// Construct the maxflow problem for a certain pair
CompactMaxFlow<double> *StrongSimulation_MC::ConstructNetwork(Pair *p)
{
  int &s1 = p->x, &s2 = p->y; // State IDs
  bool result, known_result;
  int s1_suc = row_starts[s1+1] - row_starts[s1]; // Get number of successor states
  int s2_suc = row_starts[s2+1] - row_starts[s2]; //
  CompactMaxFlow<double> *sim = new CompactMaxFlow<double>;
  
  // Create the network
  result = sim->CreateNetwork(cols, non_zeros, &rmap, row_starts[s1], s1_suc, row_starts[s2], s2_suc, known_result, 0);
  if (known_result && !result)
  {
    delete sim;
    return 0;
  }
  else assert(result);
  
  return sim;
}

// Compute strong simulation on a particular pair in the relation
bool StrongSimulation_MC::DecideStrongSimulation(Pair *p)
{
  bool r;
  CompactMaxFlow<double> *sim;
  
  if (p->x == p->y) return true;
  if (row_starts[p->x] == row_starts[p->x + 1]) return true;
  if (row_starts[p->y] == row_starts[p->y + 1]) return false;
  
  // If the simulation network does not exist, create it. If creation fails, it means that the result
  // is already known to be false; construction error will cause an assertion failure. If the network
  // exists, update it.
#ifdef OPT_CACHE_NETS
  if (!p->simulation)
  {
    if (!(p->simulation = ConstructNetwork(p))) return false;
#ifdef DEBUG
    else ++stats.num_nets_cached;
#endif//DEBUG
  }
  else
  {
    p->simulation->UpdateNetwork(&rmap, false);
#ifdef DEBUG
    ++stats.num_cache_hits;
#endif//DEBUG
  }
  sim = p->simulation;
#else//OPT_CACHE_NETS
  sim = ConstructNetwork(p);
  if (!sim) return false;
#endif//OPT_CACHE_NETS
  
  // Evaluate simulation
  r = sim->IsFlowTotal();
  
#ifdef OPT_CACHE_NETS
  if (!r)
  {
    delete p->simulation;
    p->simulation = 0;
  }
#else//OPT_CACHE_NETS
  delete sim;
#endif//OPT_CACHE_NETS

return r;
}

#if defined(OPT_PARTITION)
// Compute the first partition for the state and network partitioning schemes.
// For this purpose, states are ordered by the order defined by state_Less() and subsequent
// states that are identical under !state_Different() are assigned the same partition.
void StrongSimulation_MC::MakeFirstPartition()
{
  int n, cur_part, np, s;
  StateOrder cmp(this);
  
  order = new int[n_states];
  for (n = 0; n < n_states; ++n) order[n] = n;
  
  // The array order[] defines the permutation that puts the states in the desired order.
  // The order in the sparse matrix is left untouched.
  std::sort(order, order + n_states, cmp);
  
  np = 1;
  
  // Determine partitions
  partitions = 1;
  partition = new int[n_states];
  memset(partition, 0, sizeof(int) * n_states);
  for (n = 1; n < n_states; ++n)
  {
    if (partition[order[n - 1]] == partition[order[n]])
    {
      if (state_Different(order[n - 1], order[n]))
      {
        cur_part = partition[order[n]];
        for (s = n; s < n_states && partition[order[s]] == cur_part; ++s) partition[order[s]] = np;
        ++np;
        ++partitions;
      }
    }
  }
  delete [] order;
  order = NULL;
}
#endif

// Determine whether state s1 is less than state s2 (this is a largely arbitrary relation which helps partition the states)
bool StrongSimulation_MC::state_Less(int s1, int s2)
{
  int n, successors, l1, l2;
  
  if (row_starts[s1 + 1] - row_starts[s1] < row_starts[s2 + 1] - row_starts[s2]) return true;
  if (row_starts[s1 + 1] - row_starts[s1] > row_starts[s2 + 1] - row_starts[s2]) return false;
  
  for (n = 0, successors = row_starts[s1 + 1] - row_starts[s1]; n < successors; ++n)
  {
    //if (non_zeros[row_starts[s1] + n] < non_zeros[row_starts[s2] + n]) return true;
    //else if (non_zeros[row_starts[s1] + n] > non_zeros[row_starts[s2] + n]) return false;
    if (CompactMaxFlow<double>::_Tless(non_zeros[row_starts[s1] + n], non_zeros[row_starts[s2] + n])) return true;
    else if (CompactMaxFlow<double>::_Tless(non_zeros[row_starts[s2] + n], non_zeros[row_starts[s1] + n])) return false;
    l1 = Label(cols[row_starts[s1] + n]);
    l2 = Label(cols[row_starts[s2] + n]);
    if (l1 < l2) return true;
    else if (l1 > l2) return false;
  }
  
  return false;
}

// Determine whether two states are probabilistically different (includes labels)
bool StrongSimulation_MC::state_Different(int s1, int s2)
{
  int n, successors;
  
  if (row_starts[s1 + 1] - row_starts[s1] != row_starts[s2 + 1] - row_starts[s2]) return true;
  
  for (n = 0, successors = row_starts[s1 + 1] - row_starts[s1]; n < successors; ++n)
  {
    //if (non_zeros[row_starts[s1] + n] != non_zeros[row_starts[s2] + n]) return true;
    if (!CompactMaxFlow<double>::_Teq(non_zeros[row_starts[s1] + n], non_zeros[row_starts[s2] + n])) return true;
    if (Label(cols[row_starts[s1] + n]) != Label(cols[row_starts[s2] + n])) return false;
  }
  
  return false;
}

// Sort the successors of each state by probability (1st) and label (2nd)
void StrongSimulation_MC::SortSuccessors()
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
