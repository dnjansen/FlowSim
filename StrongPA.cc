#include "Strong.h"
#include "compactmaxflow.cc"

unsigned int StrongSimulation_PA::Simulate(ProbabilisticModel *model, std::set<std::pair<int,int> > *result)
{
  int r, iterations = 0;
  Pair *pp;
  
  if (!model || model->Type() != ProbabilisticModel::PA) return 0;
  
  // Copy sparse matrix structures
  m = (ProbabilisticAutomaton*)model;
  n_states = m->n;
  state_starts = m->state_starts;
  
  non_zeros = new double[m->nnz];
  memcpy(non_zeros, m->non_zeros, sizeof(double) * m->nnz);
  cols = new int[m->nnz];
  memcpy(cols, m->cols, sizeof(int) * m->nnz);
  actions = new int[m->na];
  memcpy(actions, m->actions, sizeof(int) * m->na);
  row_starts = new int[m->na + 1];
  memcpy(row_starts, m->row_starts, sizeof(int) * (m->na + 1));
  
  rmap.Create(1, n_states);
  
  action_masks = 0;
  
#ifdef DEBUG
  memset(&stats, 0, sizeof(stats));
  CompactMaxFlow<double>::ResetStats();
  stats.mem_model = (sizeof(double) * m->nnz) + (sizeof(int) * (m->nnz + m->na + m->na + 1));
#endif//DEBUG

  //if (optflags & OPT_PARTITION)
#ifdef OPT_PARTITION
  {
    SortSuccessors();
    MakeFirstPartition();
    
    stats.num_partitions = partitions;
    stats.mem_partition_map = partitions * partitions;
  }
#endif

  InitializeActionMasks();

  if (m->ContinuousTimeModel()) size_of_relation = BuildRelationMap_CPA();
  else size_of_relation = BuildRelationMap_PA();
  
#ifdef DEBUG
  stats.mem_relation = sizeof(Pair) * size_of_relation;
  stats.num_initial_pairs = size_of_relation;
#endif//DEBUG

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
  
#ifdef DEBUG
  stats.mem_relation_map = rmap.MemoryUsage();
  stats.num_iterations = iterations;
  stats.num_final_pairs = size_of_relation;
  stats.mem_maxflow = CompactMaxFlow<double>::global_space_peak;
  stats.num_maxflow = CompactMaxFlow<double>::global_times_invoked;
  stats.num_p_invariant_fails = CompactMaxFlow<double>::global_p_inv_fails;
  stats.num_sig_arc_fails = CompactMaxFlow<double>::global_sig_arc_fails;
  CompactMaxFlow<double>::ResetStats();
#endif//DEBUG

  //if (optflags & OPT_PARTITION)
#ifdef OPT_PARTITION
  {
    delete [] order;
    delete [] partition;
  }
#endif
  
  delete [] cols;
  delete [] non_zeros;
  delete [] actions;
  delete [] row_starts;
  delete [] action_masks;
  
  // Free relation
  pp = relation;
  while (pp)
  {
    relation = pp->next;
    delete pp;
    pp = relation;
  }
  
  return size_of_relation;
}

#ifdef WITH_VERIFIER
// Verify that a set of pairs is a simulation relation for the given model
bool StrongSimulation_PA::Verify(ProbabilisticModel *model, std::set<std::pair<int,int> > &hypothesis,
        std::set<std::pair<int,int> > *false_positives, std::set<std::pair<int,int> > *false_negatives)
{
  bool res, global_res = true, cpa;
  std::set<std::pair<int,int> >::iterator hi;
  Pair p;
  double *dist_sums;
  
  if (!model || model->Type() != ProbabilisticModel::PA) return false;
  
  // Copy sparse matrix structures
  m = (ProbabilisticAutomaton*)model;
  non_zeros = m->non_zeros;
  cols = m->cols;
  row_starts = m->row_starts;
  state_starts = m->state_starts;
  actions = m->actions;
  n_states = m->n;
  
  // For CPAs, compute distribution sums and normalize
  if (cpa = model->ContinuousTimeModel())
  {
    non_zeros = new double[m->na];
    memcpy(non_zeros, m->non_zeros, sizeof(double) * m->na);
    dist_sums = new double[m->na];
    for (int l, n = 0; n < m->na; ++n)
    {
      dist_sums[n] = 0.0;
      for (l = row_starts[n]; l < row_starts[n + 1]; ++m) dist_sums[n] += non_zeros[l];
      for (l = row_starts[n]; l < row_starts[n + 1]; ++m) non_zeros[l] /= dist_sums[n];
    }
  }
  
  rmap.Create(n_states);
  
  size_of_relation = hypothesis.size();
  
  for (hi = hypothesis.begin(); hi != hypothesis.end(); ++hi)
  {
    rmap.Set(hi->first, hi->second);
  }
  
  rmap.Commit();
  
  for (p.x = 0; p.x < n_states; ++p.x)
  {
    for (p.y = 0; p.y < n_states; ++p.y)
    {
      if (p.x == p.y) continue;
      res = (Label(p.x) == Label(p.y)
             && (cpa ? CompactMaxFlow<double>::_Tleq(dist_sums[p.x], dist_sums[p.y]) : true)
             && DecideStrongSimulation(&p));
      if (res != rmap(p.x, p.y))
      {
        global_res = false;
        if (res && !rmap(p.x, p.y) && false_negatives) false_negatives->insert(std::make_pair(p.x, p.y));
        else if (!res && rmap(p.x, p.y) && false_positives) false_positives->insert(std::make_pair(p.x, p.y));
      }
    }
  }
  
  if (cpa)
  {
    delete [] non_zeros;
    delete [] dist_sums;
  }
  
  return global_res;
}
#endif//WITH_VERIFIER

// Compute the initial relation for deterministic time PAs.
int StrongSimulation_PA::BuildRelationMap_PA()
{
  int m, n, size = 0;
  bool forward, backward;

  for (m = 0; m < n_states - 1; ++m)
  {
    for (n = m + 1; n < n_states; ++n)
    {
      if (Label(m) == Label(n))
      {
        // For L(m) = L(n), check that Act(m) \subset Act(n) and vice versa, unless
        // partitioning is activated
        forward = true;
        backward = true;

#ifndef OPT_PARTITION
        for (int i = 0; i < action_mask_pitch && (forward || backward); ++i)
        {
          if ((action_masks[(m * action_mask_pitch) + i] ^ action_masks[(n * action_mask_pitch) + i]) & ~(action_masks[(n * action_mask_pitch) + i])) forward = false;
          if ((action_masks[(n * action_mask_pitch) + i] ^ action_masks[(m * action_mask_pitch) + i]) & ~(action_masks[(m * action_mask_pitch) + i])) backward = false;
        }
#endif
        
        if (forward)  rmap.Set(m, n), ++size;
        if (backward) rmap.Set(n, m), ++size;
      }
    }
  }
  
  rmap.Commit();
  
  return size;
}

// Compute the initial relation for continuous time PAs.
int StrongSimulation_PA::BuildRelationMap_CPA()
{
  int m, n, size = 0;
  bool forward, backward;
  double *dist_sums = new double[this->m->na];
  
  // Compute sum for each distribution and normalize
  for (n = 0; n < this->m->na; ++n)
  {
    dist_sums[n] = 0.0;
    for (m = row_starts[n]; m < row_starts[n + 1]; ++m) dist_sums[n] += non_zeros[m];
    for (m = row_starts[n]; m < row_starts[n + 1]; ++m) non_zeros[m] /= dist_sums[n];
  }

  for (m = 0; m < n_states - 1; ++m)
  {
    for (n = m + 1; n < n_states; ++n)
    {
      if (Label(m) == Label(n))
      {
        // For L(m) = L(n), check that Act(m) \subset Act(n) and vice versa, unless
        // partitioning is activated, and also check that R(m,*) <= R(n,*) and vice versa
        forward =  (dist_sums[m] <= dist_sums[n]);
        backward = (dist_sums[n] <= dist_sums[m]);

#ifndef OPT_PARTITION
        for (int i = 0; i < action_mask_pitch && (forward || backward); ++i)
        {
          if ((action_masks[(m * action_mask_pitch) + i] ^ action_masks[(n * action_mask_pitch) + i]) & ~(action_masks[(n * action_mask_pitch) + i])) forward = false;
          if ((action_masks[(n * action_mask_pitch) + i] ^ action_masks[(m * action_mask_pitch) + i]) & ~(action_masks[(m * action_mask_pitch) + i])) backward = false;
        }
#endif
        
        if (forward)  rmap.Set(m, n), ++size;
        if (backward) rmap.Set(n, m), ++size;
      }
    }
  }
  
  rmap.Commit();
  
  return size;
}

// Naively iterate the relation by repeatedly invoking StrongSimulation
// on every pair.
int StrongSimulation_PA::IterateRelation(bool first)
{
  int s1, s2, new_size = size_of_relation;
  Pair p, *pp, **anchor = &relation;
  
  // If we are in the first iteration, add those pairs to the relation
  // that simulate and satisfy the initial condition. Otherwise, delete
  // the pairs that don't simulate.
  if (first)
  {
    relation = 0;
    for (s1 = 0; s1 < n_states; ++s1)
    {
      for (s2 = 0; s2 < n_states; ++s2)
      {
        if (s1 == s2 || !rmap(s1, s2)) continue;
        p.x = s1;
        p.y = s2;
        if (DecideStrongSimulation(&p))
        {
          pp = new Pair;
          pp->x = s1;
          pp->y = s2;
          pp->next = relation;
          relation = pp;
        }
        else
        {
          rmap.Clear(s1, s2);
          --new_size;
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
int StrongSimulation_PA::IterateRelation_FirstPartition()
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
        pp->next = relation;
        relation = pp;
        break;
      default: // Unknown result; execute simulation and write result to the partition map
        pair.x = s1;
        pair.y = s2;
        if (DecideStrongSimulation(&pair))
        {
          pp = new Pair;
          pp->x = s1;
          pp->y = s2;
          pp->next = relation;
          relation = pp;
          c = '+';
        }
        else
        {
          c = '-';
          --new_size;
          rmap.Clear(s1, s2);
        }
        partition_map[(partitions * partition[s1]) + partition[s2]] = c;
      }
    }
  }
  
  // Drop partition map; this is only valid for the duration of the first partition
  delete [] partition_map;
  
  // Copy the new relation map back into the main buffer
  rmap.Commit();
  
  return new_size;
}

// Compute strong simulation on a particular pair in the relation
bool StrongSimulation_PA::DecideStrongSimulation(Pair *p)
{
  int a1, a2;
  bool sim_one;
  CompactMaxFlow<double> *sim;
  
#ifdef OPT_PARTITION
  // This test can't be done in the initial relation if we are using state partitioning
  for (int i = 0; i < action_mask_pitch; ++i)
  {
    if ((action_masks[(p->x * action_mask_pitch) + i] ^ action_masks[(p->y * action_mask_pitch) + i]) & ~(action_masks[(p->y * action_mask_pitch) + i])) return false;
  }
#endif
  
  for (a1 = state_starts[p->x]; a1 < state_starts[p->x + 1]; ++a1)
  {
    sim_one = false;
    for (a2 = state_starts[p->y]; a2 < state_starts[p->y + 1]; ++a2)
    {
      if (actions[a1] == actions[a2])
      {
        if (row_starts[a1] == row_starts[a1 + 1])
        {
          sim_one = true;
          break;
        }
        
        sim = ConstructNetwork(a1, a2);
        if (sim)
        {
          sim_one = sim->IsFlowTotal();
          delete sim;
          if (sim_one) break;
        }
      }
    }
    
    if (!sim_one) return false;
  }
  
  return true;
}

// Construct the maxflow problem for a certain pair
CompactMaxFlow<double> *StrongSimulation_PA::ConstructNetwork(int a1, int a2)
{
  bool result, known_result;
  int s1_suc = row_starts[a1+1] - row_starts[a1]; // Get number of successor states
  int s2_suc = row_starts[a2+1] - row_starts[a2]; //

  // Create the network
  CompactMaxFlow<double> *simulation = new CompactMaxFlow<double>;
  //result = simulation->CreateNetwork(cols, non_zeros, &rmap, row_starts[a1], s1_suc, row_starts[a2], s2_suc, known_result, optflags);
  result = simulation->CreateNetwork(cols, non_zeros, &rmap, row_starts[a1], s1_suc, row_starts[a2], s2_suc, known_result, 0);
  if (known_result && !result)
  {
    delete simulation;
    return 0;
  }
  
  return simulation;
}

// Sort the successors of each state by action (1st) probability (2nd) and label (3rd)
void StrongSimulation_PA::SortSuccessors()
{
  int i, j, rs;
  
  int *action_order = new int[m->na + 1];
  int *successor_order = new int[m->nnz + 1];
  
  int *new_cols = new int[m->nnz];
  double *new_nz = new double[m->nnz];
  int *new_actions = new int[m->na];
  int *new_row_starts = new int[m->na + 1];
  
  ActionOrder acmp(this);
  SuccessorOrder scmp(this);
  
  for (i = 0; i < m->na; ++i) action_order[i] = i;
  for (i = 0; i < m->nnz; ++i) successor_order[i] = i;
  
  // Sort successor set for each action
  for (i = 0; i < n_states; ++i)
  {
    for (j = state_starts[i]; j < state_starts[i+1]; ++j)
    {
      if (row_starts[j+1] - row_starts[j] <= 1) continue;
      std::sort(successor_order + row_starts[j], successor_order + row_starts[j+1], scmp);
    }
  }
  
  // Apply new order
  //NOTE: Successors MUST be sorted for action sorting to work correctly.
  for (i = 0; i < m->nnz; ++i)
  {
    new_cols[i] = cols[successor_order[i]];
    new_nz[i] = non_zeros[successor_order[i]];
  }
  memcpy(cols, new_cols, sizeof(int) * m->nnz);
  memcpy(non_zeros, new_nz, sizeof(double) * m->nnz);
  
  // Sort actions for each state
  for (i = 0, rs = 0; i < n_states; ++i)
  {
    if (state_starts[i+1] - state_starts[i] == 0) continue;
    else if (state_starts[i+1] - state_starts[i] == 1)
    {
      rs += row_starts[state_starts[i] + 1] - row_starts[state_starts[i]];
      continue;
    }
    std::sort(action_order + state_starts[i], action_order + state_starts[i+1], acmp);
    for (j = state_starts[i]; j < state_starts[i+1]; ++j)
    {
      memcpy(new_cols + rs, cols + row_starts[action_order[j]], sizeof(int) * (row_starts[action_order[j] + 1] - row_starts[action_order[j]]));
      memcpy(new_nz + rs, non_zeros + row_starts[action_order[j]], sizeof(double) * (row_starts[action_order[j] + 1] - row_starts[action_order[j]]));
      rs += row_starts[action_order[j] + 1] - row_starts[action_order[j]];
    }
  }
  
  // Prepare row_starts for re-ordering
  for (i = 0; i < m->na; ++i) row_starts[m->na - i] -= row_starts[m->na - i - 1];
  
  // Create new arrays with the updated order
  new_row_starts[0] = 0;
  for (i = 0; i < m->na; ++i)
  {
    new_actions[i] = actions[action_order[i]];
    new_row_starts[i + 1] = row_starts[action_order[i] + 1];
  }
  
  // Fix the new row_starts array
  for (i = 0; i < m->na; ++i) new_row_starts[i + 1] += new_row_starts[i];
  
  delete [] successor_order;
  delete [] action_order;
  
  // Delete old arrays and copy new ones over
  delete [] cols;
  delete [] non_zeros;
  delete [] actions;
  delete [] row_starts;
  
  cols = new_cols;
  non_zeros = new_nz;
  actions = new_actions;
  row_starts = new_row_starts;
}

// Compute action masks. The action mask is a bitset where every bit corresponds to a certain
// action. If a state has at least one non-zero distribution for a certain action, the bit is
// on, otherwise it is off. This allows us to perform the test Act(s)=Act(s') in O(l) time
// where l is the number of bytes required to represent the mask.
void StrongSimulation_PA::InitializeActionMasks()
{
  int i, j;
  
  action_mask_pitch = (m->da >> 3) + ((m->da & 0x7) ? 1 : 0);
  action_masks = new unsigned char[n_states * action_mask_pitch];
  memset(action_masks, 0, n_states * action_mask_pitch);
  
  for (i = 0; i < n_states; ++i)
  {
    for (j = state_starts[i]; j < state_starts[i+1]; ++j)
    {
      action_masks[(i * action_mask_pitch) + (actions[j] >> 3)] |= (1 << (actions[j] & 0x7));
    }
  }
}

// Compute the first partition for the state and network partitioning schemes.
// For this purpose, states are ordered by the order defined by state_Less() and subsequent
// states that are identical under !state_Different() are assigned the same partition.
void StrongSimulation_PA::MakeFirstPartition()
{
  int n, cur_part, np, s;
  StateOrder cmp(this);
  
  order = new int[n_states + 1];
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
}

// Determine whether state s1 is less than state s2 (this is a largely arbitrary relation which helps partition the states)
bool StrongSimulation_PA::state_Less(int s1, int s2)
{
  int n, successors, a, a1, a2;

  if (s1 == s2) return false;
  
  if (state_starts[s1 + 1] - state_starts[s1] < state_starts[s2 + 1] - state_starts[s2]) return true;
  if (state_starts[s1 + 1] - state_starts[s1] > state_starts[s2 + 1] - state_starts[s2]) return false;
  
  for (a = 0; a < state_starts[s1 + 1] - state_starts[s1]; ++a)
  {
    a1 = state_starts[s1] + a;
    a2 = state_starts[s2] + a;
    
    if (actions[a1] < actions[a2]) return true;
    if (actions[a1] > actions[a2]) return false;
    
    if (row_starts[a1 + 1] - row_starts[a1] < row_starts[a2 + 1] - row_starts[a2]) return true;
    if (row_starts[a1 + 1] - row_starts[a1] > row_starts[a2 + 1] - row_starts[a2]) return false;
    
    for (n = 0, successors = row_starts[a1 + 1] - row_starts[a1]; n < successors; ++n)
    {
      if (CompactMaxFlow<double>::_Tless(non_zeros[row_starts[a1] + n], non_zeros[row_starts[a2] + n])) return true;
      else if (CompactMaxFlow<double>::_Tless(non_zeros[row_starts[a2] + n], non_zeros[row_starts[a1] + n])) return false;
      if (Label(cols[row_starts[a1] + n]) < Label(cols[row_starts[a2] + n])) return true;
      if (Label(cols[row_starts[a1] + n]) > Label(cols[row_starts[a2] + n])) return false;
    }
  }
  
  return false;
}

// Determine whether two states are probabilistically different (includes labels)
bool StrongSimulation_PA::state_Different(int s1, int s2)
{
  int n, successors, a, a1, a2;
  
  if (state_starts[s1 + 1] - state_starts[s1] != state_starts[s2 + 1] - state_starts[s2]) return true;
  
  for (a = 0; a < state_starts[s1 + 1] - state_starts[s1]; ++a)
  {
    a1 = state_starts[s1] + a;
    a2 = state_starts[s2] + a;
    
    if (actions[a1] != actions[a2]) return true;
    
    if (row_starts[a1 + 1] - row_starts[a1] != row_starts[a2 + 1] - row_starts[a2]) return true;
    
    for (n = 0, successors = row_starts[a1 + 1] - row_starts[a1]; n < successors; ++n)
    {
      if (!CompactMaxFlow<double>::_Teq(non_zeros[row_starts[a1] + n], non_zeros[row_starts[a2] + n])) return true;
      if (Label(cols[row_starts[a1] + n]) != Label(cols[row_starts[a2] + n])) return true;
    }
  }
  
  return false;
}
