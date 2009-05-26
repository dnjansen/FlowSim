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



#include "StrongQ.h"
#include "compactmaxflow.cc"

#include <boost/config.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/transitive_closure.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>

using namespace std;

// Compute the simulation relation
unsigned int StrongSimulation_Quotient::Simulate(ProbabilisticModel *target, set<pair<int,int> > *result)
{
  int n, m, size;
  
  // Assert that we have valid pointers
  assert(target);
  assert(label_func);
  
#ifdef QLOG
  qlog = fopen("quotientlog", "wb");
#endif
  
  // This algorithm can't handle CTMC/CPA, only DTMC/PA
  assert(!target->ContinuousTimeModel());
  
  // Get model size parameters
  model = target;
  nStates = model->States();
  nDistributions = model->Distributions();
  
  // Copy pointers to sparse matrix and set up model type abstraction
  switch (model->Type())
  {
  // Markov Chain
  case ProbabilisticModel::MC:
    row_starts = ((MarkovChain*)model)->row_starts;
    cols = ((MarkovChain*)model)->cols;
    non_zeros = ((MarkovChain*)model)->non_zeros;
    state_starts = 0;
    actions = 0;
    sim = new Simulator_MC(this);
    break;
  // Probabilistic Automaton
  case ProbabilisticModel::PA:
    row_starts = ((ProbabilisticAutomaton*)model)->row_starts;
    cols = ((ProbabilisticAutomaton*)model)->cols;
    non_zeros = ((ProbabilisticAutomaton*)model)->non_zeros;
    state_starts = ((ProbabilisticAutomaton*)model)->state_starts;
    actions = ((ProbabilisticAutomaton*)model)->actions;
    sim = new Simulator_PA(this);
    break;
  default:
    assert(false);
  }
  
  // Allocate space for quotient automaton and auxiliary structures
  // (allocate worst case space so we don't have to worry about that later)
  lifted_non_zeros = new double[model->Transitions()];
  lifted_cols = new int[model->Transitions()];
  lifted_row_starts = new int[model->Distributions() + 1];
  partition = new int[model->States() * 2];
  new_partition = partition + model->States();
  
  InitializeActionMasks();
  InitializeRelation();

#ifdef DEBUG
  memset(&stats, 0, sizeof(stats));
  CompactMaxFlow<double>::ResetStats();
  stats.mem_model = (sizeof(double) * model->Transitions())
                  + (sizeof(int) * (model->Transitions() + model->Distributions()))
                  + (sizeof(unsigned char) * action_mask_pitch * model->States());
  stats.mem_relation = 0;
  stats.num_initial_pairs = nBlocks;
  stats.mem_partition_map = 3 * nStates * sizeof(int);
#endif//DEBUG

#ifdef QLOG
  fprintf(qlog, "\n===== Initialization =====\nBlocks:");
  for (n = 0; n < nBlocks; ++n)
  {
    fprintf(qlog, " %d = {", n);
    for (int m = 0, i = 0; m < nStates; ++m)
    {
      if (partition[m] == n) fprintf(qlog, "%s%d", (i++ == 0 ? "" : ","), m);
    }
    fprintf(qlog, "}");
  }
  fprintf(qlog, "\nPartition Relation:");
  for (n = 0; n < nBlocks; ++n)
  {
    for (int m = 0; m < nBlocks; ++m)
    {
      if (rmap(n, m) && n != m) fprintf(qlog, " (%d,%d)", n, m);
    }
  }
  fprintf(qlog, "\n\n");
#endif

  // Find maximum stable partition pair
  for (iteration = 1; Iterate() != -1; ++iteration)
  {
#ifdef QLOG
    fprintf(qlog, "\n===== After Iteration %d =====\nBlocks:", iteration);
    for (n = 0; n < nBlocks; ++n)
    {
      fprintf(qlog, " %d = {", n);
      for (int m = 0, i = 0; m < nStates; ++m)
      {
        if (partition[m] == n) fprintf(qlog, "%s%d", (i++ == 0 ? "" : ","), m);
      }
      fprintf(qlog, "}");
    }
    fprintf(qlog, "\nPartition Relation:");
    for (n = 0; n < nBlocks; ++n)
    {
      for (int m = 0; m < nBlocks; ++m)
      {
        if (rmap(n, m) && n != m) fprintf(qlog, " (%d,%d)", n, m);
      }
    }
    fprintf(qlog, "\n\n");
#endif
  }
  
#ifdef QLOG
  fprintf(qlog, "\n===== After Final Iteration =====\nBlocks:");
  for (n = 0; n < nBlocks; ++n)
  {
    fprintf(qlog, " %d = {", n);
    for (int m = 0, i = 0; m < nStates; ++m)
    {
      if (partition[m] == n) fprintf(qlog, "%s%d", (i++ == 0 ? "" : ","), m);
    }
    fprintf(qlog, "}");
  }
  fprintf(qlog, "\nPartition Relation:");
  for (n = 0; n < nBlocks; ++n)
  {
    for (int m = 0; m < nBlocks; ++m)
    {
      if (rmap(n, m) && n != m) fprintf(qlog, " (%d,%d)", n, m);
    }
  }
  fprintf(qlog, "\n\n");
#endif

  // Count pairs in the actual relation and place result in the provided set
  rmap.ReportCurrent(false);
  size = 0;
  if (result) result->clear();
  for (n = 0; n < nStates; ++n)
  {
    for (m = 0; m < nStates; ++m)
    {
      if (n != m && (partition[n] == partition[m] || rmap(partition[n], partition[m])))
      {
        ++size;
        if (result) result->insert(std::make_pair(n, m));
      }
    }
  }

#ifdef DEBUG
  stats.mem_relation_map = rmap.MemoryUsage();
  stats.num_partitions = nBlocks;
  stats.num_iterations = iteration;
  stats.num_final_pairs = size;
  stats.mem_maxflow = CompactMaxFlow<double>::global_space_peak;
  stats.num_maxflow = CompactMaxFlow<double>::global_times_invoked;
  stats.num_p_invariant_fails = CompactMaxFlow<double>::global_p_inv_fails;
  stats.num_sig_arc_fails = CompactMaxFlow<double>::global_sig_arc_fails;
  stats.min_complexity = CompactMaxFlow<double>::min_complexity;
  stats.max_complexity = CompactMaxFlow<double>::max_complexity;
  CompactMaxFlow<double>::ResetStats();
#endif//DEBUG

  delete [] lifted_row_starts;
  delete [] lifted_cols;
  delete [] lifted_non_zeros;
  delete [] partition;
  delete sim;
  sigma.clear();
  
  if (action_masks) delete [] action_masks;
  
  return size;
}

// Compute action masks. The action mask is a bitset where every bit corresponds to a certain
// action. If a state has at least one non-zero distribution for a certain action, the bit is
// on, otherwise it is off. This allows us to perform the test Act(s)=Act(s') in O(l) time
// where l is the number of bytes required to represent the mask.
void StrongSimulation_Quotient::InitializeActionMasks()
{
  int i, j, a;
  
  if (model->Type() != ProbabilisticModel::PA)
  {
    action_masks = 0;
    return;
  }
  
  a = ((ProbabilisticAutomaton*)model)->da;
  action_mask_pitch = (a >> 3) + ((a & 0x7) ? 1 : 0);
  action_masks = new unsigned char[nStates * action_mask_pitch];
  memset(action_masks, 0, nStates * action_mask_pitch);
  
  for (i = 0; i < nStates; ++i)
  {
    for (j = state_starts[i], a = -1; j < state_starts[i+1]; ++j)
    {
      if (actions[j] == a) continue;
      a = actions[j];
      action_masks[(i * action_mask_pitch) + (a >> 3)] |= (1 << (a & 0x7));
    }
  }
}

// Create the pairs in the initial relation
void StrongSimulation_Quotient::InitializeRelation()
{
  int n, m, p, a, b;
  
  // Create initial partition into blocks by L(s) and Act(s)
  memset(partition, 0, sizeof(int) * nStates);
  for (n = 0, p = 0; n < nStates - 1; )
  {
    ++p;
    a = 0;
    
    for (m = n + 1; m < nStates; ++m)
    {
      if (partition[n] == partition[m] && (Label(n) != Label(m) || !Act_Equal(n, m))) partition[m] = p, a = 1;
    }
    
    while (partition[n] != p && n < nStates - 1) ++n;
  }
  
  nBlocks = p + a;
  
  rmap.Create(nBlocks, nStates);

  // Build relation map for the initial blocks we have. All pairs of blocks (B,B')
  // where L(B) = L(B') and Act(B) \subset_eq Act(B') are in the initial relation.
  sigma.resize(nBlocks, std::set<int>());
  for (n = 0; n < nStates; ++n)
  {
    sigma[partition[n]].insert(n);
    for (m = 0; m < nStates; ++m)
    {
      if (n < nBlocks && m < nBlocks)
      {
        if (n == m) rmap.Set(n, m);
        else
        {
          a = 0, b = 0;
          while (partition[a] != n) ++a;
          while (partition[b] != m) ++b;
          if (Label(a) != Label(b))
          {
            rmap.Clear(n, m);
            continue;
          }
          if (Act_Subseteq(a, b)) rmap.Set(n, m);
          else rmap.Clear(n, m);
        }
      }
    }
  }
  rmap.Commit();
  
  // Compute lifted distributions for initial partition pair
  LiftDistributions();
}

// Compute lifted distributions for current partition
void StrongSimulation_Quotient::LiftDistributions()
{
  using namespace std;
  
  int i, j, s;
  map<int,double> mu;
  map<int,double>::iterator mi;
  vector<set<int> >::iterator si;
  set<int>::iterator qi;
  
#ifdef QLOG
  fprintf(qlog, "[LIFT] Computing lifted distributions, total %d dist, %d blocks\n", nDistributions, nBlocks);
#endif
  
  // Lift distributions
  lifted_row_starts[0] = 0;
  for (i = 0; i < nDistributions; ++i) // For each distribution i...
  {
    for (j = row_starts[i]; j < row_starts[i + 1]; ++j) // For each transition j of i...
    {
      // Lift distribution i from state space to partition space
      s = partition[cols[j]];
      if ((mi = mu.find(s)) == mu.end()) mu.insert(make_pair(s, non_zeros[j]));
      else mi->second += non_zeros[j];
    }
    
    // Write lifted distribution into quotient automaton structures
    j = lifted_row_starts[i];
    for (mi = mu.begin(); mi != mu.end(); ++mi, ++j)
    {
      lifted_cols[j] = mi->first;
      lifted_non_zeros[j] = mi->second;
    }
    lifted_row_starts[i + 1] = j;
    
    mu.clear();
  }
  
  if (model->Type() == ProbabilisticModel::PA) FindCommonDistributions_PA();
  else if (model->Type() == ProbabilisticModel::MC) FindCommonDistributions_MC();
}

// Find \forall distributions in PA-type models
void StrongSimulation_Quotient::FindCommonDistributions_PA()
{
  int i, j, s;
  set<pair<int,int> > candidates, keep;
  set<pair<int,int> >::iterator mu;
  //pair<multimap<int,int>::iterator,multimap<int,int>::iterator> ci;
  vector<set<int> >::iterator si;
  set<int>::iterator qi;
  
  forall.clear();
  
#ifdef QLOG
  fprintf(qlog, "[FALL] Computing \\forall-distributions for %d blocks\n", nBlocks);
#endif
  
  // Iterate through blocks to find \forall distributions
  for (si = sigma.begin(), j = 0; si != sigma.end(); ++si, ++j)
  {
    // Put the distributions of the first state in each block into
    // the set of candidates for \forall distributions
    s = *si->begin();
    candidates.clear();
    for (i = state_starts[s]; i < state_starts[s + 1]; ++i)
    {
      candidates.insert(make_pair(actions[i], i));
    }
    
    // Iterate through the states in the block (skipping the first) and find \forall distributions
    for (qi = ++(si->begin()); qi != si->end() && candidates.size() > 0; ++qi)
    {
      // Go through distributions of state *qi and keep those that match
      // a distribution in the candidate set
      for (i = state_starts[*qi]; i < state_starts[*qi + 1]; ++i)
      {
        for (mu = candidates.begin(); mu != candidates.end(); ++mu)
        {
          if (actions[i] == mu->first && Dist_Equal(i, mu->second)) keep.insert(*mu);
        }
      }
      
      candidates = keep;
      keep.clear();
    }
    
    // Mark any remaining distributions as \forall distributions
    for (mu = candidates.begin(); mu != candidates.end(); ++mu) forall.insert(make_pair(j, *mu));
    
#ifdef QLOG
    fprintf(qlog, "[FALL] Block %d {", j);
    for (qi = si->begin(); qi != si->end(); ++qi) fprintf(qlog, "%s%d", (qi == si->begin() ? "" : ","), *qi);
    fprintf(qlog, "} >--> {");
    for (mu = candidates.begin(); mu != candidates.end(); ++mu) fprintf(qlog, "%s%d", (mu == candidates.begin() ? "" : ","), mu->second);
    fprintf(qlog, "}\n");
#endif
  }
}

// Find \forall distributions in MC-type models
void StrongSimulation_Quotient::FindCommonDistributions_MC()
{
  int j, candidate;
  vector<set<int> >::iterator si;
  set<int>::iterator qi;
  
  forall.clear();
  
  // Iterate through blocks to find \forall distributions
  for (si = sigma.begin(), j = 0; si != sigma.end(); ++si, ++j)
  {
    // Use the distribution of the first state as the candidate (only one for MC)
    candidate = *si->begin();
    
    // Iterate through the states in the block and find \forall distributions
    for (qi = si->begin(); qi != si->end(); ++qi)
    {
      if (!Dist_Equal(candidate, *qi)) break;
    }
    
    // Mark remaining distribution as \forall distribution
    if (qi == si->end()) forall.insert(make_pair(j, make_pair(0, candidate)));
  }
}

// Purge partition relation
void StrongSimulation_Quotient::PurgePartitionRelation()
{
  bool repeat;
  int i, j;
  set<pair<int,int> > relation, survivors;
  set<pair<int,int> >::iterator ri;
  
#ifdef QLOG
  int iter = 1;
  fprintf(qlog, "[PURG] Entering inner loop\n");
#endif

  // Make rmap(.,.) operate on the working copy
  rmap.ReportCurrent(true);
  
  // Iterate through entire map and put all pairs that simulate into the local set
  repeat = false;
  for (i = 0; i < nBlocks; ++i)
  {
    for (j = 0; j < nBlocks; ++j)
    {
      if (i == j || !rmap(i, j)) continue;
      
      if (!sim->qRq(i, j, false))
      {
#ifdef QLOG
        fprintf(qlog, "[PURG] Purged (%d,%d) from partition relation\n", i, j);
#endif
        rmap.Clear(i, j);
        repeat = true;
      }
      else relation.insert(make_pair(i, j));
    }
  }
  
#ifdef QLOG
  fprintf(qlog, "[PURG] End iteration 1\n");
#endif
  
  // Keep repeating inner block until nothing changes in the relation
  while (repeat)
  {
    repeat = false;
    
    // Iterate through all blocks
    for (ri = relation.begin(); ri != relation.end(); ++ri)
    {
      i = ri->first;
      j = ri->second;
      
      // Remove this pair from relation if j cannot simulate i
      if (!sim->qRq(i, j, false))
      {
#ifdef QLOG
        fprintf(qlog, "[PURG] Purged (%d,%d) from partition relation\n", i, j);
#endif
        rmap.Clear(i, j);
        repeat = true;
      }
      else survivors.insert(make_pair(i, j));
    }
    
    relation = survivors;
    survivors.clear();
#ifdef QLOG
    ++iter;
    fprintf(qlog, "[PURG] End iteration %d\n", iter);
#endif
  }
  
#ifdef QLOG
  fprintf(qlog, "[PURG] Leaving inner loop after %d iterations\n", iter);
#endif
  
  // Flush cached networks (if any) because they will be invalid in the next iteration
  sim->Flush();
}

// Improved partitioning algorithm which is applied to all iterations.
int StrongSimulation_Quotient::Iterate()
{
  using namespace boost;
  
  int i, j, nextblock, parentblocks, components;
  bool started_new, sigma_unchanged, gamma_unchanged;
  
  std::map<int,int> par;
  std::vector<std::set<int> > vertices;
  std::vector<int> block_rep, scc, scc_rev;
  std::vector<default_color_type> colormap;
  
  std::set<int>::iterator mia, mib;
  
  adjacency_list<vecS, vecS, directedS> *digraph, *digraph_closure;
  
#ifdef QLOG
  fprintf(qlog, "[REFN] Beginning refinement, %d starting blocks\n", nBlocks);
#endif

  memcpy(new_partition, partition, sizeof(int) * nStates);
  
  nextblock = nBlocks;
  parentblocks = nBlocks;
  
  // Make rmap(.,.) operate on base copy, changes have to be committed before taking effect
  rmap.ReportCurrent(false);
  
  // Initialize set of vertices Q
  vertices.resize(nBlocks, std::set<int>());
  for (i = 0; i < nBlocks; ++i) vertices[i].insert(i);
  
  // Find vertices Q (lines 5-7)
  block_rep.resize(nBlocks);
  for (i = 0; i < nStates; ++i)
  {
    started_new = false;
    for (j = i + 1; j < nStates; ++j)
    {
      if (new_partition[i] == new_partition[j])
      {
        // For two states i,j in the same block, put j into a new block if Steps(i) != Steps(j)
        if (!Steps_Equal(i, j))
        {
          if (!started_new) block_rep.push_back(0);
          
          // Put state into new block
          new_partition[j] = nextblock;
          started_new = true;
          
          // Build vertex sets for all graphs: B -> {Q, ...}
          vertices[partition[j]].insert(new_partition[j]);
        }
      }
    }
    if (started_new) ++nextblock;
    
    block_rep[new_partition[i]] = i;
  }
  
#ifdef QLOG
  fprintf(qlog, "[REFN] Initial refinement: %d refined blocks:\n[REFN]", nextblock);
  for (int n = 0; n < nextblock; ++n)
  {
    fprintf(qlog, " %d = {", n);
    for (int m = 0, i = 0; m < nStates; ++m)
    {
      if (new_partition[m] == n) fprintf(qlog, "%s%d", (i++ == 0 ? "" : ","), m);
    }
    fprintf(qlog, "}");
  }
  fprintf(qlog, "\n");
#endif
  
  // Number of working blocks
  nBlocks = nextblock;
  
  // Create graph to compute SCCs. We create one large graph. This is equivalent
  // to creating a graph for each parent block B since there will be no
  // edges between independent blocks B and B'.
  nextblock = parentblocks;
  digraph = new adjacency_list<vecS, vecS, directedS>(nBlocks);
  
  // Iterate through subblocks Q in each parent block B
  for (i = 0; i < parentblocks; ++i)
  {
    // Iterate over vertices Q in parent block B
    for (mia = vertices[i].begin(); mia != vertices[i].end(); ++mia)
    {
      for (mib = vertices[i].begin(); mib != vertices[i].end(); ++mib)
      {
        if (*mia == *mib) continue;
        
        // If sub-block Q' simulates sub-block Q in the lifted model, add edge (Q,Q') to graph
        if (sim->sRs(block_rep[*mia], block_rep[*mib], true))
        {
          add_edge(vertex(*mia, *digraph), vertex(*mib, *digraph), *digraph);
        }
      }
    }
  }
  block_rep.clear();
  
  // Compute strongly connected components, build parent relation and new partition.
  // With scc_rev we save rebuilding the graph to test Reach(Q,Q'); since we don't know
  // how the partitions are renamed by the algorithm for SCCs, we need a mapping.
  scc.resize(nBlocks, 0);
  scc_rev.resize(nBlocks, 0);
  components = strong_components(*digraph, &scc[0]);
  for (i = 0; i < nStates; ++i)
  {
    scc_rev[scc[new_partition[i]]] = new_partition[i];
    new_partition[i] = scc[new_partition[i]];
    par.insert(std::make_pair(new_partition[i], partition[i]));
  }
  
#ifdef QLOG
  fprintf(qlog, "[REFN] Got %d SCCs in graph (all blocks combined)\n", components);
#endif
  
  // Compute transitive closure
  digraph_closure = new adjacency_list<vecS, vecS, directedS>;
  transitive_closure(*digraph, *digraph_closure);
  delete digraph;
  
  nBlocks = components;
  
  // Sigma_(i+1) == Sigma_(i)
  sigma_unchanged = (nBlocks == parentblocks);

  if (!sigma_unchanged)
  {
    sigma.clear();
    sigma.resize(nBlocks, std::set<int>());
    for (i = 0; i < nStates; ++i) sigma[new_partition[i]].insert(i);
#ifdef QLOG
    fprintf(qlog, "[REFN] %d blocks in refined partition:\n[REFN]", nBlocks);
    for (int n = 0; n < nBlocks; ++n)
    {
      fprintf(qlog, " %d = {", n);
      for (int m = 0, i = 0; m < nStates; ++m)
      {
        if (new_partition[m] == n) fprintf(qlog, "%s%d", (i++ == 0 ? "" : ","), m);
      }
      fprintf(qlog, "}");
    }
    fprintf(qlog, "\n");
#endif
  }
#ifdef QLOG
  else fprintf(qlog, "[REFN] Partition sigma is the same as last iteration\n");
#endif
  
  // Initialize new partition relation (line 8)
  for (i = 0; i < nBlocks; ++i)
  {
    for (j = 0; j < nBlocks; ++j)
    {
      // Identical blocks are always in the relation
      if (i == j)
      {
        rmap.Set(i, j);
        continue;
      }
      
      if (par[i] == par[j])
      {
        // Blocks i,j have the same parent block. Determine if (i,j) in
        // relation by testing if the edge (i,j) is in the transitive closure.
        if (edge(scc_rev[i], scc_rev[j], *digraph_closure).second) rmap.Set(i, j);
        else rmap.Clear(i, j);
      }
      else
      {
        // Blocks i,j have different parent blocks. Determine if (i,j) in
        // relation by checking if (Par(i),Par(j)) was in previous relation
        if (rmap(par[i], par[j])) rmap.Set(i, j);
        else rmap.Clear(i, j);
      }
    }
  }
  
  memcpy(partition, new_partition, sizeof(int) * nStates);
  
#ifdef QLOG
  rmap.ReportCurrent(true);
  fprintf(qlog, "[REFN] New partition relation:");
  for (int n = 0; n < nBlocks; ++n)
  {
    for (int m = 0; m < nBlocks; ++m)
    {
      if (rmap(n, m) && n != m) fprintf(qlog, " (%d,%d)", n, m);
    }
  }
  fprintf(qlog, "\n");
#endif

  delete digraph_closure;
  
  // Reconstruct quotient automaton which is required by lines 9-15 and is also used at the
  // beginning of the next iteration
  LiftDistributions();

  // Purge partition relation based on updated lifted distributions (lines 9-15)
  PurgePartitionRelation();
  
  gamma_unchanged = !rmap.MapChanged();
  
  rmap.Commit();
  
  return (sigma_unchanged && gamma_unchanged ? -1 : nBlocks);
}

// For any model type, test if lifted distribution mu2 simulates mu1 under the given partition relation.
#ifdef OPT_CACHE_NETS
bool StrongSimulation_Quotient::Simulator::muRmu(int mu1, int mu2, bool no_cache)
#else//OPT_CACHE_NETS
bool StrongSimulation_Quotient::Simulator::muRmu(int mu1, int mu2, bool)
#endif//OPT_CACHE_NETS
{
  bool known_result, res;
#ifdef OPT_CACHE_NETS
  // Even if parametric maxflow is enabled, do not cache this network or look it up in the cache
  if (no_cache)
  {
    CompactMaxFlow<double> net;
    
    res = net.CreateNetwork(base->lifted_cols, base->lifted_non_zeros, &base->rmap, base->lifted_row_starts[mu1],
                  base->lifted_row_starts[mu1 + 1] - base->lifted_row_starts[mu1], base->lifted_row_starts[mu2],
                  base->lifted_row_starts[mu2 + 1] - base->lifted_row_starts[mu2], known_result, 0);
    
    if (known_result) return res;
    
    return net.IsFlowTotal();
  }
  
  CompactMaxFlow<double> *cmf;
  std::map<std::pair<int,int>,CompactMaxFlow<double>*>::iterator i;
  
  i = cache.find(std::make_pair(mu1, mu2));
  
  if (i == cache.end() || !i->second) // Network not cached, create
  {
    cmf = new CompactMaxFlow<double>;
    res = cmf->CreateNetwork(base->lifted_cols, base->lifted_non_zeros, &base->rmap, base->lifted_row_starts[mu1],
                  base->lifted_row_starts[mu1 + 1] - base->lifted_row_starts[mu1], base->lifted_row_starts[mu2],
                  base->lifted_row_starts[mu2 + 1] - base->lifted_row_starts[mu2], known_result, 0);
    
    if (!known_result) res = cmf->IsFlowTotal();
    
    if (res)
    {
      cache.insert(std::make_pair(std::make_pair(mu1, mu2), cmf));
#ifdef DEBUG
      ++base->stats.num_nets_cached;
#endif//DEBUG
      return true;
    }
    else
    {
      delete cmf;
      return false;
    }
  }
  else // Network cached, update
  {
#ifdef DEBUG
    ++base->stats.num_cache_hits;
#endif//DEBUG

    cmf = i->second;
    cmf->UpdateNetwork(&base->rmap, false);
    
    if (cmf->IsFlowTotal()) return true;
    
    delete cmf;
    cache[std::make_pair(mu1, mu2)] = 0;
    
    return false;
  }
  
  assert(false); // This point should not be reached
  return false;
#else//OPT_CACHE_NETS
  CompactMaxFlow<double> cmf;
  
  res = cmf.CreateNetwork(base->lifted_cols, base->lifted_non_zeros, &base->rmap, base->lifted_row_starts[mu1],
                 base->lifted_row_starts[mu1 + 1] - base->lifted_row_starts[mu1], base->lifted_row_starts[mu2],
                 base->lifted_row_starts[mu2 + 1] - base->lifted_row_starts[mu2], known_result, 0);
  
  if (known_result) return res;
  
  return cmf.IsFlowTotal();
#endif//OPT_CACHE_NETS
}

void StrongSimulation_Quotient::Simulator::Flush()
{
#ifdef OPT_CACHE_NETS
  std::map<std::pair<int,int>,CompactMaxFlow<double>*>::iterator i;
  
  for (i = cache.begin(); i != cache.end(); ++i) if (i->second) delete i->second;
  cache.clear();
#endif//OPT_CACHE_NETS
}

// For Markov chains, test if state s2 can simulate s1. Since every state has exactly one
// distribution in a Markov chain, this is equivalent to muRmu(s1,s2).
bool StrongSimulation_Quotient::Simulator_MC::sRs(int s1, int s2, bool no_cache)
{
  return muRmu(s1, s2, no_cache);
}

// Test if there exists an s2 -> mu2 in q2 which can simulate mu1 such that s1 -> mu1
// for all s1 in q1 (true is returned if not all states s1 in q1 have the same
// distribution mu1)
bool StrongSimulation_Quotient::Simulator_MC::qRq(int q1, int q2, bool no_cache)
{
  int pi;
  set<int>::iterator i;
  multimap<int,pair<int,int> >::iterator mi;
  
  mi = base->forall.find(q1);
  if (mi == base->forall.end()) return true;
  else pi = mi->second.second;
  
  for (i = base->sigma[q2].begin(); i != base->sigma[q2].end(); ++i)
  {
    if (muRmu(pi, *i, no_cache)) return true;
  }
  
  return false;
}

// For PAs, for every action of s1 test if there exists an action of s2 such that the
// corresponding distributions simulate under muRmu().
bool StrongSimulation_Quotient::Simulator_PA::sRs(int s1, int s2, bool no_cache)
{
  int i, j;
  for (i = base->state_starts[s1]; i < base->state_starts[s1 + 1]; ++i)
  {
    for (j = base->state_starts[s2]; j < base->state_starts[s2 + 1]; ++j)
    {
      if (base->actions[i] == base->actions[j] && muRmu(i, j, no_cache)) break;
    }
    
    if (j == base->state_starts[s2 + 1]) return false;
  }
  
  return true;
}

// For every distribution that all states in q1 have in common, see if there is a distribution
// in q2 that can simulate it
bool StrongSimulation_Quotient::Simulator_PA::qRq(int q1, int q2, bool no_cache)
{
  pair<multimap<int,pair<int,int> >::iterator,multimap<int,pair<int,int> >::iterator> mu;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> pi;
  multimap<int,int> needed;
  multimap<int,int>::iterator mmi;
  set<int>::iterator i;
  int j;
  
  // Get \forall-distributions of q1 and put them in a set, return true if there are none
  mu = base->forall.equal_range(q1);
  if (mu.first == mu.second) return true;
  for (; mu.first != mu.second; ++mu.first) needed.insert(mu.first->second);
  
  // Iterate over all states *i in q2
  for (i = base->sigma[q2].begin(); i != base->sigma[q2].end(); ++i)
  {
    // Iterate over all distributions j of each state *i in q2
    for (j = base->state_starts[*i]; j < base->state_starts[*i + 1]; ++j)
    {
      // For all \forall distributions of q1 with the same action as j, check if they are
      // simulated by j. If yes, remove them from the "needed" set. Return true if the set
      // is empty.
      pi = needed.equal_range(base->actions[j]);
      while (pi.first != pi.second)
      {
        if (muRmu(pi.first->second, j, no_cache))
        {
          mmi = pi.first;
          ++pi.first;
          needed.erase(mmi);
          if (needed.size() == 0) return true;
        }
        else ++pi.first;
      }
    }
  }
  
#ifdef QLOG
  fprintf(base->qlog, "[    ] (%d,%d) not in Gamma because of \\forall dists:", q1, q2);
  for (mmi = needed.begin(); mmi != needed.end(); ++mmi) fprintf(base->qlog, " %d", mmi->second);
  fprintf(base->qlog, "\n");
#endif
  
  return false;
}
