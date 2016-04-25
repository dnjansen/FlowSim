/*****************************************************************************/
/*!
 *   Copyright 2014-2015 David N. Jansen
 *
 *   This file is part of FlowSim.
 *
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

// The simulation algorithm in this file follows closely the publication:
// Silvia Crafa; Francesco Ranzato: Bisimulation and simulation algorithms
// on probabilistic transition systems by abstract interpretation.
// Formal Methods in System Design (2012) 40:356-376.
// DOI 10.1007/s10703-012-0147-3


#include "StrongCR.h"

using namespace std;

#define OPTIMIZE_MEMORY
#undef REPORT_MEMORY

#ifdef REPORT_MEMORY
  static clock_t clkstart;
  ssize_t old_mem_used, mem_remove, mem_deleted, mem_listener;
#endif

// Compute the simulation relation
unsigned int StrongSimulation_CR::Simulate(ProbabilisticModel *target,
            std::set<std::pair<int,int> > *result)
{
  int size;
  
  // Assert that we have valid pointers
  assert(NULL != target);
  assert(NULL != label_func);
  
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
    row_starts = dynamic_cast<MarkovChain*>(model)->row_starts;
    cols = dynamic_cast<MarkovChain*>(model)->cols;
    non_zeros = dynamic_cast<MarkovChain*>(model)->non_zeros;
    state_starts = NULL;
    actions = NULL;
    break;
  // Probabilistic Automaton
  case ProbabilisticModel::PA:
    row_starts = dynamic_cast<ProbabilisticAutomaton*>(model)->row_starts;
    cols = dynamic_cast<ProbabilisticAutomaton*>(model)->cols;
    non_zeros = dynamic_cast<ProbabilisticAutomaton*>(model)->non_zeros;
    state_starts = dynamic_cast<ProbabilisticAutomaton*>(model)->state_starts;
    actions = dynamic_cast<ProbabilisticAutomaton*>(model)->actions;
    break;
  default:
    assert(false);
  }

  // generate quotient automaton
  #ifdef DEBUG
    stats.ResetStats();
    CompactMaxFlow<double>::ResetStats();
    // stats.num_initial_pairs = nBlocks;
  #endif//DEBUG

  // Algorithm of Crafa / Ranzato
  PBis();

  // Count pairs in the actual relation and place result in the provided set
  size = 0;
  if (result) result->clear();
  for (int n = 0; n < nStates; ++n)
  {
    for (int m = 0; m < nStates; ++m)
    {
      if (n != m && rmap(n, m) )
      {
        ++size;
        if (result) result->insert(std::make_pair(n, m));
      }
    }
  }

  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"after creating result set\"\n",
                (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener);
  #endif
  for (int d = 0; d < nDistributions; ++d) {
    for (int a = 0; a < nActions; ++a) {
      if (NULL != distributions[d].in[a].Count) {
        delete [] distributions[d].in[a].Count;
        RegisterMemFree(sizeof(*distributions[d].in[a].Count) * nStates);
        distributions[d].in[a].Count = NULL;
      }
    }
    delete [] distributions[d].in;
    RegisterMemFree(sizeof(*distributions[d].in) * nActions);
    distributions[d].in = NULL;
  }
  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"after freeing Count[,,] and "
                "in[]\"\n", (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener);
  #endif
  delete [] distributions;          distributions = NULL;
  RegisterMemFree(sizeof(*distributions) * nDistributions);
  #ifdef REPORT_MEMORY
    fprintf(stderr,"%lu,%zd,%zu,%zd,%zd,%zd,\"after freeing distributions\"\n",
                (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener);
  #endif
  if (action_masks) {
    delete [] action_masks;         action_masks = NULL;
    RegisterMemFree(nStates * action_mask_pitch * sizeof(*action_masks));
    #ifdef REPORT_MEMORY
      fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"after freeing mark[,]\"\n",
                (unsigned long)(clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener);
    #endif
  }
  assert (NULL == Deleted);
  for (int n = nStates * nStates; --n >= 0; ) {
    Listener_t *temp = Listener[n];
    Listener[n] = NULL;
    while (NULL != temp) {
      Listener_t *next = temp->next;
      delete temp;
      #ifdef REPORT_MEMORY
        mem_listener -= sizeof(*temp);
      #endif
      RegisterMemFree(sizeof(*temp));
      temp = next;
    }
    #ifndef OPTIMIZE_MEMORY
      if (0 == n % (nStates * nStates / 16))
        fprintf(stderr,"%lu,%zd,%zu,%zd,%zd,%zd,\"during freeing "
                "Listener[,]\"\n", (unsigned long)(clock() - clkstart),
                CurMemUsed(), CompactMaxFlow<double>::global_space, mem_remove,
                mem_deleted, mem_listener);
    #endif
  }
  delete [] Listener;
  #ifdef REPORT_MEMORY
    mem_listener -= sizeof(*Listener) * nStates * nStates;
  #endif
  RegisterMemFree(sizeof(*Listener) * nStates * nStates);
  Listener = NULL;
  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"after freeing Listener[,]\","
                "%zd,\"is the total memory peak\"\n",
                (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener, MaxMemUsed());
  #endif

  #ifdef DEBUG
  rmap.CollectStats(&stats);
  rmap.clear_mem();
  Rdmap.CollectStats(&stats);
  Rdmap.clear_mem();
  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"after freeing R and Rd\"\n",
                (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener);
  #endif
  CompactMaxFlow<double>::CollectStats(&stats);
  stats.CollectStats();
  // stats.mem_partition_map = 0;
  // stats.mem_model = 0;
  // stats.num_partitions = nBlocks;
  stats.num_final_pairs = size;
  #endif//DEBUG

  return size;
}

// Figure 7 of Crafa/Ranzato
void StrongSimulation_CR::Initialize()
{
  // Preparation: create an array of distinct distributions

  #ifdef REPORT_MEMORY
    old_mem_used = CurMemUsed();
  #endif
  RegisterMemAlloc(sizeof(*distr_index) * model->Distributions());
  distr_index = new int[model->Distributions()];
  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"before finding distributions\","
                "%zd,\"is the size of the distribution index\"\n",
                (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener, CurMemUsed() - old_mem_used);
    old_mem_used = CurMemUsed();
  #endif
  // possibly less memory is used, but let's allocate for the worst case
  // we will call RegisterMemAlloc() when the actual memory is used.
  distributions = new Distribution[model->Distributions()];

  switch (model->Type())
  {
  case ProbabilisticModel::MC:
    nActions = 1;
    break;
  case ProbabilisticModel::PA:
    nActions = dynamic_cast<ProbabilisticAutomaton*>(model)->da;
    break;
  }
  // Unify distributions that are identical. This only works well if the
  // distributions are stored in the same order, for example if all column
  // numbers in one row are ascending. (Note that ...::Parse checks for this.)
  nDistributions = 0;
  #ifdef REPORT_MEMORY
    int p = 0;
  #endif
  for (int n = 0; n < model->Distributions(); ++n) {
    const int startn = row_starts[n], sizen = row_starts[n+1] - startn;
    for (int m = 0; ; ++m) {
      if (m >= nDistributions) {
        // The distribution n is new.
        RegisterMemAlloc(sizeof(Distribution) + sizeof(in_t) * nActions);
        // stats.mem_model += (row_starts[n+1] - row_starts[n])
        //             * (sizeof(int) + sizeof(double));
        distributions[nDistributions].row_starts = startn;
        distributions[nDistributions].row_ends = row_starts[n+1];
        distributions[nDistributions].in = new in_t[nActions];
        distr_index[n] = nDistributions++;
        #ifdef REPORT_MEMORY
          p += sizen;
        #endif
        break; // inner loop
      }
      const int startm = distributions[m].row_starts;
      if (distributions[m].row_ends - startm == sizen
                  && 0 == memcmp(&cols[startn], &cols[startm],
                              sizen * sizeof(*cols))
                  && 0 == memcmp(&non_zeros[startn], &non_zeros[startm],
                              sizen * sizeof(*non_zeros)))
      {
        // The distribution n is identical to distribution m.
        distr_index[n] = m;
        break; // inner loop
      }
    }
  }
  // nDistributions is now the number of distinct distributions.
  #ifdef REPORT_MEMORY
    fprintf(stderr, "\nunique_|Distr|    %8d\n"
                      "unique_p          %8d\n"
                      "|P|_=_m           %8d\n"
                      "number_of_arrows  %8d\n",
                nDistributions, p, model->Distributions(),
                row_starts[model->Distributions()]);
  #endif
  // distributions[] is an array of distinct distributions.
  // distr_index[] indicates, for each successor distribution (as numbered in
  // row_starts[] and row_ends[]), where in the array distributions[] it is
  // located.
  #ifdef REPORT_MEMORY
    fprintf(stderr,"%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 7, before line 2\",%zd,"
                "\"is the size of in[]\"\n",
                (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener, CurMemUsed() - old_mem_used);
    old_mem_used = CurMemUsed();
  #endif

  // Now create the in and predecessor sets
  if (model->Type() != ProbabilisticModel::PA) {
    // Lines 2-4: find for each distribution its in-actions.
    // This is trivial: each distribution has the only action as its in-action.
    // Still, we need to find the predecessors for the in-action.
    RegisterMemAlloc(nStates * sizeof(*distributions->in->pre));
    for (int n = nStates; --n >= 0; ) {
      // This loop counts backwards so that the list pre in the end will be
      // ordered the right way.
      int_list * &pre = distributions[distr_index[n]].in->pre;
      // The test in the following line is superfluous: each n will be inserted
      // exactly once into some pre-list.
      // if (NULL == pre || pre->el != n) {
        int_list *new_el = new int_list;
        new_el->el = n;
        new_el->next = pre;
        pre = new_el;
      // }
    }
  } else {
    // Lines 2-4: find for each distribution its in-actions.
    // Note that we do the calculation in a different order because we do not
    // have easy access to the in-sets of distributions. Therefore, we have to
    // run through all states and add them to the in-sets of their successors.
    // This needs time complexity m, which fits in the time complexity.
    RegisterMemAlloc(state_starts[nStates] * sizeof(*distributions->in->pre));
    for (int n = nStates; --n >= 0; ) {
      // This loop counts backwards so that the list pre in the end will be
      // ordered the right way.
      for (int m = state_starts[n]; m < state_starts[n+1]; ++m) {
        int_list * &pre = distributions[distr_index[m]].in[actions[m]].pre;
        // The test in the following line is *almost* superfluous: n will be
        // inserted into the same pre-list twice only if n has two actions with
        // the same label and the same distribution.
        // We suppress the test anyway because it makes RegisterMemAlloc()
        // simpler (see the above call, which registers all ``new int_list''
        // allocations at once).
        // if (NULL == pre || pre->el != n) {
          int_list *new_el = new int_list;
          new_el->el = n;
          new_el->next = pre;
          pre = new_el;
        // }
      }
    }
  }
  delete [] distr_index;
  RegisterMemFree(model->Distributions() * sizeof(*distr_index));
  distr_index = NULL;
  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 7, before line 5\",%zd,"
                "\"is the size of pre[] - the size of the distribution "
                "index.\"\n", (unsigned long) (clock()-clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener, CurMemUsed() - old_mem_used);
  #endif
  // Lines 5-12: Initialize Rs based on labels and actions allowed in each
  // state.
  InitializeActionMasks(); // Lines 5-8
  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 7, before line 9\"\n",
                (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener);
  #endif
  rmap.ReportCurrent(true);
  rmap.Create(nStates, 0);
  #ifdef REPORT_MEMORY
    fprintf(stderr,"%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 7, before line 9 (bis)\","
                "%zd,\"is the size of rmap\"\n",
                (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener, CurMemUsed() - old_mem_used);
  #endif

/*
  // This code would be a simpler translation of lines 9-12 than the two loops
  // below.  However, it has time complexity n^2*|Act|, which does not fit into
  // the theoretical bound.  The execution time stays the same.
  for (int n = 0; n < nStates; ++n) {
    for (int m = 0; m < nStates; ++m) {
      if (n == m || Label(n) != Label(m))
        continue;
      if (Act_Subseteq(m, n)) {
        rmap.Set(m, n);
        #ifdef DEBUG
          stats.num_initial_pairs++;
        #endif
      }
      // else rmap.Clear(m, n); -- superfluous because it's initialized to 0
    }
  }
*/
  // The following seems to be forgotten in the algorithm of Fig. 7.
  // It has time complexity n^2 and fits into the stated complexity.
  for (int n = 0; n < nStates; ++n) {
    for (int m = 0; m < nStates; ++m) {
      if (n == m || Label(n) != Label(m))
        continue;
      rmap.Set(m, n);
      #ifdef DEBUG
        stats.num_initial_pairs++;
      #endif
    }
  }

  // Lines 9-12 of Figure 7:
  if (NULL != action_masks) {
    for (int d = 0; d < nDistributions; ++d) {
      for (int a = 0; a < nActions; ++a) {
        for (int_list *x = distributions[d].in[a].pre; NULL != x; x = x->next){
          for (int y = 0; y < nStates; ++y) {
            if (!(action_masks[y * action_mask_pitch
                                        + a / (sizeof(*action_masks)*CHAR_BIT)]
                                & 1 << a % (sizeof(*action_masks)*CHAR_BIT)))
            {
              #ifdef DEBUG
                if (rmap(x->el, y)) {
                  rmap.Clear(x->el, y);
                  stats.num_initial_pairs--;
                }
              #else
                rmap.Clear(x->el, y);
              #endif
            }
          }
        }
      }
    }

    delete [] action_masks;
    RegisterMemFree(nStates * action_mask_pitch * sizeof(*action_masks));
    action_masks = NULL;
  }
  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 7, before line 13\",,"
                "\"I have just deallocated mark[,].\"\n",
                (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener);
    old_mem_used = CurMemUsed(); // after deallocation of the action masks
  #endif

  // Line 13: Initialize Rd_inv
  Rdmap.ReportCurrent(true);
  Rdmap.Create(nDistributions);
  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 7, before line 13 (bis)\","
                "%zd,\"is the size of Rdmap\"\n",
                (unsigned long)(clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener, CurMemUsed() - old_mem_used);
    old_mem_used = CurMemUsed();
  #endif
  for (int d = 0; d < nDistributions; ++d) {
    for (int e = 0; e < nDistributions; ++e) {
      if (d == e)
        continue;
      #ifdef OPTIMIZE_MEMORY
        // Only distributions that have at least one predecessor action in
        // common are relevant to be compared. This saves about half of the
        // memory in dining cryptographers example crypt4.
        // This test does, however, increase the theoretical time complexity to
        // |Distr|^2 * |Act|.
        // (We can probably implement the test in complexity |Distr| * m, which
        // would be ok.)
        int a;
        for (a = 0; a < nActions; ++a) {
          if (NULL != distributions[d].in[a].pre
                                        && NULL != distributions[e].in[a].pre)
            break;
        }
        if (!(a < nActions)) continue;
      #endif
      CompactMaxFlow<double>* cmf = Init_SMF(d, e);
      if (NULL != cmf) {
        RegisterMemAlloc(sizeof(std::pair<int,CompactMaxFlow<double>*>)
                  + MAP_OVERHEAD);
        distributions[e].Rd_inv[d] = cmf;
        Rdmap.Set(d, e);
      }
    }
    #ifdef REPORT_MEMORY
      if (0 == d % (nDistributions / 16)) {
        fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 7, in line 13\"\n",
                (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener);
      }
    #endif
  }
  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 7, before line 14\",%zd,"
                "\"is the size of Rd_inv + the relevant flow networks\"\n",
                (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener, CurMemUsed() - old_mem_used);
    old_mem_used = CurMemUsed();
  #endif
  // Lines 14-21: Initialize Count[] to zero.
  // This has time and space complexity m * n. The space complexity may or may
  // not be larger than the one stated in Crafa/Ranzato (m * |Distr|) because
  // we simplified the data structure, but it still fits into the overall space
  // complexity.
  for (int d = 0; d < nDistributions; ++d) {
    for (int a = 0; a < nActions; ++a) {
      if (NULL != distributions[d].in[a].pre) {
        RegisterMemAlloc(sizeof(*distributions[d].in[a].Count) * nStates);
        distributions[d].in[a].Count = new int[nStates];
        memset(distributions[d].in[a].Count, '\0', sizeof(int) * nStates);
      }
    }
  }
  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 7, before line 18\",%zd,"
                "\"is the size of Count[,,]\"\n",
                (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener, CurMemUsed() - old_mem_used);
    old_mem_used = CurMemUsed();
  #endif
  // Now we set up the actual counters.
  for (int d = 0; d < nDistributions; ++d) {
    for (int a = 0; a < nActions; ++a) {
      for (int_list * x = distributions[d].in[a].pre; NULL != x; x = x->next)
      {
        // Because d is not in distributions[d].Rd_inv, we have to add it
        // separately:
        distributions[d].in[a].Count[x->el]++;
        for (std::map<int,CompactMaxFlow<double>*>::const_iterator e
                    = distributions[d].Rd_inv.begin();
                    distributions[d].Rd_inv.end() != e; ++e)
        {
          if (NULL != distributions[e->first].in[a].pre) {
            distributions[e->first].in[a].Count[x->el]++;
          }
        }
      }
    }
  }
  // Initialize Listener (lines 30-32)
  // We do this before initialising Remove_a(d) because Listener needs the
  // rmap as it was before removing any pairs from the relation.
  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 7, before line 30\"\n",
                (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener);
  #endif
  assert(NULL == Listener);
  #ifdef REPORT_MEMORY
    mem_listener += sizeof(*Listener) * nStates * nStates;
  #endif
  RegisterMemAlloc(sizeof(*Listener) * nStates * nStates);
  Listener = new Listener_t *[nStates * nStates];
  for (int e = 0; e < nDistributions; ++e) {
    #ifdef OPTIMIZE_MEMORY
    // we deviate from the published algorithm here:
    // we only add those pairs of distributions that are in Rd.
    for (std::map<int,CompactMaxFlow<double>*>::const_iterator dd
                = distributions[e].Rd_inv.begin();
                distributions[e].Rd_inv.end() != dd; ++dd)
    {
      int d = dd->first;
    #else
    for (int d = 0; d < nDistributions; ++d) {
      if (d == e)
        continue;
    #endif
      for (int n = distributions[d].row_starts;
                  n < distributions[d].row_ends; ++n)
      {
        for (int m = distributions[e].row_starts;
                    m < distributions[e].row_ends; ++m)
        {
          if (cols[n] == cols[m])
            continue;
          #ifdef OPTIMIZE_MEMORY
            // we deviate from the published algorithm here:
            // we only add those pairs of states that are in R.
            if (!rmap(cols[n], cols[m]))
              continue;
          #endif
          #ifdef REPORT_MEMORY
            mem_listener += sizeof(Listener_t);
          #endif
          RegisterMemAlloc(sizeof(Listener_t));
          Listener_t * new_el = new Listener_t;
          new_el->d = d;
          new_el->e = e;
          new_el->next = Listener[cols[n] * nStates + cols[m]];
          Listener[cols[n] * nStates + cols[m]] = new_el;
        }
      }
    }
    #ifdef REPORT_MEMORY
      if (0 == e % (nDistributions / 8)) {
        fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 7, in line 31\"\n",
                (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener);
      }
    #endif
  }
  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 7, before line 22\",%zd,"
                "\"is the size of Listener[,]\"\n",
                (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener, CurMemUsed() - old_mem_used);
    old_mem_used = CurMemUsed();
  #endif
  // Initialize Remove_a(d) (lines 22-26)
  #ifdef OPTIMIZE_MEMORY
    // Because the Remove lists need an awful lot of memory, we do not
    // calculate them all before we use them. In this loop, we combine the
    // initialization of Remove with preStabilize().
    probStable = true;
    // Initialize Deleted (line 33)
    Deleted = NULL;
    // Lines 22-26: Initialize Remove
    for (int d = 0; d < nDistributions; ++d) {
      for (int a = 0; a < nActions; ++a) {
        if (NULL == distributions[d].in[a].pre) continue;
        int_list *remove_a_d = NULL;
        for (int x = 0; x < nStates; ++x) {
          if (0 == distributions[d].in[a].Count[x]) {
            // mem_remove += sizeof(*remove_a_d);
            RegisterMemAlloc(sizeof(*remove_a_d));
            int_list *new_el = new int_list;
            new_el->el = x;
            new_el->next = remove_a_d;
            remove_a_d = new_el;
          }
        }
        // Here follows some code from preStabilize():
        if (NULL != remove_a_d) {
          // line 5 of preStabilize(): for each a-predecessor of d
          for (int_list *t = distributions[d].in[a].pre; NULL != t; t=t->next)
          {
            // line 6 of preStabilize(): for each w in Remove_a(d)
            int_list *remove = remove_a_d;
            do {
              int w = remove->el;
              // line 7 of preStabilize():
              if (/* t->el == w || */ rmap(t->el, w)) {
                rmap.Clear(t->el, w);
                // The pair (t->el, w) cannot be an element of Deleted already:
                // as soon as it is added to Deleted, it is removed from rmap.
                // So we add it to the list of deleted states without further
                // check.
                #ifdef REPORT_MEMORY
                  mem_deleted += sizeof(*Deleted);
                #endif
                RegisterMemAlloc(sizeof(*Deleted));
                Deleted_t *new_el = new Deleted_t;
                new_el->t = t->el; new_el->w = w;
                new_el->next = Deleted;
                Deleted = new_el;
                probStable = false;
              }
              remove = remove->next;
            } while (NULL != remove);
          }
          // line 4 of preStabilize():
          // distributions[d].in[a].Remove = NULL;
          do {
            int_list *temp = remove_a_d;
            remove_a_d = remove_a_d->next;
            delete temp;
            // mem_remove -= sizeof(*temp);
            RegisterMemFree(sizeof(*temp));
          } while (NULL != remove_a_d);
        }
        // Here the code from preStabilize() ends.
        assert(NULL == distributions[d].in[a].Remove);
      }
      #ifdef REPORT_MEMORY
        if (0 == d % (nDistributions / 8)) {
          fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig 7, in line 22\"\n",
                (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener);
        }
      #endif
    }
    // Initialize stability flags
    // probStable = ...; -- see above
    preStable = true;
    // if (probStable)
    //   return;
  #else // OPTIMIZE_MEMORY
    preStable = true;
    // Lines 22-26: Initialize Remove
    for (int d = 0; d < nDistributions; ++d) {
      for (int a = 0; a < nActions; ++a) {
        if (NULL == distributions[d].in[a].pre)
          continue;
        assert(NULL == distributions[d].in[a].Remove);
        for (int x = 0; x < nStates; ++x) {
          if (0 == distributions[d].in[a].Count[x]) {
            #ifdef REPORT_MEMORY
              mem_remove += sizeof(*distributions[d].in[a].Remove);
            #endif
            RegisterMemAlloc(sizeof(*distributions[d].in[a].Remove));
            int_list *new_el = new int_list;
            new_el->el = x;
            new_el->next = distributions[d].in[a].Remove;
            distributions[d].in[a].Remove = new_el;
            preStable = false;
          }
        }
      }
      #ifdef REPORT_MEMORY
        if (0 == d % (nDistributions / 8)) {
          fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 7, in line 22 "
                "(bis)\"\n", (unsigned long)(clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener);
        }
      #endif
    }
    // Initialize stability flags
    probStable = true;
    // preStable = ...; -- see above
    // Initialize Deleted (line 33)
    Deleted = NULL;
    // if (preStable)
    //   return;
  #endif // OPTIMIZE_MEMORY
  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 7, at end\",%zd,\"is the "
                "size of Deleted[,]\"\n", (unsigned long) (clock() - clkstart),
                CurMemUsed(), CompactMaxFlow<double>::global_space, mem_remove,
                mem_deleted, mem_listener, CurMemUsed() - old_mem_used);
    old_mem_used = CurMemUsed();
  #endif
  return;
}

// Compute action masks. The action mask is a bitset where every bit
// corresponds to a certain action. If a state has at least one non-zero
// distribution for a certain action, the bit is on, otherwise it is off. This
// allows us to perform the test Act(s)=Act(s') in O(l) time, where l is the
// number of words required to represent the mask.
void StrongSimulation_CR::InitializeActionMasks()
{
  if (model->Type() != ProbabilisticModel::PA)
  {
    action_masks = NULL;
    action_mask_pitch = 1;
    return;
  }
  
  action_mask_pitch = (sizeof(*action_masks)*CHAR_BIT - 1 + nActions)
              / (sizeof(*action_masks)*CHAR_BIT);
  RegisterMemAlloc(nStates * action_mask_pitch * sizeof(*action_masks));
  // Line 5 of Fig. 7:
  action_masks = new unsigned int[nStates * action_mask_pitch];
  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 7, before line 6\",%zd,"
                "\"is the size of mark[,] (will be deallocated before the "
                "peak)\"\n", (unsigned long) (clock()-clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener, CurMemUsed() - old_mem_used);
    old_mem_used = CurMemUsed();
  #endif
  memset(action_masks, '\0', nStates*action_mask_pitch*sizeof(*action_masks));

/* Lines 6-8 of Fig. 7 would be, translated more literally:
  for (int d = 0; d < nDistributions; ++d) {
    for (int a = 0; a < nActions; ++a) {
      for (int_list *p = distributions[d].in[a].pre; NULL != p; p = p->next)
      {
        action_masks[p->el * action_mask_pitch
                                        + a / (sizeof(*action_masks)*CHAR_BIT)]
                  |= 1 << a % (sizeof(*action_masks)*CHAR_BIT);
      }
    }
  }
The following code has time complexity n*|Act|, which fits into the stated time
complexity, as |Act| <= m.
*/
  for (int i = 0; i < nStates; ++i)
  {
    for (int j = state_starts[i], a = -1; j < state_starts[i+1]; ++j)
    {
      if (actions[j] == a) continue;
      a = actions[j];
      action_masks[i*action_mask_pitch + a/(sizeof(*action_masks)*CHAR_BIT)]
                  |= 1 << a % (sizeof(*action_masks)*CHAR_BIT);
    }
  }
}


// Figure 6 in Crafa/Ranzato: preStabilize and probStabilize
bool StrongSimulation_CR::preStabilize() {
  // line 2:
  assert(NULL == Deleted);
  // line 3: search for a nonempty set of removed states
  for (int e = 0; e < nDistributions; ++e) {
    for (int a = 0; a < nActions; ++a) {
      if (NULL != distributions[e].in[a].Remove) {
        // line 5: for each a-predecessor of e
        for (int_list *t = distributions[e].in[a].pre; NULL != t; t = t->next)
        {
          // line 6: for each w in Remove_a(e)
          int_list *remove = distributions[e].in[a].Remove;
          do {
            int w = remove->el;
            // line 7:
            if (/* t->el == w || */ rmap(t->el, w)) {
              rmap.Clear(t->el, w);
              // The pair (t->el, w) cannot be an element of Deleted already:
              // as soon as it is added to Deleted, it is removed from rmap. So
              // we add it to the list of deleted states without further check.
              #ifdef REPORT_MEMORY
                mem_deleted += sizeof(*Deleted);
              #endif
              RegisterMemAlloc(sizeof(*Deleted));
              Deleted_t *new_el = new Deleted_t;
              new_el->t = t->el; new_el->w = w;
              new_el->next = Deleted;
              Deleted = new_el;
              probStable = false;
            }
            remove = remove->next;
          } while (NULL != remove);
        }
        // line 4:
        int_list *remove = distributions[e].in[a].Remove;
        distributions[e].in[a].Remove = NULL;
        do {
          int_list *temp = remove;
          remove = remove->next;
          delete temp;
          #ifdef REPORT_MEMORY
            mem_remove -= sizeof(*temp);
          #endif
          RegisterMemFree(sizeof(*temp));
        } while (NULL != remove);
      }
    }
    #ifdef REPORT_MEMORY
      if (0 == e % (nDistributions / 8)) {
        fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 6, in line 3 of "
                "preStabilize()\"\n", (unsigned long) (clock() - clkstart),
                CurMemUsed(), CompactMaxFlow<double>::global_space, mem_remove,
                mem_deleted, mem_listener);
      }
    #endif
  }
  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 6, before line 8 of "
                "preStabilize()\"\n", (unsigned long) (clock() - clkstart),
                CurMemUsed(), CompactMaxFlow<double>::global_space, mem_remove,
                mem_deleted, mem_listener);
  #endif
  return probStable;
}


bool StrongSimulation_CR::probStabilize() {
  #ifdef REPORT_MEMORY
    int counter = 500;
  #endif
  // line 2: for all edges in Deleted
  while (NULL != Deleted) {
    // *tw is a pair of states:
    int t = Deleted->t, w = Deleted->w;
    {
      Deleted_t *temp = Deleted;
      Deleted = Deleted->next;
      delete temp;
      #ifdef REPORT_MEMORY
        mem_deleted -= sizeof(*temp);
      #endif
      RegisterMemFree(sizeof(*temp));
    }
    // line 3: for all distribution pairs in the Listener of tw
    Listener_t *de = Listener[t*nStates + w];
    Listener[t*nStates + w] = NULL;
    while (NULL != de) {
      int d = de->d, e = de->e;
      // lines 4+5:
      if (Rdmap(d, e) && !SMF(d, e, t, w)) {
        // (It can happen that (d Rd e) has been deleted based on some other
        // Listener -- the above test ``Rdmap(d, e)'' is really required.)

        // distributions[e].Rd_inv.erase(d);--this is already done in SMF().
        // line 6:
        for (int b = 0; b < nActions; ++b) {
          if (NULL == distributions[d].in[b].pre)
            continue;
          // line 7:
          for (int_list *s = distributions[e].in[b].pre; NULL != s; s=s->next)
          {
            // lines 8 + 9:
            if (--distributions[d].in[b].Count[s->el] <= 0) {
              // line 10:
              #ifdef REPORT_MEMORY
                mem_remove += sizeof(*distributions[d].in[b].Remove);
              #endif
              RegisterMemAlloc(sizeof(*distributions[d].in[b].Remove));
              int_list *new_el = new int_list;
              new_el->el = s->el;
              new_el->next = distributions[d].in[b].Remove;
              distributions[d].in[b].Remove = new_el;
              preStable = false;
            }
          }
        }
      }
      // we can delete this element from the Listener-list.
      Listener_t *temp = de;
      de = de->next;
      delete temp;
      #ifdef REPORT_MEMORY
        mem_listener -= sizeof(*temp);
      #endif
      RegisterMemFree(sizeof(*temp));
    }
    #ifdef REPORT_MEMORY
      if (--counter <= 0) {
        fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 6, in line 3 of "
                "probStabilize()\"\n", (unsigned long) (clock() - clkstart),
                CurMemUsed(), CompactMaxFlow<double>::global_space, mem_remove,
                mem_deleted, mem_listener);
        counter = state_starts[nStates] * 4;
      }
    #endif
  }
  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 6, before line 11 of "
                "probStabilize()\"\n", (unsigned long) (clock() - clkstart),
                CurMemUsed(), CompactMaxFlow<double>::global_space, mem_remove,
                mem_deleted, mem_listener);
  #endif
  return preStable;
}

// Figure 3 in Crafa/Ranzato: PBis
void StrongSimulation_CR::PBis()
{
  #ifdef REPORT_MEMORY
    fprintf(stderr,"\n\"time (1/%lu sec)\",\"total memory\",\"maxflow memory\""
                ",\"Remove memory\",\"Deleted memory\",\"Listener memory\"\n0,"
                "%zd,%zu,%zd,%zd,%zd,\"Fig. 4, before line 2\"\n",
                (unsigned long) CLOCKS_PER_SEC, CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener);
    mem_remove = 0; mem_deleted = 0; mem_listener = 0;
    clkstart = clock();
  #endif
  Initialize();
  #ifdef DEBUG
    stats.num_iterations = 0;
  #endif
  while (!(preStable && probStable)) {
    if (!preStable) {
      probStable = preStabilize();
      preStable = true;
    }
    if (!probStable) {
      preStable = probStabilize();
      probStable = true;
    }
    #ifdef DEBUG
      stats.num_iterations++;
    #endif
  }
  #ifdef REPORT_MEMORY
    fprintf(stderr, "%lu,%zd,%zu,%zd,%zd,%zd,\"Fig. 4, before line 6\"\n",
                (unsigned long) (clock() - clkstart), CurMemUsed(),
                CompactMaxFlow<double>::global_space, mem_remove, mem_deleted,
                mem_listener);
  #endif
}

// Additional procedures
CompactMaxFlow<double>* StrongSimulation_CR::Init_SMF(int d, int e)
// The function creates an initial flow network to verify whether e can
// simulate d.  If yes, the network is stored in the cache.
{
  bool known_result, result;
  CompactMaxFlow<double>* cmf;
  CompactMaxFlow<double>::RegisterMemAlloc(sizeof(*cmf));
  cmf = new CompactMaxFlow<double>();

  result = cmf->CreateNetwork(cols, non_zeros, &rmap,
              distributions[d].row_starts,
              distributions[d].row_ends - distributions[d].row_starts,
              distributions[e].row_starts,
              distributions[e].row_ends - distributions[e].row_starts,
              known_result, 0);
  if (!known_result) {
    result = cmf->IsFlowTotal(false);
  }
  if (result) {
    #ifdef DEBUG
      stats.num_nets_cached++;
    #endif
    return cmf;
  } else {
    delete cmf;
    CompactMaxFlow<double>::RegisterMemFree(sizeof(*cmf));
    return NULL;
  }
}

bool StrongSimulation_CR::SMF(int d, int e, int t, int w)
// The function checks whether the flow network between d and e still holds if
// edge tw is deleted from it.
// If the function returns false, it also deletes the pair (d,e) from the
// relation Rd.
{
  CompactMaxFlow<double>* cmf = distributions[e].Rd_inv[d];

  cmf->DeleteArc(t, w);
  #ifdef DEBUG
    stats.num_cache_hits++;
  #endif
  if (!cmf->IsFlowTotal(false)) {
    delete cmf;
    CompactMaxFlow<double>::RegisterMemFree(sizeof(*cmf));
    distributions[e].Rd_inv.erase(d);
    RegisterMemFree(sizeof(std::pair<int,CompactMaxFlow<double>*>)
                + MAP_OVERHEAD);
    Rdmap.Clear(d, e);
    return false;
  }
  return true;
}
