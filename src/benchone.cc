/*****************************************************************************/
/*!
 *   Copyright 2009-2014 Jonathan Bogdoll, Holger Hermanns, Lijun Zhang,
 *                       David N. Jansen
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



#include <string>
#include <cmath>
#include "benchmark.h"
#include "Strong.h"
#include "StrongQ.h"
#include "StrongCR.h"


int LabelFunction(void *userdata, int s)
{
  int &n_states = *(int*)(((void**)userdata)[0]);
  int *labels = *(int**)(((void**)userdata)[1]);
  int &n_labels = *(int*)(((void**)userdata)[2]);
  
  if (s >= 0 && s < n_states) return (labels ? labels[s] : 0);
  
  return (labels ? n_labels : 1);
}

/* Reads the label file. It is not required that every state has a label.
Labels can be any string. The label file may contain comment lines (with #);
however, if a comment appears on a line with a label, it is taken to be part
of the label, so it will be regarded as a different label than one without
comment. Every line may specify a single state, a state range (n1-n2) or a list
of states (n1,n2,n3...); lists and ranges cannot be mixed. Comment lines should
not be longer than 511 characters.
(A label may actually consist of multiple words. This can be used to specify
states that satisfy more than one atomic proposition. Note, however, that
``a b'' and ``b a'' are treated as different labels.) */
const char *load_label_data(const char *fn, int states, int **assoc, int &num_labels)
{
  FILE *f = 0;
  char *buffer, *cp, *np;
  int amountread, state, n, f_state, cur_label, next_label;
  std::set<int> t_states;
  std::map<std::string,int> t_labels;
  std::string f_label;
  std::set<int>::iterator i;
  
  if (!fn || *fn == 0)
  {
    *assoc = 0;
    return 0;
  }
  else if (*fn == '!')
  {
    num_labels = (int)strtoul(fn + 1, 0, 10);
    if (num_labels == 1)
    {
      *assoc = 0;
      return 0;
    }
    else if (num_labels != 0)
    {
      *assoc = new int[states];
      for (n = 0; n < states; ++n) (*assoc)[n] = n % num_labels;
      return 0;
    }
  }
  
  f = fopen(fn, "rb");
  if (!f) return "File not readable";
  
  buffer = new char[1024];
  buffer[1023] = '\0';
  cp = buffer;
  state = 0;
  next_label = 1;
  
  *assoc = new int[states];
  for (n = 0; n < states; ++n) (*assoc)[n] = 0;
  
  amountread = fread(buffer, 1, 1023, f);
  buffer[amountread] = '\0';
    while (*cp)
    {
      switch (state)
      {
      case 0:
        /* expect whitespace, then a state number, list or range, or comment */
        while (isspace(*cp)) ++cp;
        if (*cp == '#')
        {
          np = strchr(cp, '\n');
          if (!np) cp = cp + strlen(cp) - 1;
          else cp = np + 1;
          break;
        }
        state = 1;
        break;
      case 1:
        /* expect a state number, list or range */
        if (isdigit(*cp))
        {
          t_states.clear();
          t_states.insert(strtol(cp, &cp, 10));
          if (*cp == ',') ++cp, state = 2;
          else if (*cp == '-') ++cp, state = 3;
          else if (isblank(*cp)) state = 4;
          else
          {
            fclose(f);
            delete [] buffer;
            delete [] *assoc;
            *assoc = NULL;
            return "Format error; expected list, range or single state";
          }
        }
        else
        {
          fclose(f);
          delete [] buffer;
          delete [] *assoc;
          *assoc = NULL;
          return "Format error; expected state";
        }
        break;
      case 2:
        /* expect a state number (as part of a list) */
        if (!isdigit(*cp))
        {
          fclose(f);
          delete [] buffer;
          delete [] *assoc;
          *assoc = NULL;
          return "Format error; expected state number (in list)";
        }
        t_states.insert(strtol(cp, &cp, 10));
        if (*cp == ',')
        {
          ++cp;
          break;
        }
        else if (isblank(*cp)) state = 4;
        else
        {
          fclose(f);
          delete [] buffer;
          delete [] *assoc;
          *assoc = NULL;
          return "Format error; expected comma or end-of-list";
        }
        break;
      case 3:
        /* expect a state number (as upper bound of a range) */
        if (!isdigit(*cp))
        {
          fclose(f);
          delete [] buffer;
          delete [] *assoc;
          *assoc = NULL;
          return "Format error; expected state number (in range)";
        }
        f_state = strtol(cp, &cp, 10);
        for (n = (*t_states.begin()) + 1; n <= f_state; ++n) t_states.insert(n);
        if (!isblank(*cp))
        {
          fclose(f);
          delete [] buffer;
          delete [] *assoc;
          *assoc = NULL;
          return "Format error; expected white-space after range";
        }
        state = 4;
        break;
      case 4:
        /* expect white space, then a label, then a newline */
        while (isblank(*cp)) ++cp;
        np = strchr(cp, '\n');
        if (!np)
        {
          fclose(f);
          delete [] buffer;
          delete [] *assoc;
          *assoc = NULL;
          return "Format error; expected line-break after label";
        }
        while (isspace(*np)) --np;
        f_label.assign(cp, np - cp + 1);
        if (t_labels.find(f_label) != t_labels.end()) cur_label = t_labels[f_label];
        else
        {
          t_labels[f_label] = next_label;
          cur_label = next_label++;
        }
        for (i = t_states.begin(); i != t_states.end(); ++i)
        {
          if (*i >= 0 && *i < states) (*assoc)[*i] = cur_label;
        }
        cp = np + 1;
        state = 0;
        break;
      }
      if (cp - buffer >= 512)
      {
        cp = stpcpy(buffer, cp);
        /* cp = buffer + strlen(buffer); */
        if (!feof(f)) {
          amountread = fread(cp, 1, 1023 - (cp - buffer), f);
          cp[amountread] = '\0';
        }
        cp = buffer;
      }
    }
  
  fclose(f);
  delete [] buffer;
  
  t_states.clear();
  for (n = 0; n < states; ++n) t_states.insert((*assoc)[n]);
  num_labels = t_states.size();
  
  return 0;
}

int main(int argc, char *argv[])
{
  Benchmark b;
  double epsilon;
  const char *errmsg;
  unsigned long avg, type, flags = 0;
  int n_states, n_labels, *labels;
  void *lf_info[3] = {&n_states, &labels, &n_labels};
  bool continuous_model;
  FILE *f, *out;
  ProbabilisticModel *pm;
  SimulationRelation *ss;
  SimulationStatistics stats;

  if (argc < 7 || argc > 8)
  {
    fprintf(stderr, "Usage: %s simtype mtype labels epsilon avg model [dump]\n\n", argv[0]);
    fprintf(stderr, "simtype - simulation type: must always be \"strong\"\n");
    fprintf(stderr, "mtype   - model type: dtmc, ctmc, pa, cpa\n");
    fprintf(stderr, "labels  - file to read state/label association from\n");
    fprintf(stderr, "          (or \"!n\" for n labels, assigned in turn)\n");
    fprintf(stderr, "epsilon - floating-point approximation threshold (-1 = automatic)\n");
    fprintf(stderr, "avg     - perform benchmark this many times to obtain average values\n");
    fprintf(stderr, "model   - input model filename. may be \"-\" to read from stdin\n");
    fprintf(stderr, "dump    - optionally, dump the relation in this file\n\n");
    return 0;
  }

#ifdef OPT_PARTITION
  flags |= 0x1;
#endif
#ifdef OPT_QUOTIENT
  flags |= 0x2;
#endif
#ifdef OPT_P_INVARIANT
  flags |= 0x4;
#endif
#ifdef OPT_SIGNIFICIANT_ARC
  flags |= 0x8;
#endif
#ifdef OPT_CACHE_NETS
  flags |= 0x10;
#endif
#ifdef OPT_CRAFA_RANZATO
  /* use the optimization by Crafa/Ranzato 2012: avoid redoing calculations for
     multiple states that have the same successor distribution. Implemented by
     David N. Jansen. */
  flags |= 0x20;
#  ifndef OPT_PARTITION
#    error "If OPT_CRAFA_RANZATO is set, also OPT_PARTITION must be set."
#  endif
#  ifndef OPT_CACHE_NETS
#    error "If OPT_CRAFA_RANZATO is set, also OPT_CACHE_NETS must be set."
#  endif
#  ifdef OPT_QUOTIENT
#    error "Please do not define OPT_QUOTIENT and OPT_CRAFA_RANZATO at the same time."
#  endif
#endif

  if (strcmp(argv[1], "strong"))
  {
    fprintf(stderr, "Simulation type '%s' not recognized\n", argv[1]);
    return -1;
  }

  if (!strcmp(argv[2], "dtmc")) type = 1;
  else if (!strcmp(argv[2], "ctmc")) type = 2;
  else if (!strcmp(argv[2], "pa")) type = 3;
  else if (!strcmp(argv[2], "cpa")) type = 4;
  else
  {
    fprintf(stderr, "Model type '%s' not recognized\n", argv[2]);
    return -1;
  }
  
  epsilon = strtod(argv[4], 0);
  avg = strtoul(argv[5], 0, 10);
  f = ((!strcmp(argv[6], "-")) ? stdin : fopen(argv[6], "rb"));
  if (!f)
  {
    fprintf(stderr, "Failed to open model '%s'\n", argv[6]);
    return -1;
  }
  
  if (argc >= 8)
  {
    out = fopen(argv[7], "wb");
    if (!out)
    {
      fprintf(stderr, "Failed to open output file '%s'\n", argv[7]);
      return -1;
    }
  }
  else out = 0;

#if defined(OPT_QUOTIENT) || defined(OPT_CRAFA_RANZATO)
  switch (type)
  {
  case 1:
    pm = new MarkovChain;
    continuous_model = false;
    break;
  case 2:
    pm = new MarkovChain;
    continuous_model = true;
    break;
  case 3:
    pm = new ProbabilisticAutomaton;
    continuous_model = false;
    break;
  case 4:
    pm = new ProbabilisticAutomaton;
    continuous_model = true;
    break;
  }
#  ifdef OPT_QUOTIENT
  ss = new StrongSimulation_Quotient;
#  else
  ss = new StrongSimulation_CR;
#  endif
#else//defined(OPT_QUOTIENT) || defined(OPT_CRAFA_RANZATO)
  switch (type)
  {
  case 1:
    pm = new MarkovChain;
    /*if (simtype == 1) ss = new WeakSimulation_MC;
    else */ss = new StrongSimulation_MC;
    continuous_model = false;
    break;
  case 2:
    pm = new MarkovChain;
    ss = new StrongSimulation_MC;
    continuous_model = true;
    break;
  case 3:
    pm = new ProbabilisticAutomaton;
    ss = new StrongSimulation_PA;
    continuous_model = false;
    break;
  case 4:
    pm = new ProbabilisticAutomaton;
    ss = new StrongSimulation_PA;
    continuous_model = true;
    break;
  }
#endif//defined(OPT_QUOTIENT) || defined(OPT_CRAFA_RANZATO)

  pm->Parse(f, continuous_model);
  if (strcmp(argv[5], "-")) fclose(f);

  if (epsilon < 0) epsilon = pm->PrecisionThreshold();

  if ((errmsg = load_label_data(argv[3], pm->States(), &labels, n_labels)))
  {
    fprintf(stderr, "Failed to load labels from '%s': %s\n", argv[3], errmsg);
    return -1;
  }

  n_states = pm->States();
  ss->SetLabelFunction(LabelFunction, lf_info);
  ss->SetFPPrecision(epsilon);

  b.Bench(pm, ss, 0, 1, avg, out);
  fprintf(stderr, "<%2lu>", flags);
  
  
  if (out) fclose(out);

  memcpy(&stats, b.GetStats(0), sizeof(SimulationStatistics));

  printf("%e %e %e %u %u %u %u %u %u %u %zu %zu %zu %zu %zu %u %u %u %u "
    "%.3e %.3e %.3e\n",
    b.GetUserTime(0), b.GetSystemTime(0), b.GetRealTime(0), stats.num_partitions,
    stats.num_iterations, stats.num_initial_pairs, stats.num_final_pairs,
    stats.num_maxflow,
    stats.num_p_invariant_fails, stats.num_sig_arc_fails, stats.mem_relation_map,
    stats.mem_partition_map, stats.mem_relation,
    stats.mem_maxflow, stats.mem_model, stats.min_complexity, stats.max_complexity,
    stats.num_nets_cached, stats.num_cache_hits,
    b.GetUserTimeStdev(0), b.GetSystemTimeStdev(0), b.GetRealTimeStdev(0));

  delete pm;
  delete ss;
  if (labels) delete [] labels;

  return 0;
}
