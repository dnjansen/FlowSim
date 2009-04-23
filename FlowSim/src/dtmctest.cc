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



#ifndef DEBUG
#define DEBUG
#endif

#include "modelbuilder.h"
#include "prmodel.h"
#include "Strong.h"
#include "StrongMC.cc"
#include "Weak.h"
#include "WeakMC.cc"
#include "bench.cc"

MarkovChain *GenerateModel_1(int k, double p, double q)
{
  ModelBuilder mb;
  int n;
  
  assert(k >= 1);
  assert(p >= 0 && p <= 1);
  assert(q >= 0 && q <= 1);
  
  for (n = 0; n < k; ++n)
  {
    mb.Add(n, n + 1, p);
    mb.Add(n, k + n + 2, 1.0 - p);
    mb.Add(k + n + 1, n + 1, 1.0 - q);
    mb.Add(k + n + 1, k + n + 2, q);
  }
  
  mb.Add((2 * k) + 1, 2 * (k + 1), 1.0);
  
  return mb.BuildMC((2 * k) + 3, 0);
}

MarkovChain *GenerateModel_2(int k, double p, double q)
{
  ModelBuilder mb;
  int n;
  
  assert(k >= 1);
  assert(p > .1 && p <= 1);
  assert(q > .1 && q <= 1);
  
  for (n = 0; n < k; ++n)
  {
    mb.Add(n, n + 1, p);
    mb.Add(n, k + n + 2, 0.9 - p);
    mb.Add(n, 2 * (k + 1), 0.1);
    
    mb.Add(k + n + 1, n + 1, 0.9 - q);
    mb.Add(k + n + 1, k + n + 2, q);
    mb.Add(k + n + 1, k + n + 1, 0.1);
  }
  
  mb.Add((2 * k) + 1, 2 * (k + 1), 1.0);
  mb.Add(2 * (k + 1), 2 * (k + 1), 1.0);
  
  return mb.BuildMC((2 * k) + 3, 0);
}

int LabelFunction(void *userdata, int s)
{
  int k = (int)userdata;
  
  if (s < 0 || s > (2 * k) + 2) return 2;
  
  if (s <= 2 * k && s != k) return 0;
  if (s == k || s == (2 * k) + 1) return 1;
  if (s == (2 * k) + 2) return 0;
  
  return 2;
}

void usage()
{
  fprintf(stderr, "Usage: dtmctest\n");
  fprintf(stderr, "Usage: dtmctest --export [prefix]\n");
  fprintf(stderr, "Usage: dtmctest --bench {weak|strong} avg report\n\n");
  fprintf(stderr, "If --export is given, generated models are exported instead of simulated.\n");
  fprintf(stderr, "Models are exported as (prefix)_k_p_q.dtmc with the default prefix being \"dtmctest\".\n\n");
  fprintf(stderr, "If --bench is given, generated models are benchmarked avg times\n");
  fprintf(stderr, "and the results will be stored in the file identified by report.\n\n");
  exit(-1);
}

void print_stats(FILE *out, Benchmark &bm)
{
  SimulationStatistics stats;
  memcpy(&stats, bm.GetStats(0), sizeof(stats));
  fprintf(out, "%30s: %.8f\n", "user time", bm.GetUserTime(0));
  fprintf(out, "%30s: %.8f\n", "system time", bm.GetSystemTime(0));
  fprintf(out, "%30s: %.8f\n", "real time", bm.GetRealTime(0));
  fprintf(out, "%30s: %u\n",  "num iterations", stats.num_iterations);
  fprintf(out, "%30s: %u\n",  "num initial pairs", stats.num_initial_pairs);
  fprintf(out, "%30s: %u\n",  "num final pairs", stats.num_final_pairs);
  fprintf(out, "%30s: %u\n",  "num blocks", stats.num_partitions);
  fprintf(out, "%30s: %u\n",  "num maxflow", stats.num_maxflow);
  fprintf(out, "%30s: %u\n",  "num p-invariant violated", stats.num_p_invariant_fails);
  fprintf(out, "%30s: %u\n",  "num significant arc deleted", stats.num_sig_arc_fails);
  fprintf(out, "%30s: %u\n",  "num pmf nets cached", stats.num_nets_cached);
  fprintf(out, "%30s: %u\n",  "num pmf cache hits", stats.num_cache_hits);
  fprintf(out, "%30s: %lu\n", "mem relation map", stats.mem_relation_map);
  fprintf(out, "%30s: %lu\n", "mem partition map", stats.mem_partition_map);
  fprintf(out, "%30s: %lu\n", "mem relation", stats.mem_relation);
  fprintf(out, "%30s: %lu\n", "mem maxflow", stats.mem_maxflow);
  fprintf(out, "%30s: %lu\n", "mem model", stats.mem_model);
  fprintf(out, "\n");
}

int main(int argc, char *argv[])
{
  MarkovChain *mc;
  StrongSimulation_MC sims;
  WeakSimulation_MC simw;
  int k, k0, k1, avg;
  double p, q;
  char buf[80], *pnext;
  std::set<std::pair<int,int> > result;
  std::set<std::pair<int,int> >::iterator ri;
  FILE *rl, *report;
  bool exportmodel = false, benchmodel = false, benchweak;
  const char *exportprefix = "dtmctest", *benchreport;
  Benchmark bm;
  
  if (argc > 1)
  {
    if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) usage();
    if (!strcmp(argv[1], "--export"))
    {
      if (argc > 3) usage();
      if (argc > 2) exportprefix = argv[2];
      exportmodel = true;
    }
    else if (!strcmp(argv[1], "--bench"))
    {
      if (argc < 5) usage();
      if (!strcmp(argv[2], "weak")) benchweak = true;
      else if (!strcmp(argv[2], "strong")) benchweak = false;
      else usage();
      avg = strtol(argv[3], 0, 10);
      benchreport = argv[4];
      benchmodel = true;
    }
    else usage();
  }
  
  printf("please enter range of k as two space separated values (e.g.: 1 100): ");
  fgets(&buf[0], 80, stdin);
  k0 = strtol(&buf[0], &pnext, 10);
  k1 = strtol(pnext, 0, 10);
  
  printf("please enter space separated parameters p and q (e.g.: 0.5 0.5): ");
  fgets(&buf[0], 80, stdin);
  p = strtod(&buf[0], &pnext);
  q = strtod(pnext, 0);
  
  sims.SetFPPrecision((p < q ? p : q) * 0.01);
  simw.SetFPPrecision((p < q ? p : q) * 0.01);

  if (!exportmodel)
  {
    printf("Approximation threshold: %f\n", (p < q ? p : q) * 0.01);
    rl = fopen("resultlog", "wb");
    if (!benchmodel) printf("  k   n   m  weak[iter initial final] strong[iter initial final]  k+1 2k+1\n");
  }
  if (benchmodel)
  {
    report = fopen(benchreport, "wb");
    if (!report)
    {
      fprintf(stderr, "Cannot open '%s' (benchmark report) for writing, errno=%d\n", benchreport, errno);
      return -1;
    }
  }
  for (k = k0; k <= k1; ++k)
  {
    mc = GenerateModel_2(k, p, q);
    
    if (exportmodel)
    {
      sprintf(&buf[0], "%s_%0*d_%.2f_%.2f.dtmc", exportprefix, (int)floor(::log10((double)k1)) + 1, k, p, q);
      rl = fopen(&buf[0], "wb");
      mc->Write(rl);
      fclose(rl);
      delete mc;
      continue;
    }
    
    sims.SetLabelFunction(LabelFunction, (void*)k);
    simw.SetLabelFunction(LabelFunction, (void*)k);
    
    if (benchmodel)
    {
      printf("k = %*d [", (int)floor(::log10((double)k1)) + 1, k);
      fflush(stdout);
      if (benchweak)
      {
        fprintf(rl, "=== Final   weak relation of k=%d p=%.2f q=%.2f ===\n", k, p, q);
        bm.Bench(mc, &simw, 0, 1, avg, rl);
        fprintf(rl, "\n");
        fprintf(report, "===   Weak simulation statistics of k=%*d p=%.2f q=%.2f ===\n", (int)floor(::log10((double)k1)) + 1, k, p, q);
        print_stats(report, bm);
      }
      else
      {
        fprintf(rl, "=== Final strong relation of k=%d p=%.2f q=%.2f ===\n", k, p, q);
        bm.Bench(mc, &sims, 0, 1, avg, rl);
        fprintf(rl, "\n");
        fprintf(report, "=== Strong simulation statistics of k=%*d p=%.2f q=%.2f ===\n", (int)floor(::log10((double)k1)) + 1, k, p, q);
        print_stats(report, bm);
      }
      printf("...]\n");
      continue;
    }
    
    result.clear();
    sims.Simulate(mc, 0/*&result*/);
    simw.Simulate(mc, &result);
    
    printf("%3d %3d %3d %10d %7d %5d %12d %7d %5d %5d %4d\n", k, mc->States(), mc->Transitions(),
           simw.stats.num_iterations, simw.stats.num_initial_pairs, simw.stats.num_final_pairs, 
           sims.stats.num_iterations, sims.stats.num_initial_pairs, sims.stats.num_final_pairs, k + 1, (2 * k) + 1);
    
    fprintf(rl, "==== k=%d: %u pairs ====\n", k, result.size());
    for (ri = result.begin(); ri != result.end(); ++ri) fprintf(rl, "%4d %4d\n", ri->first, ri->second);
    fprintf(rl, "\n");
    
    delete mc;
  }
  
  if (exportmodel) printf("Models exported.\n");
  if (benchmodel) fclose(report);

  return 0;
}
