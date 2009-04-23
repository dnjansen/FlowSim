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
#include "StrongMF.h"
#include "Strong.h"

#include <iostream>
#include <sys/resource.h>

MarkovChain *DTMC(int n, int a, int b);
ProbabilisticAutomaton *RandomPA(int n, double p, int a, int ma, int Ma, double r);

static void GetRuntimes(rusage *r1, rusage *r2, double *t_usr, double *t_sys)
{
  *t_usr = double((r2->ru_utime.tv_sec) - (r1->ru_utime.tv_sec)) + ((double(r2->ru_utime.tv_usec) * 0.000001) - (double(r1->ru_utime.tv_usec) * 0.000001));
  *t_sys = double((r2->ru_stime.tv_sec) - (r1->ru_stime.tv_sec)) + ((double(r2->ru_stime.tv_usec) * 0.000001) - (double(r1->ru_stime.tv_usec) * 0.000001));
}

int main___()
{
  rusage t0, t1;
  int n, m;
  //MarkovChain *sp;
  //StrongSimulation_DTMC sim_cmf;
  ProbabilisticAutomaton *pa;
  StrongSimulation_PA sim_pa;
  double t_usr[4], t_sys, ttotal;
  
  for (n = 1; n <= 10; ++n)
  {
    StrongSimulation::SetFPPrecision(1e-10);
    ttotal = 0.0;
    for (m = 0; m < 10; ++m)
    {
      //sp = DTMC(1000, 1, 1);
      pa = RandomPA(100, 0.1, 4, 1, 4, 0.2);
      getrusage(0, &t0);
      //sim_cmf.Simulate(sp, 0);
      sim_pa.Simulate(pa, 0);
      getrusage(0, &t1);
      //delete sp;
      delete pa;
      GetRuntimes(&t0, &t1, &t_usr[0], &t_sys);
      ttotal += t_usr[0];
    }
    
    printf("%2d %.8f\n", n, ttotal * 0.1);
  }
  
  return 0;
}

unsigned long simulation(MarkovChain *s)
{
  //rusage t0, t1;
  //double t_sys, t_usr;
  
  StrongSimulation::SetFPPrecision(1e-8); // Models exported by PRISM are inaccurately represented
  //getrusage(0, &t0);
  StrongSimulation_DTMC rel3;
  rel3.Simulate(s, 0);
  return 0;
  //getrusage(0, &t1);
  //GetRuntimes(&t0, &t1, &t_usr, &t_sys);
  
  //return t_usr;
}

int main__unused()
{
  MarkovChain *s;
  FILE *f;
  char filename[80];
  unsigned int n;//, m;
  int mods[] = {3,4,5,6,8};
  //double ttotal, tmin, tmax, t;
  unsigned long memstats[sizeof(mods) / sizeof(mods[0])];
  
  for (n = 0; n < sizeof(mods) / sizeof(int); ++n)
  {
    sprintf(&filename[0], "examples/leader3_%d.model", mods[n]);
    f = fopen(&filename[0], "rb");
    s = new MarkovChain();
    s->Parse(f);
    fclose(f);
    
    /*tmin = 1e+100;
    tmax = 0.0;
    ttotal = 0.0;
    for (m = 0; m < 4; ++m)
    {
      t = simulation(s);
      if (t < tmin) tmin = t;
      if (t > tmax) tmax = t;
      ttotal += t;
    }*/
    
    memstats[n] = simulation(s);
    
    delete s;
    
    //printf("%s\t%.5f\t%.5f\t%.5f\n", &filename[0], tmin, ttotal * .25, tmax);
    //fflush(stdout);
  }
  
  f = fopen("values.txt", "a+b");
  for (n = 0; n < sizeof(mods) / sizeof(int); ++n)
  {
    if (n > 0) fprintf(f, "\t");
    fprintf(f, "%.3f", double(memstats[n]) / 1024.0);
  }
  fprintf(f, "\n");
  fclose(f);
  
  return 0;
}

int main(int argc, char *argv[])
{
  /*bool DEBUG=false;
  
  int m=10,n=10;

  cout<<"argc:"<<argc<<endl;
  cout<<"argv[0]:"<< argv[0]<<endl;
  for(int i=1; i<argc; i++){
    cout<<"argv["<<i<<"]:"<<argv[i]<<endl;
    if(strcmp(argv[i], "-leftnode") == 0 || strcmp(argv[i], "-l") == 0)
      m = atoi(argv[++i]);
    else if(strcmp(argv[i], "-rightnode") == 0 || strcmp(argv[i], "-r") == 0)
      n = atoi(argv[++i]);
    else if(strcmp(argv[i], "-d")==0)
      DEBUG=true;
    else{      }
  }*/
  
  
  //return 0;
  
  if (argc == 2) stdin = freopen(argv[1], "rb", stdin);

  // sparse matrix
  //MarkovChain spa;
  ProbabilisticAutomaton spa;
  //parse the input stream. The command line is:
  // command < input.file
  spa.Parse();
  
  rusage t0, t1;
  double t_sys, t_usr;
  
  /*getrusage(0, &t0);
  StrongMF rel1(&spa);
  rel1.run();
  getrusage(0, &t1);
  GetRuntimes(&t0, &t1, &t_usr, &t_sys);
  
  printf(" MF: %.6f seconds\n", t_usr);*/
  
  StrongSimulation::SetFPPrecision(1e-8); // Models exported by PRISM are inaccurately represented
  getrusage(0, &t0);
  StrongSimulation_PA rel3;
  std::set<std::pair<int,int> > result;
  rel3.SetNumberOfLabels(3);
  rel3.SetOptimization(1);
  rel3.Simulate(&spa, &result);
  getrusage(0, &t1);
  GetRuntimes(&t0, &t1, &t_usr, &t_sys);
  
#if 0
#ifdef DEBUG
  printf("num_partitions = %u\n", rel3.stats.num_partitions);
  printf("num_iterations = %u\n", rel3.stats.num_iterations);
  printf("num_initial_pairs = %u\n", rel3.stats.num_initial_pairs);
  printf("num_final_pairs = %u\n", rel3.stats.num_final_pairs);
  printf("num_partition_trees = %u\n", rel3.stats.num_partition_trees);
  printf("num_partition_tree_clusters = %u\n", rel3.stats.num_partition_tree_clusters);
  printf("num_maxflow = %u\n", rel3.stats.num_maxflow);
  printf("\n");
  printf("mem_relation_map = %lu\n", rel3.stats.mem_relation_map);
  printf("mem_partition_map = %lu\n", rel3.stats.mem_partition_map);
  printf("mem_partition_trees = %lu\n", rel3.stats.mem_partition_trees);
  printf("mem_relation = %lu\n", rel3.stats.mem_relation);
  printf("mem_maxflow = %lu\n", rel3.stats.mem_maxflow);
#endif//DEBUG
#endif
  
  //FILE *fout = fopen("relation", "wb");
  for (std::set<std::pair<int,int> >::iterator it = result.begin(); it != result.end(); it++)
  {
    fprintf(stdout, "(%d,%d)\n", (*it).first, (*it).second);
  }
  //fclose(fout);
  
  //printf("%.5f\n", t_usr);
  
  return 0;
}
