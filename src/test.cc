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



/////////////////////////////////////////////////////////////////////////////////////////////////
// Test equivalency of original and modified hi_pr algorithm by generating random problems
// and comparing results. Also test the correctness of the parametric maxflow optimizations
// by modifying problems and checking every instance with the hi_pr algorithm.
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <set>

#include "maxflow.cc"

#define MAX_ARC_PER_NODE    10
#define HI_PR_BINARY        "/home/siberion/maxflow/hi_pr"

// Quick fix solution data structures; storage is partially redundant
struct mf_Arc
{
  int from, to, cap;
};

struct mf_Node
{
  mf_Arc in[MAX_ARC_PER_NODE], out[MAX_ARC_PER_NODE];
};

struct mf_Problem
{
  mf_Node source, sink, *set1, *set2;
  int a, b, cmin, cmax;
};

// Create a random bipartite maxflow problem where the left side has between a_min and a_max
// nodes, the right side has between b_min and b_max nodes (source and sink don't count) and
// all capacities are between cap_min and cap_max. The pointer must be freed with the
// delete_problem() function (see below). 0 may be returned if an error occurs. Use the
// MAX_ARC_PER_NODE define to tune how many arcs may go into or come out of a node (this does
// not apply to source and sink).
//
// The generation algorithm will randomly place between arcs_min and arcs_max arcs between the
// two sides of the bipartite network. This will not be possible if a_max and/or b_max are
// chosen too low. In that case, the function fails and returns 0.
mf_Problem *generate_problem(int a_min, int a_max, int b_min, int b_max, int arcs_min, int arcs_max, int cap_min, int cap_max)
{
  int a, b, x, y, n, m, arcs;
  
  std::set<int> listofarcs;
  
  if (a_max * MAX_ARC_PER_NODE < arcs_min || b_max * MAX_ARC_PER_NODE < arcs_min) return 0;
  
  mf_Problem *p = new mf_Problem;
  
  p->a = a_min + (rand() % (a_max - a_min + 1));
  p->b = b_min + (rand() % (b_max - b_min + 1));
  p->cmin = cap_min;
  p->cmax = cap_max;
  
  p->set1 = new mf_Node[p->a];
  p->set2 = new mf_Node[p->b];
  
  // Generate arcs from source to left set of nodes
  for (n = 0; n < p->a; ++n)
  {
    (p->set1 + n)->in[0].cap = cap_min + (rand() % (cap_max - cap_min + 1));
    (p->set1 + n)->in[0].from = 1;
    (p->set1 + n)->in[0].to = 2 + n;
    for (m = 1; m < MAX_ARC_PER_NODE; ++m) (p->set1 + n)->in [m].cap = 0;
    for (m = 0; m < MAX_ARC_PER_NODE; ++m) (p->set1 + n)->out[m].cap = 0;
  }
  
  // Generate arcs from right set of nodes to sink
  for (n = 0; n < p->b; ++n)
  {
    (p->set2 + n)->out[0].cap = cap_min + (rand() % (cap_max - cap_min + 1));
    (p->set2 + n)->out[0].from = 2 + p->a + n;
    (p->set2 + n)->out[0].to = 1 + p->a + p->b + 1;
    for (m = 1; m < MAX_ARC_PER_NODE; ++m) (p->set2 + n)->out[m].cap = 0;
    for (m = 0; m < MAX_ARC_PER_NODE; ++m) (p->set2 + n)->in [m].cap = 0;
  }
  
  arcs_max = ((arcs_max > (p->a * MAX_ARC_PER_NODE)) ? (p->a * MAX_ARC_PER_NODE) : arcs_max);
  arcs_max = ((arcs_max > (p->b * MAX_ARC_PER_NODE)) ? (p->b * MAX_ARC_PER_NODE) : arcs_max);
  arcs_max = ((arcs_max > (p->a * p->b)) ? (p->a * p->b) : arcs_max);
  arcs = arcs_min + (rand() % (arcs_max - arcs_min + 1));
  
  // Generate arcs in between the two sides
  for (n = 0; n < arcs; ++n)
  {
    do
    {
      do
      {
        a = rand() % p->a;
        for (x = 0; x < MAX_ARC_PER_NODE && (p->set1 + a)->out[x].cap != 0; ++x);
      }
      while ((p->set1 + a)->out[x].cap != 0);
      do
      {
        b = rand() % p->b;
        for (y = 0; y < MAX_ARC_PER_NODE && (p->set2 + b)->in[y].cap != 0; ++y);
      }
      while ((p->set2 + b)->in[y].cap != 0);
    }
    while (listofarcs.find((a * p->b) + b) != listofarcs.end());
    listofarcs.insert((a * p->b) + b);
    
    (p->set1 + a)->out[x].cap = cap_min + (rand() % (cap_max - cap_min + 1));
    (p->set1 + a)->out[x].from = 2 + a;
    (p->set1 + a)->out[x].to = 2 + p->a + b;
    (p->set2 + b)->in [y].cap = (p->set1 + a)->out[x].cap;
    (p->set2 + b)->in [y].from = 2 + a;
    (p->set2 + b)->in [y].to = 2 + p->a + b;
    
  }
  
  return p;
}

// Free memory used by a problem description
void delete_problem(mf_Problem *p)
{
  if (p == 0) return;
  delete[] p->set1;
  delete[] p->set2;
  delete p;
}

// Write the DIMACS representation of a problem to the designated stream
void write_DIMACS(mf_Problem *p, FILE *f)
{
  int a, b, arcs = p->a + p->b;
  
  for (a = 0; a < p->a; ++a)
  {
    for (b = 0; b < MAX_ARC_PER_NODE; ++b)
    {
      if ((p->set1 + a)->out[b].cap > 0) ++arcs;
    }
  }
  
  fprintf(f, "p max %d %d\n", 2 + p->a + p->b, arcs);
  fprintf(f, "n 1 s\nn %d t\n", 2 + p->a + p->b);
  
  for (a = 0; a < p->a; ++a)
  {
    fprintf(f, "a %d %d %d\n", (p->set1 + a)->in[0].from, (p->set1 + a)->in[0].to, (p->set1 + a)->in[0].cap);
    for (b = 0; b < MAX_ARC_PER_NODE; ++b)
    {
      if ((p->set1 + a)->out[b].cap > 0) fprintf(f, "a %d %d %d\n", (p->set1 + a)->out[b].from, (p->set1 + a)->out[b].to, (p->set1 + a)->out[b].cap);
    }
  }
  
  for (a = 0; a < p->b; ++a)
  {
    fprintf(f, "a %d %d %d\n", (p->set2 + a)->out[0].from, (p->set2 + a)->out[0].to, (p->set2 + a)->out[0].cap);
  }
  fprintf(f, "e\n");
}

// Use the hi_pr algorithm to solve a maxflow problem and store the result in the variable
// pointed to by *flow. The function returns false if an error occurred.
bool hi_pr(mf_Problem *p, long *flow)
{
  // Write DIMACS representation of problem into a temp file
  char tempfile[] = "/tmp/maxflow.tmp";
  FILE *t = fopen(&tempfile[0], "wb");
  write_DIMACS(p, t);
  fclose(t);
  
  // Run hi_pr on problem
  char cmdline[80];
  sprintf(&cmdline[0], "%s < %s", HI_PR_BINARY, &tempfile[0]);
  FILE *hipr = popen(&cmdline[0], "r");
  char line[80];
  
  // Read output and find the "c flow ..." line
  if (hipr)
  {
    do
    {
      fgets(&line[0], 80, hipr);
      if (!strncmp(&line[0], "c flow", 6))
      {
        *flow = strtol(&line[7], 0, 10);
        break;
      }
    }
    while (!feof(hipr));
    pclose(hipr);
    
    return (strncmp(&line[0], "c flow", 6) == 0);
  }
  
  exit(0);
  unlink(&tempfile[0]);
  return false;
}

// Modify a maxflow problem randomly, checking the result of the GenMaxflow class against
// the hi_pr algorithm in each iteration. The mods parameter specifies how many times the
// network's capacities should be randomly modified. After that, random arcs will be
// removed until the flow becomes 0.
//
// The function returns a character which indicates the result. A '.' means success, any
// other character means that there was some sort of error. See function body for details.
char compare_algorithms(mf_Problem *p, int mods)
{
  int a, b;
  ParametricMaxflow<long> mf;
  long flow1, flow2, cap;
  
  mf.SetOptimize(true);
  
  // Initialize the GenMaxflow class with the problem
  if (!mf.BeginNetwork(p->a + p->b + 2)) return '1';
  
  for (a = 0; a < p->a; ++a)
  {
    if (!mf.AddArc(mf.Source(), (p->set1 + a)->in[0].to, (p->set1 + a)->in[0].cap)) return '2';
    for (b = 0; b < MAX_ARC_PER_NODE; ++b)
    {
      if ((p->set1 + a)->out[b].cap > 0)
      {
        if (!mf.AddArc((p->set1 + a)->out[b].from, (p->set1 + a)->out[b].to, (p->set1 + a)->out[b].cap)) return '2';
      }
    }
  }
  
  for (a = 0; a < p->b; ++a)
  {
    if (!mf.AddArc((p->set2 + a)->out[0].from, mf.Sink(), (p->set2 + a)->out[0].cap)) return '2';
  }
  
  if (!mf.EndNetwork()) return '3';
  
  // Network should be classified as simple and non-parametric.
  if (!mf.IsSimple()) return 's';
  if (mf.IsParametric()) return 'p';
  
  // Check that the flow is the same as that computed by hi_pr
  flow1 = mf.Flow();
  if (!hi_pr(p, &flow2)) return '4';
  if (flow1 != flow2) return '5';
  
  do
  {
    // Make a random modification. The branches correspond to parts of the network
    switch (rand() % 3)
    //switch (1)
    {
    case 0: // Source arc
      a = rand() % p->a;
      if ((p->set1 + a)->in[0].cap == 0) continue;
      if (mods > 0)
      {
        --mods;
        cap = p->cmin + (rand() % (p->cmax - p->cmin + 1));
        if (!mf.SetArcCapacity(mf.Source(), (p->set1 + a)->in[0].to, cap)) return '\\';
        (p->set1 + a)->in[0].cap = cap;
      }
      else
      {
        if (!mf.RemoveArc(mf.Source(), (p->set1 + a)->in[0].to)) return '/';
        (p->set1 + a)->in[0].cap = 0;
      }
      break;
    case 1: // Middle arc
      a = rand() % p->a;
      for (b = MAX_ARC_PER_NODE - 1; (p->set1 + a)->out[b].cap == 0 && b != -1; --b);
      if (b == -1) continue;
      if (mods > 0)
      {
        --mods;
        cap = p->cmin + (rand() % (p->cmax - p->cmin + 1));
        if (!mf.SetArcCapacity((p->set1 + a)->out[b].from, (p->set1 + a)->out[b].to, cap)) return '\\';
        (p->set1 + a)->out[b].cap = cap;
      }
      else
      {
        if (!mf.RemoveArc((p->set1 + a)->out[b].from, (p->set1 + a)->out[b].to)) return '/';
        (p->set1 + a)->out[b].cap = 0;
      }
      break;
    case 2: // Sink arc
      a = rand() % p->b;
      if ((p->set2 + a)->out[0].cap == 0) continue;
      if (mods > 0)
      {
        --mods;
        cap = p->cmin + (rand() % (p->cmax - p->cmin + 1));
        if (!mf.SetArcCapacity((p->set2 + a)->out[0].from, mf.Sink(), cap)) return '\\';
        (p->set2 + a)->out[0].cap = cap;
      }
      else
      {
        if (!mf.RemoveArc((p->set2 + a)->out[0].from, mf.Sink())) return '/';
        (p->set2 + a)->out[0].cap = 0;
      }
      break;
    }
    
    // Verify that the flows are still the same
    flow1 = mf.Flow();
    if (!hi_pr(p, &flow2)) return '6';  
    if (flow1 != flow2)
    {
      //printf("removed %d --> %d\n", (p->set1 + a)->out[b].from, (p->set1 + a)->out[b].to);
      //write_DIMACS(p, stdout);
      //exit(0);
      return '7';
    }
  }
  while (flow2 > 0);
  
  return '.';
}

void GetRuntimes(rusage *r1, rusage *r2, double *t_usr, double *t_sys)
{
  *t_usr = double((r2->ru_utime.tv_sec * 1000) - (r1->ru_utime.tv_sec * 1000)) + ((double(r2->ru_utime.tv_usec) * 0.001) - (double(r1->ru_utime.tv_usec) * 0.001));
  *t_sys = double((r2->ru_stime.tv_sec * 1000) - (r1->ru_stime.tv_sec * 1000)) + ((double(r2->ru_stime.tv_usec) * 0.001) - (double(r1->ru_stime.tv_usec) * 0.001));
}

void test_for_runtime(mf_Problem *p, int arcs)
{
  int n, i, j, a, b;
  ParametricMaxflow<long> mf;
  
  // Initialize the GenMaxflow class with the problem
  assert(mf.BeginNetwork(p->a + p->b + 2));
  
  for (a = 0; a < p->a; ++a)
  {
    assert(mf.AddArc(mf.Source(), (p->set1 + a)->in[0].to, (p->set1 + a)->in[0].cap));
    for (b = 0; b < MAX_ARC_PER_NODE; ++b)
    {
      if ((p->set1 + a)->out[b].cap > 0)
      {
        assert(mf.AddArc((p->set1 + a)->out[b].from, (p->set1 + a)->out[b].to, (p->set1 + a)->out[b].cap));
      }
    }
  }
  
  for (a = 0; a < p->b; ++a)
  {
    assert(mf.AddArc((p->set2 + a)->out[0].from, mf.Sink(), (p->set2 + a)->out[0].cap));
  }
  
  assert(mf.EndNetwork());
  
  mf.Flow();
  
  for (i = 0, j = 0, n = 0; n < arcs; ++n, ++j)
  {
    while (j < MAX_ARC_PER_NODE && (p->set1 + i)->out[j].cap == 0) ++j;
    if (j == MAX_ARC_PER_NODE)
    {
      ++i;
      j = 0;
      continue;
    }
    mf.RemoveArc((p->set1 + i)->out[j].from, (p->set1 + i)->out[j].to);
    mf.Flow();
  }
}

int main()
{
  srand(time(0));
  mf_Problem *p;
  int n;
  char c;
  
  // Run 1000 tests. If you see only dots, everything is just fine.
  for (n = 0; n < 1000; ++n)
  {
    //p = generate_problem(2, 2, 2, 2, 2, 4, 1, 10);
    p = generate_problem(5, 20, 5, 20, 10, 60, 1, 40);
    if (!p)
    {
      printf("*\n");
      return 0;
    }
    c = compare_algorithms(p, 10);
    printf("%c", c);
    delete_problem(p);
    fflush(stdout);
    if (c == 's') break;
  }
  
  /*
  int m;
  double tusr, tsys, ttotal;
  rusage r0, r1;
  for (n = 10; n <= 200; n += 10)
  {
    ttotal = 0.0;
    for (m = 0; m < 1000; ++m)
    {
      p = generate_problem(20, 20, 20, 20, 1, 50, n, n);
      getrusage(0, &r0);
      test_for_runtime(p, n);
      getrusage(0, &r1);
      delete_problem(p);
      GetRuntimes(&r0, &r1, &tusr, &tsys);
      total += tusr + tsys;
    }
    printf("%d %.6f\n", n, ttotal);
  }*/
  
  return 0;
}

/*#include <sys/time.h>
#include <sys/resource.h>
#include <stdlib.h>
#include "Strong.h"

void GetRuntimes(rusage *r1, rusage *r2, double *t)
{
  double t_usr, t_sys;
  t_usr = double((r2->ru_utime.tv_sec * 1000) - (r1->ru_utime.tv_sec * 1000)) + ((double(r2->ru_utime.tv_usec) * 0.001) - (double(r1->ru_utime.tv_usec) * 0.001));
  t_sys = double((r2->ru_stime.tv_sec * 1000) - (r1->ru_stime.tv_sec * 1000)) + ((double(r2->ru_stime.tv_usec) * 0.001) - (double(r1->ru_stime.tv_usec) * 0.001));
  *t = t_usr + t_sys;
}

int main()
{
  int i, j, n;
  Sparse *s;
  StrongMF *s_mf;
  StrongCMF *s_cmf;
  rusage r0, r1;
  double t, t_mf, t_cmf;
  
  for (i = 10; i <= 200; i += 10)
  //i = 100;
  {
    //for (i = 10; i <= 200; i += 10)
    //i = 100;
    {
      t_mf = 0.0;
      t_cmf = 0.0;
      //for (n = 0; n < 4; ++n)
      {
        s = DTMC(i, 0.5);
        //delete s;
        
        /*getrusage(0, &r0);
        s_mf = new StrongMF(s);
        s_mf->run();
        getrusage(0, &r1);
        GetRuntimes(&r0, &r1, &t);
        t_mf += t;
        delete s_mf;*/
        
        getrusage(0, &r0);
        s_cmf = new StrongCMF(s);
        s_cmf->run();
        getrusage(0, &r1);
        GetRuntimes(&r0, &r1, &t);
        t_cmf += t;
        delete s_cmf;
        
        delete s;
        
        //fprintf(stderr, ".");
      }
      printf("%d %f %f\n", i, t_mf, t_cmf);
      fprintf(stderr, "%d %f %f\n", i, t_mf, t_cmf);
      fflush(stdout);
    }
  }
  return 0;
}
*/
