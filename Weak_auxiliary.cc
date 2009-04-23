#include <set>
#include <map>
#include <utility>
#include <vector>
#include <unistd.h>
#include <stdio.h>
#include <math.h>

#define HIPR_PROGRAM       "./hipr/hi_pr"

using namespace std;

// In the network defined by the set of arcs, contract the given vertices into a single vertex. It
// is assumed that there are no self-loops in the networks considered by this software.
void ContractVertices(const set<pair<pair<int,int>,double> > &arcs, set<pair<pair<int,int>,double> > &result, const set<int> &vertices, int id)
{
  map<int,double> in, out;
  map<int,double>::iterator j;
  set<pair<pair<int,int>,double> >::iterator i;
  bool bi, bo;
  
  // Initialize the maps with all predecessors and successors of the set of vertices
  for (i = arcs.begin(); i != arcs.end(); ++i)
  {
    bi = (vertices.find((*i).first.second) != vertices.end());
    bo = (vertices.find((*i).first.first) != vertices.end());
    
    if (bi && bo) continue;
    else if (bi && !bo) in[(*i).first.first] = 0.0;
    else if (!bi && bo) out[(*i).first.second] = 0.0;
    else result.insert(*i);
  }
  
  // Determine the sum of all ingoing and outgoing arc capacities as well as inserting
  // those arcs which remain unchanged into the result network
  for (i = arcs.begin(); i != arcs.end(); ++i)
  {
    bi = (vertices.find((*i).first.second) != vertices.end());
    bo = (vertices.find((*i).first.first) != vertices.end());
    
    if (bi && bo) continue;
    else if (bi && !bo) in[(*i).first.first] += (*i).second;
    else if (!bi && bo) out[(*i).first.second] += (*i).second;
  }
  
  // Insert all arcs going into the new contracted vertex
  for (j = in.begin(); j != in.end(); ++j)
  {
    result.insert(make_pair(make_pair(j->first, id), j->second));
  }
  
  // Inert all arcs leaving the new contracted vertex
  for (j = out.begin(); j != out.end(); ++j)
  {
    result.insert(make_pair(make_pair(id, j->first), j->second));
  }
}

// Uses the external hipr program to compute mincuts and maximum flow. This is only a
// temporary provision and will be changed in time since this implementation is obviously
// not parametric.
void AnalyzeNetwork(const set<pair<pair<int,int>,double> > &arcs, int source, int sink, double param1,
                    double param2, set<int> *max_s_mincut, set<int> *max_t_mincut, double *maxflow,
                    FILE *dbgout = 0, const char *debug = 0)
{
  char network[80];
  char buffer[160];
  FILE *f;
  int fd, tempstr_pos, nodes, n, tail, head;
  set<pair<pair<int,int>,double> >::iterator i;
  set<int> allset, lset, rset, notX1;
  map<int,int> nodeID;
  vector<int> nodeIDrev;
  bool seen_mincut_marker;
  double flow_value;
  
  nodeIDrev.clear();
  nodeIDrev.push_back(source);
  
  // Determine which nodes are on which side of the bipartite network, determine set of
  // all nodes and determine reverse node ID translation (we translate for hi_pr)
  for (i = arcs.begin(); i != arcs.end(); ++i)
  {
    if ((*i).first.first == source) lset.insert((*i).first.second);
    else if ((*i).first.second == sink) rset.insert((*i).first.first);
    
    if (i->first.first != source && allset.find(i->first.first) == allset.end()) nodeIDrev.push_back(i->first.first);
    if (i->first.second != sink && allset.find(i->first.second) == allset.end()) nodeIDrev.push_back(i->first.second);
    
    allset.insert((*i).first.first);
    allset.insert((*i).first.second);
  }
  
  nodes = allset.size();
  
  nodeIDrev.push_back(sink);
  
  // Determine forward node ID translation
  for (n = 0; n < (int)nodeIDrev.size(); ++n)
  {
    nodeID[nodeIDrev[n]] = n + 1;
  }
  
  // Generate temporary file for interaction with hipr
  strcpy(&network[0], HIPR_PROGRAM);
  strcat(&network[0], " < ");
  tempstr_pos = strlen(&network[0]);
  strcpy(&network[tempstr_pos], "/tmp/maxflow.XXXXXXXXXX");
  fd = mkstemp(&network[tempstr_pos]);
  
  // Construct hipr definition of forward network
  f = fdopen(fd, "wb");
  fprintf(f, "p max %u %u\nn %d s\nn %d t\n", nodes, arcs.size(), nodeID[source], nodeID[sink]);
  if (debug) fprintf(dbgout, "%s p max %u %u\n%s n %d s\n%s n %d t\n", debug, nodes, arcs.size(), debug, nodeID[source], debug, nodeID[sink]);
  for (i = arcs.begin(); i != arcs.end(); ++i)
  {
    tail = nodeID[i->first.first];
    head = nodeID[i->first.second];
    if (i->first.first == source) fprintf(f, "a %d %d %.20f\n", tail, head, (*i).second * param1);
    else if (i->first.second == sink) fprintf(f, "a %d %d %.20f\n", tail, head, (*i).second * param2);
    else fprintf(f, "a %d %d %.20f\n", tail, head, (*i).second);
    
    if (debug)
    {
      if (i->first.first == source) fprintf(dbgout, "%s a %d %d %.20f\n", debug, tail, head, (*i).second * param1);
      else if (i->first.second == sink) fprintf(dbgout, "%s a %d %d %.20f\n", debug, tail, head, (*i).second * param2);
      else fprintf(dbgout, "%s a %d %d %.20f\n", debug, tail, head, (*i).second);
    }
  }
  fclose(f);
  
  // Execute hipr and parse results
  f = popen(&network[0], "r");
  seen_mincut_marker = false;
  flow_value = -1.0;
  fgets(&buffer[0], 160, f);
  do
  {
    if (!strncmp(&buffer[0], "c flow:", 7)) flow_value = strtod(&buffer[7], 0);
    
    if (seen_mincut_marker)
    {
      if (buffer[2] < '0' || buffer[2] > '9') break;
      //if (sink_mincut) sink_mincut->insert(strtol(&buffer[2], 0, 10));
      allset.erase(nodeIDrev[strtol(&buffer[2], 0, 10) - 1]);
    }
    else if (!strncmp(&buffer[0], "c nodes on the sink side", 24)) seen_mincut_marker = true;
    
    fgets(&buffer[0], 160, f);
  }
  while (!feof(f));
  pclose(f);
  
  unlink(&network[tempstr_pos]);
  
  if (maxflow) *maxflow = flow_value;
  if (max_s_mincut) *max_s_mincut = allset;
  
  // Generate temporary file for interaction with hipr
  strcpy(&network[0], HIPR_PROGRAM);
  strcat(&network[0], " < ");
  tempstr_pos = strlen(&network[0]);
  strcpy(&network[tempstr_pos], "/tmp/maxflow.XXXXXXXXXX");
  fd = mkstemp(&network[tempstr_pos]);
  
  // Construct hipr definition of backward network
  f = fdopen(fd, "wb");
  fprintf(f, "p max %u %u\nn %d s\nn %d t\n", nodes, arcs.size(), nodeID[sink], nodeID[source]);
  for (i = arcs.begin(); i != arcs.end(); ++i)
  {
    tail = nodeID[i->first.first];
    head = nodeID[i->first.second];
    if (i->first.second == sink) fprintf(f, "a %d %d %.20f\n", head, tail, (*i).second * param2);
    else if (i->first.first == source) fprintf(f, "a %d %d %.20f\n", head, tail, (*i).second * param1);
    else fprintf(f, "a %d %d %.20f\n", head, tail, (*i).second);
  }
  fclose(f);
  
  // Execute hipr and parse results
  f = popen(&network[0], "r");
  if (max_t_mincut) max_t_mincut->clear();
  seen_mincut_marker = false;
  fgets(&buffer[0], 160, f);
  do
  {
    if (seen_mincut_marker)
    {
      if (buffer[2] < '0' || buffer[2] > '9') break;
      if (max_t_mincut) max_t_mincut->insert(nodeIDrev[strtol(&buffer[2], 0, 10) - 1]);
    }
    else if (!strncmp(&buffer[0], "c nodes on the sink side", 24)) seen_mincut_marker = true;
    
    fgets(&buffer[0], 160, f);
  }
  while (!feof(f));
  pclose(f);
    
  fflush(stdout);

  unlink(&network[tempstr_pos]);
}

// A private routine that is a generalization of FindSmallestBreakpoint(). Instead of finding the
// overall smallest breakpoint, it finds the smallest breakpoint to the right of a line segment
// identified by its alpha and beta coefficients.
double FindNextBreakpoint(const set<pair<pair<int,int>,double> > &arcs, int source, int sink, double alpha0, double beta0, double *alpha1, double *beta1)
{
  set<pair<pair<int,int>,double> >::const_iterator i, j;
  double capsum, lambda, alpha, beta, cut;
  
  set<int> X1;
  
  // Compute sum of all source arc capacities, the largest sink arc capacity
  for (i = arcs.begin(), capsum = 0.0, lambda = 0.0; i != arcs.end(); ++i)
  {
    if (i->first.first == source) capsum += i->second;
    if (i->first.second == sink) lambda += i->second;
  }
  for (i = arcs.begin(); i != arcs.end(); ++i)
  {
    if (i->first.second == sink && i->second < lambda) lambda = i->second;
  }
  
  lambda = (capsum / lambda) + 1;
  
  // Repeat step 1 as often as necessary
  while (true)
  {
    AnalyzeNetwork(arcs, source, sink, 1.0, lambda, 0, &X1, &cut);
    
    //printf("%f %f %f %f %f %f\n", alpha0, beta0, alpha, beta, lambda, cut);
    
    //if (alpha0 + (lambda * beta0) == cut) break;
    if (fabs(alpha0 + (lambda * beta0) - cut) < 0.0000001) break;
    
    // Compute non-constant coefficient beta
    beta = 0.0;
    for (i = arcs.begin(); i != arcs.end(); ++i)
    {
      if (i->first.second == sink)
      {
        if (X1.find(i->first.first) != X1.end()) beta += i->second;
      }
    }
    
    // Compute constant coefficient alpha
    alpha = cut - (lambda * beta);
    
    // Test of lambda is the breakpoint we are looking for
    //if (alpha0 + (lambda * beta0) == alpha + (lambda * beta)) break;
    if (fabs(alpha0 + (lambda * beta0) - alpha + (lambda * beta)) < 0.0000001) break;
    
    // Assert invariant
    assert(alpha > alpha0 && beta < beta0);
    
    // Compute next lambda
    lambda = (alpha - alpha0) / (beta0 - beta);
  }
  
  if (alpha1) *alpha1 = alpha;
  if (beta1)  *beta1  = beta;
  
  return lambda;
}

// Determines the smallest breakpoint of a given network. Uses the same provisional data
// structures as the other functions in this file.
// Assumptions: All source arcs have constant capacity. All sink arcs have parametric capacity
//              and their constant components are zero. The capacity of all in-between arcs is
//              assumed to be infinite.
double FindSmallestBreakpoint(const set<pair<pair<int,int>,double> > &arcs, int source, int sink)
{
  set<pair<pair<int,int>,double> >::const_iterator i, j;
  double capsum, lambda, beta0;
  
  set<int> X2;
  
  // Compute set of all vertices with inbound edges
  for (i = arcs.begin(); i != arcs.end(); ++i)
  {
    X2.insert(i->first.second);
  }
  
  // Compute sum of all source arc capacities, the largest sink arc capacity and
  // the slope of the line segment before the first breakpoint
  for (i = arcs.begin(), capsum = 0.0, lambda = 0.0, beta0 = 0.0; i != arcs.end(); ++i)
  {
    //if (i->first.first == source) capsum += i->second;
    if (i->first.second == sink)
    {
      lambda += i->second;
      if (X2.find(i->first.first) != X2.end()) beta0 += i->second;
    }
  }
  /*for (i = arcs.begin(); i != arcs.end(); ++i)
  {
    if (i->first.second == sink && i->second < lambda) lambda = i->second;
  }*/
  
  return FindNextBreakpoint(arcs, source, sink, 0, beta0, 0, 0);
  
  /*lambda = (capsum / lambda) + 1;
  
  // Repeat step 1 as often as necessary
  while (true)
  {
    AnalyzeNetwork(arcs, source, sink, 1.0, lambda, 0, &X2, &cut);
    
    beta = 0.0;
    for (i = arcs.begin(); i != arcs.end(); ++i)
    {
      if (i->first.second == sink)
      {
        if (X2.find(i->first.first) != X2.end()) beta += i->second;
      }
    }
    
    alpha = cut - (lambda * beta);
    
    if (lambda * beta0 == alpha + (lambda * beta)) break;
    
    lambda = alpha / (beta0 - beta);
  }
  
  return lambda;*/
}

// Determines the largest breakpoint of a given network and optionally the flow value at that
// breakpoint. Uses the same provisional data structures as the other functions in this file.
// Assumptions: All source arcs have constant capacity. All sink arcs have parametric capacity
//              and their constant components are zero. The capacity of all in-between arcs is
//              assumed to be infinite.
double FindLargestBreakpoint(const set<pair<pair<int,int>,double> > &arcs, int source, int sink, double *maxflow)
{
  set<pair<pair<int,int>,double> >::const_iterator i, j;
  double capsum, lambda, alpha, beta, beta0, alpha0, cut;
  
  set<int> X2;
  
  // Compute set of all vertices with inbound edges
  for (i = arcs.begin(); i != arcs.end(); ++i)
  {
    X2.insert(i->first.second);
  }
  
  // Compute sum of all source arc capacities, the largest sink arc capacity and
  // the slope of the line segment before the first breakpoint
  for (i = arcs.begin(), capsum = 0.0, lambda = 0.0, beta0 = 0.0; i != arcs.end(); ++i)
  {
    if (i->first.first == source) capsum += i->second;
    if (i->first.second == sink)
    {
      lambda += i->second;
      if (X2.find(i->first.first) != X2.end()) beta0 += i->second;
    }
  }
  for (i = arcs.begin(); i != arcs.end(); ++i)
  {
    if (i->first.second == sink && i->second < lambda) lambda = i->second;
  }
  
  lambda = (capsum / lambda) + 1;
  
  // Determine maximum flow value as alpha0
  AnalyzeNetwork(arcs, source, sink, 1.0, lambda, 0, 0, &alpha0);
  
  // Determine a lambda below the largest breakpoint
  lambda = alpha0 / beta0;
  
  while (true)
  {
    AnalyzeNetwork(arcs, source, sink, 1.0, lambda, 0, &X2, &cut);
    
    beta = 0.0;
    for (i = arcs.begin(); i != arcs.end(); ++i)
    {
      if (i->first.second == sink)
      {
        if (X2.find(i->first.first) != X2.end()) beta += i->second;
      }
    }
    
    alpha = cut - (lambda * beta);
    
    if (alpha + (lambda * beta) == alpha0) break;
    
    lambda = (alpha0 - alpha) / beta;
  }
  
  if (maxflow) *maxflow = alpha0;
  
  return lambda;
}

void print_cut(const char *caption, const set<int> &cut)
{
  set<int>::const_iterator i;
  printf("Cut %s:", caption);
  for (i = cut.begin(); i != cut.end(); ++i) printf(" %d", *i);
  printf("\n");
}

void print_contract(const char *caption, const set<int> &vertices, int new_id)
{
  set<int>::const_iterator i;
  printf("Contraction %s:", caption);
  for (i = vertices.begin(); i != vertices.end(); ++i) printf(" %d", *i);
  printf(" --> %d\n", new_id);
}

/*void slice(set<pair<pair<int,int>,double> > &arcs, int source, int sink, double lambda1, double lambda3, set<double> &breakpoints)
{
  double cut1, cut3, lambda2, c0, v0, c1, v1;
  set<pair<pair<int,int>,double> > temp;
  set<pair<pair<int,int>,double> >::const_iterator i;
  set<int> X2min, X2max, not_X2;
  
  printf("+++ slice %f %f\n", lambda1, lambda3);
  
  printf("\n>>> network\n");
  for (i = arcs.begin(); i != arcs.end(); ++i)
  {
    printf("   %2d %2d     %f\n", i->first.first, i->first.second, i->second);
  }
  printf("\n");
  
  // Step S1
  // Compute cut costs for ({s}, V - {s}) and (V - {t}, {t})
  for (cut1 = 0.0, cut3 = 0.0, i = arcs.begin(); i != arcs.end(); ++i)
  {
    if (i->first.first == source) cut1 += i->second;
    if (i->first.second == sink) cut3 += i->second;
  }
  
  lambda2 = cut1 / cut3;
  
  // Step S2
  AnalyzeNetwork(arcs, source, sink, 1.0, lambda2, &X2max, &X2min, 0);
  
  printf("lambda2 = %f\n", lambda2);
  
  print_cut("X2 ", X2min);
  print_cut("X2'", X2max);
  
  // Step S3
  for (i = arcs.begin(), c0 = 0.0, v0 = 0.0, c1 = 0.0, v1 = 0.0; i != arcs.end(); ++i)
  {
    assert(!(X2min.find(i->first.second) != X2min.end() && X2min.find(i->first.first) == X2min.end()));
    
    if (X2min.find(i->first.first) != X2min.end() && X2min.find(i->first.second) == X2min.end())
    {
      if (i->first.first == source) v0 += i->second;
      else c0 += i->second;
    }
    
    if (X2max.find(i->first.first) != X2max.end() && X2max.find(i->first.second) == X2max.end())
    {
      if (i->first.first == source) v1 += i->second;
      else c1 += i->second;
    }
  }
  
  if (c0 != c1 || v0 != v1)
  {
    breakpoints.insert(lambda2);
    printf("added breakpoint %f\n", lambda2);
  }
  
  // Step S4
  // Case X2 != {s}
  if ((X2min.size() != 1 || *X2min.begin() != source))
  {
    not_X2.clear();
    for (i = arcs.begin(); i != arcs.end(); ++i)
    {
      if (X2min.find(i->first.first) == X2min.end()) not_X2.insert(i->first.first);
      if (X2min.find(i->first.second) == X2min.end()) not_X2.insert(i->first.second);
    }
    
    print_contract("!X2", not_X2, ((not_X2.find(sink) == not_X2.end()) ? *not_X2.begin() : sink));
    
    temp.clear();
    ContractVertices(arcs, temp, not_X2, ((not_X2.find(sink) == not_X2.end()) ? *not_X2.begin() : sink));
    
    slice(temp, source, sink, lambda1, lambda2, breakpoints);
  }
  
  // Case !X2' != {t}
  not_X2.clear();
  for (i = arcs.begin(); i != arcs.end(); ++i)
  {
    if (X2max.find(i->first.first) == X2max.end()) not_X2.insert(i->first.first);
    if (X2max.find(i->first.second) == X2max.end()) not_X2.insert(i->first.second);
  }

  if ((not_X2.size() != 1 || *not_X2.begin() != sink))
  {
    print_contract("X2'", X2max, ((X2max.find(source) == X2max.end()) ? *X2max.begin() : source));
    
    temp.clear();
    ContractVertices(arcs, temp, X2max, ((X2max.find(source) == X2max.end()) ? *X2max.begin() : source));
    
    slice(temp, source, sink, lambda2, lambda3, breakpoints);
  }
  
  printf("--- slice %f %f\n", lambda1, lambda3);
}*/

// Determines all breakpoints of a given network. Uses the same provisional data structures as
// the other functions in this file.
// Assumptions: All source arcs have constant capacity. All sink arcs have parametric capacity
//              and their constant components are zero. The capacity of all in-between arcs is
//              assumed to be infinite.
void FindBreakpoints(const set<pair<pair<int,int>,double> > &arcs, int source, int sink, set<double> &breakpoints)
{
  /*double lambda1, lambda3;
  //set<int> X1, not_X1, X3;
  //set<pair<pair<int,int>,double> > temp, contracted;
  set<pair<pair<int,int>,double> >::const_iterator i;
  
  // Step 1: compute smallest and largest breakpoints
  lambda1 = FindSmallestBreakpoint(arcs, source, sink);
  lambda3 = FindLargestBreakpoint(arcs, source, sink, 0);
  breakpoints.insert(lambda1);
  breakpoints.insert(lambda3);
  
  if (lambda1 == lambda3) return;
  
  AnalyzeNetwork(arcs, source, sink, 1.0, lambda1, 0, &X1, &cost1);
  AnalyzeNetwork(arcs, source, sink, 1.0, lambda3, 0, &X3, &cost3);
  
  X1.clear();
  for (i = arcs.begin(); i != arcs.end(); ++i)
  {
    X1.insert(i->first.first);
  }
  
  slice(arcs, source, sink, lambda1, lambda3, X1, X3, cost1, cost3, breakpoints);*/
  
  set<pair<pair<int,int>,double> >::const_iterator i;
  double alpha0, beta0, alpha1, beta1, lambda;
  
  set<int> X2;
  
  // Compute set of all vertices with inbound edges
  for (i = arcs.begin(); i != arcs.end(); ++i)
  {
    X2.insert(i->first.second);
  }
  
  // Compute the slope of the line segment before the first breakpoint
  alpha0 = 0.0;
  for (i = arcs.begin(), beta0 = 0.0; i != arcs.end(); ++i)
  {
    if (i->first.second == sink && X2.find(i->first.first) != X2.end()) beta0 += i->second;
  }
  
  while (true)
  {
    lambda = FindNextBreakpoint(arcs, source, sink, alpha0, beta0, &alpha1, &beta1);
    breakpoints.insert(lambda);
    //if (alpha0 == alpha1) break;
    if (beta1 == 0.0) break;
    alpha0 = alpha1;
    beta0 = beta1;
  }
  
  // Compute cuts (X1,!X1) and (X3,!X3)
  /*AnalyzeNetwork(arcs, source, sink, 1.0, lambda1, &X1, 0, 0);
  AnalyzeNetwork(arcs, source, sink, 1.0, lambda3, 0, &X3, 0);
  
  // Compute !X1 from X1
  for (i = arcs.begin(); i != arcs.end(); ++i)
  {
    if (X1.find(i->first.first) == X1.end()) not_X1.insert(i->first.first);
    if (X1.find(i->first.second) == X1.end()) not_X1.insert(i->first.second);
  }
  
  // Create G' by contracting the vertices in X3 and !X1
  ContractVertices(arcs, temp, X3, source);
  ContractVertices(temp, contracted, not_X1, sink);
  temp.clear();
  
  X1.clear();
  X3.clear();
  AnalyzeNetwork(contracted, source, sink, 1.0, lambda1, &X1, &X3, 0);
  print_cut("1: S_max", X1);
  print_cut("1: T_max", X3);
  X1.clear();
  X3.clear();
  AnalyzeNetwork(contracted, source, sink, 1.0, lambda3, &X1, &X3, 0);
  print_cut("3: S_max", X1);
  print_cut("3: T_max", X3);
  
  slice(contracted, source, sink, lambda1, lambda3, breakpoints);*/
}

// Create a network with nps nodes connected to source and the same number connected to sink. Each of
// those nodes has exactly one arc connected to a node on the other side such that every node has only
// one arc connected (not counting source/sink arcs). All source arcs have the same capacity, sink arcs
// have ascending capacity. The source of the network will have ID 1. The ID of the sink is returned.
int make_network(set<pair<pair<int,int>,double> > &arcs, int nps)
{
  int n;
  
  for (n = 0; n < nps; ++n)
  {
    arcs.insert(make_pair(make_pair(1, n + 2), (double)nps));
    arcs.insert(make_pair(make_pair(n + 2, n + nps + 2), (double)(nps * nps * nps)));
    arcs.insert(make_pair(make_pair(n + nps + 2, nps + nps + 2), (double)(n + 1)));
  }
  
  return (2 * nps) + 2;
}

#if 0
int main()
{
  set<pair<pair<int,int>,double> > arcs;
  set<pair<pair<int,int>,double> >::iterator ai;
  set<double> bp;
  set<double>::iterator i;
  set<int> Xs, Xt;
  double lambda;
  int sink;
  
  // generate network with 10 nodes on each side (s + 10 + 10 + t = 22 total in the maxflow problem)
  sink = make_network(arcs, 2);
  
  for (ai = arcs.begin(); ai != arcs.end(); ++ai)
  {
    printf("%3d ==%4.0f ==> %3d\n", ai->first.first, ai->second, ai->first.second);
  }
  
  FindBreakpoints(arcs, 1, sink, bp);
  printf("Breakpoints:");
  for (i = bp.begin(); i != bp.end(); ++i) printf(" %f", *i);
  printf("\n");
  
  printf("\nMinimum cuts (X,¬X) for each breakpoint:\n");
  for (i = bp.begin(); i != bp.end(); ++i)
  {
    AnalyzeNetwork(arcs, 1, sink, 1.0, *i, &Xs, &Xt, 0);
    
    printf(">>> %f\n", *i);
    print_cut("Elements in X for max| X|", Xs);
    print_cut("Elements in X for max|¬X|", Xt);
    printf("\n");
  }
  
  return 0;
}
#endif
