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



#include "StrongMF.h"
#include "maxflow.cc"

#define epsilon 0.000000001

StrongMF::StrongMF(MarkovChain* r)
{
#ifdef DEBUG
  debug = true;
#else//DEBUG
  debug = false;
#endif//DEBUG
  
  relation = NULL;
  
  non_zeros = r->non_zeros;
  cols = r->cols;
  row_starts = r->row_starts;
  
  labels = 0;
  label_row_starts = 0;
  
  n_states = r->n;
  
  size_of_relation = 0;
  
  InitializeLabels();   // FIXME Set arbitrary labels for testing, this will be changed later
  InitializeRelation();
}

StrongMF::~StrongMF()
{
  if(relation != NULL) delete [] relation;
  
  if(labels != NULL)
    delete [] labels;
  if(label_row_starts != NULL)
    delete [] label_row_starts;
  delete[] relation_map;
}

// Initialize labels in an arbitrary way. This is for testing only.
void StrongMF::InitializeLabels()
{
  int n;
  
  n_labels = n_states;
  
  labels = new int[n_labels];
  label_row_starts = new int[n_states];
  
  for (n = 0; n < n_states; ++n)
  {
    *(label_row_starts + n) = n;
    *(labels + n) = n % 3;
  }
}

int StrongMF::BuildRelationMap(bool **ptr_relmap)
{
  std::set<int> distinct_labels;
  std::set<int>::const_iterator label;
  int n, m, s, size;
  int *label_group = new int[n_states]; // list of states which have a particular label

  bool *relmap = new bool[n_states * n_states];
    
  // Generate a set of all labels found in this problem
  for (n = 0; n < n_labels; ++n) distinct_labels.insert(*(labels + n));
  
  // Iterate over labels and generate a matrix for the relation while counting the size
  size = 0;
  memset(relmap, 0, sizeof(bool) * n_states * n_states);
  for (label = distinct_labels.begin(); label != distinct_labels.end(); label++)
  {
    memset(label_group, 0, sizeof(int) * n_states);
    for (m = 0, s = 0; m < n_labels; ++m) // Find all states which have the current label set
    {
      if (m == *(label_row_starts + s + 1)) ++s;
      if (*(labels + m) == *label)
      {
        *(label_group + s) = 1;
        if (s == n_states - 1) break;
        m = *(label_row_starts + s + 1) - 1;
      }
    }
    for (m = 0; m < n_states; ++m) // Add pairs to relation belonging to the current label
    {
      if (*(label_group + m))
      {
        for (s = m; s < n_states; ++s)
        {
          if (*(label_group + s)) // Pairs (m,s) and (s,m) share at least one label, add them.
          {
            if (!*(relmap + (m * n_states) + s))
            {
              *(relmap + (m * n_states) + s) = true;
              ++size;
            }
            if (!*(relmap + (s * n_states) + m))
            {
              *(relmap + (s * n_states) + m) = true;
              ++size;
            }
          }
        }
      }
    }
  }
  
  delete[] label_group;
  
  if (ptr_relmap) *ptr_relmap = relmap;
  else delete[] relmap;
  
  return size;  
}

/**
* Corresponds the SimulationRelation
*/
void StrongMF::run()
{
  //index of the new relations
  std::vector<int> nr;
  int iterations = 1, i;
  bool *newrel = new bool[n_states * n_states];
  
  do
  {
    nr.clear();
    
    if(debug) printf("Start of iteration %d, the relation size is %d\n", iterations, size_of_relation);
  
    memcpy(newrel, relation_map, sizeof(bool) * n_states * n_states);
    
    for(i=0; i<size_of_relation; i++)
    {
      if (StrongSimulation(relation + i)) nr.push_back(i);
      else *(newrel + ((relation + i)->x * n_states) + (relation + i)->y) = false;
    }
  
    memcpy(relation_map, newrel, sizeof(bool) * n_states * n_states);
    
    if(debug) printf("\nEnd of iteration %d, the new relation size is %d\n", iterations, nr.size());
    
    iterations++; 

    if(nr.size() == (unsigned int)size_of_relation) break;
    
    size_of_relation = nr.size();
    Pair* new_relation = new Pair[size_of_relation];
    for(int i=0; i<size_of_relation; i++)
    {
      new_relation[i] = relation[nr[i]];
    }
    
    delete [] relation;
    relation = new_relation;
  }
  while (nr.size() > 0);
  
  delete [] newrel;
  
  if (debug)
  {
    FILE *f = fopen("mf.relation", "wb");
    for(int i=0; i<size_of_relation; i++) if (relation[i].x != relation[i].y) fprintf(f, "(%d,%d)\n", relation[i].x, relation[i].y);
    fclose(f);
  }
}

//the parameter num_par is the number of partitions of the states space.
void StrongMF::InitializeRelation()
{
  int n, m, s;
  
  size_of_relation = BuildRelationMap(&relation_map);
    
  // Generate array of pairs now that we know how many pairs there are total and generate
  // maxflow problems.
  relation = new Pair[size_of_relation];
  for (n = 0, s = 0; n < n_states; ++n)
  {
    for (m = 0; m < n_states; ++m)
    {
      if (*(relation_map + (n * n_states) + m))
      {
        assert(s < size_of_relation);
        (relation + s)->x = n;
        (relation + s)->y = m;
        (relation + s)->simulation = 0;
        ++s;
      }
    }
  }
  
  assert(s == size_of_relation);
  
  delete[] labels;
  delete[] label_row_starts;
  labels = 0;
  label_row_starts = 0;
}

void StrongMF::ConstructNetwork(Pair *p, int s1, int s2, bool *relation_map)
{
  double bottom1, bottom2;
  
  int s1_suc = row_starts[s1+1] - row_starts[s1]; // Get number of successor states
  int s2_suc = row_starts[s2+1] - row_starts[s2]; //
  int n, m, nodes;
  
  ParametricMaxflow<double> *simulation = new ParametricMaxflow<double>;
  
  simulation->SetOptimize(false);
  
  nodes = 1 + s1_suc + 2 + s2_suc + 1; // source, sink, successor states plus two auxiliary states
  
  // Network must be created from scratch
  simulation->BeginNetwork(nodes);
  
  // Add arcs to the successors of s1 as well as the arcs between the two sets of successors
  for (n = 0, bottom1 = 1.0; n < s1_suc; ++n)
  {
    simulation->AddArc(simulation->Source(), cols[row_starts[s1] + n], non_zeros[row_starts[s1] + n]);
    bottom1 -= non_zeros[row_starts[s1] + n];
    
    // Add "infinite" capacity arcs between successor states that are within the relation
    for (m = 0; m < s2_suc; ++m)
    {
      if (*(relation_map + (cols[row_starts[s1] + n] * n_states) + cols[row_starts[s2] + m]) || cols[row_starts[s1] + n] == cols[row_starts[s2] + m])
      {
        simulation->AddArc(cols[row_starts[s1] + n], cols[row_starts[s2] + m] + n_states, 1.0);
      }
    }
  }
  
  // Add arcs to the successors of s2
  for (n = 0, bottom2 = 1.0; n < s2_suc; ++n)
  {
    simulation->AddArc(cols[row_starts[s2] + n] + n_states, simulation->Sink(), non_zeros[row_starts[s2] + n]);
    bottom2 -= non_zeros[row_starts[s2] + n];
  }
  
  // Set required flow. This value may become greater than 1, in which case the simulation condition is
  // unfulfillable for this network.
  //bottom1 -= bottom2;
  p->required_flow = 1.0 - bottom1;
  
  simulation->EndNetwork();
  
  p->simulation = simulation;
}

bool StrongMF::StrongSimulation(Pair *p)
{
  double flow;
  
  if (p->x == p->y) return true;
  
  if (row_starts[p->x] == row_starts[p->x + 1]) return true;
  if (row_starts[p->y] == row_starts[p->y + 1]) return false;
  
  ConstructNetwork(p, p->x, p->y, relation_map);
  
  flow = p->simulation->Flow();
  
  delete p->simulation;
  p->simulation = 0;
  
  if (flow > p->required_flow - epsilon && flow < p->required_flow + epsilon) return true;
  
  return false;
}
