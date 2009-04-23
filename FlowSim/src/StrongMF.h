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



#ifndef STRONG_MF_H
#define STRONG_MF_H

#include <set>
#include <vector>
#include "prmodel.h"
#include "maxflow.h"

// Decide strong simulation using maxflow algorithm without any optimizations
class StrongMF
{
public:
  StrongMF(MarkovChain* r);
  virtual ~StrongMF();

  void run();  // the main loop of the simulation procedure
  
protected:
  // One pair of the relation. Stores the two states being related together with their maxflow problem
  class Pair
  {
  public:
    Pair() { x=0, y=0, simulation=0; required_flow = 1.0; }
    ~Pair() {}
    int x,y;
    ParametricMaxflow<double> *simulation;
    double required_flow;
  };
  
  bool debug;
  void InitializeLabels();
  
  int BuildRelationMap(bool**);
  
  void InitializeRelation();
  bool StrongSimulation(Pair*);
  void ConstructNetwork(Pair*, int, int, bool*);

  double *non_zeros; // list of probabilities for transitions
  int *cols;         // list of successor states
  int *row_starts;   // list indices
  int n_states;      // number of states
  
  int *labels, n_labels;
  int *label_row_starts;
  
  //number of pairs in relation
  int size_of_relation;

  //The array of the relation (R).
  Pair *relation;
  
  bool *relation_map;
};

#endif//STRONG_MF_H
