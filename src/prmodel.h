/*****************************************************************************/
/*!
 *   Copyright 2009 Jonathan Bogdoll, Holger Hermanns, Lijun Zhang
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



#ifndef _PRMODEL_H_
#define _PRMODEL_H_

#include <stdio.h>

class ProbabilisticModel
{
public:
  ProbabilisticModel(bool c) { continuous = c; }
  virtual ~ProbabilisticModel() {}
  
  virtual void Parse(FILE* = 0, bool = false) = 0;
  virtual void Write(FILE* = stdout) = 0;
  
  virtual int States() = 0;
  virtual int Transitions() = 0;
  virtual int Distributions() = 0;
  virtual double PrecisionThreshold() = 0;
  
  typedef enum {MC, PA} ModelType;
  
  ModelType Type() { return mt; }
  bool ContinuousTimeModel() { return continuous; }
  
protected:
  ModelType mt;
  bool continuous;
};

class MarkovChain : public ProbabilisticModel
{
public:
  MarkovChain(bool = false);
  virtual ~MarkovChain();
  
  int n; // num states
  long nnz; // num non zeros
  
  double *non_zeros;
  int *cols;
  int *row_starts;

  void Parse(FILE* = 0, bool = false);
  void Write(FILE* = stdout);
  int States() { return n; }
  int Transitions() { return nnz; }
  int Distributions() { return n; }
  double PrecisionThreshold();

private:
  void InitMatrix();
};

class ProbabilisticAutomaton : public ProbabilisticModel
{
public:
  ProbabilisticAutomaton(bool = false);
  ~ProbabilisticAutomaton();
  
  int n; // number of states
  int na; // number of actions
  int nnz; // number of transitions
  int da; // number of different types of actions
  
  int *state_starts;
  int *actions;
  int *row_starts;
  double *non_zeros;
  int *cols;
  int *atable;
  
  int Action(int i) { return atable[i]; }
  
  void Parse(FILE* = 0, bool = false);
  void Write(FILE* = stdout);
  int States() { return n; }
  int Transitions() { return nnz; }
  int Distributions() { return na; }
  double PrecisionThreshold();
  
private:
  void InitMatrix();
};

#endif//_PRMODEL_H_
