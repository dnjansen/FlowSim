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
