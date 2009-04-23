#define WITH_VERIFIER

#include <stdio.h>
#include "prmodel.h"
#include "Strong.h"

#include "StrongMC.cc"
#include "StrongPA.cc"

void load_relation(FILE *f, std::set<std::pair<int,int> > &rel)
{
  char *buf, *p, *q;
  long size;
  int x, y;
  
  fseek(f, 0, SEEK_END);
  size = ftell(f);
  rewind(f);
  
  buf = new char[size + 1];
  fread(buf, 1, size, f);
  buf[size] = 0;
  
  p = buf, q = buf;
  while (*p)
  {
    if (*p >= '0' && *p <= '9') *(q++) = *p;
    else *(q++) = ' ';
    ++p;
  }
  *q = 0;
  
  p = buf;
  while (*p)
  {
    x = strtol(p, &p, 10);
    y = strtol(p, &p, 10);
    rel.insert(std::make_pair(x, y));
    while (*p == ' ') ++p;
  }
  
  delete [] buf;
}

int LabelFunction(void *userdata, int s)
{
  int *param = (int*)userdata;
  
  if (s >= 0 && s < param[0])
  {
    if (param[1] <= 1) return 0;
    return s % param[1];
  }
  
  return param[1];
}

int main(int argc, char *argv[])
{
  double epsilon;
  unsigned long type, labels;
  int label_data[2];
  bool continuous_model;
  FILE *f, *rel;
  ProbabilisticModel *pm;
  StrongSimulation *ss;
  std::set<std::pair<int,int> > hypothesis, false_pos, false_neg;
  std::set<std::pair<int,int> >::iterator ri;

  if (argc < 6 || argc > 7)
  {
    fprintf(stderr, "Usage: %s type labels epsilon model relation\n\n", argv[0]);
    fprintf(stderr, "Test if a given relation is a strong simulation relation on a given model.\n");
    fprintf(stderr, "type     - model type, one of dtmc, ctmc, pa, cpa\n");
    fprintf(stderr, "labels   - number of different labels to assign\n");
    fprintf(stderr, "epsilon  - floating-point approximation threshold (-1 = automatic)\n");
    fprintf(stderr, "model    - input model filename. may be \"-\" to read from stdin\n");
    fprintf(stderr, "relation - read the relation from this file\n\n");
    return 0;
  }

  if (!strcmp(argv[1], "dtmc")) type = 1;
  else if (!strcmp(argv[1], "ctmc")) type = 2;
  else if (!strcmp(argv[1], "pa")) type = 3;
  else if (!strcmp(argv[1], "cpa")) type = 4;
  else
  {
    fprintf(stderr, "Model type '%s' not recognized\n", argv[1]);
    return -1;
  }

  labels = strtoul(argv[2], 0, 10);
  epsilon = strtod(argv[3], 0);
  f = ((!strcmp(argv[4], "-")) ? stdin : fopen(argv[4], "rb"));
  if (!f)
  {
    fprintf(stderr, "Failed to open model '%s'\n", argv[4]);
    return -1;
  }
  
  rel = fopen(argv[5], "rb");
  if (!rel)
  {
    fprintf(stderr, "Failed to open relation file '%s'\n", argv[5]);
    return -1;
  }

#ifdef OPT_QUOTIENT
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
  ss = new StrongSimulation_Quotient;
#else//OPT_QUOTIENT
  switch (type)
  {
  case 1:
    pm = new MarkovChain;
    ss = new StrongSimulation_MC;
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
#endif//OPT_QUOTIENT

  pm->Parse(f, continuous_model);
  if (strcmp(argv[5], "-")) fclose(f);

  label_data[0] = pm->States();
  label_data[1] = labels;
  ss->SetLabelFunction(LabelFunction, (void*)&label_data[0]);
  if (epsilon < 0) ss->SetFPPrecision(pm->PrecisionThreshold());
  else ss->SetFPPrecision(epsilon);
  
  load_relation(rel, hypothesis);
  fclose(rel);
  
  ss->Verify(pm, hypothesis, &false_pos, &false_neg);
  
  if (false_pos.size() + false_neg.size() == 0) printf("The given relation is a strong simulation relation for the given model.\n");
  else
  {
    for (ri = false_pos.begin(); ri != false_pos.end(); ++ri) printf("-%5d %5d\n", ri->first, ri->second);
    for (ri = false_neg.begin(); ri != false_neg.end(); ++ri) printf("+%5d %5d\n", ri->first, ri->second);
  }
  
  delete pm;
  delete ss;

  return 0;
}
