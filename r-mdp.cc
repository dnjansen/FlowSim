#include <stdio.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <set>
#include "prmodel.h"

/*
 * Generate a random PA.
 * n - number of states
 * a - minimum number of successors per action (not state!)
 * b - maximum number of successors per action
 * ac- total number of different actions
 * ma- minimum of actions outgoing from a state (at least 1)
 * Ma- maximum of actions outgoing from a state
 * r - likelihood of all actions going out from a state being the same
 *     (0.0 = all actions are different, 1.0 = all actions are the same)
 * lbias - linearity bias: 0.0=uniformly random, 1.0=acyclic model
 *
 * Probabilities are uniformly distributed between successor states and
 * P(s,S)=1 for all states.
 */
ProbabilisticAutomaton *RandomPA(int n, int a, int b, int ac, int ma, int Ma, double r, double lbias)
{
  int s, i, j, k, na, da, row, action, transitions, linearchoices;
  std::vector<int> row_starts, state_starts, actions, cols;
  std::set<int> successors;
  std::set<int>::iterator si;
  ProbabilisticAutomaton *mdp = 0;
  
  if (Ma < ma || b < a || a < 0) return 0;
  
  srand(time(0));
  
  row = 0;
  action = 0;
  
  // Generate actions
  for (i = 0; i < n; ++i)
  {
    state_starts.push_back(action);
    if (Ma == ma) na = ma;
    else na = ma + (rand() % (Ma - ma + 1));
    
    for (k = 0; k < na; ++k)
    {
      if (k == 0)
      {
        da = rand() % ac;
        actions.push_back(da);
      }
      else
      {
        actions.push_back((da + int((rand() % ac) * (1.0 - r))) % ac);
      }
      ++action;
      
      // For each action, generate a set of successors
      row_starts.push_back(row);
      if (a == b) transitions = a;
      else transitions = a + (rand() % (b - a + 1));
      linearchoices = n - i - 1;
      if (transitions > linearchoices) transitions = linearchoices + (int)ceil((1.0 - lbias) * (transitions - linearchoices));
      for (j = 0; j < transitions; ++j)
      {
        if (rand() > lbias * RAND_MAX || linearchoices <= 0)
        {
          do
          {
            s = rand() % n;
          }
          while (successors.find(s) != successors.end());
          successors.insert(s);
          if (s >= i) --linearchoices;
        }
        else
        {
          if (i == n - 1) s = i;
          else do
          {
            s = i + (rand() % (n - i));
          }
          while (successors.find(s) != successors.end());
          successors.insert(s);
          --linearchoices;
        }
      }
      row += (int)successors.size();
      for (si = successors.begin(); si != successors.end(); ++si) cols.push_back(*si);
      successors.clear();
    }
  }
  
  state_starts.push_back(action);
  row_starts.push_back(row);
  
  // Initialize target MDP structure
  mdp = new ProbabilisticAutomaton;
  
  mdp->n = n;
  mdp->na = actions.size();
  mdp->nnz = row_starts[row_starts.size() - 1];
  mdp->da = ac;
  
  mdp->state_starts = new int[state_starts.size()];
  mdp->row_starts = new int[row_starts.size()];
  mdp->non_zeros = new double[cols.size()];
  mdp->cols = new int[cols.size()];
  mdp->actions = new int[actions.size()];
  mdp->atable = new int[mdp->da];
  
  for (i = 0; i < (int)state_starts.size(); ++i)
  {
    *(mdp->state_starts + i) = state_starts[i];
  }
  
  // Generate transition probabilities (uniformly distributed) for the transitions defined above
  for (i = 0, k = 0; i < (int)row_starts.size() - 1; ++i)
  {
    transitions = row_starts[i+1] - row_starts[i];
    for (j = 0; j < transitions; ++j, ++k)
    {
      *(mdp->cols + k) = cols[k];
      *(mdp->non_zeros + k) = (1.0 / transitions);
    }
  }
  
  for (i = 0; i < (int)row_starts.size(); ++i)
  {
    *(mdp->row_starts + i) = row_starts[i];
  }
  
  for (i = 0; i < (int)actions.size(); ++i)
  {
    *(mdp->actions + i) = actions[i];
  }
  
  for (i = 0; i < mdp->da; ++i) mdp->atable[i] = i;
  
  return mdp;
}

// Compile main function only if r-mdp is compiled as a standalone program and not linked elsewhere
#ifdef _STANDALONE_
int main(int argc, char *argv[])
{
  if (argc < 5 || argc > 9)
  {
    printf("Usage: r-mdp N A B NA [min-A [max-A [R [Lbias]]]]\n");
    printf("\n");
    printf("    N  number of states, >0\n");
    printf("    A  minimum number of successors per action (not state!), A >= 0\n");
    printf("    B  maximum number of successors per action, B >= A\n");
    printf("   NA  total number of different actions, >0\n");
    printf("min-A  minimum number of actions outgoing from each state, >0, default 1\n");
    printf("max-A  maximum number of actions outgoing from each state, >=min-A, default max(min-A,NA)\n");
    printf("    R  likelihood of all actions outgoing from a state being the same,\n");
    printf("       i.e. R=0.0: actions are random, R=1.0: all actions the same, default 0.0\n");
    printf("Lbias  linearity bias, 0.0=uniformly random, 1.0=acyclic model, default 0.0\n");
    printf("\n");
    return -1;
  }
  
  int n = strtol(argv[1], 0, 10);
  int a = strtol(argv[2], 0, 10);
  int b = strtol(argv[3], 0, 10);
  int na = strtol(argv[4], 0, 10);
  int mina = 1, maxa = (mina < na ? na : mina);
  double r = 0.0, lbias = 0.0;
  
  if (argc >= 6) mina = strtol(argv[5], 0, 10);
  if (argc >= 7) maxa = strtol(argv[6], 0, 10);
  if (argc >= 8) r = strtod(argv[7], 0);
  if (argc >= 9) lbias = strtod(argv[8], 0);
  
  if (n <= 0 || a < 0 || b < a || na <= 0 || mina < 0 || maxa < mina || r < 0.0 || r > 1.0 || lbias < 0.0 || lbias > 1.0)
  {
    printf("Usage: r-mdp N A B NA [min-A [max-A [R [Lbias]]]]\n");
    printf("\n");
    printf("    N  number of states, >0\n");
    printf("    A  minimum number of successors per action (not state!), A >= 0\n");
    printf("    B  maximum number of successors per action, B >= A\n");
    printf("   NA  total number of different actions, >0\n");
    printf("min-A  minimum number of actions outgoing from each state, >=0, default 1\n");
    printf("max-A  maximum number of actions outgoing from each state, >=min-A, default max(min-A,NA)\n");
    printf("    R  likelihood of all actions outgoing from a state being the same,\n");
    printf("       i.e. R=0.0: all actions different, R=1.0: all actions the same, default 0.5\n");
    printf("Lbias  linearity bias, 0.0=uniformly random, 1.0=acyclic model, default 0.0\n");
    printf("\n");
    return -1;
  }
  
  ProbabilisticAutomaton *m = RandomPA(n, a, b, na, mina, maxa, r, lbias);
  
  m->Write();
  
  delete m;
  
  return 0;
}
#endif
