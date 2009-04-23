#include "Weak.h"
#include "Weak_auxiliary.cc"

#include <set>
#include <utility>

using namespace std;

int formaterror(const char *s, int line)
{
  fprintf(stderr, "format error in line %d: %s\n", line, s);
  return -1;
}

int main(int argc, char *argv[])
{
  int nodes, source, sink, numarcs = 0, linenum, tail, head;
  double cap;
  char line[384];
  bool have_p = false, have_s = false, have_t = false;
  
  set<pair<pair<int,int>,double> > arcs;
  set<double> bp;
  set<double>::iterator i;
  set<int> Xs, Xt;
  
  linenum = 0;
  while (1)
  {
    fgets(&line[0], 384, stdin);
    if (feof(stdin) || !strcmp(&line[0], "\n")) break;
    ++linenum;
    
    if (!strncmp(&line[0], "c ", 2)) continue;
    
    if (!have_p)
    {
      if (sscanf(&line[0], "p max %d %d\n", &nodes, &numarcs) < 2) return formaterror("problem definition [p max %d %d] not found", linenum);
      have_p = true;
    }
    else if (!have_s)
    {
      if (sscanf(&line[0], "n %d s\n", &source) < 1) return formaterror("source node defintion not found", linenum);
      if (source != 1) fprintf(stderr, "line %d: warning: source node should be \"1\" but is %d\n", linenum, source);
      have_s = true;
    }
    else if (!have_t)
    {
      if (sscanf(&line[0], "n %d t\n", &sink) < 1) return formaterror("sink node defintion not found", linenum);
      if (sink != nodes) fprintf(stderr, "line %d: warning: sink node should be \"%d\" but is %d\n", linenum, nodes, sink);
      have_t = true;
    }
    else
    {
      if (numarcs == 0) return formaterror("already parsed all arcs but still got input", linenum);
      if (sscanf(&line[0], "a %d %d %le\n", &tail, &head, &cap) < 3) return formaterror("couldn't parse arc definition", linenum);
      arcs.insert(make_pair(make_pair(tail, head), cap));
      --numarcs;
    }
  }
  
  if (numarcs != 0)
  {
    fprintf(stderr, "warning: parsed less than specified number of arcs\n");
    return -1;
  }
  
  FindBreakpoints(arcs, source, sink, bp);
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
