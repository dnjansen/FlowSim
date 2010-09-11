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



#include "benchmark.h"
#include <sys/resource.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <libgen.h>
#include <vector>
#include <utility>

#include "bench.cc"
#include "compactmaxflow.cc"

char myself[128];

// Parses random model parameters from a file; the parameters are assumed
// to be space-separated and in the same order as expected by the r-dtmc tool
// plus a number of additional parameters
bool ParseRandom(FILE *f, RandomModel *prm)
{
  long pos, size;
  char *buffer, *p;
  
  pos = ftell(f);
  fseek(f, 0, SEEK_END);
  size = ftell(f) - pos;
  fseek(f, pos, SEEK_SET);
  if (size > 512) size = 512;
  
  buffer = new char[size];
  fread(buffer, 1, size, f);
  
  prm->xtarget = 0xffff;
  prm->ztarget = 0xffff;
  
  p = buffer;
  while (*p == ' ' || *p == '\t') ++p;
  if (!strncmp(p, "$x", 2)) p += 2, prm->xtarget = 0x0001;
  else if (!strncmp(p, "$z", 2)) p += 2, prm->ztarget = 0x0001;
  else prm->n = strtol(p, &p, 10);
  while (*p == ' ' || *p == '\t') ++p;
  if (!strncmp(p, "$x", 2)) p += 2, prm->xtarget = 0x0002;
  else if (!strncmp(p, "$z", 2)) p += 2, prm->ztarget = 0x0002;
  else prm->a = strtol(p, &p, 10);
  while (*p == ' ' || *p == '\t') ++p;
  if (!strncmp(p, "$x", 2)) p += 2, prm->xtarget = 0x0003;
  else if (!strncmp(p, "$z", 2)) p += 2, prm->ztarget = 0x0003;
  else prm->b = strtol(p, &p, 10);
  while (*p == ' ' || *p == '\t') ++p;
  if (!strncmp(p, "$x", 2)) p += 2, prm->xtarget = 0x0004;
  else if (!strncmp(p, "$z", 2)) p += 2, prm->ztarget = 0x0004;
  else prm->c = strtol(p, &p, 10);
  while (*p == ' ' || *p == '\t') ++p;
  if (!strncmp(p, "$x", 2)) p += 2, prm->xtarget = 0x0100;
  else if (!strncmp(p, "$z", 2)) p += 2, prm->ztarget = 0x0100;
  else prm->fb = strtod(p, &p);
  while (*p == ' ' || *p == '\t') ++p;
  if (!strncmp(p, "$x", 2)) p += 2, prm->xtarget = 0x0200;
  else if (!strncmp(p, "$z", 2)) p += 2, prm->ztarget = 0x0200;
  else prm->lb = strtod(p, &p);
  while (*p == ' ' || *p == '\t') ++p;
  if (!strncmp(p, "$x", 2)) p += 2, prm->xtarget = 0x0300;
  else if (!strncmp(p, "$z", 2)) p += 2, prm->ztarget = 0x0300;
  else prm->cb = strtod(p, &p);
  while (*p == ' ' || *p == '\t') ++p;
  if (!strncmp(p, "$x", 2)) p += 2, prm->xtarget = 0x0400;
  else if (!strncmp(p, "$z", 2)) p += 2, prm->ztarget = 0x0400;
  else prm->pb = strtod(p, &p);
  while (*p == ' ' || *p == '\t') ++p;
  if (!strncmp(p, "$x", 2)) p += 2, prm->xtarget = 0x0500;
  else if (!strncmp(p, "$z", 2)) p += 2, prm->ztarget = 0x0500;
  else prm->sb = strtod(p, &p);
  while (*p == ' ' || *p == '\t') ++p;
  prm->avg = strtoul(p, &p, 10);
  while (*p == ' ' || *p == '\t') ++p;
  prm->xstart = strtod(p, &p);
  while (*p == ' ' || *p == '\t') ++p;
  prm->xend = strtod(p, &p);
  while (*p == ' ' || *p == '\t') ++p;
  prm->zstart = strtod(p, &p);
  while (*p == ' ' || *p == '\t') ++p;
  prm->zend = strtod(p, &p);
  while (*p == ' ' || *p == '\t') ++p;
  prm->steps = strtoul(p, &p, 10);
  while (*p == ' ' || *p == '\t') ++p;
  if (!strncmp(p, "$x", 2)) p += 2, prm->xtarget = 0x0600;
  else if (!strncmp(p, "$z", 2)) p += 2, prm->ztarget = 0x0600;
  else prm->labels = strtoul(p, &p, 10);
  
  delete [] buffer;
  
  return true;
}

MarkovChain *DTMC(int, int, int, double, double, int, double, double, double);
void Write_DTMC(MarkovChain*, FILE*);

// Perform the benchmark on models generated according to rm using the given simulator.
// flagvector is an array of bitsets specifying which optimizations to use; see Strong.h.
// flagvectorsize is the number of items in flagvector. averages is the number
// of times the simulation is to be performed and averaged over.
void benchmark_random_model(RandomModel *rm, unsigned long *flagvector,
                            unsigned int flagvectorsize, unsigned int averages, 
                            double fpprecision, double *transitions, bool *tabulate,
                            double **result, bool rmap_extra, double *result_rmap, int pitch, int res0)
{
  unsigned int n, m;
  char command[384], temp[128];
  FILE *f;
  int transsum = 0;
  SimulationStatistics stats;
  double utime, stime, rtime;
  MarkovChain *mkc;

  for (n = 0; n < flagvectorsize; ++n)
  {
    for (m = 0; m < 11; ++m)
    {
      if (tabulate[m]) result[m][(n * pitch) + res0] = 0.0;
      if (tabulate[3] && rmap_extra && result_rmap) result_rmap[res0] = 0.0;
    }
  }
    
  for (n = 0; n < rm->avg; ++n)
  {
    strcpy(&temp[0], "rdbenchXXXXXX");
    f = fdopen(mkstemp(&temp[0]), "wb");
    mkc = DTMC(rm->n, rm->a, rm->b, rm->fb, rm->lb, rm->c, rm->cb, rm->pb, rm->sb);
    transsum += mkc->Transitions();
    mkc->Write(f);
    fclose(f);
    delete mkc;

    for (m = 0; m < flagvectorsize; ++m)
    {
      sprintf(&command[0], "%s/benchone_%lu strong dtmc \"!%u\" %e %u %s", &myself[0], flagvector[m], rm->labels, fpprecision, averages, &temp[0]);

      f = popen(&command[0], "r");
      if (fscanf(f, "%le %le %le %u %u %u %u %u %u %u %lu %lu %lu %lu %lu %u %u %u %u\n",
        &utime, &stime, &rtime, &stats.num_partitions,
        &stats.num_iterations, &stats.num_initial_pairs, &stats.num_final_pairs,
        &stats.num_maxflow,
        &stats.num_p_invariant_fails, &stats.num_sig_arc_fails, &stats.mem_relation_map,
        &stats.mem_partition_map, &stats.mem_relation,
        &stats.mem_maxflow, &stats.mem_model, &stats.min_complexity, &stats.max_complexity,
        &stats.num_nets_cached, &stats.num_cache_hits) < 19)
      {
        utime = -1, rtime = -1, stime = -1;
        memset(&stats, 0, sizeof(stats));
        fprintf(stderr, "Warning: Benchmark [%d %d %d %d %.3f %.3f %.3f %.3f %.3f], cfg %02lx failed\n", rm->n, rm->a, rm->b, rm->c, rm->fb, rm->lb, rm->cb, rm->pb, rm->sb, flagvector[m]);
      }
      pclose(f);

      if (tabulate[0]) result[0][(m * pitch) + res0] = utime;
      if (tabulate[1]) result[1][(m * pitch) + res0] = stime;
      if (tabulate[2]) result[2][(m * pitch) + res0] = rtime;
      if (tabulate[3])
      {
        result[3][(m * pitch) + res0] += double(stats.mem_relation + stats.mem_maxflow
                                                 + stats.mem_partition_map + stats.mem_model);
        if (!rmap_extra) result[3][(m * pitch) + res0] += double(stats.mem_relation_map);
        else result_rmap[res0] += double(stats.mem_relation_map);
      }
      if (tabulate[4]) result[4][(m * pitch) + res0] += double(stats.num_initial_pairs);
      if (tabulate[5]) result[5][(m * pitch) + res0] += double(stats.num_final_pairs);
      if (tabulate[6]) result[6][(m * pitch) + res0] += double(stats.num_partitions);
      if (tabulate[7]) result[7][(m * pitch) + res0] += double(stats.num_iterations);
      if (tabulate[8]) result[8][(m * pitch) + res0] += double(stats.num_maxflow);
      if (tabulate[9]) result[9][(m * pitch) + res0] += double(stats.num_p_invariant_fails);
      if (tabulate[10]) result[10][(m * pitch) + res0] += double(stats.num_sig_arc_fails);
      if (tabulate[11]) result[11][(m * pitch) + res0] += double(stats.min_complexity);
      if (tabulate[12]) result[12][(m * pitch) + res0] += double(stats.max_complexity);
      if (tabulate[13]) result[13][(m * pitch) + res0] += double(stats.num_nets_cached);
      if (tabulate[14]) result[14][(m * pitch) + res0] += double(stats.num_cache_hits);
    }

    unlink(&temp[0]);
  }

  if (rm->avg > 1)
  {
    for (n = 0; n < flagvectorsize; ++n)
    {
      for (m = 0; m < 11; ++m)
      {
        if (tabulate[m]) result[m][(n * pitch) + res0] /= rm->avg;
      }
    }

    if (tabulate[3] && rmap_extra && result_rmap) result_rmap[res0] /= rm->avg;
    
    *transitions = double(transsum) / rm->avg;
  }
}

int minprecision(double d)
{
  int p = 0;
  d = fabs(d);
  while (d < 1)
  {
    ++p;
    d *= 10;
  }
  return p;
}

// Print usage of the program
void usage()
{
  const char *usage = "benchmark [options] model ... \n"
"Options:\n"
"  --type, -t <type> Model type; one out of ctmc, dtmc, pa, cpa, random.\n"
"                    'random' models are scripts containing instructions\n"
"                    for generating models with certain parameters.\n"
"  --avg, -a <N>     Average times over N simulations (default 1)\n"
"  --name, -o <name> Name of benchmark (default 'benchmark')\n"
"  --labels, -n <N>  [This option is deprecated and is ignored with a warning]\n"
"  --lext, -x <ext>  Suffix which is appended to model file names and used\n"
"                    as a label specification for that model. Has no effect\n"
"                    on random models (default '.label')\n"
"  --data, -d <data> Data to tabulate. One of usertime, systemtime,\n"
"                    realtime, memory, initialsize, finalsize, partitions,\n"
"                    iterations, maxflow, pivfail, safail, all. Multiple\n"
"                    specifications are possible. (default usertime)\n"
"  --extra-rmap      For memory, show memory used by relation map in a\n"
"                    separate line (LaTeX and plain output only)\n"
"  --model-info      Add state and transition numbers in separate rows\n"
"                    (LaTeX and plain output only)\n"
"  --plot, -p <src>  Generate <name>_<data>.data, <name>_<data>.gnuplot\n"
"                    and <name>_<data>.eps. 'src' is the data source and\n"
"                    may be one of n, m, r, r+ (see below)\n"
"  --latex, -l       Generate a LaTeX table in <name>_<data>.tex\n"
"  --quiet, -q       Suppress plain-text output to stdout\n"
"  --dumprel         Dump computed simulation relation to <name>_rel.txt\n"
"  --time-unit <U>   Unit for time specs (ms, s, m, h)\n"
"  --space-unit <U>  Unit for memory (B, kB, MB, GB)\n"
"  --precision <N>   Number of decimals (default 3)\n"
"  --fp-approx <e>   Values with a difference less than e are equal\n"
"  --{x|y}rng <R>    Explicitly specify ranges in format start:end, determined\n"
"                    from data by default (gnuplot output only). Has no effect\n"
"                    in 'r' and 'r+' plotting modes.\n"
"  --opt, -O <opt>   Specify an optimization to benchmark in the format\n"
"                    opt[:title] where opt is an integer, optionally\n"
"                    followed by a colon-separated alias for the configuration.\n"
"                    See sources for details. Specify multiple times to\n"
"                    compare different configurations or \"all\" to\n"
"                    benchmark all possible configurations.\n"
"  --                Stop options scanning\n"
"\n"
"Units: By default, the average of all values is taken and the unit is\n"
"       adjusted to put the average in an adequate range.\n"
"\n"
"Data type: usertime    - processing time in user mode\n"
"           systemtime  - processing time in system mode\n"
"           realtime    - real time that elapsed\n"
"           memory      - peak amount of memory consumed at one point in time\n"
"           initialsize - # of pairs/blocks in the initial relation\n"
"           finalsize   - # of pairs in result set\n"
"           partitions  - # of state partitions found\n"
"           iterations  - # of iterations used\n"
"           maxflow     - # of times maxflow was used and not aborted\n"
"                         early due to a trivial network\n"
"           pivfail     - # of times maxflow was aborted due to P-Invariant\n"
"           safail      - # of times a significiant arc was deleted\n"
"           cache       - number of networks cached by parametric maxflow\n"
"           cachehits   - number of times a saved network was reused\n"
"           all         - tabulate all of the above\n"
"\n"
"Data source: n  - states\n"
"             m  - transitions\n"
"             r  - random models: use $x variable\n"
"             r+ - random models: use $x and $z, creates 3D plot\n"
"\n";

  fputs(usage, stderr);
  exit(0);
}

void progress_callback(unsigned int done, unsigned int total)
{
  static unsigned int progress = 0;
  unsigned int current;
  
  if (done == 0 && total == 0)
  {
    progress = 0;
    return;
  }
  
  current = (unsigned int)floor(20.0 / double(total) * double(done));
  current -= progress;
  progress += current;
  
  while (current--) fprintf(stderr, ".");
}

// Find the greatest number in a set and return the number of digits
// prior to the decimal point.
int find_max_digits(double *data, unsigned int count)
{
  int max = 0, x;
  unsigned int n;
  
  for (n = 0; n < count; ++n) if ((x = (data[n] == 0.0 ? 1 : (int)floor(::log10(fabs(data[n]))))) > max) max = x;
  
  return max + 1;
}

// Return the average over a set.
double find_average(double *data, unsigned int count)
{
  double d = 0.0;
  unsigned int n;
  
  for (n = 0; n < count; ++n) d += data[n];
  
  return d / count;
}

// Find minimum/maximum of a set of values.
double find_data_range(double *data, unsigned int count, bool find_min)
{
  double min, max = 0.0;
  unsigned int n;
  
  for (n = 0; n < count; ++n) if (data[n] > max) max = data[n];
  
  min = max;
  
  for (n = 0; n < count; ++n) if (data[n] < min) min = data[n];
  
  return (find_min ? min : max);
}

// Determine the best unit to display a time value in (assuming the value is given in seconds)
unsigned int find_best_time_unit(double avg)
{
  if (avg < 2) return 1;
  else if (avg < 100) return 2;
  else if (avg < 3600) return 3;
  else return 4;
}

// Determine the best unit to display an amount of memory in (assuming the value is given in bytes)
unsigned int find_best_space_unit(double avg)
{
  if (avg < 1024) return 1;
  else if (avg < 1048576) return 2;
  else if (avg < 1048576 * 1024) return 3;
  else return 4;
}

// Transform the unit of all values in a set assuming the values are currently in the
// base unit (seconds/bytes).
void unit_transform(double *data, unsigned int count, unsigned int to, bool space)
{
  double multipliers[2][4] = {
    {1000.0, 1.0, 1.0 / 60.0, 1.0 / 3600.0}, {1.0, 1.0 / 1024.0, 1.0 / 1048576.0, 1.0 / (1048576.0 * 1024.0)}
  }, multiplier;
  
  multiplier = multipliers[space ? 1 : 0][to - 1];
  
  for (unsigned int n = 0; n < count; ++n) data[n] *= multiplier;
}

// Return a pointer to the filename portion of a path without
// allocating a separate buffer
const char *basenameptr(const char *path)
{
  const char *slash = strrchr(path, '/');
  if (slash) return slash + 1;
  slash = strrchr(path, '\\');
  if (slash) return slash + 1;
  return path;
}

// Compute the width of the plain-text table for a particular result
int find_table_width(double *data, unsigned int rows, unsigned int cols, unsigned long *flagvector,
                     std::vector<const char*> &cfg_titles, std::vector<const char*> &files, double *extra,
                     const char *extra_title, int precision, std::vector<unsigned int> &m_n,
                     std::vector<double> &m_m, std::vector<unsigned int> &m_a)
{
  unsigned int n, m;
  int num_width, str_width, len;
  std::vector<const char*>::iterator i;
  
  // Determine column formatting width for model columns
  n = find_max_digits(data, rows * cols);
  if (extra && extra_title)
  {
    m = find_max_digits(extra, cols);
    if (m > n) num_width = m;
    else num_width = n;
  }
  else num_width = n;
  if (precision > 0) num_width += precision + 1;
  
  for (n = 0; n < cols; ++n)
  {
    if ((int)strlen(basenameptr(files[n])) > num_width) num_width = strlen(basenameptr(files[n]));
    if (m_n.size() == files.size() && (int)floor(::log10(m_n[n])) + 1 > num_width) num_width = (int)floor(::log10(m_n[n])) + 1;
    if (m_m.size() == files.size() && (int)floor(::log10(m_m[n])) + 1 > num_width) num_width = (int)floor(::log10(m_m[n])) + 1;
    if (m_a.size() == files.size() && (int)floor(::log10(m_a[n])) + 1 > num_width) num_width = (int)floor(::log10(m_a[n])) + 1;
  }
  
  // Determine column formatting width for first column
  if (extra && extra_title) str_width = strlen(extra_title);
  else str_width = 0;
  for (n = 0, i = cfg_titles.begin(); n < rows; ++n, ++i)
  {
    len = (*i ? strlen(*i) : (flagvector[n] == 0 ? 0 : (int)floor(log((double)flagvector[n]) / log(16.0))) + 3);
    if (len > str_width) str_width = len;
  }
  
  if (m_m.size() == files.size() && str_width < 6) str_width = 6;
  if (m_a.size() == files.size() && str_width < 7) str_width = 7;
  
  return str_width + (cols * (num_width + 1));
}

// Print a plain-text table of the benchmark results for one particular data type
void generate_plain(FILE *f, double *data, unsigned int rows, unsigned int cols, unsigned long *flagvector,
                   std::vector<const char*> &cfg_titles, std::vector<const char*> &files, double *extra,
                   const char *extra_title, int precision, std::vector<unsigned int> &m_n,
                   std::vector<double> &m_m, std::vector<unsigned int> &m_a)
{
  int num_width, str_width, len;
  unsigned int n, m;
  std::vector<const char*>::iterator i;
  
  // Determine column formatting width for model columns
  n = find_max_digits(data, rows * cols);
  if (extra && extra_title)
  {
    m = find_max_digits(extra, cols);
    if (m > n) num_width = m;
    else num_width = n;
  }
  else num_width = n;
  if (precision > 0) num_width += precision + 1;
  
  for (n = 0; n < cols; ++n)
  {
    if ((int)strlen(basenameptr(files[n])) > num_width) num_width = strlen(basenameptr(files[n]));
    if (m_n.size() == files.size() && (int)floor(::log10(m_n[n])) + 1 > num_width) num_width = (int)floor(::log10(m_n[n])) + 1;
    if (m_m.size() == files.size() && (int)floor(::log10(m_m[n])) + 1 > num_width) num_width = (int)floor(::log10(m_m[n])) + 1;
    if (m_a.size() == files.size() && (int)floor(::log10(m_a[n])) + 1 > num_width) num_width = (int)floor(::log10(m_a[n])) + 1;
  }
  
  // Determine column formatting width for first column
  if (extra && extra_title) str_width = strlen(extra_title);
  else str_width = 0;
  for (n = 0, i = cfg_titles.begin(); n < rows; ++n, ++i)
  {
    len = (*i ? strlen(*i) : (flagvector[n] == 0 ? 0 : (int)floor(log((double)flagvector[n]) / log(16.0))) + 3);
    if (len > str_width) str_width = len;
  }
  
  if (m_m.size() == files.size() && str_width < 6) str_width = 6;
  if (m_a.size() == files.size() && str_width < 7) str_width = 7;
  
  // Print column headers
  fprintf(f, "%*s ", str_width, "");
  for (i = files.begin(); i != files.end(); ++i) fprintf(f, "%-*s ", num_width, basenameptr(*i));
  fputs("\n", f);
  
  // Print model info if provided
  if (m_n.size() == files.size())
  {
    fprintf(f, "%-*s ", str_width, "States");
    for (n = 0; n < cols; ++n) fprintf(f, "%*u ", num_width, m_n[n]);
    fputs("\n", f);
  }
  if (m_m.size() == files.size())
  {
    fprintf(f, "%-*s ", str_width, "Trans.");
    for (n = 0; n < cols; ++n) fprintf(f, "%*u ", num_width, (unsigned int)round(m_m[n]));
    fputs("\n", f);
  }
  if (m_a.size() == files.size())
  {
    fprintf(f, "%-*s ", str_width, "Actions");
    for (n = 0; n < cols; ++n) fprintf(f, "%*u ", num_width, m_a[n]);
    fputs("\n", f);
  }
  if (m_n.size() == files.size() || m_m.size() == files.size() || m_a.size() == files.size()) fputs("\n", f);
  
  // Print extra data row if provided
  if (extra && extra_title)
  {
    fprintf(f, "%-*s ", str_width, extra_title);
    for (n = 0; n < cols; ++n) fprintf(f, "%*.*f ", num_width, precision, extra[n]);
    fputs("\n", f);
  }
  
  // Print table body
  for (n = 0, i = cfg_titles.begin(); n < rows; ++n, ++i)
  {
    if (*i) fprintf(f, "%-*s ", str_width, *i);
    else fprintf(f, "0x%-*lX ", str_width - 2, flagvector[n]);
    
    for (m = 0; m < cols; ++m)
    {
      fprintf(f, "%*.*f ", num_width, precision, data[n + (m * rows)]);
    }
    
    fputs("\n", f);
  }
}

const char *get_axis_name(unsigned int target, bool shrt)
{
  bool lc = ((target & 0x8000) != 0);
  target &= 0x7fff;
  
  switch (target)
  {
  case 0x0001: return (shrt ? (lc ? "states" : "States") : "Number of States");
  case 0x0002: return (shrt ? (lc ? "minsuc" : "Min.Suc.") : "Minimum Number of Successors");
  case 0x0003: return (shrt ? (lc ? "maxsuc" : "Max.Suc.") : "Maximum Number of Successors");
  case 0x0004: return (shrt ? (lc ? "clster" : "Clusters") : "Number of Clusters");
  case 0x0100: return (shrt ? (lc ? "fbias" : "F-Bias") : "Fanout Bias");
  case 0x0200: return (shrt ? (lc ? "lbias" : "L-Bias") : "Linearity Bias");
  case 0x0300: return (shrt ? (lc ? "cbias" : "C-Bias") : "Clustering Bias");
  case 0x0400: return (shrt ? (lc ? "pbias" : "P-Bias") : "Probability Distribution Bias");
  case 0x0500: return (shrt ? (lc ? "sbias" : "S-Bias") : "Successor Bias");
  case 0x0600: return (shrt ? (lc ? "labels" : "Labels") : "Number of Labels");
  }
  
  return "";
}

// Print a plain-text table of the benchmark results for one particular data type for random models
void gen_rnd_plain(FILE *f, RandomModel *prm, double *data, unsigned int rows, unsigned int cols,
                   unsigned long *flagvector, std::vector<const char*> &cfg_titles, double *extra,
                   const char *extra_title, int precision, std::vector<unsigned int> &m_n,
                   std::vector<double> &m_m, std::vector<unsigned int>&, const char *mtitle,
                   const char *datatype, bool dataonly = false)
{
  int num_width, len, prstepx, prstepz, xsteps, zsteps, totalwidth;
  int n, m, l, k;
  std::vector<const char*>::iterator i;
  std::vector<int> cwidth;
  const char *xtitle, *ztitle;
  
  xtitle = get_axis_name(prm->xtarget, true);
  ztitle = get_axis_name(prm->ztarget, true);
  
  // Determine width of values in table
  num_width = find_max_digits(data, rows * cols);
  if (precision > 0) num_width += precision + 1;
  
  // Determine $x column width (if any)
  if (prm->xtarget != 0xffff)
  {
    m = (fabs(prm->xend) > fabs(prm->xstart) ? find_max_digits(&prm->xend, 1) : find_max_digits(&prm->xstart, 1));
    if (prm->xtarget & 0xff) prstepx = 0;
    else prstepx = precision;
    m += (prstepx > 0 ? prstepx + 1 : 0);
    if (prm->xstart < 0.0 || prm->xend < 0.0) ++m;
    if ((int)strlen(xtitle) > m) cwidth.push_back(strlen(xtitle));
    else cwidth.push_back(m);
  }
  
  // Determine $z column width (if any)
  if (prm->ztarget != 0xffff)
  {
    m = (fabs(prm->zend) > fabs(prm->zstart) ? find_max_digits(&prm->zend, 1) : find_max_digits(&prm->zstart, 1));
    if (prm->ztarget & 0xff) prstepz = 0;
    else prstepz = precision;
    m += (prstepz > 0 ? prstepz + 1 : 0);
    if (prm->zstart < 0.0 || prm->zend < 0.0) ++m;
    if ((int)strlen(ztitle) > m) cwidth.push_back(strlen(ztitle));
    else cwidth.push_back(m);
  }
  
  // Determine model info column widths
  if (m_n.size() == rows)
  {
    len = 6;
    for (n = 0; n < (int)rows; ++n)
    {
      m = (int)floor(::log10(m_n[n])) + 1;
      if (m > len) len = m;
    }
    cwidth.push_back(len);
  }
  if (m_m.size() == rows)
  {
    len = 6;
    for (n = 0; n < (int)rows; ++n)
    {
      m = (int)floor(::log10(m_m[n])) + 1;
      if (m > len) len = m;
    }
    cwidth.push_back(len);
  }
  
  // Determine extra data column width
  if (extra && extra_title)
  {
    len = find_max_digits(extra, rows) + precision + 1;
    if ((int)strlen(extra_title) > len) cwidth.push_back(strlen(extra_title));
    else cwidth.push_back(len);
  }
  
  // Determine width for remaining data columns
  for (i = cfg_titles.begin(), n = 0; n < (int)cols; ++n, ++i)
  {
    len = (*i ? strlen(*i) : (flagvector[n] == 0 ? 0 : (int)floor(log((double)flagvector[n]) / log(16.0))) + 3);
    if (len > num_width) cwidth.push_back(len);
    else cwidth.push_back(num_width);
  }
  
  // Determine width of entire table
  totalwidth = cwidth.size() - 1;
  for (n = 0; n < (int)cwidth.size(); ++n) totalwidth += cwidth[n];
  
  // Print table header if supplied
  if (!dataonly)
  {
    for (n = 0; n < totalwidth; ++n) fputs("=", f);
    if (mtitle)
    {
      len = strlen(mtitle) + (datatype ? strlen(datatype) + 2 : 0);
      if (len > totalwidth)
      {
        for (n = 0; n < len - totalwidth; ++n) fputs("=", f);
      }
      fputs("\n", f);
      if (datatype) fprintf(f, "%s: %s\n", mtitle, datatype);
      else fprintf(f, "%s\n", mtitle);
      for (n = 0; n < totalwidth; ++n) fputs("-", f);
    }
    fputs("\n", f);
  
    // Print column headers
    n = 0;
    if (prm->xtarget != 0xffff) fprintf(f, "%-*s ", cwidth[n++], xtitle);
    if (prm->ztarget != 0xffff) fprintf(f, "%-*s ", cwidth[n++], ztitle);
    if (m_n.size() == rows) fprintf(f, "%-*s ", cwidth[n++], "States");
    if (m_m.size() == rows) fprintf(f, "%-*s ", cwidth[n++], "Trans.");
    if (extra && extra_title) fprintf(f, "%-*s ", cwidth[n++], extra_title);
    for (i = cfg_titles.begin(), m = 0; i != cfg_titles.end(); ++i, ++m)
    {
      if (*i) fprintf(f, "%-*s ", cwidth[n++], *i);
      else fprintf(f, "%#-*lx ", cwidth[n++], flagvector[m]);
    }
    fputs("\n", f);
  }
  
  // Print table body
  xsteps = (prm->xtarget == 0xffff ? 1 : prm->steps);
  zsteps = (prm->ztarget == 0xffff ? 1 : prm->steps);
  for (m = 0; m < xsteps; ++m)
  {
    for (k = 0; k < zsteps; ++k)
    {
      n = 0;
      if (prm->xtarget != 0xffff)
      {
        if (prm->ztarget == 0xffff || k == 0) fprintf(f, "%*.*f ", cwidth[n++], prstepx, prm->xstart + ((prm->xend - prm->xstart) * m / (prm->steps - 1)));
        else fprintf(f, "%*s ", cwidth[n++], "");
      }
      if (prm->ztarget != 0xffff) fprintf(f, "%*.*f ", cwidth[n++], prstepz, prm->zstart + ((prm->zend - prm->zstart) * k / (prm->steps - 1)));
      if (m_n.size() == rows) fprintf(f, "%*u ", cwidth[n++], m_n[(m * zsteps) + k]);
      if (m_m.size() == rows) fprintf(f, "%*.0f ", cwidth[n++], m_m[(m * zsteps) + k]);
      if (extra && extra_title) fprintf(f, "%*.*f ", cwidth[n++], precision, extra[(m * zsteps) + k]);
      for (l = 0; l < (int)cols; ++l)
      {
        fprintf(f, "%*.*f ", cwidth[n++], precision, data[(l * rows) + (m * zsteps) + k]);
      }
      fputs("\n", f);
    }
  }
  
  if (!dataonly)
  {
    for (n = 0; n < totalwidth; ++n) fputs("=", f);
    fputs("\n\n", f);
  }
}

// Generate a table that can be inserted directly into a LaTeX document
void generate_latex(FILE *f, double *data, unsigned int rows, unsigned int cols, unsigned long *flagvector,
                    std::vector<const char*> &cfg_titles, std::vector<const char*> &files, double *extra,
                    const char *extra_title, int precision, const char *table_title, const char *unit,
                    std::vector<unsigned int> &m_n, std::vector<double> &m_m, std::vector<unsigned int> &m_a)
{
  int num_width, str_width, len;
  unsigned int n, m;
  std::vector<const char*>::iterator i;
  
  // Determine column formatting width for model columns
  n = find_max_digits(data, rows * cols);
  if (extra && extra_title)
  {
    m = find_max_digits(extra, cols);
    if (m > n) num_width = m;
    else num_width = n;
  }
  else num_width = n;
  if (precision > 0) num_width += precision + 1;
  
  for (n = 0; n < cols; ++n)
  {
    if ((int)strlen(basenameptr(files[n])) > num_width) num_width = strlen(basenameptr(files[n]));
    if (m_n.size() == files.size() && (int)floor(::log10(m_n[n])) + 1 > num_width) num_width = (int)floor(::log10(m_n[n])) + 1;
    if (m_m.size() == files.size() && (int)floor(::log10(m_m[n])) + 1 > num_width) num_width = (int)floor(::log10(m_m[n])) + 1;
    if (m_a.size() == files.size() && (int)floor(::log10(m_a[n])) + 1 > num_width) num_width = (int)floor(::log10(m_a[n])) + 1;
  }
  
  // Determine column formatting width for first column
  if (extra && extra_title) str_width = strlen(extra_title);
  else str_width = 0;
  for (n = 0, i = cfg_titles.begin(); n < rows; ++n, ++i)
  {
    len = (*i ? strlen(*i) : (flagvector[n] == 0 ? 0 : (int)floor(log((double)flagvector[n]) / log(16.0))) + 3);
    if (len > str_width) str_width = len;
  }
  
  if (m_m.size() == files.size() && str_width < 6) str_width = 6;
  if (m_a.size() == files.size() && str_width < 7) str_width = 7;
  
  // Print LaTeX table head
  fputs("\\begin{table*}\n", f);
  fputs("\\centering\n", f);
  fputs("\\begin{tabular}{|", f);
  for (n = 0; n < cols; ++n) fputs("r|", f);
  fputs("}\n", f);
  fputs("\\hline\n", f);
  
  // Print column headers
  fprintf(f, "%*s ", str_width, "");
  for (n = 0, i = files.begin(); n < cols; ++n, ++i) fprintf(f, "& %-*s ", num_width, basenameptr(*i));
  fputs("\\\\\n", f);
  fputs("\\hline\n", f);
  
  // Print model info if provided
  if (m_n.size() == files.size())
  {
    fprintf(f, "%-*s ", str_width, "States");
    for (n = 0; n < cols; ++n) fprintf(f, "& %*u ", num_width, m_n[n]);
    fputs("\\\\\n", f);
  }
  if (m_m.size() == files.size())
  {
    fprintf(f, "%-*s ", str_width, "Trans.");
    for (n = 0; n < cols; ++n) fprintf(f, "& %*u ", num_width, (unsigned int)round(m_m[n]));
    fputs("\\\\\n", f);
  }
  if (m_a.size() == files.size())
  {
    fprintf(f, "%-*s ", str_width, "Actions");
    for (n = 0; n < cols; ++n) fprintf(f, "& %*u ", num_width, m_a[n]);
    fputs("\\\\\n", f);
  }
  if (m_n.size() == files.size() || m_m.size() == files.size() || m_a.size() == files.size()) fputs("\\hline\n", f);
  
  // Print extra data row if provided
  if (extra && extra_title)
  {
    fprintf(f, "%-*s ", str_width, extra_title);
    for (n = 0; n < cols; ++n) fprintf(f, "& %*.*f ", num_width, precision, extra[n]);
    fputs("\\\\\n\\hline\n", f);
  }
  
  // Print table body
  for (n = 0, i = cfg_titles.begin(); n < rows; ++n, ++i)
  {
    if (*i) fprintf(f, "%-*s ", str_width, *i);
    else fprintf(f, "0x%-*lX ", str_width - 2, flagvector[n]);
    for (m = 0; m < cols; ++m)
    {
      fprintf(f, "& %*.*f ", num_width, precision, data[n + (m * rows)]);
    }
    fputs("\\\\\n", f);
  }
  
  // Print LaTeX table footer
  fputs("\\hline\n", f);
  fputs("\\end{tabular}\n", f);
  fprintf(f, "\\caption{StrongSimulation benchmark. Values are in %s\\label{tab:%s}}\n", unit, table_title);
  fputs("\\end{table*}\n", f);
}

// Generate the data file for a plot; requires a list of values to be used to represent the models
// on the X axis (plot_key), typically number of states or transitions or similar
void generate_gnuplot_data(FILE *f, double *data, double *plot_key, unsigned int rows, unsigned int cols,
                           unsigned int keys, unsigned int grid = 0, bool rotate = false)
{
  unsigned int n, m;
  
  if (rotate)
  {
    for (n = 0; n < cols; ++n)
    {
      for (m = 0; m < keys; ++m)
      {
        fprintf(f, "%.10f ", plot_key[(m * cols) + n]);
      }
    
      for (m = 0; m < rows; ++m)
      {
        fprintf(f, " %.10f", data[n + (m * cols)]);
      }
    
      if (grid > 1)
      {
        if (n % grid == grid - 1) fprintf(f, "\n");
      }
      fprintf(f, "\n");
    }
  }
  else
  {
    for (n = 0; n < cols; ++n)
    {
      for (m = 0; m < keys; ++m)
      {
        fprintf(f, "%.10f ", plot_key[(m * cols) + n]);
      }
    
      for (m = 0; m < rows; ++m)
      {
        fprintf(f, " %.10f", data[m + (n * rows)]);
      }
    
      if (grid > 1)
      {
        if (n % grid == grid - 1) fprintf(f, "\n");
      }
      fprintf(f, "\n");
    }
  }
}

// Generate a gnuplot macro file for a plot
void generate_gnuplot_macro(FILE *f, const char *xlabel, const char *ylabel, const char *unit,
                            double xlow, double xhigh, double ylow, double yhigh,
                            const char *title, const char *dataname, const char *filename, unsigned long *flagvector,
                            std::vector<const char*> &cfg_titles, unsigned long flagvectorsize)
{
  unsigned int n;
  std::vector<const char*>::iterator i;
  
  fprintf(f, "set xlabel \"%s\" font \"Helvetica,24\"\n", xlabel);
  fprintf(f, "set ylabel \"%s (%s)\" font \"Helvetica,24\"\n", ylabel, unit);
  fprintf(f, "set xrange [%f:%f]\n", xlow, xhigh);
  fprintf(f, "set yrange [%f:%f]\n", ylow, yhigh);
  fprintf(f, "set data style lines\n");
  fprintf(f, "set term postscript enhanced color\n");
  fprintf(f, "set title \"%s\" font \"Helvetica,24\"\n", title);
  fprintf(f, "set output \"%s\"\n", filename);
  
  for (n = 0, i = cfg_titles.begin(); n < flagvectorsize; ++n, ++i)
  {
    if (*i) fprintf(f, "%s \"%s\" using 1:%u t \"%s\"", (n == 0 ? "plot" : "    "), dataname, n + 2, *i);
    else    fprintf(f, "%s \"%s\" using 1:%u t \"0x%lX\"", (n == 0 ? "plot" : "    "), dataname, n + 2, flagvector[n]);
    if (n < flagvectorsize - 1) fputs(", \\\n", f);
    else fputs("\n", f);
  }
}

// Begin a gnuplot macro file for a plot of random models
void rnd_2d_gnuplot_macro(FILE *f, const char *xlabel, const char *ylabel, const char *unit,
                          const char *title, const char *filename, const char *dataname, unsigned long *flagvector,
                          unsigned long flagvectorsize, std::vector<const char*> cfg_titles,
                          std::set<RandomModelPlotInfo*> modelmap)
{
  unsigned int n;
  std::vector<const char*>::iterator i;
  std::set<RandomModelPlotInfo*>::iterator j;
  char datasource[160];
  
  fprintf(f, "set xlabel \"%s\" font \"Helvetica,24\"\n", xlabel);
  fprintf(f, "set ylabel \"%s (%s)\" font \"Helvetica,24\"\n", ylabel, unit);
  fprintf(f, "set data style lines\n");
  fprintf(f, "set term postscript enhanced color\n");
  fprintf(f, "set title \"%s\" font \"Helvetica,24\"\n", title);
  fprintf(f, "set output \"%s\"\n", filename);
  
  for (j = modelmap.begin(); j != modelmap.end(); ++j)
  {
    sprintf(&datasource[0], &(*j)->datasource[0], dataname);
    for (n = 0, i = cfg_titles.begin(); n < flagvectorsize; ++n, ++i)
    {
      if (*i) fprintf(f, "%s \"%s\" using 1:%u t \"%s\"", (n == 0 ? "plot" : "    "), &datasource[0], n + 2, *i);
      else    fprintf(f, "%s \"%s\" using 1:%u t \"0x%lX\"", (n == 0 ? "plot" : "    "), &datasource[0], n + 2, flagvector[n]);
      if (n < flagvectorsize - 1) fputs(", \\\n", f);
      else fputs("\n", f);
    }
  }
}

// Begin a gnuplot macro file for a plot of random models
void rnd_3d_gnuplot_macro(FILE *f, const char *xlabel, const char *ylabel, const char *zlabel, const char *unit,
                          const char *title, const char *filename, const char *dataname, unsigned long *flagvector,
                          unsigned long flagvectorsize, std::vector<const char*> cfg_titles,
                          std::set<RandomModelPlotInfo*> modelmap)
{
  unsigned int n;
  std::vector<const char*>::iterator i;
  std::set<RandomModelPlotInfo*>::iterator j;
  char datasource[160];
  
  fprintf(f, "set xlabel \"%s\" font \"Helvetica,24\"\n", xlabel);
  fprintf(f, "set ylabel \"%s\" font \"Helvetica,24\"\n", zlabel);
  fprintf(f, "set zlabel \"%s (%s)\" font \"Helvetica,24\"\n", ylabel, unit);
  fprintf(f, "set data style lines\n");
  fprintf(f, "set pm3d\n");
  fprintf(f, "set term postscript enhanced color\n");
  fprintf(f, "set title \"%s\" font \"Helvetica,24\"\n", title);
  fprintf(f, "set output \"%s\"\n", filename);
  
  for (j = modelmap.begin(); j != modelmap.end(); ++j)
  {
    sprintf(&datasource[0], &(*j)->datasource[0], dataname);
    for (n = 0, i = cfg_titles.begin(); n < flagvectorsize; ++n, ++i)
    {
      if (*i) fprintf(f, "%s \"%s\" using 1:2:%u t \"%s\"", (n == 0 ? "splot" : "     "), &datasource[0], n + 3, *i);
      else    fprintf(f, "%s \"%s\" using 1:2:%u t \"0x%lX\"", (n == 0 ? "splot" : "     "), &datasource[0], n + 3, flagvector[n]);
      if (n < flagvectorsize - 1) fputs(", \\\n", f);
      else fputs("\n", f);
    }
  }
}

bool check_label_file(const char *model, const char *suffix, bool quiet)
{
  char *buf = new char[strlen(model) + strlen(suffix) + 1];
  FILE *f;
  bool success;
  
  strcpy(buf, model);
  strcat(buf, suffix);
  
  f = fopen(buf, "rb");
  success = (bool)f;
  if (success) fclose(f);
  
  if (!success && !quiet) fprintf(stderr, "Warning: Label specification '%s' for model '%s' can't be read.\n"
                                          "         All states will have the same label in this model.\n", buf, model);
  
  return success;
}

// Main program
int main(int argc, char *argv[])
{
  // Massive variable declaration section of doom!
  int n, m, j, k;
  char *p, modeltype[8], command[384];
  FILE *f;
  unsigned int averages = 1, rsteps, zsteps;
  const char *labelext = ".label";
  int precision = 3, total;
  const char *title = "benchmark";
  bool gen_plot_n = false, gen_plot_m = false, gen_latex = false, gen_plain = true, all_opt_cfgs = false, gen_r2d = false, gen_r3d = false;
  bool data_given = false, tabulate[15] = {false, false, false, false, false, false, false, false, false, false, false, false, false}, rmap_extra = false;
  double *result[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, *plot_key, *result_rmap = 0;
  double xlow, xhigh, ylow, yhigh, low, high, average, fpprecision = -1.0, transsum, utime, stime, rtime;
  bool xrange = false, yrange = false, model_info = false, random_model = false, dumprelation = false, dumpsuccessful = false;
  
  const char *t_units[] = {"", "ms", "s", "m", "h"};
  const char *s_units[] = {"", "Bytes", "kB", "MB", "GB"};
  unsigned int t_unit = 0, s_unit = 0;
  
  const char *datanames[15] = {"usertime", "systemtime", "realtime", "memory", "initialsize",
    "finalsize", "partitions", "iterations", "maxflow", "pivfail", "safail", "mincmplx", "maxcmplx",
    "cache", "cachehits"};
  char datasource[160], filename[160];
  
  ProbabilisticModel *pm = 0;
  RandomModel rm;
  
  unsigned long *flagvector = 0, flagvectorsize = 0;
  unsigned long allflagsvector[] = {0, 1, 2, 4, 5, 6, 16, 17, 18, 20, 21, 22, 24, 25, 26, 28, 29, 30};
  
  SimulationStatistics stats;
  
  std::vector<unsigned long> cfgs;
  std::vector<const char*> cfg_titles;
  std::vector<const char*> files;
  std::vector<const char*>::iterator i;
  std::vector<unsigned int> model_states, model_actions;
  std::vector<double> model_transitions;
  
  // If no parameters were given, print usage
  if (argc < 2) usage();
  
  strcpy(&command[0], argv[0]);
  strcpy(&myself[0], dirname(&command[0]));
  
  // Parse argument vector
  for (n = 1; n < argc; ++n)
  {
    if (!strcmp(argv[n], "--type") || !strcmp(argv[n], "-t"))
    {
      if (++n >= argc)
      {
        fprintf(stderr, "Error: Option %s requires a parameter\n", argv[n-1]);
        return 1;
      }
      if (!strcasecmp(argv[n], "dtmc")) pm = new MarkovChain;
      else if (!strcasecmp(argv[n], "ctmc")) pm = new MarkovChain;
      else if (!strcasecmp(argv[n], "pa")) pm = new ProbabilisticAutomaton;
      else if (!strcasecmp(argv[n], "cpa")) pm = new ProbabilisticAutomaton;
      else if (!strcasecmp(argv[n], "random"))
      {
        random_model = true;
        fprintf(stderr, "Info: Model type 'random' currently supports only DTMC type random models\n");
      }
      else
      {
        fprintf(stderr, "Error: Unknown model type: %s\n", argv[n]);
        return 1;
      }

      strcpy(&modeltype[0], argv[n]);
      p = &modeltype[0];
      while (*p)
      {
        if (*p >= 'A' && *p <= 'Z') *p |= 0x20;
        ++p;
      }
    }
    else if (!strcmp(argv[n], "--avg") || !strcmp(argv[n], "-a"))
    {
      if (++n >= argc)
      {
        fprintf(stderr, "Error: Option %s requires a parameter\n", argv[n-1]);
        return 1;
      }
      averages = strtoul(argv[n], 0, 0);
    }
    else if (!strcmp(argv[n], "--labels") || !strcmp(argv[n], "-n"))
    {
      if (++n >= argc)
      {
        fprintf(stderr, "Error: Option %s requires a parameter\n", argv[n-1]);
        return 1;
      }
      fprintf(stderr, "Warning: Option %s has been deprecated. Refer to the manual.\n", argv[n-1]);
    }
    else if (!strcmp(argv[n], "--lext") || !strcmp(argv[n], "-x"))
    {
      if (++n >= argc)
      {
        fprintf(stderr, "Error: Option %s requires a parameter\n", argv[n-1]);
        return 1;
      }
      labelext = argv[n];
    }
    else if (!strcmp(argv[n], "--name") || !strcmp(argv[n], "-o"))
    {
      if (++n >= argc)
      {
        fprintf(stderr, "Error: Option %s requires a parameter\n", argv[n-1]);
        return 1;
      }
      title = argv[n];
    }
    else if (!strcmp(argv[n], "--dumprel"))
    {
      dumprelation = true;
    }
    else if (!strcmp(argv[n], "--precision"))
    {
      if (++n >= argc)
      {
        fprintf(stderr, "Error: Option %s requires a parameter\n", argv[n-1]);
        return 1;
      }
      precision = strtol(argv[n], 0, 0);
    }
    else if (!strcmp(argv[n], "--fp-approx"))
    {
      if (++n >= argc)
      {
        fprintf(stderr, "Error: Option %s requires a parameter\n", argv[n-1]);
        return 1;
      }
      fpprecision = strtod(argv[n], 0);
      CompactMaxFlow<double>::precision = fpprecision;
    }
    else if (!strcmp(argv[n], "--data") || !strcmp(argv[n], "-d"))
    {
      if (++n >= argc)
      {
        fprintf(stderr, "Error: Option %s requires a parameter\n", argv[n-1]);
        return 1;
      }
      if (!strcasecmp(argv[n], "all")) for (m = 0; m < int(sizeof(tabulate) / sizeof(bool)); ++m) tabulate[m] = true;
      else
      {
        for (m = 0; m < int(sizeof(datanames) / sizeof(datanames[0])); ++m)
        {
          if (!strcasecmp(argv[n], datanames[m]))
          {
            tabulate[m] = true;
            data_given = true;
            break;
          }
        }
        if (m == sizeof(datanames) / sizeof(datanames[0])) fprintf(stderr, "Warning: Data type %s not known\n", argv[n]);
      }
    }
    else if (!strcmp(argv[n], "--time-unit"))
    {
      if (++n >= argc)
      {
        fprintf(stderr, "Error: Option %s requires a parameter\n", argv[n-1]);
        return 1;
      }
      if (!strcasecmp(argv[n], "ms")) t_unit = 1;
      else if (!strcasecmp(argv[n], "s")) t_unit = 2;
      else if (!strcasecmp(argv[n], "m")) t_unit = 3;
      else if (!strcasecmp(argv[n], "h")) t_unit = 4;
      else fprintf(stderr, "Warning: Time unit %s not known\n", argv[n]);
    }
    else if (!strcmp(argv[n], "--space-unit"))
    {
      if (++n >= argc)
      {
        fprintf(stderr, "Error: Option %s requires a parameter\n", argv[n-1]);
        return 1;
      }
      if (!strcasecmp(argv[n], "b")) s_unit = 1;
      else if (!strcasecmp(argv[n], "kb")) s_unit = 2;
      else if (!strcasecmp(argv[n], "mb")) s_unit = 3;
      else if (!strcasecmp(argv[n], "gb")) s_unit = 4;
      else fprintf(stderr, "Warning: Space unit %s not known\n", argv[n]);
    }
    else if (!strcmp(argv[n], "--plot") || !strcmp(argv[n], "-p"))
    {
      if (++n >= argc)
      {
        fprintf(stderr, "Error: Option %s requires a parameter\n", argv[n-1]);
        return 1;
      }
      if (!strcasecmp(argv[n], "n"))
      {
        if (gen_plot_m || gen_r2d || gen_r3d)
        {
          gen_plot_n = false, gen_plot_m = false, gen_r2d = false, gen_r3d = false;
          fprintf(stderr, "Warning: --plot n overrides preceding --plot option\n");
        }
        gen_plot_n = true;
      }
      else if (!strcasecmp(argv[n], "m"))
      {
        if (gen_plot_n || gen_r2d || gen_r3d)
        {
          gen_plot_n = false, gen_plot_m = false, gen_r2d = false, gen_r3d = false;
          fprintf(stderr, "Warning: --plot n overrides preceding --plot option\n");
        }
        gen_plot_m = true;
      }
      else if (!strcasecmp(argv[n], "r"))
      {
        if (gen_plot_n || gen_plot_m || gen_r3d)
        {
          gen_plot_n = false, gen_plot_m = false, gen_r2d = false, gen_r3d = false;
          fprintf(stderr, "Warning: --plot n overrides preceding --plot option\n");
        }
        gen_r2d = true;
      }
      else if (!strcasecmp(argv[n], "r+"))
      {
        if (gen_plot_n || gen_plot_m || gen_r2d)
        {
          gen_plot_n = false, gen_plot_m = false, gen_r2d = false, gen_r3d = false;
          fprintf(stderr, "Warning: --plot n overrides preceding --plot option\n");
        }
        gen_r3d = true;
      }
      else fprintf(stderr, "Warning: Plot data source '%s' not known\n", argv[n]);
    }
    else if (!strcmp(argv[n], "--latex") || !strcmp(argv[n], "-l")) gen_latex = true;
    else if (!strcmp(argv[n], "--quiet") || !strcmp(argv[n], "-q")) gen_plain = false;
    else if (!strcmp(argv[n], "--extra-rmap")) rmap_extra = true;
    else if (!strcmp(argv[n], "--model-info")) model_info = true;
    else if (!strcmp(argv[n], "--xrng") || !strcmp(argv[n], "--yrng"))
    {
      if (++n >= argc)
      {
        fprintf(stderr, "Error: Option %s requires a parameter\n", argv[n-1]);
        return 1;
      }
      low = strtod(argv[n], &p);
      if (*p != ':')
      {
        fprintf(stderr, "Error: Can't parse range %s\n", argv[n]);
        return 1;
      }
      high = strtod(p + 1, 0);
      if ((argv[n][2] == 'x' && xrange) || (argv[n][2] == 'y' && yrange))
      {
        fprintf(stderr, "Warning: %s %s overrides preceding %s\n", argv[n-1], argv[n], argv[n-1]);
      }
      if (argv[n][2] == 'x') xrange = true, xlow = low, xhigh = high;
      else if (argv[n][2] == 'y') yrange = true, ylow = low, yhigh = high;
    }
    else if (!strcmp(argv[n], "--opt") || !strcmp(argv[n], "-O"))
    {
      if (++n >= argc)
      {
        fprintf(stderr, "Error: Option %s requires a parameter\n", argv[n-1]);
        return 1;
      }
      if (!strcasecmp(argv[n], "all"))
      {
        all_opt_cfgs = true;
        if (cfgs.size() != 0) fprintf(stderr, "Warning: '--opt all' overrides all other --opt parameters\n");
      }
      else if (!all_opt_cfgs)
      {
        cfgs.push_back(strtoul(argv[n], &p, 0));
        if (*p && *p != ':')
        {
          fprintf(stderr, "Error: Can't parse optimization string %s\n", argv[n]);
          return 1;
        }
        cfg_titles.push_back(*p ? p + 1 : 0);
      }
      else fprintf(stderr, "Warning: '%s %s' ignored due to preceding '--opt all'\n", argv[n-1], argv[n]);
    }
    else if (!strcmp(argv[n], "--"))
    {
      ++n;
      break;
    }
    else if (*argv[n] != '-') break;
    else
    {
      fprintf(stderr, "Error: Unrecognized option '%s'\n", argv[n]);
      return 1;
    }
  }
  
  // Verify that we have a model type
  if (!random_model && !pm)
  {
    fprintf(stderr, "Error: No model type specified\n");
    return 1;
  }
  
  if(random_model && dumprelation)
  {
    fprintf(stderr, "Error: --dumprel cannot be used with random models\n");
    return 1;
  }
  
  if ((gen_r2d || gen_r3d) && !random_model)
  {
    fprintf(stderr, "Error: 'r' and 'r+' type plots can only be performed on random models\n");
    return 1;
  }
  
  // Default data type to tabulate
  if (!data_given && !dumprelation) tabulate[0] = true;
  
  // Warn if --extra-map and Quotient algorithm
  for (m = 0; rmap_extra && m < (int)cfgs.size(); ++m)
  {
    if (cfgs[m] & 2)
    {
      fprintf(stderr, "Warning: --extra-rmap should not be used with quotient algorithm.\n"
                      "         Results will be incorrect.\n");
      break;
    }
  }
  
  // No auto-unit in random mode warning
  if (random_model && (t_unit == 0 || s_unit == 0))
  {
    for (m = 0; m < 4; ++m)
    {
      if (tabulate[m] && t_unit == 0 && m < 3)
      {
        fprintf(stderr, "Warning: In random model mode, units are not adjusted automatically.\n"
                        "         Results will be displayed in seconds.\n");
        t_unit = 2;
        m = 3;
      }
      if (tabulate[m] && s_unit == 0 && m == 3)
      {
        fprintf(stderr, "Warning: In random model mode, units are not adjusted automatically.\n"
                        "         Results will be displayed in kilo-bytes.\n");
        s_unit = 2;
        break;
      }
    }
  }
  
  // Default optimization settings
  if (cfgs.size() == 0 && !all_opt_cfgs)
  {
    fprintf(stderr, "Warning: No optimizations specified; will benchmark with everything turned off\n");
  }
  
  // Parse input models and check if they are stochastic or not
  for (; n < argc; ++n)
  {
    f = fopen(argv[n], "rb");
    if (f)
    {
      files.push_back(argv[n]);
      
      if (random_model)
      {
        ParseRandom(f, &rm);
        if (rm.xtarget == 0xffff && rm.ztarget == 0xffff && (gen_plot_n || gen_plot_m || gen_r2d || gen_r3d))
        {
          fprintf(stderr, "Error: %s: Random models without variable parameters can't be plotted currently\n", argv[n]);
          return 1;
        }
        if ((rm.xtarget != 0xffff || rm.ztarget != 0xffff) && gen_latex)
        {
          fprintf(stderr, "Error: %s: LaTeX output is not supported for andom models with variable parameters\n", argv[n]);
          return 1;
        }
        if (rm.xtarget != 0xffff && rm.ztarget != 0xffff && gen_r2d)
        {
          fprintf(stderr, "Error: %s: Model has two variables but 2D plotting was selected\n", argv[n]);
          return 1;
        }
        if ((rm.xtarget == 0xffff || rm.ztarget == 0xffff) && gen_r3d)
        {
          fprintf(stderr, "Error: %s: Model has only one variable but 3D plotting was selected\n", argv[n]);
          return 1;
        }
        fclose(f);
        continue;
      }
      
      pm->Parse(f);
      
      fclose(f);
      
      check_label_file(argv[n], labelext, false);
    }
    else fprintf(stderr, "Warning: Model '%s' could not be opened, errno=%d\n", argv[n], errno);
  }
  
  // Make sure we got valid models above
  if (files.size() == 0)
  {
    fprintf(stderr, "Error: No models specified\n");
    return 1;
  }
  
  // Set up flagvector variable to pass to the benchmark class
  if (all_opt_cfgs)
  {
    flagvector = &allflagsvector[0];
    flagvectorsize = sizeof(allflagsvector) / sizeof(unsigned long);
    cfg_titles.clear();
    for (n = 0; n < (int)flagvectorsize; ++n) cfg_titles.push_back(0);
  }
  else if (cfgs.size() == 0)
  {
    flagvector = new unsigned long[1];
    *flagvector = 0;
    flagvectorsize = 1;
    cfg_titles.push_back("None");
  }
  else
  {
    flagvector = new unsigned long[cfgs.size()];
    flagvectorsize = cfgs.size();
    for (n = 0; n < (int)flagvectorsize; ++n) flagvector[n] = cfgs[n];
  }
  
  // If we are dealing with random models, branch out here. This should really
  // be in a separate function with all the variables passed in a struct but I don't
  // want to change the rest of the code too much.
  if (random_model)
  {
    RandomModelPlotInfo rmpi, *prmpi = new RandomModelPlotInfo[files.size()], *last_rmpi;
    last_rmpi = prmpi;
    std::map<unsigned int, std::set<RandomModelPlotInfo*> > modelmap2;
    std::map<unsigned int, std::set<RandomModelPlotInfo*> > modelmap3;
    std::map<unsigned int, std::set<RandomModelPlotInfo*> >::iterator mmi2;
    std::map<unsigned int, std::set<RandomModelPlotInfo*> >::iterator mmi3;
    
    for (i = files.begin(), m = 0; i != files.end(); ++i, ++m)
    {
      model_states.clear();
      model_transitions.clear();
      
      f = fopen(*i, "rb");
      ParseRandom(f, &rm);
      fclose(f);
      
      if (rm.xtarget != 0xffff && rm.ztarget != 0xffff)
      {
        rmpi.mdl = rm;
        sprintf(&rmpi.datasource[0], "%s_%%s_%s_%s.data", title, get_axis_name(rm.xtarget | 0x8000, true), get_axis_name(rm.ztarget | 0x8000, true));
        rmpi.id = m;
        memcpy(last_rmpi, &rmpi, sizeof(RandomModelPlotInfo));
        modelmap3[rm.xtarget | (rm.ztarget << 16)].insert(last_rmpi);
        ++last_rmpi;
      }
      else if (rm.xtarget != 0xffff && rm.ztarget == 0xffff)
      {
        rmpi.mdl = rm;
        sprintf(&rmpi.datasource[0], "%s_%%s_%s.data", title, get_axis_name(rm.xtarget | 0x8000, true));
        rmpi.id = m;
        memcpy(last_rmpi, &rmpi, sizeof(RandomModelPlotInfo));
        modelmap2[rm.xtarget].insert(last_rmpi);
        ++last_rmpi;
      }
      else if (rm.xtarget == 0xffff && rm.ztarget != 0xffff)
      {
        rmpi.mdl = rm;
        sprintf(&rmpi.datasource[0], "%s_%%s_%s.data", title, get_axis_name(rm.ztarget | 0x8000, true));
        rmpi.id = m;
        memcpy(last_rmpi, &rmpi, sizeof(RandomModelPlotInfo));
        modelmap3[rm.ztarget].insert(last_rmpi);
        ++last_rmpi;
      }
      
      zsteps = (rm.ztarget == 0xffff ? 1 : rm.steps);
      rsteps = (rm.xtarget == 0xffff ? 1 : rm.steps) * zsteps;
      
      for (n = 0; n < int(sizeof(tabulate) / sizeof(tabulate[0])); ++n)
      {
        if (tabulate[n]) result[n] = new double[flagvectorsize * rsteps];
      }
      
      if (tabulate[3] && rmap_extra) result_rmap = new double[flagvectorsize * rsteps];
      
      for (j = 0; j < (rm.xtarget == 0xffff ? 1 : (int)rm.steps); ++j)
      {
        for (k = 0; k < (rm.ztarget == 0xffff ? 1 : (int)rm.steps); ++k)
        {
          switch (rm.xtarget)
          {
          case 0x0001: rm.n = (int)(rm.xstart + ((rm.xend - rm.xstart) * j / (rm.steps - 1))); break;
          case 0x0002: rm.a = (int)(rm.xstart + ((rm.xend - rm.xstart) * j / (rm.steps - 1))); break;
          case 0x0003: rm.b = (int)(rm.xstart + ((rm.xend - rm.xstart) * j / (rm.steps - 1))); break;
          case 0x0004: rm.c = (int)(rm.xstart + ((rm.xend - rm.xstart) * j / (rm.steps - 1))); break;
          case 0x0100: rm.fb = rm.xstart + ((rm.xend - rm.xstart) * j / (rm.steps - 1)); break;
          case 0x0200: rm.lb = rm.xstart + ((rm.xend - rm.xstart) * j / (rm.steps - 1)); break;
          case 0x0300: rm.cb = rm.xstart + ((rm.xend - rm.xstart) * j / (rm.steps - 1)); break;
          case 0x0400: rm.pb = rm.xstart + ((rm.xend - rm.xstart) * j / (rm.steps - 1)); break;
          case 0x0500: rm.sb = rm.xstart + ((rm.xend - rm.xstart) * j / (rm.steps - 1)); break;
          case 0x0600: rm.labels = (unsigned int)(rm.xstart + ((rm.xend - rm.xstart) * j / (rm.steps - 1))); break;
          }
          
          switch (rm.ztarget)
          {
          case 0x0001: rm.n = (int)(rm.zstart + ((rm.zend - rm.zstart) * k / (rm.steps - 1))); break;
          case 0x0002: rm.a = (int)(rm.zstart + ((rm.zend - rm.zstart) * k / (rm.steps - 1))); break;
          case 0x0003: rm.b = (int)(rm.zstart + ((rm.zend - rm.zstart) * k / (rm.steps - 1))); break;
          case 0x0004: rm.c = (int)(rm.zstart + ((rm.zend - rm.zstart) * k / (rm.steps - 1))); break;
          case 0x0100: rm.fb = rm.zstart + ((rm.zend - rm.zstart) * k / (rm.steps - 1)); break;
          case 0x0200: rm.lb = rm.zstart + ((rm.zend - rm.zstart) * k / (rm.steps - 1)); break;
          case 0x0300: rm.cb = rm.zstart + ((rm.zend - rm.zstart) * k / (rm.steps - 1)); break;
          case 0x0400: rm.pb = rm.zstart + ((rm.zend - rm.zstart) * k / (rm.steps - 1)); break;
          case 0x0500: rm.sb = rm.zstart + ((rm.zend - rm.zstart) * k / (rm.steps - 1)); break;
          case 0x0600: rm.labels = (unsigned int)(rm.zstart + ((rm.zend - rm.zstart) * k / (rm.steps - 1))); break;
          }
          
          fprintf(stderr, "%s [%d %d %d %d %+.3f %.3f %.3f %.3f %.3f] [%d/%d]: ", *i, rm.n, rm.a, rm.b, rm.c, rm.fb, rm.lb, rm.cb, rm.pb, rm.sb, (j * (rm.ztarget == 0xffff ? 1 : rm.steps)) + k + 1, rsteps);
          benchmark_random_model(&rm, flagvector, flagvectorsize, averages, fpprecision, &transsum, tabulate, result, rmap_extra, result_rmap, rsteps, (j * zsteps) + k);
          fprintf(stderr, "\n");

          if (model_info)
          {
            model_states.push_back(rm.n);
            model_transitions.push_back(transsum);
          }
        }
      }
      fprintf(stderr, "\n");
      
      if (result[0]) unit_transform(result[0], flagvectorsize * rsteps, t_unit, false);
      if (result[1]) unit_transform(result[1], flagvectorsize * rsteps, t_unit, false);
      if (result[2]) unit_transform(result[2], flagvectorsize * rsteps, t_unit, false);
      if (result[3]) unit_transform(result[3], flagvectorsize * rsteps, s_unit, true);
      if (result_rmap) unit_transform(result_rmap, flagvectorsize * rsteps, s_unit, true);
      
      if (gen_plain)
      {
        for (n = 0; n < int(sizeof(tabulate) / sizeof(tabulate[0])); ++n)
        {
          if (tabulate[n])
          {
            gen_rnd_plain(stdout, &rm, result[n], rsteps, flagvectorsize, flagvector, cfg_titles,
                          (n == 3 ? result_rmap : 0), "Map size", (n > 3 ? 0 : precision), model_states,
                          model_transitions, model_actions, basenameptr(*i), datanames[n]);
          }
        }
      }
      
      if (gen_r2d)
      {
        plot_key = new double[rm.steps];
        
        for (n = 0; n < (int)rm.steps; ++n)
        {
          plot_key[n] = (rm.xtarget == 0xffff ? (rm.zstart + ((rm.zend - rm.zstart) * n / (rm.steps - 1))) : (rm.xstart + ((rm.xend - rm.xstart) * n / (rm.steps - 1))));
        }
        
        for (n = 0; n < int(sizeof(tabulate) / sizeof(tabulate[0])); ++n)
        {
          if (tabulate[n])
          {
            snprintf(&datasource[0], 160, "%s_%s_%s.data", title, datanames[n], get_axis_name((rm.xtarget == 0xffff ? rm.ztarget : rm.xtarget) | 0x8000, true));
            f = fopen(&datasource[0], "wb");
            gen_rnd_plain(f, &rm, result[n], rsteps, flagvectorsize, flagvector, cfg_titles,
                          (n == 3 ? result_rmap : 0), "Map size", (n > 3 ? 0 : 10), model_states,
                          model_transitions, model_actions, basenameptr(*i), datanames[n], true);
            fclose(f);
          }
        }
        
        delete [] plot_key;

        for (n = 0; n < int(sizeof(tabulate) / sizeof(tabulate[0])); ++n)
        {
          if (tabulate[n])
          {
            for (mmi2 = modelmap2.begin(); mmi2 != modelmap2.end(); ++mmi2)
            {
              snprintf(&filename[0], 160, "%s_%s_%s.gnuplot", title, datanames[n], get_axis_name(mmi2->first | 0x8000, true));
              f = fopen(&filename[0], "wb");
              snprintf(&filename[0], 160, "%s_%s_%s.eps", title, datanames[n], get_axis_name(mmi2->first | 0x8000, true));
              rnd_2d_gnuplot_macro(f, get_axis_name(mmi2->first, false), datanames[n],
                                  (n == 3 ? s_units[s_unit] : (n > 3 ? "" : t_units[t_unit])),
                                  title, &filename[0], datanames[n], flagvector, flagvectorsize, cfg_titles, mmi2->second);
              fclose(f);
            }
          }
        }
      }
      
      if (gen_r3d)
      {
        plot_key = new double[2 * rsteps];
        
        for (n = 0; n < (int)rm.steps; ++n)
        {
          for (k = 0; k < (int)rm.steps; ++k)
          {
            plot_key[(n * rm.steps) + k] = rm.xstart + ((rm.xend - rm.xstart) * n / (rm.steps - 1));
            plot_key[(rsteps) + (k * rm.steps) + n] = rm.zstart + ((rm.zend - rm.zstart) * n / (rm.steps - 1));
          }
        }
        
        for (n = 0; n < int(sizeof(tabulate) / sizeof(tabulate[0])); ++n)
        {
          if (tabulate[n])
          {
            snprintf(&datasource[0], 160, "%s_%s_%s_%s.data", title, datanames[n], get_axis_name(rm.xtarget | 0x8000, true), get_axis_name(rm.ztarget | 0x8000, true));
            f = fopen(&datasource[0], "wb");
            generate_gnuplot_data(f, result[n], plot_key, flagvectorsize, rsteps, 2, rm.steps, true);
            fclose(f);
          }
        }
        
        delete [] plot_key;

        for (n = 0; n < int(sizeof(tabulate) / sizeof(tabulate[0])); ++n)
        {
          if (tabulate[n])
          {
            for (mmi3 = modelmap3.begin(); mmi3 != modelmap3.end(); ++mmi3)
            {
              snprintf(&filename[0], 160, "%s_%s_%s_%s.gnuplot", title, datanames[n], get_axis_name((mmi3->first & 0xffff) | 0x8000, true), get_axis_name((mmi3->first >> 16) | 0x8000, true));
              f = fopen(&filename[0], "wb");
              snprintf(&filename[0], 160, "%s_%s_%s_%s.eps", title, datanames[n], get_axis_name((mmi3->first & 0xffff) | 0x8000, true), get_axis_name((mmi3->first >> 16) | 0x8000, true));
              rnd_3d_gnuplot_macro(f, get_axis_name(mmi3->first & 0xffff, false), datanames[n],
                                  get_axis_name(mmi3->first >> 16, false),
                                  (n == 3 ? s_units[s_unit] : (n > 3 ? "" : t_units[t_unit])),
                                  title, &filename[0], datanames[n], flagvector, flagvectorsize, cfg_titles, mmi3->second);
              fclose(f);
            }
          }
        }
      }
      
      for (n = 0; n < int(sizeof(tabulate) / sizeof(tabulate[0])); ++n)
      {
        if (tabulate[n]) delete [] result[n];
      }
      
      if (tabulate[3] && rmap_extra) delete [] result_rmap;
    }

    delete [] prmpi;
          
    if (!all_opt_cfgs) delete [] flagvector;
    
    return 0;
  }
  
  // Prepare result tables
  for (n = 0; n < int(sizeof(tabulate) / sizeof(tabulate[0])); ++n)
  {
    if (tabulate[n]) result[n] = new double[flagvectorsize * files.size()];
  }
  
  if (tabulate[3] && rmap_extra) result_rmap = new double[files.size()];
  
  // Go through the models and benchmark for each
  for (i = files.begin(), m = 0; i != files.end(); ++i, ++m)
  {
    f = fopen(*i, "rb");
    pm->Parse(f);
    fclose(f);
    
    fprintf(stderr, "%s: ", *i);
    
    if (model_info)
    {
      model_states.push_back(pm->States());
      model_transitions.push_back(pm->Transitions());
      if (pm->Type() == ProbabilisticModel::PA) model_actions.push_back(((ProbabilisticAutomaton*)pm)->na);
    }
    
    // Copy results into the appropriate tables
    for (n = 0; n < (int)flagvectorsize; ++n)
    {
      if (check_label_file(*i, labelext, true))
        sprintf(&command[0], "%s/benchone_%lu strong %s %s%s %e %u %s", &myself[0], flagvector[n], &modeltype[0], *i, labelext, fpprecision, averages, *i);
      else
        sprintf(&command[0], "%s/benchone_%lu strong %s \"\" %e %u %s", &myself[0], flagvector[n], &modeltype[0], fpprecision, averages, *i);
      
      if(dumprelation && !dumpsuccessful)
        sprintf(&command[strlen(&command[0])], " %s_rel.txt", title);

      f = popen(&command[0], "r");
      if (fscanf(f, "%le %le %le %u %u %u %u %u %u %u %lu %lu %lu %lu %lu %u %u %u %u\n",
        &utime, &stime, &rtime, &stats.num_partitions,
        &stats.num_iterations, &stats.num_initial_pairs, &stats.num_final_pairs,
        &stats.num_maxflow,
        &stats.num_p_invariant_fails, &stats.num_sig_arc_fails, &stats.mem_relation_map,
        &stats.mem_partition_map, &stats.mem_relation,
        &stats.mem_maxflow, &stats.mem_model, &stats.min_complexity, &stats.max_complexity,
        &stats.num_nets_cached, &stats.num_cache_hits) < 19)
      {
        utime = -1, rtime = -1, stime = -1;
        memset(&stats, 0, sizeof(stats));
        fprintf(stderr, "Warning: Benchmark of %s, cfg %02lx failed\n", *i, flagvector[n]);
      }
      else dumpsuccessful = true;
      pclose(f);

      if (tabulate[0]) result[0][n + (m * flagvectorsize)] = utime;
      if (tabulate[1]) result[1][n + (m * flagvectorsize)] = stime;
      if (tabulate[2]) result[2][n + (m * flagvectorsize)] = rtime;
      if (tabulate[3])
      {
        result[3][n + (m * flagvectorsize)] = double(stats.mem_relation + stats.mem_maxflow
                                                       + stats.mem_partition_map + stats.mem_model);
        if (!rmap_extra) result[3][n + (m * flagvectorsize)] += double(stats.mem_relation_map);
      }
      if (tabulate[4]) result[4][n + (m * flagvectorsize)] = double(stats.num_initial_pairs);
      if (tabulate[5]) result[5][n + (m * flagvectorsize)] = double(stats.num_final_pairs);
      if (tabulate[6]) result[6][n + (m * flagvectorsize)] = double(stats.num_partitions);
      if (tabulate[7]) result[7][n + (m * flagvectorsize)] = double(stats.num_iterations);
      if (tabulate[8]) result[8][n + (m * flagvectorsize)] = double(stats.num_maxflow);
      if (tabulate[9]) result[9][n + (m * flagvectorsize)] = double(stats.num_p_invariant_fails);
      if (tabulate[10]) result[10][n + (m * flagvectorsize)] = double(stats.num_sig_arc_fails);
      if (tabulate[11]) result[11][n + (m * flagvectorsize)] = double(stats.min_complexity);
      if (tabulate[12]) result[12][n + (m * flagvectorsize)] = double(stats.max_complexity);
      if (tabulate[13]) result[13][n + (m * flagvectorsize)] = double(stats.num_nets_cached);
      if (tabulate[14]) result[14][n + (m * flagvectorsize)] = double(stats.num_cache_hits);
    }
    
    if (tabulate[3] && rmap_extra) result_rmap[m] = double(stats.mem_relation_map);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  
  // Adjust time unit (determine best choice if not given)
  average = 0.0;
  total = 0;
  for (n = 0; n < 3; ++n)
  {
    if (tabulate[n])
    {
      for (m = 0; m < (int)(flagvectorsize * files.size()); ++m)
      {
        average += result[n][m];
      }
      total += flagvectorsize * files.size();
    }
  }
  if (total > 0)
  {
    average /= total;
    if (t_unit == 0) t_unit = find_best_time_unit(average);
    if (result[0]) unit_transform(result[0], flagvectorsize * files.size(), t_unit, false);
    if (result[1]) unit_transform(result[1], flagvectorsize * files.size(), t_unit, false);
    if (result[2]) unit_transform(result[2], flagvectorsize * files.size(), t_unit, false);
  }
  
  // Adjust space unit (determine best choice if not given)
  average = 0.0;
  if (tabulate[3])
  {
    for (m = 0; m < (int)(flagvectorsize * files.size()); ++m) average += result[3][m];
    average /= flagvectorsize * files.size();
    if (s_unit == 0) s_unit = find_best_space_unit(average);
    if (result[3]) unit_transform(result[3], flagvectorsize * files.size(), s_unit, true);
    if (result_rmap) unit_transform(result_rmap, files.size(), s_unit, true);
  }
  
  // Generate LaTeX output if requested
  if (gen_latex)
  {
    for (n = 0; n < int(sizeof(tabulate) / sizeof(tabulate[0])); ++n)
    {
      if (tabulate[n])
      {
        snprintf(&filename[0], 160, "%s_%s.tex", title, datanames[n]);
        f = fopen(&filename[0], "wb");
        generate_latex(f, result[n], flagvectorsize, files.size(), flagvector, cfg_titles, files,
                       (n == 3 ? result_rmap : 0), "Map size", (n > 3 ? 0 : precision), title,
                       (n == 3 ? s_units[s_unit] : (n > 3 ? "" : t_units[t_unit])), model_states,
                       model_transitions, model_actions);
        fclose(f);
      }
    }
  }
  
  // Generate plot if requested
  if (gen_plot_n || gen_plot_m)
  {
    plot_key = new double[files.size()];
    
    // Set up plot_key (see generate_plot_data() for more information)
    for (n = 0, i = files.begin(); i != files.end(); ++n, ++i)
    {
      f = fopen(*i, "rb");
      pm->Parse(f);
      fclose(f);
      plot_key[n] = (gen_plot_n ? pm->States() : pm->Transitions());
    }
    
    // Generate output files
    for (n = 0; n < int(sizeof(tabulate) / sizeof(tabulate[0])); ++n)
    {
      if (tabulate[n])
      {
        snprintf(&datasource[0], 160, "%s_%s.data", title, datanames[n]);
        f = fopen(&datasource[0], "wb");
        generate_gnuplot_data(f, result[n], plot_key, flagvectorsize, files.size(), 1);
        fclose(f);
        
        snprintf(&filename[0], 160, "%s_%s.gnuplot", title, datanames[n]);
        f = fopen(&filename[0], "wb");
        snprintf(&filename[0], 160, "%s_%s.eps", title, datanames[n]);
        generate_gnuplot_macro(f, (gen_plot_n ? "States" : "Transitions"), datanames[n],
                            (n == 3 ? s_units[s_unit] : (n > 3 ? "" : t_units[t_unit])),
                            (xrange ? xlow : find_data_range(plot_key, files.size(), true)),
                            (xrange ? xhigh : find_data_range(plot_key, files.size(), false)),
                            (yrange ? ylow : find_data_range(result[n], flagvectorsize * files.size(), true)),
                            (yrange ? yhigh : find_data_range(result[n], flagvectorsize * files.size(), false)),
                            title, &datasource[0], &filename[0], flagvector, cfg_titles, flagvectorsize);
        fclose(f);
      }
    }
    delete [] plot_key;
  }
  
  // Print plain-text table to stdout
  if (gen_plain)
  {
    int title_len, title_fill;
    for (n = 0; n < int(sizeof(tabulate) / sizeof(tabulate[0])); ++n)
    {
      if (tabulate[n])
      {
        title_len = strlen(title) + strlen(datanames[n]) + strlen(n == 3 ? s_units[s_unit] : (n > 3 ? "" : t_units[t_unit])) + 7;
        title_fill = find_table_width(result[n], flagvectorsize, files.size(), flagvector, cfg_titles, files,
                                      (n == 3 ? result_rmap : 0), "Map size", (n > 3 ? 0 : precision), model_states,
                                      model_transitions, model_actions) - title_len;
        
        if (title_fill <= 4) title_fill = 4;
        total = (title_fill / 2) + (title_fill % 2);
        while (total--) fputs("-", stdout);
        fprintf(stdout, " %s: %s (%s) ", title, datanames[n], (n == 3 ? s_units[s_unit] : (n > 3 ? "" : t_units[t_unit])));
        total = title_fill / 2;
        while (total--) fputs("-", stdout);
        fputs("\n", stdout);
        
        generate_plain(stdout, result[n], flagvectorsize, files.size(), flagvector, cfg_titles, files,
                       (n == 3 ? result_rmap : 0), "Map size", (n > 3 ? 0 : precision), model_states,
                       model_transitions, model_actions);
        
        total = title_len + title_fill;
        while (total--) fputs("-", stdout);
        fputs("\n\n", stdout);
      }
    }
  }
  
  // Clean up
  if (!all_opt_cfgs) delete [] flagvector;
  
  for (n = 0; n < int(sizeof(result) / sizeof(result[0])); ++n) if (result[n]) delete [] result[n];
  if (tabulate[3] && rmap_extra) delete [] result_rmap;
  
  //if (pss) delete pss;
  if (pm) delete pm;
  
  return 0;
}
