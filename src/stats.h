/*****************************************************************************/
/*!
 *   Copyright 2009-2014 Jonathan Bogdoll, Holger Hermanns, Lijun Zhang,
 *                       David N. Jansen
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


// The intended use of this file is as follows.
// Immediately before some memory on the heap is allocated, one should call
// RegisterMemAlloc(). (Before, because the constructor may also allocate some
// temporary storage, which might affect the global peak.)
// Immediately after some memory is freed, one should call RegisterMemFree().
// However, if it is difficult to find the size of the allocated data (e. g.
// the size of a set changes when one erases its elements), I sometimes call
// RegisterMemFree() a bit early.
// To take into account the overhead of sets, maps, and vectors, I #defined
// three constants: SET_OVERHEAD, MAP_OVERHEAD and VEC_OVERHEAD.

// At the beginning of the measurement, one calls ResetStats(). At the end,
// a call to CollectStats() will add the global space peak to the member
// mem_model of the instance. Some memory may be specialized and needs to be
// accounted for in other members of SimulationStatistics; then one should
// subtract the respective value from mem_model and add it to the correct
// member. For an example, see CompactMaxFlow<double>::CollectStats().

// However, for the standard analysis provided by the utility flowsim, it does
// not make much of a difference; only mem_relation_map will be handled
// separately (upon request), and flowsim sums the other memory amounts before
// presenting the table or graph.


#ifndef STATS_H
#define STATS_H

/* #define REPORT_MEMORY */

#ifdef DEBUG

#include <cstdio>
#include <cstring>

// __tree_node is the relevant type in my version of the g++ library.
#define SET_OVERHEAD (sizeof(std::__tree_node<int,void*>)-sizeof(int))
#define MAP_OVERHEAD (sizeof(std::__tree_node<int,void*>)-sizeof(int))
#define VEC_OVERHEAD ((size_t) 0)

  struct SimulationStatistics
  {
    unsigned int num_partitions, num_iterations, num_initial_pairs, num_final_pairs,
      num_maxflow, num_p_invariant_fails, num_sig_arc_fails, min_complexity,
      max_complexity, num_nets_cached, num_cache_hits;
    ssize_t mem_relation_map, mem_partition_map, mem_relation, mem_maxflow,
                mem_model;

    inline void ResetStats() {
      memset(this, '\0', sizeof(*this));
      mem_current = 0;
      mem_max = 0;
    }

    inline void CollectStats() {
      if (0 != mem_current)
        fprintf(stderr,"**** Memory leak: mem_current is %zd B ****\n",
                    mem_current);
      mem_model += mem_max;
      //fprintf(stderr, "At the peak, %zu B = %.3g MB was used in total.\n",
      //            mem_max, mem_max * (1/1048576.0));
      //fprintf(stderr, "Up to %zu B = %.3g MB was used for flow networks.\n",
      //            mem_maxflow, mem_maxflow * (1/1048576.0));
      //fprintf(stderr, "Up to %zu B = %.3g MB was used for relation maps.\n",
      //            mem_relation_map, mem_relation_map * (1/1048576.0));
      mem_current = 0;
      mem_max = 0;
    }

  private:
    friend void RegisterMemAlloc(size_t size);
    friend void RegisterMemFree(size_t size);
    friend ssize_t MaxMemUsed(void);
    friend ssize_t CurMemUsed(void);

    static ssize_t mem_current, mem_max;
  };

  inline void RegisterMemAlloc(size_t size)
  {
    SimulationStatistics::mem_current += size;
    if (SimulationStatistics::mem_max < SimulationStatistics::mem_current) {
#ifdef REPORT_MEMORY
      if (SimulationStatistics::mem_max >> 32
				< SimulationStatistics::mem_current >> 32)
        fprintf(stderr, "%zd GB allocated\n",
			(SimulationStatistics::mem_current + (1 << 29)) >> 30);
#endif
      SimulationStatistics::mem_max = SimulationStatistics::mem_current;
    }
  }

  inline void RegisterMemFree(size_t size)
  {
    SimulationStatistics::mem_current -= size;
  }

  inline ssize_t MaxMemUsed(void) { return SimulationStatistics::mem_max; }
  inline ssize_t CurMemUsed(void) { return SimulationStatistics::mem_current; }
#else//DEBUG
  inline void RegisterMemAlloc(size_t /*size*/) {}
  inline void RegisterMemFree(size_t /*size*/) {}
#endif//DEBUG

#endif//STATS_H
