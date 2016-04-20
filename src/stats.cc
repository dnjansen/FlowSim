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


#include "stats.h"

// definition of static variables for memory
ssize_t SimulationStatistics::mem_current, SimulationStatistics::mem_max;
