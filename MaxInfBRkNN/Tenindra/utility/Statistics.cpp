/* Copyright (C) 2015 Tenindra Abeywickrama
 *
 * This file is part of Road Network kNN Experimental Evaluation.
 *
 * Road Network kNN Experimental Evaluation is free software; you can
 * redistribute it and/or modify it under the terms of the GNU Affero 
 * General Public License as published by the Free Software Foundation; 
 * either version 3 of the License, or (at your option) any later version.
 *
 * Road Network kNN Experimental Evaluation is distributed in the hope 
 * that it will be useful, but WITHOUT ANY WARRANTY; without even the 
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE. See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public 
 * License along with Road Network kNN Experimental Evaluation; see 
 * LICENSE.txt; if not, see <http://www.gnu.org/licenses/>.
 */

#include "Statistics.h"

#include <sstream>
#include <assert.h>

Statistics::Statistics() {}

void Statistics::clear()
{
    this->statistics.clear();
}

void Statistics::initialiseStatistic(std::string desc, double value)
{
    statistics[desc] = value;
}

// Assumes it has already been intialised
void Statistics::incrementStatistic(std::string desc, double value)
{
    assert(statistics.find(desc) != statistics.end() && "Statistic to be increment has not been initialised");
    statistics[desc] += value;
}

void Statistics::mergeStatistics(Statistics& inputStatistics)
{
    for (auto it = inputStatistics.statistics.begin(); it != inputStatistics.statistics.end(); ++it) {
        if (this->statistics.find(it->first) == this->statistics.end()) {
            this->statistics[it->first] = inputStatistics.statistics[it->first];
        } else {
            this->statistics[it->first] += inputStatistics.statistics[it->first];
        }
    }
}

void Statistics::normalizeStatistics(int normalizationFactor)
{
    for (auto it = this->statistics.begin(); it != this->statistics.end(); ++it) {
        this->statistics[it->first] = it->second/normalizationFactor;
    }
}

void Statistics::populateTupleFields(Tuple& tuple, int decimalPlaces)
{
    std::stringstream ss;
    for (auto it = this->statistics.begin(); it != this->statistics.end(); ++it) {
        ss.str("");
        ss.precision(decimalPlaces);
        ss << std::fixed;
        ss << it->second;
        tuple.addSupplementaryFields(it->first,ss.str());
    }
}

double Statistics::getStatistic(std::string desc)
{
    assert(statistics.find(desc) != statistics.end() && "Statistic to return has not been initialised");
    return statistics[desc];
}
