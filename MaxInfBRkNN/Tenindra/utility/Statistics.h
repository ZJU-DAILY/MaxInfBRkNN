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

#ifndef _STATISTICS_H
#define _STATISTICS_H

#include "../tuple/Tuple.h"

#include <string>
#include <unordered_map>

class Statistics {
    
    public:
        Statistics();
        void clear();
        void initialiseStatistic(std::string desc, double value);
        void incrementStatistic(std::string desc, double value);
        void mergeStatistics(Statistics& inputStatistics);
        void normalizeStatistics(int normalizationFactor);
        void populateTupleFields(Tuple& tuple, int decimalPlaces = 3);
        double getStatistic(std::string desc);
        
        std::unordered_map<std::string,double> statistics;
        
};

#endif // _STATISTICS_H
