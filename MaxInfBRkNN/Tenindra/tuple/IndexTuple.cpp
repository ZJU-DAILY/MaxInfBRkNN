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

#include "IndexTuple.h"

#include <sstream>

IndexTuple::IndexTuple(std::string networkName, int nodes, int edges, std::string indexMethod, double constructionTimeMs, double memoryUsageMB): 
    Tuple(networkName,nodes,edges), indexMethod(indexMethod), constructionTimeMs(constructionTimeMs), memoryUsageMB(memoryUsageMB) {}
    
std::string IndexTuple::getTupleString()
{
    std::stringstream ss;
    ss.precision(3);
    ss << std::fixed;
    ss << Tuple::getTupleString() << "\t";
    ss << this->indexMethod << "\t";
    ss << this->constructionTimeMs << "\t";
    ss << this->memoryUsageMB << "\t";
    ss << this->getAdditionalFields();
    return ss.str();
}

std::string IndexTuple::getMultilineTupleString()
{
    std::stringstream ss;
    ss.precision(3);
    ss << std::fixed;
    ss << "\nIndex Construction Results: \n";
    ss << Tuple::getMultilineTupleString();
    ss << "Index Method: " << this->indexMethod << "\n";
    ss << "Construction Time: " << this->constructionTimeMs << "ms\n";
    ss << "Memory Usage: " << this->memoryUsageMB << "MB\n";
    return ss.str();
}
