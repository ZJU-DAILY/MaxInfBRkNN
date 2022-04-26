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

#include "ObjectIndexTuple.h"

#include <sstream>

ObjectIndexTuple::ObjectIndexTuple(std::string networkName, int nodes, int edges, std::string indexMethod, 
                             double constructionTimeMs, double memoryUsageMB, std::string objSetType,
                             double objSetDensity,int objSetVariable, int objSetSize): 
    Tuple(networkName,nodes,edges), indexMethod(indexMethod), constructionTimeMs(constructionTimeMs), 
    memoryUsageMB(memoryUsageMB), objSetType(objSetType), objSetDensity(objSetDensity), objSetVariable(objSetVariable), objSetSize(objSetSize) {}
    
std::string ObjectIndexTuple::getTupleString()
{
    std::stringstream ss;
    ss.precision(4);
    ss << std::fixed;
    ss << Tuple::getTupleString() << "\t";
    ss << this->indexMethod << "\t";
    ss << this->constructionTimeMs << "\t";
    ss << this->memoryUsageMB << "\t";
    ss << this->objSetType << "\t";
    ss.precision(5);
    ss << std::fixed;
    ss << this->objSetDensity << "\t";
    ss << this->objSetVariable << "\t";
    ss << this->objSetSize << "\t";
    ss << this->getAdditionalFields();
    return ss.str();
}

std::string ObjectIndexTuple::getMultilineTupleString()
{
    std::stringstream ss;
    ss.precision(4);
    ss << std::fixed;
    ss << "\nObject Index Construction Results: \n";
    ss << Tuple::getMultilineTupleString();
    ss << "Object Index Method: " << this->indexMethod << "\n";
    ss << "Construction Time: " << this->constructionTimeMs << "ms\n";
    ss << "Memory Usage: " << this->memoryUsageMB << "MB\n";
    ss.precision(5);
    ss << std::fixed;
    ss << "Object Set: Type=" << this->objSetType << ", Density=" << this->objSetDensity << ", Variable=" << this->objSetVariable << ", Size=" << this->objSetSize << "\n";
    return ss.str();
}
