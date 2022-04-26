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

#include "KnnQueryTuple.h"

#include "../utility/utility.h"
#include <sstream>

KnnQueryTuple::KnnQueryTuple(std::string networkName, int nodes, int edges, NodeID queryNodeID, std::string queryMethod, 
                             int k, double queryTimeMs, std::string objSetType, double objSetDensity, int objSetVariable, int objSetSize, 
                             std::vector<NodeID>& kNNs, std::vector<EdgeWeight>& kNNDistances): 
    Tuple(networkName,nodes,edges), queryNodeID(queryNodeID), queryMethod(queryMethod), k(k), queryTimeMs(queryTimeMs), 
    objSetType(objSetType), objSetDensity(objSetDensity), objSetVariable(objSetVariable), objSetSize(objSetSize) {

    // Create string to represented kNN result pairs (object,distance)
    this->kNNsDistancesString = utility::implodeNodesAndDistances(kNNs,kNNDistances);
}

std::string KnnQueryTuple::getTupleString()
{
    std::stringstream ss;
    ss.precision(3);
    ss << std::fixed;
    ss << Tuple::getTupleString() << "\t";
    ss << this->queryMethod << "\t";
    ss << this->k << "\t";
    ss << this->queryTimeMs << "\t";
    ss << this->objSetType << "\t";
    ss.precision(5);
    ss << std::fixed;
    ss << this->objSetDensity << "\t";
    ss << this->objSetVariable << "\t";
    ss << this->objSetSize << "\t";
    ss << this->queryNodeID << "\t";
    ss << this->getAdditionalFields();
    return ss.str();
}

std::string KnnQueryTuple::getMultilineTupleString()
{
    std::stringstream ss;
    ss.precision(3);
    ss << std::fixed;
    ss << "\nkNN Query Results: \n";
    ss << Tuple::getMultilineTupleString();
    ss << "Query Method: " << this->queryMethod << "\n";
    ss << "k: " << this->k << "\n";
    ss << "Query Time: " << this->queryTimeMs << "ms\n";
    ss.precision(5);
    ss << std::fixed;
    ss << "Object Set: Type=" << this->objSetType << ", Density=" << this->objSetDensity << ", Variable=" << this->objSetVariable << ", Size=" << this->objSetSize << "\n";
    ss << "Query Node: " << this->queryNodeID << "\n";    
    ss << "kNNs and Distances: " + this->kNNsDistancesString + "\n";
    return ss.str();
}

std::string KnnQueryTuple::getResultsString()
{
    std::stringstream ss;
    ss << this->queryNodeID << "\t";
    ss << this->queryMethod << "\t";
    ss << this->kNNsDistancesString;
    return ss.str();
}
