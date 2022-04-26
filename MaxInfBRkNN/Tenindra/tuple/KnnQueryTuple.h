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

#ifndef _KNNQUERYTUPLE_H
#define _KNNQUERYTUPLE_H

#include "Tuple.h"
#include "../common.h"

#include <string>
#include <vector>

class KnnQueryTuple: public Tuple {
    
public:
    std::string getTupleString();
    std::string getMultilineTupleString();
    KnnQueryTuple(std::string networkName, int nodes, int edges, NodeID queryNodeID, std::string queryMethod, int k,
                  double queryTimeMs, std::string objSetType, double objSetDensity, int objSetVariable, int objSetSize, 
                  std::vector<NodeID>& kNNs, std::vector<EdgeWeight>& kNNDistances);
    std::string getResultsString();
    
private:
    NodeID queryNodeID;
    std::string queryMethod;
    int k;
    double queryTimeMs;
    std::string objSetType;
    double objSetDensity;
    int objSetVariable, objSetSize;
    std::string kNNsDistancesString;
};

#endif // _KNNQUERYTUPLE_H
