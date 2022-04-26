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

#ifndef _OBJECTINDEXTUPLE_H
#define _OBJECTINDEXTUPLE_H

#include "Tuple.h"

#include <string>

class ObjectIndexTuple: public Tuple {
    
public:
    std::string getTupleString();
    std::string getMultilineTupleString();
    ObjectIndexTuple(std::string networkName, int nodes, int edges, std::string indexMethod, 
                     double constructionTimeMs, double memoryUsageMB, std::string objSetType, 
                     double objSetDensity, int objSetVariable, int objSetSize);
    
protected:
    std::string indexMethod;
    double constructionTimeMs;
    double memoryUsageMB;
    std::string objSetType;
    double objSetDensity;
    int objSetVariable;
    int objSetSize;
};

#endif // _OBJECTINDEXTUPLE_H
