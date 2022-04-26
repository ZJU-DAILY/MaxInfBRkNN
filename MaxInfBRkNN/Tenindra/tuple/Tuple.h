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

#ifndef _TUPLE_H
#define _TUPLE_H

#include <string>
#include <vector>

class Tuple {
    
public:
    virtual std::string getTupleString();
    virtual std::string getMultilineTupleString();
    virtual ~Tuple() {};
    Tuple(std::string networkName, int nodes, int edge);
    void addSupplementaryFields(std::string desc, std::string data);
    std::string getAdditionalFields();
    void setAdditionalFields(std::vector<std::string>& fields);
    
protected:
    std::string networkName;
    int nodes;
    int edges;
    std::string supplementaryData;
    std::vector<std::string> additionalFields;
};

#endif // _TUPLE_H
