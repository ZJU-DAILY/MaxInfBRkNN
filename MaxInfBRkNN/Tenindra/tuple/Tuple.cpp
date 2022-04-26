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

#include "Tuple.h"

#include <sstream>

Tuple::Tuple(std::string networkName, int nodes, int edges)
{
    this->networkName = networkName;
    this->nodes = nodes;
    this->edges = edges;
    this->supplementaryData = "none";
}

std::string Tuple::getTupleString()
{
    std::stringstream ss;
    ss.precision(3);
    ss << std::fixed;
    ss << this->networkName << "\t";
    ss << this->nodes << "\t";
    ss << this->edges << "\t";
    ss << this->supplementaryData;
    return ss.str();
}

std::string Tuple::getMultilineTupleString()
{
    std::stringstream ss;
    ss.precision(3);
    ss << std::fixed;
    ss << "Road Network: " << this->networkName << " (Nodes=" << this->nodes << ",Edges=" << this->edges << ")\n";;
    ss << "Supplementary Data: " << this->supplementaryData << "\n";
    return ss.str();
}

void Tuple::addSupplementaryFields(std::string desc, std::string data)
{
    if (this->supplementaryData == "none") {
        this->supplementaryData = "";
    }
    if (this->supplementaryData != "") {
        this->supplementaryData += ",";
    }
    this->supplementaryData += desc + "=" + data;
}

std::string Tuple::getAdditionalFields()
{
    std::stringstream ss;
    ss.precision(3);
    ss << std::fixed;
    for(std::size_t i = 0; i < this->additionalFields.size(); ++i) {
        ss << this->additionalFields[i];
        if (i != this->additionalFields.size()-1) {
            // If not last field we add a tab
            ss << "\t";
        }
    }
    return ss.str();
}

void Tuple::setAdditionalFields(std::vector<std::string>& fields)
{
    for(std::size_t i = 0; i < fields.size(); ++i) {
        this->additionalFields.push_back(fields[i]);
    }
}

