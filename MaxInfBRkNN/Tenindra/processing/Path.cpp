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

#include "Path.h"

#include <sstream>

Path::Path(NodeID sourceID, bool isShortest, EdgeWeight pathLength): sourceID(sourceID), isShortest(isShortest), pathLength(pathLength) {}

void Path::addToBeginning(NodeID nodeID)
{
    this->pathNodes.push_front(nodeID);
}

void Path::addToEnd(NodeID nodeID) {
    this->pathNodes.push_back(nodeID);
}

void Path::addToBeginning(NodeID nodeID, EdgeWeight weightToNode)
{
    this->pathNodes.push_front(nodeID);
}

void Path::addToEnd(NodeID nodeID, EdgeWeight weightToNode) {
    this->pathNodes.push_back(nodeID);
}

EdgeWeight Path::getLength()
{
    return this->pathLength;
}

bool Path::isShorestPath()
{
    return this->isShortest;
}

int Path::getNumLinks()
{
    return pathNodes.size();
}

std::string Path::getPathString()
{
    std::stringstream ss;
    ss << this->sourceID;
    for (std::size_t i = 0; i < this->pathNodes.size(); ++i) {
        ss << " -> " << this->pathNodes[i];
    }
    return ss.str();
}

void Path::setPathLength(EdgeWeight pathLength) {
    this->pathLength = pathLength;
}

