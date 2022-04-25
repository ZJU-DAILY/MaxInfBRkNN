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

#ifndef _PATH_H
#define _PATH_H

#include "../common.h"

#include <deque>

class Path {

    public:
        Path(NodeID sourceID, bool isShortest, EdgeWeight pathLength = 0);
        EdgeWeight getLength();
        bool isShorestPath();
        std::string getPathString();
        void addToBeginning(NodeID nodeID);
        void addToEnd(NodeID nodeID);
        void addToBeginning(NodeID nodeID, EdgeWeight weightToNode);
        void addToEnd(NodeID nodeID, EdgeWeight weightToNode);
        int getNumLinks();
        void setPathLength(EdgeWeight pathLength);
        
    private:
        NodeID sourceID;
        bool isShortest;
        std::deque<NodeID> pathNodes; // source is not included
        EdgeWeight pathLength = 0;
};

#endif // _PATH_H