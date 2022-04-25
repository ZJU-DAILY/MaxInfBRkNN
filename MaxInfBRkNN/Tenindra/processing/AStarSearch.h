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

#ifndef _ASTARSEARCH_H
#define _ASTARSEARCH_H

#include "../common.h"
#include "Graph.h"
#include "Path.h"

#include <vector>

struct AStarHeapElement {
    NodeID node;
    NodeID predecessor;
    EdgeWeight sourceNodeDist;
    AStarHeapElement(NodeID node, NodeID predecessor, EdgeWeight sourceNodeDist): 
        node(node), predecessor(predecessor), sourceNodeDist(sourceNodeDist) {}
};

class AStarSearch {

    public:
        Path findShortestPath(Graph_RoadNetwork& graph, NodeID source, NodeID target, std::vector<NodeID>& shortestPathTree);

        PathDistance findShortestPathDistance(Graph_RoadNetwork& graph, NodeID source, NodeID target);

};

#endif // _ASTARSEARCH_H