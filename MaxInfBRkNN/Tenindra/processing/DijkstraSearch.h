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

#ifndef _DIJKSTRASEARCH_H
#define _DIJkSTRASEARCH_H

#include "../common.h"
#include "DynamicGraph.h"
#include "Graph.h"
#include "Path.h"
#include "../queue/MinPriorityQueueWithDK.h"
#include "../queue/MinPriorityQueue.h"
#include "../queue/BinaryMinHeap.h"

#include <vector>
#include <unordered_set>
#include <unordered_map>

class DijkstraSearch {

    public:

        void findSSMTDistances(DynamicGraph& graph, NodeID source,
                               std::unordered_set<NodeID>& targetSet, 
                               std::unordered_map<NodeID,EdgeWeight>& results,
                               BinaryMinHeap<EdgeWeight,NodeID> *pqueue);

        Path findShortestPath(Graph_RoadNetwork& graph, NodeID source, NodeID target, std::vector<NodeID>& shortestPathTree);

        EdgeWeight findShortestPathDistance(Graph_RoadNetwork& graph, NodeID source, NodeID target);

        void findSSSPDistances(Graph_RoadNetwork& graph, NodeID source, std::vector<EdgeWeight>& targetDistances,
                               MinPriorityQueueWithDK<EdgeWeight,NodeID> *pqueue);

        void findSSSPDistances(Graph_RoadNetwork& graph, NodeID source, std::vector<EdgeWeight>& targetDistances,
                               MinPriorityQueue<EdgeWeight,NodeID> *pqueue);

        void findSSSPDistances(Graph_RoadNetwork& graph, NodeID source, std::vector<EdgeWeight>& targetDistances);

        void findSSMTDistances(Graph_RoadNetwork& graph, NodeID source,
                               std::unordered_set<NodeID>& targetSet, 
                               std::unordered_map<NodeID,EdgeWeight>& results,
                               BinaryMinHeap<EdgeWeight,NodeID> *pqueue);

        EdgeWeight findShortestPathDistanceSubgraph(Graph_RoadNetwork& graph, NodeID source, NodeID target, std::vector<bool>& edgeInSubgraph);
        
        void findSSMTDistancesSubgraph(Graph_RoadNetwork& graph, NodeID source,
                               std::unordered_set<NodeID>& targetSet, 
                               std::unordered_map<NodeID,EdgeWeight>& results,
                               std::vector<bool>& edgeInSubgraph);

        bool colourizeMap(Graph_RoadNetwork& graph, NodeID source, std::vector<EdgeID>& colourMap,
                          std::vector<EdgeWeight>& distances, MinPriorityQueue<EdgeWeight,NodeLinkPair> *pqueue);

        bool colourizeMap(Graph_RoadNetwork& graph, NodeID source, std::vector<EdgeID>& colourMap,
                          std::vector<EdgeWeight>& distances);

};

#endif // _DIJKSTRASEARCH_H