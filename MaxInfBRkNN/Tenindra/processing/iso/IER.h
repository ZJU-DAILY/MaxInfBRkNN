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

#ifndef _IER_H
#define _IER_H

#include "Graph.h"
#include "StaticRtree.h"
#include "Gtree.h"
#include "MortonList.h"
#include "pruned_highway_labeling.h"
#include "../common.h"
#include "../utility/Statistics.h"

#include <vector>

class IER {

    public:
        void getKNNsByDijkstra(StaticRtree& rtree, unsigned int k, NodeID queryNodeID, 
                               std::vector<NodeID>& kNNs, std::vector<EdgeWeight>& kNNDistances, 
                               Graph_RoadNetwork& graph);
        void getKNNsByDijkstraTravelTimes(StaticRtree& rtree, unsigned int k, NodeID queryNodeID, 
                                          std::vector<NodeID>& kNNs, std::vector<EdgeWeight>& kNNDistances, 
                                          Graph_RoadNetwork& graph);
        void getKNNsByGtree(StaticRtree& rtree, unsigned int k, NodeID queryNodeID, 
                               std::vector<NodeID>& kNNs, std::vector<EdgeWeight>& kNNDistances, 
                               Graph_RoadNetwork& graph, Gtree& gtree);
        void getKNNsBySILC(StaticRtree& rtree, unsigned int k, NodeID queryNodeID, 
                           std::vector<NodeID>& kNNs, std::vector<EdgeWeight>& kNNDistances, 
                           Graph_RoadNetwork& graph, SILCPathOracle& silc);
        void getKNNsByGtreeTravelTimes(StaticRtree& rtree, unsigned int k, NodeID queryNodeID, 
                               std::vector<NodeID>& kNNs, std::vector<EdgeWeight>& kNNDistances, 
                               Graph_RoadNetwork& graph, Gtree& gtree);
        void getKNNsByPHL(StaticRtree& rtree, unsigned int k, NodeID queryNodeID, 
                          std::vector<NodeID>& kNNs, std::vector<EdgeWeight>& kNNDistances, 
                          Graph_RoadNetwork& graph, PrunedHighwayLabeling& phl);
        void getKNNsByPHLTravelTimes(StaticRtree& rtree, unsigned int k, NodeID queryNodeID, 
                                     std::vector<NodeID>& kNNs, std::vector<EdgeWeight>& kNNDistances, 
                                     Graph_RoadNetwork& graph, PrunedHighwayLabeling& phl);
        Statistics stats;

};

#endif // _IER_H