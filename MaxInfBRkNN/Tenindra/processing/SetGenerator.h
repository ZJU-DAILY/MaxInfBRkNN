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

#ifndef _SETGENERATOR_H
#define _SETGENERATOR_H

#include "Graph.h"

#include <vector>

class SetGenerator {

    public:
        std::vector<NodeID> generateRandomSampleSet(Graph_RoadNetwork& graph, unsigned int numObjects, int status = -1, int outdegree = 0);
        std::vector<NodeID> generateRandomSampleSet(int numNodes, unsigned int numObjects);
        std::vector<NodeID> generateRandomClusteredSampleSet(Graph_RoadNetwork& graph, double clusterProbability, unsigned int numObjects, unsigned int clusterSize);
        void generatePartitions(Graph_RoadNetwork& graph, int numPartitions, std::vector<std::unordered_set<NodeID>>& partitionsNodeSets,
                                std::vector<std::vector<NodeID>>& partitionsNodes);
        std::vector<NodeID> generatePartitionSampleSet(Graph_RoadNetwork& graph, double hdProbability, double ldProbability, double highDensity, double lowDensity, std::vector<std::unordered_set<NodeID>>& partitionsNodeSets, std::vector<std::vector<NodeID>>& partitionsNodes, int numObjects = 0);
//         std::vector<NodeID> generateSampleSetByPartitions(Graph_RoadNetwork& graph, int clusterSize, std::vector<std::unordered_set<NodeID>>& partitionsNodeSets, std::vector<std::vector<NodeID>>& partitionsNodes);
//         std::vector<NodeID> generatePartitionBasedSampleSet(Graph_RoadNetwork& graph, int numPartitions, int clusterSize);
        std::vector<NodePair> generateRandomPairSet(Graph_RoadNetwork& graph, unsigned int numPairs);
        // Choose an approximate "centre" of the graph then divided the graph into buckets 
        // based on the distance from the centre according to a exponential formula
        // We can chooose how bucket distance are calculated (either powersof2 or equidistant)
        std::vector<std::vector<NodeID>> getNodeBuckets(Graph_RoadNetwork& graph, int numBuckets, std::string distanceType);
        // Given sets of graph node buckets, select setSize number of nodes from all buckets
        // larger than the minBucket (using the above function this guarantees all nodes
        // are at least some number of buckets away from the centre
        std::vector<NodeID> generateMinNetworkDistSampleSet(std::vector<std::vector<NodeID>>& nodeBuckets, int minBucket, int setSize);
        // Given sets of graph node buckets, select setSize number of nodes from all buckets
        // larger than the minBucket but less than max bucket (using the above function 
        // this guarantees all nodes are within some range of buckets. If the min and max 
        // bucket are equal, we select nodes exclusively from one bucket
        std::vector<NodeID> generateMinMaxNetworkDistSampleSet(std::vector<std::vector<NodeID>>& nodeBuckets, int minBucket, int maxBucket, int setSize);
        
        // Generate a random connected subgraph of size n
        std::unordered_set<NodeID> generateSubgraph(Graph_RoadNetwork& graph, Coordinate left, Coordinate bottom, Coordinate right, Coordinate top);
        std::vector<NodeID> generateSampleSubet(std::vector<NodeID> sampleSet, unsigned int numObjects);

};

#endif // _SETGENERATOR_H