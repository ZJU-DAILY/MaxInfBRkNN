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

#include "SetGenerator.h"

#include "DijkstraSearch.h"
#include "StaticRtree.h"
#include "../queue/MinPriorityQueue.h"
#include "../queue/BinaryMinHeap.h"
#include "../utility/METISWrapper.h"

#include <random>
#include <unordered_set>
#include <deque>

// Status -1 indicates indicates ignore outdegrees, status 0 indicates exclude given degree, status 1 indicates include only given outdegree
std::vector<NodeID> SetGenerator::generateRandomSampleSet(Graph_RoadNetwork& graph, unsigned int numObjects, int status, int outdegree)
{
    std::vector<bool> testedNode(graph.getNumNodes(),false);
    std::vector<NodeID> sampleSet;
    
    int n = graph.getNumNodes()-1;

    std::random_device rd;
    std::default_random_engine e1(rd());
    std::uniform_int_distribution<int> uniform_dist(0, n);    

    NodeID randomNodeID;

    // Randomly select a nodes until object set is complete
    int numEdges;
    while (sampleSet.size() < numObjects) {
        randomNodeID = uniform_dist(e1);
        if (!testedNode[randomNodeID]) {
            if (status == -1) {
                testedNode[randomNodeID] = true;
                sampleSet.push_back(randomNodeID);
            } else if (status == 0) {
                numEdges = graph.getEdgeListSize(randomNodeID)-graph.getEdgeListStartIndex(randomNodeID);
                testedNode[randomNodeID] = true; // Whether we pick this node or not we don't need to test it again
                if (numEdges != outdegree) {
                    // We exclude the given outdegree (assumes there are enough nodes
                    // with other outdegrees to satisfy numObjects)
                    sampleSet.push_back(randomNodeID);
                }
            } else if (status == 1) {
                numEdges = graph.getEdgeListSize(randomNodeID)-graph.getEdgeListStartIndex(randomNodeID);
                testedNode[randomNodeID] = true; // Whether we pick this node or not we don't need to test it again
                if (numEdges == outdegree) {
                    // We only include the given outdegree (assumes there are enough nodes
                    // with outdegree to satisfy numObjects)
                    sampleSet.push_back(randomNodeID);
                }
            } else {
                std::cerr << "Invalid status provided" << std::endl;
                std::exit(1);
            }
        }
    }
 
    return sampleSet;
}

std::vector<NodeID> SetGenerator::generateRandomSampleSet(int numNodes, unsigned int numObjects)
{
    std::unordered_set<NodeID> nodesInSet;
    std::vector<NodeID> sampleSet;
    
    int n = numNodes-1;

    std::random_device rd;
    std::default_random_engine e1(rd());
    std::uniform_int_distribution<int> uniform_dist(0, n);    

    NodeID randomNodeID;

    // Randomly select a nodes until object set is complete
    while (sampleSet.size() < numObjects) {
        randomNodeID = uniform_dist(e1);
        if (nodesInSet.find(randomNodeID) == nodesInSet.end()) {
            nodesInSet.insert(randomNodeID);
            sampleSet.push_back(randomNodeID);
        }
    }
 
    return sampleSet;
}

void SetGenerator::generatePartitions(Graph_RoadNetwork& graph, int numPartitions, std::vector<std::unordered_set<NodeID>>& partitionsNodeSets,
                                      std::vector<std::vector<NodeID>>& partitionsNodes)
{
    std::unordered_set<NodeID> subgraph = graph.getNodesIDsUset();
    std::vector<NodeID> METISIdxToNodeID;
    METISIdxToNodeID.resize(graph.getNumNodes());
    int numNodes = graph.getNumNodes();

    METISWrapper metis(graph.getNumNodes(),graph.getNumEdges(),numPartitions);

    partitionsNodes.resize(numPartitions);
    partitionsNodeSets.resize(numPartitions);
    
    if (numPartitions > 1) {
        metis.populateMETISArrays(graph,subgraph,METISIdxToNodeID,true);
        metis.partitionSubgraph(numNodes);
        for (int i = 0; i < numNodes; ++i) {
            // Create child graphs using original NodeID not the METIS idx
            partitionsNodes[metis.parts[i]].push_back(METISIdxToNodeID[i]);
            partitionsNodeSets[metis.parts[i]].insert(METISIdxToNodeID[i]);
        }
    } else {
        for (int i = 0; i < numNodes; ++i) {
            partitionsNodes[0].push_back(i);
            partitionsNodeSets[0].insert(i);
        }
    }

    // Check partition worked correctly
    for (int i = 0; i < numPartitions; ++i) {
        //std::cout << "Subgraph " << i << " Size = " << childGraphs[i].size() << std::endl;
        if (partitionsNodes[i].size() == 0) {
            std::cerr << "METIS partitioning produced an empty subgraph, use a smaller number of partitions" << std::endl;
            std::exit(1);
        }
    }
}

std::vector<NodeID> SetGenerator::generateRandomClusteredSampleSet(Graph_RoadNetwork& graph, double clusterProbability, unsigned int numObjects, unsigned int clusterSize)
{
    std::unordered_set<NodeID> nodesInSet, visited;
    std::vector<NodeID> sampleSet;
    std::deque<NodeDistancePair> queue;
    
    int n = graph.getNumNodes()-1, adjListStart, nextAdjListStart, objectsAdded;

    std::random_device rd;
    std::default_random_engine e1(rd());
    std::uniform_int_distribution<int> uniform_dist(0, n);    
    std::uniform_int_distribution<int> cluster_dist(1,clusterSize); // Vary cluster size
    std::bernoulli_distribution pick(clusterProbability); // 5% chance that a visited verterx is added cluster

    NodeID randomNodeID, adjNode, nextNode, randomClusterSize;
    EdgeWeight hops;

//     std::vector<int> dist(clusterSize+1);
//     std::vector<int> numo(clusterSize+1);
    
    // Randomly select a nodes and grow cluster around it until object set is complete
    while (sampleSet.size() < numObjects) {
        randomNodeID = uniform_dist(e1);
        randomClusterSize = cluster_dist(e1);
        if (nodesInSet.find(randomNodeID) == nodesInSet.end()) {
            queue.clear();
            visited.clear();
            queue.push_back(NodeDistancePair(randomNodeID,0));
            visited.insert(randomNodeID);
            objectsAdded = 0;
            while (queue.size() > 0) {
                NodeDistancePair nextElement = queue.front();
                nextNode = nextElement.first;
                hops = nextElement.second;
                queue.pop_front();

                if (pick(e1)) {
                    nodesInSet.insert(nextNode);
                    sampleSet.push_back(nextNode);
                    ++objectsAdded;
                    if (objectsAdded == randomClusterSize) {
                        break;
                    }
                }

                adjListStart = graph.getEdgeListStartIndex(nextNode);
                nextAdjListStart = graph.getEdgeListSize(nextNode);

                for (int j = adjListStart; j < nextAdjListStart; ++j) {
                    adjNode = graph.edges[j].first;
                    // Only grow search with objects not already in sampleSet
                    if (nodesInSet.find(adjNode) == nodesInSet.end()) {
                        if (visited.find(adjNode) == visited.end()) {
                            queue.push_back(NodeDistancePair(adjNode,hops+1));
                            visited.insert(adjNode);
                        }
                    }
                }
            }
//             ++dist[randomClusterSize];
//             numo[randomClusterSize] += objectsAdded;
        }
    }
//     for (std::size_t i = 1; i < dist.size(); ++i) {
//         std::cout << "Size " << i << ": Clusters = " << dist[i] << ", Objects = " << numo[i] << std::endl;
//     }
//     std::cout << "Total Objects = " << sampleSet.size() << std::endl << std::endl;
 
    return sampleSet;
}

std::vector<NodeID> SetGenerator::generatePartitionSampleSet(Graph_RoadNetwork& graph, double hdProbability, double ldProbability, double highDensity,
                                                             double lowDensity, std::vector<std::unordered_set< NodeID>>& partitionsNodeSets, 
                                                             std::vector<std::vector<NodeID>>& partitionsNodes, int numObjects)
{
    int numPartitions = partitionsNodeSets.size();
    std::vector<bool> isPartitionProcessed(numPartitions,false);
    std::vector<NodeID> sampleSet;
    std::unordered_set<NodeID> nodesInSet;

    std::random_device rd;
    std::default_random_engine e1(rd());
    std::uniform_int_distribution<int> uniform_dist1(0,numPartitions-1);
    
    int numHDPartitions = std::min(static_cast<int>(std::ceil(numPartitions*hdProbability)),numPartitions);
    int numLDPartitions = std::min(static_cast<int>(std::ceil(numPartitions*ldProbability)),numPartitions-numHDPartitions);
    
    int n, randomNum, randomPartition, objectsRequired, objectsAdded, partitionsProcessed = 0;
    NodeID randomNode;
    while (partitionsProcessed < numHDPartitions) {
        randomPartition = uniform_dist1(e1);
        if (!isPartitionProcessed[randomPartition]) {
            isPartitionProcessed[randomPartition] = true;
            partitionsProcessed++;
            n = partitionsNodes[randomPartition].size();
            objectsRequired = std::ceil(highDensity*n);
            std::uniform_int_distribution<int> uniform_dist2(0, n-1);
            objectsAdded = 0;
            while (objectsAdded < objectsRequired) {
                randomNum = uniform_dist2(e1);
                auto rand_it = std::next(std::begin(partitionsNodes[randomPartition]), randomNum);
                randomNode = *rand_it;
                if (nodesInSet.find(randomNode) == nodesInSet.end()) {
                    nodesInSet.insert(randomNode);
                    sampleSet.push_back(randomNode);
                    ++objectsAdded;
                }
            }            
        }
    }
    
    partitionsProcessed = 0;
    while (partitionsProcessed < numLDPartitions) {
        randomPartition = uniform_dist1(e1);
        if (!isPartitionProcessed[randomPartition]) {
            isPartitionProcessed[randomPartition] = true;
            partitionsProcessed++;
            n = partitionsNodes[randomPartition].size();
            objectsRequired = std::ceil(lowDensity*n);
            std::uniform_int_distribution<int> uniform_dist2(0, n-1);
            objectsAdded = 0;
            while (objectsAdded < objectsRequired) {
                randomNum = uniform_dist2(e1);
                auto rand_it = std::next(std::begin(partitionsNodes[randomPartition]), randomNum);
                randomNode = *rand_it;
                if (nodesInSet.find(randomNode) == nodesInSet.end()) {
                    nodesInSet.insert(randomNode);
                    sampleSet.push_back(randomNode);
                    ++objectsAdded;
                }
            }            
        }
    }
    return sampleSet;   
}

std::vector<NodePair> SetGenerator::generateRandomPairSet(Graph_RoadNetwork& graph, unsigned int numPairs)
{
    std::vector<NodeID> nodes = graph.getNodesIDsVector();
    std::vector<NodePair> samplePairSet;

    int n = graph.getNumNodes();

    std::random_device rd;
    std::default_random_engine e1(rd());
    std::uniform_int_distribution<int> uniform_dist(0, n-1);    

    NodeID nodeID1, nodeID2, randomNum;

    // Randomly select two nodes at a time until set is complete
    randomNum = uniform_dist(e1);
    if (nodes.size() > 0) {
        auto rand_it = std::next(std::begin(nodes), randomNum);
        while (samplePairSet.size() < numPairs) {
            nodeID1 = *rand_it;
            randomNum = uniform_dist(e1);
            rand_it = std::next(std::begin(nodes), randomNum);
            nodeID2 = *rand_it;
            randomNum = uniform_dist(e1);
            rand_it = std::next(std::begin(nodes), randomNum);
            samplePairSet.push_back(std::make_pair(nodeID1,nodeID2));
        }
    }
    return samplePairSet;

}

std::vector<std::vector<NodeID>> SetGenerator::getNodeBuckets(Graph_RoadNetwork& graph, int numBuckets, std::string distanceType)
{
    // Find an approximate central node of the graph
    
    std::vector<NodeID> graphNodes = graph.getNodesIDsVector();
    
    // Create Rtree on the graph nodes and then find the closest
    // node by Euclidean to the Euclidean centre of the graph
    StaticRtree rtree(8,"",1,1,graph.getNumNodes());
    rtree.bulkLoad(graphNodes,graph.coordinates);
    
    Coordinate minX, maxX, minY, maxY, centreX, centreY;
    graph.getMinMaxCoordinates(minX,maxX,minY,maxY);
    centreX = (minX+maxX)/2;
    centreY = (minY+maxY)/2;
    
    std::vector<NodeID> kNNs;
    std::vector<EuclideanDistanceSqrd> kNNDistSqrd;
    rtree.getKNNs(1,centreX,centreY,kNNs,kNNDistSqrd);
    NodeID centralNode;
    if (kNNs.size() > 0) {
        centralNode = kNNs[0];
    } else {
        std::cerr << "Could not find a nearest neighbour to the central X and Y coordinate" << std::endl;
        std::exit(1);
    }
    
    // Compute shortest path distances to all other nodes from the central node
    DijkstraSearch dijk;
    std::vector<EdgeWeight> spDistances(graph.getNumNodes());
    MinPriorityQueue<EdgeWeight,NodeID>* pqueue = new BinaryMinHeap<EdgeWeight,NodeID>();
    dijk.findSSSPDistances(graph,centralNode,spDistances,pqueue);
    
    // Find the maximum distance from the central node to create buckets
    EdgeWeight maxSPDist = 0;
    for (NodeID n = 0; n < spDistances.size(); ++n) {
        if (spDistances[n] > maxSPDist) {
            maxSPDist = spDistances[n];
        }
    }
    
    // Now compute the upper distance for each bucket
    // and then sort all the nodes into the buckets
    std::vector<EdgeWeight> maxBucketDistances(numBuckets);
    if (distanceType == "powersof2") {
        for (int i = 1; i <= numBuckets; ++i) {
            maxBucketDistances[i-1] = maxSPDist*(1/std::pow(2,numBuckets-i));
            // Note: The final distance will equal maxSPDist
        }
    } else if (distanceType == "equidistant") {
        double separation = static_cast<double>(maxSPDist)/numBuckets;
        for (int i = 1; i <= numBuckets; ++i) {
            maxBucketDistances[i-1] = separation*i;
            // Note: The final distance will equal maxSPDist
        }
    }
    
    std::vector<std::vector<NodeID>> buckets(numBuckets);
    for (NodeID n = 0; n < spDistances.size(); ++n) {
        // Find the largest possible bucket it will fit into
        // start fromthe smallest (i == 0)
        for (std::size_t i = 0; i < maxBucketDistances.size(); ++i) {
            if (spDistances[n] < maxBucketDistances[i]) {
                buckets[i].push_back(n);
                break;
            }
        }
    }
    
    delete pqueue;
    
    return buckets;
}

std::vector<NodeID> SetGenerator::generateMinNetworkDistSampleSet(std::vector<std::vector<NodeID>>& nodeBuckets, int minBucket, int setSize)
{
    //assert(minBucket < nodeBuckets.size());
    return this->generateMinMaxNetworkDistSampleSet(nodeBuckets,minBucket,nodeBuckets.size()-1,setSize);
}

std::vector<NodeID> SetGenerator::generateMinMaxNetworkDistSampleSet(std::vector<std::vector<NodeID>>& nodeBuckets, int minBucket, int maxBucket, int setSize)
{
    //assert(minBucket < nodeBUckets.size() && maxBucket < nodeBUckets.size());
    //assert(minBucket <= maxBucket);

    // Note: This is equivalent to choosing the density implied by setSize from each bucket
    std::vector<NodeID> validNodes;
    for (int i = minBucket; i <= maxBucket; ++i) {
        // We add all nodes in minBucket and greater to one vector so they all have an
        // equal chance of being selected
        validNodes.insert(validNodes.end(),nodeBuckets[i].begin(),nodeBuckets[i].end());
    }
    
    std::unordered_set<NodeID> nodesInSet;
    std::vector<NodeID> sampleSet;
    
    int n = validNodes.size();

    std::random_device rd;
    std::default_random_engine e1(rd());
    std::uniform_int_distribution<int> uniform_dist(0, n-1);

    NodeID nodeID, randomNum;

    // Randomly select a nodes until object set is complete
    randomNum = uniform_dist(e1);
    if (validNodes.size() > setSize) {
        auto rand_it = std::next(std::begin(validNodes), randomNum);
        while (sampleSet.size() < setSize) {
            nodeID = *rand_it;
            // Set will be unique
            if (nodesInSet.find(nodeID) == nodesInSet.end()) {
                nodesInSet.insert(nodeID);
                sampleSet.push_back(nodeID);
            }
            randomNum = uniform_dist(e1);
            rand_it = std::next(std::begin(validNodes), randomNum);
        }
    } else {
        std::cerr << "There are not enough nodes between the bucket(s) " << minBucket << " and " << maxBucket << " to satisfy the set size" << std::endl;
        std::cout << "Nodes Available = " << validNodes.size() << std::endl;
        std::cout << "Required Set Size = " << setSize << std::endl;
        std::exit(1);
    }
    
    return sampleSet;
}

std::unordered_set<NodeID> SetGenerator::generateSubgraph(Graph_RoadNetwork& graph, Coordinate left, Coordinate bottom, Coordinate right, Coordinate top)
{
    //assert(left <= right && bottom <= top && "Euclidean window given is invalid");
    // Retrieve all nodes in the given range by an Rtree
    std::vector<NodeID> graphNodes = graph.getNodesIDsVector();
    StaticRtree rtree(8,"",1,1,graph.getNumNodes());
    rtree.bulkLoad(graphNodes,graph.coordinates);

    Rectangle window;
    window.setDimensions(left,bottom,right,top);
    std::vector<NodeID> windowNodes = rtree.windowQuery(window);
    std::unordered_set<NodeID> sampleSet(windowNodes.begin(),windowNodes.end());
    return sampleSet;
}

std::vector< NodeID > SetGenerator::generateSampleSubet(std::vector<NodeID> sampleSet, unsigned int numObjects)
{
    std::unordered_set<NodeID> nodesInSet;
    std::vector<NodeID> subset;
    
    int n = sampleSet.size()-1;

    std::random_device rd;
    std::default_random_engine e1(rd());
    std::uniform_int_distribution<int> uniform_dist(0, n);    

    // Randomly select a nodes until object set is complete
    NodeID randomNum, randomNodeID;
    while (subset.size() < numObjects) {
        randomNum = uniform_dist(e1);
        auto rand_it = std::next(std::begin(sampleSet), randomNum);
        randomNodeID = *rand_it;
        if (nodesInSet.find(randomNodeID) == nodesInSet.end()) {
            nodesInSet.insert(randomNodeID);
            subset.push_back(randomNodeID);
        }
    }
 
    return subset;
}

