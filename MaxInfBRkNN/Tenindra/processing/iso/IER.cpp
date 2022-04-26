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

#include "IER.h"
#include "time.h"
#include "DijkstraSearch.h"
#include "../queue/BinaryMaxHeap.h"
#include "iostream"
using namespace std;

void IER::getKNNsByDijkstra(StaticRtree& rtree, unsigned int k, NodeID queryNodeID, 
                            std::vector<NodeID>& kNNs, std::vector<EdgeWeight>& kNNDistances, Graph_RoadNetwork& graph)
{
    DijkstraSearch dijkstra;
    
    // Retrieve kNN by euclidean distance
    std::vector<NodeID> euclideanKNNs;
    std::vector<EuclideanDistanceSqrd> euclideanKNNDistances;
    Coordinate queryNodeX, queryNodeY;
    graph.getCoordinates(queryNodeID,queryNodeX,queryNodeY);
    BinaryMinHeap<EuclideanDistanceSqrd,RtreeDataTuple> heap = rtree.getKNNs(k,queryNodeX,queryNodeY,euclideanKNNs,euclideanKNNDistances);
    // Note: We keep heap so that we incrementally retrieve further euclidean NNs

    // We compute the network distances to each of these
    BinaryMaxHeap<EdgeWeight,NodeID> knnCandidates;
    EdgeWeight spDist, Dk, euclidDist;
    for (std::size_t i = 0; i < euclideanKNNs.size(); ++i) {
        spDist = dijkstra.findShortestPathDistance(graph,queryNodeID,euclideanKNNs[i]);
        knnCandidates.insert(euclideanKNNs[i],spDist);
    }

    // While the euclidean distance to the next euclidean nearest neigbour
    // is smaller than than network distance to the current kth neighbour
    // we can potentially find a closer nearest neighbour. Keep searching
    // until this lower bound exceed the kth neighbour network distance.
    Dk = knnCandidates.getMaxKey();
    NodeID nextEuclidNN;
    EuclideanDistanceSqrd currEuclidDistSqrd;
    while (true) {
        if (rtree.getNextNearestNeighbour(heap,queryNodeX,queryNodeY,nextEuclidNN,currEuclidDistSqrd)) {
            // Note: If there were less than k objects to begin with then getNextNearestNeighbour
            // would return false (as the heap is empty) and we not reach here
            euclidDist = std::floor(std::sqrt(currEuclidDistSqrd)); // Floor as it's a lowerbound
            if (euclidDist < Dk) {
                spDist = dijkstra.findShortestPathDistance(graph,queryNodeID,nextEuclidNN);
                if (spDist < Dk) {
                    // Only insert if it is a better kNN candidate
                    knnCandidates.insert(nextEuclidNN,spDist);
                    knnCandidates.extractMaxElement();
                    Dk = knnCandidates.getMaxKey();
                }
            } else {
                break;
            }
        } else {
            // This mean no nearest neighbours were found (we have reported
            // all objects) so we can stop the search
            break;
        }
    }

    knnCandidates.populateKNNs(kNNs,kNNDistances);
    
}

void IER::getKNNsByDijkstraTravelTimes(StaticRtree& rtree, unsigned int k, NodeID queryNodeID, std::vector<NodeID>& kNNs, 
                                       std::vector<EdgeWeight>& kNNDistances, Graph_RoadNetwork& graph)
{
    DijkstraSearch dijkstra;
    // Retrieve kNN by euclidean distance
    std::vector<NodeID> euclideanKNNs;
    std::vector<EuclideanDistanceSqrd> euclideanKNNDistances;
    Coordinate queryNodeX, queryNodeY;
    graph.getCoordinates(queryNodeID,queryNodeX,queryNodeY);
    BinaryMinHeap<EuclideanDistanceSqrd,RtreeDataTuple> heap = rtree.getKNNs(k,queryNodeX,queryNodeY,euclideanKNNs,euclideanKNNDistances);
    // Note: We keep heap so that we incrementally retrieve further euclidean NNs

    // We compute the network distances to each of these
    BinaryMaxHeap<EdgeWeight,NodeID> knnCandidates;
    EdgeWeight spDist, Dk, minTimeByEuclid;
    for (std::size_t i = 0; i < euclideanKNNs.size(); ++i) {
        spDist = dijkstra.findShortestPathDistance(graph,queryNodeID,euclideanKNNs[i]);
        knnCandidates.insert(euclideanKNNs[i],spDist);
    }

    // While the euclidean distance to the next euclidean nearest neigbour
    // is smaller than than network distance to the current kth neighbour
    // we can potentially find a closer nearest neighbour. Keep searching
    // until this lower bound exceed the kth neighbour network distance.
    Dk = knnCandidates.getMaxKey();
    NodeID nextEuclidNN;
    EuclideanDistanceSqrd currEuclidDistSqrd;
    while (true) {
        if (rtree.getNextNearestNeighbour(heap,queryNodeX,queryNodeY,nextEuclidNN,currEuclidDistSqrd)) {
            // Note: If there were less than k objects to begin with then getNextNearestNeighbour
            // would return false (as the heap is empty) and we not reach here
            minTimeByEuclid = std::ceil(std::sqrt(currEuclidDistSqrd)/graph.getMaxGraphSpeedByEdge()); // Floor as it's a lowerbound
            if (minTimeByEuclid < Dk) {
                spDist = dijkstra.findShortestPathDistance(graph,queryNodeID,nextEuclidNN);
                if (spDist < Dk) {
                    // Only insert if it is a better kNN candidate
                    knnCandidates.insert(nextEuclidNN,spDist);
                    knnCandidates.extractMaxElement();
                    Dk = knnCandidates.getMaxKey();
                }
            } else {
                break;
            }
        } else {
            // This mean no nearest neighbours were found (we have reported
            // all objects) so we can stop the search
            break;
        }
    }

    knnCandidates.populateKNNs(kNNs,kNNDistances);
}

void IER::getKNNsByGtree(StaticRtree& rtree, unsigned int k, NodeID queryNodeID, 
                         std::vector<NodeID>& kNNs, std::vector<EdgeWeight>& kNNDistances, 
                         Graph_RoadNetwork& graph, Gtree& gtree)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("false_hit",0);
    this->stats.initialiseStatistic("false_old_candidate",0);
    this->stats.initialiseStatistic("false_new_candidate",0);
    this->stats.initialiseStatistic("spdist_calls",0);
#endif
    
    // Retrieve kNN by euclidean distance
    std::vector<NodeID> euclideanKNNs;
    std::vector<EuclideanDistanceSqrd> euclideanKNNDistances;
    Coordinate queryNodeX, queryNodeY;
    graph.getCoordinates(queryNodeID,queryNodeX,queryNodeY);
    BinaryMinHeap<EuclideanDistanceSqrd,RtreeDataTuple> heap = rtree.getKNNs(k,queryNodeX,queryNodeY,euclideanKNNs,euclideanKNNDistances);
    // Note: We keep heap so that we incrementally retrieve further euclidean NNs

    // We compute the network distances to each of these
    BinaryMaxHeap<EdgeWeight,NodeID> knnCandidates;
    EdgeWeight spDist, Dk, euclidDist;
#if !defined(UNOPTIMISED_IER)
    std::vector<bool> visited(gtree.getTreeSize(),false);
    BinaryMinHeap<EdgeWeight,NodeID> leafQueue;
    std::vector<EdgeWeight> leafVerticesDistances;
    std::unordered_set<NodeID> leafVerticesVisited;
#endif
    for (std::size_t i = 0; i < euclideanKNNs.size(); ++i) {
#if defined(COLLECT_STATISTICS)
        this->stats.incrementStatistic("spdist_calls",1);
#endif
#if defined(UNOPTIMISED_IER)
        spDist = gtree.getShortestPathDistance(graph,queryNodeID,euclideanKNNs[i]);
#else
//         spDist = gtree.getRepeatedShortestPathDistance(graph,queryNodeID,euclideanKNNs[i],visited);
        spDist = gtree.getRepeatedShortestPathDistance(graph,queryNodeID,euclideanKNNs[i],visited,leafVerticesDistances,leafQueue,leafVerticesVisited);
#endif
#if defined(COLLECT_STATISTICS)
        this->stats.mergeStatistics(gtree.stats);
#endif
        knnCandidates.insert(euclideanKNNs[i],spDist);
    }

    // While the euclidean distance to the next euclidean nearest neigbour
    // is smaller than than network distance to the current kth neighbour
    // we can potentially find a closer nearest neighbour. Keep searching
    // until this lower bound exceed the kth neighbour network distance.
    Dk = knnCandidates.getMaxKey();
    NodeID nextEuclidNN;
    EuclideanDistanceSqrd currEuclidDistSqrd;
    
    while (true) {
        if (rtree.getNextNearestNeighbour(heap,queryNodeX,queryNodeY,nextEuclidNN,currEuclidDistSqrd)) {
            // Note: If there were less than k objects to begin with then getNextNearestNeighbour
            // would return false (as the heap is empty) and we not reach here
            euclidDist = std::floor(std::sqrt(currEuclidDistSqrd)); // Floor as it's a lowerbound
            if (euclidDist < Dk) {
#if defined(COLLECT_STATISTICS)
                // This means some candidate (existing or new) was a poor choice
                this->stats.incrementStatistic("false_hit",1);
                this->stats.incrementStatistic("spdist_calls",1);
#endif
#if defined(UNOPTIMISED_IER)
                spDist = gtree.getShortestPathDistance(graph,queryNodeID,nextEuclidNN);
#else
//                 spDist = gtree.getRepeatedShortestPathDistance(graph,queryNodeID,nextEuclidNN,visited);
                spDist = gtree.getRepeatedShortestPathDistance(graph,queryNodeID,nextEuclidNN,visited,leafVerticesDistances,leafQueue,leafVerticesVisited);
#endif
#if defined(COLLECT_STATISTICS)
                this->stats.mergeStatistics(gtree.stats);
#endif
                if (spDist < Dk) {
#if defined(COLLECT_STATISTICS)
                    // This means some old candidate was a poor choice
                    this->stats.incrementStatistic("false_old_candidate",1);
#endif
                    // Only insert if it is a better kNN candidate
                    knnCandidates.insert(nextEuclidNN,spDist);
                    knnCandidates.extractMaxElement();
                    Dk = knnCandidates.getMaxKey();
                }
#if defined(COLLECT_STATISTICS)
                // This means the new candidate was the bad choice
                  else {
                    this->stats.incrementStatistic("false_new_candidate",1);
                }
#endif
            } else {
                break;
            }
        } else {
            // This mean no nearest neighbours were found (we have reported
            // all objects) so we can stop the search
            break;
        }
    }

    knnCandidates.populateKNNs(kNNs,kNNDistances);
}

void IER::getKNNsByGtreeTravelTimes(StaticRtree& rtree, unsigned int k, NodeID queryNodeID, 
                         std::vector<NodeID>& kNNs, std::vector<EdgeWeight>& kNNDistances, 
                         Graph_RoadNetwork& graph, Gtree& gtree)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("false_hit",0);
    this->stats.initialiseStatistic("false_old_candidate",0);
    this->stats.initialiseStatistic("false_new_candidate",0);
    this->stats.initialiseStatistic("spdist_calls",0);
#endif
    
    // Retrieve kNN by euclidean distance
    std::vector<NodeID> euclideanKNNs;
    std::vector<EuclideanDistanceSqrd> euclideanKNNDistances;
    Coordinate queryNodeX, queryNodeY;
    graph.getCoordinates(queryNodeID,queryNodeX,queryNodeY);
    BinaryMinHeap<EuclideanDistanceSqrd,RtreeDataTuple> heap = rtree.getKNNs(k,queryNodeX,queryNodeY,euclideanKNNs,euclideanKNNDistances);
    // Note: We keep heap so that we incrementally retrieve further euclidean NNs

    // We compute the network distances to each of these
    BinaryMaxHeap<EdgeWeight,NodeID> knnCandidates;
    EdgeWeight spDist, Dk, minTimeByEuclid;
#if !defined(UNOPTIMISED_IER)
    std::vector<bool> visited(gtree.getTreeSize(),false);
    BinaryMinHeap<EdgeWeight,NodeID> leafQueue;
    std::vector<EdgeWeight> leafVerticesDistances;
    std::unordered_set<NodeID> leafVerticesVisited;
#endif
    for (std::size_t i = 0; i < euclideanKNNs.size(); ++i) {
#if defined(COLLECT_STATISTICS)
        this->stats.incrementStatistic("spdist_calls",1);
#endif
#if defined(UNOPTIMISED_IER)
        spDist = gtree.getShortestPathDistance(graph,queryNodeID,euclideanKNNs[i]);
#else
//         spDist = gtree.getRepeatedShortestPathDistance(graph,queryNodeID,euclideanKNNs[i],visited);
        spDist = gtree.getRepeatedShortestPathDistance(graph,queryNodeID,euclideanKNNs[i],visited,leafVerticesDistances,leafQueue,leafVerticesVisited);
#endif
#if defined(COLLECT_STATISTICS)
        this->stats.mergeStatistics(gtree.stats);
#endif
        knnCandidates.insert(euclideanKNNs[i],spDist);
    }

    // While the euclidean distance to the next euclidean nearest neigbour
    // is smaller than than network distance to the current kth neighbour
    // we can potentially find a closer nearest neighbour. Keep searching
    // until this lower bound exceed the kth neighbour network distance.
    Dk = knnCandidates.getMaxKey();
    NodeID nextEuclidNN;
    EuclideanDistanceSqrd currEuclidDistSqrd;
    while (true) {
        if (rtree.getNextNearestNeighbour(heap,queryNodeX,queryNodeY,nextEuclidNN,currEuclidDistSqrd)) {
            // Note: If there were less than k objects to begin with then getNextNearestNeighbour
            // would return false (as the heap is empty) and we not reach here
            minTimeByEuclid = std::ceil(std::sqrt(currEuclidDistSqrd)/graph.getMaxGraphSpeedByEdge()); // Floor as it's a lowerbound
            if (minTimeByEuclid < Dk) {
#if defined(COLLECT_STATISTICS)
                // This means some candidate (existing or new) was a poor choice
                this->stats.incrementStatistic("false_hit",1);
                this->stats.incrementStatistic("spdist_calls",1);
#endif
#if defined(UNOPTIMISED_IER)
                spDist = gtree.getShortestPathDistance(graph,queryNodeID,nextEuclidNN);
#else
//                 spDist = gtree.getRepeatedShortestPathDistance(graph,queryNodeID,nextEuclidNN,visited);
                spDist = gtree.getRepeatedShortestPathDistance(graph,queryNodeID,nextEuclidNN,visited,leafVerticesDistances,leafQueue,leafVerticesVisited);
#endif
#if defined(COLLECT_STATISTICS)
                this->stats.mergeStatistics(gtree.stats);
#endif
                if (spDist < Dk) {
#if defined(COLLECT_STATISTICS)
                    // This means some old candidate was a poor choice
                    this->stats.incrementStatistic("false_old_candidate",1);
#endif
                    // Only insert if it is a better kNN candidate
                    knnCandidates.insert(nextEuclidNN,spDist);
                    knnCandidates.extractMaxElement();
                    Dk = knnCandidates.getMaxKey();
                }
#if defined(COLLECT_STATISTICS)
                // This means the new candidate was the bad choice
                  else {
                    this->stats.incrementStatistic("false_new_candidate",1);
                }
#endif
            } else {
                break;
            }
        } else {
            // This mean no nearest neighbours were found (we have reported
            // all objects) so we can stop the search
            break;
        }
    }

    knnCandidates.populateKNNs(kNNs,kNNDistances);
}

void IER::getKNNsByPHL(StaticRtree& rtree, unsigned int k, NodeID queryNodeID, std::vector<NodeID>& kNNs, 
                       std::vector<EdgeWeight>& kNNDistances, Graph_RoadNetwork& graph, PrunedHighwayLabeling& phl)
{
    // Retrieve kNN by euclidean distance
    std::vector<NodeID> euclideanKNNs;
    std::vector<EuclideanDistanceSqrd> euclideanKNNDistances;
    Coordinate queryNodeX, queryNodeY;
    graph.getCoordinates(queryNodeID,queryNodeX,queryNodeY);
    BinaryMinHeap<EuclideanDistanceSqrd,RtreeDataTuple> heap = rtree.getKNNs(k,queryNodeX,queryNodeY,euclideanKNNs,euclideanKNNDistances);
    // Note: We keep heap so that we incrementally retrieve further euclidean NNs

    // We compute the network distances to each of these
    BinaryMaxHeap<EdgeWeight,NodeID> knnCandidates;
    EdgeWeight spDist, Dk, euclidDist;
    for (std::size_t i = 0; i < euclideanKNNs.size(); ++i) {
        NodeID _node = euclideanKNNs[i];
        clock_t startTime, endTime;
        startTime = clock();
        spDist = phl.Query(queryNodeID,_node);
        endTime = clock();
        std::cout<<"SPDist(v"<<queryNodeID<<",v"<<_node<<")="<<spDist<<endl;   //jins
        cout<<"runtime is:"<<(double )(endTime - startTime ) / CLOCKS_PER_SEC * 1000000.0 <<" us"<< endl;
        getchar();

        knnCandidates.insert(euclideanKNNs[i],spDist);
    }

    // While the euclidean distance to the next euclidean nearest neigbour
    // is smaller than than network distance to the current kth neighbour
    // we can potentially find a closer nearest neighbour. Keep searching
    // until this lower bound exceed the kth neighbour network distance.
    Dk = knnCandidates.getMaxKey();
    NodeID nextEuclidNN;
    EuclideanDistanceSqrd currEuclidDistSqrd;
    while (true) {
        if (rtree.getNextNearestNeighbour(heap,queryNodeX,queryNodeY,nextEuclidNN,currEuclidDistSqrd)) {
            // Note: If there were less than k objects to begin with then getNextNearestNeighbour
            // would return false (as the heap is empty) and we not reach here
            euclidDist = std::floor(std::sqrt(currEuclidDistSqrd)); // Floor as it's a lowerbound
            if (euclidDist < Dk) {
                spDist = phl.Query(queryNodeID,nextEuclidNN);
                if (spDist < Dk) {
                    // Only insert if it is a better kNN candidate
                    knnCandidates.insert(nextEuclidNN,spDist);
                    knnCandidates.extractMaxElement();
                    Dk = knnCandidates.getMaxKey();
                }
            } else {
                break;
            }
        } else {
            // This mean no nearest neighbours were found (we have reported
            // all objects) so we can stop the search
            break;
        }
    }

    knnCandidates.populateKNNs(kNNs,kNNDistances);
}


void IER::getKNNsByPHLTravelTimes(StaticRtree& rtree, unsigned int k, NodeID queryNodeID, std::vector<NodeID>& kNNs, 
                                  std::vector<EdgeWeight>& kNNDistances, Graph_RoadNetwork& graph, PrunedHighwayLabeling& phl)
{
    // Retrieve kNN by euclidean distance
    std::vector<NodeID> euclideanKNNs;
    std::vector<EuclideanDistanceSqrd> euclideanKNNDistances;
    Coordinate queryNodeX, queryNodeY;
    graph.getCoordinates(queryNodeID,queryNodeX,queryNodeY);
    BinaryMinHeap<EuclideanDistanceSqrd,RtreeDataTuple> heap = rtree.getKNNs(k,queryNodeX,queryNodeY,euclideanKNNs,euclideanKNNDistances);
    // Note: We keep heap so that we incrementally retrieve further euclidean NNs

    // We compute the network distances to each of these
    BinaryMaxHeap<EdgeWeight,NodeID> knnCandidates;
    EdgeWeight spDist, Dk, minTimeByEuclid;
    for (std::size_t i = 0; i < euclideanKNNs.size(); ++i) {
        spDist = phl.Query(queryNodeID,euclideanKNNs[i]);
        knnCandidates.insert(euclideanKNNs[i],spDist);
    }

    // While the euclidean distance to the next euclidean nearest neigbour
    // is smaller than than network distance to the current kth neighbour
    // we can potentially find a closer nearest neighbour. Keep searching
    // until this lower bound exceed the kth neighbour network distance.
    Dk = knnCandidates.getMaxKey();
    NodeID nextEuclidNN;
    EuclideanDistanceSqrd currEuclidDistSqrd;
    while (true) {
        if (rtree.getNextNearestNeighbour(heap,queryNodeX,queryNodeY,nextEuclidNN,currEuclidDistSqrd)) {
            // Note: If there were less than k objects to begin with then getNextNearestNeighbour
            // would return false (as the heap is empty) and we not reach here
            minTimeByEuclid = std::ceil(std::sqrt(currEuclidDistSqrd)/graph.getMaxGraphSpeedByEdge()); // Floor as it's a lowerbound
            if (minTimeByEuclid < Dk) {
                spDist = phl.Query(queryNodeID,nextEuclidNN);
                if (spDist < Dk) {
                    // Only insert if it is a better kNN candidate
                    knnCandidates.insert(nextEuclidNN,spDist);
                    knnCandidates.extractMaxElement();
                    Dk = knnCandidates.getMaxKey();
                }
            } else {
                break;
            }
        } else {
            // This mean no nearest neighbours were found (we have reported
            // all objects) so we can stop the search
            break;
        }
    }

    knnCandidates.populateKNNs(kNNs,kNNDistances);
}


void IER::getKNNsBySILC(StaticRtree& rtree, unsigned int k, NodeID queryNodeID, 
                        std::vector<NodeID>& kNNs, std::vector<EdgeWeight>& kNNDistances, 
                        Graph_RoadNetwork& graph, SILCPathOracle& silc)
{
    // Retrieve kNN by euclidean distance
    std::vector<NodeID> euclideanKNNs;
    std::vector<EuclideanDistanceSqrd> euclideanKNNDistances;
    Coordinate queryNodeX, queryNodeY;
    graph.getCoordinates(queryNodeID,queryNodeX,queryNodeY);
    BinaryMinHeap<EuclideanDistanceSqrd,RtreeDataTuple> heap = rtree.getKNNs(k,queryNodeX,queryNodeY,euclideanKNNs,euclideanKNNDistances);
    // Note: We keep heap so that we incrementally retrieve further euclidean NNs

    // We compute the network distances to each of these
    BinaryMaxHeap<EdgeWeight,NodeID> knnCandidates;
    EdgeWeight spDist, Dk, euclidDist;
    for (std::size_t i = 0; i < euclideanKNNs.size(); ++i) {
#if defined(UNOPTIMISED_IER)
        spDist = silc.findShortestPathDistance(graph,queryNodeID,euclideanKNNs[i]);
#else
        spDist = silc.findShortestPathDistanceOptimised(graph,queryNodeID,euclideanKNNs[i]);
#endif
        knnCandidates.insert(euclideanKNNs[i],spDist);
    }

    // While the euclidean distance to the next euclidean nearest neigbour
    // is smaller than than network distance to the current kth neighbour
    // we can potentially find a closer nearest neighbour. Keep searching
    // until this lower bound exceed the kth neighbour network distance.
    Dk = knnCandidates.getMaxKey();
    NodeID nextEuclidNN;
    EuclideanDistanceSqrd currEuclidDistSqrd;
    while (true) {
        if (rtree.getNextNearestNeighbour(heap,queryNodeX,queryNodeY,nextEuclidNN,currEuclidDistSqrd)) {
            // Note: If there were less than k objects to begin with then getNextNearestNeighbour
            // would return false (as the heap is empty) and we not reach here
            euclidDist = std::floor(std::sqrt(currEuclidDistSqrd)); // Floor as it's a lowerbound
            if (euclidDist < Dk) {
#if defined(UNOPTIMISED_IER)
                spDist = silc.findShortestPathDistance(graph,queryNodeID,nextEuclidNN);
#else
                spDist = silc.findShortestPathDistanceOptimised(graph,queryNodeID,nextEuclidNN);
#endif
                if (spDist < Dk) {
                    // Only insert if it is a better kNN candidate
                    knnCandidates.insert(nextEuclidNN,spDist);
                    knnCandidates.extractMaxElement();
                    Dk = knnCandidates.getMaxKey();
                }
            } else {
                break;
            }
        } else {
            // This mean no nearest neighbours were found (we have reported
            // all objects) so we can stop the search
            break;
        }
    }

    knnCandidates.populateKNNs(kNNs,kNNDistances);
}