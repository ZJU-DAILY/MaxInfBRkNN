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

#include "INE.h"

#include "../queue/BinaryMinHeapWithDK.h"
#include "../queue/BinaryMinHeap.h"

void INE::getKNNsByDynamicGraph(DynamicGraph& graph, unsigned int k, NodeID queryNodeID, std::vector<NodeID>& kNNs, std::vector<EdgeWeight>& kNNDistances)
{
    // Note: We assumpe vectors passed to this method are empty
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("nodes_settled",0);
    this->stats.initialiseStatistic("edges_relaxed",0);
#endif

#if defined(INE_QUEUE_OPTIMISED)
    BinaryMinHeap<EdgeWeight,NodeID> pqueue;
#else
    BinaryMinHeapWithDK<EdgeWeight,NodeID> pqueue;
#endif
    EdgeWeight minDist, newDistance;
    NodeID minDistNodeID, adjNode;
#if defined(INE_VISITED_OPTIMISED)
    std::vector<bool> isNodeSettled(graph.getNumNodes(),false);
#else
    std::unordered_set<NodeID> settledNodeIDs;
#endif

    // Initialize with priority queue with query node ID
    pqueue.insert(queryNodeID,0);
    
    while (pqueue.size() > 0) {
        // Extract and remove node with smallest distance from query point
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue.getMinKey();
        minDistNodeID = pqueue.extractMinElement();
#if defined(INE_QUEUE_OPTIMISED)
    // Note: If heap is not optimised then algorithm uses decrease key 
    // so the queue will never have settled nodes in it
    #if defined(INE_VISITED_OPTIMISED)
        if (!isNodeSettled[minDistNodeID]) {
            isNodeSettled[minDistNodeID] = 1;
    #else
        if (settledNodeIDs.find(minDistNodeID) == settledNodeIDs.end()) {
            settledNodeIDs.insert(minDistNodeID);
    #endif
#else 
    #if defined(INE_VISITED_OPTIMISED)
        isNodeSettled[minDistNodeID] = 1;
    #else
        settledNodeIDs.insert(minDistNodeID);
    #endif
#endif
#if defined(COLLECT_STATISTICS)
            this->stats.incrementStatistic("nodes_settled",1);
#endif

            if (graph.isObject(minDistNodeID)) {
                // If the minimum is an object we have found a kNN
                kNNs.push_back(minDistNodeID);
                kNNDistances.push_back(minDist);
                if (kNNs.size() == k) {
                    // If this is the kth nearest neighbour object
                    // then we no longer need to inspect neighbours
                    break;
                }
            }

            // Inspect each neighbour and update pqueue using edge weights
            for (std::size_t i = 0; i < graph.nodes[minDistNodeID].adjNodes.size(); ++i) {
                adjNode = graph.nodes[minDistNodeID].adjNodes[i];
#if defined(INE_VISITED_OPTIMISED)
                if (!isNodeSettled[adjNode]) {
#else
                if (settledNodeIDs.find(adjNode) == settledNodeIDs.end()) {
#endif
                    // Only update those we haven't already settled
                    newDistance = minDist+graph.nodes[minDistNodeID].adjNodeWgts[i];
#if !defined(INE_QUEUE_OPTIMISED)
// Note: If heap is not optimised then decrease key is supported
                    if (!pqueue.contains(adjNode)) {
                        // This is the first time we have seen this node (i.e. distance infinity)
#endif
                        pqueue.insert(adjNode,newDistance);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("edges_relaxed",1);
#endif
#if !defined(INE_QUEUE_OPTIMISED)
                    } else {
                        // Only decrease if a new distance is smaller than current distance
                        if (newDistance < pqueue.getKey(adjNode)) {
                            pqueue.decreaseKey(adjNode,newDistance);
                        }
                    }
#endif
                }
            }
#if defined(INE_QUEUE_OPTIMISED)
        }
#endif
    }
    
}

void INE::getKNNs(Graph_RoadNetwork& graph, unsigned int k, NodeID queryNodeID, std::vector<NodeID>& kNNs, std::vector<EdgeWeight>& kNNDistances)
{
    // Note: We assumpe vectors passed to this method are empty
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("nodes_settled",0);
    this->stats.initialiseStatistic("edges_relaxed",0);
#endif

    BinaryMinHeap<EdgeWeight,NodeID> pqueue;
    EdgeWeight minDist, newDistance;
    NodeID minDistNodeID, adjNode;
    std::vector<bool> isNodeSettled(graph.getNumNodes(),false);
//     std::unordered_set<NodeID> settledNodeIDs;
    int adjListStart, nextAdjListStart;
    
    // Initialize with priority queue with query node ID
    pqueue.insert(queryNodeID,0);
    
    while (pqueue.size() > 0) {
        // Extract and remove node with smallest distance from query point
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue.getMinKey();
        minDistNodeID = pqueue.extractMinElement();
        if (!isNodeSettled[minDistNodeID]) {
            isNodeSettled[minDistNodeID] = 1;
//         if (settledNodeIDs.find(minDistNodeID) == settledNodeIDs.end()) {
//             settledNodeIDs.insert(minDistNodeID);
#if defined(COLLECT_STATISTICS)
            this->stats.incrementStatistic("nodes_settled",1);
#endif

            if (graph.isObject(minDistNodeID)) {
                // If the minimum is an object we have found a kNN
                kNNs.push_back(minDistNodeID);
                kNNDistances.push_back(minDist);
                if (kNNs.size() == k) {
                    // If this is the kth nearest neighbour object
                    // then we no longer need to inspect neighbours
                    break;
                }
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistNodeID);
            nextAdjListStart = graph.getEdgeListSize(minDistNodeID);
            // Note: We have to make sure we don't exceed size of graph.edges vector
            
            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjNode = graph.edges[i].first;
                if (!isNodeSettled[adjNode]) {
//                 if (settledNodeIDs.find(adjNode) == settledNodeIDs.end()) {
                    // Only update those we haven't already settled
                    newDistance = minDist + graph.edges[i].second;
                    pqueue.insert(adjNode,newDistance);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("edges_relaxed",1);
#endif
                }
            }
        }
    }
    
}

