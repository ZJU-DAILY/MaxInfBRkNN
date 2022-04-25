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

#include "DijkstraSearch.h"

#include "../queue/radix_heap.h"

// We assume priority has been cleared
void DijkstraSearch::findSSMTDistances(DynamicGraph& graph, NodeID source, 
                                       std::unordered_set<NodeID>& targetSet, 
                                       std::unordered_map<NodeID,EdgeWeight>& results,
                                       BinaryMinHeap<EdgeWeight,NodeID> *pqueue)
{
    std::vector<bool> isNodeSettled(graph.getNumNodes(),false);
//     std::unordered_set<NodeID> settledNodeIDs;
    
    EdgeWeight minDist;
    NodeID minDistNodeID, neighbourNodeID;
    
    // Initialize with priority queue with source node
    pqueue->insert(source,0);
    std::size_t targetsFound = 0;
    
    while (pqueue->size() > 0) {
        // Extract and remove node with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue->getMinKey();
        minDistNodeID = pqueue->extractMinElement();
        if (!isNodeSettled[minDistNodeID]) {
            isNodeSettled[minDistNodeID] = 1;
//         if (settledNodeIDs.find(minDistNodeID) == settledNodeIDs.end()) {
//             settledNodeIDs.insert(minDistNodeID); // Mark it as "seen" so we can avoid later

            if (targetSet.find(minDistNodeID) != targetSet.end()) {
                // We have found the SP distance to one of the target set
                // So we add it to result set and remove it from input set
                results[minDistNodeID] = minDist; // This is the final SP distance from source
                ++targetsFound;
                if (targetsFound == targetSet.size()) {
                    break;
                }
            }

            // Inspect each neighbour and update pqueue using edge weights
            for (std::size_t i = 0; i < graph.nodes[minDistNodeID].adjNodes.size(); ++i) {
                neighbourNodeID = graph.nodes[minDistNodeID].adjNodes[i];
                // Only update those we haven't already settled
                if (!isNodeSettled[neighbourNodeID]) {
//                 if (settledNodeIDs.find(neighbourNodeID) == settledNodeIDs.end()) {
                    pqueue->insert(neighbourNodeID,minDist + graph.nodes[minDistNodeID].adjNodeWgts[i]);
                }
            }
        }
    }

    //assert (targetSet.size() == targetsFound && "Dijkstra search could not find all targets in target set");
}

Path DijkstraSearch::findShortestPath(Graph_RoadNetwork& graph, NodeID source, NodeID target, std::vector<NodeID>& shortestPathTree)
{
    // We assume shortestPathTree has been resized for the size of the graph
    // However it does not need to have zero values as we overwrite them
    
    BinaryMinHeap<EdgeWeight,NodePair> pqueue;
    std::vector<bool> isNodeSettled(graph.getNumNodes(),false);
//     std::unordered_set<NodeID> settledNodeIDs;
    
    EdgeWeight minDist, distanceToTarget = 0;
    NodeID minDistNodeID, adjNode;
    int adjListStart, nextAdjListStart;
    
    // Initialize with priority queue with source node
    pqueue.insert(NodePair(source,constants::UNUSED_NODE_ID),0);
    
    while (pqueue.size() > 0) {
        // Extract and remove node with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue.getMinKey();
        NodePair minElement = pqueue.extractMinElement();
        minDistNodeID = minElement.first;
        if (!isNodeSettled[minDistNodeID]) {
            isNodeSettled[minDistNodeID] = true;
//         if (settledNodeIDs.find(minDistNodeID) == settledNodeIDs.end()) {
//             settledNodeIDs.insert(minDistNodeID); // Mark it as "seen" so we can avoid later
            shortestPathTree[minDistNodeID] = minElement.second;

            if (minDistNodeID == target) {
                // If the minimum is the target we have finished searching
                distanceToTarget = minDist;
                break;
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistNodeID);
            nextAdjListStart = graph.getEdgeListSize(minDistNodeID);

            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjNode = graph.edges[i].first;
                // Only update those we haven't already settled
                if (!isNodeSettled[adjNode]) {
//                 if (settledNodeIDs.find(adjNode) == settledNodeIDs.end()) {
                    pqueue.insert(NodePair(adjNode,minDistNodeID),minDist+graph.edges[i].second);
                }
            }
        }
    }
    
    // Retrieve and return shortest path
    Path path(source, true, distanceToTarget);
    for (NodeID currentNode = target; currentNode != source; currentNode = shortestPathTree[currentNode]) {
        path.addToBeginning(currentNode);
    }
    
    return path;

}

EdgeWeight DijkstraSearch::findShortestPathDistance(Graph_RoadNetwork& graph, NodeID source, NodeID target)
{
    BinaryMinHeap<EdgeWeight,NodeID> pqueue;
    
    std::vector<bool> isNodeSettled(graph.getNumNodes(),false);
//     std::unordered_set<NodeID> settledNodeIDs;
    
    EdgeWeight minDist, distanceToTarget = 0;
    NodeID minDistNodeID, adjNode;
    int adjListStart, nextAdjListStart;
    
    // Initialize with priority queue with source node
    pqueue.insert(source,0);
    
    while (pqueue.size() > 0) {
        // Extract and remove node with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue.getMinKey();
        minDistNodeID = pqueue.extractMinElement();
        if (!isNodeSettled[minDistNodeID]) {
            isNodeSettled[minDistNodeID] = true;
//         if (settledNodeIDs.find(minDistNodeID) == settledNodeIDs.end()) {
//             settledNodeIDs.insert(minDistNodeID); // Mark it as "seen" so we can avoid later

            if (minDistNodeID == target) {
                // If the minimum is the target we have finished searching
                distanceToTarget = minDist;
                break;
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistNodeID);
            nextAdjListStart = graph.getEdgeListSize(minDistNodeID);

            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjNode = graph.edges[i].first;
                // Only update those we haven't already settled
                if (!isNodeSettled[adjNode]) {
//                 if (settledNodeIDs.find(adjNode) == settledNodeIDs.end()) {
                    // This is the first time we have seen this node (i.e. distance infinity)
                    pqueue.insert(adjNode,minDist+graph.edges[i].second);
                }
            }
        }
    }

    return distanceToTarget;

}

// We assume priority has been cleared
void DijkstraSearch::findSSSPDistances(Graph_RoadNetwork& graph, NodeID source, std::vector<EdgeWeight>& targetDistances,
                                       MinPriorityQueueWithDK<EdgeWeight,NodeID> *pqueue)
{
    // We assume targetDistances has been resized for the size of the graph
    // However it does not need to have zero values as we overwrite them
    
    // Note: Using vector of bools is significantly faster than unordered_set.
    // Unlike the case where Dijkstra's search is between one source and
    // and one target, this is will never be less space efficient than unordered_set
    // because they will have the same number of elements. In fact std::vector<bool>
    // is more space efficient than a standard vector which is more space
    // and time efficient than a unordered_set (array vs. hash table)
    std::vector<bool> isNodeSettled(graph.getNumNodes(),false);
    int adjListStart, nextAdjListStart;
    
    EdgeWeight minDist, newDistance;
    NodeID minDistNodeID, adjNode;
    
    // Initialize with priority queue with source node
    pqueue->insert(source,0);
    
    while (pqueue->size() > 0) {
        // Extract and remove node with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue->getMinKey();
        minDistNodeID = pqueue->extractMinElement();
        isNodeSettled[minDistNodeID] = true;
        targetDistances[minDistNodeID] = minDist;

        // Inspect each neighbour and update pqueue using edge weights
        adjListStart = graph.getEdgeListStartIndex(minDistNodeID);
        nextAdjListStart = graph.getEdgeListSize(minDistNodeID);

        for (int i = adjListStart; i < nextAdjListStart; ++i) {
            adjNode = graph.edges[i].first;
            // Only update those we haven't already settled
            if (!isNodeSettled[adjNode]) {
                if (!pqueue->contains(adjNode)) {
                    // This is the first time we have seen this node (i.e. distance infinity)
                    pqueue->insert(adjNode,minDist+graph.edges[i].second);
                } else {
                    newDistance = minDist+graph.edges[i].second;
                    // Only decrease if a new distance is smaller than current distance
                    if (newDistance < pqueue->getKey(adjNode)) {
                        pqueue->decreaseKey(adjNode,newDistance);
                    }
                }
            }
        }
    }
}

// We assume priority has been cleared
void DijkstraSearch::findSSSPDistances(Graph_RoadNetwork& graph, NodeID source, std::vector<EdgeWeight>& targetDistances,
                                       MinPriorityQueue<EdgeWeight,NodeID> *pqueue)
{
    // We assume targetDistances has been resized for the size of the graph
    // However it does not need to have zero values as we overwrite them
    
    // Note: Using vector of bools is significantly faster than unordered_set.
    // Unlike the case where Dijkstra's search is between one source and
    // and one target, this is will never be less space efficient than unordered_set
    // because they will have the same number of elements. In fact std::vector<bool>
    // is more space efficient than a standard vector which is more space
    // and time efficient than a unordered_set (array vs. hash table)
    std::vector<bool> isNodeSettled(graph.getNumNodes(),false);
    int adjListStart, nextAdjListStart;
    
    EdgeWeight minDist;
    NodeID minDistNodeID, adjNode;
    
    // Initialize with priority queue with source node
    pqueue->insert(source,0);
    
    while (pqueue->size() > 0) {
        // Extract and remove node with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue->getMinKey();
        minDistNodeID = pqueue->extractMinElement();
        if (!isNodeSettled[minDistNodeID]) {
            isNodeSettled[minDistNodeID] = true;
            targetDistances[minDistNodeID] = minDist;

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistNodeID);
            nextAdjListStart = graph.getEdgeListSize(minDistNodeID);

            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjNode = graph.edges[i].first;
                // Only update those we haven't already settled
                if (!isNodeSettled[adjNode]) {
                    pqueue->insert(adjNode,minDist+graph.edges[i].second);
                }
            }
        }
    }
}

void DijkstraSearch::findSSSPDistances(Graph_RoadNetwork& graph, NodeID source, std::vector<EdgeWeight>& targetDistances)
{
    // We assume targetDistances has been resized for the size of the graph
    // However it does not need to have zero values as we overwrite them
    
    // Note: Using vector of bools is significantly faster than unordered_set.
    // Unlike the case where Dijkstra's search is between one source and
    // and one target, this is will never be less space efficient than unordered_set
    // because they will have the same number of elements. In fact std::vector<bool>
    // is more space efficient than a standard vector which is more space
    // and time efficient than a unordered_set (array vs. hash table)
    std::vector<bool> isNodeSettled(graph.getNumNodes(),false);
    int adjListStart, nextAdjListStart;
    
    radix_heap::pair_radix_heap<EdgeWeight,NodeID> pqueue;
    
    EdgeWeight minDist;
    NodeID minDistNodeID, adjNode;
    
    // Initialize with priority queue with source node
    pqueue.push(0,source);
    
    while (!pqueue.empty()) {
        // Extract and remove node with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDistNodeID = pqueue.top_value();
        minDist = pqueue.top_key();
        pqueue.pop();
        if (!isNodeSettled[minDistNodeID]) {
            isNodeSettled[minDistNodeID] = true;
            targetDistances[minDistNodeID] = minDist;

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistNodeID);
            nextAdjListStart = graph.getEdgeListSize(minDistNodeID);

            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjNode = graph.edges[i].first;
                // Only update those we haven't already settled
                if (!isNodeSettled[adjNode]) {
                    pqueue.push(minDist+graph.edges[i].second,adjNode);
                }
            }
        }
    }
}

// We assume priority has been cleared
void DijkstraSearch::findSSMTDistances(Graph_RoadNetwork& graph, NodeID source,
                                       std::unordered_set<NodeID>& targetSet, 
                                       std::unordered_map<NodeID,EdgeWeight>& results,
                                       BinaryMinHeap<EdgeWeight,NodeID> *pqueue)
{
    std::vector<bool> isNodeSettled(graph.getNumNodes(),false);
//     std::unordered_set<NodeID> settledNodeIDs;
    
    EdgeWeight minDist;
    NodeID minDistNodeID, adjNode;
    int adjListStart, nextAdjListStart;
    
    // Initialize with priority queue with source node
    pqueue->insert(source,0);
    std::size_t targetsFound = 0;
    
    while (pqueue->size() > 0) {
        // Extract and remove node with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue->getMinKey();
        minDistNodeID = pqueue->extractMinElement();
        if (!isNodeSettled[minDistNodeID]) {
            isNodeSettled[minDistNodeID] = true;
//         if (settledNodeIDs.find(minDistNodeID) == settledNodeIDs.end()) {
//             settledNodeIDs.insert(minDistNodeID); // Mark it as "seen" so we can avoid later

            if (targetSet.find(minDistNodeID) != targetSet.end()) {
                // We have found the SP distance to one of the target set
                // So we add it to result set and remove it from input set
                results[minDistNodeID] = minDist; // This is the final SP distance from source
                ++targetsFound;
                if (targetsFound == targetSet.size()) {
                    break;
                }
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistNodeID);
            nextAdjListStart = graph.getEdgeListSize(minDistNodeID);

            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjNode = graph.edges[i].first;
                // Only update those we haven't already settled
                if (!isNodeSettled[adjNode]) {
//                 if (settledNodeIDs.find(adjNode) == settledNodeIDs.end()) {
                    // This is the first time we have seen this node (i.e. distance infinity)
                    pqueue->insert(adjNode,minDist+graph.edges[i].second);
                }
            }
        }
    }

    //assert (targetSet.size() == targetsFound && "Dijkstra search could not find all targets in target set");
}

EdgeWeight DijkstraSearch::findShortestPathDistanceSubgraph(Graph_RoadNetwork& graph, NodeID source, NodeID target, std::vector<bool>& edgeInSubgraph)
{
    BinaryMinHeap<EdgeWeight,NodeID> pqueue;
//     std::vector<bool> isNodeSettled(graph.getNumNodes(),false);
    // We use unordered_set because we are only searching subgraph (no benefit 
    // from allocating std::vector<bool> for all nodes, just overhead)
    std::unordered_set<NodeID> settledNodeIDs;
    
    EdgeWeight minDist, distanceToTarget = 0;
    NodeID minDistNodeID, adjNode;
    int adjListStart, nextAdjListStart;
    
    // Initialize with priority queue with source node
    pqueue.insert(source,0);
    
    while (pqueue.size() > 0) {
        // Extract and remove node with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue.getMinKey();
        minDistNodeID = pqueue.extractMinElement();
//         if (!isNodeSettled[minDistNodeID]) {
//             isNodeSettled[minDistNodeID] = true;
        if (settledNodeIDs.find(minDistNodeID) == settledNodeIDs.end()) {
            settledNodeIDs.insert(minDistNodeID); // Mark it as "seen" so we can avoid later

            if (minDistNodeID == target) {
                // If the minimum is the target we have finished searching
                distanceToTarget = minDist;
                break;
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistNodeID);
            nextAdjListStart = graph.getEdgeListSize(minDistNodeID);

            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjNode = graph.edges[i].first;
                // Only update those we haven't already settled and within subgraph
//                 if (subgraph.find(adjNode) != subgraph.end() && !isNodeSettled[adjNode]) {
//                 if (subgraph.find(adjNode) != subgraph.end() && settledNodeIDs.find(adjNode) == settledNodeIDs.end()) {
                if (settledNodeIDs.find(adjNode) == settledNodeIDs.end() && edgeInSubgraph[i]) {
//                 if (!isNodeSettled[adjNode] && edgeInSubgraph[i]) {
                    // This is the first time we have seen this node (i.e. distance infinity)
                    pqueue.insert(adjNode,minDist+graph.edges[i].second);
                }
            }
        }
    }

    return distanceToTarget;
}

void DijkstraSearch::findSSMTDistancesSubgraph(Graph_RoadNetwork& graph, NodeID source,
                                               std::unordered_set<NodeID>& targetSet, 
                                               std::unordered_map<NodeID,EdgeWeight>& results,
                                               std::vector<bool>& edgeInSubgraph)
{
    BinaryMinHeap<EdgeWeight,NodeID> pqueue;
//     std::vector<bool> isNodeSettled(graph.getNumNodes(),false);
    // Note: We use unordered_set here because subgraph are small (they are limited by tau
    // in Gtree's case). So there is no gain.
    std::unordered_set<NodeID> settledNodeIDs;
    
    EdgeWeight minDist;
    NodeID minDistNodeID, adjNode;
    int adjListStart, nextAdjListStart;
    
    // Initialize with priority queue with source node
    pqueue.insert(source,0);
    std::size_t targetsFound = 0;
    
    while (pqueue.size() > 0) {
        // Extract and remove node with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue.getMinKey();
        minDistNodeID = pqueue.extractMinElement();
//         if (!isNodeSettled[minDistNodeID]) {
//             isNodeSettled[minDistNodeID] = 1;
        if (settledNodeIDs.find(minDistNodeID) == settledNodeIDs.end()) {
            settledNodeIDs.insert(minDistNodeID); // Mark it as "seen" so we can avoid later

            if (targetSet.find(minDistNodeID) != targetSet.end()) {
                // We have found the SP distance to one of the target set
                // So we add it to result set and remove it from input set
                results[minDistNodeID] = minDist; // This is the final SP distance from source
                ++targetsFound;
                if (targetsFound == targetSet.size()) {
                    break;
                }
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistNodeID);
            nextAdjListStart = graph.getEdgeListSize(minDistNodeID);
            
            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjNode = graph.edges[i].first;
                // Only update those we haven't already settled and within subgraph
                if (settledNodeIDs.find(adjNode) == settledNodeIDs.end() && edgeInSubgraph[i]) {
//                 if (!isNodeSettled[adjNode] && edgeInSubgraph[i]) {
//                 if (subgraph.find(adjNode) != subgraph.end() && !isNodeSettled[adjNode]) {
//                 if (subgraph.find(adjNode) != subgraph.end() && settledNodeIDs.find(adjNode) == settledNodeIDs.end()) {
                    // This is the first time we have seen this node (i.e. distance infinity)
                    pqueue.insert(adjNode,minDist+graph.edges[i].second);
                }
            }
        }
    }

    //assert (targetSet.size() == targetsFound && "Dijkstra search could not find all targets in target set");
}

// Note: This does not set a colour for the source, up to the caller to decide
// We assume priority has been cleared
bool DijkstraSearch::colourizeMap(Graph_RoadNetwork& graph, NodeID source, std::vector<EdgeID>& colourMap,
                                  std::vector<EdgeWeight>& distances, MinPriorityQueue<EdgeWeight,NodeLinkPair> *pqueue) {
    // We assume colourMap and distance have been resized for the size of the graph
    // However they do not need to have zero values as we overwrite all of them (except source)
    
    // Note: Using vector of bools is significantly faster than unordered_set.
    // Unlike the case where Dijkstra's search is between one source and
    // and one target, this is will never be less space efficient than unordered_set
    // because they will have the same number of elements. In fact std::vector<bool>
    // is more space efficient than a standard vector which is more space
    // and time efficient than a unordered_set (array vs. hash table)
    std::vector<bool> isNodeSettled(graph.getNumNodes(),false);
    int adjListStart, nextAdjListStart;
    
    EdgeWeight minDist;
    NodeID minDistNodeID, neighbourNodeID;
    
    adjListStart = graph.getEdgeListStartIndex(source);
    nextAdjListStart = graph.getEdgeListSize(source);
    int numColours = nextAdjListStart - adjListStart; // Max number of colours is number of outward edges from source
    
    // Initialize colours of the sources neighbours (as themselves)
    for (int i = adjListStart; i < nextAdjListStart; ++i) {
        colourMap[graph.edges[i].first] = static_cast<EdgeID>(i);
        pqueue->insert(NodeLinkPair(graph.edges[i].first,static_cast<EdgeID>(i-adjListStart)),graph.edges[i].second);
    }
    
    // Note: That it is possible that the shortest path to one of the adjacent nodes
    // of the source is through another adjacent node of the source, in which case
    // the colour for the original adjacent node will not be propagated to the rest 
    // of the graph because it's colour will overwritten as we will reach it through
    // the other adjacent node and update the colour as per usual
    
    while (pqueue->size() > 0) {
        // Extract and remove node with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue->getMinKey();
        NodeLinkPair minElement = pqueue->extractMinElement();
        minDistNodeID = minElement.first;
        if (!isNodeSettled[minDistNodeID]) {
            isNodeSettled[minDistNodeID] = true;
                        
            colourMap[minDistNodeID] = minElement.second;
            
            // Store shortest path distance from source
            distances[minDistNodeID] = minDist;

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistNodeID);
            nextAdjListStart = graph.getEdgeListSize(minDistNodeID);
            
            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                neighbourNodeID = graph.edges[i].first;
                // Only update those we haven't already settled
                if (!isNodeSettled[neighbourNodeID]) {
                    pqueue->insert(NodeLinkPair(neighbourNodeID,minElement.second),minDist+graph.edges[i].second);
                }
            }
        }
    }
    
    return (numColours == 1);
}

// Note: This does not set a colour for the source, up to the caller to decide
// We assume priority has been cleared
bool DijkstraSearch::colourizeMap(Graph_RoadNetwork& graph, NodeID source, std::vector<EdgeID>& colourMap,
                                  std::vector<EdgeWeight>& distances) {
//     BinaryMinHeap<EdgeWeight,NodeLinkPair> pqueue;
    radix_heap::pair_radix_heap<EdgeWeight,NodeLinkPair> pqueue;
    
    // We assume colourMap and distance have been resized for the size of the graph
    // However they do not need to have zero values as we overwrite all of them (except source)
    
    // Note: Using vector of bools is significantly faster than unordered_set.
    // Unlike the case where Dijkstra's search is between one source and
    // and one target, this is will never be less space efficient than unordered_set
    // because they will have the same number of elements. In fact std::vector<bool>
    // is more space efficient than a standard vector which is more space
    // and time efficient than a unordered_set (array vs. hash table)
    std::vector<bool> isNodeSettled(graph.getNumNodes(),false);
    int adjListStart, nextAdjListStart;
    
    EdgeWeight minDist;
    NodeID minDistNodeID, neighbourNodeID;
    
    adjListStart = graph.getEdgeListStartIndex(source);
    nextAdjListStart = graph.getEdgeListSize(source);
    int numColours = nextAdjListStart - adjListStart; // Max number of colours is number of outward edges from source
    
    // Initialize colours of the sources neighbours (as themselves)
    for (int i = adjListStart; i < nextAdjListStart; ++i) {
        colourMap[graph.edges[i].first] = static_cast<EdgeID>(i);
//         pqueue.insert(NodeLinkPair(graph.edges[i].first,static_cast<EdgeID>(i-adjListStart)),graph.edges[i].second);
        pqueue.push(graph.edges[i].second,NodeLinkPair(graph.edges[i].first,static_cast<EdgeID>(i-adjListStart)));
    }
    
    // Note: That it is possible that the shortest path to one of the adjacent nodes
    // of the source is through another adjacent node of the source, in which case
    // the colour for the original adjacent node will not be propagated to the rest 
    // of the graph because it's colour will overwritten as we will reach it through
    // the other adjacent node and update the colour as per usual
    
//     while (pqueue.size() > 0) {
    while (!pqueue.empty()) {
        // Extract and remove node with smallest distance from source
        // and mark it as "settled" so we do not inspect again
//         minDist = pqueue.getMinKey();
//         NodeLinkPair minElement = pqueue.extractMinElement();
//         minDistNodeID = minElement.first;
        minDist = pqueue.top_key();
        NodeLinkPair minElement = pqueue.top_value();
        minDistNodeID = minElement.first;
        pqueue.pop();
        if (!isNodeSettled[minDistNodeID]) {
            isNodeSettled[minDistNodeID] = true;

            colourMap[minDistNodeID] = minElement.second;
            
            // Store shortest path distance from source
            distances[minDistNodeID] = minDist;

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistNodeID);
            nextAdjListStart = graph.getEdgeListSize(minDistNodeID);
            
            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                neighbourNodeID = graph.edges[i].first;
                // Only update those we haven't already settled
                if (!isNodeSettled[neighbourNodeID]) {
//                     pqueue.insert(NodeLinkPair(neighbourNodeID,minElement.second),minDist+graph.edges[i].second);
                    pqueue.push(minDist+graph.edges[i].second,NodeLinkPair(neighbourNodeID,minElement.second));
                }
            }
        }
    }
    
    return (numColours == 1);
}