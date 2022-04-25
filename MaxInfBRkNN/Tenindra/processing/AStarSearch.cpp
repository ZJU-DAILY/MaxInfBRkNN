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

#include "AStarSearch.h"

#include "../queue/BinaryMinHeap.h"

Path AStarSearch::findShortestPath(Graph_RoadNetwork& graph, NodeID source, NodeID target,
                                   std::vector<NodeID>& shortestPathTree)
{
    // We assume shortestPathTree has been resized for the size of the graph
    // However it does not need to have zero values as we overwrite them
    
    BinaryMinHeap<EdgeWeight,AStarHeapElement> pqueue;
    std::vector<bool> isNodeSettled(graph.getNumNodes(),false);
    
//     EdgeWeight adjNodeWgt;
    EdgeWeight  minDist, sourceToAdjNodeDist, distanceToTarget = 0;
    NodeID minDistNodeID, adjNode;
    int adjListStart, adjListSize;
    
    // Initialize priority queue with source node
    EdgeWeight minSourceTargetDist = static_cast<EdgeWeight>(graph.getEuclideanDistance(source,target));
    pqueue.insert(AStarHeapElement(source,constants::UNUSED_NODE_ID,0),minSourceTargetDist);
    
    while (pqueue.size() > 0) {
        // Extract and remove node with smallest possible distance to target
        // and mark it as "settled" so we do not inspect again
        // Note: In A* search this willl still lead to a optimal solution
        // if the heuristic function we use is consistent (i.e. never overestimates)
        AStarHeapElement minElement = pqueue.extractMinElement();
        minDistNodeID = minElement.node;
        if (!isNodeSettled[minDistNodeID]) {
            isNodeSettled[minDistNodeID] = true; // Mark it as "settled" so we can avoid later
            minDist = minElement.sourceNodeDist;
            shortestPathTree[minDistNodeID] = minElement.predecessor;
            
            if (minDistNodeID == target) {
                // If the minimum is the target we have finished searching
                distanceToTarget = minDist;
                break;
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistNodeID);
            adjListSize = graph.getEdgeListSize(minDistNodeID);
            
            for (int i = adjListStart; i < adjListSize; ++i) {
                adjNode = graph.edges[i].first;
                // Only consider neighbours we haven't already settled
                if (!isNodeSettled[adjNode]) {
//                     // Heuristic Consistency Test Output
//                     EdgeWeight currentToTargetEst = graph.getEuclideanDistance(minDistNodeID,target);
//                     EdgeWeight neighbourToTargetEst = graph.getEuclideanDistance(adjNode,target);
//                     if (currentToTargetEst > adjNodeWgt+neighbourToTargetEst) {
//                         std::cout << "currentToTargetEst = " << currentToTargetEst << std::endl;
//                         std::cout << "adjNodeWgt = " << adjNodeWgt << std::endl;
//                         std::cout << "neighbourToTargetEst = " << neighbourToTargetEst << std::endl;
//                         std::cout << "Diff = " << currentToTargetEst - adjNodeWgt - neighbourToTargetEst << std::endl;
//                     }
                    //assert (currentToTargetEst <= adjNodeWgt+neighbourToTargetEst && "Heuristic function is not consistent");
                    
                    sourceToAdjNodeDist = minDist + graph.edges[i].second;
                    minSourceTargetDist = sourceToAdjNodeDist + static_cast<EdgeWeight>(graph.getEuclideanDistance(adjNode,target));
                    pqueue.insert(AStarHeapElement(adjNode,minDistNodeID,sourceToAdjNodeDist),minSourceTargetDist);
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

EdgeWeight AStarSearch::findShortestPathDistance(Graph_RoadNetwork& graph, NodeID source, NodeID target)
{
    BinaryMinHeap<EdgeWeight,NodeDistancePair> pqueue;
    std::vector<bool> isNodeSettled(graph.getNumNodes(),false);
    
//     EdgeWeight adjNodeWgt;
    EdgeWeight  minDist, sourceToAdjNodeDist, distanceToTarget = 0;
    NodeID minDistNodeID, adjNode;
    int adjListStart, adjListSize;
    
    // Initialize priority queue with source node
    EdgeWeight minSourceTargetDist = static_cast<EdgeWeight>(graph.getEuclideanDistance(source,target));
    pqueue.insert(NodeDistancePair(source,0),minSourceTargetDist);
    
    while (pqueue.size() > 0) {
        // Extract and remove node with smallest possible distance to target
        // and mark it as "settled" so we do not inspect again
        // Note: In A* search this willl still lead to a optimal solution
        // if the heuristic function we use is consistent (i.e. never overestimates)
        NodeDistancePair minElement = pqueue.extractMinElement();
        minDistNodeID = minElement.first;
        if (!isNodeSettled[minDistNodeID]) {
            isNodeSettled[minDistNodeID] = true; // Mark it as "settled" so we can avoid later
            minDist = minElement.second;
            
            if (minDistNodeID == target) {
                // If the minimum is the target we have finished searching
                distanceToTarget = minDist;
                break;
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistNodeID);
            adjListSize = graph.getEdgeListSize(minDistNodeID);
            
            for (int i = adjListStart; i < adjListSize; ++i) {
                adjNode = graph.edges[i].first;
                // Only consider neighbours we haven't already settled
                if (!isNodeSettled[adjNode]) {
//                     // Heuristic Consistency Test Output
//                     EdgeWeight currentToTargetEst = graph.getEuclideanDistance(minDistNodeID,target);
//                     EdgeWeight neighbourToTargetEst = graph.getEuclideanDistance(adjNode,target);
//                     if (currentToTargetEst > adjNodeWgt+neighbourToTargetEst) {
//                         std::cout << "currentToTargetEst = " << currentToTargetEst << std::endl;
//                         std::cout << "adjNodeWgt = " << adjNodeWgt << std::endl;
//                         std::cout << "neighbourToTargetEst = " << neighbourToTargetEst << std::endl;
//                         std::cout << "Diff = " << currentToTargetEst - adjNodeWgt - neighbourToTargetEst << std::endl;
//                     }
                    //assert (currentToTargetEst <= adjNodeWgt+neighbourToTargetEst && "Heuristic function is not consistent");
                    
                    sourceToAdjNodeDist = minDist + graph.edges[i].second;
                    minSourceTargetDist = sourceToAdjNodeDist + static_cast<EdgeWeight>(graph.getEuclideanDistance(adjNode,target));
                    pqueue.insert(NodeDistancePair(adjNode,sourceToAdjNodeDist),minSourceTargetDist);
                }
            }
        }
    }

    return distanceToTarget;
}
