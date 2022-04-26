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

#include "ROAD.h"

#include "DijkstraSearch.h"
#include "../utility/utility.h"

#include <deque>

ROAD::ROAD(std::string networkName, int numNodes, int numEdges, int fanout, int levels): 
    networkName(networkName), numNodes(numNodes), numEdges(numEdges), fanout(fanout), levels(levels) {
    this->routeOverlay.reserve(numNodes);
    this->borders.resize(numNodes,false);
}
    
std::string ROAD::getNetworkName()
{
    return this->networkName;
}

int ROAD::getNumEdges()
{
    return this->numEdges;
}

int ROAD::getNumNodes()
{
    return this->numNodes;
}

int ROAD::getFanout()
{
    return this->fanout;
}

int ROAD::getNumLevels()
{
     return this->levels;
}

int ROAD::getNumBorders()
{
    int numBorders = 0;
    for (std::size_t i = 0; i < this->rnetTree.size(); ++i) {
        numBorders += this->rnetTree[i].getNumBorders();
    }
    return numBorders;
}

int ROAD::getBorderToBorderRelationships()
{
    int numB2BRelationships = 0;
    for (std::size_t i = 0; i < this->rnetTree.size(); ++i) {
        numB2BRelationships += this->rnetTree[i].getNumBorders()*this->rnetTree[i].getNumBorders();
    }
    return numB2BRelationships;
}

// Returns the number of levels in graph as opposed to what we told it to build
int ROAD::getRealNumLevels()
{
    std::vector<std::vector<int>> treeLevelIdxs = this->getRnetIdxsByLevel();
    return treeLevelIdxs.size();
}

void ROAD::printLevels()
{
    std::vector<std::vector<int>> treeLevelIdxs = this->getRnetIdxsByLevel();
    
    for (int i = treeLevelIdxs.size()-1; i >= 0; --i) {
        std::cout << "Level " << i << ": " << treeLevelIdxs[i].size() << std::endl;
    }
}

std::vector< std::vector<int>> ROAD::getRnetIdxsByLevel()
{
    std::vector<std::vector<int>> treeLevelIdxs;
    std::vector<int> currentLevel, nextLevel;
    nextLevel.push_back(0);

    while(nextLevel.size() != 0){
        treeLevelIdxs.push_back(nextLevel);
        currentLevel.swap(nextLevel);
        nextLevel.clear();
        for (int i: currentLevel){
            for (std::size_t j = 0; j < this->rnetTree[i].children.size(); j++ ){
                nextLevel.push_back(this->rnetTree[i].children[j]);
            }
        }
    }
    return treeLevelIdxs;
}


void ROAD::buildRouteOverlay(Graph_RoadNetwork& graph)
{
    METISIdxToNodeID.resize(this->numNodes);
    METISWrapper metis(this->numNodes,this->numEdges,this->fanout);
    
    // Route overlay consists of nodes for every road network vertice
    // Note: We assume that road network vertices are numbered from 0 to this->numNodes
    for (int i = 0; i < this->numNodes; ++i) {
        RouteOverlayNode roNode(i);
        this->routeOverlay.push_back(roNode);
    }

    // Add level 0 Rnet (i.e. the whole road network)
    std::unordered_set<NodeID> subgraph = graph.getNodesIDsUset();
    RnetNode rnet(constants::ROOT_PARENT_INDEX,0);
    this->rnetTree.push_back(rnet);
    this->partitionRnet(0,0,subgraph,graph,metis);
    // Release temporary objects made for partitioning
    utility::releaseSTLCollection(this->METISIdxToNodeID);
    
    this->buildShortcutTrees(graph);
    // Release temporary objects made for shortcut tree population
//     for (std::size_t i = 0; i < this->routeOverlay.size(); ++i) {
//         utility::releaseSTLCollection(this->routeOverlay[i].rnetIdxToScTreeIdx);
//     }
}

void ROAD::buildShortcutTrees(Graph_RoadNetwork& graph)
{
    std::vector<std::vector<int>> treeLevelIdxs = this->getRnetIdxsByLevel();
    
    BinaryMinHeap<EdgeWeight,NodeStatusPair> *pqueue = new BinaryMinHeap<EdgeWeight,NodeStatusPair>();
    int rnetIdx;
    
    for (int i = treeLevelIdxs.size()-1; i >= 0; --i) {
        for (std::size_t j = 0; j < treeLevelIdxs[i].size(); ++j) {
            rnetIdx = treeLevelIdxs[i][j];
            std::unordered_map<NodeID,EdgeWeight> siblingBorderDistances;
            siblingBorderDistances.reserve(this->rnetTree[rnetIdx].borders.size());
            if (this->rnetTree[rnetIdx].isLeafRnet()) {
                for (NodeID source: this->rnetTree[rnetIdx].borders) {
                    pqueue->clear();
                    this->findSSMTDistances(graph,source,this->rnetTree[rnetIdx].bordersUset,siblingBorderDistances,pqueue);
                    this->routeOverlay[source].setShortcutStartIdx(rnetIdx,this->shortcutEdges.size());
                    for (NodeID target: this->rnetTree[rnetIdx].borders) {
                        // We check for zero distance because some shortcut are to be avoided by Lemma 4 in paper
                        if (source != target && siblingBorderDistances[target] != 0) {
                            this->addShortcut(target,siblingBorderDistances[target]);
                        }
                    }
                    this->routeOverlay[source].setShortcutEndPlusOneIdx(rnetIdx,this->shortcutEdges.size());
                }
                // These are no longer needed again (and we do not serialize them)
                utility::releaseSTLCollection(this->rnetTree[rnetIdx].borders);
                utility::releaseSTLCollection(this->rnetTree[rnetIdx].bordersUset);
            } else {
                for (NodeID source: this->rnetTree[rnetIdx].borders) {
                    pqueue->clear();
                    this->findSSMTDistancesByRnet(rnetIdx,source,this->rnetTree[rnetIdx].bordersUset,siblingBorderDistances,pqueue);
                    this->routeOverlay[source].setShortcutStartIdx(rnetIdx,this->shortcutEdges.size());
                    for (NodeID target: this->rnetTree[rnetIdx].borders) {
                        // We check for zero distance because some shortcut are to be avoided by Lemma 4 in paper
                        if (source != target && siblingBorderDistances[target] != 0) {
                            this->addShortcut(target,siblingBorderDistances[target]);
                        }
                    }
                    this->routeOverlay[source].setShortcutEndPlusOneIdx(rnetIdx,this->shortcutEdges.size());
                }
                // These are no longer needed again (and we do not serialize them)
                utility::releaseSTLCollection(this->rnetTree[rnetIdx].borders);
                utility::releaseSTLCollection(this->rnetTree[rnetIdx].bordersUset);
            }
        }
    }
    delete pqueue;
}

void ROAD::partitionRnet(int rnetIdx, int rnetLevel, std::unordered_set<NodeID>& subgraph, 
                         Graph_RoadNetwork& originalGraph, METISWrapper& metis)
{
    int currentLevel = rnetLevel + 1;
    idx_t n = subgraph.size();
    
    // Note: this->METISIdxToNodeID will map the METIS idx (which will be number from
    // 0 to the size of the subgraph, to the original node id
    metis.populateMETISArrays(originalGraph,subgraph,this->METISIdxToNodeID,true);

    metis.partitionSubgraph(n);

    // Create sets of subgraphs from partitioned parent graph
    // and add them to Rnet (and this will recurse)
    std::vector<std::unordered_set<NodeID>> childGraphs(this->fanout);
    std::vector<int> childRnetIdxs(this->fanout);
    
    for (int i = 0; i < n; ++i) {
        // Create child graphs using original NodeID not the METIS idx
        childGraphs[metis.parts[i]].insert(this->METISIdxToNodeID[i]);
    }
    
    bool isLeaf = true; // Are all children leaf nodes?
    if (currentLevel+1 < this->levels) {
        // If we can partition again then no
        isLeaf = false;
    }

    int adjListStart, adjListSize;
    for (int i = 0; i < this->fanout; ++i) {
        if (childGraphs[i].size() == 0 || childGraphs[i].size() == subgraph.size()) {
            // This happens when METIS could not partition the graph.
            // Since we cannot partioning this rnet becomes a leaf and we need to re-process
            // and do the all the operations we would have done previous if it were a leaf rnet
            this->rnetTree[rnetIdx].setLeafRnet();
            for (auto nodeIt = subgraph.begin(); nodeIt != subgraph.end(); ++nodeIt) {
                adjListStart = originalGraph.getEdgeListStartIndex(*nodeIt);
                adjListSize = originalGraph.getEdgeListSize(*nodeIt);
                for (int j = adjListStart; j < adjListSize; ++j) {
                    if (subgraph.find(originalGraph.edges[j].first) == subgraph.end()) {
                        // If it's adjacent node is not in the current child graph
                        // then this node is potentially a border
                        if (this->borders[originalGraph.edges[j].first]) {
                            // We add this rnet as leaf rnet of this node so we can later search
                            // up rnet hierarchy to find which rnets this node belongs too
                            // Note: That there could more than one leaf rnet
                            this->routeOverlay[originalGraph.edges[j].first].addLeafRnet(rnetIdx);
                        }
                    }
                }
                
                // We add this rnet as leaf rnet of this node so we can later search
                // up rnet hierarchy to find which rnets the node belongs too
                // Note: That there could more than one leaf rnet e.g. border
                this->routeOverlay[*nodeIt].addLeafRnet(rnetIdx);
                if (!this->borders[*nodeIt]) {
                    // If this a node in a leaf Rnet and it's not a border, then we need to add it's 
                    // adjacent nodes to the base node of the root shortcut tree node
                    if (!this->routeOverlay[*nodeIt].hasBaseNode(0)) {
                        // We also set create a base node in the shortcut tree for this rnet
                        // so that we can add all the adjacent nodes of this node that are
                        // contained within the rnet
                        this->routeOverlay[*nodeIt].addBaseToShortcutTree(0);
                    }
                    // Note: We actually only ever process each border and rnet combination
                    // once so it is not necessary to check hasBaseNode here
                    this->routeOverlay[*nodeIt].setBaseEdgeStartIdx(0,this->shortcutEdges.size());
                    for (int j = adjListStart; j < adjListSize; ++j) {
                        this->addShortcut(originalGraph.edges[j].first,originalGraph.edges[j].second);
                    }
                    this->routeOverlay[*nodeIt].setBaseEdgeEndPlusOneIdx(0,this->shortcutEdges.size());
                }
            }

            // For each border in a leaf Rnet we must add all it's neighbours contained within this Rnet
            // to the base shortcut tree node corresponding to this Rnet
            // Note: We do this after finding all borders as we must include any border we have added 
            // artificially (i.e. because METIS partitions by cutting edges not by selecting shared nodes)
            for (NodeID border: this->rnetTree[rnetIdx].borders) {
                if (!this->routeOverlay[border].hasBaseNode(rnetIdx)) {
                    // We also set create a base node in the shortcut tree for this rnet
                    // so that we can add all the adjacent nodes contained within the rnet
                    this->routeOverlay[border].addBaseToShortcutTree(rnetIdx);
                }
                // Note: We actually only ever process each border and rnet combination
                // once so it is not necessary to check hasBaseNode here
                adjListStart = originalGraph.getEdgeListStartIndex(border);
                adjListSize = originalGraph.getEdgeListSize(border);
                this->routeOverlay[border].setBaseEdgeStartIdx(rnetIdx,this->shortcutEdges.size());
                for (int k = adjListStart; k < adjListSize; ++k) {
                    if (subgraph.find(originalGraph.edges[k].first) != subgraph.end()
                        || this->rnetTree[rnetIdx].bordersUset.find(originalGraph.edges[k].first) != this->rnetTree[rnetIdx].bordersUset.end()) {
                        // If the adjacent node is in the current subgraph then we add it to the tree
                        // Note: We need to search borders because subgraph won't include
                        // the nodes we artificially included in this Rnet (as they were already borders)
                        // but they would have added to the borders vector while processing another subgraph
                        this->addShortcut(originalGraph.edges[k].first,originalGraph.edges[k].second);
                    }
                }
                this->routeOverlay[border].setBaseEdgeEndPlusOneIdx(rnetIdx,this->shortcutEdges.size());
            }
            return;
        }
    }
    
    // Create child Rnets for each partition
    int childRnetIdx;
    for (int i = 0; i < this->fanout; ++i) {
        assert(childGraphs[i].size() != 0 && childGraphs[i].size() != subgraph.size() && "METIS partitioning failed, should return before reaching here");
        childRnetIdx = this->rnetTree.size();
        RnetNode rnet(rnetIdx,currentLevel);
        this->rnetTree.push_back(rnet);
        childRnetIdxs[i] = childRnetIdx;
        // Make new Rnet a child of parent
        this->rnetTree[rnetIdx].addChild(childRnetIdx);
        if (isLeaf) {
            // If this is the last level then it is a leaf Rnet
            this->rnetTree[childRnetIdx].setLeafRnet();
        }
    }

    // Note: Subgraphs provided by METIS are separated on edges, i.e. there will be
    // two border nodes. So we must manually determine which of these will be the 
    // common border required for Rnets. To do this, when we find an adjacent node
    // to a node in the subgraph, and if it is a border of another subgraph we artificially
    // consider the adjacent node a node of the current subgraph. If it is not a border 
    // of the the other subgraph, then we consider the node in the current subgraph the border
    
    for (int i = 0; i < this->fanout; ++i) {
        childRnetIdx = childRnetIdxs[i];
        for (auto nodeIt = childGraphs[i].begin(); nodeIt != childGraphs[i].end(); ++nodeIt) {
            if (this->borders[*nodeIt]) {
                // If it is border at a higher level (this must necessarily be a 
                // parent rnet of the current rnet) then it is a border here too
                this->routeOverlay[*nodeIt].addRnetToShortcutTree(childRnetIdx,rnetIdx,this->rnetTree[childRnetIdx].getLevel());
                this->rnetTree[childRnetIdx].addBorder(*nodeIt);
            } 
            adjListStart = originalGraph.getEdgeListStartIndex(*nodeIt);
            adjListSize = originalGraph.getEdgeListSize(*nodeIt);
            for (int j = adjListStart; j < adjListSize; ++j) {
                if (childGraphs[i].find(originalGraph.edges[j].first) == childGraphs[i].end()) {
                    // If it's adjacent node is not in the current child graph
                    // then this node is potentially a border
                    if (this->borders[originalGraph.edges[j].first]) {
                        // If the adjacent node in another subgraph is a border already
                        // then it is also a border of this subgraph
                        // Note: This node (*nodeIt) could still be a border to a different subgraph (if
                        // we find an adjacent node not in this subgraph and is not a border already)
                        
                        // Now we need to process this adjacent node as if it were a part of this subgraph
                        // Note: This is the only opportunity to do so or it will be missing as we are
                        // only artificially considering it a node of the current subgraph
                        
                        this->rnetTree[childRnetIdx].addBorder(originalGraph.edges[j].first);
                        // Since this adjacent node is a border then this Rnet belongs in its shortcut tree
                        if (!this->routeOverlay[originalGraph.edges[j].first].hasRnet(childRnetIdx)) {
                            this->routeOverlay[originalGraph.edges[j].first].addRnetToShortcutTree(childRnetIdx,rnetIdx,this->rnetTree[childRnetIdx].getLevel());
                        }
                        if (isLeaf) {
                            // We add this rnet as leaf rnet of this node so we can later search
                            // up rnet hierarchy to find which rnets this node belongs too
                            // Note: That there could more than one leaf rnet
                            this->routeOverlay[originalGraph.edges[j].first].addLeafRnet(childRnetIdx);
                        }
                    } else if (!this->borders[*nodeIt]) {
                        // If the adjacent node is not a border but in another subgraph then this node must be a border
                        // Note: This means that this node will become a border for that subgraph - we just have not
                        // added borders for the subgraph yet
                        this->rnetTree[childRnetIdx].addBorder(*nodeIt);
                        this->borders[*nodeIt] = true;
                        // If this node is a border then this Rnet belongs in this nodes shortcut tree
                        this->routeOverlay[*nodeIt].addRnetToShortcutTree(childRnetIdx,rnetIdx,this->rnetTree[childRnetIdx].getLevel());
                    }
                }
            }
            
            if (isLeaf) {
                // We add this rnet as leaf rnet of this node so we can later search
                // up rnet hierarchy to find which rnets the node belongs too
                // Note: That there could more than one leaf rnet e.g. border
                this->routeOverlay[*nodeIt].addLeafRnet(childRnetIdx);
                if (!this->borders[*nodeIt]) {
                    // If this a node in a leaf Rnet and it's not a border, then we need to add it's 
                    // adjacent nodes to the base node of the root shortcut tree node
                    if (!this->routeOverlay[*nodeIt].hasBaseNode(0)) {
                        // We also set create a base node in the shortcut tree for this rnet
                        // so that we can add all the adjacent nodes of this node that are
                        // contained within the rnet
                        this->routeOverlay[*nodeIt].addBaseToShortcutTree(0);
                    }
                    // Note: We actually only ever process each border and rnet combination
                    // once so it is not necessary to check hasBaseNode here
                    this->routeOverlay[*nodeIt].setBaseEdgeStartIdx(0,this->shortcutEdges.size());
                    for (int j = adjListStart; j < adjListSize; ++j) {
                        this->addShortcut(originalGraph.edges[j].first,originalGraph.edges[j].second);
                    }
                    this->routeOverlay[*nodeIt].setBaseEdgeEndPlusOneIdx(0,this->shortcutEdges.size());
                }
            }
        }
        
        if (isLeaf) {
            // For each border in a leaf Rnet we must add all it's neighbours contained within this Rnet
            // to the base shortcut tree node corresponding to this Rnet
            // Note: We do this after finding all borders as we must include any border we have added 
            // artificially (i.e. because METIS partitions by cutting edges not by selecting shared nodes)
            for (NodeID border: this->rnetTree[childRnetIdx].borders) {
                if (!this->routeOverlay[border].hasBaseNode(childRnetIdx)) {
                    // We also set create a base node in the shortcut tree for this rnet
                    // so that we can add all the adjacent nodes contained within the rnet
                    this->routeOverlay[border].addBaseToShortcutTree(childRnetIdx);
                }
                // Note: We actually only ever process each border and rnet combination
                // once so it is not necessary to check hasBaseNode here
                adjListStart = originalGraph.getEdgeListStartIndex(border);
                adjListSize = originalGraph.getEdgeListSize(border);
                this->routeOverlay[border].setBaseEdgeStartIdx(childRnetIdx,this->shortcutEdges.size());
                for (int k = adjListStart; k < adjListSize; ++k) {
                    if (childGraphs[i].find(originalGraph.edges[k].first) != childGraphs[i].end()
                        || this->rnetTree[childRnetIdx].bordersUset.find(originalGraph.edges[k].first) != this->rnetTree[childRnetIdx].bordersUset.end()) {
                        // If the adjacent node is in the current subgraph then we add it to the tree
                        // Note: We need to search borders because childGraphs[i] won't include
                        // the nodes we artificially included in this Rnet (as they were already borders)
                        // but they would have added to the borders vector while processing another subgraph
                        this->addShortcut(originalGraph.edges[k].first,originalGraph.edges[k].second);
                    }
                }
                this->routeOverlay[border].setBaseEdgeEndPlusOneIdx(childRnetIdx,this->shortcutEdges.size());
            }
        }
    }
    
    // Recursively partition subgraphs further if level limit not reached
    if (!isLeaf) {
        for (int i = 0; i < this->fanout; ++i) {
            this->partitionRnet(childRnetIdxs[i],currentLevel,childGraphs[i],originalGraph,metis);
            utility::releaseSTLCollection(childGraphs[i]);
        }
    }    
}

// Note: Caller must check results map for any distance that are zero and omit adding shortcuts.
// We use 0 to shortcuts that satisfy Lemma 4 where want to discard shortcuts that have borders in the path
// (as by defintion these are guaranteed to be made up of two other shortcuts)
void ROAD::findSSMTDistancesByRnet(int rnetIdx, NodeID source, std::unordered_set<NodeID>& targetSet, 
                                   std::unordered_map<NodeID,EdgeWeight>& results,
                                   BinaryMinHeap<EdgeWeight,NodeStatusPair> *pqueue)
{
    std::vector<bool> isNodeSettled;
    isNodeSettled.resize(this->numNodes);
    
    EdgeWeight minDist, newDistance;
    NodeID minDistNodeID, neighbourNodeID;
    
    // Initialize with priority queue with source node
    pqueue->insert(NodeStatusPair(source,false),0);
    // Note: Second element in pair holds whether the shortest path
    // to the first element (the node) contains a border node
    std::size_t targetsFound = 0;
    bool minNodeIsBorder;
    
    int level = this->rnetTree[rnetIdx].getLevel();
    
    int shortcutListStart, shortcutListSize;

    while (pqueue->size() > 0 && targetsFound < targetSet.size()) {
        // Extract and remove node with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue->getMinKey();
        NodeStatusPair minElement = pqueue->extractMinElement();
        minDistNodeID = minElement.first;
        if (!isNodeSettled[minDistNodeID]) {
            isNodeSettled[minDistNodeID] = 1;

            minNodeIsBorder = false;
            if (targetSet.find(minDistNodeID) != targetSet.end()) {
                // We have found the SP distance to one of the target set
                // So we add it to result set and remove it from input set
                if (minElement.second) {
                    // In this case, the target border has another border in it's path,
                    // so we omit this shortcut (by marking it 0 - the caller must filter it)
                    // as there will be a shortcut from the other border in the path
                    results[minDistNodeID] = 0;
                } else {
                    results[minDistNodeID] = minDist; // This is the final SP distance from source
                }
                if (minDistNodeID != source) {
                    // Obviously the source will be border in every path so we ignore it
                    minNodeIsBorder = true;
                }
                ++targetsFound;
            }

            // Inspect the shortcuts of each child Rnet for the min node
            std::vector<int> children = this->routeOverlay[minDistNodeID].getShortcutsIdxsAtLevel(level+1);
            
            for (std::size_t i = 0; i < children.size(); ++i) {
                shortcutListStart = this->routeOverlay[minDistNodeID].shortcutTree[children[i]].shortcutStartIdx;
                shortcutListSize = this->routeOverlay[minDistNodeID].shortcutTree[children[i]].shortcutEndPlusOneIdx;
                
                for (int j = shortcutListStart; j < shortcutListSize; ++j) {
                    neighbourNodeID = this->shortcutEdges[j].first;
                    // Only update those we haven't already settled
                    if (!isNodeSettled[neighbourNodeID]) {
                        newDistance = minDist + this->shortcutEdges[j].second;
                        if (minNodeIsBorder || minElement.second) {
                            // Current shortest path to neighbourNodeID is through minDistNodeID
                            // If it is a border or the path to it is through then the current
                            // shortest path involves a border
                            pqueue->insert(NodeStatusPair(neighbourNodeID,true),newDistance);
                        } else {
                            // This means that current shortest path to neighbourNodeID does not
                            // involve a border so set false in case it was set true in earlier
                            pqueue->insert(NodeStatusPair(neighbourNodeID,false),newDistance);
                        }
                    }
                }                
            }
        }
    }
    
//     if (targetsFound < targetSet.size()) {
//         this->rnetTree[rnetIdx].print();
//         std::cout << "source = " << source << std::endl;
//         std::cout << "targetsFound = " << targetsFound << std::endl;
//         std::cout << "targetSet.size() = " << targetSet.size() << std::endl;
//         std::cout << "PQueue size: "  << pqueue->size() << std::endl;
//         std::cout << "Not Found Targets: ";
//         for (auto target: targetSet) {
//             bool found = false;
//             for (auto it: results) {
//                 if (it.first == target) {
//                     found = true;
//                     break;
//                 }
//             }
//             if (!found) {
//                 std::cout << target << ", ";
//                 //this->routeOverlay[target].print();
//             }
//         }
//         std::cout << std::endl;
//     }
//      assert (targetSet.size() == targetsFound && "Dijkstra search could not find all targets in target set");
}

// Note: Caller must check results map for any distance that are zero and omit adding shortcuts.
// We use 0 to shortcuts that satisfy Lemma 4 where want to discard shortcuts that have borders in the path
// (as by defintion these are guaranteed to be made up of two other shortcuts)
void ROAD::findSSMTDistances(Graph_RoadNetwork& graph, NodeID source,
                             std::unordered_set<NodeID>& targetSet, 
                             std::unordered_map<NodeID,EdgeWeight>& results,
                             BinaryMinHeap<EdgeWeight,NodeStatusPair> *pqueue)
{
    std::vector<bool> isNodeSettled;
    isNodeSettled.resize(graph.getNumNodes());
    
    EdgeWeight minDist, newDistance;
    NodeID minDistNodeID, neighbourNodeID;
    int adjListStart, adjListSize;
    
    // Initialize with priority queue with source node
    pqueue->insert(NodeStatusPair(source,false),0);
    // Note: Second element in pair holds whether the shortest path
    // to the first element (the node) contains a border node
    std::size_t targetsFound = 0;
    bool minNodeIsBorder;
    
    while (pqueue->size() > 0 && targetsFound < targetSet.size()) {
        // Extract and remove node with smallest distance from source
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue->getMinKey();
        NodeStatusPair minElement = pqueue->extractMinElement();
        minDistNodeID = minElement.first;
        if (!isNodeSettled[minDistNodeID]) {
            isNodeSettled[minDistNodeID] = 1;

            minNodeIsBorder = false;
            if (targetSet.find(minDistNodeID) != targetSet.end()) {
                // We have found the SP distance to one of the target set
                // So we add it to result set and remove it from input set
                if (minElement.second) {
                    // In this case, the target border has another border in it's path,
                    // so we omit this shortcut (by marking it 0 - the caller must filter it)
                    // as there will be a shortcut from the other border in the path
                    results[minDistNodeID] = 0;
                } else {
                    results[minDistNodeID] = minDist; // This is the final SP distance from source
                }
                if (minDistNodeID != source) {
                    // Obviously the source will be border in every path so we avoid this
                    minNodeIsBorder = true;
                }
                ++targetsFound;
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistNodeID);
            adjListSize = graph.getEdgeListSize(minDistNodeID);
            for (int i = adjListStart; i < adjListSize; ++i) {
                neighbourNodeID = graph.edges[i].first;
                // Only update those we haven't already settled
                if (!isNodeSettled[neighbourNodeID]) {
                    newDistance = minDist + graph.edges[i].second;
                    if (minNodeIsBorder || minElement.second) {
                        // Current shortest path to neighbourNodeID is through minDistNodeID
                        // If it is a border or the path to it is through then the current
                        // shortest path involves a border
                        pqueue->insert(NodeStatusPair(neighbourNodeID,true),newDistance);
                    } else {
                        // This means that current shortest path to neighbourNodeID does not
                        // involve a border so set false in case it was set true in earlier
                        pqueue->insert(NodeStatusPair(neighbourNodeID,false),newDistance);
                    }
                }
            }
        }
    }
//     assert (targetSet.size() == targetsFound && "Dijkstra search could not find all targets in target set");
}

int ROAD::getRnetTreeSize()
{
    return this->rnetTree.size();
}

const std::vector<int>& ROAD::getLeafRnets(NodeID node) const
{
    return this->routeOverlay[node].getLeafRnets();
}

int ROAD::getParentRnet(int rnetIdx)
{
    return this->rnetTree[rnetIdx].getParent();
}

void ROAD::getKNNs(AssociationDirectory& assocDir, unsigned int k, NodeID queryNode, std::vector<NodeID>& kNNs, 
                   std::vector<EdgeWeight>& kNNDistances) {
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("nodes_bypassed",0);
    this->stats.initialiseStatistic("rnets_bypassed",0);
    this->stats.initialiseStatistic("rnets_evaluated",0);
    this->stats.initialiseStatistic("shortcuts_traversed",0);
    
    // Count number of nodes in each leaf Rnet
    std::vector<bool> bypassedRnet(this->rnetTree.size(),false);
    if (this->rnetNodeCounts.size() == 0) {
        // Then we have never computed the node counts
        this->rnetNodeCounts.assign(this->rnetTree.size(),0);
        std::unordered_set<int> rnetLeafs;
        for (std::size_t i = 0; i < this->routeOverlay.size(); ++i) {
            for (std::size_t j = 0; j < this->routeOverlay[i].leafIdxs.size(); ++j) {
                rnetLeafs.insert(this->routeOverlay[i].leafIdxs[j]);
                rnetNodeCounts[this->routeOverlay[i].leafIdxs[j]]++;
            }
        }
        // Propagate up leaf Rnet counts 
        int parentRnetIdx;
        for (auto it = rnetLeafs.begin(); it != rnetLeafs.end(); ++it) {
            parentRnetIdx = this->rnetTree[*it].getParent();
            while (parentRnetIdx != -1) {
                rnetNodeCounts[parentRnetIdx] += rnetNodeCounts[*it];
                // Go up a level (until we reach root)
                parentRnetIdx = this->rnetTree[parentRnetIdx].getParent();   
            }
        }
    }
#endif
    
    BinaryMinHeap<EdgeWeight,NodeID> pqueue;
    EdgeWeight minDist;
    NodeID minDistNodeID;
    std::vector<bool> isNodeVisited(this->numNodes,false);
//     std::vector<int> stack;
    std::deque<int> stack;
    int scIdx, shortcutListStart, shortcutListSize;
    
    // Initialize with priority queue with query node ID
    pqueue.insert(queryNode,0);
    
    while (pqueue.size() > 0) {
        // Extract and remove node with smallest distance from query point
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue.getMinKey();
        minDistNodeID = pqueue.extractMinElement();

        if (!isNodeVisited[minDistNodeID]) {
            if (assocDir.isObject(minDistNodeID)) {
                // If the minimum is an object we have found a kNN
                kNNs.push_back(minDistNodeID);
                kNNDistances.push_back(minDist);
                if (kNNs.size() == k) {
                    // If this is the kth nearest neighbour object
                    // then we no longer need to inspect shortcut trees
                    break;
                }
            }
            
            // Choose Path 
            stack.clear();
            // The root of the shortcut hierarchy presents the top-level Rnet which contains whole graph
            // Each graph node will be part of at least two Rnet (and one of those is the originating Rnet).
            // But it may have more than 2 Rnets, for example if the border is shared by 3 Rnets (theoretically
            // this could be upto fanout Rnets). So for every shortcut tree node we have add all children.
            
            // Note: If this graph node is not a border then root will have a special shortcut tree "base node".
            // If it were a border however, any of it's shortcut tree node has only one base node at most, 
            // but there may be multipl base nodes in the whole tree.
            for (std::size_t i = 0; i < this->routeOverlay[minDistNodeID].shortcutTree[0].children.size(); ++i) {
                stack.push_back(this->routeOverlay[minDistNodeID].shortcutTree[0].children[i]);
            }            
            
            // Note: It is possible that this node is part border of more than two Rnets
            // at the same level so we need to keep adding shortcuts until stack is empty
            while (stack.size() > 0) {
                scIdx = stack.front();
                stack.pop_front();
                shortcutListStart = this->routeOverlay[minDistNodeID].shortcutTree[scIdx].shortcutStartIdx;
                shortcutListSize = this->routeOverlay[minDistNodeID].shortcutTree[scIdx].shortcutEndPlusOneIdx;             
                
                if (!this->routeOverlay[minDistNodeID].shortcutTree[scIdx].isBaseNode()) {
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("rnets_evaluated",1);
#endif
                    if(!assocDir.hasObject(this->routeOverlay[minDistNodeID].shortcutTree[scIdx].getRnetIdx())) {
#if defined(COLLECT_STATISTICS)
                        int rnetIdx = this->routeOverlay[minDistNodeID].shortcutTree[scIdx].getRnetIdx();
                        if (!bypassedRnet[rnetIdx]) {
                            bypassedRnet[rnetIdx] = true; // Don't count more than once
                            this->stats.incrementStatistic("rnets_bypassed",1);
                            this->stats.incrementStatistic("nodes_bypassed",this->rnetNodeCounts[rnetIdx]);
                        }
#endif
                        // If this is not a leaf Rnet, and it does not contain any nodes
                        // then we add all shortcuts to priority (we skip this Rnet essentially)
                        for (int i = shortcutListStart; i < shortcutListSize; ++i) {
                            if (!isNodeVisited[this->shortcutEdges[i].first]) {
#if defined(COLLECT_STATISTICS)
                                this->stats.incrementStatistic("shortcuts_traversed",1);
#endif
                                pqueue.insert(this->shortcutEdges[i].first,minDist+this->shortcutEdges[i].second);
                            }
                        }
                    } else {
                        // However if it does have an object, then we need to traverse
                        // down the shortcut tree to check the next level sub-Rnet
                        for (std::size_t i = 0; i < this->routeOverlay[minDistNodeID].shortcutTree[scIdx].children.size(); ++i) {
                            stack.push_back(this->routeOverlay[minDistNodeID].shortcutTree[scIdx].children[i]);
                        }
                    }
                } else {
                    // If this Rnet is a base node that means we exhausted all Rnet skipping
                    // opportunities and we add all the "shortcuts" (which are are actually
                    // some of the original edges leaving this node).
                    // Note: We say "some" of the original edges, because there maybe more than one
                    // base node in the shortcut tree which the original nodes would be split amongst
                    for (int i = shortcutListStart; i < shortcutListSize; ++i) {
                        if (!isNodeVisited[this->shortcutEdges[i].first]) {
                            pqueue.insert(this->shortcutEdges[i].first,minDist+this->shortcutEdges[i].second);
                        }
                    }
                }
            }
            
            isNodeVisited[minDistNodeID] = 1;
        }
    }
}

double ROAD::computeIndexSize()
{
    double memoryUsage = 0;
    for (std::size_t i = 0; i < this->routeOverlay.size(); ++i) {
        memoryUsage += this->routeOverlay[i].computeIndexSizeBytes();
    }
    for (std::size_t i = 0; i < this->rnetTree.size(); ++i) {
        memoryUsage += this->rnetTree[i].computeIndexSizeBytes();
    }
    memoryUsage += sizeof(NodeEdgeWeightPair)*this->shortcutEdges.size();
    memoryUsage += this->borders.size()/8;
    return memoryUsage/(1024*1024);
}

double ROAD::computeMemoryUsage()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    memoryUsage += this->networkName.size();
    for (std::size_t i = 0; i < this->routeOverlay.size(); ++i) {
        memoryUsage += this->routeOverlay[i].computeMemoryUsageBytes();
    }
    memoryUsage += sizeof(RouteOverlayNode)*(this->routeOverlay.capacity()-this->routeOverlay.size());
    for (std::size_t i = 0; i < this->rnetTree.size(); ++i) {
        memoryUsage += this->rnetTree[i].computeMemoryUsageBytes();
    }
    memoryUsage += sizeof(RnetNode)*(this->rnetTree.capacity()-this->rnetTree.size());
    memoryUsage += sizeof(NodeID)*this->METISIdxToNodeID.capacity();
    memoryUsage += sizeof(NodeEdgeWeightPair)*this->shortcutEdges.capacity();
    memoryUsage += this->borders.capacity()/8;
    return memoryUsage/(1024*1024);
}

void ROAD::addShortcut(NodeID target, EdgeWeight distance)
{
    this->shortcutEdges.push_back(NodeEdgeWeightPair(target,distance));
}

/*
 * ShortcutTreeNode
 */

ShortcutTreeNode::ShortcutTreeNode(int rnetIdx, int rnetLevel): 
     rnetIdx(rnetIdx), rnetLevel(rnetLevel), hasBase(false)
{

}

void ShortcutTreeNode::addChild(int childIdx)
{
    this->children.push_back(childIdx);
}

const std::vector<int>& ShortcutTreeNode::getChildren() const
{
    return this->children;
}

int ShortcutTreeNode::getLevel()
{
    return this->rnetLevel;
}

int ShortcutTreeNode::getBaseIdx()
{
    return this->children[0];
}

bool ShortcutTreeNode::hasBaseNode()
{
    return this->hasBase;
}

void ShortcutTreeNode::setHasBaseNode()
{
    this->hasBase = true;
}

int ShortcutTreeNode::getRnetIdx()
{
    return this->rnetIdx;
}

bool ShortcutTreeNode::isBaseNode()
{
    return this->rnetIdx == -1;
}

void ShortcutTreeNode::print()
{
    std::cout << "Rnet Index: " << this->rnetIdx << std::endl;
    std::cout << "Rnet Level: " << this->rnetLevel << std::endl;
    std::cout << "isBase: " << this->isBaseNode() << std::endl;
    std::cout << "Children (" << this->children.size() << "): ";
    for (std::size_t i = 0; i < this->children.size(); ++i) {
        if (i != 0) {
            std::cout << ", ";
        }
        std::cout << i << " => " << children[i];
    }
    std::cout << std::endl;
    std::cout << std::endl;
}

double ShortcutTreeNode::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(int)*4;
    memoryUsage += sizeof(bool);
    memoryUsage += sizeof(this->children);
    memoryUsage += sizeof(int)*this->children.size();
    return memoryUsage;
}

double ShortcutTreeNode::computeMemoryUsageBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    memoryUsage += sizeof(int)*this->children.capacity();
    return memoryUsage;
}

void ShortcutTreeNode::setShortcutStartIdx(int shortcutStartIdx)
{
    this->shortcutStartIdx = shortcutStartIdx;
}

void ShortcutTreeNode::setShortcutEndPlusOneIdx(int shortcutEndPlusOneIdx)
{
    this->shortcutEndPlusOneIdx = shortcutEndPlusOneIdx;
}

/*
 * RouteOverlayNode
 */

RouteOverlayNode::RouteOverlayNode(NodeID id)
{
    // All shortcut tree have root Rnet (i.e. whole graph at top)
    ShortcutTreeNode scTreeNode(0,0);
    this->shortcutTree.push_back(scTreeNode);
}

const std::vector<int>& RouteOverlayNode::getLeafRnets() const
{
    return this->leafIdxs;
}

void RouteOverlayNode::addRnetToShortcutTree(int rnetIdx, int parentRnetIdx, int rnetLevel)
{
    int parentScTreeIdx;
    if (this->rnetIdxToScTreeIdx.find(parentRnetIdx) != this->rnetIdxToScTreeIdx.end()) {
        // If the parent also exists in this shortcut tree then we add the rnet as a child
        parentScTreeIdx = this->rnetIdxToScTreeIdx[parentRnetIdx];
    } else {
        // Otherwise the parent is the root shortcut tree node
        parentScTreeIdx = 0;
    }
    
    int scTreeIdx = this->shortcutTree.size();
    ShortcutTreeNode scTreeNode(rnetIdx,rnetLevel);
    this->shortcutTree.push_back(scTreeNode);
    this->rnetIdxToScTreeIdx[rnetIdx] = scTreeIdx;
    this->shortcutTree[parentScTreeIdx].addChild(scTreeIdx);
}

void RouteOverlayNode::addLeafRnet(int rnetIdx)
{
    this->leafIdxs.push_back(rnetIdx);
}

void RouteOverlayNode::addBaseToShortcutTree(int rnetIdx)
{
    // Add base child for this rnet
    int parentScTreeIdx;
    if (this->rnetIdxToScTreeIdx.find(rnetIdx) != this->rnetIdxToScTreeIdx.end()) {
        // If the rnet also exists in this shortcut tree then we add the base as a child
        parentScTreeIdx = this->rnetIdxToScTreeIdx[rnetIdx];
    } else {
        // Otherwise the rnet is the root shortcut tree node
        parentScTreeIdx = 0;
    }
    int scTreeIdx = this->shortcutTree.size();
    ShortcutTreeNode scTreeNode(-1,-1);
    this->shortcutTree.push_back(scTreeNode);
    this->shortcutTree[parentScTreeIdx].addChild(scTreeIdx);
    this->shortcutTree[parentScTreeIdx].setHasBaseNode();
}

bool RouteOverlayNode::hasBaseNode(int rnetIdx)
{
    int parentScTreeIdx;
    if (rnetIdx != 0) {
        // If the rnet also exists in this shortcut tree then we add the base as a child
        parentScTreeIdx = this->rnetIdxToScTreeIdx[rnetIdx];
    } else {
        // Otherwise the rnet is the root shortcut tree node
        parentScTreeIdx = 0;
    }
    return this->shortcutTree[parentScTreeIdx].hasBaseNode();
}

const std::vector<int>& RouteOverlayNode::getShortcutTreeChildrenByRnetIdx(int rnetIdx)
{
    int scTreeIdx = this->rnetIdxToScTreeIdx[rnetIdx];
    return this->shortcutTree[scTreeIdx].getChildren();
}

const std::vector<int>& RouteOverlayNode::getShortcutTreeChildren(int scTreeIdx) const
{
    return this->shortcutTree[scTreeIdx].getChildren();
}

int RouteOverlayNode::getRnetIdx(int scTreeIdx)
{
    return this->shortcutTree[scTreeIdx].getRnetIdx();
}

bool RouteOverlayNode::hasRnet(int rnetIdx)
{
    return this->rnetIdxToScTreeIdx.find(rnetIdx) != this->rnetIdxToScTreeIdx.end();
}

const std::vector<ShortcutTreeNode>& RouteOverlayNode::getShortcutTree() const
{
    return this->shortcutTree;
}

void RouteOverlayNode::print()
{
    std::cout << "Leaf Indexes (" << this->leafIdxs.size() << "): ";
    for (std::size_t i = 0; i < this->leafIdxs.size(); ++i) {
        if (i != 0) {
            std::cout << ", ";
        }
        std::cout << i << " => " << leafIdxs[i];
    }
    std::cout << std::endl;
    std::cout << "Shortcut Tree (" << this->shortcutTree.size() << "): " << std::endl;
    for (std::size_t i = 0; i < this->shortcutTree.size(); ++i) {
        std::cout << "\nShortcut Tree Node " << i << ":";
        shortcutTree[i].print();
    }
    std::cout << std::endl;    
}

std::vector<int> RouteOverlayNode::getShortcutsIdxsAtLevel(int level)
{
    std::vector<int> levelIdxs;
    std::vector<int> stack;
    
    stack.push_back(0); // Insert root node
    while (stack.size() > 0) {
        int scIdx = stack[stack.size()-1];
        stack.pop_back();
        
        if (this->shortcutTree[scIdx].getLevel() == level) {
            levelIdxs.push_back(scIdx);
        } else if (this->shortcutTree[scIdx].getLevel() < level) {
            if (!this->hasBaseNode(this->shortcutTree[scIdx].getRnetIdx())) {
                stack.insert(stack.end(),this->shortcutTree[scIdx].children.begin(),this->shortcutTree[scIdx].children.end());
            }
        }
    }
    return levelIdxs;
}

double RouteOverlayNode::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(this->leafIdxs);
    memoryUsage += sizeof(int)*this->leafIdxs.size();
    memoryUsage += sizeof(this->shortcutTree);
    for (std::size_t i = 0; i < this->shortcutTree.size(); ++i) {
        memoryUsage += this->shortcutTree[i].computeIndexSizeBytes();
    }
    return memoryUsage;
}

double RouteOverlayNode::computeMemoryUsageBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    memoryUsage += sizeof(int)*this->leafIdxs.capacity();
    for (std::size_t i = 0; i < this->shortcutTree.size(); ++i) {
        memoryUsage += this->shortcutTree[i].computeMemoryUsageBytes();
    }
    memoryUsage += sizeof(ShortcutTreeNode)*(this->shortcutTree.capacity()-this->shortcutTree.size());
//     memoryUsage += utility::estimateUnorderedMapMemoryUsageBytes(this->rnetIdxToScTreeIdx.size(),sizeof(std::pair<int,int>),this->rnetIdxToScTreeIdx.bucket_count());    
    memoryUsage += sizeof(int)*2*this->rnetIdxToScTreeIdx.size();
    return memoryUsage;
}

void RouteOverlayNode::setShortcutStartIdx(int rnetIdx, int shortcutStartIdx)
{
    int scTreeIdx = this->rnetIdxToScTreeIdx[rnetIdx];
    this->shortcutTree[scTreeIdx].setShortcutStartIdx(shortcutStartIdx);
}

void RouteOverlayNode::setShortcutEndPlusOneIdx(int rnetIdx, int shortcutEndPlusOneIdx)
{
    int scTreeIdx = this->rnetIdxToScTreeIdx[rnetIdx];
    this->shortcutTree[scTreeIdx].setShortcutEndPlusOneIdx(shortcutEndPlusOneIdx);
}

void RouteOverlayNode::setBaseEdgeStartIdx(int rnetIdx, int shortcutStartIdx)
{
    // Assume rnetIdx is a leaf Rnet
    int parentScTreeIdx;
    if (this->rnetIdxToScTreeIdx.find(rnetIdx) != this->rnetIdxToScTreeIdx.end()) {
        // If the rnet also exists in this shortcut tree then we add the base as a child
        parentScTreeIdx = this->rnetIdxToScTreeIdx[rnetIdx];
    } else {
        // Otherwise the rnet is the root shortcut tree node
        parentScTreeIdx = 0;
    }
    int scTreeIdx = this->shortcutTree[parentScTreeIdx].getBaseIdx();
    this->shortcutTree[scTreeIdx].setShortcutStartIdx(shortcutStartIdx);
}

void RouteOverlayNode::setBaseEdgeEndPlusOneIdx(int rnetIdx, int shortcutEndPlusOneIdx)
{
    // Assume rnetIdx is a leaf Rnet
    int parentScTreeIdx;
    if (this->rnetIdxToScTreeIdx.find(rnetIdx) != this->rnetIdxToScTreeIdx.end()) {
        // If the rnet also exists in this shortcut tree then we add the base as a child
        parentScTreeIdx = this->rnetIdxToScTreeIdx[rnetIdx];
    } else {
        // Otherwise the rnet is the root shortcut tree node
        parentScTreeIdx = 0;
    }
    int scTreeIdx = this->shortcutTree[parentScTreeIdx].getBaseIdx();
    this->shortcutTree[scTreeIdx].setShortcutEndPlusOneIdx(shortcutEndPlusOneIdx);
}

/*
 * RnetNode
 */

RnetNode::RnetNode(int parentRnetIdx, int level): parentRnetIdx(parentRnetIdx), level(level), isLeaf(false), numBorders(0)
{
    
}

int RnetNode::getParent()
{
    return this->parentRnetIdx;
}

bool RnetNode::isLeafRnet()
{
    return this->isLeaf;
}

void RnetNode::setLeafRnet()
{
    this->isLeaf = true;
}

void RnetNode::addChild(int childRnetIdx)
{
    this->children.push_back(childRnetIdx);
}

int RnetNode::getLevel()
{
    return this->level;
}

void RnetNode::addBorder(NodeID node)
{
    if (!this->isBorder(node)) {
        ++this->numBorders;
        this->borders.push_back(node);
        this->bordersUset.insert(node);
    }
}

bool RnetNode::isBorder(NodeID node)
{
    return this->bordersUset.find(node) != this->bordersUset.end();
}

const std::vector<int>& RnetNode::getChildren() const
{
    return this->children;
}

const std::vector<NodeID>& RnetNode::getBorders() const
{
    return this->borders;
}

const std::unordered_set<NodeID>& RnetNode::getBordersUset() const
{
    return this->bordersUset;
}

void RnetNode::print()
{
    std::cout << "Parent Index: " << this->parentRnetIdx << std::endl;
    std::cout << "isLeaf: " << this->isLeaf << std::endl;
    std::cout << "Level: " << this->level << std::endl;
    std::cout << "Children (" << this->children.size() << "): ";
    for (std::size_t i = 0; i < this->children.size(); ++i) {
        if (i != 0) {
            std::cout << ", ";
        }
        std::cout << i << " => " << children[i];
    }
    std::cout << std::endl;
    std::cout << "Borders (" << this->borders.size() << "): ";
    for (std::size_t i = 0; i < this->borders.size(); ++i) {
        if (i != 0) {
            std::cout << ", ";
        }
        std::cout << i << " => " << borders[i];
    }
    std::cout << std::endl << std::endl;
}

double RnetNode::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(int)*2;
    memoryUsage += sizeof(bool);
    memoryUsage += sizeof(this->children);
    memoryUsage += sizeof(int)*this->children.size();
    return memoryUsage;
}

double RnetNode::computeMemoryUsageBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    memoryUsage += sizeof(int)*this->children.capacity();
    memoryUsage += sizeof(NodeID)*this->borders.capacity();
    memoryUsage += utility::estimateUnorderedMapMemoryUsageBytes(this->bordersUset.size(),sizeof(NodeID),this->bordersUset.bucket_count()); 
    return memoryUsage;
}

int RnetNode::getNumBorders()
{
    return this->numBorders;
}

/*
 * AssociationDirectory
 */

AssociationDirectory::AssociationDirectory(std::string setType, double setDensity, int setVariable, int setSize, int rnetTreeSize):
                             objSetType(setType), objSetDensity(setDensity), objSetVariable(setVariable), objSetSize(setSize)
{
    this->rnetAssociationDirectory.assign(rnetTreeSize,false);
}

std::string AssociationDirectory::getObjSetType()
{
    return this->objSetType;
}

double AssociationDirectory::getObjSetDensity()
{
    return this->objSetDensity;
}

int AssociationDirectory::getObjSetSize()
{
    return this->objSetSize;
}

int AssociationDirectory::getObjSetVariable()
{
    return this->objSetVariable;
}

void AssociationDirectory::addObject(NodeID node)
{
    this->objectSet.insert(node);
}

void AssociationDirectory::addRnetAssociation(int rnetIdx)
{
    this->rnetAssociationDirectory[rnetIdx] = true;
}

bool AssociationDirectory::hasObject(int rnetIdx)
{
    return this->rnetAssociationDirectory[rnetIdx];
}

bool AssociationDirectory::isObject(NodeID node)
{
    return this->objectSet.find(node) != this->objectSet.end();
}

double AssociationDirectory::computeIndexSize()
{
    double memoryUsage = 0;
    memoryUsage += this->rnetAssociationDirectory.size()/8; // std::vector<bool> only use 1 bit per element
    memoryUsage += sizeof(NodeID)*this->objectSet.size();
    return memoryUsage/(1024*1024);
}

double AssociationDirectory::computeMemoryUsage()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    memoryUsage += this->objSetType.size();
    memoryUsage += this->rnetAssociationDirectory.capacity()/8; // std::vector<bool> only use 1 bit per element
    memoryUsage += utility::estimateUnorderedMapMemoryUsageBytes(this->objectSet.size(),sizeof(NodeID),this->objectSet.bucket_count());
    return memoryUsage/(1024*1024);
}
