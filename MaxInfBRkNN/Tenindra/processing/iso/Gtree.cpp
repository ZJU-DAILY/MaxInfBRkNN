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

#include "Gtree.h"
#include "DynamicGraph.h"

#include "DijkstraSearch.h"
#include "../utility/StopWatch.h"
#include "../utility/utility.h"

#include <assert.h>
#include <iostream>

Gtree::Gtree(std::string networkName, int numNodes, int numEdges, int fanout, std::size_t maxLeafSize): 
    networkName(networkName), numNodes(numNodes), numEdges(numEdges), fanout(fanout), maxLeafSize(maxLeafSize) {
    this->nodeIDLeafVerticesVecIdx.assign(this->numNodes,-1);
    this->nodeIDLeafBordersVecIdx.assign(this->numNodes,-1);
    this->edgeInLeafSubgraph.assign(this->numEdges,true);
    this->nodeLeafIdxs.assign(this->numNodes,-1);
}

Gtree::Gtree() {}

std::string Gtree::getNetworkName()
{
    return this->networkName;
}

int Gtree::getNumEdges()
{
    return this->numEdges;
}

int Gtree::getNumNodes()
{
    return this->numNodes;
}

int Gtree::getFanout()
{
    return this->fanout;
}

int Gtree::getTreeSize()
{
    return this->treeNodes.size();
}

std::size_t Gtree::getMaxLeafSize()
{
     return this->maxLeafSize;
}

int Gtree::getLeafIndex(NodeID nodeID)
{
    return this->nodeLeafIdxs[nodeID];
}

int Gtree::getParentIndex(int treeIdx)
{
    return this->treeNodes[treeIdx].getParentIdx();
}

void Gtree::buildGtree(Graph_RoadNetwork& graph){
    
    this->buildTreeHierarchy(graph);
    this->computeDistanceMatrix(graph);
    this->initialiseGtreeQueryStructure();
    
}

void Gtree::buildTreeHierarchy(Graph_RoadNetwork& graph)
{
    METISIdxToNodeID.resize(this->numNodes);
    METISWrapper metis(this->numNodes,this->numEdges,this->fanout);

    std::unordered_set<NodeID> subgraph = graph.getNodesIDsUset();
    this->addNode(ROOT_PARENT_INDEX,subgraph,graph,metis);
}

void Gtree::addNode(int parentIdx, std::unordered_set<NodeID>& subgraph, Graph_RoadNetwork& originalGraph, METISWrapper& metis)
{
    // Create Gtree node and add to tree
    int treeIdx = this->treeNodes.size(), childvecIdx;
    GtreeNode treeNode(treeIdx,parentIdx,subgraph.size());
    this->treeNodes.push_back(treeNode);

    if (parentIdx != ROOT_PARENT_INDEX) {
        // Add this Gtree node as child of parent (if not root)
        childvecIdx = this->treeNodes[parentIdx].addChild(treeIdx);
        this->treeNodes[treeIdx].setParentChildIdx(childvecIdx);
        
        // Add parent's Gtree path to this nodes Gtree path then add parent
        this->treeNodes[treeIdx].addGtreePathNodes(this->treeNodes[parentIdx].gtreePath);
    }
    this->treeNodes[treeIdx].addGtreePathNode(treeIdx);
    
    // Find and add children of this node
    if (subgraph.size() > this->maxLeafSize) {
        this->addChildren(treeIdx, subgraph, originalGraph,metis);
    } else {
        // If this has less than tau graph node then we stop partitions
        this->treeNodes[treeIdx].setLeafNode();
        this->leafIdxs.push_back(treeIdx);
    }
    
    // Determine and add borders for this Gtree node
    if (parentIdx != ROOT_PARENT_INDEX) {
        // We mark the starting index of this tree nodes borders
        // in the parent's childBordersVec - we can do this as
        // we are about find all these borders and also add 
        // consecutively to the parent's childBordersVec
        this->treeNodes[parentIdx].addChildOffsetInChildBorderVec();
        int adjListStart, adjListSize, leafVerticesVecIdx, leafBordersVecIdx;
        bool isBorder;
        
        for (NodeID node: subgraph) {
            // If this node is a leaf node, we also need to map 
            // all it's graph nodes to they're Gtree leaf index
            if (this->treeNodes[treeIdx].isLeafNode()) {
                // Add each node as a vertice if it is (this
                // will include borders)
                leafVerticesVecIdx = this->treeNodes[treeIdx].addLeafVertex(node);
                this->setLeafVerticesVecIdx(node,leafVerticesVecIdx);
                this->nodeLeafIdxs[node] = treeIdx;
            }
            
            adjListStart = originalGraph.getEdgeListStartIndex(node);
            adjListSize = originalGraph.getEdgeListSize(node);
            isBorder = false;
            for (int i = adjListStart; i < adjListSize; ++i) {
                if (subgraph.find(originalGraph.edges[i].first) == subgraph.end()) {
                    // This mean we have found a neighbour of that outside this subgraph i.e. 
                    // this is node is a border so we add to list of borders of this tree node
                    // and the set of child borders for the parent
                    if (!isBorder) {
                        leafBordersVecIdx = this->treeNodes[treeIdx].addBorder(node);
                        if (this->treeNodes[treeIdx].isLeafNode()) {
                            this->setLeafBordersVecIdx(node,leafBordersVecIdx);
                        }
                        this->treeNodes[parentIdx].addChildBorder(node);
                        isBorder = true; // So that we don't add it again
                    }
                    if (this->treeNodes[treeIdx].isLeafNode()) {
                        // Since this edge leads to a node that is not in the subgraph
                        this->setEdgeNotInSubgraph(i);
                    }
                }
            }
        }

    }
    
}

void Gtree::addChildren(int parentTreeIdx, std::unordered_set<NodeID>& parentGraph, Graph_RoadNetwork& originalGraph, METISWrapper& metis) {
    // Assuming original graph is undirected graph (if not we must make undirected
    // in order for METIS to function)
    
    idx_t n = parentGraph.size();
    
    metis.populateMETISArrays(originalGraph,parentGraph,this->METISIdxToNodeID,true);

    metis.partitionSubgraph(n);    
    
    // Create sets of subgraphs from partitioned parent graph
    // and add them to Gtree (and this will recurse)
    std::vector<std::unordered_set<NodeID>> childGraphs(this->fanout);
    
    for (int i = 0; i < n; ++i) {
        //assert (nodePartitions[i] <= this->fanout && "Invalid partition assigned to graph node, too many partitions");
        // Create child graphs using original NodeID not the METIS idx
        childGraphs[metis.parts[i]].insert(this->METISIdxToNodeID[i]);
    }
    
    // Also add children of this node to Gtree
    for (std::size_t i = 0; i < childGraphs.size(); ++i) {
        this->addNode(parentTreeIdx,childGraphs[i],originalGraph,metis);
        // Note: We release child graph at this point as it is not need again
        // and this would improve overall memory usage (it won't affect vector size)
        utility::releaseSTLCollection(childGraphs[i]);
    }
}

void Gtree::computeDistanceMatrix(Graph_RoadNetwork& graph)
{
    std::vector<std::vector<int>> treeLevelIdxs = this->getTreeNodesByLevel();
    
    DijkstraSearch dijkstra;
    BinaryMinHeap<EdgeWeight,NodeID> *pqueue = new BinaryMinHeap<EdgeWeight,NodeID>();
    int currentIdx;
    DynamicGraph tempGraph(graph);
    std::unordered_set<NodeID>  *targets;
    std::vector<NodeID> adjNodes;
    std::vector<EdgeWeight> adjNodeWgts;
    std::vector<NodeID> *sourcesVec, *targetsVec;
    
    for (int i = treeLevelIdxs.size()-1; i >= 0; --i) {
        // Clear memory in unordered_map or it will continue to grow
        // Note: According to the paper, total number of borders at each level should be O(n)
        // so if we clear this map for each level then we should it's total size should be O(n)
            
        for (std::size_t j = 0; j < treeLevelIdxs[i].size(); ++j) {
            currentIdx = treeLevelIdxs[i][j];
            if (this->treeNodes[currentIdx].isLeafNode()) {
                // In a leaf we find distances from border to all leaf vertices
                sourcesVec = &this->treeNodes[currentIdx].getBorders();
                targets = &this->treeNodes[currentIdx].getLeafVerticesUset();
                targetsVec = &this->treeNodes[currentIdx].getLeafVertices();
            } else {
                // In a non-leaf we get distances from all child borders to all other child borders
                // Note: It's possible some of these distances have already been computed
                sourcesVec = &this->treeNodes[currentIdx].getChildBorders();
                targets = &this->treeNodes[currentIdx].getChildBordersUset();
                targetsVec = &this->treeNodes[currentIdx].getChildBorders();
            }
            
            std::unordered_map<NodeID,EdgeWeight> siblingBorderDistances;
            siblingBorderDistances.reserve(targetsVec->size());
            //Using single std::unordered_map
            int rowLength = targetsVec->size();
            this->treeNodes[currentIdx].matrixRowLength = rowLength;
            this->treeNodes[currentIdx].distanceMatrix.reserve(sourcesVec->size()*rowLength);
            if (this->treeNodes[currentIdx].isLeafNode()) {
                // If it is a leaf node we can search using original Graph
                // who's data structure is faster than DynamicGraph_RoadNetwork
                for (std::size_t i = 0; i < sourcesVec->size(); ++i) {
                    pqueue->clear();
                    dijkstra.findSSMTDistances(graph,(*sourcesVec)[i],(*targets),siblingBorderDistances,pqueue);
                    for (std::size_t j = 0; j < targetsVec->size(); ++j) {
                        this->treeNodes[currentIdx].distanceMatrix.push_back(siblingBorderDistances[(*targetsVec)[j]]);
                    }
                }
            } else {
                for (std::size_t i = 0; i < sourcesVec->size(); ++i) {
                    pqueue->clear();
                    dijkstra.findSSMTDistances(tempGraph,(*sourcesVec)[i],(*targets),siblingBorderDistances,pqueue);
                    for (std::size_t j = 0; j < targetsVec->size(); ++j) {
                        this->treeNodes[currentIdx].distanceMatrix.push_back(siblingBorderDistances[(*targetsVec)[j]]);
                    }
                }
            }
            
            // All future searches will be on this nodes border set (using closure property in paper)
            // because these will be the parent nodes child borders-> Therefore if we need only remove
            // unnecessary edges from these nodes (unimportant nodes will become disconnected)
            //sources = this->treeNodes[currentIdx].getBordersUset();
            sourcesVec = &this->treeNodes[currentIdx].getBorders();

            NodeID border;
            int sourceIdx, targetIdx, distMatrixRowOffset;
            for (std::size_t i = 0; i < sourcesVec->size(); ++i) {
                border = (*sourcesVec)[i];
                if (this->treeNodes[currentIdx].isLeafNode()) {
                    sourceIdx = i;
                    distMatrixRowOffset = sourceIdx*this->treeNodes[currentIdx].matrixRowLength;
                } else {
                    sourceIdx = this->treeNodes[currentIdx].getBorderIdxInChildBorderVec(i);
                    distMatrixRowOffset = sourceIdx*this->treeNodes[currentIdx].matrixRowLength;
                }
                
                // Preserve edges to outside this gtree node (i.e. subgraph)
                adjNodes.clear();
                adjNodeWgts.clear();
                for (std::size_t i = 0; i < tempGraph.nodes[border].adjNodes.size(); ++i) {
                    if (targets->find(tempGraph.nodes[border].adjNodes[i]) == targets->end()) {
                        // This check whether the adj node is within the current gtree node
                        // if it is we do not need to preserve the edge
                        adjNodes.push_back(tempGraph.nodes[border].adjNodes[i]);
                        adjNodeWgts.push_back(tempGraph.nodes[border].adjNodeWgts[i]);
                    }
                }
                // Note: That this will make the graph disconnected but this doesn't
                // matter as removing disconnected node won't change search results and
                // we are only disconnecting nodes that are not borders of the subgraph
                // (i.e. they will not be needed again at higher levels)
                
                tempGraph.nodes[border].adjNodes = std::move(adjNodes);
                tempGraph.nodes[border].adjNodeWgts = std::move(adjNodeWgts);
                for (std::size_t j = 0; j < sourcesVec->size(); ++j) {
                    targetIdx = this->treeNodes[currentIdx].getBorderIdxInChildBorderVec(j); // We will return dist matrix idx whether leaf or not
                    if (border != (*sourcesVec)[j]) {
                        tempGraph.insertImaginaryNonInvertibleEdge(border,(*sourcesVec)[j],this->treeNodes[currentIdx].distanceMatrix[distMatrixRowOffset+targetIdx]);
                    }
                }
            }
        }
    }
    
    delete pqueue;
}

void Gtree::initialiseGtreeQueryStructure()
{
    // Allocate Memory for Gtree Query
    this->sourceToTreeNodeBorderDist.resize(this->treeNodes.size());
    for (size_t i = 0; i < this->treeNodes.size(); ++i) {
        this->sourceToTreeNodeBorderDist[i].resize(this->treeNodes[i].bordersVec.size());
    }
}

void Gtree::printLevels()
{
    std::vector<std::vector<int>> treeNodeLevel = this->getTreeNodesByLevel();

    for (int i = treeNodeLevel.size()-1; i >= 0; --i) {
        std::cout << "Level " << i << ": " << treeNodeLevel[i].size() << std::endl;
    }
}

int Gtree::getNumLevels()
{
    return this->getTreeNodesByLevel().size();
}

int Gtree::getNumBorders()
{
    int numBorders = 0;
    for (std::size_t i = 0; i < this->treeNodes.size(); ++i) {
        numBorders += this->treeNodes[i].getNumBorders();
    }
    return numBorders;
}

int Gtree::getBorderToBorderRelationships()
{
    int numB2BRelationships = 0;
    for (std::size_t i = 0; i < this->treeNodes.size(); ++i) {
        numB2BRelationships += this->treeNodes[i].getNumBorders()*this->treeNodes[i].getNumBorders();
    }
    return numB2BRelationships;
}

int Gtree::getAvgPathCost(int treeIdx)
{
    int avgPathCost = 0;
    if (this->treeNodes[treeIdx].isLeafNode()) {
        avgPathCost = this->treeNodes[treeIdx].getNumBorders();
    } else {
        //assert(this->treeNodes[treeIdx].children.size() > 0); Non-leaf nodes must have at least one child
        int totalChildPathCosts = 0, totalParentChildCost = 0;
        for (std::size_t i = 0; i < this->treeNodes[treeIdx].children.size(); ++i) {
            totalParentChildCost += this->treeNodes[treeIdx].getNumBorders()*this->treeNodes[this->treeNodes[treeIdx].children[i]].getNumBorders();
            totalChildPathCosts += this->getAvgPathCost(this->treeNodes[treeIdx].children[i]);
        }
        avgPathCost += totalParentChildCost/this->treeNodes[treeIdx].children.size();
        avgPathCost += totalChildPathCosts/this->treeNodes[treeIdx].children.size();
    }
    return avgPathCost;
}

std::vector< std::vector< int > > Gtree::getTreeNodesByLevel()
{
    std::vector<std::vector<int>> treeNodeLevel;
    std::vector<int> currentLevel, nextLevel;
    nextLevel.push_back(0);

    while(nextLevel.size() != 0){
            treeNodeLevel.push_back(nextLevel);
            currentLevel.swap(nextLevel);
            nextLevel.clear();
            for (int i: currentLevel){
                for (std::size_t j = 0; j < this->treeNodes[i].children.size(); j++ ){
                    nextLevel.push_back( this->treeNodes[i].children[j] );
                }
            }
    }

    return treeNodeLevel;
}

void Gtree::printNode(int treeIdx)
{
    this->treeNodes[treeIdx].printNode();
}

EdgeWeight Gtree::getShortestPathDistance(Graph_RoadNetwork& graph, NodeID u, NodeID v)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("computations_materialized",0);
    this->stats.initialiseStatistic("computations_executed",0);
    this->stats.initialiseStatistic("computations_total",0);
#endif    
    EdgeWeight spDist = 0;
    int uLeaf = this->getLeafIndex(u);
    int vLeaf = this->getLeafIndex(v);
    
    if (uLeaf == vLeaf) {
        // This mean they are both in the same node
        spDist = this->SPDistLeaf(u,v,uLeaf,graph);
    } else {
        int LCAIdx = 0;
        std::vector<int>& uPathFromRoot = this->treeNodes[uLeaf].getGtreePathFromRoot();
        std::vector<int>& vPathFromRoot = this->treeNodes[vLeaf].getGtreePathFromRoot();

        // Since the two nodes are in different leaf nodes we can guarantee that
        // there is at least one more node (the parent of both) in the G-tree path
        // which we call the the lowest common ancestor (LCA)
        
        // Search down Gtree from root until we find first ancestor that is different
        // this means the previous ancestor is the LCA and the indexes point its children
        unsigned int i, j;
        for (i = 0, j = 0; i < uPathFromRoot.size() && j < vPathFromRoot.size(); ++i, ++j) {
            if (uPathFromRoot[i] != vPathFromRoot[j]) {
                // When i = 0 and j = 0 it is referring to the root so in that
                // case uPathFromRoot[i] does equal vPathFromRoot[j]. This means
                // when they are equal i > 0, so i-1 is safe here
                LCAIdx = uPathFromRoot[i-1];
                break;
            }
            // Note: We can guarantee that LCAIdx will be set here. The only situation
            // it would not be set is both u and v were in the same leaf but we have
            // guarnateed that is not the case in the if/else statement
        }
        
        // Now search up G-tree (from source leaf) until we reach i and j, then we 
        // search down (to target leaf) computing shortest path distance
        
        // From source to source leaf
        this->SPDistToSourceLeafNode(u,uLeaf);
        
        // From source leaf to first child of LCA
        int x = i; // This is safe, depth of tree is O(log(n)) and n is at most 24 million in US dataset
        for (int k = uPathFromRoot.size()-1; k > x; --k) {
            // Since k > x and x is at worst 0, k-1 is safe here
            this->SPDistToParentNode(uPathFromRoot[k],uPathFromRoot[k-1],false);
        }
        
        // From first child of LCA to second child of LCA
        this->SPDistToSiblingNode(uPathFromRoot[i],vPathFromRoot[j],LCAIdx,false);
    
        // From second child of LCA to target leaf
        for (std::size_t k = j; k < vPathFromRoot.size()-1; ++k) {
            // Note the size()-1 in the above condition
            this->SPDistToChildNode(vPathFromRoot[k+1],vPathFromRoot[k],false);
        }
        
        // We assume target has not been visited
        spDist = this->SPDistToLeafTarget(v,vLeaf);
    }
    
    return spDist;
}

EdgeWeight Gtree::SPDistLeaf(NodeID u, NodeID v, int leafNode, Graph_RoadNetwork& graph)
{
    // Assume u and v are in the same leaf node
    EdgeWeight borderToBorderDist, spDist = 0;
    
    spDist = this->DijkstraDist(u,v,leafNode,graph);
    borderToBorderDist = this->BorderDist(u,v,leafNode);
    if (borderToBorderDist < spDist) {
        spDist = borderToBorderDist;
    }
    
    return spDist;
}

std::unordered_map<NodeID,EdgeWeight> Gtree::DijkstraDistMultiTarget(NodeID u, std::unordered_set<NodeID>& targets, int leafNode, Graph_RoadNetwork& graph)
{
    // Assume u and all targets are in the same leaf node and number of targets > 0
    
    DijkstraSearch dijkstra;
    std::unordered_map<NodeID,EdgeWeight> targetDistances;
    
    // This is equivalent of DijkDist function in paper, except we optimise by 
    // doing multi-target search in case there are many objects in source leaf
    dijkstra.findSSMTDistancesSubgraph(graph,u,targets,targetDistances,this->edgeInLeafSubgraph);
    
    return targetDistances;
}

EdgeWeight Gtree::SPDist(NodeID u, NodeID v, std::vector<int>& gtreePath, int firstLCAChild)
{
    // Assume gtreePath size is 2 or more (i.e. u and v are not in same leaf node)
    // Assume that no results have been materialised (this for shortest path query)
    // Note: Path cannot be size 1 since that make one parent of u and v a leaf node
    // which cannot be possible as a leaf node has no children
    
    EdgeWeight spDist = 0;
    
    // Distance from query node to source leaf borders
    std::size_t i = 0;
    this->SPDistToSourceLeafNode(u,gtreePath[i]);
    
    // Distance up the G-tree hierarchy from source leaf
    for (; gtreePath[i] != firstLCAChild; ++i) {
        this->SPDistToParentNode(gtreePath[i], gtreePath[i+1],false);
    }
    
    // Distance between siblings of LCA
    int LCAIdx = this->treeNodes[gtreePath[i]].getParentIdx();
    this->SPDistToSiblingNode(gtreePath[i],gtreePath[i+1],LCAIdx,false);
    ++i;
    
    // Distance down G-tree hierarchy to target leaf
    for (; i+1 < gtreePath.size(); ++i) {
        this->SPDistToChildNode(gtreePath[i+1],gtreePath[i],false);
    }
    
    spDist = this->SPDistToLeafTarget(v,gtreePath[i]);
    
    return spDist;
}

EdgeWeight Gtree::BorderDist(NodeID u, NodeID v, int leafNode)
{
    // Assume that u and v ARE in the same leaf node
    
    EdgeWeight borderDist, minDist = 0;
    
#if defined(GTREE_STL_HASH_TABLE_DIST_MATRIX) || defined(GTREE_GOOGLE_DENSEHASH_DIST_MATRIX)
    if (this->treeNodes[leafNode].bordersVec.size() > 0) {
        // Initialise the minDist as the distance from the source to the first
        // border plus the distance from that border to the target
        minDist = this->distanceMatrix[this->treeNodes[leafNode].bordersVec[0]][u] + 
            this->distanceMatrix[this->treeNodes[leafNode].bordersVec[0]][v];
        for (std::size_t i = 0; i < this->treeNodes[leafNode].bordersVec.size(); ++i) {
            for (std::size_t j = 0; j < this->treeNodes[leafNode].bordersVec.size(); ++j) {
                // ASSUMING UNDIRECTED GRAPH (I.E. BORDER TO LEAF VERTEX DISTANCE = LEAF VERTEX TO BORDER DISTANCE
                // THIS IS BECAUSE WE HAVEN'T CALCULATED LEAF VERTEX TO BORDER DISTANCES
                if (i != j) {
                    // BorderDist return the distance between two nodes using outside vertices
                    // if i == j then we are entering and leaving through the same border
                    // which won't be smaller than the Dijkstra dist so we ignore it
                    borderDist = this->distanceMatrix[this->treeNodes[leafNode].bordersVec[i]][u] + 
                        this->distanceMatrix[this->treeNodes[leafNode].bordersVec[i]][this->treeNodes[leafNode].bordersVec[j]] + 
                        this->distanceMatrix[this->treeNodes[leafNode].bordersVec[j]][v];
                    if (borderDist < minDist) {
                        minDist = borderDist;
                    }
                }
            }
        }
    }    
#else
    int uIdx = this->getIdxInLeafVerticesVec(u);
    int vIdx = this->getIdxInLeafVerticesVec(v);
    int intermediateBorderIdx;
    if (this->treeNodes[leafNode].bordersVec.size() > 0) {
        intermediateBorderIdx = this->treeNodes[leafNode].getBorderIdxInChildBorderVec(0);
        minDist = this->treeNodes[leafNode].distanceMatrix[uIdx] 
            + this->treeNodes[leafNode].distanceMatrix[intermediateBorderIdx] 
            + this->treeNodes[leafNode].distanceMatrix[vIdx];
        for (std::size_t i = 0; i < this->treeNodes[leafNode].bordersVec.size(); ++i) {
            for (std::size_t j = 0; j < this->treeNodes[leafNode].bordersVec.size(); ++j) {
                // ASSUMING UNDIRECTED GRAPH (I.E. BORDER TO LEAF VERTEX DISTANCE = LEAF VERTEX TO BORDER DISTANCE
                // THIS IS BECAUSE WE HAVEN'T CALCULATED LEAF VERTEX TO BORDER DISTANCES
                if (i != j) {
                    // BorderDist return the distance between two nodes using outside vertices
                    // if i == j then we are entering and leaving through the same border
                    // which won't be smaller than the Dijkstra dist so we ignore it
                    intermediateBorderIdx = this->treeNodes[leafNode].getBorderIdxInChildBorderVec(j);
                    borderDist = this->treeNodes[leafNode].distanceMatrix[i*this->treeNodes[leafNode].matrixRowLength+uIdx] 
                        + this->treeNodes[leafNode].distanceMatrix[i*this->treeNodes[leafNode].matrixRowLength+intermediateBorderIdx] 
                        + this->treeNodes[leafNode].distanceMatrix[j*this->treeNodes[leafNode].matrixRowLength+vIdx];
                    if (borderDist < minDist) {
                        minDist = borderDist;
                    }
                }
            }
        }
    }
#endif

#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("computations_executed",this->treeNodes[leafNode].bordersVec.size()*this->treeNodes[leafNode].bordersVec.size());
#endif

    return minDist;
}

EdgeWeight Gtree::DijkstraDist(NodeID u, NodeID v, int leafNode, Graph_RoadNetwork& graph)
{
    // Assume that u and v ARE in the same leaf node
    
    DijkstraSearch dijk;
    
    return dijk.findShortestPathDistanceSubgraph(graph,u,v,this->edgeInLeafSubgraph);
}

// Note: The returned path does not include the LCA
std::vector<int> Gtree::getGtreePath(NodeID u, NodeID v, int& firstLCAChild)
{
    // Note: If u and v are in the same leaf node then path will be 
    // size 0 and LCA will contain leaf node
    std::vector<int> path;
    
    int uLeaf = this->getLeafIndex(u);
    int vLeaf = this->getLeafIndex(v);
    
    std::vector<int>& uPathFromRoot = this->treeNodes[uLeaf].getGtreePathFromRoot();
    std::vector<int>& vPathFromRoot = this->treeNodes[vLeaf].getGtreePathFromRoot();

    // The only case where firstLCAChild is not set is when both nodes are in the same leaf
    // Note: This will be overwritten later if it needs to be
    firstLCAChild = uPathFromRoot.back();
    // Note: Using back() is safe as uPathFromRoot will have at least root in it
    
    // Search down Gtree from root until we find first ancestor that is different
    // this means the previous ancestor is the lowest common ancestor
    //assert(uPathFromRoot[0] == vPathFromRoot[0] == 0 && "First node in path from root is not root!")
    unsigned int i, j;
    for (i = 0, j = 0; i < uPathFromRoot.size() && j < vPathFromRoot.size(); ++i, ++j) {
        if (uPathFromRoot[i] != vPathFromRoot[j]) {
            firstLCAChild = uPathFromRoot[i];
            break;
        }
    }
    
    // Note: That we do not include LCA in return path
    int x = static_cast<int>(i); // This is OK, i is greater than 0 as above
    for (int k = uPathFromRoot.size()-1; k >= x; --k) {
        // i represents the first non-common ancestor
        // so this will be the part of the path
        path.push_back(uPathFromRoot[k]);
    }

    // The path to v is in reverse order (we descend down)
    for (std::size_t k = j; k < vPathFromRoot.size(); ++k) {
        // j represents the first non-common ancestor
        // so this will be the part of the path
        path.push_back(vPathFromRoot[k]);
    }

    return path;

}

void Gtree::getKNNs(OccurenceList& occList, unsigned int k, NodeID queryNode, std::vector<NodeID>& kNNs, 
                    std::vector<EdgeWeight>& kNNDistances, Graph_RoadNetwork& graph)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("computations_materialized",0);
    this->stats.initialiseStatistic("computations_executed",0);
    this->stats.initialiseStatistic("computations_total",0);
#endif
    // Note: We store both G-tree node indexes and road network NodeIDs in the
    // same priority queue. We  map the G-tree node indexes to a NodeID that 
    // is larger than the last node in the graph (and we check whenever 
    // a dequeue happens to determine what type min element is).
    BinaryMinHeap<EdgeWeight,NodeID> pqueue;
    
    int sourceLeaf = this->getLeafIndex(queryNode);
    int Tn, prevTn, treeIdx/*, firstLCAChild*/;
    EdgeWeight Tmin, dist, minKey, borderDist;
    NodeID minNode;
    
    Tn = sourceLeaf;
    Tmin = this->SPDistToSourceLeafNode(queryNode,Tn);
    
    // Note: In order to use a priority queue with primitive types we
    // map the tree indexes as NodeID greater than the largest NodeID
    // for the current graph. Whenever we deque we can tell if the element
    // is a tree node or graph node by check if it's larger than the
    // number of nodes in the graph (since they are number 0 to n-1)
    NodeID firstMappedTreeNodeID = static_cast<NodeID>(this->numNodes);
    
    // If the source leaf node contains any object we must find them first
    if (occList.leafOccurenceList[sourceLeaf].size() > 0) {
#if defined(UNOPTIMISED_GTREE_LEAF_SEARCH)
        std::unordered_map<NodeID,EdgeWeight> leafObjectDistances = this->DijkstraDistMultiTarget(queryNode,occList.leafOccurenceSet[sourceLeaf],sourceLeaf,graph);
        for(NodeID leafObject: occList.leafOccurenceList[sourceLeaf]) {
            borderDist = this->BorderDist(queryNode,leafObject,sourceLeaf);
            if (borderDist < leafObjectDistances[leafObject]) {
                // If a shortest path can be found by leaving a border and re-entering source leaf through
                // another border then we use that distance because it is the actual shortest path distance
                pqueue.insert(leafObject,borderDist);
            } else {
                pqueue.insert(leafObject,leafObjectDistances[leafObject]);
            }
        }
#else
        if (this->getSourceLeafkNNsByINE(queryNode,k,occList.leafOccurenceSet[sourceLeaf],sourceLeaf,graph,kNNs,kNNDistances,pqueue)) {
            // If the in-leaf INE search was able to find definite k objects then we exist now
            return;
        }
#endif        
    }
    
    while (pqueue.size() > 0 || Tn != 0) {
        // UpdateT
        if (pqueue.size() == 0) {
            prevTn = Tn;
            Tn = this->getParentIndex(Tn);
            
            if (Tn != 0) {
                // We don't set Tmin if Tn is the root because this represents
                // infinite instead - it's faster to check Tn != 0 than Tmin < max int
                Tmin = this->SPDistToParentNode(prevTn,Tn);
#if defined(COLLECT_STATISTICS)
                this->stats.incrementStatistic("computations_materialized",this->getComputations(sourceLeaf,prevTn));
#endif
            }
            
            // Tn is a leaf only at start, but pqueue is size 1 at this point
            for(int childIdx: occList.nonLeafOccurenceList[Tn]) {
                if (childIdx != prevTn) {
                    // We only add children that we haven't already visited
                    dist = this->SPDistToSiblingNode(prevTn,childIdx,Tn);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("computations_materialized",this->getComputations(sourceLeaf,prevTn));
#endif
                    pqueue.insert(childIdx+this->numNodes,dist);
                }
            }
            
        }
        
        if (pqueue.size() > 0) {
            // This check is different to paper, but is necessary
            // as above UpdateT will not guarantee pqueue will have
            // size greater than 0 - in the case where the previous
            // Tn parent doesn't have any children with objects
            minKey = pqueue.getMinKey();
            minNode = pqueue.extractMinElement();
            if (Tn != 0 && minKey > Tmin) {
                // Tn == 0 is equivalent to Tmin being infinity
                // So using Tn != 0 we can avoid setting Tmin to infinity (i.e. max int value)
                // which was seen to slow down comparisons

                // UpdateT
                prevTn = Tn;
                Tn = this->getParentIndex(Tn);
                if (Tn != 0) {
                    // We don't set Tmin if Tn is the root because this represents
                    // infinite instead - it's faster to check Tn != 0 than Tmin < max int
                    Tmin = this->SPDistToParentNode(prevTn,Tn);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("computations_materialized",this->getComputations(sourceLeaf,prevTn));
#endif
                }
                
                // Tn is guaranteed to be a non-leaf here because we move to parent
                for(int childIdx: occList.nonLeafOccurenceList[Tn]) {
                    if (childIdx != prevTn) {
                        dist = this->SPDistToSiblingNode(prevTn,childIdx,Tn);
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("computations_materialized",this->getComputations(sourceLeaf,prevTn));
#endif
                        pqueue.insert(childIdx+this->numNodes,dist);
                    }
                }
                
                pqueue.insert(minNode,minKey);
            } else if (minNode < firstMappedTreeNodeID) {
                // Min element is a graph node (i.e. object)
                kNNs.push_back(minNode);
                kNNDistances.push_back(minKey);
                if (kNNs.size() == k) {
                    // If this is the kth nearest neighbour object
                    // then we no longer need to iterate
                    break;
                }
            } else {
                treeIdx = minNode-this->numNodes;
                if (this->treeNodes[treeIdx].isLeafNode()) {
                    for(NodeID leafObject: occList.leafOccurenceList[treeIdx]) {
                        dist = this->SPDistToLeafTarget(leafObject,treeIdx);
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("computations_materialized",this->getComputations(sourceLeaf,treeIdx));
#endif
                        pqueue.insert(leafObject,dist);
                    }
                } else {
                    for(int childIdx: occList.nonLeafOccurenceList[treeIdx]) {
                        dist = this->SPDistToChildNode(childIdx,treeIdx);
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("computations_materialized",this->getComputations(sourceLeaf,treeIdx));
#endif
                        pqueue.insert(childIdx+this->numNodes,dist);
                    }
                }
            }
        }
    }
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("computations_total",this->stats.getStatistic("computations_materialized")+this->stats.getStatistic("computations_executed"));
#endif
}








bool Gtree::getSourceLeafkNNsByINE(NodeID queryNode, unsigned int k, std::unordered_set<NodeID>& targets, int leafNode, Graph_RoadNetwork& graph, std::vector<NodeID>& kNNs,
                                   std::vector<EdgeWeight>& kNNDistances, BinaryMinHeap<EdgeWeight,NodeID>& pqueue)
{
    // Note: We assumpe vectors passed to this method are empty

    BinaryMinHeap<EdgeWeight,NodeID> localQueue;
    EdgeWeight minDist;
    NodeID minDistNodeID, adjNode;
//     std::vector<bool> isNodeSettled(graph.getNumNodes(),false);
//     // We use unordered_set because we are only searching subgraph (no benefit 
//     // from allocating std::vector<bool> for all nodes, just overhead)
    std::unordered_set<NodeID> settledNodeIDs;
    // We can use the fact that we store the leafVerticesVec of each NodeID
    // in the this->nodeIDLeafBordersVecIdx vector to use a vector<bool>
    // to store the visited status of leaf vertices during INE search
//     std::vector<bool> leafVertexVisited(this->treeNodes[leafNode].leafVerticesVec.size(),false);
    int adjListStart, nextAdjListStart;
    
    // Initialize with priority queue with query node ID
    localQueue.insert(queryNode,0);
    
    bool borderEncountered = false;
    unsigned int targetsFound = 0;
    
    // We also exist if all target nodes have been found
    while (localQueue.size() > 0 && targetsFound < targets.size()) {
        // Extract and remove node with smallest distance from query point
        // and mark it as "settled" so we do not inspect again
        minDist = localQueue.getMinKey();
        minDistNodeID = localQueue.extractMinElement();
//         if (!isNodeSettled[minDistNodeID]) {
//             isNodeSettled[minDistNodeID] = 1;
        if (settledNodeIDs.find(minDistNodeID) == settledNodeIDs.end()) {
            settledNodeIDs.insert(minDistNodeID);
//         if (!leafVertexVisited[this->nodeIDLeafVerticesVecIdx[minDistNodeID]]) {
//             leafVertexVisited[this->nodeIDLeafVerticesVecIdx[minDistNodeID]] = true;
            
            if (targets.find(minDistNodeID) != targets.end()) {
                ++targetsFound;
                // If the minimum is an object we have found a kNN
                if (!borderEncountered) {
                    // If we have not encounted a border so far then
                    // all objects we encounter in the source leaf
                    // cannot be bettered by a object outside the source leaf
                    // so we can count these as kNN results
                    kNNs.push_back(minDistNodeID);
                    kNNDistances.push_back(minDist);
                } else {
                    pqueue.insert(minDistNodeID,minDist);
                }
                if (kNNs.size() == k) {
                    return true;
                } else if (targetsFound == k) {
                    // We do not need to add more than k candidates within leaf
                    // because at best these are the actual kNN neighbours and at
                    // worst some of them will be replaced by closer kNNs outside this leaf
                    break;
                }
            }

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistNodeID);
            nextAdjListStart = graph.getEdgeListSize(minDistNodeID);
            // Note: We have to make sure we don't exceed size of graph.edges vector
            
            bool isBorder = false;
            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjNode = graph.edges[i].first;
                if (this->isEdgeInLeafSubgraph(i)) {
//                     if (!isNodeSettled[adjNode]) {
                    if (settledNodeIDs.find(adjNode) == settledNodeIDs.end()) {
//                     if (!leafVertexVisited[this->nodeIDLeafVerticesVecIdx[adjNode]]) {
                        // This is the first time we have seen this node (i.e. distance infinity)
                        localQueue.insert(adjNode,minDist+graph.edges[i].second);
                    }
                } else {
                    if (!isBorder) {
                        isBorder = true;
                    }
                    if (!borderEncountered) {
                        borderEncountered = true;
                    }
                }
            }
            
            // If it is a border we add each of the other borders for this leaf
            // node and the distance from this border (assuming they have not 
            // already been settled - if they are settled we have already found
            // shortest path to them so the path through this border cannot
            // be an improvement)
            if (isBorder) {
//                 int borderIdx = this->treeNodes[leafNode].getBorderIdxInBorderVec(minDistNodeID);
                int borderIdx = this->getIdxInLeafBordersVec(minDistNodeID);
                for (std::size_t i = 0; i < this->treeNodes[leafNode].bordersVec.size(); ++i) {
                    // Note since this is leaf the above functions actual gets idx in leafVerticesVec
                    // Recall that the dimension of the distance matrix is bordersVec.size()*leafVerticesVec.size()
                    // and stores the distances from borders to leaf vertices so we can use these values
                    // compute the distances between two borders in the leaf vertice
//                     if (!isNodeSettled[this->treeNodes[leafNode].bordersVec[i]]) {
                    if (settledNodeIDs.find(this->treeNodes[leafNode].bordersVec[i]) == settledNodeIDs.end()) {
//                     if (!leafVertexVisited[this->nodeIDLeafVerticesVecIdx[this->treeNodes[leafNode].bordersVec[i]]]) {
                        int targetBorderIdx = this->treeNodes[leafNode].getBorderIdxInChildBorderVec(i);
#if defined(GTREE_STL_HASH_TABLE_DIST_MATRIX) || defined(GTREE_GOOGLE_DENSEHASH_DIST_MATRIX)
                        localQueue.insert(this->treeNodes[leafNode].bordersVec[i],minDist+this->distanceMatrix[minDistNodeID][this->treeNodes[leafNode].bordersVec[i]]);
#else
                        localQueue.insert(this->treeNodes[leafNode].bordersVec[i],minDist+this->treeNodes[leafNode].distanceMatrix[borderIdx*this->treeNodes[leafNode].matrixRowLength+targetBorderIdx]);
#endif
                    }
                }
#if defined(COLLECT_STATISTICS)
                this->stats.incrementStatistic("computations_executed",this->treeNodes[leafNode].bordersVec.size());
#endif
            }
        }
    }
    return false;
}

EdgeWeight Gtree::SPDistToSourceLeafNode(NodeID u, int sourceLeafIdx)
{
    EdgeWeight spDist = 0, sourceToNextBorderDist;
#if defined(GTREE_STL_HASH_TABLE_DIST_MATRIX) || defined(GTREE_GOOGLE_DENSEHASH_DIST_MATRIX)
    if (this->treeNodes[sourceLeafIdx].bordersVec.size() > 0) {
        spDist = this->distanceMatrix[this->treeNodes[sourceLeafIdx].bordersVec[0]][u];
        this->sourceToTreeNodeBorderDist[sourceLeafIdx][0] = spDist;
        for (std::size_t i = 1; i < this->treeNodes[sourceLeafIdx].bordersVec.size(); ++i) {
            sourceToNextBorderDist = this->distanceMatrix[this->treeNodes[sourceLeafIdx].bordersVec[i]][u];
            // ASSUMING DIRECTED GRAPH (I.E. BORDER TO LEAF VERTEX DISTANCE = LEAF VERTEX TO BORDER DISTANCE
            // THIS IS BECAUSE WE HAVEN'T CALCULATED LEAF VERTEX TO BORDER DISTANCES
            this->sourceToTreeNodeBorderDist[sourceLeafIdx][i] = sourceToNextBorderDist;
            if (sourceToNextBorderDist < spDist) {
                spDist = sourceToNextBorderDist;
            }
        }
    }
#else
    int uIdx = this->getIdxInLeafVerticesVec(u);
    if (this->treeNodes[sourceLeafIdx].bordersVec.size() > 0) {
        spDist = this->treeNodes[sourceLeafIdx].distanceMatrix[uIdx];
        this->sourceToTreeNodeBorderDist[sourceLeafIdx][0] = spDist;
        for (std::size_t i = 1; i < this->treeNodes[sourceLeafIdx].bordersVec.size(); ++i) {
            sourceToNextBorderDist = this->treeNodes[sourceLeafIdx].distanceMatrix[i*this->treeNodes[sourceLeafIdx].matrixRowLength+uIdx];
            // ASSUMING DIRECTED GRAPH (I.E. BORDER TO LEAF VERTEX DISTANCE = LEAF VERTEX TO BORDER DISTANCE
            // THIS IS BECAUSE WE HAVEN'T CALCULATED LEAF VERTEX TO BORDER DISTANCES
            this->sourceToTreeNodeBorderDist[sourceLeafIdx][i] = sourceToNextBorderDist;
            if (sourceToNextBorderDist < spDist) {
                spDist = sourceToNextBorderDist;
            }
        }
    }
#endif
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("computations_executed",this->treeNodes[sourceLeafIdx].bordersVec.size());
#endif
    return spDist;
}

EdgeWeight Gtree::SPDistToParentNode(int childTreeIdx, int parentTreeIdx, bool computeSPDist)
{
    EdgeWeight spDist = 0, sourceToNextBorderDist, sourceToChildBorderDist, sourceToParentBorderDist;
    
    // Note: We assume we have already calculated distance to childTreeIdx from source
    
#if defined(GTREE_STL_HASH_TABLE_DIST_MATRIX) || defined(GTREE_GOOGLE_DENSEHASH_DIST_MATRIX)
    if (this->treeNodes[childTreeIdx].bordersVec.size() > 0) {
        sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[childTreeIdx][0];
        for (std::size_t k = 0; k < this->treeNodes[parentTreeIdx].bordersVec.size(); ++k) {
            sourceToParentBorderDist = sourceToChildBorderDist + 
                this->distanceMatrix[this->treeNodes[childTreeIdx].bordersVec[0]][this->treeNodes[parentTreeIdx].bordersVec[k]];
            // I.e. this is the first time we are computing distances to parent's border set so we initialise
            this->sourceToTreeNodeBorderDist[parentTreeIdx][k] = sourceToParentBorderDist;
        }
        for (std::size_t j = 1; j < this->treeNodes[childTreeIdx].bordersVec.size(); ++j) {
            sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[childTreeIdx][j];
            for (std::size_t k = 0; k < this->treeNodes[parentTreeIdx].bordersVec.size(); ++k) {
                sourceToParentBorderDist = sourceToChildBorderDist + 
                    this->distanceMatrix[this->treeNodes[childTreeIdx].bordersVec[j]][this->treeNodes[parentTreeIdx].bordersVec[k]];
                if (sourceToParentBorderDist < this->sourceToTreeNodeBorderDist[parentTreeIdx][k]) {
                    this->sourceToTreeNodeBorderDist[parentTreeIdx][k] = sourceToParentBorderDist;
                }
            }
        }
    }
#else
    int childPos, parentBorderIdx, childBorderOffset, childBorderIdx;

    childPos = this->treeNodes[childTreeIdx].getParentChildIdx();
    childBorderOffset = this->treeNodes[parentTreeIdx].getChildOffsetInChildBorderVec(childPos);
    
    if (this->treeNodes[childTreeIdx].bordersVec.size() > 0) {
        childBorderIdx = childBorderOffset;
        sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[childTreeIdx][0];
        for (std::size_t k = 0; k < this->treeNodes[parentTreeIdx].bordersVec.size(); ++k) {
            parentBorderIdx = this->treeNodes[parentTreeIdx].getBorderIdxInChildBorderVec(k);
            sourceToParentBorderDist = sourceToChildBorderDist + this->treeNodes[parentTreeIdx].distanceMatrix[childBorderIdx*this->treeNodes[parentTreeIdx].matrixRowLength+parentBorderIdx];
            // I.e. this is the first time we are computing distances to parent's border set so we initialise
            this->sourceToTreeNodeBorderDist[parentTreeIdx][k] = sourceToParentBorderDist;
        }
        for (std::size_t j = 1; j < this->treeNodes[childTreeIdx].bordersVec.size(); ++j) {
            childBorderIdx = childBorderOffset + j;
            sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[childTreeIdx][j];
            for (std::size_t k = 0; k < this->treeNodes[parentTreeIdx].bordersVec.size(); ++k) {
                parentBorderIdx = this->treeNodes[parentTreeIdx].getBorderIdxInChildBorderVec(k);
                sourceToParentBorderDist = sourceToChildBorderDist + this->treeNodes[parentTreeIdx].distanceMatrix[childBorderIdx*this->treeNodes[parentTreeIdx].matrixRowLength+parentBorderIdx];
                if (sourceToParentBorderDist < this->sourceToTreeNodeBorderDist[parentTreeIdx][k]) {
                    this->sourceToTreeNodeBorderDist[parentTreeIdx][k] = sourceToParentBorderDist;
                }
            }
        }
    }
#endif
    
    if (this->sourceToTreeNodeBorderDist[parentTreeIdx].size() > 0 && computeSPDist) {
        spDist = this->sourceToTreeNodeBorderDist[parentTreeIdx][0];
        for (std::size_t i = 1; i < this->sourceToTreeNodeBorderDist[parentTreeIdx].size(); ++i) {
            sourceToNextBorderDist = this->sourceToTreeNodeBorderDist[parentTreeIdx][i];
            if (sourceToNextBorderDist < spDist) {
                spDist = sourceToNextBorderDist;
            }
        }
    }

#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("computations_executed",this->treeNodes[childTreeIdx].bordersVec.size()*this->treeNodes[parentTreeIdx].bordersVec.size());
#endif
    return spDist;
}

EdgeWeight Gtree::SPDistToSiblingNode(int firstLCAChildIdx, int targetLCAChildIdx, int LCAIdx, bool computeSPDist)
{
    EdgeWeight spDist = 0, sourceToNextBorderDist, sourceToChildBorderDist, sourceToBorderDist;

#if defined(GTREE_STL_HASH_TABLE_DIST_MATRIX) || defined(GTREE_GOOGLE_DENSEHASH_DIST_MATRIX)
    if (this->treeNodes[firstLCAChildIdx].bordersVec.size() > 0) {
        sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[firstLCAChildIdx][0];
        for (std::size_t k = 0; k < this->treeNodes[targetLCAChildIdx].bordersVec.size(); ++k) {
            sourceToBorderDist = sourceToChildBorderDist + 
                this->distanceMatrix[this->treeNodes[firstLCAChildIdx].bordersVec[0]][this->treeNodes[targetLCAChildIdx].bordersVec[k]];
            // I.e. this is the first time we are computing distances to tthe target set
            this->sourceToTreeNodeBorderDist[targetLCAChildIdx][k] = sourceToBorderDist;
        }
        for (std::size_t j = 1; j < this->treeNodes[firstLCAChildIdx].bordersVec.size(); ++j) {
            sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[firstLCAChildIdx][j];
            for (std::size_t k = 0; k < this->treeNodes[targetLCAChildIdx].bordersVec.size(); ++k) {
                sourceToBorderDist = sourceToChildBorderDist + 
                    this->distanceMatrix[this->treeNodes[firstLCAChildIdx].bordersVec[j]][this->treeNodes[targetLCAChildIdx].bordersVec[k]];
                if (sourceToBorderDist < this->sourceToTreeNodeBorderDist[targetLCAChildIdx][k]) {
                    this->sourceToTreeNodeBorderDist[targetLCAChildIdx][k] = sourceToBorderDist;
                }
            }
        }
    }
#else
    int childPos, childBorderOffset, childBorderIdx;
    int targetChildPos, targetBorderOffset, targetBorderIdx;
    
    childPos = this->treeNodes[firstLCAChildIdx].getParentChildIdx();
    childBorderOffset = this->treeNodes[LCAIdx].getChildOffsetInChildBorderVec(childPos);
    targetChildPos = this->treeNodes[targetLCAChildIdx].getParentChildIdx();
    targetBorderOffset = this->treeNodes[LCAIdx].getChildOffsetInChildBorderVec(targetChildPos);
    
    if (this->treeNodes[firstLCAChildIdx].bordersVec.size() > 0) {
        childBorderIdx = childBorderOffset;
        sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[firstLCAChildIdx][0];
        for (std::size_t k = 0; k < this->treeNodes[targetLCAChildIdx].bordersVec.size(); ++k) {
            targetBorderIdx = targetBorderOffset + k;
            sourceToBorderDist = sourceToChildBorderDist + this->treeNodes[LCAIdx].distanceMatrix[childBorderIdx*this->treeNodes[LCAIdx].matrixRowLength+targetBorderIdx];
            // I.e. this is the first time we are computing distances to tthe target set
            this->sourceToTreeNodeBorderDist[targetLCAChildIdx][k] = sourceToBorderDist;
        }
        for (std::size_t j = 1; j < this->treeNodes[firstLCAChildIdx].bordersVec.size(); ++j) {
            childBorderIdx = childBorderOffset + j;
            sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[firstLCAChildIdx][j];
            for (std::size_t k = 0; k < this->treeNodes[targetLCAChildIdx].bordersVec.size(); ++k) {
                targetBorderIdx = targetBorderOffset + k;
                sourceToBorderDist = sourceToChildBorderDist + this->treeNodes[LCAIdx].distanceMatrix[childBorderIdx*this->treeNodes[LCAIdx].matrixRowLength+targetBorderIdx];
                if (sourceToBorderDist < this->sourceToTreeNodeBorderDist[targetLCAChildIdx][k]) {
                    this->sourceToTreeNodeBorderDist[targetLCAChildIdx][k] = sourceToBorderDist;
                }
            }
        }
    }
#endif

    if (this->sourceToTreeNodeBorderDist[targetLCAChildIdx].size() > 0 && computeSPDist) {
        spDist = this->sourceToTreeNodeBorderDist[targetLCAChildIdx][0];
        for (std::size_t i = 1; i < this->sourceToTreeNodeBorderDist[targetLCAChildIdx].size(); ++i) {
            sourceToNextBorderDist = this->sourceToTreeNodeBorderDist[targetLCAChildIdx][i];
            if (sourceToNextBorderDist < spDist) {
                spDist = sourceToNextBorderDist;
            }
        }
    }
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("computations_executed",this->treeNodes[firstLCAChildIdx].bordersVec.size()*this->treeNodes[targetLCAChildIdx].bordersVec.size());
#endif

    return spDist;
}

EdgeWeight Gtree::SPDistToChildNode(int childTreeIdx, int parentTreeIdx, bool computeSPDist)
{
    EdgeWeight spDist = 0, sourceToNextBorderDist, sourceToChildBorderDist, sourceToBorderDist;
    
#if defined(GTREE_STL_HASH_TABLE_DIST_MATRIX) || defined(GTREE_GOOGLE_DENSEHASH_DIST_MATRIX)
    if (this->treeNodes[parentTreeIdx].bordersVec.size() > 0) {
        sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[parentTreeIdx][0];
        for (std::size_t k = 0; k < this->treeNodes[childTreeIdx].bordersVec.size(); ++k) {
            sourceToBorderDist = sourceToChildBorderDist + 
                this->distanceMatrix[this->treeNodes[parentTreeIdx].bordersVec[0]][this->treeNodes[childTreeIdx].bordersVec[k]];
            // I.e. this is the first time we are computing distances to tthe target set
            this->sourceToTreeNodeBorderDist[childTreeIdx][k] = sourceToBorderDist;
        }
        for (std::size_t j = 1; j < this->treeNodes[parentTreeIdx].bordersVec.size(); ++j) {
            sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[parentTreeIdx][j];
            for (std::size_t k = 0; k < this->treeNodes[childTreeIdx].bordersVec.size(); ++k) {
                sourceToBorderDist = sourceToChildBorderDist + 
                    this->distanceMatrix[this->treeNodes[parentTreeIdx].bordersVec[j]][this->treeNodes[childTreeIdx].bordersVec[k]];
                if (sourceToBorderDist < this->sourceToTreeNodeBorderDist[childTreeIdx][k]) {
                    this->sourceToTreeNodeBorderDist[childTreeIdx][k] = sourceToBorderDist;
                }
            }
        }
    }
#else
    int childPos, parentBorderIdx, childBorderOffset, childBorderIdx;
    childPos = this->treeNodes[childTreeIdx].getParentChildIdx();
    childBorderOffset = this->treeNodes[parentTreeIdx].getChildOffsetInChildBorderVec(childPos);

    if (this->treeNodes[parentTreeIdx].bordersVec.size() > 0) {
        parentBorderIdx = this->treeNodes[parentTreeIdx].getBorderIdxInChildBorderVec(0);
        sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[parentTreeIdx][0];
        for (std::size_t k = 0; k < this->treeNodes[childTreeIdx].bordersVec.size(); ++k) {
            childBorderIdx = childBorderOffset + k;
            sourceToBorderDist = sourceToChildBorderDist + this->treeNodes[parentTreeIdx].distanceMatrix[parentBorderIdx*this->treeNodes[parentTreeIdx].matrixRowLength+childBorderIdx];
            // I.e. this is the first time we are computing distances to tthe target set
            this->sourceToTreeNodeBorderDist[childTreeIdx][k] = sourceToBorderDist;
        }
        for (std::size_t j = 1; j < this->treeNodes[parentTreeIdx].bordersVec.size(); ++j) {
            parentBorderIdx = this->treeNodes[parentTreeIdx].getBorderIdxInChildBorderVec(j);
            sourceToChildBorderDist = this->sourceToTreeNodeBorderDist[parentTreeIdx][j];
            for (std::size_t k = 0; k < this->treeNodes[childTreeIdx].bordersVec.size(); ++k) {
                childBorderIdx = childBorderOffset + k;
                sourceToBorderDist = sourceToChildBorderDist + this->treeNodes[parentTreeIdx].distanceMatrix[parentBorderIdx*this->treeNodes[parentTreeIdx].matrixRowLength+childBorderIdx];
                if (sourceToBorderDist < this->sourceToTreeNodeBorderDist[childTreeIdx][k]) {
                    this->sourceToTreeNodeBorderDist[childTreeIdx][k] = sourceToBorderDist;
                }
            }
        }
    }
#endif

    if (this->sourceToTreeNodeBorderDist[childTreeIdx].size() > 0 && computeSPDist) {
        spDist = this->sourceToTreeNodeBorderDist[childTreeIdx][0];
        for (std::size_t i = 1; i < this->sourceToTreeNodeBorderDist[childTreeIdx].size(); ++i) {
            sourceToNextBorderDist = this->sourceToTreeNodeBorderDist[childTreeIdx][i];
            if (sourceToNextBorderDist < spDist) {
                spDist = sourceToNextBorderDist;
            }
        }
    }
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("computations_executed",this->treeNodes[parentTreeIdx].bordersVec.size()*this->treeNodes[childTreeIdx].bordersVec.size());
#endif
    
    return spDist;
    
}

EdgeWeight Gtree::SPDistToLeafTarget(NodeID target, int leafIdx)
{
    EdgeWeight spDist = 0, candidateSourceToTargetDist;
    
    // Now find the border which has the shortest distance to the target in the leaf node
#if defined(GTREE_STL_HASH_TABLE_DIST_MATRIX) || defined(GTREE_GOOGLE_DENSEHASH_DIST_MATRIX)
    if (this->treeNodes[leafIdx].bordersVec.size() > 0) {
        spDist = this->sourceToTreeNodeBorderDist[leafIdx][0] + this->distanceMatrix[this->treeNodes[leafIdx].bordersVec[0]][target];
        for (std::size_t i = 0; i < this->treeNodes[leafIdx].bordersVec.size(); ++i) {
            candidateSourceToTargetDist = this->sourceToTreeNodeBorderDist[leafIdx][i] + this->distanceMatrix[this->treeNodes[leafIdx].bordersVec[i]][target];
            if (candidateSourceToTargetDist < spDist) {
                spDist = candidateSourceToTargetDist;
            }
        }
    }
#else
    int vIdx = this->getIdxInLeafVerticesVec(target);

    if (this->treeNodes[leafIdx].bordersVec.size() > 0) {
        spDist = this->sourceToTreeNodeBorderDist[leafIdx][0] + this->treeNodes[leafIdx].distanceMatrix[vIdx];
        for (std::size_t i = 0; i < this->treeNodes[leafIdx].bordersVec.size(); ++i) {
            candidateSourceToTargetDist = this->sourceToTreeNodeBorderDist[leafIdx][i] + this->treeNodes[leafIdx].distanceMatrix[i*this->treeNodes[leafIdx].matrixRowLength+vIdx];
            if (candidateSourceToTargetDist < spDist) {
                spDist = candidateSourceToTargetDist;
            }
        }
    }
#endif
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("computations_executed",this->treeNodes[leafIdx].bordersVec.size());
#endif
    
    return spDist;
}

int Gtree::getIdxInLeafVerticesVec(NodeID u)
{
    //assert(this->nodeIDLeafVerticesVecIdx[u] != -1 && "This node was never added to a leaf Gtree node");
    return this->nodeIDLeafVerticesVecIdx[u];
}

void Gtree::setLeafVerticesVecIdx(NodeID u, int leafVerticesVecIdx)
{
    this->nodeIDLeafVerticesVecIdx[u] = leafVerticesVecIdx;
}

int Gtree::getIdxInLeafBordersVec(NodeID u)
{
    return this->nodeIDLeafBordersVecIdx[u];
}

void Gtree::setLeafBordersVecIdx(NodeID u, int leafBordersVecIdx)
{
    this->nodeIDLeafBordersVecIdx[u] = leafBordersVecIdx;
}

bool Gtree::isEdgeInLeafSubgraph(NodeID edge)
{
    return this->edgeInLeafSubgraph[edge];
}

void Gtree::setEdgeNotInSubgraph(NodeID edge)
{
    this->edgeInLeafSubgraph[edge] = false;
}

double Gtree::computeIndexSize()
{
    double memoryUsage = 0;
    for (std::size_t i = 0; i < this->treeNodes.size(); ++i) {
        memoryUsage += this->treeNodes[i].computeIndexSizeBytes();
    }
    memoryUsage += sizeof(int)*this->leafIdxs.size();
    memoryUsage += sizeof(std::vector<EdgeWeight>)*this->sourceToTreeNodeBorderDist.size();
    for (std::size_t i = 0; i < sourceToTreeNodeBorderDist.size(); ++i) {
        memoryUsage += sizeof(EdgeWeight)*this->sourceToTreeNodeBorderDist[i].size();
    }
    memoryUsage += sizeof(int)*this->nodeIDLeafVerticesVecIdx.size();
    memoryUsage += sizeof(int)*this->nodeIDLeafBordersVecIdx.size();
    memoryUsage += this->edgeInLeafSubgraph.size()/8; // std::vector<bool> only use 1 bit per element
    memoryUsage += sizeof(int)*this->nodeLeafIdxs.size();
    return memoryUsage/(1024*1024);
}

double Gtree::computeMemoryUsage()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    for (std::size_t i = 0; i < this->treeNodes.size(); ++i) {
        memoryUsage += this->treeNodes[i].computeMemoryUsageBytes();
    }
    memoryUsage += sizeof(GtreeNode)*(this->treeNodes.capacity()-this->treeNodes.size());
    memoryUsage += sizeof(int)*this->leafIdxs.capacity();
    memoryUsage += sizeof(std::vector<EdgeWeight>)*this->sourceToTreeNodeBorderDist.capacity();
    for (std::size_t i = 0; i < sourceToTreeNodeBorderDist.size(); ++i) {
        memoryUsage += sizeof(EdgeWeight)*this->sourceToTreeNodeBorderDist[i].capacity();
    }
    memoryUsage += sizeof(int)*this->nodeIDLeafVerticesVecIdx.capacity();
    memoryUsage += sizeof(int)*this->nodeIDLeafBordersVecIdx.capacity();
    memoryUsage += this->edgeInLeafSubgraph.capacity()/8; // std::vector<bool> only use 1 bit per element
    memoryUsage += sizeof(int)*this->nodeLeafIdxs.capacity();
    memoryUsage += sizeof(NodeID)*this->METISIdxToNodeID.capacity();
    memoryUsage += this->networkName.size();
#if defined(GTREE_STL_HASH_TABLE_DIST_MATRIX) || defined(GTREE_GOOGLE_DENSEHASH_DIST_MATRIX)
    memoryUsage += utility::estimateUnorderedMapMemoryUsageBytes(this->distanceMatrix.size(),sizeof(std::pair<NodeID,std::unordered_map<NodeID,EdgeWeight>>),this->distanceMatrix.bucket_count());
    for (auto it = this->distanceMatrix.begin(); it != this->distanceMatrix.end(); ++it) {
        memoryUsage += sizeof(it->second)+utility::estimateUnorderedMapMemoryUsageBytes(it->second.size(),sizeof(std::pair<NodeID,EdgeWeight>),it->second.bucket_count());
    }
#endif
    return memoryUsage/(1024*1024);
}

double Gtree::computeDistanceMatrixMemoryUsage()
{
    double memoryUsage = 0;
#if defined(GTREE_STL_HASH_TABLE_DIST_MATRIX) || defined(GTREE_GOOGLE_DENSEHASH_DIST_MATRIX)
    memoryUsage += utility::estimateUnorderedMapMemoryUsageBytes(this->distanceMatrix.size(),sizeof(std::pair<NodeID,std::unordered_map<NodeID,EdgeWeight>>),this->distanceMatrix.bucket_count());
    for (auto it = this->distanceMatrix.begin(); it != this->distanceMatrix.end(); ++it) {
        memoryUsage += sizeof(it->second)+utility::estimateUnorderedMapMemoryUsageBytes(it->second.size(),sizeof(std::pair<NodeID,EdgeWeight>),it->second.bucket_count());
    }
#else
    for (std::size_t i = 0; i < this->treeNodes.size(); ++i) {
        memoryUsage += this->treeNodes[i].computeDistanceMatrixMemoryUsageBytes();
    }
#endif
    return memoryUsage/(1024*1024);
}

EdgeWeight Gtree::getRepeatedShortestPathDistance(Graph_RoadNetwork& graph, NodeID u, NodeID v, std::vector<bool>& visited, std::vector<EdgeWeight>& leafVertexDistances,
                                                             BinaryMinHeap<EdgeWeight,NodeID>& pqueue, std::unordered_set<NodeID>& leafVertexVisited)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("computations_materialized",0);
    this->stats.initialiseStatistic("computations_executed",0);
    this->stats.initialiseStatistic("computations_total",0);
#endif
    EdgeWeight spDist = 0;
    int uLeaf = this->getLeafIndex(u);
    int vLeaf = this->getLeafIndex(v);
    
    if (uLeaf == vLeaf) {
        // This mean they are both in the same node
//         spDist = this->SPDistLeaf(u,v,uLeaf,graph);
        if (leafVertexDistances.size() == 0) {
            // This means repeated leaf search has not be initialised yet
            leafVertexDistances.resize(this->treeNodes[uLeaf].leafVerticesVec.size(),0);
//             leafVertexVisited.resize(this->treeNodes[uLeaf].leafVerticesVec.size(),false);
            pqueue.insert(u,0); // Insert source node into queue for leaf search
        }
        spDist = this->getRepeatedSourceLeafShortestPathDistance(graph,uLeaf,v,leafVertexDistances,pqueue,leafVertexVisited);
    } else {
        int LCAIdx = 0;
        std::vector<int>& uPathFromRoot = this->treeNodes[uLeaf].getGtreePathFromRoot();
        std::vector<int>& vPathFromRoot = this->treeNodes[vLeaf].getGtreePathFromRoot();

        // Since the two nodes are in different leaf nodes we can guarantee that
        // there is at least one more node (the parent of both) in the G-tree path
        // which we call the the lowest common ancestor (LCA)
        
        // Search down Gtree from root until we find first ancestor that is different
        // this means the previous ancestor is the LCA and the indexes point its children
        unsigned int i, j;
        for (i = 0, j = 0; i < uPathFromRoot.size() && j < vPathFromRoot.size(); ++i, ++j) {
            if (uPathFromRoot[i] != vPathFromRoot[j]) {
                // When i = 0 and j = 0 it is referring to the root so in that
                // case uPathFromRoot[i] does equal vPathFromRoot[j]. This means
                // when they are equal i > 0, so i-1 is safe here
                LCAIdx = uPathFromRoot[i-1];
                break;
            }
            // Note: We can guarantee that LCAIdx will be set here. The only situation
            // it would not be set is both u and v were in the same leaf but we have
            // guarnateed that is not the case in the if/else statement
        }
        
        // Now search up G-tree (from source leaf) until we reach i and j, then we 
        // search down (to target leaf) computing shortest path distance
        
        // From source to source leaf
        if (!visited[uLeaf]) {
            this->SPDistToSourceLeafNode(u,uLeaf);
            visited[uLeaf] = true;
        }
        
        // From source leaf to first child of LCA
        int x = i; // This is safe, depth of tree is O(log(n)) and n is at most 24 million in US dataset
        for (int k = uPathFromRoot.size()-1; k > x; --k) {
            // Since k > x and x is at worst 0, k-1 is safe here
            if (!visited[uPathFromRoot[k-1]]) {
                this->SPDistToParentNode(uPathFromRoot[k],uPathFromRoot[k-1],false);
                visited[uPathFromRoot[k-1]] = true;
            }
        }
        
        // From first child of LCA to second child of LCA
        if (!visited[vPathFromRoot[j]]) {
            this->SPDistToSiblingNode(uPathFromRoot[i],vPathFromRoot[j],LCAIdx,false);
            visited[vPathFromRoot[j]] = true;
        }

        // From second child of LCA to target leaf
        for (std::size_t k = j; k < vPathFromRoot.size()-1; ++k) {
            // Note the size()-1 in the above condition
            if (!visited[vPathFromRoot[k+1]]) {
                this->SPDistToChildNode(vPathFromRoot[k+1],vPathFromRoot[k],false);
                visited[vPathFromRoot[k+1]] = true;
            }
        }
        
        // We assume target has not been visited
        spDist = this->SPDistToLeafTarget(v,vLeaf);
    }
    
    return spDist;
}

EdgeWeight Gtree::getRepeatedSourceLeafShortestPathDistance(Graph_RoadNetwork& graph, int leafNode, NodeID v, std::vector<EdgeWeight>& leafVertexDistances,
                                                            BinaryMinHeap<EdgeWeight,NodeID>& pqueue, std::unordered_set<NodeID>& leafVertexVisited)
{
    // Assume vectors have been resized for number of vertices in source leaf
    // and query vertex has been inserted into pqueue
    //assert(this->getLeafIndex()[v] != leafNode && "Target is not in source leaf node!");
    
    if (leafVertexVisited.find(v) != leafVertexVisited.end()) {
        // If it's been already visited in the Dijkstra's search
        // we can return the distance immediately
        return leafVertexDistances[this->nodeIDLeafVerticesVecIdx[v]];
    }
    
    EdgeWeight minDist;
    NodeID minDistNodeID, adjNode;
    int adjListStart, nextAdjListStart;
    
    // We also exist if all target nodes have been found
    while (pqueue.size() > 0) {
        // Extract and remove node with smallest distance from query point
        // and mark it as "settled" so we do not inspect again
        minDist = pqueue.getMinKey();
        minDistNodeID = pqueue.extractMinElement();
        if (leafVertexVisited.find(minDistNodeID) == leafVertexVisited.end()) {
            if (minDistNodeID == v) {
                // Re-insert it back into queue so that we can dequeue
                // and add neighbours in repeated search
                pqueue.insert(minDistNodeID,minDist);
                return minDist;
            }

            leafVertexVisited.insert(minDistNodeID);
            leafVertexDistances[this->nodeIDLeafVerticesVecIdx[minDistNodeID]] = minDist;

            // Inspect each neighbour and update pqueue using edge weights
            adjListStart = graph.getEdgeListStartIndex(minDistNodeID);
            nextAdjListStart = graph.getEdgeListSize(minDistNodeID);
            // Note: We have to make sure we don't exceed size of graph.edges vector
            
            bool isBorder = false;
            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjNode = graph.edges[i].first;
                if (this->isEdgeInLeafSubgraph(i)) {
                    if (leafVertexVisited.find(adjNode) == leafVertexVisited.end()) {
                        // This is the first time we have seen this node (i.e. distance infinity)
                        pqueue.insert(adjNode,minDist+graph.edges[i].second);
                    }
                } else {
                    if (!isBorder) {
                        isBorder = true;
                    }
                }
            }
            
            // If it is a border we add each of the other borders for this leaf
            // node and the distance from this border (assuming they have not 
            // already been settled - if they are settled we have already found
            // shortest path to them so the path through this border cannot
            // be an improvement)
            if (isBorder) {
                int borderIdx = this->getIdxInLeafBordersVec(minDistNodeID);
                for (std::size_t i = 0; i < this->treeNodes[leafNode].bordersVec.size(); ++i) {
                    // Note since this is leaf the above functions actual gets idx in leafVerticesVec
                    // Recall that the dimension of the distance matrix is bordersVec.size()*leafVerticesVec.size()
                    // and stores the distances from borders to leaf vertices so we can use these values
                    // compute the distances between two borders in the leaf vertice
                    if (leafVertexVisited.find(this->treeNodes[leafNode].bordersVec[i]) == leafVertexVisited.end()) {
                        int targetBorderIdx = this->treeNodes[leafNode].getBorderIdxInChildBorderVec(i);
#if defined(GTREE_STL_HASH_TABLE_DIST_MATRIX) || defined(GTREE_GOOGLE_DENSEHASH_DIST_MATRIX)
                        pqueue.insert(this->treeNodes[leafNode].bordersVec[i],minDist+this->distanceMatrix[minDistNodeID][this->treeNodes[leafNode].bordersVec[i]]);
#else
                        pqueue.insert(this->treeNodes[leafNode].bordersVec[i],minDist+this->treeNodes[leafNode].distanceMatrix[borderIdx*this->treeNodes[leafNode].matrixRowLength+targetBorderIdx]);
#endif
                    }
                }
#if defined(COLLECT_STATISTICS)
                this->stats.incrementStatistic("computations_executed",this->treeNodes[leafNode].bordersVec.size());
#endif
            }
            
        }
    }
    assert(false && "Could not find target in G-tree source leaf node!");
    return 0;
}


void Gtree::populateUnorderedMapDistanceMatrix()
{
#if defined(GTREE_GOOGLE_DENSEHASH_DIST_MATRIX)
    this->distanceMatrix.set_empty_key(constants::UNUSED_NODE_ID); // Required by google::dense_hash_map
#endif
    
    std::vector<std::vector<int>> treeNodeLevel = this->getTreeNodesByLevel();
    
    for (std::size_t i = 0; i < treeNodeLevel.size(); ++i) {
        for (std::size_t j = 0; j < treeNodeLevel[i].size(); ++j) {
            int treeNodeIdx = treeNodeLevel[i][j];
            if (!this->treeNodes[treeNodeIdx].isLeafNode()) {
                // For each pair of children we retrieve the border-to-border distances
                for (std::size_t k = 0; k < this->treeNodes[treeNodeIdx].children.size(); ++k) {
                    int sourceChild = this->treeNodes[treeNodeIdx].children[k];
                    int sourceChildOffset = this->treeNodes[treeNodeIdx].getChildOffsetInChildBorderVec(k);
                    // For each child's border find the distances to each border of every other child
                    for (std::size_t l = 0; l < this->treeNodes[sourceChild].bordersVec.size(); ++l) {
                        int sourceBorderIdx = sourceChildOffset + l; // Position of source border in parent's distance matrix
#if defined(GTREE_GOOGLE_DENSEHASH_DIST_MATRIX)
                        if (this->distanceMatrix.find(this->treeNodes[sourceChild].bordersVec[l]) == this->distanceMatrix.end()) {
                            this->distanceMatrix[this->treeNodes[sourceChild].bordersVec[l]].set_empty_key(constants::UNUSED_NODE_ID);
                        }
#endif
                        for (std::size_t m = 0; m < this->treeNodes[treeNodeIdx].children.size(); ++m) {
                            int targetChild = this->treeNodes[treeNodeIdx].children[m];
                            int targetChildOffset = this->treeNodes[treeNodeIdx].getChildOffsetInChildBorderVec(m);
                            for (std::size_t n = 0; n < this->treeNodes[targetChild].bordersVec.size(); ++n) {
                                int targetBorderIdx = targetChildOffset + n;
                                this->distanceMatrix[this->treeNodes[sourceChild].bordersVec[l]][this->treeNodes[targetChild].bordersVec[n]] = 
                                    this->treeNodes[treeNodeIdx].distanceMatrix[sourceBorderIdx*this->treeNodes[treeNodeIdx].matrixRowLength+targetBorderIdx];
                            }
                        }
                    }
                }
            } else {
                // For leaf nodes we retrieve the distance from each border to each leaf vertice
                for (std::size_t k = 0; k < this->treeNodes[treeNodeIdx].bordersVec.size(); ++k) {
#if defined(GTREE_GOOGLE_DENSEHASH_DIST_MATRIX)
                    if (this->distanceMatrix.find(this->treeNodes[treeNodeIdx].bordersVec[k]) == this->distanceMatrix.end()) {
                        this->distanceMatrix[this->treeNodes[treeNodeIdx].bordersVec[k]].set_empty_key(constants::UNUSED_NODE_ID);
                    }
#endif
                    for (std::size_t l = 0; l < this->treeNodes[treeNodeIdx].leafVerticesVec.size(); ++l) {
                        this->distanceMatrix[this->treeNodes[treeNodeIdx].bordersVec[k]][this->treeNodes[treeNodeIdx].leafVerticesVec[l]] = 
                            this->treeNodes[treeNodeIdx].distanceMatrix[k*this->treeNodes[treeNodeIdx].matrixRowLength+l];
                    }
                }
            }
        }
    }
}

int Gtree::getComputations(int leafIdx, int targetIdx)
{
    int firstLCAChild;
    std::vector<int> gtreePath = this->getGtreePath(leafIdx,targetIdx,firstLCAChild);
    
    int numComputations = 0;
    std::size_t x = 0;
    
    if (gtreePath.size() > 0) {
        assert (gtreePath.size() >= 2 && "Invalid path - must contain at least two nodes e.g. two leaves or leaf and parent");
        numComputations += this->treeNodes[gtreePath[0]].getNumBorders();
        
        // Add computation up the path until we reach firstLCAChild
        for (x = 1; x < gtreePath.size() && gtreePath[x-1] != firstLCAChild; ++x) {
            numComputations += this->treeNodes[gtreePath[x-1]].bordersVec.size()*this->treeNodes[gtreePath[x]].bordersVec.size();
        }
        
        // Add computation between the two children of the LCA
        // Note: If the parent of leafIdx is the firstLCAChild then this loop does nothing
        // but x is now smaller than the gtreePath (because targetIdx is in the same path
        // as the leafIdx) then the following code does nothing
        if (x < gtreePath.size()) {
            numComputations += this->treeNodes[gtreePath[x-1]].bordersVec.size()*this->treeNodes[gtreePath[x]].bordersVec.size();
            ++x;
        }
        
        // Add computation down the path until we reach the end
        for (; x < gtreePath.size(); ++x) {
            numComputations += this->treeNodes[gtreePath[x-1]].bordersVec.size()*this->treeNodes[gtreePath[x]].bordersVec.size();
        }
    }

    return numComputations;
}

/*
 * GtreeNode
 */

GtreeNode::GtreeNode(int treeIdx, int parentIdx, int numVertices): 
    treeIdx(treeIdx), parentIdx(parentIdx), numVertices(numVertices), isLeaf(false) {}

GtreeNode::GtreeNode() {}
    
bool GtreeNode::isLeafNode()
{
    return this->isLeaf;
}

void GtreeNode::setLeafNode()
{
    this->isLeaf = true;
}

int GtreeNode::getParentIdx()
{
    return this->parentIdx;
}

int GtreeNode::getTreeIdx()
{
    return this->treeIdx;
}

int GtreeNode::getNumVertices()
{
    return this->numVertices;
}

int GtreeNode::getParentChildIdx()
{
    return this->parentChildIdx;
}

void GtreeNode::setParentChildIdx(int idx)
{
    this->parentChildIdx = idx;
}

int GtreeNode::addChild(int treeIdx)
{
    this->children.push_back(treeIdx);
    return this->children.size()-1;
}

std::vector<int>& GtreeNode::getChildren()
{
    return this->children;
}

int GtreeNode::addBorder(NodeID node)
{
    // Borders are only ever added to G-tree node once
    this->bordersVec.push_back(node);
    int bordersVecIdx = this->bordersVec.size()-1;
    if (!this->isLeaf) {
        // If this is not a leaf, it's childBordersVec will be full
        // so we want to remember the array index for this border in
        // childBordersVec
        this->borderOffsetsInChildBorderVec.push_back(childBorderToChildBorderVecIdx[node]);
    } else {
        this->borderOffsetsInChildBorderVec.push_back(leafVerticeToLeafVerticeVecIdx[node]);
    }
    return bordersVecIdx;
}

std::vector<NodeID>& GtreeNode::getBorders()
{
    return this->bordersVec;
}

int GtreeNode::getNumBorders()
{
    return this->bordersVec.size();
}

int GtreeNode::addLeafVertex(NodeID node)
{
    this->leafVerticesUset.insert(node);
    this->leafVerticesVec.push_back(node);
    int leafVerticesVecIdx = this->leafVerticesVec.size()-1;
    this->leafVerticeToLeafVerticeVecIdx[node] = leafVerticesVecIdx;
    return leafVerticesVecIdx;
}

bool GtreeNode::isLeafVertex(NodeID node)
{
    // If iterator is not equal to end, border exists
    return this->leafVerticesUset.find(node) != this->leafVerticesUset.end();
}

std::vector<NodeID>& GtreeNode::getLeafVertices()
{
    return this->leafVerticesVec;
}

std::unordered_set<NodeID >& GtreeNode::getLeafVerticesUset()
{
    return this->leafVerticesUset;
}

void GtreeNode::addChildBorder(NodeID node)
{
    // Child borders are only ever added once as they only belong to one child
    this->childBordersUset.insert(node);
    this->childBordersVec.push_back(node);
    this->childBorderToChildBorderVecIdx[node] = this->childBordersVec.size()-1;
}

bool GtreeNode::isChildBorder(NodeID node)
{
    // If iterator is not equal to end, border exists
    return this->childBordersUset.find(node) != this->childBordersUset.end();
}

std::vector<NodeID>& GtreeNode::getChildBorders()
{
    return this->childBordersVec;
}

std::unordered_set<NodeID>& GtreeNode::getChildBordersUset()
{
    return this->childBordersUset;
}

void GtreeNode::addGtreePathNode(int nodeIdx)
{
    this->gtreePath.push_back(nodeIdx);
}

void GtreeNode::addGtreePathNodes(const std::vector<int>& parentPath)
{
    this->gtreePath.insert(this->gtreePath.end(),parentPath.begin(),parentPath.end());
}

std::vector<int>& GtreeNode::getGtreePathFromRoot()
{
    return this->gtreePath;
}

std::vector<EdgeWeight>& GtreeNode::getDistanceMatrix()
{
    return this->distanceMatrix;
}

int GtreeNode::getBorderIdxInChildBorderVec(int borderIdx)
{
    return this->borderOffsetsInChildBorderVec[borderIdx];
}

int GtreeNode::getChildOffsetInChildBorderVec(int childIdx)
{
    return this->childOffsetsInChildBorderVec[childIdx];
}

void GtreeNode::addChildOffsetInChildBorderVec()
{
    this->childOffsetsInChildBorderVec.push_back(this->childBordersVec.size());
}

void GtreeNode::printNode()
{
    std::cout << "Index: " << this->treeIdx << std::endl;
    std::cout << "Parent Index: " << this->parentIdx << std::endl;
    std::cout << "Parent's Child Vec Index: " << this->parentChildIdx << std::endl;
    std::cout << "isLeaf: " << this->isLeaf << std::endl  << std::endl;
    std::cout << "Children (" << this->children.size() << "): ";
    for (std::size_t i = 0; i < this->children.size(); ++i) {
        if (i != 0) {
            std::cout << ", ";
        }
        std::cout << i << " => " << children[i];
    }
    std::cout << std::endl << std::endl;
    std::cout << "Borders (" << this->bordersVec.size() << "): ";
    for (std::size_t i = 0; i < this->bordersVec.size(); ++i) {
        if (i != 0) {
            std::cout << ", ";
        }
        std::cout << i << " => " << bordersVec[i];
    }
    std::cout << std::endl << std::endl;
    std::cout << "Child Borders (" << this->childBordersVec.size() << "): ";
    for (std::size_t i = 0; i < this->childBordersVec.size(); ++i) {
        if (i != 0) {
            std::cout << ", ";
        }
        std::cout << i << " => " << childBordersVec[i];
    }
    std::cout << std::endl << std::endl;
    std::cout << "Border Offsets In Child Border Vec (" << this->borderOffsetsInChildBorderVec.size() << "): ";
    for (std::size_t i = 0; i < this->borderOffsetsInChildBorderVec.size(); ++i) {
        if (i != 0) {
            std::cout << ", ";
        }
        std::cout << i << " => " << borderOffsetsInChildBorderVec[i];
    }
    std::cout << std::endl << std::endl;
    std::cout << "Child Offsets In Child Border Vec (" << this->childOffsetsInChildBorderVec.size() << "): ";
    for (std::size_t i = 0; i < this->childOffsetsInChildBorderVec.size(); ++i) {
        if (i != 0) {
            std::cout << ", ";
        }
        std::cout << i << " => " << childOffsetsInChildBorderVec[i];
    }
    std::cout << std::endl << std::endl;
    std::cout << "Leaf Vertices (" << this->leafVerticesVec.size() << "): ";
    for (std::size_t i = 0; i < this->leafVerticesVec.size(); ++i) {
        if (i != 0) {
            std::cout << ", ";
        }
        std::cout << i << " => " << leafVerticesVec[i];
    }
    std::cout << std::endl << std::endl;
    std::cout << "Distance Matrix (" << this->distanceMatrix.size() << "): ";
    for (std::size_t i = 0; i < this->distanceMatrix.size(); ++i) {
        if (i != 0) {
            std::cout << ", ";
        }
        std::cout << i << " => " << distanceMatrix[i];
    }
    std::cout << std::endl << std::endl;

}

double GtreeNode::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(int)*5;
    memoryUsage += sizeof(bool);
    memoryUsage += sizeof(this->children);
    memoryUsage += sizeof(this->bordersVec);
    memoryUsage += sizeof(this->leafVerticesVec);
    memoryUsage += sizeof(this->gtreePath);
    memoryUsage += sizeof(this->distanceMatrix);
    memoryUsage += sizeof(this->borderOffsetsInChildBorderVec);
    memoryUsage += sizeof(this->childOffsetsInChildBorderVec);
    memoryUsage += sizeof(int)*this->children.size();
    memoryUsage += sizeof(NodeID)*this->bordersVec.size();
    memoryUsage += sizeof(NodeID)*this->leafVerticesVec.size();
    memoryUsage += sizeof(int)*this->gtreePath.size();
    memoryUsage += sizeof(EdgeWeight)*this->distanceMatrix.size();
    memoryUsage += sizeof(int)*this->borderOffsetsInChildBorderVec.size();
    memoryUsage += sizeof(int)*this->childOffsetsInChildBorderVec.size();
    return memoryUsage;
}

double GtreeNode::computeMemoryUsageBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    memoryUsage += sizeof(int)*this->children.capacity();
    memoryUsage += sizeof(NodeID)*this->bordersVec.capacity();
    memoryUsage += sizeof(NodeID)*this->leafVerticesVec.capacity();
    memoryUsage += sizeof(NodeID)*this->childBordersVec.capacity();
    memoryUsage += utility::estimateUnorderedMapMemoryUsageBytes(this->childBordersUset.size(),sizeof(NodeID),this->childBordersUset.bucket_count());
    memoryUsage += utility::estimateUnorderedMapMemoryUsageBytes(this->leafVerticesUset.size(),sizeof(NodeID),this->leafVerticesUset.bucket_count());
    memoryUsage += sizeof(int)*this->gtreePath.capacity();
    memoryUsage += sizeof(EdgeWeight)*this->distanceMatrix.capacity();
    memoryUsage += sizeof(int)*this->borderOffsetsInChildBorderVec.capacity();
    memoryUsage += sizeof(int)*this->childOffsetsInChildBorderVec.capacity();
    memoryUsage += utility::estimateUnorderedMapMemoryUsageBytes(this->childBorderToChildBorderVecIdx.size(),sizeof(std::pair<NodeID,int>),this->childBorderToChildBorderVecIdx.bucket_count());
    memoryUsage += utility::estimateUnorderedMapMemoryUsageBytes(this->leafVerticeToLeafVerticeVecIdx.size(),sizeof(std::pair<NodeID,int>),this->leafVerticeToLeafVerticeVecIdx.bucket_count());
    return memoryUsage;
}

double GtreeNode::computeDistanceMatrixMemoryUsageBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(this->distanceMatrix);
    memoryUsage += sizeof(EdgeWeight)*this->distanceMatrix.capacity();
    return memoryUsage;
}

/*
 * OccurenceList
 */

OccurenceList::OccurenceList(std::string setType, double setDensity, int setVariable, int setSize):
                             objSetType(setType), objSetDensity(setDensity), objSetVariable(setVariable), objSetSize(setSize)
{}

OccurenceList::OccurenceList() {}

void OccurenceList::addLeafOccurence(int leafIdx, NodeID object)
{
    if (this->leafOccurenceSet.find(leafIdx) == this->leafOccurenceSet.end()) {
        std::unordered_set<NodeID> occurenceSet;
        occurenceSet.insert(object);
        std::vector<NodeID> occurenceList;
        occurenceList.push_back(object);
        this->leafOccurenceSet[leafIdx] = occurenceSet;
        this->leafOccurenceList[leafIdx] = occurenceList;
    } else {
        if (this->leafOccurenceSet[leafIdx].find(object) == this->leafOccurenceSet[leafIdx].end()) {
            this->leafOccurenceSet[leafIdx].insert(object);
            this->leafOccurenceList[leafIdx].push_back(object);
        }
    }
}

void OccurenceList::addParentOccurence(int parentIdx, int childIdx)
{
    if (this->nonLeafOccurenceSet.find(parentIdx) == this->nonLeafOccurenceSet.end()) {
        std::unordered_set<int> occurenceSet;
        occurenceSet.insert(childIdx);
        this->nonLeafOccurenceSet[parentIdx] = occurenceSet;
        std::vector<int> occurenceList;
        occurenceList.push_back(childIdx);
        this->nonLeafOccurenceList[parentIdx] = occurenceList;
    } else {
        if (this->nonLeafOccurenceSet[parentIdx].find(childIdx) == this->nonLeafOccurenceSet[parentIdx].end()) {
            this->nonLeafOccurenceSet[parentIdx].insert(childIdx);
            this->nonLeafOccurenceList[parentIdx].push_back(childIdx);
        }
    }
}

const std::vector<NodeID>& OccurenceList::getLeafObjects(int leafIdx)
{
    return this->leafOccurenceList[leafIdx];
}

const std::unordered_set<NodeID>& OccurenceList::getLeafObjectsSet(int leafIdx)
{
    return this->leafOccurenceSet[leafIdx];
}

const std::vector<int>& OccurenceList::getNonLeafOccurenceList(int treeIdx)
{
    return this->nonLeafOccurenceList[treeIdx];
}

std::string OccurenceList::getObjSetType()
{
    return this->objSetType;
}

double OccurenceList::getObjSetDensity()
{
    return this->objSetDensity;
}

int OccurenceList::getObjSetSize()
{
    return this->objSetSize;
}

int OccurenceList::getObjSetVariable()
{
    return this->objSetVariable;
}

double OccurenceList::computeIndexSize()
{
    double memoryUsage = 0;
    for (auto it = this->leafOccurenceList.begin(); it != this->leafOccurenceList.end(); ++it) {
        memoryUsage += sizeof(int)+sizeof(*it)+sizeof(NodeID)*it->second.size();
    }
    for (auto it = this->nonLeafOccurenceList.begin(); it != this->nonLeafOccurenceList.end(); ++it) {
        memoryUsage += sizeof(int)+sizeof(*it)+sizeof(int)*it->second.size();
    }
    return memoryUsage/(1024*1024);
    // Note: We don't include leafOccurenceSet even though it is serialized, because
    // it contains the same data as leafOccurenceList and can easily reconstructed
    // from that instead. But we do consider this part of memory usage below.
}

double OccurenceList::computeMemoryUsage()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    memoryUsage += this->objSetType.size();
    memoryUsage += utility::estimateUnorderedMapMemoryUsageBytes(this->leafOccurenceList.size(),sizeof(std::pair<int,std::vector<NodeID>>),this->leafOccurenceList.bucket_count());
    for (auto it = this->leafOccurenceList.begin(); it != this->leafOccurenceList.end(); ++it) {
        memoryUsage += sizeof(NodeID)*it->second.capacity();
    }
    memoryUsage += utility::estimateUnorderedMapMemoryUsageBytes(this->nonLeafOccurenceList.size(),sizeof(std::pair<int,std::vector<int>>),this->nonLeafOccurenceList.bucket_count());
    for (auto it = this->nonLeafOccurenceList.begin(); it != this->nonLeafOccurenceList.end(); ++it) {
        memoryUsage += sizeof(int)*it->second.capacity();
    }
    memoryUsage += utility::estimateUnorderedMapMemoryUsageBytes(this->leafOccurenceSet.size(),sizeof(std::pair<int,std::unordered_set<NodeID>>),this->leafOccurenceSet.bucket_count());
    for (auto it = this->leafOccurenceSet.begin(); it != this->leafOccurenceSet.end(); ++it) {
        memoryUsage += utility::estimateUnorderedMapMemoryUsageBytes(it->second.size(),sizeof(NodeID),it->second.bucket_count());
    }
    memoryUsage += utility::estimateUnorderedMapMemoryUsageBytes(this->nonLeafOccurenceSet.size(),sizeof(std::pair<int,std::unordered_set<int>>),this->nonLeafOccurenceSet.bucket_count());
    for (auto it = this->nonLeafOccurenceSet.begin(); it != this->nonLeafOccurenceSet.end(); ++it) {
        memoryUsage += utility::estimateUnorderedMapMemoryUsageBytes(it->second.size(),sizeof(int),it->second.bucket_count());
    }
    return memoryUsage/(1024*1024);
}
