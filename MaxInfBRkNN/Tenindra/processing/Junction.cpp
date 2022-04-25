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

#include "Junction.h"

#include <deque>
#include <tuple>

Junction::Junction(std::string networkName, int numNodes, int numEdges): 
    networkName(networkName), numNodes(numNodes), numEdges(numEdges) {
    
}
    
std::string Junction::getNetworkName()
{
    return this->networkName;
}

int Junction::getNumEdges()
{
    return this->numEdges;
}

int Junction::getNumNodes()
{
    return this->numNodes;
}

void Junction::buildJunction(Graph_RoadNetwork& graph)
{
    this->isNoThruRoadNode.assign(graph.getNumNodes(),false);
    this->noThruRoadID.assign(graph.getNumNodes(),-1);
    this->nextJunctionNodeAndDistance.assign(graph.getNumEdges(),IntDistancePair(-2,0));
    std::vector<bool> isSettledJunction(graph.getNumNodes());
    // Tuple: 0: Node, Source Junction, Source Edge, Distance from Junction, Previous Node
    std::deque<std::tuple<NodeID,NodeID,NodeID,EdgeWeight,NodeID>> queue; // Pair of node, source junction and edge (as index in edge vector) taken from source juction, and previous node

    int adjListStart, nextAdjListStart, outDegree;
    NodeID currNode, adjNode, prevNode, sourceEdge, sourceJuction, nextNode, sourceNode;
    EdgeWeight adjNodeWgt, sourceJuctionDist;
    unsigned int numIntraJuctionRoadNodes = 0, numNoThruRoadNodes = 0;

    // Find a root node that has more than 3 edges (a junction)
    std::unordered_set<NodeID> visitedNodes;
    queue.push_back(std::make_tuple(0,0,0,0,0)); // Insert node 0
    while (queue.size() > 0) {
        currNode = std::get<0>(queue.front());
        queue.pop_front();

        if (visitedNodes.find(currNode) == visitedNodes.end()) {
            visitedNodes.insert(currNode);
            
            adjListStart = graph.getEdgeListStartIndex(currNode);
            nextAdjListStart = graph.getEdgeListSize(currNode);
            outDegree = nextAdjListStart-adjListStart;
            
            // If this node has more than 2 neighbours then it is a juction and we have found a root
            if (outDegree > 2) {
                break;
            }

            // Add each neighbour to queue
            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                adjNode = graph.edges[i].first;
                queue.push_back(std::make_tuple(adjNode,0,0,0,0));
            }
        }
    }
    queue.clear();
    //visitedNodes.clear();

    isSettledJunction[currNode] = true;
    adjListStart = graph.getEdgeListStartIndex(currNode);
    nextAdjListStart = graph.getEdgeListSize(currNode);
    outDegree = nextAdjListStart - adjListStart;
    assert(outDegree > 2 && "Root node is not a junction");
    for (int i = adjListStart; i < nextAdjListStart; ++i) {
        adjNode = graph.edges[i].first;
        adjNodeWgt =  graph.edges[i].second;
        queue.push_back(std::make_tuple(adjNode,currNode,i,adjNodeWgt,currNode));
    }
    
    while (queue.size() > 0) {
        currNode = std::get<0>(queue.front());
        sourceJuction = std::get<1>(queue.front());
        sourceEdge = std::get<2>(queue.front()); // First edge (from junction)
        sourceJuctionDist = std::get<3>(queue.front()); 
        prevNode = std::get<4>(queue.front());
        queue.pop_front();
        
        adjListStart = graph.getEdgeListStartIndex(currNode);
        nextAdjListStart = graph.getEdgeListSize(currNode);
        outDegree = nextAdjListStart - adjListStart;
        
        if (outDegree > 2) {
            // If this is the first time we have visited this junction
            // we need add all its edge make sure we propagate whole graph
            // Note: It is is possible we may reach here again as there maybe
            // incident edges that we have processed yet
            if (!isSettledJunction[currNode]) {
                isSettledJunction[currNode] = true;
                for (int i = adjListStart; i < nextAdjListStart; ++i) {
                    // Only update those we haven't already settled
                    adjNode = graph.edges[i].first;
                    adjNodeWgt =  graph.edges[i].second;
                    queue.push_back(std::make_tuple(adjNode,currNode,i,adjNodeWgt,currNode));
                }
            }

            // This is junction node so we mark all edges leading to it from the source juction
            // as with this juction until we reach this node
            this->nextJunctionNodeAndDistance[sourceEdge].first = currNode;
            this->nextJunctionNodeAndDistance[sourceEdge].second = sourceJuctionDist; // Includes weight of first edge
            nextNode = graph.edges[sourceEdge].first;
            sourceNode = nextNode;
            prevNode = sourceJuction;
            sourceJuctionDist -= graph.edges[sourceEdge].second;
            numIntraJuctionRoadNodes = 1;
            while (sourceNode != currNode) {
                ++numIntraJuctionRoadNodes;
                adjListStart = graph.getEdgeListStartIndex(sourceNode);
                nextAdjListStart = graph.getEdgeListSize(sourceNode);
                outDegree = nextAdjListStart-adjListStart;
                assert(outDegree == 2 && "We arrived here from an incorrect source juction");
                for (int i = adjListStart; i < nextAdjListStart; ++i) {
                    adjNode = graph.edges[i].first;
                    adjNodeWgt =  graph.edges[i].second;
                    if (adjNode != prevNode) {
                        // Only traverse back towards the juction we came from
                        // i.e. don't take the edge back to the prev node
                        this->nextJunctionNodeAndDistance[i].first = currNode;
                        this->nextJunctionNodeAndDistance[i].second = sourceJuctionDist;
                        nextNode = adjNode;
                        sourceJuctionDist -= adjNodeWgt;
                    }
                    // Note: The other edge (i.e. leading from currNode will be handled separately
                    // when we search outwards from the junction
                }
                prevNode = sourceNode;
                sourceNode = nextNode;
            }
            if (this->intraJuctionRoadLengths.size() <= numIntraJuctionRoadNodes) {
                // Note: resize doesn't change the values of existing values (unlike assign)
                this->intraJuctionRoadLengths.resize(numIntraJuctionRoadNodes+1,0);
            }
            this->intraJuctionRoadLengths[numIntraJuctionRoadNodes]++;
        } else if (outDegree == 2) {
            for (int i = adjListStart; i < nextAdjListStart; ++i) {
                // Only traverse away from the juction we came from
                // i.e. don't take the edge back to the prev node
                adjNode = graph.edges[i].first;
                adjNodeWgt =  graph.edges[i].second;
                if (adjNode != prevNode) {
                    // Note: We can't just avoid visited node because this node
                    // may visited again from the target junction back to source junction
                    // when the roles reversed
                    queue.push_back(std::make_tuple(adjNode,sourceJuction,sourceEdge,sourceJuctionDist+adjNodeWgt,currNode));
                }
            }
        } else /*if (outDegree == 1)*/ {
            assert(outDegree == 1);
            // This is deadend node so we mark all edges leading to it from source juction
            // as deadend edges and mark edges leading away from it with source juction
            // until we reach the source juction itself

            // Mark the source junction edge we followed as a deadend edge
            this->nextJunctionNodeAndDistance[sourceEdge].first = -1;
            
            this->isNoThruRoadNode[currNode] = true;
            this->noThruRoadID[currNode] = this->numNoThruRoads; // We give each node an ID that corresponds to which no thru road it belongs on
            this->nextJunctionNodeAndDistance[adjListStart].first = sourceJuction;
            this->nextJunctionNodeAndDistance[adjListStart].second = sourceJuctionDist; // Assumes undirected graph
            nextNode = prevNode; // == graph.edges[adjListStart].first
            prevNode = currNode;
            currNode = nextNode;
            sourceJuctionDist -= graph.edges[adjListStart].second;
            numNoThruRoadNodes = 1;
            while (currNode != sourceJuction) {
                ++numNoThruRoadNodes;
                this->isNoThruRoadNode[currNode] = true;
                this->noThruRoadID[currNode] = this->numNoThruRoads;
                adjListStart = graph.getEdgeListStartIndex(currNode);
                nextAdjListStart = graph.getEdgeListSize(currNode);
                outDegree = nextAdjListStart-adjListStart;
                assert(outDegree == 2 && "We arrived here from an incorrect source juction");
                for (int i = adjListStart; i < nextAdjListStart; ++i) {
                    adjNode = graph.edges[i].first;
                    adjNodeWgt =  graph.edges[i].second;
                    if (adjNode != prevNode) {
                        // Only traverse back towards the juction we came from
                        // i.e. don't take the edge back to the prev node
                        this->nextJunctionNodeAndDistance[i].first = sourceJuction;
                        this->nextJunctionNodeAndDistance[i].second = sourceJuctionDist; // Assumes undirected graph
                        nextNode = adjNode;
                        sourceJuctionDist -= adjNodeWgt;
                    } else {
                        // We mark the edge towards the dead as leading to a deadend
                        this->nextJunctionNodeAndDistance[i].first = -1;
                    }
                }
                prevNode = currNode;
                currNode = nextNode;
            }
            ++this->numNoThruRoads;
            if (this->noThruRoadLengths.size() <= numNoThruRoadNodes) {
                // Note: resize doesn't change the values of existing values (unlike assign)
                this->noThruRoadLengths.resize(numNoThruRoadNodes+1,0);
            }
            this->numTotalNoThruRoadNodes += numNoThruRoadNodes;
            this->noThruRoadLengths[numNoThruRoadNodes]++;            
        }
    }
    
    // Check there are -2 values in nextJunctionNode (means edge never set)
    for (std::size_t i = 0; i < this->numEdges; ++i) {
        if (this->nextJunctionNodeAndDistance[i].first == -2) {
            std::cout << "Edge " << i << " with target " << graph.edges[i].first << " was not set" << std::endl;
        }
        assert(this->nextJunctionNodeAndDistance[i].first != -2 && "Unset edge found");
    }
    
    // Note: To extend this to directed graphs we need to consider incident edges e.g. a 
    // deadend might be a one-way street or 2-outdegree street might be a one-way Y split
}

double Junction::computeIndexSize()
{
    double memoryUsageBytes = 0;
    memoryUsageBytes += this->isNoThruRoadNode.size()/8;
    memoryUsageBytes += sizeof(int)*this->noThruRoadID.size();
    memoryUsageBytes += sizeof(IntDistancePair)*this->nextJunctionNodeAndDistance.size();
    return memoryUsageBytes/(1024*1024);
}

double Junction::computeMemoryUsage()
{
    double memoryUsageBytes = 0;
    memoryUsageBytes += sizeof(*this);
    memoryUsageBytes += this->networkName.size();
    memoryUsageBytes += this->isNoThruRoadNode.capacity()/8;
    memoryUsageBytes += sizeof(int)*this->noThruRoadID.capacity();
    memoryUsageBytes += sizeof(IntDistancePair)*this->nextJunctionNodeAndDistance.capacity();
    return memoryUsageBytes/(1024*1024);
}

void Junction::verifyJunctionsAndDistances(Graph_RoadNetwork& graph)
{
    int outDegree, endNodeOutdegree;
    int junctionID;
    NodeID endNode;
    EdgeWeight junctionDist, distanceToEnd;
    bool isOnNoThruRoad;
    int adjListStart, nextAdjListStart;
    
    for (NodeID i = 0; i < this->numNodes; ++i) {
        // Verify junction or deadend status
        isOnNoThruRoad = this->isNoThruRoadNode[i];
        adjListStart = graph.getEdgeListStartIndex(i);
        nextAdjListStart = graph.getEdgeListSize(i);
        outDegree = nextAdjListStart - adjListStart;
        for (int j = adjListStart; j < nextAdjListStart; ++j) {
            junctionID = this->nextJunctionNodeAndDistance[j].first;
            junctionDist = this->nextJunctionNodeAndDistance[j].second;
            if (junctionID == -1 && outDegree <= 2) {
                assert(isOnNoThruRoad && "Node is not a no thru road node without being a junction");
            } 
            
            endNode = graph.findRoadEnd(i,j,endNodeOutdegree,distanceToEnd);
            if (junctionID == -1) {
                assert(endNodeOutdegree == 1 && "Final node on the road is not an deadend");
            } else {
                assert(static_cast<NodeID>(junctionID) == endNode && "Stored junction is incorrect");
                assert(distanceToEnd == junctionDist && "Junction distance is incorrect");
            }
        }
    }
    std::cout << "Junction index verification completed!" << std::endl;
}

void Junction::printGraphStats()
{
    std::cout << "\nJunctions Statistics for Network " << this->networkName << ":" << std::endl;
    std::cout << "Number of No Thru Roads = " << this->numNoThruRoads << std::endl;
    std::cout << "Number of Nodes in No Thru Roads = " << this->numTotalNoThruRoadNodes << std::endl;
    int totalNoThruRoadsEdges = 0, totalIntraJunctionRoadsEdges = 0;
    std::cout << "Intra Junction Road Lengths Frequencies:" << std::endl;
    for (std::size_t i = 0; i < this->intraJuctionRoadLengths.size(); ++i) {
        if (this->intraJuctionRoadLengths[i] != 0) {
            std::cout << "\tLength = " << i << ": " << this->intraJuctionRoadLengths[i] << " (Total Edges = " << this->intraJuctionRoadLengths[i]*i << ")" << std::endl;
        }
        totalIntraJunctionRoadsEdges += this->intraJuctionRoadLengths[i]*i;
    }
    std::cout << "No Thru Road Lengths Frequencies:" << std::endl;
    for (std::size_t i = 0; i < this->noThruRoadLengths.size(); ++i) {
        if (this->noThruRoadLengths[i] != 0) {
            std::cout << "\tLength = " << i << ": " << this->noThruRoadLengths[i] << " (Total Edges = " << this->noThruRoadLengths[i]*i << ")" << std::endl;
        }
        totalNoThruRoadsEdges += this->noThruRoadLengths[i]*i;
    }
    std::cout << "Total Edges Captured = " << totalIntraJunctionRoadsEdges+totalNoThruRoadsEdges*2 << std::endl;
}
