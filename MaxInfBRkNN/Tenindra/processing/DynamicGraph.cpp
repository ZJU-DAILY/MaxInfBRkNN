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

#include "DynamicGraph.h"
#include "../utility/utility.h"

#include <fstream>
#include <iostream>
#include <assert.h>
#include <cmath>

DynamicGraph::DynamicGraph(const DynamicGraph& obj)
{
    this->networkName = obj.networkName;
    this->numNodes = obj.numNodes;
    this->numEdges = obj.numEdges;
    for (std::size_t i = 0; i < obj.nodes.size(); ++i) {
        Node_DG nodeCopy(obj.nodes[i]);
        this->nodes.push_back(nodeCopy);
    }
}

// Create an updateable graph from non-updateable graph
DynamicGraph::DynamicGraph(Graph_RoadNetwork& obj)
{
    this->networkName = obj.getNetworkName();
    this->numNodes = obj.getNumNodes();
    this->numEdges = obj.getNumEdges();
    for (std::size_t i = 0; i < obj.coordinates.size(); ++i) {
        this->insertNode(i,obj.coordinates[i].first,obj.coordinates[i].second);
    }

    int adjListStart, adjListSize;
    NodeID numNodes = static_cast<NodeID>(this->numNodes);
    for (NodeID i = 0; i < numNodes; ++i) {
        adjListStart = obj.getEdgeListStartIndex(i);
        adjListSize = obj.getEdgeListSize(i);
        
        for (int j = adjListStart; j < adjListSize; ++j) {
            this->insertEdge(i,obj.edges[j].first,obj.edges[j].second);
        }        
    }    
}

void DynamicGraph::insertNode(NodeID nodeID, Coordinate x, Coordinate y)
{
    assert(nodeID == this->nodes.size() && "Node not inserted in order");
    
    Node_DG newNode(nodeID,x,y);
    this->nodes.push_back(newNode);
    
    // Store the minimum and maximum X and Y coordinates values
    // This is needed for quadtree contruction of this graph
    if (nodeID == 0) {
        // This is the first node, so we need to initialize min/max coordinates
        this->minXCoord = x;
        this->maxXCoord = x;
        this->minYCoord = y;
        this->maxYCoord = y;
        this->minXNode = nodeID;
        this->maxXNode = nodeID;
        this->minYNode = nodeID;
        this->maxYNode= nodeID;
    } else {
        if (x < this->minXCoord) {
            this->minXCoord = x;
            this->minXNode = nodeID;
        }
        if (x > this->maxXCoord) {
            this->maxXCoord = x;
            this->maxXNode = nodeID;
        }
        if (y < this->minYCoord) {
            this->minYCoord = y;
            this->minYNode = nodeID;
        }
        if (y > this->maxYCoord) {
            this->maxYCoord = y;
            this->maxYNode = nodeID;
        }
    }
}

void DynamicGraph::insertEdge(NodeID sourceID, NodeID targetID, EdgeWeight edgeWeight) {
    assert(!this->nodes[sourceID].edgeExists(targetID) && "Parallel edge detected");
    assert(sourceID != targetID && "Self-loop detected");
    if (this->type == constants::DISTANCE_WEIGHTS) {
        EdgeWeight euclidDist = std::ceil(this->getEuclideanDistance(sourceID,targetID));
        if (edgeWeight < euclidDist) {
            std::cout << this->networkName << std::endl;
            std::cout << sourceID << std::endl;
            std::cout << targetID << std::endl;
            std::cout << edgeWeight << std::endl;
            std::cout << euclidDist << std::endl;
        }
        assert(edgeWeight >= euclidDist && "Edge weight smaller than Euclidean distance between neighbours");
    }    
    assert(edgeWeight != 0 && "Zero-edge weight detected");
    
    this->nodes[sourceID].adjNodes.push_back(targetID);
    this->nodes[targetID].invAdjNodes.push_back(sourceID);
    this->nodes[sourceID].adjNodeWgts.push_back(edgeWeight);
    this->nodes[targetID].invAdjNodeWgts.push_back(edgeWeight);
    assert(this->nodes[sourceID].adjNodes.size() < (constants::MAX_EDGES-1) && "Too many neighbours - this will break EMPTY_COLOUR and SPECIAL_COLOUR constants used by SILC");
}

void DynamicGraph::insertImaginaryNonInvertibleEdge(NodeID sourceID, NodeID targetID, EdgeWeight edgeWeight) {
    // This function is for the benefit of G-tree (and any other method that
    // relies on adding new edges to make Dijkstra computation faster). In this
    // we do not need inverted edges.
    this->nodes[sourceID].adjNodes.push_back(targetID);
    this->nodes[sourceID].adjNodeWgts.push_back(edgeWeight);
}

std::string DynamicGraph::getNetworkName()
{
    return this->networkName;
}

int DynamicGraph::getNumNodes() {
    return this->numNodes;
}

int DynamicGraph::getNumEdges()
{
    return this->numEdges;
}

int DynamicGraph::countNodes() {
    return this->nodes.size();
}

int DynamicGraph::countEdges() {
    int numEdges = 0;
    for (auto it = this->nodes.begin(); it != this->nodes.end(); ++it) {
        numEdges += it->adjNodes.size();
    }
    return numEdges;
}

void DynamicGraph::parseGraphFile(std::string grFilePath, std::string coFilePath)
{
    NodeID sourceID, targetID;
    EdgeWeight weight;
    std::string name, type;
    int nodes, edges;
    Coordinate x, y;

    // Parse Coordinates File
    // We expect all coordinates to be numbered from 0 to n-1 and are supplied in order
    std::ifstream coordFile(coFilePath, std::ios::in);
    if (coordFile.is_open()) {
        coordFile >> name >> type >> nodes;
        this->networkName = name;
        this->numNodes = nodes;
            
        std::cout << "Now parsing " << nodes << " nodes for graph " << this->networkName << std::endl;

        while (coordFile >> sourceID >> x >> y)  {
           this->insertNode(sourceID,x,y);
        }
    } else {
        std::cerr << "Cannot open coordinates file " << coFilePath << std::endl;
        exit(1);
    }    
    coordFile.close();
    
    // Parse Graph File
    std::ifstream graphFS(grFilePath, std::ios::in);
    if (graphFS.is_open()) {
        // First line is graph information
        graphFS >> name >> type >> nodes >> edges;
        this->numEdges = edges;
        this->type = type; // Edges file sets type not coordinates
            
        std::cout << "Now parsing " << edges << " edges for graph " << this->networkName << std::endl;

        while (graphFS >> sourceID >> targetID >> weight)  {
            this->insertEdge(sourceID,targetID,weight);
        }
    } else {
        std::cerr << "Cannot open graph file " << grFilePath << std::endl;
        exit(1);
    }
    graphFS.close();
    
    std::cout << "\nSuccessfully parsed graph " << this->networkName << " with " << this->numNodes
                << " nodes and " << this->numEdges << " edges!" << std::endl;

    int countedEdges = this->countEdges();
    int countedNodes = this->countNodes();
    //std::cout << "\nSuccessfully parsed graph " << this->networkName << " with " << countedNodes
    //            << " nodes and " << countedEdges << " edges!" << std::endl;
    assert(countedEdges == this->getNumEdges() && "Number of edges added not equal to expected number");
    assert(countedNodes == this->getNumNodes() && "Number of nodes added not equal to expected number");
}

EuclideanDistance DynamicGraph::getEuclideanDistance(NodeID node1, NodeID node2)
{
    Coordinate x1 = this->nodes[node1].x, y1 = this->nodes[node1].y;
    Coordinate x2 = this->nodes[node2].x, y2 = this->nodes[node2].y;
    EuclideanDistance dist = std::sqrt(std::pow(x2-x1,2) + std::pow(y2-y1,2));
    return dist;
}

std::vector<NodeID> DynamicGraph::getNodesIDsVector()
{
    std::vector<NodeID> nodeIDs;
    for (NodeID i = 0; i < this->nodes.size(); ++i) {
        nodeIDs.push_back(i);
    }
    return nodeIDs;
}

std::unordered_set<NodeID> DynamicGraph::getNodesIDsUset()
{
    std::unordered_set<NodeID> nodeIDs;
    for (NodeID i = 0; i < this->nodes.size(); ++i) {
        nodeIDs.insert(i);
    }
    return nodeIDs;
}

const std::vector<NodeID>& DynamicGraph::getAdjNodes(NodeID node) const
{
    return this->nodes[node].adjNodes;
}

const std::vector<EdgeWeight>& DynamicGraph::getAdjNodeWgts(NodeID node) const
{
    return this->nodes[node].adjNodeWgts;
}

const std::vector<NodeID>& DynamicGraph::getInvAdjNodes(NodeID node) const
{
    return this->nodes[node].invAdjNodes;
}

const std::vector<EdgeWeight>& DynamicGraph::getInvAdjNodeWgts(NodeID node) const
{
    return this->nodes[node].invAdjNodeWgts;
}

void DynamicGraph::parseObjectFile(std::string objectSetFilePath, std::string& setType, double& setDensity, int& setVariable, int& setSize)
{
    this->nodeObjects.assign(this->numNodes,false);
    
    NodeID objectID;
    std::ifstream objSetFS(objectSetFilePath, std::ios::in);
    std::string networkName;
    
    if (objSetFS.is_open()) {
        // First line is object set information
        objSetFS >> networkName >> setType >> setDensity >> setSize;
        
        if (networkName != this->getNetworkName()) {
            std::cerr << "This object set was not generated for current graph " 
                << "(generated for " << networkName << ")!" << std::endl;
            exit(1);
        }
        
        while (objSetFS >> objectID)  {
            this->nodeObjects[objectID] = true;
        }
    } else {
        std::cerr << "Cannot open object set file " << objectSetFilePath << std::endl;
        exit(1);
    }
    
    objSetFS.close();
}

void DynamicGraph::parseObjectSet(std::vector<NodeID>& set)
{
    this->nodeObjects.assign(this->numNodes,false);
    for(std::size_t i = 0; i < set.size(); ++i) {
        this->nodeObjects[set[i]] = true;
    }
}

void DynamicGraph::resetAllObjects()
{
    this->nodeObjects.assign(this->numNodes,false);
}

bool DynamicGraph::isObject(NodeID nodeID)
{
    return this->nodeObjects[nodeID];
}


bool DynamicGraph::isConnected(NodeID node1, NodeID node2)
{
    // Since edge lists have few edges in road networks
    // this should be not and not much more inefficient
    // than using unordered_map to do reverse lookup
    for (auto it = this->nodes[node1].adjNodes.begin(); it != this->nodes[node1].adjNodes.end(); ++it) {
        if (*it == node2) {
            return true;
        }
    }
    return false;
}

bool DynamicGraph::isUndirectedGraph()
{
    bool undirected = true;
    for (auto it = this->nodes.begin(); it != this->nodes.end(); ++it) {
        for (auto adjNodeIt = it->adjNodes.begin(); adjNodeIt != it->adjNodes.end(); ++adjNodeIt) {
            if (!this->isConnected(*adjNodeIt,it->id)) {
                undirected = false;
                break;
            }
        }
        if (!undirected) {
            break;
        }
    }
    return undirected;

}

bool DynamicGraph::isConnectedGraph()
{
    bool connected = true;
    for (auto it = this->nodes.begin(); it != this->nodes.end(); ++it) {
        if (it->adjNodes.size() == 0) {
            connected  = false;
        }
    }
    return connected;
}

Coordinate DynamicGraph::getMaxCoordinateRange()
{
    // Return the largest range in coordinate values
    Coordinate xRange = maxXCoord-minXCoord, yRange = maxYCoord-minYCoord;
    return (xRange > yRange) ? xRange : yRange;
}

void DynamicGraph::getMinMaxCoordinates(Coordinate &minX, Coordinate &maxX, Coordinate &minY, Coordinate &maxY) {
    minX = this->minXCoord;
    maxX = this->maxXCoord;
    minY = this->minYCoord;
    maxY = this->maxYCoord;
}

EdgeWeight DynamicGraph::getEdgeWeight(NodeID u, NodeID v) {
    return this->nodes[u].getEdgeWeight(v);
}

EdgeWeight DynamicGraph::getNeighbourAndEdgeWeightByPosition(NodeID u, EdgeID pos, EdgeWeight& edgeWeight) {
    return this->nodes[u].getNeighbourAndEdgeWeightByPosition(pos,edgeWeight);
}

void DynamicGraph::getCoordinates(NodeID u, Coordinate &x, Coordinate &y) {
    this->nodes[u].getCoordinates(x,y);
}

void DynamicGraph::getTranslatedCoordinates(NodeID u, Coordinate &x, Coordinate &y, int xTranslation, int yTranslation) {
    // xTranslation and yTranslation represent the translation of 
    // the origin to keep all x and y values positive
    
    this->nodes[u].getCoordinates(x,y);
    x = x - xTranslation;
    y = y - yTranslation;
}

std::vector<int> DynamicGraph::countOutdegreeFrequencies()
{
    std::vector<int> frequencies;
    unsigned int outdegree;
    for (NodeID i = 0; i < this->nodes.size(); ++i) {
        outdegree = this->nodes[i].adjNodes.size();
        if (frequencies.size() <= outdegree) {
            frequencies.resize(outdegree+1,0);
        }
        ++frequencies[outdegree];
    }
    return frequencies;
}

void DynamicGraph::printGraphStats()
{
    std::cout << "\nStatistics for Network " << this->networkName << ":" << std::endl;
    std::cout << "Edge Type = " << this->type << ", Nodes = " <<this->numNodes << ", Edges = " << this->numEdges << std::endl;
    std::cout << "Max Coordinate Range = " << this->getMaxCoordinateRange() << std::endl;
    std::vector<int> outdegreeFrequencies = this->countOutdegreeFrequencies();
    std::cout << "\nOutdegree Frequencies:" << std::endl;
    for (std::size_t i = 0; i < outdegreeFrequencies.size(); ++i) {
        std::cout << "Outdegree of " << i << ": " << outdegreeFrequencies[i] << std::endl;
    }
    std::cout << "\nConnected = " << (this->isConnectedGraph() ? "Yes" : "No") << std::endl;
    std::cout << "Undirected = " << (this->isUndirectedGraph() ? "Yes" : "No") << std::endl;
    std::cout << "Min (X,Y) = " << "(" << minXCoord << "," << minYCoord << ")" << std::endl;
    std::cout << "Max (X,Y) = " << "(" << maxXCoord << "," << maxYCoord << ")" << std::endl;
    std::cout << "Min X Node= " << minXNode << ", Max X Node= " << maxXNode << std::endl;
    std::cout << "Min Y Node= " << minYNode << ", Max Y Node= " << maxYNode << std::endl;
}

double DynamicGraph::computeIndexSize()
{
    double memoryUsageBytes = 0;
    for (std::size_t i = 0; i < this->nodes.size(); ++i) {
        memoryUsageBytes += this->nodes[i].computeIndexSizeBytes();
    }
    return memoryUsageBytes/(1024*1024);
}

double DynamicGraph::computeINEIndexSize()
{
    double memoryUsageBytes = 0;
    for (std::size_t i = 0; i < this->nodes.size(); ++i) {
        memoryUsageBytes += this->nodes[i].computeINEIndexSizeBytes();
    }
    return memoryUsageBytes/(1024*1024);
}

double DynamicGraph::computeObjectSetIndexSize()
{
    double memoryUsageBytes = 0;
    memoryUsageBytes += this->nodeObjects.size()/8; // std::vector<bool> only use 1 bit per element
    return memoryUsageBytes/(1024*1024);
}

double DynamicGraph::computeMemoryUsage()
{
    double memoryUsageBytes = 0;
    memoryUsageBytes += sizeof(*this);
    memoryUsageBytes += this->networkName.size();
    memoryUsageBytes += this->type.size();
    for (std::size_t i = 0; i < this->nodes.size(); ++i) {
        memoryUsageBytes += this->nodes[i].computeMemoryUsageBytes();
    }
    memoryUsageBytes += sizeof(Node_DG)*(this->nodes.capacity()-this->nodes.size());
    memoryUsageBytes += this->nodeObjects.capacity()/8; // std::vector<bool> only use 1 bit per element
    return memoryUsageBytes/(1024*1024);
}

/*
 * DynamicGraph::Node
 */

Node_DG::Node_DG(const Node_DG& obj)
{
    this->id = obj.id;
    this->x = obj.x;
    this->y = obj.y;
    this->adjNodes = obj.adjNodes;
    this->adjNodeWgts = obj.adjNodeWgts;
}

const std::vector<NodeID> &Node_DG::getAdjNodes() const {
    return this->adjNodes;
}

const std::vector< EdgeWeight > &Node_DG::getAdjNodeWgts() const {
    return this->adjNodeWgts;
}

EdgeWeight Node_DG::getEdgeWeight(NodeID v) {
    for (std::size_t i = 0; i < this->adjNodes.size(); ++i) {
        if (this->adjNodes[i] == v) {
            return adjNodeWgts[i];
        }
    }
    std::cerr << "There is no edge between two nodes" << std::endl;
    exit(1);
    return 0;    
}

NodeID Node_DG::getNeighbourAndEdgeWeightByPosition(EdgeID pos, EdgeWeight& edgeWeight)
{
    assert (pos < adjNodes.size() && "No such edge exists");
    edgeWeight = adjNodeWgts[pos];
    return adjNodes[pos];
}

void Node_DG::getCoordinates(Coordinate &x, Coordinate &y) {
    x = this->x;
    y = this->y;
}

bool Node_DG::edgeExists(NodeID neighbour)
{
    for (std::size_t i = 0; i < this->adjNodes.size(); ++i) {
        if (this->adjNodes[i] == neighbour) {
            return true;
        }
    }
    // Edge not found so we return false
    return false;
}

double Node_DG::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(NodeID);
    memoryUsage += sizeof(Coordinate)*2;
    memoryUsage += sizeof(this->adjNodes);
    memoryUsage += sizeof(this->adjNodeWgts);
    memoryUsage += sizeof(this->invAdjNodes);
    memoryUsage += sizeof(this->invAdjNodeWgts);
    memoryUsage += sizeof(NodeID)*this->adjNodes.size();
    memoryUsage += sizeof(EdgeWeight)*this->adjNodeWgts.size();
    memoryUsage += sizeof(NodeID)*this->invAdjNodes.size();
    memoryUsage += sizeof(EdgeWeight)*this->invAdjNodeWgts.size();
    return memoryUsage;
}

double Node_DG::computeINEIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(NodeID);
    memoryUsage += sizeof(this->adjNodes);
    memoryUsage += sizeof(this->adjNodeWgts);
    memoryUsage += sizeof(NodeID)*this->adjNodes.size();
    memoryUsage += sizeof(EdgeWeight)*this->adjNodeWgts.size();
    return memoryUsage;
}

double Node_DG::computeMemoryUsageBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    memoryUsage += sizeof(NodeID)*this->adjNodes.capacity();
    memoryUsage += sizeof(EdgeWeight)*this->adjNodeWgts.capacity();
    memoryUsage += sizeof(NodeID)*this->invAdjNodes.capacity();
    memoryUsage += sizeof(EdgeWeight)*this->invAdjNodeWgts.capacity();
    return memoryUsage;
}
