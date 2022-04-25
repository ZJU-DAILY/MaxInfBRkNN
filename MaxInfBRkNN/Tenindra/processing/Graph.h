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

#ifndef _GRAPH_H
#define _GRAPH_H

#include "../common.h"

#include <vector>
#include <unordered_set>
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <iostream>
#include <sstream>
#include "../../source_config.h"
#include  <map>

using namespace std;
/*
 * Directed Graph Stored In-Memory
 */



class Graph_RoadNetwork {
    
    public:
        Graph_RoadNetwork() {}
        std::string getNetworkName();
        int getNumNodes();
        int getNumEdges();
        std::string getEdgeType();
        void setGraphData(std::string name, int numNodes, int numEdges, std::string edgeType);
        void insertNode(NodeID nodeID, Coordinate x, Coordinate y) {
            assert(nodeID == this->coordinates.size() && "Node not inserted in order");

            this->coordinates.push_back(CoordinatePair(x,y));
            this->numNodes +=1;
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
                this->maxYNode = nodeID;
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

        // We assume edges are added in order (i.e. all edges for sourceID are added together
        void insertEdge (NodeID sourceID, NodeID targetID, EdgeWeight edgeWeight) {
            assert(sourceID != targetID && "Self-loop detected");
            ///std::cout << "insert edge("<<sourceID<<","<<targetID <<"), dist_weight="<< edgeWeight << std:: endl;
            assert(edgeWeight != 0 && "Zero-edge weight detected");

            if (sourceID == this->firstEdgeIndex.size()) {
                // Then this is the first edge encountered
                this->firstEdgeIndex.push_back(this->edges.size());
            } else {
                assert(sourceID == this->firstEdgeIndex.size()-1 && "Edge not inserted in order");
                // Check if edge has already been
                int adjListStart = this->firstEdgeIndex[sourceID], adjListSize = this->edges.size();
                for (int i = adjListStart; i < adjListSize; ++i) {
                    assert(targetID != this->edges[i].first && "Parallel edge detected");
                }
            }
            this->edges.push_back(NodeEdgeWeightPair(targetID,edgeWeight));
            this->numEdges += 1;
        }
        //parse road network file into binary data file
        void parseGraphFile(std::string grFilePath, std::string coFilePath){
            NodeID sourceID, targetID;
            EdgeWeight weight;
            std::string name, type;
            unsigned int nodes, edges;
            Coordinate x, y;
            this->networkName = dataset_name;
            // Parse Coordinates File
            // We expect all coordinates to be numbered from 0 to n-1 and are supplied in order
            std::ifstream coordFile(coFilePath, std::ios::in);
            if (coordFile.is_open()) {
                //coordFile >> name >> type >> nodes;
                //this->networkName = name;
                //this->numNodes = nodes;

                std::cout << "Now parsing " << nodes << " nodes for graph " << this->networkName << std::endl;
                std::string line;
                while (getline(coordFile, line))  {
                    istringstream tt(line);
#ifndef Large_DATASET
                    double x_value, y_value;
                    tt >> sourceID >> x_value >> y_value;
                    x = x_value * scale_coordinate;
                    y = y_value * scale_coordinate;
                    this->insertNode(sourceID,x,y);
#else
                    tt >> sourceID >> x >> y;
                    this->insertNode(sourceID,x,y);
#endif
                }
            } else {
                std::cerr << "Cannot open coordinates file " << coFilePath << std::endl;
                exit(1);
            }
            coordFile.close();

            // Parse Graph File
            // We expect all edges to be supplied in order of the source ID (from 0 to n-1)
            std::ifstream graphFS(grFilePath, std::ios::in);
            if (graphFS.is_open()) {
                // First line is graph information
                //graphFS >> name >> type >> nodes >> edges;
                //this->numEdges = edges;
                //this->type = type; // Edges file sets type not coordinates

                std::cout << "Now parsing " << edges << " edges for graph " << this->networkName << std::endl;
                string line;
                while (getline(graphFS, line))  {
                    istringstream ss(line);
#ifdef Large_DATASET
                    float weight_read; int edge_id;
                    ss >> edge_id >> sourceID >> targetID >> weight_read;
                    int weight_scaled = (int )(weight_read*scale_edge);
                    this->insertEdge(sourceID,targetID,weight_scaled);
#else
                    float weight_read; int edge_id;
                    ss >> edge_id >> sourceID >> targetID >> weight_read;
                    int weight_scaled = (int )(weight_read*scale_edge);
                    this->insertEdge(sourceID,targetID,weight_scaled);
#endif
                }
                // The last element firstEdgeIndex is the size of the edges vector
                this->firstEdgeIndex.push_back(this->edges.size());
            } else {
                std::cerr << "Cannot open graph file " << grFilePath << std::endl;
                exit(1);
            }
            graphFS.close();

            std::cout << "\nSuccessfully parsed graph " << this->networkName << " with " << this->numNodes
                      << " nodes and " << this->numEdges << " edges!" << std::endl;
            //getchar();
            std::cout << "this->edges.size="<<this->edges.size()<<endl;
            std::cout << "this->getNumEdges="<<this->getNumEdges()<<endl;
            std::cout << "firstEdgeIndex.size="<<this->firstEdgeIndex.size()<<endl;
            std::cout << "this->getNumNodes="<<this->getNumNodes()<<endl;
            std::cout << "this->coordinates.size="<<this->coordinates.size()<<endl;

            //assert(this->edges.size() == static_cast<unsigned int>(this->getNumEdges()) && "Number of edges added not equal to expected number");
            assert(this->firstEdgeIndex.size() == static_cast<unsigned int>(this->getNumNodes()+1) && "Number of nodes added not equal to expected number");
            assert(this->coordinates.size() == static_cast<unsigned int>(this->getNumNodes()) && "Number of coordinate pairs added not equal to expected number");
        }




        std::vector<NodeID> getNodesIDsVector();
        std::unordered_set<NodeID> getNodesIDsUset();
        Coordinate getMaxCoordinateRange();
        void getMinMaxCoordinates(Coordinate& minX, Coordinate& maxX, Coordinate& minY, Coordinate& maxY);
        NodeID getNeighbourAndEdgeWeightByPosition(NodeID u, EdgeID pos, EdgeWeight& edgeWeight);
        NodeID getNeighbourAndEdgeWeightByPosition(NodeID u, EdgeID pos, EdgeWeight& edgeWeight, NodeID& edgeIndex);
        void getCoordinates(NodeID u, Coordinate& x, Coordinate& y);
        void getTranslatedCoordinates(NodeID u, Coordinate& x, Coordinate& y, int xTranslation, int yTranslation);
        EuclideanDistance getEuclideanDistance(NodeID node1, NodeID node2);
        int getEdgeListStartIndex(NodeID node);
        int getEdgeListSize(NodeID node);
        
        // Graph Validation/Statistics Functions
        bool isConnected(NodeID node1, NodeID node2);
        bool isUndirectedGraph();
        bool isConnectedGraph();
        bool isContiguous();
        std::vector<int> countOutdegreeFrequencies();
        void printGraphStats();
        NodeID findRoadEnd(NodeID sourceNode, NodeID edgeIndex, int& outDegree, EdgeWeight& distanceToEnd);
        double getMaxGraphSpeedByEdge();
        double getMinGraphSpeedByEdge();
        
        // Convert Graph to Other Formats and Output
        void outputToTSVFile(std::string outputFilePath);
        void outputDIMACSFiles(std::string grFilePath, std::string coFilePath);
        void outputToDDSGFile(std::string outputFilePath);
        void outputZeroedCoordinatesToFile(std::string outputFilePath);
        Graph_RoadNetwork createSubGraph(std::string subgraphName, std::unordered_set<NodeID>& subgraphNodes);
        
        // INE Functions
        void parseObjectFile(std::string objectSetFilePath, std::string& setType, double& setDensity, int& setVariable, int& setSize);
        void parseObjectSet(std::vector<NodeID>& set);
        void resetAllObjects();
        bool isObject(NodeID nodeID);
        
        // Memory Usage Functions
        double computeIndexSize();
        double computeINEIndexSize();
        double computeMemoryUsage();
        
        // For each node in the graph (i.e. from 0 to numNodes), this vector
        // has the index of it's first edge in the edges graph. It also has one  The last
        // additional element with the size of the edges vector (so if the graph has
        // n nodes then this vector has n+1 elements).
        // Note: This is to stop us from falling off the edges vector when iterating.
        std::vector<unsigned int> firstEdgeIndex; 
        std::vector<NodeEdgeWeightPair> edges;
        std::vector<CoordinatePair> coordinates; // x is first, y is second (in pair)
        
    private:
        friend class boost::serialization::access;
        
        std::string networkName = "Unknown";
        unsigned int numNodes = 0;
        unsigned int numEdges = 0;
        std::string type = "dis";
        
        Coordinate minXCoord = 0, maxXCoord = 0, minYCoord = 0, maxYCoord = 0;
        NodeID minXNode = 0, maxXNode = 0, minYNode = 0, maxYNode = 0;
        double maxGraphSpeedByEdge = 0, minGraphSpeedByEdge = std::numeric_limits<double>::max();

        // Non-Serialized Members
        std::vector<bool> nodeObjects; // Used at run-time only
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->networkName;
            ar & this->numNodes;
            ar & this->numEdges;
            ar & this->type;
            ar & this->firstEdgeIndex;
            ar & this->edges;
            ar & this->coordinates;
            ar & this->minXCoord;
            ar & this->maxXCoord;
            ar & this->minYCoord;
            ar & this->maxYCoord;
            ar & this->minXNode;
            ar & this->maxXNode;
            ar & this->minYNode;
            ar & this->maxYNode;
            //ar & this->maxGraphSpeedByEdge;
            //ar & this->minGraphSpeedByEdge;
        }    
        
};

#endif // _GRAPH_H
