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

#ifndef _SIMPLEQUADTREE_H
#define _SIMPLEQUADTREE_H

#include "Quadtree.h"

class SimpleQuadtreeNode {

    public:
        SimpleQuadtreeNode() {};
        SimpleQuadtreeNode(MortonCode& parentCode, Direction dir);
        SimpleQuadtreeNode(bool isRoot, bool isLeaf);
        ~SimpleQuadtreeNode();
        bool isRootNode();
        bool isLeafNode();
        void setLeaf(bool isLeaf);
        int getNumNodes(); // This wil get number > 0 even for non-leaf nodes
        const std::vector<NodeID>& getNodes() const;
        NodeID getFirstNode();
        void addNode(NodeID node);
        void pushDownNodes(std::vector<MortonCode>& allMortonCodes);
        void pushDownNodes(std::vector<MortonCode>& allMortonCodes, std::unordered_map<NodeID,int>& nodeToIdxMap);
        void pushDownNode(NodeID node, MortonCode& mortonCode);
        bool contains(const MortonCode& point);
        const MortonCode& getMortonCode() const;
        void getBlockCoordinates(Coordinate& x, Coordinate& y, int& width);
        void setBlockCoordinates(Coordinate x, Coordinate y, int width);
        void computeAndSetBlockCoordinates(int xTranslation, int yTranslation, int maxRangeExp);
        double computeIndexSizeBytes();
        double computeMemoryUsageBytes();
        int getTotalRegions();
        
        Coordinate graphX, graphY; // Coordinates of bottom left corner of region in original graph coordinates
        int width;
        
        // Child quadtree nodes (if non leaf)
        SimpleQuadtreeNode *sw, *se, *nw, *ne;
        
        MortonCode mortonCode;
        std::vector<NodeID> nodes; // Empty for non-leaf nodes

    private:
        friend class boost::serialization::access;
        
        bool isRoot;
        bool isLeaf;
        int numNodes;
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->isRoot;
            ar & this->isLeaf;
            ar & this->sw;
            ar & this->se;
            ar & this->nw;
            ar & this->ne;
            ar & this->mortonCode;
            ar & this->nodes;
            ar & this->numNodes;
            ar & this->graphX;
            ar & this->graphY;
            ar & this->width;
        }

};

class SimpleQuadtree {

    public:
        SimpleQuadtree() {};
        SimpleQuadtree(int maxLeafItems);
        ~SimpleQuadtree();
        void buildFromGraph(Graph_RoadNetwork& graph);
        void buildFromPointSet(Graph_RoadNetwork& graph, std::vector<NodeID>& pointSet, std::string setType, double setDensity, int setVariable, int setSize);
        int getMaxLeafItems();
        int getMaxRange();
        int getMaxRangeExp();
        int getXTranslation();
        int getYTranslation();
        void insert(NodeID node);
        void split(SimpleQuadtreeNode* treeNode);
        const std::vector<MortonCode>& getMortonCodes() const; 
        void populateCoordinates(NodeID node, Coordinate& x, Coordinate& y);
        EuclideanDistance getEuclideanDistance(NodeID u, NodeID v);
        double getObjSetDensity();
        int getObjSetSize();
        int getObjSetVariable();
        std::string getObjSetType();
        double computeIndexSize();
        double computeMemoryUsage();
        int getTotalRegions();
        
        SimpleQuadtreeNode *rootNode = NULL;

        // Non-Serialized Members
        // We store a relative x and y coordinate for every node
        std::vector<Coordinate> relativeXCoordinates;
        std::vector<Coordinate> relativeYCoordinates;
        std::vector<MortonCode> mortonCodes;
        std::unordered_map<NodeID,int> nodeToIdxMap;
        
    private:
        friend class boost::serialization::access;
        
        int maxLeafItems;
        Coordinate maxRange;
        int maxRangeExp; // Min maxRangeExp such that 2^maxRangeExp > maxRange
        int xTranslation;
        int yTranslation;
        bool isPointSet;
        std::string objSetType;
        double objSetDensity;
        int objSetSize, objSetVariable;

        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->maxLeafItems;
            ar & this->maxRange;
            ar & this->maxRangeExp;
            ar & this->xTranslation;
            ar & this->yTranslation;
            ar & this->isPointSet;
            ar & this->objSetType;
            ar & this->objSetDensity;
            ar & this->objSetVariable;
            ar & this->objSetSize;
            //ar & this->relativeXCoordinates;
            //ar & this->relativeYCoordinates;
            //ar & this->mortonCodes;
            //ar & this->nodeToIdxMap;
            ar & this->rootNode;
        }
};

#endif // _SIMPLEQUADTREE_H