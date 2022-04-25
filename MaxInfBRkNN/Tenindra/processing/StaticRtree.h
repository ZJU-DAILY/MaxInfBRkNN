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

#ifndef _STATICTREE_H
#define _STATICTREE_H

#include "../queue/BinaryMinHeap.h"
#include "../common.h"

#include <vector>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>

/* 
 * 2D In-Memory Static Rtree Implementation for Point Data
 */

class StaticRtreeNode;

struct RtreeDataTuple {
    bool isNode;
    StaticRtreeNode* nodePtr;
    NodeID objectID;
    RtreeDataTuple(StaticRtreeNode* nodePtr) : 
        isNode(true), nodePtr(nodePtr) {}
    RtreeDataTuple(NodeID objectID) : 
        isNode(false), objectID(objectID) {}
};

class Rectangle {

    public:
        Rectangle() {}
        void setDimensions(Coordinate left, Coordinate bottom, Coordinate right, Coordinate top);
        void setDimensions(Rectangle& rect);
        void setCentroid();
        void expand(Rectangle& rect);
        EuclideanDistanceSqrd getMinDistSqrd(Coordinate x, Coordinate y);
        EuclideanDistanceSqrd getMinMaxDistSqrd(Coordinate x, Coordinate y);
        EuclideanDistanceSqrd getObjectDistSqrd(Coordinate x, Coordinate y);
        void printRectangle();
        double computeIndexSizeBytes();
        double computeMemoryUsageBytes();
        bool intersects(Rectangle& rect);
        
        Coordinate left;
        Coordinate bottom;
        Coordinate right;
        Coordinate top;
        Coordinate centreX;
        Coordinate centreY;

    private:
        friend class boost::serialization::access;
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->left;
            ar & this->bottom;
            ar & this->right;
            ar & this->top;
            ar & this->centreX;
            ar & this->centreY;
        }
};

class StaticRtreeEntry {

    public:
        StaticRtreeEntry() {}
        void setEntryData(CoordinatePair& ptCoords, StaticRtreeNode* rtreeNode, NodeID id = 0);
        void setEntryData(Rectangle rect, StaticRtreeNode* rtreeNode, NodeID id = 0);
        static bool compareX(const StaticRtreeEntry& entry1, const StaticRtreeEntry& entry2);
        static bool compareY(const StaticRtreeEntry& entry1, const StaticRtreeEntry& entry2);
        void printEntry();
        double computeIndexSizeBytes();
        double computeMemoryUsageBytes();
        
        Rectangle mbr;
        StaticRtreeNode* rtreeNode;
        NodeID id; // Object ID for entries of leaf node

    private:
        friend class boost::serialization::access;

        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->mbr;
            ar & this->rtreeNode;
            ar & this->id;
        }
};


class StaticRtreeNode {

    public:
        StaticRtreeNode() {}
        StaticRtreeNode(int level);
        void addEntry(StaticRtreeEntry entry);
        bool isLeaf();
        void printNode();
        double computeIndexSizeBytes();
        double computeMemoryUsageBytes();
        
        std::vector<StaticRtreeEntry> entries;
        int level; // 0 for leaf nodes

    private:
        friend class boost::serialization::access;
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->entries;
            ar & this->level;
        }
};


class StaticRtree {

    public:
        StaticRtree() {}
        StaticRtree(int branchFactor);
        StaticRtree(int branchFactor, std::string setType, double setDensity, int setVariable, int setSize);
        ~StaticRtree();
        void bulkLoad(std::vector<NodeID>& objectIDs, std::vector<CoordinatePair>& objectCoords);
        StaticRtreeNode* buildLevel(int level, std::vector<StaticRtreeEntry>& levelEntries);
        void printTree();
        BinaryMinHeap<EuclideanDistanceSqrd,RtreeDataTuple> getKNNs(unsigned int k, Coordinate queryPointX, Coordinate queryPointY, 
                                                                    std::vector<NodeID>& kNNs, std::vector<EuclideanDistanceSqrd>& kNNDistancesSqrd);
        bool getNextNearestNeighbour(BinaryMinHeap<EuclideanDistanceSqrd,RtreeDataTuple>& heap, Coordinate queryPointX, 
                                     Coordinate queryPointY, NodeID& nn, EuclideanDistanceSqrd& nnDistSqrd);
        std::vector<NodeID> windowQuery(Rectangle& rect);
        std::vector<StaticRtreeNode*> getAllNodePtrs();
        double getObjSetDensity();
        int getObjSetSize();
        int getObjSetVariable();
        std::string getObjSetType();
        int getBranchFactor();
        double computeIndexSize();
        double computeMemoryUsage();
        
    private:
        friend class boost::serialization::access;
        
        StaticRtreeNode* root;
        int branchFactor;
        std::string objSetType;
        double objSetDensity;
        int objSetVariable,  objSetSize;

        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->root;
            ar & this->branchFactor;
            ar & this->objSetType;
            ar & this->objSetDensity;
            ar & this->objSetVariable;
            ar & this->objSetSize;
        }
};

#endif // _STATICTREE_H