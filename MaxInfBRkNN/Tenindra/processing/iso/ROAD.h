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

#ifndef _ROAD_H
#define _ROAD_H

#include "Graph.h"
#include "../queue/BinaryMinHeap.h"
#include "../utility/METISWrapper.h"
#include "../utility/Statistics.h"

#include <string>
#include <unordered_map>
#include <map>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_set.hpp>

class AssociationDirectory {

    public:
        AssociationDirectory(std::string setType, double setDensity, int setVariable, int setSize, int rnetTreeSize);
        AssociationDirectory() {};
        double getObjSetDensity();
        int getObjSetSize();
        int getObjSetVariable();
        std::string getObjSetType();
        bool isObject(NodeID node);
        bool hasObject(int rnetIdx);
        void addObject(NodeID node);
        void addRnetAssociation(int rnetIdx);
        double computeIndexSize();
        double computeMemoryUsage();
        
    private:
        friend class boost::serialization::access;
        
        std::string objSetType;
        double objSetDensity;
        int objSetVariable, objSetSize;
        
        std::unordered_set<NodeID> objectSet;
        std::vector<bool> rnetAssociationDirectory;
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->objSetType;
            ar & this->objSetDensity;
            ar & this->objSetVariable;
            ar & this->objSetSize;
            ar & this->objectSet;
            ar & this->rnetAssociationDirectory;
        }
    
};

class ShortcutTreeNode {
    public:
        ShortcutTreeNode(int rnetIdx, int rnetLevel);
        ShortcutTreeNode() {}
        void addChild(int childIdx);
        const std::vector<int>& getChildren() const;
        bool hasBaseNode();
        int getBaseIdx();
        void setHasBaseNode();
        int getRnetIdx();
        void print();
        int getLevel();
        bool isBaseNode();
        double computeIndexSizeBytes();
        double computeMemoryUsageBytes();
        void setShortcutStartIdx(int shortcutStartIdx);
        void setShortcutEndPlusOneIdx(int shortcutEndPlusOneIdx);
        
        int shortcutStartIdx;
        int shortcutEndPlusOneIdx;
        std::vector<int> children;
        // Note: When the ShortcutTreeNode is a leaf Rnet then it's child 
        // will be a base ShortcutTreeNode with the border's original edges 
        // contained in the leaf Rnet
        
    private:
        friend class boost::serialization::access;
        
        int rnetIdx; // Note: This must be -1 for base nodes
        int rnetLevel; // Note: This must be -1 for base nodes
        bool hasBase;
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->rnetIdx;
            ar & this->rnetLevel;
            ar & this->hasBase;
            ar & this->children;
            ar & this->shortcutStartIdx;
            ar & this->shortcutEndPlusOneIdx;
        }
};

class RouteOverlayNode {
    public:
        RouteOverlayNode(NodeID id);
        RouteOverlayNode() {}
        void addRnetToShortcutTree(int rnetIdx, int parentRnetIdx, int rnetLevel);
        void addLeafRnet(int rnetIdx);
        void addBaseToShortcutTree(int rnetIdx);
        bool hasBaseNode(int rnetIdx);
        void print();
        const std::vector<int>& getLeafRnets() const;
        const std::vector<int>& getShortcutTreeChildrenByRnetIdx(int rnetIdx);
        const std::vector<int>& getShortcutTreeChildren(int scTreeIdx) const;
        std::vector<int> getShortcutsIdxsAtLevel(int level);
        int getRnetIdx(int scTreeIdx);
        bool hasRnet(int rnetIdx);
        const std::vector<ShortcutTreeNode>& getShortcutTree() const;
        double computeIndexSizeBytes();
        double computeMemoryUsageBytes();
        void setShortcutStartIdx(int rnetIdx, int shortcutStartIdx);
        void setShortcutEndPlusOneIdx(int rnetIdx, int shortcutEndPlusOneIdx);
        void setBaseEdgeStartIdx(int leafRnetIdx, int shortcutStartIdx);
        void setBaseEdgeEndPlusOneIdx(int leafRnetIdx, int shortcutEndPlusOneIdx);
        
        std::vector<int> leafIdxs; // This is a vector because leaf Rnet border could belong to upto fanout leaf Rnets
        std::vector<ShortcutTreeNode> shortcutTree;
        
        // Non-Serialized Members
        // When building add a new rnet to the shortcut tree we need to know 
        // the shortcut tree index of the parent rnet to allow traversal
        // Note: This doesn't have to be serialized, we search tree from
        // root depth-first when querying
        //std::unordered_map<int,int> rnetIdxToScTreeIdx;
        std::map<int,int> rnetIdxToScTreeIdx;
        // Note: std::map is faster for small numbers of elements

    private:
        friend class boost::serialization::access;
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->leafIdxs;
            ar & this->shortcutTree;
        }                 
};

class RnetNode {
    public:
        RnetNode(int parentRnetIdx, int level);
        RnetNode() {}
        bool isLeafRnet();
        void setLeafRnet();
        void addChild(int childRnetIdx);
        int getLevel();
        void addBorder(NodeID node);
        int getParent();
        const std::vector<int>& getChildren() const;
        const std::vector<NodeID>& getBorders() const;
        const std::unordered_set<NodeID>& getBordersUset() const;
        void print();
        bool isBorder(NodeID node);
        double computeIndexSizeBytes();
        double computeMemoryUsageBytes();
        int getNumBorders();
        
        std::vector<int> children;
        
        // Non-Serialized Members
        std::vector<NodeID> borders;
        std::unordered_set<NodeID> bordersUset;

    private:
        friend class boost::serialization::access;
        
        int parentRnetIdx;
        int level;
        bool isLeaf;
        int numBorders;

        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->parentRnetIdx;
            ar & this->level;
            ar & this->isLeaf;
            ar & this->children;
            //ar & this->borders;
            //ar & this->bordersUset;
        }                 
};

class ROAD {

    public:
        ROAD(std::string networkName, int numNodes, int numEdges, int fanout, int levels);
        ROAD() {};
        void buildRouteOverlay(Graph_RoadNetwork& graph);
        void buildShortcutTrees(Graph_RoadNetwork& graph);
        std::string getNetworkName();
        int getNumNodes();
        int getNumEdges();
        int getFanout();
        int getNumLevels();
        int getNumBorders();
        int getBorderToBorderRelationships();
        void printLevels();
        int getRealNumLevels();
        int getRnetTreeSize();
        std::vector<std::vector<int>> getRnetIdxsByLevel();
        const std::vector<int>& getLeafRnets(NodeID node) const;
        int getParentRnet(int rnetIdx);
        void partitionRnet(int parentIdx, int parentLevel, std::unordered_set<NodeID>& subgraph, Graph_RoadNetwork& originalGraph, METISWrapper& metis);
        void findSSMTDistances(Graph_RoadNetwork& graph, NodeID source,
                               std::unordered_set<NodeID>& targetSet, 
                               std::unordered_map<NodeID,EdgeWeight>& results,
                               BinaryMinHeap<EdgeWeight,NodeStatusPair> *pqueue);
        void findSSMTDistancesByRnet(int rnetIdx, NodeID source, std::unordered_set<NodeID>& targetSet, 
                                     std::unordered_map<NodeID,EdgeWeight>& results,
                                     BinaryMinHeap<EdgeWeight,NodeStatusPair> *pqueue);
        void getKNNs(AssociationDirectory& assocDir, unsigned int k, NodeID queryNodeID, std::vector<NodeID>& kNNs, 
                     std::vector<EdgeWeight>& kNNDistances);
        double computeIndexSize();
        double computeMemoryUsage();
        void addShortcut(NodeID target, EdgeWeight distance);
        
        std::vector<RouteOverlayNode> routeOverlay;
        std::vector<RnetNode> rnetTree;
        std::vector<NodeEdgeWeightPair> shortcutEdges;
        std::vector<bool> borders;
        Statistics stats;
        
    private:
        friend class boost::serialization::access;
        
        std::string networkName;
        int numNodes;
        int numEdges;
        int fanout;
        int levels;
        
        // Non-Serialized Members
        std::vector<NodeID> METISIdxToNodeID;
#if defined(COLLECT_STATISTICS)
        std::vector<int> rnetNodeCounts;
#endif
    
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->networkName;
            ar & this->numNodes;
            ar & this->numEdges;
            ar & this->fanout;
            ar & this->levels;
            ar & this->routeOverlay;
            ar & this->rnetTree;
            ar & this->shortcutEdges;
            ar & this->borders;
        }    
};


#endif // _ROAD_H