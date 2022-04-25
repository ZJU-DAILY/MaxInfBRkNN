/* Original Work Copyright (C) 2012 Lingkun Wu, Xiaokui Xiao, Dingxiong Deng, Gao Cong, Andy Diwen Zhu, Shuigeng Zhou
 * Modified Work Copyright (C) 2015 Tenindra Abeywickrama
 *
 * This file is part of Road Network kNN Experimental Evaluation.
 * 
 * The following file is a derivative work of the code from 
 * http://sourceforge.net/projects/ntu-sp-exp/ developed for the paper below.
 * The authors have requested users of this code to cite the paper.
 * 
 * Lingkun Wu, Xiaokui Xiao, Dingxiong Deng, Gao Cong, Andy Diwen Zhu, Shuigeng Zhou
 * Shortest Path and Distance Queries on Road Networks: An Experimental Evaluation
 * PVLDB 5(5): 406-417 (2012)
 *
 * That work is in turn a derivative work of the code from
 * Contraction Hierarchies (http://algo2.iti.kit.edu/english/routeplanning.php),
 * licensed under AGPLv3. Thus this work is also licensed under that license.
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

#ifndef _QUADTREE_H
#define _QUADTREE_H

#include "Graph.h"
#include "../common.h"

#include <vector>
#include <unordered_map>
#include <deque>
#include <limits>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>

enum Direction{
    SW, SE, NW, NE
};

namespace quadtree_constants {
    const int NUM_QUADRANTS = 4;
    const int MAX_COORDINATE_LENGTH = 32; // Since it is a signed int
    const int MAX_MORTON_LENGTH = 64; // Since it is a signed int
    const int QUADRANT_MASK_POSITIONS = 63; // Needs 2 digits in 64-bit z-number for child's quadrant
                                            // so range is 0 to 62 (64 - 2) so 63 combinations
    const EdgeID EMPTY_COLOUR = static_cast<EdgeID>(constants::MAX_EDGES);
    // Note: We check in graph parsing that the number of adjacent nodes
    // for any node is smaller MAX_EDGES-1 to ensure SPECIAL_COLOUR won't clash
    const EdgeID SPECIAL_COLOUR = static_cast<EdgeID>(constants::MAX_EDGES-1);

    const unsigned long long B[] = {0x5555555555555555, 0x3333333333333333, 0x0F0F0F0F0F0F0F0F, 
                                    0x00FF00FF00FF00FF, 0x0000FFFF0000FFFF, 0x00000000FFFFFFFF};
    const unsigned long long S[] = {1ULL, 2ULL, 4ULL, 8ULL, 16ULL};

    const unsigned long long SW_MC = 0ULL; // 00
    const unsigned long long SE_MC = 1ULL; // 01
    const unsigned long long NW_MC = 2ULL; // 10
    const unsigned long long NE_MC = 3ULL; // 11

    // The following four arrays contain the all the possible positions of
    // quadrant codes within a 64-bit z-number
    const unsigned long long ARRAY_SW_BIT_MASK[QUADRANT_MASK_POSITIONS] = {
        SW_MC << 62, SW_MC << 61, SW_MC << 60, SW_MC << 59, SW_MC << 58, 
        SW_MC << 57, SW_MC << 56, SW_MC << 55, SW_MC << 54, SW_MC << 53, 
        SW_MC << 52, SW_MC << 51, SW_MC << 50, SW_MC << 49, SW_MC << 48,
        SW_MC << 47, SW_MC << 46, SW_MC << 45, SW_MC << 44, SW_MC << 43, 
        SW_MC << 42, SW_MC << 41, SW_MC << 40, SW_MC << 39, SW_MC << 38, 
        SW_MC << 37, SW_MC << 36, SW_MC << 35, SW_MC << 34, SW_MC << 33, 
        SW_MC << 32, SW_MC << 31, SW_MC << 30, SW_MC << 29, SW_MC << 28, 
        SW_MC << 27, SW_MC << 26, SW_MC << 25, SW_MC << 24, SW_MC << 23, 
        SW_MC << 22, SW_MC << 21, SW_MC << 20, SW_MC << 19, SW_MC << 18, 
        SW_MC << 17, SW_MC << 16, SW_MC << 15, SW_MC << 14, SW_MC << 13, 
        SW_MC << 12, SW_MC << 11, SW_MC << 10, SW_MC << 9, SW_MC << 8, 
        SW_MC << 7, SW_MC << 6, SW_MC << 5, SW_MC << 4, SW_MC << 3, 
        SW_MC << 2, SW_MC << 1, SW_MC << 0
    };

    const unsigned long long ARRAY_SE_BIT_MASK[QUADRANT_MASK_POSITIONS] = {
        SE_MC << 62, SE_MC << 61, SE_MC << 60, SE_MC << 59, SE_MC << 58, 
        SE_MC << 57, SE_MC << 56, SE_MC << 55, SE_MC << 54, SE_MC << 53, 
        SE_MC << 52, SE_MC << 51, SE_MC << 50, SE_MC << 49, SE_MC << 48,
        SE_MC << 47, SE_MC << 46, SE_MC << 45, SE_MC << 44, SE_MC << 43, 
        SE_MC << 42, SE_MC << 41, SE_MC << 40, SE_MC << 39, SE_MC << 38, 
        SE_MC << 37, SE_MC << 36, SE_MC << 35, SE_MC << 34, SE_MC << 33, 
        SE_MC << 32, SE_MC << 31, SE_MC << 30, SE_MC << 29, SE_MC << 28, 
        SE_MC << 27, SE_MC << 26, SE_MC << 25, SE_MC << 24, SE_MC << 23, 
        SE_MC << 22, SE_MC << 21, SE_MC << 20, SE_MC << 19, SE_MC << 18, 
        SE_MC << 17, SE_MC << 16, SE_MC << 15, SE_MC << 14, SE_MC << 13, 
        SE_MC << 12, SE_MC << 11, SE_MC << 10, SE_MC << 9, SE_MC << 8, 
        SE_MC << 7, SE_MC << 6, SE_MC << 5, SE_MC << 4, SE_MC << 3, 
        SE_MC << 2, SE_MC << 1, SE_MC << 0
    };

    const unsigned long long ARRAY_NW_BIT_MASK[QUADRANT_MASK_POSITIONS] = {
        NW_MC << 62, NW_MC << 61, NW_MC << 60, NW_MC << 59, NW_MC << 58, 
        NW_MC << 57, NW_MC << 56, NW_MC << 55, NW_MC << 54, NW_MC << 53, 
        NW_MC << 52, NW_MC << 51, NW_MC << 50, NW_MC << 49, NW_MC << 48,
        NW_MC << 47, NW_MC << 46, NW_MC << 45, NW_MC << 44, NW_MC << 43, 
        NW_MC << 42, NW_MC << 41, NW_MC << 40, NW_MC << 39, NW_MC << 38, 
        NW_MC << 37, NW_MC << 36, NW_MC << 35, NW_MC << 34, NW_MC << 33, 
        NW_MC << 32, NW_MC << 31, NW_MC << 30, NW_MC << 29, NW_MC << 28, 
        NW_MC << 27, NW_MC << 26, NW_MC << 25, NW_MC << 24, NW_MC << 23, 
        NW_MC << 22, NW_MC << 21, NW_MC << 20, NW_MC << 19, NW_MC << 18, 
        NW_MC << 17, NW_MC << 16, NW_MC << 15, NW_MC << 14, NW_MC << 13, 
        NW_MC << 12, NW_MC << 11, NW_MC << 10, NW_MC << 9, NW_MC << 8, 
        NW_MC << 7, NW_MC << 6, NW_MC << 5, NW_MC << 4, NW_MC << 3, 
        NW_MC << 2, NW_MC << 1, NW_MC << 0
    };

    const unsigned long long ARRAY_NE_BIT_MASK[QUADRANT_MASK_POSITIONS] = {
        NE_MC << 62, NE_MC << 61, NE_MC << 60, NE_MC << 59, NE_MC << 58, 
        NE_MC << 57, NE_MC << 56, NE_MC << 55, NE_MC << 54, NE_MC << 53, 
        NE_MC << 52, NE_MC << 51, NE_MC << 50, NE_MC << 49, NE_MC << 48,
        NE_MC << 47, NE_MC << 46, NE_MC << 45, NE_MC << 44, NE_MC << 43, 
        NE_MC << 42, NE_MC << 41, NE_MC << 40, NE_MC << 39, NE_MC << 38, 
        NE_MC << 37, NE_MC << 36, NE_MC << 35, NE_MC << 34, NE_MC << 33, 
        NE_MC << 32, NE_MC << 31, NE_MC << 30, NE_MC << 29, NE_MC << 28, 
        NE_MC << 27, NE_MC << 26, NE_MC << 25, NE_MC << 24, NE_MC << 23, 
        NE_MC << 22, NE_MC << 21, NE_MC << 20, NE_MC << 19, NE_MC << 18, 
        NE_MC << 17, NE_MC << 16, NE_MC << 15, NE_MC << 14, NE_MC << 13, 
        NE_MC << 12, NE_MC << 11, NE_MC << 10, NE_MC << 9, NE_MC << 8, 
        NE_MC << 7, NE_MC << 6, NE_MC << 5, NE_MC << 4, NE_MC << 3, 
        NE_MC << 2, NE_MC << 1, NE_MC << 0
    };
}

class MortonCode {

    public:
        MortonCode() {};
        MortonCode(Coordinate x, Coordinate y, int maxRangeExp); // For points
        MortonCode(Coordinate x, Coordinate y, int width, int maxRangeExp); // For blocks
        MortonCode(MortonNumber zNumber, MortonLength length);
        void setCode(const MortonCode& parentBlock, Direction dir);
        MortonNumber getZNumber() const;
        EdgeID getLength() const;
        bool contains(const MortonCode& point) const;
        int compareMortonCode(const MortonCode& point) const;
        Direction getQuadrant();
        MortonCode getParentCode();
        Direction getNextQuadrant(Direction quadrant);
        Direction getPrevQuadrant(Direction quadrant);
        void decompose(int maxRangeExp, Coordinate& x, Coordinate& y, int& width);
        bool operator==(const MortonCode& code) const;
        bool operator<(const MortonCode& code) const;
        bool operator>(const MortonCode& code) const;
        void printCode();
        MortonCode getNextMergeCode();
        MortonCode getPrevMergeCode();
        bool isFirstMergeCode();
        bool isLastMergeCode();
        Direction getNextQuadrant();
        Direction getPrevQuadrant();
        MortonNumber getChildZNumber(MortonNumber parentZNumber, Direction childQuadrant);
        double computeIndexSizeBytes();
        double computeMemoryUsageBytes();
        
    private:
        friend class boost::serialization::access;
        
        MortonNumber zNumber = 0;
        EdgeID len = 0;
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->zNumber;
            ar & this->len;
        }
};

class QuadtreeNode {

    public:
        QuadtreeNode() {};
        QuadtreeNode(MortonCode& parentCode, Direction dir, QuadtreeNode* prevLeaf = NULL, QuadtreeNode* nextLeaf = NULL);
        QuadtreeNode(bool isRoot, bool isLeaf);
        ~QuadtreeNode();
        void setNextLeaf(QuadtreeNode* nextLeaf);
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
        void copyDownNodesCheckColour(NodeID source, const std::vector<MortonCode>& allMortonCodes, std::vector<EdgeID>& colourMap, 
                                      std::vector<DistanceRatio>& distRatios, std::vector<QuadtreeNode*>& treeNodes, 
                                      std::deque<bool>& sameColourStatus, std::vector<EdgeID>& colours, 
                                      std::vector<DistanceRatio>& minBlockRatios, std::vector<DistanceRatio>& maxBlockRatios,
                                      std::vector<std::vector<DistanceRatio>>& childDistRatios);
        bool contains(const MortonCode& point);
        const MortonCode& getMortonCode() const;
        double computeIndexSizeBytes();
        double computeMemoryUsageBytes();
        
        // Child quadtree nodes (if non leaf)
        QuadtreeNode *sw, *se, *nw, *ne;
        
        // Doubly-linked list of leaf nodes
        QuadtreeNode *prevLeaf, *nextLeaf;
        
        MortonCode mortonCode;
        std::vector<NodeID> nodes; // Only full for leaf nodes

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
            ar & this->prevLeaf;
            ar & this->nextLeaf;
            ar & this->mortonCode;
            ar & this->nodes;
            ar & this->numNodes;
        }

};

class Quadtree {

    public:
        Quadtree() {};
        Quadtree(int maxLeafItems);
        ~Quadtree();
        void buildFromGraph(Graph_RoadNetwork& graph);
        int getMaxLeafItems();
        int getMaxRange();
        int getMaxRangeExp();
        int getXTranslation();
        int getYTranslation();
        void insert(NodeID node);
        void split(QuadtreeNode* treeNode);
        int countLeafNodes(int& totalNodes, int& nonEmptyNodes);
        const std::vector<MortonCode>& getMortonCodes() const; 
        void populateCoordinates(NodeID node, Coordinate& x, Coordinate& y);
        EuclideanDistance getEuclideanDistance(NodeID u, NodeID v);
        double computeIndexSize();
        double computeMemoryUsage();
        
        QuadtreeNode *rootNode;
        QuadtreeNode *firstLeaf, *lastLeaf;

        // We store a relative x and y coordinate for every node
        std::vector<Coordinate> relativeXCoordinates;
        std::vector<Coordinate> relativeYCoordinates;
        std::vector<MortonCode> mortonCodes;
        
        // Non-Serialized Members
        std::unordered_map<NodeID,int> nodeToIdxMap;
        
    private:
        friend class boost::serialization::access;
        
        int maxLeafItems;
        Coordinate maxRange;
        int maxRangeExp; // Min maxRangeExp such that 2^maxRangeExp > maxRange
        int xTranslation;
        int yTranslation;
        bool isPointSet;

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
            ar & this->relativeXCoordinates;
            ar & this->relativeYCoordinates;
            ar & this->mortonCodes;
            ar & this->rootNode;
            ar & this->firstLeaf;
            ar & this->lastLeaf;
        }
};

#endif // _QUADTREE_H