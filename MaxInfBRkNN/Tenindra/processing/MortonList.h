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

#ifndef _MORTONLIST_H
#define _MORTONLIST_H

#include "Graph.h"
#include "Quadtree.h"
#include "SimpleQuadtree.h"
#include "Junction.h"
#include "Path.h"
#include "StaticRtree.h"
#include "../queue/BinaryMinHeap.h"
#include "../queue/BinaryMaxHeapWithDK.h"
#include "../common.h"
#include "../utility/Statistics.h"

#include <list>
#include <bitset>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>

struct DataTuple {
    DistanceBound lowerBound;
    DistanceBound upperBound;
    NodeID intermediate;
    EdgeWeight distance; // From source to intermediate
    bool isBlock;
    SimpleQuadtreeNode* region;
    NodeID point;
    DataTuple(DistanceBound lowerBound, DistanceBound upperBound, NodeID intermediate, EdgeWeight distance, SimpleQuadtreeNode* region) : 
        lowerBound(lowerBound), upperBound(upperBound), intermediate(intermediate), distance(distance), isBlock(true), region(region) {}
    DataTuple(DistanceBound lowerBound, DistanceBound upperBound, NodeID intermediate, EdgeWeight distance, NodeID point) : 
        lowerBound(lowerBound), upperBound(upperBound), intermediate(intermediate), distance(distance), isBlock(false), point(point) {}
};

class MortonBlock {

    public:
        MortonBlock() {};
        MortonBlock(MortonCode code, EdgeID link);
        MortonBlock(MortonNumber zNumber, EdgeID mortonLen, EdgeID link);
        EdgeID getLink();
        void setDistanceInterval(DistanceRatio minDistRatio, DistanceRatio maxDistRatio);
        DistanceRatio getMinDistRatio();
        DistanceRatio getMaxDistRatio();
        void printBlock();
        void decomposeAndGetRatios(int maxRangeExp, Coordinate& x, Coordinate& y, int& width, DistanceRatio& min, DistanceRatio& max);
        EdgeID getLinkAndRatios(DistanceRatio& min, DistanceRatio& max);
        MortonNumber getMortonNumber();
        EdgeID getMortonLength();
        int compareMortonCode(const MortonCode& code) const;
        MortonBlock createMergeBlock(EdgeID mergeColour);
        MortonCode getMortonCode();
        MortonCode getNextMergeCode();
        MortonCode getPrevMergeCode();
        bool isFirstMergeBlock();
        bool isLastMergeBlock();
        Direction getNextQuadrant();
        Direction getPrevQuadrant();
        MortonNumber getChildZNumber(MortonNumber parentZNumber, Direction childQuadrant);
        void decompose(int maxRangeExp, Coordinate& x, Coordinate& y, int& width);
        bool operator==(const MortonCode& code) const;
        double computeIndexSizeBytes();
        double computeMemoryUsageBytes();
        
    private:
        friend class boost::serialization::access;
        
        MortonNumber zNumber;
        EdgeID len; // zNumber length (i.e. used bits)
        EdgeID link;
        DistanceRatio minDistRatio;
        DistanceRatio maxDistRatio;
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->zNumber;
            ar & this->len;
            ar & this->link;
            ar & this->minDistRatio;
            ar & this->maxDistRatio;
        }        

};

class MortonList {

    public:
        MortonList() {};
        void Mortonize(NodeID source, Quadtree& qt, std::vector<EdgeID>& colourMap, std::vector<EdgeWeight>& SSSPDistances, bool singleColour);
        void mergeAndSearchBackwards(std::vector<MortonBlock>& mortonList, bool& mergeColourSet, EdgeID& mergeColour, std::vector<MortonBlock>& mergeCandidates, MortonCode& nextMergeCode);
        void addAndClearMergeCandidates(std::vector<MortonBlock>& mortonList, bool& mergeColourSet, EdgeID& mergeColour, std::vector<MortonBlock>& mergeCandidates);
        EdgeID getLink(const MortonCode& targetCode);
        EdgeID getLinkAndRatios(const MortonCode& targetCode, DistanceRatio& minRatio, DistanceRatio& maxRatio);
        std::size_t getBlockIdxContainingTarget(const MortonCode& targetCode);
        int getNumBlocks();
        bool isListSorted();
        bool isListMergeable(int& count);
        bool areBlockIntervalsValid(NodeID source, Quadtree& qt, std::vector<EdgeWeight>& SSSPDistances);
        void splitAndAddToList(NodeID source, Quadtree& qt, QuadtreeNode* splitNode, std::vector<EdgeID>& colourMap, 
                               std::vector<MortonBlock>& mortonList, std::vector<DistanceRatio>& nodeDistRatios);
        std::size_t getBlockIdxContainingTargetByBST(const MortonCode& targetCode);
        void getIntervalDistance(Coordinate x, Coordinate y, Coordinate regionLeft, Coordinate regionBottom, int regionWidth, 
                                 int maxRangeExp, int maxRange, DistanceBound& lowerBound, DistanceBound& upperBound);
        void rangeSearch(std::size_t parent, MortonCode& minCode, MortonCode& maxCode, int maxRangeExp, Coordinate x, Coordinate y, Coordinate regionLeft, Coordinate regionBottom, 
                         Coordinate regionRight, Coordinate regionTop, DistanceBound& minLowerBound, DistanceBound& maxUpperBound, bool& initialized);
        std::size_t getLeftChild(std::size_t parent);
        std::size_t getRightChild(std::size_t parent);
        bool isTreeSorted(int pos);
        bool areSubtreesValid(int pos);
        bool isLeftSubtreeValid(int leftChild, MortonCode& parentCode);
        bool isRightSubtreeValid(int rightChild, MortonCode& parentCode);
        bool unionBoundsWithCandidateBlock(int blockIdx, int maxRangeExp, Coordinate x, Coordinate y, 
                                           Coordinate regionLeft, Coordinate regionBottom, Coordinate regionRight, Coordinate regionTop, 
                                           DistanceBound& minLowerBound, DistanceBound& maxUpperBound, bool& initialized);
        void computeNewSearchBoundaries(MortonCode& minCode, MortonCode& maxCode, MortonBlock& dividingBlock, MortonNumber& leftMaxZ, MortonNumber& rightMinZ, int maxRangeExp);
        void loadPattern(unsigned int firstVal, unsigned int repeatVal, int pos, std::bitset<64>& zNumber, int maxRangeExp);
        bool checkLeftChildOverlapping(int parentIdx, int maxRangeExp, Coordinate regionLeft, Coordinate regionBottom, Coordinate regionRight, Coordinate regionTop);
        bool checkRightChildOverlapping(int parentIdx, int maxRangeExp, Coordinate regionLeft, Coordinate regionBottom, Coordinate regionRight, Coordinate regionTop);
        
        // Distance Browsing Versions of Functions
        bool unionBoundsWithCandidateBlockForDistBrws(int blockIdx, int maxRangeExp, Coordinate x, Coordinate y, 
                                                    Coordinate regionLeft, Coordinate regionBottom, Coordinate regionRight, Coordinate regionTop,
                                                    DistanceBound& minLowerBound, bool& initialized);   
        void rangeSearchForDistBrws(std::size_t parent, MortonCode& minCode, MortonCode& maxCode, int maxRangeExp, Coordinate x, Coordinate y, Coordinate regionLeft, Coordinate regionBottom, 
                         Coordinate regionRight, Coordinate regionTop, DistanceBound& minLowerBound, bool& initialized);
        void getIntervalDistanceForDistBrws(Coordinate x, Coordinate y, Coordinate regionLeft, Coordinate regionBottom, int regionWidth, 
                                 int maxRangeExp, int maxRange, DistanceBound& lowerBound);
        void printBST();
        void printMortonList();
        double computeIndexSizeBytes();
        double computeMemoryUsageBytes();

    private:
        friend class boost::serialization::access;
        
        std::vector<MortonBlock> binarySearchTreeMortonList;
        
        // Non-Serialized
        std::vector<MortonBlock> mortonList;
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->binarySearchTreeMortonList;
        }

};

class SILCPathOracle {
 
    public:
        SILCPathOracle() {};
        void buildPathOracle(Quadtree& qt, Graph_RoadNetwork& graph);
        Path findShortestPath(Graph_RoadNetwork& graph, NodeID source, NodeID target);
        EdgeWeight findShortestPathDistance(Graph_RoadNetwork& graph, NodeID source, NodeID target);
        Path findShortestPathOptimised(Graph_RoadNetwork& graph, NodeID source, NodeID target);
        EdgeWeight findShortestPathDistanceOptimised(Graph_RoadNetwork& graph, NodeID source, NodeID target);
        NodeID nextInPath(Graph_RoadNetwork& graph, NodeID intermediate, MortonCode& targetCode, EdgeWeight& edgeWeight);
        void refineDistance(Graph_RoadNetwork& graph, NodeID target, MortonCode& targetCode,
                            NodeID& next, EdgeWeight& sourceToNextDistance,
                            DistanceBound &lowerBound, DistanceBound &upperBound);
        void optimisedRefineDistance(Graph_RoadNetwork& graph, NodeID target, MortonCode& targetCode, 
                            NodeID& next, EdgeWeight& sourceToNextDistance,
                            DistanceBound &lowerBound, DistanceBound &upperBound, Junction& junc);
        void intervalDistance(NodeID source, SimpleQuadtree& qt, SimpleQuadtreeNode* region, 
                              DistanceBound& lowerBound, DistanceBound& upperBound);
        void intervalDistanceForDistBrws(NodeID source, SimpleQuadtree& qt, SimpleQuadtreeNode* region, 
                              DistanceBound& lowerBound);
        std::string getNetworkName();
        int getNumNodes();
        int getNumEdges();
        void getKNNs(SimpleQuadtree& objectHierarchy, unsigned int k, NodeID queryNode, std::vector<NodeID>& kNNs, 
                     std::vector<EdgeWeight>& kNNDistances, Graph_RoadNetwork& graph);
        void getKNNsByDistanceBrowsing(SimpleQuadtree& objectHierarchy, unsigned int k, NodeID queryNode, std::vector<NodeID>& kNNs, 
                                       std::vector<EdgeWeight>& kNNDistances, Graph_RoadNetwork& graph);
        void getKNNsByOptimisedSILC(SimpleQuadtree& objectHierarchy, unsigned int k, NodeID queryNode, std::vector<NodeID>& kNNs, 
                                    std::vector<EdgeWeight>& kNNDistances, Graph_RoadNetwork& graph, Junction& junc);
        void getKNNsByOptimisedDistanceBrowsing(SimpleQuadtree& objectHierarchy, unsigned int k, NodeID queryNode, std::vector<NodeID>& kNNs, 
                                                std::vector<EdgeWeight>& kNNDistances, Graph_RoadNetwork& graph, Junction& junc);
        void getKNNsByDistanceBrowsingViaRtree(StaticRtree& rtree, unsigned int k, NodeID queryNode, std::vector<NodeID>& kNNs, 
                                    std::vector<EdgeWeight>& kNNDistances, Graph_RoadNetwork& graph, Junction& junc);
        MortonCode convertToMortonCode(Graph_RoadNetwork& graph, NodeID node);
        void findRealKNNDistances(Graph_RoadNetwork& graph, NodeID queryNode, std::vector<NodeID>& kNNs, std::vector<EdgeWeight>& kNNDistances);
        double computeIndexSize();
        double computeMemoryUsage();
        
        std::vector<Coordinate> relativeXCoordinates;
        std::vector<Coordinate> relativeYCoordinates;
        std::vector<MortonCode> mortonCodes;
        Statistics stats;
        
    private:
        friend class boost::serialization::access;
        
        std::vector<MortonList> pathOracle;
        std::string networkName;
        unsigned int numNodes;
        unsigned int numEdges;
        Coordinate maxRange;
        int maxRangeExp;        
        int xTranslation;
        int yTranslation;
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->pathOracle;
            ar & this->networkName;
            ar & this->numNodes;
            ar & this->numEdges;
            ar & this->maxRange;
            ar & this->maxRangeExp;
            ar & this->xTranslation;
            ar & this->yTranslation;
            ar & this->relativeXCoordinates;
            ar & this->relativeYCoordinates;
            ar & this->mortonCodes;
        }           

};

#endif // _MORTONLIST_H