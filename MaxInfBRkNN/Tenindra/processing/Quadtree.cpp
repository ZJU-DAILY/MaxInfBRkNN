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

#include "Quadtree.h"

#include "../utility/geometry.h"
#include "../utility/utility.h"

#include <bitset>
#include <cmath>

Quadtree::Quadtree(int maxLeafItems): maxLeafItems(maxLeafItems), maxRange(0), xTranslation(0), yTranslation(0), 
    isPointSet(false) {}

Quadtree::~Quadtree()
{
    // Recursively delete all quadtree nodes
    delete rootNode;
}

void Quadtree::buildFromGraph(Graph_RoadNetwork &graph) {
    
    Coordinate maxCoordinateRange = graph.getMaxCoordinateRange();

    // Find number bits needed to represent largest coordinate
    this->maxRange = 1; // 2^0
    this->maxRangeExp = 0;
    // The biggest number that maxRangeExp bits can reprsent is actually
    // maxRange-1. So if maxCoordinateRange is equal to maxRange, then at
    // maxRangeExp does have no bits to represent maxCoordinateRange so
    // we must continue loop in this case and hence the <= sign
    while (maxRange <= maxCoordinateRange) {
        this->maxRange = this->maxRange << 1;
        this->maxRangeExp++;
    }
    this->maxRange -= 1; 
    // Worst case there will be a node with x and y coordinates with 
    // maxRangeExp bits. This means when we interleave x and y bits,
    // in the worst case the max morton code length is 2*maxrangeExp
    
    this->rootNode = new QuadtreeNode(true,true);
    this->firstLeaf = this->rootNode;
    this->lastLeaf = this->rootNode;
    
    // We adjust the coordinate of the node with the minimum X (resp. Y) 
    // coordinate such that it has a coordinate of 0 and all other X (resp. Y)
    // coordinates are adjusted all adjust relative to this 
    
    NodeID numNodes = graph.getNumNodes();
    Coordinate x, y, minX, maxX, minY, maxY;
    graph.getMinMaxCoordinates(minX,maxX,minY,maxY);
    this->xTranslation = minX;
    this->yTranslation = minY;
    // Reserve space for all nodes
    this->mortonCodes.reserve(numNodes);
    this->relativeXCoordinates.reserve(numNodes);
    this->relativeYCoordinates.reserve(numNodes);
    for (NodeID node = 0; node < numNodes; ++node) {
        graph.getCoordinates(node,x,y);
        // We translate all X and Y coordinates such that (0,0) is the most
        // SW point and coordinates increase in the N and E directions
        x = x - this->xTranslation; // Translate right, e.g. -100 -(-100)
        y = y - this->yTranslation; 
        // Note: In order to have (0,0) be in the SW corner we adjust the X 
        // and Y coordinates. If coordinate is negative translate right
        // and if coordinate is position when translate left.
        //assert(x >= 0 && "Adjusted X coordinate is not positive");
        //assert(y >= 0 && "Adjusted Y coordinate is not positive");

        // Note: We store coordinates in the quadtree so that we don't require graph
        // again. We also pre-compute morton codes to make it easier to point location.
        // Both of these help improve SILC index building efficiency
        MortonCode code(x,y,maxRangeExp);
        this->mortonCodes.push_back(code);
        this->relativeXCoordinates.push_back(x);
        this->relativeYCoordinates.push_back(y);
        this->insert(node);
    }

//     // Debug Output
//     std::cout << "X Translation = " << this->xTranslation << std::endl;
//     std::cout << "Y Translation = " << this->yTranslation << std::endl;
//     std::cout << "Actual Range = " << maxCoordinateRange << std::endl;
//     std::cout << "Power of 2 Range = " << maxRange << std::endl;
//     std::cout << "Exponent = " << maxRangeExp << std::endl;
//     std::cout << "X Coords Size = " << this->mortonCodes.size() << std::endl;
//     std::cout << "Y Coords Size = " << this->relativeYCoordinates.size() << std::endl;
//     std::cout << "Morton Code Size = " << this->mortonCodes.size() << std::endl;
//     int vertexesInLeafNodes = 0, nonEmptyLeafNodes = 0;
//     std::cout << "Leaf Nodes = " << this->countLeafNodes(vertexesInLeafNodes, nonEmptyLeafNodes) << std::endl;
//     std::cout << "Nodes in Quadtree Leaf Nodes = " << vertexesInLeafNodes << std::endl;    
//     std::cout << "Non-Empty Leaf Nodes = " << nonEmptyLeafNodes << std::endl;
}

void Quadtree::insert(NodeID node)
{
    QuadtreeNode* currentNode = this->rootNode;
    
    int index;
    if (!this->isPointSet) {
        // If this is not a point set then we don't need to map 
        // between the internal index and the node ID
        index = node;
    } else {
        index = this->nodeToIdxMap[node];
    }

    // Find leaf node in which this node is contained
    while(!currentNode->isLeafNode()) {
        // Determine which of current (non-leaf) nodes children contains node
        if (currentNode->sw->contains(this->mortonCodes[index])) {
            currentNode = currentNode->sw;
        } else if (currentNode->se->contains(this->mortonCodes[index])) {
            currentNode = currentNode->se;
        } else if (currentNode->nw->contains(this->mortonCodes[index])) {
            currentNode = currentNode->nw;
        } else if (currentNode->ne->contains(this->mortonCodes[index])) {
            currentNode = currentNode->ne;
        } else {
            std::cerr << "No child quadrant of current node contains point" << std::endl;
            exit(1);
        }
    }
    
    currentNode->addNode(node);
    if (currentNode->getNumNodes() > this->maxLeafItems) {
        this->split(currentNode);
    }
}

void Quadtree::split(QuadtreeNode* treeNode)
{
    // Mark node is being no longer a leaf node
    treeNode->setLeaf(false);
    
    // Create and update tree node pointers for split
    treeNode->sw = new QuadtreeNode(treeNode->mortonCode,SW,treeNode->prevLeaf);
    treeNode->se = new QuadtreeNode(treeNode->mortonCode,SE,treeNode->sw);
    treeNode->nw = new QuadtreeNode(treeNode->mortonCode,NW,treeNode->se);
    treeNode->ne = new QuadtreeNode(treeNode->mortonCode,NE,treeNode->nw,treeNode->nextLeaf);
    
    // Note: We don't know the next leaf nodes until all of them have created
    treeNode->sw->setNextLeaf(treeNode->se);
    treeNode->se->setNextLeaf(treeNode->nw);
    treeNode->nw->setNextLeaf(treeNode->ne);
    
    if (treeNode->prevLeaf != NULL) {
        treeNode->prevLeaf->nextLeaf = treeNode->sw;
        treeNode->prevLeaf = NULL;
    }
    if (treeNode->nextLeaf != NULL) {
        treeNode->nextLeaf->prevLeaf = treeNode->ne;
        treeNode->nextLeaf = NULL;
    }
    
    if (!this->isPointSet) {
        // If this is not a point set then we don't need to map 
        // between the internal index and the node ID
        treeNode->pushDownNodes(this->mortonCodes);
    } else {
        treeNode->pushDownNodes(this->mortonCodes,this->nodeToIdxMap);
    }
    
    // If this node was the first or last leaf 
    // we need update these with the child nodes
    if (treeNode == this->firstLeaf) {
        this->firstLeaf = treeNode->sw;
    }
    if (treeNode == this->lastLeaf) {
        this->lastLeaf = treeNode->ne;
    }
    
    // It is possible that all of items were placed in a single quadrant
    // so we must check and recursively split those as well
    if (treeNode->sw->getNumNodes() > this->maxLeafItems) {
        this->split(treeNode->sw);
    }
    if (treeNode->se->getNumNodes() > this->maxLeafItems) {
        this->split(treeNode->se);
    }
    if (treeNode->nw->getNumNodes() > this->maxLeafItems) {
        this->split(treeNode->nw);
    }
    if (treeNode->ne->getNumNodes() > this->maxLeafItems) {
        this->split(treeNode->ne);
    }
}

int Quadtree::countLeafNodes(int& totalNodes, int& nonEmptyNodes) {
    totalNodes = 0;
    nonEmptyNodes = 0;
    int count = 0;
    QuadtreeNode* currentNode = this->firstLeaf;
    while (currentNode != NULL) {
        totalNodes += currentNode->getNumNodes();
        if (currentNode->getNumNodes() > 0) {
            ++nonEmptyNodes;
        }
        currentNode = currentNode->nextLeaf;
        ++count;
    }
    totalNodes = 0;
    nonEmptyNodes = 0;
    count = 0;
    currentNode = this->lastLeaf;
    while (currentNode != NULL) {
        totalNodes += currentNode->getNumNodes();
        if (currentNode->getNumNodes() > 0) {
            ++nonEmptyNodes;
        }
        currentNode = currentNode->prevLeaf;
        ++count;
    }
    return count;
}

int Quadtree::getMaxLeafItems()
{
    return this->maxLeafItems;
}

int Quadtree::getMaxRange() {
    return this->maxRange;
}

int Quadtree::getMaxRangeExp() {
    return this->maxRangeExp;
}

int Quadtree::getXTranslation() {
    return this->xTranslation;
}

int Quadtree::getYTranslation() {
    return this->yTranslation;
}

void Quadtree::populateCoordinates(NodeID node, Coordinate &x, Coordinate &y) {
    // Note: We cannot use this function if object has been deserialized
    assert(this->relativeXCoordinates.size() > 0 && "Relative X and Y coordinates do not exist (they are not serialized)");
    x = this->relativeXCoordinates[node];
    y = this->relativeYCoordinates[node];
}

EuclideanDistance Quadtree::getEuclideanDistance(NodeID s, NodeID t) {
    // Note: We cannot use this function if object has been deserialized
    assert(this->relativeXCoordinates.size() > 0 && "Relative X and Y coordinates do not exist (they are not serialized)");
    return geometry::getEuclideanDist(this->relativeXCoordinates[s],this->relativeYCoordinates[s],
                                      this->relativeXCoordinates[t],this->relativeYCoordinates[t]);
}

const std::vector<MortonCode> &Quadtree::getMortonCodes() const {
    return this->mortonCodes;
}

double Quadtree::computeIndexSize()
{
    double memoryUsage = 0;
    std::deque<QuadtreeNode*> stack;
    stack.push_back(this->rootNode);
    while (stack.size() > 0) {
        QuadtreeNode* currNode = stack.front();
        stack.pop_front();
        memoryUsage += currNode->computeIndexSizeBytes();
        if (!currNode->isLeafNode()) {
            stack.push_back(currNode->sw);
            stack.push_back(currNode->se);
            stack.push_back(currNode->nw);
            stack.push_back(currNode->ne);
        }
    }
    memoryUsage += sizeof(Coordinate)*relativeXCoordinates.size();
    memoryUsage += sizeof(Coordinate)*relativeYCoordinates.size();
    for (std::size_t i = 0; i < this->mortonCodes.size(); ++i) {
        memoryUsage += this->mortonCodes[i].computeIndexSizeBytes();
    }
    return memoryUsage/(1024*1024);
}

double Quadtree::computeMemoryUsage() {
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    std::deque<QuadtreeNode*> stack;
    stack.push_back(this->rootNode);
    while (stack.size() > 0) {
        QuadtreeNode* currNode = stack.front();
        stack.pop_front();
        memoryUsage += currNode->computeMemoryUsageBytes();
        if (!currNode->isLeafNode()) {
            stack.push_back(currNode->sw);
            stack.push_back(currNode->se);
            stack.push_back(currNode->nw);
            stack.push_back(currNode->ne);
        }
    }
    memoryUsage += sizeof(Coordinate)*relativeXCoordinates.capacity();
    memoryUsage += sizeof(Coordinate)*relativeYCoordinates.capacity();
    for (std::size_t i = 0; i < this->mortonCodes.size(); ++i) {
        memoryUsage += this->mortonCodes[i].computeMemoryUsageBytes();
    }
    memoryUsage += sizeof(MortonCode)*(this->mortonCodes.capacity()-this->mortonCodes.size());
    memoryUsage += utility::estimateUnorderedMapMemoryUsageBytes(this->nodeToIdxMap.size(),sizeof(std::pair<NodeID,int>),this->nodeToIdxMap.bucket_count());
    return memoryUsage/(1024*1024);
}

/*
 * QuadtreeNode
 */

QuadtreeNode::QuadtreeNode(MortonCode& parentCode, Direction dir, QuadtreeNode* prevLeaf /* = NULL */, QuadtreeNode* nextLeaf /* = NULL */)
{
    this->mortonCode.setCode(parentCode,dir);
    this->isRoot = false; // Cannot be root because it has a parent
    this->isLeaf = true;
    this->sw = this->se = this->nw = this->ne = NULL;
    this->prevLeaf = prevLeaf;
    this->nextLeaf = nextLeaf;
    this->numNodes = 0;
}

QuadtreeNode::QuadtreeNode(bool isRoot, bool isLeaf)
{
    this->isRoot = isRoot;
    this->isLeaf = isLeaf;
    this->sw = this->se = this->nw = this->ne = NULL;
    this->prevLeaf = NULL;
    this->nextLeaf = NULL;
    this->numNodes = 0;
}

QuadtreeNode::~QuadtreeNode() {
    // Recursively delete all quadtree nodes
    delete this->sw;
    delete this->se;
    delete this->nw;
    delete this->ne;
}

void QuadtreeNode::setNextLeaf(QuadtreeNode *nextLeaf) {
    this->nextLeaf = nextLeaf;
}

bool QuadtreeNode::isRootNode() {
    return this->isRoot;
}

int QuadtreeNode::getNumNodes()
{
    return this->numNodes;
}

bool QuadtreeNode::isLeafNode() {
    return this->isLeaf;
}

void QuadtreeNode::setLeaf(bool isLeaf) {
    this->isLeaf = isLeaf;
}

const std::vector<NodeID>& QuadtreeNode::getNodes() const {
    return this->nodes;
}

NodeID QuadtreeNode::getFirstNode()
{
    return this->nodes[0];
}

void QuadtreeNode::addNode(NodeID node) {
    this->numNodes++;
    this->nodes.push_back(node);
}

void QuadtreeNode::pushDownNodes(std::vector<MortonCode>& allMortonCodes) {
    for (std::size_t i = 0; i < this->nodes.size(); ++i) {
        this->pushDownNode(this->nodes[i],allMortonCodes[this->nodes[i]]);
    }
    this->nodes.clear();
}

void QuadtreeNode::pushDownNodes(std::vector<MortonCode>& allMortonCodes, std::unordered_map<NodeID,int>& nodeToIdxMap) {
    int index;
    for (std::size_t i = 0; i < this->nodes.size(); ++i) {
        index = nodeToIdxMap[this->nodes[i]];
        this->pushDownNode(this->nodes[i],allMortonCodes[index]);
    }
    this->nodes.clear();
}

void QuadtreeNode::pushDownNode(NodeID node, MortonCode& mortonCode) {
    if (this->sw->contains(mortonCode)) {
        this->sw->addNode(node);
    } else if (this->se->contains(mortonCode)) {
        this->se->addNode(node);
    } else if (this->nw->contains(mortonCode)) {
        this->nw->addNode(node);
    } else if (this->ne->contains(mortonCode)) {
        this->ne->addNode(node);
    } else {
        std::cerr << "Cannot push down, no child quadrant contains point" << std::endl;
        exit(1);
    }
}

void QuadtreeNode::copyDownNodesCheckColour(NodeID source, const std::vector<MortonCode>& allMortonCodes, std::vector<EdgeID>& colourMap, 
                                      std::vector<DistanceRatio>& distRatios, std::vector<QuadtreeNode*>& treeNodes, 
                                      std::deque<bool>& sameColourStatus, std::vector<EdgeID>& colours, 
                                      std::vector<DistanceRatio>& minBlockRatios, std::vector<DistanceRatio>& maxBlockRatios,
                                      std::vector<std::vector<DistanceRatio>>& childDistRatios) {
    //assert(distRatios.size() == this->nodes.size() && "Incorrect number of ratios");
    //assert(treeNodes.size() == 4 && "Incorrect number of tree nodes");
    
    std::deque<bool> empty;
    std::deque<bool> ratioIntialised;
    
    sameColourStatus.resize(4);
    colours.resize(4);
    minBlockRatios.resize(4);
    maxBlockRatios.resize(4);
    childDistRatios.resize(4);
    for (std::size_t i = 0; i < 4; ++i) {
        sameColourStatus[i] = true;
        colours[i] = quadtree_constants::EMPTY_COLOUR;
        empty.push_back(true);
        ratioIntialised.push_back(false);
    }
    
    //bool added;
    for (std::size_t i = 0; i < this->nodes.size(); ++i) {
        //added = false;
        for (std::size_t j = 0; j < 4; ++j) {
            if (treeNodes[j]->contains(allMortonCodes[this->nodes[i]])) {
                //added = true;
                treeNodes[j]->addNode(this->nodes[i]);
                
                childDistRatios[j].push_back(distRatios[i]);
                
                if (!ratioIntialised[j]) {
                    minBlockRatios[j] = distRatios[i];
                    maxBlockRatios[j] = distRatios[i];
                    ratioIntialised[j] = true;
                } else {
                    if (distRatios[i] < minBlockRatios[j]) {
                        minBlockRatios[j] = distRatios[i];
                    } 
                    if (distRatios[i] > maxBlockRatios[j]) {
                        maxBlockRatios[j] = distRatios[i];
                    }
                }
                
                if (sameColourStatus[j]) {
                    if (empty[j]) {
                        empty[j] = false;
                        colours[j] = colourMap[this->nodes[i]];
                    } else if (colours[j] != colourMap[this->nodes[i]] && colourMap[nodes[i]] != quadtree_constants::SPECIAL_COLOUR) {
                        // Note: SPECIAL_COLOUR represents the source. Since we only give the source
                        // a different colour to ensure that it's MortonBlock is included in the list
                        // (i.e. to cover case where it is an object). But we don't care what colour
                        // it has (never need to retrieve it) so if there is only one other colour in this
                        // quadtree leaf block then we consider the source having the same colour
                        if (colours[j] == quadtree_constants::SPECIAL_COLOUR) {
                            // Another the possibility is that the colours are different
                            // because the source was encountered first, in this case
                            // we update the colour of this block (and as above do not
                            // consider the source having a different colour
                            colours[j] = colourMap[nodes[i]];
                        } else {
                            sameColourStatus[j] = false;
                        }
                        //assert(colours[j] != quadtree_constants::EMPTY_COLOUR && "Valid colour has not been set");
                    }
                }
            }
        }
        //assert ((this->nodes[i] == source || added) && "Valid (i.e. non-source) node was not added to a child");
    }
}

bool QuadtreeNode::contains(const MortonCode& point) {
    return this->mortonCode.contains(point);
}

const MortonCode& QuadtreeNode::getMortonCode() const {
    return this->mortonCode;
}

double QuadtreeNode::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(QuadtreeNode*)*6;
    memoryUsage += sizeof(bool)*2;
    memoryUsage += sizeof(int);
    memoryUsage += sizeof(this->nodes);
    memoryUsage += this->mortonCode.computeIndexSizeBytes();
    memoryUsage += sizeof(NodeID)*this->nodes.size();
    return memoryUsage;
}

double QuadtreeNode::computeMemoryUsageBytes() {
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    memoryUsage += sizeof(NodeID)*this->nodes.capacity();
    return memoryUsage;
}

/*
 * MortonCode
 */

// For points (width is 1, so len represents whole maxRangeExp)
MortonCode::MortonCode(Coordinate x, Coordinate y, int maxRangeExp)
{
    // Need to check that max length can fit in this MortonCode
    //assert(maxRangeExp <= quadtree_constants::MAX_MORTON_LENGTH/2 && "Morton code will not fit in 64-bit zNumber");
    //assert(x >= 0 && y >= 0 && "Coordinates are not both positive");
    
    MortonNumber mx = x, my = y;

    mx = (mx | (mx << quadtree_constants::S[4])) & quadtree_constants::B[4];
    mx = (mx | (mx << quadtree_constants::S[3])) & quadtree_constants::B[3];
    mx = (mx | (mx << quadtree_constants::S[2])) & quadtree_constants::B[2];
    mx = (mx | (mx << quadtree_constants::S[1])) & quadtree_constants::B[1];
    mx = (mx | (mx << quadtree_constants::S[0])) & quadtree_constants::B[0];

    my = (my | (my << quadtree_constants::S[4])) & quadtree_constants::B[4];
    my = (my | (my << quadtree_constants::S[3])) & quadtree_constants::B[3];
    my = (my | (my << quadtree_constants::S[2])) & quadtree_constants::B[2];
    my = (my | (my << quadtree_constants::S[1])) & quadtree_constants::B[1];
    my = (my | (my << quadtree_constants::S[0])) & quadtree_constants::B[0];
    
    this->len = maxRangeExp << 1; // * 2
    this->zNumber = mx | (my << 1);
    this->zNumber <<= (quadtree_constants::MAX_MORTON_LENGTH - this->len);
    // Note: We can use this->len here instead of maxRangeExp << 1
    // because this is a point, the length will always be maxRangeExp
}

// For blocks (width is arbitrary and length will depend on width)
MortonCode::MortonCode(Coordinate x, Coordinate y, int width, int maxRangeExp)
{
    // Determine len of the morton code that represents this block
    int widthPower = 0;
    int widthCopy = width;
    while (widthCopy >>= 1) {
        ++widthPower;
    }

    MortonNumber mx = x, my = y;

    mx = (mx | (mx << quadtree_constants::S[4])) & quadtree_constants::B[4];
    mx = (mx | (mx << quadtree_constants::S[3])) & quadtree_constants::B[3];
    mx = (mx | (mx << quadtree_constants::S[2])) & quadtree_constants::B[2];
    mx = (mx | (mx << quadtree_constants::S[1])) & quadtree_constants::B[1];
    mx = (mx | (mx << quadtree_constants::S[0])) & quadtree_constants::B[0];

    my = (my | (my << quadtree_constants::S[4])) & quadtree_constants::B[4];
    my = (my | (my << quadtree_constants::S[3])) & quadtree_constants::B[3];
    my = (my | (my << quadtree_constants::S[2])) & quadtree_constants::B[2];
    my = (my | (my << quadtree_constants::S[1])) & quadtree_constants::B[1];
    my = (my | (my << quadtree_constants::S[0])) & quadtree_constants::B[0];
    
    this->len = (maxRangeExp-widthPower) << 1; // * 2
    this->zNumber = mx | (my << 1);
    this->zNumber <<= (quadtree_constants::MAX_MORTON_LENGTH - (maxRangeExp << 1));
}

MortonCode::MortonCode(MortonNumber zNumber, MortonLength length): zNumber(zNumber), len(length) {}

void MortonCode::setCode(const MortonCode& parentBlock, Direction dir)
{
    MortonNumber parentZNumber = parentBlock.getZNumber();
    EdgeID parentLen = parentBlock.getLength();
    //assert (parentLen+2 <= quadtree_constants::MAX_MORTON_LENGTH && "No room to split block");
    if (dir == SW) {
        this->zNumber = parentZNumber | quadtree_constants::ARRAY_SW_BIT_MASK[parentLen];
        // Note: This is equivalent to SW_MC << 64 - 2 - parentLen
        // I.e. we are the morton code for this quadrant to the end
        // of the parent block's morton code
    } else if (dir == SE) {
        this->zNumber = parentZNumber | quadtree_constants::ARRAY_SE_BIT_MASK[parentLen];
    } else if (dir == NW) {
        this->zNumber = parentZNumber | quadtree_constants::ARRAY_NW_BIT_MASK[parentLen];
    } else if (dir == NE) {
        this->zNumber = parentZNumber | quadtree_constants::ARRAY_NE_BIT_MASK[parentLen];
    } else {
        std::cerr << "Invalid direction provided" << std::endl;
        exit(1);
    }
    this->len = parentLen + 2;
}

MortonNumber MortonCode::getZNumber() const
{
    return this->zNumber;
}

EdgeID MortonCode::getLength() const
{
    return this->len;
}

bool MortonCode::contains(const MortonCode& point) const
{
    //assert (point.getLength() >= this->len && "Input code is shorter than this code");
    MortonNumber testNum = point.getZNumber();
    // Clear least significant bits that are not used by this code
    testNum = testNum >> (quadtree_constants::MAX_MORTON_LENGTH - this->len);
    testNum = testNum << (quadtree_constants::MAX_MORTON_LENGTH - this->len);
    // If these two numbers are the same we return true because
    // the test code is prefixed by this code
    return (this->zNumber == testNum);
}

/*
 * If point code is smaller than this morton code, return -1
 * If point code is contained in this morton code, return 0
 * If point code is larger than this morton code, return 1
 */
int MortonCode::compareMortonCode(const MortonCode& point) const {
    // We assume morton code is shorter than input
    //assert (point.len >= this->len && "Input code is shorter than this code");
    MortonNumber testNum = point.getZNumber();
    // Clear least significant bits that are not used by this code
    testNum = testNum >> (quadtree_constants::MAX_MORTON_LENGTH - this->len);
    testNum = testNum << (quadtree_constants::MAX_MORTON_LENGTH - this->len);
    
    if (testNum < this->zNumber) {
        return -1;
    } else if (testNum > this->zNumber) {
        return 1;
    } else /* if (testNum == this->zNumber */ {
        return 0;
    }
}

void MortonCode::decompose(int maxRangeExp, Coordinate& x, Coordinate& y, int& width)
{
    MortonNumber unshiftedZ = this->zNumber >> (quadtree_constants::MAX_MORTON_LENGTH - (maxRangeExp << 1)); 
    // Note: We can't use this->length here because we don't know if it is point or block
    
    MortonNumber mx = unshiftedZ & quadtree_constants::B[0], my = (unshiftedZ >> 1) & quadtree_constants::B[0];

    mx = (mx | (mx >> quadtree_constants::S[0])) & quadtree_constants::B[1];
    mx = (mx | (mx >> quadtree_constants::S[1])) & quadtree_constants::B[2];
    mx = (mx | (mx >> quadtree_constants::S[2])) & quadtree_constants::B[3];
    mx = (mx | (mx >> quadtree_constants::S[3])) & quadtree_constants::B[4];
    mx = (mx | (mx >> quadtree_constants::S[4])) & quadtree_constants::B[5];

    my = (my | (my >> quadtree_constants::S[0])) & quadtree_constants::B[1];
    my = (my | (my >> quadtree_constants::S[1])) & quadtree_constants::B[2];
    my = (my | (my >> quadtree_constants::S[2])) & quadtree_constants::B[3];
    my = (my | (my >> quadtree_constants::S[3])) & quadtree_constants::B[4];
    my = (my | (my >> quadtree_constants::S[4])) & quadtree_constants::B[5];
    
    x = mx;
    y = my;
    
    unsigned int widthPower = maxRangeExp - (this->len >> 1); // No precision will be lost because this->len is always multiple of two
    width = 1;
    width = width << widthPower;
}

bool MortonCode::operator==(const MortonCode& code) const  {
    if (this->len == code.getLength() && this->zNumber == code.getZNumber()) {
        return true;
    } else {
        return false;
    }
}

bool MortonCode::operator<(const MortonCode &code) const {
    if (this->zNumber < code.getZNumber()) {
        return true;
    } else {
        return false;
    }
}

bool MortonCode::operator>(const MortonCode &code) const  {
    if (this->zNumber > code.getZNumber()) {
        return true;
    } else {
        return false;
    }
}

Direction MortonCode::getQuadrant() {
    // We assume the morton code is for a block not a point
    MortonNumber smallestQuadrant = this->zNumber << (this->len - 2);
    smallestQuadrant = smallestQuadrant >> (quadtree_constants::MAX_MORTON_LENGTH - 2);
    if (smallestQuadrant == quadtree_constants::SW_MC) {
        return SW;
    } else if (smallestQuadrant == quadtree_constants::SE_MC) {
        return SE;
    } else if (smallestQuadrant == quadtree_constants::NW_MC) {
        return NW;
    } else if (smallestQuadrant == quadtree_constants::NE_MC) {
        return NE;
    } else {
        std::cerr << "Invalid quadrant from getQuadrant" << std::endl;
        exit(1);
    }
}

MortonCode MortonCode::getParentCode() {
    MortonNumber parentNumber = 0;
    if (this->len > 2) {
        // We must shift by number less than MAX_MORTON_LENGTH otherwise it's undefined in C++
        // (but in our case shifting by length is equivalent to setting to 0)
        parentNumber = this->zNumber >> (quadtree_constants::MAX_MORTON_LENGTH - this->len + 2);
        parentNumber = parentNumber << (quadtree_constants::MAX_MORTON_LENGTH - this->len + 2);
    }
    MortonCode parentCode(parentNumber,this->len-2);
    return parentCode;
}

Direction MortonCode::getNextQuadrant(Direction quadrant) {
    if (quadrant == SW) {
        return SE;
    } else if (quadrant == SE) {
        return NW;
    } else if (quadrant == NW) {
        return NE;
    } else if (quadrant == NE) {
        return SW;
    } else {
        return quadrant;
    }
}

Direction MortonCode::getPrevQuadrant(Direction quadrant) {
    if (quadrant == SW) {
        return NE;
    } else if (quadrant == SE) {
        return SW;
    } else if (quadrant == NW) {
        return SE;
    } else if (quadrant == NE) {
        return NW;
    } else {
        return quadrant;
    }
}

void MortonCode::printCode() {
     std::cout << "Morton Code = " << std::bitset<64>(this->zNumber) << " (Len=" << static_cast<int>(this->len) << ")" << std::endl;
}

// This should not be called on a NE block (check before using isLastMergeCode)
MortonCode MortonCode::getNextMergeCode()
{
    Direction nextQuadrant = this->getNextQuadrant();
    //assert (nextQuadrant != this->getQuadrant() && nextQuadrant > this->getQuadrant() && "Invalid input quadrant");
    MortonNumber parentNumber = 0;
    if (this->len > 2) {
        // We must shift by number less than MAX_MORTON_LENGTH otherwise it's undefined in C++
        // (but in our case shifting by length is equivalent to setting to 0)
        parentNumber = this->zNumber >> (quadtree_constants::MAX_MORTON_LENGTH - this->len + 2);
        parentNumber = parentNumber << (quadtree_constants::MAX_MORTON_LENGTH - this->len + 2);
    }    
    MortonNumber nextMergeNum = this->getChildZNumber(parentNumber,nextQuadrant);
    MortonCode nextCode(nextMergeNum,this->len);
    return nextCode;
}

// This should not be called on a SW block (check before using isFirstMergeCode)
MortonCode MortonCode::getPrevMergeCode()
{
    Direction prevQuadrant = this->getPrevQuadrant();
    //assert (prevQuadrant != this->getQuadrant() && prevQuadrant < this->getQuadrant() && "Invalid input quadrant");
    MortonNumber parentNumber = 0;
    if (this->len > 2) {
        // We must shift by number less than MAX_MORTON_LENGTH otherwise it's undefined in C++
        // (but in our case shifting by length is equivalent to setting to 0)
        parentNumber = this->zNumber >> (quadtree_constants::MAX_MORTON_LENGTH - this->len + 2);
        parentNumber = parentNumber << (quadtree_constants::MAX_MORTON_LENGTH - this->len + 2);
    }
    MortonNumber prevMergeNum = this->getChildZNumber(parentNumber,prevQuadrant);
    MortonCode prevCode(prevMergeNum,this->len);
    return prevCode;
}

bool MortonCode::isFirstMergeCode()
{
    if (this->len >= 2) {
        // If length is zero it cannot be merged with anything so we return false
        MortonNumber quadrantDigits = this->zNumber << (this->len - 2);
        quadrantDigits = quadrantDigits >> (quadtree_constants::MAX_MORTON_LENGTH - 2);
        if (quadrantDigits == quadtree_constants::SW_MC) {
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

bool MortonCode::isLastMergeCode()
{
    if (this->len >= 2) {
        // If length is zero it cannot be merged with anything so we return false
        MortonNumber quadrantDigits = this->zNumber << (this->len - 2);
        quadrantDigits = quadrantDigits >> (quadtree_constants::MAX_MORTON_LENGTH - 2);
        if (quadrantDigits == quadtree_constants::NE_MC) {
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

Direction MortonCode::getNextQuadrant()
{
    // This code quadrant is the last two occupied digits
    MortonNumber quadrantDigits = this->zNumber << (this->len - 2);
    quadrantDigits = quadrantDigits >> (quadtree_constants::MAX_MORTON_LENGTH - 2);
    if (quadrantDigits == quadtree_constants::SW_MC) {
        return SE;
    } else if (quadrantDigits == quadtree_constants::SE_MC) {
        return NW;
    } else if (quadrantDigits == quadtree_constants::NW_MC) {
        return NE;
    } else if (quadrantDigits == quadtree_constants::NE_MC) {
        std::cerr << "Current quadrant is NE - cannot getNextQuadrant" << std::endl;
        exit(1);
    } else {
        std::cerr << "Invalid quadrant from getQuadrant" << std::endl;
        exit(1);
    }
}

Direction MortonCode::getPrevQuadrant()
{
    // This code quadrant is the last two occupied digits
    MortonNumber quadrantDigits = this->zNumber << (this->len - 2);
    quadrantDigits = quadrantDigits >> (quadtree_constants::MAX_MORTON_LENGTH - 2);
    if (quadrantDigits == quadtree_constants::SW_MC) {
        std::cerr << "Current quadrant is SW - cannot getNextQuadrant" << std::endl;
        exit(1);
    } else if (quadrantDigits == quadtree_constants::SE_MC) {
        return SW;
    } else if (quadrantDigits == quadtree_constants::NW_MC) {
        return SE;
    } else if (quadrantDigits == quadtree_constants::NE_MC) {
        return NW;
    } else {
        std::cerr << "Invalid quadrant from getQuadrant" << std::endl;
        exit(1);
    }
}

MortonNumber MortonCode::getChildZNumber(MortonNumber parentZNumber, Direction childQuadrant)
{
    MortonNumber childZNumber = 0;
    if (childQuadrant == SW) {
       childZNumber = parentZNumber | quadtree_constants::ARRAY_SW_BIT_MASK[this->len-2];
        // Note: This is equivalent to SW_MC << 64 - 2 - parentLen
        // I.e. we are the morton code for this quadrant to the end
        // of the parent block's morton code
    } else if (childQuadrant == SE) {
        childZNumber = parentZNumber | quadtree_constants::ARRAY_SE_BIT_MASK[this->len-2];
    } else if (childQuadrant == NW) {
        childZNumber = parentZNumber | quadtree_constants::ARRAY_NW_BIT_MASK[this->len-2];
    } else if (childQuadrant == NE) {
        childZNumber = parentZNumber | quadtree_constants::ARRAY_NE_BIT_MASK[this->len-2];
    } else {
        std::cerr << "Invalid direction provided" << std::endl;
        exit(1);
    }
    return childZNumber;
}

double MortonCode::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(this->zNumber);
    memoryUsage += sizeof(this->len);
    return memoryUsage;
}

double MortonCode::computeMemoryUsageBytes() {
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    return memoryUsage;
}
