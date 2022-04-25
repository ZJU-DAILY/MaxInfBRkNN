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

#include "SimpleQuadtree.h"

#include "../utility/geometry.h"
#include "../utility/utility.h"

#include <cmath>

SimpleQuadtree::SimpleQuadtree(int maxLeafItems): maxLeafItems(maxLeafItems), maxRange(0), xTranslation(0), yTranslation(0), 
    isPointSet(false), objSetType(""), objSetDensity(0), objSetSize(0) {}

SimpleQuadtree::~SimpleQuadtree()
{
    // Recursively delete all quadtree nodes
    delete rootNode;
}

void SimpleQuadtree::buildFromGraph(Graph_RoadNetwork &graph) {
    
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
    
    // We adjust the coordinate of the node with the minimum X (resp. Y) 
    // coordinate such that it has a coordinate of 0 and all other X (resp. Y)
    // coordinates are adjusted all adjust relative to this 
    
    NodeID numNodes = graph.getNumNodes();
    Coordinate x, y, minX, maxX, minY, maxY;
    graph.getMinMaxCoordinates(minX,maxX,minY,maxY);
    this->xTranslation = minX;
    this->yTranslation = minY;
    
    this->rootNode = new SimpleQuadtreeNode(true,true);
    this->rootNode->setBlockCoordinates(this->xTranslation,this->yTranslation,this->maxRange+1);
    
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
//     std::cout << "Nodes in SimpleQuadtree Leaf Nodes = " << vertexesInLeafNodes << std::endl;    
//     std::cout << "Non-Empty Leaf Nodes = " << nonEmptyLeafNodes << std::endl;    
}

void SimpleQuadtree::buildFromPointSet(Graph_RoadNetwork &graph, std::vector<NodeID>& pointSet, std::string setType, double setDensity, int setVariable, int setSize) {
    this->objSetType = setType;
    this->objSetDensity = setDensity;
    this->objSetVariable = setVariable;
    this->objSetSize = setSize;
    this->isPointSet = true;
    
    // Compute maximum coordinate range to see how many bits are required to
    // capture morton and how much to translate X and Y coordinates
    Coordinate x, y, minX, maxX, minY, maxY;
    if (pointSet.size() > 0) {
        graph.getCoordinates(pointSet[0],x,y);
        minX = x;
        maxX = x;
        minY = y;
        maxY = y;
        for (std::size_t i = 1; i < pointSet.size(); ++i) {
            graph.getCoordinates(pointSet[i],x,y);
            if (x < minX) {
                minX = x;
            }
            if (x > maxX) {
                maxX = x;
            }
            if (y < minY) {
                minY = y;
            }
            if (y > maxY) {
                maxY = y;
            }
        }
    } else {
        std::cerr << "Empty object set provided!" << std::endl;
        exit(1);
    }
    
    Coordinate xRange = maxX-minX, yRange = maxY-minY;
    Coordinate maxCoordinateRange = (xRange > yRange) ? xRange : yRange;
    this->xTranslation = minX;
    this->yTranslation = minY;
    
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

    // Debug Output
//     std::cout << "minX = " << minX << std::endl;
//     std::cout << "maxX = " << maxX << std::endl;
//     std::cout << "minY = " << minY << std::endl;
//     std::cout << "maxY = " << maxY << std::endl;
//     std::cout << "MinX Node= " << minXNode << std::endl;
//     std::cout << "MaxX Node= " << maxXNode << std::endl;
//     std::cout << "MinY Node= " << minYNode << std::endl;
//     std::cout << "MaxY Node= " << maxYNode << std::endl;
//     std::cout << "xRange = " << xRange << std::endl;
//     std::cout << "yRange = " << yRange << std::endl;
    
    this->rootNode = new SimpleQuadtreeNode(true,true);
    this->rootNode->setBlockCoordinates(this->xTranslation,this->yTranslation,this->maxRange+1);

    // Reserve space for all nodes
    this->mortonCodes.reserve(this->objSetSize);
    this->relativeXCoordinates.reserve(this->objSetSize);
    this->relativeYCoordinates.reserve(this->objSetSize);

    for (std::size_t i = 0; i < pointSet.size(); ++i) {
        // Since we do not have all nodes in the graph we need a way to map nodes to
        // their index in the respective data
        this->nodeToIdxMap[pointSet[i]] = i;

        // We translate all X and Y coordinates such that (0,0) is the most
        // SW point and coordinates increase in the N and E directions
        // i.e. all X and Y coordinates are positive
        graph.getTranslatedCoordinates(pointSet[i],x,y,this->xTranslation,this->yTranslation);
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
        this->insert(pointSet[i]);
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
//     std::cout << "Nodes in SimpleQuadtree Leaf Nodes = " << vertexesInLeafNodes << std::endl;    
//     std::cout << "Non-Empty Leaf Nodes = " << nonEmptyLeafNodes << std::endl;
}

void SimpleQuadtree::insert(NodeID node)
{
    SimpleQuadtreeNode* currentNode = this->rootNode;
    
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

void SimpleQuadtree::split(SimpleQuadtreeNode* treeNode)
{
    // Mark node is being no longer a leaf node
    treeNode->setLeaf(false);
    
    // Create and update tree node pointers for split
    treeNode->sw = new SimpleQuadtreeNode(treeNode->mortonCode,SW);
    treeNode->sw->computeAndSetBlockCoordinates(this->xTranslation,this->yTranslation,this->maxRangeExp);
    treeNode->se = new SimpleQuadtreeNode(treeNode->mortonCode,SE);
    treeNode->se->computeAndSetBlockCoordinates(this->xTranslation,this->yTranslation,this->maxRangeExp);
    treeNode->nw = new SimpleQuadtreeNode(treeNode->mortonCode,NW);
    treeNode->nw->computeAndSetBlockCoordinates(this->xTranslation,this->yTranslation,this->maxRangeExp);
    treeNode->ne = new SimpleQuadtreeNode(treeNode->mortonCode,NE);
    treeNode->ne->computeAndSetBlockCoordinates(this->xTranslation,this->yTranslation,this->maxRangeExp);
    
    if (!this->isPointSet) {
        // If this is not a point set then we don't need to map 
        // between the internal index and the node ID
        treeNode->pushDownNodes(this->mortonCodes);
    } else {
        treeNode->pushDownNodes(this->mortonCodes,this->nodeToIdxMap);
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

int SimpleQuadtree::getMaxLeafItems()
{
    return this->maxLeafItems;
}

int SimpleQuadtree::getMaxRange() {
    return this->maxRange;
}

int SimpleQuadtree::getMaxRangeExp() {
    return this->maxRangeExp;
}

int SimpleQuadtree::getXTranslation() {
    return this->xTranslation;
}

int SimpleQuadtree::getYTranslation() {
    return this->yTranslation;
}

std::string SimpleQuadtree::getObjSetType()
{
    return this->objSetType;
}

double SimpleQuadtree::getObjSetDensity()
{
    return this->objSetDensity;
}

int SimpleQuadtree::getObjSetSize()
{
    return this->objSetSize;
}

int SimpleQuadtree::getObjSetVariable()
{
    return this->objSetVariable;
}

void SimpleQuadtree::populateCoordinates(NodeID node, Coordinate &x, Coordinate &y) {
    // Note: We cannot use this function if object has been deserialized
    assert(this->relativeXCoordinates.size() > 0 && "Relative X and Y coordinates do not exist (they are not serialized)");
    x = this->relativeXCoordinates[node];
    y = this->relativeYCoordinates[node];
}

EuclideanDistance SimpleQuadtree::getEuclideanDistance(NodeID s, NodeID t) {
    // Note: We cannot use this function if object has been deserialized
    assert(this->relativeXCoordinates.size() > 0 && "Relative X and Y coordinates do not exist (they are not serialized)");
    return geometry::getEuclideanDist(this->relativeXCoordinates[s],this->relativeYCoordinates[s],
                                      this->relativeXCoordinates[t],this->relativeYCoordinates[t]);
}

const std::vector<MortonCode> &SimpleQuadtree::getMortonCodes() const {
    return this->mortonCodes;
}

double SimpleQuadtree::computeIndexSize()
{
    double memoryUsage = 0;
    std::deque<SimpleQuadtreeNode*> stack;
    stack.push_back(this->rootNode);
    while (stack.size() > 0) {
        SimpleQuadtreeNode* currNode = stack.front();
        stack.pop_front();
        memoryUsage += currNode->computeIndexSizeBytes();
        if (!currNode->isLeafNode()) {
            stack.push_back(currNode->sw);
            stack.push_back(currNode->se);
            stack.push_back(currNode->nw);
            stack.push_back(currNode->ne);
        }
    }
    return memoryUsage/(1024*1024);
}

double SimpleQuadtree::computeMemoryUsage() {
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    std::deque<SimpleQuadtreeNode*> stack;
    stack.push_back(this->rootNode);
    while (stack.size() > 0) {
        SimpleQuadtreeNode* currNode = stack.front();
        stack.pop_front();
        memoryUsage += currNode->computeMemoryUsageBytes();
        if (!currNode->isLeafNode()) {
            stack.push_back(currNode->sw);
            stack.push_back(currNode->se);
            stack.push_back(currNode->nw);
            stack.push_back(currNode->ne);
        }
    }
    memoryUsage += this->objSetType.size();
    memoryUsage += sizeof(Coordinate)*relativeXCoordinates.capacity();
    memoryUsage += sizeof(Coordinate)*relativeYCoordinates.capacity();
    for (std::size_t i = 0; i < this->mortonCodes.size(); ++i) {
        memoryUsage += this->mortonCodes[i].computeMemoryUsageBytes();
    }
    memoryUsage += sizeof(MortonCode)*(this->mortonCodes.capacity()-this->mortonCodes.size());
    memoryUsage += utility::estimateUnorderedMapMemoryUsageBytes(this->nodeToIdxMap.size(),sizeof(std::pair<NodeID,int>),this->nodeToIdxMap.bucket_count());
    return memoryUsage/(1024*1024);
}

int SimpleQuadtree::getTotalRegions()
{
    int totalRegions = 0;
    if (this->rootNode != NULL) {
        totalRegions = this->rootNode->getTotalRegions();
    }
    return totalRegions;
}


/*
 * SimpleQuadtreeNode
 */

SimpleQuadtreeNode::SimpleQuadtreeNode(MortonCode& parentCode, Direction dir)
{
    this->mortonCode.setCode(parentCode,dir);
    this->isRoot = false; // Cannot be root because it has a parent
    this->isLeaf = true;
    this->sw = this->se = this->nw = this->ne = NULL;
    this->numNodes = 0;
}

SimpleQuadtreeNode::SimpleQuadtreeNode(bool isRoot, bool isLeaf)
{
    this->isRoot = isRoot;
    this->isLeaf = isLeaf;
    this->sw = this->se = this->nw = this->ne = NULL;
    this->numNodes = 0;
}

SimpleQuadtreeNode::~SimpleQuadtreeNode() {
    // Recursively delete all quadtree nodes
    delete this->sw;
    delete this->se;
    delete this->nw;
    delete this->ne;
}

bool SimpleQuadtreeNode::isRootNode() {
    return this->isRoot;
}

int SimpleQuadtreeNode::getNumNodes()
{
    return this->numNodes;
}

bool SimpleQuadtreeNode::isLeafNode() {
    return this->isLeaf;
}

void SimpleQuadtreeNode::setLeaf(bool isLeaf) {
    this->isLeaf = isLeaf;
}

const std::vector<NodeID>& SimpleQuadtreeNode::getNodes() const {
    return this->nodes;
}

NodeID SimpleQuadtreeNode::getFirstNode()
{
    return this->nodes[0];
}

void SimpleQuadtreeNode::addNode(NodeID node) {
    this->numNodes++;
    this->nodes.push_back(node);
}

void SimpleQuadtreeNode::pushDownNodes(std::vector<MortonCode>& allMortonCodes) {
    for (std::size_t i = 0; i < this->nodes.size(); ++i) {
        this->pushDownNode(this->nodes[i],allMortonCodes[this->nodes[i]]);
    }
    this->nodes.clear();
}

void SimpleQuadtreeNode::pushDownNodes(std::vector<MortonCode>& allMortonCodes, std::unordered_map<NodeID,int>& nodeToIdxMap) {
    int index;
    for (std::size_t i = 0; i < this->nodes.size(); ++i) {
        index = nodeToIdxMap[this->nodes[i]];
        this->pushDownNode(this->nodes[i],allMortonCodes[index]);
    }
    this->nodes.clear();
}

void SimpleQuadtreeNode::pushDownNode(NodeID node, MortonCode& mortonCode) {
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

bool SimpleQuadtreeNode::contains(const MortonCode& point) {
    return this->mortonCode.contains(point);
}

const MortonCode& SimpleQuadtreeNode::getMortonCode() const {
    return this->mortonCode;
}

void SimpleQuadtreeNode::getBlockCoordinates(Coordinate& left, Coordinate& bottom, int& width)
{
    left = this->graphX;
    bottom = this->graphY;
    width = this->width;
}

void SimpleQuadtreeNode::computeAndSetBlockCoordinates(int xTranslation, int yTranslation, int maxRangeExp)
{
    this->mortonCode.decompose(maxRangeExp,this->graphX,this->graphY,this->width);
    this->graphX += xTranslation;
    this->graphY += yTranslation;
}

void SimpleQuadtreeNode::setBlockCoordinates(Coordinate x, Coordinate y, int width)
{
    this->graphX = x;
    this->graphY = y;
    this->width = width;
}

double SimpleQuadtreeNode::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(SimpleQuadtreeNode*)*4;
    memoryUsage += sizeof(bool)*2;
    memoryUsage += sizeof(Coordinate)*2;
    memoryUsage += sizeof(int)*2;
    memoryUsage += sizeof(this->nodes);
    memoryUsage += this->mortonCode.computeIndexSizeBytes();
    memoryUsage += sizeof(NodeID)*this->nodes.size();
    return memoryUsage;
}

double SimpleQuadtreeNode::computeMemoryUsageBytes() {
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    memoryUsage += sizeof(NodeID)*this->nodes.capacity();
    return memoryUsage;
}

int SimpleQuadtreeNode::getTotalRegions()
{
    int totalRegions = 1;
    if (this->sw != NULL) {
        totalRegions += sw->getTotalRegions();
    }
    if (this->se != NULL) {
        totalRegions += se->getTotalRegions();
    }
    if (this->nw != NULL) {
        totalRegions += nw->getTotalRegions();
    }
    if (this->ne != NULL) {
        totalRegions += ne->getTotalRegions();
    }
    return totalRegions;
}
