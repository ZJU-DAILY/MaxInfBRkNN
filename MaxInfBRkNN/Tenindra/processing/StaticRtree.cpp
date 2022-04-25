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

#include "StaticRtree.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <unordered_map>

StaticRtree::StaticRtree(int branchFactor): 
    branchFactor(branchFactor), objSetType(""), objSetDensity(0), objSetSize(0)
{
    root = NULL;
}

StaticRtree::StaticRtree(int branchFactor, std::string setType, double setDensity, int setVariable, int setSize): 
    branchFactor(branchFactor), objSetType(setType), objSetDensity(setDensity), objSetVariable(setVariable), objSetSize(setSize)
{
    root = NULL;
}

StaticRtree::~StaticRtree()
{
    // Collect all pointers to nodes and delete them
    std::vector<StaticRtreeNode*> allNodePtrs = this->getAllNodePtrs();
    
    for (std::size_t i = 0; i < allNodePtrs.size(); ++i) {
        // Note: This will try to delete child pointers in leaf nodes
        // (which are null) but this is fine as they are simply no-ops
        delete allNodePtrs[i];
    }
}

void StaticRtree::bulkLoad(std::vector<NodeID>& objectIDs, std::vector<CoordinatePair>& objectCoords)
{
    // The algorithm used to bulkload object into the R-tree is based on:
    // Leutenegger, Scott T., Mario Lopez, and Jeffrey Edgington
    // STR: A simple and efficient algorithm for R-tree packing
    // Proceedings of the 13th IEEE International Conference on Data
    // Engineering (ICDE), pp. 497-506, 1997
    
    // Create a (leaf) entry for each input object
    std::size_t numEntries = objectIDs.size();
    std::vector<StaticRtreeEntry> leafLevelEntries(numEntries);
    for (std::size_t i = 0; i < numEntries; ++i) {
        leafLevelEntries[i].setEntryData(objectCoords[i],NULL,objectIDs[i]);
    }
    
    // Recursively build all levels (starting at leaf)
    this->root = this->buildLevel(0, leafLevelEntries);
    
}

StaticRtreeNode* StaticRtree::buildLevel(int level, std::vector<StaticRtreeEntry>& levelEntries)
{
    int numNodes = std::ceil(static_cast<double>(levelEntries.size())/this->branchFactor);
    int numSlices = std::ceil(std::sqrt(numNodes));
//     std::cout << "\n\nProcessing level " << level << std::endl;
//     std::cout << "numNodes = " << numNodes << std::endl;
//     std::cout << "numSlices = " << numSlices << std::endl;
//     std::cout << "branchFactor = " << this->branchFactor << std::endl;
//     std::cout << "levelEntries.size() = " << levelEntries.size() << std::endl;
    
    // Order entries by their X coordinate and divide into numSlices slices
    // Note: Each slice contain numSlices*this->branchFactor entries
    std::sort(levelEntries.begin(),levelEntries.end(),StaticRtreeEntry::compareX);
    
    // Order each slice by the Y coordinate and create node out of every 
    // this->branchFactor group of nodes (except the last may have fewer)
    std::vector<StaticRtreeNode*> nodes;
    std::vector<Rectangle> rectangles;
    std::vector<StaticRtreeEntry>::iterator sliceStartIt, sliceEndIt, nodeStartIt, nodeEndIt;
    std::size_t nextSliceStartIdx = 0, nextGroupStartIdx = 0;
    for (std::size_t i = 0; i < levelEntries.size(); ) {
        // Sort the next group of this->branchFactor entries by their Y value
//         std::cout << "\nProcessing New Slice" << std::endl;
//         std::cout << "i = " << i << std::endl;
        // This is essentially considering one of the numSlices partitions
        // and we sort it according to the y value
        sliceStartIt = levelEntries.begin();
        std::advance(sliceStartIt,i);
        int maxEntriesInSlice = numSlices*this->branchFactor;
//         std::cout << "maxEntriesInSlice = " << maxEntriesInSlice << std::endl;
        if (i + maxEntriesInSlice < levelEntries.size()) {
            sliceEndIt = levelEntries.begin();
            std::advance(sliceEndIt,i + maxEntriesInSlice); // Note: We advanced by i first (for next slice)
            nextSliceStartIdx = i + maxEntriesInSlice;
        } else {
            // This mean we have reached the last slices and there
            //  were not be enough to fill the slice completely
            sliceEndIt = levelEntries.end();
            nextSliceStartIdx = levelEntries.size();
        }
//         std::cout << "Actual Entries in slices = " << nextSliceStartIdx-i << std::endl;
        std::sort(sliceStartIt,sliceEndIt,StaticRtreeEntry::compareY);
        
        // Create node for each group of entries in the current slice
        for (std::size_t j = i; j < nextSliceStartIdx;) {
            StaticRtreeNode* rtreeNode = new StaticRtreeNode(level);
            Rectangle directoryRect;
            int maxEntriesInGroup = this->branchFactor;
            nodeStartIt = levelEntries.begin();
            std::advance(nodeStartIt,j);
            if (j + maxEntriesInGroup < nextSliceStartIdx) {
                nodeEndIt = levelEntries.begin();
                std::advance(nodeEndIt,j + maxEntriesInGroup); // Note: We advanced by j first (for next group in slice)
                nextGroupStartIdx = j + maxEntriesInGroup;
            } else {
                // This means we reached the last group of entries in this slice and 
                // there were not enough of them for to fill the node completely
                nodeEndIt = sliceEndIt;
                nextGroupStartIdx = nextSliceStartIdx;
            }
            assert (nodeStartIt != nodeEndIt);
            if (nodeStartIt != nodeEndIt) {
//                 std::cout << "Add First Entry To A Node" << std::endl;
//                 nodeStartIt->mbr.printRectangle();
                // Initialise directory rectangle to cover first entry
                directoryRect.setDimensions(nodeStartIt->mbr);
                rtreeNode->addEntry(std::move(*nodeStartIt));
                ++nodeStartIt;
                for (; nodeStartIt != nodeEndIt; ++nodeStartIt) {
//                     std::cout << "Add New Entry To A Node" << std::endl;
//                     nodeStartIt->mbr.printRectangle();
                    directoryRect.expand(nodeStartIt->mbr);
                    rtreeNode->addEntry(std::move(*nodeStartIt));
                }
            } else {
                std::cout << "This shouldn't happen!" << std::endl;
                std::exit(1);
            }
            nodes.push_back(rtreeNode);
            rectangles.push_back(directoryRect);
//             std::cout << "nextGroupStartIdx = " << nextGroupStartIdx << std::endl;
            j = nextGroupStartIdx;
        }
        i = nextSliceStartIdx;
    }
    
//     std::cout << "\nNum nodes created = " << nodes.size() << std::endl;
    
    StaticRtreeNode* rootNode = NULL;
    std::size_t numEntries = nodes.size();
    if (numEntries == 1) {
        // This mean all entries were added to the root
        // and we can stop building the tree
        rootNode = nodes[0];
    } else {
        // Create the entries for each node and recursively build next level
        std::vector<StaticRtreeEntry> nextLevelEntries(numEntries);
        for (std::size_t i = 0; i < numEntries; ++i) {
            nextLevelEntries[i].setEntryData(std::move(rectangles[i]),nodes[i]);
        }
        rootNode = this->buildLevel(level+1,nextLevelEntries);
    }
    
    return rootNode;
}

BinaryMinHeap<EuclideanDistanceSqrd,RtreeDataTuple> StaticRtree::getKNNs(unsigned int k, Coordinate queryPointX, Coordinate queryPointY, 
                                                                         std::vector<NodeID>& kNNs, std::vector<double>& kNNDistancesSqrd)
{
    BinaryMinHeap<EuclideanDistanceSqrd,RtreeDataTuple> heap;
    EuclideanDistanceSqrd minDist;
    
    heap.insert(RtreeDataTuple(this->root),0);
    while (heap.size() > 0) {
        minDist = heap.getMinKey();
        RtreeDataTuple minElement = heap.extractMinElement();
        if (minElement.isNode) {
            if(minElement.nodePtr->isLeaf()) {
                // Entries in leaf node are objects, compute Euclidean dist to it
                for (std::size_t i = 0; i < minElement.nodePtr->entries.size(); ++i) {
                    minDist = minElement.nodePtr->entries[i].mbr.getObjectDistSqrd(queryPointX,queryPointY);
                    heap.insert(RtreeDataTuple(minElement.nodePtr->entries[i].id),minDist);
                }
            } else {
                // Entries in non-leaf nodes are other tree nodes, compute mindist
                for (std::size_t i = 0; i < minElement.nodePtr->entries.size(); ++i) {
                    minDist = minElement.nodePtr->entries[i].mbr.getMinDistSqrd(queryPointX,queryPointY);
                    heap.insert(RtreeDataTuple(minElement.nodePtr->entries[i].rtreeNode),minDist);
                }
            }
        } else {
            // Then we have found object and the minimum distance to it is the actual distance
            // and since there is no closer node or object, it is the next nearest neighbour
            kNNs.push_back(minElement.objectID);
            kNNDistancesSqrd.push_back(minDist);
            if (kNNs.size() == k) {
                break;
            }
        }
    }
    
    return heap;
}

bool StaticRtree::getNextNearestNeighbour(BinaryMinHeap<EuclideanDistanceSqrd,RtreeDataTuple>& heap, Coordinate queryPointX, 
                                          Coordinate queryPointY, NodeID& nn, EuclideanDistanceSqrd& nnDistSqrd)
{
    EuclideanDistanceSqrd minDist;
    while (heap.size() > 0) {
        minDist = heap.getMinKey();
        RtreeDataTuple minElement = heap.extractMinElement();
        if (minElement.isNode) {
            if(minElement.nodePtr->isLeaf()) {
                // Entries in leaf node are objects, compute Euclidean dist to it
                for (std::size_t i = 0; i < minElement.nodePtr->entries.size(); ++i) {
                    minDist = minElement.nodePtr->entries[i].mbr.getObjectDistSqrd(queryPointX,queryPointY);
                    heap.insert(RtreeDataTuple(minElement.nodePtr->entries[i].id),minDist);
                }
            } else {
                // Entries in non-leaf nodes are other tree nodes, compute mindist
                for (std::size_t i = 0; i < minElement.nodePtr->entries.size(); ++i) {
                    minDist = minElement.nodePtr->entries[i].mbr.getMinDistSqrd(queryPointX,queryPointY);
                    heap.insert(RtreeDataTuple(minElement.nodePtr->entries[i].rtreeNode),minDist);
                }
            }
        } else {
            // Then we have found object and the minimum distance to it is the actual distance
            // and since there is no closer node or object, it is the next nearest neighbour
            nn = minElement.objectID;
            nnDistSqrd = minDist;
            return true;
        }
    }
    // If we reach it means there are no objects to be found left
    return false;
}

std::vector<NodeID> StaticRtree::windowQuery(Rectangle& rect)
{
    std::vector<NodeID> windowsPoints;
    std::vector<StaticRtreeNode*> queue;
    
    queue.push_back(this->root);
    while (queue.size() > 0) {
        StaticRtreeNode* currentNode = queue.back();
        queue.pop_back();
        if(currentNode->isLeaf()) {
            for (std::size_t i = 0; i < currentNode->entries.size(); ++i) {
                if (currentNode->entries[i].mbr.intersects(rect)) {
                    windowsPoints.push_back(currentNode->entries[i].id);
                }
            }
        } else {
            for (std::size_t i = 0; i < currentNode->entries.size(); ++i) {
                if (currentNode->entries[i].mbr.intersects(rect)) {
                    queue.push_back(currentNode->entries[i].rtreeNode);
                }
            }
        }
    }
    return windowsPoints;
}

std::vector<StaticRtreeNode*> StaticRtree::getAllNodePtrs()
{
    std::vector<StaticRtreeNode*> allNodePtrs;
    std::vector<StaticRtreeNode*> currentLevel, nextLevel;
    
    nextLevel.push_back(this->root);
    while(nextLevel.size() != 0){
        allNodePtrs.insert(allNodePtrs.end(), nextLevel.begin(), nextLevel.end());
        currentLevel.swap(nextLevel);
        nextLevel.clear();
        for (std::size_t i = 0; i < currentLevel.size(); ++i) {
            if (currentLevel[i] != NULL) {
                for (std::size_t j = 0; j < currentLevel[i]->entries.size(); ++j) {
                    nextLevel.push_back(currentLevel[i]->entries[j].rtreeNode);
                }
            }
        }
    }
    return allNodePtrs;
}

void StaticRtree::printTree()
{
    std::unordered_map<int,int> levelTotals;
    std::vector<std::vector<StaticRtreeNode*>> treeNodePointersByLevel;
    std::vector<StaticRtreeNode*> currentLevel, nextLevel;
    
    nextLevel.push_back(this->root);
    while(nextLevel.size() != 0){
        treeNodePointersByLevel.push_back(nextLevel);
        currentLevel.swap(nextLevel);
        nextLevel.clear();
        for (std::size_t i = 0; i < currentLevel.size(); ++i) {
            if (currentLevel[i] != NULL) {
                if (levelTotals.find(currentLevel[i]->level) != levelTotals.end()) {
                    levelTotals[currentLevel[i]->level]++;
                } else {
                    levelTotals[currentLevel[i]->level] = 1;
                }
                for (std::size_t j = 0; j < currentLevel[i]->entries.size(); ++j) {
                    nextLevel.push_back(currentLevel[i]->entries[j].rtreeNode);
                }
            } else {
                // Then this is an actual data object
                if (levelTotals.find(-1) != levelTotals.end()) {
                    levelTotals[-1]++;
                } else {
                    levelTotals[-1] = 1;
                }
            }
        }
    }
    
    int numLevels = levelTotals.size();
    for (int i = -1; i < numLevels-1; ++i) {
        if (i == -1) {
            std::cout << "Data Objects: " << levelTotals[i] << std::endl;
        } else {
            std::cout << "Level " << i << ": " << levelTotals[i] << std::endl;
        }
    }
    
//     for (int i = treeNodePointersByLevel.size()-2; i >= 0; --i) {
//         /*int levelSize = 0;
//         for (std::size_t j = 0; j < treeNodeLevel[i].size(); ++j) {
//             levelSize += this->treeNodes[i].getNumVertices();
//         }*/
//         std::cout << "\nLevel " << static_cast<int>(treeNodePointersByLevel.size())-2-i << " Nodes: " << std::endl;
//         for (std::size_t j = 0; j < treeNodePointersByLevel[i].size(); ++j) {
//             std::cout << "Printing Node " << j << std::endl;
//             treeNodePointersByLevel[i][j]->printNode();
//             std::cout << std::endl;
//         }
//         std::cout << std::endl;
//     }
}

std::string StaticRtree::getObjSetType()
{
    return this->objSetType;
}

double StaticRtree::getObjSetDensity()
{
    return this->objSetDensity;
}

int StaticRtree::getObjSetSize()
{
    return this->objSetSize;
}

int StaticRtree::getObjSetVariable()
{
    return this->objSetVariable;
}

int StaticRtree::getBranchFactor()
{
    return this->branchFactor;
}

double StaticRtree::computeIndexSize()
{
    double memoryUsageBytes = 0;
    std::vector<StaticRtreeNode*> allNodePtrs = this->getAllNodePtrs();
    for (std::size_t i = 0; i < allNodePtrs.size(); ++i) {
        if (allNodePtrs[i] != NULL) {
            memoryUsageBytes += allNodePtrs[i]->computeIndexSizeBytes();;
        }
    }
    return memoryUsageBytes/(1024*1024);
}

double StaticRtree::computeMemoryUsage()
{
    double memoryUsageBytes = 0;
    memoryUsageBytes += sizeof(*this);
    std::vector<StaticRtreeNode*> allNodePtrs = this->getAllNodePtrs();
    for (std::size_t i = 0; i < allNodePtrs.size(); ++i) {
        if (allNodePtrs[i] != NULL) {
            memoryUsageBytes += allNodePtrs[i]->computeMemoryUsageBytes();;
        }
    }
    return memoryUsageBytes/(1024*1024);
}

/*
 * StaticRtreeNode
 */

StaticRtreeNode::StaticRtreeNode(int level): level(level)
{

}

void StaticRtreeNode::addEntry(StaticRtreeEntry entry)
{
    this->entries.push_back(entry);

}

bool StaticRtreeNode::isLeaf()
{
    return this->level == 0;
}

void StaticRtreeNode::printNode()
{
    for (std::size_t i = 0; i < this->entries.size(); ++i) {
        this->entries[i].printEntry();
    }
}

double StaticRtreeNode::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    for (std::size_t i = 0; i < this->entries.size(); ++i) {
        memoryUsage += this->entries[i].computeIndexSizeBytes();
    }
    memoryUsage += sizeof(int);
    return memoryUsage;

}

double StaticRtreeNode::computeMemoryUsageBytes()
{
    double memoryUsage = 0;
    for (std::size_t i = 0; i < this->entries.size(); ++i) {
        memoryUsage += this->entries[i].computeMemoryUsageBytes();
    }
    memoryUsage += sizeof(*this);
    return memoryUsage;
}

/*
 * StaticRtreeEntry
 */

void StaticRtreeEntry::setEntryData(CoordinatePair& ptCoords, StaticRtreeNode* rtreeNode, NodeID id)
{
    this->mbr.setDimensions(ptCoords.first,ptCoords.second,ptCoords.first,ptCoords.second);
    this->mbr.setCentroid();
    this->rtreeNode = rtreeNode;
    this->id = id;
}

void StaticRtreeEntry::setEntryData(Rectangle rect, StaticRtreeNode* rtreeNode, NodeID id)
{
    this->mbr = rect;
    this->mbr.setCentroid(); 
    // Note: We assume centroid  has not been set (which allows us to 
    // delay computing it until all rectangle expansions are complete)
    this->rtreeNode = rtreeNode;
    this->id = id;
}

bool StaticRtreeEntry::compareX(const StaticRtreeEntry& entry1, const StaticRtreeEntry& entry2)
{
    return entry1.mbr.centreX < entry2.mbr.centreX;
}

bool StaticRtreeEntry::compareY(const StaticRtreeEntry& entry1, const StaticRtreeEntry& entry2)
{
    return entry1.mbr.centreY < entry2.mbr.centreY;
}

void StaticRtreeEntry::printEntry()
{
    std::cout << "Rectangle: ";
    this->mbr.printRectangle();
    std::cout << "ID = " << this->id << std::endl;
}

double StaticRtreeEntry::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(StaticRtreeNode*);
    memoryUsage += sizeof(NodeID);
    memoryUsage += this->mbr.computeIndexSizeBytes();
    return memoryUsage;

}

double StaticRtreeEntry::computeMemoryUsageBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    return memoryUsage;
}

/* 
 * Rectangle
 */

void Rectangle::setDimensions(Coordinate left, Coordinate bottom, Coordinate right, Coordinate top)
{
    this->left = left;
    this->bottom = bottom;
    this->right = right;
    this->top = top;
}

void Rectangle::setDimensions(Rectangle& rect)
{
    this->left = rect.left;
    this->bottom = rect.bottom;
    this->right = rect.right;
    this->top = rect.top;
}

void Rectangle::setCentroid()
{
    this->centreX = (this->left+this->right)/2;
    this->centreY = (this->bottom+this->top)/2;
}

void Rectangle::expand(Rectangle& rect)
{
    // Expand the current rectangle to all cover given rectangle
    this->left = this->left > rect.left ? rect.left : this->left;
    this->bottom = this->bottom > rect.bottom ? rect.bottom : this->bottom;
    this->right = this->right < rect.right ? rect.right : this->right;
    this->top = this->top < rect.top ? rect.top : this->top;
}

EuclideanDistanceSqrd Rectangle::getMinDistSqrd(Coordinate x, Coordinate y)
{
    Coordinate nearestX;
    if (x < this->left) {
        nearestX = this->left;
    } else if (x > this->right) {
        nearestX = this->right;
    } else {
        nearestX = x;
    }
    
    Coordinate nearestY;
    if (y < this->bottom) {
        nearestY = this->bottom;
    } else if (y > this->top) {
        nearestY = this->top;
    } else {
        nearestY = y;
    }
    
    return std::pow(nearestX-x,2) + std::pow(nearestY-y,2);
}

EuclideanDistanceSqrd Rectangle::getMinMaxDistSqrd(Coordinate x, Coordinate y)
{
    Coordinate furthestX, furthestY, nearestX, nearestY;

    // Get the corner of the rectangle that is further from the pointers
    // Note: This must necessarily be one of the corners
    EuclideanDistance midX = left + (right-left)/2, midY = bottom + (top-bottom)/2;

    if (x >= midX) {
        furthestX = left;
    } else {
        furthestX = right;
    }
    
    if (y >= midY) {
        furthestY = bottom;
    } else {
        furthestY = top;
    }
    
    // Compute the distance to this further point
    EuclideanDistanceSqrd maxDistInXSqrd, maxDistInYSqrd, maxDistSqrd, minMaxDistSqrd, currMinMaxDistSqrd;
    maxDistInXSqrd = std::pow(furthestX-x,2);
    maxDistInYSqrd = std::pow(furthestY-y,2);
    maxDistSqrd = maxDistInXSqrd + maxDistInYSqrd;
    
    // Remove the x (and then y) component of the distance to the furthest point
    // and replace it with the x (and then y component) of the nearest point
    // in the x (and then y) dimension) to calculate minmaxdist
    if (x <= midX) {
        nearestX = left;
    } else {
        nearestX = right;
    }
    minMaxDistSqrd = maxDistSqrd - maxDistInXSqrd + std::pow(nearestX-x,2);
    
    if (y <= midX) {
        nearestY = bottom;
    } else {
        nearestY = top;
    }
    currMinMaxDistSqrd = maxDistSqrd - maxDistInYSqrd + std::pow(nearestY-y,2);
    
    if (currMinMaxDistSqrd < minMaxDistSqrd) {
        minMaxDistSqrd = currMinMaxDistSqrd;
    }
    
    return minMaxDistSqrd;
}

// This assumes that the rectage is a point (therefore both corners are the same point)
EuclideanDistanceSqrd Rectangle::getObjectDistSqrd(Coordinate x, Coordinate y)
{
    return std::pow(left-x,2) + std::pow(bottom-y,2);
}

bool Rectangle::intersects(Rectangle& rect)
{
    if (left < rect.right && right > rect.left &&
        bottom < rect.top && top > rect.bottom) {
        return true;
    } else {
        return false;
    }
}

void Rectangle::printRectangle()
{
    std::cout << "LB = (" << this->left << "," << this->bottom << "), RT = (" << this->right << "," << this->top << ")" << std::endl;
}

double Rectangle::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(Coordinate)*6;
    return memoryUsage;
}

double Rectangle::computeMemoryUsageBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    return memoryUsage;
}
