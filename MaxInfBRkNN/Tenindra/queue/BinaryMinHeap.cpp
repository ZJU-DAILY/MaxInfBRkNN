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

#include "BinaryMinHeap.h"

#include "../processing/MortonList.h"
#include "../processing/StaticRtree.h"
#include "../processing/AStarSearch.h"

#include <iostream>

// Note: Wherever possible we try to use move semantics on Element
// since it can be object or complex types. We expect Key to be 
// primitive type so we do not force move semantics on that

template <typename Key, typename Element>
void BinaryMinHeap<Key,Element>::init(int n)
{
    this->bHeap.reserve(n);
}

template <typename Key, typename Element>
void BinaryMinHeap<Key,Element>::insert(Element element, Key key)
{
    // Add element at the end of the heap and siftUp 
    // until it is in the correct position
    this->bHeap.push_back(ElementKey(std::move(element),key));
    this->siftUp(bHeap.size()-1);
}

template <typename Key, typename Element>
void BinaryMinHeap<Key,Element>::clear()
{
    this->bHeap.clear();
}

template <typename Key, typename Element>
void BinaryMinHeap<Key,Element>::print()
{
    std::cout << "Not implemented yet" << std::endl;
}

template <typename Key, typename Element>
unsigned int BinaryMinHeap<Key,Element>::size()
{
    return this->bHeap.size();
}

template <typename Key, typename Element>
Key BinaryMinHeap<Key,Element>::getMinKey()
{
    return this->bHeap[0].key;
}

template <typename Key, typename Element>
Element BinaryMinHeap<Key,Element>::extractMinElement()
{
    Element minElement = std::move(this->bHeap[0].element);
    // Move the last element to the root
    int lastHeapIdx = this->bHeap.size()-1;
    if (lastHeapIdx != 0) {
        // siftDown the new root element until it is in the correct position
        ElementKey lastHeapPair = std::move(this->bHeap[lastHeapIdx]);
        this->bHeap.pop_back();
        this->bHeap[0] = lastHeapPair;
        this->siftDown(0);
    } else {
        this->bHeap.pop_back();
    }
    return minElement;
}

/*
 * Private Methods
 */

template <typename Key, typename Element>
void BinaryMinHeap<Key,Element>::siftUp(int heapIndex)
{
    int parentIndex = this->getParentIndex(heapIndex);
    // If parent is the larger than the element we are sifting we swap
    // and it is also still in the heap (i.e. index >= 0)
    // Note: That if the child is smaller than the parent, it must be smaller
    // than the other child (the other child conforms to the heap property)
    if (heapIndex != 0 && bHeap[parentIndex].key > bHeap[heapIndex].key) {
        this->swapKeysAndElements(parentIndex,heapIndex);
        this->siftUp(parentIndex);
    }
}

template <typename Key, typename Element>
void BinaryMinHeap<Key,Element>::siftDown(int heapIndex)
{
    int leftChildIndex = this->getLeftChildIndex(heapIndex);
    int rightChildIndex = this->getRightChildIndex(heapIndex);
    int heapSize =  this->bHeap.size();
    int smallestChildIndex = -1;

    // Note: In practice rightChildIndex is usually smaller than heapSize
    // because we siftDown from the start of queue (it is only larger at
    // the end of the queue). So by structuring our code in this way
    // we can ensure we take advantage of any branch prediction benefits.
    if (rightChildIndex < heapSize) {
        // If the parent has two children pick the smallest
        if (this->bHeap[rightChildIndex].key < this->bHeap[leftChildIndex].key) {
            smallestChildIndex = rightChildIndex;
        } else {
            smallestChildIndex = leftChildIndex;
        }
    } else if (leftChildIndex < heapSize) {
        // The parent has only one child (it must be smallest)
        smallestChildIndex = leftChildIndex;
    } 
    
    // If the parent has at least 1 child and it is smaller than the parent we swap
    if (smallestChildIndex != -1 && this->bHeap[heapIndex].key > this->bHeap[smallestChildIndex].key) {
        this->swapKeysAndElements(smallestChildIndex,heapIndex);            
        this->siftDown(smallestChildIndex); // Continue sifting down the original parent (now swapped with child)
    }
}

template <typename Key, typename Element>
void BinaryMinHeap<Key,Element>::swapKeysAndElements(int indexA, int indexB)
{
    using std::swap; // Use argument-dependent lookup to find best swap
    swap(this->bHeap[indexA],this->bHeap[indexB]);
}

template <typename Key, typename Element>
int BinaryMinHeap<Key,Element>::getLeftChildIndex(int parentIndex)
{
    return parentIndex*2 + 1;
}

template <typename Key, typename Element>
int BinaryMinHeap<Key,Element>::getRightChildIndex(int parentIndex)
{
    return parentIndex*2 + 2;
}

template <typename Key, typename Element>
int BinaryMinHeap<Key,Element>::getParentIndex(int childIndex)
{
    // Note: static_case will floor() the expression
    return (childIndex-1)/2; // In C++ integer division is floored so we don't need to call std::floor or static_cast
}

template class BinaryMinHeap<DistanceBound,DataTuple>;
template class BinaryMinHeap<EuclideanDistanceSqrd,RtreeDataTuple>;
template class BinaryMinHeap<EdgeWeight,NodeID>;
template class BinaryMinHeap<LongPathDistance,NodeID>;
template class BinaryMinHeap<EdgeWeight,NodeLinkPair>;
template class BinaryMinHeap<EdgeWeight,NodeStatusPair>;
template class BinaryMinHeap<EdgeWeight,NodePair>;
template class BinaryMinHeap<EdgeWeight,AStarHeapElement>;