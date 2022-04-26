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

#include "BinaryMinHeapWithDK.h"

#include <iostream>

template <typename Key, typename Element>
void BinaryMinHeapWithDK<Key,Element>::init(int n)
{
    this->bHeap.reserve(n);
    this->dataElementToHeapIndex.reserve(n);
}

template <typename Key, typename Element>
void BinaryMinHeapWithDK<Key,Element>::insert(Element element, Key key)
{
    // Add element at the end of the heap and siftUp 
    // until it is in the correct position
    this->bHeap.push_back(std::make_pair(std::move(element),key));
    int lastHeapIdx = bHeap.size()-1;
    this->dataElementToHeapIndex[element] = lastHeapIdx;
    this->siftUp(lastHeapIdx);    
}

template <typename Key, typename Element>
Key BinaryMinHeapWithDK<Key,Element>::getKey(Element element)
{
    // Note: We assume caller has used contains and verified the
    // element actually exists in the queue
    return this->bHeap[this->dataElementToHeapIndex[element]].second;
}

template <typename Key, typename Element>
void BinaryMinHeapWithDK<Key,Element>::decreaseKey(Element element, Key newKey)
{
    // Note: We assume caller has used getKey and checked new key is
    // smaller than than current key
    int heapIndex = this->dataElementToHeapIndex[element];
    // Update key and siftUp until element is in correct position
    // Note: This is valid because we will only make it smaller
    // and so it can only go up the binary heap not down
    this->bHeap[heapIndex].second = newKey;
    this->siftUp(heapIndex);    
}

template <typename Key, typename Element>
bool BinaryMinHeapWithDK<Key,Element>::contains(Element element)
{
    if (this->dataElementToHeapIndex.find(element) != this->dataElementToHeapIndex.end()) {
        return true;
    } else {
        return false;
    }
}

template <typename Key, typename Element>
void BinaryMinHeapWithDK<Key,Element>::clear()
{
    this->bHeap.clear();
    this->dataElementToHeapIndex.clear();    
}

template <typename Key, typename Element>
void BinaryMinHeapWithDK<Key,Element>::print()
{
    std::cout << "Not implemented yet" << std::endl;
}

template <typename Key, typename Element>
unsigned int BinaryMinHeapWithDK<Key,Element>::size()
{
    return this->bHeap.size();
}

template <typename Key, typename Element>
Key BinaryMinHeapWithDK<Key,Element>::getMinKey()
{
    return this->bHeap[0].second;
}

template <typename Key, typename Element>
Element BinaryMinHeapWithDK<Key,Element>::extractMinElement()
{
    Element minElement = std::move(this->bHeap[0].first);
    // Move the last element to the root 
    int lastHeapIdx = this->bHeap.size()-1;
    if (lastHeapIdx != 0) {
        // siftDown the new root element until it is in the correct position
        std::pair<Element,Key> lastHeapPair = std::move(this->bHeap[lastHeapIdx]);
        this->bHeap.pop_back();
        this->bHeap[0] = lastHeapPair;
        this->dataElementToHeapIndex[this->bHeap[0].first] = 0;
        this->dataElementToHeapIndex.erase(minElement);
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
void BinaryMinHeapWithDK<Key,Element>::siftUp(int heapIndex)
{
    // If parent is the larger than the element we are sifting we swap
    // and it is also still in the heap (i.e. index >= 0)
    // Note: That if the child is smaller than the parent, it must be smaller
    // than the other child (the other child conforms to the heap property)
    int parentIndex = this->getParentIndex(heapIndex);
    if (parentIndex >= 0 && bHeap[parentIndex].second > bHeap[heapIndex].second) {
        this->swapKeysAndElements(parentIndex,heapIndex);
        this->siftUp(parentIndex);
    }
}

template <typename Key, typename Element>
void BinaryMinHeapWithDK<Key,Element>::siftDown(int heapIndex)
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
        if (this->bHeap[rightChildIndex].second < this->bHeap[leftChildIndex].second) {
            smallestChildIndex = rightChildIndex;
        } else {
            smallestChildIndex = leftChildIndex;
        }
    } else if (leftChildIndex < heapSize) {
        // The parent has only one child (it must be smallest)
        smallestChildIndex = leftChildIndex;
    } 
    
    // If the parent has at least 1 child and it is smaller than the parent we swap
    if (smallestChildIndex != -1 && this->bHeap[heapIndex].second > this->bHeap[smallestChildIndex].second) {
        this->swapKeysAndElements(smallestChildIndex,heapIndex);            
        this->siftDown(smallestChildIndex); // Continue sifting down the original parent (now swapped with child)
    }
}

template <typename Key, typename Element>
void BinaryMinHeapWithDK<Key,Element>::swapKeysAndElements(int indexA, int indexB)
{
    using std::swap; // Use argument-dependent lookup to find best swap
    swap(this->bHeap[indexA],this->bHeap[indexB]);
    // Update mapping of Elements to their corresponding heap index
    this->dataElementToHeapIndex[this->bHeap[indexB].first] = indexB;
    this->dataElementToHeapIndex[this->bHeap[indexA].first] = indexA;
}

template <typename Key, typename Element>
int BinaryMinHeapWithDK<Key,Element>::getLeftChildIndex(int parentIndex)
{
    return parentIndex*2 + 1;
}

template <typename Key, typename Element>
int BinaryMinHeapWithDK<Key,Element>::getRightChildIndex(int parentIndex)
{
    return parentIndex*2 + 2;
}

template <typename Key, typename Element>
int BinaryMinHeapWithDK<Key,Element>::getParentIndex(int childIndex)
{
    // Note: static_case will floor() the expression
    return (childIndex-1)/2;  // In C++ integer division is floored so we don't need to call std::floor or static_cast
}

template class BinaryMinHeapWithDK<EdgeWeight,NodeID>;