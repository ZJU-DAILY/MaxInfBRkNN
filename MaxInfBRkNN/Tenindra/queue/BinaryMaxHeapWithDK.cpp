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

#include "BinaryMaxHeapWithDK.h"

#include <iostream>
#include <algorithm>

template <typename Key, typename Element>
void BinaryMaxHeapWithDK<Key,Element>::init(int n)
{
    this->bHeap.reserve(n);
    this->dataElementToHeapIndex.reserve(n);
}

template <typename Key, typename Element>
void BinaryMaxHeapWithDK<Key,Element>::insert(Element element, Key key)
{
    // Add element at the end of the heap and siftUp 
    // until it is in the correct position
    this->bHeap.push_back(std::make_pair(std::move(element),key));
    int lastHeapIdx = bHeap.size()-1;
    this->dataElementToHeapIndex[element] = lastHeapIdx;
    this->siftUp(lastHeapIdx);
}

template <typename Key, typename Element>
Key BinaryMaxHeapWithDK<Key,Element>::getKey(Element element)
{
    // Note: We assume caller has used contains and verified the
    // element actually exists in the queue
    return this->bHeap[this->dataElementToHeapIndex[element]].second;
}

template <typename Key, typename Element>
void BinaryMaxHeapWithDK<Key,Element>::increaseKey(Element element, Key newKey)
{
    // Note: We assume caller has used getKey and checked new key is
    // larger than than current key
    int heapIndex = this->dataElementToHeapIndex[element];
    this->bHeap[heapIndex].second = newKey;
    this->siftUp(heapIndex);
}

template <typename Key, typename Element>
void BinaryMaxHeapWithDK<Key,Element>::decreaseKey(Element element, Key newKey)
{
    // Note: We assume caller has used getKey and checked new key is
    // smaller than than current key
    int heapIndex = this->dataElementToHeapIndex[element];
    this->bHeap[heapIndex].second = newKey;
    this->siftDown(heapIndex);
}

template <typename Key, typename Element>
bool BinaryMaxHeapWithDK<Key,Element>::contains(Element element)
{
    if (this->dataElementToHeapIndex.find(element) != this->dataElementToHeapIndex.end()) {
        return true;
    } else {
        return false;
    }
}

template <typename Key, typename Element>
void BinaryMaxHeapWithDK<Key,Element>::clear()
{
    this->bHeap.clear();
    this->dataElementToHeapIndex.clear();    
}

template <typename Key, typename Element>
void BinaryMaxHeapWithDK<Key,Element>::print()
{
    std::cout << "Priority Queue = [";
    for (std::size_t i = 0; i < bHeap.size(); ++i) {
        if (i != 0) {
            std::cout << ", ";
            if (i % 5 == 0) {
                std::cout << std::endl << "                  ";
            }
        }
        std::cout << "(" << bHeap[i].first << "," << bHeap[i].second << ")";
    }
    std::cout << "]" << std::endl;
}

template <typename Key, typename Element>
unsigned int BinaryMaxHeapWithDK<Key,Element>::size()
{
    return this->bHeap.size();
}

template <typename Key, typename Element>
Key BinaryMaxHeapWithDK<Key,Element>::getMaxKey()
{
    return this->bHeap[0].second;
}

template <typename Key, typename Element>
Element BinaryMaxHeapWithDK<Key,Element>::extractMaxElement()
{
    // Move the last element to the root
    Element maxElement = std::move(this->bHeap[0].first);
    int lastHeapIdx = this->bHeap.size()-1;
    std::pair<Element,Key> lastHeapPair = std::move(this->bHeap[lastHeapIdx]);
    this->bHeap.pop_back();
    this->bHeap[0] = lastHeapPair;
    this->dataElementToHeapIndex[this->bHeap[0].first] = 0;
    this->dataElementToHeapIndex.erase(maxElement);
    
    // siftDown the new root element until it is in the correct position
    this->siftDown(0);
    
    return maxElement;      
}

template <typename Key, typename Element>
void BinaryMaxHeapWithDK<Key,Element>::deleteElement(Element element)
{
    // Note: We assume caller has used contains and verified their
    // element actually exists in the queue

    // Move the last element to the given element index (overwriting
    // delete element in process)
    int elementIdx = this->dataElementToHeapIndex[element];
    int lastHeapIdx = this->bHeap.size()-1;
    std::pair<Element,Key> lastHeapPair = std::move(this->bHeap[lastHeapIdx]);
    this->bHeap.pop_back();
    this->bHeap[elementIdx] = lastHeapPair;
    this->dataElementToHeapIndex[this->bHeap[elementIdx].first] = elementIdx;
    this->dataElementToHeapIndex.erase(element);
    
    // We now re-heapify by moving swapped element to its (new) correct position
    // But if the swapped element is larger than parent, then we need to sift up
    int parentIdx = this->getParentIndex(elementIdx);
    if (this->bHeap[parentIdx].second < this->bHeap[elementIdx].second) {
        this->siftUp(elementIdx);
    } else {
        this->siftDown(elementIdx);
    }    
}

// This may be destructive depending on implementation to below
template <typename Key, typename Element>
void BinaryMaxHeapWithDK<Key,Element>::populateKNNs(std::vector<Element>& elements, std::vector<Key>& keys)
{
    // Sort heap key value in increasing order and then copy elements one by one
    std::sort(this->bHeap.begin(),this->bHeap.end(),SortPair());
    for (std::size_t i = 0; i < this->bHeap.size(); ++i) {
        elements.push_back(this->bHeap[i].first);
        keys.push_back(this->bHeap[i].second);
    }
}

/*
 * Private Methods
 */

template <typename Key, typename Element>
void BinaryMaxHeapWithDK<Key,Element>::siftUp(int heapIndex)
{
    int parentIndex = this->getParentIndex(heapIndex);
    
    // If parent is the smaller than the element we are sifting we swap
    // and it is also still in the heap (i.e. index >= 0)
    // Note: That if the child is larger than the parent, it must be larger
    // than the other child (the other child conforms to the heap property)
    if (parentIndex >= 0 && bHeap[parentIndex].second < bHeap[heapIndex].second) {
        this->swapKeysAndElements(parentIndex,heapIndex);
        this->siftUp(parentIndex);
    }
}

template <typename Key, typename Element>
void BinaryMaxHeapWithDK<Key,Element>::siftDown(int heapIndex)
{
    int leftChildIndex = this->getLeftChildIndex(heapIndex);
    int rightChildIndex = this->getRightChildIndex(heapIndex);
    int heapSize =  this->bHeap.size();
    int largestChildIdx = -1;

    if (rightChildIndex < heapSize) {
        // If the parent has two children pick the largest
        if (this->bHeap[rightChildIndex].second > this->bHeap[leftChildIndex].second) {
            largestChildIdx = rightChildIndex;
        } else {
            largestChildIdx = leftChildIndex;
        }
    } else if (leftChildIndex < heapSize) {
        // The parent has only one child (it must be the largest)
        largestChildIdx = leftChildIndex;
    } 
    
    // If the parent has at least 1 child and it is larger than the parent we swap
    if (largestChildIdx != -1 && this->bHeap[heapIndex].second < this->bHeap[largestChildIdx].second) {
        this->swapKeysAndElements(largestChildIdx,heapIndex);            
        this->siftDown(largestChildIdx); // Continue sifting down the original parent (now swapped with child)
    }
}

template <typename Key, typename Element>
void BinaryMaxHeapWithDK<Key,Element>::swapKeysAndElements(int indexA, int indexB)
{
    using std::swap; // Use argument-dependent lookup to find best swap
    swap(this->bHeap[indexA],this->bHeap[indexB]);
    // Update mapping of Elements to their corresponding heap index
    this->dataElementToHeapIndex[this->bHeap[indexB].first] = indexB;
    this->dataElementToHeapIndex[this->bHeap[indexA].first] = indexA;
}

template <typename Key, typename Element>
int BinaryMaxHeapWithDK<Key,Element>::getLeftChildIndex(int parentIndex)
{
    return parentIndex*2 + 1;
}

template <typename Key, typename Element>
int BinaryMaxHeapWithDK<Key,Element>::getRightChildIndex(int parentIndex)
{
    return parentIndex*2 + 2;
}

template <typename Key, typename Element>
int BinaryMaxHeapWithDK<Key,Element>::getParentIndex(int childIndex)
{
    // Note: static_case will floor() the expression
    return static_cast<int>((childIndex-1)/2);
}

template class BinaryMaxHeapWithDK<DistanceBound,NodeID>;