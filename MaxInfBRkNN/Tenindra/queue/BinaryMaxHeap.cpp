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

#include "BinaryMaxHeap.h"

#include <iostream>
#include <algorithm>

template <typename Key, typename Element>
void BinaryMaxHeap<Key,Element>::init(int n)
{
    this->bHeap.reserve(n);
}

template <typename Key, typename Element>
void BinaryMaxHeap<Key,Element>::insert(Element element, Key key)
{
    // Add element at the end of the heap and siftUp 
    // until it is in the correct position
    this->bHeap.push_back(ElementKey(std::move(element),key));
    int lastHeapIdx = bHeap.size()-1;
    this->siftUp(lastHeapIdx);
}

template <typename Key, typename Element>
void BinaryMaxHeap<Key,Element>::clear()
{
    // Add element at the end of the heap and siftUp 
    // until it is in the correct position
    this->bHeap.clear();
}

template <typename Key, typename Element>
void BinaryMaxHeap<Key,Element>::print()
{
    std::cout << "Priority Queue = [";
    for (std::size_t i = 0; i < bHeap.size(); ++i) {
        if (i != 0) {
            std::cout << ", ";
            if (i % 5 == 0) {
                std::cout << std::endl << "                  ";
            }
        }
        std::cout << "(" << bHeap[i].element << "," << bHeap[i].key << ")";
    }
    std::cout << "]" << std::endl;
}

template <typename Key, typename Element>
unsigned int BinaryMaxHeap<Key,Element>::size()
{
    return this->bHeap.size();
}

template <typename Key, typename Element>
Key BinaryMaxHeap<Key,Element>::getMaxKey()
{
    return this->bHeap[0].key;
}

template <typename Key, typename Element>
Element BinaryMaxHeap<Key,Element>::extractMaxElement()
{
    // Move the last element to the root
    Element maxElement = std::move(this->bHeap[0].element);
    int lastHeapIdx = this->bHeap.size()-1;
    ElementKey lastHeapPair = std::move(this->bHeap[lastHeapIdx]);
    this->bHeap.pop_back();
    this->bHeap[0] = lastHeapPair;
    
    // siftDown the new root element until it is in the correct position
    this->siftDown(0);
    
    return maxElement;      
}

template <typename Key, typename Element>
void BinaryMaxHeap<Key,Element>::deleteElement(Element element)
{
    // We delete the element if it can be found (the caller should
    // avoid trying to delete if heap key is greater than max key)
    bool found = false;
    int elementIdx = 0;
    for (std::size_t i = 0; i < this->bHeap.size(); ++i) {
        if (this->bHeap[i].element == element) {
            found = true;
            elementIdx = i;
            break;
        }
    }
    
    if (found) {
        // Move the last element to the given element index (overwriting
        // delete element in process)
        int lastHeapIdx = this->bHeap.size()-1;
        ElementKey lastHeapPair = std::move(this->bHeap[lastHeapIdx]);
        this->bHeap.pop_back();
        this->bHeap[elementIdx] = lastHeapPair;
        
        int parentIdx = this->getParentIndex(elementIdx);
        if (this->bHeap[parentIdx].key < this->bHeap[elementIdx].key) {
            this->siftUp(elementIdx);
        } else {
            this->siftDown(elementIdx);
        }
    }
}

// This may be destructive depending on implementation to below
template <typename Key, typename Element>
void BinaryMaxHeap<Key,Element>::populateKNNs(std::vector<Element>& elements, std::vector<Key>& keys)
{
    // Sort heap key value in increasing order and then copy elements one by one
    std::sort(this->bHeap.begin(),this->bHeap.end());
    for (std::size_t i = 0; i < this->bHeap.size(); ++i) {
        elements.push_back(this->bHeap[i].element);
        keys.push_back(this->bHeap[i].key);
    }
}

/*
 * Private Methods
 */

template <typename Key, typename Element>
void BinaryMaxHeap<Key,Element>::siftUp(int heapIndex)
{
    int parentIndex = this->getParentIndex(heapIndex);
    
    // If parent is the smaller than the element we are sifting we swap
    // and it is also still in the heap (i.e. index >= 0)
    // Note: That if the child is larger than the parent, it must be larger
    // than the other child (the other child conforms to the heap property)
    if (parentIndex >= 0 && bHeap[parentIndex].key < bHeap[heapIndex].key) {
        this->swapKeysAndElements(parentIndex,heapIndex);
        this->siftUp(parentIndex);
    }
}

template <typename Key, typename Element>
void BinaryMaxHeap<Key,Element>::siftDown(int heapIndex)
{
    int leftChildIndex = this->getLeftChildIndex(heapIndex);
    int rightChildIndex = this->getRightChildIndex(heapIndex);
    int heapSize =  this->bHeap.size();
    int largestChildIdx = -1;

    if (rightChildIndex < heapSize) {
        // If the parent has two children pick the largest
        if (this->bHeap[rightChildIndex].key > this->bHeap[leftChildIndex].key) {
            largestChildIdx = rightChildIndex;
        } else {
            largestChildIdx = leftChildIndex;
        }
    } else if (leftChildIndex < heapSize) {
        // The parent has only one child (it must be the largest)
        largestChildIdx = leftChildIndex;
    } 
    
    // If the parent has at least 1 child and it is larger than the parent we swap
    if (largestChildIdx != -1 && this->bHeap[heapIndex].key < this->bHeap[largestChildIdx].key) {
        this->swapKeysAndElements(largestChildIdx,heapIndex);            
        this->siftDown(largestChildIdx); // Continue sifting down the original parent (now swapped with child)
    }
}

template <typename Key, typename Element>
void BinaryMaxHeap<Key,Element>::swapKeysAndElements(int indexA, int indexB)
{
    // Update mapping of Elements to their corresponding heap index
    using std::swap; // Use argument-dependent lookup to find best swap
    swap(this->bHeap[indexA],this->bHeap[indexB]);
}

template <typename Key, typename Element>
int BinaryMaxHeap<Key,Element>::getLeftChildIndex(int parentIndex)
{
    return parentIndex*2 + 1;
}

template <typename Key, typename Element>
int BinaryMaxHeap<Key,Element>::getRightChildIndex(int parentIndex)
{
    return parentIndex*2 + 2;
}

template <typename Key, typename Element>
int BinaryMaxHeap<Key,Element>::getParentIndex(int childIndex)
{
    // Note: static_case will floor() the expression
    return static_cast<int>((childIndex-1)/2);
}

template class BinaryMaxHeap<DistanceBound,NodeID>;