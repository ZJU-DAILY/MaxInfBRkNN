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

#ifndef _BINARYMINHEAPWITHDK_H
#define _BINARYMINHEAPWITHDK_H

#include "MinPriorityQueueWithDK.h"

#include <unordered_map>
#include <vector>

// Can only be used for primitive types for Element and Key
// (because hash function to unordered_map is not custom)

template <typename Key, typename Element>
class BinaryMinHeapWithDK: public MinPriorityQueueWithDK<Key,Element> {

    public:
        void init(int n);
        void insert(Element element, Key key);
        Key getKey(Element element);
        void decreaseKey(Element element, Key key);
        bool contains(Element element);
        unsigned int size();
        void clear();
        void print();
        Key getMinKey();
        Element extractMinElement();

    private:
        void siftUp(int heapIndex);
        void siftDown(int heapIndex);
        int getLeftChildIndex(int parentIndex);
        int getRightChildIndex(int parentIndex);
        int getParentIndex(int childIndex);
        void swapKeysAndElements(int indexA, int indexB);
        
//         std::vector<Key> heap;
//         std::vector<Element> dataElements;
        std::vector<std::pair<Element,Key>> bHeap;
        std::unordered_map<Element,int> dataElementToHeapIndex;

};

#endif // _BINARYMINHEAPWITHDK_H