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

#ifndef _BINARYHEAP_H
#define _BINARYHEAP_H

#include "MinPriorityQueue.h"

#include <vector>

// Can only be used for primitive types for Key or complex types
// which have defined comparison operators

template <typename Key, typename Element>
class BinaryMinHeap: public MinPriorityQueue<Key,Element> {

    public:
        void init(int n);
        void insert(Element element, Key key);
        unsigned int size();
        void clear();
        Key getMinKey();
        Element extractMinElement();
        void print();

    private:
        void siftUp(int heapIndex);
        void siftDown(int heapIndex);
        int getLeftChildIndex(int parentIndex);
        int getRightChildIndex(int parentIndex);
        int getParentIndex(int childIndex);
        void swapKeysAndElements(int indexA, int indexB);
        struct ElementKey {
            Element element;
            Key key;
            ElementKey(Element _element, Key _key): element(_element), key(_key) {}
            inline bool operator<(ElementKey const& rhs) const { return key > rhs.key; }
        };
        
        std::vector<ElementKey> bHeap;

};

#endif // _BINARYHEAP_H