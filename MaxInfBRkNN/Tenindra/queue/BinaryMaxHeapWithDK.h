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

#ifndef _BINARYMAXHEAPWITHDK_H
#define _BINARYMAXHEAPWITHDK_H

#include "../common.h"

#include <unordered_map>
#include <vector>

// Can only be used for primitive types for Element and Key
// (because hash function to unordered_map is not custom)

template <typename Key, typename Element>
class BinaryMaxHeapWithDK {

    public:
        void init(int n);
        void insert(Element element, Key key);
        Key getKey(Element element);
        void increaseKey(Element element, Key key);
        void decreaseKey(Element element, Key key);
        bool contains(Element element);
        unsigned int size();
        void clear();
        void print();
        Key getMaxKey();
        Element extractMaxElement();
        void deleteElement(Element element);
        // Populate kNNs in increasing order
        void populateKNNs(std::vector<Element>& elements, std::vector<Key>& keys);
        
    private:
        void siftUp(int heapIndex);
        void siftDown(int heapIndex);
        int getLeftChildIndex(int parentIndex);
        int getRightChildIndex(int parentIndex);
        int getParentIndex(int childIndex);
        void swapKeysAndElements(int indexA, int indexB);
        class SortPair {
            public:
                // Sort in increasing order (for use with populating kNNs)
                int operator()(const std::pair<Element,Key>& pair1, const std::pair<Element,Key>& pair2)
                {
                    if (pair1.second == pair2.second) {
                        return pair1.first < pair2.first;
                    } else {
                        return pair1.second < pair2.second;
                    }
                }    
        };
        
        std::vector<std::pair<Element,Key>> bHeap;
        std::unordered_map<Element,int> dataElementToHeapIndex;

};

#endif // _BINARYMAXHEAPWITHDK_H