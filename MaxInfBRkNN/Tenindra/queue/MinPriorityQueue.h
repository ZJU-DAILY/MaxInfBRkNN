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

#ifndef _MINPRIORITYQUEUE_H
#define _MINPRIORITYQUEUE_H

#include "../common.h"

template <typename Key, typename Element>
class MinPriorityQueue {

    public:
        virtual ~MinPriorityQueue() {}
        // Queue should be able to allocate memory it's implementation needs
        // (e.g. reserve calls for vectors) in advance so that reallocation
        // can be avoided. Implement init() to call in constuctor().
        virtual void init(int n) = 0;
        virtual void insert(Element element, Key key) = 0;
        virtual Element extractMinElement() = 0;
        virtual Key getMinKey() = 0;
        virtual unsigned int size() = 0;
        virtual void clear() = 0;
        virtual void print() = 0;
};

#endif // _MINPRIORITYQUEUE_H