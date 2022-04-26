/* Original Work Copyright (C) 2005, 2006, 2007, 2008 Robert Geisberger, Dominik Schultes, Peter Sanders, Universitaet Karlsruhe (TH)
 * Modified Work Copyright (C) 2012 Lingkun Wu, Xiaokui Xiao, Dingxiong Deng, Gao Cong, Andy Diwen Zhu, Shuigeng Zhou
 *
 * This file is part of Road Network kNN Experimental Evaluation.
 * 
 * The following file is a derivative work of the code from 
 * http://sourceforge.net/projects/ntu-sp-exp/ developed for the paper below.
 * The authors have requested users of this code to cite the paper.
 * 
 * Lingkun Wu, Xiaokui Xiao, Dingxiong Deng, Gao Cong, Andy Diwen Zhu, Shuigeng Zhou
 * Shortest Path and Distance Queries on Road Networks: An Experimental Evaluation
 * PVLDB 5(5): 406-417 (2012)
 *
 * That work is in turn a derivative work of the code from
 * Contraction Hierarchies (http://algo2.iti.kit.edu/english/routeplanning.php),
 * licensed under AGPLv3. Thus this work is also licensed under that license.
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

#ifndef BINARYHEAP_H
#define BINARYHEAP_H

#include <vector>

	/**
	 * In case of the (FIFO)BinaryHeap, we distinguish between an ExternalKey,
	 * which is provided by the application (which uses the heap) and an internal Key,
	 * which is used inside the heap. The ExternalKey is always a part of the internal Key.
	 * A "key extractor" extracts the ExternalKey out of the internal Key.
	 * 
	 * In the simple case of a BinaryHeap without the FIFO property, the internal Key just
	 * corresponds to the ExternalKey. Thus, the SimpleKeyExtractor just returns the
	 * given internal Key as the ExternalKey.
	 * @param ExternalKey the type of the external key = the type of the internal key
	 */
template < typename ExternalKey >
class SimpleKeyExtractor
{
public:
    /** Returns the given internal Key as the ExternalKey. */
    static ExternalKey key(ExternalKey k) {return k;}
};


	/**
	 * Represents an element of a BinaryHeap.
	 * @param ExternalKey the type of the external key (see SimpleKeyExtractor)
	 * @param MetaExtKey class that provides data about the external key, e.g. the max value
	 * @param Data the type of the application-specific data that is associated with this element
	 * @param Count the type of the counter that counts the number of heap elements / 
	 *              the number of heap operations
	 */
template < typename ExternalKey,
           typename MetaExtKey,
           typename Data,
           typename Count >
class BinaryHeapElement
{
public:
    /** Constructs a DUMMY element. */
    BinaryHeapElement() : _key(MetaExtKey::MAX_VALUE), _index(0) {}
            
    /**
     * Constructor.
     * @param key key of the new element
     * @param index index of the new element within the heap
     */
    BinaryHeapElement(ExternalKey key, Count index) : _key(key), _index(index) {}
    
    
    /** Returns the external key of this element. */
    ExternalKey key() const {return _key;}
    
    /** Sets the external key of this element. */
    void key(ExternalKey newKey) {_key = newKey;}

    bool isDummy() const {return key() == MetaExtKey::MAX_VALUE;}
    
    
    /** Returns a reference to the application-specific data object of this element. */
    Data& data() {return _data;}
    
    /** Returns a reference to the application-specific data object of this element. */
    const Data& data() const {return _data;}
    
	void SetData(const Data& data) { _data = data; }
    
    /** Marks that this element has been deleted from the heap. */
    void markDeleted() {index(0);}
    
    /** Returns true iff this element has been deleted from the heap. */
    bool hasBeenDeleted() const {return (_index==0);}

    
    /** Returns this element's index within the heap. */
    Count index() const {return _index;}
    
    /** Sets this element's index within the heap. */
    void index(Count newIndex) {_index = newIndex;}
    
private:
    /** The external key of this element. */
    ExternalKey _key;    
    /** The application-specific data that is associated with this element. */
    Data _data;
    /**
     * This element's index within the heap.
     * 0 is used to indicate that this element has been deleted from the heap.
     */
    Count _index;
};

	/**
	 * Represents a binary heap.
	 * @param ExternalKey the type of the external key (see SimpleKeyExtractor)
	 * @param MetaExtKey class that provides data about the external key, e.g. the max value
	 * @param Data the type of the application-specific data that is associated with this element
	 * @param Count the type of the counter that counts the number of heap elements / 
	 *              the number of heap operations
	 * @param Key the type of the internal key (see SimpleKeyExtractor)
	 * @param KeyExtractor the "key extractor" that is used (see SimpleKeyExtractor)
	 */
template < typename ExternalKey,
           typename MetaExtKey,
           typename Data,
           typename Count,
           typename Key = ExternalKey,
           typename KeyExtractor = SimpleKeyExtractor<ExternalKey> >
class BinaryHeap
{
public:
    /** The type of the elements stored in this heap. */
    typedef BinaryHeapElement<ExternalKey, MetaExtKey, Data, Count> PQElement;
    
    typedef Data PQData;
    
    /** Constructor. */
    BinaryHeap() {
		//reserve(1000);
        insertDummy();
        _heap.push_back( IndexKey() ); // add a dummy element
        _heap[0].second = MetaExtKey::MIN_VALUE;
    }

    /** Returns the size of this heap (= the number of elements). */
    Count size() const { 
        return _heap.size()-1; // subtract dummy element (index 0)
    }
    
    void reserve( size_t n )
    {
        _elements.reserve(n+1);
        _heap.reserve(n+1);
    }
                            
    /** Returns true iff the heap is empty. */
    bool empty() const { return (size() == 0); }

    /** 
     * Returns the (external) key of the minimum element in this heap
     * or the maximum value of the ExternalKey data type iff the heap is empty.
     */
    ExternalKey min() const {   
        if (empty()) return MetaExtKey::MAX_VALUE;
        return KeyExtractor::key(_heap[1].second);
    }

    /**
     * Creates a new element with the given (internal) key and inserts it into the heap.
     * @return the index of the new element.
     */
    Count insert(Key key) {
        Count element = _elements.size();
        Count index = _heap.size();
        _elements.push_back( PQElement(KeyExtractor::key(key), index) );
       _heap.push_back(IndexKey(element,key));
        upheap(index);        
        return element;
    }
    
	//Count insert(Key key, Data data) {
	//	Count element = _elements.size();
	//	Count index = _heap.size();
	//	_elements.push_back( PQElement(KeyExtractor::key(key), index) );
	//	_heap.push_back(IndexKey(element,key));
	//	_elements[index].SetData(data);
	//	upheap(index);        
	//	return element;
	//}

    /** Replaces a dummy element and inserts the actual element into the heap. */
    void insert(Key key, Count elementIndex) {          
        PQElement& element = _elements[elementIndex];
        Count index = _heap.size();
        element.key(KeyExtractor::key(key));        
        element.index(index);
        _heap.push_back(IndexKey(elementIndex,key));
        upheap(index);
    }
    
    /** Adds a dummy element, which is NOT inserted into the heap, yet. */
    void insertDummy() {        
        _elements.push_back( PQElement() );
    }

    /**
     * Initialises the _elements vector so that it contains the given number of dummy elements.
     * This is useful when a search is performed after another search has not cleaned up.
     * In particular, this is used when the main phase of the landmark query is executed
     * so that the shortest path tree of the initial phase can be kept in order to be able
     * to compose the complete path at the end.
     */
    void initDummies(Count number) {
        assert( empty() );
        _elements.resize(number);
    }
    
    bool isDummy(Count element) const {
        return (_elements[element].isDummy());
    }
    
    /** Returns a reference to the elements. */
    vector<PQElement>& elements() {return _elements;}
    
    /** Returns a reference to the elements. */
    const vector<PQElement>& elements() const {return _elements;}

    /**
     * Returns the index of the minimum element in this heap
     * or 0 iff the heap is empty.
     */
    Count minElement() const {
        if (empty()) return 0;
        return _heap[1].first;
    }

	Data minElementData() const {
		if(empty()) return 0;
		return _elements[_heap[1].first].data();
	}

    /**
     * Deletes the minimum element in this heap.
     * Precondition: The heap is not empty.
     * @return the index of the deleted element
     */
    Count deleteMin() {
        assert( ! empty() );

        const Count element = _heap[1].first;
        _heap[1] = _heap[_heap.size()-1];
        Count index = 1;
        const Count droppingElement = _heap[index].first;
        _heap.pop_back();
        if (size() > 1) {
            // downheap:
            // Move the element at the top downwards
            // until the heap condition is restored.            
            
            Key k = _heap[index].second;
            Count nextIndex = 2*index;
            while (nextIndex < _heap.size()) {
                nextIndex += 
                    ((nextIndex+1 < _heap.size()) &&
                    (_heap[nextIndex].second > _heap[nextIndex+1].second));
                
                assert( _elements[_heap[nextIndex].first].key() == KeyExtractor::key(_heap[nextIndex].second) );
                
                if (k <= _heap[nextIndex].second) break;
                
                _heap[index] = _heap[nextIndex];
                _elements[_heap[nextIndex].first].index(index);
                index = nextIndex;
                nextIndex *= 2;
            }
            _heap[index].first = droppingElement;
            _heap[index].second = k;
            
            // end of downheap
        }
        _elements[droppingElement].index(index);
        
        _elements[element].markDeleted();
        COUNTING( counter.incDouble(COUNT_DEL_MIN) );
        return element;
    }

    /**
     * Deletes an arbitrary element in this heap.
     * Precondition: The heap is not empty.
     * @return the index of the deleted element
     */
    Count deleteArbitrary() {
        assert( ! empty() );         
        Count element = _heap.back().first;
        _heap.pop_back();
        _elements[element].markDeleted();
        return element;
    }

    /**
     * Increases the key of the given element:
     * sets the key to the given value.
     */
    void increaseKey(Count element, Key newKey) {
        Count index = _elements[element].index();
        assert( index < _heap.size() );
        assert( _heap[index].first == element );
        
        _elements[element].key(KeyExtractor::key(newKey));
        _heap[index].second = newKey;
        downheap(index);
        COUNTING( counter.incDouble(COUNT_INCR_KEY) );
    }

	 /**
     * Delete a given element:
     */
    void deleteElement(Count element) {
		if(element > size())
		{
			return;
		}
        Count index = _elements[element].index();
        assert( index < _heap.size() );
        assert( _heap[index].first == element );
        
		Count last = _heap.back().first;

		_elements[element].key(_elements[last].key());
		_elements[element].SetData(_elements[last].data());
		printf("delete %d, %d, after %d %d\n", 
			_elements[last].key(), _elements[last].data(),
			_elements[element].key(), _elements[element].data());
	    _heap[index].second = _heap.back().second;
	    _heap.pop_back();
	    _elements[last].markDeleted();
		
		downheap(index);
    }

    /**
     * Decreases the key of the given element:
     * sets the key to the given value.
     */
    void decreaseKey(Count element, Key newKey) {
        Count index = _elements[element].index();
        assert( index < _heap.size() );
        assert( _heap[index].first == element );
        
        _elements[element].key(KeyExtractor::key(newKey));
        _heap[index].second = newKey;
        upheap(index);
        COUNTING( counter.incDouble(COUNT_DECR_KEY) );
    }
    
    /**
     * Updates the key of the given element:
     * sets the key to the given value.
     * The key can be larger or smaller than the previous key.
     */
    void updateKey(Count element, Key newKey) {
        Count index = _elements[element].index();
        assert( index < _heap.size() );
        assert( _heap[index].first == element );
        Key oldKey = _heap[index].second;

        _elements[element].key(KeyExtractor::key(newKey));
        _heap[index].second = newKey;
        if (newKey < oldKey)
        {        
            upheap(index);
            COUNTING( counter.incDouble(COUNT_DECR_KEY) );
        }
        else if (newKey > oldKey)
        {
            downheap(index);
            COUNTING( counter.incDouble(COUNT_INCR_KEY) );
        }
    }
	NodeID elementsSize()
	{
		return _elements.size();
	}

    void clear() {
        _elements.clear();
		_heap.clear();
		_heap.push_back( IndexKey() ); // add a dummy element
		_heap[0].second = MetaExtKey::MIN_VALUE;
        insertDummy();
    }
    
	~BinaryHeap()
	{
		//_elements.clear();
		//_heap.clear();
		vector <PQElement>().swap(_elements);
		vector <IndexKey>().swap(_heap);
	}

    /**
     * For debugging purpose.
     * @return _heap is still a binary heap.
     */
    bool checkHeapProperty() 
    {
        bool result = true;
        for ( unsigned int i = 2; i < _heap.size(); i++ )
        {
            if (!(_heap[i/2].second <= _heap[i].second))
            {
                cerr << "_heap[" << (i/2) << "] = " << _heap[i/2].second;
                cerr << " > " << _heap[i].second << " = _heap[" << i << "]" << endl;
                result = false;
                break;
            }            
        }
        return result;
    }
    
    unsigned int totalHeapSize() const
	{
		unsigned int pairSize = sizeof(Count) + sizeof(Key) ;
		unsigned int  eSize = sizeof(PQElement) * _elements.size();
		unsigned int hSize = pairSize * _heap.size();
		return eSize + hSize;
	}
private:
    typedef pair<Count, Key> IndexKey;
    
    /**
     * The elements of this heap.
     * The order corresponds to the order of insertion and
     * is not changed by any heap operation. That implies
     * that an index of this vector can be used as a pointer
     * to the corresponding element.
     * It is possible that this vector contains holes (i.e.
     * dummy elements), namely for each element that has been
     * inserted into the other pqueue (for the other search direction)
     * but not to this one.
     * The first element (index 0) is a dummy element so that the index 0
     * can be used to mark elements that have not been inserted into the heap.
     */
    vector<PQElement> _elements;
    
    /** 
     * "Pointers" (first) to the elements of this heap in the right heap order.
     * In addition (second), the internal key of the corresponding element.
     * The first element (index 0) is a dummy element so that the index 0
     * can be used to mark elements that have been deleted from the heap.
     */
    vector<IndexKey> _heap;
    
    
    

    /** 
     * Move the element with the given index upwards
     * until the heap condition is restored.
     */
    void upheap(Count index) {
        Count risingElement = _heap[index].first;
        Key k = _heap[index].second;
        while (_heap[index / 2].second > k) {
            assert( index > 1 );
            assert( _elements[_heap[index / 2].first].key() == KeyExtractor::key(_heap[index / 2].second) );
            _heap[index] = _heap[index / 2];
            _elements[_heap[index].first].index(index);
            index /= 2;
        }
        _heap[index].first = risingElement;
        _heap[index].second = k;
        _elements[risingElement].index(index);
    }

    /** 
     * Move the element with the given index downwards.
     * until the heap condition is restored. 
     */
    void downheap(Count index) {
        Count descendingElement = _heap[index].first;
        Key k = _heap[index].second;        
        Count maxIndex;
        if (2*index < _heap.size() && _heap[2*index].second < k) 
        {
            maxIndex = 2*index;
        }
        else
        {
            maxIndex = index;
        }
        if ((2*index + 1) < _heap.size() && _heap[2*index+1].second < _heap[maxIndex].second) 
        {
            maxIndex = 2*index + 1;   
        }
        while (maxIndex != index) {
            assert( index >= 1 );
            assert( _elements[_heap[maxIndex].first].key() == KeyExtractor::key(_heap[maxIndex].second) );
            _heap[index] = _heap[maxIndex];
            _elements[_heap[index].first].index(index);
            index = maxIndex;
            _heap[index].second = k;

            if (2*index < _heap.size() && _heap[2*index].second < k) 
            {
                maxIndex = 2*index;
            }
            else
            {
                maxIndex = index;
            }
            if ((2*index + 1) < _heap.size() && _heap[2*index+1].second < _heap[maxIndex].second) 
            {
                maxIndex = 2*index + 1;   
            }
        }
        _heap[index].first = descendingElement;
        _elements[descendingElement].index(index);
    }

};



#endif //BINARYHEAP_H


//	/*
//	 * This binary heap is based on priority queue in stl
//	 *
//	 */
//#ifndef BINARYHEAP_H
//#define BINARYHEAP_H
//
//#include <vector>
//
//	/**
//	 * In case of the (FIFO)BinaryHeap, we distinguish between an ExternalKey,
//	 * which is provided by the application (which uses the heap) and an internal Key,
//	 * which is used inside the heap. The ExternalKey is always a part of the internal Key.
//	 * A "key extractor" extracts the ExternalKey out of the internal Key.
//	 * 
//	 * In the simple case of a BinaryHeap without the FIFO property, the internal Key just
//	 * corresponds to the ExternalKey. Thus, the SimpleKeyExtractor just returns the
//	 * given internal Key as the ExternalKey.
//	 * @param ExternalKey the type of the external key = the type of the internal key
//	 */
//template < typename ExternalKey >
//class SimpleKeyExtractor
//{
//public:
//    /** Returns the given internal Key as the ExternalKey. */
//    static ExternalKey key(ExternalKey k) {return k;}
//};
//
//
//	/**
//	 * Represents an element of a BinaryHeap.
//	 * @param ExternalKey the type of the external key (see SimpleKeyExtractor)
//	 * @param MetaExtKey class that provides data about the external key, e.g. the max value
//	 * @param Data the type of the application-specific data that is associated with this element
//	 * @param Count the type of the counter that counts the number of heap elements / 
//	 *              the number of heap operations
//	 */
//template < typename ExternalKey,
//           typename MetaExtKey,
//           typename Data,
//           typename Count >
//class BinaryHeapElement
//{
//public:
//    /** Constructs a DUMMY element. */
//    BinaryHeapElement() : _key(MetaExtKey::MAX_VALUE), _index(0) {}
//            
//    /**
//     * Constructor.
//     * @param key key of the new element
//     * @param index index of the new element within the heap
//     */
//    BinaryHeapElement(ExternalKey key, Count index) : _key(key), _index(index) {}
//    
//    
//    /** Returns the external key of this element. */
//    ExternalKey key() const {return _key;}
//    
//    /** Sets the external key of this element. */
//    void key(ExternalKey newKey) {_key = newKey;}
//
//    bool isDummy() const {return key() == MetaExtKey::MAX_VALUE;}
//    
//    
//    /** Returns a reference to the application-specific data object of this element. */
//    Data& data() {return _data;}
//    
//    /** Returns a reference to the application-specific data object of this element. */
//    const Data& data() const {return _data;}
//    
//    
//    /** Marks that this element has been deleted from the heap. */
//    void markDeleted() {index(0);}
//    
//    /** Returns true iff this element has been deleted from the heap. */
//    bool hasBeenDeleted() const {return (_index==0);}
//
//    
//    /** Returns this element's index within the heap. */
//    Count index() const {return _index;}
//    
//    /** Sets this element's index within the heap. */
//    void index(Count newIndex) {_index = newIndex;}
//    
//private:
//    /** The external key of this element. */
//    ExternalKey _key;    
//    /** The application-specific data that is associated with this element. */
//    Data _data;
//    /**
//     * This element's index within the heap.
//     * 0 is used to indicate that this element has been deleted from the heap.
//     */
//    Count _index;
//};
//
//
//struct IndexKey 
//{
//	NodeID _index;
//	EdgeWeight _distance;
//	NodeID _id;
//
//	IndexKey(NodeID index = 0, EdgeWeight weight = 0, NodeID id = 0):
//	_index(index), _distance(weight), _id(id)
//	{
//
//	}
//	bool operator == (const IndexKey &a){
//		return (_id == a._id) && (_distance == a._distance);
//	}
//};
//
//struct IndexKeyComparision{
//	bool operator() (const IndexKey& left, const IndexKey& right)
//	{
//		if (left._distance == right._distance)
//			return left._index  > right._index;
//		else
//			return (left._distance > right._distance);
//	}
//};
//
//template<class T, class Compare>
//class PQV {
//	vector<T> v;
//	Compare comp;
//public:
//	// Don't need to call make_heap(); it's empty:
//	PQV(Compare cmp = Compare()) : comp(cmp) {}
//	void push(const T& x) {
//		// Put it at the end:
//		v.push_back(x);
//		// Re-adjust the heap:
//		push_heap(v.begin(), v.end(), comp);
//	}
//	void pop() {
//		// Move the top element to the last position:
//		pop_heap(v.begin(), v.end(), comp);
//		// Remove that element:
//		v.pop_back();
//	}
//	const T& top() { return v.front(); }
//	const T& back(){ return v.back(); }
//	void pop_back()
//	{
//		v.pop_back();
//	}
//	bool empty() const { return v.empty(); }
//	int size() const { return v.size(); }
//	bool update(T &t, EdgeWeight &newKey) 
//	{
//		std::vector<T>::iterator it = find(v.begin(), v.end(), t);
//		EdgeWeight oldKey = t._distance;
//		if (it != v.end()){
//			it->_distance = newKey;
//			if (oldKey > newKey)
//				make_heap(v.begin(), it + 1, comp);
//			else
//				make_heap(it, v.end(), comp);
//			return true;
//		}
//		else{
//			return false;
//		}
//	}
//	
//	bool contains(const T& x) const {
//		return (find(v.begin(), v.end(), x) == v.end());
//	}
//	typedef vector<T> TVec;
//	TVec vector() {
//		TVec r(v.begin(), v.end());
//		// It's already a heap
//		sort_heap(r.begin(), r.end(), comp);
//		// Put it into priority-queue order:
//		reverse(r.begin(), r.end());
//		return r;
//	}
//
//};
//	/**
//	 * Represents a binary heap.
//	 * @param ExternalKey the type of the external key (see SimpleKeyExtractor)
//	 * @param MetaExtKey class that provides data about the external key, e.g. the max value
//	 * @param Data the type of the application-specific data that is associated with this element
//	 * @param Count the type of the counter that counts the number of heap elements / 
//	 *              the number of heap operations
//	 * @param Key the type of the internal key (see SimpleKeyExtractor)
//	 * @param KeyExtractor the "key extractor" that is used (see SimpleKeyExtractor)
//	 */
//template < typename ExternalKey,
//           typename MetaExtKey,
//           typename Data,
//           typename Count,
//           typename Key = ExternalKey,
//           typename KeyExtractor = SimpleKeyExtractor<ExternalKey> >
//class BinaryHeap
//{
//public:
//    /** The type of the elements stored in this heap. */
//    typedef BinaryHeapElement<ExternalKey, MetaExtKey, Data, Count> PQElement;
//    
//    typedef Data PQData;
//    
//	typedef PQV <IndexKey, IndexKeyComparision> PQHeap;
//    /** Constructor. */
//    BinaryHeap() {
//		//reserve(1000);
//        insertDummy();
//        //_heap.push_back( IndexKey() ); // add a dummy element
//        //_heap[0].second = MetaExtKey::MIN_VALUE;
//    }
//
//    /** Returns the size of this heap (= the number of elements). */
//    Count size() const { 
//        //return _heap.size()-1; // subtract dummy element (index 0)
//		return _heap.size();
//    }
//    
//    void reserve( size_t n )
//    {
//        _elements.reserve(n+1);
//        //_heap.reserve(n+1);
//    }
//                            
//    /** Returns true iff the heap is empty. */
//    bool empty() const { return (size() == 0); }
//
//    /** 
//     * Returns the (external) key of the minimum element in this heap
//     * or the maximum value of the ExternalKey data type iff the heap is empty.
//     */
//    ExternalKey min() {   
//        //if (empty()) return MetaExtKey::MAX_VALUE;
//        //return KeyExtractor::key(_heap[1].second);
//		if (empty()) return MetaExtKey::MAX_VALUE;
//		return KeyExtractor::key(_heap.top()._distance);
//
//    }
//
//    /**
//     * Creates a new element with the given (internal) key and inserts it into the heap.
//     * @return the index of the new element.
//     */
//    Count insert(Key key, NodeID id) {
//        /*Count element = _elements.size();
//        Count index = _heap.size();
//        _elements.push_back( PQElement(KeyExtractor::key(key), index) );
//        _heap.push_back(IndexKey(element,key));
//        upheap(index);        
//        return element;*/
//
//		Count element = _elements.size();
//		_elements.push_back( PQElement(KeyExtractor::key(key), 0));
//		_heap.push(IndexKey(element, key, id));
//		return element;
//    }
//    
//    /** Replaces a dummy element and inserts the actual element into the heap. */
//    void insert(Key key, Count elementIndex, NodeID id) {          
//        /*PQElement& element = _elements[elementIndex];
//        Count index = _heap.size();
//        element.key(KeyExtractor::key(key));        
//        element.index(index);
//        _heap.push_back(IndexKey(elementIndex,key));
//        upheap(index);*/
//
//		PQElement& element = _elements[elementIndex];
//		element.key(KeyExtractor::key(key));        
//		_heap.push(IndexKey(elementIndex, key, id));
//    }
//    
//    /** Adds a dummy element, which is NOT inserted into the heap, yet. */
//    void insertDummy() {        
//        _elements.push_back( PQElement() );
//    }
//
//    /**
//     * Initialises the _elements vector so that it contains the given number of dummy elements.
//     * This is useful when a search is performed after another search has not cleaned up.
//     * In particular, this is used when the main phase of the landmark query is executed
//     * so that the shortest path tree of the initial phase can be kept in order to be able
//     * to compose the complete path at the end.
//     */
//    void initDummies(Count number) {
//        assert( empty() );
//        _elements.resize(number);
//    }
//    
//    bool isDummy(Count element) const {
//        return (_elements[element].isDummy());
//    }
//    
//    /** Returns a reference to the elements. */
//    vector<PQElement>& elements() {return _elements;}
//    
//    /** Returns a reference to the elements. */
//    const vector<PQElement>& elements() const {return _elements;}
//
//    /**
//     * Returns the index of the minimum element in this heap
//     * or 0 iff the heap is empty.
//     */
//    Count minElement() const {
//       /* if (empty()) return 0;
//        return _heap[1].first;*/
//		if (empty()) return 0;
//		return _heap.top()._index;
//    }
//
//    /**
//     * Deletes the minimum element in this heap.
//     * Precondition: The heap is not empty.
//     * @return the index of the deleted element
//     */
//    Count deleteMin() {
//        //assert( ! empty() );
//        //const Count element = _heap[1].first;
//        //_heap[1] = _heap[_heap.size()-1];
//        //Count index = 1;
//        //const Count droppingElement = _heap[index].first;
//        //_heap.pop_back();
//        //if (size() > 1) {
//        //    // downheap:
//        //    // Move the element at the top downwards
//        //    // until the heap condition is restored.            
//        //    
//        //    Key k = _heap[index].second;
//        //    Count nextIndex = 2*index;
//        //    while (nextIndex < _heap.size()) {
//        //        nextIndex += 
//        //            ((nextIndex+1 < _heap.size()) &&
//        //            (_heap[nextIndex].second > _heap[nextIndex+1].second));
//        //        
//        //        assert( _elements[_heap[nextIndex].first].key() == KeyExtractor::key(_heap[nextIndex].second) );
//        //        
//        //        if (k <= _heap[nextIndex].second) break;
//        //        
//        //        _heap[index] = _heap[nextIndex];
//        //        _elements[_heap[nextIndex].first].index(index);
//        //        index = nextIndex;
//        //        nextIndex *= 2;
//        //    }
//        //    _heap[index].first = droppingElement;
//        //    _heap[index].second = k;
//        //    
//        //    // end of downheap
//        //}
//        //_elements[droppingElement].index(index);
//        //
//        //_elements[element].markDeleted();
//        //COUNTING( counter.incDouble(COUNT_DEL_MIN) );
//        //return element;
//		
//		const Count element = _heap.top()._index;
//		_heap.pop();
//		_elements[element].markDeleted();
//		return element;
//
//    }
//
//    /**
//     * Deletes an arbitrary element in this heap.
//     * Precondition: The heap is not empty.
//     * @return the index of the deleted element
//     */
//    Count deleteArbitrary() {
//        /*assert( ! empty() );         
//        Count element = _heap.back().first;
//        _heap.pop_back();
//        _elements[element].markDeleted();
//        return element;*/
//		Count element = _heap.back()._index;
//		_heap.pop_back();
//		_elements[element].markDeleted();
//		return element;
//		//return 0;
//		
//    }
//
//    /**
//     * Increases the key of the given element:
//     * sets the key to the given value.
//     */
//    void increaseKey(Count element, Key newKey) {
//       /* Count index = _elements[element].index();
//        assert( index < _heap.size() );
//        assert( _heap[index].first == element );
//        
//        _elements[element].key(KeyExtractor::key(newKey));
//        _heap[index].second = newKey;
//        downheap(index);
//        COUNTING( counter.incDouble(COUNT_INCR_KEY) );*/
//		Count index = _elements[element].index();
//
//		Key oldKey = _elements[element].key();
//		_elements[element].key(KeyExtractor::key(newKey));
//		NodeID id = _elements[element].data().nodeID();
//		IndexKey oldIndexKey(index, oldKey, id);
//		_heap.update(oldIndexKey, newKey);
//    }
//    
//    /**
//     * Decreases the key of the given element:
//     * sets the key to the given value.
//     */
//    void decreaseKey(Count element, Key newKey) {
//        /*Count index = _elements[element].index();
//        assert( index < _heap.size() );
//        assert( _heap[index].first == element );
//        
//        _elements[element].key(KeyExtractor::key(newKey));
//        _heap[index].second = newKey;
//        upheap(index);
//        COUNTING( counter.incDouble(COUNT_DECR_KEY) );*/
//
//		Count index = _elements[element].index();
//		assert( index < _heap.size() );
//		assert( _heap[index].first == element );
//
//		Key oldKey = _elements[element].key();
//		_elements[element].key(KeyExtractor::key(newKey));
//		NodeID id = _elements[element].data().nodeID();
//		IndexKey oldIndexKey(index, oldKey, id);
//		_heap.update(oldIndexKey, newKey);
//		COUNTING( counter.incDouble(COUNT_DECR_KEY) );
//    }
//    
//    /**
//     * Updates the key of the given element:
//     * sets the key to the given value.
//     * The key can be larger or smaller than the previous key.
//     */
//    void updateKey(Count element, Key newKey) {
//        /*Count index = _elements[element].index();
//        assert( index < _heap.size() );
//        assert( _heap[index].first == element );
//        Key oldKey = _heap[index].second;
//
//        _elements[element].key(KeyExtractor::key(newKey));
//        _heap[index].second = newKey;
//        if (newKey < oldKey)
//        {        
//            upheap(index);
//            COUNTING( counter.incDouble(COUNT_DECR_KEY) );
//        }
//        else if (newKey > oldKey)
//        {
//            downheap(index);
//            COUNTING( counter.incDouble(COUNT_INCR_KEY) );
//        }*/
//		Count index = _elements[element].index();
//
//		Key oldKey = _elements[element].key();
//		_elements[element].key(KeyExtractor::key(newKey));
//		NodeID id = _elements[element].data().nodeID();
//		IndexKey oldIndexKey(index, oldKey, id);
//		_heap.update(oldIndexKey, newKey);
//		
//    }
//	NodeID elementsSize()
//	{
//		return _elements.size();
//	}
//
//    void clear() {
//        _elements.clear();
//		//_heap.clear();
//        insertDummy();
//    }
//    
//	~BinaryHeap()
//	{
//		//_elements.clear();
//		//_heap.clear();
//		vector <PQElement>().swap(_elements);
//		//vector <IndexKey>().swap(_heap);
//	}
//
//    
//    unsigned int totalHeapSize() const
//	{
//		unsigned int pairSize = sizeof(Count) + sizeof(Key) ;
//		unsigned int  eSize = sizeof(PQElement) * _elements.size();
//		unsigned int hSize = pairSize * _heap.size();
//		return eSize + hSize;
//	}
//private:
//    //typedef pair<Count, Key> IndexKey;
//    
//    /**
//     * The elements of this heap.
//     * The order corresponds to the order of insertion and
//     * is not changed by any heap operation. That implies
//     * that an index of this vector can be used as a pointer
//     * to the corresponding element.
//     * It is possible that this vector contains holes (i.e.
//     * dummy elements), namely for each element that has been
//     * inserted into the other pqueue (for the other search direction)
//     * but not to this one.
//     * The first element (index 0) is a dummy element so that the index 0
//     * can be used to mark elements that have not been inserted into the heap.
//     */
//    vector<PQElement> _elements;
//    
//    /** 
//     * "Pointers" (first) to the elements of this heap in the right heap order.
//     * In addition (second), the internal key of the corresponding element.
//     * The first element (index 0) is a dummy element so that the index 0
//     * can be used to mark elements that have been deleted from the heap.
//     */
//    PQHeap _heap;
//    
//
//};
//
//
//
//#endif // BINARYHEAP_H
