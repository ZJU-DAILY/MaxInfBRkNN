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

#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <functional>
#include <set>
#include <map>
#include <unordered_map>
//#include "../../hr_time.h"
#include "../../types.h"

using namespace std;

// #include <hash_map>
// #include <hash_set>

/** 
 * Type of a node ID.
 * Also used for everything that is bounded by the number of nodes,
 * e.g. the number of components.
 */
//typedef int NodeID;
// typedef unsigned int NodeID; // Exported to types.h

/** Type of an edge ID. */
typedef unsigned int OtherEdgeID;

/** Type of a x or y coordinate. */
// typedef double CoordinateType; // Exported to types.h

/** Type of a source/target node pair. */
typedef pair<NodeID, NodeID> stPair;
typedef vector<stPair> stPairs;

/** 
 * Special node ID which is normally not used
 * so that it can be used as flag.
 */

static const NodeID SPECIAL_NODEID = __INT_MAX__ * 2U + 1;

static const OtherEdgeID MAX_EDGEID = __INT_MAX__ * 2U;

/** Type of a level ID. */
typedef int LevelID;

#include "edge.h"
#include "path.h"


/**
 * Represents a matrix that can store a result of
 * a many-to-many computation.
 */
template <typename value_type, typename size_type = unsigned int>
class Matrix
{
    friend ostream& operator<<( ostream& os, const Matrix& matrix ) {
		os << matrix.noOfRows() << " * " << matrix.noOfCols() << endl;
        for (size_type i = 0; i < matrix.noOfEntries(); i++) {
            if (i > 0) {
                if (i % matrix.noOfCols() == 0) os << endl; else os << " ";
            }
            os << matrix._data[i];
        }
		os << endl;
        return os;
    }
    
public:
    bool operator== (const Matrix<value_type, size_type>& matrix) const {
        if (noOfRows() != matrix.noOfRows()) return false;
        if (noOfCols() != matrix.noOfCols()) return false;
        return (_data == matrix._data);
    }

    /** Constructor. */
	Matrix()
	{

	}
    Matrix(const size_type r, const size_type c) : _rows(r), _cols(c), _data(r*c) {}

	void serialize(ostream& out) {	
		//VERBOSE( cout << "serialize matrix: " << _rows << " * " << _cols << endl );
		PrimitiveSerializer <size_type>::serialize(out, _rows);
		PrimitiveSerializer <size_type>::serialize(out, _cols);
		VectorSerializer< value_type, NodeID>::serialize(out, _data);
		//write all the matrix
		//VERBOSE( cout << "done." << endl );
	}
	/** Deserializes the graph from the given stream. */
	void deserialize(istream& in) {
		//VERBOSE( cout << "deserialize matrix: " << endl );
		PrimitiveSerializer <size_type>::deserialize(in, _rows);
		PrimitiveSerializer <size_type>::deserialize(in, _cols);
		// nodes
		VectorSerializer< value_type, NodeID >::deserialize(in, _data);
		//VERBOSE( cout << "done. " << _rows << " * " << _cols << endl );
	} 


	void setRowAndCol(size_type r, size_type c)
	{
		_rows = r;
		_cols = c;
	}


    /** Inits all entries to a given value. */
    void init(const value_type& v) {
		_data.clear();
		_data.insert( _data.end(), _rows * _cols, v );
    }

    size_type noOfRows() const {return _rows;}
    size_type noOfCols() const {return _cols;}
    size_type noOfEntries() const {return _data.size();}

    value_type value(const size_type r, const size_type c) const 
	{
		return _data[r * _cols + c];
	}

    /**
    * Sets a specified entry to the given value.
   
    */
    void set(const size_type r, const size_type c, const value_type& v) {
        _data[r * _cols + c] = v;
    }

    /**
    * Sets a specified entry to the given value
   
    */
    void improve(const size_type r, const size_type c, const value_type& v) {
        improve(r * _cols + c, v);
    }

	void improveNojudge(const size_type r, const size_type c, const value_type& v) {
		_data[r * _cols + c] = v;
		//improve(index(r, c), v);
	}

    /**
    * Sets a specified entry to the given value
   
    */
    void improve(const size_type i, const value_type& v) {
        assert( i < noOfEntries() );
        if (v < _data[i]) _data[i] = v;
    }

    /**
     * Returns the index of a specified entry.
 
     */
    size_type index(const size_type r, const size_type c) const {
        assert( r < noOfRows() );
        assert( c < noOfCols() );
        assert( r * _cols + c < noOfEntries() );
        return r * _cols + c;
    }

	size_type totalItems()
	{
		return _data.size();
	}

	void clear(){
		_data.clear();
	}
    


private:
    /** The number of rows. */
   size_type _rows;

    /** The number of columns. */
    size_type _cols;

    /**
     * All entries.
     */
    vector<value_type> _data;
};


// template <>
// struct hash_compare<pair<NodeID, NodeID> > {
// 	// the mean number of elements per bucket
// 	static const size_t bucket_size = 4;
// 	// the minimum number of buckets
// 	static const size_t min_buckets = 8;
// 	// hash function
// 	size_t operator()(const pair<NodeID, NodeID>& key) const 
// 	{
// 		//return key.first * 31 + key.second;
// 		return hash_value(key.first) + hash_value(key.second);
// 	}
// 	// compare function
// 	bool operator()(const pair<NodeID, NodeID>& lhs, const pair<NodeID, NodeID>& rhs) const 
// 	{
// 		return lhs < rhs;
// 	}
// 
// };

// typedef unordered_map <pair<NodeID, NodeID>, EdgeWeight, Hash_Key, Equal_Key> LocalTableType;
// typedef hash_map <pair<NodeID, NodeID>, EdgeWeight> LocalTableType;
#endif // GRAPH_H
