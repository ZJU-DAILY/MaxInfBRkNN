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

#ifndef SEARCHSPACES_H
#define SEARCHSPACES_H


/**
 * Represents forward and backward search spaces.
 * Used by the preprocessing of transit-node routing
 * and for many-to-many computations.
 * In both cases, we first collect backward search spaces
 * for several targets; then, we rearrange the search spaces;
 * and finally, we repeatedly obtain a forward search space that is
 * 'intersected' with the backward search spaces.
 * This class supports all the required methods.
 * There are some tricks to save memory, namely a difference encoding
 * and distinguishing between entries that fit into shorts and entries
 * that fit into integers. Function template parameters allow to switch
 * these tricks on or off.
 */
class SearchSpaces
{
private:
    
    /**
    * While a search from some node x is performed, at each encountered node u
    * we add a SearchSpaceEntry containing the origin x and the distance from x to u.
    * Instead of storing x, we store (x - u), which often needs less bits.
    */
    template <typename SSENodeID_type, typename SSEEdgeWeight_type>
    class SearchSpaceEntry
    {
    public:
        typedef SSENodeID_type SSENodeID;
        typedef SSEEdgeWeight_type SSEEdgeWeight;

        /** Default constructor. */
        SearchSpaceEntry() {}

        /**
       * Constructor.
       * @param u the ID of the encountered node; ignored;
       *          only to be compatible with CompleteSearchSpaceEntry
       * @param o the origin x of the search minus u
       * @param d the distance from the origin to the node u
       */
        SearchSpaceEntry(NodeID u, SSENodeID o, SSEEdgeWeight d) : _origin(o), _dist(d) {}

        /**
       * Returns the origin of the search.
       * Only valid if the difference encoding has been used!
       * @param u the ID of the node that this entry belongs to
       */
        NodeID origin(NodeID u) const {return (u + _origin);}

        /**
       * Returns the origin of the search.
       * Only valid if NO difference encoding has been used!
       */
        NodeID origin() const {return _origin;}

        EdgeWeight dist() const {return _dist;}
    
        void serialize(ostream& out) const {
            writePrimitive(out, _origin);
            writePrimitive(out, _dist);
        }
    
        void deserialize(istream& in) const {
            readPrimitive(in, _origin);
            readPrimitive(in, _dist);
        }
        
    private:
        SSENodeID _origin;
        SSEEdgeWeight _dist;
    };

    /**
    * Extending SearchSpaceEntry, stores explictely the ID of the node
    * that the entry belongs to.
    */
    template <typename SSENodeID_type, typename SSEEdgeWeight_type>
    class CompleteSearchSpaceEntry : public SearchSpaceEntry<SSENodeID_type, SSEEdgeWeight_type>
    {
    private:
        typedef SearchSpaceEntry<SSENodeID_type, SSEEdgeWeight_type> Super;
        
    public:
        typedef SSENodeID_type SSENodeID;
        typedef SSEEdgeWeight_type SSEEdgeWeight;
    
        /** Used to group a vector of entries by node ID. */
        bool operator< (const CompleteSearchSpaceEntry& e) const {
            return (nodeID() < e.nodeID());
        }

        /** Default constructor. */
        CompleteSearchSpaceEntry() {}

        /**
       * Constructor.
       * @param u the ID of the node that the entry belongs to
       * @param o the origin x of the search minus u
       * @param d the distance from the origin to the node u
       */
        CompleteSearchSpaceEntry(NodeID u, SSENodeID o, SSEEdgeWeight d)
        : Super(u, o, d), _nodeID(u) {}
    
        NodeID nodeID() const {return _nodeID;}
        
    private:
        NodeID _nodeID;
    };

    /**
    * After a search from all relevant nodes has been performed, all search
    * space entries that have been obtained are grouped and stored in an
    * IndexedSearchSpace, which allows to directly access all entries for
    * a specific node.
    */
    template <typename Entry>
    class IndexedSearchSpace
    {
    public:
        typedef typename vector<Entry>::const_iterator const_iterator;

        /** Default constructor. */
        IndexedSearchSpace() {}

        /**
       * Copy constructor.
       * Used to copy a dynamic to a static indexed search space.
       */
        template <typename EntryIn>
        IndexedSearchSpace(const IndexedSearchSpace<EntryIn>& iss) : _index(iss.index()) {
            _searchSp.resize(iss.searchSp().size());
            copy(iss.searchSp().begin(), iss.searchSp().end(), _searchSp.begin());
        }
        
        void serialize(ostream& out) const {
            VectorSerializer< Entry, NodeID, ComplexSerializer<Entry> >::serialize(out, _searchSp);
            VectorSerializer< NodeID, NodeID >::serialize(out, _index);
        }
    
        void deserialize(istream& in) {
            VectorSerializer< Entry, NodeID, ComplexSerializer<Entry> >::deserialize(in, _searchSp);
            VectorSerializer< NodeID, NodeID >::deserialize(in, _index);
        }
    
        void freeMemory() {
            vector<Entry>().swap(_searchSp);
            vector<NodeID>().swap(_index);
        }

        /** Returns a const_iterator to the first entry for the given node. */
        const_iterator begin(NodeID u) const {
            assert( u < _index.size() );
            return _searchSp.begin() + _index[u];
        }

        /**
       * Returns a const_iterator to the entry after the last entry of the
       * given node.
       */
        const_iterator end(NodeID u) const {
            return begin(u+1);
        }
    
        NodeID size() const {
            return _searchSp.size();
        }

        const vector<Entry>& searchSp() const {return _searchSp;}

        const vector<NodeID>& index() const {return _index;}
        
    protected:
        /** Contains all search space entries. */
        vector<Entry> _searchSp;

        /** Contains the indices of the first entries in the _searchSp vector. */
        vector<NodeID> _index;
    };

    /**
    * During the search process all search space entries are just added
    * in an unordered way to a DynamicIndexedSearchSpace. After a sequence
    * of 'add's, we have to call 'sort', which transforms this data structure
    * such that the methods of the super class IndexedSearchSpace can be used.
    */
    template <typename Entry>
    class DynamicIndexedSearchSpace : public IndexedSearchSpace<Entry>
    {
    private:
        typedef IndexedSearchSpace<Entry> Super;
        typedef typename Entry::SSENodeID SSENodeID;
        typedef typename Entry::SSEEdgeWeight SSEEdgeWeight;
        using Super::_searchSp;
        using Super::_index;
        
    public:
        /**
       * Adds a new entry.
       * @param u the ID of the node that the entry belongs to
       * @param o the origin x of the search minus u
       * @param d the distance from the origin to the node u
       */
        void add(NodeID u, SSENodeID o, SSEEdgeWeight d) {
            // our own growth strategy for the _searchSp vector:
            // if it is already very big, we prefer a linear growth to an exponential growth
            // advantage: better chances that we don't exceed the available memory
            const unsigned int memoryThreshold = 100000000;
            if ((_searchSp.size() == _searchSp.capacity()) && (_searchSp.size() >= memoryThreshold)) {
                _searchSp.reserve(_searchSp.size() + memoryThreshold);
            }

            // add the entry
            _searchSp.push_back(Entry(u, o, d));
        }

        /**
       * AFTER the search space has been filled (i.e., NO further 'add' operations
       * will follow!), we can sort the search space (group it by node ID).
       * After that, we may use the access methods of the super class.
       */
        void sort(NodeID n) {
            // group by node ID
            std::sort( _searchSp.begin(), _searchSp.end() );

            // build index
            _index.resize(n + 1);
            NodeID oldID = 0;
            _index[0] = 0;
            for (NodeID i = 0; i < _searchSp.size(); i++) {
                NodeID currentID = _searchSp[i].nodeID();
                if (currentID != oldID) {
                    for (NodeID j = oldID + 1; j <= currentID; j++) {
                        _index[j] = i;
                    }
                    oldID = currentID;
                }
            }
            for (NodeID j = oldID + 1; j <= n; j++) {
                _index[j] = _searchSp.size();
            }

            // add dummy element (where the last end index points to)
            _searchSp.push_back(Entry());
        }
    };


public:
    // note: we have two instances of each type of search space, one based on 'int's and
    //       one based on 'short's. Whenever possible we store an entry in the 'short' version;
    //       if the values of an entry are too big, we store it in the 'int' version.
    typedef DynamicIndexedSearchSpace< CompleteSearchSpaceEntry<int, EdgeWeight> > DISSInt;
    typedef DynamicIndexedSearchSpace< CompleteSearchSpaceEntry<short, short> > DISSShort;
    typedef IndexedSearchSpace< SearchSpaceEntry<int, EdgeWeight> > ISSInt;
    typedef IndexedSearchSpace< SearchSpaceEntry<short, short> > ISSShort;
    
    void setCurrentNode(const NodeID cN) {_currentNode = cN;}
    
    /** Sorts the backward search spaces and swaps them out to hard disk. */
    void swapOutSearchSpacesBW(const string filename, const NodeID noOfNodes) {
        VERBOSE( cout << "  search space sizes: " << _searchSpacesBwDynInt.size()
                      << " + " << _searchSpacesBwDynShort.size() << endl );
        
        ofstream out(filename.c_str());
        _searchSpacesBwDynInt.sort(noOfNodes);
        _searchSpacesBwDynInt.serialize(out);
        _searchSpacesBwDynInt.freeMemory();
        
        _searchSpacesBwDynShort.sort(noOfNodes);
        _searchSpacesBwDynShort.serialize(out);
        _searchSpacesBwDynShort.freeMemory();
    }

    /**
    * Reads the backward search spaces from hard disk.
    * A sequence of a swap out and a swap in can be used to transform
    * a dynamic search space to a (static) search space, which needs less
    * memory and which contains its elements in an ordered fashion.
    */
    void swapInSearchSpacesBW(const string filename) {
        ifstream in(filename.c_str());
        _searchSpacesBwInt.deserialize(in);
        _searchSpacesBwShort.deserialize(in);

        VERBOSE( cout << "  search space sizes: " << _searchSpacesBwInt.size()
                      << " + " << _searchSpacesBwShort.size() << endl );
    }

    
    
    // *********************************************************
    // ****************** called from Dijkstra *****************
    // *********************************************************


    /**
    * Adds an entry to the backward search space of the current node.
    * @param u the encountered node
    * @param dist the backward distance from '_currentNode' to u
    */
    template <bool useDiffEncoding, bool allowShorts>
    void addToSearchSpaceBW(NodeID u, EdgeWeight dist) {
        // compute '_currentNode' minus u since this value often
        // can be stored using less bits
        int diff = _currentNode;
        if (useDiffEncoding) diff -= u;
        
        if (allowShorts && (dist < 0x10000) && (diff >= -0x8000) && (diff < 0x8000)) {
            // entry fits in 'short's
            _searchSpacesBwDynShort.add(u, diff, dist);
        }
        else {
            // entry only fits in 'int's
            _searchSpacesBwDynInt.add(u, diff, dist);
        }
    }

    /**
    * Adds an entry to the forward search space of the current node s.
    * @param u the encountered node
    * @param dist the distance from s to u
    */
    void addToSearchSpaceFW(NodeID u, EdgeWeight dist) {
        _searchSpaceFW.push_back(Edge(u, dist, false, false, false));
        // note: the flags (false, false, false) have NO meaning
    }


    
protected:
    // Backward search spaces (of all relevant nodes).
    // We have a dynamic and a static variant.
    // The dynamic variant is used while the search spaces are determined.
    // Then, the dynamic variant is tranformed to the static variant,
    // which is used to access the data.
    // We have two instances of each type of search space, one based on 'int's and
    // one based on 'short's. Whenever possible we store an entry in the 'short' version;
    // if the values of an entry are too big, we store it in the 'int' version.
    DISSInt _searchSpacesBwDynInt;
    DISSShort _searchSpacesBwDynShort;
    ISSInt _searchSpacesBwInt;
    ISSShort _searchSpacesBwShort;

    /** Forward search space (of the current node). */
    vector<Edge> _searchSpaceFW;

    NodeID _currentNode;
};

#endif // SEARCHSPACES_H
