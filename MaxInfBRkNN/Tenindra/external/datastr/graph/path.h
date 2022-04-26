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

#ifndef PATH_H
#define PATH_H

#include <list>

using namespace std;

/**
 * Represents a path in a graph.
 * Two types of path objects can be dealt with.
 * Type A: Normally, a path is a list of x+1 nodes and x edge weights
 * and shortcut flags.
 * Type B: In addition to the lists of Type A, a list of x edge IDs is managed.
 * It is important to ensure that both types are not mixed up.
 */
class Path
{
    /** Outputs data about the path for debugging purposes. */ 
    friend ostream& operator<<( ostream& os, const Path& p ) {
        if (p._notConnected) {
            os << "NOT CONNECTED" << endl;
            return os;
        }
        if (p.empty()) {
            os << "EMPTY" << endl;
            return os;
        }
        assert( p.checkInvariant() );
        os << setprecision(20);
        os << p.length() << ": " << p.node(0);
        for (OtherEdgeID e = 0; e < p.noOfEdges(); e++)
            os << " --(" << p._weights[e] << ")-> " << p.node(e+1); 
        os << endl;
        return os;
    }
    
public:
    /** Returns true iff the given path is identical to this path. */
    bool operator== (const Path& p) const {
        if (_notConnected != p._notConnected) return false;
        if (_notConnected) return true;
        return ((_nodes == p._nodes) && (_weights == p._weights));
    }

    /** Compares two paths by length (=sum of weights). */
    bool operator< (const Path& p) const {
        return (length() < p.length());
    }

    /** 
     * Constructor.
     * Creates an empty path.
     */
    Path() {clear();}
    
    /**
     * Constructor.
     * Creates a path with a given start node.
     */
    Path(NodeID id) {
        clear();
        addFirstNode(id);
    }

    /** Clears this path (so that it becomes an empty path). */
    void clear() {
        _nodes.clear();
        _weights.clear();
        _edges.clear();
        _shortcuts.clear();
        _length = 0;
        _notConnected = false;
    }

    /** Adds a start node to this (empty) path. */
    void addFirstNode(NodeID id) {
        _nodes.push_back(id);
        assert( checkInvariant() );
    }

    /** 
     * Adds an edge to this path. 
     * Note: the source of the edge is the formerly last node in the path.
     * @param id the target of the edge
     * @param w the weight of the edge
     * @param shc the shortcut flag of the edge
     */
    void add(NodeID id, EdgeWeight w, bool shc) {
        _nodes.push_back(id);
        _weights.push_back(w);
        _shortcuts.push_back(shc);
        _length += w;
        assert( checkInvariant() );
    }

    /** 
     * Adds an edge to this path.
     * Note: the source of the edge is the formerly last node in the path.
     * @param nodeID the target of the edge
     * @param w the weight of the edge
     * @param edgeID the ID of the edge
     * @param shc the shortcut flag of the edge
     */
    void add(NodeID nodeID, EdgeWeight w, OtherEdgeID edgeID, bool shc) {
        add(nodeID, w, shc);
        _edges.push_back(edgeID);
    }

    /**
     * Adds a path to this path.
     * Note: the first node is NOT added; we assume that it is
     * identical to the formerly last node of this path.
     */
    void add(const Path& path) {
        _nodes.insert( _nodes.end(), path._nodes.begin() + 1, path._nodes.end() );
        _weights.insert( _weights.end(), path._weights.begin(), path._weights.end() );
        _edges.insert( _edges.end(), path._edges.begin(), path._edges.end() );
        _shortcuts.insert( _shortcuts.end(), path._shortcuts.begin(), path._shortcuts.end() );
        _length += path.length();
        assert( checkInvariant() );
    }

    /**
     * Adds sequences of node ids and edge weights to this path.
     * Note: The shortcut flags are filled with dummy entries.
     * @param beginNodes begin of the node sequence
     * @param endNodes end of the node sequence
     * @param beginWeights begin of the weights sequence
     * @param endWeights end of the weights sequence
     * @param w sum of all weights in the weights sequence
     */
    void add(vector<NodeID>::const_iterator beginNodes, vector<NodeID>::const_iterator endNodes,
             vector<EdgeWeight>::const_iterator beginWeights, vector<EdgeWeight>::const_iterator endWeights,
             EdgeWeight w) {
        _nodes.insert( _nodes.end(), beginNodes, endNodes );
        _weights.insert( _weights.end(), beginWeights, endWeights );
        _shortcuts.resize( _weights.size() );
        _length += w;
        assert( checkInvariant() );
    }

    /**
     * Reverses this path. (The first node becomes the last one and vice versa.)
     * Note that this might lead to unwanted results if edge IDs are given since
     * only the order of the edge IDs is reversed; the IDs of the reverse edges
     * are not determined.
     */
    void reverse() {
        std::reverse( _nodes.begin(), _nodes.end() );
        std::reverse( _weights.begin(), _weights.end() );
        std::reverse( _edges.begin(), _edges.end() );
        std::reverse( _shortcuts.begin(), _shortcuts.end() );
    }

    /** 
     * When this path is used to store a path between two nodes s and t,
     * this flag can be set to indicate that s and t aren't connected.
     */
    void setNotConnected() {_notConnected = true;}

    /** Returns the number of nodes that belong to this path. */
    NodeID noOfNodes() const {return _nodes.size();}
    
    /** Returns the number of edges that belong to this path. */
    OtherEdgeID noOfEdges() const {
        assert( _weights.size() == _shortcuts.size() );
        return _weights.size();
    }

    /** Returns true iff this path is empty. */
    bool empty() const {return noOfNodes() == 0;}
    
    /** Returns the length (sum of weights) of this path. */
    EdgeWeight length() const {
        // special case: length "infinity" means "not connected"
        if (_notConnected) return Weight::MAX_VALUE;
 
        return _length;
    }

    /** Returns the node with the given index within this path. */
    NodeID node(NodeID index) const {
        assert( index < _nodes.size() );
        return _nodes[index];
    }

    NodeID firstNode() const {
        assert( ! empty() );
        return node(0);
    }

    NodeID lastNode() const {
        assert( ! empty() );
        return node(noOfNodes() - 1);
    }

    EdgeWeight weight(OtherEdgeID index) const {
        assert( index < _weights.size() );
        return _weights[index];
    }

    OtherEdgeID edge(OtherEdgeID index) const {
        assert( _edges.size() == _weights.size() );
        assert( index < _edges.size() );
        return _edges[index];
    }

    void setEdge(OtherEdgeID index, OtherEdgeID newID) {
	assert( index < _edges.size() );
	_edges[index] = newID;
    }
    
    bool isShortcut(OtherEdgeID index) const {
        assert( index < _shortcuts.size() );
        return _shortcuts[index];
    }

    bool isNotConnected() const {return _notConnected;}
    

    /**
     * Computes the differences between this path and a given path.
     * @param p the path that this path is compared to
     * @param diff1 a subpath of this path that begins just before the first 
     *              occurence of a difference and ends just after the last
     * @param diff2 a subpath of p that begins just before the first 
     *              occurence of a difference and ends just after the last
     */
    void diff(const Path& p, Path& diff1, Path& diff2) const {
        diff1.clear();
        diff2.clear();

        if (_notConnected != p._notConnected) {
            // completely different
            diff1 = (*this);
            diff2 = p;
            return;
        }

        if (*this == p) return; // completely equal
        
        OtherEdgeID i; // moves in both paths forwards
        OtherEdgeID j; // moves in this path backwards
        OtherEdgeID k; // moves in the path p backwards
        // look for the first difference:
        // move forwards until a difference occurs
        for (i=0; 
               (i<noOfEdges()) && (i<p.noOfEdges()) && 
               (node(i)==p.node(i)) && (_weights[i]==p._weights[i]);
             i++) ;
        // look for the last difference:
        // move backwards until a difference occurs
        for (j = noOfEdges(), k = p.noOfEdges();
             (j > 0) && (k > 0) && (node(j) == p.node(k)) && (_weights[j-1] == p._weights[k-1]);
             j--, k--) ;

        // if applicable, take one step back
        // because the output should begin just BEFORE the first difference
        // and end just AFTER the last difference
        if (i > 0) i--;
        if (j < noOfEdges()) j++;
        if (k < p.noOfEdges()) k++;
        
        subpath( diff1, i, j);   // copy the appropriate subpath of this path to diff1
        p.subpath( diff2, i, k); // copy the appropriate subpath of p to diff2
    }

    /**
     * Exports all edges of this path to the given stream using the given color
     * (in order to be able to draw the path).
     */
    void exportEdges(ostream& out, int color) {
        for (OtherEdgeID e = 0; e < noOfEdges(); e++)
            out << node(e) << " " << node(e+1) << " " << _weights[e] << " " << color << endl;
    }
    
private:
    /** The nodes of this path. */
    vector<NodeID> _nodes;
    
    /** The weights of the edges of this path. */
    vector<EdgeWeight> _weights;

    /** The IDs of the edges. Optional, i.e., may be left empty. */
    vector<OtherEdgeID> _edges;

    /** For each edge in this path, a shortcut flag. */
    vector<bool> _shortcuts;
    
    /** The length (=sum of edge weights) of this path. */
    EdgeWeight _length;
    
    /** 
     * A flag that indicates that two nodes are not connected
     * in case that this path should represent a connection
     * between those two nodes.
     */
    bool _notConnected;

    
    /**
     * Checks the invariant that the number of nodes of a non-empty
     * path is always equal to the number of edges plus 1.
     */
    bool checkInvariant() const {
        return (noOfNodes() == noOfEdges() + 1);
    }

    /** Adds the specified subpath of this path to the given path p. */
    void subpath(Path& p, NodeID begin, NodeID end) const {
        p.addFirstNode(node(begin));
        for (NodeID u = begin+1; u <= end; u++) p.add( node(u), _weights[u-1], _shortcuts[u-1] );
    }
    
};


/** My own path structure, used for highway dimension computation*/
class EasyPath{
	/** Outputs data about the path for debugging purposes. */ 
friend ostream& operator<<( ostream& os, EasyPath& p ) {
	os  << p.noOfNodes() << " nodes in all  " << p.node(0);
	for (NodeID e = 1; e < p.noOfNodes(); e++)
		os << " --> " << p.node(e); 
	os << endl;
	return os;
}
public:
	EasyPath(Path& p)
	{
		NodeID size = p.noOfNodes();
		_enodes.resize(size);
		for (NodeID i = 0; i < size; i++){
			_enodes[i] = p.node(i);
		}
	}
	
	EasyPath(list <NodeID> &li)
	{
		_enodes.resize(li.size());
		list <NodeID>::iterator it = li.begin();
		int i = 0;
		while (it != li.end()){
			_enodes[i] = *it;
			it++;
			i++;
		}
	}
	void setPath(Path& p)
	{
		NodeID size = p.noOfNodes();
		_enodes.resize(size);
		for (NodeID i = 0; i < size; i++){
			_enodes[i] = p.node(i);
		}
	}

	void setPath(EasyPath& p)
	{
		NodeID size = p.noOfNodes();
		_enodes.resize(size);
		for (NodeID i = 0; i < size; i++){
			_enodes[i] = p.node(i);
		}
	}
	NodeID noOfNodes() const{
		return _enodes.size();
	}
	NodeID node(NodeID index){
		return _enodes[index];
	}

	 EasyPath() {clear();}
    
    /**
     * Constructor.
     * Creates a path with a given start node.
     */
    EasyPath(NodeID id) {
        clear();
        addFirstNode(id);
    }

    /** Clears this path (so that it becomes an empty path). */
    void clear() {
        _enodes.clear();
		_notConnected = false;
    }

    /** Adds a start node to this (empty) path. */
    void addFirstNode(NodeID id) {
        _enodes.push_back(id);
    }
	void add(NodeID id) {
		_enodes.push_back(id);
	}

	/**
     * Adds a path to this path.
     * Note: the first node is NOT added; we assume that it is
     * identical to the formerly last node of this path.
     */
    void add(const EasyPath& path) {
        _enodes.insert( _enodes.end(), path._enodes.begin() + 1, path._enodes.end() );
    }

	/**
     * Adds a path to this path.
     * Note: the first node is also added
     */
    void combine(const EasyPath& path) {
        _enodes.insert( _enodes.end(), path._enodes.begin(), path._enodes.end() );
    }
	 /**
     * Reverses this path. (The first node becomes the last one and vice versa.)
     * Note that this might lead to unwanted results if edge IDs are given since
     * only the order of the edge IDs is reversed; the IDs of the reverse edges
     * are not determined.
     */
    void reverse() {
        std::reverse( _enodes.begin(), _enodes.end() );
    }

	NodeID firstNode() const {
		assert( ! empty() );
		return _enodes[0];
	}

	NodeID lastNode() const {
		assert( ! empty() );
		return _enodes[noOfNodes() - 1];
	}
	
	void setNotConnected()
	{
		_notConnected = true;
	}
	NodeID memorySize(){
		return _enodes.size() * sizeof(NodeID);
	}
	void printtofile(FILE *fp)
	{
		for (NodeID i = 0; i < _enodes.size(); i++){
			fprintf(fp,"%d ", _enodes[i] );
		}
		fprintf(fp, "\n");
	}
	~EasyPath()
	{
		_enodes.clear();
	}
private:
	vector <NodeID> _enodes;
	bool _notConnected;
};


/** Represents a collection of paths. */
class Paths : public vector<Path>
{
public:
    /**
     * Compares two collections of paths.
     * @param os the stream that the result of the comparison/the differences are written to
     * @param p the path collection that this collection is compared to
     * @param verbose write details about the differences (true) / only a summary (false)
     */
    void diff(ostream& os, const Paths& p, bool verbose = false) const {
        if (size() != p.size()) {
            os << "Not the same number of paths." << endl;
            return;
        }

        NodeID countDifferentPaths = 0;
        NodeID countDifferentLengths = 0;
        Path diff1, diff2; // store the differences of one comparison of two paths
        vector< pair<Path, Path> > diffs; // stores the differences of all comparisons,
                                          // all entries are unique (duplicates aren't stored)
        
        // compare path i in this collection with path i in the collection p
        for (NodeID i = 0; i < size(); i++) {
            (*this)[i].diff( p[i], diff1, diff2 ); // perform comparison
            if ((! diff1.empty()) || (! diff2.empty())) { // if both paths aren't identical
                countDifferentPaths++;
                if (diff1.length() != diff2.length()) { // if the lengths of both paths differ
                    countDifferentLengths++;
                    if (verbose && (! p[i].empty())) {
                        os << "Different path lengths: " << p[i].node(0) << " -> "
                           << p[i].node(p[i].noOfNodes()-1)
                           << ": " << (*this)[i].length() << " != " << p[i].length() << endl;
                    }
                }
                pair<Path, Path> newDiff( diff1, diff2 );
                bool unique = true;
                // check whether the current differences pair is unique (or a duplicate)
                for (NodeID j = 0; j < diffs.size(); j++)
                    if (diffs[j] == newDiff) {
                        unique = false;
                        break;
                    }
                if (unique) diffs.push_back(newDiff);
            }
        }
        // write summary
        os << countDifferentPaths << " out of " << size() << " pairs of paths differ. ("
           << diffs.size() << " distinct differences)" << endl;

        if (countDifferentLengths == 0) {            
            os << "For each pair of paths, the lengths are equal." << endl << endl;
        }
        else {
            os << endl << "!!! WARNING !!! The path lengths differ in "
               << countDifferentLengths << " cases. !!! WARNING !!!" << endl << endl;
        }
        
        if (verbose) {
            // write details
            for (NodeID i = 0; i < diffs.size(); i++) {
                os << diffs[i].first << diffs[i].second << endl;
            }
        }
    }

    /** 
     * Exports all edges of all paths of this collection to the given stream 
     * using the given color (in order to be able to draw the paths).
     */
    void exportEdges(ostream& out, int color) {
        for (NodeID i = 0; i < size(); i++) (*this)[i].exportEdges(out, color);
    }
};


#endif // PATH_H
