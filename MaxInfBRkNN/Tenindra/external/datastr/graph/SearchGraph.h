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
 
#ifndef _DATASTR_GRAPH_SEARCHGRAPH_H
#define _DATASTR_GRAPH_SEARCHGRAPH_H

#include <math.h>
#include "UpdateableGraph.h"

namespace datastr { namespace graph {
    
/**
 * Search graph node order. A good node order can speedup the query.
 */
enum SearchGraphNodeOrder
{
    // original node id
    SGNO_ORIGINAL,
    
    // level id (distinct level for each node)
    SGNO_LEVEL,
    
    // Special order, large partitions according to level
    // order according to original node id within level.
    SGNO_LOGLEVEL_ORIGINAL,
};

/**
 * Search graph for queries. It stores only the necessary information
 * for queries, the edges incident to an node with the other edge in the same original
 * higher level. Its interface is closly related to the Static/DynamicGraph.
 * Internally, an adjacency array is kept.
 */  
class SearchGraph
{
public:
    /** Represents a Node in the SearchGraph. */
    class SearchNode
    {
    public:
        /** Default Constructor. */
        SearchNode() : _firstLevelEdge(0), _pqElement(0) {}

        /**
         * Constructor.
         */
        SearchNode(const OtherEdgeID fLE) {
            _firstLevelEdge = fLE;
        }

        /** Returns the index of the first edge leaving this node. */
        OtherEdgeID firstLevelEdge() const { return _firstLevelEdge; }
        /** Sets the index of the first edge leaving this node. */
        void setFirstLevelEdge(const OtherEdgeID fLE) { _firstLevelEdge = fLE; }

        /** Returns the index of the corresponding elements in the pqueues. */
        NodeID pqElement() const {return _pqElement & ((1U<<31)-1);}
    
        /** Sets the index of the corresponding elements in the pqueues. */
        void pqElement(const NodeID pqE) 
        {
            assert( pqE < ((unsigned int)1<<31) );
            _pqElement = (pqE | (_pqElement & (1<<31)));
            assert( pqElement() == pqE );
        }

        /** Returns wheter the node is in the core (the topmost level) */
        bool isInCore() const { return (_pqElement & (1U<<31)) != 0; }
        /** Sets wheter the node is in the core (the topmost level) */
        void setInCore(bool inCore)
        {
            _pqElement = (pqElement() | (inCore << 31));
        }

        /** 
         * Serialize the node. 
         * Warning: it depends on the node representation in the main memory.
         */
        void serialize(ostream& out) const {
			//cout << this->_firstLevelEdge << " " << this->_pqElement << endl;
            out.write((char*)this,sizeof(SearchNode)/sizeof(char));
        }

        /** 
         * Deserialize the node. 
         * Warning: it depends on the node representation in the main memory.
         */
        void deserialize(istream& in) {
            in.read((char*)this,sizeof(SearchNode)/sizeof(char));
        }

        /** Not used. Provides same interface as the UpdateableGraph. */
        bool isTarget() const { assert(false); return false; }

        
    private:
        OtherEdgeID _firstLevelEdge;
        NodeID _pqElement;
    };

    /** The node type used in this graph. */
    typedef SearchNode MyNode;

    
public:
    /** Constructor. Builds the graph. */
    SearchGraph(UpdateableGraph* updGraph, const SearchGraphNodeOrder nodeOrder) {
        construct(updGraph, nodeOrder);
    }
    
    /** Constructor. Builds the graph. */
    SearchGraph(istream& in)
    {
        deserialize(in);
    }


    /** Returns the number of nodes. */
    NodeID noOfNodes() const {return (_nodes.size()-1) /*substract dummy node*/ ;}

    /** 
     * Returns the number of edges that memory has been allocated for. 
     * For a search graph that is the sam as the number of nodes since
     * there are no holes in the edge array.
     */
    OtherEdgeID noOfEdges() const {return _edges.size();}

    /** Returns the index of the first edge of u. */
    OtherEdgeID firstEdge(const NodeID u) const {
        assert( false );
        return SPECIAL_NODEID;
    }

    /** Returns the index of the first edge of u. */
    OtherEdgeID firstLevelEdge(const NodeID u) const {
        assert( u < noOfNodes() );
        return _nodes[u].firstLevelEdge();
    }

    /** Returns the index+1 of the last edge of u. */
    OtherEdgeID lastEdge(const NodeID u) const {
        assert( u < noOfNodes() );
        return _nodes[u+1].firstLevelEdge();
    }
    
    /**
    * Returns the number of edges stored at node u
    * irrespective of the direction.
    */
    OtherEdgeID degree(const NodeID u) const {
        assert( u < noOfNodes() );
        return lastEdge(u)-firstEdge(u);
    }

    /**
    * Returns the number of existing edges.
    */
    OtherEdgeID noOfExistingEdges() const {
        return noOfEdges();
    }

    const MyNode& node(const NodeID u) const {
        assert( u < noOfNodes() );
        return _nodes[u];
    }

    MyNode& node(const NodeID u) {
        assert( u < noOfNodes() );
        return _nodes[u];
    }

    const Edge& edge(const OtherEdgeID e) const {
        assert( e < _edges.size() );
        return _edges[e];
    }

    Edge& edge(const OtherEdgeID e) {
        assert( e < _edges.size() );
        return _edges[e];
    }

    /**
    * The nodes can have internally a differnt
    * node id to increase the performance (cache).
    * Returns for a given original NodeID, the
    * internally used NodeID.
    */
    NodeID mapExtToIntNodeID(const NodeID ext) const {
        assert( ext < _mapExtToIntNodeIDs.size() );
        return _mapExtToIntNodeIDs[ext];
    }

    /** Serializes the graph to the given stream. */
    void serialize(ostream& out) {
        VERBOSE( cout << "datastr::graph::SearchGraph::serialize " << noOfNodes()
                      << " " << noOfEdges() << endl );

        // nodes
        VectorSerializer< SearchNode, NodeID, ComplexSerializer<SearchNode> >::serialize(out, _nodes);

        // edges
        VectorSerializer< Edge, OtherEdgeID, ComplexSerializer<Edge> >::serialize(out, _edges);

        // node-id mapping        
        VectorSerializer< NodeID, NodeID >::serialize(out, _mapExtToIntNodeIDs);

        VERBOSE( cout << "done." << endl );
        VERBOSE( printMemoryUsage(cout) );
    }

    /** Deserializes the graph from the given stream. */
    void deserialize(istream& in) {
//         VERBOSE( cout << "datastr::graph::SearchGraph::deserialize " << flush );

        // nodes
        VectorSerializer< SearchNode, NodeID, ComplexSerializer<SearchNode> >::deserialize(in, _nodes);
        
        // edges
        VectorSerializer< Edge, NodeID, ComplexSerializer<Edge> >::deserialize(in, _edges);

        // node-id mapping
        VectorSerializer< NodeID, NodeID >::deserialize(in, _mapExtToIntNodeIDs);

//         VERBOSE( cout << noOfNodes() << " " << noOfEdges() << endl; )
//         VERBOSE( cout << "done." << endl );
    }        

private:
    vector<MyNode> _nodes;

    vector<Edge> _edges;

    /** Maps original ('external') node IDs to the IDs used internally. */
    vector<NodeID> _mapExtToIntNodeIDs;


    /** Sort after the second element of a pair. */
    template < class T1, class T2 >
    struct CompareSecond : public binary_function<pair<T1, T2>, pair<T1, T2>, bool>
    {
        bool operator()(const pair<T1, T2>& a, const pair<T1, T2>& b)
        {
            return (a.second < b.second);
        }
    };
    
    
    /**
    * Construct the SearchGraph from an UpdateableGraph. Only edges having their
    * other incident edge in the same or higher level are kept because only these
    * edges are required for a correct query algorithm. Also a new, more cache-efficient
    * enumeration of the nodes (node ids) can be computed.
    * This nodeOrder has no influence on the node order used for contraction hierarchy construction.
    */
    void construct(UpdateableGraph* updGraph, const SearchGraphNodeOrder nodeOrder = SGNO_ORIGINAL) {
        LevelID levels = (LevelID)(log((double)updGraph->noOfNodes())/log((double)2))-4;
        assert( levels > 0 );
        switch (nodeOrder)
        {
            case SGNO_ORIGINAL:
                VERBOSE( cout << "Node-order: original node-id" << endl; )
                break;
            case SGNO_LEVEL:
                VERBOSE( cout << "Node-order: level" << endl; )
                break;
            case SGNO_LOGLEVEL_ORIGINAL:
                VERBOSE( cout << "Node-order: " << levels << " levels + original node-id" << endl; )
                break;
            default:
                VERBOSE( cout << "Node-order: unknown (" << nodeOrder << ")" << endl; )
                break;
        }
            

        // Store the pair <node,level of node> in an array. 
        vector< pair<LevelID, NodeID> > levelNodes(updGraph->noOfNodes());
        for (NodeID u = 0; u < updGraph->noOfNodes(); u++) {
            levelNodes[u].first = updGraph->node(u).level();
            levelNodes[u].second = u;
        }

        if (nodeOrder == SGNO_LEVEL || nodeOrder == SGNO_LOGLEVEL_ORIGINAL )
        {
            // sort levels after level-id
            sort( levelNodes.begin(), levelNodes.end() );
    
    
            // If a new node order is requested.
            if ( nodeOrder == SGNO_LOGLEVEL_ORIGINAL )
            {
                // partition nodes into levels, level 0 size = n/2, level 1 size = n/4,...
                // within the partition sort by original node id
                NodeID oldBorder = 0;
                NodeID step = updGraph->noOfNodes()/2;
                NodeID border = step;
                LevelID partition = 0;
                CompareSecond<LevelID,NodeID> comp;
                
                while (border <= levelNodes.size() && partition < levels)
                {
                    LevelID exactLevel = levelNodes[border-1].first;
                    
                    // do not split existing levels
                    while (border < levelNodes.size() && levelNodes[border].first == exactLevel) border++;
        
                    // sort by node id within the partition            
                    sort( levelNodes.begin() + oldBorder, levelNodes.begin() + border, comp );
        
                    if (border == levelNodes.size()) break;
                    partition++;
                    oldBorder = border;
        
                    step /= 2;
                    if (step > 0)
                    {                
                        border = border + step;
                        if ( border >levelNodes.size() ) border = levelNodes.size();
                    }
                    else
                    {
                        border = levelNodes.size();
                    }
                }
            }

        }
        
        // node-id mapping
        _mapExtToIntNodeIDs.resize(updGraph->noOfNodes());
        for (NodeID i = 0; i < levelNodes.size(); i++) {
            _mapExtToIntNodeIDs[levelNodes[i].second] = i;
        }
                

        // build adjacency array
        _nodes.resize(levelNodes.size()+1);
        NodeID noOfNodesInCore = 0;
        for ( NodeID i = 0; i < levelNodes.size(); i++ )
        {
            OtherEdgeID iFirstEdge = _edges.size();
            _nodes[i].setFirstLevelEdge(iFirstEdge);
            
            // a node is in the core level if the level == n
            _nodes[i].setInCore(levelNodes[i].first == (LevelID)updGraph->noOfNodes());
            if (_nodes[i].isInCore()) noOfNodesInCore++;
            
            NodeID u = levelNodes[i].second; 
            OtherEdgeID lastEdge = updGraph->lastEdge(u);
            for (OtherEdgeID e = updGraph->firstLevelEdge(u); e < lastEdge; e++) {
                const Edge& edge = updGraph->edge(e);
                
                // do not add witness shortcuts
                if ( edge.type() == EDGE_TYPE_WITNESS_SHORTCUT ) continue;
                    
                assert( updGraph->node(edge.target()).level() >= updGraph->node(u).level() );
                assert( edge.target() < updGraph->noOfNodes() );
                
                _edges.push_back(Edge(mapExtToIntNodeID(edge.target()),
                                      edge.weight(),
                                      edge.type(),
                                      edge.isDirected(0),
                                      edge.isDirected(1),
                                      CH_EXPAND(edge.isShortcut() ? mapExtToIntNodeID(edge.shortcutMiddle()) : ) SPECIAL_NODEID,  // warning: dangerous syntax
                                      edge.shortcutEdge1(),
                                      edge.shortcutEdge2()));

            }
        }
        // guard border node at the end because if
        // lastEdge(n-1) is accessed, firstLevelEdge(n) is returned.
        _nodes[levelNodes.size()].setFirstLevelEdge(_edges.size()); 
        
        VERBOSE( printMemoryUsage(cout); )
        VERBOSE( cout << "#core nodes: " << noOfNodesInCore << endl; )

    }

    /** Used for debugging purposes. Prints all edges (u,v). */
    void printAllEdges(ostream& out, const NodeID u) const {
        out << "all edges from " << u << ": ";
        bool first = true;
        OtherEdgeID lastE = lastEdge(u);
        for (OtherEdgeID e = firstEdge(u); e < lastE; e++) {
            if (! first) out << " | ";
            first = false;

            out << edge(e) << " [" << edge(e).type()  << "]";
        }
        out << endl;
    }

public:
    /** Prints the total memory usage of this UpdateableGraph. */
    void printMemoryUsage(ostream& out) const {
        out << "** memory usage on hard disk [MB (bytes per node)] **" << endl;
        
        unsigned int memoryTotal = 0;
        memoryTotal += printMemoryUsage(out, "edges", noOfEdges() * sizeof(Edge));
        memoryTotal += printMemoryUsage(out, "node data", _nodes.size() * sizeof(SearchNode));
        memoryTotal += printMemoryUsage(out, "node ID mapping", _mapExtToIntNodeIDs.size() * sizeof(NodeID));
        
        printMemoryUsage(out, "TOTAL", memoryTotal);
        out << endl;
    }

    /** Prints the memory usage of one particular data structure of this UpdateableGraph. */
    unsigned int printMemoryUsage(ostream& out, const string descr, const unsigned int mem) const {
        double megaBytes = (double)((mem / (double)(1024*1024)) + 0.5);
        double bytesPerNode = (double)((mem / (double)noOfNodes()) + 0.5);
        out << "   " << descr.c_str() << ": " << megaBytes << " (" << bytesPerNode << ")" << endl;
        return mem;
    }
    
    double getIndexSize() {
        double memoryTotal = 0;
        memoryTotal += noOfEdges() * sizeof(Edge);
        memoryTotal += _nodes.size() * sizeof(SearchNode);
        memoryTotal += _mapExtToIntNodeIDs.size() * sizeof(NodeID);
        return memoryTotal/(1024*1024);
    }

};

} } // namespace

#endif // _DATASTR_GRAPH_SEARCHGRAPH_H
