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

#ifndef PQNODE_H
#define PQNODE_H


#include "../pqueue/binaryHeap.h"


/** 
 * Encapsulates data that is associated with a node that is in the
 * pqueue, only stores NodeID
 */
class SimplePQueueNode
{
public:
    /**
     * Initializes this object.
     * @param nodeID the id of the node this object is associated with
     */
    void init(NodeID nodeID) {
        _nodeID = nodeID;
    }

    /** Returns the id of the node this object is associated with. */ 
    NodeID nodeID() const {return _nodeID;}
    
private:
    /** The id of the node this object is associated with. */
    NodeID _nodeID; // possibility to save memory: use parentEdge.target()
};

/** 
 * Encapsulates data that is associated with a node that is in the
 * pqueue (for a certain search direction).
 */
class PQueueNode
{
public:
    /**
     * Initializes this object.
     * @param nodeID the id of the node this object is associated with
     */
    void init(NodeID nodeID) {
        _nodeID = nodeID;
    }

    /**
     * Marks that this node is the start node of a search.
     * This is done by setting the parent pointer so
     * that a self-loop is created.
     */
    void setStartNode() { _parentNode = _nodeID; }

    bool isStartNode() const {return (_parentNode == _nodeID);}

    /**
     * This method is never used. However, the compiler is not able to realize this fact.
     * @see PQueueNodeConstruction::updateParent(...)
     */
    void updateParent(NodeID parentNodeID, OtherEdgeID parentOtherEdgeID, EdgeWeight parentEdgeWeight, EdgeWeight parentKey,
		      const PQueueNode& parentData, EdgeWeight dH, NodeID next) {
	assert(false);
    }

    /**
     * Set a new parent for this node.
     * @param parentNode the parent node
     * @param parentEdge the edge that leads from the parent to this node     
     */
    void updateParent(NodeID parentNode, OtherEdgeID parentEdge) {
        _parentNode = parentNode;
        _parentEdge = parentEdge;
    }

    /** Returns the id of the node this object is associated with. */ 
    NodeID nodeID() const {return _nodeID;}
    /** Returns the id of the parent node. */
    NodeID parentNode() const {return _parentNode;}
    /** Returns the id of the edge that leads from the parent to this node. */
    OtherEdgeID parentEdge() const {return _parentEdge;}

    
    /**
    * The stalling concept is only used in case of the HNR query.
    * In PQueueNodeDynQuery, this method is overwritten appropriately.
    * In all other cases (in particular during a HH query), there are
    * NO stalled nodes.
    * @see DijkstraTemplate::obtainRelevantSearchSpace
    */
    bool stalled() const {return false;}
    
    
    // The following dummy methods are required to provide the same interface
    // in all PQueueNode subclasses so that the compiler does not complain.
    // These dummy methods are never called.
    bool isCovered() const { assert(false); return false; }
    void unstall() { assert(false); }
    EdgeWeight stallKey() const { assert(false); return 0; }
    void stallKey(EdgeWeight k) { assert(false); }
    NodeID hwyNodeCount() const { assert(false); return 0; }
    void setHwyNodeCount(NodeID c) { assert(false); }
    void incHwyNodeCount() { assert(false); }
    void setMaverick() { assert(false); }
    bool isActive() const { assert( false ); return false; }
    void setActive(bool act) { assert( false ); }
    void setLevel(LevelID level, EdgeWeight dist) { assert( false ); }
    LevelID level() const { assert( false ); return 0; }
    EdgeWeight distanceToLevelBorder() const { assert( false ); return 0; }
    EdgeWeight potential() const { assert( false ); return 0; }
    void setPotential(EdgeWeight pot) { assert( false ); }
    void next(const NodeID n) { assert(false); }
    NodeID next() const { assert(false); return 0; }
    bool isTarget() const { assert( false ); return false; }
    void setTarget(bool value) { assert( false ); }
    bool isEliminated() const { assert(false); return false; }
    void setEliminated(bool value) { assert(false); }
    unsigned int hops() const { assert( false ); return 0; }
    void setHops(unsigned int v) { assert( false ); }    
    
private:
    /** The id of the node this object is associated with. */
    NodeID _nodeID; // possibility to save memory: use parentEdge.target()
    /** The id of the parent node. */
    NodeID _parentNode;
    /** The id of the edge that leads from the parent to this node. */
    OtherEdgeID _parentEdge;
};


/**
 * Local search node. Stores several flags and a hop-counter
 * that is used to stop the search erlier.
 */
class PQueueNodeLocalSearch : public PQueueNode
{
public:
    PQueueNodeLocalSearch() 
    : _target(false), 
      _eliminated(false), 
      _hops(0) {}
    
    /** If node is a target of the local search */
    bool isTarget() const { return _target; }
    void setTarget(bool value) { _target = value; }

    
    /* If node is already elminated / contracted. */
    bool isEliminated() const { return _eliminated; }
    void setEliminated(bool value) { _eliminated = value; }
    
    /* Hops (edges) to the source node. */
    unsigned int hops() const { return _hops; }
    void setHops(unsigned int v) { _hops = v; }    
private:
    bool _target:1;
    bool _eliminated:1;
    unsigned int _hops:30;
};


/**
 * Encapsulates data that is associated with a node that is in the
 * pqueue (for a certain search direction) during a highway search.
 */
class PQueueNodeHwySearch : public PQueueNode
{
public:
    PQueueNodeHwySearch() : _potential(0), _active(true) {}
    
    /**
     * Sets the search level and the distance to the border
     * of the applicable neighbourhood.
     */
    void setLevel(LevelID level, EdgeWeight dist) {
        _level = level;
        _distToLevelBorder = dist;
    }

    void setPotential(EdgeWeight pot) {
        _potential = pot;
    }

    /** Returns the search level. */
    LevelID level() const {return _level;}
    
    /** Returns the distance to the border of the applicable neighbourhood. */
    EdgeWeight distanceToLevelBorder() const {return _distToLevelBorder;}

    EdgeWeight potential() const {return _potential;}

    bool isActive() const { return _active; }
    
    void setActive(bool act) { _active = act; }
    
private:
    /** The search level. */
    LevelID _level;
    
    /** 
     * The distance to the border of the applicable neighbourhood.
     * When this border is crossed, the current search level is left.
     */
    EdgeWeight _distToLevelBorder;

    /**
     * Used to cache the potential of this node during the main phase of a landmark query to
     * avoid multiple calls of the potential function.
     */
    EdgeWeight _potential;

    /**
     * A flag used to perform "active/passive searches".
     * Used for transit node routing: determine middle access points.
     */
    bool _active;
};


/**
 * Encapsulates data that is associated with a node that is in the
 * pqueue (for a certain search direction) during the construction
 * of the dynamic method.
 */
class PQueueNodeDynConstr : public PQueueNode
{
private:
    typedef unsigned short CounterType;
    
public:
    PQueueNodeDynConstr() : _hwyNodeCounter(0), _maverick(false) {}

    CounterType hwyNodeCount() const { return _hwyNodeCounter; }
    
    void setHwyNodeCount(const CounterType c) {
        _hwyNodeCounter = c;
        _maverick = false;
    }
    
    void incHwyNodeCount() { _hwyNodeCounter++; }

    bool isActive() const { return ((hwyNodeCount() == 0) && (!_maverick)); }

    bool isCovered() const { return (hwyNodeCount() > 0); }

    void setMaverick() { _maverick = true; }
    
private:
    /**
    * Counts the number of nodes that belong to the overlay graph
    * on the path from the root to this node. Used for the stall-
    * in-advance technique. A node with a counter value >= 1 is
    * covered. A node with a counter value of 1 whose parent is
    * not covered is a covering node.
    */
    CounterType _hwyNodeCounter;

    bool _maverick;
};


/**
 * Encapsulates data that is associated with a node that is in the
 * pqueue (for a certain search direction) during a query of the
 * dynamic method.
 */
class PQueueNodeDynQuery : public PQueueNode
{
public:
    PQueueNodeDynQuery() {unstall();}
    
    bool stalled() const { return (_stallKey != Weight::MAX_VALUE); }
    void unstall() { _stallKey = Weight::MAX_VALUE; }

    EdgeWeight stallKey() const { return _stallKey; }
    void stallKey(EdgeWeight k) { _stallKey = k; }
    
private:
    /**
    * If this node u has been stalled, some path from the source to u
    * has been found that is shorter than the key of u. The length of
    * this path is stored in the 'stallKey'. If the node has not been
    * stalled, the stallKey is 'infinity'.
    */
    EdgeWeight _stallKey;
};


/**
 * Encapsulates data that is associated with a node that is in the
 * pqueue during the construction.
 */
class PQueueNodeConstruction : public PQueueNode
{
private:
    // we rely on MAVERICK > ACTIVE > PASSIVE; thus, don't change the order
    static const char PASSIVE = 0;
    static const char ACTIVE = 1;
    static const char MAVERICK = 2;
    
public:
    /** Constructor. */
    PQueueNodeConstruction() :
        _distToBorderOFs1(Weight::SIGNED_MAX_VALUE), _dist_s0u(Weight::MAX_VALUE), _next(0), _status(ACTIVE), _slackDefined(false) {}

    /**
     * Set a new parent for this node.
     * @param parentNodeID the parent node
     * @param parentOtherEdgeID the edge that leads from the parent to this node
     * @param parentEdgeWeight the weight of the edge from the parent to this node
     * @param parentKey the distance from source to parent
     * @param parentData a reference to the PQueueNodeConstruction object of the parent node
     * @param dH the neighbourhood radius of this node
     * @param next the index of the next parent or 0 if there is only a single parent
     */
    void updateParent(NodeID parentNodeID, OtherEdgeID parentOtherEdgeID, EdgeWeight parentEdgeWeight, EdgeWeight parentKey,
		      const PQueueNodeConstruction& parentData, EdgeWeight dH, NodeID next) {
	assert( dH <= Weight::SIGNED_MAX_VALUE );
	assert( parentEdgeWeight <= Weight::SIGNED_MAX_VALUE );
        
        // if this node is s_1 (i.e., its parent has key 0), then the distance to the neighbourhood border of s_1
        // is set to the neighbourhood radius of s_1; otherwise, we adapt the distance stored at the parent node
	SignedEdgeWeight currentDist = (parentKey == 0) ? dH : (parentData.distToBorderOFs1() - parentEdgeWeight);
	    
	if (next == 0) { // if there is only a single tentative parent...
            // set the single parent
	    PQueueNode::updateParent( parentNodeID, parentOtherEdgeID );
	    distToBorderOFs1( currentDist );
	    // inherit data from the parent
	    dist_s0u( parentData.dist_s0u() );
	    _status = parentData._status;	    
	}
	else { // if there are multiple tentative parents...
            // do not just inherit data from the new parent, but take the maxima of the already stored data
            // and the data stored at the new parent
	    distToBorderOFs1( max(distToBorderOFs1(), currentDist) );
	    dist_s0u( max(dist_s0u(), parentData.dist_s0u()) );
	    _status = max(_status, parentData._status);
	}
	setNext(next);
    }

    /**
     * This method is never used. However, the compiler is not able to realize this fact.
     * @see PQueueNode::updateParent(...)
     */
    void updateParent(NodeID parentNode, OtherEdgeID parentEdge) {
	assert(false);
    }

    /** Returns true iff the neighbourhood of s_1 has been left. */
    bool nbrhOFs1Left() const {return (_distToBorderOFs1 < 0);}

    /** 
     * Returns true iff the distance from s_0 to u has been defined.
     * @see _dist_s0u
     */
    bool defined_dist_s0u() const {return (_dist_s0u != Weight::MAX_VALUE);}
    /**
     * Returns the distance from s_0 to u.
     * @see _dist_s0u
     */
    EdgeWeight dist_s0u() const {return _dist_s0u;}
    /**
     * Sets the distance from s_0 to u.
     * @see _dist_s0u
     */
    void dist_s0u(EdgeWeight d) {_dist_s0u = d;}

        

    /**
     * Initializes the slack.
     * @param sl the dH value of the corresponding node
     * @see _slack
     */
    void initSlack(EdgeWeight sl) {
	assert( sl <= Weight::SIGNED_MAX_VALUE );
	SignedEdgeWeight signedSlack = (SignedEdgeWeight)sl;
	
        if (! _slackDefined) {
            // leaf
	    setSlack(signedSlack);
        }
        else {
            // inner node
            if (signedSlack < slack()) setSlack(signedSlack);
        }
    }
    
    /**
     * Updates the slack.
     * @param child the child from which the update is invoked
     * @param dist the distance from the child to this node
     * @return true iff the slack of this node GETS negative
     *         (note: false if the slack IS already negative, but this
     *          child does not cause a negative value)
     */
    bool updateSlack(const PQueueNodeConstruction& child, EdgeWeight dist) {
	assert( dist <= Weight::SIGNED_MAX_VALUE );
	const SignedEdgeWeight diff = child.slack() - (SignedEdgeWeight)dist;
        assert( (diff < 0) || (child.slack() >= (SignedEdgeWeight)dist) );
        if (diff < 0) {
            // the child causes a negative value of this slack
	    setSlack(diff);
            return true;
        }
        
        // if this slack is already negative, do nothing
        if ((_slackDefined) && (slack() < 0)) return false;
        
        // if this slack hasn't been defined yet
        // or the child causes a decrease of this slack...        
        if ((!_slackDefined) || (diff < slack())) setSlack(diff); // update this slack
        return false;
    }

    /** 
     * Deactivates this node iff the abort criterion is fulfilled. 
     * @param dist the distance from s_0 to p
     * @param dH the value d_H(p)
     */
    void abort(EdgeWeight dist, EdgeWeight dH) {        
        if (defined_dist_s0u()) { // d(s_0, u) has to be defined; otherwise, we can't abort
            assert( dist >= dist_s0u() );
            EdgeWeight w = dist - dist_s0u(); // compute d(u,p)            
            if (w > dH) _status = PASSIVE; // deactivate
        }
    }

    /** Returns true iff this node is active. */
    bool isActive() const { return (_status >= ACTIVE); }

    /** Returns true iff this node is a maverick. */
    bool isMaverick() const { return (_status == MAVERICK); }

    /** Marks that this node is a maverick. */
    void setMaverick() {
        // note: only an active node can be a maverick
        if (isActive()) _status = MAVERICK;
    }

    /** Sets the index of the next parent. */
    void setNext(const NodeID n) {
	_next = (unsigned short)n;
	assert( next() == n );
    }

    /** Returns the index of the next parent. */
    NodeID next() const {return _next;}
    
private:    
    /**
     * The distance to the border of the neighbourhood of s_1.
     * Also used to store the slack during the backward evaluation phase
     * of the construction, where the distance to the neighbourhood border
     * is no longer needed.
     */
    SignedEdgeWeight _distToBorderOFs1;
    
    /** 
     * The distance from s_0 to u.
     * u is the direct predecessor of v, which is the last node within the
     * neighbourhood of s_1.
     */
    EdgeWeight _dist_s0u;
    
    /** The index of the next parent. 0 if there is no next parent. */
    unsigned short _next;

    /** The status of the node (passive, active, maverick). */
    char _status;

    /** True iff the slack has already been defined. */
    bool _slackDefined;
    

    
    /** Returns the distance to the border of the neighbourhood of s_1. */
    SignedEdgeWeight distToBorderOFs1() const {return _distToBorderOFs1;}

    /** Sets the distance to the border of the neighbourhood of s_1. */
    void distToBorderOFs1(const SignedEdgeWeight d) {_distToBorderOFs1 = d;}

    /** Sets the slack. */
    void setSlack(const SignedEdgeWeight s) {
        // we use _distToBorderOFs1 to store the slack in order to save memory
	_distToBorderOFs1 = s;
	_slackDefined = true;
    }

    /** Returns the slack. */
    SignedEdgeWeight slack() const {
	assert( _slackDefined );
	return _distToBorderOFs1;
    }
};

/** Priority Queue that is used by the original version of Dijkstra's algorithm. */
typedef BinaryHeap< EdgeWeight, Weight, PQueueNode, NodeID > NormalPQueue;
/** Priority Queue that is used by the multilevel query algorithm ("highway search"). */
typedef BinaryHeap< EdgeWeight, Weight, PQueueNodeHwySearch, NodeID > HwyPQueue;
/** Priority Queue that is used during the construction. */
typedef BinaryHeap< EdgeWeight, Weight, PQueueNodeConstruction, NodeID > ConstrPQueue;

typedef BinaryHeap< EdgeWeight, Weight, PQueueNodeDynConstr, NodeID > DynConstrPQueue;
typedef BinaryHeap< EdgeWeight, Weight, PQueueNodeDynQuery, NodeID > DynQueryPQueue;

/** Priority Queue that is used during the node contraction of CH for the local searches. */
typedef BinaryHeap< EdgeWeight, Weight, PQueueNodeLocalSearch, NodeID > LocalSearchPQueue;

typedef NormalPQueue::PQElement NormalPQElement;
typedef HwyPQueue::PQElement HwyPQElement;
typedef ConstrPQueue::PQElement ConstrPQElement;
typedef DynConstrPQueue::PQElement DynConstrPQElement;
typedef DynQueryPQueue::PQElement DynQueryPQElement;

#endif // PQNODE_H
