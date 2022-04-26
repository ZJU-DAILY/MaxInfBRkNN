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

#ifndef EDGE_H
#define EDGE_H

#include "../../io/serialize.h"

/** Normal edge. */
const static NodeID EDGE_TYPE_ORIGINAL = 0;
/** Normal shortcut edge. */
const static NodeID EDGE_TYPE_SHORTCUT = 1;
/** Shortcut edge representing a witness path. Such shortcuts
 *  are used during construction if possible witnesses are known
 *  in advance.
 */
const static NodeID EDGE_TYPE_WITNESS_SHORTCUT = 2;

/** 
 * Represents a weighted edge assuming the source node is known,
 * i.e. only the target node is stored.
 * Base class for both EdgeCHFast
 */
template <NodeID TARGET_SHIFT>
class EdgeTemplate
{
    /** Outputs data about the edge for debugging purposes. */
    friend std::ostream& operator<<( std::ostream& os, const EdgeTemplate &edge ) {
        os << "(";
        if (edge.isDirected(1)) os << "<";
        os << "-";
        if (edge.isDirected(0)) os << ">";
        os << edge.target() << ", " << edge.weight() << ")";
        return os;
    }
    
public:
    /**
     * Returns true iff the given edge and this edge are equal with respect
     * to target and direction, but NOT necessarily weight.
     * Note that edge1 == edge2 does not exclude edge1 "<" edge2.
     */
    bool operator== (const EdgeTemplate& e) const {
        return ((target() == e.target()) &&
                (isDirected(0) == e.isDirected(0)) && (isDirected(1) == e.isDirected(1)));
    }
    
    /** 
     * Compares two edges first by target, in case of equality by weight.
     * Returns false if both edges are identical.
     */
    bool operator< (const EdgeTemplate& e) const {
        if (target() == e.target()) return (weight() < e.weight());
        return (target() < e.target());
    }

    /** Mark (during construction) that this edge belongs to the highway network. */
    void setHighwayEdge() {setFlag(HWY);}
    
    /** Reset the flag that this edge belongs to the highway network. */
    void unsetHighwayEdge() {unsetFlag(HWY);}

    NodeID target() const {return (_target >> TARGET_SHIFT);}
    EdgeWeight weight() const {return _weight;}

    bool isHighwayEdge() const {return (isFlag(HWY));}

    /** Returns true iff this edge is open in the given direction. */
    bool isDirected(int searchID) const {return (isFlag(direction(searchID)));}
    
    /** Returns true iff this edge is open in both directions (two-way). */
    bool isBidirected() const {return (isDirected(0) && isDirected(1));}

    bool isClosed() const {return ( (!isDirected(0)) && (!isDirected(1)) );}
    
    void makeTwoWay() {
        setFlag(direction(0));
        setFlag(direction(1));
    }
    
    void makeOneWay(int searchID) {
        setFlag(direction(searchID));
        unsetFlag(direction(1-searchID));
    }

    void makeClosed() {
        unsetFlag(direction(0));
        unsetFlag(direction(1));
    }

    /**
    * Ensures that this edge is at least open in the same
    * directions as the given edge.
    */
    void makeOpen(const EdgeTemplate& e) {
        if (e.isDirected(0)) setFlag(direction(0));
        if (e.isDirected(1)) setFlag(direction(1));
    }

    /**
    * Ensures that this edge is at least open in the opposite
    * directions as the given edge.
    */
    void makeOpenReversed(const EdgeTemplate& e) {
        if (e.isDirected(0)) setFlag(direction(1));
        if (e.isDirected(1)) setFlag(direction(0));
    }

    /**
    * Ensures that this edge is not open for any direction
    * that the given edge is open.
    */
    void makeClosed(const EdgeTemplate& e) {
        if (e.isDirected(0)) unsetFlag(direction(0));
        if (e.isDirected(1)) unsetFlag(direction(1));
    }

    /** Serializes this edge to a given stream. */
    void serialize(ostream& out) const {
        writePrimitive(out, _target);
        writePrimitive(out, _weight);
    }

    /** Deserializes this edge from a given stream. */
    void deserialize(istream& in) {
        readPrimitive(in, _target);
        readPrimitive(in, _weight);
    }
    
protected:
    // counting from LSB (1) to MSB (32):
    // bit 1 and 2
    NodeID direction(int searchID) const {
        assert( (searchID >= 0) && (searchID < 2) );
        return searchID+1;
    }
    // bit 3
    static const NodeID HWY = 4;
    
    bool isFlag(NodeID flag) const {return ((_target & flag) != 0);}
    void setFlag(NodeID flag) {_target |= flag;}
    void unsetFlag(NodeID flag) {_target &= (~flag);}
    
    void setTarget(const NodeID t) 
    {
        _target = t << TARGET_SHIFT; 
        // treat SPECIAL_NODEID for unused edges different
        assert( t == SPECIAL_NODEID || target() == t );
    }
    
    /**
     * The target node.
     * In addition, some flags.
     */
    NodeID _target;
    
    /** The edge weight. */
    EdgeWeight _weight;
};


/** 
 * Represents a weighted edge assuming the source node is known,
 * i.e. only the target node is stored.
 * Used for Contraction Hierachies.
 *
 * There are two different edge classes for CH, EdgeCHFast and EdgeCH.
 * As the name says, EdgeCHFast is faster than EdgeCH but also more compicated.
 * The difference is how to store multiple values in one integer:
 *   EdgeCHFast: manually write the code to perform shift and mask operations
 *   EdgeCH: let the compiler generate this code. This is done by
 *           adding the number of bytes to each variable with a colon,
 *           e.g. NodeID _type:2; requires only 2 bytes.
 *
 */
class EdgeCHFast : public EdgeTemplate<4>
{
public:
    /** Default constructor. target and weight are set to 0.*/
    EdgeCHFast() {
        _target = 0;
        _weight = 0;
        assert( isClosed() );
        setShortcutOriginalEdgeCount(1);
    }
    
    /** 
     * Constructor. 
     * @param target the target node
     * @param weight the weight
     * @param type edge type (see constants EDGE_TYPE_...)
     * @param forward from source to target open
     * @param backward from target to source open
     * @param shortcutMiddle Not used, see EdgeCHExpand.
     * @param shortcutEdge1 Not used, see EdgeCHExpand.
     * @param shortcutEdge2 Not used, see EdgeCHExpand.
     * @param shortcutEdgeOriginalEdgeCount used for a priority term during CH node ordering
     */
    EdgeCHFast(NodeID target, EdgeWeight weight, NodeID type, bool forward, bool backward,
        NodeID shortcutMiddle = SPECIAL_NODEID, OtherEdgeID shortcutEdge1 = SPECIAL_NODEID, 
        OtherEdgeID shortcutEdge2 = SPECIAL_NODEID, OtherEdgeID shortcutEdgeOriginalEdgeCount = 1) {
        setTarget(target);
        if (forward) setFlag(direction(0));
        if (backward) setFlag(direction(1));
        setType(type);
        
        _weight = weight;
        setShortcutOriginalEdgeCount(shortcutEdgeOriginalEdgeCount);
    }

    /** Constructs the reverse edge of the given edge. */
    EdgeCHFast(NodeID target, const EdgeCHFast& reverseEdge) {
        setTarget(target);
        if (reverseEdge.isDirected(1)) setFlag(direction(0));
        if (reverseEdge.isDirected(0)) setFlag(direction(1));
        setType(reverseEdge.type());
        
        _weight = reverseEdge.weight();
        copyShortcutInfo(reverseEdge, true);
    }

    /** Weight of the edge */
    void setWeight(EdgeWeight w) {_weight = w;}

    /** Is a shortcut edge. */
    bool isShortcut() const {
        return (type() != 0);
    }

    /** Type of the edge, see EDGE_TYPE_... constanst above. */
    NodeID type() const {return ((_target >> TYPE_SHIFT) & 3);}
    void setType(NodeID k) {
        assert( k < 4 );
        _target -= ( type() << TYPE_SHIFT );
        _target += (k << TYPE_SHIFT);
        assert( type() == k );
    }

    bool isIdentical(const EdgeCHFast& e) const {
        return ((*this == e) && (weight() == e.weight()) && (type() == e.type()));
    }

    bool isReverse(NodeID u, const EdgeCHFast& e) const {
        return ((target() == u) && (weight() == e.weight()) &&
                (isDirected(0) == e.isDirected(1)) && (isDirected(1) == e.isDirected(0)) &&
                (type() == e.type()));
    }

    /** 
     *  Methods required to unpack shortcuts. Not supported.
     *  EdgeCHExpand supports them (switch USE_CH_EXPAND on in config.h)
     */
    void setShortcutMiddle(NodeID middle) { }
    NodeID shortcutMiddle() const { return SPECIAL_NODEID; }
    void setShortcutEdge1(OtherEdgeID e) { }
    NodeID shortcutEdge1() const { return shortcutEdgeLimit(); }
    void setShortcutEdge2(OtherEdgeID e) { }
    NodeID shortcutEdge2() const { return shortcutEdgeLimit(); }
    NodeID shortcutEdgeLimit() const { return SPECIAL_NODEID; }
    void copyShortcutInfo(const EdgeCHFast& e, bool reverse) 
    {
        #ifdef COUNT_SHORTCUT_ORIGINAL_EDGES
        _shortcutOriginalEdgeCount = e.shortcutOriginalEdgeCount();
        #endif
    }

    /** 
     * Counter for the number of original edges a shortcut represents.
     * Used as priority term during node ordering.
     * (switch COUNT_SHORTCUT_ORIGINAL_EDGES on in config.h)
     */
    void setShortcutOriginalEdgeCount(const OtherEdgeID count)
    {
        #ifdef COUNT_SHORTCUT_ORIGINAL_EDGES
        _shortcutOriginalEdgeCount = count;           
        #endif
    }
    OtherEdgeID shortcutOriginalEdgeCount() const
    {
        #ifdef COUNT_SHORTCUT_ORIGINAL_EDGES
        return _shortcutOriginalEdgeCount;
        #else
        return SPECIAL_NODEID;
        #endif
    }
    
protected:
    // counting from LSB (1) to MSB (32):
    // 1,2: direction
    // 3,4: type
    //
    static const NodeID TYPE_SHIFT = 2;
    
    #ifdef COUNT_SHORTCUT_ORIGINAL_EDGES
    OtherEdgeID _shortcutOriginalEdgeCount;
    #endif
};


/** 
 * Represents a weighted edge assuming the source node is known,
 * i.e. only the target node is stored.
 * Used for Contraction Hierarchies.
 * Despite the name, EdgeCHExpand seems to be faster
 * than EdgeCHExpandFast. This class is currently Not
 * used.
 */
class EdgeCHExpandFast : public EdgeCHFast
{
public:
    /** Default constructor. target and weight are set to 0.*/
    EdgeCHExpandFast()
        : EdgeCHFast() {
        _shortcut = 0;
    }
    
    /** 
     * Constructor.
     * @param target the target node
     * @param weight the weight
     * @param type edge type (see constants EDGE_TYPE_...)
     * @param forward from source to target open
     * @param backward from target to source open
     * @param shortcutMiddle middle nod of the path a shortcut represents
     * @param shortcutEdge1 relative index into the edge array of the first edge of the path
     * @param shortcutEdge2 relative index into the edge array of the second edge of the path
     * @param shortcutEdgeOriginalEdgeCount used for a priority term during CH node ordering
     */
    EdgeCHExpandFast(NodeID target, EdgeWeight weight, NodeID type, bool forward, bool backward, 
        NodeID shortcutMiddle = SPECIAL_NODEID, OtherEdgeID shortcutEdge1 = SPECIAL_NODEID, 
        OtherEdgeID shortcutEdge2 = SPECIAL_NODEID, OtherEdgeID shortcutOriginalEdgeCount = 1) :
        EdgeCHFast(target,weight,type,forward,backward) {
        setShortcutMiddle(shortcutMiddle);
        setShortcutEdge1(shortcutEdge1);
        setShortcutEdge2(shortcutEdge2);
        setShortcutOriginalEdgeCount(shortcutOriginalEdgeCount);
    }

    /** Constructs the reverse edge of the given edge. */
    EdgeCHExpandFast(NodeID target, const EdgeCHExpandFast& reverseEdge) :
        EdgeCHFast(target, reverseEdge) {
        copyShortcutInfo(reverseEdge, true);
    }
    EdgeWeight weight() const {return _weight & WEIGHT_MASK;}
    void setWeight(EdgeWeight w) {
        assert( w < (1 << (32-SHORTCUT_SHIFT)) );
        _weight = (_weight & ~WEIGHT_MASK) | w;
        assert( w == weight() );
    }
    

    // ***
    // To expand a path with a recursive unpacking routine,
    // three values are used, the middle node (shortcutMiddle)
    // and two relative indices into the edge array of the middle
    // node (shortcutEdge1, shortcutEdge2) that form the path
    // the shortcut represents.
    // ***

    void setShortcutMiddle(NodeID middle) { 
        assert( middle == SPECIAL_NODEID || middle < (1 << (32-SHORTCUT_SHIFT)) );
        _shortcut = (_shortcut & SHORTCUT_EDGE_LIMIT) | (middle << SHORTCUT_SHIFT); 
        assert( middle == SPECIAL_NODEID || shortcutMiddle() == middle );
    }
    NodeID shortcutMiddle() const { return _shortcut >> SHORTCUT_SHIFT; }
    
    void setShortcutEdge1(OtherEdgeID shortcutEdge1) 
    {
        if (shortcutEdge1 < SHORTCUT_EDGE_LIMIT)
        {
            _weight = weight() | (shortcutEdge1 << (32-SHORTCUT_SHIFT));
            assert( this->shortcutEdge1() == shortcutEdge1 );
        }
        else
        {
            _weight = _weight | (SHORTCUT_EDGE_LIMIT << (32-SHORTCUT_SHIFT));
            assert( this->shortcutEdge1() == shortcutEdgeLimit() );
        }
    }
    OtherEdgeID shortcutEdge1() const { return _weight >> (32-SHORTCUT_SHIFT); }

    void setShortcutEdge2(OtherEdgeID shortcutEdge2) 
    { 
        if (shortcutEdge2 < SHORTCUT_EDGE_LIMIT)
        {
            _shortcut = (_shortcut & ~SHORTCUT_EDGE_LIMIT) | shortcutEdge2;
            assert( this->shortcutEdge2() == shortcutEdge2 );
        }
        else
        {
            _shortcut = _shortcut | SHORTCUT_EDGE_LIMIT;
            assert( this->shortcutEdge2() == shortcutEdgeLimit() );
        }
    }
    OtherEdgeID shortcutEdge2() const { return _shortcut & SHORTCUT_EDGE_LIMIT; }
    
    OtherEdgeID shortcutEdgeLimit() const { return SHORTCUT_EDGE_LIMIT; }

    void copyShortcutInfo(const EdgeCHExpandFast& e, bool reverse) 
    { 
        setShortcutMiddle(e.shortcutMiddle());
        if ( true || !reverse )
        {
            setShortcutEdge1(e.shortcutEdge1());
            setShortcutEdge2(e.shortcutEdge2());
        }
        else
        {
            setShortcutEdge1(e.shortcutEdge2());
            setShortcutEdge2(e.shortcutEdge1());
        }
        #ifdef COUNT_SHORTCUT_ORIGINAL_EDGES
        _shortcutOriginalEdgeCount = e.shortcutOriginalEdgeCount();           
        #endif
    }
    
    /** Serializes this edge to a given stream. */
    void serialize(ostream& out) const {
        // not tested
        out.write((char*)this,sizeof(EdgeCHExpandFast)/sizeof(char));
    }

    /** Deserializes this edge from a given stream. */
    void deserialize(istream& in) {
        // not tested
        in.read((char*)this,sizeof(EdgeCHExpandFast)/sizeof(char));
    }
    
    
    
protected:
    OtherEdgeID _shortcut;
    static const NodeID SHORTCUT_SHIFT = 6;
    static const NodeID SHORTCUT_EDGE_LIMIT = (1 << SHORTCUT_SHIFT) - 1;
    static const NodeID WEIGHT_MASK = (1 << (32-SHORTCUT_SHIFT)) - 1;
};

/** 
 * Represents a weighted edge assuming the source node is known,
 * i.e. only the target node is stored.
 * Used for Contraction Hierachies.
 *
 * There are two different edge classes for CH, EdgeCHFast and EdgeCH.
 * See EdgeCHFast for the differences.
 */
class EdgeCH
{

    /** Outputs data about the edge for debugging purposes. */
    friend std::ostream& operator<<( std::ostream& os, const EdgeCH &edge ) {
        os << "(";
        if (edge.isDirected(1)) os << "<";
        os << "-";
        if (edge.isDirected(0)) os << ">";
        os << edge.target() << ", " << edge.weight() << ")";
        return os;
    }

public:
    /** Default constructor. target and weight are set to 0.*/
    EdgeCH() {
        _target = 0;
        _weight = 0;
        _type = 0;
        _flags = 0;
        setShortcutOriginalEdgeCount(1);
        assert( isClosed() );
    }
    
    /** 
     * Constructor.
     * @param target the target node
     * @param weight the weight
     * @param type edge type (see constants EDGE_TYPE_...)
     * @param forward from source to target open
     * @param backward from target to source open
     * @param shortcutMiddle Not used, see EdgeCHExpand.
     * @param shortcutEdge1 Not used, see EdgeCHExpand.
     * @param shortcutEdge2 Not used, see EdgeCHExpand.
     * @param shortcutEdgeOriginalEdgeCount used for a priority term during CH node ordering
     */
    EdgeCH(NodeID target, EdgeWeight weight, NodeID type, bool forward, bool backward, 
        NodeID shortcutMiddle = SPECIAL_NODEID, OtherEdgeID shortcutEdge1 = SPECIAL_NODEID, 
        OtherEdgeID shortcutEdge2 = SPECIAL_NODEID, OtherEdgeID shortcutEdgeOriginalEdgeCount = 1) {
        setTarget(target);
        _flags = 0;
        if (forward) setFlag(direction(0));
        if (backward) setFlag(direction(1));
        setType(type);
        _weight = weight;
        setShortcutOriginalEdgeCount(shortcutEdgeOriginalEdgeCount);
        
    }

    /** Constructs the reverse edge of the given edge. */
    EdgeCH(NodeID target, const EdgeCH& reverseEdge) {
        setTarget(target);
        _flags = 0;
        makeOpenReversed(reverseEdge);
        setType(reverseEdge.type());
        _weight = reverseEdge.weight();
        copyShortcutInfo(reverseEdge, true);
    }


    /**
     * Returns true iff the given edge and this edge are equal with respect
     * to target and direction, but NOT necessarily weight.
     * Note that edge1 == edge2 does not exclude edge1 "<" edge2.
     */
    bool operator== (const EdgeCH& e) const {
        return ((target() == e.target()) &&
                (isDirected(0) == e.isDirected(0)) && (isDirected(1) == e.isDirected(1)));
    }
    
    /** 
     * Compares two edges first by target, in case of equality by weight.
     * Returns false if both edges are identical.
     */
    bool operator< (const EdgeCH& e) const {
        if (target() == e.target()) return (weight() < e.weight());
        return (target() < e.target());
    }

    NodeID target() const {return _target; }
    EdgeWeight weight() const {return _weight; }

    /** Returns true iff this edge is open in the given direction. */
    bool isDirected(int searchID) const {return (isFlag(direction(searchID)));}
    
    /** Returns true iff this edge is open in both directions (two-way). */
    bool isBidirected() const {return (isDirected(0) && isDirected(1));}

    bool isClosed() const {return ( (!isDirected(0)) && (!isDirected(1)) );}
    
    bool isHighwayEdge() const { assert( false ); return false; }
    void setHighwayEdge() { assert(false);}
    
    void makeTwoWay() {
        setFlag(direction(0));
        setFlag(direction(1));
    }
    
    void makeOneWay(int searchID) {
        setFlag(direction(searchID));
        unsetFlag(direction(1-searchID));
    }

    void makeClosed() {
        unsetFlag(direction(0));
        unsetFlag(direction(1));
    }

    /**
    * Ensures that this edge is at least open in the same
    * directions as the given edge.
    */
    void makeOpen(const EdgeCH& e) {
        if (e.isDirected(0)) setFlag(direction(0));
        if (e.isDirected(1)) setFlag(direction(1));
    }

    /**
    * Ensures that this edge is at least open in the opposite
    * directions as the given edge.
    */
    void makeOpenReversed(const EdgeCH& e) {
        if (e.isDirected(0)) setFlag(direction(1));
        if (e.isDirected(1)) setFlag(direction(0));
    }

    /**
    * Ensures that this edge is not open for any direction
    * that the given edge is open.
    */
    void makeClosed(const EdgeCH& e) {
        if (e.isDirected(0)) unsetFlag(direction(0));
        if (e.isDirected(1)) unsetFlag(direction(1));
    }

    /** Serializes this edge to a given stream. */
    void serialize(ostream& out) const {
        // not tested
        out.write((char*)this,sizeof(EdgeCH)/sizeof(char));
    }

    /** Deserializes this edge from a given stream. */
    void deserialize(istream& in) {
        // not tested
        in.read((char*)this,sizeof(EdgeCH)/sizeof(char));
    }

    void setWeight(EdgeWeight w) {
        _weight = w;
    }

    NodeID type() const {return _type;}
    
    void setType(NodeID k) {
        _type = k;
    }
    
    bool isShortcut() const { return type() != EDGE_TYPE_ORIGINAL; }

    bool isIdentical(const EdgeCH& e) const {
        return ((*this == e) && (weight() == e.weight()) && (type() == e.type()));
    }

    bool isReverse(NodeID u, const EdgeCH& e) const {
        return ((target() == u) && (weight() == e.weight()) &&
                (isDirected(0) == e.isDirected(1)) && (isDirected(1) == e.isDirected(0)) &&
                (type() == e.type()));
    }
    
    void setShortcutMiddle(NodeID middle) { }
    NodeID shortcutMiddle() const { return SPECIAL_NODEID; }

    void setShortcutEdge1(OtherEdgeID e) { }
    OtherEdgeID shortcutEdge1() const { return SPECIAL_NODEID; }

    void setShortcutEdge2(OtherEdgeID e) { }
    OtherEdgeID shortcutEdge2() const { return SPECIAL_NODEID; }
    OtherEdgeID shortcutEdgeLimit() const { return SPECIAL_NODEID; }
    
    void copyShortcutInfo(const EdgeCH& e, bool reverse)
    {
        #ifdef COUNT_SHORTCUT_ORIGINAL_EDGES
        _shortcutOriginalEdgeCount = e.shortcutOriginalEdgeCount();           
        #endif
    }
    
    void setShortcutOriginalEdgeCount(const OtherEdgeID count)
    {
        #ifdef COUNT_SHORTCUT_ORIGINAL_EDGES
        _shortcutOriginalEdgeCount = count;           
        #endif
    }
    
    OtherEdgeID shortcutOriginalEdgeCount() const
    {
        #ifdef COUNT_SHORTCUT_ORIGINAL_EDGES
        return _shortcutOriginalEdgeCount;
        #else
        return SPECIAL_NODEID;
        #endif
    }
    
protected:
    // The colon behind the variable denotes the number of bytes
    // the variable uses in the main memory.
    NodeID _target:26;
    OtherEdgeID _shortcutEdge1:6;
    EdgeWeight _weight:28;
    NodeID _type:2;
    NodeID _flags:2;

    #ifdef COUNT_SHORTCUT_ORIGINAL_EDGES
    OtherEdgeID _shortcutOriginalEdgeCount;
    #endif

    // counting from LSB (1) to MSB (32):
    // bit 1 and 2
    NodeID direction(int searchID) const {
        assert( (searchID >= 0) && (searchID < 2) );
        return searchID+1;
    }
    
    bool isFlag(NodeID flag) const {return ((_flags & flag) != 0);}
    void setFlag(NodeID flag) {_flags |= flag;}
    void unsetFlag(NodeID flag) {_flags &= (~flag);}
    
    void setTarget(const NodeID t) {
        _target = t;
        assert( t == SPECIAL_NODEID || _target == t );
    }
        
};


/** 
 * Represents a weighted edge assuming the source node is known,
 * i.e. only the target node is stored.
 * Used for Contraction Hierarchies.
 *
 * Along with the information normally stored with an EdgeCH
 * three values are stored to expand shortcut edges during
 * a recurisve unpacking routine. Remember that each shortcut
 * represents a path consiting of exactly two edges. To reconstruct
 * this path, we need the middle node of this path (shortcutMiddle).
 * To speedup the unpacking, we also store the relative indices
 * to the two edges of the path into the edge array of the middle node.
 * See DikjstraCH::pathTo() for more details to path expansion.
 */
class EdgeCHExpand : public EdgeCH
{
public:
    /** Default constructor. target and weight are set to 0.*/
    EdgeCHExpand()
        : EdgeCH() {
        _shortcutMiddle = 0;
        _shortcutEdge1 = 0;
        _shortcutEdge2 = 0;
		_bMarked = true;
		_iEdgeBound = 0;
    }
    
    /** 
     * Constructor.
     * @param target the target node
     * @param weight the weight
     * @param type edge type (see constants EDGE_TYPE_...)
     * @param forward from source to target open
     * @param backward from target to source open
     * @param shortcutMiddle middle nod of the path a shortcut represents
     * @param shortcutEdge1 relative index into the edge array of the first edge of the path
     * @param shortcutEdge2 relative index into the edge array of the second edge of the path
     * @param shortcutEdgeOriginalEdgeCount used for a priority term during CH node ordering
     */
    EdgeCHExpand(NodeID target, EdgeWeight weight, NodeID type, bool forward, bool backward, 
        NodeID shortcutMiddle = SPECIAL_NODEID, OtherEdgeID shortcutEdge1 = SPECIAL_NODEID, 
        OtherEdgeID shortcutEdge2 = SPECIAL_NODEID, OtherEdgeID shortcutOriginalEdgeCount = 1) :
        EdgeCH(target,weight,type,forward,backward) {
        setShortcutMiddle(shortcutMiddle);
        if (shortcutEdge1 < SHORTCUT_EDGE_LIMIT)
        {
            _shortcutEdge1 = shortcutEdge1;
        }
        else
        {
            _shortcutEdge1 = SHORTCUT_EDGE_LIMIT;
        }
        if (shortcutEdge2 < SHORTCUT_EDGE_LIMIT)
        {
            _shortcutEdge2 = shortcutEdge2;
        }
        else
        {
            _shortcutEdge2 = SHORTCUT_EDGE_LIMIT;
        }
		_bMarked = true;
		_iEdgeBound = 0;
        setShortcutOriginalEdgeCount(shortcutOriginalEdgeCount);
    }

    /** Constructs the reverse edge of the given edge. */
    EdgeCHExpand(NodeID target, const EdgeCHExpand& reverseEdge) :
        EdgeCH(target, reverseEdge) {
        copyShortcutInfo(reverseEdge, true);
		_bMarked = reverseEdge._bMarked;
		_iEdgeBound = reverseEdge._iEdgeBound;
    }

    // ***
    // To expand a path with a recursive unpacking routine,
    // three values are used, the middle node (shortcutMiddle)
    // and two relative indices into the edge array of the middle
    // node (shortcutEdge1, shortcutEdge2) that form the path
    // the shortcut represents.
    // ***

    void setShortcutMiddle(NodeID middle) { 
        _shortcutMiddle = middle; 
        assert( middle == SPECIAL_NODEID || _shortcutMiddle == middle );
    }
    NodeID shortcutMiddle() const { return _shortcutMiddle; }
    
    void setShortcutEdge1(OtherEdgeID shortcutEdge1) 
    {
        if (shortcutEdge1 < SHORTCUT_EDGE_LIMIT)
        {
            _shortcutEdge1 = shortcutEdge1;
        }
        else
        {
            _shortcutEdge1 = SHORTCUT_EDGE_LIMIT;
        }
    }
    OtherEdgeID shortcutEdge1() const { return _shortcutEdge1; }

    void setShortcutEdge2(OtherEdgeID shortcutEdge2) 
    { 
        if (shortcutEdge2 < SHORTCUT_EDGE_LIMIT)
        {
            _shortcutEdge2 = shortcutEdge2;
        }
        else
        {
            _shortcutEdge2 = SHORTCUT_EDGE_LIMIT;
        }
    }
    OtherEdgeID shortcutEdge2() const { return _shortcutEdge2; }
    
    OtherEdgeID shortcutEdgeLimit() const { return SHORTCUT_EDGE_LIMIT; }

    void copyShortcutInfo(const EdgeCHExpand& e, bool reverse) 
    { 
        setShortcutMiddle(e.shortcutMiddle());
        if ( true || !reverse )
        {
            _shortcutEdge1 = e.shortcutEdge1();
            _shortcutEdge2 = e.shortcutEdge2();
        }
        else
        {
            _shortcutEdge1 = e.shortcutEdge2();
            _shortcutEdge2 = e.shortcutEdge1();
        }
        #ifdef COUNT_SHORTCUT_ORIGINAL_EDGES
        _shortcutOriginalEdgeCount = e.shortcutOriginalEdgeCount();           
        #endif
		_bMarked = e._bMarked;
		_iEdgeBound = e._iEdgeBound;
		_iReserveOtherEdgeID = e._iReserveOtherEdgeID;
    }
    
    /** Serializes this edge to a given stream. */
    void serialize(ostream& out) const {
        // not tested
        out.write((char*)this,16/sizeof(char));
    }

    /** Deserializes this edge from a given stream. */
    void deserialize(istream& in) {
        // not tested
        in.read((char*)this,16/sizeof(char));
		_bMarked = true;
		_iEdgeBound = 0;
    }
    
protected:
    // The colon behind the variable denotes the number of bytes
    // the variable uses in the main memory.
    OtherEdgeID _shortcutEdge2:6;
    NodeID _shortcutMiddle:26;
    static const NodeID SHORTCUT_EDGE_LIMIT = (1 << 6) - 1;

public:
	bool _bMarked;
	EdgeWeight _iEdgeBound;
	OtherEdgeID _iReserveOtherEdgeID;
};

// Here, we decide which edge to use.
// Globally declaring the Edge class is unsatisfying and
// legacy of Dominik Schultes code. 
#ifdef USE_CH_EXPAND
typedef EdgeCHExpand Edge;
#else
typedef EdgeCHFast Edge;
#endif


/** Represents a weighted edge including the source node. */
class CompleteEdge : public Edge
{
    /** Outputs data about the edge for debugging purposes. */
    friend std::ostream& operator<<( std::ostream& os, const CompleteEdge &edge ) {
        os << "(" << edge.source();
        
        if (edge.isDirected(1)) os << "<";
        os << "-";
        if (edge.isDirected(0)) os << ">";
        
        os << edge.target() << ", " << edge.weight() << ")";
        return os;
    }
    
public:
    /** 
     * Compares two edges first by source, in case of equality by target, and
     * in case of equality by weight.
     * If all these values are identical, returns true iff this edge is two-way
     * and the given edge is at most one-way.
     */
    bool operator< (const CompleteEdge& e) const {
        if (source() == e.source()) {
            if (target() == e.target()) {
                if (weight() == e.weight()) {
                    return (isDirected(0) && isDirected(1) &&
                            ((! e.isDirected(0)) || (! e.isDirected(1))));
                }
                return (weight() < e.weight());
            }
            return (target() < e.target());
        }
        return (source() < e.source());
    }

    /** Constructor. */
    CompleteEdge() : Edge(), _source(0) {}
    
    /** Constructor. */
    CompleteEdge(NodeID s, NodeID t, EdgeWeight w, bool shortcut, bool forward, bool backward) : 
            Edge(t, w, shortcut, forward, backward), _source(s) {}
    
    /** Constructor. */
    CompleteEdge(const Edge& e, NodeID s) : Edge(e), _source(s) {}

    NodeID source() const {return _source;}

    void serialize(ostream& out) const {
        Edge::serialize(out);
        writePrimitive(out, _source);
    }

    void deserialize(istream& in) {
        Edge::deserialize(in);
        readPrimitive(in, _source);
    }

private:
    NodeID _source;
};

/** Represents an edge in the search space. */
class SearchSpaceEdge
{
public:
    SearchSpaceEdge()
    : _source(0), _target(0), _edgeID(0), _keyOfTarget(0), _searchLevel(0), _core(true) {}
    
    SearchSpaceEdge(NodeID s, NodeID t, OtherEdgeID e, EdgeWeight k, LevelID lev = 0, bool c = true)
    : _source(s), _target(t), _edgeID(e), _keyOfTarget(k), _searchLevel(lev), _core(c) {}

    void overwriteSearchLevel(LevelID lev) {_searchLevel = lev;}
    
    NodeID source() const {return _source;}
    NodeID target() const {return _target;}
    OtherEdgeID edgeID() const {return _edgeID;}
    EdgeWeight keyOfTarget() const {return _keyOfTarget;}
    LevelID searchLevel() const {return _searchLevel;}
    bool core() const {return _core;}
    
private:
    NodeID _source;
    NodeID _target;
    OtherEdgeID _edgeID;
    EdgeWeight _keyOfTarget;
    LevelID _searchLevel;
    bool _core;
};

/** Represents an easy edge used in TNR, SILC  and PCP */
/*
	An edge contains a source node and a target node
 */
class EasyEdge{
		friend ostream& operator <<(ostream& os, const EasyEdge &e)
		{
			os << " ( " << e._source << " , " << e._target << " ) " << endl;
			return os;
		}
	public:
		EasyEdge(){

		}
		EasyEdge(NodeID s, NodeID t):_source(s), _target(t)
		{

		}

		NodeID source() const
		{
			return _source;
		}
		NodeID target() const
		{
			return _target;
		}

		void setSource(NodeID s){
			_source = s;
		}
		void setTarget(NodeID t){
			_target = t;
		}
		/** 
		 * Compares two edges first by source, in case of equality by target
		 */
		bool operator< (const EasyEdge& e) const {
			if (source() == e.source()) {
				return (target() < e.target());
			}
			return (source() < e.source());
		}
		private:
			NodeID _source;
			NodeID _target;

	};
#endif 
