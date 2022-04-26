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

#ifndef MANYTOMANY_H
#define MANYTOMANY_H

#include "../datastr/graph/graph.h"
#include "../datastr/graph/SearchGraph.h"

typedef datastr::graph::SearchGraph TransitGraph;

#include "../stats/utils.h"
#include "../processing/DijkstraCH.h"
#include "../datastr/searchSpaces.h"

#include <unordered_set>

/**
 * Provides methods to perform many-to-many computations.
 */
template <typename Graph, typename DijkstraFW, typename DijkstraBW, bool performBucketScans = true>
class ManyToMany : public SearchSpaces
{
public:

	const static int MAX_TREE_NODE_NUM = 10000;
    /** Constructor. */
    ManyToMany(Graph *const g, const LevelID earlyStopLevel) : _g(g), _dFW(g), _dBW(g), _earlyStopLevel(earlyStopLevel) {
        // early stop level currently not supported
        //_dFW.setEarlyStopLevel(earlyStopLevel);
        //_dBW.setEarlyStopLevel(earlyStopLevel);
    }
	ManyToMany()
	{

	}
	void setManyToMany(Graph *const g, const LevelID earlyStopLevel)
	{
		_g = g;
		_dFW.setDij(g);
		_dBW.setDij(g);
		_earlyStopLevel = earlyStopLevel;
	}

	
	//
    /**
     * Reference implementation.
     * Performs point-to-point queries in order to fill the distance table.
     */
    void computeMatrixNaive(const vector<NodeID>& sources, const vector<NodeID>& targets, Matrix<EdgeWeight>& result) {
        assert( sources.size() == result.noOfRows() );
        assert( targets.size() == result.noOfCols() );

        VERBOSE( cout << "computing reference solution ..." << endl );

        VERBOSE( Percent progress(sources.size()) );
        for (NodeID u = 0; u < sources.size(); u++) {
            VERBOSE( progress.printStatus(u) );
            for (NodeID v = 0; v < targets.size(); v++) {
                const NodeID s = sources[u];
                const NodeID t = targets[v];
                const EdgeWeight w = _dFW.bidirSearch(s, t);
                _dFW.clear();
                result.set(u, v, w);
            }
        }

        VERBOSE( cout << "done." << endl );
    }

   
	/**
     * Efficient implementation.
     * Applies many-to-many algorithm in order to fill the distance table.
     */
    void computeMatrix(const vector<NodeID>& sources, const vector<NodeID>& targets, Matrix<EdgeWeight>& result) {

      // VERBOSE( cout << "computing " << sources.size() << " x " << targets.size() << " table ..." << endl );

        const double start = timestamp();

        // init table
        result.init(Weight::MAX_VALUE);

        // backward search
#ifdef DEBUG
        VERBOSE( cout << "backward search" << endl );
        VERBOSE( Percent progress(targets.size()) );
#endif
        for (NodeID v = 0; v < targets.size(); v++) {
#ifdef DEBUG
            VERBOSE( progress.printStatus(v) );
#endif
            const NodeID t = targets[v];
            setCurrentNode(v);
            _dBW.bidirSearch(SPECIAL_NODEID, t);
            _dBW.obtainRelevantSearchSpace(*this);
            _dBW.clear();
        }
        double elapsedTime = timestamp() - start;
#ifdef DEBUG
        VERBOSE( cout << elapsedTime << " ms" << endl );       
        VERBOSE( cout << "  search space sizes: int " << _searchSpacesBwDynInt.size()
                      << " + short " << _searchSpacesBwDynShort.size() << endl );
#endif

        // only the 'int'-search space is used
        assert( _searchSpacesBwDynShort.size() == 0 );
        
        // sort and copy
#ifdef DEBUG
        VERBOSE( cout << "sort and copy" << endl );
#endif
        _searchSpacesBwDynInt.sort(_g->noOfNodes());
        _searchSpacesBwInt = _searchSpacesBwDynInt;
        elapsedTime = timestamp() - start;
#ifdef DEBUG
        VERBOSE( cout << elapsedTime << " ms" << endl );
        // forward search
        VERBOSE( cout << "forward search" << endl );
        VERBOSE( progress.reinit(sources.size()) );
#endif
        COUNTING( unsigned long long bucketScans = 0 );
	    //COUNTING( unsigned long long bucketScansTop = 0 );
        for (NodeID u = 0; u < sources.size(); u++) {
#ifdef DEBUG
            VERBOSE( progress.printStatus(u) );
#endif
            _dFW.bidirSearch(sources[u], SPECIAL_NODEID);
            _dFW.obtainRelevantSearchSpace(*this);
            _dFW.clear();

            if (performBucketScans) {
                const NodeID matrixIndexOffset = result.index(u, 0);
                for (NodeID i = 0; i < _searchSpaceFW.size(); i++) {
                    const NodeID via = _searchSpaceFW[i].target();
                    const EdgeWeight distFW = _searchSpaceFW[i].weight();
                
                    const ISSInt::const_iterator endInt = _searchSpacesBwInt.end(via);
                    for (ISSInt::const_iterator it = _searchSpacesBwInt.begin(via); it != endInt; it++) {
						//cout << it->origin() << " ";
                        result.improve(matrixIndexOffset + it->origin(), distFW + it->dist());
                        COUNTING( bucketScans++ );
			            //COUNTING( if (_g->node(via).level() >= _earlyStopLevel) bucketScansTop++ );
                    }
					//cout << endl;
                }
            }
            
            _searchSpaceFW.clear();
        }
        elapsedTime = timestamp() - start;
#ifdef DEBUG
        VERBOSE( cout << elapsedTime << " ms" << endl; )
#endif 
	    //COUNTING( cout << "bucket scans (top): " << bucketScansTop << endl );
    }


	
	/**
     * Efficient implementation.
     * Applies many-to-many algorithm in order to fill the distance table, also stores the path.
	 * by dingxiong based on computeMatrix
     */
    void computeMatrixStorePath(const vector<NodeID>& sources, const vector<NodeID>& targets, Matrix<EdgeWeight>& result) {

        //VERBOSE( cout << "computing " << sources.size() << " x " << targets.size() << " table ..." << endl );
        const double start = timestamp();

        // init table
        result.init(Weight::MAX_VALUE);
		
		//init map
		for (NodeID i = 0; i < sources.size(); i++)
			_sourcesMapToID.insert(make_pair(sources[i], i));
		for (NodeID i = 0; i < targets.size(); i++)
			_targetsMapToID.insert(make_pair(targets[i], i));

		//init rooted tree and the intermediate node 
		_viaNode.resize(sources.size());
		for (NodeID i = 0; i < sources.size(); i++)
			_viaNode[i].resize(targets.size());
		_rootedTreeBW.resize(targets.size());
		_rootedTreeFW.resize(sources.size());
		
		//_rootedTreeBW = new DynQueryPQueue[targets.size()];
		//_rootedTreeFW = new DynQueryPQueue[sources.size()];
		//cout << _rootedTreeBW[0]._hSize << " " <<_rootedTreeBW[0]._eSize << endl;

        // backward search, and also store every backward tree
       // VERBOSE( cout << "backward search" << endl );
       // VERBOSE( Percent progress(targets.size()) );
        for (NodeID v = 0; v < targets.size(); v++) {
           // VERBOSE( progress.printStatus(v) );
            const NodeID t = targets[v];
            setCurrentNode(v);
			//DijkstraBW *dBW = new DijkstraBW(_g);
            //_dBW.bidirSearch(SPECIAL_NODEID, t);

			//revised by ddx
			_dBW.bidirSearchWithTree(SPECIAL_NODEID, t, 1, _rootedTreeBW[v]);
			/*cout << endl;
			cout << "copy heap: "<< _rootedTreeBW[v].elementsSize() << endl;
			for (NodeID i = 0; i < _rootedTreeBW[v].elementsSize(); i++){
				cout << i << " " << _rootedTreeBW[v].getElements(i).data().nodeID() << " ";
			}*/
			//revised end
			_dBW.obtainRelevantSearchSpace(*this);
            _dBW.clear();
        }
        double elapsedTime = timestamp() - start;
       
       // VERBOSE( cout << elapsedTime << " ms" << endl );
       

       // VERBOSE( cout << "  search space sizes: int " << _searchSpacesBwDynInt.size()
        //              << " + short " << _searchSpacesBwDynShort.size() << endl );

        // only the 'int'-search space is used
        assert( _searchSpacesBwDynShort.size() == 0 );
        
        // sort and copy
       // VERBOSE( cout << "sort and copy" << endl );
        _searchSpacesBwDynInt.sort(_g->noOfNodes());
        _searchSpacesBwInt = _searchSpacesBwDynInt;
        elapsedTime = timestamp() - start;
     //   VERBOSE( cout << elapsedTime << " ms" << endl );

        // forward search
    //    VERBOSE( cout << "forward search" << endl );
    //    VERBOSE( progress.reinit(sources.size()) );
        COUNTING( unsigned long long bucketScans = 0 );
	    //COUNTING( unsigned long long bucketScansTop = 0 );
        for (NodeID u = 0; u < sources.size(); u++) {
          //  VERBOSE( progress.printStatus(u) );
			//DijkstraFW *dFW = new DijkstraFW(_g);
			// after new a DijkstraFW, remember to delete it

            //_dFW.bidirSearch(sources[u], SPECIAL_NODEID);
			/*revised by ddx*/
			_dFW.bidirSearchWithTree(sources[u], SPECIAL_NODEID, 0, _rootedTreeFW[u]);
			/*revised end*/
			
			_dFW.obtainRelevantSearchSpace(*this);
            _dFW.clear();

            if (performBucketScans) {
                const NodeID matrixIndexOffset = result.index(u, 0);
                for (NodeID i = 0; i < _searchSpaceFW.size(); i++) {
                    const NodeID via = _searchSpaceFW[i].target();
                    const EdgeWeight distFW = _searchSpaceFW[i].weight();
                
                    const ISSInt::const_iterator endInt = _searchSpacesBwInt.end(via);
                    for (ISSInt::const_iterator it = _searchSpacesBwInt.begin(via); it != endInt; it++) {
						if (distFW + it->dist() < result.value(u, it->origin())){
							_viaNode[u][it->origin()] = via;
							result.improveNojudge(u, it->origin(), distFW + it->dist());
						}
                        COUNTING( bucketScans++ );
			            //COUNTING( if (_g->node(via).level() >= _earlyStopLevel) bucketScansTop++ );
                    }
                }
            }
            
            _searchSpaceFW.clear();
        }
		/*for (NodeID i = 0; i < sources.size(); i++){
			for (NodeID j = 0; j < targets.size(); j++){
				cout << _viaNode[i][j] << " ";
			}
			cout << endl;
		}*/
        elapsedTime = timestamp() - start;
       // VERBOSE( cout << elapsedTime << " ms" << endl; )
			// build the map from the pqQueue
		buildPQIndex(sources, targets);
	    //COUNTING( cout << "bucket scans (top): " << bucketScansTop << endl );
    }
    
	void buildPQIndex(const vector <NodeID> &sources, const vector <NodeID> &targets)
	{
		_targetsPQIndex.resize(targets.size());
		_sourcesPQIndex.resize(sources.size());
		for (NodeID i = 0; i < targets.size(); i++){
			//VERBOSE(cout << _rootedTreeBW[i].elementsSize() << endl;)
			for (NodeID j = 0; j < _rootedTreeBW[i].elementsSize(); j++){
				if (!_rootedTreeBW[i].elements()[j].isDummy()){
					//cout << _rootedTreeBW[i].getElements(j).data().nodeID() << endl;
					_targetsPQIndex[i].insert(make_pair(_rootedTreeBW[i].elements()[j].data().nodeID(), j));
				}	
			}
		}
		for (NodeID i = 0; i < sources.size(); i++){
			//VERBOSE(cout << _rootedTreeFW[i].elements().size() << endl;)
			for (NodeID j = 0; j < _rootedTreeFW[i].elementsSize(); j++){
				if (!_rootedTreeFW[i].elements()[j].isDummy()){
					_sourcesPQIndex[i].insert(make_pair(_rootedTreeFW[i].elements()[j].data().nodeID(), j));
				}	
			}
		}

	}
	
	/**
	 * find the common vertices in one to many algorithms
	 * by dingxiong
	 */
	bool findOneToMany_CommonVertices(const vector <NodeID> &sources, const vector <NodeID> &targets, 
		unordered_set <NodeID> &prevVertexSet, unordered_set <NodeID> &curVertexSet)
	{
		NodeID size = targets.size();
		NodeID bef, nex;
		NodeID i;
		Path path;
		unordered_set <NodeID>::iterator nodeit;
		bool expand = true;
		for (i = 0; i < size; i++){
			pathToInMany(path, targets[i], _viaNode[0][i], 1, false, expand);
			EasyPath epath(path);		
			bef = epath.node(0);
			if (prevVertexSet.empty()){
				curVertexSet.insert(bef);
			}
			else{
				nodeit = prevVertexSet.find(bef);
				if (nodeit != prevVertexSet.end())
					curVertexSet.insert(bef);
			}
			for (NodeID i = 1; i < epath.noOfNodes(); i++){
				nex = epath.node(i);
				if (prevVertexSet.empty()){
					curVertexSet.insert(nex);
				}
				else{
					nodeit = prevVertexSet.find(nex);
					if (nodeit != prevVertexSet.end())
						curVertexSet.insert(nex);
				}
			}
		}
		for (i = 0; i < size; i++){
			pathToInMany(path, sources[0], _viaNode[i], 0, true, expand);
			EasyPath epath(path);		
			bef = epath.node(0);
			if (prevVertexSet.empty()){
				curVertexSet.insert(bef);
			}
			else{
				nodeit = prevVertexSet.find(bef);
				if (nodeit != prevVertexSet.end())
					curVertexSet.insert(bef);
			}
			for (NodeID i = 1; i < epath.noOfNodes(); i++){
				nex = epath.node(i);
				if (prevVertexSet.empty()){
					curVertexSet.insert(nex);
				}
				else{
					nodeit = prevVertexSet.find(nex);
					if (nodeit != prevVertexSet.end())
						curVertexSet.insert(nex);
				}
			}
			if (curVertexSet.empty() == true)
				return false;
		}
		return true;
	}

  
	/**
    * Determines the shortest path from the source node s of the search 
    * that has been performed last to the given node t.
    * @param path a reference to the path object that accepts the result
    * @param t the target node (not relevant if searchID = -1)
    * @param searchID the search direction that should be considered; special value -1 is used
    *                 to indicate that the result of a bidirectional search should be considered
    * @param forward true/false = return the path from s to t / from t to s
    * @param expand flag: expand the contracted path (i.e., return the expanded path)
    *                 This is a recursive unpacking routine that replaces each shortcut
    *                 by its two originating edges.
    * @return a reference to the given path object that accepts the result
     */
    Path& pathToInMany(Path& path, NodeID source, NodeID target, int searchID = 0, bool forward = true, bool expand = false) {
        
		path.clear();
        if (searchID == -1) { // special case: return result of the bidirectional search
            if (source == target) {
                path.setNotConnected();
                return path;
            }

			NodeID sid = _sourcesMapToID[source];
			NodeID tid = _targetsMapToID[target];
			NodeID via = _viaNode[sid][tid];


            
			if ( expand )
            {
                // The query algorithm is bidirectional. We have two
                // partial paths, one of the forward search and one of the
                // backward search.
                // Get forward path reversed, not expanded because
                // this is the natural way using the parent pointers.
                Path pathForwardReverse;
                pathToInMany(pathForwardReverse, source, via, 0, false, false);
				//VERBOSE(cout << pathForwardReverse << endl;)
               
				// Expand forward path, this writes the expanded
                // path in a second object and implicitely reverses
                // the path.
                expandPath<true,false,true>(pathForwardReverse, path);
                assert( pathForwardReverse.length() == path.length() );
				//VERBOSE(cout << path << endl;)

                // The backward path will be processed in code below
                // that is used for the unidirectional case.
                // We just need to set these variables accordingly.
                searchID = 1;
				source = target;
                target = via;
                forward = !forward;

            }
            else
            {

                // first part of the path from s to the intermediate node 'via'
                pathToInMany(path, source, via, 0, true, expand);

                // second part of the path from the intermediate node 'via' to t
                Path path2;
                pathToInMany(path2, target, via, 1, false, expand);

                // combine both parts
                path.add(path2);

                // We have the path forwards.
                if (! forward) path.reverse(); // reverse, if required

                return path;
            }

        }
        else
        {
            // add first node to path, all following nodes are added
            // with an edge
            path.addFirstNode(target);
        }

        if ( expand )
        {
            // simple expand routine, because shortcut 
            // from parent->child is original parent->middle->child
            
            // If we found a shortcut, we replace it by the two originating
            // edges. One edge is put on a stack and the second one is processed
            // in the next pass of the loop. _expandPathStack is this stack.
            assert( _expandPathStack.empty() );
            NodeID parent = SPECIAL_NODEID;
            OtherEdgeID e = SPECIAL_NODEID;
            
            // If we found a shortcut, we replace it by the two originating
            // edges. One edge is put on a stack and the second one is processed
            // in the next pass of the loop. processSecondHalf indicates this.
            bool processSecondHalf = false;
			//cout << source << " " << target << endl; 
            
            // Main loop, will run until all edges are replaced by shorcut edges.
            while (true) 
            {
                // This variable determines which of the two edges of a shortcut is
                // accessed first in the edge array. Intuitively, the edge that
                // is processed first should be accessed first. But tests indicate
                // otherwise.
                //bool shortcutParent1Child2 = !forward;
                
                // To process an edge, we need the edge id of this edge
                // and the parent in the Dijkstra search graph to ensure
                // a correct order of the edges in the unpacked path.
                // The current position of the processing is marked by node t.
                // Either in the last pass of the loop an edge was specified
                // (processSecondHalf == true)...
                if ( !processSecondHalf )
                {
                    // ... or we take the next edge from the path found by the query algorithm if 
                    // the stack is empty ...
                    if ( _expandPathStack.empty() )
                    {
						NodeID tindex = _targetsMapToID[source];
						map <NodeID, NodeID>::iterator mit = _targetsPQIndex[tindex].find(target);
                        if ( mit == _targetsPQIndex[tindex].end() ) {
                            path.setNotConnected();
                            return path;
                        }
                        // determine parent
                        parent = _rootedTreeBW[tindex].elements()[mit->second].data().parentNode();
                        if (parent == target) break;
                        e = _rootedTreeBW[tindex].elements()[mit->second].data().parentEdge();
                    }
                    // ... otherwise we take the edge on top of the stack.
                    else
                    {
                        parent = _expandPathStack.top().first;
                        e = _expandPathStack.top().second;
                        _expandPathStack.pop();
						
                        //shortcutParent1Child2 = !shortcutParent1Child2;
                    }
                }
				
                // process edge
                expandEdgeMany(path, target, parent, e, processSecondHalf, forward);
				
				
            }
            assert( _expandPathStack.empty() );
        }
        
        // No path expansion requested, only return path found by query algorithm.
        else
        {
            // determine path backwards, starting with t
            while (true) {
				NodeID sindex = _sourcesMapToID[source];
				//DynQueryPQueue dq = _rootedTreeFW[sindex];
				map <NodeID, NodeID>::iterator mit = _sourcesPQIndex[sindex].find(target);
                if (mit == _sourcesPQIndex[sindex].end()) {
                    path.setNotConnected();
                    return path;
                }
                // determine parent
                NodeID parent = _rootedTreeFW[sindex].elements()[mit->second].data().parentNode();
                if (parent == target) break;
                target = parent;
                OtherEdgeID e = _rootedTreeFW[sindex].elements()[mit->second].data().parentEdge();

                path.add( target, _g->edge(e).weight(), e, _g->edge(e).isShortcut() );
            }
        }

        // We have the path backwards.
        // reverse, if required
        // (Note that the expander cannot deal with a backward path since the
        //  included edge IDs would not correspond to this direction.)
        if (forward) path.reverse();
        return path;
    }

	
	/**
     * Expands a given path. This is almost the same path expansion routine as in pathTo() but
     * reads the path that should be expanded from a Path object instead of the Dijkstra search tree.
     * This is used for the forward search since an extracted path would be read from the Dijkstra
     * search three reversed, meaning from the target (or via node) to the source node. We do not want
     * to reverse the expanded path (100 times longer than the path found by the query algorithm). 
     * So we store the reversed path found by the query algorithm, process it backwards to get an expanded
     * path forwards.
     * @param reverse should be expanded path be reversed to the input path
     * @param forward is the input path forward from source node to target node (or via node)
     * @param newPath should the pathExpanded return path be cleared and newly initalized
     * @param path a reference to the original path object 
     * @param pathExpanded a reference to the path object that accepts the expanded result
     */
    template<bool reverse, bool forward, bool newPath>
    void expandPath(const Path& path, Path& pathExpanded) 
    {
        if (newPath) pathExpanded.clear();
        // simple expand routine, because shortcut from parent->child is original parent->middle->child
        // assertion: edge-array is sorted by target
        if (path.isNotConnected())
        {
            pathExpanded.setNotConnected();
            return;
        }

        assert( !path.empty() );

        // The current position of the processing is marked by node t.
        // The current index in the input path is variable index.
        NodeID t;
        NodeID index;
        if ( !reverse )
        {
            t = path.firstNode();
            index = 0;
        }
        {
            t = path.lastNode();
            index = path.noOfEdges() - 1;
        }
        if (newPath) pathExpanded.addFirstNode(t);


        // If we found a shortcut, we replace it by the two originating
        // edges. One edge is put on a stack and the second one is processed
        // in the next pass of the loop. _expandPathStack is this stack.
        assert( _expandPathStack.empty() );
        NodeID parent = SPECIAL_NODEID;
        OtherEdgeID e = SPECIAL_NODEID;
        
        
        
        // If we found a shortcut, we replace it by the two originating
        // edges. One edge is put on a stack and the second one is processed
        // in the next pass of the loop. processSecondHalf indicates this.
        bool processSecondHalf = false;
        
        // Main loop, will run until all edges are processed.
        // These can be on the stack or on the input path.
        while (!(_expandPathStack.empty() &&
                (( reverse && (index == SPECIAL_NODEID) )
                || ( !reverse && index >= path.noOfEdges() )))) 
        {
            // This variable determines which of the two edges of a shortcut is
            // accessed first in the edge array. Intuitively, the edge that
            // is processed first should be accessed first. But tests indicate
            // otherwise.
            //bool shortcutParent1Child2 = !forward;

            // To process an edge, we need the edge id of this edge
            // and the parent in the Dijkstra search graph to ensure
            // a correct order of the edges in the unpacked path.
            // The current position of the processing is marked by node t.
            // Either in the last pass of the loop an edge was specified
            // (processSecondHalf == true)...
            if ( !processSecondHalf )
            {
                // ... or we take the next edge from the path found by the query algorithm if 
                // the stack is empty ...
                if ( _expandPathStack.empty() )
                {
                    assert( index != SPECIAL_NODEID );
                    if ( !reverse )
                    {
                        e = path.edge(index);
                        index++;
                        parent = path.node(index);
                    }
                    else
                    {
                        e = path.edge(index);
                        parent = path.node(index);
                        if ( index > 0 ) index--;
                        else index = SPECIAL_NODEID;
                        //shortcutParent1Child2 = !shortcutParent1Child2;
                    }
                }
                // ... otherwise we take the edge on top of the stack.
                else
                {
                    parent = _expandPathStack.top().first;
                    e = _expandPathStack.top().second;
                    _expandPathStack.pop();
                    //shortcutParent1Child2 = !shortcutParent1Child2;
                }
            }
            
            // process edge
            expandEdgeMany(pathExpanded, t, parent, e, processSecondHalf, forward );
        }
        assert( _expandPathStack.empty() );
    }
	
	/**
     * This subroutine is common to pathTo() and expandPath(), give an edge by its edge id e,
     * it is either a shortcut and gets expanded or not and will be appended to path.
     * A shortcut is expanded into two edges, the first is put on the stack, the second
     * is prepared to be processed next (processSecondHalf flag).
     * The variables have the same name as the ones used in pathTo() and expandPath().
     */
    void expandEdgeMany(Path& path, NodeID& t, NodeID& parent, OtherEdgeID& e, bool& processSecondHalf, const bool forward)
    {
        assert( parent != SPECIAL_NODEID );
        assert( e != SPECIAL_NODEID );
        const Edge& edge = _g->edge(e);

        // A shortcut parent->child is split into two edges parent->middle->child
        if ( edge.isShortcut() )
        {
            NodeID middle = edge.shortcutMiddle();
            OtherEdgeID firstEdge = _g->firstLevelEdge(middle);
            
            // edge.shortcutEdge1() and edge.shortcutEdge2() are a relative index
            // into the adjacency array of the middle node. They are valid if
            // they are != edge.shortcutEdgeLimit()
            OtherEdgeID eParent, eChild;
            if ( true /*|| shortcutParent1Child2*/ )
            {
                eParent = edge.shortcutEdge1();
                eChild = edge.shortcutEdge2();
            }
            else
            {
                eParent = edge.shortcutEdge2();
                eChild = edge.shortcutEdge1();
            }

            // eParent and eChild are valid relative indices
            if ( eParent != edge.shortcutEdgeLimit() && eChild != edge.shortcutEdgeLimit() )
            {
                eParent += firstEdge;
                eChild  += firstEdge;
                // Switch child <-> parent if required, the edges in the
                // expanded path need to be in the right order.
                if ( _g->edge(eChild).target() == parent )
                {
                    OtherEdgeID dummy = eParent;
                    eParent = eChild;
                    eChild = dummy;   
                }
            }

            // at least on of eParent or eChild is invalid, need to scan through the OtherEdgeID
            // array of the middle node to find the required edges. This can happen if there 
            // are more than shortcutEdgeLimit() edges incident to a node.
            else 
            {
                if ( eParent == edge.shortcutEdgeLimit() ) 
                {
                    eParent = SPECIAL_NODEID;
                    if ( eChild == edge.shortcutEdgeLimit() ) eChild = SPECIAL_NODEID;
                    else
                    {
                        eChild += firstEdge;
                        // Switch child <-> parent if required, the edges in the
                        // expanded path need to be in the right order.
                        if ( _g->edge(eChild).target() == parent )
                        {
                            eParent = eChild;
                            eChild = SPECIAL_NODEID;
                        }
                        // Relative index to the child node was wrong.
                        else if ( _g->edge(eChild).target() != t )
                        {
                            eChild = SPECIAL_NODEID;
                        }
                    }

                }
                else 
                {
                    eParent += firstEdge;
                    // Switch child <-> parent if required, the edges in the
                    // expanded path need to be in the right order.
                    if ( _g->edge(eParent).target() == t )
                    {
                        eChild = eParent;
                        eParent = SPECIAL_NODEID;
                    }
                    // Relative index to the parent node was wrong.
                    else if ( _g->edge(eParent).target() != parent )
                    {
                        eParent = SPECIAL_NODEID;
                    }
                }     

                NodeID directionParent = forward ? 1 : 0;
                NodeID directionChild  = 1-directionParent;
                OtherEdgeID lastEdge = _g->lastEdge(middle);
                // scan through all edges, starting at shortcutEdgeLimit since otherwise
                // the relative index would be valid.
                for ( OtherEdgeID eMiddle = firstEdge + edge.shortcutEdgeLimit(); eMiddle < lastEdge; eMiddle++ )
                {
                    const Edge& edgeMiddle = _g->edge(eMiddle);
                    if ( edgeMiddle.isDirected( directionParent ) && edgeMiddle.target() == parent )
                    {
                        eParent = eMiddle;
                        if ( eChild != SPECIAL_NODEID ) break;
                    }
                    if ( edgeMiddle.isDirected( directionChild ) && edgeMiddle.target() == t )
                    {
                        eChild = eMiddle;
                        if ( eParent != SPECIAL_NODEID ) break;
                    }
                }
            }
            //cout << middle << " " << parent << endl;
            // The first of the two edges is put on the stack ...
            _expandPathStack.push( make_pair( parent, eParent ) );
            // .. the other edge one is process in the next pass of the loop. 
            processSecondHalf = true;
            parent = middle;
            e = eChild;
        }
        
        // Edge was no shortcut, can be added to the expanded path
        else
        {
            processSecondHalf = false;
            t = parent;
            path.add( t, edge.weight(), e, edge.isShortcut() );
        }
    }

	
	void insert(NodeID *&arrays, int pos, int v, int length){
		int k=0;
		for(k=length;k>pos;k--)
			arrays[k]=arrays[k-1];
		arrays[pos]=v;
	}


	int binarysearch(NodeID *&arrays, int l, int h, NodeID v){
		int cur=l;
		int low=l;
		int high=h;
		while(low<high){
			cur=low+((high-low)/2);
			if(arrays[cur]==v)
				return cur;
			else{
				if(arrays[cur]>v)
					high=cur;
				else
					low=cur+1;
			}
		}
		return -low-1;
	}


	void clear()
	{
		for (NodeID i = 0; i < _viaNode.size(); i++)
			_viaNode[i].clear();
		_viaNode.clear();
		_rootedTreeBW.clear();
		_rootedTreeFW.clear();
			
		_sourcesMapToID.clear();
		_targetsMapToID.clear();

		for (NodeID i = 0; i < _sourcesPQIndex.size(); i++)
			_sourcesPQIndex[i].clear();
		_sourcesPQIndex.clear();

		for (NodeID i = 0; i < _sourcesPQIndex.size(); i++)
			_targetsPQIndex[i].clear();
		_targetsPQIndex.clear();
	}
	
	//used for level two transit node distance tables building
	//by dingxiong
	void computeMatrixForLocalTransitNode(const vector<NodeID>& sources, const vector<NodeID>& targets, int ascale, CoordinateType cellLength, CoordinateType range, 
		NodeID** &levelTwoLocalMapping, EdgeWeight** &levelTwoTransitDistTable, map <NodeID, NodeID> &reverseMapping, int localNum,
		vector <CoordinateType> &xcord, vector <CoordinateType> &ycord)
	{

		// VERBOSE( cout << "computing " << sources.size() << " x " << targets.size() << " table ..." << endl );

		const double start = timestamp();			

		// backward search

		for (NodeID v = 0; v < targets.size(); v++) {

			const NodeID t = targets[v];
			setCurrentNode(v);
			_dBW.bidirSearch(SPECIAL_NODEID, t);
			_dBW.obtainRelevantSearchSpace(*this);
			_dBW.clear();
		}
		double elapsedTime = timestamp() - start;

		// only the 'int'-search space is used
		assert( _searchSpacesBwDynShort.size() == 0 );

		// sort and copy

		_searchSpacesBwDynInt.sort(_g->noOfNodes());
		_searchSpacesBwInt = _searchSpacesBwDynInt;
		elapsedTime = timestamp() - start;

		COUNTING( unsigned long long bucketScans = 0 );
		
		//COUNTING( unsigned long long bucketScansTop = 0 );
		for (NodeID u = 0; u < sources.size(); u++) {
			Matrix <EdgeWeight> oneRow(1, targets.size());
			oneRow.init(Weight::MAX_VALUE);
			if ( u % 1000 == 0)
				cout << "locap mapping " << u << endl;
			_dFW.bidirSearch(sources[u], SPECIAL_NODEID);
			_dFW.obtainRelevantSearchSpace(*this);
			_dFW.clear();

			if (performBucketScans) {
				//const NodeID matrixIndexOffset = result.index(u, 0);
				for (NodeID i = 0; i < _searchSpaceFW.size(); i++) {
					const NodeID via = _searchSpaceFW[i].target();
					const EdgeWeight distFW = _searchSpaceFW[i].weight();

					const ISSInt::const_iterator endInt = _searchSpacesBwInt.end(via);
					for (ISSInt::const_iterator it = _searchSpacesBwInt.begin(via); it != endInt; it++) {
						oneRow.improve(it->origin(), distFW + it->dist());
						//result.improve(matrixIndexOffset + it->origin(), distFW + it->dist());
						COUNTING( bucketScans++ );	
					}				
				}
			}
			_searchSpaceFW.clear();
			
			int x1, x2, y1, y2;
			int absx, absy;
			// for local transit computing
			x1 = (NodeID) (xcord[reverseMapping[u]] / cellLength);	//Ӧ�ü���ţ�earnestwu
			y1 = (NodeID) (ycord[reverseMapping[u]] / cellLength);	//Ӧ�ü���ţ�earnestwu
			if (xcord[reverseMapping[u]] == range){
				x1 = ascale - 1;
			}
			if (ycord[reverseMapping[u]] == range)
				y1 = ascale - 1;
			for (NodeID j = 0; j < sources.size(); j++){
				x2 = (NodeID) (xcord[reverseMapping[j]] / cellLength);	//Ӧ�ü���ţ�earnestwu
				y2 = (NodeID) (ycord[reverseMapping[j]] / cellLength);	//Ӧ�ü���ţ�earnestwu
				if (xcord[reverseMapping[j]] == range){
					x2 = ascale - 1;
				}
				if (ycord[reverseMapping[j]] == range)
					y2 = ascale - 1;
				if ((x1 - x2) > 0)
					absx = x1 - x2;
				else
					absx = x2 - x1;
				if ((y1 - y2) > 0)
					absy = y1 - y2;
				else
					absy = y2 - y1;
			    
				if (!(absx > 5 || absy > 5)){
					
					int low = 2;
					int high = levelTwoLocalMapping[u][1];
					int pos = binarysearch(levelTwoLocalMapping[u], low, high, reverseMapping[j]);
					if (pos < 0){
						if (high >= (int)levelTwoLocalMapping[u][0]){
							levelTwoLocalMapping[u] = (NodeID *)realloc(levelTwoLocalMapping[u], (high + localNum) * sizeof (NodeID));
							levelTwoLocalMapping[u][0] += localNum;
							levelTwoTransitDistTable[u] = (EdgeWeight *)realloc(levelTwoTransitDistTable[u], (high + localNum) * sizeof(EdgeWeight));
							levelTwoTransitDistTable[u][0] += localNum;
						}
						insert(levelTwoLocalMapping[u], -pos-1, reverseMapping[j], high);
						if (u != j)
							insert(levelTwoTransitDistTable[u], -pos-1, oneRow.value(0, j), high);
						else
							insert(levelTwoTransitDistTable[u], -pos-1, 0, high);
						levelTwoLocalMapping[u][1]++;
						levelTwoTransitDistTable[u][1]++;
						// insert to the distance table
					}
				}
			}
			
		}
		//matrixout.close();
		elapsedTime = timestamp() - start;
	}	

	//destructor
	~ManyToMany()
	{
	};


private:
    Graph * _g;
    DijkstraFW _dFW;
    DijkstraBW _dBW;
    LevelID _earlyStopLevel;
	
	stack< pair<NodeID,OtherEdgeID> > _expandPathStack;

	//below are used for output the many to many paths between sources and targets
	//revised by dingxiong
	vector <vector <NodeID>> _viaNode;
	vector <DynQueryPQueue> _rootedTreeBW;
	vector <DynQueryPQueue> _rootedTreeFW;
	
	//sources and target maps to the internal ID
	map <NodeID, NodeID> _sourcesMapToID;
	map <NodeID, NodeID> _targetsMapToID;

	vector <map <NodeID, NodeID> > _sourcesPQIndex;
	vector <map <NodeID, NodeID> > _targetsPQIndex;
	

};

#endif // MANYTOMANY_H
