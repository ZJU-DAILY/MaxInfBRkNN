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

#ifndef _PROCESSING_CONSTRUCTCH_H
#define _PROCESSING_CONSTRUCTCH_H

#include "../EliminationWeight.h"
#include "DijkstraCH.h"

// Intervals for verbose output.
static const unsigned int PROGRESS_NODES_INIT = 100000;
static const unsigned int PROGRESS_NODES_ELIMINATE = 100000;
static const unsigned int PROGRESS_NODES_CONSTRUCT = 100000;

static const unsigned int SAVE_STATS_NONE = 0;
static const unsigned int SAVE_STATS_CONTRACT = 1;
static const unsigned int SAVE_STATS_ALL = 2;

typedef double BetweennessValue;
typedef unsigned int ReachValue;

namespace processing {
    
    /**
     * Constructs a Contraction Hieararchy (CH). This is the main class for preprocessing.
     * Either a CH is created from scratch including a new node order using a priority 
     * queue, see createHierarchy(), or a hierarchy construction with a given node order
     * is performed using constructHierarchy().
     * @param Graph graph class
     * @param saveShortcutsWitnesses should shortcuts and witnesses be written to file? 
     *                               note: local edge reduction and one hop backward search
     *                                     currently not supported if TRUE
     * @param saveShortcutsText save shortcuts to text file, requested by D.Delling
     * @param saveStats save stats to file (slows down node ordering)
     * @param oneHopBackwardSearch on k-hop-search, do (k-1)-hop forward-search with dijktra and 1-hop backward (edge array scan)
     * @param onTheFlyWitnessCheck if true, do not create "witness-shortcuts" before 
     *        the contraction and read the witnesses directly from file during contraction.
     */
    template 
        <
        typename Graph,
        bool saveShortcutsWitnesses = false,
        bool saveShortcutsText = false,
        unsigned int saveStats = SAVE_STATS_NONE,
        bool oneHopBackwardSearch = true,
        bool onTheFlyWitnessCheck = false
        >
	class ConstructCH {
    public:

        enum ConstructPhase
        {
            PHASE_NODEORDER_INIT,             // initialize priority queue with all nodes
            PHASE_NODEORDER_ELIMINATE,        // eliminate nodes by priority queue
            PHASE_CONSTRUCT,                  // eliminate nodes by previously known node order
        };

        /**
         * Represents parameters used for node ordering. It specifies the coefficients for 
         * the linear combination of the priority terms, the settled nodes and hop limits to the 
         * local searches and other parameters.
         * Those parameters are usually set in ../command/NodeOrder.h
         */
        struct WeightCalculation
        {
            WeightCalculation() : 
                edgeDiffMult(0), newEdgesMult(0), delNeighbMult(0), searchSpaceMult(0), 
                betweennessAdd(0), reachAdd(0), voronoiMult(0), searchPathHopBorderMult(0), 
                searchPathHopBorderOriginalEdges(false), shortcutOriginalEdgeSumMult(0),
                maxSettledApprox(0), maxSettledElim(0),  lazyUpdateRecalcLimit(0), updateHops(false), 
                localReduceEdges(!saveShortcutsWitnesses)
                
            {
            }
            // coefficient for edge difference
            int edgeDiffMult;
            // coefficient for number of new edges
            int newEdgesMult;
            // coefficient for deleted neighbors
            int delNeighbMult;
            // coefficient for search space size of local searches for contraction
            double searchSpaceMult;
            // coefficient for relative betweenness
            int betweennessAdd;
            // coefficient for relative reach
            int reachAdd;
            // coefficient for Voronoi region size
            int voronoiMult;
            // coefficient for upper bound on search paths hops
            int searchPathHopBorderMult;
            // in combination with searchPathHopBorderMult does calculate an upper bound on the costs
            // for a time-dependend search, requested by Veit Batz
            bool searchPathHopBorderOriginalEdges;
            // coefficient for sum of original edges new shortcuts represent
            int shortcutOriginalEdgeSumMult;            
            // limit of settled nodes in local searches during weight calculation
            unsigned int maxSettledApprox;
            // limit of settled nodes in local searches during contraction
            unsigned int maxSettledElim;
            // hop-/degree-limits
            vector<double> maxHops;
            // lazy updates, parameter also specifies the check intervals
            // for complete recalculation of priority queue
            unsigned int lazyUpdateRecalcLimit;
            // perform update on all affected nodes after contraction
            // this only works with hop-limits and is really slow
            // should only be used for testing
            bool updateHops;
            // wheter to use local edge reduction (recommended)
            bool localReduceEdges;
        };

        /**
         * Represents parameters used for hierarchy construction. It specifies 
         * the settled nodes and hop limits to the 
         * local searches and other parameters.
         * Those parameters are usually set in ../command/Construct.h
         */
        struct ContractParameters
        {
            ContractParameters() : 
                maxSettledElim(0), localReduceEdges(true), coreSize(0)
            { }
            unsigned int maxSettledElim;
            vector<double> maxHops;
            bool localReduceEdges;
            NodeID coreSize;
        };

        ///////////////////////////////////////////////////////////////////////
        /**
         * Constructor.
         * @param graph Graph object, will be modified so that it becomes a CH.
         */
        ConstructCH(Graph* graph) 
            :
            _graph(graph),
            _localDijkstra(graph),
            _dVoronoi(graph),
            _updateBfs(graph->noOfNodes(), false),
            _testDijkstra(graph)
        {
        }
        
        ///////////////////////////////////////////////////////////////////////
        ~ConstructCH()
        {
            // cleanup objects that are stored in pointers
        }


        ///////////////////////////////////////////////////////////////////////
        /** 
         * Main method to create a hierarchy including node order.
         * This method is the public method that should be called.
         * @param weightCalc node ordering and contraction parameters, see class
         *        WeightCalc above.
         * @param betweenness betweenness values of each node, can be NULL
         * @param reach reach values of each node, can be NULL
         * @param statsFile filename prefix for statistic files
         * @param testShortestPath perform time-consuming tests, only for debugging purpose
         */
        void createHierarchy(WeightCalculation weightCalc, 
            vector<BetweennessValue>* betweenness, 
            vector<ReachValue>* reach,
            string statsFile = "", bool testShortestPaths = false)
        {
            // initalize global variables
            _weightCalc = weightCalc;
            _betweenness = betweenness;
            _reach = reach;
            _statsFile = statsFile;
            _testShortestPaths = testShortestPaths;
            _localReduceEdges = weightCalc.localReduceEdges;

            assert( !saveShortcutsWitnesses || (!_localReduceEdges && !oneHopBackwardSearch) );

            // perform time-consuming tests during the hierarchy creation
            // used for debugging purposes
            if (_testShortestPaths) initTestShortestPaths();
                
            // initialize priority queue with all nodes
            initPQueue();

			eliminateByPQueue();
        }
        
        
        /**
         * Main method to create a hierarchy if the node order is already given
         * by readLevels().
         * @param contractParams contraction parameters, see class ContractParameters above.
         */
        void constructHierarchy(ContractParameters& contractParams)
        {
            _contractParams = contractParams;
            _localReduceEdges = contractParams.localReduceEdges;
            constructByLevel();
        }
        
        /**
         * Main method to create a hierarchy if the node order is already given
         * by readLevels() and witness paths and shortcuts are stored.
         * This is only a preliminary version supporting witness paths.
         * @param contractParams contraction parameters, see class ContractParameters above.
         * @param inShortcuts input stream to read shortcut edges
         * @param inWitnesses input stream to read witness paths
         */
        void constructHierarchyWithWitnesses(ContractParameters& contractParams, istream& inShortcuts, istream& inWitnesses)
        {
            _contractParams = contractParams;
            constructByLevelRetrieve(inShortcuts, inWitnesses);
        }
        
        
        /**
         * Write exact levels (node order) to file from 0 to n-1.
         * The levels are stored in the graph datastructure and
         * are calculated during the hierarchy creation by priority queue,
         * see createHierarchy().
         * Format: (binary)
         *         no of nodes (unsigned int)
         *         per node: level (unsigned int)
         */
        void writeLevels(ostream& out)
        {
             // first write no of nodes as unsigned int
             unsigned int buffer = _graph->noOfNodes();
             out.write((char*)&buffer, sizeof(unsigned int)/sizeof(char));
             // write levels to file as unsigned int
             for (NodeID node = 0; node < _graph->noOfNodes(); node++) {
                 buffer = pqData(node).level;
                 assert( buffer != UINT_MAX );
                 out.write((char*)&buffer, sizeof(unsigned int)/sizeof(char));
             }
        }
        
        /**
         * Read exact levels (node order) from file from 0 to n-1.
         * The levels will be stored in a separate array (_levelNodeList) and will
         * be stored in the graph datastructure not until the contraction
         * of the node, constructHierarchy().
         * Format: (binary)
         *         no of nodes (unsigned int)
         *         per node: level (unsigned int)
         */
        void readLevels(ifstream& in)
        {
            // read levels from stream
            unsigned int buffer;
            in.read((char*)&buffer, sizeof(unsigned int)/sizeof(char));
            assert( buffer == _graph->noOfNodes() );

            _levelNodeList.clear();
            _levelNodeList.reserve(_graph->noOfNodes());
            for (NodeID u = 0; u < _graph->noOfNodes() ; u++) {
                in.read((char*)&buffer, sizeof(unsigned int)/sizeof(char));
                _levelNodeList.push_back( make_pair( buffer, u ) );
            }
            sort( _levelNodeList.begin(), _levelNodeList.end() );
            
            _noOfNodes = _graph->noOfNodes();
            _noOfEdges = _graph->noOfExistingEdges();
            
        }
        
        bool savesShortcutsWitnesses()
        {
            return saveShortcutsWitnesses;   
        }
        
        /**
         *  Output stream to write shortcuts, only used if template parameter
         *  saveShortcutsWitnesses is enabled.
         */
        void storeShortcuts(ostream& out)
        {
            _storeShortcuts.store(out);
        }

        /**
         *  Output stream to write witness paths, only used if template parameter
         *  saveShortcutsWitnesses is enabled.
         */
        void storeWitnesses(ostream& out)
        {
            _storeWitnesses.store(out);
        }
        
        /**
         *  Returns the maximum average degree occured during the last hierarchy creation.
         */
        double lastMaxAvgDegree() 
        {
            return _maxAvgDegree;
        }

      
		/**
         * Cleanum after hierarchy construction.
         */        
        void clear()
    	{
    	    _levelNodeList.clear();
    	}
    	
        /**
         * Partition a file with exact levels (total node order)
         * into given partitions, first partition contains highest levels.
         * This is used to convert a CH to a hierarchy that can be used by
         * Dominik Schultes HNR and TNR code.
         * @param in input hcn-file
         * @param out output hcn-file
         * @param levelsStr comma-separated list of partition/level sizes
         * @param a if this parameter is != 0, then the a most unimportant nodes will get a 
         *          common level and above them has each node its own level.
         */
        // Note: Command::CreateVector is commented out so code below does not work
        // But this function is not used anywhere so we are commenting this out as well
//         static void partitionExactInLevels(istream& in, ostream& out, const string& levelsStr, unsigned int a)
//         {
//             unsigned int buffer;
//             in.read((char*)&buffer, sizeof(unsigned int)/sizeof(char));
//             
//             vector<NodeID> nodeByLevel(buffer, SPECIAL_NODEID); // reversed, index 0 is highest level
//             vector<unsigned int> levelByNode(buffer, UINT_MAX);
//             
//             if (a == 0)
//             {
//                 for ( unsigned int i = 0; i < levelByNode.size(); i++ )
//                 {
//                     in.read((char*)&buffer, sizeof(unsigned int)/sizeof(char));
//                     assert( buffer != UINT_MAX );
//                     assert( buffer < nodeByLevel.size() );
//                     assert( nodeByLevel[nodeByLevel.size() - 1 - buffer] == SPECIAL_NODEID );
//                     nodeByLevel[nodeByLevel.size() - 1 - buffer] = i;
//                 }
// 
//                 // parse levels str
//                 vector<unsigned int> levels;
//                 Command::createVector(levelsStr,levels,buffer);
//                 VERBOSE( cout << "levels "; copy(levels.begin(),levels.end(), 
//                 ostream_iterator<NodeID>(cout, ",")); cout << endl );
//                 int level = levels.size();
//                 NodeID levelCount = 0;
//                 NodeID levelLimit = level > 0 ? levels[0] : 0;
//     
//                 for ( vector<NodeID>::iterator iter = nodeByLevel.begin(); iter != nodeByLevel.end(); iter++ )
//                 {
//     
//                     levelByNode[*iter] = level;
//     
//                     levelCount++;
//                     if (levelCount == levelLimit) {
//                         level--;
//                         levelCount = 0;
//                         levelLimit = level > 0 ? levels[levels.size()-level] : 0;
//                     }
//                 }
//             }
//             else
//             {
// 
//                 for ( unsigned int i = 0; i < levelByNode.size(); i++ )
//                 {
//                     in.read((char*)&buffer, sizeof(unsigned int)/sizeof(char));
//                     assert( buffer != UINT_MAX );
//                     assert( buffer < nodeByLevel.size() );
//                     assert( nodeByLevel[buffer] == SPECIAL_NODEID );
//                     nodeByLevel[buffer] = i;
//                 }
// 
//                 for ( NodeID i = 0; i < a; i++ )
//                 {
//                     levelByNode[nodeByLevel[i]] = 0;   
//                 }
//                 for ( NodeID i = a; i < levelByNode.size(); i++ )
//                 {
//                     levelByNode[nodeByLevel[i]] = i-(a-1);   
//                 }
//             }
// 
//              // first write no of nodes as unsigned int
//              buffer = levelByNode.size();
//              out.write((char*)&buffer, sizeof(unsigned int)/sizeof(char));
//             // write levels to file as unsigned int
//              for (NodeID node = 0; node < levelByNode.size(); node++) {
//                 assert( levelByNode[node] < UINT_MAX );
//                 buffer = levelByNode[node];
//                 out.write((char*)&buffer, sizeof(unsigned int)/sizeof(char));
//              }
//         }
    	
        
        
    private:


        /**
         *  Data structure stored for each node in the priority queue used to
         *  determine the node order. Note that this information is even accessible
         *  after a node has been removed from the priority queue.
         */
        struct PQueueNodeElimination
        {
            PQueueNodeElimination(): 
                   level(UINT_MAX), searchSpace(0), deletedNeighbors(0), 
                   searchPathHopBorder(0),
                   voronoiOwner(SPECIAL_NODEID), voronoiNextBorderNode(SPECIAL_NODEID), voronoiNumber(0)
                {}
                
            // Level of the node. This is initially UINT_MAX and will be set 
            // to the calculated level directly after the contraction 
            // of the node.
            unsigned int level;
            
            // Size of the search spaces of the local searches. used
            // as priority term to speedup the contraction.
            unsigned int searchSpace;
            
            // deleted/contracted neighbors counter
            unsigned int deletedNeighbors;
            
            // Upper bound on search paths found by the query algorithm.
            // Datatype double is used since this bound is modified for
            // time-dependend routing (parameter "-T")
            double searchPathHopBorder;
            
            // Those three variables are used for the priority term
            // Voronoi regions. They allow to manage the Voronoi regions
            // of the remaining nodes.
            // The voronoi region that this node belongs to.
            NodeID voronoiOwner;
            // The nodes in a voronoi region are stored in a linked list.
            // The list is terminated by SPECIAL_NODEID.
            NodeID voronoiNextBorderNode;
            // Contains the number of nodes in the voronoi region if the nodes
            // is not contracted. Otherwise, if the node is already contracted,
            // it contains the distance to the voronoiOwner.
            unsigned int voronoiNumber;

            /** Wheter the node is eliminated, see unsigned int level. */
            bool isEliminated()
            {
                return level != UINT_MAX;
            }
        };

        /**
         * Datastructure to temporarely store new shortcut edges. An additional
         * structure is required since the Edge class does not store the source
         * node. See ../datastr/graph/edge.h for the Edge class.
         */
        struct NewEdge
        {
            ///////////////////////////////////////////////////////////////////
            NewEdge (const NodeID source, const NodeID target, const EdgeWeight weight, const bool isBidirected,
                     const NodeID shortcutMiddle, const OtherEdgeID shortcutParentEdge, const OtherEdgeID shortcutChildEdge,
                     const OtherEdgeID shortcutOriginalEdgeCount)
                : source(source), edge(target, weight, EDGE_TYPE_SHORTCUT, true, isBidirected, 
                         shortcutMiddle, shortcutParentEdge, shortcutChildEdge, shortcutOriginalEdgeCount)
            { }
            NodeID source;
            Edge edge;
        };
        
        /** 
         * Store shortcuts in a file (binary format), those are used to speedup the hierarchy construction
         * in a second step. (only preliminary)
         */
        class StoreShortcuts
        {
            public:
            StoreShortcuts(): _out(NULL), _lastNode(SPECIAL_NODEID), _lastIn(SPECIAL_NODEID) { }
            /** Output stream to store shortcuts to. */
            void store(ostream& out) { _out = &out; }
            /** 
             * Write a shortcut to file, it omits the (middle) node and the input node if
             * it is the same as the previous shortcut. This is achieved using flags.
             */
            void addShortcut(const NodeID node, const NodeID in, const NodeID out, bool bidir)
            {
                if (_out != NULL )
                {
                    if (node != _lastNode)
                    {
                        Data data(TYPE_NODE, node);
                        _out->write((char*)&data,sizeof(Data)/sizeof(char));
                        //cout << "write TYPE_NODE, " << node << endl;
                        _lastNode = node;
                        _lastIn = SPECIAL_NODEID;
                    }
                    if (in != _lastIn)
                    {
                        Data data(TYPE_IN, in);
                        _out->write((char*)&data,sizeof(Data)/sizeof(char));
                        //cout << "write TYPE_IN, " << in << endl;
                        _lastIn = in;
                    }
                    if (bidir)
                    {
                        Data data(TYPE_OUT_BIDIRECTIONAL, out);
                        _out->write((char*)&data,sizeof(Data)/sizeof(char));
                        //cout << "write TYPE_OUT_BIDIRECTIONAL, " << out << endl;
                    }
                    else
                    {
                        Data data(TYPE_OUT_UNIDIRECTIONAL, out);
                        _out->write((char*)&data,sizeof(Data)/sizeof(char));
                        //cout << "write TYPE_OUT_UNIDIRECTIONAL, " << out << endl;
                        //_out->flush();
                        //exit(1);
                    }
                }
            }
                
            private:
            // Node id is the contracted node.
            static const unsigned int TYPE_NODE = 0;
            // Node id is a source node, as contracted node the previously 
            // written contracted node is used.
            static const unsigned int TYPE_IN = 1;
            // Node id is the target node, as contracted node and source
            // node the previously written nodes are used.
            // Distinguish between unidirectional and bidirectional shortcuts.
            static const unsigned int TYPE_OUT_UNIDIRECTIONAL = 2;
            static const unsigned int TYPE_OUT_BIDIRECTIONAL = 3;

            /** 
             * Data structure, it can store a node id plus a 2-bit flag that 
             * allows a more space efficient storage of shortcuts, see addShortcut()
             */ 
            struct Data
            {
                Data(unsigned int type, unsigned int id): type(type), id(id) { }
                unsigned int type:2;
                unsigned int id:30;
            };
            ostream* _out;   
            NodeID _lastNode;
            NodeID _lastIn;
        };
        
        /** 
         * Retrieve shortcuts from file (binary format), those are used to speedup the hierarchy construction.
         * (only preliminary)
         */
        class RetrieveShortcuts
        {
            public:
            /** @param in Input stream to retrieve shortcuts from. */
            RetrieveShortcuts(istream& in)
            : _in(in), _lastNode(SPECIAL_NODEID), _lastIn(SPECIAL_NODEID) 
            {
                _preEof = _in.eof();
                read();
            }

            /**
             * The shortcuts are stored in the order that they are needed during construction,
             * see StoreShortcuts.
             * "Node" is the contracted node, "In" and "Out" are the two incident nodes of the
             * shortcut edge.
             * To extract the shortcuts, three nested loops are required.
             * while (nextNode(node)) { while nextIn(in)) { while nextOut(out)) { ... }}}
             */
            bool nextNode(NodeID& node)
            {
                if ( _eof ) return false;
                assert( _data.type == TYPE_NODE );
                node = _data.id;
                read();
                return true;
            }
            bool nextIn(NodeID& in)
            {
                if ( _eof || _data.type != TYPE_IN ) return false;
                assert( _data.type == TYPE_IN );
                in = _data.id;
                read();
                return true;
            }
            bool nextOut(NodeID& out, bool& bidir)
            {
                if ( _eof || (_data.type != TYPE_OUT_UNIDIRECTIONAL && _data.type != TYPE_OUT_BIDIRECTIONAL) ) return false;
                assert( _data.type == TYPE_OUT_UNIDIRECTIONAL || _data.type == TYPE_OUT_BIDIRECTIONAL );
                out = _data.id;
                bidir = (_data.type == TYPE_OUT_BIDIRECTIONAL);
                read();
                return true;
            }
                
            private:
            // Node id is the contracted node.
            static const unsigned int TYPE_NODE = 0;
            // Node id is a source node, as contracted node the previously 
            // read contracted node is used.
            static const unsigned int TYPE_IN = 1;
            // Node id is the target node, as contracted node and source
            // node the previously read nodes are used.
            // Distinguish between unidirectional and bidirectional shortcuts.
            static const unsigned int TYPE_OUT_UNIDIRECTIONAL = 2;
            static const unsigned int TYPE_OUT_BIDIRECTIONAL = 3;
            /** 
             * Data structure, it can store a node id plus a 2-bit flag that 
             * allows a more space efficient storage of shortcuts.
             */ 
            struct Data
            {
                unsigned int type:2;
                unsigned int id:30;
            };
            istream& _in;   
            NodeID _lastNode;
            NodeID _lastIn;
            Data _data;
            bool _eof;
            bool _preEof;
            void read()
            {
                if ( !_preEof )
                {
                    assert( !_in.eof() );
                    _in.read((char*)&_data, sizeof(Data)/sizeof(char));
                }
                _eof = _preEof;
                _preEof = _in.eof();
            }
        };
        
        /** 
         * Store witness paths in a file (binary format), those are used to speedup the hierarchy construction
         * in a second step, see Construct.h. (only preliminary)
         */
        class StoreWitnesses
        {
            public:
            StoreWitnesses(): _out(NULL), _lastIn(SPECIAL_NODEID), _lastNode(SPECIAL_NODEID) { }
            /** Output stream to store witness paths to. */
            void store(ostream& out) { _out = &out; }
            /** 
             * Write the initial information of a witness path to file, it omits the (middle) node and the input node if
             * it is the same as the previous path. This is achieved using flags.
             */
            void beginWitness(const NodeID node, const NodeID in, const NodeID out)
            {
                if (_out != NULL )
                {
                    if (node != _lastNode)
                    {
                        Data data1(TYPE_NODE, node);
                        _out->write((char*)&data1,sizeof(Data)/sizeof(char));
                        _lastNode = node;
                        Data data2(TYPE_IN, in);
                        _out->write((char*)&data2,sizeof(Data)/sizeof(char));
                        _lastIn = in;
                    }
                    else if (in != _lastIn)
                    {
                        Data data(TYPE_IN, in);
                        _out->write((char*)&data,sizeof(Data)/sizeof(char));
                        _lastIn = in;
                    }
                    Data data(TYPE_OUT, out);
                    _out->write((char*)&data,sizeof(Data)/sizeof(char));
                }
            }
            
            /** 
             * Store the next step (node) of the path.
             */
            void nextStep(const NodeID step)
            {
                if (_out != NULL )
                {
                    Data data(TYPE_STEP, step);
                    _out->write((char*)&data,sizeof(Data)/sizeof(char));
                }
            }
                
            private:
            // Node id is the contracted node.
            static const unsigned int TYPE_NODE = 0;
            // Node id is a source node, as contracted node the previously 
            // written contracted node is used.
            static const unsigned int TYPE_IN = 1;
            // Node id is the target node, as contracted node and source
            // node the previously written nodes are used.
            static const unsigned int TYPE_OUT = 2;
            // Node id is the next node on the witness path having the contracted node, source
            // and target node the previously written nodes used.
            static const unsigned int TYPE_STEP = 3;
            /** 
             * Data structure, it can store a node id plus a 2-bit flag that 
             * allows a more space efficient storage of witness paths, see beginWitness()
             */ 
            struct Data
            {
                Data(unsigned int type, unsigned int id): type(type), id(id) { }
                unsigned int type:2;
                unsigned int id:30;
            };
            ostream* _out;   
            NodeID _lastIn;
            NodeID _lastNode;
        };
        
        /** 
         * Retrieve witness paths from file (binary format), those are used to speedup the hierarchy construction.
         * (only preliminary)
         */
        class RetrieveWitnesses
        {
            public:
            /** @param in input stream to retrieve witnesses from. */
            RetrieveWitnesses(istream& in)
            : _in(in), _lastNode(SPECIAL_NODEID), _lastIn(SPECIAL_NODEID), _lastOut(SPECIAL_NODEID)
            {
                _preEof = _in.eof();
                read();
            }
            
            /**
             * The witness paths are stored in the order that they are needed during construction,
             * see StoreWitnesses.
             * nextWitness() returns the contracted node and the source and target
             * node of the witness paths. The nodes within the path are returned by nextStep().
             * Since we do not know in advance if a witness exists, we can
             * "peek" for the next witness without moving the "file pointer" forward. (simply speaking)
             * If we only peeked, we need to perform a read() operation before we can read
             * the nodes within the witness paths with nextStep().
             */
            bool nextWitness(NodeID& node, NodeID& in, NodeID& out, bool peek = false)
            {
                if ( _eof ) return false;
                if (_data.type == TYPE_NODE)
                {
                    _lastNode = _data.id;
                    read();
                    assert( !_eof );
                }
                if (_data.type == TYPE_IN)
                {
                    _lastIn = _data.id;
                    read();
                    assert( !_eof );
                }
                assert( _data.type == TYPE_OUT );
                node = _lastNode;
                in = _lastIn;
                out = _data.id;
                if ( !peek ) read();
                return true;
            }
            
            /**
             * Returns the nodes within the witness paths.
             * @return true as long as there is a step, used in a while-loop.
             */
            bool nextStep(NodeID& step)
            {
                if ( _eof || _data.type != TYPE_STEP ) return false;
                assert( _data.type == TYPE_STEP );
                step = _data.id;
                read();
                return true;
            }

            /**
             * Perform a read operation. Required if we "peeked"
             * the next witness path with nextWitness().
             */
            void read()
            {
                if ( !_preEof )
                {
                    assert( !_in.eof() );
                    _in.read((char*)&_data, sizeof(Data)/sizeof(char));   
                }
                _eof = _preEof;
                _preEof = _in.eof();
            }
            private:
            // Node id is the contracted node.
            static const unsigned int TYPE_NODE = 0;
            // Node id is a source node, as contracted node the previously 
            // read contracted node is used.
            static const unsigned int TYPE_IN = 1;
            // Node id is the target node, as contracted node and source
            // node the previously read nodes are used.
            static const unsigned int TYPE_OUT = 2;
            // Node id is the next node on the witness path having the contracted node, source
            // and target node the previously read nodes used.
            static const unsigned int TYPE_STEP = 3;
            /** 
             * Data structure, it can store a node id plus a 2-bit flag that 
             * allows a more space efficient storage of witness paths.
             */ 
            struct Data
            {
                unsigned int type:2;
                unsigned int id:30;
            };
            istream& _in;   
            NodeID _lastNode;
            NodeID _lastIn;
            NodeID _lastOut;
            Data _data;
            bool _eof;
            bool _preEof;
        };
        
        
    
        // graph datastructure, usually a UpdateableGraph is used, see ../datastr/graph/UpdateableGraph.h
        Graph* _graph;
        
        // Contains the coefficients for the linear combination of the priority and the local search and hop limits.
        WeightCalculation _weightCalc;

        // Contains the local search and hop limits used during hierarchy construction if node order is given.
        ContractParameters _contractParams;
        
        // Priority queue used for node order calculation.
        BinaryHeap< EliminationWeight::Type, EliminationWeight, PQueueNodeElimination, NodeID > _pqElimination;
            
        // Dijkstra algorithm used for local searches during contraction, see DijkstraCH.h
        LocalDijkstraContract _localDijkstra;

        // number of edges in the remaining (not contracted) graph
        OtherEdgeID _noOfEdges;
        // number of nodes in the remaining (not contracted) graph
        NodeID _noOfNodes;
        
        // For debugging purpose only, perform time-consuming test during hierarchy creation.
        bool _testShortestPaths;
        
        // Store betweeness value, node id as index.
        vector<BetweennessValue>* _betweenness;
        // Store reach value, node id as index.
        vector<ReachValue>* _reach;
        
        // Array for new edges, used in processNode().
        vector<NewEdge> _newEdges;
        
        // Filename prefix for statistic files.
        string _statsFile;
        // Output filestream for main statistic file that contains a line for each
        // contracted node. Template parameter saveStats needs to be enabled.
        ofstream _stats;
        
        // Output filestream, each line starts with the contracted node
        // followed by the incident nodes of all shortcuts. Template parameter 
        // saveShortcutsText needs to be enabled.
        ofstream _shortcutsText;
        
        // Output filestream, used to calculate the size of the stored witnesses
        // if stored as trees of shortcuts. Only for testing.
        ofstream _tree;
        
        // Objects to store shortcuts and witness paths to file. Template parameter
        // saveShortcutsWitnesses needs to be enabled.
        StoreShortcuts _storeShortcuts;
        StoreWitnesses _storeWitnesses;
        
        // Counts the number of witnesses (value) with given length (index). those
        // values will be written to file if template parameter saveStats is enabled.
        vector<unsigned int> _witnessCounter;
        // Same as _witnessCounter, but does regard (as index) the number of settled
        // nodes during local search until a witness is found.
        vector<unsigned int> _witnessSearchSpaceCounter;
        
        // Current level during contraction.
        unsigned int _currentLevel;
        
        // Switch deciding wheter lazy updates are enabled.
        bool _lazyUpdate;
        
        // Current hop-limit for local searches. 0 = no hop-limit.
        unsigned int _maxHops;
        // Current degree-limit that will trigger a hop-limit change.
        double _maxHopsDegreeLimit;
        // The hop-/degree-limits are stored in a vector in the class WeightCalculation. This
        // variable points to the index of the current hop-limit.
        unsigned int _maxHopsIndex;
        // Use local edge reduction.
        bool _localReduceEdges;
        
        // Dijkstra object to update Voronoi regions. The Voronoi region of the contracted
        // node is distributed among the neighboring regions based on shortest paths.
        DijkstraUpdateVoronoi _dVoronoi;
        // BFS queue, if updateHops is enabled, update all nodes within hop-limit.
        vector<bool> _updateBfs;
        // maximum average degree of all remaining graphs, used for statistics
        double _maxAvgDegree;

        // List of pairs <level,node id> sorted by level ascending.
        // This list is filled by readLevels() and used during hierarchy construction
        // since the level of a node is stored in the graph datastructure not
        // until the node is contracted.
        vector< pair<LevelID,NodeID> > _levelNodeList;

        // Used for debugging purpose. Perform time-consuming tests during
        // hierarchy creation. Theses parameters specify the number of tests etc.
        static const unsigned int _noOfTestCases = 10;
        stPairs _testRuns;
        vector<Path> _testPaths;
        processing::DijkstraCH<datastr::graph::UpdateableGraph, DynQueryPQueue, 2, true> _testDijkstra;
    
        /**
         *  Maps node id to index into priority queue _pqElimination.
         *  Since the priority queue is initalized by ascending node id,
         *  the mapping is simple.
         */
        inline NodeID pqNodeToIndex(NodeID node)
        {
            return node + 1;
        }
        
        /**
         *  Maps index of the priority queue _pqElimination to the node id.
         *  Since the priority queue is initalized by ascending node id,
         *  the mapping is simple.
         */
        inline NodeID pqIndexToNode(NodeID index)
        {
            assert( index >= 1 );
            return index - 1;
        }
        
        /**
         * Returns the data stored with each node (element) index
         * the priority queue. Note that this data is still available
         * even a node has been deleted from the pqueue.
         */
        PQueueNodeElimination& pqData(NodeID node)
        {
            assert( pqNodeToIndex(node) < _pqElimination.elements().size() );
            return _pqElimination.elements()[pqNodeToIndex(node)].data();   
        }

        /**
         * Insert all nodes into priority queue.
         */
        void initPQueue() {
            VERBOSE( cout << "Initialize elimination weights..." << endl );
            VERBOSE( double timeStart = timestamp(); );
            VERBOSE( double timeLast = timestamp(); );
            
            // initalize statistical variables
            _noOfEdges = 0;
            _noOfNodes = 0;
            
            // Prepare hop-limits, since staged hop-limits are used that change
            // after reaching certain average degrees, some overhead is required.
            initHopLimit(_weightCalc.maxHops);

			//cout <<  "maxHops " << _maxHops << " maxdegreelimit " << _maxHopsDegreeLimit << endl;

            // start with level 0, the most unimportant level
            _currentLevel = 0;
            
            // enable lazy update if the check interval > 0
            _lazyUpdate = _weightCalc.lazyUpdateRecalcLimit > 0;
            
            // need to init witness array for fast mtm 2 hops search            
            // if a 2-hop search is in the list of hop-limits
            initPossibleWitnessesMTM(_weightCalc.maxHops);

            // initialize memory for priority queue, important for large graphs otherwise the memory
            // becomes fragmented and a bad_malloc occurs.
            _pqElimination.reserve( _graph->noOfNodes() );
            
            VERBOSE( Percent percent( _graph->noOfNodes() ) );
            
            // Nodes are inserted into the priority queue by ascending node id.
            // This allows a simple mapping from node id to index of the priority queue
            // element.
			//fstream fp("weight.txt", ios::out);
            for ( NodeID node = 0; node < _graph->noOfNodes(); node++ )
            {
                VERBOSE( 
                    if (node > 0 && node % PROGRESS_NODES_INIT == 0) 
                    {
                        double now = timestamp();
                        cout << node << " nodes, " << (now-timeStart) << " mseconds, " ;
                        cout << (now-timeLast) << " mseconds last, ";
                        cout << ((((now-timeStart) / (node)) * _graph->noOfNodes())-(now-timeStart)) << " mseconds remaining";
                        cout << endl;
                        timeLast = now;
                    }
                )
                    
                VERBOSE_CONTRACT( cout << "node " << node << flush;)
                

                // Calculate elimination weight/priority. Some priority terms are returned
                // that need to be stored for each node for later priority updates.
                VERBOSE_CONTRACT( double time1 = timestamp(); )
                unsigned int searchSpace;
                int edgeDiff;
                unsigned int newEdges;
                EliminationWeight::Type weight = calculateEliminationWeight<PHASE_NODEORDER_INIT>(node, &searchSpace, &edgeDiff, &newEdges);
				
			
				//fp << node << edgeDiff << " " << weight << endl;
			
                    
                // insert node into priority queue
                VERBOSE_CONTRACT( time1 = timestamp() - time1; )
                VERBOSE_CONTRACT( cout << " calc " << fixed << setprecision(3) << time1 * 1000 << flush; )
                _pqElimination.insert(weight);
                //assert( index == pqNodeToIndex(node) );
                
                // store priority terms (node information) in the data structure assigned
                // to each node, required for later priority updates.
                PQueueNodeElimination& data = pqData(node);
                data.searchSpace = searchSpace;
                data.voronoiOwner = node;
                assert( data.deletedNeighbors == 0 );

                // count number of remaining nodes and edges
                _noOfNodes++;
                _noOfEdges += _graph->lastEdge(node) - _graph->firstEdge(node);
                VERBOSE_CONTRACT( cout << " edges " << _noOfEdges << flush; )

                VERBOSE_CONTRACT( cout << endl; )
                VERBOSE( percent.printStatus(node) );
            }
            
            VERBOSE( cout << "#nodes: " << _noOfNodes << " / #edges: " << _noOfEdges << endl; )
            assert( _pqElimination.checkHeapProperty() );

        }
        

        
        /** 
         * Calculates the elimination weight (priority) of a remaining node in the current
         * state of contraction. The priority is a linear combination of several 
         * priority terms.
         * There are several return values (call-by-reference) that are
         * returned, if != NULL.
         * @param phase phase of hierarchy creation, during PHASE_NODEORDER_INIT no accessible
         *        to the priority queue is allowed since not all nodes have been inserted.
         */
        template < ConstructPhase phase >
        EliminationWeight::Type calculateEliminationWeight(const NodeID node, 
            unsigned int* searchSpaceResult = NULL, int* edgeDiffResult = NULL, 
            unsigned int* newEdgesResult = NULL)
        {
            // variables representing several priority terms
            unsigned int searchSpace = 0;
            int edgeDiff = 0;
            unsigned int newEdges = 0;
            unsigned int inDegree = 0;
            unsigned int outDegree = 0;
            unsigned int deletedNeighbors = 0;
            unsigned int shortcutOriginalEdgeSum = 0;
            if ( phase != PHASE_NODEORDER_INIT )
            {
                deletedNeighbors = pqData(node).deletedNeighbors;
            }
            
            // Perform a simulated contraction of node to calculate several priority terms. This
            // yields e.g. the edge difference. The second template parameter, here true, specifies
            // the simulation.
            processNode<phase,true>(node,&searchSpace, &edgeDiff, &newEdges, &inDegree, &outDegree, 
                &shortcutOriginalEdgeSum);
            
            if (searchSpaceResult != NULL) *searchSpaceResult = searchSpace;
            if (edgeDiffResult != NULL) *edgeDiffResult = edgeDiff;
            if (newEdgesResult != NULL) *newEdgesResult = newEdges;
                
            // Result will contain the priority.
            // Priority terms: edge difference, number of new edges, deleted neighbors, search space of local searches
            EliminationWeight::Type result = ((EliminationWeight::Type)edgeDiff*_weightCalc.edgeDiffMult
                 + (EliminationWeight::Type)newEdges*_weightCalc.newEdgesMult
                 + (EliminationWeight::Type)deletedNeighbors*_weightCalc.delNeighbMult
                 + ((EliminationWeight::Type)searchSpace*_weightCalc.searchSpaceMult));
			
			//cout << node << " " << edgeDiff << " " << searchSpace << endl;
            // size of Voronoi region, extract square root
            if ( _weightCalc.voronoiMult != 0)
            {
                unsigned int voronoiNumber = 0;
                if ( phase != PHASE_NODEORDER_INIT )
                {
                    voronoiNumber = pqData(node).voronoiNumber;
                }
                result += sqrt(((EliminationWeight::Type)voronoiNumber)*_weightCalc.voronoiMult);
            }
                    
            // Upper bound on search paths found by query algorithm. There is a command-line 
            // argument (-T) that changes the behaviour to time-dependend queries. In this case,
            // the log10(+1) is used since the priority term eventually exceeds even a long for
            // the PTV Western Europe road network.
            if ( _weightCalc.searchPathHopBorderMult != 0) 
            {
                if ( phase != PHASE_NODEORDER_INIT )
                {
                    if ( !_weightCalc.searchPathHopBorderOriginalEdges )
                    {
                        result += ((EliminationWeight::Type)pqData(node).searchPathHopBorder*_weightCalc.searchPathHopBorderMult);
                    }
                    else
                    {
                        result += _weightCalc.searchPathHopBorderMult*log10(1+((EliminationWeight::Type)pqData(node).searchPathHopBorder));
                    }
                }
            }
                    
            // Relative betweenness, the fraction of remaining nodes with smaller betweenness
            // is the priority term in [0,1] and thus needs an approriate wheight.
            if (_weightCalc.betweennessAdd != 0)
            {
                result += calculateSmallerNeighbors(node, _betweenness)*_weightCalc.betweennessAdd;
            }
            // Relative reach, the fraction of remaining nodes with smaller reach
            // is the priority term in [0,1] and thus needs an approriate wheight.
            if (_weightCalc.reachAdd != 0)
            {
                result += calculateSmallerNeighbors(node, _reach)*_weightCalc.reachAdd;
            }
                
            // count the original edges, the shortcuts added during the simulated contraction represent
            if (_weightCalc.shortcutOriginalEdgeSumMult != 0)
            {
                result += _weightCalc.shortcutOriginalEdgeSumMult*sqrt((EliminationWeight::Type)shortcutOriginalEdgeSum);
            }
                            
            return result;
        }
        
        /**
         * Calculates the fraction of neighbors that have smaller value in vector values than unsigned.
         * Used to calculate the priority term for relative betweenness and reach.
         * @param u node-id
         * @param values vector of values (index node-id)
         */
        template < typename T >
        double calculateSmallerNeighbors(NodeID u, const vector<T>* values)
        {
            OtherEdgeID lastEdge = _graph->lastEdge(u);
            vector<NodeID> neighbors;
            for ( OtherEdgeID e = _graph->firstLevelEdge(u); e < lastEdge; e++ )
            {
                neighbors.push_back(_graph->edge(e).target());
            }
            sort( neighbors.begin(), neighbors.end() );
            NodeID previous = SPECIAL_NODEID;
            unsigned int noOfNeighbors = 0;
            unsigned int noOfSmallNeighbors = 0;
            T myValue = (*values)[u];
            for ( vector<NodeID>::const_iterator iter = neighbors.begin(); iter != neighbors.end(); iter++ )
            {
                if ( *iter == previous ) continue;
                previous = *iter;
                noOfNeighbors++;
                if ( (*values)[*iter] < myValue ) noOfSmallNeighbors++;
            }
            return noOfNeighbors == 0 ? 0 : (double)noOfSmallNeighbors/noOfNeighbors;
        }
        
       
        /** 
         * Calculates the elimination weight (priority) of node and updates the priority queue accordingly.
         */
        EliminationWeight::Type updateEliminationWeight(NodeID node, unsigned int* searchSpaceResult = NULL,
            int* edgeDiffResult = NULL, unsigned int* newEdgesResult = NULL)
        {

            unsigned int searchSpace = 0;
            int edgeDiff = 0 ;
            unsigned int newEdges = 0;
            
            // calculate elimination weight (priority), this only happens during the elimiation phase
            // (PHASE_NODEORDER_ELIMINATION), the other nodeorder phase is the initalization phase.
            EliminationWeight::Type newWeight = calculateEliminationWeight<PHASE_NODEORDER_ELIMINATE>(node, 
                &searchSpace, &edgeDiff, &newEdges);
                
            // update the priority queue key and data
            _pqElimination.updateKey( pqNodeToIndex(node), newWeight);
            pqData(node).searchSpace = searchSpace;

            // return some priority terms by call-by-reference
            if (searchSpaceResult != NULL) *searchSpaceResult = searchSpace;
            if (edgeDiffResult != NULL) *edgeDiffResult = edgeDiff;
            if (newEdgesResult != NULL) *newEdgesResult = newEdges;
            return newWeight;
        }
        

        /** 
         * Update the priorities of all remaining nodes in the priority queue.
         * This happens after a hop-limit change or too many lazy update within
         * the check interval.
         * 
         */
        void updatePQueue()
        {
            VERBOSE( cout << "Recalculate all remaining weights..." << endl; )
            unsigned int i = 0;
            VERBOSE( double timeStart = timestamp(); );
            VERBOSE( double timeLast = timestamp(); );
            VERBOSE( Percent percent(_noOfNodes) );
            for ( NodeID v = 0; v < _graph->noOfNodes(); v++ )
            {
                if (! pqData(v).isEliminated() )
                {
                    
                    VERBOSE( 
                        if (i > 0 && i % PROGRESS_NODES_INIT == 0) 
                        {
                            double now = timestamp();
                            cout << i << " nodes, " << (now-timeStart) << " seconds, " ;
                            cout << (now-timeLast) << " last, ";
                            cout << ((((now-timeStart) / (i)) * _noOfNodes)-(now-timeStart)) << " remaining";
                            cout << endl;
                            timeLast = now;
                        }
                    )

                    updateEliminationWeight(v);
                    
                    VERBOSE( percent.printStatus(i) );
                    i++;
                }
            }
        }
        
        /** 
         * Main phase: node contraction by priority queue.
         * Always the node with the lowest priority is contracted.
         */
        void eliminateByPQueue()
        {
            VERBOSE( cout << "Eliminatate by priority-queue..." << endl );

            // Save statistics to file
            if (_statsFile != "")
            {
                if (saveStats != SAVE_STATS_NONE) 
                {
                    _stats.open((_statsFile+".contract").c_str());
                    if (!_stats.is_open()) { cerr << "Cannot write to " << (_statsFile+".contract") << endl; exit(1); }

                    if ( saveStats == SAVE_STATS_ALL )
                    {
                        _tree.open((_statsFile+".tree").c_str());
                        if (!_tree.is_open()) { cerr << "Cannot write to " << (_statsFile+".tree") << endl; exit(1); }
                    }
                }

                if (saveShortcutsText) 
                {
                    _shortcutsText.open((_statsFile+".shortcuts").c_str());
                    if (!_shortcutsText.is_open()) { cerr << "Cannot write to " << (_statsFile+".shortcuts") << endl; exit(1); }
                }

            }

            VERBOSE( double timeStart = timestamp(); );
            VERBOSE( double timeLast = timestamp(); );

            unsigned int i = 0;
            VERBOSE( Percent percent(_noOfNodes) );

            // lazy updates counter
            NodeID lazyUpdateCounter = 0;
            // lazy updates counter at the beginnig of the last check interval
            // if there are too many lazy updates during a check interval
            // the whole priority queue is updated
            NodeID lastLazyUpdateCounter = 0;
            
			/*ofstream congraph("congraph.txt", ios::out);*/


            // Main loop: in each loop run, exactly one node is eliminated
            //            so this loop is executed noOfNodes() times.
            while (_pqElimination.min() != EliminationWeight::MAX_VALUE) {

                VERBOSE( 
                    if (i > 0 && _noOfNodes > 0 && i % PROGRESS_NODES_ELIMINATE == 0) 
                    {
                        double now = timestamp();
                        cout << i << " nodes, " << (now-timeStart) << " seconds, " ;
                        cout << (now-timeLast) << " last, ";
                        cout << (((now-timeStart) / (i)) * _noOfNodes) << " remaining, ";
                        cout << "avg degree: " << ((double)_noOfEdges/_noOfNodes);
                        cout << " max hops: " << _maxHops;
                        cout << endl;
                        timeLast = now;
                    }
                )
                
                // lazy update: update min element, only remove it if it is still the min element
                if (_lazyUpdate)
                {

                    // if there where more than one lazy update for each node in the last check
                    // interval (on average), the update the whole priority queue
                    if (i % _weightCalc.lazyUpdateRecalcLimit == 0)
                    {
                        if (lazyUpdateCounter-lastLazyUpdateCounter > _weightCalc.lazyUpdateRecalcLimit)
                        {
                            VERBOSE( cout << "There were too many lazy updates, updating whole priority queue." << endl; )
                            updatePQueue();
                        }
                        lastLazyUpdateCounter = lazyUpdateCounter;
                    }
                    
                    NodeID index = _pqElimination.minElement();
                    NodeID node = pqIndexToNode(index);
                    
                    // recalculate the elimination weight (priority) of the currently
                    // topmost (most unimportant) node, if it does not equals the
                    // stored priority, repeat this step with the now topmost node
                    EliminationWeight::Type oldWeight = _pqElimination.min();
                    EliminationWeight::Type newWeight = updateEliminationWeight(node);
					//cout << "oldweight: " << oldWeight << " newweight: " << newWeight  << " " << (oldWeight == newWeight) << endl;
                    while ( oldWeight != newWeight )
                    {
                        VERBOSE_CONTRACT( cout << "resubmit " << oldWeight << " != " << newWeight << endl; )
                        lazyUpdateCounter++;
                        VERBOSE( if (lazyUpdateCounter % PROGRESS_NODES_ELIMINATE == 0) cout << lazyUpdateCounter << " lazy updates" << endl; )
                            
                        // If the topmost node remains the same but had a different priority than
                        // stored, another update is not necessary. However the lazyUpdateCounter
                        // should be increased, so this condition is not checked at the entrance 
                        // of the loop.
                        if ( index == _pqElimination.minElement() ) break;
                        index = _pqElimination.minElement();
                        node = pqIndexToNode(index);
                        oldWeight = _pqElimination.min();
                        newWeight = updateEliminationWeight(node);
                    }
                }
                
                VERBOSE_CONTRACT( cout << "[" << i << "]"; )
                if (saveStats != SAVE_STATS_NONE && _stats.is_open())
                {
                    _stats << fixed << setprecision(0) << _pqElimination.min() << " ";
                }
                    
                double time1, time2;
                if (saveStats != SAVE_STATS_NONE) time1 = timestamp();
                
                // remove topmost node with lowest priority from priority queue
                NodeID index = _pqElimination.deleteMin();
                NodeID node = pqIndexToNode(index);

				//cout << node << endl;
                
                if (saveStats != SAVE_STATS_NONE) time1 = timestamp() - time1;
                VERBOSE_CONTRACT( cout << " node " << node << flush; )
                VERBOSE_CONTRACT( if (saveStats != SAVE_STATS_NONE) cout << " delete " << fixed << setprecision(3) << time1 * 1000 << flush; )
                if (saveStats != SAVE_STATS_NONE) time2 = timestamp();

                
                unsigned int searchSpace;
                int edgeDiff;
                unsigned int newEdges;
                unsigned int inDegree;
                unsigned int outDegree;
                unsigned int voronoiNumber = pqData(node).voronoiNumber;
                
                // process node, meaning contraction (adds shortcuts),
                // update of level of node and update of priority of neighbors
                processNode<PHASE_NODEORDER_ELIMINATE,false>(node, &searchSpace, &edgeDiff, &newEdges, &inDegree, &outDegree);
			

                if (saveStats != SAVE_STATS_NONE) time2 = timestamp() - time2;
                VERBOSE_CONTRACT( if (saveStats != SAVE_STATS_NONE) cout << " eliminate " << fixed << setprecision(3) << time2 * 1000 << flush; )

                // staged hop-limits: if the average degree exceeds
                // a limit, a new hop-limit is specified and the whole
                // priority queue gets updated.
                if ( updateHopLimit(_weightCalc.maxHops) )
                {
                    VERBOSE( cout << "Switch to max hops " << _maxHops << " after " << (i+1) << " nodes, " << (timestamp()-timeStart) << " seconds" << endl; )
                    updatePQueue();
                }


                // debugging: Save statistics for each contracted node including some
                // priority terms, tempalte parmeter saveStats needs to be enabled.
                if (saveStats != SAVE_STATS_NONE && _stats.is_open())
                {
                    _stats << node << " ";
                    _stats << fixed << setprecision(0) << _noOfNodes << " " << _noOfEdges;
                    _stats << " " << setprecision(8) << (time1+time2) << setprecision(0) << " " << edgeDiff << " " << newEdges;
                    _stats << " " << searchSpace << " " << pqData(node).searchSpace;
                    _stats << " " << pqData(node).deletedNeighbors;
                    _stats << " " << inDegree << " " << outDegree;
                    _stats << " " << voronoiNumber << " " << pqData(node).searchPathHopBorder;
                    _stats << endl;
                }

                _currentLevel++;

                VERBOSE_CONTRACT( cout << endl; )

                // debugging: Perform time consuming tests to ensure correctness
                // of contraction.
                if (_testShortestPaths )
                {
                    assert( _graph->checkReverseGraphExists() );
                    assert( _pqElimination.checkHeapProperty() );
                    testShortestPaths();
                }
                VERBOSE( percent.printStatus(i) );
                i++;
            }

            assert( _graph->checkReverseGraphExists() );
            
            if (_stats.is_open()) _stats.close();
            if (_shortcutsText.is_open()) _shortcutsText.close();
            if (_tree.is_open()) _tree.close();
            
            // statistics: Store witness counter to file. Those counters
            // are in two global vectors and have been filled in processNode().
            if (saveStats == SAVE_STATS_ALL)
            {
                ofstream out((_statsFile+".witness-counter").c_str());
                if (!out.is_open()) { cerr << "Cannot write to " << (_statsFile+".witness-counter") << endl; exit(1); }
                for ( vector<unsigned int>::const_iterator iter = _witnessCounter.begin(); iter != _witnessCounter.end(); iter++ )
                {
                    out << *iter << endl;
                }
                out.close();
                out.open((_statsFile+".witness-search-space-counter").c_str());
                if (!out.is_open()) { cerr << "Cannot write to " << (_statsFile+".witness-search-space-counter") << endl; exit(1); }
                for ( vector<unsigned int>::const_iterator iter = _witnessSearchSpaceCounter.begin(); iter != _witnessSearchSpaceCounter.end(); iter++ )
                {
                    out << *iter << endl;
                }
                out.close();
            }

            VERBOSE( cout << "#edges: " << _graph->noOfExistingEdges() << " / " << _graph->noOfEdges() << endl );
        }
        
        
        /** 
         * Create a contraction hierarchy if the node order (levels) is already given.
         * Basically this is one loop where each node is contracted, ascending
         * by importance (level).
         */
        void constructByLevel()
        {
            VERBOSE( cout << "Eliminate by given level..." << endl );

            // insert all nodes in max level
            // so local search will only regard not already processed nodes
            // this is already done in importGraphListOfEdgesUpdateable
            
            VERBOSE( double timeStart = timestamp(); );
            VERBOSE( double timeLast = timestamp(); );

            VERBOSE( Percent percent(_graph->noOfNodes() ) );
            
            _currentLevel = 0;

            // Prepare hop-limits, since staged hop-limits are used that change
            // after reaching certain average degrees, some overhead is required.
            initHopLimit(_contractParams.maxHops);
            
            // need to init wittness array for fast mtm 2 hops search            
            // if a 2-hop search is in the list of hop-limits
            initPossibleWitnessesMTM(_contractParams.maxHops);

            NodeID k = 0;
            while ( k < _graph->noOfNodes()  && _noOfNodes > _contractParams.coreSize) {
                VERBOSE( 
                    if (k > 0 && k % PROGRESS_NODES_CONSTRUCT == 0) 
                    {
                        double now = timestamp();
                        cout << k << " nodes, " << (now-timeStart) << " seconds, " ;
                        cout << (now-timeLast) << " last, ";
                        cout << ((((now-timeStart) / (k)) * _graph->noOfNodes())-(now-timeStart)) << " remaining, ";
                        cout << "avg degree: " << ((double)_noOfEdges/_noOfNodes);
                        cout << " max hops: " << _maxHops;
                        cout << endl;
                        timeLast = now;
                    }
                )
                
                // contract nodes by level.                
                assert( k < _levelNodeList.size() );
                NodeID u = _levelNodeList[k].second;
                assert( u < _graph->noOfNodes() );
                VERBOSE_CONTRACT( cout << "[" << k << "]" << flush; )
                VERBOSE_CONTRACT( cout << " node " << u << flush; )
                VERBOSE_CONTRACT( double time1 = timestamp(); )
                
                // *** Main step: Node contraction ***
                processNode<PHASE_CONSTRUCT, false>(u);
                
                VERBOSE_CONTRACT( time1 = timestamp() - time1; )
                VERBOSE_CONTRACT( cout << " eliminate " << fixed << setprecision(3) << time1 * 1000 << flush; )
                VERBOSE_CONTRACT( cout << endl; )

                // staged hop-limits: if the average degree exceeds
                // a limit, a new hop-limit is specified.
                if ( updateHopLimit(_contractParams.maxHops) )
                {
                    VERBOSE( cout << "Switch to max hops " << _maxHops << " after " << (k+1) << " nodes, avg degree: " << ((double)_noOfEdges/_noOfNodes) << ", " << (timestamp()-timeStart) << " seconds" << endl; )
                }
                
                k++;
                
                _currentLevel++;
                
                VERBOSE( percent.printStatus(k) );
            }

            assert( _graph->checkReverseGraphExists() );
            VERBOSE( cout << "#edges: " << _graph->noOfExistingEdges() << " / " << _graph->noOfEdges() << endl );
            VERBOSE( cout << "max avg degree: " << _maxAvgDegree << endl; )
        }
        
        /** 
         * Create a contraction hierarchy if the node order (levels) is already given
         * and additionally all shortcuts and witness paths used during the node ordering.
         * The idea was to use this information to create a hierarchy for a graphs
         * with slightly changed edge weights. But this implementation did not works
         * that well and is only preliminary.
         * There are currently two algorithms implemented:
         *   1) Add shortcuts and add a "witness shortcut" edge for each witness path. Then perform
         *      the contraction, the witness path can then be "found" by a 1-hop search.
         *      The shortcuts are also added.
         *   2) Use the stored witness paths and shortcuts on the fly during contraction
         *      to omit local searches if a valid witness exists.
         */
        void constructByLevelRetrieve(istream& inShortcuts, istream& inWitnesses)
        {
            VERBOSE( cout << "Eliminatate by given level..." << endl );
            VERBOSE( cout << "Shortcuts and Witnesses are also given." << endl );

            // insert all nodes in max level
            // so local search will only regard not already processed nodes
            // this is already done in importGraphListOfEdgesUpdateable
            
            VERBOSE( double timeStart = timestamp(); );
            VERBOSE( double timeLast = timestamp(); );

            VERBOSE( double now; )

    	    RetrieveWitnesses rWitnesses(inWitnesses);

            // First Algorithm:
            //    Add shortcuts and add a "witness shortcut" edge for each witness path. Then perform
            //    the contraction, the witness path can then be "found" by a 1-hop search.
            if ( !onTheFlyWitnessCheck )
            {

                VERBOSE( cout << "Add shortcuts..." << endl );
                addShortcuts(inShortcuts);
                VERBOSE( now = timestamp(); )
                VERBOSE( cout << (now-timeStart) << " seconds, " << (now-timeLast) << " last" << endl; )
                VERBOSE( timeLast = now; )
    
                VERBOSE( cout << "Add witness-shortcuts..." << endl );
                addWitnessShortcuts(rWitnesses);
                VERBOSE( now = timestamp(); )
                VERBOSE( cout << (now-timeStart) << " seconds, " << (now-timeLast) << " last" << endl; )
                VERBOSE( timeLast = now; )
                VERBOSE( cout << "Validate witnesses..." << endl );
            }   
            else
            {             
                VERBOSE( cout << "Check witnesses..." << endl );
            }

            // two hop search for new shortcuts            
            _possibleWitnesses.clear();
            _possibleWitnesses.resize(_graph->noOfNodes(), Weight::MAX_VALUE);

            VERBOSE( Percent percent(_graph->noOfNodes() - 2) );
            
            _currentLevel = 0;
            
            NodeID k = 0;
            while ( k < _graph->noOfNodes() ) {
                VERBOSE( 
                    if (k > 0 && k % PROGRESS_NODES_CONSTRUCT == 0) 
                    {
                        now = timestamp();
                        cout << k << " nodes, " << (now-timeStart) << " seconds, " ;
                        cout << (now-timeLast) << " last, ";
                        cout << ((((now-timeStart) / (k)) * _graph->noOfNodes())-(now-timeStart)) << " remaining, ";
                        cout << "avg degree: " << ((double)_noOfEdges/_noOfNodes);
                        cout << endl;
                        timeLast = now;
                    }
                )
                
                NodeID u = _levelNodeList[k].second;
                VERBOSE_CONTRACT( cout << "[" << k << "]" << flush; )
                VERBOSE_CONTRACT( cout << " node " << u << flush; )
                VERBOSE_CONTRACT( double time1 = timestamp(); )
                
                
                // Shortcuts and "witness shortcuts" are already added to the graph.
                // Only a limited local search is necessary.
                if ( !onTheFlyWitnessCheck )
                {
                    validateWitnesses(u);
                }
                // Perform contraction but try to avoid local searches and use 
                // the stored information (shortcuts, witnesses) instead.
                else
                {
                    checkWitnesses(u, rWitnesses);
                }
                
                _graph->node(u).setLevel(_currentLevel);
                _graph->changeNodeLevelOnlyReverseEdges(u);
                VERBOSE_CONTRACT( time1 = timestamp() - time1; )
                VERBOSE_CONTRACT( cout << " eliminate " << fixed << setprecision(3) << time1 * 1000 << flush; )
                VERBOSE_CONTRACT( cout << endl; )

                k++;
                _currentLevel++;
                
                VERBOSE( percent.printStatus(k) );
            }
            assert( _graph->checkReverseGraphExists() );
            
            VERBOSE( cout << "#edges: " << _graph->noOfExistingEdges() << " / " << _graph->noOfEdges() << endl );
            
        }

        /**
         * Initialize staged hop-limits. The currently active hop-limit
         * is stored in _maxHops. The update of the hop-limit, triggered
         * by average-degree-limits, is performed in updateHopLimit().
         */
        void initHopLimit(const vector<double>& maxHops)
        {
            _maxAvgDegree = 0;
            _maxHops = 0;
            _maxHopsIndex = 0;
            _maxHopsDegreeLimit = 0;
            if ( _maxHopsIndex < maxHops.size() )
            {
                _maxHops = (unsigned int)maxHops[_maxHopsIndex];
                _maxHopsIndex++;
                
                if ( _maxHopsIndex < maxHops.size() )
                {
                    _maxHopsDegreeLimit = maxHops[_maxHopsIndex];
                    _maxHopsIndex++;
                }
            }
        }
        
        /**
         * Performs a hop-limit switch based on average-degree-limits. It also keeps track of the 
         * maximum average degree during contraction.
         * @return if a switch occured.
         */
        bool updateHopLimit(const vector<double>& maxHops)
        {
            // staged hop-limits: if the average degree exceeds
            // a limit, a new hop-limit is specified.
            if (_noOfNodes > 0)
            {
                double currentAvgDegree = (double)_noOfEdges/_noOfNodes;
                if (_maxAvgDegree < currentAvgDegree) _maxAvgDegree = currentAvgDegree;

                if (_maxHopsIndex < maxHops.size() && currentAvgDegree >= _maxHopsDegreeLimit)
                {
                    _maxHops = (unsigned int)maxHops[_maxHopsIndex];
                    _maxHopsIndex++;
                    
                    if ( _maxHopsIndex < maxHops.size() )
                    {
                        _maxHopsDegreeLimit = maxHops[_maxHopsIndex];
                        _maxHopsIndex++;
                    }
                    return true;
                }
            }
            return false;
        }

        // ***
        // The next structure, variables and procedures with MTM in their name
        // are for the local 2-hop search for witnesses. They implement a simple
        // many-to-many algorithm, see the diploma thesis for a theoretical
        // description.
        // ***
        
        /**
        * The bucket entries of a node are stored as a linked list in an array (vector).
        * It contains the node and the distance from this node. The target node
        * is known as the owner of the bucket.
        */
        struct LinkedNodeMTM
        {
            NodeID nodeID;
            EdgeWeight weight;
            NodeID next; // next index
        };

        // array of buckets, stored as linked lists
        // It is filled in initBucketsMTM() and read
        // in findPossibleWitnessesMTM(). The first
        // index to the bucket b(x) of node x
        // is stored in the pqElement attribute
        // of the node in the graph datastructure.
        vector<LinkedNodeMTM> _linkedListMTM;

        // index is node id, value is distance to this node id
        // It is filled in findPossibleWitnessesMTM() and concerns
        // the node v that is specified in the call of this procedure.
        vector<EdgeWeight> _possibleWitnesses;

        /**
         * Initalize possible witnesses array.
         */
        void initPossibleWitnessesMTM(const vector<double>& maxHops)
        {
            for ( unsigned int i = 0; i < maxHops.size(); i += 2 )
            {
                if ( maxHops[i] == 2 )
                {
                    _possibleWitnesses.clear();
                    _possibleWitnesses.resize(_graph->noOfNodes(), Weight::MAX_VALUE);
                    break;
                }
            }
        }

        /**
        * First step of the 2-hop many-to-many search. The incoming edges (x,v)
        * from each node v that is incident to a outgoing edge (node,v) of the 
        * currently processed node are scanned and the distance d(x,v) is stored
        * in the bucket b(x) along with v. The currently processed node (node)
        * is ignored.
        * @param node currently processed node
        * @param node firstEdge start index into edge array for node
        * @param node lastEdge last index+1 into edge array for node
        */
        void initBucketsMTM(const NodeID node, const OtherEdgeID firstEdge, const OtherEdgeID lastEdge)
        {
            assert( _linkedListMTM.empty() );

            // one step backward search from all targets of the outgoing edges
            for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
            {
                const Edge& edgeOut = _graph->edge(eOut);
                // only outgoing edges of the currently processed node
                if (  !edgeOut.isDirected(0) ) continue;
                NodeID v = edgeOut.target();

                // scan all incoming edges (x,v) of node v and store
                // the distance to v in the bucket b(x) of x.
                OtherEdgeID vLastEdge = _graph->lastEdge(v);
                for ( OtherEdgeID e = _graph->firstLevelEdge(v); e < vLastEdge; e++ )
                {
                    const Edge& edge = _graph->edge(e);
                    if ( !edge.isDirected(1) || edge.target() == node ) continue;
                    _linkedListMTM.push_back( LinkedNodeMTM() );
                    LinkedNodeMTM& link = _linkedListMTM.back();
                    link.nodeID = v;
                    link.weight = edge.weight();
                    link.next = _graph->node(edge.target()).pqElement();
                    _graph->node(edge.target()).pqElement(_linkedListMTM.size());
                }
            }
        }

        /**
        * Second step of the 2-hop many-to-many search. The buckets
        * of the 1-hop backward search are now filled. We now perform
        * a 1-hop forward search (equals a edge array scan) from a
        * node v that is incident to an incoming edge (v,node). We use
        * the buckets to find the length of 2-hop witnesses.
        * The scan of the bucket b(v) gives the length 1-hop witnesses.
        * Since we are interested in the shortest witness to each node,
        * we simply use a vector to store the shortest length to each
        * reached node. The currently processed node (node) is ignored.
        * @param node currently processed node
        * @param v node incident to an incoming edge (v,node)
        */
        void findPossibleWitnessesMTM(const NodeID node, const NodeID v)
        {
            assert(_possibleWitnesses.size() >= _graph->noOfNodes());
            // Check wheter the array storing the lengths of the
            // witnesses is clean.
            #ifndef NDEBUG
            OtherEdgeID firstEdge = _graph->firstLevelEdge(node);
            OtherEdgeID lastEdge = _graph->lastEdge(node);
            for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
            {
                const Edge& edgeOut = _graph->edge(eOut);
                // only outgoing edges
                if (  !edgeOut.isDirected(0) ) continue;
                assert( _possibleWitnesses[edgeOut.target()] == Weight::MAX_VALUE );
            }
            #endif

            // Scan the bucket b(v) to find the length
            // of one-hop witnesses.
            NodeID index = _graph->node(v).pqElement();
            while (index != 0)
            {
                assert( index <= _linkedListMTM.size() );
                LinkedNodeMTM& link = _linkedListMTM[index-1];
                if (_possibleWitnesses[link.nodeID] > link.weight)
                {
                    _possibleWitnesses[link.nodeID] = link.weight;
                }
                index = link.next;
            }

            // Perform a 1-hop forward search (edge array scan)
            // and scan the buckets of the reached nodes
            // to find the length of 2-hop witnesses.
            OtherEdgeID vLastEdge = _graph->lastEdge(v);
            for ( OtherEdgeID e = _graph->firstLevelEdge(v); e < vLastEdge; e++ )
            {
                const Edge& edge = _graph->edge(e);
                if ( !edge.isDirected(0) || edge.target() == node ) continue;
                index = _graph->node(edge.target()).pqElement();
                while ( index != 0 )
                {
                    assert( index <= _linkedListMTM.size() );
                    LinkedNodeMTM& link = _linkedListMTM[index-1];
                    if (link.nodeID != v)
                    {
                        EdgeWeight newWeight = link.weight + edge.weight();
                        if (_possibleWitnesses[link.nodeID] > newWeight)
                        {
                            _possibleWitnesses[link.nodeID] = newWeight;
                        }
                    }
                    index = link.next;
                }
            }
        }

        /**
         * Cleanup step of the 2-hop many-to-many search. The bucket entries
         * are cleared and the start pointers stored in the graph data structure
         * to the buckets are removed. Similar to a clear() in the DijkstraCH class.
         * @param node currently processed node
         * @param node firstEdge start index into edge array for node
         * @param node lastEdge last index+1 into edge array for node
         */
        void clearMTM(const NodeID node, const OtherEdgeID firstEdge, const OtherEdgeID lastEdge)
        {
            for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
            {
                const Edge& edgeOut = _graph->edge(eOut);
                if ( !edgeOut.isDirected(0) ) continue;
                NodeID v = edgeOut.target();
                OtherEdgeID vLastEdge = _graph->lastEdge(v);
                for ( OtherEdgeID e = _graph->firstLevelEdge(v); e < vLastEdge; e++ )
                {
                    const Edge& edge = _graph->edge(e);
                    if ( !edge.isDirected(1) || edge.target() == node ) continue;
                    _graph->node(edge.target()).pqElement(0);
                }
            }
            _linkedListMTM.clear();                    
        }
        
        /**
         * Stores a new shortcut edge in the buffer for new edges. These adges are added to the graph
         * after all shortcuts for the contraction of a node are found.
         */
        void prepareNewShortcutEdge(const NodeID node, const OtherEdgeID firstEdge, const OtherEdgeID lastEdge, 
            const OtherEdgeID eIn, const Edge& edgeIn, const OtherEdgeID eOut, const Edge& edgeOut)
        {
            // On bidirectional shortcuts, only insert one bidirectional
            // shortcut instead of two unidirectional ones. The check
            // uses the fact that the adjacency array of node is sorted.
            bool bidir = edgeIn.isBidirected() && edgeOut.isBidirected();
            bool b = !bidir || edgeIn.target() < edgeOut.target(); // adjacency array sorted
            if (!b)
            {
                b = true;
                for ( typename vector<NewEdge>::iterator iter = _newEdges.begin(); 
                      b && iter != _newEdges.end(); iter++ )
                {
                    if (iter->source == edgeOut.target() 
                        && iter->edge.target() == edgeIn.target()
                        && iter->edge.weight() == (edgeIn.weight() + edgeOut.weight()))
                    {
                        iter->edge.makeTwoWay();
                        b = false; // second shortcut not necessary
                    }
                }
            }
            // Add shortcut edge. All shortcut edges required for the
            // contraction of node are saved in an array and later added
            // at once.
            if (b)
            {
                // Prepare shortcut edge, information of path expansion are stored
                // if supported by the edge datastructure. And also the
                // number of original edges a shortcut represent, is stored
                // using the attribute shortcutOriginalEdgeCount.
                _newEdges.push_back( NewEdge( edgeIn.target(), edgeOut.target(), 
                        edgeIn.weight() + edgeOut.weight(), false /* unidirectional */,
                        node, eIn-firstEdge, eOut-firstEdge,
                        edgeIn.shortcutOriginalEdgeCount() + edgeOut.shortcutOriginalEdgeCount() ) );
            }
        }
        
        /**
         * Perform a local edge reduction using the information stored in the Dijkstra object
         * used for local searches. Only edges (v,x) that are not on any shortest path
         * since their weight is larger than the shortest path distance between v and x, 
         * are removed.
         */
        void reduceEdgesLocal(const NodeID node, const Edge& edgeIn)
        {
            // reduce edges, process edges descending because otherwise it 
            // will cause problems at removal (rearrangeing of edges)
            OtherEdgeID rE = _graph->lastEdge(edgeIn.target());
            if (rE > 0)
            {
                // Scan through all outgoing edges of the start node of the local search.
                // If the incident node is settled and the distance found by the local
                // search is smaller than the edge weight, remove the edge.
                // Note: the edges are processed descending, this makes it easier to
                // directly remove the edges without causing index problems.
                OtherEdgeID rFirstEdge = _graph->firstLevelEdge(edgeIn.target());
                for ( rE--; rE >= rFirstEdge; rE-- )
                {
                    Edge& rEdge = _graph->edge(rE);
                    if ( rEdge.isDirected(0) && rEdge.target() != node )
                    {
                        if (_localDijkstra.isSettled(0, rEdge.target() ) && _localDijkstra.distanceTo(rEdge.target(), 0) < rEdge.weight())
                        {
                            OtherEdgeID rReverseE = _graph->reverseLevelEdge(edgeIn.target(), rE);
                            // Special case: need to consider bidirectional flags. We can unly
                            // remove unidirectional edges, but bidirectional edges becomes
                            // one-way. Also the reverse edge in the graph data structure
                            // needs to be considered since each edge is stored in the 
                            // adjacency array of both incident nodes.
                            if ( rEdge.isBidirected() )
                            {
                                rEdge.makeOneWay(1);
                                _graph->edge(rReverseE).makeOneWay(0);
                            }
                            else
                            {
                                NodeID reverseNodeID = rEdge.target();
                                OtherEdgeID reverseOtherEdgeID = _graph->reverseLevelEdge( edgeIn.target(), rE );
                                _graph->removeEdge(edgeIn.target(), rE);
                                _graph->removeEdge(reverseNodeID, reverseOtherEdgeID);
                                _noOfEdges -= 2;
                            }
                        }
                    }
                    if (rE == 0) break;                                
                }
            }
        }
         
        /**
         * Store a witness path to file to use it in a later hiearchy construction.
         * Also some statistics for witnesses are accumulated if a template
         * parameter is enabled.
         */   
        void storeWitness(const NodeID node, const Edge& edgeIn, const Edge& edgeOut)
        {
            if ( (saveStats == SAVE_STATS_ALL || saveShortcutsWitnesses) )
            {
                
                NodeID current = edgeOut.target();
                const vector<NodeID>& settledNodes = _localDijkstra.settledNodes(0);
                
                if (saveStats == SAVE_STATS_ALL)
                {
                    NodeID index = 0;
                    while (settledNodes[index] != current)
                    {
                        index++;
                        assert( index < settledNodes.size() );
                    }
                    if (_witnessSearchSpaceCounter.size() <= index) 
                    {
                        _witnessSearchSpaceCounter.resize(index, 0);
                    }
                    assert( index > 0 );
                    _witnessSearchSpaceCounter[index-1]++;
                }
                
                if (saveStats == SAVE_STATS_ALL || saveShortcutsWitnesses)
                {
                    unsigned int counter = 0;
                    while (current != edgeIn.target())
                    {
                        NodeID next = _localDijkstra.parentOf(current, 0);
                        if (_tree.is_open())
                        {
                            _tree << edgeIn.target() << " " << next << " " << current << endl;
                        }
                        if ( saveShortcutsWitnesses && next != edgeIn.target() )
                        {
                            // do not save witnesses of hop length 1
                            if ( current == edgeOut.target() )
                            {
                                _storeWitnesses.beginWitness(node, edgeIn.target(),edgeOut.target());
                            }
                            _storeWitnesses.nextStep(next);
                        }
                        current = next;
                        counter++;
                    }
                    if (saveStats == SAVE_STATS_ALL)
                    {
                        if (_witnessCounter.size() <= counter) 
                        {
                            _witnessCounter.resize(counter, 0);
                        }
                        assert( counter > 0 );
                        _witnessCounter[counter-1]++;
                    }
                }
            }
        }
        
        /**
         * After the contraction of node, we usually update all neighbors of it.
         * This is only a heuristic but works quite well. Before we can actually
         * recalculate the elimination weight (priority), we need to update some
         * attributes 
         *   - the deleted neighbor counter
         *   - the upper limit on search path lengts
         *   - distribute the Voronoi region of the contracted nodes among the
         *     neighboring regions.
         * For testing: In case of hop-limits, we can update all affected nodes.
         */
        void updateAfterContraction(const NodeID node, const OtherEdgeID firstEdge, const OtherEdgeID lastEdge)
        {
            // Update the deleted neibhor counter and the upper bound on search paths hops.
            // Need to do this before update of elimination weights.
            NodeID previous = SPECIAL_NODEID;
            if (!_weightCalc.searchPathHopBorderOriginalEdges)
            {
                // possibly new search path hops = upper bound of current node + 1 
                double newsearchPathHopBorder = pqData(node).searchPathHopBorder + 1;
                assert( newsearchPathHopBorder >= pqData(node).searchPathHopBorder );
                for ( OtherEdgeID e = firstEdge; e < lastEdge; e++ )
                {
                    const Edge& edge = _graph->edge(e);
                    if (edge.target() == previous) continue;
                    previous = edge.target();

                    pqData(edge.target()).deletedNeighbors++;
                    
                    // increase upper search path hop bound, if necessary
                    if (pqData(edge.target()).searchPathHopBorder < newsearchPathHopBorder)
                    {
                        pqData(edge.target()).searchPathHopBorder = newsearchPathHopBorder;
                    }
                }
            }
            
            // This is an upper bound for the costs of the search paths by a time dependend search.
            // Requested by Veit Batz.
            else
            {
                for ( OtherEdgeID e = firstEdge; e < lastEdge; e++ )
                {
                    const Edge& edge = _graph->edge(e);
                    if (edge.target() == previous) continue;
                    previous = edge.target();
                    // we ignore the direction of the edge and use just the first edge
                    double newsearchPathHopBorder = 2*pqData(node).searchPathHopBorder + edge.shortcutOriginalEdgeCount();
                    assert( newsearchPathHopBorder >= pqData(node).searchPathHopBorder );

                    pqData(edge.target()).deletedNeighbors++;
                    if (pqData(edge.target()).searchPathHopBorder < newsearchPathHopBorder)
                    {
                        pqData(edge.target()).searchPathHopBorder = newsearchPathHopBorder;
                    }
                }
            }
            
            // Update the Voronoi regions, meaning distribute the nodes in the Voronoi region of
            // the current node among the the neighboring Voronoi regions.
            if (_weightCalc.voronoiMult != 0 && _noOfNodes > 1)
            {
                distributeVoronoiRegion(node);
            }

            VERBOSE_CONTRACT( cout << " (" << flush; )
            
            // Update weights of all neigbors.
            if ( !_weightCalc.updateHops )
            {
                // do not update a neighbor twice if there is a separate
                // incoming and outgoing edge to this neighbor.
                // The edges in the edge array are sorted by the node id
                // of the incident node.
                previous = SPECIAL_NODEID; 
                for ( OtherEdgeID e = firstEdge; e < lastEdge; e++ )
                {
                    const Edge& edge = _graph->edge(e);
                    if (edge.target() == previous) continue;
                    previous = edge.target();

                    VERBOSE_CONTRACT( cout << " node " << edge.target() << flush; )
                    updateEliminationWeight(edge.target());
                }
            }
            // Test case: update all nodes within the hop-limit.
            // This is quite time consuming and shows only little
            // impact on the resulting hierarchy.
            // A BFS is used to find all nodes, the direction of the
            // edges is ignored.
            else
            {
                assert( _maxHops > 0 );
                stack< pair<NodeID, NodeID> > bfs;
                vector<NodeID> updatedNodes;
                bfs.push( make_pair( node, 0 ) );
                NodeID maxHops = _maxHops;
                // Special case: 1-hop search can affect nodes 2 hops away
                if ( maxHops == 1 ) maxHops = 2;
                while ( !bfs.empty() )
                {
                    NodeID v = bfs.top().first;
                    NodeID hops = bfs.top().second + 1;
                    bfs.pop();
                    
                    NodeID vLastEdge = _graph->lastEdge(v);
                    for ( OtherEdgeID e = _graph->firstLevelEdge(v); e < vLastEdge; e++ )
                    {
                        const Edge& edge = _graph->edge(e);
                        if ( !_updateBfs[edge.target()] )
                        {
                            VERBOSE_CONTRACT( cout << " node " << edge.target() << flush; )
                            unsigned int singleSearchSpace = 0;
                            updateEliminationWeight(edge.target(), &singleSearchSpace);
                            //searchSpace += singleSearchSpace;
                            updatedNodes.push_back(edge.target());
                            _updateBfs[edge.target()] = true;
                            
                            if ( hops < maxHops )
                            {
                                bfs.push( make_pair( edge.target(), hops ) );   
                            }
                        }
                    }
                }
                for ( vector<NodeID>::const_iterator iter = updatedNodes.begin(); iter != updatedNodes.end(); iter++ )
                {
                    _updateBfs[*iter] = false;   
                }
            }
        }
            

        /**
        * Main procedure to contract a single node. Contracting a node means
        * adding shortcut edges so that shortest path distances are preserved
        * after the node is removed. The necessity of shortcuts is checked
        * by searching for witness paths that show an alternative path not
        * including node. There are currently four local searches with significant
        * implemenation differences:
        *  - 1-hop search
        *  - 2-hop search using many-to-many
        *  - local Dijkstra search with 1-hop backward search
        *  - local Dijkstra search
        * Their implementation differs many in the search for witnesses but
        * their decision to add a shortcut is always based on the distance
        * found between nodes incident to incoming and outgoing edges. And the
        * search always ignores the currently processed node.
        * @param phase phase, PHASE_NODEORDER_INIT, PHASE_NODEORDER_ELIMINATE, 
        *              PHASE_CONSTRUCT possible, during PHASE_NODEORDER_INIT
        *              is no access to the priority queue allowed since it just gets filled.
        * @param simulateOnly only simulate contraction for weight calculation, see 
        *                     calculateEliminationWeight()
        * @param node currently processed node
        * The additional parameters are return values, that are filled 
        * if != NULL (call-by-reference).
        */
        template < ConstructPhase phase, bool simulateOnly >
        void processNode(const NodeID node, unsigned int* searchSpaceResult = NULL, 
            int* edgeDiffResult = NULL, unsigned int* newEdgesResult = NULL, 
            unsigned int* inDegreeResult = NULL, unsigned int* outDegreeResult = NULL,
            unsigned int* shortcutOriginalEdgeSum = NULL)
        {
            // range in the edge array, for each node the edges are
            // partitioned into edges where the second incident edge,
            // called target independent of the direction of the edge,
            // is already contracted/eliminated or not.
            //      firstEdge
            //          .. 
            //        (edges to contracted nodes) 
            //          .. 
            //      firstLevelEdge
            //          .. 
            //        (edges to remaining nodes) 
            //          .. 
            //      lastEdge
            OtherEdgeID firstEdge = _graph->firstLevelEdge(node); 
            OtherEdgeID lastEdge = _graph->lastEdge(node);
			//if (node == 0)
			

            // We sort the edges so we can process them more efficiently. This
            // is e.g. used to add bidirectional flags to shortcuts, if
            // they exists in both directions with the same weight.
            if (firstEdge < lastEdge)
            {
                _graph->sortEdges(firstEdge, lastEdge);
            }

			//ofstream congraph("graph.txt", ios::out);
			//congraph << _graph->noOfNodes() << " "  << _graph->noOfEdges() << endl;
		
			/*if (phase == PHASE_NODEORDER_ELIMINATE)
			{
				cout << endl;
				cout << "first : last " << firstEdge << " " << lastEdge << endl;
				for ( OtherEdgeID e = firstEdge; e < lastEdge; e++ )
				{
					const Edge& edge = _graph->edge(e);
					int dir = 3;
					if (edge.isDirected(0)) dir -= 2;
					if (edge.isDirected(1)) dir -= 1;
					cout << node << " " << edge.target() << " " << edge.weight() << " " << dir << endl;
				}
			}*/
		

            if ( phase == PHASE_NODEORDER_INIT || phase == PHASE_NODEORDER_ELIMINATE )
            {
                // calculate the incoming degree and outgoing degree
                // of the node (always wrt to the remainig graph)
                // used for weight (priority) calculation
                if (inDegreeResult != NULL || outDegreeResult != NULL)
                {
                    unsigned int inDegree = 0;
                    unsigned int outDegree = 0;
                    for (OtherEdgeID edgeID  = firstEdge; edgeID < lastEdge; edgeID++) {
                        const Edge& edge = _graph->edge(edgeID);
                        if (edge.isDirected(0)) outDegree++;
                        if (edge.isDirected(1)) inDegree++;
                    }
                    if ( inDegreeResult    != NULL ) *inDegreeResult    = inDegree;
                    if ( outDegreeResult   != NULL ) *outDegreeResult   = outDegree;
                }
            }

            // the edge difference is implemented as the difference
            // in the number of entries in the edge array.
            // Since each remaining edge is still stored twice,
            // at both incident edges, multiply by 2.
            // Initially, the edges incident to the currently
            // processed nodes count negative. Each required shortcut
            // increases the difference by 2.
            int edgeDiff = -2*(lastEdge-firstEdge);
            // New edge counter, similar to edgeDiff.
            unsigned int newEdgesCounter = 0;
            // Search space, a priority term. Depends on the
            // implementation of the local searches.
            size_t searchSpace = 0;

            // Only act if the current node has remainig edges.
            if (firstEdge < lastEdge)
            {

                // debug: Check for parallel edges, these must not exist.
                // Usually, this is ensured by UpdateableGraph::addShortcutEdge().
                #ifndef NDEBUG
                NodeID previousIn = SPECIAL_NODEID;
                NodeID previousOut = SPECIAL_NODEID;
                for ( OtherEdgeID e = firstEdge; e < lastEdge; e++ )
                {
                    const Edge& edge = _graph->edge(e);
                    if ( edge.isDirected(1) )
                    {
                        if ( previousIn == edge.target() )
                        {
                            cout << "node " << node << endl;
                            for (OtherEdgeID edgeID  = firstEdge; edgeID < lastEdge; edgeID++)
                            {
                                cout << _graph->edge(edgeID) << endl;
                            }
                        }
                        assert( previousIn != edge.target() );
                        previousIn = edge.target();
                    }
                    if ( edge.isDirected(0) )
                    {
                        if ( previousOut == edge.target() )
                        {
                            cout << "node " << node << endl;
                            for (OtherEdgeID edgeID  = firstEdge; edgeID < lastEdge; edgeID++)
                            {
                                cout << _graph->edge(edgeID) << endl;
                            }
                        }
                        assert( previousOut != edge.target() );
                        previousOut = edge.target();
                    }
                }
                #endif

                // ***
                // The next if blocks specify the local searches because
                // there are for siginficantly different implementations:
                //   - 1-hop search
                //   - 2-hop search using many-to-many
                //   - local Dijkstra search with 1-hop backward search
                //   - local Dijkstra search
                // ***
                
                // 1-hop search
                // ------------
                // This is a simple scan of the adjacency array.
                if (_maxHops == 1)
                {
                    // The search is implemented as a triply nested loop.
                    // First loop:  Scan edge array of node for nodes incident to 
                    //              incoming edges (v,node)
                    // Second loop: Scan edge array of node for nodes incident to
                    //              outgoing edges (node,w)
                    // Third loop:  Scan edge array of v to find an edge to w.
                    // A witness path is found, if its length is not longer
                    // than the path <v,node,w>.
                    for ( OtherEdgeID eIn = firstEdge; eIn < lastEdge; eIn++ )
                    {
                        const Edge& edgeIn = _graph->edge(eIn);
                        if ( !edgeIn.isDirected(1) ) continue;

                        OtherEdgeID inFirstEdge = _graph->firstLevelEdge(edgeIn.target());
                        OtherEdgeID inLastEdge = _graph->lastEdge(edgeIn.target());
                        
                        if ( phase == PHASE_NODEORDER_INIT || phase == PHASE_NODEORDER_ELIMINATE )
                        {
                            searchSpace += inLastEdge-inFirstEdge;
                        }

                        for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
                        {
                            const Edge& edgeOut = _graph->edge(eOut);
                            if ( !edgeOut.isDirected(0) ) continue;

                            // A path to the start node is not necessary.
                            if ( edgeIn.target() == edgeOut.target() ) continue;

                            bool foundWitness = false;
                            for ( OtherEdgeID e = inFirstEdge; e < inLastEdge; e++ )
                            {
                                const Edge& edge = _graph->edge(e);
                                if ( edge.isDirected(0) && edge.target() == edgeOut.target() )
                                {
                                    if ( edge.weight() <= (edgeIn.weight() + edgeOut.weight()) )
                                    {
                                        foundWitness = true;
                                    }
                                    break;
                                }
                            }

                            // If no witness path is found, after the contraction of node,
                            // shortest paths distance would possibly increase. The countermeasure
                            // is to add a shortcut edge with the length of the shortest path.
                            // Since the search is limted, sometimes shortcuts are added, 
                            // that are not necessary. But they do not invalidate the correctness.
                            if ( !foundWitness )
                            {
                                prepareNewShortcutEdge( node, firstEdge, lastEdge, eIn, edgeIn, eOut, edgeOut );
                            }

                            if ( phase == PHASE_NODEORDER_INIT || phase == PHASE_NODEORDER_ELIMINATE )
                            {
                                // Count each scan of the adjacency array as 1 search space unit.
                                // This works fine. The search space is not so relevant for
                                // the 1-hop search since at least in roat networks the degree
                                // of nodes is small.
                                searchSpace++;
                            }
                        }
                    }
                }

                // 2-hop search
                // ------------
                // The 2-hop search is implemented as a simplified many-to-many search
                // from all nodes incident to incoming edges to all nodes incident to outgoing
                // edges. It uses additional procedures.
                else if (_maxHops == 2)
                {
                    // 1-hop backward search from nodes incident to outgoing edges to fill buckets.
                    initBucketsMTM(node, firstEdge, lastEdge);

                    // For each node incident to an incoming edge, perform a
                    // 1-hop forward search using the previously initalized buckets
                    // to find distances (ommiting current node (node)) to
                    // nodes 1 or 2 hops away. The distances are stored
                    // in the array _possibleWitnesses, index is node id.
                    for ( OtherEdgeID eIn = firstEdge; eIn < lastEdge; eIn++ )
                    {
                        const Edge& edgeIn = _graph->edge(eIn);
                        if ( !edgeIn.isDirected(1) ) continue;

                        // Fill _possibleWitnesses vector.
                        findPossibleWitnessesMTM(node, edgeIn.target());
                        
                        if ( phase == PHASE_NODEORDER_INIT || phase == PHASE_NODEORDER_ELIMINATE )
                        {
                            searchSpace += _graph->lastEdge(edgeIn.target())-_graph->lastEdge(edgeIn.target());
                        }

                        // Check the distance to each node incident to an outgoing edge.
                        for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
                        {
                            const Edge& edgeOut = _graph->edge(eOut);
                            if ( !edgeOut.isDirected(0) ) continue;

                            // A path to the start node is not necessary.
                            if ( edgeIn.target() == edgeOut.target() ) continue;

                            // Get length of shortest witness.
                            // witnessWeight == Weight::MAX_VALUE if no witness exists
                            EdgeWeight witnessWeight = _possibleWitnesses[edgeOut.target()];
                            _possibleWitnesses[edgeOut.target()] = Weight::MAX_VALUE;

                            // If shortest path would have increased length, a shortcut edge may be necessary.
                            // Since the search is limted, sometimes shortcuts are added, 
                            // that are not necessary. But they do not invalidate the correctness.
                            if ( witnessWeight > (edgeIn.weight() + edgeOut.weight()))
                            {
                                prepareNewShortcutEdge( node, firstEdge, lastEdge, eIn, edgeIn, eOut, edgeOut );
                            }
                        }
                    }
                    
                    if ( phase == PHASE_NODEORDER_INIT || phase == PHASE_NODEORDER_ELIMINATE )                    
                    {
                        searchSpace += _linkedListMTM.size();
                    }
                    
                    // cleanup the used mtm data structures
                    clearMTM(node, firstEdge, lastEdge);
                }

                // local dijkstra search, with 1-hop backward search
                // -------------------------------------------------
                // The 1-hop backward search is only used if hop-limits are specified.
                // The current hop-limit is stored in _maxHops.
                else if ( oneHopBackwardSearch && _maxHops > 0 )
                {
                    
                    // Prepare target-flags for the local search. The search is stopped
                    // if all targets are settled. These target flags need to be removed
                    // after the local searches.
                    NodeID noOfTargets = 0;
                    for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
                    {
                        const Edge& edgeOut = _graph->edge(eOut);
                        if ( !edgeOut.isDirected(0) ) continue;

                        // Since we perform an additional 1-hop-backward search,
                        // we do use the nodes adjacent to the outgoing node
                        // as targets.
                        OtherEdgeID outFirstEdge = _graph->firstLevelEdge(edgeOut.target());
                        OtherEdgeID outLastEdge = _graph->lastEdge(edgeOut.target());
                        for ( OtherEdgeID e = outFirstEdge; e < outLastEdge; e++ )
                        {
                            const Edge& edge = _graph->edge(e);
                            if ( !edge.isDirected(1) || edge.target() == node ) continue;
                            if (!_graph->node(edge.target()).isTarget())
                            {
                                _graph->node(edge.target()).setTarget(true);
                                noOfTargets++;
                            }
                        }
                    }
    
                        
                    // A local search starting at each node incident to an incoming edge
                    // is performed. The distances to the nodes incident to the outgoing edges
                    // are used to decide the necessity of shortcut edges.
                    for ( OtherEdgeID eIn = firstEdge; eIn < lastEdge; eIn++ )
                    {
                        const Edge& edgeIn = _graph->edge(eIn);
                        
                        // need a incoming edge
                        if ( !edgeIn.isDirected(1) ) continue;
    
                        // This is the second stop criterion for the local search. We only want to find
                        // witness paths that are at most as long as the path via the currently processed node.
                        // If we exceed this distance, we can stop the search.
                        // The maximum distance is the length of the incoming edge plus 
                        // the maximum length of an outgoing edge minus the minimum incoming
                        // edge of the node incident to this outgoing edge.
                        // maxOutDist is this length without the length of the incoming edge.
                        int maxOutDist = 0;
                        for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
                        {
                            const Edge& edgeOut = _graph->edge(eOut);
                            if ( !edgeOut.isDirected(0) ) continue;
                                
                            // A path to the start node is not necessary.
                            if ( edgeIn.target() == edgeOut.target() ) continue;
                                
                            // one hop backward search
                            EdgeWeight minBackDist = Weight::MAX_VALUE;
                            OtherEdgeID outFirstEdge = _graph->firstLevelEdge(edgeOut.target());
                            OtherEdgeID outLastEdge = _graph->lastEdge(edgeOut.target());
                            for ( OtherEdgeID e = outFirstEdge; e < outLastEdge; e++ )
                            {
                                const Edge& edge = _graph->edge(e);
                                if ( !edge.isDirected(1) || edge.target() == node ) continue;
                                if ( minBackDist > edge.weight() ) minBackDist = edge.weight();
                            }
                            
                            // Special case: the node incident to the outgoing edge has no ingoing edges 
                            // except for the one from "node". Look out for problems converting unsigned int to int.
                            if (minBackDist < Weight::MAX_VALUE)
                            {
                                if ( maxOutDist < ((int)edgeOut.weight() - (int)minBackDist)) maxOutDist = ((int)edgeOut.weight() - (int)minBackDist);
                            }
                        }
                        
                        // maximum number of settled nodes during local search, 0 = infinite.
                        // This limit can differ between weight calculation (simulateOnly=true) 
                        // and actual contraction.
                        unsigned int maxSettled;
                        if ( phase == PHASE_NODEORDER_INIT || phase == PHASE_NODEORDER_ELIMINATE )
                        {
                            maxSettled = _weightCalc.maxSettledElim;
                            if ( simulateOnly )
                            {
                                maxSettled = _weightCalc.maxSettledApprox;
                            }
                        }
                        else
                        {
                            maxSettled = _contractParams.maxSettledElim;    
                        }
                        
                        // Hop-limit, decrease the hop-limit by 1 since an additional
                        // 1-hop backward search is performed.
                        unsigned int maxHops = _maxHops;
                        if (maxHops > 0) maxHops--;
                        
                        // *** Perform local search. ***
                        // Only perform local search if there is hope to find a witness.
                        // This check is necessary since maxOutDist may be negative.
                        if ( (maxOutDist + (int)edgeIn.weight()) > 0 )
                        {
                            _localDijkstra.searchWithoutTarget(
                                edgeIn.target(),
                                0 /*forward*/,
                                node /*ignore node that will be contracted*/,
                                maxOutDist + edgeIn.weight() /* maximum distance of the search */,
                                noOfTargets,
                                maxSettled /* max settled nodes 0 == inf */,
                                maxHops);
                                
                            if ( phase == PHASE_NODEORDER_INIT || phase == PHASE_NODEORDER_ELIMINATE )
                            {
                                // add up search space sizes of the local searches for weight calculation
                                searchSpace += _localDijkstra.noOfSettledNodes();
                            }
                        }
                        
                        
                        // Local edge reduction uses the results of the local search to remove
                        // edges that are not on any shortest path.
                        if (_localReduceEdges && (phase == PHASE_NODEORDER_ELIMINATE || phase == PHASE_CONSTRUCT) && !simulateOnly)
                        {
                            reduceEdgesLocal( node, edgeIn );
                        }
                        
                        // We finished the local search from start node v = edgeIn.target(). Now we
                        // want to find witness paths to each node w incident to an outgoing edge (node,w)
                        // Either the local search settled w or we perform a 1-hop backward search from w.
                        for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
                        {
                            const Edge& edgeOut = _graph->edge(eOut);
                            if ( !edgeOut.isDirected(0) ) continue;

                            // A path to the start node is not necessary.
                            if ( edgeIn.target() == edgeOut.target() ) continue;

                            bool foundWitness = false;
                            
                            // First check if the local search reached edgeOut.target() with a 
                            // _maxHops - 1 witness path. If there is no such witness, do there
                            // 1-hop backward search.
                            if ( !_localDijkstra.isSettled(0, edgeOut.target()) 
                                || _localDijkstra.distanceTo(edgeOut.target(), 0) > (edgeIn.weight() + edgeOut.weight())
                                || _localDijkstra.parentOf(edgeOut.target(), 0) == node)
                                
                            {
                                // Scan through the edges (x,w) incoming to the node w incident to the outgoing edge
                                // (node,w). (w = edgeOut.target()) to perform a 1-hop backward search. 
                                // Then we check if x is settled by the local search with to find additional
                                // witness paths.
                                OtherEdgeID outFirstEdge = _graph->firstLevelEdge(edgeOut.target());
                                OtherEdgeID outLastEdge = _graph->lastEdge(edgeOut.target());
                                for ( OtherEdgeID e = outFirstEdge; e < outLastEdge; e++ )
                                {
                                    const Edge& edge = _graph->edge(e);
                                    if ( !edge.isDirected(1) || edge.target() == node ) continue;
                                        
                                    if (_localDijkstra.isSettled(0, edge.target()) 
                                        && _localDijkstra.distanceTo(edge.target(), 0) + edge.weight() <= (edgeIn.weight() + edgeOut.weight())
                                        && _localDijkstra.parentOf(edge.target(), 0) != node)
                                    {
                                        foundWitness = true;
                                        break;
                                    }
                                }
                            }
                            else
                            {
                                foundWitness = true;
                            }
    
                            
                            // If shortest path would have increased length, a shortcut edge may be necessary.
                            // Since the search is limted, sometimes shortcuts are added, 
                            // that are not necessary. But they do not invalidate the correctness.
                            if ( !foundWitness )
                            {
                                prepareNewShortcutEdge( node, firstEdge, lastEdge, eIn, edgeIn, eOut, edgeOut );
                            }
                        }    
                        _localDijkstra.clear();
                        
                    }
                    
                    // Remove target-flags that are previously set since these are
                    // stored in the global graph datastructure.
                    for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
                    {
                        const Edge& edgeOut = _graph->edge(eOut);
                        if ( !edgeOut.isDirected(0) ) continue;
                            
                        // one hop backward search
                        OtherEdgeID outFirstEdge = _graph->firstLevelEdge(edgeOut.target());
                        OtherEdgeID outLastEdge = _graph->lastEdge(edgeOut.target());
                        for ( OtherEdgeID e = outFirstEdge; e < outLastEdge; e++ )
                        {
                            const Edge& edge = _graph->edge(e);
                            if ( !edge.isDirected(1) || edge.target() == node ) continue;
                            if (_graph->node(edge.target()).isTarget())
                            {
                                _graph->node(edge.target()).setTarget(false);
                                noOfTargets--;
                            }
                        }
                        
                    }
                    
                    // In case of edge reduction, it is possible that an edge (v,w) has been removed
                    // with (v,node) an incoming edge and (node,w) an outgoing edge. Because of the
                    // 1-hop backward search from node w, the node v becomes a target. So it is 
                    // necessary to check nodes incident to ingoing edges for target flags, too.
                    if (_localReduceEdges && (phase == PHASE_NODEORDER_ELIMINATE || phase == PHASE_CONSTRUCT) && !simulateOnly)
                    {
                        for ( OtherEdgeID eIn = firstEdge; eIn < lastEdge; eIn++ )
                        {
                            const Edge& edgeIn = _graph->edge(eIn);
                            if ( !edgeIn.isDirected(1) ) continue;
                            if (_graph->node(edgeIn.target()).isTarget())
                            {
                                _graph->node(edgeIn.target()).setTarget(false);
                                noOfTargets--;
                            }
                        }
                    }
                    assert( noOfTargets == 0 );
                }
                
                // simple local search with dijkstra
                // ---------------------------------
                // This is the basic local search without any speedups.
                else
                {
					//cout << " perform " << endl;
                    // Prepare target-flags for the local search. The search is stopped
                    // if all targets are settled. These target flags need to be removed
                    // after the local searches.
                    NodeID noOfTargets = 0;
                    for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
                    {
                        const Edge& edgeOut = _graph->edge(eOut);
                        if ( !edgeOut.isDirected(0) ) continue;
                        if (!_graph->node(edgeOut.target()).isTarget())
                        {
                            _graph->node(edgeOut.target()).setTarget(true);
                            noOfTargets++;
                        }
                    }
    
                        
                    // A local search starting at each node incident to an incoming edge
                    // is performed. The distances to the nodes incident to the outgoing edges
                    // are used to decide the necessity of shortcut edges.
                    for ( OtherEdgeID eIn = firstEdge; eIn < lastEdge; eIn++ )
                    {
                        const Edge& edgeIn = _graph->edge(eIn);
                        if ( !edgeIn.isDirected(1) ) continue;
    
                        // This is the second stop criterion for the local search. We only want to find
                        // witness paths that are at most as long as the path via the currently processed node.
                        // If we exceed this distance, we can stop the search.
                        // The maximum distance is the length of the incoming edge plus 
                        // the maximum length of an outgoing edge.
                        // maxOutDist is this length without the length of the incoming edge.
                        EdgeWeight maxOutDist = 0;
                        for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
                        {
                            const Edge& edgeOut = _graph->edge(eOut);
                            if ( !edgeOut.isDirected(0) ) continue;
                                
                            // A path to the start node is not necessary.
                            if ( edgeIn.target() == edgeOut.target() ) continue;
                                
                            if (maxOutDist < edgeOut.weight()) maxOutDist = edgeOut.weight();
                        }
                        
                        // We assume that each edge has positive weight. If maxOutDist == 0,
                        // then there exists no outgoing edge (or only one leading to the
                        // start node), we do not need a local search and also no shortcuts.
                        if (maxOutDist == 0) continue;
                            
                        // maximum number of settled nodes during local search, 0 = infinite.
                        // This limit can differ between weight calculation (simulateOnly=true) 
                        // and actual contraction.
                        unsigned int maxSettled;
                        if ( phase == PHASE_NODEORDER_INIT || phase == PHASE_NODEORDER_ELIMINATE )
                        {
                            maxSettled = _weightCalc.maxSettledElim;
                            if ( simulateOnly )
                            {
                                maxSettled = _weightCalc.maxSettledApprox;
                            }
                        }
                        else
                        {
                            maxSettled = _contractParams.maxSettledElim;    
                        }
                        
                        // *** Perform local search. ***
                        _localDijkstra.searchWithoutTarget(
                            edgeIn.target(),
                            0 /*forward*/,
                            node /*ignore node that will be eliminated*/,
                            maxOutDist + edgeIn.weight() /* maximum distance of the search */,
                            noOfTargets,
                            maxSettled /* max settled nodes 0 == inf */,
                            _maxHops);
                            
                        if ( phase == PHASE_NODEORDER_INIT || phase == PHASE_NODEORDER_ELIMINATE )
                        {
                            searchSpace += _localDijkstra.noOfSettledNodes();
                        }
                        
                        // Local edge reduction uses the results of the local search to remove
                        // edges that are not on any shortest path.
                        if (_localReduceEdges && (phase == PHASE_NODEORDER_ELIMINATE || phase == PHASE_CONSTRUCT) && !simulateOnly)
                        {
                            reduceEdgesLocal( node, edgeIn );  
                        }
                        
                        // The local search starting at v ( = edgeIn.target() ) is finished.
                        // Now we check for each node w incident to an outging edge (node,w) ( = edgeOut )
                        // if the path <v,node,w> is shorter than the shortest path distance v -> w ignoring
                        // node. If so, a shortcut edge is necessary to save shortest paths distances.
                        // Since the local search may be limited, sometimes shortcuts are added, 
                        // that are not necessary. But they do not invalidate the correctness.
                        for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
                        {
                            const Edge& edgeOut = _graph->edge(eOut);
                            if ( !edgeOut.isDirected(0) ) continue;

                            // no self loops
                            if ( edgeIn.target() == edgeOut.target() ) continue;
    
                            // Compare length of path <v,node,w> to distance between v -> w ignoring node.
                            if ( !_localDijkstra.isSettled(0, edgeOut.target()) 
                                || _localDijkstra.distanceTo(edgeOut.target(), 0) > (edgeIn.weight() + edgeOut.weight())
                                || _localDijkstra.parentOf(edgeOut.target(), 0) == node)
                            {
                                
                                prepareNewShortcutEdge( node, firstEdge, lastEdge, eIn, edgeIn, eOut, edgeOut );
                            }
                            
                            // If a witness exists, we may want to store the witness to use This
                            // information in a later hieararchy construction.
                            else
                            {
                                if ( phase == PHASE_NODEORDER_ELIMINATE && !simulateOnly )
                                {
                                    storeWitness( node, edgeIn, edgeOut );
                                }
                            }
                        }    
                        _localDijkstra.clear();
                        
                    }
                    
                    // Remove previously added target flags from the global graph datastructure.
                    // If we do not, they manipulate subsequent local searches.
                    for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
                    {
                        const Edge& edgeOut = _graph->edge(eOut);
                        if ( !edgeOut.isDirected(0) ) continue;
                        _graph->node(edgeOut.target()).setTarget(false);
                    }
                }
                
            }
            
            VERBOSE_CONTRACT( cout << " search space " << searchSpace << flush; )
            VERBOSE_CONTRACT( cout << " old " << 2*(lastEdge-firstEdge) << flush; )
            VERBOSE_CONTRACT( cout << " new " << 2*_newEdges.size() << flush; )

            // Save shortcuts to text file, starting witn node, then the triple (v,w,weight) for each shortcut (v,w).
            if ( saveShortcutsText && phase == PHASE_NODEORDER_ELIMINATE && !simulateOnly && _shortcutsText.is_open() )
            {
                _shortcutsText << node;
            }

            // Now we add necessary shortcuts that are stored in a buffer (_newEdges).
            // Also some statistics and outputs are performed.
            unsigned int shortcutOriginalEdgeSumTemp = 0;
            for ( typename vector<NewEdge>::iterator iter = _newEdges.begin(); iter != _newEdges.end(); iter++ )
            {
                // *** Add shortcuts ***
                // The addShorcutEdge...() function returns the edge difference in terms of
                // used entries in the edge array. It returns -2, 0 or 2.
                int diff;
                if ( simulateOnly )
                {
                    diff = _graph->addShortcutEdgeSimulate(iter->source, iter->edge);
                }
                else
                {
                    diff = _graph->addShortcutEdge(iter->source, iter->edge);
                    
                    // Output of shortcut edges that can be used in a later hierarchy construction.
                    if ( phase == PHASE_NODEORDER_ELIMINATE && saveShortcutsWitnesses)
                    {
                        _storeShortcuts.addShortcut(node, iter->source, iter->edge.target(), iter->edge.isBidirected());
                    }
                }
                
                edgeDiff += diff;

                if ( phase == PHASE_NODEORDER_INIT || phase == PHASE_NODEORDER_ELIMINATE )
                {
                    // Used for weight calculation.
                    newEdgesCounter += diff;
                    
                    // Priority term: Count sum of original edges the new shortcuts represent.
                    // Count bidirectional edges twice.
                    shortcutOriginalEdgeSumTemp += iter->edge.shortcutOriginalEdgeCount();
                    if (iter->edge.isBidirected()) shortcutOriginalEdgeSumTemp += iter->edge.shortcutOriginalEdgeCount();
                }

                // output shorcuts, one line per contracted node
                if ( saveShortcutsText && phase == PHASE_NODEORDER_ELIMINATE && !simulateOnly && _shortcutsText.is_open() )
                {
                    _shortcutsText << " " << iter->source << " " << iter->edge.target() << " " << iter->edge.weight();
                    if (iter->edge.isBidirected()) _shortcutsText << " " << iter->edge.target() << " " << iter->source << " " << iter->edge.weight();
                }
                
            }
            
            // Clear buffer of new edges.
            _newEdges.clear();
            

            // Update new level of contracted node. In case of node ordering Also
            // update attributes and neighbors.
			//ofstream congraph("levelgraph.txt", ios::out);
            if ( !simulateOnly && (phase == PHASE_NODEORDER_ELIMINATE || phase == PHASE_CONSTRUCT) )
            {
				//cout << "currentLevel: " << _currentLevel << endl;
				
				

				//for ( NodeID u = 0; u < _graph->noOfNodes(); u++ )
				//{
				//	//congraph << graph->edge(u) << endl;
				//	OtherEdgeID lastEdge = _graph->lastEdge(_graph->mapExtToIntNodeID(u));
				//	cout << _graph->firstEdge(u) << " " << lastEdge << endl;
				//	for ( OtherEdgeID e = _graph->firstLevelEdge(_graph->mapExtToIntNodeID(u)); e < lastEdge; e++ )
				//	{
				//		const Edge& edge = _graph->edge(e);
				//		int dir = 3;
				//		if (edge.isDirected(0)) dir -= 2;
				//		if (edge.isDirected(1)) dir -= 1;
				//		congraph << u << " " << edge.target() << " " << edge.weight() << " " << dir << endl;
				//	}
				//}

                // before updating neighbors set current node into right level
                _graph->node(node).setLevel(_currentLevel);
                _graph->changeNodeLevelOnlyReverseEdges(node);	
				                
                _noOfEdges += edgeDiff;
                _noOfNodes--;
                
                if ( phase == PHASE_NODEORDER_ELIMINATE )
                {
                    pqData(node).level = _currentLevel;
                    VERBOSE_CONTRACT( cout << " edges " << _noOfEdges << flush; )
                    updateAfterContraction( node, firstEdge, lastEdge );
                    VERBOSE_CONTRACT( cout << " )" << flush; )
                }

				//for ( NodeID u = 0; u < _graph->noOfNodes(); u++ )
				//{
				//	//congraph << graph->edge(u) << endl;
				//	OtherEdgeID lastEdge = _graph->lastEdge(_graph->mapExtToIntNodeID(u));
				//	cout << _graph->firstEdge(u) << " " << lastEdge << endl;
				//	for ( OtherEdgeID e = _graph->firstLevelEdge(_graph->mapExtToIntNodeID(u)); e < lastEdge; e++ )
				//	{
				//		const Edge& edge = _graph->edge(e);
				//		int dir = 3;
				//		if (edge.isDirected(0)) dir -= 2;
				//		if (edge.isDirected(1)) dir -= 1;
				//		cout << u << " " << edge.target() << " " << edge.weight() << " " << dir << endl;
				//	}
				//}
            }

            if ( phase == PHASE_NODEORDER_INIT || phase == PHASE_NODEORDER_ELIMINATE )
            {
                if ( saveShortcutsText && phase == PHASE_NODEORDER_ELIMINATE && !simulateOnly && _shortcutsText.is_open() )
                {
                    _shortcutsText << endl;
                }
                // Return attributes/priority terms, if requested.
                if ( shortcutOriginalEdgeSum != NULL ) *shortcutOriginalEdgeSum = shortcutOriginalEdgeSumTemp;
                if ( searchSpaceResult != NULL ) *searchSpaceResult = searchSpace;
                if ( edgeDiffResult    != NULL ) *edgeDiffResult    = edgeDiff;
                if ( newEdgesResult    != NULL ) *newEdgesResult    = newEdgesCounter;
            }
            
        }
        
        /** 
         * Distribute the Voronoi region R(node) of the current node among the neighboring voronoi regions.
         * This is done using a modified Dijkstra algorithm:
         *   Init: For each node x in R(node), that has an incoming edge (y,x) with
         *   y in R(z), z != node, add node x to the priority queue of Dijkstra algorithm
         *   with distance d(z,y)+w(y,x) and possibly new owner z.
         *   Then we always extract the node with the shortest distance and add it to the
         *   Voronoi region of its new owner. If during the edge relaxtion an updateAfterContraction
         *   is successful, also the possibly new owner is updated.
         */
        void distributeVoronoiRegion(const NodeID node)
        {
            VERBOSE_CONTRACT_VORONOI( cout << "distributeVornoiRegion " << node << " voronoiNumber " << pqData(node).voronoiNumber << endl; )
            // The nodes of the Voronoi region are stored in a single linked list that
            // is terminated by SPECIAL_NODEID. The DijkstraCH class provides
            // only low level methods to update/insert nodes and extract The
            // node with the lowest distance. The main logic is in this subroutine.
            NodeID current = node;
            while ( current != SPECIAL_NODEID )
            {
                VERBOSE_CONTRACT_VORONOI( cout << "current " << current << endl; )
                
                // Under all nodes that are incident to an incoming node of the currently
                // regarded node in R(node), and that are not in R(node), take the node
                // that is on the shortest path from its owner to the currently regarded node.
                EdgeWeight minDist = Weight::MAX_VALUE;
                NodeID minDistNode = SPECIAL_NODEID;
                OtherEdgeID lastEdge = _graph->lastEdge(current);
                for ( OtherEdgeID e = _graph->firstEdge(current); e < lastEdge; e++ )
                {
                    const Edge& edge = _graph->edge(e);
                    if ( !edge.isDirected(1) ) continue;
                    const PQueueNodeElimination& data = pqData(edge.target());
                    // Ignore nodes in Voronoi region R(node).
                    if ( data.voronoiOwner == node ) continue;
                    // Special case: because the graph is directed, not all nodes can be distributed.
                    // Such nodes that cannot be distributed are identified by the fact that their
                    // Voronoi region owner is already eliminated (contracted). These nodes
                    // must not be regarded for the Vornoi region distribution.
                    if ( pqData(data.voronoiOwner).isEliminated() ) continue;
                    EdgeWeight dist = data.voronoiNumber + edge.weight();
                    VERBOSE_CONTRACT_VORONOI( cout << "from " << edge.target() << " dist " << dist << endl; )
                    if ( dist < minDist )
                    {
                        minDist = dist;
                        minDistNode = edge.target();
                    }
                }
                // Special case: It is possibly, that no neighboring Voronoi region exists that can reach 
                // the currently regarded node. In this case, the node is ignored.
                if ( minDistNode != SPECIAL_NODEID )
                {
                    VERBOSE_CONTRACT_VORONOI( cout << "insert " << current << " dist " << minDist << " parent " << minDistNode << endl; )
                    _dVoronoi.insertNode( current, minDist, minDistNode );
                }
                current = pqData(current).voronoiNextBorderNode;
            }
            
            NodeID v;
            EdgeWeight dist;
            NodeID parent;
            // *** Process nodes in the priority queue. Main work loop ***
            // searchNext() returns the node with the lowest distance in the priority queue
            // by call-by-reference, along with the distance and the parent. the
            // return value of the function indicates whetere the priority queue is not empty.
            while ( _dVoronoi.searchNext(v,dist,parent) )
            {
                // Node v is assigned to the Voronoi region of its parent.
                PQueueNodeElimination& data = pqData(v);
                data.voronoiOwner = pqData(parent).voronoiOwner;
                data.voronoiNumber = dist;
                PQueueNodeElimination& dataOwner = pqData(data.voronoiOwner);
                
                VERBOSE_CONTRACT_VORONOI( cout << "node " << v << " dist " << dist << " parent " << parent << " new owner " << data.voronoiOwner << endl; )

                assert( !dataOwner.isEliminated() );
               
                // Update stats of voronoi owner
                dataOwner.voronoiNumber++;
                
                // Add v to linked list of voronoi Region of new owner.
                data.voronoiNextBorderNode = dataOwner.voronoiNextBorderNode;
                dataOwner.voronoiNextBorderNode = v;

                // Relax edges, but only into voronoi region of "node"
                OtherEdgeID lastEdge = _graph->lastEdge(v);
                for ( OtherEdgeID e = _graph->firstEdge(v); e < lastEdge; e++ )
                {
                    const Edge& edge = _graph->edge(e);
                    if ( !edge.isDirected(0) || pqData(edge.target()).voronoiOwner != node ) continue;
                    EdgeWeight newDist = dist + edge.weight();
                    _dVoronoi.updateNode( edge.target(), newDist, v );
                }
            }
            _dVoronoi.clear();            
        }

        // ***
        // Methods for storing witness paths and shortcut edges.
        // ***

    	/** 
    	 * Adds listes shortcut edges, works even if original graph has changed edge weights
    	 * or bidir edges split to unidir edges. But does not work if the edges a shortcut
    	 * represent are not there.
    	 */
    	void addShortcuts(istream& ins)
    	{
    	    RetrieveShortcuts r(ins);
    	    NodeID node, in, out;
    	    bool bidir;
    	    while (r.nextNode(node))
    	    {
    	        OtherEdgeID firstEdge =_graph->firstLevelEdge(node);
    	        OtherEdgeID lastEdge = _graph->lastEdge(node);
    	        while (r.nextIn(in))
    	        {
    	            bool bIn = false;  // in edge found?
    	            for ( OtherEdgeID eIn = firstEdge; eIn < lastEdge; eIn++ )
    	            {
    	                const Edge& edgeIn = _graph->edge(eIn);
    	                if (edgeIn.isDirected(1) && edgeIn.target() == in)
    	                {
            	            while (r.nextOut(out, bidir))
            	            {
            	                bool bOut = false;  // out edge found?
            	                for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
            	                {
            	                    const Edge& edgeOut = _graph->edge(eOut);
            	                    if (edgeOut.isDirected(0) && edgeOut.target() == out)
            	                    {
            	                        if (!bidir || (edgeIn.isBidirected() && edgeOut.isBidirected()) )
            	                        {
            	                            _noOfEdges += _graph->addShortcutEdge(in, Edge(out, edgeIn.weight()+edgeOut.weight(), EDGE_TYPE_SHORTCUT, true, bidir, node));
            	                        }
            	                        
            	                        // shortcut was bidirected but inEdge and outEdge are not both bidirected
            	                        // try to find edges for shortcut in opposite direction
            	                        else
            	                        {
            	                            EdgeWeight weight = edgeIn.weight()+edgeOut.weight();
            	                            EdgeWeight weight2 = 0;
                                            bool bIn2 = false;
                                            bool bOut2 = false;
            	                            if (edgeIn.isBidirected())
            	                            {
            	                                bIn2 = true;
            	                                weight2 += edgeIn.weight();    
            	                            }
            	                            else
            	                            {
            	                                for ( OtherEdgeID eIn2 = firstEdge; eIn2 < lastEdge; eIn2++ )
    	                                        {
    	                                            const Edge& edgeIn2 = _graph->edge(eIn2);
    	                                            if (edgeIn2.isDirected(0) && edgeIn2.target() != in)
    	                                            {
    	                                                bIn2 = true;
    	                                                weight2 += edgeIn2.weight();
    	                                                break;
    	                                            }
    	                                        }
    	                                        if (!bIn2)
    	                                        {
    	                                            cerr << "Did not find an in edge: " << in << " <- " << node << endl;
    	                                        }
    	                                    }
    	                                    
    	                                    if (edgeOut.isBidirected())
    	                                    {
    	                                        bOut2 = true;
    	                                        weight2 += edgeOut.weight();
    	                                    }
            	                            else
            	                            {
            	                                for ( OtherEdgeID eOut2 = firstEdge; eOut2 < lastEdge; eOut2++ )
    	                                        {
    	                                            const Edge& edgeOut2 = _graph->edge(eOut2);
    	                                            if (edgeOut2.isDirected(1) && edgeOut2.target() != out)
    	                                            {
    	                                                bOut2 = true;
    	                                                weight2 += edgeOut2.weight();
    	                                                break;
    	                                            }
    	                                        }
    	                                        if (!bOut2)
    	                                        {
    	                                            cerr << "Did not find an out edge: " << node << " <- " << out << endl;
    	                                        }
    	                                    }
    
                                            // if weight for both directions is the same, add bidir edge
                                            if (bIn2 && bOut2 && weight == weight2)
                                            {
        	                                    _noOfEdges += _graph->addShortcutEdge(in, Edge(out, weight, EDGE_TYPE_SHORTCUT, true, true, node));
        	                                }
        	                                // else add two edges with different weight
        	                                else
        	                                {
        	                                    _noOfEdges += _graph->addShortcutEdge(in, Edge(out, weight, EDGE_TYPE_SHORTCUT, true, false, node));
        	                                    
        	                                    // only add second shortcut if all edges have been found
        	                                    if (bIn2 && bOut2)
        	                                    {
        	                                        _noOfEdges += _graph->addShortcutEdge(out, Edge(in, weight2, EDGE_TYPE_SHORTCUT, true, false, node));
        	                                    }
        	                                }
            	                        }   
                	                    bOut = true;
                	                    break;
            	                    }
            	                }
                	            if ( !bOut )
                	            {
                	                cout << "Did not find an out edge: " << node << " -> " << out << endl;
                	            }
            	            }
            	            bIn = true;
            	            break;
            	        }
    	            }
    	            if ( !bIn )
    	            {
    	                cout << "Did not find an in edge: " << in << " -> " << node << endl;
    	            }
                }
            }   
    	}
    	
    	/**
    	 * Calculates weight of witnesses and adds witness shortcuts
    	 * works even if original graph has changed edge weights.
    	 * Does not work if edges on the witness path have been deleted.
    	 */
    	void addWitnessShortcuts(RetrieveWitnesses& r)
    	{
    	    NodeID node, in, out;
    	    while (r.nextWitness(node, in, out))
    	    {
    	        NodeID current = out;
    	        NodeID next;
    	        EdgeWeight weight = 0;
  	            bool bNext = true;
    	        while ( r.nextStep(next) )
    	        {
    	            if ( bNext )
    	            {
        	            bNext = false;
        	            OtherEdgeID lastEdge = _graph->lastEdge(current);
        	            for ( OtherEdgeID e = _graph->firstLevelEdge(current); e < lastEdge; e++ )
        	            {
        	                const Edge& edge = _graph->edge(e);
        	                if ( edge.isDirected(1) && edge.target() == next )
        	                {
        	                    bNext = true;
        	                    weight += edge.weight();
        	                    break;
        	                }
        	            }
        	            current = next;
        	        }
    	        }
    	        
    	        if ( bNext )
    	        {
        	        // find edge to in
    	            OtherEdgeID lastEdge = _graph->lastEdge(current);
    	            bNext = false;
    	            for ( OtherEdgeID e = _graph->firstLevelEdge(current); e < lastEdge; e++ )
    	            {
    	                const Edge& edge = _graph->edge(e);
    	                if ( edge.isDirected(1) && edge.target() == in )
    	                {
    	                    bNext = true;
    	                    weight += edge.weight();
    	                    break;
    	                }
    	            }
                }
                
                // only add shortcut edge if all nodes of the path have been found
                if ( bNext )
                {
    	            // add a shortcut edge for this witness marked as shortcut 
    	            // attention: edge.isShortcut() means edge is a witness-shortcut
    	            _noOfEdges += _graph->addShortcutEdge(in, Edge(out, weight, EDGE_TYPE_WITNESS_SHORTCUT, true, false));
    	            //_graph->addEdge(in, Edge(out, weight, true, true, false));
    	            //_graph->addEdge(out, Edge(in, weight, true, false, true));
    	        }
    	        
    	        if ( !bNext )
    	        {
    	            cout << "missing edges for witness" << endl;
    	        }
    	        
    	    }
    	}

    	void validateWitnesses(const NodeID node)
    	{
    	    OtherEdgeID firstEdge = _graph->firstLevelEdge(node);
    	    OtherEdgeID lastEdge = _graph->lastEdge(node);
    	    bool bInitBucketsMTM = false;
    	    
    	    
    	    // sort (OPEN QUESTION: does this always help)
    	    if (firstEdge < lastEdge)
    	    {
    	        vector<NewEdge> newEdges;
    	        _graph->sortEdges(firstEdge, lastEdge);

                for ( OtherEdgeID eIn = firstEdge; eIn < lastEdge; eIn++ )
                {
                    const Edge& edgeIn = _graph->edge(eIn);
                    if ( !edgeIn.isDirected(1) || edgeIn.isShortcut() ) continue;
                        
                    OtherEdgeID inFirstEdge = _graph->firstLevelEdge(edgeIn.target());
                    OtherEdgeID inLastEdge = _graph->lastEdge(edgeIn.target());
                    
                    bool bFindPossibleWitnessesMTM = false;
                        
                    for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
                    {
                        const Edge& edgeOut = _graph->edge(eOut);
                        if ( edgeIn.target() == edgeOut.target() || !edgeOut.isDirected(0) || edgeOut.isShortcut() ) continue;
                        
                        bool bValid = false;
                        for ( OtherEdgeID e = inFirstEdge; e < inLastEdge; e++ )
                        {
                            const Edge& edge = _graph->edge(e);
                            if ( edge.isDirected(0) && edge.target() == edgeOut.target() )
                            {
                                if ( edge.weight() <= (edgeIn.weight() + edgeOut.weight()) )
                                {
                                    bValid = true;
                                    break;
                                }
                                //break;
                            }
                        }
                        
                        // if there is no valid shortcut, create one
                        if ( !bValid )
                        {
                            // make a simplified 2 hop search for witnesses
                            if ( !bInitBucketsMTM )
                            {
                                initBucketsMTM(node, firstEdge, lastEdge);
                                bInitBucketsMTM = true;
                            }
                            if ( !bFindPossibleWitnessesMTM )
                            {
                                findPossibleWitnessesMTM(node, edgeIn.target());                            
                                bFindPossibleWitnessesMTM = true;
                            }
                            
                            // if no two-hop witnesses exists, add shortcut
                            if (_possibleWitnesses[edgeOut.target()] > (edgeIn.weight() + edgeOut.weight()))
                            {
                                cout << "witness " << edgeIn.target() << " -> " << edgeOut.target() << " no longer valid, create shortcut" << endl;
                                prepareNewShortcutEdge( node, firstEdge, lastEdge, eIn, edgeIn, eOut, edgeOut );
                            }
                        }
                    }
                    
                    // clear _possibleWitnesses array
                    if ( bFindPossibleWitnessesMTM )
                    {
                        for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
                        {
                            const Edge& edgeOut = _graph->edge(eOut);
                            if ( !edgeOut.isDirected(0) ) continue;
                                
                            // no self loops
                            if ( edgeIn.target() == edgeOut.target() ) continue;
        
                            _possibleWitnesses[edgeOut.target()] = Weight::MAX_VALUE;
                        }
                    }
                    
                }
                
                if ( bInitBucketsMTM )
                {
                    // cleanup
                    clearMTM(node, firstEdge, lastEdge);
                }

                _noOfEdges -= 2*(lastEdge - firstEdge);
                for ( typename vector<NewEdge>::const_iterator iter = _newEdges.begin(); iter != _newEdges.end(); iter++ )
                {
                    _noOfEdges += _graph->addShortcutEdge(iter->source, iter->edge);
                }
                _newEdges.clear();
        	}
        	_noOfNodes--;
        }

    	void checkWitnesses(const NodeID node, RetrieveWitnesses& r)
    	{
    	    OtherEdgeID firstEdge = _graph->firstLevelEdge(node);
    	    OtherEdgeID lastEdge = _graph->lastEdge(node);
    	    bool bInitBucketsMTM = false;
    	    
    	    if (firstEdge < lastEdge)
    	    {
    	        // sort (OPEN QUESTION: does this always help)
    	        _graph->sortEdges(firstEdge, lastEdge);

                for ( OtherEdgeID eIn = firstEdge; eIn < lastEdge; eIn++ )
                {
                    const Edge& edgeIn = _graph->edge(eIn);
                    if ( !edgeIn.isDirected(1) || edgeIn.isShortcut() ) continue;
                        
                    OtherEdgeID inFirstEdge = _graph->firstLevelEdge(edgeIn.target());
                    OtherEdgeID inLastEdge = _graph->lastEdge(edgeIn.target());
                    
                    bool bFindPossibleWitnessesMTM = false;
                        
                    for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
                    {
                        const Edge& edgeOut = _graph->edge(eOut);
                        if ( edgeIn.target() == edgeOut.target() || !edgeOut.isDirected(0) || edgeOut.isShortcut() ) continue;

                        bool shortcutNecessary = true;

                        // check for stored witness
                        NodeID wNode, wIn, wOut;
                        if ( r.nextWitness(wNode, wIn, wOut, true) )
                        {
                            // witness available, validate witness
                            if ( wNode == node && wIn == edgeIn.target() && wOut == edgeOut.target() )
                            {
                                r.read(); // no peek
                                
                                NodeID wCurrent = wOut;
                                NodeID wNext;
                                EdgeWeight wWeight = 0;
                                bool bNext = true;
                                while ( r.nextStep(wNext) )
                                {
                                    if ( bNext )
                                    {
                                        OtherEdgeID wLastEdge = _graph->lastEdge(wCurrent);
                                        bNext = false;
                                        for ( OtherEdgeID e = _graph->firstLevelEdge(wCurrent); e < wLastEdge; e++ )
                                        {
                                            const Edge& edge = _graph->edge(e);
                                            if ( edge.isDirected(1) && edge.target() == wNext )
                                            {
                                                bNext = true;
                                                wWeight += edge.weight();
                                                break;
                                            }
                                        }
                                        wCurrent = wNext;
                                    }
                                }

                                if ( bNext )
                                {                                
                                    // find edge to in
                                    OtherEdgeID wLastEdge = _graph->lastEdge(wCurrent);
                                    bNext = false;
                                    for ( OtherEdgeID e = _graph->firstLevelEdge(wCurrent); e < wLastEdge; e++ )
                                    {
                                        const Edge& edge = _graph->edge(e);
                                        if ( edge.isDirected(1) && edge.target() == wIn )
                                        {
                                            bNext = true;
                                            wWeight += edge.weight();
                                            break;
                                        }
                                    }
                                }
                                
                                // if witness is short enough, no shortcut necessary
                                if ( bNext && (wWeight <= (edgeIn.weight() + edgeOut.weight())) )
                                {
                                    shortcutNecessary = false;
                                }
                            }
                        }
                        
                        // if a shortcut is necessary, do a one hop search
                        if ( shortcutNecessary )
                        {
                            for ( OtherEdgeID e = inFirstEdge; e < inLastEdge; e++ )
                            {
                                const Edge& edge = _graph->edge(e);
                                if ( edge.isDirected(0) && edge.target() == edgeOut.target() )
                                {
                                    if ( edge.weight() <= (edgeIn.weight() + edgeOut.weight()) )
                                    {
                                        shortcutNecessary = false;
                                        break;
                                    }
                                    //break;
                                }
                            }
                        }
                        
                        // if there is no valid shortcut, create one
                        if ( shortcutNecessary )
                        {
                            // make a simplified 2 hop search for witnesses
                            if ( !bInitBucketsMTM )
                            {
                                initBucketsMTM(node, firstEdge, lastEdge);
                                bInitBucketsMTM = true;
                            }
                            if ( !bFindPossibleWitnessesMTM )
                            {
                                findPossibleWitnessesMTM(node, edgeIn.target());                            
                                bFindPossibleWitnessesMTM = true;
                            }
                            
                            // if no two-hop witnesses exists, add shortcut
                            if (_possibleWitnesses[edgeOut.target()] > (edgeIn.weight() + edgeOut.weight()))
                            {
                                prepareNewShortcutEdge( node, firstEdge, lastEdge, eIn, edgeIn, eOut, edgeOut );
                            }
                        }
                    }
                    
                    // clear _possibleWitnesses array
                    if ( bFindPossibleWitnessesMTM )
                    {
                        for ( OtherEdgeID eOut = firstEdge; eOut < lastEdge; eOut++ )
                        {
                            const Edge& edgeOut = _graph->edge(eOut);
                            if ( !edgeOut.isDirected(0) ) continue;
                                
                            // no self loops
                            if ( edgeIn.target() == edgeOut.target() ) continue;
        
                            _possibleWitnesses[edgeOut.target()] = Weight::MAX_VALUE;
                        }
                    }
                    
                }
                
                if ( bInitBucketsMTM )
                {
                    // cleanup
                    clearMTM(node, firstEdge, lastEdge);
                }
                
                _noOfEdges -= 2*(lastEdge - firstEdge);
                for ( typename vector<NewEdge>::const_iterator iter = _newEdges.begin(); iter != _newEdges.end(); iter++ )
                {
                    _noOfEdges += _graph->addShortcutEdge(iter->source, iter->edge);
                }
                _newEdges.clear();
        	}
        	_noOfNodes--;
        }

        /** 
         * Debugging routines: check whetere shortest paths are preserved during contraction.
         * These are time consuming and should only be used for debugging and on small graphs.
         * And they also contain some custom code like special node ids where errors occured.
         */
        void initTestShortestPaths()
        {
            srand(25);
            for (NodeID x = 0; x < _noOfTestCases; x++) {
                NodeID s = randomNodeID(_graph->noOfNodes());
                NodeID t = randomNodeID(_graph->noOfNodes());
                _testRuns.push_back( stPair(s, t) );
            }
            _testRuns.push_back( stPair(24066, 24537) );

            VERBOSE_CONTRACT( cout << "initTestShortestPaths" << endl; )
            _testPaths.resize(_testRuns.size(), Weight::MAX_VALUE);
            
            for (NodeID x = 0; x < _testRuns.size(); x++) {
                _testDijkstra.bidirSearch(_testRuns[x].first, _testRuns[x].second);
                _testDijkstra.pathTo(_testPaths[x],_testRuns[x].second,-1);
                VERBOSE_CONTRACT( cout << _testPaths[x] << endl; )
                _testDijkstra.clear();
            }
        }

        /** 
         * Debugging routines: check whetere shortest paths are preserved during contraction.
         * These are time consuming and should only be used for debugging and on small graphs.
         * And they also contain some custom code like special node ids where errors occured.
         */
        void testShortestPaths()
        {
            VERBOSE_CONTRACT( cout << "testShortestPaths" << endl; )
            for (NodeID x = 0; x < _testRuns.size(); x++) {
                _testDijkstra.bidirSearch(_testRuns[x].first, _testRuns[x].second);
                Path postPath;
                _testDijkstra.pathTo(postPath,_testRuns[x].second,-1);
                _testDijkstra.clear();
                Path& prePath = _testPaths[x];
                
                if (prePath.length() != postPath.length())
                {
                    cerr << "source " << _testRuns[x].first << " level " << _graph->node(_testRuns[x].first).level() << " target " << _testRuns[x].second << " level " << _graph->node(_testRuns[x].second).level() << " dist pre " << prePath.length() << " post " << postPath.length() << endl;
                    cerr << "PRE " << prePath << endl;
                    cerr << "POST " << postPath << endl;
                    exit(1);
                }
            }
        }
	};
} // namespace

#endif // _PROCESSING_CONSTRUCTCH_H
