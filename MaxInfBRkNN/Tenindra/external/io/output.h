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

#ifndef OUTPUT_H
#define OUTPUT_H

#include "../stats/utils.h"

/**
 * Export search graph to text file in ddsg-format.
 */
template< typename Graph >
void exportSearchGraph(ostream& out, const Graph* g)
{
    out << "d" << endl;
    OtherEdgeID m = 0;
    vector<NodeID> mapIntToExtNodeID(g->noOfNodes());
    for ( NodeID u = 0; u < g->noOfNodes(); u++ )
    {
        m += g->lastEdge(g->mapExtToIntNodeID(u)) - g->firstLevelEdge(g->mapExtToIntNodeID(u));
        mapIntToExtNodeID[g->mapExtToIntNodeID(u)] = u;
    }
    
    out << g->noOfNodes() << " " << m << endl;
    for ( NodeID u = 0; u < g->noOfNodes(); u++ )
    {
        OtherEdgeID lastEdge = g->lastEdge(g->mapExtToIntNodeID(u));
        for ( OtherEdgeID e = g->firstLevelEdge(g->mapExtToIntNodeID(u)); e < lastEdge; e++ )
        {
            const Edge& edge = g->edge(e);
            int dir = 3;
            if (edge.isDirected(0)) dir -= 2;
            if (edge.isDirected(1)) dir -= 1;
            out << u << " " << mapIntToExtNodeID[edge.target()] << " " << edge.weight() << " " << dir << endl;
        }
    }
}

/**
 * Export a contraction hierarchy in a binary format, see ../docu/chDocu.html.
 * supported graph classes: datastr::graph::UpdateableGraph
 */
template< typename Graph>
void exportContractionHierarchy(ostream& out, const Graph* g)
{

    bool expand = false;
    CH_EXPAND( expand = true );
    if (!expand)
    {
        cerr << "No shortcut information stored into edges, cannot proceed with export!" << endl;
        cerr << "Please #define USE_CH_EXPAND in config.h" << endl;
        exit(1);
    }
    
    double time = timestamp();

    unsigned int version = 1;
    unsigned int n = g->noOfNodes();
    unsigned int m1 = 0;
    unsigned int m2 = 0;

    // edge flags    
    const unsigned int CH_EDGE_FLAG_FORWARD  = 1 << 0;
    const unsigned int CH_EDGE_FLAG_BACKWARD = 1 << 1;
    const unsigned int CH_EDGE_FLAG_SHORTCUT = 1 << 2;
    
    for ( NodeID u = 0; u < g->noOfNodes(); u++ )
    {
        OtherEdgeID lastEdge = g->lastEdge(u);
        for ( OtherEdgeID e = g->firstLevelEdge(u); e < lastEdge; e++ )
        {
            const Edge& edge = g->edge(e);
            if ( edge.isShortcut() )
            {
                m2++;
            }
            else
            {
                m1++;
            }
        }
    }
    
    VERBOSE( cout << "Contraction Hierarchy file format version " << version << endl; )
    VERBOSE( cout << n << " nodes" << endl; )
    VERBOSE( cout << m1 << " original edges" << endl; )
    VERBOSE( cout << m2 << " shortcut edges" << endl; )
    
    // === BEGIN OF EXPORT ===
    
    // header
    // ------
    // id
    out << "CH\r\n";
    // version
    writePrimitive( out, version );
    // n (number of nodes)
    writePrimitive( out, n );
    // m1 ( number of original (=non shortcut) edges )
    writePrimitive( out, m1 );
    // m2 ( number of shortcut edges )
    writePrimitive( out, m2 );
    
    // levels for nodes 0..(n-1)
    // -------------------------

    VERBOSE( cout << "Export levels..." << endl; );
    VERBOSE( Percent percent1( n ); )
    VERBOSE( unsigned int i1 = 0; )
    unsigned int level;
    for ( NodeID u = 0; u < n; u++ )
    {
        level = g->node(u).level();
        writePrimitive( out, level );
        VERBOSE( percent1.printStatus( i1++ ); )
    }
    
    // orignal (=non shortcut) edges
    // -----------------------------
    VERBOSE( cout << "Export original edges..." << endl; );
    VERBOSE( Percent percent2( m1 ); )
    VERBOSE( unsigned int i2 = 0; )
    for ( NodeID u = 0; u < g->noOfNodes(); u++ )
    {
        OtherEdgeID lastEdge = g->lastEdge(u);
        for ( OtherEdgeID e = g->firstLevelEdge(u); e < lastEdge; e++ )
        {
            const Edge& edge = g->edge(e);
            if ( !edge.isShortcut() )
            {
                unsigned int source = u;
                unsigned int target = edge.target();
                unsigned int weight = edge.weight();
                
                // edge flags
                unsigned int flags = 0;
                if ( edge.isDirected(0) ) flags |= CH_EDGE_FLAG_FORWARD;
                if ( edge.isDirected(1) ) flags |= CH_EDGE_FLAG_BACKWARD;
                
                // source node
                writePrimitive( out, source );
                // target node
                writePrimitive( out, target );
                // weight
                writePrimitive( out, weight );
                // flags
                writePrimitive( out, flags );
                
                VERBOSE( percent2.printStatus( i2++ ); )       
            }
        }
    }
    
    // shortcut edges
    // --------------
    VERBOSE( cout << "Export shortcut edges..." << endl; );
    VERBOSE( Percent percent3( m2 ); )
    VERBOSE( unsigned int i3 = 0; )
    for ( NodeID u = 0; u < g->noOfNodes(); u++ )
    {
        OtherEdgeID lastEdge = g->lastEdge(u);
        for ( OtherEdgeID e = g->firstLevelEdge(u); e < lastEdge; e++ )
        {
            const Edge& edge = g->edge(e);
            if ( edge.isShortcut() )
            {
                unsigned int source = u;
                unsigned int target = edge.target();
                unsigned int weight = edge.weight();
                unsigned int middle = edge.shortcutMiddle();
                
                // edge flags
                unsigned int flags = CH_EDGE_FLAG_SHORTCUT;
                if ( edge.isDirected(0) ) flags |= CH_EDGE_FLAG_FORWARD;
                if ( edge.isDirected(1) ) flags |= CH_EDGE_FLAG_BACKWARD;
                
                // source node
                writePrimitive( out, source );
                // target node
                writePrimitive( out, target );
                // weight
                writePrimitive( out, weight );
                // flags
                writePrimitive( out, flags );
                // middle node of shortcut
                writePrimitive( out, middle );
                
                VERBOSE( percent3.printStatus( i3++ ); )       
            }
        }
    }
    unsigned int end = 0x12345678;
    writePrimitive( out, end );    
    // === END OF EXPORT
    
    time = timestamp() - time;
    VERBOSE( cout << "Export done, took " << time << " seconds." << endl; )
    
}

#endif // OUTPUT_H
