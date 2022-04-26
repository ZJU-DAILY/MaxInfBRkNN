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

#ifndef CREATEGRAPH_H
#define CREATEGRAPH_H


#include <cmath>
//#include <cstdio>

#include "../datastr/graph/graph.h"
#include "../datastr/graph/UpdateableGraph.h"


/** The factor that the weights from the input are multiplied with. */
double weightFactor = 1;

/** Returns a random node identifier between 0 and n-1. */
inline NodeID randomNodeID(NodeID n) {
    return (NodeID)(rand() / (double)(RAND_MAX+1.0) * n);
}

/** Returns a random edge weight. */
inline EdgeWeight randomEdgeWeight() {
    const EdgeWeight maxEdgeWeight = (EdgeWeight)1e2;
    return (EdgeWeight)(rand() / (double)(RAND_MAX+1.0) * maxEdgeWeight);
}

/**
 * Removes parallel edges and self loops from the input and writes
 * the unique edges to the output.
 */
void removeDuplicates(vector<CompleteEdge>& input, vector<CompleteEdge>& output) {
    // sort the edge list first by source, then by target and weight
    sort(input.begin(), input.end());
    
    output.reserve( input.size() );

    VERBOSE( OtherEdgeID selfLoops = 0 );
    VERBOSE( OtherEdgeID parallelEdges = 0 );
    VERBOSE( OtherEdgeID weightDependsOnDir = 0 );
    VERBOSE( OtherEdgeID mergeTwoOneWays = 0 );
    VERBOSE( OtherEdgeID downgradeToOneWay = 0 );

    CompleteEdge lastEdge(0,0,0,false,true,true);
    for (OtherEdgeID i=0; i<input.size(); i++) { // scan the original edge list
        CompleteEdge currentEdge = input[i];
        if (currentEdge.source() == currentEdge.target()) {
            // skip self loop
            VERBOSE( selfLoops++ );
            VERBOSE2(cout << "self-loop " << currentEdge << " skipped." << endl); 
        }
        else {
            if ((currentEdge.source() != lastEdge.source()) ||
                (currentEdge.target() != lastEdge.target())) {
                
                lastEdge = currentEdge;
                // add unique edge to the new edge list (that doesn't contain duplicates)
                output.push_back( currentEdge );
            }
            else {
                assert( i > 0 );
                // Does the current edge point in another direction
                // than the previous edge ?
                if (((currentEdge.isDirected(0)) && (! lastEdge.isDirected(0))) ||
                    ((currentEdge.isDirected(1)) && (! lastEdge.isDirected(1)))) {
                    
                    if (currentEdge.weight() == lastEdge.weight()) {
                        // same weight
                        // -> sorting order implies that both edges are one-way and
                        //    complementary
                        // -> don't add the current edge, but make the last edge two-way
                        assert( (currentEdge.isDirected(0) != lastEdge.isDirected(0)) &&
                                (currentEdge.isDirected(1) != lastEdge.isDirected(1)) &&
                                (currentEdge.isDirected(0) != currentEdge.isDirected(1)) );
                        output.back().makeTwoWay();
                        lastEdge.makeTwoWay();
                        VERBOSE( mergeTwoOneWays++ );
                        continue;
                    }
                    
                    if (currentEdge.isBidirected()) {
                        // Since the last (one-way) edge is shorter, 
                        // the current (two-way) edge is only useful 
                        // for the opposite direction
                        assert( currentEdge.weight() > lastEdge.weight() );                        
                        if (lastEdge.isDirected(0)) currentEdge.makeOneWay(1);
                        if (lastEdge.isDirected(1)) currentEdge.makeOneWay(0);
                        VERBOSE( downgradeToOneWay++ );
                    }
                    
                    // set both direction flags to indicate that no further
                    // edge is needed for the current s-t-pair
                    lastEdge.makeTwoWay();
                    
                    // add unique edge to the new edge list
                    output.push_back( currentEdge );
                    VERBOSE( weightDependsOnDir++ );
                    VERBOSE2(cout << "parallel edge " << currentEdge
                             << " kept; predecessor: " << input[i-1] << endl);
                }
                else {
                    // skip parallel edge
                    // note: since parallel edges are sorted by weight,
                    //       it is ensured that the shorter one "survives"
                    VERBOSE( parallelEdges++ );
                    VERBOSE2(cout << "parallel edge " << currentEdge
                             << " skipped; predecessor: " << input[i-1] << endl);
                }
            }           
        }
    }

    // print statistics
    VERBOSE( cout << "Start with " << input.size() << " edges." << endl );
    VERBOSE( if (selfLoops > 0) cout << selfLoops << " self-loop(s) skipped." << endl );
    VERBOSE( if (parallelEdges > 0) cout << parallelEdges << " parallel edge(s) skipped." << endl );
    VERBOSE( if (weightDependsOnDir > 0) cout << weightDependsOnDir 
             << " parallel edge(s) kept (different directions -> different weights)." << endl );
    VERBOSE( if (downgradeToOneWay > 0) cout << "  " << downgradeToOneWay 
             << " two-way edge(s) downgraded to one-way "
             << "(since the other direction is covered by a shorter edge)." << endl );
    VERBOSE( if (mergeTwoOneWays > 0) cout << mergeTwoOneWays 
             << " one-way edge(s) merged with an opposite edge of the same length." << endl );
}





/**
 * Reads a list of edges from the given stream and adds them to the given vector.
 * @return the number of nodes
 */
NodeID readGraphFromStream(istream &in, const bool adaptWeight, const bool addIsolatedNode,
                           vector<CompleteEdge>& edgeList) {
    NodeID n; // number of nodes
    OtherEdgeID m; // number of edges
    NodeID source, target;
    int dir; // direction (0 = open, 1 = forward, 2 = backward, 3 = closed) [adopted from PTV]
    char dummy;
    
    // directed ?
    bool directed = false;
    if (in.peek() == 'd') {
        directed = true;
        in >> dummy;
    }

    in >> n >> m;

    VERBOSE( cout << "Import Graph: n = " << n << ", m = " << m << endl );
    
    if (adaptWeight) {
        // read the edge list a first time 
        // in order to determine the maximum weight in the input
        
        double weight; // the current weight read from the input
        double maxWeight = 0; // the maximum weight
        double weightSum = 0; // the sum of all weights
        
        // if false, no overflow (w.r.t. the maximum integer) occured 
        //           when computing the sum of all weights
        // if true, there might be an overflow
        bool overflow = false;
                               

        for (OtherEdgeID i=0; i<m; i++) { // scan the edge list
            in >> source >> target >> weight;
            if (directed) in >> dir;
            if (weight > maxWeight) maxWeight = weight; // update the maximum weight
            
            if ((weightSum > Weight::MAX_INTEGER / 2) || 
                 (weight > Weight::MAX_INTEGER / 2)) overflow = true; 
            weightSum += weight; // update the sum of weights
        }
        
        in.seekg(0); // re-read the input
        if (directed) in >> dummy;
        in >> n >> m;
        
        // the allowed weight per edge;
        // this upper bound is based on the fact that a shortest path always consists
        // of less than n edges
		double allowedWeight = (double)(Weight::MAX_INTEGER / n);
        //double allowedWeight = trunc(Weight::MAX_INTEGER / n);
        
        // if multiplied by this factor, all edge weights (incl. the maximum one)
        // are less than or equal to the allowed weight
        weightFactor = allowedWeight / maxWeight;

        VERBOSE( cout << setprecision(20) << "  adapt weight: allowed weight = " << allowedWeight
                 << ", max weight in input = " << maxWeight << ", factor = " << weightFactor <<  endl );

        if (! overflow) {
            // if there was no overflow, we can try to use the sum of weights
            // in order to find a less restrict bound for the edge weights
            
            // if the edge weights are multiplied by this factor,
            // the sum of all edge weights does not exceed the maximum integer;
            // note that a shortest path does not contain any edge twice
            double weightFactor2 = Weight::MAX_INTEGER / weightSum;
            
            // if the second factor is better than the first one, use the second one;
            // "better" means that we exhaust the given range more carefully so that
            // the loss of precision by this initial rounding is less serious
            if (weightFactor2 > weightFactor) {
                weightFactor = weightFactor2;

                VERBOSE( cout << "                sum of weights = " << weightSum 
                              << ", factor = " << weightFactor2 << endl );
            }
        }       
    }
    
    VERBOSE( OtherEdgeID onewayStreets = 0 );
    edgeList.reserve(2*m); // each edge will be added twice (once for each direction)
	
	/*graph format is one-way edge, revised by ddx*/
	//edgeList.reserve(m);
        
    for (OtherEdgeID i=0; i<m; i++) {
        EdgeWeight weight;
        
        in >> source >> target;

        if (adaptWeight) {
            double tmp;
            in >> tmp;
            weight = static_cast<EdgeWeight>( (int)(tmp * weightFactor) );
        }
        else {
            in >> weight;
        }
        
        dir = 0;
        if (directed) in >> dir;
        bool forward = true;
        bool backward = true;
        if (dir == 1) backward = false;
        if (dir == 2) forward = false;
        VERBOSE( if ((dir == 1) || (dir == 2)) onewayStreets++ );
        // note: all CLOSED roads (dir = 3) are considered as OPEN !

        CompleteEdge edge1(source, target, weight, false, forward, backward);
		//commentted by ddx
        CompleteEdge edge2(target, source, weight, false, backward, forward);
        edgeList.push_back(edge1);
		// commentted by ddx
        edgeList.push_back(edge2);
        
        // Check wheter end of line is reached
        if ( in.peek() !=  '\n' && in.peek() != '\r' )
        {
          cerr << "After edge number " << i << " is no newline." << endl;
          cerr << "Input format error, aborting." << endl;
          exit(1);
        }

    }
    
    VERBOSE( cout << onewayStreets << " oneway street(s)." << endl );

    if (addIsolatedNode) n++;
    return n;
}
    
/**
 * Imports a list of edges from a given input stream and
 * builds an UpdateableGraph. The node levels are read from
 * the specified file 'filenameHwyNodes', otherwise the levels are set to n.
 * Note: To use the graph to create a contraction hierarchy (../processing/ConstructCH.h)
 *       initially all nodes must have level n in the graph. The node order is specified
 *       separately.
 */
datastr::graph::UpdateableGraph* importGraphListOfEdgesUpdateable(istream &inGraph, const bool adaptWeight,
                                        const bool addIsolatedNode,  string filenameHwyNodes) {
    vector<CompleteEdge> edgeList;
    const NodeID n = readGraphFromStream(inGraph, adaptWeight, addIsolatedNode, edgeList);
    vector<CompleteEdge> edgeListWithoutDuplicates;
    removeDuplicates(edgeList, edgeListWithoutDuplicates);
    
    // Free edgeList because it is no longer needed
    edgeList.clear();
    vector<CompleteEdge>().swap(edgeList);
    
	/*vector <CompleteEdge>::iterator it = edgeList.begin();
	fstream edgefile("edge.txt", ios::out);
	edgefile << edgeList.size() << endl;
	for (int i= 0; it != edgeList.end(); it++, i++)
		edgefile << i << " " << *it << endl;*/

    // default level is for all n nodes is n
    vector<LevelID> nodeLevels(n, n);
    //if (filenameHwyNodes.c_str() != "")
    //{
    //    VERBOSE( cout << "open file '" << filenameHwyNodes.c_str() << "' to read hwy node levels" << endl );
    //    ifstream in(filenameHwyNodes.c_str());
    //    if (! in) {
    //        VERBOSE( cout << "failed" << endl );
    //        exit(-1);
    //    }
    //    NodeID n2;
    //    LevelID max = 0;
    //    readPrimitive(in, n2);
    //    assert( n == n2 );
    //    for (NodeID u = 0; u < n; u++) {
    //        LevelID k;
    //        readPrimitive(in, k);
    //        if (k > max) max = k;
    //        nodeLevels[u] = k;
    //    }
    //    VERBOSE( cout << "max hwy node level = " << max << endl );
    //}

    return new datastr::graph::UpdateableGraph(edgeListWithoutDuplicates, nodeLevels);
}

/**
 * Imports a list of edges from a given input stream and
 * builds an UpdateableGraph. The node levels are read from
 * the specified file 'filenameHwyNodes', otherwise the levels are set to n.
 * Note: To use the graph to create a contraction hierarchy (../processing/ConstructCH.h)
 *       initially all nodes must have level n in the graph. The node order is specified
 *       separately.
 */
datastr::graph::UpdateableGraph* importGraphListOfEdgesUpdateable(istream &inGraph, const bool adaptWeight,
                                        const bool addIsolatedNode) {
    vector<CompleteEdge> edgeList;
    const NodeID n = readGraphFromStream(inGraph, adaptWeight, addIsolatedNode, edgeList);
    vector<CompleteEdge> edgeListWithoutDuplicates;
    removeDuplicates(edgeList, edgeListWithoutDuplicates);

    return new datastr::graph::UpdateableGraph(edgeListWithoutDuplicates, n);
}

datastr::graph::UpdateableGraph* importGraphListOfEdgesUpdateable(istream &inGraph, const bool adaptWeight,
                                        const bool addIsolatedNode, fstream &cordfile) {
    vector<CompleteEdge> edgeList;
    const NodeID n = readGraphFromStream(inGraph, adaptWeight, addIsolatedNode, edgeList);
    vector<CompleteEdge> edgeListWithoutDuplicates;
    removeDuplicates(edgeList, edgeListWithoutDuplicates);

    return new datastr::graph::UpdateableGraph(edgeListWithoutDuplicates, n);
}

datastr::graph::UpdateableGraph* importGraphListOfEdgesUpdateable(istream &inGraph, const bool adaptWeight,
                                        const bool addIsolatedNode, string filenameHwyNodes, 
										fstream &cordfile, vector<CoordinateType> &xcord, vector <CoordinateType> &ycord,
										CoordinateType &xrange, CoordinateType &yrange) {
    vector<CompleteEdge> edgeList;
    const NodeID n = readGraphFromStream(inGraph, adaptWeight, addIsolatedNode, edgeList);
    vector<CompleteEdge> edgeListWithoutDuplicates;
    removeDuplicates(edgeList, edgeListWithoutDuplicates);
    
	/*vector <CompleteEdge>::iterator it = edgeList.begin();
	fstream edgefile("edge.txt", ios::out);
	edgefile << edgeList.size() << endl;
	for (int i= 0; it != edgeList.end(); it++, i++)
		edgefile << i << " " << *it << endl;*/

    // default level is for all n nodes is n
    vector<LevelID> nodeLevels(n, n);

    if (filenameHwyNodes != "")
    {
        VERBOSE( cout << "open file '" << filenameHwyNodes.c_str() << "' to read hwy node levels" << endl );
        ifstream in(filenameHwyNodes.c_str());
        if (! in) {
            VERBOSE( cout << "failed" << endl );
            exit(-1);
        }
        NodeID n2;
        LevelID max = 0;
        readPrimitive(in, n2);
        assert( n == n2 );
        for (NodeID u = 0; u < n; u++) {
            LevelID k;
            readPrimitive(in, k);
            if (k > max) max = k;
            nodeLevels[u] = k;
        }
        VERBOSE( cout << "max hwy node level = " << max << endl );
    }

    return new datastr::graph::UpdateableGraph(edgeListWithoutDuplicates, nodeLevels, cordfile, xcord, ycord, xrange, yrange);
}

#endif // CREATEGRAPH_H
