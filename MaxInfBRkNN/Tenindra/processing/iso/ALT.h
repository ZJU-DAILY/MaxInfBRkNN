/* Copyright (C) 2015 Tenindra Abeywickrama
 *
 * This file is part of Road Network kNN Experimental Evaluation.
 *
 * Road Network kNN Experimental Evaluation is free software; you can
 * redistribute it and/or modify it under the terms of the GNU Affero 
 * General Public License as published by the Free Software Foundation; 
 * either version 3 of the License, or (at your option) any later version.
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

#ifndef _ALT_H
#define _ALT_H

#include "Graph.h"
#include "Path.h"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>

enum LANDMARK_TYPE
{
    RANDOM
};

class ALT {

    public:
        ALT(std::string networkName, int numNodes, int numEdges);
        ALT() {};
        void buildALT(Graph_RoadNetwork& graph, LANDMARK_TYPE landmarkType, unsigned int numLandmarks);
        std::string getNetworkName();
        unsigned int getNumNodes();
        unsigned int getNumEdges();
        double computeIndexSize();
        EdgeWeight getLowerBound(NodeID s, NodeID t);
        void getLowerAndUpperBound(NodeID s, NodeID t, EdgeWeight& lb, EdgeWeight& ub);
        Path findShortestPath(Graph_RoadNetwork& graph, NodeID source, NodeID target, std::vector<NodeID>& shortestPathTree);
        PathDistance findShortestPathDistance(Graph_RoadNetwork& graph, NodeID source, NodeID target);
        LANDMARK_TYPE getLandmarkType(std::string selectionMethod, bool& success);
        
    private:
        friend class boost::serialization::access;
        
        std::string networkName;
        unsigned int numNodes;
        unsigned int numEdges;
        unsigned int numLandmarks;
        std::vector<NodeID> landmarks;
        std::vector<int> vertexFromLandmarkDistances;
        
        // Non-Serialized Members
    
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & networkName;
            ar & numNodes;
            ar & numEdges;
            ar & numLandmarks;
            ar & landmarks;
            ar & vertexFromLandmarkDistances;
        }    
};


#endif // _ALT_H