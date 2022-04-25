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

#ifndef _JUNCTION_H
#define _JUNCTION_H

#include "Graph.h"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp> 

class Junction {

    public:
        Junction(std::string networkName, int numNodes, int numEdges);
        Junction() {};
        void buildJunction(Graph_RoadNetwork& graph);
        std::string getNetworkName();
        int getNumNodes();
        int getNumEdges();
        double computeIndexSize();
        double computeMemoryUsage();
        void verifyJunctionsAndDistances(Graph_RoadNetwork& graph);
        void printGraphStats();
        
        std::vector<bool> isNoThruRoadNode;
        std::vector<int> noThruRoadID; // Indicates which no thru road each node belongs on (-1 if it doesn't belong on one)
        std::vector<IntDistancePair> nextJunctionNodeAndDistance;
        // Note: Size of vector is number of edges and distance 
        // component includes that edges weight
        
        // Statistics Collection
        std::vector<int> intraJuctionRoadLengths;
        std::vector<int> noThruRoadLengths;
        unsigned int numTotalNoThruRoadNodes = 0, numNoThruRoads = 0;        

private:
        friend class boost::serialization::access;
        
        std::string networkName;
        unsigned int numNodes;
        unsigned int numEdges;
        
        // Non-Serialized Members

        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->networkName;
            ar & this->numNodes;
            ar & this->numEdges;
            ar & this->isNoThruRoadNode;
            ar & this->noThruRoadID;
            ar & this->nextJunctionNodeAndDistance;
        }    
};


#endif // _JUNCTION_H