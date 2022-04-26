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

#include "serialization.h"

#include "../processing/Graph.h"
#include "../processing/DynamicGraph.h"
//#include "../processing/Gtree.h"
//#include "../processing/ROAD.h"
#include "../processing/MortonList.h"
#include "../processing/Quadtree.h"
#include "../processing/SimpleQuadtree.h"
#include "../processing/Junction.h"
#include "../processing/StaticRtree.h"
//#include "../processing/ALT.h"
#include "../processing/V2PNVD.h"
#include "StopWatch.h"

#include <iostream>
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

namespace serialization
{
    template <typename T>
    T getIndexFromBinaryFile(std::string idxFilePath) {
        T index;
        //StopWatch sw;
      
        std::ifstream binaryIndexFS(idxFilePath,std::ios::in | std::ios::binary);
        if (binaryIndexFS.is_open()) {
            boost::archive::binary_iarchive ia(binaryIndexFS);
            //sw.start();
            ia >> index;
            //sw.stop();
        } else {
            std::cerr << "Cannot open binary index file " << idxFilePath << std::endl;
            exit(1);
        }
        binaryIndexFS.close();
        
        //std::cout << "Binary index deserialized in " << sw.getTimeMs() << "ms" << std::endl;        
        
        return index;
    }
    
    template <typename T>
    void populateIndexFromBinaryFile(std::string idxFilePath, T& index) {
        //StopWatch sw;
      
        std::ifstream binaryIndexFS(idxFilePath,std::ios::in | std::ios::binary);
        if (binaryIndexFS.is_open()) {
            boost::archive::binary_iarchive ia(binaryIndexFS);
            //sw.start();
            ia >> index;
            //sw.stop();
        } else {
            std::cerr << "Cannot open binary index file " << idxFilePath << std::endl;
            exit(1);
        }
        binaryIndexFS.close();
        
        //std::cout << "Binary index deserialized in " << sw.getTimeMs() << "ms" << std::endl;        
    }
    
    template <typename T>
    void outputIndexToBinaryFile(T& index, std::string idxFilePath) {
        //StopWatch sw;

        std::ofstream binaryIndexFS(idxFilePath,std::ios::out | std::ios::binary | std::ios::trunc);
        {
            boost::archive::binary_oarchive oa(binaryIndexFS);
            //sw.start();
            oa << index;
            //sw.stop();
        }
        binaryIndexFS.close();

        //std::cout << "Binary index serialized in " << sw.getTimeMs() << "ms" << std::endl;        
    }
    
    // Instatiate Templates
    //below this for  hash index of nvd and v2p
    template V2PNVDHashIndex getIndexFromBinaryFile<V2PNVDHashIndex>(std::string idxFilePath);
    template void outputIndexToBinaryFile<V2PNVDHashIndex>(V2PNVDHashIndex& index, std::string idxFilePath);
    template HybridPOINNHashIndex getIndexFromBinaryFile<HybridPOINNHashIndex>(std::string idxFilePath);
    template void outputIndexToBinaryFile<HybridPOINNHashIndex>(HybridPOINNHashIndex& index, std::string idxFilePath);

    template SocialGraphBinary getIndexFromBinaryFile<SocialGraphBinary>(std::string dataFilePath);
    template void outputIndexToBinaryFile<SocialGraphBinary>(SocialGraphBinary& data, std::string dataFilePath);

    template CheckInBinary getIndexFromBinaryFile<CheckInBinary>(std::string dataFilePath);
    template void outputIndexToBinaryFile<CheckInBinary>(CheckInBinary& data, std::string dataFilePath);


    template Graph_RoadNetwork getIndexFromBinaryFile<Graph_RoadNetwork>(std::string idxFilePath);
    template void outputIndexToBinaryFile<Graph_RoadNetwork>(Graph_RoadNetwork& index, std::string idxFilePath);
    template DynamicGraph getIndexFromBinaryFile<DynamicGraph>(std::string idxFilePath);
    template void outputIndexToBinaryFile<DynamicGraph>(DynamicGraph& index, std::string idxFilePath);
   /* template Gtree getIndexFromBinaryFile<Gtree>(std::string idxFilePath);
    template void outputIndexToBinaryFile<Gtree>(Gtree& index, std::string idxFilePath);
    template OccurenceList getIndexFromBinaryFile<OccurenceList>(std::string idxFilePath);
    template void outputIndexToBinaryFile<OccurenceList>(OccurenceList& index, std::string idxFilePath);*/
    /*template ROAD getIndexFromBinaryFile<ROAD>(std::string idxFilePath);
    template void outputIndexToBinaryFile<ROAD>(ROAD& index, std::string idxFilePath);
    template AssociationDirectory getIndexFromBinaryFile<AssociationDirectory>(std::string idxFilePath);
    template void outputIndexToBinaryFile<AssociationDirectory>(AssociationDirectory& index, std::string idxFilePath);*/
    template Quadtree getIndexFromBinaryFile<Quadtree>(std::string idxFilePath);
    template void outputIndexToBinaryFile<Quadtree>(Quadtree& index, std::string idxFilePath);
    template SILCPathOracle getIndexFromBinaryFile<SILCPathOracle>(std::string idxFilePath);
    template void outputIndexToBinaryFile<SILCPathOracle>(SILCPathOracle& index, std::string idxFilePath);
    template SimpleQuadtree getIndexFromBinaryFile<SimpleQuadtree>(std::string idxFilePath);
    template void outputIndexToBinaryFile<SimpleQuadtree>(SimpleQuadtree& index, std::string idxFilePath);
    template Junction getIndexFromBinaryFile<Junction>(std::string idxFilePath);
    template void outputIndexToBinaryFile<Junction>(Junction& index, std::string idxFilePath);
    template StaticRtree getIndexFromBinaryFile<StaticRtree>(std::string idxFilePath);
    template void outputIndexToBinaryFile<StaticRtree>(StaticRtree& index, std::string idxFilePath);
    /*template ALT getIndexFromBinaryFile<ALT>(std::string idxFilePath);
    template void outputIndexToBinaryFile<ALT>(ALT& index, std::string idxFilePath);*/
    
    template void populateIndexFromBinaryFile<SILCPathOracle>(std::string idxFilePath, SILCPathOracle& index);
    template void populateIndexFromBinaryFile<Junction>(std::string idxFilePath, Junction& index);
    
}