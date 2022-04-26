/* Original Work Copyright (C) 2012 Lingkun Wu, Xiaokui Xiao, Dingxiong Deng, Gao Cong, Andy Diwen Zhu, Shuigeng Zhou
 * Modified Work Copyright (C) 2015 Tenindra Abeywickrama
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

#include "ShortestPathWrapper.h"
#include "INE.h"

#include "StaticRtree.h"
#include "../tuple/KnnQueryTuple.h"
#include "../utility/StopWatch.h"
#include "../utility/utility.h"
#include "../utility/serialization.h"

#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

#include "../external/config.h"
#include "../external/stats/utils.h"
#include "../external/io/createGraph.h"
#include "../external/io/output.h"
// #include "../external/io/coordinateio.h"
#include "../external/processing/ConstructCH.h"
#include "../external/processing/DijkstraCH.h"
#include "../external/processing/TransitNode.h"
#include "../external/datastr/graph/SearchGraph.h"
#include "../external/io/CoordinateLoader.h"

#include <iostream>
#include <fstream>

typedef processing::ConstructCH<datastr::graph::UpdateableGraph, false, false, SAVE_STATS_NONE> ProcessingNodeOrder;

void ShortestPathWrapper::buildIntermediaryGraph(std::string graphFilePath, std::string coordinateFilePath, std::string bgrIntFilePath, std::string bcoIntFilePath)
{
    datastr::graph::UpdateableGraph* graph;
    std::ifstream grFile(graphFilePath, std::ios::in);
    if (grFile.is_open()) {
         graph = importGraphListOfEdgesUpdateable(grFile, false, false, "");
    } else {
        std::cerr << "Cannot open text graph file " << graphFilePath << std::endl;
        std::exit(1);
    }
    grFile.close();
    
    std::ofstream bgrFile(bgrIntFilePath, std::ios::out | std::ios::binary);
    if (bgrFile.is_open()) {
         graph->serialize(bgrFile);
    } else {
        std::cerr << "Cannot write to binary graph file " << bgrIntFilePath << std::endl;
        std::exit(1);
    }
    bgrFile.close();

    delete graph;

    std::vector<CoordinateType> g_stXcord, g_stYcord;
    CoordinateLoader::readCoordinate(coordinateFilePath,g_stXcord,g_stYcord);
    CoordinateLoader::serializeCoordinate(bcoIntFilePath,g_stXcord,g_stYcord);
}

void ShortestPathWrapper::buildCH(std::string bgrIntFilePath, std::string chFilePath, std::string noFilePath)
{
    datastr::graph::UpdateableGraph* graph;
    std::ifstream bgrFile(bgrIntFilePath, std::ios::in);
    if (bgrFile.is_open()) {
         graph = new datastr::graph::UpdateableGraph(bgrFile);
    } else {
        std::cerr << "Cannot open binary graph file " << bgrIntFilePath << std::endl;
        std::exit(1);
    }
    bgrFile.close();
    
    // The stringstream ss consists of a "-" separated list of all parameters and can be used as prefix
    // for filenames if "x" is specified as parameter value.
    std::stringstream ss;
    ss << "choutput/exact-"; 

    // Initialize pseudorandom number generator allow reproduction of results.
    // Currently, no random data is used during node ordering but it does not harm
    // either to specifiy the random seed.   
    std::srand(22);
    
    // Object for node ordering.
    ProcessingNodeOrder c(graph);
            
    // Prepare output files for witness and shortcut information that can 
    // be used in a later hierarchy construction step. 
    // To export the wheter witness information, a template parameter
    // in the ProcessingNodeOrder class needs to be set. There is
    // currently node switch in ../config.h
    std::ofstream outShortcuts;
    std::ofstream outWitnesses;
    if ( c.savesShortcutsWitnesses() )
    {                    
        outShortcuts.open((ss.str()+".shortcuts").c_str());
        if (!outShortcuts.is_open()) { 
            std::cerr << "Cannot write to " << (ss.str()+".shortcuts") << std::endl; 
            exit(1); 
        }
        c.storeShortcuts(outShortcuts);
        outWitnesses.open((ss.str()+".witnesses").c_str());
        if (!outWitnesses.is_open()) { 
            std::cerr << "Cannot write to " << (ss.str()+".witnesses") << std::endl; 
            exit(1);
        }
        c.storeWitnesses(outWitnesses);
    }
    
    /*********************************************
     * Node ordering and hierarchy construction. *
     *********************************************/
            
    /**
        we use the aggressive parameter setting "EVSQL" in the preprocessing for CH,
        and have not use betweenness value and reach value
        commented by dingxiong
        */

    // stores most of the relevant parameters for node ordering
    // especially the coefficients for the linear combination of 
    // the priority terms and the limits to the local searches
    ProcessingNodeOrder::WeightCalculation calc; 

    // coefficient for edge difference
    calc.edgeDiffMult = 190;
                
    // coefficient for search space size of local searches for contraction
    calc.searchSpaceMult = 1;

    // coefficient for Voronoi region size
    calc.voronoiMult = 60;
    
    // coefficient for upper bound on search paths hops
    calc.searchPathHopBorderMult = 145;


    // limit of settled nodes in local searches during weight calculation
    calc.maxSettledApprox = 1000;

    // limit of settled nodes in local searches during contraction
    calc.maxSettledElim = 1000;

    // lazy updates, parameter also specifies the check intervals
    // for complete recalculation of priority queue
    calc.lazyUpdateRecalcLimit = 1000;

    // perform update on all affected nodes after contraction
    // this only works with hop-limits and is really slow
    // should only be used for testing
    calc.updateHops = false;
    // omit local edge reduction
    calc.localReduceEdges = !calc.localReduceEdges;
    
    /*
        below is  another parameters option which can be combined with above setting,
        if a specified parameter is needed, just uncomment these relative codes
        commented by dingxiong
    */

    // coefficient for deleted neighbors
    //calc.delNeighbMult = 120; 

    // coefficient for number of new edges
    //calc.newEdgesMult = 100;
    
    // This is the main step, a node order is created including a contraction hierarchy.
    StopWatch sw;
    sw.start();
    c.createHierarchy(calc, NULL, NULL);
    sw.pause();

    // Now write the node order to a file.
    std::ofstream noFile(noFilePath);
    if (noFile.is_open()) { 
        c.writeLevels(noFile);
    } else {
        std::cerr << "Cannot write to " << noFilePath << std::endl; 
        std::exit(1);
    }        
    noFile.close();
    if (outShortcuts.is_open()) {
        outShortcuts.close();
    }
    if (outWitnesses.is_open()) {
        outWitnesses.close();
    }
        
    // output search-graph to binary format (depends on switches in config.h)
    sw.resume();
    datastr::graph::SearchGraph* searchGraph = new datastr::graph::SearchGraph(graph, datastr::graph::SGNO_ORIGINAL);
    sw.stop();
    
    constructionTimeMs += sw.getTimeMs();
    memoryUsageMB = searchGraph->getIndexSize();
    
    std::ofstream chFile(chFilePath, std::ios::out | std::ios::binary);
    if (chFile.is_open()) { 
        searchGraph->serialize(chFile);
    } else {
        std::cerr << "Cannot write to " << chFilePath << std::endl; 
        std::exit(1);
    }        
    chFile.close();

    delete graph;
    delete searchGraph;
}

void ShortestPathWrapper::buildTNR(std::string bgrIntFilePath, std::string bcoIntFilePath, std::string chFilePath, std::string tnrFilePath, int ascale)
{
    datastr::graph::UpdateableGraph* graph;
    std::ifstream bgrFile(bgrIntFilePath, std::ios::in | std::ios::binary);
    if (bgrFile.is_open()) {
         graph = new datastr::graph::UpdateableGraph(bgrFile);
    } else {
        std::cerr << "Cannot open binary graph file " << bgrIntFilePath << std::endl;
        std::exit(1);
    }
    bgrFile.close();
    
    datastr::graph::SearchGraph* searchGraph;
    std::ifstream sgrFile(chFilePath, std::ios::in | std::ios::binary);
    if (sgrFile.is_open()) {
         searchGraph = new datastr::graph::SearchGraph(sgrFile);
    } else {
        std::cerr << "Cannot open search graph file " << chFilePath << std::endl;
        std::exit(1);
    }
    sgrFile.close();
    
    std::vector<CoordinateType> g_stXcord, g_stYcord;
    CoordinateLoader::deserializeCoordinate(bcoIntFilePath,g_stXcord,g_stYcord);
    CoordinateType g_Range = CoordinateLoader::getCoordinateRange(g_stXcord,g_stYcord);
    
//     std::ifstream bcoFile(bcoIntFilePath, std::ios::in | std::ios::binary);
//     if (bcoFile.is_open()) {
//         deserializeCoordinate(bcoFile, g_stXcord, g_stYcord);
//         getCoordinateRange(g_stXcord, g_stYcord);
//     } else {
//         std::cerr << "Cannot open search graph file " << bcoIntFilePath << std::endl;
//         std::exit(1);
//     }
//     bcoFile.close();

    /***********************
     * Build transit nodes *
     ***********************/

    //set temp output file
    stringstream ss;
    std::string transitFile, accessNodeFile, boundFile;
    ss << "output/transitnum_a=" << ascale << ".ch";
    transitFile = ss.str();
    ss.clear();
    ss << "output/accessnode_a==" << ascale << ".ch";
    accessNodeFile = ss.str();
    ss.clear();
    ss << "output/transitnode_BAY_t_test_a=" << ascale << ".ch";
    boundFile = ss.str();
    ss.clear();
    fstream transitnum(transitFile, std::ios::out | std::ios::binary);
    fstream accessnode(accessNodeFile, std::ios::out);
    fstream bound(boundFile, std::ios::out);

    StopWatch sw;
    sw.start();
    TransitNodeTest tnr(searchGraph,graph, ascale, g_Range);
    tnr.computeAccessNodeForEdge(graph, g_stXcord, g_stYcord, accessnode);
    tnr.computeAllCellTransitNode(g_stXcord, g_stYcord, bound);
    tnr.processEachNodeToCell(g_stXcord, g_stYcord);
    tnr.buildDistanceTable();
    sw.stop();

    constructionTimeMs += sw.getTimeMs();
    std::ofstream tnrFile(tnrFilePath, std::ios::out | std::ios::binary);
    if (tnrFile.is_open()) { 
        tnr.serialize(tnrFile);
    } else {
        std::cerr << "Cannot write to " << tnrFilePath << std::endl; 
        std::exit(1);
    }        
    tnrFile.close();
    memoryUsageMB = searchGraph->getIndexSize() + graph->getIndexSize() + getFileSizeMB(tnrFilePath);
    // Note: TNR query also needs to coordinate vectors (for locating a vertex's cell)
    memoryUsageMB += sizeof(CoordinateType)*(g_stXcord.size()+g_stYcord.size())/(1024*1024);
    
    std::cout << "Index Size = " << memoryUsageMB << "MB" << std::endl;
    
    delete graph;
    delete searchGraph;
}

// Get construction time of last constructed index (but only for CH and TNR)
double ShortestPathWrapper::getConstructionTimeMs()
{
    return constructionTimeMs;
}

// Get memory usage of last constructed index (but only for CH and TNR)
double ShortestPathWrapper::getIndexSizeMB()
{
    return memoryUsageMB;
}

double ShortestPathWrapper::getFileSizeMB(string filePath)
{
    double fileSizeBytes = 0;
    std::ifstream fs(filePath,std::ios::ate | std::ios::binary);
    fileSizeBytes = fs.tellg();
    return fileSizeBytes/(1024*1024);
}

void ShortestPathWrapper::runIERQueries(Graph_RoadNetwork& myGraph, std::string method, std::string bgrIntFilePath, std::string bcoIntFilePath, std::string chFilePath, std::string tnrFilePath,
                                        std::vector<int>& branchFactors, std::vector<NodeID>& queryNodes, std::vector<int>& kValues, std::size_t numSets,
                                        std::vector<double> objDensities, std::vector<std::string> objTypes, std::vector<int> objVariable,
                                        std::string filePathPrefix, std::string statsOutputFile, bool verifyKNN, std::vector<std::string> specialFields) {
    
    datastr::graph::UpdateableGraph* graph;
    std::ifstream bgrFile(bgrIntFilePath, std::ios::in | std::ios::binary);
    if (bgrFile.is_open()) {
         graph = new datastr::graph::UpdateableGraph(bgrFile);
    } else {
        std::cerr << "Cannot open binary graph file " << bgrIntFilePath << std::endl;
        std::exit(1);
    }
    bgrFile.close();
    
    datastr::graph::SearchGraph* searchGraph;
    std::ifstream sgrFile(chFilePath, std::ios::in | std::ios::binary);
    if (sgrFile.is_open()) {
         searchGraph = new datastr::graph::SearchGraph(sgrFile);
    } else {
        std::cerr << "Cannot open search graph file " << chFilePath << std::endl;
        std::exit(1);
    }
    sgrFile.close();
    
    std::vector<CoordinateType> g_stXcord, g_stYcord;
    CoordinateLoader::deserializeCoordinate(bcoIntFilePath,g_stXcord,g_stYcord);

    TransitNodeTest tnr;
    processing::DijkstraCH<datastr::graph::SearchGraph, NormalPQueue, 2, true> dijkstraCHTest;
    
    std::string desc = "";
    if (method == constants::IER_TNR_KNN_QUERY) {
        std::ifstream tnrFile(tnrFilePath, std::ios::in | std::ios::binary);
        if (!tnrFile.is_open()) {
            std::cerr << "Cannot open TNR index file " << tnrFilePath << std::endl;
            std::exit(1);
        }
        tnr.loadIndex(searchGraph, graph, tnrFile);
        desc = "using CH+TNR";
    } else {
        dijkstraCHTest.loadGraph(searchGraph);
        desc = "using CH";
    }

    std::vector<NodeID> kNNs, ineKNNs;
    std::vector<EdgeWeight> kNNDistances, ineKNNDistances;

    std::string objSetType;
    double objSetDensity;
    int objSetSize, objSetVariable;
    
    std::vector<std::string> parameterKeys, parameterValues;
    
    StopWatch sw;
    double totalQueryTime;
    int totalQueries = numSets*queryNodes.size();
    int totalObjects;
    
    INE ine;
    std::string message;
    Coordinate queryNodeX, queryNodeY;
    
    for (std::size_t i = 0; i < objTypes.size(); ++i) {
        for (std::size_t m = 0; m < branchFactors.size(); ++m) {
            parameterKeys.clear();
            parameterValues.clear();
            parameterKeys.push_back("branchfactor");
            parameterValues.push_back(std::to_string(branchFactors[m]));
            for (std::size_t j = 0; j < objDensities.size(); ++j) {
                for (std::size_t n = 0; n < objVariable.size(); ++n) {
                    for (std::size_t k = 0; k < kValues.size(); ++k) {
                        totalQueryTime = 0;
                        totalObjects = 0;
                        for (std::size_t l = 0; l < numSets; ++l) {
                            std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(myGraph.getNetworkName(),constants::OBJ_IDX_RTREE,objTypes[i],objDensities[j],objVariable[n],l,parameterKeys,parameterValues);
                            StaticRtree rtree = serialization::getIndexFromBinaryFile<StaticRtree>(objIdxFilePath);
                            totalObjects += rtree.getObjSetSize();
                            
                            if (verifyKNN) {
                                myGraph.resetAllObjects();
                                std::string objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(myGraph.getNetworkName(),objTypes[i],objDensities[j],objVariable[n],l);
                                myGraph.parseObjectFile(objSetFile,objSetType,objSetDensity,objSetVariable,objSetSize);
                            }
                            
                            for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                                kNNs.clear();
                                kNNDistances.clear();
                                kNNs.reserve(kValues[k]);
                                kNNDistances.reserve(kValues[k]);
                                if (method == constants::IER_TNR_KNN_QUERY) {
                                    if (myGraph.getEdgeType() == constants::DISTANCE_WEIGHTS) {
                                        sw.reset();
                                        sw.start();
                                        myGraph.getCoordinates(*queryNodeIt,queryNodeX,queryNodeY);
                                        tnr.getKNNs(rtree,kValues[k],*queryNodeIt,kNNs,kNNDistances,queryNodeX,queryNodeY,g_stXcord,g_stYcord);
                                        sw.stop();
                                    } else {
                                        //assert(myGraph.getEdgeType() == constants::TIME_WEIGHTS);
                                        sw.reset();
                                        sw.start();
                                        myGraph.getCoordinates(*queryNodeIt,queryNodeX,queryNodeY);
                                        tnr.getKNNsByTravelTimes(rtree,kValues[k],*queryNodeIt,kNNs,kNNDistances,queryNodeX,queryNodeY,g_stXcord,g_stYcord,myGraph.getMaxGraphSpeedByEdge());
                                        sw.stop();
                                    }
                                } else {
                                    dijkstraCHTest.clear();
                                    if (myGraph.getEdgeType() == constants::DISTANCE_WEIGHTS) {
                                        sw.reset();
                                        sw.start();
                                        myGraph.getCoordinates(*queryNodeIt,queryNodeX,queryNodeY);
                                        dijkstraCHTest.getKNNs(rtree,kValues[k],*queryNodeIt,kNNs,kNNDistances,queryNodeX,queryNodeY);
                                        sw.stop();
                                    } else {
                                        //assert(myGraph.getEdgeType() == constants::TIME_WEIGHTS);
                                        sw.reset();
                                        sw.start();
                                        myGraph.getCoordinates(*queryNodeIt,queryNodeX,queryNodeY);
                                        dijkstraCHTest.getKNNsByTravelTimes(rtree,kValues[k],*queryNodeIt,kNNs,kNNDistances,queryNodeX,queryNodeY,myGraph.getMaxGraphSpeedByEdge());
                                        sw.stop();
                                    }                                    
                                }
                                totalQueryTime += sw.getTimeUs();
                                if (verifyKNN) {
                                    ineKNNs.clear();
                                    ineKNNDistances.clear();
                                    ine.getKNNs(myGraph,kValues[k],*queryNodeIt,ineKNNs,ineKNNDistances);
                                    // Note: We are marking Distance Browsing kNN as unsorted until equal max priority problem in algorithm is resolved
                                    if (!utility::verifyKNN(ineKNNs,ineKNNDistances,kNNs,kNNDistances,false,kValues[k],message,true)) {
                                        std::cout << "Verfication failed for IER by " << desc << " on object index " << objIdxFilePath << " for query node " << *queryNodeIt << " with k = " << kValues[k] << std::endl;
                                        std::cout << "Message: " << message << std::endl;
                                    }
                                }
                            }
                        }

                        double ierQueryTimeMs = totalQueryTime/totalQueries;

                        // Collect stats and return to output to file
                        kNNs.clear(); // Clear so that we don't pass last executed queries results to stats tuple
                        kNNDistances.clear();
                        KnnQueryTuple stats(myGraph.getNetworkName(),myGraph.getNumNodes(),myGraph.getNumEdges(),totalQueries,method,
                                            kValues[k],ierQueryTimeMs,objTypes[i],objDensities[j],objVariable[n],static_cast<int>(totalObjects/numSets),kNNs,kNNDistances);
                        stats.setAdditionalFields(specialFields);
                        stats.addSupplementaryFields("rtree_branch_factor",std::to_string(branchFactors[m]));
                        std::ofstream statsFS(statsOutputFile, std::ios::out | std::ios::app);
                        if (statsFS.is_open()) {
                            statsFS << stats.getTupleString() << std::endl;
                        } else {
                            std::cerr << "Cannot open stats file " << statsOutputFile << std::endl;
                        }                        
                    }
                }
            }
        }
    }
    std::cout << "IER by " << desc << " kNN queries successfully executed for " << myGraph.getNetworkName() << std::endl;
    if (verifyKNN) {
        std::cout << "IER by " << desc << " kNN verification completed for " << myGraph.getNetworkName() << std::endl;
    }
    delete graph;
    delete searchGraph;
}

void ShortestPathWrapper::runRealWorldIERQueries(Graph_RoadNetwork& myGraph, string bgrIntFilePath, string bcoIntFilePath,
                                                 string chFilePath, string tnrFilePath, std::vector<int>& branchFactors, 
                                                 std::vector<NodeID>& queryNodes, std::vector<int>& kValues, string filePathPrefix, 
                                                 std::string statsOutputFile, std::vector<string>& fileNames, std::vector<string>& setNames)
{
    datastr::graph::UpdateableGraph* graph;
    std::ifstream bgrFile(bgrIntFilePath, std::ios::in | std::ios::binary);
    if (bgrFile.is_open()) {
         graph = new datastr::graph::UpdateableGraph(bgrFile);
    } else {
        std::cerr << "Cannot open binary graph file " << bgrIntFilePath << std::endl;
        std::exit(1);
    }
    bgrFile.close();
    
    datastr::graph::SearchGraph* searchGraph;
    std::ifstream sgrFile(chFilePath, std::ios::in | std::ios::binary);
    if (sgrFile.is_open()) {
         searchGraph = new datastr::graph::SearchGraph(sgrFile);
    } else {
        std::cerr << "Cannot open search graph file " << chFilePath << std::endl;
        std::exit(1);
    }
    sgrFile.close();
    
    std::vector<CoordinateType> g_stXcord, g_stYcord;
    CoordinateLoader::deserializeCoordinate(bcoIntFilePath,g_stXcord,g_stYcord);

    TransitNodeTest tnr;
    processing::DijkstraCH<datastr::graph::SearchGraph, NormalPQueue, 2, true> dijkstraCHTest;
    

    std::ifstream tnrFile(tnrFilePath, std::ios::in | std::ios::binary);
    if (!tnrFile.is_open()) {
        std::cerr << "Cannot open TNR index file " << tnrFilePath << std::endl;
        std::exit(1);
    }
    tnr.loadIndex(searchGraph, graph, tnrFile);
    
    std::vector<NodeID> kNNs;
    std::vector<EdgeWeight> kNNDistances;    
    std::vector<std::string> parameterKeys, parameterValues;
    StopWatch sw;
    double totalQueryTime;
    int totalQueries = queryNodes.size();
    int totalObjects;
    Coordinate queryNodeX, queryNodeY;

    for (std::size_t i = 0; i < kValues.size(); ++i) {
        kNNs.reserve(kValues[i]);
        kNNDistances.reserve(kValues[i]);
        for (std::size_t j = 0; j < fileNames.size(); ++j) {
            std::string objIdxFilePath = filePathPrefix + "/real_world_pois/" + myGraph.getNetworkName() + "/" +fileNames[j] + "." + constants::OBJ_IDX_RTREE;
            StaticRtree rtree = serialization::getIndexFromBinaryFile<StaticRtree>(objIdxFilePath);
            
            if (myGraph.getEdgeType() == constants::DISTANCE_WEIGHTS) {
                sw.reset();
                sw.start();
                for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                    kNNs.clear();
                    kNNDistances.clear();
                    myGraph.getCoordinates(*queryNodeIt,queryNodeX,queryNodeY);
                    tnr.getKNNs(rtree,kValues[i],*queryNodeIt,kNNs,kNNDistances,queryNodeX,queryNodeY,g_stXcord,g_stYcord);
                }
                sw.stop();
                totalQueryTime = sw.getTimeUs();
            } else {
                sw.reset();
                sw.start();
                for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                    kNNs.clear();
                    kNNDistances.clear();
                    myGraph.getCoordinates(*queryNodeIt,queryNodeX,queryNodeY);
                    tnr.getKNNsByTravelTimes(rtree,kValues[i],*queryNodeIt,kNNs,kNNDistances,queryNodeX,queryNodeY,g_stXcord,g_stYcord,myGraph.getMaxGraphSpeedByEdge());
                }
                sw.stop();
                totalQueryTime = sw.getTimeUs();
            }
            
            double queryTimeMs = totalQueryTime/totalQueries;
            
            // Collect stats and return to output to file
            kNNs.clear(); // Clear so that we don't pass last executed queries results to stats tuple
            kNNDistances.clear();
            KnnQueryTuple stats(myGraph.getNetworkName(),myGraph.getNumNodes(),myGraph.getNumEdges(),totalQueries,constants::IER_TNR_KNN_QUERY,
                                kValues[i],queryTimeMs,setNames[j],rtree.getObjSetDensity(),rtree.getObjSetVariable(),rtree.getObjSetSize(),kNNs,kNNDistances);
            stats.addSupplementaryFields("rtree_branch_factor",std::to_string(rtree.getBranchFactor()));
            std::ofstream statsFS(statsOutputFile, std::ios::out | std::ios::app);
            if (statsFS.is_open()) {
                statsFS << stats.getTupleString() << std::endl;
            } else {
                std::cerr << "Cannot open stats file " << statsOutputFile << std::endl;
            }             
        }
    }

    dijkstraCHTest.loadGraph(searchGraph);    

    for (std::size_t i = 0; i < kValues.size(); ++i) {
        kNNs.reserve(kValues[i]);
        kNNDistances.reserve(kValues[i]);
        for (std::size_t j = 0; j < fileNames.size(); ++j) {
            std::string objIdxFilePath = filePathPrefix + "/real_world_pois/" + myGraph.getNetworkName() + "/" +fileNames[j] + "." + constants::OBJ_IDX_RTREE;
            StaticRtree rtree = serialization::getIndexFromBinaryFile<StaticRtree>(objIdxFilePath);
            
            if (myGraph.getEdgeType() == constants::DISTANCE_WEIGHTS) {
                sw.reset();
                sw.start();
                for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                    kNNs.clear();
                    kNNDistances.clear();
                    dijkstraCHTest.clear();
                    myGraph.getCoordinates(*queryNodeIt,queryNodeX,queryNodeY);
                    dijkstraCHTest.getKNNs(rtree,kValues[i],*queryNodeIt,kNNs,kNNDistances,queryNodeX,queryNodeY);
                }
                sw.stop();
                totalQueryTime = sw.getTimeUs();
            } else {
                sw.reset();
                sw.start();
                for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                    kNNs.clear();
                    kNNDistances.clear();
                    dijkstraCHTest.clear();
                    myGraph.getCoordinates(*queryNodeIt,queryNodeX,queryNodeY);
                    dijkstraCHTest.getKNNsByTravelTimes(rtree,kValues[i],*queryNodeIt,kNNs,kNNDistances,queryNodeX,queryNodeY,myGraph.getMaxGraphSpeedByEdge());
                }
                sw.stop();
                totalQueryTime = sw.getTimeUs();
            }
            
            double queryTimeMs = totalQueryTime/totalQueries;
            
            // Collect stats and return to output to file
            kNNs.clear(); // Clear so that we don't pass last executed queries results to stats tuple
            kNNDistances.clear();
            KnnQueryTuple stats(myGraph.getNetworkName(),myGraph.getNumNodes(),myGraph.getNumEdges(),totalQueries,constants::IER_CH_KNN_QUERY,
                                kValues[i],queryTimeMs,setNames[j],rtree.getObjSetDensity(),rtree.getObjSetVariable(),rtree.getObjSetSize(),kNNs,kNNDistances);
            stats.addSupplementaryFields("rtree_branch_factor",std::to_string(rtree.getBranchFactor()));
            std::ofstream statsFS(statsOutputFile, std::ios::out | std::ios::app);
            if (statsFS.is_open()) {
                statsFS << stats.getTupleString() << std::endl;
            } else {
                std::cerr << "Cannot open stats file " << statsOutputFile << std::endl;
            }             
        }
    }
    delete graph;
    delete searchGraph;
}
