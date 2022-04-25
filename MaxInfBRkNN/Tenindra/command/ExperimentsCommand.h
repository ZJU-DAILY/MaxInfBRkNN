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

#ifndef _EXPERIMENTSCOMMAND_H
#define _EXPERIMENTSCOMMAND_H

#include "Command.h"
#include "ExperimentsCommand.h"

#include "../processing/Graph.h"
#include "../processing/DynamicGraph.h"
//#include "../processing/Gtree.h"
//#include "../processing/ROAD.h"
//#include "../processing/MortonList.h"
//#include "../processing/INE.h"
//#include "../processing/IER.h"
#include "../processing/SetGenerator.h"
//#include "../processing/ALT.h"
//#include "../processing/ShortestPathWrapper.h"
#include "../tuple/IndexTuple.h"
#include "../tuple/ObjectIndexTuple.h"
#include "../tuple/KnnQueryTuple.h"
#include "../common.h"
#include "../utility/StopWatch.h"
#include "../utility/Statistics.h"
#include "../utility/utility.h"
#include "../utility/serialization.h"

#include <unordered_map>


class ExperimentsCommand: public Command {

public:
    void execute(int argc, char* argv[]){


    }


    void buildingPHL(){
        cout<<"begin buildingPHL..."<<endl;
        TEST_START
        std::string experiment = "";
        std::string bgrFilePath = "";


        std::string parameters = "";
        unsigned int numSets = 0;
        std::string method = "";
        unsigned int numPoints = 0;

        //路网图数据二进制文件目录
        stringstream binaryGraphFilePath;
        binaryGraphFilePath<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name<<".bin";
        bgrFilePath = binaryGraphFilePath.str();

        //索引构建参数
        parameters = "gtree=0,road=0,silc=0,phl=1,ch=0,tnr=0,alt=1,gtree_fanout=4,gtree_maxleafsize=64,road_fanout=4,road_levels=7,silc_maxquadtreeleafsize=1000,tnr_gridsize=128,alt_numlandmarks=16,alt_landmarktype=random";

        //路网索引构建处理相关日志信息输出目录
        stringstream statsPath;
        statsPath<<DATASET_HOME<<dataset_name<<"/road/"<<"stats/index_stats.txt";
        std::string statsOutputFile = statsPath.str();

        //索引相关文件输出根目录
        stringstream roadDataPath;
        roadDataPath<<DATASET_HOME<<dataset_name<<"/";
        std::string filePathPrefix  = roadDataPath.str();

        this -> buildPhlIndexes(bgrFilePath,parameters,filePathPrefix,statsOutputFile);
        TEST_END
        //TEST_DURA_PRINT("************PHL build")


    }

    void showCommandUsage(std::string programName){
        std::cerr << "Usage: " << programName << " -c " + constants::EXPERIMENTS_CMD + " -e <experimental step>\n\n"
                  << "Steps:\n"
                  << utility::getFormattedUsageString(constants::EXP_BUILD_INDEXES,"1. Build all indexes") + "\n"
                  << utility::getFormattedUsageString(constants::EXP_BUILD_OBJ_INDEXES,"2. Generate object sets and build object indexes") + "\n"
                  << utility::getFormattedUsageString(constants::EXP_RUN_KNN,"3. Run all standard kNN query experiments") + "\n"
                  << utility::getFormattedUsageString(constants::EXP_RUN_KNN_OPTIMIZATIONS,"4. Run individual method kNN query experiments") + "\n"
                  << utility::getFormattedUsageString(constants::EXP_RUN_KNN_RW_POI,"5. Build object indexes for real-world POIs") + "\n"
                  << utility::getFormattedUsageString(constants::EXP_RUN_KNN_RW_POI,"6. Run kNN querying on real-world POIs experiments") + "\n";
    }
    void showPhaseUsage(std::string method, std::string programName){
        if (method == constants::EXP_BUILD_INDEXES ) {
            std::cerr << "Usage: " << programName << " -c " + constants::EXPERIMENTS_CMD
                      << " -e " + constants::EXP_BUILD_INDEXES + " -g <binary graph file>\n"
                      << "-p <parameter key-value pairs> -f <index output file path prefix> -s <stats output file>\n";
        } else if (method == constants::EXP_BUILD_OBJ_INDEXES ) {
            std::cerr << "Usage: " << programName << " -c " + constants::EXPERIMENTS_CMD
                      << " -e " + constants::EXP_BUILD_OBJ_INDEXES + " -g <binary graph file>\n"
                      << "-p <parameter key-value pairs> -n <num sets> -d <list of object set densities or partitions>\n"
                      << "-t <list of object types>  -v <list of some object variable> -f <index output file path prefix> -s <stats output file>\n";
        } else if (method == constants::EXP_RUN_KNN ) {
            std::cerr << "Usage: " << programName << " -c " + constants::EXPERIMENTS_CMD
                      << " -e " + constants::EXP_RUN_KNN + " -g <binary graph file>\n"
                      << "-q <query node file> -k <k values> -p <parameters> -n <num sets> -d <list of object set densities or partitions>\n"
                      << "-t <list of object types>  -v <list of some object variable> -f <index output file path prefix> -s <stats output file>\n";
        } else if (method == constants::EXP_RUN_KNN_OPTIMIZATIONS ) {
            std::cerr << "Usage: " << programName << " -c " + constants::EXPERIMENTS_CMD
                      << " -e " + constants::EXP_RUN_KNN_OPTIMIZATIONS + " -m <method>\n"
                      << "-g <binary graph file> -q <query node file> -k <k values> -p <parameters>\n"
                      << "-n <num sets> -d <list of object densities> -t <list of object types>\n"
                      << " -v <list of some object variable> -f <index output file path prefix> -s <stats output file>\n";
        } else if (method == constants::EXP_RUN_KNN_RW_POI ) {
            std::cerr << "Usage: " << programName << " -c " + constants::EXPERIMENTS_CMD
                      << " -e " + constants::EXP_RUN_KNN_RW_POI + " -g <binary graph file>\n"
                      << "-q <query node file> -k <k values> -p <parameters> -r <rw POI set list file>\n"
                      << "-f <index output file path prefix> -s <stats output file>\n";
        } else if (method == constants::EXP_BUILD_RW_POI_OBJ_INDEXES ) {
            std::cerr << "Usage: " << programName << " -c " + constants::EXPERIMENTS_CMD
                      << " -e " + constants::EXP_RUN_KNN_RW_POI + " -g <binary graph file>\n"
                      << "-p <parameters> -r <rw POI set list file>\n"
                      << "-f <index output file path prefix> -s <stats output file>\n";
        } else {
            std::cerr << "Invalid experiment phase!" << std::endl;
            this->showCommandUsage(programName);
        }
    }

private:
    //void buildIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string statsOutputFile);
    void buildExternalIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string statsOutputFile){
        std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);

        bool buildPHLIdx = parameterMap["phl"] == "1";
        bool buildCHIdx = parameterMap["ch"] == "0";
        bool buildTNRIdx = parameterMap["tnr"] == "0";

        // Find all additional fields we need to add to index stats tuples
        std::vector<std::string> specialFields;
        int field = 0;
        while (true) {
            std::string key = "special_field_" + std::to_string(field);
            if (parameterMap.find(key) != parameterMap.end()) {
                specialFields.push_back(parameterMap[key]);
                ++field;
            } else {
                break;
            }
        }

        // Put all of these into functions (except input validation) so indexes get released after serialization
        std::string idxOutputFile = "";
        std::vector<std::string> parameterKeys, parameterValues;

        // Some additional file paths to create intermediary and additional files for some indexes
        std::string graphDataOutputFile = "", coordDataOutputFile = "", bgrIntFile = "", bcoIntFile = "", chFile = "";

        std::string networkName = "";
        int numNodes = 0, numEdges = 0;
        if (buildCHIdx || buildPHLIdx || buildTNRIdx) {
            // Use scope to destroy graph data structure to have extra memory for methods that don't need it
            Graph_RoadNetwork tempGraph = serialization::getIndexFromBinaryFile<Graph_RoadNetwork>(bgrFileName);
            networkName = tempGraph.getNetworkName();
            numNodes = tempGraph.getNumNodes();
            numEdges = tempGraph.getNumEdges();

            if (buildPHLIdx) {
                parameterKeys.clear();
                parameterValues.clear();
                graphDataOutputFile = filePathPrefix + "/data/" + utility::constructIndexFileName(networkName,"tsv",parameterKeys,parameterValues);
                tempGraph.outputToTSVFile(graphDataOutputFile);
            }
        }

        if (buildPHLIdx) {
            parameterKeys.clear();
            parameterValues.clear();
            graphDataOutputFile = filePathPrefix + "/data/" + utility::constructIndexFileName(networkName,"tsv",parameterKeys,parameterValues);
            idxOutputFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(networkName,constants::IDX_PHL_CMD,parameterKeys,parameterValues);
            this->buildPHL(networkName,numNodes,numEdges,graphDataOutputFile,idxOutputFile,statsOutputFile,specialFields);

            // Delete intermediary text file as it is not needed
            std::remove(graphDataOutputFile.c_str());
        }


        std::cout << "External index building complete!" << std::endl;
    }
    void buildPhlIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string statsOutputFile){

        std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);

        bool buildPHLIdx = parameterMap["phl"] == "1";
        bool buildCHIdx = parameterMap["ch"] == "0";
        bool buildTNRIdx = parameterMap["tnr"] == "0";

        // Find all additional fields we need to add to index stats tuples
        std::vector<std::string> specialFields;
        int field = 0;
        while (true) {
            std::string key = "special_field_" + std::to_string(field);
            if (parameterMap.find(key) != parameterMap.end()) {
                specialFields.push_back(parameterMap[key]);
                ++field;
            } else {
                break;
            }
        }

        // Put all of these into functions (except input validation) so indexes get released after serialization
        std::string idxOutputFile = "";
        std::vector<std::string> parameterKeys, parameterValues;

        // Some additional file paths to create intermediary and additional files for some indexes
        std::string graphDataOutputFile = "", coordDataOutputFile = "", bgrIntFile = "", bcoIntFile = "", chFile = "";

        std::string networkName = "";
        int numNodes = 0, numEdges = 0;
        if (buildCHIdx || buildPHLIdx || buildTNRIdx) {
            // Use scope to destroy graph data structure to have extra memory for methods that don't need it
            Graph_RoadNetwork tempGraph = serialization::getIndexFromBinaryFile<Graph_RoadNetwork>(bgrFileName);
            networkName = tempGraph.getNetworkName();
            numNodes = tempGraph.getNumNodes();
            numEdges = tempGraph.getNumEdges();

            if (buildPHLIdx) {
                parameterKeys.clear();
                parameterValues.clear();
                graphDataOutputFile = filePathPrefix + "road/" + utility::constructIndexFileName(networkName,"tsv",parameterKeys,parameterValues);
                cout<<graphDataOutputFile<<endl;
                tempGraph.outputToTSVFile(graphDataOutputFile);
            }
        }

        if (buildPHLIdx) {
            parameterKeys.clear();
            parameterValues.clear();
            graphDataOutputFile = filePathPrefix + "road/" + utility::constructIndexFileName(networkName,"tsv",parameterKeys,parameterValues);
            idxOutputFile = filePathPrefix + "road/" + utility::constructIndexFileName(networkName,constants::IDX_PHL_CMD,parameterKeys,parameterValues);
            cout<<"graphDataOutputFile="<<graphDataOutputFile<<endl;
            cout<<"idxOutputFile="<<idxOutputFile<<endl;
            this->buildPHL(networkName,numNodes,numEdges,graphDataOutputFile,idxOutputFile,statsOutputFile,specialFields);

            // Delete intermediary text file as it is not needed
            std::remove(graphDataOutputFile.c_str());
        }


        std::cout << "Phl index building complete!" << std::endl;
    }



    std::unordered_map<std::string,std::string> getParameters(std::string parameters)
    {
        std::unordered_map<std::string,std::string> parameterMap;
        std::vector<std::string> pairs = utility::splitByDelim(parameters,',');

        for (std::size_t i = 0; i < pairs.size(); ++i) {
            std::vector<std::string> pair = utility::splitByDelim(pairs[i],'=');
            if (pair.size() == 2) {
//             std::cout << "Key = " << pair[0] << std::endl;
//             std::cout << "Value = " << pair[1] << std::endl;
                parameterMap[pair[0]] = pair[1];
            } else {
                std::cerr << "Invalid key-value pair in parameter string" << std::endl;
                exit(1);
            }
        }
        return parameterMap;
    }

    void buildPHL(std::string networkName, int numNodes, int numEdges, std::string dataOutputFile, std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields)
    {
        PrunedHighwayLabeling phl;
        phl.ConstructLabel(dataOutputFile.c_str());

        double processingTimeMs = phl.getConstructionTime()*1000; // convert to ms
        double memoryUsage = phl.computeIndexSize();

        IndexTuple stats(networkName,numNodes,numEdges,constants::IDX_PHL_CMD,processingTimeMs,memoryUsage);
        stats.setAdditionalFields(specialFields);
        this->outputCommandStats(statsOutputFile,stats.getTupleString());

        phl.StoreLabel(idxOutputFile.c_str());

        std::cout << "PHL index successfully created!" << std::endl;
    }





    /*void runSingleMethodQueries(std::string method, std::string bgrFileName,
                                                std::string queryNodeFile, std::string kValues, std::string parameters,
                                                std::size_t numSets, std::string objDensities, std::string objTypes, std::string objVariable,
                                                std::string filePathPrefix, std::string statsOutputFile)
    {
        Graph_RoadNetwork graph = serialization::getIndexFromBinaryFile<Graph_RoadNetwork>(bgrFileName);
        std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);

        std::vector<NodeID> queryNodes = utility::getPointSetFromFile(queryNodeFile);

        bool verifykNN = parameterMap["verify"] == "1";

        // Find all additional fields we need to add to kNN stats tuples
        std::vector<std::string> specialFields;
        int field = 0;
        while (true) {
            std::string key = "special_field_" + std::to_string(field);
            if (parameterMap.find(key) != parameterMap.end()) {
                specialFields.push_back(parameterMap[key]);
                ++field;
            } else {
                break;
            }
        }

        std::vector<std::string> strObjDensitiesVec = utility::splitByDelim(objDensities,',');
        std::vector<double> objDensitiesVec;
        for(std::size_t i = 0; i < strObjDensitiesVec.size(); ++i) {
            double density = std::stod(strObjDensitiesVec[i]);
            if (density > 0) {
                objDensitiesVec.push_back(density);
            } else {
                std::cerr << "Invalid density in list provided!\n";
                exit(1);
            }
        }
        std::vector<std::string> objTypesVec = utility::splitByDelim(objTypes,',');
        std::vector<std::string> strKValuesVec = utility::splitByDelim(kValues,',');
        std::vector<int> kValuesVec;
        for(std::size_t i = 0; i < strKValuesVec.size(); ++i) {
            int k = std::stoi(strKValuesVec[i]);
            if (k > 0) {
                kValuesVec.push_back(k);
            } else {
                std::cerr << "Invalid k value in list provided!\n";
                exit(1);
            }
        }
        std::vector<std::string> strObjVariables = utility::splitByDelim(objVariable,',');
        std::vector<int> objVariableVec;
        for(std::size_t i = 0; i < strObjVariables.size(); ++i) {
            int variable = std::stoi(strObjVariables[i]);
            if (variable > 0) {
                objVariableVec.push_back(variable);
            } else {
                std::cerr << "Invalid variable in list provided (must be greater than zero)!\n";
                exit(1);
            }
        }
        if (objDensitiesVec.size() == 0 || objTypesVec.size() == 0 || objVariableVec.size() == 0 || kValuesVec.size() == 0) {
            std::cerr << "Not enough densities or types provided!\n";
            exit(1);
        }

        std::string gtreeIdxFile, roadIdxFile, silcIdxFile, juncIdxFile, phlIdxFile, altIdxFile;
        std::vector<std::string> parameterKeys, parameterValues;

        if (method == constants::IER_PHL_KNN_QUERY) {
            std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
            if (branchFactors.size() == 0) {
                std::cerr << "No R-tree branch factors provided " << method << "!\n";
                exit(1);
            }
            parameterKeys.clear();
            parameterValues.clear();
            phlIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_PHL_CMD,parameterKeys,parameterValues);
            if (utility::fileExists(phlIdxFile)) {
                this->runIERQueries(graph,constants::PHL_SPDIST_QUERY,phlIdxFile,branchFactors,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,specialFields);
            } else {
                std::cerr << "PHL index not built for " << graph.getNetworkName() << ", unable to execute IER-PHL\n";
            }
        }else {
            std::cerr << "Could not recognise method - check kNN query command" << std::endl;
            exit(1);
        }

    }*/


    void runDistanceQueriesByPHL(std::string method, std::string bgrFileName,
                                 std::string filePathPrefix, std::string statsOutputFile)
    {
        Graph_RoadNetwork graph = serialization::getIndexFromBinaryFile<Graph_RoadNetwork>(bgrFileName);


        if (true) {


            std::string phlIdxFile = filePathPrefix + "/indexes/" + graph.getNetworkName() + "." + constants::IDX_PHL_CMD;
            //utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_PHL_CMD,parameterKeys,parameterValues);


            if (utility::fileExists(phlIdxFile)) {
                //this->runIERQueries(graph,constants::PHL_SPDIST_QUERY,phlIdxFile,branchFactors,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,specialFields);
                PrunedHighwayLabeling phl;
                //phl.LoadLabel(phlIdxFile)
                phl.LoadLabel(phlIdxFile.c_str());
                int dist = phl.Query(10, 1588);
                std::cout<<"SPDist(10, 1588)="<<dist<<std::endl;
            } else {
                std::cerr << "PHL index not built for " << graph.getNetworkName() << ", unable to execute IER-PHL\n";
            }
        }else {
            std::cerr << "Could not recognise method - check kNN query command" << std::endl;
            exit(1);
        }

    }

};

#endif // _EXPERIMENTSCOMMAND_H