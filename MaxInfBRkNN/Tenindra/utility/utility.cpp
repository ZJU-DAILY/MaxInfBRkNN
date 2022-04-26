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

#include "utility.h"

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <assert.h>
#include <unordered_set>
#include <algorithm>

namespace utility
{
    // kNN1 is the expected results, and kNN2 is the results being tested
    bool verifyKNN(std::vector<NodeID>& kNN1, std::vector<EdgeWeight>& kNNDistances1, 
                   std::vector<NodeID>& kNN2, std::vector<EdgeWeight>& kNNDistances2,
                   bool unsorted, unsigned int k, std::string& message, bool verifyDistances)
    {
        if (kNN1.size() != kNNDistances1.size()) {
            message = "Baseline number of kNNs does not match number of distances";
            return false;
        }
        
        if (kNN1.size() < kNN2.size()) {
            message = "Too many kNN results";
            return false;
        }
        
        if (kNN1.size() > kNN2.size()) {
            message = "Too few kNN results";
            return false;
        }

        if (kNN2.size() != kNNDistances2.size()) {
            message = "Number of kNNs does not match number of distances";
            return false;
        }
        
        std::vector<NodeDistancePair> kNN1Combined = utility::combineNodesAndDistanceVectors(kNN1,kNNDistances1);
        std::vector<NodeDistancePair> kNN2Combined = utility::combineNodesAndDistanceVectors(kNN2,kNNDistances2);

        if (unsorted) {
            std::sort(kNN1Combined.begin(),kNN1Combined.end(), utility::compareNodeDistancePair);
            std::sort(kNN2Combined.begin(),kNN2Combined.end(), utility::compareNodeDistancePair);
        }
        
        for (std::size_t i = 0; i < kNN1Combined.size(); ++i) {
            if (kNN1Combined[i].first != kNN2Combined[i].first) {
                if (kNN1Combined[i].second != kNN2Combined[i].second) {
                    message = "Incorrect kNN found: " + std::to_string(kNN2Combined[i].first) + " should be " + std::to_string(kNN1Combined[i].first);
                    return false;
                }
            } else if (verifyDistances) {
                if (kNN1Combined[i].second != kNN2Combined[i].second) {
                    message = "Incorrect kNN distance found for " + std::to_string(kNN2Combined[i].first) + ": distance " + 
                        std::to_string(kNN2Combined[i].second) + " should be " + std::to_string(kNN1Combined[i].second);
                    return false;
                }
            }
        }
        
        return true;
    }
    
    std::vector<NodeDistancePair> combineNodesAndDistanceVectors(std::vector<NodeID>& nodes, std::vector<EdgeWeight>& distances)
    {
        std::vector<NodeDistancePair> nodesAndDistances;
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            nodesAndDistances.push_back(NodeDistancePair(nodes[i],distances[i]));
        }
        return nodesAndDistances;
    }
    
    bool compareNodeDistancePair(const NodeDistancePair& pair1, const NodeDistancePair& pair2)
    {
        if (pair1.second < pair2.second) {
            return true;
        } else if (pair1.second == pair2.second) {
            // If equal use the first to break tie
            if (pair1.first < pair2.first) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }
    
    std::string getFormattedUsageString(std::string commandName, std::string description) {

        // <--leading white space--><--name--><--mid white space--><--description-->
        // <--2 chars--><--15 chars max--><--mid white space-->

        assert(commandName.length()+1 < constants::MAX_COMMAND_NAME_SIZE && "Command name is too long to format");
        
        int leadingSpaceSize = 2;
        int middleSpaceSize = leadingSpaceSize + constants::MAX_COMMAND_NAME_SIZE - commandName.size();
        
        std::string leadingSpace(leadingSpaceSize,' '); // "  "
        std::string middleSpace(middleSpaceSize,' '); // "  "
        
        return leadingSpace + commandName + middleSpace + description;
    }
    
    std::string implodeNodesAndDistances(std::vector<NodeID> kNNs, std::vector<NodeID> kNNDistances)
    {
        std::string implodedString = "", pair = "";
        for (std::size_t i = 0; i < kNNs.size(); ++i) {
            if (i != 0) {
                implodedString += ",";
            }
            pair = "[" + std::to_string(kNNs[i]) + "," + std::to_string(kNNDistances[i]) + "]";
            implodedString += pair;
        }
        return implodedString;
    }
    
    std::vector<NodeID> getPointSetFromFile(std::string queryPointsFilePath)
    {
        std::vector<NodeID> queryNodesIDs;
        NodeID queryNodeID;
        std::string networkName, setType; 
        double setDensity; 
        int setSize, setVariable;

        std::ifstream queryPointFS(queryPointsFilePath, std::ios::in);
        if (queryPointFS.is_open()) {
            // First line is query point set information
            queryPointFS >> networkName >> setType >> setDensity >> setVariable >> setSize;

            while (queryPointFS >> queryNodeID)  {
                queryNodesIDs.push_back(queryNodeID);
            }
        } else {
            std::cerr << "Cannot open point set file " << queryPointsFilePath << std::endl;
            exit(1);
        }
        queryPointFS.close();
        return queryNodesIDs;
    }
    
    std::vector<NodeID> getPointSetFromFile(std::string queryPointsFilePath, std::string& setType, double& setDensity, int& setSize, int& setVariable)
    {
        std::vector<NodeID> queryNodesIDs;
        NodeID queryNodeID;
        std::string networkName; 

        std::ifstream queryPointFS(queryPointsFilePath, std::ios::in);
        if (queryPointFS.is_open()) {
            // First line is query point set information
            queryPointFS >> networkName >> setType >> setDensity >> setVariable >> setSize;

            while (queryPointFS >> queryNodeID)  {
                queryNodesIDs.push_back(queryNodeID);
            }
        } else {
            std::cerr << "Cannot open point set file " << queryPointsFilePath << std::endl;
            exit(1);
        }
        queryPointFS.close();
        return queryNodesIDs;
    }
    
    std::unordered_set<NodeID> getPointSetFromFileUset(std::string queryPointsFilePath, std::string& setType, double& setDensity, int& setSize, int& setVariable)
    {
        std::unordered_set<NodeID> queryNodesIDs;
        NodeID queryNodeID;
        std::string networkName; 

        std::ifstream queryPointFS(queryPointsFilePath, std::ios::in);
        if (queryPointFS.is_open()) {
            // First line is query point set information
            queryPointFS >> networkName >> setType >> setDensity >> setVariable >> setSize;

            while (queryPointFS >> queryNodeID)  {
                queryNodesIDs.insert(queryNodeID);
            }
        } else {
            std::cerr << "Cannot open point set file " << queryPointsFilePath << std::endl;
            exit(1);
        }
        queryPointFS.close();
        return queryNodesIDs;
    }
    
    std::vector<bool> populateGraphObjectStatusFromFile(int graphNumNodes, std::string queryPointsFilePath, std::string& setType, double& setDensity, int& setSize)
    {
        std::vector<bool> queryNodesIDs(graphNumNodes,false);
        NodeID queryNodeID;
        std::string networkName; 

        std::ifstream queryPointFS(queryPointsFilePath, std::ios::in);
        if (queryPointFS.is_open()) {
            // First line is query point set information
            queryPointFS >> networkName >> setType >> setDensity >> setSize;

            while (queryPointFS >> queryNodeID)  {
                queryNodesIDs[queryNodeID] = true;
            }
        } else {
            std::cerr << "Cannot open point set file " << queryPointsFilePath << std::endl;
            exit(1);
        }
        queryPointFS.close();
        return queryNodesIDs;
    }
    
    std::vector<NodePair> getPointPairSetFromFile(std::string queryPointsFilePath, std::string& setType, unsigned int& numPairs)
    {
        std::vector<NodePair> nodePairs;
        NodeID sourceID, targetID;
        std::string networkName, uselessField; 

        std::ifstream queryPointFS(queryPointsFilePath, std::ios::in);
        if (queryPointFS.is_open()) {
            // First line is query point set information
            queryPointFS >> networkName >> setType >> numPairs;

            while (queryPointFS >> sourceID >> targetID)  {
                nodePairs.push_back(std::make_pair(sourceID,targetID));
            }
        } else {
            std::cerr << "Cannot open query point file " << queryPointsFilePath << std::endl;
            exit(1);
        }
        queryPointFS.close();
        return nodePairs;
    }
    
    std::string constructIndexFileName(std::string networkName, std::string indexType, std::vector<std::string>& parameterNames, std::vector<std::string>& parameterValues)
    {
        std::stringstream fileName;
        fileName << networkName;
        for (std::size_t i = 0; i < parameterNames.size(); ++i) {
            fileName << "_" << parameterNames[i] << "=" << parameterValues[i];
        }
        fileName << "." << indexType; // Index type is the extension
        return fileName.str();
    }

    std::string constructObjectIndexFileName(std::string networkName, std::string indexType, std::string setType, double setDensity, int setVariable, int setIdx, std::vector<std::string>& parameterNames, std::vector<std::string>& parameterValues)
    {
        std::stringstream fileName;
        fileName << networkName << "_" ;
        fileName << setType << "_" ;
        fileName << setDensity << "_" ;
        fileName << setVariable << "_" ;
        if (setDensity == 1 && setType == constants::RAND_OBJ_SET) {
            // If density is one we only need to create one object index
            fileName << std::setfill('0') << std::setw(3) << 0;
        } else {
            fileName << std::setfill('0') << std::setw(3) << setIdx;
        }
        for (std::size_t i = 0; i < parameterNames.size(); ++i) {
            fileName << "_" << parameterNames[i] << "=" << parameterValues[i];
        }
        fileName << "." << indexType; // Index type is the extension
        return fileName.str();
    }

    std::string constructObjsectSetFileName(std::string networkName, std::string setType, double setDensity, int setVariable, int setIdx)
    {
        std::stringstream fileName;
        fileName << networkName << "_" ;
        fileName << setType << "_" ;
        fileName << setDensity << "_" ;
        fileName << setVariable << "_" ;
        if (setDensity == 1 && setType == constants::RAND_OBJ_SET) {
            // If density is one we only need to create one object set
            fileName << std::setfill('0') << std::setw(3) << 0;
        } else {
            fileName << std::setfill('0') << std::setw(3) << setIdx;
        }
        fileName << ".txt";
        return fileName.str();
    }

    std::string constructQueryPointSetFileName(std::string networkName, std::string setType, std::string setSize, int setIdx)
    {
        std::stringstream fileName;
        fileName << networkName << "_" ;
        fileName << setType << "_" ;
        fileName << setSize << "_" ;
        fileName << std::setfill('0') << std::setw(3) << setIdx;
        fileName << ".txt";
        return fileName.str();
    }

    std::vector<std::string> splitByDelim(std::string input, char delimeter)
    {
        std::vector<std::string> columns;
        std::stringstream ss(input);
        std::string column;
        while (std::getline(ss, column, delimeter)) {
            columns.push_back(column);
        }
        return columns;
    }
    
    std::vector<int> getIngetersFromStringList(std::string input, char delimeter, int minValue)
    {
        std::vector<std::string> strValues = utility::splitByDelim(input,delimeter);
        std::vector<int> intValues;
        for(std::size_t i = 0; i < strValues.size(); ++i) {
            int value = std::stoi(strValues[i]);
            if (value > minValue) {
                intValues.push_back(value);
            } else {
                std::cerr << "Value in string is less than minimum value of " << minValue << std::endl;
                exit(1);    
            }
        }
        return intValues;
    }

    void writeSampleSet(std::string outputFile, std::string networkName, std::string setType, 
                                                double setDensity, int setSize, int setVariable, std::vector<NodeID>& sampleSet)
    {
        std::ofstream outputFS(outputFile, std::ios::out | std::ios::trunc);
        if (outputFS.is_open()) {
            // Write header line with object set details
            outputFS << networkName << " " << setType << " " 
                << setDensity << " " << setVariable << " " << setSize << std::endl;
            
            // Write one object NodeID per line
            for (auto it = sampleSet.begin(); it != sampleSet.end(); ++it) {
                outputFS << *it << std::endl;
            }
        } else {
            std::cerr << "Cannot open output file " << outputFile << std::endl;
        }
        outputFS.close();
    }

    void writePairSet(std::string outputFile, std::string networkName, std::string setType, 
                                                unsigned int setSize, std::vector<NodePair>& pairSet)
    {
        std::ofstream outputFS(outputFile, std::ios::out | std::ios::trunc);
        if (outputFS.is_open()) {
            // Write header line with apir set details
            outputFS << networkName << " " << setType << " " << setSize << std::endl;
            
            // Write source ID followed by target ID per line
            for (auto it = pairSet.begin(); it != pairSet.end(); ++it) {
                outputFS << it->first << " " << it->second << std::endl;
            }
        } else {
            std::cerr << "Cannot open output file " << outputFile << std::endl;
        }
        outputFS.close();
    }
    
    void flushCache()
    {
        int size = 2000000; // 16MB in total to flush cache
        std::vector<double> bigVector(size,1);
        volatile double sum = 0;
        for (double element: bigVector) {
            sum += element;
        }
    }
    
    double estimateUnorderedMapMemoryUsageBytes(int numElements, int elementSize, int numBuckets)
    {
        double memoryUsage = 0;
        memoryUsage += numElements*elementSize;
        memoryUsage += numBuckets*((sizeof(std::size_t)+sizeof(void*)));
        memoryUsage += 0.5*numElements*(sizeof(void*)+elementSize);
        return memoryUsage;
    }

    bool fileExists(const std::string& fileName)
    {
        std::ifstream file(fileName);
        if (file.good()) {
            file.close();
            return true;
        } else {
            file.close();
            return false;
        }
    }

}