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

#ifndef _UTILITY_H
#define _UTILITY_H

#include "../common.h"

#include <vector>
#include <unordered_set>

namespace utility
{
    std::string getFormattedUsageString(std::string commandName, std::string description);
    std::string implodeNodesAndDistances(std::vector<NodeID> kNNs,std::vector<NodeID> kNNDistances);
    std::vector<NodeID> getPointSetFromFile(std::string queryPointsFilePath);
    std::vector<NodeID> getPointSetFromFile(std::string queryPointsFilePath, std::string& setType, double& setDensity, int& setSize, int& setVariable);
    std::unordered_set<NodeID> getPointSetFromFileUset(std::string queryPointsFilePath, std::string& setType, double& setDensity, int& setSize, int& setVariable);
    std::vector<bool> populateGraphObjectStatusFromFile(int graphNumNodes, std::string queryPointsFilePath, std::string& setType, double& setDensity, int& setSize);
    std::vector<NodePair> getPointPairSetFromFile(std::string queryPointsFilePath, std::string& setType, unsigned int& numPairs);
    std::string constructIndexFileName(std::string networkName, std::string indexType, std::vector<std::string>& parameterNames, std::vector<std::string>& parameterValues);
    std::string constructObjectIndexFileName(std::string networkName, std::string indexType, std::string setType, 
                                             double setDensity, int setVariable, int setIdx, std::vector<std::string>& parameterNames, 
                                             std::vector<std::string>& parameterValues);
    std::string constructObjsectSetFileName(std::string networkName, std::string setType, double setDensity, int setVariable, int setIdx);
    std::string constructQueryPointSetFileName(std::string networkName, std::string setType, std::string setSize, int setIdx);
    std::vector<std::string> splitByDelim(std::string input, char delimeter);
    std::vector<int> getIngetersFromStringList(std::string input, char delimeter, int minValue);
    void writeSampleSet(std::string outputFile, std::string networkName, std::string setType, 
                        double setDensity, int setSize, int setVariable, std::vector<NodeID>& sampleSet);
    void writePairSet(std::string outputFile, std::string networkName, std::string setType, 
                        unsigned int setSize, std::vector<NodePair>& pairSet);
    bool verifyKNN(std::vector<NodeID>& kNN1, std::vector<EdgeWeight>& kNNDistances1, 
                   std::vector<NodeID>& kNN2, std::vector<EdgeWeight>& kNNDistances2,
                   bool unsorted, unsigned int k, std::string& message, bool verifyDistances = false);
    bool compareNodeDistancePair(const NodeDistancePair& pair1, const NodeDistancePair& pair2);
    std::vector<NodeDistancePair> combineNodesAndDistanceVectors(std::vector<NodeID>& nodes, std::vector<EdgeWeight>& distances);
    void flushCache();
    double estimateUnorderedMapMemoryUsageBytes(int numElements, int elementSize, int numBuckets);
    bool fileExists(const std::string& fileName);
    
    template <typename STLCollection>
    void releaseSTLCollection(STLCollection& collection) {
        STLCollection tmp;
        collection.swap(tmp);
    }
}

#endif // _UTILITY_H