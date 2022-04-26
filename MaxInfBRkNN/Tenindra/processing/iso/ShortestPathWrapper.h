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

#ifndef _SHORTEST_PATH_WRAPPER
#define _SHORTEST_PATH_WRAPPER

#include <string>
#include <vector>

#include "Graph.h"
#include "../common.h"

class ShortestPathWrapper {

    public:
        void buildIntermediaryGraph(std::string graphFilePath, std::string coordinateFilePath, std::string bgrIntFilePath, std::string bcoIntFilePath);
        void buildCH(std::string bgrIntFilePath, std::string chFilePath, std::string noFilePath);
        void buildTNR(std::string bgrIntFilePath, std::string bcoIntFilePath, std::string chFilePath, std::string tnrFilePath, int ascale);
        double getConstructionTimeMs();
        double getIndexSizeMB();
        void runIERQueries(Graph_RoadNetwork& graph, std::string method, std::string bgrIntFilePath, std::string bcoIntFilePath, std::string chFilePath, std::string tnrFilePath,
                           std::vector<int>& branchFactors, std::vector<NodeID>& queryNodes, std::vector<int>& kValues, std::size_t numSets,
                           std::vector<double> objDensities, std::vector<std::string> objTypes, std::vector<int> objVariable,
                           std::string filePathPrefix, std::string statsOutputFile, bool verifyKNN, std::vector<std::string> specialFields = {});
        void runRealWorldIERQueries(Graph_RoadNetwork& graph, std::string bgrIntFilePath, std::string bcoIntFilePath, std::string chFilePath, std::string tnrFilePath,
                           std::vector<int>& branchFactors, std::vector<NodeID>& queryNodes, std::vector<int>& kValues, 
                           std::string filePathPrefix, std::string statsOutputFile, std::vector<std::string>& fileNames, std::vector<std::string>& setNames);
    private:
        double getFileSizeMB(std::string filePath);
        
        double constructionTimeMs = 0;
        double memoryUsageMB = 0;
        
};


#endif // _SHORTEST_PATH_WRAPPER