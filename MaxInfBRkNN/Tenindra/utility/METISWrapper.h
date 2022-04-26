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

#ifndef _METISWRAPPER_H
#define _METISWRAPPER_H

#include "../processing/DynamicGraph.h"
#include "../processing/Graph.h"

#include <metis.h>
#include <unordered_set>
#include <vector>

class METISWrapper {

    public:
        METISWrapper(int maxNumNodes, int maxNumEdges, int numParts);
        ~METISWrapper();
        void initializeMETISOptions();
        void setNumCuts(idx_t ncuts);
        void setNumIterations(idx_t niter);
        idx_t partitionSubgraph(idx_t& nvtxs);
        idx_t partitionSubgraphNonContiguous(idx_t& nvtxs);
        void populateMETISArrays(DynamicGraph& graph, std::unordered_set<NodeID>& subgraphNodes, std::vector<NodeID>& METISIdxToNodeID, bool setEdgeWeightsToOne = false);
        void populateMETISArrays(Graph_RoadNetwork& graph, std::unordered_set<NodeID>& subgraphNodes, std::vector<NodeID>& METISIdxToNodeID, bool setEdgeWeightsToOne = false);
        void updateNumPartitions(int numParts);
        idx_t *parts;
         
    private:
        idx_t options[METIS_NOPTIONS];
        idx_t ncons;
        idx_t nparts;
        idx_t *xadj;
        idx_t *adjncy;
        idx_t *adjwgt;
        int maxNumNodes;
        int maxNumEdges;
};

#endif // _METISWRAPPER_H