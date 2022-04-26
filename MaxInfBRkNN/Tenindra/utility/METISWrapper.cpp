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

#include "METISWrapper.h"

#include <unordered_map>

METISWrapper::METISWrapper(int maxNumNodes, int maxNumEdges, int numParts): 
    maxNumNodes(maxNumNodes), maxNumEdges(maxNumEdges)
{
    this->initializeMETISOptions();
    this->nparts = static_cast<idx_t>(numParts);
}

METISWrapper::~METISWrapper()
{
    delete[] xadj;
    delete[] parts;
    delete[] adjncy;
    delete[] adjwgt;    
}

void METISWrapper::initializeMETISOptions()
{
    // Configure METIS Partitioning Options
    METIS_SetDefaultOptions(this->options);
    this->options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY; // Partition graph k-ways (this will be contiguous)
    this->options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL; // Minimise actual number of borders (slower than cut)
//     this->options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // Minimise edge cut (only estimated of total communication)
    this->options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
    this->options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_RANDOM;
    this->options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM;
    this->options[METIS_OPTION_UFACTOR] = 500;
    this->options[METIS_OPTION_CONTIG] = 1;
    this->options[METIS_OPTION_NUMBERING] = 0;
    this->ncons = 1; // Number of balancing constraints (at least 1)
    
    this->xadj = new idx_t[this->maxNumNodes+1];
    this->parts = new idx_t[this->maxNumNodes];
    this->adjncy = new idx_t[this->maxNumEdges];
    this->adjwgt = new idx_t[this->maxNumEdges];
}

void METISWrapper::setNumCuts(idx_t ncuts)
{
    this->options[METIS_OPTION_NCUTS] = ncuts;
}

void METISWrapper::setNumIterations(idx_t niter)
{
    this->options[METIS_OPTION_NITER] = niter;
}

void METISWrapper::updateNumPartitions(int numParts)
{
    this->nparts = static_cast<idx_t>(numParts);
}

idx_t METISWrapper::partitionSubgraph(idx_t& nvtxs)
{
    idx_t objVal;
    
    /*int status = */METIS_PartGraphKway(
                                     &nvtxs, 
                                     &this->ncons, 
                                     this->xadj, 
                                     this->adjncy, 
                                     NULL, 
                                     NULL, 
                                     this->adjwgt,
                                     &this->nparts, 
                                     NULL, 
                                     NULL, 
                                     this->options, 
                                     &objVal, 
                                     this->parts);
    
    //assert(status == METIS_OK && "METIS partitioning failed");
    
    return objVal;
    
}

idx_t METISWrapper::partitionSubgraphNonContiguous(idx_t& nvtxs)
{
    idx_t objVal;
    
    this->options[METIS_OPTION_CONTIG] = 0;
    
    /*int status = */METIS_PartGraphKway(
                                     &nvtxs, 
                                     &this->ncons, 
                                     this->xadj, 
                                     this->adjncy, 
                                     NULL, 
                                     NULL, 
                                     this->adjwgt,
                                     &this->nparts, 
                                     NULL, 
                                     NULL, 
                                     this->options, 
                                     &objVal, 
                                     this->parts);

    this->options[METIS_OPTION_CONTIG] = 1;
    
    //assert(status == METIS_OK && "METIS partitioning failed");
    
    return objVal;
    
}

void METISWrapper::populateMETISArrays(DynamicGraph& graph, std::unordered_set<NodeID>& subgraphNodes, std::vector<NodeID>& METISIdxToNodeID, bool setEdgeWeightsToOne)
{
    // Assume subgraph size + 1 is equal to size of xadj

    /*
     * Convert Graph to Compressed Storage Format (required by METIS)
     */
    std::unordered_map<NodeID, int> nodeIDToMetisIdx;
    int xadj_pos = 0, adjncy_pos = 0;
    
    for (auto nodeIt = subgraphNodes.begin(); nodeIt != subgraphNodes.end(); ++nodeIt, ++xadj_pos) {
        this->xadj[xadj_pos] = adjncy_pos;
        METISIdxToNodeID[xadj_pos] = *nodeIt; // So we can retrieve subgraph NodeIDs after partitioning from Metis index
        nodeIDToMetisIdx[*nodeIt] = xadj_pos; // So we can convert neighbour NodeIDs to their Metis index
        
        const std::vector<NodeID>& adjNodes = graph.getAdjNodes(*nodeIt);
        const std::vector<EdgeWeight>& adjNodeWgts = graph.getAdjNodeWgts(*nodeIt);
        for (std::size_t i = 0; i < adjNodes.size(); ++i) {
            if (subgraphNodes.find(adjNodes[i]) != subgraphNodes.end()) {
                // We only consider neighbours that are in the subgraph we are partitiong
                this->adjncy[adjncy_pos] = adjNodes[i]; // Note this must to be converted to MetisIdx after
                if (!setEdgeWeightsToOne) {
                    this->adjwgt[adjncy_pos] = adjNodeWgts[i];
                } else {
                    this->adjwgt[adjncy_pos] = 1;
                }
                ++adjncy_pos;
            }
        }
    }
    
    // Must set final position of xadj as end of adjncy
    // I.e. the index after the end of neighbours for the last node
    this->xadj[xadj_pos] = adjncy_pos;

    // Converted neighbour NodeIDs to corresponding MetisIdx
    for (int i = 0; i < adjncy_pos /* == adjncy.size() */; ++i) {
        adjncy[i] = nodeIDToMetisIdx[adjncy[i]];
    }

}

void METISWrapper::populateMETISArrays(Graph_RoadNetwork& graph, std::unordered_set<NodeID>& subgraphNodes, std::vector<NodeID>& METISIdxToNodeID, bool setEdgeWeightsToOne)
{
    // Assume subgraph size + 1 is equal to size of xadj

    /*
     * Convert Graph to Compressed Storage Format (required by METIS)
     */
    std::unordered_map<NodeID, int> nodeIDToMetisIdx;
    int xadj_pos = 0, adjncy_pos = 0;
    int adjListStart, nextAdjListStart;

    for (auto nodeIt = subgraphNodes.begin(); nodeIt != subgraphNodes.end(); ++nodeIt, ++xadj_pos) {
        this->xadj[xadj_pos] = adjncy_pos;
        METISIdxToNodeID[xadj_pos] = *nodeIt; // So we can retrieve subgraph NodeIDs after partitioning from Metis index
        nodeIDToMetisIdx[*nodeIt] = xadj_pos; // So we can convert neighbour NodeIDs to their Metis index
        
        adjListStart = graph.getEdgeListStartIndex(*nodeIt);
        nextAdjListStart = graph.getEdgeListSize(*nodeIt);
        for (int i = adjListStart; i < nextAdjListStart; ++i) {
            if (subgraphNodes.find(graph.edges[i].first) != subgraphNodes.end()) {
                // We only consider neighbours that are in the subgraph we are partitiong
                this->adjncy[adjncy_pos] = graph.edges[i].first; // Note this must to be converted to MetisIdx after
                if (!setEdgeWeightsToOne) {
                    this->adjwgt[adjncy_pos] = graph.edges[i].second;
                } else {
                    this->adjwgt[adjncy_pos] = 1;
                }
                ++adjncy_pos;
            }            
        }
    }
    
    // Must set final position of xadj as end of adjncy
    // I.e. the index after the end of neighbours for the last node
    this->xadj[xadj_pos] = adjncy_pos;

    // Converted neighbour NodeIDs to corresponding MetisIdx
    for (int i = 0; i < adjncy_pos /* == adjncy.size() */; ++i) {
        adjncy[i] = nodeIDToMetisIdx[adjncy[i]];
    }

}
