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

#include "MortonList.h"

#include "DijkstraSearch.h"
#include "../queue/MinPriorityQueue.h"
#include "../queue/BinaryMinHeap.h"
#include "../utility/StopWatch.h"
#include "../utility/geometry.h"
#include "../utility/utility.h"

#include <tuple>
#include <bitset>
#include <iomanip>

/*
 * SILCPathOracle
 */

void SILCPathOracle::buildPathOracle(Quadtree &qt, Graph_RoadNetwork &graph) {
    this->networkName = graph.getNetworkName();
    this->numNodes = graph.getNumNodes();
    this->numEdges = graph.getNumEdges();
    this->maxRange = qt.getMaxRange();
    this->maxRangeExp = qt.getMaxRangeExp();
    this->xTranslation = qt.getXTranslation();
    this->yTranslation = qt.getYTranslation();

    this->relativeXCoordinates = qt.relativeXCoordinates;
    this->relativeYCoordinates = qt.relativeYCoordinates;
    this->mortonCodes = qt.mortonCodes;
    
    // Parallesed Version using OpenMP
    // Note: Quadtree is used by Mortonize operation as read-only
    // operations, so it is fine to shared amongst threads
    DijkstraSearch dijkstra;

    // We know the path oracle will have a MortonList for each node
    this->pathOracle.resize(this->numNodes);
    
    StopWatch sw;
    double dijkstraTime = 0, mortonListTime = 0;
    int totalMortonBlocks = 0;

    #pragma omp parallel for
    for (NodeID node = 0; node < this->numNodes; ++node) {
        // We set the origin to a special colour to preserve it's Morton block
        // as colourizeMap will not give it a colour        
        // Note: colourMap and SSSPDistances will be overwritten each time
        // so we don't waste time clearing them (pqueue is different)
        std::vector<EdgeID> colourMap(this->numNodes);
        std::vector<EdgeWeight> SSSPDistances(this->numNodes);
        bool singleColour = dijkstra.colourizeMap(graph, node, colourMap, SSSPDistances);
        colourMap[node] = quadtree_constants::SPECIAL_COLOUR;
        this->pathOracle[node].Mortonize(node,qt,colourMap,SSSPDistances,singleColour);
        if (node % 50000 == 0) {
            std::cout << "Created Morton Lists for " << node << " nodes" << std::endl;
        }
    }
    
//     int totalMortonBlocks = 0;
//     for (NodeID node = 0; node < this->numNodes; ++node) {
//         totalMortonBlocks += this->pathOracle[node].getNumBlocks();
//     }
//     std::cout << "Total number of Morton blocks created = " << totalMortonBlocks << std::endl;
    
}

std::string SILCPathOracle::getNetworkName() {
    return this->networkName;
}

int SILCPathOracle::getNumNodes() {
    return this->numNodes;
}

int SILCPathOracle::getNumEdges() {
    return this->numEdges;
}

NodeID SILCPathOracle::nextInPath(Graph_RoadNetwork& graph, NodeID intermediate, MortonCode& targetCode, EdgeWeight& edgeWeight) {
    EdgeID pos = this->pathOracle[intermediate].getLink(targetCode);
//     assert(pos != quadtree_constants::SPECIAL_COLOUR && "We are calling nextInPath without checking if intermediate == target");
    NodeID next = graph.getNeighbourAndEdgeWeightByPosition(intermediate,pos,edgeWeight);
    return next;
}

// Note: nextPos is actually the position of the node after next
// because we retrieved while retrieving the ratios for next
void SILCPathOracle::refineDistance(Graph_RoadNetwork& graph, NodeID target, MortonCode& targetCode,
                                    NodeID& next, EdgeWeight& sourceToNextDistance,
                                    DistanceBound &lowerBound, DistanceBound &upperBound) {
    // Note: We should explicitly check that next is not the target
    // before calling refineDistance. Or guarantee that it we will
    // not call refineDistance in that case, because we will 
    // fully refine lowerBound and upperBound as soon as we find
    // the new next node to be the target
    
    DistanceRatio nextMinRatio, nextMaxRatio;

    if (next == target) {
        // Note: In this case nextPost is invalid (it will be invalid because we do not update it)
        // but this is fine because next is also not updated (will always equal target)
        lowerBound = sourceToNextDistance;
        upperBound = sourceToNextDistance;
    } else {
#if defined(COLLECT_STATISTICS)
        this->stats.incrementStatistic("total_refinements",1);
#endif
        EdgeID pos = this->pathOracle[next].getLinkAndRatios(targetCode,nextMinRatio,nextMaxRatio);
        
        EuclideanDistance nextToTargetSpatialDist = graph.getEuclideanDistance(next,target);
        // Note: If next == target the above will be zero and bounds will be the same and equal 
        // to actual network distance from source to target
        // min/max ratio x ds(next,target) + dn(source,intermediate) + w(intermediate,next)
        DistanceBound newLowerBound = std::floor(nextMinRatio*std::floor(nextToTargetSpatialDist)) + sourceToNextDistance;
        DistanceBound newUpperBound = std::ceil(nextMaxRatio*std::ceil(nextToTargetSpatialDist)) + sourceToNextDistance;

        if (lowerBound < newLowerBound) {
            lowerBound = newLowerBound;
        }
        if (upperBound > newUpperBound) {
            upperBound = newUpperBound;
        }
        
//         // DEBUG ONLY
//         // Verify interval is valid by checking shortest path distance to target from next
//         EdgeWeight actualDistance = sourceToNextDistance + this->findShortestPathDistance(graph,next,target);
//         if (lowerBound > actualDistance) {
//             std::cerr << "Lower bound is larger than actual shortest path distance!" << std::endl;
//             std::cerr << "Lower Bound = " << lowerBound << std::endl;
//             std::cerr << "Actual Distance = " << actualDistance << std::endl;
//         }
//         if (upperBound < actualDistance) {
//             std::cerr << "Upper bound is smaller than actual shortest path distance!" << std::endl;
//             std::cerr << "Upper Bound = " << upperBound << std::endl;
//             std::cerr << "Actual Distance = " << actualDistance << std::endl;
//         }

        // While retrieving the min and max ratios for the current "next" node
        // we also retrieved the EdgeID of the node from the current "next"
        // that is subsequent node in the shortest path to the target. So
        // we update next to the this new node and the sourceToNextDistance
        // to match it (i.e. we add the edge weight from current "next")
        
        // This is the next intermediate node, we retrieve it now while we 
        // are retrieving the min and max ratios. We can then use it in the
        // next refinement (this save the nextInPath call in the paper's 
        // version of the algorithm)
        EdgeWeight nextWeight;
        next = graph.getNeighbourAndEdgeWeightByPosition(next,pos,nextWeight);        
        sourceToNextDistance = sourceToNextDistance + nextWeight;
    }
}

// Note: nextPos is actually the position of the node after next
// because we retrieved while retrieving the ratios for next
void SILCPathOracle::optimisedRefineDistance(Graph_RoadNetwork& graph, NodeID target, MortonCode& targetCode,
                                             NodeID& next, EdgeWeight& sourceToNextDistance,
                                             DistanceBound &lowerBound, DistanceBound &upperBound, Junction& junc) {
    // Note: We should explicitly check that next is not the target
    // before calling refineDistance. Or guarantee that it we will
    // not call refineDistance in that case, because we will 
    // fully refine lowerBound and upperBound as soon as we find
    // the new next node to be the target
    
    DistanceRatio nextMinRatio, nextMaxRatio;

    if (next == target) {
        // Note: In this case nextPost is invalid (it will be invalid because we do not update it)
        // but this is fine because next is also not updated (will always equal target)
        lowerBound = sourceToNextDistance;
        upperBound = sourceToNextDistance;
    } else {
#if defined(COLLECT_STATISTICS)
        this->stats.incrementStatistic("total_refinements",1);
#endif
        EdgeID pos = this->pathOracle[next].getLinkAndRatios(targetCode,nextMinRatio,nextMaxRatio);
        
        EuclideanDistance nextToTargetSpatialDist = graph.getEuclideanDistance(next,target);
        // Note: If next == target the above will be zero and bounds will be the same and equal 
        // to actual network distance from source to target
        // min/max ratio x ds(next,target) + dn(source,intermediate) + w(intermediate,next)
        DistanceBound newLowerBound = std::floor(nextMinRatio*std::floor(nextToTargetSpatialDist)) + sourceToNextDistance;
        DistanceBound newUpperBound = std::ceil(nextMaxRatio*std::ceil(nextToTargetSpatialDist)) + sourceToNextDistance;

        if (lowerBound < newLowerBound) {
            lowerBound = newLowerBound;
        }
        if (upperBound > newUpperBound) {
            upperBound = newUpperBound;
        }
        
//         // DEBUG ONLY
//         // Verify interval is valid by checking shortest path distance to target from next
//         EdgeWeight actualDistance = sourceToNextDistance + this->findShortestPathDistance(graph,next,target);
//         if (lowerBound > actualDistance) {
//             std::cerr << "Lower bound is larger than actual shortest path distance!" << std::endl;
//             std::cerr << "Lower Bound = " << lowerBound << std::endl;
//             std::cerr << "Actual Distance = " << actualDistance << std::endl;
//         }
//         if (upperBound < actualDistance) {
//             std::cerr << "Upper bound is smaller than actual shortest path distance!" << std::endl;
//             std::cerr << "Upper Bound = " << upperBound << std::endl;
//             std::cerr << "Actual Distance = " << actualDistance << std::endl;
//         }

        // While retrieving the min and max ratios for the current "next" node
        // we also retrieved the EdgeID of the node from the current "next"
        // that is subsequent node in the shortest path to the target. So
        // we update next to the this new node and the sourceToNextDistance
        // to match it (i.e. we add the edge weight from current "next")
        
        // This is the next intermediate node, we retrieve it now while we 
        // are retrieving the min and max ratios. We can then use it in the
        // next refinement (this save the nextInPath call in the paper's 
        // version of the algorithm)
        EdgeWeight nextWeight, sourceToNewNextDistance;
        NodeID edgeToNextIndex, prevNode = next, adjNode; // edgeToNextIndex is the edge we follow to reach new next node
        next = graph.getNeighbourAndEdgeWeightByPosition(next,pos,nextWeight,edgeToNextIndex);
        sourceToNextDistance = sourceToNextDistance + nextWeight;
        
        // See if we can use junction index to move next node along a bit further
        int outdegree, adjListStart, nextAdjListStart;
        if (junc.isNoThruRoadNode[next]) {
            // If the next node is a no through road node then we search opposite from the 
            // direction we came, which either towards a target or towards a junction
            // so do not assume we do not start on a non-no thru road node
            adjListStart = graph.getEdgeListStartIndex(next);
            nextAdjListStart = graph.getEdgeListSize(next);
            outdegree = nextAdjListStart - adjListStart;
            while (next != target && outdegree != 2) {
                if (graph.edges[adjListStart].first != prevNode) {
                    adjNode = graph.edges[adjListStart].first;
                    sourceToNextDistance += graph.edges[adjListStart].second;
                } else {
                    adjNode = graph.edges[adjListStart+1].first; 
                    sourceToNextDistance += graph.edges[adjListStart].second;
                }
                prevNode = next;
                next = adjNode;
                adjListStart = graph.getEdgeListStartIndex(next);                
                nextAdjListStart = graph.getEdgeListSize(next);
                outdegree = nextAdjListStart - adjListStart;
            }
            return;
        } else {
            adjListStart = graph.getEdgeListStartIndex(next);
            nextAdjListStart = graph.getEdgeListSize(next);
            outdegree = nextAdjListStart - adjListStart;
            if (outdegree == 2) {
                // As long as the distance to the next junction is smaller than the lowerbound
                // then we can skip this "road" and go straight to the next junction
                // Note: outdegree cannot be 1 because then next would have to be target if we followed some edge
                sourceToNewNextDistance = sourceToNextDistance - nextWeight + junc.nextJunctionNodeAndDistance[edgeToNextIndex].second;
                // Note: The distance to the junction is from the previous node, so we subtract nextWeight
                if (lowerBound >= sourceToNewNextDistance) {
                    // Here the new next node is the junction from the current next node
                    // because the lower bound is larger than the distance to it (or 
                    // the junction is the target)
                    next = junc.nextJunctionNodeAndDistance[edgeToNextIndex].first;
                    sourceToNextDistance = sourceToNewNextDistance;
                } else {
                    // If the distance to junction is greater than or equal to the lowerbound
                    // this means we can incrementally move towards the junction
                    if (upperBound <= sourceToNewNextDistance) {
                        sourceToNewNextDistance = sourceToNextDistance; // Reset this to be distance to actual new next
                        // This means the target is on the road between next and junction
                        while (next != target) {
                            if (graph.edges[adjListStart].first != prevNode) {
                                adjNode = graph.edges[adjListStart].first;
                                sourceToNewNextDistance = sourceToNewNextDistance + graph.edges[adjListStart].second;
                            } else {
                                adjNode = graph.edges[adjListStart+1].first; 
                                sourceToNewNextDistance = sourceToNewNextDistance + graph.edges[adjListStart+1].second;
                            }
                            prevNode = next;
                            next = adjNode;
                            sourceToNextDistance = sourceToNewNextDistance;
                            adjListStart = graph.getEdgeListStartIndex(next);
                        }
                    } else {
                        NodeID junctionNode = junc.nextJunctionNodeAndDistance[edgeToNextIndex].first;
                        sourceToNewNextDistance = sourceToNextDistance; // Reset this to be distance to actual new next
                        // This means either the target could be on the road between next and junction
                        // or it is past the junction so we incrementally search towards junction 
                        // stopping if we find target
                        while (next != target && next != junctionNode) {
                            if (graph.edges[adjListStart].first != prevNode) {
                                adjNode = graph.edges[adjListStart].first;
                                sourceToNewNextDistance = sourceToNewNextDistance + graph.edges[adjListStart].second;
                            } else {
                                adjNode = graph.edges[adjListStart+1].first; 
                                sourceToNewNextDistance = sourceToNewNextDistance + graph.edges[adjListStart+1].second;
                            }
                            prevNode = next;
                            next = adjNode;
                            sourceToNextDistance = sourceToNewNextDistance;
                            adjListStart = graph.getEdgeListStartIndex(next);
                        }
                    }
                    // Either next == target or next == junctionNode. In either case
                    // sourceToNextDistance could improve the lowerBound. E.g. if next 
                    // is the junctionNode then we already the lowerBound is smaller
                    // than the distance to the junctionNode
                    if (lowerBound < sourceToNextDistance) {
                        lowerBound = sourceToNextDistance;
                    }
                    // If we did reach target then we 
                    if (next == target && sourceToNextDistance > upperBound) {
                        upperBound = sourceToNextDistance;
                    }
                    // Note: The above will help either or both SILC or Distance Browsing 
                    // kNN because they use LB or UP (respectively) to differentiate neighbours
                }
            }
        }
    }
}

Path SILCPathOracle::findShortestPath(Graph_RoadNetwork& graph, NodeID source, NodeID target) {
    Path path(source,true,0);
 
    NodeID intermediate;
    EdgeWeight distance = 0, weight;
    intermediate = source;
    while (intermediate != target) {
        intermediate = this->nextInPath(graph,intermediate,this->mortonCodes[target],weight);
        distance += weight;
        path.addToEnd(intermediate,weight);
    }
    path.setPathLength(distance);
    
    return path;
}

EdgeWeight SILCPathOracle::findShortestPathDistance(Graph_RoadNetwork& graph, NodeID source, NodeID target) {
    NodeID intermediate;
    EdgeWeight distance = 0, weight;
    intermediate = source;
    while (intermediate != target) {
        intermediate = this->nextInPath(graph,intermediate,this->mortonCodes[target],weight);
        distance += weight;
    }
    
    return distance;
}

Path SILCPathOracle::findShortestPathOptimised(Graph_RoadNetwork& graph, NodeID source, NodeID target) {
    Path path(source,true,0);

    NodeID intermediate, prev;
    EdgeWeight distance = 0, weight;
    int adjListStart, nextAdjListStart;
    intermediate = source;
    if (intermediate != target) {
        // We need to handle case where source node has only two edges
        // in this case we need to do a nextInPath call so that we no
        // which direction to head in
        prev = intermediate;
        intermediate = this->nextInPath(graph,intermediate,this->mortonCodes[target],weight);
        distance += weight;
        while (intermediate != target) {
            adjListStart = graph.getEdgeListStartIndex(intermediate);
            nextAdjListStart = graph.getEdgeListSize(intermediate);
            if (nextAdjListStart-adjListStart == 2) {
                if (graph.edges[adjListStart].first != prev) {
                    prev = intermediate;
                    intermediate = graph.edges[adjListStart].first;
                    weight = graph.edges[adjListStart].second;
                } else {
                    prev = intermediate;
                    intermediate = graph.edges[adjListStart+1].first;
                    weight = graph.edges[adjListStart+1].second;
                }
            } else {
                prev = intermediate;
                intermediate = this->nextInPath(graph,intermediate,this->mortonCodes[target],weight);
            }
            distance += weight;
            path.addToEnd(intermediate,weight);
        }
    }

    return path;
}

EdgeWeight SILCPathOracle::findShortestPathDistanceOptimised(Graph_RoadNetwork& graph, NodeID source, NodeID target)
{
    NodeID intermediate, prev;
    EdgeWeight distance = 0, weight;
    int adjListStart, nextAdjListStart;
    intermediate = source;
    if (intermediate != target) {
        // We need to handle case where source node has only two edges
        // in this case we need to do a nextInPath call so that we no
        // which direction to head in
        prev = intermediate;
        intermediate = this->nextInPath(graph,intermediate,this->mortonCodes[target],weight);
        distance += weight;
        while (intermediate != target) {
            adjListStart = graph.getEdgeListStartIndex(intermediate);
            nextAdjListStart = graph.getEdgeListSize(intermediate);
            if (nextAdjListStart-adjListStart == 2) {
                if (graph.edges[adjListStart].first != prev) {
                    prev = intermediate;
                    intermediate = graph.edges[adjListStart].first;
                    distance += graph.edges[adjListStart].second;
                } else {
                    prev = intermediate;
                    intermediate = graph.edges[adjListStart+1].first;
                    distance += graph.edges[adjListStart+1].second;
                }
            } else {
                prev = intermediate;
                intermediate = this->nextInPath(graph,intermediate,this->mortonCodes[target],weight);
                distance += weight;
            }
        }
    }

    return distance;
}

void SILCPathOracle::getKNNs(SimpleQuadtree& objectHierarchy, unsigned int k, NodeID queryNode, std::vector<NodeID>& kNNs, 
                             std::vector<EdgeWeight>& kNNDistances, Graph_RoadNetwork& graph)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("total_regions",objectHierarchy.getTotalRegions());
    this->stats.initialiseStatistic("regions_inserted",0);
    this->stats.initialiseStatistic("objects_inserted",0);
    this->stats.initialiseStatistic("total_refinements",0);
#endif
    DistanceBound lowerBound, upperBound;
    BinaryMinHeap<DistanceBound,DataTuple> pqueue;
    NodeID object, next;
    EdgeWeight sourceToNextDistance;

    // Initialise with entire region encompassing object hierarchy
    this->intervalDistance(queryNode,objectHierarchy,objectHierarchy.rootNode,lowerBound,upperBound);
    pqueue.insert(DataTuple(lowerBound,upperBound,queryNode,0,objectHierarchy.rootNode),lowerBound);
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("regions_inserted",1);
#endif
    
    while(pqueue.size() > 0) {
        DataTuple minElement = pqueue.extractMinElement();
        if (minElement.isBlock) {
            if (!minElement.region->isLeafNode()) {
                // If it is not a leaf node then we add each child which contains objects
                if (minElement.region->sw->getNumNodes() > 0) {
                    // We also do not want to add any empty blocks which can exist since object 
                    // hierarchy is not a linear Quadtree (obviously it will not have children)
                    this->intervalDistance(minElement.intermediate,objectHierarchy,minElement.region->sw,lowerBound,upperBound);
                    //assert(upperBound >= lowerBound && "Upper bound smaller than lower bound");
                    pqueue.insert(DataTuple(lowerBound,upperBound,minElement.intermediate,minElement.distance,minElement.region->sw),lowerBound);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("regions_inserted",1);
#endif
                    // Note: In case of queue elements that are regions minElement.intermediate 
                    // will always be equal to queryNode and minElement.distance will be 0
                }
                if (minElement.region->se->getNumNodes() > 0) {
                    this->intervalDistance(minElement.intermediate,objectHierarchy,minElement.region->se,lowerBound,upperBound);
                    //assert(upperBound >= lowerBound && "Upper bound smaller than lower bound");
                    pqueue.insert(DataTuple(lowerBound,upperBound,minElement.intermediate,0,minElement.region->se),lowerBound);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("regions_inserted",1);
#endif
                }
                if (minElement.region->nw->getNumNodes() > 0) {
                    this->intervalDistance(minElement.intermediate,objectHierarchy,minElement.region->nw,lowerBound,upperBound);
                    //assert(upperBound >= lowerBound && "Upper bound smaller than lower bound");
                    pqueue.insert(DataTuple(lowerBound,upperBound,minElement.intermediate,0,minElement.region->nw),lowerBound);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("regions_inserted",1);
#endif
                }
                if (minElement.region->ne->getNumNodes() > 0) {
                    this->intervalDistance(minElement.intermediate,objectHierarchy,minElement.region->ne,lowerBound,upperBound);
                    //assert(upperBound >= lowerBound && "Upper bound smaller than lower bound");
                    pqueue.insert(DataTuple(lowerBound,upperBound,minElement.intermediate,0,minElement.region->ne),lowerBound);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("regions_inserted",1);
#endif
                }
            } else {
                SimpleQuadtreeNode *leafRegion = minElement.region;
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("objects_inserted",leafRegion->nodes.size());
#endif
                for (std::size_t i = 0; i < leafRegion->nodes.size(); ++i) {
                    object = leafRegion->nodes[i];
                    next = minElement.intermediate; // This will be the queryNode
                    sourceToNextDistance = 0; // next is source, so source to source dist is 0
                    lowerBound = minElement.lowerBound; // We inherit the dist bounds to the intersecting region containing element
                    upperBound = minElement.upperBound;
                    // We provide the source (queryNode) as the next node in the path from the source to the
                    // to the object and the lowerbound and upperbound inherited from the intersecting region
                    // so as far as refineDistance is concerned, all input values are valid and it will update
                    // next and sourceToNextDistance to be the first edge and weight in the path the object
                    // unless the target is the query node (i.e. object == next) in which we won't refine 
                    // further anyway.                    
                    this->refineDistance(graph,object,this->mortonCodes[object],next,sourceToNextDistance,lowerBound,upperBound);
                    pqueue.insert(DataTuple(lowerBound,upperBound,next,sourceToNextDistance,object),lowerBound);
                }
            }
        } else {
            if (minElement.upperBound > pqueue.getMinKey()) {
                // Note: We change the defintion of INTERSECTS from the paper. In the paper
                // we should reach here even if minElement.upperBound == pqueue.getMinKey()
                // However in this case, pqueue.getMinKey() provides the lower bound of the
                // next smallest element. So if minElement.upperBound == lower bound of the
                // next smallest element, then the next smallest element cannot be a closer kNN
                // to the query node than the current minElement.point. So it is safe to 
                // report the minElement.point as a kNN. Refining further would be a waste.
                lowerBound = minElement.lowerBound;
                upperBound = minElement.upperBound;
                // Note: If we reach here then minElement.intermediate cannot equal minElement.point
                // because this can only happen when an object has been fully refined. In that case
                // minElement.lowerBound == minElement.upperBound, and it would be smaller than
                // or equal the lower bound of any other element in the queue and we would report 
                // as a result. Assert above checks this is the case.
                next = minElement.intermediate;
                sourceToNextDistance = minElement.distance;
                //assert(minElement.upperBound != minElement.lowerBound && "All fully refined object should be reported as results");
                this->refineDistance(graph,minElement.point,this->mortonCodes[minElement.point],next,sourceToNextDistance,lowerBound,upperBound);
                //assert(upperBound >= lowerBound && "Upper bound smaller than lower bound");
                pqueue.insert(DataTuple(lowerBound,upperBound,next,sourceToNextDistance,minElement.point),lowerBound);
            } else {
                //assert(minElement.upperBound >= minElement.lowerBound && "Upper bound bigger than lower bound");
                kNNs.push_back(minElement.point);
                kNNDistances.push_back(minElement.lowerBound);
                if (kNNs.size() == k) {
                    break;
                }
            }
        }
    }
}

void SILCPathOracle::getKNNsByDistanceBrowsing(SimpleQuadtree& objectHierarchy, unsigned int k, NodeID queryNode, std::vector<NodeID>& kNNs, 
                                               std::vector<EdgeWeight>& kNNDistances, Graph_RoadNetwork& graph)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("total_regions",objectHierarchy.getTotalRegions());
    this->stats.initialiseStatistic("regions_inserted",0);
    this->stats.initialiseStatistic("objects_inserted",0);
    this->stats.initialiseStatistic("total_refinements",0);
#endif
    DistanceBound lowerBound, upperBound;
    BinaryMinHeap<DistanceBound,DataTuple> pqueue;
    BinaryMaxHeapWithDK<DistanceBound,NodeID> knnCandidates;
    knnCandidates.init(k);
    NodeID object, next;
    EdgeWeight sourceToNextDistance;

    DistanceBound Dk = 0;
    bool DkInfinity = true;
    
    // Initialise with entire region encompassing object hierarchy
    this->intervalDistanceForDistBrws(queryNode,objectHierarchy,objectHierarchy.rootNode,lowerBound);
    pqueue.insert(DataTuple(lowerBound,0,queryNode,0,objectHierarchy.rootNode),lowerBound);
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("regions_inserted",1);
#endif
    
    while(pqueue.size() > 0) {
        DataTuple minElement = pqueue.extractMinElement();
        if (!minElement.isBlock) {
            // Min element is an object
            if (!DkInfinity && minElement.lowerBound >= Dk) {
                // If Dk is infinity, then lowerbound cannot be greater
                break;
            } else if (minElement.upperBound > pqueue.getMinKey() || (minElement.upperBound == pqueue.getMinKey() && minElement.lowerBound != minElement.upperBound)) {
                // The >= condition in the above statement (accordining to the original paper) has been changed 
                // to a > condition (and the additional condition handles an edge-case which can cause misorderings)
                // because it potentially leads to an infinite loop. This happens because if minElement.upperBound
                // == pqueue.getMinKey() then we will try refine it, even if it has been fully refined.
                // As long as the lowerbound <= Dk we will re-insert into pqueue. Since it is fully 
                // refined, it will still be the minmum element in pqueue and we have an infinite loop.
                // The second condition is added because we need to fully refine two objects with the
                // same upperbound otherwise they will remain in kNNCandidates in the wrong order.
                
                if (DkInfinity || minElement.upperBound <= Dk) {
                    // We try to remove element, if upperbound is smaller than or equal to Dk
                    // then it could already in knnCandidates
                    if (knnCandidates.contains(minElement.point)) {
                        knnCandidates.deleteElement(minElement.point);
                    }
                }
                
                // Note: If we reach here then minElement.intermediate cannot equal minElement.point
                // because this can only happen when an object has been fully refined. In that case
                // minElement.lowerBound == minElement.upperBound (see refineDistance), so since it's
                // the minElement in the queue we must also have minElement.upperBound <= pqueue.getMinKey()
                // and therefore it would not reach here. Note that in the case of the second part of the 
                // condition it excludes element that have been fully refined (would equal Dk in that case).
                next = minElement.intermediate;
                sourceToNextDistance = minElement.distance;
                lowerBound = minElement.lowerBound;
                upperBound = minElement.upperBound;
                this->refineDistance(graph,minElement.point,this->mortonCodes[minElement.point],next,sourceToNextDistance,lowerBound,upperBound);
                // Note: refineDistance will update next and sourceToNextDistance to the next intermediate node
                // and it will also retrieve loweBound and upperBound for the current intermediate node (next)

                if (DkInfinity || upperBound <= Dk) {
                    knnCandidates.insert(minElement.point,upperBound);
                    if (knnCandidates.size() == k) {
                        Dk = knnCandidates.getMaxKey();
                        DkInfinity = false;
                    } else if (knnCandidates.size() > k) {
                        // Note: knnCandidates cannot have more than k elements
                        // because each time we insert we check if size > k and 
                        // remove the largest element if yes
                        knnCandidates.extractMaxElement();
                        Dk = knnCandidates.getMaxKey();
                        //assert(!DkInfinity && "Dk cannot be infinity because k candidates have already been found");
                    }
                }
                
                // We have that Dk >= pqueue.getMinKey. This because the element associated with Dk must be in the queue
                // otherwise the algorithm will stop (because if it is not in pqueue then that means Dk < pqueue.getMinKey).
                // But if the element associated with Dk has been fully refined then lowerBound == Dk below, and this element
                // would not be re-insert into the queue acording to the original algorithm in the paper. Then for the next 
                // dequeued (which is the last candidate for the kNN will not be considered and we may potentially lose a 
                // kNN result). To fix this we change the < to be <= to below.
                if (DkInfinity || lowerBound <= Dk) {
                    pqueue.insert(DataTuple(lowerBound,upperBound,next,sourceToNextDistance,minElement.point),lowerBound);
                }
            } else {
            }
        } else if (!DkInfinity && minElement.lowerBound >= Dk) {
            // Min element is not an object and it cannot have an object that can be a kNN
            // and since it is the min element no other element can be or have an kNN
            break;
        } else {
            // Min element is a block
            if (minElement.region->isLeafNode()) {
                // Min element is a leaf block
                SimpleQuadtreeNode *leafRegion = minElement.region;
                for (std::size_t i = 0; i < leafRegion->nodes.size(); ++i) {
                    object = leafRegion->nodes[i];
                    next = minElement.intermediate; // This will be the queryNode
                    sourceToNextDistance = 0; // next is source, so source to source dist is 0
                    lowerBound = minElement.lowerBound;
                    // Note: Can't inherit upper bound from intersecting region containing element
                    // because we don't set this for any object hierarchy region in dist brws 
                    upperBound = std::numeric_limits<EdgeWeight>::max(); 
                    // We provide the source (queryNode) as the next node in the path from the source to the
                    // to the object and the lowerbound and upperbound inherited from the intersecting region
                    // so as far as refineDistance is concerned, all input values are valid and it will update
                    // next and sourceToNextDistance to be the first edge and weight in the path to the object
                    // unless the target is the query node (i.e. object == next) in which we won't refine 
                    // further anyway.
                    this->refineDistance(graph,object,this->mortonCodes[object],next,sourceToNextDistance,lowerBound,upperBound);
                    if (DkInfinity || lowerBound < Dk) {
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("objects_inserted",1);
#endif
                        pqueue.insert(DataTuple(lowerBound,upperBound,next,sourceToNextDistance,object),lowerBound);
                        if (DkInfinity || upperBound < Dk) {
                            knnCandidates.insert(object,upperBound);
                            if (knnCandidates.size() == k) {
                                Dk = knnCandidates.getMaxKey();
                                DkInfinity = false;
                            } else if (knnCandidates.size() > k) {
                                // Note: knnCandidates cannot have more than k elements
                                // because each time we insert we check if size > k and 
                                // remove the largest element if yes
                                knnCandidates.extractMaxElement();
                                Dk = knnCandidates.getMaxKey();
                                //assert(!DkInfinity && "Dk cannot be infinity because k candidates have already been found");
                            }
                        }
                    }
                }
            } else {
                // Min element is a non-leaf block
                if (minElement.region->sw->getNumNodes() > 0) {
                    // We also do not want to add any empty blocks which can exist since object 
                    // hierarchy is not a linear Quadtree (obviously it will not have children)
                    this->intervalDistanceForDistBrws(minElement.intermediate,objectHierarchy,minElement.region->sw,lowerBound);
                    if (DkInfinity || lowerBound < Dk) {
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("regions_inserted",1);
#endif
                        pqueue.insert(DataTuple(lowerBound,0,minElement.intermediate,minElement.distance,minElement.region->sw),lowerBound);
                    }
                    // Note: In case of queue elements that are regions minElement.intermediate 
                    // will always be equal to queryNode and minElement.distance will be 0
                }
                if (minElement.region->se->getNumNodes() > 0) {
                    this->intervalDistanceForDistBrws(minElement.intermediate,objectHierarchy,minElement.region->se,lowerBound);
                    if (DkInfinity || lowerBound < Dk) {
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("regions_inserted",1);
#endif
                        pqueue.insert(DataTuple(lowerBound,0,minElement.intermediate,0,minElement.region->se),lowerBound);
                    }
                }
                if (minElement.region->nw->getNumNodes() > 0) {
                    this->intervalDistanceForDistBrws(minElement.intermediate,objectHierarchy,minElement.region->nw,lowerBound);
                    if (DkInfinity || lowerBound < Dk) {
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("regions_inserted",1);
#endif
                        pqueue.insert(DataTuple(lowerBound,0,minElement.intermediate,0,minElement.region->nw),lowerBound);
                    }
                }
                if (minElement.region->ne->getNumNodes() > 0) {
                    this->intervalDistanceForDistBrws(minElement.intermediate,objectHierarchy,minElement.region->ne,lowerBound);
                    if (DkInfinity || lowerBound < Dk) {
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("regions_inserted",1);
#endif
                        pqueue.insert(DataTuple(lowerBound,0,minElement.intermediate,0,minElement.region->ne),lowerBound);
                    }
                }
            }
        }
    }
    
    // Note: It's possible that knnCandidates size is less than k if there are
    // less than k objects in the whole road network
    knnCandidates.populateKNNs(kNNs,kNNDistances);
}

void SILCPathOracle::getKNNsByOptimisedSILC(SimpleQuadtree& objectHierarchy, unsigned int k, NodeID queryNode, std::vector<NodeID>& kNNs, 
                                            std::vector<EdgeWeight>& kNNDistances, Graph_RoadNetwork& graph, Junction& junc)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("total_regions",objectHierarchy.getTotalRegions());
    this->stats.initialiseStatistic("regions_inserted",0);
    this->stats.initialiseStatistic("objects_inserted",0);
    this->stats.initialiseStatistic("total_refinements",0);
#endif
    DistanceBound lowerBound, upperBound;
    BinaryMinHeap<DistanceBound,DataTuple> pqueue;
    NodeID object, next;
    EdgeWeight sourceToNextDistance;
    
    // Initialise with entire region encompassing object hierarchy
    this->intervalDistance(queryNode,objectHierarchy,objectHierarchy.rootNode,lowerBound,upperBound);
    pqueue.insert(DataTuple(lowerBound,upperBound,queryNode,0,objectHierarchy.rootNode),lowerBound);
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("regions_inserted",1);
#endif
    
    while(pqueue.size() > 0) {
        DataTuple minElement = pqueue.extractMinElement();
        if (minElement.isBlock) {
            if (!minElement.region->isLeafNode()) {
                // If it is not a leaf node then we add each child which contains objects
                if (minElement.region->sw->getNumNodes() > 0) {
                    // We also do not want to add any empty blocks which can exist since object 
                    // hierarchy is not a linear Quadtree (obviously it will not have children)
                    this->intervalDistance(minElement.intermediate,objectHierarchy,minElement.region->sw,lowerBound,upperBound);
                    //assert(upperBound >= lowerBound && "Upper bound smaller than lower bound");
                    pqueue.insert(DataTuple(lowerBound,upperBound,minElement.intermediate,minElement.distance,minElement.region->sw),lowerBound);
                    // Note: In case of queue elements that are regions minElement.intermediate 
                    // will always be equal to junctionNode and minElement.distance will be 0
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("regions_inserted",1);
#endif
                }
                if (minElement.region->se->getNumNodes() > 0) {
                    this->intervalDistance(minElement.intermediate,objectHierarchy,minElement.region->se,lowerBound,upperBound);
                    //assert(upperBound >= lowerBound && "Upper bound smaller than lower bound");
                    pqueue.insert(DataTuple(lowerBound,upperBound,minElement.intermediate,minElement.distance,minElement.region->se),lowerBound);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("regions_inserted",1);
#endif
                }
                if (minElement.region->nw->getNumNodes() > 0) {
                    this->intervalDistance(minElement.intermediate,objectHierarchy,minElement.region->nw,lowerBound,upperBound);
                    //assert(upperBound >= lowerBound && "Upper bound smaller than lower bound");
                    pqueue.insert(DataTuple(lowerBound,upperBound,minElement.intermediate,minElement.distance,minElement.region->nw),lowerBound);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("regions_inserted",1);
#endif
                }
                if (minElement.region->ne->getNumNodes() > 0) {
                    this->intervalDistance(minElement.intermediate,objectHierarchy,minElement.region->ne,lowerBound,upperBound);
                    //assert(upperBound >= lowerBound && "Upper bound smaller than lower bound");
                    pqueue.insert(DataTuple(lowerBound,upperBound,minElement.intermediate,minElement.distance,minElement.region->ne),lowerBound);
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("regions_inserted",1);
#endif
                }
            } else {
                SimpleQuadtreeNode *leafRegion = minElement.region;
                for (std::size_t i = 0; i < leafRegion->nodes.size(); ++i) {
                    object = leafRegion->nodes[i];
#if defined(COLLECT_STATISTICS)
                    this->stats.incrementStatistic("objects_inserted",1);
#endif
                    next = minElement.intermediate; // This will be the queryNode
                    sourceToNextDistance = 0; // next is source, so source to source dist is 0
                    lowerBound = std::floor(graph.getEuclideanDistance(next,object));
                    if (lowerBound < minElement.lowerBound) {
                        // Note: The minElement.lowerbound is currently the shortest lower bound for the block.
                        // I.e. this is shortest possible path to any vertex in the block. This may in fact
                        // be higher than the Euclidean Distance to the vertex, so we check before setting.
                        lowerBound = minElement.lowerBound; // We inherit the dist bounds to the intersecting region containing element
                    }
                    upperBound = minElement.upperBound;
                    // We provide the source (junctionNode) as the next node in the path from the source to the
                    // to the object and the lowerbound and upperbound inherited from the intersecting region
                    // so as far as refineDistance is concerned, all input values are valid and it will update
                    // next and sourceToNextDistance to be the first edge and weight in the path the object
                    // unless the target is the query node (i.e. object == next) in which we won't refine 
                    // further anyway.
                    pqueue.insert(DataTuple(lowerBound,upperBound,next,sourceToNextDistance,object),lowerBound);
                    // Note: Unlike the original algorithm, we avoid the initial refinement, as it is not a trivial 
                    // cost and the Euclidean distance lowerbound or the lower bound inherited from the leaf block
                    // (whichever is larger) may be enough to avoid the object in future (at worst we will have to
                    // refine it later, so we essentially delay refinement as long as possible).
                }
            }
        } else {
            if (minElement.upperBound > pqueue.getMinKey()) {
                // Note: We change the defintion of INTERSECTS from the paper. In the paper
                // we should reach here even if minElement.upperBound == pqueue.getMinKey()
                // However in this case, pqueue.getMinKey() provides the lower bound of the
                // next smallest element. So if minElement.upperBound == lower bound of the
                // next smallest element, then the next smallest element cannot be a closer kNN
                // to the query node than the current minElement.point. So it is safe to 
                // report the minElement.point as a kNN. Refining further would be a waste.
                lowerBound = minElement.lowerBound;
                upperBound = minElement.upperBound;
                // Note: If we reach here then minElement.intermediate cannot equal minElement.point
                // because this can only happen when an object has been fully refined. In that case
                // minElement.lowerBound == minElement.upperBound, and it would be smaller than
                // or equal the lower bound of any other element in the queue and we would report 
                // as a result. Assert above checks this is the case.
                next = minElement.intermediate;
                sourceToNextDistance = minElement.distance;
                //assert(minElement.upperBound != minElement.lowerBound && "All fully refined object should be reported as results");
#if defined(NO_CHAIN_SKIPPING)
                this->refineDistance(graph,minElement.point,this->mortonCodes[minElement.point],next,sourceToNextDistance,lowerBound,upperBound);
#else
                this->optimisedRefineDistance(graph,minElement.point,this->mortonCodes[minElement.point],next,sourceToNextDistance,lowerBound,upperBound,junc);
#endif
                //assert(upperBound >= lowerBound && "Upper bound smaller than lower bound");
                pqueue.insert(DataTuple(lowerBound,upperBound,next,sourceToNextDistance,minElement.point),lowerBound);
            } else {
                //assert(minElement.upperBound >= minElement.lowerBound && "Upper bound bigger than lower bound");
                kNNs.push_back(minElement.point);
                kNNDistances.push_back(minElement.lowerBound);
                if (kNNs.size() == k) {
                    break;
                }                
            }
        }
    }
}

void SILCPathOracle::getKNNsByOptimisedDistanceBrowsing(SimpleQuadtree& objectHierarchy, unsigned int k, NodeID queryNode, std::vector<NodeID>& kNNs, 
                                                        std::vector<EdgeWeight>& kNNDistances, Graph_RoadNetwork& graph, Junction& junc)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("total_regions",objectHierarchy.getTotalRegions());
    this->stats.initialiseStatistic("regions_inserted",0);
    this->stats.initialiseStatistic("objects_inserted",0);
    this->stats.initialiseStatistic("total_refinements",0);
#endif
    DistanceBound lowerBound, upperBound;
    BinaryMinHeap<DistanceBound,DataTuple> pqueue;
    BinaryMaxHeapWithDK<DistanceBound,NodeID> knnCandidates;
    knnCandidates.init(k);
    NodeID object, next, iK = k;
    EdgeWeight sourceToNextDistance;

    DistanceBound Dk = 0;
    bool DkInfinity = true;

    // Initialise with entire region encompassing object hierarchy
    this->intervalDistance(queryNode,objectHierarchy,objectHierarchy.rootNode,lowerBound,upperBound);
    pqueue.insert(DataTuple(lowerBound,upperBound,queryNode,0,objectHierarchy.rootNode),lowerBound);
#if defined(COLLECT_STATISTICS)
    this->stats.incrementStatistic("regions_inserted",1);
#endif
    
    while(pqueue.size() > 0) {
        DataTuple minElement = pqueue.extractMinElement();
        if (!minElement.isBlock) {
            // Min element is an object
            if (!DkInfinity && minElement.lowerBound >= Dk) {
                // If Dk is infinity, then lowerbound cannot be greater
                break;
            } else if (minElement.upperBound > pqueue.getMinKey() || (minElement.upperBound == pqueue.getMinKey() && minElement.lowerBound != minElement.upperBound)) {
                // The >= condition in the above statement (accordining to the original paper) has been changed 
                // to a > condition (and the additional condition handles an edge-case which can cause misorderings)
                // because it potentially leads to an infinite loop. This happens because if minElement.upperBound
                // == pqueue.getMinKey() then we will try refine it, even if it has been fully refined.
                // As long as the lowerbound <= Dk we will re-insert into pqueue. Since it is fully 
                // refined, it will still be the minmum element in pqueue and we have an infinite loop.
                // The second condition is added because we need to fully refine two objects with the
                // same upperbound otherwise they will remain in kNNCandidates in the wrong order.
                
                if (DkInfinity || minElement.upperBound <= Dk) {
                    // We try to remove element, if upperbound is smaller than or equal to Dk
                    // then it could already in knnCandidates
                    if (knnCandidates.contains(minElement.point)) {
                        knnCandidates.deleteElement(minElement.point);
                    }
                }
                
                // Note: If we reach here then minElement.intermediate cannot equal minElement.point
                // because this can only happen when an object has been fully refined. In that case
                // minElement.lowerBound == minElement.upperBound (see refineDistance), so since it's
                // the minElement in the queue we must also have minElement.upperBound <= pqueue.getMinKey()
                // and therefore it would not reach here. Note that in the case of the second part of the 
                // condition it excludes element that have been fully refined (would equal Dk in that case).
                next = minElement.intermediate;
                sourceToNextDistance = minElement.distance;
                lowerBound = minElement.lowerBound;
                upperBound = minElement.upperBound;
#if defined(NO_CHAIN_SKIPPING)
                this->refineDistance(graph,minElement.point,this->mortonCodes[minElement.point],next,sourceToNextDistance,lowerBound,upperBound);
#else
                this->optimisedRefineDistance(graph,minElement.point,this->mortonCodes[minElement.point],next,sourceToNextDistance,lowerBound,upperBound,junc);
#endif
                // Note: refineDistance will update next and sourceToNextDistance to the next intermediate node
                // and it will also retrieve loweBound and upperBound for the current intermediate node (next)

                if (DkInfinity || upperBound <= Dk) {
                    knnCandidates.insert(minElement.point,upperBound);
                    if (knnCandidates.size() == k) {
                        Dk = knnCandidates.getMaxKey();
                        DkInfinity = false;
                    } else if (knnCandidates.size() > k) {
                        // Note: knnCandidates cannot have more than k elements
                        // because each time we insert we check if size > k and 
                        // remove the largest element if yes
                        knnCandidates.extractMaxElement();
                        Dk = knnCandidates.getMaxKey();
                        //assert(!DkInfinity && "Dk cannot be infinity because k candidates have already been found");
                    }
                }
                // We have that Dk >= pqueue.getMinKey. This because the element associated with Dk must be in the queue
                // otherwise the algorithm will stop (because if it is not in pqueue then that means Dk < pqueue.getMinKey).
                // But if the element associated with Dk has been fully refined then lowerBound == Dk below, and this element
                // would not be re-insert into the queue acording to the original algorithm in the paper. Then for the next 
                // dequeued (which is the last candidate for the kNN will not be considered and we may potentially lose a 
                // kNN result). To fix this we change the < to be <= to below.             
                if (DkInfinity || lowerBound <= Dk) {
                    pqueue.insert(DataTuple(lowerBound,upperBound,next,sourceToNextDistance,minElement.point),lowerBound);
                }
            } else {
                // This indicates that an object has been refine enough that it's distance interval does
                // not intersect with any other distance interval (i.e. it is kNN result). So it is
                // implicitly dropped from priority queue.
            }
        } else if (!DkInfinity && minElement.lowerBound >= Dk) {
            // Min element is not an object and it cannot have an object that can be a kNN
            // and since it is the min element no other element can be or have an kNN
            break;
        } else {
            // Min element is a block
            if (minElement.region->isLeafNode()) {
                // Min element is a leaf block
                SimpleQuadtreeNode *leafRegion = minElement.region;
                for (std::size_t i = 0; i < leafRegion->nodes.size(); ++i) {
                    object = leafRegion->nodes[i];
                    next = minElement.intermediate; // This will be the queryNode
                    sourceToNextDistance = 0; // next is source, so source to source dist is 0
                    // We provide the source (queryNode) as the next node in the path from the source to the
                    // to the object and the lowerbound and upperbound inherited from the intersecting region
                    // so as far as refineDistance is concerned, all input values are valid and it will update
                    // next and sourceToNextDistance to be the first edge and weight in the path to the object
                    // unless the target is the query node (i.e. object == next) in which we won't refine 
                    // further anyway.
                    if (!DkInfinity) {
                        // Since calculating the Euclidean distance is somewhat expensive
                        // we only bother doing this when Dk has already been found
                        // Otherwise it is not helpful as we can't use it to prune anything
                        lowerBound = std::floor(graph.getEuclideanDistance(next,object));
                        if (lowerBound < minElement.lowerBound) {
                            lowerBound = minElement.lowerBound; // We inherit the dist bounds to the intersecting region containing element
                        }
                        // Note: The Euclidean distance can be more accurate than the minElement.lowerBound
                        // e.g. when there is only one block in the quadtree (i.e., for single-edge query nodes)
                    } else {
                        lowerBound = minElement.lowerBound;
                    }
                    
                    // Note: Unlike the original algorithm, we avoid refineDistance as long possible
                    // by using the Euclidean distance or the lower bound inherited from the leaf block
                    // (whichever is larger) to prune objects that cannot be kNN

                    if (DkInfinity || lowerBound < Dk) {
                        // Unlike the original algorithm we also compute the upperBound to regions
                        // so we can also inherit it for this object
                        upperBound = minElement.upperBound; // Must intialise this to max
#if defined(NO_CHAIN_SKIPPING)
                        this->refineDistance(graph,object,this->mortonCodes[object],next,sourceToNextDistance,lowerBound,upperBound);
#else
                        this->optimisedRefineDistance(graph,object,this->mortonCodes[object],next,sourceToNextDistance,lowerBound,upperBound,junc);
#endif
                        if (DkInfinity || lowerBound < Dk) {
#if defined(COLLECT_STATISTICS)
                            this->stats.incrementStatistic("objects_inserted",1);
#endif
                            // We check lowerbound again as refinement may have improved it
                            pqueue.insert(DataTuple(lowerBound,upperBound,next,sourceToNextDistance,object),lowerBound);
                            if (DkInfinity || upperBound < Dk) {
                                knnCandidates.insert(object,upperBound);
                                if (knnCandidates.size() == k) {
                                    Dk = knnCandidates.getMaxKey();
                                    DkInfinity = false;
                                } else if (knnCandidates.size() > k) {
                                    // Note: knnCandidates cannot have more than k elements
                                    // because each time we insert we check if size > k and 
                                    // remove the largest element if yes
                                    knnCandidates.extractMaxElement();
                                    Dk = knnCandidates.getMaxKey();
                                    //assert(!DkInfinity && "Dk cannot be infinity because k candidates have already been found");
                                }
                            }
                        }
                    }
                }
            } else {
                // Min element is a non-leaf block
                if (minElement.region->sw->getNumNodes() > 0) {
                    // We also do not want to add any empty blocks which can exist since object 
                    // hierarchy is not a linear Quadtree (obviously it will not have children)
                    this->intervalDistance(minElement.intermediate,objectHierarchy,minElement.region->sw,lowerBound,upperBound);
                    if (DkInfinity || lowerBound < Dk) {
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("regions_inserted",1);
#endif
                        pqueue.insert(DataTuple(lowerBound,upperBound,minElement.intermediate,minElement.distance,minElement.region->sw),lowerBound);
                        if (minElement.region->sw->getNumNodes() >= iK && (DkInfinity || upperBound < Dk)) {
                            // If the region contains more than k objects then we can get a estimate on the upper
                            // bound on the kth candidate. If this improves the current Dk then we can use it
                            // to potentially prune other regions using this region
                            Dk = upperBound;
                        }
                    }
                    // Note: In case of queue elements that are regions minElement.intermediate 
                    // will always be equal to junctionNode and minElement.distance will be 0
                }
                if (minElement.region->se->getNumNodes() > 0) {
                    this->intervalDistance(minElement.intermediate,objectHierarchy,minElement.region->se,lowerBound,upperBound);
                    if (DkInfinity || lowerBound < Dk) {
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("regions_inserted",1);
#endif
                        pqueue.insert(DataTuple(lowerBound,upperBound,minElement.intermediate,0,minElement.region->se),lowerBound);
                        if (minElement.region->se->getNumNodes() >= iK && (DkInfinity || upperBound < Dk)) {
                            Dk = upperBound;
                        }
                    }
                }
                if (minElement.region->nw->getNumNodes() > 0) {
                    this->intervalDistance(minElement.intermediate,objectHierarchy,minElement.region->nw,lowerBound,upperBound);
                    if (DkInfinity || lowerBound < Dk) {
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("regions_inserted",1);
#endif
                        pqueue.insert(DataTuple(lowerBound,upperBound,minElement.intermediate,0,minElement.region->nw),lowerBound);
                        if (minElement.region->nw->getNumNodes() >= iK && (DkInfinity || upperBound < Dk)) {
                            Dk = upperBound;
                        }
                    }
                }
                if (minElement.region->ne->getNumNodes() > 0) {
                    this->intervalDistance(minElement.intermediate,objectHierarchy,minElement.region->ne,lowerBound,upperBound);
                    if (DkInfinity || lowerBound < Dk) {
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("regions_inserted",1);
#endif
                        pqueue.insert(DataTuple(lowerBound,upperBound,minElement.intermediate,0,minElement.region->ne),lowerBound);
                        if (minElement.region->ne->getNumNodes() >= iK && (DkInfinity || upperBound < Dk)) {
                            Dk = upperBound;
                        }
                    }
                }
            }
        }
    }
    
    // Note: It's possible that knnCandidates size is less than k if there are
    // less than k objects in the whole road network
    knnCandidates.populateKNNs(kNNs,kNNDistances);
}

void SILCPathOracle::getKNNsByDistanceBrowsingViaRtree(StaticRtree& rtree, unsigned int k, NodeID queryNode, std::vector<NodeID>& kNNs, 
                                                       std::vector<EdgeWeight>& kNNDistances, Graph_RoadNetwork& graph, Junction& junc)
{
#if defined(COLLECT_STATISTICS)
    this->stats.clear();
    this->stats.initialiseStatistic("total_regions",0);
    this->stats.initialiseStatistic("regions_inserted",0);
    this->stats.initialiseStatistic("objects_inserted",0);
    this->stats.initialiseStatistic("total_refinements",0);
#endif
    std::vector<NodeID> euclideanKNNs;
    std::vector<EuclideanDistanceSqrd> euclideanKNNDistances;
    Coordinate queryNodeX, queryNodeY;
    graph.getCoordinates(queryNode,queryNodeX,queryNodeY);
    // Note: We retrieve the k+1 Euclidean nearest neighbours as we need at least 2 elements in queue
    BinaryMinHeap<EuclideanDistanceSqrd,RtreeDataTuple> heap = rtree.getKNNs(k,queryNodeX,queryNodeY,euclideanKNNs,euclideanKNNDistances);
    // Note: We keep heap so that we incrementally retrieve further euclidean NNs

    DistanceBound lowerBound, upperBound, Dk = 0, minRtreeDist;
    BinaryMinHeap<DistanceBound,DataTuple> pqueue;
    NodeID next, nextEuclidNN;
    EuclideanDistanceSqrd currEuclidDistSqrd;
    EdgeWeight sourceToNextDistance;
    
    BinaryMaxHeapWithDK<EdgeWeight,NodeID> knnCandidates;
    for (std::size_t i = 0; i < euclideanKNNs.size(); ++i) {
        next = queryNode; // This will be the queryNode
        sourceToNextDistance = 0; // next is source, so source to source dist is 0
        lowerBound = 0;
        upperBound = std::numeric_limits<DistanceBound>::max();
#if defined(NO_CHAIN_SKIPPING)
        this->refineDistance(graph,euclideanKNNs[i],this->mortonCodes[euclideanKNNs[i]],next,sourceToNextDistance,lowerBound,upperBound);
#else
        this->optimisedRefineDistance(graph,euclideanKNNs[i],this->mortonCodes[euclideanKNNs[i]],next,sourceToNextDistance,lowerBound,upperBound,junc);
#endif
        knnCandidates.insert(euclideanKNNs[i],upperBound);
        pqueue.insert(DataTuple(lowerBound,upperBound,next,sourceToNextDistance,euclideanKNNs[i]),lowerBound);
#if defined(COLLECT_STATISTICS)
        this->stats.incrementStatistic("objects_inserted",1);
#endif
    }

    Dk = knnCandidates.getMaxKey();
    // Note: It's possible that there are less than k objects in the whole set 
    // in which we don't need to use Dk to prune, so setting it here is safe
    
    // Insert the rest of the results in Rtree as a pseudo-block using a dummy pointer
    NodeID specialNode = std::numeric_limits<NodeID>::max();
    if (heap.size() > 0) {
        minRtreeDist = std::floor(std::sqrt(heap.getMinKey()));
        pqueue.insert(DataTuple(minRtreeDist,0,queryNode,0,specialNode),minRtreeDist);
    }
    
    while(pqueue.size() > 0) {
        DataTuple minElement = pqueue.extractMinElement();
        if (minElement.lowerBound >= Dk) {
            break;
        } else if (minElement.point != specialNode) {
            if (minElement.lowerBound >= Dk) {
                break;
            } else if (minElement.upperBound > pqueue.getMinKey() || (minElement.upperBound == pqueue.getMinKey() && minElement.lowerBound != minElement.upperBound)) {
                // The >= condition in the above statement (accordining to the original paper) has been changed 
                // to a > condition (and the additional condition handles an edge-case which can cause misorderings)
                // because it potentially leads to an infinite loop. This happens because if minElement.upperBound
                // == pqueue.getMinKey() then we will try refine it, even if it has been fully refined.
                // As long as the lowerbound <= Dk we will re-insert into pqueue. Since it is fully 
                // refined, it will still be the minmum element in pqueue and we have an infinite loop.
                // The second condition is added because we need to fully refine two objects with the
                // same upperbound otherwise they will remain in kNNCandidates in the wrong order.

                // Note: We don't need to delete element as in original paper because we decrease key instead 
//                 if (minElement.upperBound <= Dk) {
//                     // We try to remove element, if upperbound is smaller than or equal to Dk
//                     // then it could already in knnCandidates
//                     if (knnCandidates.contains(minElement.point)) {
//                         knnCandidates.deleteElement(minElement.point);
//                     }
//                 }
                
                // Note: If we reach here then minElement.intermediate cannot equal minElement.point
                // because this can only happen when an object has been fully refined. In that case
                // minElement.lowerBound == minElement.upperBound (see refineDistance), so since it's
                // the minElement in the queue we must also have minElement.upperBound <= pqueue.getMinKey()
                // and therefore it would not reach here. Note that in the case of the second part of the 
                // condition it excludes element that have been fully refined (would equal Dk in that case).
                next = minElement.intermediate;
                sourceToNextDistance = minElement.distance;
                lowerBound = minElement.lowerBound;
                upperBound = minElement.upperBound;
#if defined(NO_CHAIN_SKIPPING)
                this->refineDistance(graph,minElement.point,this->mortonCodes[minElement.point],next,sourceToNextDistance,lowerBound,upperBound);
#else
                this->optimisedRefineDistance(graph,minElement.point,this->mortonCodes[minElement.point],next,sourceToNextDistance,lowerBound,upperBound,junc);
#endif
                // Note: refineDistance will update next and sourceToNextDistance to the next intermediate node
                // and it will also retrieve loweBound and upperBound for the current intermediate node (next)

                if (upperBound <= Dk) {
                    if (knnCandidates.contains(minElement.point)) {
                        // The upperBound, at worst, remains the same. It only be improved
                        // by refinement so, we know this is a decrease key operation. Since
                        // we forced to use heap that supports decrease (to avoid re-inserting
                        // the same object twice0 we use it here to decrease instead of delete
                        // and re-insert as in the original algorithm
                        knnCandidates.decreaseKey(minElement.point,upperBound);
                    } else {
                        knnCandidates.insert(minElement.point,upperBound);
                    }
                    if (knnCandidates.size() == k) {
                        Dk = knnCandidates.getMaxKey();
                    } else if (knnCandidates.size() > k) {
                        // Note: knnCandidates cannot have more than k elements
                        // because each time we insert we check if size > k and 
                        // remove the largest element if yes
                        knnCandidates.extractMaxElement();
                        Dk = knnCandidates.getMaxKey();
                    }
                }
                
                // We have that Dk >= pqueue.getMinKey. This because the element associated with Dk must be in the queue
                // otherwise the algorithm will stop (because if it is not in pqueue then that means Dk < pqueue.getMinKey).
                // But if the element associated with Dk has been fully refined then lowerBound == Dk below, and this element
                // would not be re-insert into the queue acording to the original algorithm in the paper. Then for the next 
                // dequeued (which is the last candidate for the kNN) will not be considered and we may potentially lose a 
                // kNN result. To fix this we change the < to be <= to below.
                if (lowerBound <= Dk) {
                    pqueue.insert(DataTuple(lowerBound,upperBound,next,sourceToNextDistance,minElement.point),lowerBound);
                }
            }
        } else {
            // Note: If the heap size is zero we don't bother doing the next NN search
            if (heap.size() > 0 && rtree.getNextNearestNeighbour(heap,queryNodeX,queryNodeY,nextEuclidNN,currEuclidDistSqrd)) {
                // Note: If there were less than k objects to begin with then getNextNearestNeighbour
                // would return false (as the heap is empty) and we not reach here
                next = queryNode; // This will be the queryNode
                sourceToNextDistance = 0; // next is source, so source to source dist is 0
                lowerBound = std::floor(std::sqrt(currEuclidDistSqrd));
                if (lowerBound < Dk) {
                    // We use the Euclidean distance from the source to the target to decide
                    // whether it is worth refining the distance or not
                    upperBound = std::numeric_limits<DistanceBound>::max();
#if defined(NO_CHAIN_SKIPPING)
                    this->refineDistance(graph,nextEuclidNN,this->mortonCodes[nextEuclidNN],next,sourceToNextDistance,lowerBound,upperBound);
#else
                    this->optimisedRefineDistance(graph,nextEuclidNN,this->mortonCodes[nextEuclidNN],next,sourceToNextDistance,lowerBound,upperBound,junc);
#endif
                    if (lowerBound < Dk) {
                        pqueue.insert(DataTuple(lowerBound,upperBound,next,sourceToNextDistance,nextEuclidNN),lowerBound);
#if defined(COLLECT_STATISTICS)
                        this->stats.incrementStatistic("objects_inserted",1);
#endif
                        if (upperBound < Dk) {
                            knnCandidates.insert(nextEuclidNN,upperBound);
                            if (knnCandidates.size() == k) {
                                Dk = knnCandidates.getMaxKey();
                            } else if (knnCandidates.size() > k) {
                                // Note: knnCandidates cannot have more than k elements
                                // because each time we insert we check if size > k and 
                                // remove the largest element if yes
                                knnCandidates.extractMaxElement();
                                Dk = knnCandidates.getMaxKey();
                            }
                        }
                    }
                }
                if (heap.size() > 0) {
                    // If there anymore potential results then we re-insert pseudoblock
                    minRtreeDist = std::floor(std::sqrt(heap.getMinKey()));
                    pqueue.insert(DataTuple(minRtreeDist,0,queryNode,0,specialNode),minRtreeDist);
                }
            }
        }
    }

    // Note: It's possible that knnCandidates size is less than k if there are
    // less than k objects in the whole road network
    knnCandidates.populateKNNs(kNNs,kNNDistances);
}

MortonCode SILCPathOracle::convertToMortonCode(Graph_RoadNetwork& graph, NodeID node)
{
    Coordinate x, y;
    graph.getCoordinates(node,x,y);
    // We translate all X and Y coordinates such that (0,0) is the most
    // SW point and coordinates increase in the N and E directions
    x = x - this->xTranslation; // Translate right, e.g. -100 -(-100)
    y = y - this->yTranslation; 
    // Note: In order to have (0,0) be in the SW corner we adjust the X 
    // and Y coordinates. If coordinate is negative translate right
    // and if coordinate is position when translate left.
    MortonCode targetCode(x,y,this->maxRangeExp);
    return targetCode;
}


void SILCPathOracle::intervalDistance(NodeID source, SimpleQuadtree& qt, SimpleQuadtreeNode* region, 
                                      DistanceBound& lowerBound, DistanceBound& upperBound)
{
    // The region's morton code will be relative to a different origin, so we translate it
    // This is possible because the underlying space is partitioned by the same units.
    // The minimum width of 1 representing a point is the same for both representation 
    // because we are choosing objects from the same set of points. Origin (0,0) in either 
    // quadtree may not be registration however, so there will be a translation.
    
    Coordinate regionX, regionY;
    int origWidth;
    region->getBlockCoordinates(regionX,regionY,origWidth);
    regionX -= this->xTranslation;
    regionY -= this->yTranslation;

    // We also need compute the relative coordinates of the source since we moved origin
    // when building Morton lists so that all coordinates began at (0,0) and increased
    this->pathOracle[source].getIntervalDistance(this->relativeXCoordinates[source],this->relativeYCoordinates[source],regionX,regionY,origWidth,this->maxRangeExp,this->maxRange,lowerBound,upperBound);
}

void SILCPathOracle::intervalDistanceForDistBrws(NodeID source, SimpleQuadtree& qt, SimpleQuadtreeNode* region, 
                                                 DistanceBound& lowerBound)
{
    // The region's morton code will be relative to a different origin, so we translate it
    // This is possible because the underlying space is partitioned by the same units.
    // The minimum width of 1 representing a point is the same for both representation 
    // because we are choosing objects from the same set of points. Origin (0,0) in either 
    // quadtree may not be registration however, so there will be a translation.
    
    Coordinate regionX, regionY;
    int origWidth;
    region->getBlockCoordinates(regionX,regionY,origWidth);
    regionX -= this->xTranslation;
    regionY -= this->yTranslation;

    // We also need compute the relative coordinates of the source since we moved origin
    // when building Morton lists so that all coordinates began at (0,0) and increased
    this->pathOracle[source].getIntervalDistanceForDistBrws(this->relativeXCoordinates[source],this->relativeYCoordinates[source],regionX,regionY,origWidth,this->maxRangeExp,this->maxRange,lowerBound);
}

void SILCPathOracle::findRealKNNDistances(Graph_RoadNetwork& graph, NodeID queryNode, std::vector<NodeID>& kNNs, std::vector<EdgeWeight>& kNNDistances)
{
    for (std::size_t i = 0; i < kNNs.size(); ++i) {
        kNNDistances[i] = this->findShortestPathDistance(graph,queryNode,kNNs[i]);
    }
}

double SILCPathOracle::computeIndexSize()
{
    double memoryUsage = 0;
    for (std::size_t i = 0; i < this->pathOracle.size(); ++i) {
        memoryUsage += this->pathOracle[i].computeIndexSizeBytes();
    }
    memoryUsage += sizeof(Coordinate)*relativeXCoordinates.size();
    memoryUsage += sizeof(Coordinate)*relativeYCoordinates.size();
    for (std::size_t i = 0; i < this->mortonCodes.size(); ++i) {
        memoryUsage += this->mortonCodes[i].computeIndexSizeBytes();
    }    
    return memoryUsage/(1024*1024);
}

double SILCPathOracle::computeMemoryUsage() {
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    memoryUsage += this->networkName.size();
    for (std::size_t i = 0; i < this->pathOracle.size(); ++i) {
        memoryUsage += this->pathOracle[i].computeMemoryUsageBytes();
    }
    memoryUsage += sizeof(MortonList)*(this->pathOracle.capacity()-this->pathOracle.size());
    memoryUsage += sizeof(Coordinate)*relativeXCoordinates.capacity();
    memoryUsage += sizeof(Coordinate)*relativeYCoordinates.capacity();
    for (std::size_t i = 0; i < this->mortonCodes.size(); ++i) {
        memoryUsage += this->mortonCodes[i].computeMemoryUsageBytes();
    }
    memoryUsage += sizeof(MortonCode)*(this->mortonCodes.capacity()-this->mortonCodes.size());
    return memoryUsage/(1024*1024);
}

/*
 * MortonList
 */

EdgeID MortonList::getLink(const MortonCode& targetCode) {
    // Binary search morton list to find block containing target
    int listIndex = this->getBlockIdxContainingTargetByBST(targetCode);
    return this->binarySearchTreeMortonList[listIndex].getLink();
}

EdgeID MortonList::getLinkAndRatios(const MortonCode& targetCode, DistanceRatio& minRatio, DistanceRatio& maxRatio) {
    // Binary search morton list to find block containing target
    int listIndex = this->getBlockIdxContainingTargetByBST(targetCode);
    return this->binarySearchTreeMortonList[listIndex].getLinkAndRatios(minRatio,maxRatio);
}

std::size_t MortonList::getBlockIdxContainingTarget(const MortonCode& targetCode) {
    bool found = false;
    int left = 0, right = this->mortonList.size()-1, mid, status;
    while(left <= right) {
        mid = (right+left)/2;
        status = this->mortonList[mid].compareMortonCode(targetCode);
        if (status == 0) {
            // The target is contained within this morton block
            found = true;
            break;
        } else if (status == 1) {
            // Target is larger than middle
            left = mid+1;
        } else /* if (status == = -1)*/ {
            // Target is smaller than middle
            right = mid-1;
        }
    }
    if (found) {
        return mid;
    } else {
        std::cerr << "getBlockIdxContainingTarget() could not find containing block of target morton code in morton list" << std::endl;
        exit(1);
    }
}

std::size_t MortonList::getBlockIdxContainingTargetByBST(const MortonCode& targetCode) {
    bool found = false;
    std::size_t parent = 0; // Start at root element
    int status;
    while(parent < this->binarySearchTreeMortonList.size()) {
        status = this->binarySearchTreeMortonList[parent].compareMortonCode(targetCode);
        if (status == 0) {
            // The target is contained within this morton block
            found = true;
            break;
        } else if (status == 1) {
            // Target is larger so it is in the right branch
            parent = this->getRightChild(parent);
            //assert (parent < this->binarySearchTreeMortonList.size() && "Fell off tree from right");
        } else /* if (status == = -1)*/ {
            // Target is smaller is it is in the left branch
            parent = this->getLeftChild(parent);
            //assert (parent < this->binarySearchTreeMortonList.size() && "Fell off tree left");
        }
    }
    if (found) {
        return parent;
    } else {
        std::cerr << "getBlockIdxContainingTargetByBST() could not find containing block of target morton code in morton list" << std::endl;
        exit(1);
    }
}

void MortonList::getIntervalDistance(Coordinate x, Coordinate y, Coordinate regionLeft, Coordinate regionBottom, int regionWidth, 
                                     int maxRangeExp, int maxRange, DistanceBound& lowerBound, DistanceBound& upperBound)
{
    // Store the coordinates of the corners of the region we are checking
    Coordinate regionRight, regionTop;
    geometry::getRightTopCorner(regionLeft,regionBottom,regionWidth, regionRight, regionTop);

    bool initialized = false;
    
    // Limit range to start/end blocks
    
    // Since MortonCode can only fit about maxRange numbers so we need to track
    // case where max in range is essentially infinity separately
    Coordinate containedRegionRight = regionRight-1, containedRegionTop = regionTop-1;
    if (regionRight > maxRange) {
        // Recall that maxRange is from 0 to 2^maxRangeExp-1 so > is OK
        // i.e. maxRange is a possible value
        containedRegionRight = maxRange;
    }
    if (regionTop > maxRange) {
        containedRegionTop = maxRange;
    }

    // Note: The NE point of the block would not be contained in any 
    // We want the NE point of the region that would still be contained
    // within a Morton block (decrement by 1 for the point width of 1)
    MortonCode regionSW(regionLeft,regionBottom,maxRangeExp); // Most SW point in region
    MortonCode regionNE(containedRegionRight,containedRegionTop,maxRangeExp);
    
    this->rangeSearch(0,regionSW,regionNE,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,lowerBound,upperBound,initialized);
    if (!initialized) {
        std::cerr << "No intersection for object hierarchy region was found in getIntervalDistance()" << std::endl;
    }
//     assert(initialized && "No intersection for object hierarchy region was found");
    
//     // Check all morton block for intersection
//     
//     for (std::size_t i = 0; i < binarySearchTreeMortonList.size(); ++i) {
//         this->unionBoundsWithCandidateBlock(i,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,lowerBound,upperBound,initialized);
//     }
//     if (!initialized) {
//         std::cerr << "No intersection for object hierarchy region was found in getIntervalDistance()" << std::endl;
//     }
// //     assert(initialized && "No intersection for object hierarchy region was found");

}

bool MortonList::unionBoundsWithCandidateBlock(int blockIdx, int maxRangeExp, Coordinate x, Coordinate y, 
                                               Coordinate regionLeft, Coordinate regionBottom, Coordinate regionRight, Coordinate regionTop,
                                               DistanceBound& minLowerBound, DistanceBound& maxUpperBound, bool& initialized)
{
    // This method assumes that we don't know if provided block is intersecting (just a candidate)
    
    DistanceBound newLowerBound, newUpperBound;
    Coordinate blockLeft, blockBottom, blockRight, blockTop;
    int blockWidth;
    
    // Store the coordinates of the overlapping rectangle between them
    Coordinate overlapLeft, overlapBottom, overlapRight, overlapTop;
    DistanceRatio min, max;
    EuclideanDistance minDist, maxDist;
    
    this->binarySearchTreeMortonList[blockIdx].decomposeAndGetRatios(maxRangeExp,blockLeft,blockBottom,blockWidth,min,max);
    geometry::getRightTopCorner(blockLeft,blockBottom,blockWidth,blockRight,blockTop);
    
    bool isOverlapping = geometry::getIntersectingRectangle(blockLeft, blockBottom, blockRight, blockTop, 
                                           regionLeft, regionBottom, regionRight, regionTop, 
                                           overlapLeft, overlapBottom, overlapRight, overlapTop);
    
    // If overlapping we find the intersection rectangle then
    if (isOverlapping) {
        minDist = geometry::getMinDistToRectangle(x, y, overlapLeft, overlapBottom, overlapRight, overlapTop);
        maxDist = geometry::getMaxDistToRectangle(x, y, overlapLeft, overlapBottom, overlapRight, overlapTop);
        newLowerBound = std::floor(min*minDist);
        newUpperBound = std::ceil(max*maxDist);
        if (!initialized) {
            initialized = true;
            minLowerBound = newLowerBound;
            maxUpperBound = newUpperBound;
        } else {
            if (minLowerBound > newLowerBound) {
                minLowerBound = newLowerBound;
            }
            if (maxUpperBound < newUpperBound) {
                maxUpperBound = newUpperBound;
            }
        }
    } /*else {
        std::cerr << "Range search found non-overlapping block" << std::endl;
        exit(1);
    }*/
    
    return isOverlapping;
}

void MortonList::rangeSearch(std::size_t parent, MortonCode& minCode, MortonCode& maxCode, int maxRangeExp, Coordinate x, Coordinate y, 
                             Coordinate regionLeft, Coordinate regionBottom, Coordinate regionRight, Coordinate regionTop, 
                             DistanceBound& minLowerBound, DistanceBound& maxUpperBound, bool& initialized)
{
    if (parent < this->binarySearchTreeMortonList.size()) {
        int status = this->binarySearchTreeMortonList[parent].compareMortonCode(minCode);
        if (status >= 0) {
            // This means minCode is larger than or equal to parent's code
            // The SW point of the region is contained or larger than parent Morton block
            // so we have found the smallest Morton block in range and hence
            // do not need to follow the left branch as they will be smaller
            // than the parent block and hence smaller than minCode
            if (status == 0) {
                // If SW point is contained, then it is obviously an intersecting block
                // so we add it to the list of candidates
                this->unionBoundsWithCandidateBlock(parent,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,minLowerBound,maxUpperBound,initialized);
            }
            this->rangeSearch(this->getRightChild(parent),minCode,maxCode,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,minLowerBound,maxUpperBound,initialized);
        } else {
            // This means minCode is smaller than parent's code and we could potentially
            // follow both left and right branches of parent
            // Note: It's OK to use else because 
            status = this->binarySearchTreeMortonList[parent].compareMortonCode(maxCode);
            if (status <= 0) {
                // This means maxCode is smaller than or equal to parent's code
                // The NE point of the region is contained or smaller than parent Morton block
                // so we have found the largest Morton block in range and hence
                // do not need to follow the right branch as they will be larger
                // than the parent block and hence larger than maxCode
                if (status == 0) {
                    // If NE point is contained, then it is obviously an intersecting block
                    // so we add it to the list of candidates
                    this->unionBoundsWithCandidateBlock(parent,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,minLowerBound,maxUpperBound,initialized);
                }
                this->rangeSearch(this->getLeftChild(parent),minCode,maxCode,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,minLowerBound,maxUpperBound,initialized);
            } else {
                // This is reached if both NE and SW point are within but not contained
                // in any Morton block, so it is a candidate to be in the range
                // and we must also traverse both branches
                this->unionBoundsWithCandidateBlock(parent,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,minLowerBound,maxUpperBound,initialized);
                this->rangeSearch(this->getLeftChild(parent),minCode,maxCode,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,minLowerBound,maxUpperBound,initialized);
                this->rangeSearch(this->getRightChild(parent),minCode,maxCode,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,minLowerBound,maxUpperBound,initialized);
            }
        }
    }
}

bool MortonList::checkLeftChildOverlapping(int parentIdx, int maxRangeExp, Coordinate regionLeft, Coordinate regionBottom, Coordinate regionRight, Coordinate regionTop)
{
    bool isOverlapping = false;
    Coordinate blockLeft, blockBottom, blockRight, blockTop;
    int blockWidth;
    DistanceRatio min, max;
    Coordinate overlapLeft, overlapBottom, overlapRight, overlapTop;
    unsigned int leftChildIdx = this->getLeftChild(parentIdx);
    if (leftChildIdx < this->binarySearchTreeMortonList.size()) {
        this->binarySearchTreeMortonList[leftChildIdx].decomposeAndGetRatios(maxRangeExp,blockLeft,blockBottom,blockWidth,min,max);
        geometry::getRightTopCorner(blockLeft,blockBottom,blockWidth,blockRight,blockTop);
        
        isOverlapping = geometry::getIntersectingRectangle(blockLeft, blockBottom, blockRight, blockTop, 
                                            regionLeft, regionBottom, regionRight, regionTop, 
                                            overlapLeft, overlapBottom, overlapRight, overlapTop);    
    }
    return isOverlapping;
}

bool MortonList::checkRightChildOverlapping(int parentIdx, int maxRangeExp, Coordinate regionLeft, Coordinate regionBottom, Coordinate regionRight, Coordinate regionTop)
{
    bool isOverlapping = false;
    Coordinate blockLeft, blockBottom, blockRight, blockTop;
    int blockWidth;
    DistanceRatio min, max;
    Coordinate overlapLeft, overlapBottom, overlapRight, overlapTop;
    unsigned int rightChildIdx = this->getRightChild(parentIdx);
    if (rightChildIdx < this->binarySearchTreeMortonList.size()) {
        this->binarySearchTreeMortonList[rightChildIdx].decomposeAndGetRatios(maxRangeExp,blockLeft,blockBottom,blockWidth,min,max);
        geometry::getRightTopCorner(blockLeft,blockBottom,blockWidth,blockRight,blockTop);
        
        isOverlapping = geometry::getIntersectingRectangle(blockLeft, blockBottom, blockRight, blockTop, 
                                            regionLeft, regionBottom, regionRight, regionTop, 
                                            overlapLeft, overlapBottom, overlapRight, overlapTop);    
    }

    return isOverlapping;
}

void MortonList::computeNewSearchBoundaries(MortonCode& minCode, MortonCode& maxCode, MortonBlock& dividingBlock, MortonNumber& leftMaxZ, MortonNumber& rightMinZ, int maxRangeExp)
{
    std::bitset<64> min(minCode.getZNumber());
    std::bitset<64> max(maxCode.getZNumber());
    std::bitset<64> div(dividingBlock.getMortonNumber());
    std::bitset<64> leftMax; // Default construction initializes with zeroes
    std::bitset<64> rightMin;
    
    bool finished = false;
    for (int i = quadtree_constants::MAX_MORTON_LENGTH-1; i >= quadtree_constants::MAX_MORTON_LENGTH-2*maxRangeExp && !finished; --i) {
        /*if (div[i] == 0 && min[i] == 0 && max[i] == 0) {
            // No action, continue scanning
        } else */ if (div[i] == 0 && min[i] == 0 && max[i] == 1){
            rightMin = min;
            this->loadPattern(1,0,i,rightMin,maxRangeExp);
            this->loadPattern(0,1,i,max,maxRangeExp);
        } /*else if (div[i] == 0 && min[i] == 1 && max[i] == 0){
            // This case is not possible because min <= max
        } */else if (div[i] == 0 && min[i] == 1 && max[i] == 1){
            rightMinZ = min.to_ulong();
            leftMaxZ = leftMax.to_ulong();
            finished = true;
        } else if (div[i] == 1 && min[i] == 0 && max[i] == 0){
            leftMaxZ = max.to_ulong();
            rightMinZ = rightMin.to_ulong();
            finished = true;
        } else if (div[i] == 1 && min[i] == 0 && max[i] == 1){
            leftMax = max;
            this->loadPattern(0,1,i,leftMax,maxRangeExp);
            this->loadPattern(1,0,i,min,maxRangeExp);
        } /*else if (div[i] == 1 && min[i] == 1 && max[i] == 0){
            // This case is not possible because min <= max
        } else if (div[i] == 1 && min[i] == 1 && max[i] == 1){
            // No action, continue scanning
        }*/
    }
}

void MortonList::loadPattern(unsigned int firstVal, unsigned int repeatVal, int pos, std::bitset<64>& zNumber, int maxRangeExp)
{
    //assert(repeatVal == 0 || repeatVal == 1);
    //assert(firstVal == 0 || firstVal == 1);
    //assert(pos >= 0);
    zNumber[pos] = firstVal;
    for (int i = pos-2; i >= quadtree_constants::MAX_MORTON_LENGTH-2*maxRangeExp;) {
        // Note: We only set every second bit because we are making changes for
        // one dimension only (whether x or why depends on pos)
        zNumber[i] = repeatVal;
        i -= 2;
    }
}

bool MortonList::unionBoundsWithCandidateBlockForDistBrws(int blockIdx, int maxRangeExp, Coordinate x, Coordinate y, 
                                               Coordinate regionLeft, Coordinate regionBottom, Coordinate regionRight, Coordinate regionTop,
                                               DistanceBound& minLowerBound, bool& initialized)
{
    // This method assumes that we don't know if provided block is intersecting (just a candidate)
    
    DistanceBound newLowerBound;
    Coordinate blockLeft, blockBottom, blockRight, blockTop;
    int blockWidth;
    
    // Store the coordinates of the overlapping rectangle between them
    Coordinate overlapLeft, overlapBottom, overlapRight, overlapTop;
    DistanceRatio min, max;
    EuclideanDistance minDist;
    
    this->binarySearchTreeMortonList[blockIdx].decomposeAndGetRatios(maxRangeExp,blockLeft,blockBottom,blockWidth,min,max);
    geometry::getRightTopCorner(blockLeft,blockBottom,blockWidth,blockRight,blockTop);
    
    bool isOverlapping = geometry::getIntersectingRectangle(blockLeft, blockBottom, blockRight, blockTop, 
                                           regionLeft, regionBottom, regionRight, regionTop, 
                                           overlapLeft, overlapBottom, overlapRight, overlapTop);
    
    // If overlapping we find the intersection rectangle then
    if (isOverlapping) {
        minDist = geometry::getMinDistToRectangle(x, y, overlapLeft, overlapBottom, overlapRight, overlapTop);
        newLowerBound = std::floor(min*minDist);
        if (!initialized) {
            initialized = true;
            minLowerBound = newLowerBound;
        } else {
            if (minLowerBound > newLowerBound) {
                minLowerBound = newLowerBound;
            }
        }
    } /*else {
        std::cerr << "Range search found non-overlapping block" << std::endl;
        exit(1);
    }*/
    
    return isOverlapping;
}

void MortonList::rangeSearchForDistBrws(std::size_t parent, MortonCode& minCode, MortonCode& maxCode, int maxRangeExp, Coordinate x, Coordinate y, 
                             Coordinate regionLeft, Coordinate regionBottom, Coordinate regionRight, Coordinate regionTop, 
                             DistanceBound& minLowerBound, bool& initialized)
{
    if (parent < this->binarySearchTreeMortonList.size()) {
        int status = this->binarySearchTreeMortonList[parent].compareMortonCode(minCode);
        if (status >= 0) {
            // This means minCode is larger than or equal to parent's code
            // The SW point of the region is contained or larger than parent Morton block
            // so we have found the smallest Morton block in range and hence
            // do not need to follow the left branch as they will be smaller
            // than the parent block and hence smaller than minCode
            if (status == 0) {
                // If SW point is contained, then it is obviously an intersecting block
                // so we add it to the list of candidates
                this->unionBoundsWithCandidateBlockForDistBrws(parent,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,minLowerBound,initialized);
            }
            this->rangeSearchForDistBrws(this->getRightChild(parent),minCode,maxCode,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,minLowerBound,initialized);
        } else {
            // This means minCode is smaller than parent's code and we could potentially
            // follow both left and right branches of parent
            // Note: It's OK to use else because 
            status = this->binarySearchTreeMortonList[parent].compareMortonCode(maxCode);
            if (status <= 0) {
                // This means maxCode is smaller than or equal to parent's code
                // The NE point of the region is contained or smaller than parent Morton block
                // so we have found the largest Morton block in range and hence
                // do not need to follow the right branch as they will be larger
                // than the parent block and hence larger than maxCode
                if (status == 0) {
                    // If NE point is contained, then it is obviously an intersecting block
                    // so we add it to the list of candidates
                    this->unionBoundsWithCandidateBlockForDistBrws(parent,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,minLowerBound,initialized);
                }
                this->rangeSearchForDistBrws(this->getLeftChild(parent),minCode,maxCode,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,minLowerBound,initialized);
            } else {
                // This is reached if both NE and SW point are within but not contained
                // in any Morton block, so it is a candidate to be in the range
                // and we must also traverse both branches
                this->unionBoundsWithCandidateBlockForDistBrws(parent,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,minLowerBound,initialized);
                this->rangeSearchForDistBrws(this->getLeftChild(parent),minCode,maxCode,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,minLowerBound,initialized);
                this->rangeSearchForDistBrws(this->getRightChild(parent),minCode,maxCode,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,minLowerBound,initialized);
            }
        }
    }
}

void MortonList::getIntervalDistanceForDistBrws(Coordinate x, Coordinate y, Coordinate regionLeft, Coordinate regionBottom, int regionWidth, 
                                     int maxRangeExp, int maxRange, DistanceBound& lowerBound)
{
    // Store the coordinates of the corners of the region we are checking
    Coordinate regionRight, regionTop;
    geometry::getRightTopCorner(regionLeft,regionBottom,regionWidth, regionRight, regionTop);

    // Limit range to start/end blocks
    
    // Since MortonCode can only fit about maxRange numbers so we need to track
    // case where max in range is essentially infinity separately
    Coordinate containedRegionRight = regionRight-1, containedRegionTop = regionTop-1;
    if (regionRight > maxRange) {
        // Recall that maxRange is from 0 to 2^maxRangeExp-1 so > is OK
        // i.e. maxRange is a possible value
        containedRegionRight = maxRange;
    }
    if (regionTop > maxRange) {
        containedRegionTop = maxRange;
    }
    // Note: The NE point of the block would not be contained in any 
    // We want the NE point of the region that would still be contained
    // within a Morton block (decrement by 1 for the point width of 1)
    MortonCode regionSW(regionLeft,regionBottom,maxRangeExp); // Most SW point in region
    MortonCode regionNE(containedRegionRight,containedRegionTop,maxRangeExp);
    
    bool initialized = false;
    this->rangeSearchForDistBrws(0,regionSW,regionNE,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,lowerBound,initialized);
    if (!initialized) {
        std::cerr << "No intersection for object hierarchy region was found in getIntervalDistanceForDistBrws()" << std::endl;
    }
//     assert(initialized && "No intersection for object hierarchy region was found");

//     for (std::size_t i = 0; i < binarySearchTreeMortonList.size(); ++i) {
//         this->unionBoundsWithCandidateBlockForDistBrws(i,maxRangeExp,x,y,regionLeft,regionBottom,regionRight,regionTop,lowerBound,initialized);
//     }
//     if (!initialized) {
//         std::cerr << "No intersection for object hierarchy region was found in getIntervalDistanceForDistBrws()" << std::endl;
//     }

}

std::size_t MortonList::getLeftChild(std::size_t parent)
{
    return parent*2 + 1;
}

std::size_t MortonList::getRightChild(std::size_t parent)
{
    return parent*2 + 2;
}

bool MortonList::isTreeSorted(int pos)
{
    bool sorted = true;
    std::size_t leftChild = this->getLeftChild(pos);
    MortonCode parentCode, leftCode, rightCode;
    if (leftChild < this->binarySearchTreeMortonList.size()) {
        parentCode = this->binarySearchTreeMortonList[pos].getMortonCode();
        leftCode = this->binarySearchTreeMortonList[leftChild].getMortonCode();
        if (leftCode > parentCode) {
            std::cout << "Left child at " << leftChild << " greater than parent at " << pos << std::endl;
            sorted = false;
        } else {
            sorted = this->isTreeSorted(leftChild);
        }
        // We don't have to check right child if left child is bigger than vector
        // because left child always appears before right child
        std::size_t rightChild = this->getRightChild(pos);
        if (rightChild < this->binarySearchTreeMortonList.size()) {
            rightCode = this->binarySearchTreeMortonList[rightChild].getMortonCode();
            if (rightCode < parentCode) {
                std::cout << "Right child at " << rightChild << " smaller than parent at " << pos << std::endl;
                sorted = false;
            } else {
                sorted = this->isTreeSorted(rightChild);
            }
        } 
    }
    return sorted;
}

// Check left and right subtree confirm to binary search tree property
// i.e. all nodes on left are small and all nodes on right are larger
bool MortonList::areSubtreesValid(int pos)
{
    bool isValid = true;
    MortonCode parentCode, leftCode, rightCode;
    std::size_t leftChild = this->getLeftChild(pos);
    if (leftChild < this->binarySearchTreeMortonList.size()) {
        parentCode = this->binarySearchTreeMortonList[pos].getMortonCode();
        leftCode = this->binarySearchTreeMortonList[leftChild].getMortonCode();
        if (leftCode > parentCode) {
            std::cout << "Left child at " << leftChild << " greater than parent at " << pos << std::endl;
            isValid = false;
        } else {
            isValid = this->isLeftSubtreeValid(leftChild,parentCode);
        }
        // We don't have to check right child if left child is bigger than vector
        // because left child always appears before right child
        std::size_t rightChild = this->getRightChild(pos);
        if (rightChild < this->binarySearchTreeMortonList.size()) {
            rightCode = this->binarySearchTreeMortonList[rightChild].getMortonCode();
            if (rightCode < parentCode) {
                std::cout << "Right child at " << rightChild << " smaller than parent at " << pos << std::endl;
                isValid = false;
            } else {
                isValid = this->isRightSubtreeValid(rightChild,parentCode);
            }
        } 
    }
    return isValid;
}

// All nodes on the left subtree must be smaller than parent
bool MortonList::isLeftSubtreeValid(int subtreeRoot, MortonCode& parentCode)
{
    bool isValid = true;
    MortonCode leftCode, rightCode;
    std::size_t leftChild = this->getLeftChild(subtreeRoot);
    if (leftChild < this->binarySearchTreeMortonList.size()) {
        leftCode = this->binarySearchTreeMortonList[leftChild].getMortonCode();
        if (leftCode > parentCode) {
            std::cout << "Tree node in left subtree of root parent at " << leftChild << " greater than root parent" << std::endl;
            parentCode.printCode();
            leftCode.printCode();
            isValid = false;
        } else {
            isValid = this->isLeftSubtreeValid(leftChild,parentCode);
        }
        // We don't have to check right child if left child is bigger than vector
        // because left child always appears before right child
        std::size_t rightChild = this->getRightChild(subtreeRoot);
        if (rightChild < this->binarySearchTreeMortonList.size()) {
            rightCode = this->binarySearchTreeMortonList[rightChild].getMortonCode();
            if (rightCode > parentCode) {
                std::cout << "Tree node in left subtree of root parent at " << rightChild << " greater than root parent" << std::endl;
                parentCode.printCode();
                leftCode.printCode();
                isValid = false;
            } else {
                isValid = this->isLeftSubtreeValid(rightChild,parentCode);
            }
        } 
    }
    return isValid;
}

// All nodes on the right subtree must be larger than parent
bool MortonList::isRightSubtreeValid(int subtreeRoot, MortonCode& parentCode)
{
    bool isValid = true;
    MortonCode leftCode, rightCode;
    std::size_t leftChild = this->getLeftChild(subtreeRoot);
    if (leftChild < this->binarySearchTreeMortonList.size()) {
        leftCode = this->binarySearchTreeMortonList[leftChild].getMortonCode();
        if (leftCode < parentCode) {
            std::cout << "Tree node in right subtree of root parent at " << leftChild << " smaller than root parent" << std::endl;
            parentCode.printCode();
            leftCode.printCode();
            isValid = false;
        } else {
            isValid = this->isRightSubtreeValid(leftChild,parentCode);
        }
        // We don't have to check right child if left child is bigger than vector
        // because left child always appears before right child
        std::size_t rightChild = this->getRightChild(subtreeRoot);
        if (rightChild < this->binarySearchTreeMortonList.size()) {
            rightCode = this->binarySearchTreeMortonList[rightChild].getMortonCode();
            if (rightCode < parentCode) {
                std::cout << "Tree node in right subtree of root parent at " << rightChild << " smaller than root parent" << std::endl;
                parentCode.printCode();
                leftCode.printCode();
                isValid = false;
            } else {
                isValid = this->isRightSubtreeValid(rightChild,parentCode);
            }
        } 
    }
    return isValid;
}

void MortonList::Mortonize(NodeID source, Quadtree &qt, std::vector<EdgeID> &colourMap, std::vector<EdgeWeight>& SSSPDistances, bool singleColour) {
    
    QuadtreeNode* currNode = qt.firstLeaf;
//     assert (currNode->getMortonCode().getQuadrant() == SW && "First leaf is not SW block");
    
    EdgeID colour, mergeColour;
    bool sameColour, mergeColourSet, empty;
    
    mergeColour = quadtree_constants::EMPTY_COLOUR;
    mergeColourSet = false;
    std::vector<MortonBlock> mergeCandidates;
    std::vector<MortonBlock> tempList;
    MortonCode nextMergeCode;

    std::vector<DistanceRatio> distRatios;
    DistanceRatio currRatio;
    DistanceRatio minBlockRatio = 0, maxBlockRatio = 0;
    bool intervalInitialised = false;
    
    // We know there will be at most four merge candidates
    mergeCandidates.reserve(quadtree_constants::NUM_QUADRANTS);
    
    while (currNode != NULL) {
        // Note: If quadtree node has no nodes then it will be empty and have colour quadtree_constants::EMPTY_COLOUR
        sameColour = true;
        empty = true;
        colour = quadtree_constants::EMPTY_COLOUR;
        const std::vector<NodeID>& nodes = currNode->getNodes();
        minBlockRatio = 0;
        maxBlockRatio = 0;
        intervalInitialised = false;
        if (nodes.size() > 0) {
            // Since we are iterating over all nodes in this node take opportunity to calculate
            // spatial distances and keep track of min and max intervals
            distRatios.clear();
            
            for (std::size_t i = 0; i < nodes.size(); ++i) {
                if (nodes[i] != source) {
                    currRatio = SSSPDistances[nodes[i]]/qt.getEuclideanDistance(source,nodes[i]);
                    // Note: Since qt.getEuclideanDistance returns type EuclideanDistance (i.e. double/float)
                    // the division will be floating point division
                    //assert (currRatio >= 1 && "Network distance is smaller than euclidean distance");
                } else {
                    currRatio = 0;
                }
                
                if (nodes[i] != source || !singleColour) {
                    // If the current node is the source and there is only a single colour
                    // (i.e. the source only has one edge), then we do not want to update
                    // the rations using currRatio. In this case there will be a single block
                    // and it will have a min and max ratio of 0, and SILC-based kNN methods
                    // will not be able prune Object Hierarchy region using these values.
                    // However if there is more than one colour, then we need to set the 
                    // ratio of the block containing the source to zero, to ensure that the
                    // block covering the source is retrieved first (to ensure correctness
                    // in the case that the source is an object).
                    if (!intervalInitialised) { // Can use i == 0 here instead (previously we couldn't because source was considered empty)
                        minBlockRatio = currRatio;
                        maxBlockRatio = currRatio;
                        intervalInitialised = true;
                    } else {
                        if (currRatio < minBlockRatio) {
                            minBlockRatio = currRatio;
                        }
                        if (currRatio > maxBlockRatio) {
                            maxBlockRatio = currRatio;
                        }
                    }
                    distRatios.push_back(currRatio);
                }
                
                if (empty) {
                    empty = false;
                    colour = colourMap[nodes[i]];
                    // Note: This mean we will also add source to the Morton list
                } else if (colourMap[nodes[i]] != colour && colourMap[nodes[i]] != quadtree_constants::SPECIAL_COLOUR) {
                    // Note: SPECIAL_COLOUR represents the source. Since we only give the source
                    // a different colour to ensure that it's MortonBlock is included in the list
                    // (i.e. to cover case where it is an object). But we don't care what colour
                    // it has (never need to retrieve it) so if there is only one other colour in this
                    // quadtree leaf block then we consider the source having the same colour
                    if (colour == quadtree_constants::SPECIAL_COLOUR) {
                        // Another the possibility is that the colours are different
                        // because the source was encountered first, in this case
                        // we update the colour of this block (and as above do not
                        // consider the source having a different colour
                        colour = colourMap[nodes[i]];
                    } else {
                        sameColour = false;
                    }
                    //assert(colour != quadtree_constants::EMPTY_COLOUR && "Valid colour has not been set");
                }
            }
        }
        //assert(intervalInitialised && "Block interval was never initialised!");
       
        if (!sameColour) {
            this->addAndClearMergeCandidates(tempList,mergeColourSet,mergeColour,mergeCandidates);
            
            // If the colour is different we need to split this quadnode (recursively)
            //this->splitAndAddToList(source,qt,currNode,colourMap,tempList);
            this->splitAndAddToList(source,qt,currNode,colourMap,tempList,distRatios);
        } else {
            MortonBlock block(currNode->mortonCode,colour);
            block.setDistanceInterval(minBlockRatio,maxBlockRatio);            
            bool isFirstMergeCode = currNode->mortonCode.isFirstMergeCode();
            if (mergeCandidates.size() == 0 && isFirstMergeCode) {
                // There are no merge candidates and this is SW so it is a candidate
                // and we set the next expected sibling morton block (in this case SE)
                mergeCandidates.push_back(block);
                if (!mergeColourSet && !empty) { // mergeColourSet won't ever be true here - could remove it
                    // We can't set merge colour until we find a merge candidate that is not empty
                    mergeColour = colour;
                    mergeColourSet = true;
                }

                nextMergeCode = currNode->mortonCode.getNextMergeCode();
            } else if (mergeCandidates.size() != 0 && nextMergeCode == currNode->mortonCode) {
                // We found the sibling morton block we were expecting
                // Note: SPECIAL_COLOUR represents the source. Since we only give the source
                // a different colour to ensure that it's MortonBlock is included in the list
                // (i.e. to cover case where it is an object). But we don't care what colour
                // it has (never need to retrieve it) so we allow merging it with any colour
                if ((mergeColourSet && colour == mergeColour) || !mergeColourSet || empty || 
                    colour == quadtree_constants::SPECIAL_COLOUR || mergeColour == quadtree_constants::SPECIAL_COLOUR) {
                    if (mergeColourSet && mergeColour == quadtree_constants::SPECIAL_COLOUR && !empty) {
                        // If we are merging the source node's block (i.e. with SPECIAL_COLOUR)
                        // with another block that has a valid colour (i.e. not empty) we need
                        // to update it (source node's block is now guaranteed to be in list
                        // so we don't need to preserve it's colour)
                        mergeColour = colour;
                    }

                    // If colour is quadtree_constants::EMPTY_COLOUR then empty will be true
                    // so we will reach here for empty nodes
                    mergeCandidates.push_back(block);
                    if (!mergeColourSet && !empty) {
                        // We can't set merge colour until we find a merge candidate that is not empty
                        mergeColour = colour;
                        mergeColourSet = true;
                    }
                    
                    if (mergeCandidates.size() == 4) {
                        // If we have 4 candidates (which will be 4 sibling quandrants) then we can merge
                        // (this will be recursive in case merge makes another merge possible)
                        this->mergeAndSearchBackwards(tempList,mergeColourSet,mergeColour,mergeCandidates,nextMergeCode);
                    } else {
                        // We update the next expected sibling morton block
                        nextMergeCode = currNode->mortonCode.getNextMergeCode();
                    }
                } else {
                    // If colour is not merge colour and block is not empty then we cannot continue trying to merge current candidates
                    this->addAndClearMergeCandidates(tempList,mergeColourSet,mergeColour,mergeCandidates);
                    
                    // Current node cannot be a candidate if because it is a sibling quadrant
                    tempList.push_back(block);
                }
                
            } else {
                // Covers case where mergeCandidates() == 0 && quadrant != SW
                // Merge candidates is empty because the earlier quadrants have been eliminated already
                // (e.g. because we found two differently coloured blocks)
                
                // Covers case where mergeCandidates() != 0 nextMergeCode != code
                // Merge candidates cannot be merged because an unexpected morton block was found
                // It is still possible this unexpected block might be merged if it is SW block 
                // (in fact it must be a SW block)
                    
                this->addAndClearMergeCandidates(tempList,mergeColourSet,mergeColour,mergeCandidates);
                    
                if (!isFirstMergeCode) {
                    // Current block cannot be a candidate if it is not SW block
                    tempList.push_back(block);
                } else {
                    mergeCandidates.push_back(block);
                    if (!mergeColourSet && !empty) { // mergeColourSet won't ever be true here - could remove it
                        // We can't set merge colour until we find a merge candidate that is not empty
                        mergeColour = colour;
                        mergeColourSet = true;
                    }

                    // We update the next expected sibling morton block (in this case it will be SE)
                    nextMergeCode = currNode->mortonCode.getNextMergeCode();
                }
            }
        }
        
        currNode = currNode->nextLeaf;
    }
    
    // Build morton list from temporary list of morton blocks and remove
    // any remanining empty blocks
    for (auto it = tempList.begin(); it != tempList.end(); ++it) {
        if (it->getLink() != quadtree_constants::EMPTY_COLOUR) {
            this->mortonList.push_back(std::move(*it));
        }
    }
    
    // We know that binary search tree will contain all the elements of list
    this->binarySearchTreeMortonList.reserve(this->mortonList.size());
    // Build binary tree from linear quadtree Morton list
    // Note: We exploit the fact that Morton list is sorted
    std::vector<std::tuple<int,int,int>> currentLevel, nextLevel;
    unsigned int left = 0, right = this->mortonList.size()-1, mid, pos, size, largestPowerOf2, leftChildPos, rightChildPos;
    nextLevel.push_back(std::make_tuple(0,left,right));
    while (nextLevel.size() > 0) {
        currentLevel = nextLevel;
        nextLevel.clear();
        for (auto it = currentLevel.begin(); it != currentLevel.end(); ++it) {
            pos = std::get<0>(*it);
            left = std::get<1>(*it);
            right = std::get<2>(*it);
            size = right - left + 1; // right and left inclusive
            if (left == right) {
                //assert(static_cast<unsigned int>(pos) == this->binarySearchTreeMortonList.size());
                this->binarySearchTreeMortonList.push_back(std::move(this->mortonList[left]));
            } else if (left < right) {
                // We are building a complete binary search tree for space efficiency
                // and ease of inserts. This means that all levels are full up until
                // the final level where the elements are as far left as possible.
                
                // If n is the size of a binary search then the number elements on the left
                // hand side of the root is 2^x - 1 when the left handside is full (where 2^x
                // is the largest power of 2 that is less than n). However if the left-hand
                // side is not full the number of elements is n - 1 - (2^(x-1) - 1). This
                // is the total number of elements minus the root and minus the total number 
                // of elements in the right hand. Note that the right hand is always full, 
                // except the bottom level (which is why we use (x-1). Recall that we
                // always which fill from the left so since in thise case the left-hand side
                // is not full, the right hand side must be empty. n - 1 - (2^(x-1) - 1 can be
                // simplified to n - 2^(x-1), if n is full this expression will always greater
                // than or equal to 2^x - 1, because we are not considering the right-hand side
                // elements that may be present when the left-hand side is full. So we can
                // simply choose the minimum of the two and obtain the correct number of elments
                // on the left hand side.
                
                largestPowerOf2 = geometry::maxPower2(size); // Less than size
                mid = left + std::min(largestPowerOf2-1,size-(largestPowerOf2 >> 1));

                //assert(static_cast<unsigned int>(pos) == this->binarySearchTreeMortonList.size());
                this->binarySearchTreeMortonList.push_back(std::move(this->mortonList[mid]));

                leftChildPos = pos*2 + 1;
                nextLevel.push_back(std::make_tuple(leftChildPos,left,mid-1));
                if (size > 1) {
                    rightChildPos = pos*2 + 2;
                    nextLevel.push_back(std::make_tuple(rightChildPos,mid+1,right));
                } else {
                    // This means we have no more children to add because we reached
                    // the end of input sorted array (except left child of current)
                    //(this->binarySearchTreeMortonList.size() == this->mortonList.size()-1);
                }
            }
        }
    }
    
    assert(this->binarySearchTreeMortonList.size() == this->mortonList.size() && "Not all morton list elements were added to BST");
    
    // We no longer need mortonList (it will not be serialized) so we release memoryUsage
    utility::releaseSTLCollection(this->mortonList);
    
}

void MortonList::addAndClearMergeCandidates(std::vector<MortonBlock>& mortonList, bool& mergeColourSet, EdgeID& mergeColour, std::vector<MortonBlock>& mergeCandidates) {
    for (std::size_t i = 0; i < mergeCandidates.size(); ++i) {
        mortonList.push_back(mergeCandidates[i]);
    }
    mergeCandidates.clear();
    mergeColourSet = false;
    mergeColour = quadtree_constants::EMPTY_COLOUR;
}

void MortonList::mergeAndSearchBackwards(std::vector<MortonBlock>& mortonList, bool& mergeColourSet, EdgeID& mergeColour, std::vector<MortonBlock>& mergeCandidates, MortonCode& nextMergeCode) {
    // We assume numMergeCandidates is 4 (i.e. merge is possible)
    //assert (mergeColour != quadtree_constants::EMPTY_COLOUR && "Invalid merge colour set");
    
    // Create merged code (i.e. parent code) using any of the candidates
    // because all candiates will have the same parent block
    MortonBlock mergedBlock = mergeCandidates[0].createMergeBlock(mergeColour);

    // Set min and max block ratio of merged block by choosing
    // the smallest min and largest max ratio (resp.) of the child blocks
    // Note: We can do this since of the min network distance to a point 
    // merged would also be in the min network to one of the child
    // blocks (i.e. the one containing that point).
    //assert(mergeCandidates.size() == 4 && "Not enough blocks to merge");
    bool intervalInitialised = false;
    DistanceRatio minBlockRatio = 0, maxBlockRatio = 0;
    for (std::size_t i = 0; i < mergeCandidates.size(); ++i) {
        if (mergeCandidates[i].getLink() != quadtree_constants::EMPTY_COLOUR) {
            // Note: We still need to exclude empty blocks
            if (!intervalInitialised) {
                minBlockRatio = mergeCandidates[i].getMinDistRatio();
                maxBlockRatio = mergeCandidates[i].getMaxDistRatio();
                intervalInitialised = true;
            } else {
                if (mergeCandidates[i].getMinDistRatio() < minBlockRatio) {
                    minBlockRatio = mergeCandidates[i].getMinDistRatio();
                }
                if (mergeCandidates[i].getMaxDistRatio() > maxBlockRatio) {
                    maxBlockRatio = mergeCandidates[i].getMaxDistRatio();
                }
            }
        }
    }
    //assert(intervalInitialised); // Also makes sure at least one block is non-empty
    mergedBlock.setDistanceInterval(minBlockRatio,maxBlockRatio);
    mergeCandidates.clear(); // No longer need candidates as we final merged block

    // But it's possible that this merged block can now be merge with
    // morton blocks before or after it, so remember to merge further
    if (mergedBlock.isFirstMergeBlock()) {
        // If merged block is SW, we need only check forwards for possible further merges using merged block.
        // We set the next expected sibling morton block (in this case SE) and we return from this function.
        nextMergeCode = mergedBlock.getNextMergeCode();
        mergeCandidates.push_back(mergedBlock);
        
        // Note: We leave mergeColourSet and mergeColour to their
        // current values (will be true for merged block)
    } else if (mergedBlock.getMortonLength() != 0) {
        // If we have not merge all blocks (merged block's Morton len != 0), then we search backwards to merge further

        // If merged block is anything other quadrants, we need to check backwards, pop off any 
        // blocks that satisfy merge constraints, add them to the merge candidates and then either
        // return to forward search or (if we have 4 candidates) merge and search backwards again
        
        // We set the previous block to the expected sibling of the merged block
        // Note: Since mergeQuadrant is at least SW we are guaranteed at least one 
        // backward search iteration 
        MortonCode prevMergeCode = mergedBlock.getPrevMergeCode();
        
        bool searchBackwards = true;
        
        while (searchBackwards) {
            // Note: We can safely using std::vector::back() because there must be
            // a preceeding MortonBlock in the list if the current block is not
            // the first merging bloc (i.e. SW) - we haven't removed empty blocks yet
            MortonBlock& candidateBlock = mortonList.back();
            EdgeID candidateColour = candidateBlock.getLink();
            // Note: SPECIAL_COLOUR represents the source. Since we only give the source
            // a different colour to ensure that it's MortonBlock is included in the list
            // (i.e. to cover case where it is an object). But we don't care what colour
            // it has (never need to retrieve it) so we allow merging it with any colour
            if (candidateBlock == prevMergeCode && (candidateColour == mergeColour || candidateColour == quadtree_constants::EMPTY_COLOUR ||
                candidateColour == quadtree_constants::SPECIAL_COLOUR || mergeColour == quadtree_constants::SPECIAL_COLOUR)) {
                if (mergeColour == quadtree_constants::SPECIAL_COLOUR && candidateColour != quadtree_constants::EMPTY_COLOUR) {
                    // If we are merging the source node's block (i.e. with SPECIAL_COLOUR)
                    // with another block that has a valid colour (i.e. not empty)
                    // to update it (source node's block is now guaranteed to be in list
                    // so we don't need to preserve it's colour)
                    mergeColour = candidateColour;
                }                
                // If the previous morton block is what we are expecting and it is has the same colour it is merge candidate
                // Note: We don't have to check if mergeColour is not EMPTY_COLOUR because at least block is non-empty and has been used to mergeColour
                //assert(mergeColour != quadtree_constants::EMPTY_COLOUR);
                if (candidateBlock.isFirstMergeBlock()) {
                    // Stop search once we get to first quadrant of parent block
                    mergeCandidates.insert(mergeCandidates.begin(),candidateBlock);
                    mortonList.pop_back();
                    searchBackwards = false;
                    mergeCandidates.push_back(mergedBlock);
                    if (mergeCandidates.size() == 4) {
                        // If we have 4 candidates (which will be 4 sibling quandrants) then we can merge
                        // (this will be recursive in case merge makes another merge possible)
                        // Note: This will happen if merged block was NE block (i.e. last block)
                        this->mergeAndSearchBackwards(mortonList,mergeColourSet,mergeColour,mergeCandidates,nextMergeCode);
                    } else {
                        // We return from this function and continue forward search for merge candidates
                        // because have found all expected blocks during backward search but it is not enough
                        // to do a merge (meaning we need to continue forward search)
                        // We update the next expected sibling morton block based on the merged block
                        nextMergeCode = mergeCandidates.back().getNextMergeCode();
                    }
                } else {
                    // We set the previous block to the expected sibling of the merged block
                    prevMergeCode = candidateBlock.getPrevMergeCode();
                    mergeCandidates.insert(mergeCandidates.begin(),candidateBlock);
                    mortonList.pop_back();
                }
            } else {
                // This means one of the previous morton block has children and cannot be merged so
                // there is no chance that this merged block will be merged again
                searchBackwards = false;
                
                this->addAndClearMergeCandidates(mortonList,mergeColourSet,mergeColour,mergeCandidates);
                    
                // Merged block node cannot be a candidate because it is a sibling of blocks that are not candidates
                mortonList.push_back(mergedBlock);
            }
        }
    } else {
        mortonList.push_back(mergedBlock);
    }
}

void MortonList::splitAndAddToList(NodeID source, Quadtree& qt, QuadtreeNode* splitNode, std::vector<EdgeID>& colourMap, 
                                   std::vector<MortonBlock>& mortonList, std::vector<DistanceRatio>& distRatios) {
    // Create child block for population and checking
    std::vector<QuadtreeNode*> treeNodes(4);
    treeNodes[0] = new QuadtreeNode(splitNode->mortonCode,SW);
    treeNodes[1] = new QuadtreeNode(splitNode->mortonCode,SE);
    treeNodes[2] = new QuadtreeNode(splitNode->mortonCode,NW);
    treeNodes[3] = new QuadtreeNode(splitNode->mortonCode,NE);
    
    // Create variables to be populated with child block data
    std::deque<bool> sameColourStatus;
    std::vector<EdgeID> colours;
    std::vector<DistanceRatio> maxBlockRatios;
    std::vector<DistanceRatio> minBlockRatios;
    std::vector<std::vector<DistanceRatio>> childDistRatios;
    
    splitNode->copyDownNodesCheckColour(source,qt.mortonCodes,colourMap,distRatios,treeNodes,
                                        sameColourStatus,colours,minBlockRatios,maxBlockRatios,
                                        childDistRatios);
    
    for (std::size_t i = 0; i < 4; ++i) {
        if (sameColourStatus[i]) {
            if (colours[i] != quadtree_constants::EMPTY_COLOUR) {
                MortonBlock block(treeNodes[i]->mortonCode,colours[i]);
                //assert(minBlockRatios[i] != 0 && maxBlockRatios[i] != 0 && "Copy down function did not set correct ratios for block");
                block.setDistanceInterval(minBlockRatios[i],maxBlockRatios[i]);
                mortonList.push_back(block);
            } else {
                // There's no need to add these empty blocks to the Morton list
                // as they would removed anyway and since these were created been
                // by splitting they would never be mergeable anyway
                //MortonBlock block(treeNodes[i]->mortonCode,colours[i]);
                //mortonList.push_back(block);
            }
        } else {
            //assert(colours[i] != quadtree_constants::EMPTY_COLOUR && "Valid colour has not bee set");
            this->splitAndAddToList(source,qt,treeNodes[i],colourMap,mortonList,childDistRatios[i]);
        }
    }

    // Release temporary split quadtree nodes
    for (std::size_t i = 0; i < treeNodes.size(); ++i) {
        delete treeNodes[i];
    }
}

int MortonList::getNumBlocks() {
    return this->binarySearchTreeMortonList.size();
}

void MortonList::printBST()
{
    for (std::size_t i = 0; i < this->binarySearchTreeMortonList.size(); ++i) {
        MortonCode code = this->binarySearchTreeMortonList[i].getMortonCode();
        std::cout << std::setfill('0') << std::setw(3)  << i << ": ";
        code.printCode();
        std::cout << "Link: " << static_cast<unsigned int>(this->binarySearchTreeMortonList[i].getLink()) << std::endl;        
    }
}

void MortonList::printMortonList()
{
    for (std::size_t i = 0; i < this->mortonList.size(); ++i) {
        MortonCode code = this->mortonList[i].getMortonCode();
        std::cout << std::setfill('0') << std::setw(3)  << i << ": ";
        code.printCode();
        std::cout << "Link: " << static_cast<unsigned int>(this->mortonList[i].getLink()) << std::endl;        
    }
}

bool MortonList::isListSorted() {
    assert(this->mortonList.size() > 0 && "Linear Morton list has been released!");

    if (this->mortonList.size() > 1) {
//         std::cout << "\nMorton Block " << 0 << ": " << std::endl;        
//         this->mortonList[0].printBlock();
        for (std::size_t i = 1; i < this->mortonList.size(); ++i) {
//             std::cout << "\nMorton Block " << i << ": " << std::endl;        
//             this->mortonList[i].printBlock();
            MortonCode code1 = this->mortonList[i-1].getMortonCode();
            MortonCode code2 = this->mortonList[i].getMortonCode();
            if (code2 < code1) {
                code1.printCode();
                code2.printCode();
                return false;
            }
        }
    }
    return true;
}

bool MortonList::isListMergeable(int& count) {
    std::vector<MortonBlock> mergeCandidates;
    
    EdgeID mergeColour;
    bool mergeColourSet, mergeable = false;
    
    mergeColour = quadtree_constants::EMPTY_COLOUR;
    mergeColourSet = false;
    MortonCode nextMergeCode;
    
    assert(this->mortonList.size() > 0 && "Linear Morton list has been released!");

    for (std::size_t i = 0; i < this->mortonList.size(); ++i) {
        bool isFirstMergeBlock = this->mortonList[i].isFirstMergeBlock();
        EdgeID colour = this->mortonList[i].getLink();

        if (mergeCandidates.size() == 0 && isFirstMergeBlock) {
            // There are no merge candidates and this is SW so it is a candidate
            // and we set the next expected sibling morton block (in this case SE)
            mergeCandidates.push_back(this->mortonList[i]);
            if (!mergeColourSet && colour != quadtree_constants::EMPTY_COLOUR) { // mergeColourSet won't ever be true here - we could remove it
                // We can't set merge colour until we find a merge candidate that is not empty
                mergeColour = colour;
                mergeColourSet = true;
            }
            nextMergeCode = this->mortonList[i].getNextMergeCode();
//             std::cout << "Current Block" << std::endl;
//             this->mortonList[i].getMortonCode().printCode();
//             std::cout << "Next Expected Block" << std::endl;
//             nextMergeCode.printCode();
        } else if (mergeCandidates.size() != 0 && this->mortonList[i] == nextMergeCode) {
            // We found the sibling morton block we were expecting
            if ((mergeColourSet && colour == mergeColour) || !mergeColourSet || colour == quadtree_constants::EMPTY_COLOUR 
                || colour == quadtree_constants::SPECIAL_COLOUR || mergeColour == quadtree_constants::SPECIAL_COLOUR) {
                if (mergeColourSet && mergeColour == quadtree_constants::SPECIAL_COLOUR && colour != quadtree_constants::EMPTY_COLOUR) {
                    mergeColour = colour;
                }
                
                mergeCandidates.push_back(this->mortonList[i]);
                if (!mergeColourSet && colour != quadtree_constants::EMPTY_COLOUR) {
                    // We can't set merge colour until we find a merge candidate that is not empty
                    mergeColour = colour;
                    mergeColourSet = true;
                }
                
                if (mergeCandidates.size() == 4) {
                    mergeCandidates.clear();
                    mergeColourSet = false;
                    mergeable = true;
                    ++count;
                } else {
                    // We update the next expected sibling morton block
                    nextMergeCode = this->mortonList[i].getNextMergeCode();
//                     std::cout << "Current Block" << std::endl;
//                     this->mortonList[i].getMortonCode().printCode();
//                     std::cout << "Next Expected Block" << std::endl;
//                     nextMergeCode.printCode();
                }
            } else {
//                 std::cout << "\nRejected Candidates 1" << std::endl;
//                 for(std::size_t i = 0; i < mergeCandidates.size(); ++i) {
//                     MortonCode canCode = mergeCandidates[i].getMortonCode();
//                     canCode.printCode();
//                 }
//                 std::cout << "Breaking Candidate" << std::endl;
//                 this->mortonList[i].getMortonCode().printCode();
//                 std::cout << "Breaking Colour = " << static_cast<unsigned int>(colour) << std::endl;
//                 std::cout << "Expected Colour = " << static_cast<unsigned int>(mergeColour) << std::endl;
//                 std::cout << "Merge Colour Set = " << mergeColourSet << std::endl << std::endl; 
                
                // If colour is not merge colour and block is not empty then we cannot continue trying to merge current candidates
                mergeCandidates.clear();
                mergeColourSet = false;
                
                // Current node cannot be a candidate if because it is a sibling quadrant
            }
            
        } else {
            // Covers case where mergeCandidates() == 0 && quadrant != SW
            // Merge candidates is empty because the earlier quadrants have been eliminated already
            
            // Covers case where mergeCandidates() != 0 nextMergeCode != code
            // Merge candidates cannot be merged because an unexpected morton block was found
            // It is still possible this unexpected block might be merged if it is SW block 
            // (in fact it must be a SW block)
            
//             std::cout << "\nRejected Candidates 2" << std::endl;
//             for(std::size_t i = 0; i < mergeCandidates.size(); ++i) {
//                 MortonCode canCode = mergeCandidates[i].getMortonCode();
//                 canCode.printCode();
//             }
//             std::cout << "Breaking Candidate" << std::endl;
//             this->mortonList[i].getMortonCode().printCode();
//             std::cout << "Expected Candidate" << std::endl;
//             nextMergeCode.printCode();

            mergeCandidates.clear();
            mergeColourSet = false;
                
            if (!isFirstMergeBlock) {
                // Current block cannot be a candidate if it is not SW block
            } else {
                mergeCandidates.push_back(this->mortonList[i]);
                if (!mergeColourSet && colour != quadtree_constants::EMPTY_COLOUR) { // mergeColourSet won't ever be true here - we could remove it
                    // We can't set merge colour until we find a merge candidate that is not empty
                    mergeColour = colour;
                    mergeColourSet = true;
                }

                // We update the next expected sibling morton block (in this case it will be SE)
//                 nextMergeCode = this->mortonList[i].getNextMergeCode();
//                 std::cout << "Current Block" << std::endl;
//                 this->mortonList[i].getMortonCode().printCode();
//                 std::cout << "Next Expected Block" << std::endl;
//                 nextMergeCode.printCode();
            }
        }
    }
    return mergeable;
}

bool MortonList::areBlockIntervalsValid(NodeID source, Quadtree& qt, std::vector<EdgeWeight>& SSSPDistances) {
    bool validIntervals = true;
    DistanceRatio minRatio, maxRatio;
    for (std::size_t i = 0; i < this->binarySearchTreeMortonList.size(); ++i) {
        minRatio = this->binarySearchTreeMortonList[i].getMinDistRatio();
        maxRatio = this->binarySearchTreeMortonList[i].getMaxDistRatio();
        if (minRatio == 0 && this->binarySearchTreeMortonList[i].getLink() != quadtree_constants::SPECIAL_COLOUR) {
            std::cout << "Minimum ratio for block " << i << " is 0" << std::endl;
            validIntervals = false;
        }
        if (maxRatio == 0 && this->binarySearchTreeMortonList[i].getLink() != quadtree_constants::SPECIAL_COLOUR) {
            std::cout << "Maxmimum ratio for block " << i << " is 0" << std::endl;
            validIntervals = false;
        }
        if (minRatio > maxRatio) {
            std::cout << "Mininum ratio is larger than maximum ratio for block " << i << std::endl;
            validIntervals = false;
        }
    }
    
    Coordinate x, y;
    EuclideanDistance spatialDist;
    DistanceBound lowerBound, upperBound;
    for (std::size_t i = 0; i < SSSPDistances.size(); ++i) {
        if (i != source) {
            qt.populateCoordinates(i,x,y);
            MortonCode targetCode(x,y,qt.getMaxRangeExp());
            int idx = this->getBlockIdxContainingTargetByBST(targetCode);
            spatialDist = qt.getEuclideanDistance(source,i);
            minRatio = this->binarySearchTreeMortonList[idx].getMinDistRatio();
            maxRatio = this->binarySearchTreeMortonList[idx].getMaxDistRatio();
            lowerBound = std::floor(spatialDist*minRatio);
            upperBound = std::ceil(spatialDist*maxRatio);
            if (SSSPDistances[i] < lowerBound || SSSPDistances[i] > upperBound) {
                std::cout << "Interval using block " << idx << " for target " << i << " from source " << source << " is incorrect" << std::endl;
                std::cout << "Expected Distance " << SSSPDistances[i] << std::endl;
                std::cout << "Lowerbound " << lowerBound << std::endl;
                std::cout << "minRatio " << minRatio << std::endl;
                std::cout << "Upperbound " << upperBound << std::endl;
                std::cout << "maxRatio " << maxRatio << std::endl;
                targetCode.printCode();
                validIntervals = false;
            }
        }
    }
    
    return validIntervals;
}

double MortonList::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(this->binarySearchTreeMortonList);
    memoryUsage += sizeof(this->mortonList);
    for (std::size_t i = 0; i < this->binarySearchTreeMortonList.size(); ++i) {
        memoryUsage += this->binarySearchTreeMortonList[i].computeIndexSizeBytes();
    }
    return memoryUsage;
}

double MortonList::computeMemoryUsageBytes() {
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    for (std::size_t i = 0; i < this->binarySearchTreeMortonList.size(); ++i) {
        memoryUsage += this->binarySearchTreeMortonList[i].computeMemoryUsageBytes();
    }
    memoryUsage += sizeof(MortonBlock)*(this->binarySearchTreeMortonList.capacity()-this->binarySearchTreeMortonList.size());
    for (std::size_t i = 0; i < this->mortonList.size(); ++i) {
        memoryUsage += this->mortonList[i].computeMemoryUsageBytes();
    }
    memoryUsage += sizeof(MortonBlock)*(this->mortonList.capacity()-this->mortonList.size());
    return memoryUsage;
}

/*
 * MortonBlock
 */

// MortonBlock::MortonBlock(MortonCode code, EdgeID link): mortonCode(code), link(link), minDistRatio(0), maxDistRatio(0) {}
MortonBlock::MortonBlock(MortonCode code, EdgeID link): zNumber(code.getZNumber()), len(code.getLength()), link(link), minDistRatio(0), maxDistRatio(0) {}

MortonBlock::MortonBlock(MortonNumber zNumber, EdgeID mortonLen, EdgeID link): zNumber(zNumber), len(mortonLen), link(link) {}

void MortonBlock::setDistanceInterval(DistanceRatio minDistRatio, DistanceRatio maxDistRatio) {
    this->minDistRatio = minDistRatio;
    this->maxDistRatio = maxDistRatio;
}

EdgeID MortonBlock::getLink() {
    return this->link;
}

DistanceRatio MortonBlock::getMinDistRatio() {
    return this->minDistRatio;
}

DistanceRatio MortonBlock::getMaxDistRatio() {
    return this->maxDistRatio;
}

EdgeID MortonBlock::getLinkAndRatios(DistanceRatio& min, DistanceRatio& max)
{
    min = this->minDistRatio;
    max = this->maxDistRatio;
    return this->link;
}

void MortonBlock::printBlock() {
    std::cout << "Morton Code = " << std::bitset<64>(this->zNumber) << " (Len=" << static_cast<int>(this->len) << ")" << std::endl;
    std::cout << "Link: " << static_cast<unsigned int>(this->link) << std::endl;
}

MortonNumber MortonBlock::getMortonNumber()
{
    return this->zNumber;
}

EdgeID MortonBlock::getMortonLength()
{
    return this->len;
}

int MortonBlock::compareMortonCode(const MortonCode& code) const {
    // We assume morton code is shorter than input
    //assert (code.len >= this->len && "Input code is shorter than this code");
    MortonNumber testNum = 0;
    if (this->len > 0) {
        // We must shift by number less than MAX_MORTON_LENGTH otherwise it's undefined in C++
        // (but in our case shifting by length is equivalent to setting to 0)
        testNum = code.getZNumber();
        // Clear least significant bits that are not used by this code
        testNum = testNum >> (quadtree_constants::MAX_MORTON_LENGTH - this->len);
        testNum = testNum << (quadtree_constants::MAX_MORTON_LENGTH - this->len);
    }
    
    if (testNum < this->zNumber) {
        return -1;
    } else if (testNum > this->zNumber) {
        return 1;
    } else /* if (testNum == this->zNumber */ {
        return 0;
    }    
}

// Note: We use mergeColour instead of this block link colour as this block may be empty
MortonBlock MortonBlock::createMergeBlock(EdgeID mergeColour)
{
    //assert (this->len >= 2);
    MortonNumber parentNumber = 0;
    if (this->len > 2) {
        // We must shift by number less than MAX_MORTON_LENGTH otherwise it's undefined in C++
        // (but in our case shifting by length is equivalent to setting to 0)
        parentNumber = this->zNumber >> (quadtree_constants::MAX_MORTON_LENGTH - this->len + 2);
        parentNumber = parentNumber << (quadtree_constants::MAX_MORTON_LENGTH - this->len + 2);
    }
    MortonBlock parentBlock(parentNumber,this->len-2,mergeColour);
    // Note: Parent block has len-2 because the last digits represetend the child quadrant
    return parentBlock;
}

MortonCode MortonBlock::getMortonCode() {
    return MortonCode(this->zNumber,this->len);
}

void MortonBlock::decomposeAndGetRatios(int maxRangeExp, Coordinate& x, Coordinate& y, int& width, DistanceRatio& min, DistanceRatio& max)
{
    this->decompose(maxRangeExp,x,y,width);
    min = this->minDistRatio;
    max = this->maxDistRatio;
}

// This should not be called on a NE block (check before using isLastMergeCode)
MortonCode MortonBlock::getNextMergeCode()
{
    //assert (this->len >= 2);
    Direction nextQuadrant = this->getNextQuadrant();
    //assert (nextQuadrant != this->getQuadrant() && nextQuadrant > this->getQuadrant() && "Invalid input quadrant");
    MortonNumber parentNumber = 0;
    if (this->len > 2) {
        // We must shift by number less than MAX_MORTON_LENGTH otherwise it's undefined in C++
        // (but in our case shifting by length is equivalent to setting to 0)
        parentNumber = this->zNumber >> (quadtree_constants::MAX_MORTON_LENGTH - this->len + 2);
        parentNumber = parentNumber << (quadtree_constants::MAX_MORTON_LENGTH - this->len + 2);
    }
//     std::cout << "Original ZNumber = " << std::bitset<64>(this->zNumber) << std::endl;
//     std::cout << "Parent ZNumber = " << std::bitset<64>(parentNumber) << std::endl;
    MortonNumber nextMergeNum = this->getChildZNumber(parentNumber,nextQuadrant);
//     std::cout << "Next Merge ZNumber = " << std::bitset<64>(nextMergeNum) << std::endl;
    MortonCode nextCode(nextMergeNum,this->len);
    return nextCode;
}

// This should not be called on a SW block (check before using isFirstMergeCode)
MortonCode MortonBlock::getPrevMergeCode()
{
    //assert (this->len >= 2);
    Direction prevQuadrant = this->getPrevQuadrant();
    //assert (prevQuadrant != this->getQuadrant() && prevQuadrant < this->getQuadrant() && "Invalid input quadrant");
    MortonNumber parentNumber = 0;
    if (this->len > 2) {
        // We must shift by number less than MAX_MORTON_LENGTH otherwise it's undefined in C++
        // (but in our case shifting by length is equivalent to setting to 0)
        parentNumber = this->zNumber >> (quadtree_constants::MAX_MORTON_LENGTH - this->len + 2);
        parentNumber = parentNumber << (quadtree_constants::MAX_MORTON_LENGTH - this->len + 2);
    }
    MortonNumber prevMergeNum = this->getChildZNumber(parentNumber,prevQuadrant);
    MortonCode prevCode(prevMergeNum,this->len);
    return prevCode;
}

bool MortonBlock::isFirstMergeBlock()
{
    if (this->len >= 2) {
        // If length is zero it cannot be merged with anything so we return false
        MortonNumber quadrantDigits = this->zNumber << (this->len - 2);
        quadrantDigits = quadrantDigits >> (quadtree_constants::MAX_MORTON_LENGTH - 2);
        if (quadrantDigits == quadtree_constants::SW_MC) {
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

bool MortonBlock::isLastMergeBlock()
{
    if (this->len >= 2) {
        // If length is zero it cannot be merged with anything so we return false
        MortonNumber quadrantDigits = this->zNumber << (this->len - 2);
        quadrantDigits = quadrantDigits >> (quadtree_constants::MAX_MORTON_LENGTH - 2);
        if (quadrantDigits == quadtree_constants::NE_MC) {
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

Direction MortonBlock::getNextQuadrant()
{
    //assert (this->len >= 2);
    // This code quadrant is the last two occupied digits
    MortonNumber quadrantDigits = this->zNumber << (this->len - 2);
    quadrantDigits = quadrantDigits >> (quadtree_constants::MAX_MORTON_LENGTH - 2);
    if (quadrantDigits == quadtree_constants::SW_MC) {
        return SE;
    } else if (quadrantDigits == quadtree_constants::SE_MC) {
        return NW;
    } else if (quadrantDigits == quadtree_constants::NW_MC) {
        return NE;
    } else if (quadrantDigits == quadtree_constants::NE_MC) {
        std::cerr << "Current quadrant is NE - cannot getNextQuadrant" << std::endl;
        exit(1);
    } else {
        std::cerr << "Invalid quadrant from getQuadrant" << std::endl;
        exit(1);
    }
}

Direction MortonBlock::getPrevQuadrant()
{
    //assert (this->len >= 2);
    // This code quadrant is the last two occupied digits
    MortonNumber quadrantDigits = this->zNumber << (this->len - 2);
    quadrantDigits = quadrantDigits >> (quadtree_constants::MAX_MORTON_LENGTH - 2);
    if (quadrantDigits == quadtree_constants::SW_MC) {
        std::cerr << "Current quadrant is SW - cannot getPrevQuadrant" << std::endl;
        exit(1);
    } else if (quadrantDigits == quadtree_constants::SE_MC) {
        return SW;
    } else if (quadrantDigits == quadtree_constants::NW_MC) {
        return SE;
    } else if (quadrantDigits == quadtree_constants::NE_MC) {
        return NW;
    } else {
        std::cerr << "Invalid quadrant from getQuadrant" << std::endl;
        exit(1);
    }
}

MortonNumber MortonBlock::getChildZNumber(MortonNumber parentZNumber, Direction childQuadrant)
{
    //assert (this->len >= 2);
    MortonNumber childZNumber = 0;
    if (childQuadrant == SW) {
       childZNumber = parentZNumber | quadtree_constants::ARRAY_SW_BIT_MASK[this->len-2];
        // Note: This is equivalent to SW_MC << 64 - 2 - parentLen
        // I.e. we are the morton code for this quadrant to the end
        // of the parent block's morton code
    } else if (childQuadrant == SE) {
        childZNumber = parentZNumber | quadtree_constants::ARRAY_SE_BIT_MASK[this->len-2];
    } else if (childQuadrant == NW) {
        childZNumber = parentZNumber | quadtree_constants::ARRAY_NW_BIT_MASK[this->len-2];
    } else if (childQuadrant == NE) {
        childZNumber = parentZNumber | quadtree_constants::ARRAY_NE_BIT_MASK[this->len-2];
    } else {
        std::cerr << "Invalid direction provided" << std::endl;
        exit(1);
    }
    return childZNumber;
}

void MortonBlock::decompose(int maxRangeExp, Coordinate& x, Coordinate& y, int& width)
{
    MortonNumber unshiftedZ = this->zNumber >> (quadtree_constants::MAX_MORTON_LENGTH - (maxRangeExp << 1)); 
    // Note: We can't use this->length here because we don't know if it is point or block
    
    MortonNumber mx = unshiftedZ & quadtree_constants::B[0], my = (unshiftedZ >> 1) & quadtree_constants::B[0];

    mx = (mx | (mx >> quadtree_constants::S[0])) & quadtree_constants::B[1];
    mx = (mx | (mx >> quadtree_constants::S[1])) & quadtree_constants::B[2];
    mx = (mx | (mx >> quadtree_constants::S[2])) & quadtree_constants::B[3];
    mx = (mx | (mx >> quadtree_constants::S[3])) & quadtree_constants::B[4];
    mx = (mx | (mx >> quadtree_constants::S[4])) & quadtree_constants::B[5];

    my = (my | (my >> quadtree_constants::S[0])) & quadtree_constants::B[1];
    my = (my | (my >> quadtree_constants::S[1])) & quadtree_constants::B[2];
    my = (my | (my >> quadtree_constants::S[2])) & quadtree_constants::B[3];
    my = (my | (my >> quadtree_constants::S[3])) & quadtree_constants::B[4];
    my = (my | (my >> quadtree_constants::S[4])) & quadtree_constants::B[5];
    
    x = mx;
    y = my;
    
    unsigned int widthPower = maxRangeExp - (this->len >> 1); // No precision will be lost because this->len is always multiple of two
    width = 1;
    width = width << widthPower;
}

bool MortonBlock::operator==(const MortonCode& code) const  {
    if (this->len == code.getLength() && this->zNumber == code.getZNumber()) {
        return true;
    } else {
        return false;
    }
}

double MortonBlock::computeIndexSizeBytes()
{
    double memoryUsage = 0;
    memoryUsage += sizeof(this->link);
    memoryUsage += sizeof(this->minDistRatio);
    memoryUsage += sizeof(this->maxDistRatio);
    memoryUsage += sizeof(this->zNumber);
    memoryUsage += sizeof(this->len);
    return memoryUsage;
}

double MortonBlock::computeMemoryUsageBytes() {
    double memoryUsage = 0;
    memoryUsage += sizeof(*this);
    return memoryUsage;
}
