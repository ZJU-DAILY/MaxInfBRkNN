
#include "TransformInputCommand.h"





void TransformInputCommand::showCommandUsage(std::string programName)
{
    std::cerr << "Usage: " << programName << " -c " + constants::TRANSFORM_INPUT_CMD + " -m <query method>\n\n"
              << "Methods Options:\n"
              << utility::getFormattedUsageString(TRANSFORM_DIMACS,"Process DIMACS files") + "\n"
              << utility::getFormattedUsageString(TRANSFORM_TPQ,"Process TPQ files") + "\n"
              << utility::getFormattedUsageString(TRANSFORM_OSM_DIMACS,"Process DIMACS files created from OSM") + "\n";
}

void TransformInputCommand::showMethodUsage(std::string method, std::string programName)
{
    if (method == TRANSFORM_DIMACS ) {
        std::cerr << "Usage: " << programName << " -c " + constants::TRANSFORM_INPUT_CMD
                  << " -m " + TRANSFORM_DIMACS + " -e <graph edge file>"
                  << " -n <graph node file> \n-o <prefix for output files>\n";
    } else if (method == TRANSFORM_TPQ ) {
        std::cerr << "Usage: " << programName << " -c " + constants::TRANSFORM_INPUT_CMD
                  << " -m " + TRANSFORM_TPQ + " -e <graph edge file>"
                  << " -n <graph node file> \n-o <prefix for output files> -w <weight inflate factor>"
                  << " -c <coordinate inflate factor>\n";
    } else if (method == TRANSFORM_OSM_DIMACS ) {
        std::cerr << "Usage: " << programName << " -c " + constants::TRANSFORM_INPUT_CMD
                  << " -m " + TRANSFORM_OSM_DIMACS + " -e <graph edge file>"
                  << " -n <graph node file> \n-o <prefix for output files>-r <region name> -s <sub-region name>\n";
    } else {
        this->showCommandUsage(programName);
    }
}

void TransformInputCommand::addNode(NodeID node, Coordinate x, Coordinate y)
{
    if (node == nodes.size()) {
        this->nodes.push_back(node);
        this->coordinates.push_back(CoordinatePair(x,y));
    } else {
        std::cerr << "Nodes are not in order" << std::endl;
        std::exit(1);
    }
}

void TransformInputCommand::addEdge(NodeID source, NodeID target, EdgeWeight weight)
{
    ///std::cout<<"attempt to insert ("<<source<<","<<target<<")="<<weight<<endl;
    ///getchar();
    if (weight > 0) {
        this->neighbours[source].push_back(NodeEdgeWeightPair(target,weight));
    } else {
        std::cerr << "Edge weight is zero" << std::endl;
        std::exit(1);
    }
}

bool TransformInputCommand::checkValidNode(NodeID node)
{
    if (node >= this->nodes.size()) {
        std::cerr << "Node (" << node << ") does not exist (too large)" << std::endl;
        exit(1);
    }
    return true;
}

bool TransformInputCommand::isValidEdge(NodeID source, NodeID target, int rawWeight, EdgeWeight weight)
{
    // We assume weight has been corrected for Euclidean distance
    
    bool isValidEdge = true;
    if (source == target) {
        // Remove self-loops
        ++this->numSelfLoops;
        isValidEdge = false;
    }
    
    if (isValidEdge && rawWeight < 0) {
        // Remove negative edges
        ++this->numNegativeEdgeWeights;
        isValidEdge = false;
    }
    
    // Remove parallel edges (i.e. multiple edges to same target)
    if (isValidEdge) {
        // Since weight it's positive if we reach here we cast to EdgeWeight
        for (std::size_t i = 0; i < this->neighbours[source].size(); ++i) {
            if (this->neighbours[source][i].first == target) {
                ++this->numParallelEdges;
                isValidEdge = false;
                break;
            }
        }
    }
    return isValidEdge;
}

EdgeWeight TransformInputCommand::correctEdgeWeightByEuclidDist(NodeID source, NodeID target, EdgeWeight weight)
{
    // For distance-based edge weights, edge weight cannot be 
    // smaller than euclidean distance between source and target
    if (this->edgeType == constants::DISTANCE_WEIGHTS) {
        // We ceil to make sure edge weight is always bigger than euclidean distance
        EdgeWeight minWeight = std::ceil(geometry::getEuclideanDist(this->coordinates[source].first,this->coordinates[source].second,
                                                                    this->coordinates[target].first,this->coordinates[target].second));
        if (minWeight == 0) {
            std::cerr << "There is an edge between two nodes that have same coordinates!" << std::endl;
            std::cerr << "Node " << source << " has the same coordinates as node " << target << std::endl;
            exit(1);
        }
        if (weight < minWeight) {
            this->weightAddedByEuclideanCorrections += (minWeight-weight);
            ++this->numEuclideanCorrections;
            if (weight*10 < minWeight) {
                //std::cout << "Corrected edge between node " << source << " and node " << target << " by over a factor of 10" << std::endl;
                this->numCorrectionsOverFactorTen++;
            }
            return minWeight;
        }
    }
    return weight;
}

void TransformInputCommand::fixShortestPathOverflow()
{
    StopWatch sw;
    sw.start();
    LongPathDistance estMaxSPDist = this->computePseudoDiameter(100);
    sw.stop();
    //std::cout << "Pseudo Diameter Calculated in " << sw.getTimeMs() << "ms" << std::endl;
    //std::cout << "Pseudo Diameter = " << estMaxSPDist << std::endl;
    
    // We adjust weights so that maximum possible shortest (i.e. one involving)
    // all edges fit in EdgeWeight (at it's smallest it is 32-bits)
    LongPathDistance maxPossibleSPDist = std::numeric_limits<EdgeWeight>::max();
    
    if (estMaxSPDist > maxPossibleSPDist) {
        this->wasDeflated = true;
        int deflateShift = 1;
        while ((estMaxSPDist >> deflateShift) > maxPossibleSPDist) {
            ++deflateShift;
        }
        this->deflationFactor = 1 << deflateShift;

        // We need to deflate X and Y coordinates first so Euclidean distance will match new weights
        for (std::size_t i = 0; i < this->nodes.size(); ++i) {
            this->coordinates[i].first = std::ceil(this->coordinates[i].first >> deflateShift);
            this->coordinates[i].second = std::ceil(this->coordinates[i].second>> deflateShift);
        }
        
        // Adjust all edge weights by deflation factor (but check Euclidean distances to be sure)
        unsigned long long newTotalEdgeWeight = 0;
        for (std::size_t i = 0; i < this->nodes.size(); ++i) {
            for (std::size_t j = 0; j < this->neighbours[i].size(); ++j) {
                this->neighbours[i][j].second = this->neighbours[i][j].second;
                if (this->neighbours[i][j].second == 0) {
                    std::cerr << "Deflation has caused a zero edge weight!" << std::endl;
                    exit(1);
                }
                newTotalEdgeWeight += this->neighbours[i][j].second;
                EdgeWeight minWeight = std::ceil(geometry::getEuclideanDist(this->coordinates[i].first,this->coordinates[i].second,
                                                                            this->coordinates[this->neighbours[i][j].first].first,
                                                                            this->coordinates[this->neighbours[i][j].first].second));
                // Note: It's possible that if a Coordinate collision was caused by deflation the minWeight is
                // zero. That's fine, it will not affect the edge weight and we will remove non-unique Nodes later
                if (this->neighbours[i][j].second < minWeight) {
                    this->neighbours[i][j].second = minWeight;
                }
            }
        }
        
        this->oldTotalEdgeWeight = this->totalEdgeWeight;
        this->totalEdgeWeight = newTotalEdgeWeight;
    }
}

LongPathDistance TransformInputCommand::computePseudoDiameter(int numRounds)
{
    std::vector<bool> isNodeSettled(this->numNodes,false);
    BinaryMinHeap<LongPathDistance,NodeID> pqueue;
    
    LongPathDistance minDist, maxSPDist = 0, prevMaxSPDist = 0;
    NodeID minDistNodeID, adjNode, maxSPDistNode, maxSPDistNodeDegree = std::numeric_limits<NodeID>::max();
    
    // Initialize with priority queue with source node
    pqueue.insert(0,0);
    int i;
    maxSPDistNode = 0;
    for (i = 0; (i < numRounds || numRounds == -1) && pqueue.size() != 0; ++i) {
        maxSPDist = 0;
        while (pqueue.size() > 0) {
            // Extract and remove node with smallest distance from source
            // and mark it as "settled" so we do not inspect again
            minDist = pqueue.getMinKey();
            minDistNodeID = pqueue.extractMinElement();
            if (!isNodeSettled[minDistNodeID]) {
                isNodeSettled[minDistNodeID] = true;
                
                if (maxSPDist < minDist) {
                    maxSPDist = minDist;
                    maxSPDistNode = minDistNodeID;
                    maxSPDistNodeDegree = this->neighbours[minDistNodeID].size();
                } else if (maxSPDist == minDist && maxSPDistNodeDegree > this->neighbours[minDistNodeID].size()) {
                    // If the distances are equal we choose the node with the smaller degree
                    maxSPDistNode = minDistNodeID;
                    maxSPDistNodeDegree = this->neighbours[minDistNodeID].size();
                }

                for (std::size_t i = 0; i < this->neighbours[minDistNodeID].size(); ++i) {
                    adjNode = this->neighbours[minDistNodeID][i].first;
                    // Only update those we haven't already settled
                    if (!isNodeSettled[adjNode]) {
                        pqueue.insert(adjNode,minDist+this->neighbours[minDistNodeID][i].second);
                    }
                }
            }
        }    
        
        isNodeSettled.assign(this->numNodes,false);
        pqueue.clear();
        
        // Make the node with maximum shortest path the next node
        // unless we did not improve the max shortest path distance
        // (which case we stop irrespective of number of rounds)
        //std::cout << "New Pseudo Diameter Estimate = " << maxSPDist << std::endl;            
        if (prevMaxSPDist < maxSPDist) {
            pqueue.insert(maxSPDistNode,0);
            prevMaxSPDist = maxSPDist;
        }
    }
    //std::cout << "PseudoDiameter Calculation Rounds = " << i << std::endl;
    return prevMaxSPDist;
}

void TransformInputCommand::computeCoordinateRanges()
{
    Coordinate minX, minY, maxX, maxY;  
    
    minX = this->coordinates[0].first;
    maxX = this->coordinates[0].first;
    minY = this->coordinates[0].second;
    maxY = this->coordinates[0].second;    
    for (std::size_t i = 1; i < nodes.size(); ++i) {
        if (this->coordinates[i].first < minX) {
            minX = this->coordinates[i].first;
        }
        if (this->coordinates[i].first > maxX) {
            maxX = this->coordinates[i].first;
        }
        if (this->coordinates[i].second < minY) {
            minY = coordinates[i].second;
        }
        if (this->coordinates[i].second > maxY) {
            maxY = this->coordinates[i].second;
        }
    }

    // Take opportunity to note x and y translations required to make
    // origin (0,0) and hence use Morton numbers to check uniqueness later
    this->xTranslation = minX;
    this->yTranslation = minY;

}

void TransformInputCommand::computeCoordinateTranslations(std::vector<TransformInputCommand::Point>& coordinates, int& minX, int& minY)
{
    Coordinate testMinX = coordinates[0].x;
    Coordinate testMinY = coordinates[0].y;  
    for (std::size_t i = 1; i < nodes.size(); ++i) {
        if (coordinates[i].x < testMinX) {
            testMinX = coordinates[i].x;
        }
        if (coordinates[i].y < testMinY) {
            testMinY = coordinates[i].y;
        }
    }

    // Take opportunity to note x and y translations required to make
    // origin (0,0) and hence use Morton numbers to check uniqueness later
    minX = testMinX;
    minY = testMinY;
}

void TransformInputCommand::checkNonUniqueCoordinates()
{
    int nonUniqueCoordinates = 0;
    std::unordered_set<MortonNumber> zNumbers;    
    for (NodeID i = 0; i < this->nodes.size(); ++i) {
        MortonCode code(this->coordinates[i].first-this->xTranslation,this->coordinates[i].second-this->yTranslation,31);
        MortonNumber zNumber = code.getZNumber();
        if (zNumbers.find(zNumber) == zNumbers.end()) {
            zNumbers.insert(code.getZNumber());
        } else {
//             MortonCode testCode(zNumber,1);
//             Coordinate testX, testY;
//             int width;
//             testCode.decompose(31,testX,testY,width);
//             std::cout << "Node " << i << " has non-unique coordinates (" << testX+this->xTranslation << "," << testY+this->yTranslation << ")" << std::endl;
            this->nonUniqueNodes.insert(i);
            ++nonUniqueCoordinates;
        }
    }
    if (nonUniqueCoordinates > 0) {
        std::cout << "Found " << nonUniqueCoordinates << " non-unique coordinates" << std::endl;
    }
}

void TransformInputCommand::removeNonUniqueNodes()
{
    std::cerr << "This function doesn't work - removing nodes currently may disconnect entire subgraphs but we don't remove them" << std::endl;
    if (this->nonUniqueNodes.size() == 0) {
        // Nothing to do if there are no non-unique nodes
        return;
    }
    std::unordered_map<NodeID,NodeID> newToOldNodeID;
    std::unordered_map<NodeID,NodeID> oldToNewNodeID;
    std::vector<NodeID> newNodes;
    std::vector<std::vector<NodeEdgeWeightPair>> newNeighbours;
    std::vector<CoordinatePair> newCoordinates;
    NodeID newID = 0;
    for (NodeID i = 0; i < this->nodes.size(); ++i) {
        if (this->nonUniqueNodes.find(i) == this->nonUniqueNodes.end()) {
            // This node should not be removed so we give it a new renumbered ID
            newToOldNodeID[newID] = i;
            oldToNewNodeID[i] = newID;
            newNodes.push_back(newID);
            newCoordinates.push_back(std::move(this->coordinates[i]));
            ++newID;
            assert(newID == newNodes.size() && "Re-numbering has failed");
        } else {
            this->numNonUniqueNodesRemoved++;
            // This node's edges will be implicitly removed
            this->numNonUniqueNodeEdgesRemoved += this->neighbours[i].size();
        }
    }
    this->numNodes = newNodes.size();

    // We copy the old adjacency list for each node, excluding
    // neighbours that are non-unique (removed above)
    
    // Note: This might result in some nodes becoming disconnected, 
    // so we later call removeDisconnectedNodes() to remove them
    NodeID oldNodeID, neighbourNewID;
    this->numEdges = 0;
    for (NodeID i = 0; i < newNodes.size(); ++i) {
        assert(newToOldNodeID.find(i) != newToOldNodeID.end() && "New to old NodeID mapping doesn't exist");
        oldNodeID = newToOldNodeID[i];
        std::vector<NodeEdgeWeightPair> newAdjacencyList;
        for (std::size_t j = 0; j < this->neighbours[oldNodeID].size(); ++j) {
            if (this->nonUniqueNodes.find(this->neighbours[oldNodeID][j].first) == this->nonUniqueNodes.end()) {
                assert(oldToNewNodeID.find(this->neighbours[oldNodeID][j].first) != oldToNewNodeID.end() && "Old to new NodeID mapping doesn't exist");
                neighbourNewID = oldToNewNodeID[this->neighbours[oldNodeID][j].first];
                newAdjacencyList.push_back(NodeEdgeWeightPair(neighbourNewID,this->neighbours[oldNodeID][j].second));
            } else {
                this->numNonUniqueNodeEdgesRemoved++;
            }
        }
        this->numEdges += newAdjacencyList.size();
        newNeighbours.push_back(newAdjacencyList);
    }
    this->nodes = newNodes;
    this->neighbours = newNeighbours;
    this->coordinates = newCoordinates;
    assert(this->nodes.size() == this->neighbours.size());
    assert(this->neighbours.size() == this->coordinates.size());
}

void TransformInputCommand::removeIsolatedNodes()
{
    // Note: Removing disconnected nodes will not cause more disconnected
    // to be created (unlike removing non-unique) because by definition a
    // disconnected node is not present in any other nodes adjacency list
    std::unordered_map<NodeID,NodeID> newToOldNodeID;
    std::unordered_map<NodeID,NodeID> oldToNewNodeID;
    std::vector<NodeID> newNodes;
    std::vector<std::vector<NodeEdgeWeightPair>> newNeighbours;
    std::vector<CoordinatePair> newCoordinates;
    std::unordered_set<NodeID> disconnectedNodes;
    
    NodeID newID = 0;
    for (NodeID i = 0; i < this->nodes.size(); ++i) {
        if (this->neighbours[i].size() != 0) {
            // This is not a disconnected node
            newToOldNodeID[newID] = i;
            oldToNewNodeID[i] = newID;
            newNodes.push_back(newID);
            newCoordinates.push_back(std::move(this->coordinates[i]));
            ++newID;
            assert(newID == newNodes.size() && "Re-numbering has failed");
        } else {
            disconnectedNodes.insert(i);
            this->numDisconnectedNodesRemoved++;
            // Note: We do not need to count how many edges are removed
            // because we have already checked this node has no edges
        }
    }

    if (disconnectedNodes.size() == 0) {
        // Nothing to do, no disconnected nodes found
        return;
    }
    
    this->numNodes = newNodes.size();
    NodeID oldNodeID, neighbourNewID;
    for (NodeID i = 0; i < newNodes.size(); ++i) {
        assert(newToOldNodeID.find(i) != newToOldNodeID.end() && "New to old NodeID mapping doesn't exist");
        oldNodeID = newToOldNodeID[i];
        std::vector<NodeEdgeWeightPair> newAdjacencyList;
        for (std::size_t j = 0; j < this->neighbours[oldNodeID].size();  ++j) {
            if (disconnectedNodes.find(this->neighbours[oldNodeID][j].first) == disconnectedNodes.end()) {
                assert(oldToNewNodeID.find(this->neighbours[oldNodeID][j].first) != oldToNewNodeID.end() && "Old to new NodeID mapping doesn't exist");
                neighbourNewID = oldToNewNodeID[this->neighbours[oldNodeID][j].first];
                newAdjacencyList.push_back(NodeEdgeWeightPair(neighbourNewID,this->neighbours[oldNodeID][j].second));
            } else {
                assert(false && "Cannot reach here because disconnected node, by defintion, should not appear in any other node's adjacency list");
            }
        }
        newNeighbours.push_back(newAdjacencyList);
    }
    this->nodes = newNodes;
    this->neighbours = newNeighbours;
    this->coordinates = newCoordinates;
    assert(this->nodes.size() == this->neighbours.size());
    assert(this->neighbours.size() == this->coordinates.size());
}

void TransformInputCommand::removeIsolatedRegions()
{

    std::vector<std::unordered_set<NodeID>> nodeRegions;
    NodeID currentNode, neighbourNodeID;
    int unassignedNodes = this->numNodes;
    NodeID unassignedNode = 0;
    
    std::vector<bool> isNodeFound(this->numNodes,false);
    std::deque<NodeID> pqueue;
    std::unordered_set<NodeID> currentRegion;
    
    std::size_t biggestRegionSize = 0;
    
    while (unassignedNodes > 0) {
        // Clear previous iteration data
        this->numRegionsFound++;
        pqueue.clear();
        currentRegion.clear();
        
        pqueue.push_back(unassignedNode);
        while (pqueue.size() > 0) {
            currentNode = pqueue.front();
            pqueue.pop_front();
            if (!isNodeFound[currentNode]) {
                currentRegion.insert(currentNode);
                isNodeFound[currentNode] = true;
                unassignedNodes--;

                for (std::size_t i = 0; i < neighbours[currentNode].size(); ++i) {
                    neighbourNodeID = neighbours[currentNode][i].first;
                    // Only add those we haven't already seen
                    if (!isNodeFound[neighbourNodeID]) {
                        pqueue.push_back(neighbourNodeID);
                    }
                }
            }
        }
        
        if (currentRegion.size() > biggestRegionSize) {
            biggestRegionSize = currentRegion.size();
        }
        
        nodeRegions.push_back(std::move(currentRegion));
        
        if (unassignedNodes > 0) {
            // Find a new unassigned node (as root not of new BFS search)
            // (i.e. it is not present any of the region found so far)
            bool nodeFound = false;
            for (NodeID node = 0; node < this->numNodes; ++node) {
                if (!isNodeFound[node]) {
                    nodeFound = true;
                    unassignedNode = node;
                    break;
                }
            }
            assert(nodeFound && "There are unassigned nodes, but we could not find any one of them!");
        }
    
    }
    
    // Remove all isolated regions except the biggest region
    bool biggestRegionPreserved = false; // In case there are two regions with same biggest size
    for (std::size_t i = 0; i < nodeRegions.size(); ++i) {
        if (nodeRegions[i].size() != biggestRegionSize || biggestRegionPreserved) {
            if (nodeRegions[i].size() == biggestRegionSize) {
                biggestRegionPreserved = true;
            }
            for (NodeID node: nodeRegions[i]) {
                // Clear it's adjacency list so it will be deleted as isolated node
                // This is OK because we know all of the edge are only within this region
                this->numIsolatedRegionNodesRemoved++;
                this->numIsolatedRegionEdgesRemoved += this->neighbours[node].size();
                this->numEdges -= this->neighbours[node].size(); // this->numNodes will be decremented later
                this->neighbours[node].clear();
            }
        }
    }
    
}



void TransformInputCommand::printTransformationStatistics()
{
    std::cout << "\nGraph " << this->networkName << " processed with " << this->numNodes << " nodes and " << this->numEdges << " edges" << std::endl;
    std::cout << "Total edge weight = " << this->totalEdgeWeight << std::endl;
    std::cout << "Average edge weight = " << static_cast<double>(this->totalEdgeWeight)/this->numEdges << std::endl;

    std::cout << "\nStage 1: Validate Input File Graph" << std::endl;
    std::cout << "Number of self-loops removed = " << this->numSelfLoops << std::endl;
    std::cout << "Number of parallel edges removed = " << this->numParallelEdges << std::endl;        
    std::cout << "Number of negative weight edges removed = " << this->numNegativeEdgeWeights << std::endl;        
    std::cout << "Number of edges smaller than Euclidean distance corrected = " << this->numEuclideanCorrections << std::endl;        
    std::cout << "Number of edges smaller than Euclidean distance divided by 10 corrected = " << this->numCorrectionsOverFactorTen << std::endl;        
    if (this->numEuclideanCorrections != 0) {
        std::cout << "Average euclidean correction = " << static_cast<double>(this->weightAddedByEuclideanCorrections)/this->numEuclideanCorrections << std::endl;
    }

    std::cout << "\nStage 2: Correct Possible Shortest Path Distance Overflows" << std::endl;
    if (this->wasDeflated) {
        std::cout << "Sum of edge weights exceeded max possible shortest path distance value" << std::endl;
        std::cout << "Coordinate and EdgeWeight values were adjusted by deflation factor of " << this->deflationFactor << std::endl;
        std::cout << "Sum of edge weights was reduced from " << this->oldTotalEdgeWeight << " to " << this->totalEdgeWeight << std::endl;
    } else {
        std::cout << "No deflation performed, total of all edge weights below maximum possible shortest path distance value" << std::endl;
    }

    std::cout << "\nStage 3: Remove Non-Unique Nodes" << std::endl;
    std::cout << "Number of non-unique nodes removed = " << this->numNonUniqueNodesRemoved << std::endl;
    std::cout << "Number of non-unique node edges removed = " << this->numNonUniqueNodeEdgesRemoved << std::endl;

    std::cout << "\nStage 4: Disconnect Isolated Region Nodes" << std::endl;
    std::cout << "Number of isolated regions found  = " << this->numRegionsFound-1 << std::endl; 
    std::cout << "Number of isolated region nodes made disconnected (removed later as below) = " << this->numIsolatedRegionNodesRemoved << std::endl;
    std::cout << "Number of isolated region edges removed = " << this->numIsolatedRegionEdgesRemoved << std::endl;
    
    std::cout << "\nStage 5: Remove Disconnected Nodes" << std::endl;
    std::cout << "Number of disconnected nodes removed = " << this->numDisconnectedNodesRemoved << std::endl; 
    std::cout << "Number of disconnected node edges removed = " << this->numDisconnectedNodeEdgesRemoved << std::endl; 

}


