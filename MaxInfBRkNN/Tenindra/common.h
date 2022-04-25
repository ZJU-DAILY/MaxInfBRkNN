
#ifndef _COMMON_H
#define _COMMON_H

#include <string>



typedef unsigned int ScaledWeight;
typedef unsigned int EdgeWeight;
typedef unsigned int NodeID;
typedef unsigned int EntryKey;
typedef uint8_t EdgeID; // Adjacency List Index
typedef long Coordinate;
typedef unsigned int CoordinateRange;
typedef double EuclideanDistance;
typedef double EuclideanDistanceSqrd;
typedef unsigned int DistanceBound;
typedef float DistanceRatio;
typedef unsigned long long MortonNumber;
typedef uint8_t MortonLength;
typedef EdgeWeight PathDistance;
typedef unsigned long long LongPathDistance;
typedef std::pair<NodeID,EdgeWeight> NodeDistancePair;
typedef std::pair<NodeID,EdgeWeight> NodeEdgeWeightPair;
typedef std::pair<NodeID,bool> NodeStatusPair;
typedef std::pair<NodeID,NodeID> NodePair;
typedef std::pair<NodeID,uint8_t> NodeLinkPair;
typedef std::pair<Coordinate,Coordinate> CoordinatePair; // x is first, y is second
typedef std::pair<int,EdgeWeight> IntDistancePair;
typedef std::pair<int,int> IntIntPair;


namespace constants
{
    // Command Names (must be less than constants::MAX_COMMAND_NAME_SIZE)
    std::string const IDX_GRAPH_CMD = "binary_graph";
    std::string const IDX_GTREE_CMD = "gtree";
    std::string const IDX_DYNAMICGRAPH_CMD = "dynamic_graph";
    std::string const IDX_ROUTEOVERLAY_CMD = "route_overlay";
    std::string const IDX_SILC_CMD = "silc";
    std::string const IDX_JUNC_CMD = "junc";
    std::string const IDX_PHL_CMD = "phl";
    std::string const IDX_ALT_CMD = "alt";
    std::string const IDX_CH_CMD = "ch";
    std::string const IDX_TNR_CMD = "tnr";
    std::string const NETWORK_IDX_CMD = "idx";
    std::string const QUERY_KNN_CMD = "knn";
    std::string const OBJ_IDX_CMD = "obj_idx";
    std::string const SAMPLE_SET_CMD = "sample_set";
    std::string const TRANSFORM_INPUT_CMD = "transform";
    std::string const EXPERIMENTS_CMD = "exp";

    // Point/Object Set Types
    std::string const RAND_OBJ_SET = "random";
    std::string const PARTITION_OBJ_SET = "partition";
    std::string const CLUSTER_OBJ_SET = "cluster";
    std::string const MINND_OBJ_SET = "min_nd";
    std::string const MINMAXND_OBJ_SET = "minmax_nd";
    std::string const RAND_PAIRS_SET = "random_pairs";
    std::string const POI_OBJECT_SET = "poi";

    // Object Index Types
    std::string const OBJ_IDX_GTREE = "occ_list";
    std::string const OBJ_IDX_ROAD = "assoc_dir";
    std::string const OBJ_IDX_QUADTREE = "quadtree";
    std::string const OBJ_IDX_RTREE = "rtree";

    // Shortest Path Distance Query Methods
    std::string const DIJKSTRA_SPDIST_QUERY = "dijkstra";
    std::string const ASTAR_SPDIST_QUERY = "a_star";
    std::string const GTREE_SPDIST_QUERY = "gtree";
    std::string const SILC_SPDIST_QUERY = "silc";
    std::string const PHL_SPDIST_QUERY = "phl";
    std::string const ALT_SPDIST_QUERY = "alt";
    std::string const CH_SPDIST_QUERY = "ch";
    std::string const TNR_SPDIST_QUERY = "tnr";

    // kNN Query Methods
    std::string const INE_KNN_QUERY = "ine";
    std::string const GTREE_KNN_QUERY = "gtree";
    std::string const ROAD_KNN_QUERY = "road";
    std::string const SILC_KNN_QUERY = "silc";
    std::string const OPT_SILC_KNN_QUERY = "opt_silc";
    std::string const DISTBRWS_KNN_QUERY = "dist_brws";
    std::string const OPT_DISTBRWS_KNN_QUERY = "opt_db";
    std::string const DB_RTREE_KNN_QUERY = "db_rtree";
    std::string const IER_DIJKSTRA_KNN_QUERY = "ier_dijk";
    std::string const IER_GTREE_KNN_QUERY = "ier_gtree";
    std::string const IER_SILC_KNN_QUERY = "ier_silc";
    std::string const IER_PHL_KNN_QUERY = "ier_phl";
    std::string const IER_CH_KNN_QUERY = "ier_ch";
    std::string const IER_TNR_KNN_QUERY = "ier_tnr";

    // Steps for Experiments Command
    std::string const EXP_BUILD_INDEXES = "indexes";
    std::string const EXP_BUILD_OBJ_INDEXES = "objects";
    std::string const EXP_BUILD_RW_POI_OBJ_INDEXES = "rw_poi";
    std::string const EXP_RUN_KNN = "knn";
    std::string const EXP_RUN_KNN_OPTIMIZATIONS = "knn_special";
    std::string const EXP_RUN_SPD_PHL = "spd_phl";
    std::string const EXP_RUN_KNN_RW_POI = "knn_rw_poi";
    std::string const EXP_RUN_SP = "sp";
    std::string const EXP_RUN_SP_SINGLE = "sp_single";
        
    // Constants
    int const MAX_COMMAND_NAME_SIZE = 15;
    const unsigned int MAX_EDGES = 255; // Maximum number of edges allow due to use as link in SILC
    const int ROOT_PARENT_INDEX = -1;
    const int THRU_ROAD_ID = -1;
    const NodeID UNUSED_NODE_ID = 100000000; // Largest NodeID is size largst graph (i.e. US with 25 million nodes)
    
    // Other String Constants
    std::string const DISTANCE_WEIGHTS = "d";
    std::string const TIME_WEIGHTS = "t";
    std::string const NODE_EXT = "cnode";
    std::string const EDGE_EXT = "cedge";
}

#endif // _COMMON_H

