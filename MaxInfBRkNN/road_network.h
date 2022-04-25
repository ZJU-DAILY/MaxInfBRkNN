//
// Created by jins on 10/31/19.
//

#ifndef MAXINFRKGSKQ_1_1_ROAD_NETWORK_H
#include <unordered_map>
#include <map>
#include <set>
#include <deque>
#include <stack>
//#include "gim_tree.h"

//路网信息
struct edge {
    int Ni, Nj;
    float dist;
    vector<POI> pts;
    vector<User> usrs;
    set<int> OUnionKeys;
    int poiAddress =-1111;
    int poiIdxAddress =-1111;
    int okeyAddress = -1;
};
typedef map<int, edge> EdgeMapType;
typedef map<int, map<int, float> > LeafMsType;



//distance bound between nodes
typedef map<int,map<int,float>> DistanceBoundMap;
//DistanceBoundMap MaxDisBound;
//DistanceBoundMap MinDisBound;

float **MaxDisBound;
float **MinDisBound;



typedef vector<vector<int> > AdjType;
typedef vector<vector<float> > AdjEdgeWeight;

EdgeMapType EdgeMap;    // key: i0*NodeNum+j0
LeafMsType LeafMs; //key: user_leaf poi_leaf
AdjType adjList;
AdjEdgeWeight adjWList;


#define MAXINFRKGSKQ_1_1_ROAD_NETWORK_H

#endif //MAXINFRKGSKQ_1_1_ROAD_NETWORK_H
