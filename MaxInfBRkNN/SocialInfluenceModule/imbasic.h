//
// Created by jin on 20-2-26.
//

#ifndef MAXINFBRGSTQ_IMBASIC_H
#define MAXINFBRGSTQ_IMBASIC_H

#include <unordered_map>
#include <vector>

enum PropModels {
    IC, LT,
};

typedef  PropModels  InfluModel;
//typedef map<int,set<int>>& AdjListSocial;
typedef unordered_map<int,vector<int>>& AdjListSocial;



#endif //MAXINFBRGSTQ_IMBASIC_H
