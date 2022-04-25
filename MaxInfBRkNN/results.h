//
// Created by jins on 5/18/20.
//

#ifndef MAXINFBRGSTKNN_RESULTS_H

#include <map>

using  namespace std;

//below is for the data structure of elements in TkGSKQ results
class ResultDetail{
public :
    int usr_id;
    double social;
    double distance;
    double textual;
    double score;
    double rk;
    ResultDetail(int u, double a,double b,double c, double d, double value):usr_id(u), social(a),distance(b),textual(c),score(d),rk(value) {
    }
    void printResultDetail(){

        //printf("u %d, top-k score =%lf, score(o,u)=%lf\n", usr_id,rk,gsk_score);
    }
    bool operator<  (const ResultDetail &a) const
    {
        return a.usr_id > usr_id;
    }

};



//多查询结果集合
typedef  map<int, vector<ResultDetail>> BRkGSKQResults;

typedef map<int, vector<ResultDetail>> BatchResults;

//仅仅存放 未经user层filtering处理的所有candidate users
typedef map<int, vector<ResultDetail>> FilterResults;

typedef map<int, vector<ResultDetail>> MQueryResults;


typedef  map<int, map<int,double>> PotentialUpper;

typedef pair<BRkGSKQResults,PotentialUpper> IncrementalResults;

//typedef map<int, priority_queue<ResultCurrent>> BatchResultsQueue;



#define MAXINFBRGSTKNN_RESULTS_H

#endif //MAXINFBRGSTKNN_RESULTS_H
