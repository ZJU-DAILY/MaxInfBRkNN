//
// Created by jins on 3/7/19.
//

#ifndef MAXBRGSTKNN_UTILITY_H

#include "querySet.h"
#include <vector>
#include <deque>
#include <map>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "results.h"

using namespace std;

#define FastArray	vector
#define FastList	deque
#define BitStore	vector<bool>

#define DEFAULT_CACHESIZE   4096   //4M
#define DEFAULT_BLOCKLENGTH	4096


//各类文本相似度度量函数
bool textCover(int term_id, vector<int> Pkey){
    float cnt = 0;
    for (size_t i = 0; i < Pkey.size(); i++) {
        if (Pkey[i] == term_id) {
            return  true;

        }
    }
    return false;
}

set<int> obtain_itersection_jins(set<int>& large, set<int>& small){
    set<int> inter;

    set<int>::iterator iter;
    set<int>::iterator iter2;
    int maxInSmall = -1;
    for(iter2=small.begin();iter2!=small.end();iter2++){
        int value = *iter2;
        maxInSmall = max(maxInSmall,value);
    }
    //printSetElements(small);


    for(iter=large.begin();iter!=large.end();iter++){
        int e = *iter;
        for(iter2=small.begin();iter2!=small.end();iter2++){
            int e2 = *iter2;
            if(e2==e){
                inter.insert(e);
                break;
            }
        }
        if(e>maxInSmall) break;
    }
    return  inter;
}


// used for compute Jaccard
struct triple {
    int id, first, second;
};

// see Jiaheng Lu's paper, stand for : key | uj | ij | 'uj | 'ij
struct fTuple {
    int id, first, second, third, fourth;
};

class maxList {
public:
    int id;
    float score;

    maxList(int a, float b) : id(a), score(b) {}

    void printresult() {
        printf("id: %d score is %f\n", id, score);
    }

    bool operator<(const maxList &a) const {
        return a.score > score;
    }
};

class updateList {
public:
    int id, checkNum;

    updateList(int a, int b) : id(a), checkNum(b) {}

    void printresult() {
        printf("id: %d checkNum is %d\n", id, checkNum);
    }

    bool operator<(const updateList &a) const {
        return a.checkNum > checkNum;
    }
};


// Jaccard similarity : measure the similarity of two keyword vectors
float jaccard(vector<int> Qkey, vector<int> Pkey) {
    float a = 0, b = 0, c = 0;
    vector<triple> triVec;
    for (size_t i = 0; i < Qkey.size(); i++) {
        bool flag = false;
        for (size_t j = 0; j < triVec.size(); j++) {
            if (Qkey[i] == triVec[j].id) {
                triVec[j].first++;
                flag = true;
                break;
            }
        }
        if (false == flag) {
            triple triTmp;
            triTmp.id = Qkey[i];
            triTmp.first = 1;
            triTmp.second = 0;
            triVec.push_back(triTmp);
        }
    }

    for (size_t i = 0; i < Pkey.size(); i++) {
        bool flag = false;
        for (size_t j = 0; j < triVec.size(); j++) {
            if (Pkey[i] == triVec[j].id) {
                triVec[j].second++;
                flag = true;
                break;
            }
        }
        if (false == flag) {
            triple triTmp;
            triTmp.id = Pkey[i];
            triTmp.first = 0;
            triTmp.second = 1;
            triVec.push_back(triTmp);
        }
    }
    for (size_t i = 0; i < triVec.size(); i++) {
        //cout<<triVec[i].id<<" "<<triVec[i].first<<' '<<triVec[i].second<<endl;
        a += triVec[i].first * triVec[i].second;
        b += triVec[i].first * triVec[i].first;
        c += triVec[i].second * triVec[i].second;
    }
    return a / (b + c - a);
}




vector<pair<int, int> > vecToPair(vector<int> keywords) {
    vector<pair<int, int> > result;
    for (size_t i = 0; i < keywords.size(); i++) {
        bool ifFind = false;
        for (size_t j = 0; j < result.size(); j++) {
            if (keywords[i] == result[j].first) {
                result[j].second++;
                ifFind = true;
                break;
            }
        }
        if (!ifFind) {
            result.push_back(make_pair(keywords[i], 1));
        }
    }
    return result;
}

int testVecPair()
{
    vector<int> vec;
    vec.push_back(0);   vec.push_back(0);   vec.push_back(1);   vec.push_back(1);   vec.push_back(1);
    vec.push_back(2);   vec.push_back(3);   vec.push_back(4);   vec.push_back(4);   vec.push_back(4);
    vector<pair<int,int>> vecPair = vecToPair(vec);
    for(size_t i=0;i<vecPair.size();i++)
    {
        cout<<vecPair[i].first<<' '<<vecPair[i].second<<endl;
    }
    return 0;
}

// compute max min simT bound
vector<fTuple>
fTuCompute(vector<pair<int, int> > userUKey, vector<pair<int, int> > userIKey, vector<pair<int, int> > objectUKey,
           vector<pair<int, int> > objectIKey) {
    vector<fTuple> fTuVec;

    for (size_t i = 0; i < userUKey.size(); i++) {
        bool flag = false;
        for (size_t j = 0; j < fTuVec.size(); j++) {
            if (userUKey[i].first == fTuVec[j].id) {
                fTuVec[j].first++;
                flag = true;
                break;
            }
        }
        if (false == flag) {
            fTuple fTuTmp;
            fTuTmp.id = userUKey[i].first;
            fTuTmp.first = userUKey[i].second;
            fTuTmp.second = 0;
            fTuTmp.third = 0;
            fTuTmp.fourth = 0;
            fTuVec.push_back(fTuTmp);
        }
    }

    for (size_t i = 0; i < userIKey.size(); i++) {
        bool flag = false;
        for (size_t j = 0; j < fTuVec.size(); j++) {
            if (userIKey[i].first == fTuVec[j].id) {
                fTuVec[j].second++;
                flag = true;
                break;
            }
        }
        if (false == flag) {
            fTuple fTuTmp;
            fTuTmp.id = userIKey[i].first;
            fTuTmp.first = 0;
            fTuTmp.second = userIKey[i].second;
            fTuTmp.third = 0;
            fTuTmp.fourth = 0;
            fTuVec.push_back(fTuTmp);
        }
    }

    for (size_t i = 0; i < objectUKey.size(); i++) {
        bool flag = false;
        for (size_t j = 0; j < fTuVec.size(); j++) {
            if (objectUKey[i].first == fTuVec[j].id) {
                fTuVec[j].third++;
                flag = true;
                break;
            }
        }
        if (false == flag) {
            fTuple fTuTmp;
            fTuTmp.id = objectUKey[i].first;
            fTuTmp.first = 0;
            fTuTmp.second = 0;
            fTuTmp.third = objectUKey[i].second;
            fTuTmp.fourth = 0;
            fTuVec.push_back(fTuTmp);
        }
    }

    for (size_t i = 0; i < objectIKey.size(); i++) {
        bool flag = false;
        for (size_t j = 0; j < fTuVec.size(); j++) {
            if (objectIKey[i].first == fTuVec[j].id) {
                fTuVec[j].fourth++;
                flag = true;
                break;
            }
        }
        if (false == flag) {
            fTuple fTuTmp;
            fTuTmp.id = objectIKey[i].first;
            fTuTmp.first = 0;
            fTuTmp.second = 0;
            fTuTmp.third = 0;
            fTuTmp.fourth = objectIKey[i].second;
            fTuVec.push_back(fTuTmp);
        }
    }
    return fTuVec;
}

// compute maxSimT score of U and O
float maxSimT(vector<fTuple> fTuVec) {
    float a = 0, b = 0, c = 0;

    for (size_t i = 0; i < fTuVec.size(); i++) {
        //cout<<triVec[i].id<<" "<<triVec[i].first<<' '<<triVec[i].second<<" "<<triVec[i].third<<' '<<triVec[i].fourth<<endl;
        if (fTuVec[i].second > fTuVec[i].third) {
            a += fTuVec[i].second * fTuVec[i].third;
            b += fTuVec[i].second * fTuVec[i].second;
            c += fTuVec[i].third * fTuVec[i].third;
        } else {
            if (fTuVec[i].first < fTuVec[i].fourth) {
                a += fTuVec[i].first * fTuVec[i].fourth;
                b += fTuVec[i].first * fTuVec[i].first;
                c += fTuVec[i].fourth * fTuVec[i].fourth;
            } else {
                if (fTuVec[i].fourth <= fTuVec[i].first && fTuVec[i].first <= fTuVec[i].third) {
                    a += fTuVec[i].first * fTuVec[i].first;
                    b += fTuVec[i].first * fTuVec[i].first;
                    c += fTuVec[i].first * fTuVec[i].first;
                } else {
                    a += fTuVec[i].third * fTuVec[i].third;
                    b += fTuVec[i].third * fTuVec[i].third;
                    c += fTuVec[i].third * fTuVec[i].third;
                }
            }
        }
    }
    return a / (b + c - a);
}

// compute minSimT score of U and O
float minSimT(vector<fTuple> fTuVec) {
    float a = 0, b = 0, c = 0;

    for (size_t i = 0; i < fTuVec.size(); i++) {
        //cout<<fTuVec[i].id<<" "<<fTuVec[i].first<<' '<<fTuVec[i].second<<" "<<fTuVec[i].third<<' '<<fTuVec[i].fourth<<endl;
        if (fTuVec[i].first * fTuVec[i].second >= fTuVec[i].third * fTuVec[i].fourth) {
            a += fTuVec[i].first * fTuVec[i].fourth;
            b += fTuVec[i].first * fTuVec[i].first;
            c += fTuVec[i].fourth * fTuVec[i].fourth;
        } else {
            a += fTuVec[i].second * fTuVec[i].third;
            b += fTuVec[i].second * fTuVec[i].second;
            c += fTuVec[i].third * fTuVec[i].third;
        }
    }
    return a / (b + c - a);
}


vector<vector<int> > conbinationSet(int n, int r)
{
    vector<vector<int> > conbinations;

    std::vector<bool> v(n);
    std::fill(v.end() - r, v.end(), true);

    do {
        vector<int> vecTmp;
        for (int i = 0; i < n; ++i) {
            if (v[i]) {
                vecTmp.push_back(i);
            }
        }
        conbinations.push_back(vecTmp);
    } while (std::next_permutation(v.begin(), v.end()));
    return conbinations;
}




#define _MIN(a, b) (((a) < (b))? (a) : (b)  )
#define _MAX(a, b) (((a) > (b))? (a) : (b)  )


inline int getKey(int i, int j) {
    int i0 = _MIN(i, j), j0 = _MAX(i, j);    // be careful of the order in other places
    return (i0 * VertexNum + j0);
}

inline void breakKey(int key, int &Ni, int &Nj) {
    Ni = key / VertexNum;
    Nj = key % VertexNum;
}


class QEntry {
public:
    int id;
    bool isvertex;
    int lca_pos;
    int dis;
    double score;
    double textual_score; double social_score;

    QEntry(int o, float flag, int pos, int d, double s, double t, double ss) : id(o), isvertex(flag), lca_pos(pos), dis (d),score(s){
        textual_score = t; social_score = ss;
    }
    QEntry(int o, float flag, int pos, int d, double s) : id(o), isvertex(flag), lca_pos(pos), dis (d),score(s){
    }
    ~QEntry(){
        //cout<<"nodeDistType被释放"<<endl;
    }

    void printRlt() {
        printf("adj: %d dist: %f\n", id, score);
    }

    bool operator<(const QEntry &a) const {
        return a.score > score;
    }
};


class CheckEntry{
public:
    int id;
    bool isUser;
    int b_th;
    User user;
    int node_id;
    set<int> keys_cover;
    float dist;


    CheckEntry(int check_id, int n_id, set<int> _keys, float d, int _b_th) : id(check_id), node_id(n_id), dist (d), b_th(_b_th){
        keys_cover = _keys;
        isUser = false;
    }

    CheckEntry(int check_id, int n_id, float d, User _user) : id(check_id), node_id(n_id), dist (d){
        user = _user;
        isUser = true;
    }

    ~CheckEntry(){
        //cout<<"nodeDistType被释放"<<endl;
    }

    void printRlt() {
        if(isUser)
            printf("u%d dist: %f\n", id, dist);
        else
            printf("border: v%d , 所在node: n %d, dist: %f\n", id, node_id, dist);
    }

    bool operator<(const CheckEntry &a) const {
        return a.dist < dist;
    }
};


class BatchCheckEntry{
public:
    int id;
    bool isUser;
    int b_th;
    User user;
    int node_id;
    set<int> keys_cover;
    float dist;
    int p_id;


    BatchCheckEntry(int check_id, int _p_id, int n_id, set<int> _keys, float d, int _b_th) : id(check_id), node_id(n_id), dist (d), b_th(_b_th){
        keys_cover = _keys;
        isUser = false;
        p_id = _p_id;
    }

    BatchCheckEntry(int check_id, int _p_id, int n_id, float d, User _user) : id(check_id), node_id(n_id), dist (d){
        user = _user;
        isUser = true;
        p_id = _p_id;
    }

    ~BatchCheckEntry(){
        //cout<<"nodeDistType被释放"<<endl;
    }

    void printRlt() {
        if(isUser)
            printf("u%d to p%d, dist(u,p): %f\n", id,p_id, dist);
        else
            printf("border: v%d , 所在node: n %d, to p%d, dist(b,p): %f\n", id, node_id, p_id, dist);
    }

    bool operator<(const BatchCheckEntry &a) const {
        return a.dist < dist;
    }
};


class BatchVerifyEntry {
public:
    int u_id;
    double score;
    double rk_current;
    vector<int> key;
    User u;

    BatchVerifyEntry(int id, double s, double r, User user) : u_id(id), score(s) {
        rk_current = r; u = user;
    }
    ~BatchVerifyEntry(){
        //cout<<"Result被释放"<<endl;
    }
    void printResult() {
        printf("id: %d score is %f\n", u_id, score);
        printf("key: ");

    }

    bool operator<(const BatchVerifyEntry &a) const {
        return a.u_id < u_id;
    }
};


class BatchCardinalityFisrt {
public:
    int usr_id;
    double max_score=0;
    double current_rk=0;
    priority_queue<ResultLargerFisrt> Lu;
    double gsk_score;
    User u;

    BatchCardinalityFisrt (int id, double lcl_rk, priority_queue<ResultLargerFisrt>& queue, User user)  {
        usr_id = id;
        u = user;
        current_rk = lcl_rk;
        Lu = queue;
        if(queue.top().score>0)
            max_score = queue.top().score;
    }

    

    ~BatchCardinalityFisrt (){
        //cout<<"Result被释放"<<endl;
    }
    void printResult() {

    }

    bool operator<(const BatchCardinalityFisrt  &a) const {
        return a.max_score > max_score;
    }
};


class POIHighestNode{

public:
    int p_id;
    int node_highest;

    POIHighestNode(int _p, int _n): p_id(_p), node_highest(_n){

    }
    ~POIHighestNode(){

    }

    void printRlt() {
        printf("p%d 's current highest is n %d\n", p_id, node_highest);

    }

    bool operator<(const POIHighestNode &a) const {
        return a.p_id < p_id;
    }


};


class SBEntry{
public:
    int p_id;

    float dist;

    double score_bound;


    SBEntry(int check_id, float _d, double _score) : p_id(check_id), dist(_d), score_bound (_score){

    }


    ~SBEntry(){
        //cout<<"nodeDistType被释放"<<endl;
    }

    void printRlt() {
        printf("p%d , score_lower: %f, dist: %f\n", p_id, score_bound,dist);
    }

    bool operator<(const SBEntry &a) const { //小的排前面
        return a.score_bound < score_bound;
    }
};


class GlobalSBEntry{
public:
    int term_id;

    int term_th;

    double sk_bound;


    GlobalSBEntry(int check_id, float _th, double _score) : term_id(check_id), term_th(_th), sk_bound(_score){

    }


    ~GlobalSBEntry(){
        //cout<<"nodeDistType被释放"<<endl;
    }


    bool operator<(const GlobalSBEntry &a) const {  //大的排前面
        return a.sk_bound > sk_bound;
    }
};



class nodeDistType {
public:
    int ID;
    float dist;
    bool isPOI;
    int inVertex;//(ni)
    int outVertex; //(Nj)
    float dis; //(distance to the outVertex)
    double ts_score;
    vector<int> keywords;
    //该拓展点为路网顶点
    nodeDistType(int a, float b) : ID(a), dist(b) {
        isPOI = false;
    }
    //该拓展点为POI
    nodeDistType(bool _isPOI, int a, float b, int _Ni, int _Nj,float _dis,double _score, vector<int> keys) : ID(a), dist(b) {
        isPOI = _isPOI;
        inVertex = _Ni;
        outVertex = _Nj;
        dis = _dis;
        ts_score = _score; keywords = keys;
    }
    ~nodeDistType(){
        //cout<<"nodeDistType被释放"<<endl;
    }

    void printRlt() {
        printf("adj: %d dist: %f\n", ID, dist);
    }

    void printPOI() {
        printf("p %d invertex: %d outVertex: %d\n", ID, inVertex,outVertex);
    }

    bool operator<(const nodeDistType &a) const {
        return a.dist < dist;
    }
};

enum NextExpansion{OutSide, Inside};
class nodeDistType_NVD {
public:
    int vertex_id;
    int poi_id;
    float dist;
    int term_id;
    bool isBorder = false;
    int b_th = -1;
    int  nextExpan = -1;


    //该拓展点中间条目信息
    nodeDistType_NVD(int _v, int _p, float cur_dis) {
        vertex_id = _v;
        poi_id = _p;
        dist = cur_dis;
    }

    nodeDistType_NVD(int v, int _p, float cur_dis, int _key) {
        poi_id = _p;
        term_id = _key;
        dist = cur_dis;
    }

    nodeDistType_NVD(int _v, int _p, float cur_dis, bool _isborder, int _expan) {
        vertex_id = _v;
        poi_id = _p;
        dist = cur_dis;
        isBorder = _isborder;
        nextExpan = _expan;
    }


    ~nodeDistType_NVD(){
        //cout<<"nodeDistType被释放"<<endl;
    }

    void printRlt() {
        printf("v: %d, dist: %f, poi:\n", vertex_id, dist, poi_id);
    }



    bool operator<(const nodeDistType_NVD &a) const {
        return a.dist < dist;
    }
};

class Entry_NVD {
public:
    int vertex_id;
    int poi_id;
    int dist_scale;
    float dist;


    //该拓展点中间条目信息
    Entry_NVD(int _v, int _p, float _dist, int cur_dis_scale) {
        vertex_id = _v;
        poi_id = _p;
        dist = _dist;
        dist_scale = cur_dis_scale;
    }



    bool operator<(const Entry_NVD &a) const {
        return a.dist_scale < dist_scale;
    }
};



class nodeDistType_POI {
public:
    float dist;
    int poi_id;
    vector<int> keywords;

    nodeDistType_POI(int _pid, float _dist, vector<int> keys) {
        poi_id = _pid;
        dist = _dist;
        keywords = keys;
    }



    ~nodeDistType_POI(){
        //cout<<"nodeDistType被释放"<<endl;
    }

    void printRlt() {
        printf("p: %d, dist: %f:\n", poi_id, dist);
    }



    bool operator<(const nodeDistType_POI &a) const {
        return a.dist < dist;
    }
};

class GlobalEntry_Pseudo {
public:

    float dist;
    int term_id;
    float score;

    GlobalEntry_Pseudo(int _tid, float cur_dist, float _score) {
        term_id = _tid;
        dist = cur_dist;
        score = _score;
    }



    ~GlobalEntry_Pseudo(){
        //cout<<"nodeDistType被释放"<<endl;
    }



    bool operator<(const GlobalEntry_Pseudo &a) const {
        return a.score > score;
    }
};


struct minmaxDis {
    float maxDis, minDis;
};

query transformToQ(User usr){
    query q;
    q.id= usr.id;
    q.dist = usr.dist;
    q.dis = usr.dis;
    q.Ni = usr.Ni;
    q.Nj = usr.Nj;
    q.keywords.assign(usr.keywords.begin(),usr.keywords.end());
    return q;

}



void printQueryKeywords(vector<int>  queryKeywords){
    //cout << "current queryKeywords:{" ;
    cout <<"{";
    for (size_t j = 0; j < queryKeywords.size(); j++) {
        cout << queryKeywords[j]<<",";
    }
    cout <<"}"<<endl;

}

void printCandidateSeeds(vector<int> seeds){
    cout <<"{";
    for (size_t j = 0; j < seeds.size(); j++) {
        cout << seeds[j]<<",";
        if(seeds[j]==30) {
            cout<<"find 30!"<<endl;
            getchar();
        }
    }
    cout <<"}"<<endl;
}

void printCandidateUser(vector<int> usrs){
    cout <<"{";
    for (size_t j = 0; j < usrs.size(); j++) {
        cout << usrs[j]<<",";
    }
    cout <<"}"<<endl;
}

void printSeeds(vector<int> seeds){
    cout <<"{";
    for (size_t j = 0; j < seeds.size(); j++) {
        cout << seeds[j]<<",";
    }
    cout <<"}"<<endl;
}

class Entry
{
public:
    int term;
    int count;
    Entry(int a, float b): term(a),count(b){}
    void printroad()
    {
        printf("term: %d count: %d\n",term,count);
    }
    bool operator<  (const Entry &a) const
    {
        return a.count > count;
    }
};

class EntryUsr
{
public:
    int id;
    float distance;
    EntryUsr(int a, float b): id(a),distance(b){}
    void printroad()
    {
        printf("u_id: %d distance: %d\n",id,distance);
    }
    bool operator<  (const EntryUsr &a) const
    {
        return a.id > id;
    }
};



class LCLEntry2
{
public:
    int id;
    int count;
    double score;
    LCLEntry2(int a, int b, double s): id(a),count(b),score(s){}
    void printLCLEntry2()
    {
        printf("id: %d, count: %d, score: %lf\n",id,count,score);
    }
    bool operator<  (const LCLEntry2 &entry) const
    {
        return entry.score < score;
    }
};

class Cardinarlity
{
public:
    int term;
    int cardinality;
    Cardinarlity(int a, float b): term(a),cardinality(b){}
    void printCardinarlity()
    {
        printf("term: %d count: %d\n",term,cardinality);
    }
    bool operator <  (const Cardinarlity &a) const
    {
        return a.cardinality > cardinality;
    }
};
class Margine
{
public:
    int u_id;
    double margine;
    Margine(int a, double b): u_id(a),margine(b){}
    void printCardinarlity()
    {
        printf("u: %d margine: %d\n",u_id,margine);
    }
    bool operator <  (const Margine &a) const
    {
        return a.margine > margine;
    }
};


struct Pair{
    int term;
    int usr;
};

class PairCardinarlity
{
public:
    int  pairCode;
    Pair pair;
    int cardinality;
    PairCardinarlity(Pair p, int code, int num): cardinality(num){
        pair.term = p.term;
        pair.usr = p.usr;
        pairCode = code;
    }
    void printCardinarlity()
    {
        printf("pair: t%d, u%d , count: %d\n",pair.term,pair.usr,cardinality);
    }
    bool operator <  (const PairCardinarlity &a) const
    {
        return a.cardinality > cardinality;
    }
};

class VerifDetail{
public :
        int usr_id;
        double social;
        double distance;
        double textual;
        double score;
        double rk;
        VerifDetail(int u, double a,double b,double c, double d, double value):usr_id(u), social(a),distance(b),textual(c),score(d),rk(value) {
        }
        void printVerifDetail(){

            //printf("u %d, top-k score =%lf, score(o,u)=%lf\n", usr_id,rk,score);
        }

    bool operator <  (const VerifDetail &a) const
    {
        return a.usr_id > usr_id;
    }

};




class KeySetCardinarlity
{
public:
    vector<int> keywords;
    vector<ResultDetail> elements;
    KeySetCardinarlity(vector<int> keys, vector<ResultDetail> usrs){
        keywords = keys;
        elements = usrs;
    }
    void printKeySetCardinarlity()
    {
        //printf("term: %d count: %f\n",term,cardinality);
    }
    bool operator <  (const KeySetCardinarlity &a) const
    {
        return a.elements.size() > elements.size();
    }
};


void printResults(vector<ResultDetail> results){
    priority_queue <ResultDetail> EQ;
    for(ResultDetail r: results){
        EQ.push(r);
    }
    while(!EQ.empty()){
        ResultDetail e = EQ.top();
        int usr_id = e.usr_id;
        double score = e.score;
        double rk = e.rk;
        //cout<<"u"<<usr_id<<" gsk(u,o)="<<score<<",topk="<<rk<<endl;
        cout<<"u"<<usr_id<<",";
        EQ.pop();
    }
    cout<<endl;

}

void printBatchRkGSKQResults(vector<int> stores, map<int,vector<ResultDetail>> results){
    for(int store: stores){
        priority_queue <ResultDetail> EQ;
        for(ResultDetail r: results[store]){
            EQ.push(r);
        }
        cout<<"store "<<store<<"的RkGSKQ, 结果个数= "<<results[store].size()<<"，内容："<<endl;
        while(!EQ.empty()){
            ResultDetail e = EQ.top();
            int usr_id = e.usr_id;
            double score = e.score;
            double rk = e.rk;
            //cout<<",u"<<usr_id<<" gsk(u,o)="<<score<<",topk="<<rk<<endl;
            cout<<"u"<<usr_id<<",";
            EQ.pop();
        }
        cout<<endl;
    }
    //cout<<endl;

}

void printMaxInfPOIInfo(vector<int> stores, map<int,vector<ResultDetail>> results){
    int optimal = 0; int optimal_id;
    for(int store: stores){
        priority_queue <ResultDetail> EQ;
        cout<<"store "<<store<<" inf = "<<results[store].size()<<",";//<<"，内容："<<endl;
        /*for(ResultDetail r: results[store]){
            EQ.push(r);
        }
        while(!EQ.empty()){
            ResultDetail e = EQ.top();
            int usr_id = e.usr_id;
            double score = e.score;
            double rk = e.rk;
            //cout<<",u"<<usr_id<<" gsk(u,o)="<<score<<",topk="<<rk<<endl;
            cout<<"u"<<usr_id<<",";
            EQ.pop();
        }
        cout<<endl;*/
        if(results[store].size()>optimal){
            optimal_id = store;
            optimal = results[store].size();
        }
    }

    cout<<"-----------------------------------------------------------------------"<<endl;
    cout<<"store "<<optimal_id<<" has the largest influence="<<results[optimal_id].size()<<"， 具体为：{";
    priority_queue <ResultDetail> EQ;
    for(ResultDetail r: results[optimal_id]){
           EQ.push(r);
    }
    while(!EQ.empty()){
           ResultDetail e = EQ.top();
           int usr_id = e.usr_id;
           double score = e.score;
           double rk = e.rk;
           //cout<<",u"<<usr_id<<" gsk(u,o)="<<score<<",topk="<<rk<<endl;
           cout<<"u"<<usr_id<<",";
           EQ.pop();
    }
    cout<<"}"<<endl;

}



void printElements(vector<int> elements){
    priority_queue <Entry> EQ;
    if(elements.size()>735376){
        cout<<elements.size()<<"elements.size()>5000"<<endl;
        ///getchar();
    } else{
        cout<<"elements.size="<<elements.size()<<":";
    }

    for(int ele: elements){
        Entry e(0,ele);
        EQ.push(e);
    }

    while(!EQ.empty()){
        int id = EQ.top().count;
        cout<<id<<",";
        EQ.pop();
    }
    cout<<endl;

}

void outputSeedResults(vector<int> elements, ofstream& output){
    priority_queue <Entry> EQ;
    if(elements.size()>5000){
        cout<<elements.size()<<"elements.size()>5000"<<endl;
        getchar();
    } else{
        cout<<"elements.size="<<elements.size()<<":";
    }

    for(int ele: elements){
        Entry e(0,ele);
        EQ.push(e);
    }

    while(!EQ.empty()){
        int id = EQ.top().count;
        cout<<id<<",";
        output<<"p"<<id<<",";
        EQ.pop();
    }
    cout<<endl; output<<endl;

}


void printPotentialUsers(vector<int> elements){
    priority_queue <Entry> EQ;
    if(elements.size()>5000){
        cout<<elements.size()<<"elements.size()>5000"<<endl;
        getchar();
    } else{
        cout<<"elements.size="<<elements.size()<<endl;
    }

    for(int ele: elements){
        Entry e(0,ele);
        EQ.push(e);
    }

    while(!EQ.empty()){
        int id = EQ.top().count;
        cout<<"u"<<id<<",";
        EQ.pop();
    }
    cout<<endl;

}



void printSetElements(set<int> elements){
    set<int>::iterator iter= elements.begin();
    while(iter!=elements.end()){
        cout<<*iter<<",";
        iter++;
    }

    cout<<endl;
}

void printSeedSet(set<int> elements){
    cout<<"seed size="<<elements.size()<<endl;
    set<int>::iterator iter= elements.begin();
    while(iter!=elements.end()){
        cout<<*iter<<",";
        iter++;
    }
    cout<<endl;
}

void outPutSeedSet(set<int> elements, ofstream& output){
    cout<<"seed size="<<elements.size()<<endl;
    output<<"seed size="<<elements.size()<<endl;
    set<int>::iterator iter= elements.begin();
    while(iter!=elements.end()){
        cout<<*iter<<",";
        output<<*iter<<",";
        iter++;
    }
    cout<<endl;
    output<<endl;
}


vector<int>  extractInvertedList( map <int,vector<int>> list){
    //cout<<"word count="<<list.size()<<endl;
    priority_queue <Entry> EQ;

    map<int, vector<int>>:: iterator iter;

    for(iter = list.begin(); iter!= list.end(); iter++){
        int term_id = iter ->first;
        vector<int> post_list = iter ->second;
        Entry e = Entry(term_id,post_list.size());
        EQ.push(e);
    }
    priority_queue <Entry> EQ2;
    vector<int> topfrequence;
    //int n = topk;
    while(!EQ.empty()){
        Entry e = EQ.top();
        topfrequence.push_back(e.term);
        EQ.pop();
        int term_id = e.term;
        int count = e.count;
        //cout<<term_id<<"frequence= "<<count<<", ";
        vector<int> post_list = list[term_id];
        /*for(int id: post_list)
            cout<<"u"<<id<<",";
        cout<<endl; */
    }
    cout<<endl;
    return topfrequence;


}

vector<int>  extract_Seed_Range( vector<int> seeds, int loc, double maxDist){
    //cout<<"word count="<<list.size()<<endl;
    vector<int> seeds_Range;
    for(int u: seeds){
        double distance = getDistance(loc,Users[u]);
        if(distance < maxDist) seeds_Range.push_back(u);
    }
    return seeds_Range;


}

vector<int>  extractTopkInvertedList( map <int,vector<int>> list, int topk){
    //cout<<"word count="<<list.size()<<endl;
    priority_queue <Entry> EQ;

    map<int, vector<int>>:: iterator iter;

    for(iter = list.begin(); iter!= list.end(); iter++){
        int term_id = iter ->first;
        vector<int> post_list = iter ->second;
        Entry e = Entry(term_id,post_list.size());
        EQ.push(e);
    }
    priority_queue <Entry> EQ2;
    vector<int> topfrequence;
    int n = topk;
    while(n){
        if(EQ.empty()) break;
        Entry e = EQ.top();
        topfrequence.push_back(e.term);
        EQ.pop();
        int term_id = e.term;
        int count = e.count;
        //cout<<term_id<<"frequence= "<<count<<", ";
        vector<int> post_list = list[term_id];
        /*for(int id: post_list)
            cout<<"u"<<id<<",";
        cout<<endl; */
        n--;

    }
    cout<<endl;
    return topfrequence;


}

vector<int>  extractInvertedListElements( map <int,vector<int>> list){
    //cout<<"word count="<<list.size()<<endl;
    vector<int> elements;

    map<int, vector<int>>:: iterator iter;

    for(iter = list.begin(); iter!= list.end(); iter++){
        int term_id = iter ->first;
        vector<int> post_list = iter ->second;
        for(int e:post_list)
            elements.push_back(e);
    }

    return elements;
}



//new add
typedef map<string,string> ConfigType;

typedef vector<Result> TopkResults;

typedef  map<int, vector<int>> InvList;


set<int> obtain_itersection(set<int>& set1, set<int>& set2){
    set<int> usr_UKey = set1;
    set<int> queryKey = set2;
    set<int> inter_Key;
    set_intersection(usr_UKey.begin(),usr_UKey.end(), queryKey.begin(),queryKey.end(),inserter(inter_Key,inter_Key.begin()));
    return  inter_Key;
}



class TopkQueryResult{
public:
    double topkScore;
    vector<Result> topkElements;

    TopkQueryResult(double score, priority_queue<Result> content) {
        topkScore = score;
        while(content.size()>0){
            Result r = content.top();
            topkElements.push_back(r);
            content.pop();
        }
    }

    TopkQueryResult(double score) {
        topkScore = score;
    }

    ~TopkQueryResult(){
        //cout<<"TopkQueryResult被释放"<<endl;
    }
    void printResults(int Qk){
        //cout<< "user"<< 28 << "'s top-k score="<<topkScore<<",results:"<<endl;
        //cout<< "user"<< i << "'s top-k results:"<<endl;
        int j =0;
        for(Result r: topkElements){
            cout<<"<top-"<<(Qk-j)<<",p"<<r.id<<","<<r.score<<">,";
            //cout<<"<top-"<<(Qk-j)<<",p"<<r.id<<","<<r.score<<">,";
            j++;
        }
        cout<<endl;

    }
    void printResults(int u_id, int Qk){
        cout<< "user"<< u_id << "'s top-k score="<<topkScore<<",results:"<<endl;
        //cout<< "user"<< i << "'s top-k results:"<<endl;
        int j =0;
        for(Result r: topkElements){
            cout<<"<top-"<<(Qk-j)<<",p"<<r.id<<","<<r.score<<">,";
            cout<<endl;
            //cout<<"<top-"<<(Qk-j)<<",p"<<r.id<<","<<r.score<<">,";
            j++;
        }
        cout<<endl;

    }
};

class TopkQueryCurrentResult{
public:
    double topkScore;
    vector<Result> topkElements;
    TopkQueryCurrentResult(float score) {
        topkScore = score;
    }
    TopkQueryCurrentResult(float score, priority_queue<Result> content) {
        topkScore = score;
        while(!content.empty()){
            Result r = content.top();
            topkElements.push_back(r);
            content.pop();
        }
    }
    ~TopkQueryCurrentResult(){
        //cout<<"TopkQueryResult被释放"<<endl;
    }
    void printResults(int Qk){
        //cout<< "user"<< 28 << "'s top-k score="<<topkScore<<",results:"<<endl;
        //cout<< "user"<< i << "'s top-k results:"<<endl;
        int j =0;
        for(Result r: topkElements){
            cout<<"<top-"<<(Qk-j)<<",o"<<r.id<<",score="<<r.score
                <<",r.dist="<<r.dist<<">,keywords:";
            printElements(r.key);
            j++;
        }
        //cout<<endl;

    }
};



//the data structure storeing the results of RkGSKQ
class SingleResults{   //凝练用户集结果，区别FilterResult过滤中间结果
public :
    map <int,vector<int>> invCandidateKeys;
    map <int,vector<int>> invCandidateSeeds;
    vector<ResultDetail> candidateUsr;
    SingleResults( map <int,vector<int>> a, map <int,vector<int>> b, vector<ResultDetail> c){
        invCandidateKeys = a;
        invCandidateSeeds = b;
        candidateUsr = c;
    }
    void printFilterContents(){
        cout<<"candidate usrs:";
        printResults(candidateUsr);
        cout<<"candidate keywords:";
        printQueryKeywords(extractInvertedList(invCandidateKeys));
        cout<<"candidate seed users:";
        printCandidateSeeds(extractInvertedList(invCandidateSeeds));

    }

};



void printTopkResultsQueue(priority_queue<Result> results){
    int i=results.size();
    while(results.size()>0){
        Result topk_element = results.top();
        results.pop();
        cout<<"<top"<<i<<":o"<<topk_element.id<<",score="<<topk_element.score<<", dist="<<topk_element.dist<<",keyword:";
            //"ts_score"<<topk_element.<<">";
        //<<topk_element.id<<" 位于leaf n"<<Nodes[topk_element.Ni].gtreepath.back()<<
        printElements(topk_element.key);
        i--;
    }
    //cout<<endl;
}

void printTopkResultsVector(TopkQueryCurrentResult& results){
    int i = 1;
    for( Result rr: results.topkElements){
        int o_id = rr.id;
        float score = rr.score;
        cout<<"<top"<<i<<":o"<<o_id<<",score="<<score<<">";
        i++;
    }
    cout<<endl;
}





#define MAXBRGSTKNN_UTILITY_H
#endif //MAXBRGSTKNN_UTILITY_H
