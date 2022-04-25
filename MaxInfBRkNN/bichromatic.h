//
// Created by jins on 10/31/19.
//

#ifndef MAXINFRKGSKQ_1_1_BICHROMATIC_H
#include <unordered_map>
#include <map>
#include <set>
#include <deque>
#include <stack>

typedef int NODE_ID;
typedef int TERM_ID;
typedef pair<NODE_ID,TERM_ID> Related_Pair;

class Top_result{
public:
    int id, Ni, Nj;
    double score;
    vector<int> key;

    Top_result(int i, double s) : id(i), score(s) {
        //key = keyInput;
    }
    ~Top_result(){
        //cout<<"Result2被释放"<<endl;
    }
    void print_topk_result() {
        printf("u_id: %d , score is %f\n", id, score);
        printf("key: ");
        for (size_t i = 0; i < key.size(); i++) {
            cout << key[i] << " ";
        }
        cout << endl;
    }

    bool operator<(const Top_result &a) const {
        return a.score < score;
    }
};


class Result {
public:
    int id, Ni, Nj;
    float dis, dist;double score;
    float textual_score; float social_score;
    vector<int> key;

    Result(int a, int b, int c, float e, float f, float s, vector<int> keyInput) : id(a), Ni(b), Nj(c), dis(e), dist(f),
                                                                                   score(s) {
        key = keyInput;

    }
    Result(int a, int b, int c, float e, float f, float s, float ts, float ss, vector<int> keyInput) : id(a), Ni(b), Nj(c), dis(e), dist(f),
                                                                                   score(s) {
        key = keyInput;
        textual_score = ts;
        social_score = ss;
    }
    ~Result(){
        //cout<<"Result被释放"<<endl;
    }
    void printResult() {
        printf("id: %d Ni: %d Nj: %d dis: %f dist: %f gsk_score is %f (ts=%f, social=%f) \n", id, Ni, Nj, dis, dist,
               score,textual_score,social_score);
        printf("key: ");
        for (size_t i = 0; i < key.size(); i++) {
            cout << key[i] << " ";
        }
        cout << endl;
    }

    bool operator<(const Result &a) const {  //小的在前
        return a.score < score;
    }
};

class ResultCurrent {
public:
    int o_id;
    double score;
    double score_upper;
    double relevance;
    double influence;

    ResultCurrent (int a, double f) : o_id(a), score(f) {

    }
    ResultCurrent (int a, double f, double u) : o_id(a), score(f) , score_upper(u){

    }
    ResultCurrent (int a, double f, double u,double r, double i) : o_id(a), score(f) , score_upper(u),relevance(r), influence(i){

    }
    ~ResultCurrent (){
        //cout<<"Result被释放"<<endl;
    }
    void printResult() {

    }

    bool operator<(const ResultCurrent  &a) const {
        return a.score < score;
    }
};

class ResultLargerFisrt {
public:
    int o_id;
    double score;
    double relevance;
    double influence;

    ResultLargerFisrt (int a, double f) : o_id(a), score(f) {

    }

    ResultLargerFisrt(int a, double f, double r, double i) : o_id(a), score(f) , relevance(r), influence(i){

    }
    ~ResultLargerFisrt (){
        //cout<<"Result被释放"<<endl;
    }
    void printResult() {

    }

    bool operator<(const ResultLargerFisrt  &a) const {
        return a.score > score;
    }
};

class CardinalityFisrt {
public:
    int usr_id;
    double max_score=0;
    double current_rk=0;
    priority_queue<ResultLargerFisrt> Lu;

    CardinalityFisrt (int id, double lcl_rt, priority_queue<ResultLargerFisrt>& queue)  {
        usr_id = id;
        current_rk = lcl_rt;
        Lu = queue;
        if(queue.top().score>0)
            max_score = queue.top().score;
    }

    ~CardinalityFisrt (){
        //cout<<"Result被释放"<<endl;
    }
    void printResult() {

    }

    bool operator<(const CardinalityFisrt  &a) const {
        return a.max_score > max_score;
    }
};

class CardinalityEntry {
public:
    int id;
    int size;


    CardinalityEntry (int _id, int _size)  {
       id = _id;
       size = _size;
    }

    ~CardinalityEntry (){
        //cout<<"Result被释放"<<endl;
    }
    void printResult() {

    }

    bool operator<(const CardinalityEntry  &a) const {
        return a.size > size;
    }
};



class VerifyEntry {
public:
    int u_id;
    double score;
    double rk_current;
    vector<int> key;

    VerifyEntry(int id, double s, double r) : u_id(id), score(s) {
        rk_current = r;
    }
    ~VerifyEntry(){
        //cout<<"Result被释放"<<endl;
    }
    void printResult() {
        printf("id: %d score is %f\n", u_id, score);
        printf("key: ");

    }

    bool operator<(const VerifyEntry &a) const {
        return a.u_id < u_id;
    }
};




class ULeafRefinePriorEntry {
public:
    int leaf_id;
    int priority;
    vector<int> interKey;   //该Ui与查询关键词集的交集

    set<int> usr_withinLeaf;
    set<int> ossocialte_query;

    ULeafRefinePriorEntry(int u_leaf, int r, vector<int> keys,set<int> us, set<int> os) : leaf_id(u_leaf), priority(r) {
        interKey = keys;
        usr_withinLeaf = us;
        ossocialte_query = os;
    }
    ~ULeafRefinePriorEntry(){
        //cout<<"Result被释放"<<endl;
    }
    void printResult() {
        printf("u_leaf: %d, 优先级权重: %d, 相关用户个数： %d, 关联查询对象数： %d\n", leaf_id, priority, usr_withinLeaf.size(),ossocialte_query.size());

    }

    bool operator<(const ULeafRefinePriorEntry &a) const {
        return a.priority > priority;
    }
};

typedef priority_queue<ULeafRefinePriorEntry> Priority_ULeafRefinement;




//双色体对象数据结构
struct User {
    int id;
    int Ni;
    int Nj;
    float dist;
    float dis;
    vector<int> keywords;
    set <int> keywordSet;
    vector<int> friends;
    bool topK = false;
    bool verified = false;
    double topkScore_current = -1.0;
    double topkScore_Final = -1.0;

    priority_queue<Top_result> resultFinal;
    priority_queue<Top_result> current_Results;
    bool isVisited = false;
    double u_lcl=0;
    priority_queue<ResultCurrent> u_LCL;   //下界列表
    int addressDataOrgByLeaf = 0;
};
struct POI {
    int id=-100;
    int Ni=-1;
    int Nj=-1;
    int category; //0 stand for user and 1 stand for object
    float dist;
    float dis;
    vector<int> keywords;
    set <int> keywordSet;
    vector<int> check_ins;  ////签到用户id
    float maxPopular = -1.0;
    bool isVisited = false;
    int  addressDataOrgByLeaf = 0;
};

void printPOIInfo(POI& p){
    cout<<"poi"<<p.id<<", Ni="<<p.Ni<<",Nj="<<p.Nj;
    printf(", keywords(size=%d): ",p.keywords.size());
    for(int term:p.keywords)
        cout<<term<<",";
    cout<<endl;
    printf(", check-in(size=%d): ",p.check_ins.size());
    for(int checkin_uid:p.check_ins)
        cout<<checkin_uid<<",";
    cout<<endl;
}


void printUsrInfo(User& u){
    cout<<"u"<<u.id<<", Ni="<<u.Ni<<",Nj="<<u.Nj;
    cout<<", keywords: ";
    for(int term:u.keywords)
        cout<<term<<",";
    cout<<endl;
    /*cout<<"friends: ";
    for(int f:u.friends)
        cout<<"u"<<f<<",";
    cout<<endl;*/
}




double maxIDFScore = 0;


struct matrix {
    int x;
    int y;
    float value;
};


typedef unordered_map<long, float> socialType;
typedef unordered_map<int, vector<int> > friType;
typedef unordered_map<int, vector<pair<int, int>>> checkType;
typedef unordered_map<int, int> userLeafType;
typedef unordered_map<int, vector<int>> userType;
typedef unordered_map<int, float> maxSocialScoreType;

socialType social; // key: user*poi_num + poi
//userType userMap,poiMap; // save user and poi information, key is id
friType friendMap; // key: userId, then it's friends' id
checkType checkMap; // key: poiId, then User,checkNum
checkType usercheckMap; // key: userId, then poi,checkNum
userLeafType userLeafMap;
userType userMap;
maxSocialScoreType maxSocialScore; // key: userId, then the max social score of user
vector<matrix> nei_poi;



#define MAXINFRKGSKQ_1_1_BICHROMATIC_H

#endif //MAXINFRKGSKQ_1_1_BICHROMATIC_H
