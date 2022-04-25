#ifndef MAXBRGSTKNN_GTREE_H
#define MAXBRGSTKNN_GTREE_H


#include <stdio.h>
#include <metis.h>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <unordered_map>
#include <map>
#include <set>
#include <deque>
#include <stack>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include "config.h"
#include "bichromatic.h"
#include "commonStruct.h"
#include "omp.h"

using namespace std;

#define maximum  99999.9
#define minimum  0.0
struct timeval pv2;
long long ps2, pe2;
#define PRE_START gettimeofday( &pv2, NULL ); ps2 = pv2.tv_sec * 1000000 + pv2.tv_usec;
#define PRE_END gettimeofday( &pv2, NULL ); pe2 = pv2.tv_sec * 1000000 + pv2.tv_usec;
#define PRE_TICK_PRINT(T) printf("%s runtime: %lld (us)\r\n", (#T), (pe2 - ps2) );

//maxinf总共用时
struct timeval ov;
long long os, oe;


//MRkGSKQ求解
struct timeval mv;
long long ms, me;



//总共用时
struct timeval tv;
long long ts, te;

//加载数据用时
//总共用时
struct timeval lv;
long long ls, le;

//group用时
struct timeval gv;
long long gs, ge;
//batch用时
struct timeval bv;
long long bs, be;
//filter用时
struct timeval fv;
long long fs, fe;
struct timeval pv;
long long ps, pe;
//refine用时
struct timeval rv;
long long rs, re;

struct timeval testv;
long long test_s, test_e;


#define OVERALL_START gettimeofday( &ov, NULL ); os = ov.tv_sec * 1000000 + ov.tv_usec;
#define OVERALL_END gettimeofday( &ov, NULL ); oe = ov.tv_sec * 1000000 + ov.tv_usec;
#define OVERALL_PRINT(T) printf("%s runtime: %lld (ms)\r\n", (#T), (oe - os)/1000 );

#define MRkGSKQ_START gettimeofday( &mv, NULL ); ms = mv.tv_sec * 1000000 + mv.tv_usec;
#define MRkGSKQ_END gettimeofday( &mv, NULL ); me = mv.tv_sec * 1000000 + mv.tv_usec;
#define MRkGSKQ_PRINT(T) printf("%s runtime: %lld (ms)\r\n", (#T), (me - ms)/1000 );


#define TEST_START gettimeofday( &testv, NULL ); test_s = testv.tv_sec * 1000000 + testv.tv_usec;
#define TEST_END gettimeofday( &testv, NULL ); test_e = testv.tv_sec * 1000000 + testv.tv_usec;
#define TEST_DURA_PRINT(T) printf("%s runtime: %lld (sec)\r\n", (#T), (test_e - test_s)/1000000 );


#define TIME_TICK_START gettimeofday( &tv, NULL ); ts = tv.tv_sec * 1000000 + tv.tv_usec;
#define TIME_TICK_END gettimeofday( &tv, NULL ); te = tv.tv_sec * 1000000 + tv.tv_usec;
#define TIME_TICK_PRINT(T) printf("%s runtime: %lld (ms)\r\n", (#T), (te - ts)/1000 );
#define TIME_TICK_PRINT_us(T) printf("%s runtime: %lld (us)\r\n", (#T), (te - ts) );




#define Group_START gettimeofday( &gv, NULL ); gs = gv.tv_sec * 1000000 + gv.tv_usec;
#define Group_END gettimeofday( &gv, NULL ); ge = gv.tv_sec * 1000000 + gv.tv_usec;
#define Group_PRINT(T) printf("%s runtime: %lld (ms)\r\n", (#T), (ge - gs)/1000 );

#define Batch_START gettimeofday( &bv, NULL ); bs = bv.tv_sec * 1000000 + bv.tv_usec;
#define Batch_END gettimeofday( &bv, NULL); be = bv.tv_sec * 1000000 + bv.tv_usec;
#define Batch_PRINT(T) printf("%s runtime: %lld (ms)\r\n", (#T), (be - bs)/1000 );
#define Filter_START gettimeofday( &fv, NULL ); fs = fv.tv_sec * 1000000 + fv.tv_usec;
#define PAUSE_START gettimeofday( &pv, NULL ); ps = pv.tv_sec * 1000000 + pv.tv_usec;
#define PAUSE_END gettimeofday( &pv, NULL ); pe = pv.tv_sec * 1000000 + pv.tv_usec;
//#define Filter_START gettimeofday( &pv, NULL ); fs = fv.tv_sec * 1000000 + fv.tv_usec;
#define Filter_END gettimeofday( &fv, NULL ); fe = fv.tv_sec * 1000000 + fv.tv_usec;
#define Filter_PRINT(T) printf("%s runtime: %lld (ms)\r\n", (#T), (fe - fs)/1000 );
#define Filter_PAUSE_ PRINT(T) printf("%s runtime: %lld (ms)\r\n", (#T), (fe - pe + ps -fs)/1000 );
#define Refine_START gettimeofday( &rv, NULL ); rs = rv.tv_sec * 1000000 + rv.tv_usec;
#define Refine_END gettimeofday( &rv, NULL ); re = rv.tv_sec * 1000000 + rv.tv_usec;
#define Refine_PRINT(T) printf("%s runtime: %lld (ms)\r\n", (#T), (re - rs)/1000 );




//-------------------用于GIM-tree上的剪枝
// this is for pruning
class LCLEntry
{
public:
    int id;
    int count;
    double score;
    LCLEntry(int a, double s, int c): id(a),score(s),count(c){}
    void printLCLEntry()
    {
        printf("id: %d, count: %d, score: %lf\n",id,count,score);
    }
    bool operator<  (const LCLEntry &entry) const
    {
        return entry.score < score;
    }

};

class UCLEntry
{
public:
    int id;
    int count;
    double score;
    UCLEntry(int a, double s, int c): id(a),score(s),count(c){}
    void printLCLEntry()
    {
        printf("id: %d, count: %d, score: %lf\n",id,count,score);
    }
    bool operator<  (const UCLEntry &entry) const
    {
        return entry.score > score;
    }
};


typedef struct {
    priority_queue<LCLEntry> LCL;
    priority_queue<UCLEntry> UCL;
    double topscore = -1;
}LCLResult;

void printLCLResultsInfo(LCLResult& lcl){
    LCLResult _lcl = lcl;
    while (_lcl.LCL.size()>0) {
        LCLEntry e = _lcl.LCL.top();
        e.printLCLEntry();
        _lcl.LCL.pop();
    }
}


typedef struct {
    priority_queue<UCLEntry> UCL;
    double topscore = 0;
}UCLResult;


typedef vector<LCLResult>  MultiLCLResult;

map<int,map<int,float>> SPDCache;



//路网顶点信息（包含其所在Gtree node的id）
typedef struct {
    double x, y;
    vector<int> adjnodes;
    vector<int> adjweight;
    bool isborder;
    vector<int> gtreepath; // this is used to do sub-graph locating
} Node;

//Gtree节点
typedef struct {
    vector<int> borders;
    vector<int> children;
    bool isleaf;
    //若为leaf node则
    int poiInvIndexBaseAddress;  //叶节点对应的兴趣点倒排索引数据块内容的基地址  //jordan


    vector<int> leafnodes;
    int father;
// ----- min dis -----
    vector<int> union_borders; // for non leaf node
    vector<int> mind; // min dis, row by row of union_borders
// ----- for pre query init, OCCURENCE LIST in paper -----
    vector<int> nonleafinvlist;
    vector<int> leafinvlist;
    vector<int> up_pos;
    vector<int> current_pos;

    // add new figure
    vector<float> maxUserScore;
    vector<float> minUserScore;

    vector<int> userUKey;
    vector<int> userIKey;

    vector<int> poiSet;
    vector<int> userSet;



    //add new feature by jin
    int minPopularPOI_Count =10000;
    float maxEdgeDis = minimum;
    float minEdgeDis = maximum;
    float minDis = maximum;
    float maxDis = 0;
    set<int> unionFriends;
    set<int> unionFriends_load;
    set<int> userUKeySet;
    set<int> objectUKeySet;
    set<int> node_poiCnt;

    set<int> vetexSet;
    double LCLRt;
    map<int,double> termLCLRt;
    map<int,LCLResult> termLCLDetail;
    map<int,LCLResult> termParentLCLDetail;
    map<set<int>,MultiLCLResult> cacheMultiLCLDetail;

    LCLResult  nodeLCL;
    LCLResult  parentLCL;

    MultiLCLResult nodeLCL_multi;   // this is for leaf node

    bool hasLCLRt = false;
    //add for batch RkGSKQ
    set<int> child_hasStores;
    set<int> stores;
    map<int,set<int>> related_queryNodes;   //(term_id,{node1...,node2})

    set<int> ossociate_queryNodes;   //(term_id,{node1...,node2})

    //set<int> occurrence_list;   //({node1(poi1)...,node2(poi2)})
    unordered_map<int, set<int>> inverted_list_o;    //(term_id, {entry1, entry2, entry3})
    unordered_map<int, set<int>> inverted_list_u;

    //用于记录该节点中覆盖某关键词（term)的所有对象的集合
    unordered_map<int, set<int>> term_object_Map;    //(term_id, {o_id1, o_id2, o_id3...})
    unordered_map<int, int> term_object_countMap;   //(term_id, related object size)    // this is for lower counting list computation
    unordered_map<int, set<int>> term_usr_Map;    //(term_id, {u_id1, u_id2, u_id3...}) // facilitate aggressive pruning strategy 1
    unordered_map<int, int> term_usr_countMap;   //(term_id, related user size) // for aggressive pruning strategy 1
    unordered_map<int, set<int>> term_usrLeaf_Map;  //(term_id, {leaf_id1, leaf_id2, leaf_id3...}) // for aggressive pruning strategy 1



    map<int, int> term_poiCnt;
    map<int, int> term_usrCnt;

    //term_id对应的的 poi child entries （兴趣点）
    unordered_map<int, vector<int>> invertedAddress_list_o;  // term_id, {child1,....}

    //term_id对应的的 user child entries （用户）
    unordered_map<int, vector<int>> invertedAddress_list_u;  // term_id, {child1,....}



    // vector[0] 为 该 object node在 term_id对应的 inverted list中的内部偏移位置
    map<int, vector<int>> term_innerAddr_OMap;      //term_id,  vector[0] inner address 表
    // vector[0] 为 该 usr node在 term_id 对应的...
    map<int, vector<int>> term_innerAddr_UMap;


    int maxSocialCount = -1;
    int minSocialCount = INT32_MAX;
    //从磁盘中读取的数据缓存
    int maxSocialCount_load = -1;
    int minSocialCount_load = INT32_MAX;


} TreeNode;

void printONodeInfo(TreeNode& n, int n_id){
    cout<<"onode"<<n_id<<", isleaf="<<n.isleaf<<",father="<<n.father<<endl;
    cout<<"oVocabulary=";
    for(int keyword: n.objectUKeySet)
        cout<<"<t"<<keyword<<","<<n.term_innerAddr_OMap[keyword][0]<<">,";
    cout<<endl;

}

void printUNodeInfo(TreeNode& n, int n_id){
    cout<<"unode"<<n_id<<", isleaf="<<n.isleaf<<",father="<<n.father<<endl;
    for(int keyword: n.userUKeySet)
        cout<<"<t"<<keyword<<","<<n.term_innerAddr_UMap[keyword][0]<<">,";
    cout<<endl;

}



//双色体对象集合
vector<POI> POIs;
vector<User> Users;
vector<float> termWeight_Map;

//更换数据集后，需修改
//int UserSize =1001;
//int POISize = 13467;

//统计user、poi的个数
int usrCnt = 0;
int poiCnt = 0;

//倒排索引信息
map<int, vector<int>> invListOfUser;
//term_ID, vector<poi_id>,  the inverted list of each keyword
map<int, vector<int>> invListOfPOI;



//数据集社交信息
unordered_map<int,vector<int>> friendshipMap;   // user,  follower

unordered_map<int,vector<int>> followedMap;   // user,  influencer

unordered_map<int,unordered_map<int,int>> friendShipTable;
//vector<vector<int>> friendShipTable;
//社交网络图
Graph graph;



//数据集信息索引
//indexing usrer textual information
map<int,map<int,vector<int>>> leafUsrInv;  // term_id, map<leaf_id, vector<usr_id>>
map<int,set<int>> usrTerm_leafSet;  //term_id(from usr), vector<leaf_id>
// indexing node textual information of users
map<int,map<int,set<int>>> usrNodeTermChild; // node_id, map<term_id, vector<child_id>>

//indexing poi textual information
map<int,set<int>> term_poiSetInv;
map<int,map<int,vector<int>>> leafPoiInv;  // term_id, map<leaf_id, vector<poi_id>>

map<int,set<int>> poiTerm_leafSet; // term_id(from poi), set<leaf_id>


// indexing node textual information of pois
map<int,map<int,set<int>>> poiNodeTermChild; // node_id, map<term_id, vector<child_id>>


unordered_map<int,vector<int>> userCheckInIDList;  //usr_id, vector<poi_id>
unordered_map<int,vector<int>> userCheckInCountList;  //usr_id, vector<count>
unordered_map<int,vector<pair<int,int>>> userCheckInfoList;  //usr_id, {pair<poi_id, count>....}

unordered_map<int,vector<int>> poiCheckInIDList;  //poi_id, vector<usr_id>
unordered_map<int,vector<int>> poiCheckInCountList;  //poi_id, vector<count>
unordered_map<int,vector<pair<int,int>>> poiCheckInfoList;  //poi_id, {pair<user_id, count>....}



map<int,map<int,int>> userCheckInMap;  // usr_id,  map<poi_id,checkin-count>

map<int,map<int,int>> poiCheckInMap;  // poi_id,  map<usr_id,checkin-count>

map<int,double> usrMaxSS;       // usr_id, 该用户的朋友在所有兴趣点中的最大签到次数

map<int,int> usrOptimalPOI;



map<int,map<int,double>> LCLMap;  // usrNode_id,  map<term_id,Rt>

//for social bound
map<int,map<int,int>> term_leaf_MinPopular;    // term_id, map<leaf_id, minPopular_Count>
map<int,map<int,int>> node_term_MinPopular;    // node_id, map<term_id, minPopular_Count>


// for texual bound of leaf node
map<int,map<int,set<int>>> term_leaf_KeySet;    // term_id, map<leaf_id, set<term_id>>




int noe; // number of edges
vector<Node> Nodes;
bool IsBorder[VertexNum];
vector<TreeNode> GTree;
vector<TreeNode> GTreeCache;
int id2idx_vertexInLeaf[VertexNum];


// use for metis
// idx_t = int64_t / real_t = double
idx_t nvtxs; // |vertices|
idx_t ncon; // number of weight per vertex
idx_t *xadj; // array of adjacency of indices
idx_t *adjncy; // array of adjacency nodes
idx_t *vwgt; // array of weight of nodes
idx_t *adjwgt; // array of weight of edges in adjncy
idx_t nparts; // number of parts to partition
idx_t objval; // edge cut for partitioning solution
idx_t *part; // array of partition vector
idx_t options[METIS_NOPTIONS]; // option array

// METIS setting options
void options_setting() {
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY; // _RB
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // _VOL
    options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM; // _RM
    options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_RANDOM; // _GROW _EDGE _NODE
    options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM; // _GREEDY _SEP2SIDED _SEP1SIDED

    options[METIS_OPTION_UFACTOR] = 500;
    // options[METIS_OPTION_MINCONN];
    options[METIS_OPTION_CONTIG] = 1;
    // options[METIS_OPTION_SEED];
    options[METIS_OPTION_NUMBERING] = 0;
    // options[METIS_OPTION_DBGLVL] = 0;
}


// input init
void load_road_input() {
    FILE *fin;
    // load node
    std::string node_path = getRoadInputPath(NODE);
    printf("LOADING NODE...");
    fin = fopen(node_path.c_str(), "r");
    int nid;
    double x, y;
    while (fscanf(fin, "%d %lf %lf", &nid, &x, &y) == 3) {
        Node node = {x, y};
        Nodes.push_back(node);
    }
    fclose(fin);
    printf("COMPLETE. NODE_COUNT=%d\n", (int) Nodes.size());

    // load edge
    printf("LOADING EDGE...");

    std::string edge_path = getRoadInputPath(EDGE);
    fin = fopen(edge_path.c_str(), "r");
    int eid=0;
    int snid=0, enid=0;
    int weight=0;
    //float distance=0;
    noe = 0;
    while (fscanf(fin, "%d %d %d %d", &eid, &snid, &enid, &weight) == 4) {
        noe++;
        //iweight = (int) (weight * WEIGHT_INFLATE_FACTOR);
        Nodes[snid].adjnodes.push_back(enid);
        Nodes[snid].adjweight.push_back(weight);
        Nodes[enid].adjnodes.push_back(snid);
        Nodes[enid].adjweight.push_back(weight);
    }
    fclose(fin);
    printf("COMPLETE.\n");
}

void init() {
    load_road_input();
    //cout<<"******"<<endl;
    options_setting();
}

void finalize() {
    delete xadj;
    delete adjncy;
    delete adjwgt;
    delete part;
}

// dijkstra search, used for single-source shortest path search WITHIN one gtree leaf node!
// input: s = source node
//        cands = candidate node list
//        graph = search graph(this can be set to subgraph)
vector<int> dijkstra_candidate(int s, vector<int> &cands, vector<Node> &graph) {
    // init
    set<int> todo;
    todo.clear();
    todo.insert(cands.begin(), cands.end());

    unordered_map<int, int> result1;
    result1.clear();
    set<int> visited;
    visited.clear();
    unordered_map<int, int> q;
    q.clear();
    q[s] = 0;

    // start
    int _min, minpos, adjnode, weight;
    while (!todo.empty() && !q.empty()) {
        _min = -1;
        for (unordered_map<int, int>::iterator it = q.begin(); it != q.end(); it++) {
            if (_min == -1) {
                minpos = it->first;
                _min = it->second;
            } else {
                if (it->second < _min) {
                    _min = it->second;
                    minpos = it->first;
                }
            }
        }

        // put min to result1, add to visited
        result1[minpos] = _min;
        visited.insert(minpos);
        q.erase(minpos);

        if (todo.find(minpos) != todo.end()) {
            todo.erase(minpos);
        }

        // expand
        for (int i = 0; i < graph[minpos].adjnodes.size(); i++) {
            adjnode = graph[minpos].adjnodes[i];
            if (visited.find(adjnode) != visited.end()) {
                continue;
            }
            weight = graph[minpos].adjweight[i];

            if (q.find(adjnode) != q.end()) {
                if (_min + weight < q[adjnode]) {
                    q[adjnode] = _min + weight;
                }
            } else {
                q[adjnode] = _min + weight;
            }

        }
    }

    // output
    vector<int> output;
    for (int i = 0; i < cands.size(); i++) {
        output.push_back(result1[cands[i]]);
    }

    // return
    return output;
}


/*----------------------------构建Gtree索引部分的功能函数(不正确)--------------------------------------*/


// transform original data format to that suitable for METIS
void data_transform_init(set<int> &nset) {
    // nvtxs, ncon
    nvtxs = nset.size();
    ncon = 1;

    xadj = new idx_t[nset.size() + 1];
    adjncy = new idx_t[noe * 2];
    adjwgt = new idx_t[noe * 2];


    int xadj_pos = 1;
    int xadj_accum = 0;
    int adjncy_pos = 0;

    // xadj, adjncy, adjwgt
    unordered_map<int, int> nodemap;
    nodemap.clear();

    xadj[0] = 0;
    int i = 0;
    for (set<int>::iterator it = nset.begin(); it != nset.end(); it++, i++) {
        // init node map
        nodemap[*it] = i;

        int nid = *it;
        int fanout = Nodes[nid].adjnodes.size();
        for (int j = 0; j < fanout; j++) {
            int enid = Nodes[nid].adjnodes[j];
            // ensure edges within
            if (nset.find(enid) != nset.end()) {
                xadj_accum++;

                adjncy[adjncy_pos] = enid;
                adjwgt[adjncy_pos] = Nodes[nid].adjweight[j];
                adjncy_pos++;
            }
        }
        xadj[xadj_pos++] = xadj_accum;
    }

    // adjust nodes number started by 0
    for (int i = 0; i < adjncy_pos; i++) {
        adjncy[i] = nodemap[adjncy[i]];
    }

    // adjwgt -> 1
    if (ADJWEIGHT_SET_TO_ALL_ONE) {
        for (int i = 0; i < adjncy_pos; i++) {
            adjwgt[i] = 1;
        }
    }

    // nparts
    nparts = PARTITION_PART;

    // part
    part = new idx_t[nset.size()];
}

// graph partition
// input: nset = a set of node id
// output: <node, node belong to partition id>
unordered_map<int, int> graph_partition(set<int> &nset) {
    unordered_map<int, int> result;

    // transform data to metis
    data_transform_init(nset);

    // partition, result -> part
    // k way partition
    METIS_PartGraphKway(
            &nvtxs,
            &ncon,
            xadj,
            adjncy,
            NULL,
            NULL,
            adjwgt,
            &nparts,
            NULL,
            NULL,
            options,
            &objval,
            part
    );

    // push to result
    result.clear();
    int i = 0;
    for (set<int>::iterator it = nset.begin(); it != nset.end(); it++, i++) {
        result[*it] = part[i];
    }

    // finalize
    finalize();

    return result;
}

// init status struct
typedef struct {
    int tnid; // tree node id
    set<int> nset; // node set
} Status;


int generateObject() {
    //stringstream ss;
    //ss<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name<<".object";
    string object_path = getRoadInputPath(OBJECTS);
    ofstream outfile;
    outfile.open(object_path.c_str());
    for (int i = 0; i < VertexNum; i++) {
        if (i != VertexNum - 1)
            outfile << i << " " << i << endl;
        else
            outfile << i << " " << i;
    }
    outfile.close();

    return 0;
}


// gtree construction
void buildTreeStructure() {
    // init root
    TreeNode root;
    root.isleaf = false;
    root.father = -1;
    GTree.push_back(root);
    printf("Gtree has been added with root node %d\n",GTree.size());

    // init stack
    stack<Status> buildstack;
    Status rootstatus;
    rootstatus.tnid = 0;
    rootstatus.nset.clear();
    for (int i = 0; i < Nodes.size(); i++) {
        rootstatus.nset.insert(i);
    }
    buildstack.push(rootstatus); //push_back status of root

    // start to build
    unordered_map<int, int> presult;
    set<int> childset[PARTITION_PART];


    while (buildstack.size() > 0) {
        // pop top status
        Status current = buildstack.top();
        buildstack.pop();

        // update gtreepath
        for (set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++) {
            Nodes[*it].gtreepath.push_back(current.tnid);
        }

        // check cardinality
        if (current.nset.size() <= LEAF_CAP) {
            // build leaf node
            GTree[current.tnid].isleaf = true;
            GTree[current.tnid].leafnodes.clear();
            for (set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++) {
                GTree[current.tnid].leafnodes.push_back(*it);
            }
            continue;
        }

        // partition
//		printf("PARTITIONING...NID=%d...SIZE=%d...", current.tnid, (int)current.nset.size() );
        presult = graph_partition(current.nset);
//		printf("COMPLETE.\n");

        // construct child node set
        for (int i = 0; i < PARTITION_PART; i++) {
            childset[i].clear();
        }
        int slot;
        for (set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++) {
            slot = presult[*it];
            childset[slot].insert(*it);
        }

        // generate child tree nodes
        int childpos;
        for (int i = 0; i < PARTITION_PART; i++) {
            TreeNode tnode;
            tnode.isleaf = false;
            tnode.father = current.tnid;

            // insert to GTree first
            GTree.push_back(tnode);
            if(GTree.size()%10000==0)
                printf("Gtree has been added with node %d\n",GTree.size());
            childpos = GTree.size() - 1;
            GTree[current.tnid].children.push_back(childpos);

            // calculate border nodes
            GTree[childpos].borders.clear();
            for (set<int>::iterator it = childset[i].begin(); it != childset[i].end(); it++) {

                bool isborder = false;
                for (int j = 0; j < Nodes[*it].adjnodes.size(); j++) {
                    if (childset[i].find(Nodes[*it].adjnodes[j]) == childset[i].end()) {
                        isborder = true;
                        break;
                    }
                }
                if (isborder) {
                    GTree[childpos].borders.push_back(*it);
                    // update globally
                    Nodes[*it].isborder = true;
                }
            }

            // add to stack
            Status ongoingstatus;
            ongoingstatus.tnid = childpos;
            ongoingstatus.nset = childset[i];
            buildstack.push(ongoingstatus);

        }

    }
    printf("********************Gtree has %d nodes \n",GTree.size());


}

// dump gtree index to file
void gtree_save() {
    // FILE_GTREE
    //stringstream ss;
    //ss<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name<<".gtree";
    string filePath = getRoadInputPath(GTREE);

    FILE *fout = fopen(filePath.c_str(), "wb");
    int *buf = new int[Nodes.size()];
    for (int i = 0; i < GTree.size(); i++) {
        // borders
        int count_borders = GTree[i].borders.size();
        fwrite(&count_borders, sizeof(int), 1, fout);
        copy(GTree[i].borders.begin(), GTree[i].borders.end(), buf);
        fwrite(buf, sizeof(int), count_borders, fout);
        // children
        int count_children = GTree[i].children.size();
        fwrite(&count_children, sizeof(int), 1, fout);
        copy(GTree[i].children.begin(), GTree[i].children.end(), buf);
        fwrite(buf, sizeof(int), count_children, fout);
        // isleaf
        fwrite(&GTree[i].isleaf, sizeof(bool), 1, fout);
        /*if(GTree[i].isleaf){
            cout<<"node"<<i<<"is leaf"<<endl;
        }*/
        // leafnodes
        int count_leafnodes = GTree[i].leafnodes.size();
        fwrite(&count_leafnodes, sizeof(int), 1, fout);
        copy(GTree[i].leafnodes.begin(), GTree[i].leafnodes.end(), buf);
        fwrite(buf, sizeof(int), count_leafnodes, fout);
        // father
        fwrite(&GTree[i].father, sizeof(int), 1, fout);
    }
    fclose(fout);

    // FILE_NODES_GTREE_PATH
    //stringstream ss2;
    //ss2<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name<<".paths";
    filePath = getRoadInputPath(GPATHS);
    fout = fopen(filePath.c_str(), "wb");
    for (int i = 0; i < Nodes.size(); i++) {
        int count = Nodes[i].gtreepath.size();
        fwrite(&count, sizeof(int), 1, fout);
        copy(Nodes[i].gtreepath.begin(), Nodes[i].gtreepath.end(), buf);
        fwrite(buf, sizeof(int), count, fout);
    }
    fclose(fout);
    delete[] buf;
}

// calculate the distance matrix, algorithm shown in section 5.2 of paper
void hierarchy_shortest_path_calculation() {
    // level traversal
    vector<vector<int> > treenodelevel;

    vector<int> current;
    current.clear();
    current.push_back(0);
    treenodelevel.push_back(current);

    vector<int> mid;
    while (current.size() != 0) {
        mid = current;
        current.clear();
        for (int i = 0; i < mid.size(); i++) {
            for (int j = 0; j < GTree[mid[i]].children.size(); j++) {
                current.push_back(GTree[mid[i]].children[j]);
            }
        }
        if (current.size() == 0) break;
        treenodelevel.push_back(current);
    }

    // bottom up calculation
    // temp graph
    vector<Node> graph;
    graph = Nodes;
    vector<int> cands;
    vector<int> result;
    unordered_map<int, unordered_map<int, int> > vertex_pairs;

    // do dijkstra
    int s, t, tn, nid, cid, weight;
    vector<int> tnodes, tweight;
    set<int> nset;

    for (int i = treenodelevel.size() - 1; i >= 0; i--) {
        for (int j = 0; j < treenodelevel[i].size(); j++) {
            tn = treenodelevel[i][j];

            cands.clear();
            if (GTree[tn].isleaf) {
                // cands = leafnodes
                cands = GTree[tn].leafnodes;
                // union borders = borders;
                GTree[tn].union_borders = GTree[tn].borders;
            } else {
                nset.clear();
                for (int k = 0; k < GTree[tn].children.size(); k++) {
                    cid = GTree[tn].children[k];
                    nset.insert(GTree[cid].borders.begin(), GTree[cid].borders.end());
                }
                // union borders = cands;

                cands.clear();
                for (set<int>::iterator it = nset.begin(); it != nset.end(); it++) {
                    cands.push_back(*it);
                }
                GTree[tn].union_borders = cands;
            }

            // start to do min dis
            vertex_pairs.clear();

            // for each border, do min dis
            int cc = 0;

            for (int k = 0; k < GTree[tn].union_borders.size(); k++) {
                //printf("DIJKSTRA...LEAF=%d BORDER=%d\n", tn, GTree[tn].union_borders[k] );
                result = dijkstra_candidate(GTree[tn].union_borders[k], cands, graph);
                //printf("DIJKSTRA...END\n");

                // save to map
                for (int p = 0; p < result.size(); p++) {
                    GTree[tn].mind.push_back(result[p]);
                    vertex_pairs[GTree[tn].union_borders[k]][cands[p]] = result[p];
                }
            }

            // IMPORTANT! after all border finished, degenerate graph
            // first, remove inward edges
            for (int k = 0; k < GTree[tn].borders.size(); k++) {
                s = GTree[tn].borders[k];
                tnodes.clear();
                tweight.clear();
                for (int p = 0; p < graph[s].adjnodes.size(); p++) {
                    nid = graph[s].adjnodes[p];
                    weight = graph[s].adjweight[p];
                    // if adj node in same tree node

                    if (graph[nid].gtreepath.size() <= i || graph[nid].gtreepath[i] != tn) {
                        // only leave those useful
                        tnodes.push_back(nid);
                        tweight.push_back(weight);

                    }
                }
                // cut it
                graph[s].adjnodes = tnodes;
                graph[s].adjweight = tweight;
            }
            // second, add inter connected edges
            for (int k = 0; k < GTree[tn].borders.size(); k++) {
                for (int p = 0; p < GTree[tn].borders.size(); p++) {
                    if (k == p) continue;
                    s = GTree[tn].borders[k];
                    t = GTree[tn].borders[p];
                    graph[s].adjnodes.push_back(t);
                    graph[s].adjweight.push_back(vertex_pairs[s][t]);
                }
            }
        }
    }
}

// dump distance matrix into file
void hierarchy_shortest_path_save() {
    //stringstream ss;
    //ss<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name<<".minds";
    string filePath = getRoadInputPath(MINDS); //ss.str();


    FILE *fout = fopen(filePath.c_str(), "wb");
    int *buf;
    int count;
    for (int i = 0; i < GTree.size(); i++) {
        // union borders
        count = GTree[i].union_borders.size();
        fwrite(&count, sizeof(int), 1, fout);
        buf = new int[count];
        copy(GTree[i].union_borders.begin(), GTree[i].union_borders.end(), buf);
        fwrite(buf, sizeof(int), count, fout);
        delete[] buf;
        // mind
        count = GTree[i].mind.size();
        fwrite(&count, sizeof(int), 1, fout);
        buf = new int[count];
        copy(GTree[i].mind.begin(), GTree[i].mind.end(), buf);
        fwrite(buf, sizeof(int), count, fout);
        delete[] buf;
    }
    fclose(fout);
}


void buildGTreeIndex(){
    cout<<"begin buildGTreeIndex..."<<endl;

    TEST_START
    TIME_TICK_START
    init();
    TIME_TICK_END
    TIME_TICK_PRINT("INIT")

    // gtree_build
    TIME_TICK_START
    buildTreeStructure();
    TIME_TICK_END
    TIME_TICK_PRINT("BUILD")

    // dump gtree
    gtree_save();

    // calculate distance matrix
    TIME_TICK_START
    hierarchy_shortest_path_calculation();
    TIME_TICK_END
    TIME_TICK_PRINT("MIND")

    // dump distance matrix
    hierarchy_shortest_path_save();

    TIME_TICK_START
    generateObject();
    TIME_TICK_END
    TIME_TICK_PRINT("OBJECT")
    TEST_END
    //TEST_DURA_PRINT("************G-Tree build ")

    int leafCnt = 0;

    for(size_t i = 0;i<GTree.size();i++)
    {
        if(GTree[i].isleaf == true)
        {
            leafCnt++;
        }
    }

    cout<<" leaf cnt is "<<leafCnt<<endl;

    cout<<"buildGTreeIndex Complete!"<<endl;

}

/*----------------------------加载Gtree索引功能函数--------------------------------------*/

void gtree_load() {
    // FILE_GTREE
    string gtreePath = getRoadInputPath(GTREE);
    FILE *fin = fopen(gtreePath.c_str(), "rb");
    int *buf = new int[Nodes.size()];
    int count_borders, count_children, count_leafnodes;
    bool isleaf;
    int father;

    // clear gtree
    GTree.clear();
    int n_id = 0;
    while (fread(&count_borders, sizeof(int), 1, fin)) {
        TreeNode tn;
        // borders
        tn.borders.clear();
        fread(buf, sizeof(int), count_borders, fin);
        for (int i = 0; i < count_borders; i++) {
            tn.borders.push_back(buf[i]);
        }
        // children
        fread(&count_children, sizeof(int), 1, fin);
        fread(buf, sizeof(int), count_children, fin);
        for (int i = 0; i < count_children; i++) {
            tn.children.push_back(buf[i]);
        }
        // isleaf
        fread(&isleaf, sizeof(bool), 1, fin);
        tn.isleaf = isleaf;
        /*if(tn.isleaf){
            cout<<"find node"<<n_id<<" is a leaf node!"<<endl;
        }*/
        // leafnodes
        fread(&count_leafnodes, sizeof(int), 1, fin);  //count_leafnodes是空的
        int num = count_leafnodes;
        fread(buf, sizeof(int), count_leafnodes, fin);
        for (int i = 0; i < count_leafnodes; i++) {  //记录vertex在leafnodes数组中的idx
            int v = buf[i]; int idx = i;
            tn.leafnodes.push_back(v); id2idx_vertexInLeaf[v]= i;

        }
        // father
        fread(&father, sizeof(int), 1, fin);
        tn.father = father;

        GTree.push_back(tn);
        n_id++;
    }
    fclose(fin);

    // FILE_NODES_GTREE_PATH
    int count;  //level
    fin = fopen(getRoadInputPath(GPATHS).c_str(), "rb");
    int pos = 0;
    while (fread(&count, sizeof(int), 1, fin)) {
        fread(buf, sizeof(int), count, fin);
        // clear gtreepath
        Nodes[pos].gtreepath.clear();
        //cout<<"构建v"<<pos<<"的 node gtree_path"<<endl;
        for (int i = 0; i < count; i++) {
            int node_id = buf[i];
            Nodes[pos].gtreepath.push_back(node_id);
        }
       /* if(pos == 23550){
            cout<<"now v23550 is actually in node"<<Nodes[pos].gtreepath.back()<<endl;
        }*/
        // pos increase
        pos++;
    }
    fclose(fin);
    //getchar();


    delete[] buf;
}


// load distance matrix from file
void hierarchy_shortest_path_load() {
    FILE *fin = fopen(getRoadInputPath(MINDS).c_str(), "rb");
    int *buf;
    int count, pos = 0;
    while (fread(&count, sizeof(int), 1, fin)) {
        // union borders
        buf = new int[count];
        fread(buf, sizeof(int), count, fin);
        GTree[pos].union_borders.clear();
        for (int i = 0; i < count; i++) {
            GTree[pos].union_borders.push_back(buf[i]);
        }
        delete[] buf;
        // mind
        fread(&count, sizeof(int), 1, fin);
        buf = new int[count];
        fread(buf, sizeof(int), count, fin);
        GTree[pos].mind.clear();
        for (int i = 0; i < count; i++) {
            GTree[pos].mind.push_back(buf[i]);
        }
        pos++;
        delete[] buf;
    }
    fclose(fin);
}

// output size of Node's border
void output_borderSize() {
    ofstream out("data/sizeOfBorder");
    for (size_t i = 0; i < GTree.size(); i++) {
        out << " GTree Node " << i << "'s Border size is " << GTree[i].borders.size() << endl;
    }
    out.close();
}

// before query, we have to set OCCURENCE LIST etc.
// this is done only ONCE for a given set of objects.
void pre_query() {
    // first clear all
    for (int i = 0; i < GTree.size(); i++) {
        GTree[i].nonleafinvlist.clear();
        GTree[i].leafinvlist.clear();
    }

    // read object list
    vector<int> o;
    o.clear();

    //char file_object[100];
    //sprintf(file_object, "%s/%s.", getRoadInputPath(HOME),dataset_name);
    FILE *fin = fopen(getRoadInputPath(OBJECTS).c_str(), "r");
    int oid, id;
    while (fscanf(fin, "%d %d", &oid, &id) == 2) {
        o.push_back(oid);
    }
    fclose(fin);

    // set OCCURENCE LIST
    for (int i = 0; i < o.size(); i++) {
        int current = Nodes[o[i]].gtreepath.back();
        // add leaf inv list
        int pos = lower_bound(GTree[current].leafnodes.begin(), GTree[current].leafnodes.end(), o[i]) -
                  GTree[current].leafnodes.begin();
        GTree[current].leafinvlist.push_back(pos);
        // recursive
        int child;
        while (current != -1) {
            child = current;
            current = GTree[current].father;
            if (current == -1) break;
            if (find(GTree[current].nonleafinvlist.begin(), GTree[current].nonleafinvlist.end(), child) ==
                GTree[current].nonleafinvlist.end()) {
                GTree[current].nonleafinvlist.push_back(child);
            }
        }
    }

    // up_pos & current_pos(used for quickly locating parent & child nodes)
    unordered_map<int, int> pos_map;
    for (int i = 1; i < GTree.size(); i++) {
        GTree[i].current_pos.clear();
        GTree[i].up_pos.clear();

        // current_pos
        pos_map.clear();
        for (int j = 0; j < GTree[i].union_borders.size(); j++) {
            pos_map[GTree[i].union_borders[j]] = j;
        }
        for (int j = 0; j < GTree[i].borders.size(); j++) {
            GTree[i].current_pos.push_back(pos_map[GTree[i].borders[j]]);
        }
        // up_pos
        pos_map.clear();
        for (int j = 0; j < GTree[GTree[i].father].union_borders.size(); j++) {
            pos_map[GTree[GTree[i].father].union_borders[j]] = j;
        }
        for (int j = 0; j < GTree[i].borders.size(); j++) {
            GTree[i].up_pos.push_back(pos_map[GTree[i].borders[j]]);
        }
    }
}

/*----------------------------加载Gtree索引功能函数--------------------------------------*/


//compute max dist and min dist of two Rnet


// init search node
typedef struct {
    int id;
    bool isvertex;
    int lca_pos;
    int dis;
    double score;
    //bool isobject;
} Status_query;

typedef struct {
    int id;
    bool isobject;
    int lca_pos;
    double score;
} Topk_entry;


struct Status_query_comp {
    bool operator()(const Status_query &l, const Status_query &r) {
        return l.dis > r.dis;
        //return 1000000000*l.score < 1000000000*r.score;
    }
};

struct Status_query_max_comp {
    bool operator()(const Status_query &l, const Status_query &r) {
        //return l.dis > r.dis;
        return 1.0/r.score < 1.0/l.score;
        //return 1.0/l.score > 1.0/r.score;   //大的排后头
    }
};

typedef struct {
    int id;
    int dis;
} ResultSet;

typedef struct {
    int id;
    double score;
} TopResult;

int initGtree() {
    LOAD_START
    init();
    // load gtree index
    gtree_load();

    // load distance matrix
    hierarchy_shortest_path_load();

    pre_query();// pre_que

    LOAD_END
    printf("load Gtree Node Success! NodeNum=%d ", GTree.size()); LOAD_PRINT()
}





// ----- CORE PART -----
// knn search
// input: locid = query location, node id
//        K = top-K
// output: a vector of ResultSet, each is a tuple (node id, shortest path), ranked by shortest path distance from query location
vector<ResultSet> knn_query(int locid, int K) {
    // init priority queue & result set
    vector<Status_query> pq;
    pq.clear();
    vector<ResultSet> rstset;
    rstset.clear();

    // init upstream
    unordered_map<int, vector<int> > itm; // intermediate answer, tree node -> array
    itm.clear();
    int tn, cid, posa, posb, _min, dis;
    for (int i = Nodes[locid].gtreepath.size() - 1; i > 0; i--) {   //自底向上
        tn = Nodes[locid].gtreepath[i];
        itm[tn].clear();

        if (GTree[tn].isleaf) {
            posa = lower_bound(GTree[tn].leafnodes.begin(), GTree[tn].leafnodes.end(), locid) -
                   GTree[tn].leafnodes.begin();

            for (int j = 0; j < GTree[tn].borders.size(); j++) {
                itm[tn].push_back(GTree[tn].mind[j * GTree[tn].leafnodes.size() + posa]);
            }
        } else {
            cid = Nodes[locid].gtreepath[i + 1];
            for (int j = 0; j < GTree[tn].borders.size(); j++) {
                _min = -1;
                posa = GTree[tn].current_pos[j];
                for (int k = 0; k < GTree[cid].borders.size(); k++) {
                    posb = GTree[cid].up_pos[k];
                    dis = itm[cid][k] + GTree[tn].mind[posa * GTree[tn].union_borders.size() + posb];
                    // get min
                    if (_min == -1) {
                        _min = dis;
                    } else {
                        if (dis < _min) {
                            _min = dis;
                        }
                    }
                }
                // update
                itm[tn].push_back(_min);      //自底向上记录 loc 到其所在 叶节点，以及叶节点的祖先节点的最小距离
            }
        }

    }

    // do search
    Status_query rootstatus = {0, false, 0, 0};
    pq.push_back(rootstatus);
    make_heap(pq.begin(), pq.end(), Status_query_comp());

    vector<int> cands, result;
    int child, son, allmin, vertex;

    while (pq.size() > 0 && rstset.size() < K) {
        Status_query top = pq[0];
        pop_heap(pq.begin(), pq.end(), Status_query_comp());
        pq.pop_back();

        if (top.isvertex) {  //当前最优为某个路网顶点
            ResultSet rs = {top.id, top.dis};
            rstset.push_back(rs);
        } else {           //当前最优为某个叶节点
            if (GTree[top.id].isleaf) {
                // inner of leaf node, do dijkstra
                if (top.id == Nodes[locid].gtreepath[top.lca_pos]) {

                    cands.clear();
                    for (int i = 0; i < GTree[top.id].leafinvlist.size(); i++) {
                        cands.push_back(GTree[top.id].leafnodes[GTree[top.id].leafinvlist[i]]);
                    }
                    result = dijkstra_candidate(locid, cands, Nodes);
                    for (int i = 0; i < cands.size(); i++) {
                        Status_query status = {cands[i], true, top.lca_pos, result[i]};
                        pq.push_back(status);
                        push_heap(pq.begin(), pq.end(), Status_query_comp());
                    }

                }

                    // else do
                else {          //是叶节点，但与loc不在同一子图
                    for (int i = 0; i < GTree[top.id].leafinvlist.size(); i++) {    //获取当前最优叶节点的顶点倒排表
                        posa = GTree[top.id].leafinvlist[i];
                        vertex = GTree[top.id].leafnodes[posa];
                        allmin = -1;

                        for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                            //计算lco 经过 叶节点中第 k个 border 到达 top.id叶节点内第i个顶点的距离
                            dis = itm[top.id][k] + GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa];
                            if (allmin == -1) {
                                allmin = dis;
                            } else {
                                if (dis < allmin) {
                                    allmin = dis;   // 提取loc 经过所有叶节点top.id中的boder, 而到达第i个顶点的最短距离
                                }
                            }

                        }

                        Status_query status = {vertex, true, top.lca_pos, allmin};
                        pq.push_back(status);
                        push_heap(pq.begin(), pq.end(), Status_query_comp());

                    }
                }
            } else {    //非叶节点
                for (int i = 0; i < GTree[top.id].nonleafinvlist.size(); i++) {
                    child = GTree[top.id].nonleafinvlist[i];
                    son = Nodes[locid].gtreepath[top.lca_pos + 1];
                    // on gtreepath
                    if (child == son) {
                        Status_query status = {child, false, top.lca_pos + 1, 0};
                        pq.push_back(status);
                        push_heap(pq.begin(), pq.end(), Status_query_comp());
                    }
                    // brothers
                    else if (GTree[child].father == GTree[son].father) {
                        itm[child].clear();
                        allmin = -1;

                        for (int j = 0; j < GTree[child].borders.size(); j++) {   //把到当前最优节点的各个child的距离算出，把child加入
                            _min = -1;
                            posa = GTree[child].up_pos[j];
                            for (int k = 0; k < GTree[son].borders.size(); k++) {
                                posb = GTree[son].up_pos[k];
                                dis = itm[son][k] +
                                      GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                                if (_min == -1) {
                                    _min = dis;
                                } else {
                                    if (dis < _min) {
                                        _min = dis;
                                    }
                                }
                            }
                            itm[child].push_back(_min);  //保存到child的距离
                            // update all min
                            if (allmin == -1) {
                                allmin = _min;
                            } else if (_min < allmin) {
                                allmin = _min;
                            }
                        }
                        Status_query status = {child, false, top.lca_pos, allmin};
                        pq.push_back(status);
                        push_heap(pq.begin(), pq.end(), Status_query_comp());
                    }
                        // downstream
                    else {
                        itm[child].clear();
                        allmin = -1;

                        for (int j = 0; j < GTree[child].borders.size(); j++) {
                            _min = -1;
                            posa = GTree[child].up_pos[j];
                            for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                posb = GTree[top.id].current_pos[k];
                                dis = itm[top.id][k] +
                                      GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                                if (_min == -1) {
                                    _min = dis;
                                } else {
                                    if (dis < _min) {
                                        _min = dis;
                                    }
                                }
                            }
                            itm[child].push_back(_min);
                            // update all min
                            if (allmin == -1) {
                                allmin = _min;
                            } else if (_min < allmin) {
                                allmin = _min;
                            }
                        }
                        Status_query status = {child, false, top.lca_pos, allmin};
                        pq.push_back(status);
                        push_heap(pq.begin(), pq.end(), Status_query_comp());
                    }
                }
            }

        }
    }

    /* ----- return id list -----
    vector<int> rst;
    for ( int i = 0; i < rstset.size(); i++ ){
        rst.push_back( rstset[i].id) ;
    }

    return rst;
    */

    return rstset;
}







int testGtree() {

    init();
    // load gtree index
    gtree_load();

    // load distance matrix
    hierarchy_shortest_path_load();

    pre_query();// pre_que

    vector<ResultSet> result;
    result = knn_query(2, VertexNum);
    float maxDist = (float) result.back().dis / WEIGHT_INFLATE_FACTOR;
    cout << " max dist is " << maxDist << endl;
    //for (size_t i = 0; i < result.size(); i++) {
    //    cout << result[i].id << " " << result[i].dis << endl;
    //}
    return 0;
}


//SPSP query
float SPSP_old(int Ni,int Nj)
{
    if(SPDCache[Ni][Nj]!=0) return SPDCache[Ni][Nj];
    if(Nodes[Ni].gtreepath.back() == Nodes[Nj].gtreepath.back())
    {
        int NumofNode = VertexNum;
        adjacency_list_t adjacency_list(NumofNode);
        int GtreeNode = Nodes[Ni].gtreepath.back();
        for(int i=0;i<GTree[GtreeNode].leafnodes.size();i++)
        {
            int leafNode_Ni = GTree[GtreeNode].leafnodes[i];
            for(int j=0;j<Nodes[leafNode_Ni].adjnodes.size();j++)
            {
                int leafNode_Nj = Nodes[leafNode_Ni].adjnodes[j];
                float dist = Nodes[leafNode_Ni].adjweight[j];
                adjacency_list[leafNode_Ni].push_back(neighbor(leafNode_Nj, dist));
                adjacency_list[leafNode_Nj].push_back(neighbor(leafNode_Ni, dist));
            }
        }
        std::vector<weight_t> min_distance;

        DijkstraComputePaths(Ni, adjacency_list, min_distance);
        SPDCache[Ni][Nj] = min_distance[Nj]; //by jins
        return min_distance[Nj];
    }
    else
    {
        map<int,vector<int> > level_node;
        vector<int> SizeOfNode;
        SizeOfNode.push_back(1);
        int _start;
        for(int i = 0;i < Nodes[Ni].gtreepath.size();i++)
        {
            if(Nodes[Ni].gtreepath[i] == Nodes[Nj].gtreepath[i] && Nodes[Ni].gtreepath[i+1] != Nodes[Nj].gtreepath[i+1] )
            {
                //cout<<"the Nodes is "<<Nodes[Ni].gtreepath[i]<<endl;
                _start = i;
                break;
            }
        }
        int level_count = 0;
        //cout<<"_start is "<<_start<<endl;

        for(int i = Nodes[Ni].gtreepath.size() - 1;i > _start;i--)
        {
            int sizeofN = SizeOfNode[level_count] + GTree[Nodes[Ni].gtreepath[i]].borders.size();
            SizeOfNode.push_back(sizeofN);
            //cout<<"border is ";
            for(int j =0;j<GTree[Nodes[Ni].gtreepath[i]].borders.size();j++)
            {
                level_node[level_count].push_back(GTree[Nodes[Ni].gtreepath[i]].borders[j]);
                //cout<<GTree[Nodes[Ni].gtreepath[i]].borders[j]<<" ";
            }
            level_count++;
            //cout<<endl;
        }

        for(int i = _start + 1;i < Nodes[Nj].gtreepath.size();i++)
        {
            int sizeofN = SizeOfNode[level_count] + GTree[Nodes[Nj].gtreepath[i]].borders.size();
            SizeOfNode.push_back(sizeofN);
            //cout<<"border is ";
            for(int j =0;j<GTree[Nodes[Nj].gtreepath[i]].borders.size();j++)
            {
                level_node[level_count].push_back(GTree[Nodes[Nj].gtreepath[i]].borders[j]);
                //cout<<GTree[Nodes[Nj].gtreepath[i]].borders[j]<<" ";
            }
            level_count++;
            //cout<<endl;
        }
        int sizeofN = SizeOfNode[level_count] + 1;
        SizeOfNode.push_back(sizeofN);
        //cout<<"Node_num is "<<SizeOfNode.back()<<endl;
        //cout<<"level_count is "<<level_count<<endl;

        int Node_Num = SizeOfNode.back();
        vector<int> dist;
        vector<int> path;
        int **a;
        a = new int*[Node_Num];
        for(int i = 0; i != Node_Num; ++i)
        {
            a[i] = new int[Node_Num];
        }

        for(int i=0;i<Node_Num;i++)
        {
            for(int j=0;j<Node_Num;j++)
            {
                a[i][j] = 99999999;
            }
        }
        vector<int> Node_list;
        int _flag = 0;

        for(int i = Nodes[Ni].gtreepath.size() - 1;i > _start;i--)
        {
            int entry_id = Nodes[Ni].gtreepath[i];
            if(GTree[entry_id].isleaf == true)
            {
                for(int j=0;j<GTree[entry_id].borders.size();j++)
                {
                    int pos_j = lower_bound( GTree[entry_id].leafnodes.begin(), GTree[entry_id].leafnodes.end(), Ni ) - GTree[entry_id].leafnodes.begin();
                    //cout<<" Ni ="<<GTree[entry_id].leafnodes[pos_j]<<endl;
                    int distOfNtoN = GTree[entry_id].mind[j * GTree[entry_id].leafnodes.size() + pos_j];
                    //cout<<"Ni to "<<GTree[entry_id].borders[j]<<" dist is "<<distOfNtoN<<endl;
                    a[0][SizeOfNode[_flag] + j] = distOfNtoN;
                }
            }
            else
            {
                int highentry_id = Nodes[Ni].gtreepath[i+1];
                for(int j=0;j<GTree[highentry_id].borders.size();j++)
                {
                    int pos_i = lower_bound( GTree[entry_id].union_borders.begin(), GTree[entry_id].union_borders.end(), GTree[highentry_id].borders[j] ) - GTree[entry_id].union_borders.begin();
                    for(int k=0;k<GTree[entry_id].borders.size();k++)
                    {
                        int pos_j = lower_bound( GTree[entry_id].union_borders.begin(), GTree[entry_id].union_borders.end(), GTree[entry_id].borders[k] ) - GTree[entry_id].union_borders.begin();
                        int distOfNtoN = GTree[entry_id].mind[pos_i*GTree[entry_id].union_borders.size() + pos_j];
                        a[SizeOfNode[_flag] + j][SizeOfNode[_flag + 1] + k] = distOfNtoN;
                    }
                }
                _flag++;
            }
            //cout<<entry_id<<" ";
        }
        if(true)
        {
            int highentry_id = Nodes[Ni].gtreepath[_start + 1];
            int entry_id = Nodes[Nj].gtreepath[_start + 1];
            int start_id= Nodes[Ni].gtreepath[_start];

            for(int j=0;j<GTree[highentry_id].borders.size();j++)
            {
                int pos_i = lower_bound( GTree[start_id].union_borders.begin(), GTree[start_id].union_borders.end(), GTree[highentry_id].borders[j] ) - GTree[start_id].union_borders.begin();
                for(int k=0;k<GTree[entry_id].borders.size();k++)
                {
                    int pos_j = lower_bound( GTree[start_id].union_borders.begin(), GTree[start_id].union_borders.end(), GTree[entry_id].borders[k] ) - GTree[start_id].union_borders.begin();
                    int distOfNtoN = GTree[start_id].mind[pos_i*GTree[start_id].union_borders.size() + pos_j];
                    a[SizeOfNode[_flag] + j][SizeOfNode[_flag + 1] + k] = distOfNtoN;
                }
            }
            _flag++;
        }



        for(int i = _start + 1 ;i < Nodes[Nj].gtreepath.size();i++)
        {
            int entry_id = Nodes[Nj].gtreepath[i];
            if(GTree[entry_id].isleaf == true)
            {
                for(int j=0;j<GTree[entry_id].borders.size();j++)
                {
                    int pos_j = lower_bound( GTree[entry_id].leafnodes.begin(), GTree[entry_id].leafnodes.end(), Nj ) - GTree[entry_id].leafnodes.begin();
                    //cout<<" Ni ="<<GTree[entry_id].leafnodes[pos_j]<<endl;
                    int distOfNtoN = GTree[entry_id].mind[j * GTree[entry_id].leafnodes.size() + pos_j];
                    //cout<<"Ni to "<<GTree[entry_id].borders[j]<<" dist is "<<distOfNtoN<<endl;
                    a[SizeOfNode[_flag] + j][Node_Num - 1] = distOfNtoN;
                }
            }
            else
            {
                int highentry_id = Nodes[Nj].gtreepath[i+1];
                for(int j=0;j<GTree[entry_id].borders.size();j++)
                {
                    int pos_i = lower_bound( GTree[entry_id].union_borders.begin(), GTree[entry_id].union_borders.end(), GTree[entry_id].borders[j] ) - GTree[entry_id].union_borders.begin();
                    for(int k=0;k<GTree[highentry_id].borders.size();k++)
                    {
                        int pos_j = lower_bound( GTree[entry_id].union_borders.begin(), GTree[entry_id].union_borders.end(), GTree[highentry_id].borders[k] ) - GTree[entry_id].union_borders.begin();
                        int distOfNtoN = GTree[entry_id].mind[pos_i*GTree[entry_id].union_borders.size() + pos_j];
                        a[SizeOfNode[_flag] + j][SizeOfNode[_flag + 1] + k] = distOfNtoN;
                    }
                }
                _flag++;
            }
            //cout<<entry_id<<" ";
        }

        //fpath(a,VertexNum,dist);
        int i,j,k;
        dist.push_back(0);
        int max_Num = 99999999;
        for(i=1;i<Node_Num;i++){
            k = max_Num;
            for(j=0;j<i;j++){
                if(a[j][i] != max_Num)
                    if((dist[j]+a[j][i])<k)
                        k=dist[j]+a[j][i];
            }
            dist.push_back(k);
        }


        for(int i = 0; i != Node_Num; i++)
        {
            delete[] a[i];
        }
        delete[] a;

        SPDCache[Ni][Nj] = dist[dist.size()-1]; //by jins
        return dist[dist.size()-1];

    }
}

float SPSP(int Ni,int Nj)
{
    if(SPDCache[Ni][Nj]!=0) return SPDCache[Ni][Nj];
    if(Nodes[Ni].gtreepath.back() == Nodes[Nj].gtreepath.back())
    {
        int NumofNode = VertexNum;
        adjacency_list_t adjacency_list(NumofNode);
        int GtreeNode = Nodes[Ni].gtreepath.back();
        for(int i=0;i<GTree[GtreeNode].leafnodes.size();i++)
        {
            int leafNode_Ni = GTree[GtreeNode].leafnodes[i];
            for(int j=0;j<Nodes[leafNode_Ni].adjnodes.size();j++)
            {
                int leafNode_Nj = Nodes[leafNode_Ni].adjnodes[j];
                float dist = Nodes[leafNode_Ni].adjweight[j];
                adjacency_list[leafNode_Ni].push_back(neighbor(leafNode_Nj, dist));
                adjacency_list[leafNode_Nj].push_back(neighbor(leafNode_Ni, dist));
            }
        }
        std::vector<weight_t> min_distance;

        DijkstraComputePaths(Ni, adjacency_list, min_distance);
        SPDCache[Ni][Nj] = min_distance[Nj]; //by jins
        return min_distance[Nj];
    }
    else
    {
        map<int,vector<int> > level_node;
        vector<int> SizeOfNode;
        SizeOfNode.push_back(1);
        int _start;
        for(int i = 0;i < Nodes[Ni].gtreepath.size();i++)
        {
            if(Nodes[Ni].gtreepath[i] == Nodes[Nj].gtreepath[i] && Nodes[Ni].gtreepath[i+1] != Nodes[Nj].gtreepath[i+1] )
            {
                //cout<<"the Nodes is "<<Nodes[Ni].gtreepath[i]<<endl;
                _start = i;
                break;
            }
        }
        int level_count = 0;
        //cout<<"_start is "<<_start<<endl;

        for(int i = Nodes[Ni].gtreepath.size() - 1;i > _start;i--)
        {
            int sizeofN = SizeOfNode[level_count] + GTree[Nodes[Ni].gtreepath[i]].borders.size();
            SizeOfNode.push_back(sizeofN);
            //cout<<"border is ";
            for(int j =0;j<GTree[Nodes[Ni].gtreepath[i]].borders.size();j++)
            {
                level_node[level_count].push_back(GTree[Nodes[Ni].gtreepath[i]].borders[j]);
                //cout<<GTree[Nodes[Ni].gtreepath[i]].borders[j]<<" ";
            }
            level_count++;
            //cout<<endl;
        }

        for(int i = _start + 1;i < Nodes[Nj].gtreepath.size();i++)
        {
            int sizeofN = SizeOfNode[level_count] + GTree[Nodes[Nj].gtreepath[i]].borders.size();
            SizeOfNode.push_back(sizeofN);
            //cout<<"border is ";
            for(int j =0;j<GTree[Nodes[Nj].gtreepath[i]].borders.size();j++)
            {
                level_node[level_count].push_back(GTree[Nodes[Nj].gtreepath[i]].borders[j]);
                //cout<<GTree[Nodes[Nj].gtreepath[i]].borders[j]<<" ";
            }
            level_count++;
            //cout<<endl;
        }
        int sizeofN = SizeOfNode[level_count] + 1;
        SizeOfNode.push_back(sizeofN);
        //cout<<"Node_num is "<<SizeOfNode.back()<<endl;
        //cout<<"level_count is "<<level_count<<endl;

        int Node_Num = SizeOfNode.back();
        vector<int> dist;
        vector<int> path;
        int **a;
        a = new int*[Node_Num];
        for(int i = 0; i != Node_Num; ++i)
        {
            a[i] = new int[Node_Num];
        }

        for(int i=0;i<Node_Num;i++)
        {
            for(int j=0;j<Node_Num;j++)
            {
                a[i][j] = 99999999;
            }
        }
        vector<int> Node_list;
        int _flag = 0;

        for(int i = Nodes[Ni].gtreepath.size() - 1;i > _start;i--)
        {
            int entry_id = Nodes[Ni].gtreepath[i];
            if(GTree[entry_id].isleaf == true)
            {
                for(int j=0;j<GTree[entry_id].borders.size();j++)
                {
                    int pos_j = lower_bound( GTree[entry_id].leafnodes.begin(), GTree[entry_id].leafnodes.end(), Ni ) - GTree[entry_id].leafnodes.begin();
                    //cout<<" Ni ="<<GTree[entry_id].leafnodes[pos_j]<<endl;
                    int distOfNtoN = GTree[entry_id].mind[j * GTree[entry_id].leafnodes.size() + pos_j];
                    //cout<<"Ni to "<<GTree[entry_id].borders[j]<<" dist is "<<distOfNtoN<<endl;
                    a[0][SizeOfNode[_flag] + j] = distOfNtoN;
                }
            }
            else
            {
                int highentry_id = Nodes[Ni].gtreepath[i+1];
                for(int j=0;j<GTree[highentry_id].borders.size();j++)
                {
                    int pos_i = lower_bound( GTree[entry_id].union_borders.begin(), GTree[entry_id].union_borders.end(), GTree[highentry_id].borders[j] ) - GTree[entry_id].union_borders.begin();
                    for(int k=0;k<GTree[entry_id].borders.size();k++)
                    {
                        int pos_j = lower_bound( GTree[entry_id].union_borders.begin(), GTree[entry_id].union_borders.end(), GTree[entry_id].borders[k] ) - GTree[entry_id].union_borders.begin();
                        int distOfNtoN = GTree[entry_id].mind[pos_i*GTree[entry_id].union_borders.size() + pos_j];
                        a[SizeOfNode[_flag] + j][SizeOfNode[_flag + 1] + k] = distOfNtoN;
                    }
                }
                _flag++;
            }
            //cout<<entry_id<<" ";
        }
        if(true)
        {
            int highentry_id = Nodes[Ni].gtreepath[_start + 1];
            int entry_id = Nodes[Nj].gtreepath[_start + 1];
            int start_id= Nodes[Ni].gtreepath[_start];

            for(int j=0;j<GTree[highentry_id].borders.size();j++)
            {
                int pos_i = lower_bound( GTree[start_id].union_borders.begin(), GTree[start_id].union_borders.end(), GTree[highentry_id].borders[j] ) - GTree[start_id].union_borders.begin();
                for(int k=0;k<GTree[entry_id].borders.size();k++)
                {
                    int pos_j = lower_bound( GTree[start_id].union_borders.begin(), GTree[start_id].union_borders.end(), GTree[entry_id].borders[k] ) - GTree[start_id].union_borders.begin();
                    int distOfNtoN = GTree[start_id].mind[pos_i*GTree[start_id].union_borders.size() + pos_j];
                    a[SizeOfNode[_flag] + j][SizeOfNode[_flag + 1] + k] = distOfNtoN;
                }
            }
            _flag++;
        }



        for(int i = _start + 1 ;i < Nodes[Nj].gtreepath.size();i++)
        {
            int entry_id = Nodes[Nj].gtreepath[i];
            if(GTree[entry_id].isleaf == true)
            {
                for(int j=0;j<GTree[entry_id].borders.size();j++)
                {
                    int pos_j = lower_bound( GTree[entry_id].leafnodes.begin(), GTree[entry_id].leafnodes.end(), Nj ) - GTree[entry_id].leafnodes.begin();
                    //cout<<" Ni ="<<GTree[entry_id].leafnodes[pos_j]<<endl;
                    int distOfNtoN = GTree[entry_id].mind[j * GTree[entry_id].leafnodes.size() + pos_j];
                    //cout<<"Ni to "<<GTree[entry_id].borders[j]<<" dist is "<<distOfNtoN<<endl;
                    a[SizeOfNode[_flag] + j][Node_Num - 1] = distOfNtoN;
                }
            }
            else
            {
                int highentry_id = Nodes[Nj].gtreepath[i+1];
                for(int j=0;j<GTree[entry_id].borders.size();j++)
                {
                    int pos_i = lower_bound( GTree[entry_id].union_borders.begin(), GTree[entry_id].union_borders.end(), GTree[entry_id].borders[j] ) - GTree[entry_id].union_borders.begin();
                    for(int k=0;k<GTree[highentry_id].borders.size();k++)
                    {
                        int pos_j = lower_bound( GTree[entry_id].union_borders.begin(), GTree[entry_id].union_borders.end(), GTree[highentry_id].borders[k] ) - GTree[entry_id].union_borders.begin();
                        int distOfNtoN = GTree[entry_id].mind[pos_i*GTree[entry_id].union_borders.size() + pos_j];
                        a[SizeOfNode[_flag] + j][SizeOfNode[_flag + 1] + k] = distOfNtoN;
                    }
                }
                _flag++;
            }
            //cout<<entry_id<<" ";
        }

        //fpath(a,VertexNum,dist);
        int i,j,k;
        dist.push_back(0);
        int max_Num = 99999999;
        for(i=1;i<Node_Num;i++){
            k = max_Num;
            for(j=0;j<i;j++){
                if(a[j][i] != max_Num)
                    if((dist[j]+a[j][i])<k)
                        k=dist[j]+a[j][i];
            }
            dist.push_back(k);
        }


        for(int i = 0; i != Node_Num; i++)
        {
            delete[] a[i];
        }
        delete[] a;

        SPDCache[Ni][Nj] = dist[dist.size()-1]; //by jins
        return dist[dist.size()-1];

    }
}


double getDistance(int loc, User usr){
    int Ni = min(usr.Ni,usr.Nj);
    int Nj = max(usr.Ni,usr.Nj);
    double distance;
    double d1 = SPSP(loc,Ni)/WEIGHT_INFLATE_FACTOR +usr.dis;
    double d2 = SPSP(loc,Nj)/WEIGHT_INFLATE_FACTOR + (usr.dist - usr.dis);

    if(loc == Ni)
        return usr.dis;
    if(loc == Nj)
        return usr.dist-usr.dis;

    return min(d1,d2);


}

double getDistance(int loc, POI poi){
    int Ni = min(poi.Ni,poi.Nj);
    int Nj = max(poi.Ni,poi.Nj);
    double distance;
    double d1 = SPSP(loc,Ni)/WEIGHT_INFLATE_FACTOR +poi.dis;
    double d2 = SPSP(loc,Nj)/WEIGHT_INFLATE_FACTOR + (poi.dist - poi.dis);

    if(loc == Ni)
        return poi.dis;
    if(loc == Nj)
        return poi.dist-poi.dis;

    return min(d1,d2);


}



float usrToPOIDistance(User usr,POI poi){
    float distance =0;
    if(poi.Ni == usr.Ni && poi.Nj == usr.Nj){
        if (usr.dis > poi.dis)
            distance = usr.dis - poi.dis;
        else
            distance = poi.dis - usr.dis;
        //cout<<"同边"<<endl;
        return distance;
    }
    if(poi.Ni == usr.Ni && poi.Nj != usr.Nj){
        distance = usr.dis+poi.dis;
        //cout<<"邻边"<<endl;
        return distance;
    }

    double d1 = getDistance(poi.Ni,usr) + poi.dis;
    double d2 = getDistance(poi.Nj,usr)+poi.dist-poi.dis;
    //cout<<"spsp"<<endl;
    //if(usr.id==625&&poi.id==57) cout<<d1<<","<<d2<<endl;
    return min(d1,d2);
}



float usrToPOIDistance_u_Ni(User usr,POI poi){
    float distance =0;
    if(poi.Ni == usr.Ni && poi.Nj == usr.Nj){
        if (usr.dis > poi.dis)
            distance = usr.dis - poi.dis;
        else
            distance = poi.dis - usr.dis;
        //cout<<"同边"<<endl;
        return distance;
    }
    if(poi.Ni == usr.Ni && poi.Nj != usr.Nj){
        distance = usr.dis+poi.dis;
        //cout<<"邻边"<<endl;
        return distance;
    }

    double d1 = getDistance(usr.Ni,poi) + usr.dis;
    //double d2 = getDistance(poi.Nj,usr)+poi.dist-poi.dis;
    //cout<<"spsp"<<endl;
    //if(usr.id==625&&poi.id==57) cout<<d1<<","<<d2<<endl;
    return d1; //min(d1,d2);
}



double getDistance(User usr, POI poi){
    int u_Ni = min(usr.Ni,usr.Nj);
    int u_Nj = max(usr.Ni,usr.Nj);
    int p_Ni = min(poi.Ni,poi.Nj);
    int p_Nj = max(poi.Ni,poi.Nj);
    if(u_Ni==p_Ni && u_Nj==p_Nj)
        return fabs(poi.dis-usr.dis);


    double distance;
    double d1 = SPSP(poi.Ni,usr.Ni)/1000+poi.dis+usr.dis;
    double d2 = SPSP(poi.Ni,usr.Nj)/1000+poi.dis+usr.dist-usr.dis;
    double di = min(d1,d2);

    double d3 = SPSP(poi.Nj,usr.Ni)/1000+poi.dist-poi.dis+usr.dis;
    double d4 = SPSP(poi.Nj,usr.Nj)/1000+poi.dist-poi.dis+usr.dist-usr.dis;
    double dj = min(d3,d4);
    return min(di,dj);

}

double getDistance_upper_Gtree(User user, POI poi){

    double distance = SPSP(user.Ni,poi.Nj)/1.0+user.dis+poi.dis;
    return  distance;

}




struct MinMaxD{
    int min;
    int max;

};

struct BorderPair{
    int b1;
    int b2;
};

BorderPair **BorderPairMatrix;


MinMaxD getMinMaxDistance(int loc, int Node_id){
    MinMaxD mm;
    /*
    //验证 loc 是否在 Node_id所在的子图中
    if(
    for(int i=0; i<Nodes[loc].gtreepath.size();i++)
        if(Nodes[loc].gtreepath[i] == Node_id) {
            mm.min = knn_query(loc, 1)[0].dis / WEIGHT_INFLATE_FACTOR;
            mm.max = ;
            return  mm

        }
        */
    //loc 不在 Node_id所在的子图中
    vector<int> bs = GTree[Node_id].borders;
    float maxDis = 0;
    float minDis = 1000000;
    for(int i = 0;i<bs.size();i++){
        int border_id = bs[i];
        double distance;
        //if(SPDMap[loc][border_id]!=NULL)  distance = SPDMap[loc][border_id];
        distance = SPSP(loc, border_id)/ WEIGHT_INFLATE_FACTOR;
        if(distance < minDis) minDis = distance;
        if(distance > maxDis) maxDis = distance;
    }
    mm.min = minDis;
    mm.max = maxDis;
    if(minDis==0){  //loc即为Node的border
        if(Nodes[loc].gtreepath.back() == Node_id )
            mm.min = GTree[Node_id].minDis;
    }
    return  mm;

}


MinMaxD getMinMaxDistanceOfLeaf(int Node1, int Node2){
    MinMaxD mm;
    //loc 不在 Node_id所在的子图中
    vector<int> bs1 = GTree[Node1].borders;
    vector<int> bs2 = GTree[Node2].borders;
    float maxDis = 0;
    float minDis = 1000000;
    for(int i = 0;i<bs1.size();i++){
        int b1 = bs1[i];
        for(int i = 0;i<bs2.size();i++){
            int b2 = bs2[i];
            double distance = SPSP(b1, b2)/ WEIGHT_INFLATE_FACTOR;
            if(distance < minDis) minDis = distance;
            if(distance > maxDis) maxDis = distance;
        }
    }
    mm.min = minDis;
    mm.max = maxDis+GTree[Node1].maxEdgeDis+GTree[Node2].maxEdgeDis;
    return  mm;

}


MinMaxD getMinMaxDistanceOfNode(int loc, int Node_id){
    MinMaxD mm;

    //验证 loc 是否在 Node_id所在的子图中
    for(int i=0; i<Nodes[loc].gtreepath.size();i++){
        if(Nodes[loc].gtreepath[i] == Node_id) { //loc 在 Node_id所在的子图中
            int leaf_Of_loc = Nodes[loc].gtreepath.back();
            float mindis = GTree[leaf_Of_loc].minDis;
            mm.min = mindis;
            vector<int> bs = GTree[Node_id].borders;
            float maxDis = 0;
            float minDis = 1000000;
            for(int i = 0;i<bs.size();i++){
                int border_id = bs[i];
                double distance;
                //if(SPDMap[loc][border_id]!=NULL)  distance = SPDMap[loc][border_id];
                distance = SPSP(loc, border_id)/ WEIGHT_INFLATE_FACTOR;
                if(distance < minDis) minDis = distance;
                if(distance > maxDis) maxDis = distance;
            }
            mm.max = maxDis;
            return  mm;
        }
    }
    //loc 不在 Node_id所在的子图中
    vector<int> bs = GTree[Node_id].borders;
    float maxDis = 0;
    float minDis = 1000000;
    for(int i = 0;i<bs.size();i++){
        int border_id = bs[i];
        double distance;
        //if(SPDMap[loc][border_id]!=NULL)  distance = SPDMap[loc][border_id];
        distance = SPSP(loc, border_id)/ WEIGHT_INFLATE_FACTOR;
        if(distance < minDis) minDis = distance;
        if(distance > maxDis) maxDis = distance;
    }
    mm.min = minDis;
    mm.max = maxDis;
    return  mm;

}


MinMaxD getMinMaxDistanceBetweenNode(int poiNode, int usrNode){
    MinMaxD mm;

    //验证 两节点子图是否为同一子图
    if(poiNode == usrNode){

        vector<int> bs1 = GTree[poiNode].borders;
        vector<int> bs2 = GTree[usrNode].borders;
        float maxDis = 0; float maxDis2 = 0;
        float minDis = 1000000;
        for(int i = 0;i<bs2.size();i++){
            int b1 = bs2[i];
            for(int i = 0;i<bs2.size();i++){
                int b2 = bs2[i];
                float distance = SPSP(b1, b2)/ WEIGHT_INFLATE_FACTOR;
                if(distance > maxDis) maxDis = distance;
            }
        }
        mm.min = 0;
        mm.max = maxDis;
        return mm;
    }

    //验证 两节点的子图是否存在包含关系
    int current = poiNode; int father;
    while(current!=0) {
        father = GTree[current].father;
        if(father == usrNode){  //poiNode在usrNode的子图中
            mm.min =  GTree[poiNode].minDis;
            vector<int> bs1 = GTree[poiNode].borders;
            vector<int> bs2 = GTree[usrNode].borders;
            float maxDis = 0; float maxDis2 = 0;
            float minDis = 1000000;
            for(int i = 0;i<bs2.size();i++){
                int b1 = bs2[i];
                for(int i = 0;i<bs2.size();i++){
                    int b2 = bs2[i];
                    float distance = SPSP(b1, b2)/ WEIGHT_INFLATE_FACTOR;
                    if(distance < minDis) minDis = distance;
                    if(distance > maxDis) maxDis = distance;
                }
            }
            mm.max = maxDis+GTree[poiNode].maxEdgeDis;
            return mm;

        }
        current = father;
    }


    // poiNode 与 usrNode 非子图包含
    vector<int> bs1 = GTree[poiNode].borders;
    vector<int> bs2 = GTree[usrNode].borders;
    float maxDis = 0;
    float minDis = 1000000;
    for(int i = 0;i<bs1.size();i++){
        int b1 = bs1[i];
        for(int i = 0;i<bs2.size();i++){
            int b2 = bs2[i];
            double distance = SPSP(b1, b2)/ WEIGHT_INFLATE_FACTOR;
            if(distance < minDis) minDis = distance;
            if(distance > maxDis) maxDis = distance;
        }
    }

    float maxDis1 = 0;
    for(int i = 0;i<bs1.size();i++){
        int b1 = bs1[i];
        for(int i = 0;i<bs1.size();i++){
            int b2 = bs1[i];
            double distance = SPSP(b1, b2)/ WEIGHT_INFLATE_FACTOR;
            if(distance > maxDis1) maxDis1 = distance;
        }
    }

    float maxDis2 = 0;
    for(int i = 0;i<bs2.size();i++){
        int b1 = bs2[i];
        for(int i = 0;i<bs2.size();i++){
            int b2 = bs2[i];
            double distance = SPSP(b1, b2)/ WEIGHT_INFLATE_FACTOR;
            if(distance > maxDis2) maxDis2 = distance;
        }
    }

    mm.min = minDis;
    mm.max = maxDis+maxDis1+maxDis2;
    return  mm;

}

MinMaxD getMinMaxDistanceBetweenNode2(int Node1, int Node2){
    MinMaxD mm;
    float maxDis = 0;
    float minDis = 10000000;
    for(int v1:GTree[Node1].vetexSet){
        for(int v2: GTree[Node2].vetexSet){
            float distance = SPSP(v1, v2)/ WEIGHT_INFLATE_FACTOR;
            if(distance > maxDis) maxDis = distance;
            if(distance < minDis) minDis = distance;
        }
    }

    mm.min = minDis;
    mm.max = maxDis;
    return  mm;

}

//jins new
MinMaxD  getMinMaxDistanceWithinLeaf_slow(int leaf){
    MinMaxD mm;
    float maxDis = 0;
    float minDis = 10000000;

    vector<int> vertexes;
    for(int v :GTree[leaf].vetexSet)
        vertexes.push_back(v);
    //omp_set_num_threads(4);
//#pragma omp parallel for
    for(int i=0; i<vertexes.size();i++){
        int s = vertexes[i];
        vector<int> cands = vertexes;
        vector<int> distance_results;
        distance_results = dijkstra_candidate(s,cands,Nodes);
        for(int distance: distance_results){
            float distance_normalized = distance / 1000.0;
            maxDis = max(maxDis,distance_normalized);
        }
    }

    mm.min = 0;
    mm.max = maxDis;
    return  mm;
}


// jins new
MinMaxD getMinMaxDistanceBetweenLeaves_slow(int leaf1, int leaf2){
    MinMaxD mm;
    vector<int> bs1 = GTree[leaf1].borders;
    vector<int> bs2 = GTree[leaf2].borders;
    int b1 = -1; int b2 = -1;
    float maxDis = 0;
    float minDis = 9999999999;

    for(int _b1:bs1){
        for(int _b2: bs2){
            double distance = SPSP(_b1, _b2)/ WEIGHT_INFLATE_FACTOR;
            if(distance < minDis) {
                minDis = distance; b1 = _b1; b2 = _b2;
            }
        }
    }

    float maxDis1 = 0;
    vector<int> distance_results;
    vector<int> cands1;
    for(int v :GTree[leaf1].vetexSet)
        cands1.push_back(v);

    distance_results = dijkstra_candidate(b1,cands1,Nodes);

    for(int  distance1: distance_results){

        float distance_normalized1 = distance1/1.0;
        maxDis1 = max(distance_normalized1,maxDis1);

    }

    float maxDis2 = 0;
    vector<int> cands2;
    for(int v :GTree[leaf2].vetexSet)
        cands2.push_back(v);
    distance_results = dijkstra_candidate(b2,cands2,Nodes);

    for(int  distance2: distance_results){

        float distance_normalized2 = distance2/1.0;
        maxDis2 = max(distance_normalized2,maxDis2);

    }


    mm.min = minDis;
    mm.max = minDis+maxDis1+maxDis2;
    return  mm;

}



struct LeafDistanceBoundInfo{
    struct BorderPair borderPair;
    MinMaxD disPair;
};
LeafDistanceBoundInfo getMinMaxDistanceBetweenLeaves_BP(int leaf1, int leaf2){
    LeafDistanceBoundInfo boundInfo;
    MinMaxD mm;  BorderPair borderPair;
    vector<int> bs1 = GTree[leaf1].borders;
    vector<int> bs2 = GTree[leaf2].borders;
    int b1 = -1; int b2 = -1;
    float maxDis = 0;
    float minDis = 9999999999;

    for(int _b1:bs1){
        for(int _b2: bs2){
            double distance = SPSP(_b1, _b2)/ WEIGHT_INFLATE_FACTOR;
            if(distance < minDis) {
                minDis = distance; b1 = _b1; b2 = _b2;
                borderPair.b1 = _b1; borderPair.b2 = _b2;

            }
        }
    }

    float maxDis1 = 0;
    vector<int> distance_results;
    vector<int> cands1;
    for(int v :GTree[leaf1].vetexSet)
        cands1.push_back(v);

    distance_results = dijkstra_candidate(b1,cands1,Nodes);

    for(int  distance1: distance_results){

        float distance_normalized1 = distance1/1.0;
        maxDis1 = max(distance_normalized1,maxDis1);

    }

    float maxDis2 = 0;
    vector<int> cands2;
    for(int v :GTree[leaf2].vetexSet)
        cands2.push_back(v);
    distance_results = dijkstra_candidate(b2,cands2,Nodes);

    for(int  distance2: distance_results){

        float distance_normalized2 = distance2/1.0;
        maxDis2 = max(distance_normalized2,maxDis2);

    }
    float _tmp1 = maxDis1;
    float _tmp2 = maxDis2;

    mm.min = minDis;
    mm.max = minDis+maxDis1+maxDis2;

    boundInfo.borderPair = borderPair;
    boundInfo.disPair = mm;

    return  boundInfo;

}




class GIM_Tree{
public:
    //双色体对象(用户、兴趣点)

    vector<User> gimTree_bichromatic_Users;
    vector<POI> gimTree_bichromatic_POIs;
    //社交链接信息
    map<int,vector<int>> gimTree_social_map;
    map<int,map<int,int>> gimTree_social_table;
    Graph gimTree_graph;//整张社交网络图，其中改图中用户节点个数远大于双色体对象中用户个数
    //签到数据信息
    map<int,vector<int>> gimTree_userCheckInList;  //usr_id, vector<poi_id>
    map<int,vector<int>> gimTree_poiCheckInList;
    map<int,map<int,int>> gimTree_userCheckInMap;
    map<int,map<int,int>> gimTree_poiCheckInMap;

    //文本信息，倒排索引
    map<int, vector<int>> gimTree_invListOfUser;
    map<int, vector<int>> gimTree_invListOfPOI;
    //g-Tree索引信息
    vector<Node> gimTree_Nodes;
    vector<TreeNode> gimTree_GTree;
    //路网信息
    AdjType gimTree_adjList;

    //g-Tree各用户（节点）的关联信息
    map<int,map<int,vector<int>>> gimTree_leafUsrInv;
    map<int,set<int>> gimTree_usrTerm_leafSet;
    map<int,map<int,int>> gimTree_term_leaf_MinPopular;    // term_id, map<leaf_id, minPopular_Count>
    map<int,map<int,int>> gimTree_node_term_MinPopular;
    map<int,map<int,set<int>>> gimTree_term_leaf_KeySet;
    map<int,map<int,set<int>>> gimTree_usrNodeTermChild;
    //g-Tree各兴趣点（节点）的关联信息
    map<int,map<int,vector<int>>> gimTree_leafPoiInv;
    map<int,set<int>> gimTree_poiTerm_leafSet;
    map<int,map<int,set<int>>> gimTree_poiNodeTermChild;
    GIM_Tree(){

    }
    ~GIM_Tree(){

    }
};




// this is for distance bound pre-computation
int findLCA(int leaf1, int leaf2){
    int LCA_id =-1;
    int current_1 = GTree[leaf1].father;
    int current_2 = GTree[leaf2].father;
    if(current_1==current_2) return  current_1;
    int root =0;
    map<int,bool> flag_map;
    while(current_1!= root){
        flag_map[current_1]=true;
        current_1 = GTree[current_1].father;
    }
    while(current_2!= root){
        if(flag_map[current_2]==true){
            return  current_2;
        }
        current_2 = GTree[current_2].father;
    }
    return  root;

}


//废弃
MinMaxD  getMinMaxDistanceWithinLeaf_gtree_error(int leaf){
    MinMaxD mm;
    int max_border_vertex = -1;
    map<int, vector<int>> itm;  // vertex_id, border_id_idx, <...>
    for(int v1:GTree[leaf].vetexSet){
        int locid = v1;
        int posa = lower_bound(GTree[leaf].leafnodes.begin(), GTree[leaf].leafnodes.end(), locid) -
               GTree[leaf].leafnodes.begin();

        for (int j = 0; j < GTree[leaf].borders.size(); j++) {
            itm[v1].push_back(GTree[leaf].mind[j * GTree[leaf].leafnodes.size() + posa]);
        }
    }
    for(int v1:GTree[leaf].vetexSet){
        for(int v2:GTree[leaf].vetexSet){
            int maxDis = 0;
            int minDis = -1;
            for(int i=0; i<itm[v1].size();i++){
                int dis_sum = itm[v1][i]+itm[v2][i];
                if(minDis ==-1) minDis = dis_sum;
                else minDis = min(minDis,dis_sum);
            }
            max_border_vertex = max(minDis,max_border_vertex);

        }
    }

    mm.min = 0;
    mm.max = max_border_vertex;
    return  mm;
}


MinMaxD  getMinMaxDistanceWithinLeaf_gtree(int leaf){
    MinMaxD mm;
    int max_border_vertex = -1;
    map<int, vector<int>> itm;  // vertex_id, border_id_idx, <...>
    for(int v1:GTree[leaf].vetexSet){
        int locid = v1;
        int posa = lower_bound(GTree[leaf].leafnodes.begin(), GTree[leaf].leafnodes.end(), locid) -
                   GTree[leaf].leafnodes.begin();

        for (int j = 0; j < GTree[leaf].borders.size(); j++) {
            itm[v1].push_back(GTree[leaf].mind[j * GTree[leaf].leafnodes.size() + posa]);
        }
    }
    for(int v1:GTree[leaf].vetexSet){
        //for(int v2:GTree[leaf].vetexSet){
            int maxDis = 0;
            int minDis = -1;
            for(int i=0; i<itm[v1].size();i++){
                int Dis = itm[v1][i];//+itm[v2][i];
                //if(minDis ==-1) minDis = dis_sum;
                // minDis = min(minDis,dis_sum);
                max_border_vertex = max(Dis,max_border_vertex);
            }


        //}
    }

    mm.min = 0;
    mm.max = max_border_vertex*2;
    return  mm;
}


MinMaxD getMinMaxDistanceBetweenLeaves_gtree(int leaf1, int leaf2, map<int, unordered_map<int, vector<int>>>& itm_global, map<int,int>& leaf_in){
    int maxDis = 0;
    int minDis = -1;
    int LCA = findLCA(leaf1, leaf2);
    //cout<<"LCA=n"<<LCA<<endl;
    int min_dis_border = -1;
    MinMaxD mm;
    for(int border1: GTree[leaf1].borders){
        for(int border2: GTree[leaf2].borders){
            int i1 = itm_global[border1][LCA].size();
            int i2 = itm_global[border2][LCA].size();
            if(i1!=i2){
                cout<<"有错误！"<<endl;
                mm.max=-1;
                return mm;
            }
            //计算border与border间的最小距离
            for(int i=0;i<i1;i++){
                int sum_dis = itm_global[border1][LCA][i]+itm_global[border2][LCA][i];
                if(min_dis_border==-1) min_dis_border = sum_dis;
                else{
                    min_dis_border = min(min_dis_border,sum_dis);
                }
            }
            //更新最大边界距离
            maxDis = max(maxDis, min_dis_border);
            //更新最小边界距离
            if(minDis==-1)  minDis = min_dis_border;
            else minDis = min(minDis,min_dis_border);
        }
    }

    mm.max = (maxDis+leaf_in[leaf1]+leaf_in[leaf2]) / WEIGHT_INFLATE_FACTOR;
    mm.min = minDis / WEIGHT_INFLATE_FACTOR;
    return mm;
}





//dibu

#endif

