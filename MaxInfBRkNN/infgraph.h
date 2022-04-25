#ifndef Ingraph_H
#define Ingraph_H
#include <utility>
#include <algorithm>
using namespace std;

typedef pair<double,int> dipair;

typedef struct {
    double local_inf;
    double outer_inf;
} Inf_Pair;


#include <queue>
#include <utility>
// /SFMT/dSFMT/dSFMT.c
#include "SFMT/dSFMT/dSFMT.h"
#include "./SocialInfluenceModule/IMM/argument.h"
#include "./SocialInfluenceModule/IMM/discrete_rrset.h"
//#include "./SocialInfluenceModule/IMM/iheap.h"

//#include "./SocialInfluenceModule/IMM/head.h"
#include "./SocialInfluenceModule/IMM/graph.h"

struct CompareBySecond {
    bool operator()(pair<int, int> a, pair<int, int> b)
    {
        return a.second < b.second;
    }
};

struct CompareBySecondDouble{
    bool operator()(pair<int, double> a, pair<int, double> b)
    {
        return a.second < b.second;
    }
};


enum RRSETS{ R1, R2, R3};

class InfGraph  {  //public SocialGraph

private:
    vector<bool> visit;
    vector<int> visit_mark;


public:

    SocialGraph* social_graph = NULL;

    long long influence_upper = -1;
    long long outer_influence_only_upper = -1;
    double local_influence_upper = -1;
    double outer_influence_upper = -1;
    double overall_influence_upper = -1;
    Inf_Pair overall_influence_upperDetails;
    long long influence_upper2 = -1;
    vector<vector<int>> hyperG;
    vector<vector<int>> hyperGT;  //reverse

    vector<vector<int>> hyperGVld;
    vector<vector<int>> hyperGTVld;  //reverse


    vector<vector<int>> hyperGPre;
    vector<vector<int>> hyperGTPre;  //for pre-sampled RIS


    //for poi influence
    vector<vector<int>> hyperG_POI;
    vector<vector<int>> hyperGT_POI;  //reverse
    vector<vector<int>> hyperG_POIVld;
    vector<vector<int>> hyperGT_POIVld;  //reverse

    //for outer range influence of poi
    vector<vector<int>> hyperG_POI_OUT;
    vector<vector<int>> hyperGT_POI_OUT;  //reverse
    vector<vector<int>> hyperG_POIVld_OUT;
    vector<vector<int>> hyperGT_POIVld_OUT;  //reverse


    vector<vector<int>> hyperG_Upper_POI;
    vector<vector<int>> hyperGT_Upper_POI;  //reverse


    InfGraph(int nodes, int pois, SocialGraph* graph){
        social_graph = graph;
        int num = nodes;
        init_hyper_graph_poi(pois);
        if(social_graph->n == nodes){
            visit = vector<bool> (nodes);
            visit_mark = vector<int> (nodes);
        } else{
            cout<<"n!=nodes"<<endl;
            exit(-1);
        }

    }

    InfGraph(int nodes, SocialGraph* graph){
        social_graph = graph;
        int num = nodes;
        init_hyper_graph();
        if(social_graph->n==nodes){
            visit = vector<bool> (nodes);
            visit_mark = vector<int> (nodes);
        } else{
            cout<<"n!=nodes"<<endl;
            exit(-1);
        }

    }

    void freeMemory(){
        cout<<"free Memeory"<<endl;
        cout<<"hyperG_POI num="<<hyperG_POI.capacity()<<endl;
        vector<vector<int>> _temp;
        _temp.swap(hyperG_POI);
        cout<<"1_temp num="<<_temp.capacity()<<endl;
        //getchar();
        //delete _temp;


        vector<vector<int>>(). swap(hyperGT_POI);

        //cout<<"2_temp num="<<_temp.capacity()<<endl;

        vector<vector<int>>().swap(hyperG_POIVld);

        //cout<<"3_temp num="<<_temp.capacity()<<endl;

        vector<vector<int>>().swap(hyperGT_POIVld);

        //cout<<"4_temp num="<<_temp.capacity()<<endl;

        //getchar();



    }


    ~InfGraph(){
        printf("~begin InfGraph()\n");


    }

    void add_edge(int u,int v,double prob){
        social_graph->add_edge(u,v,prob);
    }

    void setRearchable(int u,int p){
        social_graph->setRearchable(u,p);
    }

    void init_hyper_graph(){
        hyperG.clear();
        for (int i = 0; i < social_graph->n; i++)
            hyperG.push_back(vector<int>());
        hyperGT.clear();
    }

    void init_hyper_graph_vld(){
        hyperGVld.clear();
        for (int i = 0; i < social_graph->n; i++)
            hyperGVld.push_back(vector<int>());
        hyperGTVld.clear();
    }


// this is for poi
    void init_hyper_graph_poi(int p){
        hyperG_POI.clear();
        for (int i = 0; i < p; i++)
            hyperG_POI.push_back(vector<int>());
        hyperGT_POI.clear();
    }

    void init_Two_hyper_graph_poi(int p){
        hyperG_POI.clear();
        hyperG_POIVld.clear();
        for (int i = 0; i < p; i++){
            hyperG_POI.push_back(vector<int>());
            hyperG_POIVld.push_back(vector<int>());
        }
        hyperGT_POI.clear();
        hyperGT_POIVld.clear();
    }

    void init_Two_hyper_graph_poi_out(int p){
        hyperG_POI_OUT.clear();
        hyperG_POIVld_OUT.clear();
        for (int i = 0; i < p; i++){
            hyperG_POI_OUT.push_back(vector<int>());
            hyperG_POIVld_OUT.push_back(vector<int>());
        }
        hyperGT_POI_OUT.clear();
        hyperGT_POIVld_OUT.clear();
    }

    void init_Three_hyper_graph_poi(int p, int p_u){
        hyperG_POI.clear();
        hyperG_POIVld.clear();
        hyperG_Upper_POI.clear();
        for (int i = 0; i < p; i++){
            hyperG_POI.push_back(vector<int>());
            hyperG_POIVld.push_back(vector<int>());
        }
        for(int j = 0; j< p_u; j++){
            hyperG_Upper_POI.push_back(vector<int>());
        }
        hyperGT_POI.clear();
        hyperGT_POIVld.clear();
        hyperGT_Upper_POI.clear();
    }





    int get_RR_sets_size(){
        return  hyperGT.size();
    }

    int get_RR_POI_sets_size(){
        return  hyperGT_POI.size();
    }

    int get_RR_POI_sets_outer_size(){
        return  hyperGT_POI_OUT.size();
    }



    int get_RR_POIUpper_sets_size(){  // this is for RR set of upper constraint graph
        return  hyperGT_Upper_POI.size();
    }


    /*
 * BFS starting from one node
 */
    deque<int> q;
    vector<int> seedSet;
    vector<int> seedSet2;  //superset
    vector<int> poiSet;
    vector<int> poiSet2;   //superset

    void printPOISelectionInfo(){
        for(int p_id:poiSet)
            cout<<"p"<<p_id<<",";
        cout<<endl;
    }

    int BuildHypergraphNode(int uStart, int hyperiiid, RRSETS which)
    {

        int n_visit_edge = 1;

        if(which==R1)
            hyperGT[hyperiiid].push_back(uStart);
        else if(which==R2)
            hyperGTVld[hyperiiid].push_back(uStart);


        int n_visit_mark = 0;

        q.clear();
        q.push_back(uStart);
        //ASSERT(n_visit_mark < n);
        visit_mark[n_visit_mark] = uStart;
#ifdef RIS_TRACK
        if(uStart > visit.size()){
            cout<<uStart<<"> visit.size()="<<visit.size()<<endl;
            getchar();
        }
#endif
        n_visit_mark++;
        visit[uStart] = true;
        while (!q.empty())
        {

            int expand = q.front();
            q.pop_front();

            if (true)  //influModel == IC
            {
                int i = expand;
#ifdef RIS_TRACK
                cout<<"expand u"<<expand<<" has "<<social_graph->in_neighbors[i].size()<<" friends, q.size="<<q.size()<< endl;
#endif
                //i = 44063;


                for (int j = 0; j < (int)social_graph->in_neighbors[i].size(); j++)
                {
                    //int u=expand;
                    int v = social_graph->in_neighbors[i][j];  //zheng
                    //cout<<"v"<<v;
                    n_visit_edge++;
                    //double randDouble = 0.03;// dsfmt_gv_genrand_uint32_range(1001)/1000.0;
                    double randDouble = dsfmt_gv_genrand_uint32_range(1001)/1000.0;
#ifdef RIS_TRACK
                    cout<<" randDouble= "<<randDouble<<endl;
                    cout<<" social_graph->probAR["<<i<<"]["<<j<<"]="<<social_graph->probAR[i][j]<<endl;
#endif
                    if (randDouble > social_graph->probAR[i][j])
                        continue;
                    if (visit[v])
                        continue;
                    if (!visit[v])
                    {
                        //ASSERT(n_visit_mark < n);
#ifdef RIS_TRACK
                        cout<<"n_visit_mark="<<n_visit_mark<<endl;
#endif
                        visit_mark[n_visit_mark] = v;
                        n_visit_mark++;
                        visit[v] = true;
#ifdef RIS_TRACK
                        if(v > visit.size()){
                            cout<<"v"<<v<<" > visit.size()="<<visit.size()<<endl;
                            getchar();
                        }
#endif


                    }
                    q.push_back(v);
                    //ASSERT((int)hyperGT.size() > hyperiiid);
#ifdef RIS_TRACK
                    //cout<<"hyperGT["<<hyperiiid<<"] success!"<<endl;
                    //cout<<"hyperGT.size()="<<hyperGT.size()<<", hyperiiid="<<hyperiiid<<endl;
#endif
                    if(which==R1){
                        hyperGT[hyperiiid].push_back(v);
#ifdef RIS_TRACK
                        cout<<"push u"<<v<<" into RR"<<hyperiiid<<endl;
#endif
                    }

                    else if(which==R2){
                        hyperGTVld[hyperiiid].push_back(v);
#ifdef RIS_TRACK
                        cout<<"push u"<<v<<" into RR_vld"<<hyperiiid<<endl;
#endif
                    }

#ifdef RIS_TRACK
                    cout<<"hyperGT["<<hyperiiid<<"] success!"<<endl;
#endif

                }
            }
        }
#ifdef RIS_TRACK
        cout<<"复位"<<endl;
#endif
        if(n_visit_mark>0){
#ifdef RIS_TRACK
            cout<<"n_visit_mark="<<n_visit_mark<<endl;
#endif
            for (int i = 0; i < n_visit_mark; i++){
                //cout<<"i="<<i<<endl;
                //cout<<"try to reflag u"<<visit_mark[i]<<endl;
                visit[visit_mark[i]] = false;
                //cout<<"reflag u"<<visit_mark[i]<<"success"<<endl;
            }
        }
#ifdef RIS_TRACK
        cout<<"end!"<<endl;
#endif

        return n_visit_edge;
    }




    int BuildHypergraphPOI(int u_idx, int hyperiiid, RRSETS which)
    {
        //getchar();
        //int uStart0 = reachable_list[u_idx];   //确定第u_idx个 rearchable user的id

        int uStart = social_graph->idx_usr_table[u_idx];

        //cout<<"get reachable_list["<<u_idx<<"]"<<"=u"<<uStart0<<endl;

        //cout<<"get idx_usr_table["<<u_idx<<"]"<<"=u"<<uStart<<endl;

        int n_visit_edge = 1;


        int n_visit_mark = 0;

        q.clear();
        q.push_back(uStart);
        //ASSERT(n_visit_mark < n);
        visit_mark[n_visit_mark++] = uStart;
        visit[uStart] = true;
        while (!q.empty())
        {

            int expand = q.front();
            q.pop_front();
            if (true)  //influModel == IC
            {
                int i = expand;
                //cout<<"expand from usr_id="<<i<<endl;
                //先检测该user是否有反向影响力(u-p)可达 poi
                if(social_graph->reachable_poi_table[i].size()>0){
                    //cout<<"u"<<i<<"is a potential user "<<endl;
                    for(int p: social_graph->reachable_poi_table[i]){
                        //cout<<"hyperiiid="<<hyperiiid<<",p_id= "<<p<<endl;
                        //cout<<"size="<<hyperGT_POI.size()<<endl;
                        if(which==R1)
                            hyperGT_POI[hyperiiid].push_back(p);
                        else if(which == R2)
                            hyperGT_POIVld[hyperiiid].push_back(p);

                    }

                }
                //其次沿反向 social link 向外延伸
                for (int j = 0; j < (int)social_graph->in_neighbors[i].size(); j++)
                {

                    //int u=expand;
                    int v = social_graph->in_neighbors[i][j];  //zheng

                    double randDouble = dsfmt_gv_genrand_uint32_range(1001)/1000.0;
                    if (randDouble > social_graph->probAR[i][j])
                        continue;
                    if (visit[v])
                        continue;
                    if (!visit[v])
                    {
                        visit_mark[n_visit_mark++] = v;
                        visit[v] = true;
                        //cout<<"to u"<<v<<endl;
                    }
                    q.push_back(v);
                }
            }
        }
        for (int i = 0; i < n_visit_mark; i++)
            visit[visit_mark[i]] = false;
        return n_visit_edge;
    }


    int BuildHypergraphPOI_Outer(int u_idx, int hyperiiid, RRSETS which)
    {
        //int uStart0 = reachable_list[u_idx];   //确定第u_idx个 rearchable user的id

        int uStart = social_graph->idx_usr_table[u_idx];


        int n_visit_edge = 1;


        int n_visit_mark = 0;

        q.clear();
        deque<int> q_local;q_local.clear();
        //q.push_back(uStart);
        //ASSERT(n_visit_mark < n);
        visit_mark[n_visit_mark++] = uStart;
        visit[uStart] = true;

        //先在影响力领域范围内进行一步 BFS
        for (int j = 0; j < (int)social_graph->in_neighbors[uStart].size(); j++)
        {

            //int u=expand;
            int v = social_graph->in_neighbors[uStart][j];  //zheng
            double randDouble = dsfmt_gv_genrand_uint32_range(1001)/1000.0;
            if (randDouble > social_graph->probAR[uStart][j])
                continue;
            if (visit[v])
                continue;
            if (!visit[v])
            {
                visit_mark[n_visit_mark++] = v;
                visit[v] = true;
                //cout<<"to u"<<v<<endl;
            }
            q_local.push_back(v);

        }
        //在影响力领域范围内进行二步 BFS
        while (!q_local.empty()){
            int i = q_local.front();
            q_local.pop_front();
            //其次沿反向 social link 向外延伸
            for (int j = 0; j < (int)social_graph->in_neighbors[i].size(); j++)
            {

                //int u=expand;
                int v = social_graph->in_neighbors[i][j];  //zheng

                double randDouble = dsfmt_gv_genrand_uint32_range(1001)/1000.0;
                if (randDouble > social_graph->probAR[i][j])
                    continue;
                if (visit[v])
                    continue;
                if (!visit[v])
                {
                    visit_mark[n_visit_mark++] = v;
                    visit[v] = true;
                    //cout<<"to u"<<v<<endl;
                }
                q.push_back(v);
            }
        }

        //进行外部影响力范围中的BFS
        while (!q.empty()) {

            int expand = q.front();
            q.pop_front();
            if (true)  //influModel == IC
            {
                int i = expand;
                //cout<<"expand from usr_id="<<i<<endl;
                //先检测该user是否有反向影响力(u-p)可达 poi
                if(social_graph->reachable_poi_table[i].size()>0){
                    //cout<<"u"<<i<<"is a potential user "<<endl;
                    for(int p: social_graph->reachable_poi_table[i]){
                        //if(i==uStart) continue;  //领域直接可达poi, 该情况忽略
                        //cout<<"hyperiiid="<<hyperiiid<<",p_id= "<<p<<endl;
                        //cout<<"size="<<hyperGT_POI.size()<<endl;
                        if(which==R1)
                            hyperGT_POI_OUT[hyperiiid].push_back(p);
                        else if(which == R2)
                            hyperGT_POIVld_OUT[hyperiiid].push_back(p);

                    }

                }
                //其次沿反向 social link 向外延伸
                for (int j = 0; j < (int)social_graph->in_neighbors[i].size(); j++)
                {

                    //int u=expand;
                    int v = social_graph->in_neighbors[i][j];  //zheng

                    double randDouble = dsfmt_gv_genrand_uint32_range(1001)/1000.0;
                    if (randDouble > social_graph->probAR[i][j])
                        continue;
                    if (visit[v])
                        continue;
                    if (!visit[v])
                    {
                        visit_mark[n_visit_mark++] = v;
                        visit[v] = true;
                        //cout<<"to u"<<v<<endl;
                    }
                    q.push_back(v);
                }
            }
        }
        for (int i = 0; i < n_visit_mark; i++)
            visit[visit_mark[i]] = false;
        return n_visit_edge;
    }


    int BuildHypergraphPOI_Upper(int u_idx, int hyperiiid)
    {
        //int uStart = reachable_list_upper[u_idx];   //确定第u_idx个 rearchable user的id

        int uStart = social_graph->idx_usr_table_upper[u_idx];


        int n_visit_edge = 1;


        int n_visit_mark = 0;

        q.clear();
        q.push_back(uStart);
        //ASSERT(n_visit_mark < n);
        visit_mark[n_visit_mark++] = uStart;
        visit[uStart] = true;
        while (!q.empty())
        {

            int expand = q.front();
            q.pop_front();
            if (true)  //influModel == IC
            {
                int i = expand;
                //先检测该user是否有反向影响力(u-p)可达 poi
                if(social_graph->reachable_poi_table_upper[i].size()>0){
                    for(int p: social_graph->reachable_poi_table_upper[i]){
                        //if(which==R1)
                        hyperGT_Upper_POI[hyperiiid].push_back(p);
                        //else if(which == R2)
                        //hyperGT_POIVld[hyperiiid].push_back(p);
                    }
                }
                //其次沿反向 social link 向外延伸
                for (int j = 0; j < (int)social_graph->in_neighbors[i].size(); j++)
                {
                    //int u=expand;
                    int v = social_graph->in_neighbors[i][j];  //zheng
                    double randDouble = dsfmt_gv_genrand_uint32_range(1001)/1000.0;
                    if (randDouble > social_graph->probAR[i][j])
                        continue;
                    if (visit[v])
                        continue;
                    if (!visit[v])
                    {
                        visit_mark[n_visit_mark++] = v;
                        visit[v] = true;
                    }
                    q.push_back(v);
                }
            }
        }
        for (int i = 0; i < n_visit_mark; i++)
            visit[visit_mark[i]] = false;
        return n_visit_edge;
    }




    void build_hyper_graph_r(int R, const Argument & arg)
    {

        int prevSize = hyperGT.size();
        while ((int)hyperGT.size() <= R)
            hyperGT.push_back( vector<int>() );


        vector<int> random_number;
        for (int i = 0; i < R; i++)
        {
            int random = dsfmt_gv_genrand_uint32_range(social_graph->n);  //gai
            random_number.push_back(random);
        }
#ifdef RIS_TRACK
        cout<<"social_graph->n="<<social_graph->n<<endl;
        cout<<"2"<<endl;
#endif

        //trying BFS start from same node

        for (int i = prevSize; i < R; i++)
        {
#ifdef RIS_TRACK
            cout<<"build rr set"<<i<<"for u"<<random_number[i]<<endl;
#endif
            BuildHypergraphNode(random_number[i], i, R1);
#ifdef RIS_TRACK
            cout<<"build succeed!"<<i<<endl;
#endif

        }

        int totAddedElement = 0;
        for (int i = prevSize; i < R; i++)
        {
            for (int t : hyperGT[i])
            {
                hyperG[t].push_back(i);
                //hyperG.addElement(t, i);
                totAddedElement++;
            }
        }
    }


    void build_hyper_graph_vld_r(int R, const Argument & arg)
    {

        int prevSize = hyperGTVld.size();
        while ((int)hyperGTVld.size() <= R)
            hyperGTVld.push_back( vector<int>() );
        //cout<<"1"<<endl;


        vector<int> random_number;
        for (int i = 0; i < R; i++)
        {
            int random = dsfmt_gv_genrand_uint32_range(social_graph->n);  //gai
            random_number.push_back(random);
        }
        //cout<<"social_graph->n="<<social_graph->n<<endl;
        //cout<<"2"<<endl;


        //trying BFS start from same node

        for (int i = prevSize; i < R; i++)
        {
            //cout<<"build rr set"<<i<<"for u"<<random_number[i]<<endl;
            BuildHypergraphNode(random_number[i], i, R2);
            //cout<<"build succeed!"<<i<<endl;


        }
        //cout<<"3"<<endl;

        int totAddedElement = 0;
        for (int i = prevSize; i < R; i++)
        {
            for (int t : hyperGTVld[i])
            {
                hyperGVld[t].push_back(i);
                //hyperG.addElement(t, i);
                totAddedElement++;
            }
        }
        //cout<<"4"<<endl;
    }


    void build_hyper_graph_r_poi(int R, const Argument & arg)
    {

        int prevSize = hyperGT_POI.size();
        while ((int)hyperGT_POI.size() <= R)
            hyperGT_POI.push_back(vector<int>() );


        vector<int> random_number;
        for (int i = 0; i < R; i++)
        {
            int random = dsfmt_gv_genrand_uint32_range(social_graph->getNc());  //gai
            random_number.push_back(random);
        }

        //trying BFS start from same node

        for (int i = prevSize; i < R; i++)
        {
            int u_id = social_graph->idx_usr_table[random_number[i]] ;

            BuildHypergraphPOI(u_id, i,R1);  //i:第几个 rearchable_user

        }

        int totAddedElement = 0;
        for (int i = prevSize; i < R; i++)  //根据rr set 构建 fr set(正向采样集)
        {
            for (int p : hyperGT_POI[i])  //这里的p为兴趣点的id
            {
                int p_idx = social_graph->poi_idx_table[p];
                hyperG_POI[p_idx].push_back(i);  //key是 poi_id, value 是 R id
                //hyperG.addElement(t, i);
                totAddedElement++;
            }
        }
    }


    void build_hyper_graph_vld_r_poi(int R, const Argument & arg)
    {

        int prevSize = hyperGT_POIVld.size();
        while ((int)hyperGT_POIVld.size() <= R)
            hyperGT_POIVld.push_back(vector<int>() );


        vector<int> random_number;
        for (int i = 0; i < R; i++)
        {
            int random = dsfmt_gv_genrand_uint32_range(social_graph->getNc());  //gai
            random_number.push_back(random);
        }

        //trying BFS start from same node

        for (int i = prevSize; i < R; i++)
        {
            int u_id = social_graph->idx_usr_table[random_number[i]] ;

            BuildHypergraphPOI(u_id, i,R2);  //i:第几个 rearchable_user

        }

        int totAddedElement = 0;
        for (int i = prevSize; i < R; i++)  //根据rr set 构建 fr set(正向采样集)
        {
            for (int p : hyperGT_POIVld[i])  //这里的p为兴趣点的id
            {
                int p_idx = social_graph->poi_idx_table[p];
                hyperG_POIVld[p_idx].push_back(i);  //key是 poi_id, value 是 R id
                //hyperG.addElement(t, i);
                totAddedElement++;
            }
        }
    }

/*---------------------below this for outer range influence scenario--------------------------------*/

    void build_hyper_graph_r_poi_outer(int R, const Argument & arg)
    {

        int prevSize = hyperGT_POI_OUT.size();
        while ((int)hyperGT_POI_OUT.size() <= R)
            hyperGT_POI_OUT.push_back(vector<int>() );


        vector<int> random_number;
        for (int i = 0; i < R; i++)
        {
            int random = dsfmt_gv_genrand_uint32_range(social_graph->getNc());  //gai
            random_number.push_back(random);
        }

        //trying BFS start from same node

        for (int i = prevSize; i < R; i++)
        {
            int u_id = social_graph->idx_usr_table[random_number[i]] ;

            BuildHypergraphPOI_Outer(u_id, i,R1);  //i:第几个 rearchable_user

        }

        int totAddedElement = 0;
        for (int i = prevSize; i < R; i++)  //根据rr set 构建 fr set(正向采样集)
        {
            for (int p : hyperGT_POI_OUT[i])  //这里的p为兴趣点的id
            {
                int p_idx = social_graph->poi_idx_table[p];
                hyperG_POI_OUT[p_idx].push_back(i);  //key是 poi_id, value 是 R id
                //hyperG.addElement(t, i);
                totAddedElement++;
            }
        }
    }

    void build_hyper_graph_vld_r_poi_outer(int R, const Argument & arg)
    {

        int prevSize = hyperGT_POIVld_OUT.size();
        while ((int)hyperGT_POIVld_OUT.size() <= R)
            hyperGT_POIVld_OUT.push_back(vector<int>() );


        vector<int> random_number;
        for (int i = 0; i < R; i++)
        {
            int random = dsfmt_gv_genrand_uint32_range(social_graph->getNc());  //gai
            random_number.push_back(random);
        }

        //trying BFS start from same node

        for (int i = prevSize; i < R; i++)
        {
            int u_id = social_graph->idx_usr_table[random_number[i]] ;

            BuildHypergraphPOI_Outer(u_id, i,R2);  //i:第几个 rearchable_user

        }

        int totAddedElement = 0;
        for (int i = prevSize; i < R; i++)  //根据rr set 构建 fr set(正向采样集)
        {
            for (int p : hyperGT_POIVld_OUT[i])  //这里的p为兴趣点的id
            {
                int p_idx = social_graph->poi_idx_table[p];
                hyperG_POIVld_OUT[p_idx].push_back(i);  //key是 poi_id, value 是 R id
                //hyperG.addElement(t, i);
                totAddedElement++;
            }
        }
    }
/*---------------------below this for outer range influence scenario--------------------------------*/

    void build_hyper_graph_upper_r_poi(int R, const Argument & arg)
    {

        int prevSize = hyperGT_Upper_POI.size();
        while ((int)hyperGT_Upper_POI.size() <= R)
            hyperGT_Upper_POI.push_back(vector<int>() );


        vector<int> random_number;
        for (int i = 0; i < R; i++)
        {
            int random = dsfmt_gv_genrand_uint32_range(social_graph->getNc_upper());  //gai
            random_number.push_back(random);
        }

        //trying BFS start from same node

        for (int i = prevSize; i < R; i++)
        {
            int u_id = social_graph->idx_usr_table_upper[random_number[i]] ;

            BuildHypergraphPOI_Upper(u_id,i);  //i:第几个 rearchable_user

        }

        int totAddedElement = 0;
        for (int i = prevSize; i < R; i++)  //根据rr set 构建 fr set(正向采样集)
        {
            for (int p : hyperGT_Upper_POI[i])  //这里的p为兴趣点的id
            {
                int p_idx = social_graph->poi_idx_table_upper[p];
                hyperG_Upper_POI[p_idx].push_back(i);  //key是 poi_id, value 是 R id
                //hyperG.addElement(t, i);
                totAddedElement++;
            }
        }
    }


    // for incremental

    void build_hyper_graph_r_poi_from_scratch(int R, const Argument & arg)
    {

        hyperGT_POI.clear();
        int prevSize = hyperGT_POI.size();
        while ((int)hyperGT_POI.size() <= R)
            hyperGT_POI.push_back(vector<int>() );


        vector<int> random_number;
        for (int i = 0; i < R; i++)
        {
            int random = dsfmt_gv_genrand_uint32_range(social_graph->getNc());  //gai
            int size2 = social_graph->reachable_usrs.size();
            if(size2!=social_graph->nc){
                cout<<"error!"<<endl;
                getchar();
            }

            random_number.push_back(random);
        }

        //trying BFS start from same node

        for (int i = prevSize; i < R; i++)
        {
            int idx = random_number[i];
            int u_id = social_graph->idx_usr_table[idx] ;
            //cout<<"i="<<i<<",in_neighbors.size()="<<in_neighbors.size()<<",reachable_usrs size="<<reachable_usrs.size()<<endl;
            //cout<<"idx="<<idx<<",build RR sets for u"<<u_id<<"(the "<<i<<"th RR set)"<<endl;

            //BuildHypergraphPOI(u_id, i,R1);  //i:第几个 rearchable_user
            BuildHypergraphPOI(idx, i,R1);  //i:第几个 rearchable_user

        }

        int totAddedElement = 0;
        for (int i = prevSize; i < R; i++)  //根据rr set 构建 fr set(正向采样集)
        {
            for (int p : hyperGT_POI[i])  //这里的p为兴趣点的id
            {
                int p_idx = social_graph->poi_idx_table[p];
                hyperG_POI[p_idx].push_back(i);  //key是 poi_id, value 是 R id
                //hyperG.addElement(t, i);
                totAddedElement++;
            }
        }
    }


    void build_hyper_graph_vld_r_poi_from_scratch(int R, const Argument & arg)
    {

        hyperGT_POIVld.clear();
        int prevSize = hyperGT_POIVld.size();
        while ((int)hyperGT_POIVld.size() <= R)
            hyperGT_POIVld.push_back(vector<int>() );


        vector<int> random_number;
        for (int i = 0; i < R; i++)
        {
            int random = dsfmt_gv_genrand_uint32_range(social_graph->getNc());  //gai
            random_number.push_back(random);
        }

        //trying BFS start from same node

        for (int i = prevSize; i < R; i++)
        {
            int u_id = social_graph->idx_usr_table[random_number[i]] ;

            BuildHypergraphPOI(u_id, i,R2);  //i:第几个 rearchable_user

        }

        int totAddedElement = 0;
        for (int i = prevSize; i < R; i++)  //根据rr set 构建 fr set(正向采样集)
        {
            for (int p : hyperGT_POIVld[i])  //这里的p为兴趣点的id
            {
                int p_idx = social_graph->poi_idx_table[p];
                hyperG_POIVld[p_idx].push_back(i);  //key是 poi_id, value 是 R id
                //hyperG.addElement(t, i);
                totAddedElement++;
            }
        }
    }

    void build_hyper_graph_upper_r_poi_from_scratch(int R, const Argument & arg)
    {

        hyperGT_Upper_POI.clear();
        int prevSize = hyperGT_Upper_POI.size();
        while ((int)hyperGT_Upper_POI.size() <= R)
            hyperGT_Upper_POI.push_back(vector<int>() );


        vector<int> random_number;
        for (int i = 0; i < R; i++)
        {
            int random = dsfmt_gv_genrand_uint32_range(social_graph->getNc_upper());  //gai
            random_number.push_back(random);
        }

        //trying BFS start from same node

        for (int i = prevSize; i < R; i++)
        {
            int u_id = social_graph->idx_usr_table_upper[random_number[i]] ;

            BuildHypergraphPOI_Upper(u_id,i);  //i:第几个 rearchable_user

        }

        int totAddedElement = 0;
        for (int i = prevSize; i < R; i++)  //根据rr set 构建 fr set(正向采样集)
        {
            for (int p : hyperGT_Upper_POI[i])  //这里的p为兴趣点的id
            {
                int p_idx = social_graph->poi_idx_table_upper[p];
                hyperG_Upper_POI[p_idx].push_back(i);  //key是 poi_id, value 是 R id
                //hyperG.addElement(t, i);
                totAddedElement++;
            }
        }
    }






    pair<double,double> build_seedset2(int k)    // inf,  Cover_OPT_UPPERBOUND
    {

        priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap;
        vector<int>coverage(social_graph->n, 0);

        for (int i = 0; i < social_graph->n; i++)
        {
            pair<int, int>tep(make_pair(i, (int)hyperG[i].size()));  //i: node_id,  hyperG: FRset{ rr set_id...}
            heap.push(tep);
            coverage[i] = (int)hyperG[i].size();
        }

        int maxInd;

        long long influence = 0;
        long long numEdge = hyperGT.size();

        // check if an edge is removed
        vector<bool> edgeMark(numEdge, false);
        // check if an node is remained in the heap
        vector<bool> nodeMark(social_graph->n + 1, true);

        seedSet.clear();
        while ((int)seedSet.size()<k)
        {
            pair<int, int>ele = heap.top();
            heap.pop();
            if (ele.second > coverage[ele.first])
            {
                ele.second = coverage[ele.first];
                heap.push(ele);
                continue;
            }

            maxInd = ele.first;
            vector<int>e = hyperG[maxInd];
            influence += coverage[maxInd];
            seedSet.push_back(maxInd);
            nodeMark[maxInd] = false;

            for (unsigned int j = 0; j < e.size(); ++j){
                if (edgeMark[e[j]])continue;

                vector<int>nList = hyperGT[e[j]];
                for (unsigned int l = 0; l < nList.size(); ++l){
                    if (nodeMark[nList[l]])coverage[nList[l]]--;
                }
                edgeMark[e[j]] = true;
            }
        }
        double inf = 1.0*influence / hyperGT.size();
        double cover_upper = 0.0;
        pair<int, int> ratio_pair (make_pair(inf, cover_upper));
        return ratio_pair;
    }

    double build_seedset(int k)    // inf,  Cover_OPT_UPPERBOUND
    {

        priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap;
        vector<int>coverage(social_graph->n, 0);

        for (int i = 0; i < social_graph->n; i++)
        {
            pair<int, int>tep(make_pair(i, (int)hyperG[i].size()));  //i: node_id,  hyperG: FRset{ rr set_id...}
            heap.push(tep);
            coverage[i] = (int)hyperG[i].size();
        }

        int maxInd;

        long long influence = 0;
        long long numEdge = hyperGT.size();

        // check if an edge is removed
        vector<bool> edgeMark(numEdge, false);
        // check if an node is remained in the heap
        vector<bool> nodeMark(social_graph->n + 1, true);

        seedSet.clear();
        while ((int)seedSet.size()<k)
        {
            pair<int, int>ele = heap.top();
            heap.pop();
            if (ele.second > coverage[ele.first])
            {
                ele.second = coverage[ele.first];
                heap.push(ele);
                continue;
            }

            maxInd = ele.first;
            vector<int>e = hyperG[maxInd];
            influence += coverage[maxInd];
            seedSet.push_back(maxInd);
            nodeMark[maxInd] = false;

            for (unsigned int j = 0; j < e.size(); ++j){
                if (edgeMark[e[j]])continue;

                vector<int>nList = hyperGT[e[j]];
                for (unsigned int l = 0; l < nList.size(); ++l){
                    if (nodeMark[nList[l]])coverage[nList[l]]--;
                }
                edgeMark[e[j]] = true;
            }
        }
        double inf = 1.0*influence / hyperGT.size();
        return inf;
    }

    int check_RR_Cover(vector<int> seeds){
        int rrset_size = get_RR_sets_size();
        std:: vector<bool> vecBoolVst = std::vector<bool>(rrset_size);
        //
        for(int u: seeds){
            for(auto rr_id: hyperGVld[u]) // mark each rr set in FR_set(u)
                vecBoolVst[rr_id] = true;
        }
        //count the number of marked rr set
        return std::count(vecBoolVst.begin(), vecBoolVst.end(),true);
    }

    int check_RR_POI_Cover(vector<int> pois){
        int rrset_size = get_RR_POI_sets_size();
        std:: vector<bool> vecBoolVst = std::vector<bool>(rrset_size);
        //
        for(int p: pois){
            int p_idx = social_graph->poi_idx_table[p];
            for(auto rr_id: hyperG_POIVld[p_idx]) // mark each rr set in FR_set(u)
                vecBoolVst[rr_id] = true;
        }
        //count the number of marked rr set
        return std::count(vecBoolVst.begin(), vecBoolVst.end(),true);
    }



    double build_seedset_OPIMC(int k)
    {

        priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap;
        int n = social_graph->n;
        vector<int>coverage(n, 0);

        for (int i = 0; i < n; i++)
        {
            pair<int, int>tep(make_pair(i, (int)hyperG[i].size()));
            heap.push(tep);
            coverage[i] = (int)hyperG[i].size();
        }

        int maxInd;

        long long influence = 0;
        long long numEdge = hyperGT.size();

        // check if an edge is removed
        vector<bool> edgeMark(numEdge, false);
        // check if an node is remained in the heap
        vector<bool> nodeMark(n + 1, true);

        seedSet.clear();
        while ((int)seedSet.size()<k)
        {
            pair<int, int>ele = heap.top();
            heap.pop();
            if (ele.second > coverage[ele.first])
            {
                ele.second = coverage[ele.first];
                heap.push(ele);
                continue;
            }

            maxInd = ele.first;
            vector<int>e = hyperG[maxInd];
            influence += coverage[maxInd];
            seedSet.push_back(maxInd);
            nodeMark[maxInd] = false;

            for (unsigned int j = 0; j < e.size(); ++j){
                if (edgeMark[e[j]])continue;

                vector<int>nList = hyperGT[e[j]];
                for (unsigned int l = 0; l < nList.size(); ++l){
                    if (nodeMark[nList[l]])coverage[nList[l]]--;
                }
                edgeMark[e[j]] = true;
            }
        }
        return 1.0*influence / hyperGT.size();
    }


    double do_poiSelection(int k)
    {

        priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap;
        vector<int>coverage(social_graph->np, 0);

        int np = social_graph->np;
        //所有（兴趣点）待选单元按 coverage到的采样点数目，从大到小排序
        for (int i = 0; i < np; i++)
        {
            pair<int, int>tep(make_pair(i, (int)hyperG_POI[i].size()));  //key:i(idx), cover_elements(采样点id)
            heap.push(tep);
            int size0 = (int)hyperG_POI[i].size();
            coverage[i] = size0;
        }

        int maxInd;

        long long influence = 0;
        long long numEdge = hyperGT_POI.size();  //采样个数

        // check if an edge is removed  对应采样点 （usr)
        vector<bool> edgeMark(numEdge, false);
        // check if an poi is remained in the heap 对应 兴趣点 （poi）
        vector<bool> poiMark(np + 1, true);

        poiSet.clear();

        while ((int)poiSet.size()<k)
        {
            pair<int, int>ele = heap.top();
            heap.pop();
            if (ele.second > coverage[ele.first])
            {
                ele.second = coverage[ele.first]; //更新
                heap.push(ele); //reinsert
                continue;  //next iteration
            }

            maxInd = ele.first;
            vector<int>e = hyperG_POI[maxInd];  //maxInd(poi_id)  , e={rr_id...}
            influence += coverage[maxInd];

            int p_id = social_graph->idx_poi_table[maxInd];
            poiSet.push_back(p_id);

            poiMark[maxInd] = false;

            for (unsigned int j = 0; j < e.size(); ++j){
                if (edgeMark[e[j]])continue;   //e[j]为rr set id , 查看是否已处理过

                vector<int>pList = hyperGT_POI[e[j]];
                for (unsigned int l = 0; l < pList.size(); ++l){
                    int p_idx = social_graph->poi_idx_table[pList[l]];
                    if (poiMark[p_idx])coverage[p_idx]--;
                }
                edgeMark[e[j]] = true;
            }
        }
        return 1.0*influence / hyperGT_POI.size();
    }


    double do_poiSelection_UB(int k)  // add Cover_UpperBound evaluation
    {

        priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap;
        vector<int>coverage(social_graph->getNp(), 0);
        //vector<bool> seed_flag(np,false);

        //所有（兴趣点）待选单元按 coverage到的采样点数目，从大到小排序
        for (int i = 0; i < social_graph->getNp(); i++)
        {
            pair<int, int>tep(make_pair(i, (int)hyperG_POI[i].size()));  //key:i(idx), cover_elements(采样点id)
            heap.push(tep);
            int size0 = (int)hyperG_POI[i].size();
            coverage[i] = size0;
        }

        int maxInd;

        long long influence = 0;
        influence_upper = -1;
        long long numEdge = hyperGT_POI.size();  //采样个数

        // check if an edge is removed  对应采样点 （usr)
        vector<bool> edgeMark(numEdge, false);
        // check if an poi is remained in the heap 对应 兴趣点 （poi）
        vector<bool> poiMark(social_graph->np + 1, true);

        poiSet.clear();



        while ((int)poiSet.size()<k)
        {
            pair<int, int>ele = heap.top();
            heap.pop();
            if (ele.second > coverage[ele.first])   // first:id, second:marginal gain
            {
                ele.second = coverage[ele.first]; //有变动，需更新
                heap.push(ele); //reinsert
                continue;  //next iteration
            }

            maxInd = ele.first;
            vector<int>e = hyperG_POI[maxInd];  //maxInd(poi_id)  , e={rr_id...}
            influence += coverage[maxInd];

            int p_id = social_graph->idx_poi_table[maxInd];
            poiSet.push_back(p_id);


            poiMark[maxInd] = false;

            bool cover_updated = false;
            //update the marginal gain for each poi under current POI set
            for (unsigned int j = 0; j < e.size(); ++j){
                if (edgeMark[e[j]])continue;   //e[j]为rr set id , 查看是否已处理过

                vector<int>pList = hyperGT_POI[e[j]];
                for (unsigned int l = 0; l < pList.size(); ++l){
                    int p_idx = social_graph->poi_idx_table[pList[l]];
                    if (poiMark[p_idx]){
                        coverage[p_idx]--;
                        cover_updated = true;
                    }
                }
                edgeMark[e[j]] = true;
            }

            if(true){   // marginal gain ranking has changed
                //update cover upper
                priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap2;
                int np = social_graph->np;
                for (int i = 0; i < np; i++)
                {
                    if(poiMark[i]==false) continue;
                    //把marginal gain 有变化的 元祖信息重新加入 优先队列中
                    pair<int, int>tep2(make_pair(i, (int)coverage[i]));  //key:i(idx), cover_elements(采样点id)
                    heap2.push(tep2);
                }

                long long maxMC =0;
                for(int i=0;i<k;i++){
                    pair<int, int>ele = heap2.top();
                    heap2.pop();
                    int add = ele.second;
                    maxMC+= ele.second;
                }
                if(influence_upper==-1){
                    influence_upper = influence+maxMC;
                }
                else if (influence_upper>(influence+maxMC)){
                    long long tmp = influence+maxMC;
                    influence_upper = influence+maxMC;
                }

            }



        }
        int _size= hyperGT_POI.size();
        //cout<<"-----size="<<_size<<endl;
        //cout<<"influence(cover) = "<<influence<<endl;
        //cout<<"return influence="<<influence<<"/"<<",hyperGT_POI.size="<<hyperGT_POI.size()<<endl;
        return 1.0*influence / hyperGT_POI.size();
    }


    double do_poiSelection_UB_OUT(int k)  // add Cover_UpperBound evaluation
    {
        outer_influence_only_upper = -1;
        priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap;
        vector<int>coverage(social_graph->getNp(), 0);
        //vector<bool> seed_flag(np,false);

        //所有（兴趣点）待选单元按 coverage到的采样点数目，从大到小排序
        for (int i = 0; i < social_graph->getNp(); i++)
        {
            pair<int, int>tep(make_pair(i, (int)hyperG_POI_OUT[i].size()));  //key:i(idx), cover_elements(采样点id)
            heap.push(tep);
            int size0 = (int)hyperG_POI_OUT[i].size();
            coverage[i] = size0;
        }

        int maxInd;

        long long influence = 0;
        influence_upper = -1;
        long long numEdge = hyperGT_POI_OUT.size();  //采样个数

        // check if an edge is removed  对应采样点 （usr)
        vector<bool> edgeMark(numEdge, false);
        // check if an poi is remained in the heap 对应 兴趣点 （poi）
        vector<bool> poiMark(social_graph->np + 1, true);

        poiSet.clear();



        while ((int)poiSet.size()<k)
        {
            pair<int, int>ele = heap.top();
            heap.pop();
            if (ele.second > coverage[ele.first])   // first:id, second:marginal gain
            {
                ele.second = coverage[ele.first]; //有变动，需更新
                heap.push(ele); //reinsert
                continue;  //next iteration
            }

            maxInd = ele.first;
            vector<int>e = hyperG_POI_OUT[maxInd];  //maxInd(poi_id)  , e={rr_id...}
            influence += coverage[maxInd];

            int p_id = social_graph->idx_poi_table[maxInd];
            poiSet.push_back(p_id);


            poiMark[maxInd] = false;

            bool cover_updated = false;
            //update the marginal gain for each poi under current POI set
            for (unsigned int j = 0; j < e.size(); ++j){
                if (edgeMark[e[j]])continue;   //e[j]为rr set id , 查看是否已处理过

                vector<int>pList = hyperGT_POI_OUT[e[j]];
                for (unsigned int l = 0; l < pList.size(); ++l){
                    int p_idx = social_graph->poi_idx_table[pList[l]];
                    if (poiMark[p_idx]){
                        coverage[p_idx]--;
                        cover_updated = true;
                    }
                }
                edgeMark[e[j]] = true;
            }

            if(true){   // marginal gain ranking has changed
                //update cover upper
                priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap2;
                int np = social_graph->np;
                for (int i = 0; i < np; i++)
                {
                    if(poiMark[i]==false) continue;
                    //把marginal gain 有变化的 元祖信息重新加入 优先队列中
                    pair<int, int>tep2(make_pair(i, (int)coverage[i]));  //key:i(idx), cover_elements(采样点id)
                    heap2.push(tep2);
                }

                long long maxMC =0;
                for(int i=0;i<k;i++){
                    pair<int, int>ele = heap2.top();
                    heap2.pop();
                    int add = ele.second;
                    maxMC+= ele.second;
                }
                if(outer_influence_only_upper==-1){
                    outer_influence_only_upper = influence+maxMC;
                }
                else if (outer_influence_only_upper>(influence+maxMC)){
                    long long tmp = influence+maxMC;
                    outer_influence_only_upper = influence+maxMC;
                }

            }



        }
        int _size= hyperGT_POI.size();
        //cout<<"-----size="<<_size<<endl;
        //cout<<"influence(cover) = "<<influence<<endl;
        //cout<<"return influence="<<influence<<"/"<<",hyperGT_POI.size="<<hyperGT_POI.size()<<endl;
        return 1.0*influence / hyperGT_POI.size();
    }



    double compute_usersetLocalInf0 (set<int> total_users){
        double _tmp_inf = 0;
        set<int> hop_neighbors;
        map<int,vector<double>> hop_EdgeWeightMap;
        for(int usr_id: total_users){
            _tmp_inf +=1;  //首先计入该seed_usr_id的影响力
            int seed_usr_idx = social_graph-> usr_idx_table[usr_id];
            for(int i=0;i<social_graph->out_neighbors[seed_usr_idx].size();i++){  //找出所有oneHop_usr，
                int oneHop_usr = social_graph->out_neighbors[seed_usr_idx][i];
                double activation_weight = social_graph->probAF[seed_usr_idx][i];
                hop_neighbors.insert(oneHop_usr);
                //将每条可达到oneHop_usr的（正向）边上的weight加入oneHop_EdgeWeightMap
                hop_EdgeWeightMap[oneHop_usr].push_back(activation_weight);
            }
        }
        for(int oneHop_usr: hop_neighbors){
            double bigPai = 1.0;
            for(double weight: hop_EdgeWeightMap[oneHop_usr]){
                bigPai = bigPai * (1-weight);
            }
            double influence_ratio_usr = 1.0 - bigPai;
            _tmp_inf += influence_ratio_usr;
            //pai_seed2v[i][oneHop_usr] = influence_ratio_usr;
        }
        return  _tmp_inf;
    }

#define LOCAL_INF1
    double compute_usersetLocalInf (set<int> total_users){
        double _tmp_inf = 0;
        set<int> hop_neighbors;
        map<int,vector<double>> hop_EdgeWeightMap;
        for(int usr_id: total_users){
            _tmp_inf +=1;  //首先计入该seed_usr_id的影响力
            //int seed_usr_idx = social_graph-> usr_idx_table[usr_id];
            for(int i=0;i<social_graph->out_neighbors[usr_id].size();i++){  //找出所有oneHop_usr，
                int oneHop_usr = social_graph->out_neighbors[usr_id][i];
                double activation_weight = social_graph->probAF[usr_id][i];
                hop_neighbors.insert(oneHop_usr);
                //将每条可达到oneHop_usr的（正向）边上的weight加入oneHop_EdgeWeightMap
                hop_EdgeWeightMap[oneHop_usr].push_back(activation_weight);
            }
        }
        for(int oneHop_usr: hop_neighbors){
            double bigPai = 1.0;
            for(double weight: hop_EdgeWeightMap[oneHop_usr]){
                bigPai = bigPai * (1-weight);
            }
            double influence_ratio_usr = 1.0 - bigPai;
            _tmp_inf += influence_ratio_usr;
            //pai_seed2v[i][oneHop_usr] = influence_ratio_usr;
        }
        return  _tmp_inf;
    }



    double do_poiSelection_UB_ByOneHop(int k)  // the algorithm using One Hop based approach
    {
        local_influence_upper = -1;
        double total_inf = 0;
        priority_queue<pair<int, double>, vector<pair<int, double>>, CompareBySecond>heap;
        vector<double>coverage(social_graph->getNp(), 0);
        //vector<bool> seed_flag(np,false);

        //所有（兴趣点） (的potential user集合) 待选单元按 其 One hop neighbor influence region 中的影响力进行 从大到小排序
        int poi_size = social_graph->getNp();
        map<int,map<int,double>> pai_seed2v;  // <poi_id, <user_id, ratio>


        for (int i = 0; i < poi_size; i++)
        {
            //提取该poi的potential user的one hop neighbor
            int poi_id = social_graph-> idx_poi_table[i];
            vector<int> poi_potential_user = social_graph->br_usr_table[poi_id];
            set<int> oneHop_neighbors; //poi_potential_user :找到与所有seeds集中用户 1 度相邻的user
            map<int,vector<double>> oneHop_EdgeWeightMap; //map: one hop user, <与该one hop user一度相邻的所有 seed node...>
            double local_influence_poi = 0.0;
            for(int seed_usr_id: poi_potential_user){
                local_influence_poi +=1;  //首先计入该seed_usr_id的影响力
                //int seed_usr_idx = social_graph-> usr_idx_table[seed_usr_id];
                for(int i=0;i<social_graph->out_neighbors[seed_usr_id].size();i++){  //找出所有oneHop_usr，
                    int oneHop_usr = social_graph -> out_neighbors[seed_usr_id][i];
                    double activation_weight = social_graph -> probAF[seed_usr_id][i];
                    oneHop_neighbors.insert(oneHop_usr);
                    //将每条可达到oneHop_usr的（正向）边上的weight加入oneHop_EdgeWeightMap
                    oneHop_EdgeWeightMap[oneHop_usr].push_back(activation_weight);
                }
            }
            //根据oneHop_neighbors、oneHop_EdgeWeightMap 计算 该poi的 local influence:   pai(PU(p), u) = 1 - bigPai(1-edge_weight)

            for(int oneHop_usr: oneHop_neighbors){
                double bigPai = 1.0;
                for(double weight: oneHop_EdgeWeightMap[oneHop_usr]){
                    bigPai = bigPai * (1-weight);
                }
                double influence_ratio_usr = 1.0 - bigPai;
                local_influence_poi += influence_ratio_usr;
                pai_seed2v[i][oneHop_usr] = influence_ratio_usr;
            }

            pair<int, double>tep(make_pair(i, local_influence_poi));  //key:i(idx), cover_elements(采样点id)
            heap.push(tep);
            coverage[i] = local_influence_poi;
        }

        int maxInd;

        double influence = 0;
        influence_upper = -1;
        long long numEdge = hyperGT_POI.size();  //采样个数

        // check if an edge is removed  对应采样点 （usr)
        vector<bool> edgeMark(numEdge, false);
        // check if an poi is remained in the heap 对应 兴趣点 （poi）
        vector<bool> poiMark(social_graph->np + 1, true);
        //vector<bool> usrSelectedMark(social_graph->nc + 1, false);
        bool* usrSelectedMark = new bool[social_graph->n];
        for(int i=0;i<=social_graph->nc;i++){
            usrSelectedMark[i] = false;
        }
        vector<int> selected_usrs;

        poiSet.clear();


        while ((int)poiSet.size()<k)
        {
            pair<int, double>ele = heap.top();
            heap.pop();
#ifdef LOCAL_INF
            cout<<"ele.first="<<ele.first<<endl;
            int _poi = social_graph->idx_poi_table[ele.first];
            cout<<"对应poi"<<_poi<<endl;
#endif
            if (ele.second > coverage[ele.first])   // first:id, second:marginal gain
            {
#ifdef LOCAL_INF
                cout<<"ele.first="<<ele.first<<"对应poi"<<_poi<<"的marginal gain已失效！"<<endl;
#endif
                ele.second = coverage[ele.first]; //有变动，需更新
                heap.push(ele); //reinsert
                continue;  //next iteration
            }

            maxInd = ele.first;
            //vector<int>e = hyperG_POI[maxInd];  //maxInd(poi_id)  , e={rr_id...}, 首先得到该poi对应了哪些 one Hop neighbor
            //cout<<"maxInd:"<<maxInd<<endl;
            int p_id = social_graph->idx_poi_table[maxInd];
            double delta_inf = coverage[maxInd];
#ifdef LOCAL_INF
            cout<<"select 第"<<poiSet.size()+1<<"个 poi:o"<<p_id;
            cout<<",delta_inf= "<<delta_inf<<endl;

#endif
            vector<int> potential_user_poiOPT = social_graph->br_usr_table[p_id];  //<poi_id,  {its potential user list}>
            influence += delta_inf;
            total_inf += delta_inf;

            poiSet.push_back(p_id);
            poiMark[maxInd] = false;  //标记该poi已选！
            for(int selected_u: potential_user_poiOPT){
                usrSelectedMark[selected_u] = true;
                selected_usrs.push_back(selected_u);
            }

            bool cover_updated = false;
            /*--(更新)--------根据当前poi selection情况，更新其他（有coverage的） poi 的 marginal gain--------*/
            /*-----------------1. 获得与当前POI seed集，local influence region 有重合的 poi----------------*/
            set<int> need_updated_pois;
            for (unsigned int j = 0; j < potential_user_poiOPT.size(); ++j){
                int u_id = potential_user_poiOPT[j];
                for (unsigned int p_th = 0; p_th < social_graph->reachable_poi_table[u_id].size(); ++p_th){
                    //首先更新与当前opt_poi有potential users重合的 那部分poi 的 marginal gain
                    int p_id = social_graph->reachable_poi_table[u_id][p_th];
                    //int p_idx = social_graph->poi_idx_table[p_id];
                    need_updated_pois.insert(p_id);
                }
            }
#ifdef LOCAL_INF
            printf("1. 获得与当前POI seed集，local influence region 有重合的 poi 成功！\n");
#endif
            /*-----------------2. 跟新coverage_poiSet中各个poi的marginal gain----------------*/
            //2.1 先对当前 seed的inf 进行计算： Gain (Sp)
            set<int> current_users;
            for(int u: selected_usrs) current_users.insert(u);
            double _inf0 = compute_usersetLocalInf(current_users);
            for(int poi_id: need_updated_pois){   // marginal gain ranking has changed
                set<int> total_users = current_users;
                int p_idx = social_graph->poi_idx_table[poi_id];
                if(poiMark[p_idx]==false) continue; //如果是已选为seed的poi,则跳过！
                /*--------------evaluate marginal gain for each users in the  local influence region of this poi-------------*/
                //2.2 对Gain (Sp + ps) 进行计算 (此时考虑纳入 ps)
                for(int u_id: social_graph->br_usr_table[poi_id])
                    total_users.insert(u_id);
                double _inf_combine = compute_usersetLocalInf(total_users);
                double _inf_delta = _inf_combine - _inf0;
                coverage[p_idx] = _inf_delta;
                //将更新marginal gain后的poi 信息条目加入优先队列表
                pair<int, double>tep2(make_pair(p_idx, (double)coverage[p_idx]));  //key:i(idx), cover_elements(采样点id)
                heap.push(tep2);
            }

#ifdef LOCAL_INF
            printf("2. 重新计算各个poi的marginal gain...成功！\n");
#endif

            if(true){   // marginal gain ranking has changed
                //update cover upper
                priority_queue<pair<int, double>, vector<pair<int, double>>, CompareBySecond>heap2;
                int np = social_graph->np;
                for (int i = 0; i < np; i++)
                {
                    if(poiMark[i]==false) continue;
                    //把marginal gain 有变化的 元祖信息重新加入 优先队列中
                    pair<int, int>tep2(make_pair(i, coverage[i]));  //key:i(idx), cover_elements(采样点id)
                    heap2.push(tep2);
                }
#ifdef LOCAL_INF
                printf("3. 根据更新后的情况，重新排序poi 成功！\n");
#endif

                double maxMC =0;
                for(int i=0;i<k;i++){
                    pair<int, double>ele = heap2.top();
                    heap2.pop();
                    int add = ele.second;
                    maxMC+= ele.second;
                }
#ifdef LOCAL_INF
                printf("计算累计的maxMC 成功！\n");
#endif
                if(local_influence_upper ==-1){
                    local_influence_upper = influence+maxMC;
                }
                else if (local_influence_upper>(influence+maxMC)){
                    long long tmp = influence+maxMC;
                    local_influence_upper = influence+maxMC;
                }
#ifdef LOCAL_INF
                printf("更新 local_influence_upper 成功！\n");

#endif


            }


        }
        delete []usrSelectedMark;

        return total_inf;
    }

    double do_poiSelection_Heuristic(int k)  // the algorithm using One Hop based approach
    {
        local_influence_upper = -1;
        double total_inf = 0;
        priority_queue<pair<int, double>, vector<pair<int, double>>, CompareBySecond>heap;
        vector<double>coverage(social_graph->getNp(), 0);
        //vector<bool> seed_flag(np,false);

        //所有（兴趣点） (的potential user集合) 待选单元按 其 One hop neighbor influence region 中的影响力进行 从大到小排序
        int poi_size = social_graph->getNp();
        map<int,map<int,double>> pai_seed2v;  // <poi_id, <user_id, ratio>


        for (int i = 0; i < poi_size; i++)
        {
            //提取该poi的potential user的one hop neighbor
            int poi_id = social_graph-> idx_poi_table[i];

            vector<int> poi_potential_user = social_graph->br_usr_table[poi_id];
            set<int> oneHop_neighbors; //poi_potential_user :找到与seed user集合中 1 度相邻的user
            map<int,vector<double>> oneHop_EdgeWeightMap; //map: one hop user, <与该one hop user一度相邻的所有 seed node...>
            double local_influence_poi = 0.0;
#ifdef LOCAL_INF
            printf("o%d has %d 个 potential users!\n",poi_id,poi_potential_user.size());
#endif
            for(int seed_usr_id: poi_potential_user){
                ////cout<<"seed_usr_id: poi_potential_user： u"<<seed_usr_id<<endl;
                local_influence_poi +=1;  //首先计入该seed_usr_id的影响力
                //int seed_usr_idx = social_graph-> usr_idx_table[seed_usr_id];
                ////printf("u%d has %d个one hop neighbor!,具体为：\n",seed_usr_id,social_graph->out_neighbors[seed_usr_id].size());
                if(!social_graph->out_neighbors[seed_usr_id].size()>0) continue;
                for(int i=0;i<social_graph->out_neighbors[seed_usr_id].size();i++){  //找出所有oneHop_usr，
                    int oneHop_usr = social_graph -> out_neighbors[seed_usr_id][i];
                    ////printf("oneHop neighbor:u%d\n",oneHop_usr);
                    double activation_weight = social_graph -> probAF[seed_usr_id][i];
                    oneHop_neighbors.insert(oneHop_usr);
                    //将每条可达到oneHop_usr的（正向）边上的weight加入oneHop_EdgeWeightMap
                    oneHop_EdgeWeightMap[oneHop_usr].push_back(activation_weight);
                }
            }
            //根据oneHop_neighbors、oneHop_EdgeWeightMap 计算 该poi的 local influence:   pai(PU(p), u) = 1 - bigPai(1-edge_weight)

            for(int oneHop_usr: oneHop_neighbors){
                double bigPai = 1.0;
                for(double weight: oneHop_EdgeWeightMap[oneHop_usr]){
                    bigPai = bigPai * (1-weight);
                }
                double influence_ratio_usr = 1.0 - bigPai;
                local_influence_poi += influence_ratio_usr;
                pai_seed2v[i][oneHop_usr] = influence_ratio_usr;
            }

            pair<int, double>tep(make_pair(i, local_influence_poi));  //key:i(idx), cover_elements(采样点id)
            heap.push(tep);
            //cout<<"heap push p_th"<<i<<",local_influence_poi="<<local_influence_poi<<endl;
            double _tmp_inf = local_influence_poi;
            coverage[i] = local_influence_poi;
#ifdef LOCAL_INF
            ////cout<<"_tmp_inf="<<_tmp_inf<<endl;
#endif
        }
        //getchar();

        int maxInd;

        double influence = 0;
        influence_upper = -1;

        // check if an poi is remained in the heap 对应 兴趣点 （poi）
        vector<bool> poiMark(social_graph->np + 1, true);
        //vector<bool> usrSelectedMark(social_graph->nc + 1, false);
        int _size = social_graph->n;
        bool* usrSelectedMark = new bool[social_graph->n];
        for(int i=0;i<=social_graph->n;i++){
            usrSelectedMark[i] = false;
        }
        vector<int> selected_usrs;

        poiSet.clear();


        while ((int)poiSet.size()<k)
        {
            pair<int, double>ele = heap.top();
            heap.pop();
#ifdef LOCAL_INF
            cout<<"ele.first="<<ele.first<<endl;
            int _poi = social_graph->idx_poi_table[ele.first];
            cout<<"对应poi"<<_poi<<endl;
#endif
            if (ele.second > coverage[ele.first])   // first:id, second:marginal gain
            {
#ifdef LOCAL_INF
                cout<<"ele.first="<<ele.first<<"对应poi"<<_poi<<"的marginal gain已失效！"<<endl;
#endif
                ele.second = coverage[ele.first]; //有变动，需更新
                heap.push(ele); //reinsert
                continue;  //next iteration
            }

            maxInd = ele.first;
            if(poiMark[maxInd]==false)
                continue;
            //vector<int>e = hyperG_POI[maxInd];  //maxInd(poi_id)  , e={rr_id...}, 首先得到该poi对应了哪些 one Hop neighbor
            cout<<"maxInd:"<<maxInd<<endl;
            int p_id = social_graph->idx_poi_table[maxInd];

            vector<int> potential_user_poiOPT = social_graph->br_usr_table[p_id];  //<poi_id,  {its potential user list}>
            float delta_inf = coverage[maxInd];
            influence += delta_inf;
            total_inf += delta_inf;
            poiSet.push_back(p_id);
            poiMark[maxInd] = false;  //标记该poi已选！
#ifdef LOCAL_INF
            cout<<"select 第"<<poiSet.size()<<"个 poi:o"<<p_id<<endl;
            cout<<"delta_inf="<<delta_inf<<endl;
            cout<<"total_inf="<<total_inf<<endl;
#endif


            for(int selected_u: potential_user_poiOPT){
                usrSelectedMark[selected_u] = true;
                selected_usrs.push_back(selected_u);
            }

            bool cover_updated = false;
            /*--(更新)--------根据当前poi selection情况，更新其他（有coverage的） poi 的 marginal gain--------*/
            /*-----------------1. 获得与当前POI seed集，local influence region 有重合的 poi----------------*/
            set<int> need_updated_pois;
            for (unsigned int j = 0; j < potential_user_poiOPT.size(); ++j){
                int u_id = potential_user_poiOPT[j];
                for (unsigned int p_th = 0; p_th < social_graph->reachable_poi_table[u_id].size(); ++p_th){
                    //首先更新与当前opt_poi有potential users重合的 那部分poi 的 marginal gain
                    int p_id = social_graph->reachable_poi_table[u_id][p_th];
                    //int p_idx = social_graph->poi_idx_table[p_id];
                    need_updated_pois.insert(p_id);
                }
            }
#ifdef LOCAL_INF
            printf("1. 获得与当前POI seed集，local influence region 有重合的 poi 成功！\n");
#endif
            /*-----------------2. 跟新coverage_poiSet中各个poi的marginal gain----------------*/
            //2.1 先对当前 seed的inf 进行计算： Gain (Sp)
            set<int> current_users;
            for(int u: selected_usrs) current_users.insert(u);
            double _inf0 = compute_usersetLocalInf(current_users);
            for(int poi_id: need_updated_pois){   // marginal gain ranking has changed
                set<int> total_users = current_users;
                int p_idx = social_graph->poi_idx_table[poi_id];
                if(poiMark[p_idx]==false) continue; //如果是已选为seed的poi,则跳过！
                /*--------------evaluate marginal gain for each users in the  local influence region of this poi-------------*/

                //2.2 对Gain (Sp + ps) 进行计算 (此时考虑纳入 ps)
                for(int u_id: social_graph->br_usr_table[poi_id])
                    total_users.insert(u_id);
                double _inf_combine = compute_usersetLocalInf(total_users);
                double _inf_delta = _inf_combine - _inf0;
                coverage[p_idx] = _inf_delta;
                //将更新marginal gain后的poi 信息条目加入优先队列表
                pair<int, double>tep2(make_pair(p_idx, (double)coverage[p_idx]));  //key:i(idx), cover_elements(采样点id)
                heap.push(tep2);
            }

#ifdef LOCAL_INF
            printf("2. 重新计算各个poi的marginal gain...成功！\n");
#endif

            if(true){   // marginal gain ranking has changed
                //update cover upper
                priority_queue<pair<int, double>, vector<pair<int, double>>, CompareBySecond>heap2;
                int np = social_graph->np;
                for (int i = 0; i < np; i++)
                {
                    if(poiMark[i]==false) continue;
                    //把marginal gain 有变化的 元祖信息重新加入 优先队列中
                    pair<int, int>tep2(make_pair(i, coverage[i]));  //key:i(idx), cover_elements(采样点id)
                    heap2.push(tep2);
                }
#ifdef LOCAL_INF
                printf("3. 根据更新后的情况，重新排序poi 成功！\n");
#endif

                double maxMC =0;
                for(int i=0;i<k;i++){
                    pair<int, double>ele = heap2.top();
                    heap2.pop();
                    int add = ele.second;
                    maxMC+= ele.second;
                }
#ifdef LOCAL_INF
                printf("计算累计的maxMC 成功！\n");
#endif
                if(local_influence_upper ==-1){
                    local_influence_upper = influence+maxMC;
                }
                else if (local_influence_upper>(influence+maxMC)){
                    long long tmp = influence+maxMC;
                    local_influence_upper = influence+maxMC;
                }
#ifdef LOCAL_INF
                printf("更新 local_influence_upper 成功！\n");
#endif


            }


        }

        delete []usrSelectedMark;

        return total_inf;
    }

    double do_poiSelection_Heuristic_upper(int k)  // the algorithm using One Hop based approach
    {
        local_influence_upper = -1;
        double total_inf = 0;
        priority_queue<pair<int, double>, vector<pair<int, double>>, CompareBySecond>heap;
        vector<double>coverage(social_graph->getNp_upper(), 0);
        //vector<bool> seed_flag(np,false);

        //所有（兴趣点） (的potential user集合) 待选单元按 其 One hop neighbor influence region 中的影响力进行 从大到小排序
        int poi_size = social_graph->getNp_upper();
        map<int,map<int,double>> pai_seed2v;  // <poi_id, <user_id, ratio>


        for (int i = 0; i < poi_size; i++)
        {
            //提取该poi的potential user的one hop neighbor
            int poi_id = social_graph-> idx_poi_table_upper[i];
            vector<int> poi_potential_user = social_graph->br_usr_table_upper[poi_id];
            set<int> oneHop_neighbors; //poi_potential_user :找到与所有seeds集中用户 1 度相邻的user
            map<int,vector<double>> oneHop_EdgeWeightMap; //map: one hop user, <与该one hop user一度相邻的所有 seed node...>
            double local_influence_poi = 0.0;
            for(int seed_usr_id: poi_potential_user){
                local_influence_poi +=1;  //首先计入该seed_usr_id的影响力
                //int seed_usr_idx = social_graph-> usr_idx_table_upper[seed_usr_id];
                for(int i=0;i<social_graph->out_neighbors[seed_usr_id].size();i++){  //找出所有oneHop_usr，
                    int oneHop_usr = social_graph -> out_neighbors[seed_usr_id][i];
                    double activation_weight = social_graph -> probAF[seed_usr_id][i];
                    oneHop_neighbors.insert(oneHop_usr);
                    //将每条可达到oneHop_usr的（正向）边上的weight加入oneHop_EdgeWeightMap
                    oneHop_EdgeWeightMap[oneHop_usr].push_back(activation_weight);
                }
            }
            //根据oneHop_neighbors、oneHop_EdgeWeightMap 计算 该poi的 local influence:   pai(PU(p), u) = 1 - bigPai(1-edge_weight)

            for(int oneHop_usr: oneHop_neighbors){
                double bigPai = 1.0;
                for(double weight: oneHop_EdgeWeightMap[oneHop_usr]){
                    bigPai = bigPai * (1-weight);
                }
                double influence_ratio_usr = 1.0 - bigPai;
                local_influence_poi += influence_ratio_usr;
                pai_seed2v[i][oneHop_usr] = influence_ratio_usr;
            }

            pair<int, double>tep(make_pair(i, local_influence_poi));  //key:i(idx), cover_elements(采样点id)
            heap.push(tep);
            coverage[i] = local_influence_poi;
        }

        int maxInd;

        double influence = 0;
        influence_upper = -1;

        // check if an poi is remained in the heap 对应 兴趣点 （poi）
        vector<bool> poiMark(social_graph->np_upper + 1, true);
        //vector<bool> usrSelectedMark(social_graph->nc + 1, false);
        int _size = social_graph->n;
        bool* usrSelectedMark = new bool[social_graph->n];
        for(int i=0;i<=social_graph->n;i++){
            usrSelectedMark[i] = false;
        }
        vector<int> selected_usrs;

        poiSet2.clear();


        while ((int)poiSet2.size()<k)
        {
            pair<int, double>ele = heap.top();
            heap.pop();
#ifdef LOCAL_INF
            cout<<"ele.first="<<ele.first<<endl;
            int _poi = social_graph->idx_poi_table[ele.first];
            cout<<"对应poi"<<_poi<<endl;
#endif
            if (ele.second > coverage[ele.first])   // first:id, second:marginal gain
            {
#ifdef LOCAL_INF
                cout<<"ele.first="<<ele.first<<"对应poi"<<_poi<<"的marginal gain已失效！"<<endl;
#endif
                ele.second = coverage[ele.first]; //有变动，需更新
                heap.push(ele); //reinsert
                continue;  //next iteration
            }

            maxInd = ele.first;
            if(poiMark[maxInd] == false)
                continue;
            //cout<<"maxInd:"<<maxInd<<endl;
            int p_id = social_graph->idx_poi_table_upper[maxInd];
#ifdef LOCAL_INF
            cout<<"select 第"<<poiSet2.size()+1<<"个 poi:o"<<p_id<<endl;

#endif
            vector<int> potential_user_poiOPT = social_graph->br_usr_table_upper[p_id];  //<poi_id,  {its potential user list}>
            influence += coverage[maxInd];
            total_inf += coverage[maxInd];

            poiSet2.push_back(p_id);
            poiMark[maxInd] = false;  //标记该poi已选！
            for(int selected_u: potential_user_poiOPT){
                usrSelectedMark[selected_u] = true;
                selected_usrs.push_back(selected_u);
            }

            bool cover_updated = false;
            /*--(更新)--------根据当前poi selection情况，更新其他（有coverage的） poi 的 marginal gain--------*/
            /*-----------------1. 获得与当前POI seed集，local influence region 有重合的 poi----------------*/
            set<int> need_updated_pois;
            for (unsigned int j = 0; j < potential_user_poiOPT.size(); ++j){
                int u_id = potential_user_poiOPT[j];
                for (unsigned int p_th = 0; p_th < social_graph->reachable_poi_table_upper[u_id].size(); ++p_th){
                    //首先更新与当前opt_poi有potential users重合的 那部分poi 的 marginal gain
                    int p_id = social_graph->reachable_poi_table_upper[u_id][p_th];
                    //int p_idx = social_graph->poi_idx_table[p_id];
                    need_updated_pois.insert(p_id);
                }
            }
#ifdef LOCAL_INF
            printf("1. 获得与当前POI seed集，local influence region 有重合的 poi 成功！\n");
#endif
            /*-----------------2. 跟新coverage_poiSet中各个poi的marginal gain----------------*/
            //2.1 先对当前 seed的inf 进行计算： Gain (Sp)
            set<int> current_users;
            for(int u: selected_usrs) current_users.insert(u);
            double _inf0 = compute_usersetLocalInf(current_users);
            for(int poi_id: need_updated_pois){   // marginal gain ranking has changed
                set<int> total_users = current_users;
                int p_idx = social_graph->poi_idx_table_upper[poi_id];
                if(poiMark[p_idx]==false) continue; //如果是已选为seed的poi,则跳过！
                /*--------------evaluate marginal gain for each users in the  local influence region of this poi-------------*/

                //2.2 对Gain (Sp + ps) 进行计算 (此时考虑纳入 ps)
                for(int u_id: social_graph->br_usr_table_upper[poi_id])
                    total_users.insert(u_id);
                double _inf_combine = compute_usersetLocalInf(total_users);
                double _inf_delta = _inf_combine - _inf0;
                coverage[p_idx] = _inf_delta;
                //将更新marginal gain后的poi 信息条目加入优先队列表
                pair<int, double>tep2(make_pair(p_idx, (double)coverage[p_idx]));  //key:i(idx), cover_elements(采样点id)
                heap.push(tep2);
            }

#ifdef LOCAL_INF
            printf("2. 重新计算各个poi的marginal gain...成功！\n");
#endif

            if(true){   // marginal gain ranking has changed
                //update cover upper
                priority_queue<pair<int, double>, vector<pair<int, double>>, CompareBySecond>heap2;
                int np_upper = social_graph->np_upper;
                for (int i = 0; i < np_upper; i++)
                {
                    if(poiMark[i]==false) continue;
                    //把marginal gain 有变化的 元祖信息重新加入 优先队列中
                    pair<int, int>tep2(make_pair(i, coverage[i]));  //key:i(idx), cover_elements(采样点id)
                    heap2.push(tep2);
                }
#ifdef LOCAL_INF
                printf("3. 根据更新后的情况，重新排序poi 成功！\n");
#endif

                double maxMC =0;
                for(int i=0;i<k;i++){
                    pair<int, double>ele = heap2.top();
                    heap2.pop();
                    int add = ele.second;
                    maxMC+= ele.second;
                }
#ifdef LOCAL_INF
                printf("计算累计的maxMC 成功！\n");
#endif
                if(local_influence_upper ==-1){
                    local_influence_upper = influence+maxMC;
                }
                else if (local_influence_upper>(influence+maxMC)){
                    long long tmp = influence+maxMC;
                    local_influence_upper = influence+maxMC;
                }
#ifdef LOCAL_INF
                printf("更新 local_influence_upper 成功！\n");
#endif


            }


        }
        delete []usrSelectedMark;

        return total_inf;
    }


    Inf_Pair do_poiSelection_UB_ByHybrid(int k)  // add Cover_UpperBound evaluation
    {
        Inf_Pair inf_pair;
        inf_pair.outer_inf=0.0; inf_pair.local_inf=0.0;
        overall_influence_upper =-1;
        priority_queue<pair<int, double>, vector<pair<int, double>>, CompareBySecond>heap;
        vector<double>coverage(social_graph->getNp(), 0.0);
        vector<double> coverage_local(social_graph->getNp(), 0.0);
        vector<double> coverage_outer(social_graph->getNp(), 0.0);
        vector<double> estimated_outer(social_graph->getNp(), 0.0);
        double ratio = social_graph->getNc()/hyperGT_POI_OUT.size();
        //vector<bool> seed_flag(np,false);

        //所有（兴趣点）待选单元按 coverage到的采样点数目，从大到小排序
        for (int i = 0; i < social_graph->getNp(); i++)
        {
            int poi_id = social_graph-> idx_poi_table[i];
            set<int> poi_potential_user;
            for(int u_id: social_graph->br_usr_table[poi_id]) poi_potential_user.insert(u_id);
            double local_inf = compute_usersetLocalInf(poi_potential_user);
            int cover_count = hyperG_POI_OUT[i].size();
            double outer_inf = ratio * cover_count; ////注意这里是外部影响力评估值， outer_influence_estimated = nc/RRset_num * coverage(R1,Ps)
            double overall_inf = local_inf + outer_inf;

            pair<int, double>tep(make_pair(i, overall_inf));  //key:i(idx), cover_elements(采样点id)
            heap.push(tep);
            coverage[i] = overall_inf;
            coverage_local[i] = local_inf;
            coverage_outer[i] = cover_count;
            estimated_outer[i] = outer_inf;
            //cout<<"initialize p"<<i<<",local_inf="<<local_inf<<", outer_inf="<<outer_inf<<endl;
        }

        int maxInd;

        double influence = 0;
        double influence_local =0;
        double influence_outer =0;
        influence_upper = -1;
        long long numEdge = hyperGT_POI_OUT.size();  //采样个数

        // check if an edge is removed  对应采样点 （usr)
        vector<bool> edgeMark(numEdge, false);
        // check if an poi is remained in the heap 对应 兴趣点 （poi）
        vector<bool> poiMark(social_graph->np + 1, true);
        //vector<bool> usrSelectedMark(social_graph->nc + 1, false);
        bool* usrSelectedMark = new bool[social_graph->n];
        for(int i=0;i<=social_graph->nc;i++){
            usrSelectedMark[i] = false;
        }
        vector<int> selected_usrs;

        poiSet.clear();


        while ((int)poiSet.size()<k)
        {
            pair<int, double>ele = heap.top();
            heap.pop();
            if (ele.second > coverage[ele.first])   // first:id, second:marginal gain
            {
                ele.second = coverage[ele.first]; //有变动，需更新
                heap.push(ele); //reinsert
                continue;  //next iteration
            }

            maxInd = ele.first;
            int p_id = social_graph->idx_poi_table[maxInd];  //当前最优poi, 加入POS seeds
            vector<int>e = hyperG_POI_OUT[maxInd];  //maxInd(poi_id)  , e={rr_id...}
            influence += coverage[maxInd];
            influence_local += coverage_local[maxInd];
            influence_outer += estimated_outer[maxInd];
            inf_pair.local_inf += coverage_local[maxInd];
            inf_pair.outer_inf += estimated_outer[maxInd];
            //cout<<"select poi"<<p_id;
            //cout<<",its coverage_local["<<maxInd<<"]="<<coverage_local[maxInd]<<endl;


            vector<int> potential_user_poiOPT = social_graph->br_usr_table[p_id];  //<poi_id,  {its potential user list}>

            poiSet.push_back(p_id);
            poiMark[maxInd] = false;
            //vector<int> potential_user_poiOPT = social_graph->br_usr_table[p_id];  //<poi_id,  {its potential user list}>
            for(int selected_u: potential_user_poiOPT){
                usrSelectedMark[selected_u] = true;
                selected_usrs.push_back(selected_u);
            }
            if((int)poiSet.size()==k){
                double inf_upper = local_influence_upper;
                //cout<<"local inf="<<inf_pair.local_inf<<", outer_inf="<<inf_pair.outer_inf<<endl;
                //getchar();
                return inf_pair;
            }

            bool cover_updated = false;
            ////更新 outer 影响力增益情况：update the marginal gain of outer influence for each poi under current POI seed set
            set<int> need_updated_pois_outer;
            for (unsigned int j = 0; j < e.size(); ++j){
                if (edgeMark[e[j]])continue;   //e[j]为rr set id , 查看是否已处理过

                vector<int>pList = hyperGT_POI_OUT[e[j]];
                for (unsigned int l = 0; l < pList.size(); ++l){
                    int p_idx = social_graph->poi_idx_table[pList[l]];
                    if (poiMark[p_idx]){
                        need_updated_pois_outer.insert(p_idx);
                        estimated_outer[p_idx] = estimated_outer[p_idx]-ratio*1;  //joke, 重要！ 注意有ratio加持
                        coverage_outer[p_idx]--;
                        coverage[p_idx]= coverage[p_idx] - ratio*1 ;  //joke, 重要！ 注意有ratio加持
                        cover_updated = true;
                    }
                }
                edgeMark[e[j]] = true;
            }

            ////更新 local 影响力增益情况：update the marginal gain of local influence for each poi under current POI seed set
            /*-----------------1. 获得与当前POI seed集，local influence region 有重合的 poi----------------*/
            set<int> need_updated_pois;  set<int> need_updated_pois_local;
            for (unsigned int j = 0; j < potential_user_poiOPT.size(); ++j){
                int u_id = potential_user_poiOPT[j];
                for (unsigned int p_th = 0; p_th < social_graph->reachable_poi_table[u_id].size(); ++p_th){
                    int p_id = social_graph->reachable_poi_table[u_id][p_th];
                    int p_idx = social_graph->poi_idx_table[p_id];
                    need_updated_pois.insert(p_id);
                    need_updated_pois_local.insert(p_idx);
                }
            }
            /*-----------------2. 跟新coverage_poiSet中各个poi的marginal gain的 local部分----------------*/
            //2.1 先对当前 seed的inf 进行计算： Gain (Sp)
            set<int> total_users;
            for(int u: selected_usrs) total_users.insert(u);
            double _inf0 = compute_usersetLocalInf(total_users);
            for(int poi_id: need_updated_pois){   // marginal gain ranking has changed

                int p_idx = social_graph->poi_idx_table[poi_id];
                if(poiMark[p_idx]==false) continue; //如果是已选为seed的poi,则跳过！
                /*--------------evaluate marginal gain for each users in the  local influence region of this poi-------------*/

                //2.2 对Gain (Sp + ps) 进行计算 (此时考虑纳入 ps)
                for(int u_id: social_graph->br_usr_table[poi_id]) total_users.insert(u_id);
                double _inf_combine = compute_usersetLocalInf(total_users);
                double _inf_delta =  _inf_combine - _inf0;
                coverage_local[p_idx] = _inf_delta;
                coverage[p_idx] = estimated_outer[p_idx]+ coverage_local[p_idx];  //注意：marginal gain: 外部影响力估计值.gain + 局部影响力真实值.gain
                //将更新marginal gain后的poi 信息条目加入优先队列表
            }


            if(true){   // marginal gain ranking has changed
                //update cover upper
                priority_queue<pair<int, double>, vector<pair<int, double>>, CompareBySecond> heap_overall;
                //priority_queue<pair<int, double>, vector<pair<int, double>>, CompareBySecond>heap_outer;
                int np = social_graph->np;
                for (int i = 0; i < np; i++)
                {
                    if(poiMark[i]==false) continue;
                    //把marginal gain 有变化的 元祖信息重新加入 优先队列中
                    pair<int, double>tep2(make_pair(i, coverage[i]));   //key:i(idx), cover_elements(采样点id)
                    heap_overall.push(tep2);

                }
                double maxMC = 0;
                double relatedMC_local = 0;
                double relatedMC_outer = 0;
                for(int i=0;i<k;i++){
                    pair<int, double>ele = heap_overall.top();
                    heap_overall.pop();
                    int p_idx = ele.first;
                    double add = ele.second;
                    maxMC+= ele.second;
                    relatedMC_local += coverage_local[p_idx];
                    relatedMC_outer += coverage_outer[p_idx];
                }
                if(overall_influence_upper==-1){
                    overall_influence_upper = influence+maxMC;
                    overall_influence_upperDetails.local_inf = relatedMC_local;
                    overall_influence_upperDetails.outer_inf = relatedMC_outer;
                }
                else if (overall_influence_upper>(influence+maxMC)){
                    double tmp = influence + maxMC;
                    overall_influence_upper = influence+maxMC;    //joke
                    overall_influence_upperDetails.local_inf = relatedMC_local;
                    overall_influence_upperDetails.outer_inf = relatedMC_outer;
                }

            }


        }

        int _size= hyperGT_POI.size();
        delete [] usrSelectedMark;
        //cout<<"-----size="<<_size<<endl;
        //cout<<"influence(cover) = "<<influence<<endl;
        //cout<<"return influence="<<influence<<"/"<<",hyperGT_POI.size="<<hyperGT_POI.size()<<endl;
        overall_influence_upper =-1;
        return inf_pair;
    }



    int check_RR_POI_Cover_OUT(vector<int> pois){
        int rrset_size = get_RR_POI_sets_outer_size();
        std:: vector<bool> vecBoolVst = std::vector<bool>(rrset_size);
        //
        for(int p: pois){
            int p_idx = social_graph->poi_idx_table[p];
            for(auto rr_id: hyperG_POIVld_OUT[p_idx]) // mark each rr set in FR_set(u)
                vecBoolVst[rr_id] = true;
        }
        //count the number of marked rr set
        return std::count(vecBoolVst.begin(), vecBoolVst.end(),true);
    }



    double do_poiSelection_UB_Upper(int k)  // add Cover_UpperBound evaluation
    {

        priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap;
        vector<int>coverage(social_graph->getNp_upper(), 0);
        //vector<bool> seed_flag(np,false);

        //所有（兴趣点）待选单元按 coverage到的采样点数目，从大到小排序
        for (int i = 0; i < social_graph->getNp_upper(); i++)
        {
            int _size = (int)hyperG_Upper_POI[i].size();
            pair<int, int>tep(make_pair(i, _size));  //key:i(idx), cover_elements(采样点id)
            heap.push(tep);
            int size0 = (int)hyperG_Upper_POI[i].size();
            coverage[i] = size0;
        }

        int maxInd;

        long long influence = 0;
        influence_upper2 = -1;
        long long numEdge = hyperGT_Upper_POI.size();  //采样个数

        // check if an edge is removed  对应采样点 （usr)
        vector<bool> edgeMark(numEdge, false);
        // check if an poi is remained in the heap 对应 兴趣点 （poi）
        vector<bool> poiMark(social_graph->np_upper + 1, true);

        poiSet2.clear();
        influence_upper2=-1;

        while ((int)poiSet2.size()<k)
        {
            pair<int, int>ele = heap.top();
            heap.pop();
            if (ele.second > coverage[ele.first])   // first:id, second:marginal gain
            {
                ele.second = coverage[ele.first]; //更新
                heap.push(ele); //reinsert
                continue;  //next iteration
            }

            maxInd = ele.first;
            vector<int>e = hyperG_Upper_POI[maxInd];  //maxInd(poi_id)  , e={rr_id...}
            influence += coverage[maxInd];

            int p_id = social_graph->idx_poi_table_upper[maxInd];
            poiSet2.push_back(p_id);


            poiMark[maxInd] = false;

            bool cover_updated = false;
            //update the marginal gain for each poi under current POI set
            for (unsigned int j = 0; j < e.size(); ++j){
                if (edgeMark[e[j]])continue;   //e[j]为rr set id , 查看是否已处理过

                vector<int>pList = hyperGT_Upper_POI[e[j]];
                for (unsigned int l = 0; l < pList.size(); ++l){
                    int p_idx = social_graph->poi_idx_table_upper[pList[l]];
                    if (poiMark[p_idx]){
                        coverage[p_idx]--;
                        cover_updated = true;
                    }
                }
                edgeMark[e[j]] = true;
            }

            if(true){   // marginal gain ranking has changed
                //update cover upper
                priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap2;
                int np_upper = social_graph->np_upper;
                for (int i = 0; i < np_upper; i++)
                {
                    if(poiMark[i]==false) continue;
                    pair<int, int>tep2(make_pair(i, (int)coverage[i]));  //key:i(idx), cover_elements(采样点id)
                    heap2.push(tep2);
                }

                long long maxMC =0;
                for(int i=0;i<k;i++){
                    pair<int, int>ele = heap2.top();
                    heap2.pop();
                    maxMC+= ele.second;
                }
                if(influence_upper2==-1){
                    influence_upper2 = influence+maxMC;
                }
                else if (influence_upper2>(influence+maxMC)){
                    influence_upper2 = influence+maxMC;
                }
            }



        }
        return 1.0*influence / hyperGT_Upper_POI.size();
    }




    double MC_Evaluation(int r)
    {
        cout<<"MC_Evaluation for "<<r<<"rounds..."<<endl;
        double inf_MC = 0;
        for(int i=0;i<r;i++){
            int n_visit_mark = 0;
            q.clear();
            for(int u: seedSet){
                q.push_back(u);
                visit_mark[n_visit_mark++] = u;
                visit[u] = true;
            }

            //ASSERT(n_visit_mark < n);

            while (!q.empty())
            {

                int expand = q.front();
                q.pop_front();
                if (true)  //influModel == IC
                {
                    int i = expand;
                    for (int j = 0; j < (int)social_graph->in_neighbors[i].size(); j++)
                    {
                        //int u=expand;
                        int v = social_graph->out_neighbors[i][j];
                        double randDouble =
                                dsfmt_gv_genrand_uint32_range(1001)/1000.0;
                        double prob = social_graph->probAF[i][j];
                        //cout<<"prob("<<i<<","<<j<<")="<<prob<<endl;
                        //if(i%100==0)
                        //getchar();
                        if (randDouble > prob)
                            continue;
                        if (visit[v])
                            continue;
                        if (!visit[v])
                        {
                            visit_mark[n_visit_mark++] = v;
                            visit[v] = true;
                        }
                        q.push_back(v);
                        //ASSERT((int)hyperGT.size() > hyperiiid);
                        //hyperGT[hyperiiid].push_back(v);
                    }
                }
            }
            // fu wei
            for (int i = 0; i < n_visit_mark; i++)
                visit[visit_mark[i]] = false;

            inf_MC += n_visit_mark;
        }
        //average
        inf_MC = (double) inf_MC/r;
        cout<<"inf_MC = "<<inf_MC<<endl;

        return inf_MC;
    }


    double MC_Evaluation_GivenSeeds(int r, set<int> seedSet)
    {
        cout<<"MC_Evaluation_GivenSeeds for "<<r<<"rounds..."<<endl;
        double inf_MC = 0;
        for(int i=0;i<r;i++){
            int n_visit_mark = 0;
            q.clear();
            for(int u : seedSet){
                q.push_back(u);
                visit_mark[n_visit_mark++] = u;
                visit[u] = true;
            }

            //ASSERT(n_visit_mark < n);

            while (!q.empty())
            {

                int expand = q.front();
                q.pop_front();
                if (true)  //influModel == IC
                {
                    int i = expand;
                    for (int j = 0; j < (int)social_graph->in_neighbors[i].size(); j++)
                    {
                        //int u=expand;
                        int v = social_graph->out_neighbors[i][j];
                        double randDouble =
                                dsfmt_gv_genrand_uint32_range(1001)/1000.0;
                        double prob = social_graph->probAF[i][j];
                        //cout<<"prob("<<i<<","<<j<<")="<<prob<<endl;
                        //if(i%100==0)
                        //getchar();
                        if (randDouble > prob)
                            continue;
                        if (visit[v])
                            continue;
                        if (!visit[v])
                        {
                            visit_mark[n_visit_mark++] = v;
                            visit[v] = true;
                            q.push_back(v);
                        }

                        //ASSERT((int)hyperGT.size() > hyperiiid);
                        //hyperGT[hyperiiid].push_back(v);
                    }
                }
            }
            // fu wei
            for (int i = 0; i < n_visit_mark; i++)
                visit[visit_mark[i]] = false;

            inf_MC += n_visit_mark;
        }
        //average
        inf_MC = (double) inf_MC/r;
        cout<<"inf_MC = "<<inf_MC<<endl;

        return inf_MC;
    }


    double MC_Evaluation_POI(int r)
    {
#ifdef MC_LOG
        cout<<"MC_Evaluation_POI for "<<r<<"rounds..."<<endl;
#endif
        double inf_MC = 0;
        //set<int> _pois;   //75,57,155,154,153
        //_pois.insert(57);_pois.insert(75); _pois.insert(155);_pois.insert(154);_pois.insert(153);
        set<int> real_seeds;
        for(int p: poiSet){  //for(int p: _pois){//
            for(int u: social_graph->br_usr_table[p]){
                real_seeds.insert(u);
            }
        }
        cout<<"Ps的潜在用户集合大小："<<real_seeds.size()<<"，具体为："<<endl;
        for(int u:real_seeds){
            cout<<"u"<<u<<",";
        }
        cout<<endl;
        for(int i= 1;i<=r;i++){
            int n_visit_mark = 0;
            deque<int> q2;
            q2.clear();
            //不同的poi可能有相同的potential users, 这里要注意！!!!!!
            for(int u: real_seeds){
                visit_mark[n_visit_mark++] = u;
                visit[u] = true;
                q2.push_back(u);
            }

            //ASSERT(n_visit_mark < n);

            while (!q2.empty())
            {

                int expand = q2.front();
                q2.pop_front();
                if (true)  //influModel == IC
                {
                    int i = expand;

                    for (int j = 0; j <  social_graph->out_neighbors[i].size(); j++)
                    {
                        //int u=expand;
                        int v = social_graph->out_neighbors[i][j];
                        //cout<<i<<"的第"<<j<<"个neighbor为"<<v<<endl;
                        if(v<0) continue;

                        if (!visit[v])
                        {
                            double randDouble = ((float)(rand() % 1001))/(float)1000;
                            //dsfmt_gv_genrand_uint32_range(1001)/1000.0;
                            double prob = social_graph->probAF[i][j];
                            //printf("probAF[%d][%d]=%f\n",i,v,prob);
                            if(prob >= randDouble){
                                visit_mark[n_visit_mark++] = v;
                                //cout<<"set visit["<<v<<"] as true...";
                                visit[v] = true;
                                //cout<<"successfully!"<<endl;
                                q2.push_back(v);
                                //cout<<"push "<<v<<endl;
                            }
                        }

                        //ASSERT((int)hyperGT.size() > hyperiiid);
                        //hyperGT[hyperiiid].push_back(v);
                    }
                }
            }
            // fu wei
            //cout<<"round "<<i<<" fu wei...";
            for (int i = 0; i < n_visit_mark; i++)
                visit[visit_mark[i]] = false;
            //cout<<"COMPELETE"<<endl;

            inf_MC += n_visit_mark;
#ifdef MC_LOG
            if(i%2000==0){
                cout<<"round "<<i<<", inf="<<(float)inf_MC/(float)i<<endl;
            }
#endif

        }
        //average
        inf_MC = (float) inf_MC/(float)r;
#ifdef MC_LOG
        cout<<"the number of attracted users = "<< real_seeds.size()<<"具体为："<<endl;
        printRealSeeds(real_seeds);
#endif
        cout<<"the overall influence spread= "<<inf_MC<<endl;
        //cout<<"inf_MC(POI) = "<<inf_MC<<endl;

        return inf_MC;
    }




    void printRealSeeds(set<int> _seeds){
        cout<<"{ " ;
        for(int u: _seeds){
            cout<<"u"<<u<<",";
        }
        cout<<"};"<<endl;
    }



    void pre_sampling_build_hyper_graph(int R, const Argument & arg)
    {

        int prevSize = hyperGTVld.size();
        while ((int)hyperGTVld.size() <= R)
            hyperGTVld.push_back( vector<int>() );
        //cout<<"1"<<endl;


        vector<int> random_number;
        for (int i = 0; i < R; i++)
        {
            int random = dsfmt_gv_genrand_uint32_range(social_graph->n);  //gai
            random_number.push_back(random);
        }
        //cout<<"social_graph->n="<<social_graph->n<<endl;
        //cout<<"2"<<endl;


        //trying BFS start from same node

        for (int i = 0; i < R; i++)
        {
            //cout<<"build rr set"<<i<<"for u"<<random_number[i]<<endl;
            BuildHypergraphNode(random_number[i], i, R2);
            //cout<<"build succeed!"<<i<<endl;


        }
        //cout<<"3"<<endl;

        int totAddedElement = 0;
        for (int i = prevSize; i < R; i++)
        {
            for (int t : hyperGTVld[i])
            {
                hyperGVld[t].push_back(i);
                //hyperG.addElement(t, i);
                totAddedElement++;
            }
        }
        //cout<<"4"<<endl;
    }


    void pre_sampling_FIS(int R) {
        cout << "pre_sampling Monte carlo simulation sketch for " << R << " times!" << endl;
        int social_nodeSize = social_graph->in_neighbors.size();
        vector<MCSketch> MCSketchSet;
        //先对MC simulation 初始化
        for (int round = 0; round < R; round++) {
            MCSketch instance;
            for (int u = 0; u < social_nodeSize; u++) {
                vector<int> u_neighbor;
                instance.push_back(u_neighbor);

            }
            MCSketchSet.push_back(instance);
        }
        //每轮次下，对社交网络影响力图上的每条边都进行采样，总共R个轮次
        for (int round = 0; round < R; round++) {
            for (int u = 0; u < social_nodeSize; u++) {
                int i = u;
                for (int j = 0; j < (int) social_graph->out_neighbors[i].size(); j++) {
                    int v = social_graph->out_neighbors[i][j];
                    double randDouble =
                            dsfmt_gv_genrand_uint32_range(1001) / 1000.0;
                    double prob = social_graph->probAF[i][j];

                    if (randDouble > prob) {
                        continue;
                    } else {
                        MCSketchSet[round][u].push_back(v);
                    }
                }
            }
        }

        //对每个user在R 个MC采样实例下的（双跳）节点激活情况进行统计
        for (int u = 0; u < social_graph->out_neighbors.size(); u++) {
            vector<int>  activated_Entry;
            for (int r = 0; r < R; r++) {
                //1. 进行one hop 邻居统计
                int n_visit_mark = 0; q.clear();
                int _size = MCSketchSet[r][u].size();
                for(int k=0;k<_size;k++){
                    int uk= MCSketchSet[r][u][k];
                    if (visit[uk])
                        continue;
                    q.push_back(uk);
                    visit_mark[n_visit_mark++] = uk;
                    visit[u] = true;
                    int key = social_nodeSize*r + uk;  //RR_key = <instance_id, u_id>
                    activated_Entry.push_back(key);
                }
                //2. 对后续hop neighbor进行统计
                while (!q.empty())
                {
                    int expand = q.front();
                    q.pop_front();
                    int i = expand;
                    for (int j = 0; j < MCSketchSet[r][i].size(); j++) {
                        int uj = MCSketchSet[r][i][j];
                        if (visit[uj])
                            continue;
                        if (!visit[uj]) {
                            visit_mark[n_visit_mark++] = uj;
                            visit[uj] = true;
#ifdef Global_Inf
                            q.push_back(uj);
#endif
                            int key = social_nodeSize*r + uj;  //RR_key = <instance_id, u_id>
                            activated_Entry.push_back(key);
                        }
                    }
                }
                for (int i = 0; i < n_visit_mark; i++){
                    int actived_u = visit_mark[i];
                    visit[actived_u] = false;
                }


            }//50次采样完毕！
            cout<<"u"<<u<<" 's activated node id:"<<endl;
            for(int id:activated_Entry){
                int instance_id = id /social_nodeSize;
                int u_id = id % social_nodeSize;
                //cout<<id<<",";
                cout<<"<u"<<u_id<<",i"<<instance_id<<">,";
            }
            cout<<endl;
            cout<<"u"<<u<<" 's activation node size = "<<activated_Entry.size()<<endl;
            getchar();
        }

    }

};



#endif




