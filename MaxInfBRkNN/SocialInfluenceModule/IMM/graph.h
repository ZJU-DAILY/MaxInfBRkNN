#define HEAD_INFO
//#define HEAD_TRACE
#define DISCRETE


//#include "../../SFMT/dSFMT/dSFMT.c"

#include "../imbasic.h"
using namespace std;
typedef double (*pf)(int, int);
typedef vector<vector<int>> MCSketch;


class SocialGraph {
public:
    // this is the data structure storing Social Network Graph
    int n, m, k;
    vector<int> inDeg;
    vector<vector<int>> in_neighbors;   //  in neighbor 邻居
    vector<vector<double>> probAR;  //反向激活概率


    vector<vector<int>> out_neighbors;   //  out neighbor 邻居
    vector<vector<double>> probAF;  //正向激活概率

    vector<bool> hasnode;

    /*-------for MaxInfBRGST----------------*/
    //this is for constraint graph
    int nc; int ns;
    int np; //有影响力（反近邻用户）的POI个数
    vector<bool> reachable;
    set<int> reachable_usrs;
    vector<int> reachable_list;
    map<int,int> usr_idx_table;
    unordered_map<int,int> idx_usr_table;

    set<int> inf_pois;
    vector<int> poi_list;
    map<int, int> idx_poi_table; // idx(in poi_list)  poi_id,
    map<int, int> poi_idx_table; // poi_id,  idx(in poi_list)

    set<int> reachable_source;  //u_idreachable_poi_table
    map<int, vector<int>> reachable_poi_table;   // usr_id, <poi_id...>
    map<int, vector<int>> br_usr_table;   // poi_id, <usr_id...>
    //----------------------------------------------------------------

    //this is for upper constraint graph
    int nc_upper; int ns_upper;
    int np_upper; //有影响力（反近邻用户）的POI个数
    vector<bool> reachable_upper;
    set<int> reachable_usrs_upper;
    vector<int> reachable_list_upper;
    map<int,int> usr_idx_table_upper;
    map<int,int> idx_usr_table_upper;

    set<int> inf_pois_upper;
    vector<int> poi_list_upper;
    map<int, int> idx_poi_table_upper; // idx(in poi_list)  poi_id,
    map<int, int> poi_idx_table_upper; // poi_id,  idx(in poi_list)

    set<int> reachable_source_upper;  //u_idreachable_poi_table
    map<int, vector<int>> reachable_poi_table_upper;   // usr_id, <poi_id...>
    map<int, vector<int>> br_usr_table_upper;   // poi_id, <usr_id...>
    //----------------------------------------------------------------




   InfluModel influModel;
   //jins
   SocialGraph(int vn){
        this -> n = vn;
        //m = en;

        for(int i=0;i< this->n;i++)  //为图结构初始化内存空间
        {
            in_neighbors.push_back(vector<int>());    //（反向）邻接表
            out_neighbors.push_back(vector<int>());  //正向 邻接表
            hasnode.push_back(false);     //

            reachable.push_back(false);   // 受限可达图，顶点标记向量
            reachable_upper.push_back(false);   // 受限可达图，顶点标记向量

            probAR.push_back(vector<double>()); //反向 概率表
            probAF.push_back(vector<double>());//正向 概率表
            inDeg.push_back(0);
        }
    }

   ~SocialGraph(){
       cout<<"~begin SocialGraph()"<<endl;
   }

    void add_edge(int a, int b, double p)
    {
        probAR[b].push_back(p);
        in_neighbors[b].push_back(a);
        inDeg[b]++;

        probAF[a].push_back(p);
        out_neighbors[a].push_back(b);

        hasnode[a] = true;
        hasnode[b] = true;
    }
    // 双色体反近邻usr_id,   poi_id
    void setRearchable(int u_id, int p_id){
       reachable[u_id] = true;
       reachable_poi_table[u_id].push_back(p_id);
       br_usr_table[p_id].push_back(u_id);  // poi 对应的反近邻用户
       reachable_usrs.insert(u_id);
       reachable_source.insert(u_id);
       inf_pois.insert(p_id);
   }

    void setRearchable_Upper(int u_id, int p_id){
        reachable_upper[u_id] = true;
        reachable_poi_table_upper[u_id].push_back(p_id);
        br_usr_table_upper[p_id].push_back(u_id);  // poi 对应的反近邻用户
        reachable_usrs_upper.insert(u_id);
        reachable_source_upper.insert(u_id);
        inf_pois_upper.insert(p_id);
    }

   void rearchInfo_clear(){
        for(int i = 0; i<reachable.size();i++){
           reachable[i]= false;
        }
        reachable_poi_table.clear();
        br_usr_table.clear(); // poi 对应的反近邻用户
        reachable_usrs.clear();
        reachable_source.clear();
        inf_pois.clear();




   }

    void rearchInfoUpper_clear(){
        for(int i = 0; i<reachable_upper.size();i++){
            reachable_upper[i]= false;
        }
        reachable_poi_table_upper.clear();
        br_usr_table_upper.clear();
        reachable_usrs_upper.clear();
        reachable_source_upper.clear();
        inf_pois_upper.clear();

    }



    int getNc(){
       nc = reachable_usrs.size();
       return nc;
   }
   int  getNs(){
       ns = reachable_source.size();
       return ns;
   }

    int getNp(){
        np = inf_pois.size();
        return np;
    }


    int getNc_upper(){
        nc_upper = reachable_usrs_upper.size();
        return nc_upper;
    }
    int  getNs_upper(){
        ns_upper = reachable_source_upper.size();
        return ns_upper;
    }

    int getNp_upper(){
        np_upper = inf_pois_upper.size();
        return np_upper;
    }


    //从双色体反近邻用户出发，对社交网络图进行遍历，直到访问不到新的节点
    void traFromSource(){  //obtain all the reachable user in the constraint graph
       cout<<"begin traFromSource ("<<reachable_source.size()<<"users)"<<endl;
       deque<int> q;
       q.clear();
       reachable_usrs.clear();
       reachable_list.clear();
       usr_idx_table.clear();
       idx_usr_table.clear();
       poi_list.clear();
       idx_poi_table.clear(); //下一步对应do_poiselection
       poi_idx_table.clear();

       //先将 source user (反近邻用户)加入 q
       //set<int> visited;
       for(int u: reachable_source){
           q.push_back(u);
           reachable_usrs.insert(u);
       }
       for(int i= 0; i<reachable.size();i++){
           reachable[i]= false;
       }


       //开始BFS遍历
       while(!q.empty()){
           int i = q.front();
           q.pop_front();
           for (int j = 0; j < (int)out_neighbors[i].size(); j++){
               int v = out_neighbors[i][j];
               if(reachable[v]== true)
                   continue;
               else{
                   reachable[v]= true;  //标记可达
                   reachable_usrs.insert(v); //将该用户加入受限可达图顶点集
                   q.push_back(v);
               }
           }
       }

       for(int u:reachable_usrs){
           int idx = reachable_list.size();
           reachable_list.push_back(u);
           usr_idx_table[u]= idx;
           idx_usr_table[idx] = u;
       }
       for(int p: inf_pois){
           int idx = poi_list.size();
           poi_list.push_back(p);
           idx_poi_table[idx] = p; //下一步对应do_poiselection
           poi_idx_table[p]= idx;
       }

       nc = reachable_list.size();
       np = inf_pois.size();
       cout<<"Gc有影响力兴趣点个数："<<np <<endl;
       cout<<"受限可达图顶点集个数："<< nc<<endl;
       cout<<"reachable_list中元素个数："<<reachable_list.size()<<endl;

   }

    void initialPOISelectionInfo(){  //用于heuristic算法
        poi_list.clear();
        idx_poi_table.clear();poi_idx_table.clear();
        for(int p: inf_pois){
            int idx = poi_list.size();
            poi_list.push_back(p);
            idx_poi_table[idx] = p; //下一步对应do_poiselection
            poi_idx_table[p]= idx;
        }

        np = inf_pois.size();
        cout<<"Gc有影响力兴趣点个数："<<np <<endl;

    }

    void initialPOISelectionInfo_Upper(){  //用于heuristic算法
        poi_list_upper.clear();
        idx_poi_table_upper.clear();poi_idx_table_upper.clear();
        for(int p: inf_pois_upper){
            int idx = poi_list_upper.size();
            poi_list_upper.push_back(p);
            idx_poi_table_upper[idx] = p; //下一步对应do_poiselection
            poi_idx_table_upper[p]= idx;
        }

        np_upper = inf_pois_upper.size();
        cout<<"Gc+有影响力兴趣点个数："<<np_upper <<endl;

    }


    void traFromSource_Upper(){  //obtain all the reachable user in the constraint graph
        cout<<"begin traFromSource_Upper ("<<reachable_source_upper.size()<<"users)"<<endl;
        deque<int> q;
        q.clear();
        reachable_usrs_upper.clear();
        reachable_list_upper.clear();
        usr_idx_table_upper.clear();
        idx_usr_table_upper.clear();
        poi_list_upper.clear();
        idx_poi_table_upper.clear(); //下一步对应do_poiselection
        poi_idx_table_upper.clear();


        for(int u: reachable_source_upper){
            q.push_back(u);
            reachable_usrs_upper.insert(u);
        }

        for(int i= 0; i<reachable_upper.size();i++){
            reachable_upper[i]= false;
        }

        //开始BFS遍历
        while(!q.empty()){
            int i = q.front();
            q.pop_front();
            for (int j = 0; j < (int)out_neighbors[i].size(); j++){
                int v = out_neighbors[i][j];
                if(reachable_upper[v]== true)
                    continue;
                else{
                    reachable_upper[v]= true;  //标记可达
                    reachable_usrs_upper.insert(v); //将该用户加入受限可达图顶点集
                    q.push_back(v);
                }
            }
        }
        for(int u:reachable_usrs_upper){
            int idx = reachable_list_upper.size();
            reachable_list_upper.push_back(u);
            usr_idx_table_upper[u]= idx;
            idx_usr_table_upper[idx] = u;
        }
        for(int p: inf_pois){
            int idx = poi_list_upper.size();
            poi_list_upper.push_back(p);
            idx_poi_table_upper[idx] = p; //下一步对应do_poiselection
            poi_idx_table_upper[p]= idx;
        }

        nc_upper = reachable_list_upper.size();
        np_upper = inf_pois_upper.size();
        cout<<"Gc_upper 有影响力兴趣点个数："<<np_upper <<endl;
        cout<<"受限Upper可达图顶点集个数："<< nc_upper<<endl;
        cout<<"reachable_list_upper中元素个数："<<reachable_list_upper.size()<<endl;

    }





    void transferGraph(AdjListSocial& wholeSocialLinkMap){
        cout<<"In function: transferGraph..."<<endl;
        //unordered_map<int,set<int>>::iterator it;
        //it = wholeSocialLinkMap.begin();
        int edgeNum =0;
        for(auto it= wholeSocialLinkMap.begin();it!=wholeSocialLinkMap.end();++it){
            int u2 =(int) (it->first);
            vector<int> infriendsOfU2 = it->second;
            int size =0; size= infriendsOfU2.size();
            if(size>0){
                double prob = 1.0/ infriendsOfU2.size();
                for(int u1: infriendsOfU2){
                    edgeNum++;
                    add_edge(u1,u2,prob);
                }
                it ++;
            }else{
                it ++;
                continue;
            }
        }

        cout<<"graph 导入完毕！"<<endl;
    }
};


/*class ConstraintGraph{
    int nc, mc;

    vector<vector<int>> cgT;   //
    vector<vector<double>> cprobAR;
    vector<vector<int>> cout_neighbors;   //
    vector<vector<double>> cprobAF;
};*/


