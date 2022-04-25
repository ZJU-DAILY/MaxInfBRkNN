//
// Created by jins on 10/31/19.
//

#ifndef MAXINFRKGSKQ_1_1_QUERY_BASIC_H
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <queue>
#include <sstream>
#include "querySet.h"
#include "loadInfo.h"
#include "utility.cc"
#include "functional.h"
#include "gim_tree.h"


#include "Tenindra/queue/BinaryMinHeap.h"
#include "timer.h"

using namespace std;

#define maxNum  1000000
#define Print_Debug 1;


float poiDistance[poi_num];



/*------------------------------GSK评分上下界求解函数----------------------------*/
//兴趣点与一组用户的GSK评分上界计算
double getMaxScoreOfNodeOnce_loc(int loc, int node_id,set<int> keyword_Set, set<int> checkin_usr_Set, float a){
    double min_Distance = getMinMaxDistanceOfNode(loc,node_id).min;
    //cout<<min_Distance<<endl;
    double max_relevance = 0;
    double max_influence ;
    double max_social_score;
    set<int> unionUserKey = GTree[node_id].userUKeySet;
    set<int> interKey;

    max_relevance = tfIDF(keyword_Set, unionUserKey);



    int  minPopular = 999999;
    //获取所有用户最佳兴趣点签到数目中的最小值
    for(int key:keyword_Set){

        if(node_term_MinPopular[node_id][key]!=NULL){
            minPopular = min(minPopular,node_term_MinPopular[node_id][key]);
        }
    }

    //cout<<"leaf"<<leaf_id<<", leaste POI count="<<minPopular<<endl;
    set<int> unionFriends = GTree[node_id].unionFriends;
    set<int> inter;
    set_intersection(checkin_usr_Set.begin(),checkin_usr_Set.end(),unionFriends.begin(),unionFriends.end(),inserter(inter,inter.begin()));
    //cout<<inter.size()<<endl;

    if(inter.size()==0)
        max_social_score = 0;
    else if(minPopular == 0)
        max_social_score = 1;
    else{
        double rate = inter.size()*1.0/minPopular;
        max_social_score = min(1.0,rate);
    }
    if(max_social_score == 0)
        max_influence =1.0;
    else
        max_influence =1.0 + pow(a,max_social_score);

    double max_Score = max_relevance * max_influence / (1.0+min_Distance);
    return  max_Score;
}



double getMaxScoreOfNodeOnce_poi(POI poi, int node_id,set<int> keyword_Set, set<int> checkin_usr_Set, float alpha, float a){
    double min_Distance = getMinMaxDistanceOfNode(poi.Nj,node_id).min;
    //cout<<min_Distance<<endl;
    double max_relevance = 0;
    double max_influence ;
    set<int> unionUserKey = GTree[node_id].userUKeySet;
    set<int> interKey;

    max_relevance = tfIDF(unionUserKey, keyword_Set);
    if(max_relevance==0)
        return 0.0;
    set<int> unionFriends = GTree[node_id].unionFriends;
    set<int> inter;
    inter = obtain_itersection_jins(unionFriends,checkin_usr_Set);

    int _size = inter.size();
    int min_popular = GTree[node_id].minSocialCount;
    //判断社交关联性
    double max_ratio = 1.0 ;  //这里可以优化
    if(inter.size()==0)
        max_ratio = 0;
    else if(min_popular == INT32_MAX)
        max_ratio = 1.0;
    else{
        double rate = inter.size()*1.0/min_popular;
        max_ratio = min(1.0,rate);
    }

    //cout<<"leaf"<<leaf_id<<", leaste POI count="<<minPopular<<endl;
    if(max_ratio == 0)
        max_influence =1.0;
    else
        max_influence =1.0 + pow(a,max_ratio);


    double max_social_textual_score = alpha*max_influence + (1.0-alpha)*max_relevance;

    double max_Score = max_social_textual_score / (1.0 + min_Distance);
    return  max_Score;
}

double getMaxScoreOfNodeOnce_poi_InterKey(POI poi, int node_id,set<int> inter_KeySet, set<int> checkin_usr_Set, float alpha, float a){

    double min_Distance = getMinMaxDistanceOfNode(poi.Nj,node_id).min;
    //cout<<min_Distance<<endl;
    double max_relevance = 0;
    double max_influence ;

    max_relevance = tfIDF_termSet(inter_KeySet);
    if(max_relevance==0)
        return 0.0;
    set<int> unionFriends;
#ifndef DiskAccess
    unionFriends = GTree[node_id].unionFriends;
#else
    unionFriends = GTree[node_id].unionFriends_load;
#endif
    set<int> inter;
    inter = obtain_itersection_jins(unionFriends,checkin_usr_Set);
    int _size = inter.size();
#ifndef DiskAccess
    int min_popular = GTree[node_id].minSocialCount;
#else
    int min_popular = GTree[node_id].minSocialCount_load;
#endif
    //判断社交关联性
    double max_ratio = 1.0 ;  //这里可以优化
    if(inter.size()==0)
        max_ratio = 0;
    else if(min_popular == INT32_MAX)
        max_ratio = 1.0;
    else{
        double rate = inter.size()*1.0/min_popular;
        max_ratio = min(1.0,rate);
    }

    //cout<<"leaf"<<leaf_id<<", leaste POI count="<<minPopular<<endl;
    if(max_ratio == 0)
        max_influence =1.0;
    else
        max_influence =1.0 + pow(a,max_ratio);


    double max_social_textual_score = alpha*max_influence + (1.0-alpha)*max_relevance;

    double max_Score = max_social_textual_score / (1.0 + min_Distance);
    return  max_Score;
}


double getMaxScoreOfLeafOnce(int loc, int leaf_id, double distance, set<int> keyword_Set, set<int> checkin_usr_Set, float a){

    double min_Distance = distance;
    //cout<<min_Distance<<endl;
    double max_relevance= 0;
    double max_influence ;
    double max_social_score;
    int  minPopular = 999999;
    set<int> unionUserKey = GTree[leaf_id].userUKeySet;
    //set<int> interKey;
    //set_intersection(keyword_Set.begin(),keyword_Set.end(),unionUserKey.begin(),unionUserKey.end(),inserter(interKey,interKey.begin()));
    double tfIDF_score = 0;

    max_relevance = tfIDF(unionUserKey,keyword_Set);
    //cout<< max_relevance<<endl;
    //getchar();
    if(max_relevance == 0)
        return  0.0;

    //cout<<"leaf"<<leaf_id<<", leaste POI count="<<minPopular<<endl;
    set<int> unionFriends = GTree[leaf_id].unionFriends;
    set<int> inter;
    set_intersection(checkin_usr_Set.begin(),checkin_usr_Set.end(),unionFriends.begin(),unionFriends.end(),inserter(inter,inter.begin()));
    //cout<<inter.size()<<endl;
    if(inter.size()==0)
        max_social_score = 0;
    else if(minPopular == 0)
        max_social_score = 1;
    else{
        double rate = inter.size()*1.0/minPopular;
        max_social_score = min(1.0,rate);
    }
    if(max_social_score == 0)
        max_influence =1;
    else
        max_influence =1 + pow(a,max_social_score);

    double max_Score = max_relevance * max_influence / (1+min_Distance);  //避免分母为零
    return  max_Score;
}

double getMaxScoreOfLeaf_InterKey(int loc, int leaf_id, double distance, set<int> inter_KeySet, set<int> checkin_usr_Set, float alpha, float a){

    double min_Distance = distance;
    //cout<<min_Distance<<endl;
    double max_relevance= 0;
    double max_influence =1.0 ;
    double max_social_score;
    int  minPopular = -1;
#ifndef DiskAccess
    minPopular = GTree[leaf_id].minSocialCount;
#else
    minPopular = GTree[leaf_id].minSocialCount_load;
#endif

    max_relevance = tfIDF_termSet(inter_KeySet);

    if(max_relevance == 0)
        return  0.0;

    //cout<<"leaf"<<leaf_id<<", leaste POI count="<<minPopular<<endl;
#ifndef DiskAccess
    set<int> unionFriends = GTree[leaf_id].unionFriends;
#else
    set<int> unionFriends = GTree[leaf_id].unionFriends_load;
#endif
    set<int> inter;
    inter = obtain_itersection_jins(unionFriends,checkin_usr_Set);



    int _size = inter.size();
    if(inter.size()==0)
        max_social_score = 0;
    else if(minPopular == INT32_MAX)
        max_social_score = 1;
    else{
        double rate = inter.size()*1.0/minPopular;
        max_social_score = min(1.0,rate);
    }
    if(max_social_score == 0)
        max_influence =1;
    else
        max_influence =1 + pow(a,max_social_score);


    double max_Score =
            (alpha * max_influence + (1-alpha)*max_relevance) / (1+min_Distance);  //避免分母为零
    return  max_Score;
}





/*-----------------------------求解用LCL与UCL的功能函数---------------------------*/

//usr节点继承上层节点的LCL结果，在此基础上进一步精简
LCLResult getUpperNode_LCL(int usr_node, int term_id, vector<int> related_poi_nodeSet, map<int,set<int>> leaf_poiCnt, int Qk, float a, double alpha, double upper_boud){
    priority_queue<LCLEntry> LCL;
    priority_queue<UCLEntry> UCL;

    int cnt=0;int cnt2 =0;
    int upper_cnt = 0;
    int count = 0;
    //extract relatd poi_leaves

    map<int,int> poiLeaf_include;
#ifdef TRACK
    cout<<"computing lcl for u node"<<usr_node<<" according to t"<<term_id<<endl;
#endif
    // access the leaves in related_poi_leafSetget to get the LCL of usr leaf
    for(int node:related_poi_nodeSet){  //首先该usr_node 的关键词情况是怎样的？

        double min_Score = getGSKScoreO2U_Lower(a, alpha, node, usr_node, term_id);
        //double max_Score = getGSKScoreO2U_Upper(a, node, usr_node);

        //double min_Score = 0.0;
        //if(poiLeaf_include[poi_leaf]==1) continue;
        //UCL.push(UCLEntry(node,max_Score, leaf_poiCnt[node].size())); //加入上界列表的排序，复杂度可能多log(n)
        if(cnt<Qk){
            LCL.push(LCLEntry(node, min_Score, leaf_poiCnt[node].size()));
            cnt += leaf_poiCnt[node].size();
        }
        else{
            //LCLEntry e =  LCL.top();
            while(min_Score>LCL.top().score && cnt > Qk){
                LCLEntry e =  LCL.top();
                LCL.pop();
                cnt -= e.count;
            }
            LCL.push(LCLEntry(node,min_Score, leaf_poiCnt[node].size()));
            cnt += leaf_poiCnt[node].size();
            //if(upper_boud < LCL.top().score && cnt > Qk){
            if(0 < LCL.top().score && cnt > Qk){
                LCLResult lclresult;
                lclresult.LCL = LCL;
                lclresult.UCL = UCL;
                lclresult.topscore = LCL.top().score;
                return lclresult;
            }
        }
    }
    LCLResult lclresult;
    lclresult.LCL = LCL;
    lclresult.UCL = UCL;
    lclresult.topscore = LCL.top().score;
    return lclresult;

}

//内存操作
LCLResult getUpperNode_termLCL(int usr_node, int term_id, map<int,vector<int>>& related_poi_nodeSet, int Qk, float a, double alpha, double upper_bound){
    priority_queue<LCLEntry> LCL;
    priority_queue<UCLEntry> UCL;

    int cnt=0;int cnt2 =0;
    int upper_cnt = 0;
    int count = 0;
    //extract relatd poi_leaves

    map<int,int> poiLeaf_include;
#ifdef  TRACK
    cout<<"computing lcl for u node"<<usr_node<<" according to t"<<term_id<<endl;
#endif

    // access the node in term related_poi_nodeSet to get the LCL of usr node
    for(int node:related_poi_nodeSet[term_id]){  //首先该usr_node 的关键词情况是怎样的？

        double min_Score = getGSKScoreO2U_Lower(a, alpha, node, usr_node, term_id);

        int _size = GTree[node].term_object_Map[term_id].size();
        if(cnt<Qk){
            LCL.push(LCLEntry(node, min_Score, _size));
            cnt += _size;
        }
        else{
            //LCLEntry e =  LCL.top();
            while(min_Score>LCL.top().score && cnt > Qk){
                LCLEntry e =  LCL.top();
                LCL.pop();
                cnt -= e.count;
            }
            LCL.push(LCLEntry(node,min_Score, _size));
            cnt += _size;
            float _tmp = LCL.top().score;
            if(upper_bound /1.0 < LCL.top().score && cnt > Qk){
                //if(0 < LCL.top().score && cnt > Qk){
                LCLResult lclresult;
                lclresult.LCL = LCL;
                lclresult.UCL = UCL;
                lclresult.topscore = LCL.top().score;
                return lclresult;
            }
        }
    }
    LCLResult lclresult;
    lclresult.LCL = LCL;
    lclresult.UCL = UCL;
    lclresult.topscore = LCL.top().score;
    return lclresult;

}

// 外存操作
LCLResult getUpperNode_termLCL_Disk_Small(int usr_node, int term_id, map<int,vector<int>>& related_poi_nodeSet, int Qk, float a, double alpha, double upper_bound){
    priority_queue<LCLEntry> LCL;
    priority_queue<UCLEntry> UCL;


    int cnt=0;int cnt2 =0;
    int upper_cnt = 0;
    int count = 0;
    //extract relatd poi_leaves

    map<int,int> poiLeaf_include;
#ifdef  LCL_LOG
    cout<<"computing lcl for u node"<<usr_node<<" according to t"<<term_id<<endl;
#endif

    // access the node in term related_poi_nodeSet to get the LCL of usr node
    for(int node:related_poi_nodeSet[term_id]){  //首先该usr_node 的关键词情况是怎样的？

        //cout<<"read onode "<<node<<endl;
        TreeNode oNode = getGIMTreeNodeData(node,OnlyO);

        //printONodeInfo(oNode,node);

        int _size = oNode.term_object_countMap[term_id];
        //cout<<"o node"<<node<<", object_count="<<_size<<endl;

        double min_Score = getGSKScoreO2U_Lower(a, alpha, node, usr_node, term_id);
        //cout<<"min_score="<<min_Score<<endl;


        if(cnt<Qk){
            LCL.push(LCLEntry(node, min_Score, _size));
            cnt += _size;
        }
        else{
            //LCLEntry e =  LCL.top();
            //cout<<"cnt="<<cnt<<endl;
            //cout<<"LCL.top().score="<<LCL.top().score<<endl;
            while(min_Score>LCL.top().score && cnt > Qk){
                LCLEntry e =  LCL.top();
                //cout<<"e.count="<<e.count<<endl;
                LCL.pop();
                //cout<<"LCL.pop()"<<endl;
                cnt -= e.count;
                //cout<<"cnt - e.count ="<<cnt<<endl;
            }
            LCL.push(LCLEntry(node,min_Score, _size));
            //cout<<"push LCL:"<<node<<",min_socre="<<min_Score<<",size="<<_size<<endl;

            cnt += _size;
            float _tmp = LCL.top().score;
            //cout<<"_tmp="<<_tmp<<endl;
            //cout<<"LCL.size="<<LCL.size()<<endl;
            //cout<<"upper_bound"<<upper_bound<<endl;
            if(upper_bound /1.0 < LCL.top().score && cnt > Qk){
                //cout<<"提前终止条件"<<endl;
                LCLResult lclresult;
                lclresult.LCL = LCL;
                lclresult.UCL = UCL;
                lclresult.topscore = LCL.top().score;
                return lclresult;
            }
        }
    }
    LCLResult lclresult;
    lclresult.LCL = LCL;
    lclresult.UCL = UCL;
    lclresult.topscore = LCL.top().score;
    return lclresult;

}


bool update_LCL( int usr_node, int onode, TreeNode& oNodeData, priority_queue<LCLEntry>& LCL, int term_id, int& cnt, int Qk,int a, double alpha, double upper_bound){
    int _size = oNodeData.term_object_countMap[term_id];
    //cout<<"o node"<<node<<", object_count="<<_size<<endl;
    double min_Score = getGSKScoreO2U_Lower(a, alpha, onode, usr_node, term_id);
    //cout<<"min_score="<<min_Score<<endl;

    if(cnt<Qk){
        LCL.push(LCLEntry(onode, min_Score, _size));
        //LCL.push(LCLEntry(node,min_Score, _size));
        cnt += _size;
    }
    else {
        //LCLEntry e =  LCL.top();
        //cout<<"cnt="<<cnt<<endl;
        //cout<<"LCL.top().score="<<LCL.top().score<<endl;
        while(LCL.size()>0 && min_Score>LCL.top().score && cnt > Qk){
            LCLEntry e =  LCL.top();
            //cout<<"e.count="<<e.count<<endl;
            LCL.pop();
            //cout<<"LCL.pop()"<<endl;
            cnt -= e.count;
            //cout<<"cnt - e.count ="<<cnt<<endl;
        }
        LCL.push(LCLEntry(onode,min_Score, _size));
        //cout<<"push LCL:"<<node<<",min_socre="<<min_Score<<",size="<<_size<<endl;

        cnt += _size;
        float _tmp = LCL.top().score;
        //cout<<"_tmp="<<_tmp<<endl;
        //cout<<"LCL.size="<<LCL.size()<<endl;
        //cout<<"upper_bound"<<upper_bound<<endl;
        if(upper_bound /1.0 < LCL.top().score && cnt > Qk){
            //cout<<"提前终止条件"<<endl;
            return true;
        }
    }
    return  false;
}

bool check_range(int term_id, int usr_node, int father, map<int,int>& parent_child){
    float _min = 9999999999999;
    int rate;
#ifdef LV
    rate = 1;
#else
    rate = 5;
#endif
    vector<int> related_entry = getObjectTermRelatedEntry(term_id,father);
    for(int sib : related_entry){ //GTree[father].children
        if(sib == parent_child[father]) continue;
        float dis = MinDisBound[usr_node][sib]; // MaxDisBound
        _min = min(_min, dis);

    }
    return (rate*_min>EXPANSION_MAX);
}

LCLResult getUpperNode_termLCL_Disk_Large(int usr_node, int term_id, map<int,vector<int>>& related_poi_nodeSet, int Qk, float a, double alpha, double upper_bound){
    priority_queue<LCLEntry> LCL;
    priority_queue<UCLEntry> UCL;


    int cnt=0;int cnt2 =0;
    int upper_cnt = 0;
    int count = 0;
    //extract relatd poi_leaves

    map<int,int> poiLeaf_include;
#ifdef  LCL_LOG
    cout<<"computing lcl for u node"<<usr_node<<" according to t"<<term_id<<endl;
#endif
    int node_num = related_poi_nodeSet[term_id].size();
    if(node_num >=64){//false
        //cout<<"too many node to be evaluate!, cnt="<<_cnt<<endl;
        TreeNode oNodeData_same = getGIMTreeNodeData(usr_node,OnlyO);
        if(oNodeData_same.term_object_countMap.count(term_id)>0){  //先把相同节点内的object set 加入
            bool flag = false;
            flag = update_LCL(usr_node, usr_node, oNodeData_same, LCL, term_id, cnt, Qk,a,alpha,upper_bound);
            if(flag){
                LCLResult lclresult;
                lclresult.LCL = LCL;
                lclresult.UCL = UCL;
                lclresult.topscore = LCL.top().score;
                return lclresult;
            }
        }


        int current = GTree[usr_node].father; int root =0; int delta_level =1;
        map<int, int> parent_child;  // <father_id, child_id>
        parent_child[current] = usr_node;
        while(current!=root){
            TreeNode oNodeData = getGIMTreeNodeData(current,OnlyO);
            bool relate_flag = (oNodeData.term_object_countMap.count(term_id)>0);
            if(relate_flag && check_range(term_id, usr_node,current, parent_child))  //local expansion and evaluation
                break;
            else{
                int pre = current;
                current = GTree[current].father;
                parent_child[current] = pre;
                delta_level++;
            }
        }
        vector<int> onode_list; onode_list.push_back(current);
        while(delta_level>0){
            vector<int> onode_list_next;
            for(int node: onode_list){
                if(GTree[node].isleaf){
                    onode_list_next.push_back(node);  //将叶节点重新加入队列
                }
                else{
                    vector<int> related = getObjectTermRelatedEntry(term_id,node);
                    onode_list_next.insert(onode_list_next.end(), related.begin(),related.end());
                }

            }
            onode_list_next.swap(onode_list);
            delta_level--;
        }
        for(int onode: onode_list){
            TreeNode oNodeData = getGIMTreeNodeData(onode,OnlyO);
            //printONodeInfo(oNode,node);
            bool flag = false;
            flag = update_LCL(usr_node, onode, oNodeData, LCL, term_id, cnt, Qk,a,alpha,upper_bound);
            if(flag){
                LCLResult lclresult;
                lclresult.LCL = LCL;
                lclresult.UCL = UCL;
                lclresult.topscore = LCL.top().score;
                return lclresult;
            }
        }

    }

    else{
        // access the node in term related_poi_nodeSet to get the LCL of usr node
        for(int onode:related_poi_nodeSet[term_id]){  //首先该usr_node 的关键词情况是怎样的？

            //cout<<"read onode "<<node<<endl;
            TreeNode oNodeData = getGIMTreeNodeData(onode,OnlyO);

            //printONodeInfo(oNode,node);
            bool flag = false;
            flag = update_LCL(usr_node, onode, oNodeData, LCL, term_id, cnt, Qk,a,alpha,upper_bound);
            if(flag){
                LCLResult lclresult;
                lclresult.LCL = LCL;
                lclresult.UCL = UCL;
                lclresult.topscore = LCL.top().score;
                return lclresult;
            }

        }
    }


    LCLResult lclresult;
    lclresult.LCL = LCL;
    lclresult.UCL = UCL;
    lclresult.topscore = LCL.top().score;
    return lclresult;

}


LCLResult getUpperNode_termLCL_Disk(int usr_node, int term_id, map<int,vector<int>>& related_poi_nodeSet, int Qk, float a, double alpha, double upper_bound){
#ifdef LV
    return  getUpperNode_termLCL_Disk_Large(usr_node, term_id, related_poi_nodeSet, Qk, a, alpha, upper_bound);
#else
    return  getUpperNode_termLCL_Disk_Large(usr_node, term_id, related_poi_nodeSet, Qk, a, alpha, upper_bound);
#endif
}




LCLResult getLeafNode_LCL(int usr_leaf, int term, set<int> related_poi_leafSet, map<int,set<int>> leaf_poiCnt, int Qk, float a,double alpha,double upper_boud){
    priority_queue<LCLEntry> LCL;
    priority_queue<LCLEntry2> LCL2;
    int cnt =0;
    int count = 0;

    //extract relatd poi_leaves

    map<int,int> poiLeaf_include;


    int leafAccess = 0;
    for(int leaf:related_poi_leafSet){
        leafAccess++;
        double max_Distance = MaxDisBound[usr_leaf][leaf];//getMinMaxDistanceBetweenNode(usr_leaf,leaf).max;//getMinMaxDistanceOfLeaf(usr_leaf,leaf).max; //100000;//
        //cout<<min_Distance<<endl;
        double min_relevance = 0;
        //double relevance = textRelevance(usr.keywords,qKey);
        min_relevance = tfIDF_term(term);

        double min_influence = 1.0;

        double min_Score = min_relevance * min_influence / (1+max_Distance);

        if(cnt<Qk){
            LCL.push(LCLEntry(leaf,min_Score, leaf_poiCnt[leaf].size()));
            cnt += leaf_poiCnt[leaf].size();
        }
        else{
            while(min_Score > LCL.top().score && cnt > Qk){
                LCLEntry e =  LCL.top();
                LCL.pop();
                cnt -= e.count;
            }
            LCL.push(LCLEntry(leaf, min_Score, leaf_poiCnt[leaf].size()));
            cnt += leaf_poiCnt[leaf].size();
            if(upper_boud < LCL.top().score && cnt > Qk){
                //cout<<"early, leafAccess:"<<leafAccess<<endl;
                LCLResult lclresult;
                lclresult.LCL = LCL;
                lclresult.topscore = LCL.top().score;
                return lclresult;
            }

        }
    }
    //cout<<"final, leafAccess:"<<leafAccess<<endl;
    LCLResult lclresult;
    lclresult.LCL = LCL;
    lclresult.topscore = LCL.top().score;
    return lclresult;

}



LCLResult getLeafNode_LCL_object(int usr_leaf, int term_id,  priority_queue<LCLEntry> LCL,  map<int,set<int>> node_poiCnt, int Qk, float a,double alpha,double gsk){
    priority_queue<LCLEntry> leaf_lcl;
    //priority_queue<ResultCurrent> u_lcl;
    int cnt =0;
    int count = 0;

    map<int,int> poiLeaf_include;

    double score = 0;
    while(LCL.size()>0){
        LCLEntry le = LCL.top();
        LCL.pop();
        int n = le.id;
        //node lcl 表记录的是用户节点id
        for(int o_id : node_poiCnt[n]){

            double score_lower = getGSKScoreo2U_Lower(a, alpha, o_id, usr_leaf,term_id);
            leaf_lcl.push(LCLEntry(o_id,score_lower,1));
            if(leaf_lcl.size() > Qk){
                leaf_lcl.pop();
                score  = leaf_lcl.top().score;
                if(score > gsk){
                    LCLResult lclresult;
                    lclresult.LCL = leaf_lcl;
                    lclresult.topscore = score;
                    return lclresult;
                }
            }
        }
    }
    LCLResult lclresult;
    lclresult.LCL = leaf_lcl;
    lclresult.topscore = leaf_lcl.top().score;
    return lclresult;

}

LCLResult getLeafNode_UCL_object(int usr_leaf, int term_id,  priority_queue<LCLEntry> LCL,  map<int,set<int>> node_poiCnt, int Qk, float a,double alpha, double gsk){
    priority_queue<LCLEntry> leaf_lcl;
    priority_queue<UCLEntry> leaf_ucl;
    //priority_queue<ResultCurrent> u_lcl;
    int cnt =0;
    int count = 0;

    //extract relatd poi_leaves

    map<int,int> poiLeaf_include;

    int size = LCL.size();
    while(LCL.size()>0){
        LCLEntry le = LCL.top();
        LCL.pop();
        int n = le.id;
        //node lcl 表记录的是用户节点id
        for(int o_id : node_poiCnt[n]){
            //double score = 0;
            //double score = getGSKScoreO2U(a, o_id, usr_id);  //这里最好用下界：路网距离的下界为欧式距离
            double score_lower = getGSKScoreo2U_Lower(a, alpha, o_id, usr_leaf,term_id);
            double score_upper = getGSKScoreo2U_Upper(a, alpha, o_id, usr_leaf);
            leaf_ucl.push(UCLEntry(o_id,score_upper,1));

            if(leaf_lcl.size()<Qk)
                leaf_lcl.push(LCLEntry(o_id,score_lower,1));
            else{
                leaf_lcl.push(LCLEntry(o_id,score_lower,1));
                leaf_lcl.pop();
                LCLEntry rc = leaf_lcl.top();
                if(rc.score > gsk){
                    LCLResult lclresult;
                    lclresult.LCL = leaf_lcl;
                    lclresult.topscore = rc.score;
                    return lclresult;
                }
            }
        }
    }
    LCLResult lclresult;
    lclresult.LCL = leaf_lcl;
    lclresult.UCL = leaf_ucl;
    lclresult.topscore = leaf_lcl.top().score;
    return lclresult;

}



LCLResult getLeafNode_TopkSDijkstra(int leaf_id, int term_id, int Qk,float a, float alpha,float score_bound) {

    int nodeCnt = 0;
    float Max_D = 25000;

    bool vis[VertexNum];
    bool rearch[VertexNum];
    bool poi_rearch[poi_num];
    float pointval[VertexNum];

    memset(vis, false, sizeof(vis));
    memset(rearch, false, sizeof(rearch));
    float maxDisValue = 9999999.9;
    //memset(pointval, maxDisValue, sizeof(float)*sizeof(pointval));
    //for (int i = 0; i < VertexNum; ++i)
    //pointval[i] = maxNum;

    priority_queue<nodeDistType> tmpPoint;
    priority_queue<LCLEntry> lcl;
    LCLResult result;
    float maxDisInsideLeaf = MaxDisBound[leaf_id][leaf_id];
    int leaf_border = GTree[leaf_id].borders[0];


    int qNi = leaf_border;
    float qDis = maxDisInsideLeaf;
    //cout<<maxDisInsideLeaf<<endl;
    vector<int> qKey;
    qKey.push_back(term_id);

    // expand the edge from Ni, Nj
    pointval[qNi] = qDis;

    nodeDistType tmp1(qNi, qDis);
    tmpPoint.push(tmp1);
    double  optimal_text_score = textRelevance(qKey, qKey);
    while (!tmpPoint.empty()) {
        nodeDistType tmp = tmpPoint.top();
        tmpPoint.pop();
        int tmpNi = tmp.ID;
        if (vis[tmpNi] == true)
            continue;
        nodeCnt++;
        vis[tmpNi] = true;

        vector<int> tmpAdj = adjList[tmpNi];
        for (size_t i = 0; i < tmpAdj.size(); i++) {  // expand each edge
            int tmpNj = tmpAdj[i];
            if (vis[tmpNj])
                continue;

            edge e = EdgeMap[getKey(tmpNi, tmpNj)];
            int Ni = e.Ni;
            float edgeDist = e.dist;
            vector<POI> poiList = e.pts;
            for (size_t j = 0; j < poiList.size(); j++) { // check each poi on the edge
                POI tmpPOI = poiList[j];

                //tmpPOI = POIs[60];

                if (tmpPOI.category == 0)
                    continue;

                int pId = tmpPOI.id;

                float pDis = tmpPOI.dis;
                vector<int> pKey = tmpPOI.keywords;
                //这里待加入社交

                float distance = -1;

                bool relevance = textCover(term_id, pKey);

                //必须文本相关
                if (!relevance) continue;

                double  optimal_text_score = textRelevance(qKey, qKey);

                // user and poi are on the same edge
                if ((tmpNi == qNi) || (tmpNj == qNi)) {
                    if (qDis > pDis)
                        distance = qDis - pDis;
                    else
                        distance = pDis - qDis;
                }
                    //user and poi are on the different edge
                else {
                    if (tmpNi == e.Ni)  //同一个点
                        distance = pointval[tmpNi] + pDis;
                    else //边的另一个点
                        distance = pointval[tmpNi] + edgeDist - pDis;

                }

                double simD = distance;

                float simT = 0;
                simT = tfIDF_term(term_id);

                //double pow(double x,double y) //求x的y次方
                float simS;
                //double social_score = getSocialScore(qId,pId);
                double social_score = 0;
                /*if(social_score == 0)
                    simS = 1.0;
                else
                    simS = 1.0 + pow(a,social_score);*/
                if(social_score > 0.0)
                    simS = 1.0 + pow(a,social_score);
                else
                    simS = 1.0;
                double social_textual_score = alpha*simS + (1-alpha)*simT;
                double score = social_textual_score / (1+simD);

                LCLEntry le(pId, score, 1);
                //cout<<"qId="<<qId<<"pId="<<pId<<endl;
                if (lcl.size() < Qk) {
                    lcl.push(le);
                    //cout<<"simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;
                }
                else {
                    LCLEntry _tmpNode = lcl.top();
                    float tmpScore = _tmpNode.score;

                    if (tmpScore < score) {
                        lcl.pop();
                        lcl.push(_tmpNode);
                        //cout<<"simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;
                    }
                    if(lcl.top().score > score_bound){//当前已有更优的k个结果，拓展可提前终止

                        result.LCL = lcl;
                        result.topscore = lcl.top().score;
                        return result;  //注意此时的score只是current score
                    }
                }
            }
            //首次拓展当前拓展顶点
            if(rearch[tmpNj]==false){
                pointval[tmpNj] = maxDisValue;
            }
            if (e.dist + tmp.dist < pointval[tmpNj]) {
                pointval[tmpNj] = e.dist + tmp.dist;
                nodeDistType tmp2(tmpNj, edgeDist + tmp.dist);
                //cout << tmp2.num << ' ' << tmp2.dist << endl;
                tmpPoint.push(tmp2);
                rearch[tmpNj]=true;
            }
        }
        //路网边拓展提前终止条件
        if (lcl.size() == Qk) {
            float tmpScore = lcl.top().score;
            float currentDistance = pointval[tmpNi];
            double optimal_social_textual_score = alpha*(a+1.0)+(1-alpha)*optimal_text_score;
            double bound = optimal_social_textual_score /tmpScore - 1;               //提前拓展终止条件出错
            if (currentDistance > bound){
                //cout<<"当前拓展终止"<<endl;
                break;
            }

        }
    }
    result.LCL = lcl;
    result.topscore = lcl.top().score;
    return result;
}



//jordan
//内存操作
MultiLCLResult getLeafNode_LCL_Dijkstra_multi_memory(int leaf_id, vector<int> Ukeys, int Qk,float a, float alpha,float score_bound) {

    int nodeCnt = 0;
    //float Max_D = 25000;


    set<int> Ukeyset;
    for(int keyword: Ukeys)
        Ukeyset.insert(keyword);


    bool vis[VertexNum];
    bool rearch[VertexNum];
    bool poi_rearch[poi_num];
    float pointval[VertexNum];

    memset(vis, false, sizeof(vis));
    memset(rearch, false, sizeof(rearch));
    float maxDisValue = 9999999.9;
    //memset(pointval, maxDisValue, sizeof(float)*sizeof(pointval));
    //for (int i = 0; i < VertexNum; ++i)
    //pointval[i] = maxNum;

    priority_queue<nodeDistType> tmpPoint;

    vector<priority_queue<LCLEntry>>  multi_Q;
    MultiLCLResult lcl_multi;
    set<int> stop_flag;


    for(int i=0;i<Ukeys.size();i++){
        priority_queue<LCLEntry> lcl_queue;
        multi_Q.push_back(lcl_queue);

        LCLResult lclresult;
        lcl_multi.push_back(lclresult);
    }
    priority_queue<LCLEntry>& lcl_q = multi_Q[0];

    //kingmura
    vector<vector<int>> qKey_multi;
    vector<double> optimal_text_score_list;
    for(int term: Ukeys){
        vector<int> qKey2;
        qKey2.push_back(term);
        qKey_multi.push_back(qKey2);

        double  optimal_text_score = textRelevance(qKey2, qKey2);

        optimal_text_score_list.push_back(optimal_text_score);
    }


    //选择usr leaf 节点内部最大距离
    float maxDisInsideLeaf = MaxDisBound[leaf_id][leaf_id];
    //选择该usr leaf 节点上的任意一个 border vertex进行 dj expansion
    int leaf_border = GTree[leaf_id].borders[0];


    int qNi = leaf_border;
    float qDis = maxDisInsideLeaf;
    //cout<<maxDisInsideLeaf<<endl;


    // expand the edge from Ni, Nj
    pointval[qNi] = qDis;

    nodeDistType tmp1(qNi, qDis);
    tmpPoint.push(tmp1);




    while (!tmpPoint.empty()) {
        nodeDistType tmp = tmpPoint.top();
        tmpPoint.pop();
        int tmpNi = tmp.ID;
        if (vis[tmpNi] == true)
            continue;
        nodeCnt++;
        vis[tmpNi] = true;

        vector<int> tmpAdj = adjList[tmpNi];

        //当前关键词 weight  idx
        double  optimal_text_score = -1;
        int opitmal_term_idx = -1;
        set<int>  term_idx;

        //int AdjGrpAddr; int AdjListSize;
        //AdjGrpAddr = getAdjListGrpAddr(tmpNi);

        int AdjListSize= tmpAdj.size();



        for (size_t j = 0; j < AdjListSize; j++) {  // expand each edge
            int tmpNj = tmpAdj[j];

            //getVarE(ADJNODE_A, Ref(tmpNj),AdjGrpAddr,j);

            if (vis[tmpNj])
                continue;

            edge e = EdgeMap[getKey(tmpNi, tmpNj)];
            int Ni = e.Ni;
            float edgeDist = e.dist;

            //读取Okey psudo doc 的地址,并赋予okeyAddr
            //int okeyAddr = -1;
            //getVarE(OKEY_A, Ref(okeyAddr),AdjGrpAddr,j); //AdjGrpAddr为Ni顶点在邻接表文件中的数据块基地址
            //根据okeyAddr读 edge上的 okeywords
            set<int> OKeyset;
            OKeyset = e.OUnionKeys;


            set<int> Inter;
            Inter = obtain_itersection_jins(Ukeyset,OKeyset);  //ala


            if(Inter.size() > 0){  //_relevance>0.0

                //读取边上poi集的地址(成功）
                //int PoiGrpAddress =-1; getVarE(PTKEY_A, Ref(PoiGrpAddress),AdjGrpAddr,i);  //读取边上poi集的地址(成功）
                //读取边上poi集的兴趣点个数
                int poi_size = 0;  int adj_th = j;
                poi_size = e.pts.size();
                //cout<<"poiList size= "<<poiList.size()<<endl;
                if(poi_size>0){//(PoiGrpAddress !=-1){

                    for (size_t p_th = 0; p_th < poi_size; p_th++) { // check each poi on the edge
                        /*int pId = 0;    //------------取得第n个poi的id
                        int poi_indx_address = getPoiIdxAddr(PoiIdxAddress,j);  //索引表所在的地址
                        int poi_address;
                        getPOIAddr(&poi_address, poi_indx_address);  //获取poi数据所在地址
                        */

                        //int idx_baseAddr = getPOIIdxAddrOfEdge(tmpNi,adj_th);
                        POI p = e.pts[p_th];

                        //POI p = getPOIDATA(poi_address);

                        int pId = p.id;
                        float pDis = p.dis;   //------------取得第n个poi的距离Ni的距离
                        vector<int> pKey = p.keywords;//tmpPOI.keywords;
                        //这里待加入社交

                        float distance = -1;

                        bool relevance = false;
                        for(int i=0;i<qKey_multi.size();i++){
                            //对每个关键词进行检查
                            vector<int> _qKey= qKey_multi[i];
                            relevance = textCover(_qKey[0], pKey);
                            if(relevance) break;
                        }

                        //必须文本相关
                        if (!relevance) continue;

                        optimal_text_score = -1;
                        opitmal_term_idx = -1;
                        //vector

                        term_idx.clear();

                        for(int i=0;i<qKey_multi.size();i++){
                            //对每个关键词进行检查
                            double _score = textRelevance(qKey_multi[i], pKey);
                            if(_score>0){
                                term_idx.insert(i);
                            }
                        }

                        //textRelevance(qKey, qKey);

                        // user and poi are on the same edge
                        if ((tmpNi == qNi) || (tmpNj == qNi)) {
                            if (qDis > pDis)
                                distance = qDis - pDis;
                            else
                                distance = pDis - qDis;
                        }
                            //user and poi are on the different edge
                        else {
                            if (tmpNi == e.Ni)  //同一个点
                                distance = pointval[tmpNi] + pDis;
                            else //边的另一个点
                                distance = pointval[tmpNi] + edgeDist - pDis;

                        }

                        double simD = distance;

                        float simT = 0;
                        for(int idx: term_idx){
                            int term_id = qKey_multi[idx][0];
                            simT = tfIDF_term(term_id);

                            float simS;

                            double social_score = 0;

                            if(social_score > 0.0)
                                simS = 1.0 + pow(a,social_score);
                            else
                                simS = 1.0;
                            double social_textual_score = alpha*simS + (1-alpha)*simT;
                            double score = social_textual_score / (1+simD);

                            LCLEntry le(pId, score, 1);

                            if (multi_Q[idx].size() < Qk) {
                                multi_Q[idx].push(le);
                                //cout<<"simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;
                            }
                            else {
                                LCLEntry _tmpNode = multi_Q[idx].top();
                                float tmpScore = _tmpNode.score;

                                if (tmpScore < score) {
                                    multi_Q[idx].pop();
                                    multi_Q[idx].push(_tmpNode);
                                    //cout<<"simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;
                                }

                            }
                        }

                    }
                }

            }
            else{
                //cout<<"there is no related object located on this edage"<<endl;
            }

            //首次拓展当前拓展顶点
            if(rearch[tmpNj]==false){
                pointval[tmpNj] = maxDisValue;
            }
            if (e.dist + tmp.dist < pointval[tmpNj]) {
                pointval[tmpNj] = e.dist + tmp.dist;
                nodeDistType tmp2(tmpNj, edgeDist + tmp.dist);
                //cout << tmp2.num << ' ' << tmp2.dist << endl;
                tmpPoint.push(tmp2);
                rearch[tmpNj]=true;
            }
        }
        //路网边拓展提前终止条件
        float currentDistance = pointval[tmpNi];

        for(int i = 0;i<Ukeys.size();i++){

            if (multi_Q[i].size()==Qk) {

                stop_flag.insert(i);//jordan
                if(stop_flag.size()==Ukeys.size()&&currentDistance> EXPANSION_MAX) {

                    for(int i=0;i<Ukeys.size();i++){
                        LCLResult result;
                        result.LCL = multi_Q[i];
                        if(multi_Q[i].size()==Qk)
                            result.topscore = multi_Q[i].top().score;
                        else{
                            result.topscore = 0;
                            cout<<"lcl for"<<Ukeys[i]<<"less than"<<Qk<<endl;
                        }

                        lcl_multi[i]= result;

                    }
                    //cout<<"lcl is full break"<<endl;
                    return lcl_multi;

                }

            }

        }


    }//end while

    for(int i=0;i<Ukeys.size();i++){
        LCLResult result;
        result.LCL = multi_Q[i];
        if(multi_Q[i].size()==Qk)
            result.topscore = multi_Q[i].top().score;
        else{
            result.topscore = 0;
            //cout<<"lcl for"<<Ukeys[i]<<"less than"<<Qk<<endl;
        }

        lcl_multi[i]= result;

    }

    return lcl_multi;
}

// 外存操作
MultiLCLResult getLeafNode_LCL_Dijkstra_multi_disk_0617(int leaf_id, vector<int> Ukeys, int Qk,float a, float alpha,float score_bound) {

    int nodeCnt = 0;
    //float Max_D = 25000;


    set<int> Ukeyset;
    for(int keyword: Ukeys)
        Ukeyset.insert(keyword);


    bool vis[VertexNum];
    bool rearch[VertexNum];
    bool poi_rearch[poi_num];
    float pointval[VertexNum];

    memset(vis, false, sizeof(vis));
    memset(rearch, false, sizeof(rearch));
    float maxDisValue = 9999999.9;
    //memset(pointval, maxDisValue, sizeof(float)*sizeof(pointval));
    //for (int i = 0; i < VertexNum; ++i)
    //pointval[i] = maxNum;

    priority_queue<nodeDistType> tmpPoint;

    vector<priority_queue<LCLEntry>>  multi_Q;
    MultiLCLResult lcl_multi;
    set<int> stop_flag;

    map<int, int > idx_termIdMap;
    map<int, int > id_termTHMap;
    for(int i=0;i<Ukeys.size();i++){
        int term = Ukeys[i];
        idx_termIdMap[i]=term;
        id_termTHMap[term]= i;
        priority_queue<LCLEntry> lcl_queue;
        multi_Q.push_back(lcl_queue);

        LCLResult lclresult;
        lcl_multi.push_back(lclresult);
    }
    priority_queue<LCLEntry>& lcl_q = multi_Q[0];

    //kingmura
    vector<vector<int>> qKey_multi;
    vector<double> optimal_text_score_list;
    for(int term: Ukeys){
        vector<int> qKey2;
        qKey2.push_back(term);
        qKey_multi.push_back(qKey2);

        double  optimal_text_score = textRelevance(qKey2, qKey2);

        optimal_text_score_list.push_back(optimal_text_score);
    }


    //选择usr leaf 节点内部最大距离
    float maxDisInsideLeaf = MaxDisBound[leaf_id][leaf_id];
    //选择该usr leaf 节点上的任意一个 border vertex进行 dj expansion
    int leaf_border = GTree[leaf_id].borders[0];


    int qNi = leaf_border;
    float qDis = maxDisInsideLeaf;
    //cout<<maxDisInsideLeaf<<endl;


    // expand the edge from Ni, Nj
    pointval[qNi] = qDis;

    nodeDistType tmp1(qNi, qDis);
    tmpPoint.push(tmp1);


    while (!tmpPoint.empty()) {
        nodeDistType tmp = tmpPoint.top();
        tmpPoint.pop();
        int tmpNi = tmp.ID;
        if (vis[tmpNi] == true)
            continue;
        nodeCnt++;
        vis[tmpNi] = true;

        vector<int> tmpAdj = adjList[tmpNi];

        //当前关键词 weight  idx
        double  optimal_text_score = -1;
        int opitmal_term_idx = -1;
        set<int>  term_idx;

        int AdjGrpAddr; int AdjListSize;
        AdjGrpAddr = getAdjListGrpAddr(tmpNi);

        getFixedF(SIZE_A,Ref(AdjListSize), AdjGrpAddr);


        for (size_t j = 0; j < AdjListSize; j++) {  // expand each edge
            int tmpNj;

            getVarE(ADJNODE_A, Ref(tmpNj),AdjGrpAddr,j);

            if (vis[tmpNj])
                continue;

            edge e = EdgeMap[getKey(tmpNi, tmpNj)];
            int Ni = e.Ni;
            float edgeDist = e.dist;


            set<int> OKeyset;
            //OKeyset = getOKeyEdge(okeyAddr);
            OKeyset = getOKeyEdge(tmpNi,j);


            set<int> Inter;
            Inter = obtain_itersection_jins(Ukeyset,OKeyset);  //ala


            if(Inter.size() > 0){  //_relevance>0.0

                //读取边上poi集的地址(成功）
                //int PoiGrpAddress =-1; getVarE(PTKEY_A, Ref(PoiGrpAddress),AdjGrpAddr,i);  //读取边上poi集的地址(成功）
                //读取边上poi集的兴趣点个数
                int poi_size = 0;  int adj_th = j;
                poi_size = getPOISizeOfEdge(tmpNi,adj_th);
                //cout<<"poiList size= "<<poiList.size()<<endl;
                if(poi_size>0){//(PoiGrpAddress !=-1){
                    //读取poi索引位置
                    /*
                    int PoiIdxAddress =-1; getVarE(PIDX_A, Ref(PoiIdxAddress),AdjGrpAddr,i);
                    poi_size = -3;
                    int Ni=-1; int Nj=-1; float Ddist= 0;
                    //读取该兴趣点信息，首先度边信息：Ni,Nj, Dist, poiGrpAddress--------------------------------
                    getFixedF(NI_P,Ref(Ni),PoiGrpAddress);    //
                    getFixedF(NJ_P,Ref(Nj),PoiGrpAddress);
                    //cout <<" Ni="<<Ni<<",Nj="<<Nj<<endl;
                    getFixedF(DIST_P,Ref(Ddist),PoiGrpAddress); //读取边上poi个数
                    //cout<<"Ddist="<<Ddist<<endl;
                    getFixedF(SIZE_P,Ref(poi_size),PoiGrpAddress); //读取边上poi个数
                    //cout<<"poi_size="<<poi_size<<endl;*/
                    for (size_t p_th = 0; p_th < poi_size; p_th++) { // check each poi on the edge
                        /*int pId = 0;    //------------取得第n个poi的id
                        int poi_indx_address = getPoiIdxAddr(PoiIdxAddress,j);  //索引表所在的地址
                        int poi_address;
                        getPOIAddr(&poi_address, poi_indx_address);  //获取poi数据所在地址
                        */

                        //int idx_baseAddr = getPOIIdxAddrOfEdge(tmpNi,adj_th);
                        POI p = getPOIDataOnEdge(tmpNi,adj_th, p_th);

                        //POI p = getPOIDATA(poi_address);

                        int pId = p.id;
                        float pDis = p.dis;   //------------取得第n个poi的距离Ni的距离
                        vector<int> pKey = p.keywords;//tmpPOI.keywords;
                        //这里待加入社交

                        float distance = -1;

                        bool relevance = false;
                        for(int i=0;i<qKey_multi.size();i++){
                            //对每个关键词进行检查
                            vector<int> _qKey= qKey_multi[i];
                            relevance = textCover(_qKey[0], pKey);
                            if(relevance) break;
                        }

                        //必须文本相关
                        if (!relevance) continue;

                        optimal_text_score = -1;
                        opitmal_term_idx = -1;
                        //vector

                        term_idx.clear();

                        for(int i=0;i<qKey_multi.size();i++){
                            //对每个关键词进行检查
                            double _score = textRelevance(qKey_multi[i], pKey);
                            if(_score>0){
                                term_idx.insert(i);
                            }
                        }

                        //textRelevance(qKey, qKey);

                        // user and poi are on the same edge
                        if ((tmpNi == qNi) || (tmpNj == qNi)) {
                            if (qDis > pDis)
                                distance = qDis - pDis;
                            else
                                distance = pDis - qDis;
                        }
                            //user and poi are on the different edge
                        else {
                            if (tmpNi == e.Ni)  //同一个点
                                distance = pointval[tmpNi] + pDis;
                            else //边的另一个点
                                distance = pointval[tmpNi] + edgeDist - pDis;

                        }

                        double simD = distance;

                        float simT = 0;
                        for(int idx: term_idx){
                            int term_id = qKey_multi[idx][0];
                            simT = tfIDF_term(term_id);

                            float simS;

                            double social_score = 0;

                            if(social_score > 0.0)
                                simS = 1.0 + pow(a,social_score);
                            else
                                simS = 1.0;
                            double social_textual_score = alpha*simS + (1-alpha)*simT;
                            double score = social_textual_score / (1+simD);

                            LCLEntry le(pId, score, 1);

                            if (multi_Q[idx].size() < Qk) {
                                multi_Q[idx].push(le);
                                //cout<<"simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;
                            }
                            else if(multi_Q[idx].size() == Qk){
                                LCLEntry _tmpNode = multi_Q[idx].top();
                                multi_Q[idx].push(le);
                                multi_Q[idx].pop();

                            }
                        }

                    }
                }

            }
            else{
                //cout<<"there is no related object located on this edage"<<endl;
            }

            //首次拓展当前拓展顶点
            if(rearch[tmpNj]==false){
                pointval[tmpNj] = maxDisValue;
            }
            if (e.dist + tmp.dist < pointval[tmpNj]) {
                pointval[tmpNj] = e.dist + tmp.dist;
                nodeDistType tmp2(tmpNj, edgeDist + tmp.dist);
                //cout << tmp2.num << ' ' << tmp2.dist << endl;
                tmpPoint.push(tmp2);
                rearch[tmpNj]=true;
            }
        }
        //路网边拓展提前终止条件
        float currentDistance = pointval[tmpNi];

        for(int i = 0;i<Ukeys.size();i++){

            if (multi_Q[i].size()==Qk) {

                stop_flag.insert(i);//jordan
                if(stop_flag.size()==Ukeys.size()&&currentDistance> EXPANSION_MAX) {

                    for(int i=0;i<Ukeys.size();i++){
                        LCLResult result;
                        result.LCL = multi_Q[i];
                        if(multi_Q[i].size()==Qk)
                            result.topscore = multi_Q[i].top().score;
                        else{
                            result.topscore = 0;
                            cout<<"lcl for"<<Ukeys[i]<<"less than"<<Qk<<endl;
                        }

                        lcl_multi[i]= result;

                    }
                    //cout<<"lcl is full break"<<endl;
                    return lcl_multi;

                }

            }

        }


    }//end while

    for(int i=0;i<Ukeys.size();i++){
        LCLResult result;
        result.LCL = multi_Q[i];
        if(multi_Q[i].size()==Qk)
            result.topscore = multi_Q[i].top().score;
        else{
            result.topscore = 0;
            //cout<<"lcl for"<<Ukeys[i]<<"less than"<<Qk<<endl;
        }

        lcl_multi[i]= result;

    }

    return lcl_multi;
}


MultiLCLResult getLeafNode_LCL_Dijkstra_multi_disk(int leaf_id, vector<int> Ukeys, int Qk,float a, float alpha,float score_bound) {

    int nodeCnt = 0;
    //float Max_D = 25000;


    set<int> Ukeyset;
    for(int keyword: Ukeys)
        Ukeyset.insert(keyword);


    bool vis[VertexNum];
    bool rearch[VertexNum];
    bool poi_rearch[poi_num];
    float pointval[VertexNum];

    memset(vis, false, sizeof(vis));
    memset(rearch, false, sizeof(rearch));
    float maxDisValue = 9999999.9;
    //memset(pointval, maxDisValue, sizeof(float)*sizeof(pointval));
    //for (int i = 0; i < VertexNum; ++i)
    //pointval[i] = maxNum;

    priority_queue<nodeDistType> tmpPoint;

    vector<priority_queue<LCLEntry>>  multi_Q;
    MultiLCLResult lcl_multi;
    set<int> stop_flag;

    map<int, int > idx_termIdMap;  //关键词索引表
    map<int, int > id_termTHMap;
    for(int i=0;i<Ukeys.size();i++){
        int term = Ukeys[i];
        idx_termIdMap[i]=term;
        id_termTHMap[term]= i;
        priority_queue<LCLEntry> lcl_queue;
        multi_Q.push_back(lcl_queue);

        LCLResult lclresult;
        lcl_multi.push_back(lclresult);
    }
    priority_queue<LCLEntry>& lcl_q = multi_Q[0];

    //kingmura
    vector<vector<int>> qKey_multi;
    vector<double> optimal_text_score_list;
    for(int term: Ukeys){
        vector<int> qKey2;
        qKey2.push_back(term);
        qKey_multi.push_back(qKey2);

        double  optimal_text_score = textRelevance(qKey2, qKey2);

        optimal_text_score_list.push_back(optimal_text_score);
    }


    //选择usr leaf 节点内部最大距离
    float maxDisInsideLeaf = MaxDisBound[leaf_id][leaf_id];
    //选择该usr leaf 节点上的任意一个 border vertex进行 dj expansion
    int leaf_border = GTree[leaf_id].borders[0];


    int qNi = leaf_border;
    float qDis = maxDisInsideLeaf;
    //cout<<maxDisInsideLeaf<<endl;


    // expand the edge from Ni, Nj
    pointval[qNi] = qDis;

    nodeDistType tmp1(qNi, qDis);
    tmpPoint.push(tmp1);


    while (!tmpPoint.empty()) {
        nodeDistType tmp = tmpPoint.top();
        tmpPoint.pop();
        int tmpNi = tmp.ID;
        if (vis[tmpNi] == true)
            continue;
        nodeCnt++;
        vis[tmpNi] = true;

        vector<int> tmpAdj = adjList[tmpNi];

        //当前关键词 weight  idx
        double  optimal_text_score = -1;
        int opitmal_term_idx = -1;
        set<int>  term_idx;

        int AdjGrpAddr; int AdjListSize;
        AdjGrpAddr = getAdjListGrpAddr(tmpNi);

        getFixedF(SIZE_A,Ref(AdjListSize), AdjGrpAddr);


        for (size_t j = 0; j < AdjListSize; j++) {  // expand each edge
            int tmpNj;

            getVarE(ADJNODE_A, Ref(tmpNj),AdjGrpAddr,j);

            if (vis[tmpNj])
                continue;

            edge e = EdgeMap[getKey(tmpNi, tmpNj)];
            int Ni = e.Ni;
            float edgeDist = e.dist;


            set<int> OKeyset;
            //OKeyset = getOKeyEdge(okeyAddr);
            OKeyset = getOKeyEdge(tmpNi,j);


            set<int> Inter;
            Inter = obtain_itersection_jins(Ukeyset,OKeyset);  //ala


            if(Inter.size() > 0){  //_relevance>0.0

                //读取边上poi集的地址(成功）
                //int PoiGrpAddress =-1; getVarE(PTKEY_A, Ref(PoiGrpAddress),AdjGrpAddr,i);  //读取边上poi集的地址(成功）
                //读取边上poi集的兴趣点个数
                int poi_size = 0;  int adj_th = j;
                poi_size = getPOISizeOfEdge(tmpNi,adj_th);
                //cout<<"poiList size= "<<poiList.size()<<endl;
                if(poi_size>0){//(PoiGrpAddress !=-1){

                    for (size_t p_th = 0; p_th < poi_size; p_th++) { // check each poi on the edge

                        POI p = getPOIDataOnEdge(tmpNi,adj_th, p_th);


                        int pId = p.id;
                        float pDis = p.dis;   //------------取得第n个poi的距离Ni的距离
                        vector<int> pKey = p.keywords;//tmpPOI.keywords;
                        //这里待加入社交

                        float distance = -1;

                        bool relevance = false;
                        set<int> interKeyPOI = obtain_itersection_jins(Inter, p.keywordSet);


                        //必须文本相关
                        if (!interKeyPOI.size()>0) continue;

                        optimal_text_score = -1;
                        opitmal_term_idx = -1;
                        //vector

                        term_idx.clear();


                        //textRelevance(qKey, qKey);

                        // user and poi are on the same edge
                        if ((tmpNi == qNi) || (tmpNj == qNi)) {
                            if (qDis > pDis)
                                distance = qDis - pDis;
                            else
                                distance = pDis - qDis;
                        }
                            //user and poi are on the different edge
                        else {
                            if (tmpNi == e.Ni)  //同一个点
                                distance = pointval[tmpNi] + pDis;
                            else //边的另一个点
                                distance = pointval[tmpNi] + edgeDist - pDis;

                        }

                        double simD = distance;

                        float simT = 0;
                        for(int term_id: interKeyPOI){

                            int term_idx = id_termTHMap[term_id];
                            simT = tfIDF_term(term_id);

                            float simS;

                            double social_score = 0;

                            if(social_score > 0.0)
                                simS = 1.0 + pow(a,social_score);
                            else
                                simS = 1.0;
                            double social_textual_score = alpha*simS + (1-alpha)*simT;
                            double score = social_textual_score / (1+simD);

                            LCLEntry le(pId, score, 1);

                            if (multi_Q[term_idx].size() < Qk) {
                                multi_Q[term_idx].push(le);
                                //cout<<"simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;
                            }
                            else if(multi_Q[term_idx].size() == Qk){
                                int size = multi_Q[term_idx].size();
                                //cout<<"Qk="<<Qk<<",size="<<size<<endl;
                                LCLEntry _tmpNode = multi_Q[term_idx].top();
                                multi_Q[term_idx].push(le);
                                multi_Q[term_idx].pop();

                            }
                        }

                    }
                }

            }
            else{
                //cout<<"there is no related object located on this edage"<<endl;
            }

            //首次拓展当前拓展顶点
            if(rearch[tmpNj]==false){
                pointval[tmpNj] = maxDisValue;
            }
            if (e.dist + tmp.dist < pointval[tmpNj]) {
                pointval[tmpNj] = e.dist + tmp.dist;
                nodeDistType tmp2(tmpNj, edgeDist + tmp.dist);
                //cout << tmp2.num << ' ' << tmp2.dist << endl;
                tmpPoint.push(tmp2);
                rearch[tmpNj]=true;
            }
        }
        //路网边拓展提前终止条件
        float currentDistance = pointval[tmpNi];

        for(int i = 0;i<Ukeys.size();i++){

            if (multi_Q[i].size()==Qk) {

                stop_flag.insert(i);//jordan
                if(stop_flag.size()==Ukeys.size()&&currentDistance> EXPANSION_MAX) {

                    for(int i=0;i<Ukeys.size();i++){
                        LCLResult result;
                        result.LCL = multi_Q[i];
                        if(multi_Q[i].size()==Qk)
                            result.topscore = multi_Q[i].top().score;
                        else{
                            result.topscore = 0;
                            cout<<"lcl for"<<Ukeys[i]<<"less than"<<Qk<<endl;
                        }

                        lcl_multi[i]= result;

                    }
                    //cout<<"lcl is full break"<<endl;
                    return lcl_multi;

                }

            }

        }


    }//end while

    for(int i=0;i<Ukeys.size();i++){
        LCLResult result;
        result.LCL = multi_Q[i];
        if(multi_Q[i].size()==Qk)
            result.topscore = multi_Q[i].top().score;
        else{
            result.topscore = 0;
            //cout<<"lcl for"<<Ukeys[i]<<"less than"<<Qk<<endl;
        }

        lcl_multi[i]= result;

    }

    return lcl_multi;
}



MultiLCLResult getLeafNode_LCL_Dijkstra_multi(int leaf_id, vector<int> Ukeys, int Qk,float a, float alpha,float score_bound){
#ifndef DiskAccess
    return getLeafNode_LCL_Dijkstra_multi_memory(leaf_id, Ukeys, Qk,a, alpha,score_bound);
#else
    return getLeafNode_LCL_Dijkstra_multi_disk(leaf_id, Ukeys, Qk,a, alpha,score_bound);
#endif
}




//jordan_leaf_lcl_bottom

void UpdateR_usrLeaf_Memory(int usr_leaf, int locid, int a, double alpha, double score_bound, int& node_highest,int& node_pos, vector<int> keywords, unordered_map<int, vector<int> >& itm, priority_queue<QEntry>& Queue, double& score_max, bool& eary_terminate){

#ifdef  TRACK
    cout<<"表示当前范围内已全部遍历， 需扩大遍历子图范围"<<endl;
#endif
    //表示当前范围内已全部遍历， 需扩大遍历子图范围
    int posa, posb, _min, dis;
    int root_id =0;
    int node_highest_pre = node_highest;
    int second_highest =0;
    //int node_highest = GTree[node_highest].father;
    int current = node_highest; int current_pre = current;
    set<int> entry_set;
    //
    vector<int> Ukeys = keywords;

    //cout<<"从"<<current_pre<<"向上寻找"<<endl;
    while(current != root_id){
        // 获得该中间节点下包含用户关键词的所有child
        current_pre = current;
        current = GTree[current].father;
        node_pos--;   //最高层节点层数往上一层
        for (int term: keywords) {
            if (GTree[current].inverted_list_o[term].size() > 0) {
                for (int child: GTree[current].inverted_list_o[term]) {
                    entry_set.insert(child);
                    //block_num++;
                }
            }
        }
        //父节点就有关键词
        if(entry_set.size()>0){
            node_highest = current;
            second_highest = current_pre;
            //cout<<"到"<<current<<",发现其有关键词！"<<endl;
            break;
        }
        else{

            //继续往上寻找关键词LCA
            //cout<<"到"<<current<<endl;
        }

    }
    //将扩大的查询范围子图加入
    for(int child: entry_set){
        if(child == current_pre) {
            //cout<<"跳过"<<child<<endl;
            //getchar();
            continue; //跳过已经遍历过的节点
        }
        //int son = Nodes[locid].gtreepath.back()(错误的！)
        //int node_pos = 0;
        int son = Nodes[locid].gtreepath[node_pos + 1];             //重要！容易出错
        // on gtreepath
        if (child == son) {
            double min_Dist =0;
            double score = getGSKScoreusrLeaf2O_Upper (a, alpha, min_Dist, Ukeys, GTree[child].objectUKeySet);

            if(score > score_bound){
                QEntry entry(child, false, node_pos+1, min_Dist,score);
                Queue.push(entry);
                //cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
            }

        }
            // brothers
        else if (GTree[child].father == GTree[son].father) {
            itm[child].clear();
            int allmin = -1;

            for (int j = 0; j < GTree[child].borders.size(); j++) {   //把到当前最优节点的各个child的距离算出，把child加入
                _min = -1;
                posa = GTree[child].up_pos[j];
                for (int k = 0; k < GTree[son].borders.size(); k++) {   //第k个border
                    posb = GTree[son].up_pos[k];
                    dis = itm[son][k] + GTree[current].mind[posb * GTree[current].union_borders.size() + posa];  //有问题
                    //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                    double dist = dis;
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
            double distance = allmin/1.0;



            double score = getGSKScoreusrLeaf2O_Upper(a, alpha, distance,Ukeys,
                                                      GTree[child].objectUKeySet);

            if(score > score_bound){
                QEntry entry(child, false, node_pos, allmin, score);
                Queue.push(entry);
                //cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
            }
        }
            // downstream (也按brother的方式来)
        else {
            itm[child].clear();
            int allmin = -1;

            for (int j = 0; j < GTree[child].borders.size(); j++) {
                _min = -1;
                posa = GTree[child].up_pos[j];
                for (int k = 0; k < GTree[second_highest].borders.size(); k++) {  //重要！ 容易出错！
                    posb = GTree[second_highest].current_pos[k];
                    //int right_end = second_highest;
                    dis = itm[second_highest][k] + GTree[current].mind[posb * GTree[current].union_borders.size() + posa];
                    //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
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
            double distance = allmin/1.0;
            double score = getGSKScoreusrLeaf2O_Upper(a, alpha, distance,Ukeys,
                                                      GTree[child].objectUKeySet);
            //getGSKScoreu2O_Upper(u_id, a, alpha, distance, GTree[child].objectUKeySet);

            if(score > score_bound){
                QEntry entry(child, false, node_pos, allmin, score);
                Queue.push(entry);
                //cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
            }

        }

    }
    //cout<<"update后的 topscore="<<Queue.top().score<<"score_k="<<score_max<<endl;
    if(Queue.size()>0){
        if(score_max > Queue.top().score || Queue.top().dis>EXPANSION_MAX){
            //priority_queue<QEntry> emptyQ;
            //swap(emptyQ,Queue);     // 把当前队列清空!
            eary_terminate = true;
        }
    }


}

//add expansion early termination , the value of the score is the priority
//内存操作

MultiLCLResult getLeafNode_LCL_bottom2up_multi_op_Memory(int leaf_id, vector<int> Ukeys,int K, float a, float alpha, double score_bound) {


    priority_queue<Result> resultFinal;
    priority_queue<QEntry> Queue;

    vector<priority_queue<LCLEntry>>  multi_Q;
    MultiLCLResult lcl_multi;
    set<int> stop_flag;


    int usr_leaf = leaf_id;
    for(int i=0;i<Ukeys.size();i++){
        priority_queue<LCLEntry> lcl_queue;
        multi_Q.push_back(lcl_queue);

        LCLResult lclresult;
        lcl_multi.push_back(lclresult);
    }
    priority_queue<LCLEntry>& lcl_q = multi_Q[0];

    //kingmura
    vector<vector<int>> qKey_multi;
    vector<double> optimal_text_score_list;
    for(int term: Ukeys){
        vector<int> qKey2;
        qKey2.push_back(term);
        qKey_multi.push_back(qKey2);

        double  optimal_text_score = textRelevance(qKey2, qKey2);

        optimal_text_score_list.push_back(optimal_text_score);
    }



    //获取leaf信息
    //选择usr leaf 节点内部最大距离
    float maxDisInsideLeaf = MaxDisBound[leaf_id][leaf_id];
    //选择该usr leaf 节点上的任意一个 border vertex进行 dj expansion

    int leaf_border = GTree[leaf_id].borders[0];


    //路网距离初始化
    int locid = leaf_border;
    double u_dis = maxDisInsideLeaf;
    vector<int> keywords = Ukeys;
    int node_highest_expansion;

    int u_leaf = Nodes[locid].gtreepath.back();

    // 计算 此 usr leaf 的border 到 他的所有祖先节点的border的距离
    unordered_map<int, vector<int> > itm; // intermediate answer, tree node -> array
    itm.clear();
    int cid, posa, posb, _min, dis;
    for (int i = Nodes[locid].gtreepath.size() - 1; i > 0; i--) {   //自底向上
        int tn = Nodes[locid].gtreepath[i];
        itm[tn].clear();
        //block_num++;

        if (GTree[tn].isleaf) {
            posa = lower_bound(GTree[tn].leafnodes.begin(), GTree[tn].leafnodes.end(), locid) -
                   GTree[tn].leafnodes.begin();

            for (int j = 0; j < GTree[tn].borders.size(); j++) {
                itm[tn].push_back(GTree[tn].mind[j * GTree[tn].leafnodes.size() + posa]);
            }
        } else {
            cid = Nodes[locid].gtreepath[i + 1];
            double overall_min=-1.0;
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
                itm[tn].push_back(_min);      //自底向上记录 loc 到其所在 叶节点，以及叶节点的祖先节点的(each border)最小距离
                if(overall_min == -1 || overall_min>_min)
                    overall_min= _min;

            }
            //record the highest node that within the expansion range
            if(overall_min /1.0> 2*EXPANSION_MAX){
                //if(overall_min > EXPANSION_MAX){
                node_highest_expansion = tn;
                break;
            }

        }

    }

#ifdef  TRACK
    //cout<<"完成计算 此 usr leaf 的border 到 他的所有祖先节点的border的距离"<<endl;
#endif
    //先把 u_leaf压入最大优先队列中
    int root_id = 0;
    int node_current = u_leaf;
    int node_highest = node_current;
    double score_max = -1;  double score_k = -1;
    vector<int> score_k_multi;
    for(int i=0;i<Ukeys.size();i++){
        score_k_multi.push_back(-1);
    }


    double min_Dist = maxDisInsideLeaf; int node_pos = Nodes[locid].gtreepath.size()-1;
    double score_leaf = getGSKScoreusrLeaf2O_Lower(a, alpha, min_Dist, Ukeys, GTree[u_leaf].objectUKeySet);
    QEntry entry(u_leaf,false,node_pos,0, score_leaf);
    Queue.push(entry);
#ifdef  TRACK
    //cout<<"score_leaf="<<score_leaf<<endl;
    // 自底向上进行节点遍历迭代
    //cout<<"开始自底向上进行节点遍历迭代"<<endl;
#endif

    bool eary_terminate = false;
    while(Queue.size() > 0 || node_highest != root_id){//node_highest != root_id){  //Queue.size() > 0 ||

        //当前范围内已全部遍历， 需更新扩大路网图的遍历范围
        if(Queue.empty()){
            if(node_highest==node_highest_expansion)
                eary_terminate = true;  //jordan
            //cout<<"当前范围内已全部遍历， 需更新扩大路网图的遍历范围"<<endl;
            //判断当前是否所有关键词的lcl表中都已有足够中间结果
            bool flag = false;
            for(int i = 0;i<Ukeys.size();i++){
                if (multi_Q[i].size()<K){
                    flag = false;
                    break;
                }
                else{
                    flag = true;
                }

            }
            if(flag && eary_terminate) break;
            UpdateR_usrLeaf_Memory(usr_leaf, locid, a,  alpha, score_bound, node_highest, node_pos, keywords, itm, Queue, score_max, eary_terminate);
        }


        //当前范围内遍历
        //cout<<"node_highest="<<node_highest<<endl;
        QEntry top = Queue.top();

        double topk_score = top.score;
        int current_dist = top.dis;
        //if(current_dist/WEIGHT_INFLATE_FACTOR>EXPANSION_MAX) break;

        if(score_k > topk_score) {   // 当前队列中已无更优元素
            //priority_queue<QEntry> emptyQ;
            //swap(emptyQ,Queue);     // 把当前队列清空!
            if (node_highest!= root_id) { //如果还未到最高层，那就再拓展一波
                UpdateR_usrLeaf_Memory(usr_leaf, locid, a,  alpha, score_bound, node_highest, node_pos, keywords, itm, Queue, score_max,
                                       eary_terminate);
            }
            else { //提前终止条件  //jordan,kingmura
                if(Queue.size()>0){
                    top = Queue.top();
                    bool flag = false;
                    for(int i=0;i<score_k_multi.size();i++){
                        if(score_k_multi[i] > top.score){   //提前终止条件！
                            if(top.isvertex){
                                //cout<<"提前终止条件，出队列，o="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                            }
                            else if(GTree[top.id].isleaf){
                                //cout<<"提前终止条件，出队列，leaf="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                            } else{
                                //cout<<"提前终止条件，出队列，n="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                            }
                            flag = true;
                            //getchar();
                        }
                        else{
                            flag = false;
                            break;
                        }
                    }
                    if(flag) {
                        cout<<"flag break"<<endl;
                        break; // break while loop
                    }

                }

            }
            if(eary_terminate) {
                cout<<"eary_ break"<<endl;
                break; // break while loop
            }


        }
        //当队列内还有更优的元素时
        if(Queue.size()>0){
            top = Queue.top();
            Queue.pop();
            if (top.isvertex) {  //当前最优为某个路网上的兴趣点（顶点）
                //cout<<"出队列，o="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;

                set<int> term_idx;
                term_idx.clear();
                int poi = top.id;
                //对兴趣点所覆盖的关键词进行检查
                for(int i=0;i<qKey_multi.size();i++){
                    double _score = textRelevance(qKey_multi[i], POIs[poi].keywords);
                    if(_score>0){
                        term_idx.insert(i);
                    }
                }

                for(int idx: term_idx) {
                    int term_id = qKey_multi[idx][0];
                    if(multi_Q[idx].size() < K){
                        int poi = top.id;
                        //Result tmpRlt(poi, POIs[poi].Ni, POIs[poi].Nj, POIs[poi].dis, POIs[poi].dist, top.score, POIs[poi].keywords);
                        //resultFinal.push(tmpRlt);

                        LCLEntry le(poi, top.score, 1);

                        multi_Q[idx].push(le);
                        //cout<<"simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;

                        if(multi_Q[idx].size() == K) {
                            score_k_multi[idx] = multi_Q[idx].top().score;
                            if(score_k_multi[idx]>score_max)
                                score_max = score_k_multi[idx];

                        }
                    }
                    else{  // update lcl[idx],  score_k[idx]
                        int poi = top.id;

                        LCLEntry le(poi, top.score, 1);
                        multi_Q[idx].push(le);
                        multi_Q[idx].pop();

                        score_k_multi[idx] = multi_Q[idx].top().score;
                        if(score_k_multi[idx]>score_max)
                            score_max = score_k_multi[idx];

                        //cout<<"更新score_k="<<score_k<<endl;
                    }

                }

            }
            else if (GTree[top.id].isleaf) {           //当前最优为某个叶节点
                //cout<<"出队列，leaf="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
                if (top.id == Nodes[locid].gtreepath[top.lca_pos]) {  //该叶节点为usr所在的叶节点
                    vector<int> cands;
                    vector<int> result;
                    set<int> object_set;
                    // 获得该叶节点下包含usr_leaf关键词的所有poi集合
                    for (int term: Ukeys) {
                        if (GTree[top.id].inverted_list_o[term].size() > 0) {
                            for (int poi: GTree[top.id].inverted_list_o[term]) {
                                object_set.insert(poi);
                            }
                        }
                    }
                    //获得user对这些poi的评分
                    vector<int> objects;
                    for (int poi : object_set) {
                        int Ni = POIs[poi].Ni;
                        int Nj = POIs[poi].Nj;
                        cands.push_back(Ni);
                        cands.push_back(Nj);
                        //objects.push_back(object_id);
                        objects.push_back(poi);
                        objects.push_back(poi);
                    }
                    result = dijkstra_candidate(locid, cands, Nodes);
                    //将叶节点内的各poi(顶点)插入
                    for (int i = 0; i < cands.size(); i += 2) {
                        int o_id = objects[i];

                        double distance = result[i] / 1000;
                        distance = distance + u_dis + POIs[o_id].dis;
                        double distance2 = result[i + 1] / 1000;
                        distance2 = distance2 + u_dis + POIs[o_id].dist - POIs[o_id].dis;
                        double distance_min = min(distance, distance2);
                        //double score = getGSKScoreu2O_Upper(u_id, a, distance, node_uKey_Map[top.id]);



                        double score = getGSKScoreusrLeaf2o_Upper(a, alpha, distance_min,Ukeys,POIs[o_id].keywords);
                        if(score>0.0){
                            //cout<<"score="<<score<<endl;
                        } else{
                            continue;
                        }

                        if (score > score_bound) {
                            QEntry entry(o_id, true, top.lca_pos, result[i], score);
                            Queue.push(entry);
                            //cout<<"加入o"<<o_id<<"distance="<<distance_min<<",score="<<score<<endl;

                        }

                    }
                } else {  // 是叶节点，但与loc不在同一子图
                    set<int> object_set;
                    // 获得该叶节点下包含用户关键词的所有poi
                    for (int term: Ukeys) {
                        if (GTree[top.id].inverted_list_o[term].size() > 0) {
                            for (int poi: GTree[top.id].inverted_list_o[term]) {
                                object_set.insert(poi);
                            }
                        }
                    }
                    for (int id: object_set) {
                        int o_id = id;

                        //else cout<<"发现其他"<<endl;

                        int Ni = POIs[o_id].Ni;
                        int Nj = POIs[o_id].Nj;
                        double score = 0; int allmin = -1; double distance;
                        for (int i = 0; i < GTree[top.id].leafinvlist.size(); i++) {    //获取当前最优叶节点的顶点倒排表
                            int posa = GTree[top.id].leafinvlist[i];
                            int vertex = GTree[top.id].leafnodes[posa];
                            if (vertex == Ni || vertex == Nj) {
                                allmin = -1;
                                for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                    //计算loc 经过 叶节点中第 k个 border 到达 top.id叶节点内第i个顶点的距离
                                    dis = itm[top.id][k] + GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa];
                                    //cout<<"dis="<<dis<<",itm[top.id][k]="<<itm[top.id][k]<<",边界间距离="<<GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa]<<endl;
                                    if (allmin == -1) {
                                        allmin = dis;
                                    } else {
                                        if (dis < allmin) {
                                            allmin = dis;   // 提取loc 经过所有叶节点top.id中的boder, 而到达第i个顶点的最短距离
                                        }
                                    }

                                }
                                //把所有的allmin
                                distance = allmin / 1000;
                                distance = distance + u_dis + POIs[o_id].dis;

                                score = max(score, getGSKScoreusrLeaf2o_Upper(a, alpha,
                                                                              distance,Ukeys,POIs[o_id].keywords));

                            }

                        }

                        if (score > score_bound) {
                            QEntry entry(o_id, true, top.lca_pos, distance, score);
                            Queue.push(entry);
                            //cout<<"加入o"<<o_id<<"distance="<<distance<<",score="<<score<<endl;

                        }


                    }

                }

            }
            else {    //当前最优为某个中间节点
                //cout<<"出队列，n="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
                set<int> entry_set;
                // 获得该中间节点下包含用户关键词的所有child
                for (int term: Ukeys) {
                    if (GTree[top.id].inverted_list_o[term].size() > 0) {
                        for (int child: GTree[top.id].inverted_list_o[term]) {
                            entry_set.insert(child);
                        }
                    }
                }
                //将各子条目加入优先队列
                for(int child: entry_set){

                    int son = Nodes[locid].gtreepath[top.lca_pos + 1];
                    // on gtreepath
                    if (child == son) {
                        double _min_Dist =0;
                        double score = getGSKScoreusrLeaf2O_Upper(a, alpha, _min_Dist, Ukeys, GTree[child].objectUKeySet);
                        //cout<<"child_score="<<score<<endl;


                        if(score > score_bound){
                            QEntry entry(child, false, top.lca_pos + 1, min_Dist, score);
                            Queue.push(entry);

                        }

                    }
                        // brothers
                    else if (GTree[child].father == GTree[son].father) {
                        itm[child].clear();
                        int allmin = -1;

                        for (int j = 0; j < GTree[child].borders.size(); j++) {   //把到当前最优节点的各个child的距离算出，把child加入
                            _min = -1;
                            posa = GTree[child].up_pos[j];
                            for (int k = 0; k < GTree[son].borders.size(); k++) {   //第k个border
                                posb = GTree[son].up_pos[k];
                                dis = itm[son][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                                double dist = dis;
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
                        double distance = allmin/1.0;
                        double score = getGSKScoreusrLeaf2O_Upper(a, alpha, distance, Ukeys, GTree[child].objectUKeySet);

                        if(score > score_bound){
                            QEntry entry(child, false, top.lca_pos, allmin, score);
                            Queue.push(entry);

                        }
                    }
                        // downstream
                    else {
                        itm[child].clear();
                        int allmin = -1;

                        for (int j = 0; j < GTree[child].borders.size(); j++) {
                            _min = -1;
                            posa = GTree[child].up_pos[j];
                            for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                posb = GTree[top.id].current_pos[k];
                                dis = itm[top.id][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
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
                        double distance = allmin/1.0;
                        double score = getGSKScoreusrLeaf2O_Upper(a, alpha, distance, Ukeys, GTree[child].objectUKeySet);

                        if(score > score_bound){
                            QEntry entry(child, false, top.lca_pos, allmin, score);
                            Queue.push(entry);
                        }

                    }

                }


            }


        }



    }//while

    for(int i=0;i<Ukeys.size();i++){
        LCLResult result;
        result.LCL = multi_Q[i];
        if(multi_Q[i].size()==K)
            result.topscore = multi_Q[i].top().score;
        else{
            result.topscore = 0;
            //cout<<"lcl for"<<Ukeys[i]<<"less than"<<K<<endl;
        }

        lcl_multi[i]= result;

    }

    return lcl_multi;
}

// 外存操作
void UpdateR_usrLeaf_Disk(int usr_leaf, int locid, int a, double alpha, double score_bound, int& node_highest,int& node_pos, vector<int> keywords, unordered_map<int, vector<int> >& itm, priority_queue<QEntry>& Queue, double& score_max, bool& eary_terminate){
#ifdef TRACKLEAF
    cout<<"表示当前范围内已全部遍历， 需扩大遍历子图范围"<<endl;
#endif
    int posa, posb, _min, dis;
    int root_id =0;
    int node_highest_pre = node_highest;
    int second_highest =0;
    //int node_highest = GTree[node_highest].father;
    int current = node_highest; int current_pre = current;
    set<int> entry_set;
    //
    vector<int> Ukeys = keywords;

    //cout<<"从"<<current_pre<<"向上寻找"<<endl;
    while(current != root_id){
        // 获得该中间节点下包含用户关键词的所有child
        current_pre = current;
        TreeNode oNode = getGIMTreeNodeData(current,OnlyO);

        //读取父节点
        current = GTree[current].father;


        node_pos--;   //最高层节点层数往上一层
        for (int term: keywords) {
            vector<int>  related_entry_o = getObjectTermRelatedEntry(term, current);
            int _size = related_entry_o.size();
            if (related_entry_o.size() > 0) {
                for (int child: related_entry_o) {

                    entry_set.insert(child);
#ifdef TRACKLEAF
                    cout<<"加入child o entry"<<child<<endl;
#endif
                }
            }
        }
        //父节点就有关键词
        if(entry_set.size()>0){
            node_highest = current;
            second_highest = current_pre;
            //cout<<"到"<<current<<",发现其有关键词！"<<endl;
            break;
        }
        else{

            //继续往上寻找关键词LCA
            //cout<<"到"<<current<<endl;
        }

    }
    //将扩大的查询范围子图加入
    for(int child: entry_set){
        if(child == current_pre) {
            //cout<<"跳过"<<child<<endl;
            //getchar();
            continue; //跳过已经遍历过的节点
        }
        TreeNode child_node = getGIMTreeNodeData(child,OnlyO);

        int son = Nodes[locid].gtreepath[node_pos + 1];             //重要！容易出错
        //TreeNode son_node = getGTreeNodeData(getGTreeNodeAddr(son));

        // on gtreepath
        if (child == son) {
            double min_Dist =0;

            //printSetElements(child_node.objectUKeySet);

            double score = getGSKScoreusrLeaf2O_Upper (a, alpha, min_Dist, Ukeys, child_node.objectUKeySet);

            if(score > score_bound){
                QEntry entry(child, false, node_pos+1, min_Dist,score);
                Queue.push(entry);
#ifdef TRACKLEAF
                cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
#endif
            }

        }
            // brothers
        else if (child_node.father == GTree[son].father) {


            itm[child].clear();
            int allmin = -1;

            for (int j = 0; j < child_node.borders.size(); j++) {   //把到当前最优节点的各个child的距离算出，把child加入
                _min = -1;
                posa = GTree[child].up_pos[j];
                for (int k = 0; k < GTree[son].borders.size(); k++) {   //第k个border
                    posb = GTree[son].up_pos[k];
                    dis = itm[son][k] + GTree[current].mind[posb * GTree[current].union_borders.size() + posa];
                    //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                    double dist = dis;
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
            double distance = allmin/1.0;



            double score = getGSKScoreusrLeaf2O_Upper(a, alpha, distance,Ukeys,
                                                      child_node.objectUKeySet);

            if(score > score_bound){
                QEntry entry(child, false, node_pos, allmin, score);
                Queue.push(entry);
#ifdef TRACKLEAF
                cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
#endif
            }
        }
            // downstream (也按brother的方式来)
        else {
            itm[child].clear();
            int allmin = -1;

            for (int j = 0; j < GTree[child].borders.size(); j++) {
                _min = -1;
                posa = GTree[child].up_pos[j];
                for (int k = 0; k < GTree[second_highest].borders.size(); k++) {  //重要！ 容易出错！
                    posb = GTree[second_highest].current_pos[k];
                    //int right_end = second_highest;
                    dis = itm[second_highest][k] + GTree[current].mind[posb * GTree[current].union_borders.size() + posa];
                    //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
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
            double distance = allmin/1.0;
            double score = getGSKScoreusrLeaf2O_Upper(a, alpha, distance,Ukeys,
                                                      child_node.objectUKeySet);
            //getGSKScoreu2O_Upper(u_id, a, alpha, distance, GTree[child].objectUKeySet);

            if(score > score_bound){
                QEntry entry(child, false, node_pos, allmin, score);
                Queue.push(entry);
                //cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
#ifdef TRACKLEAF
                cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
#endif
            }

        }

    }
    //cout<<"update后的 topscore="<<Queue.top().score<<"score_k="<<score_max<<endl;
    if(Queue.size()>0){
        if(score_max > Queue.top().score || Queue.top().dis>EXPANSION_MAX ){
            //priority_queue<QEntry> emptyQ;
            //swap(emptyQ,Queue);     // 把当前队列清空!
            eary_terminate = true;
        }
    }


}


//外存操作
MultiLCLResult getLeafNode_LCL_bottom2up_multi_op_Disk(int leaf_id, vector<int> Ukeys,int K, float a, float alpha, double score_bound) {


    priority_queue<Result> resultFinal;
    priority_queue<QEntry> Queue;

    vector<priority_queue<LCLEntry>>  multi_Q;
    MultiLCLResult lcl_multi;
    set<int> stop_flag;


    int usr_leaf = leaf_id;
    for(int i=0;i<Ukeys.size();i++){
        priority_queue<LCLEntry> lcl_queue;
        multi_Q.push_back(lcl_queue);

        LCLResult lclresult;
        lcl_multi.push_back(lclresult);
    }
    priority_queue<LCLEntry>& lcl_q = multi_Q[0];

    //kingmura
    vector<vector<int>> qKey_multi;
    vector<double> optimal_text_score_list;
    for(int term: Ukeys){
        vector<int> qKey2;
        qKey2.push_back(term);
        qKey_multi.push_back(qKey2);

        double  optimal_text_score = textRelevance(qKey2, qKey2);

        optimal_text_score_list.push_back(optimal_text_score);
    }



    //获取leaf信息
    //选择usr leaf 节点内部最大距离
    float maxDisInsideLeaf = MaxDisBound[leaf_id][leaf_id];
    //选择该usr leaf 节点上的任意一个 border vertex进行 dj expansion
    int leaf_border = GTree[leaf_id].borders[0];


    //路网距离初始化
    int locid = leaf_border;
    double u_dis = maxDisInsideLeaf;
    vector<int> keywords = Ukeys;
    int node_highest_expansion;

    int u_leaf = Nodes[locid].gtreepath.back();

    // 计算 此 usr leaf 的border 到 所有节点的border的距离
    unordered_map<int, vector<int> > itm; // intermediate answer, tree node -> array
    itm.clear();
    int cid, posa, posb, _min, dis;
    for (int i = Nodes[locid].gtreepath.size() - 1; i > 0; i--) {   //自底向上
        int tn = Nodes[locid].gtreepath[i];
        itm[tn].clear();
        //block_num++;

        if (GTree[tn].isleaf) {
            posa = lower_bound(GTree[tn].leafnodes.begin(), GTree[tn].leafnodes.end(), locid) -
                   GTree[tn].leafnodes.begin();

            for (int j = 0; j < GTree[tn].borders.size(); j++) {
                itm[tn].push_back(GTree[tn].mind[j * GTree[tn].leafnodes.size() + posa]);
            }
        } else {
            cid = Nodes[locid].gtreepath[i + 1];
            double overall_min=-1.0;
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
                itm[tn].push_back(_min);      //自底向上记录 loc 到其所在 叶节点，以及叶节点的祖先节点的(each border)最小距离
                if(overall_min == -1 || overall_min>_min)
                    overall_min= _min;

            }
            //record the highest node that within the expansion range
            double max_expansion = EXPANSION_MAX;
#ifdef LV
            max_expansion = EXPANSION_MAX;
#else
            max_expansion = 100*EXPANSION_MAX;  //大数据集上，由于关键词分布稀疏，需要拓展更广的范围
#endif
            if(overall_min /1.0> 2*max_expansion){
                //if(overall_min > EXPANSION_MAX){
                node_highest_expansion = tn;
                break;
            }

        }

    }

#ifdef  TRACKLEAF
    cout<<"完成计算此usr_leaf"<<leaf_id<<"的border到所有祖先节点的border的距离"<<endl;
#endif

    //先把 u_leaf压入最大优先队列中
    int root_id = 0;
    int node_current = u_leaf;
    int node_highest = node_current;
    double score_max = -1;  double score_k = -1;
    vector<int> score_k_multi;
    for(int i=0;i<Ukeys.size();i++){
        score_k_multi.push_back(-1);
    }


    double min_Dist = maxDisInsideLeaf; int node_pos = Nodes[locid].gtreepath.size()-1;
    TreeNode usr_leafNode = getGIMTreeNodeData(usr_leaf, OnlyO);
    double score_leaf = getGSKScoreusrLeaf2O_Lower(a, alpha, min_Dist, Ukeys, usr_leafNode.objectUKeySet);
    QEntry entry(u_leaf,false,node_pos,0, score_leaf);
    Queue.push(entry);
#ifdef TRACKLEAF
    // 自底向上进行节点遍历迭代
    cout<<"score_leaf="<<score_leaf<<endl;
    cout<<"开始自底向上进行节点遍历迭代"<<endl;
#endif
    bool eary_terminate = false;
    while(Queue.size() > 0 || node_highest != root_id){//node_highest != root_id){  //Queue.size() > 0 ||

        //当前范围内已全部遍历， 需更新扩大路网图的遍历范围
        if(Queue.empty()){
            if(node_highest==node_highest_expansion)
                eary_terminate = true;  //jordan
            //cout<<"当前范围内已全部遍历， 需更新扩大路网图的遍历范围"<<endl;
            //判断当前是否所有关键词的lcl表中都已有足够中间结果
            bool flag = false;
            for(int i = 0;i<Ukeys.size();i++){
                if (multi_Q[i].size()<K){
                    flag = false;
                    break;
                }
                else{
                    flag = true;
                }

            }
            if(flag && eary_terminate) break; //gangcai
            UpdateR_usrLeaf_Disk(usr_leaf, locid, a,  alpha, score_bound, node_highest, node_pos, keywords, itm, Queue, score_max, eary_terminate);
        }


        //当前范围内遍历
        //cout<<"node_highest="<<node_highest<<endl;
        QEntry top = Queue.top();

        double topk_score = top.score;
        int current_dist = top.dis;
        //if(current_dist/WEIGHT_INFLATE_FACTOR>EXPANSION_MAX) break;

        if(score_k > topk_score) {   // 当前队列中已无更优元素
            //priority_queue<QEntry> emptyQ;
            //swap(emptyQ,Queue);     // 把当前队列清空!
            if (node_highest!= root_id) { //如果还未到最高层，那就再拓展一波
                UpdateR_usrLeaf_Disk(usr_leaf, locid, a,  alpha, score_bound, node_highest, node_pos, keywords, itm, Queue, score_max,
                                     eary_terminate);
            }
            else { //提前终止条件  //jordan,kingmura
                if(Queue.size()>0){
                    top = Queue.top();
                    bool flag = false;
                    for(int i=0;i<score_k_multi.size();i++){
                        if(score_k_multi[i] > top.score){   //提前终止条件！
                            if(top.isvertex){
                                //cout<<"提前终止条件，出队列，o="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                            }
                            else if(GTree[top.id].isleaf){
                                //cout<<"提前终止条件，出队列，leaf="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                            } else{
                                //cout<<"提前终止条件，出队列，n="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                            }
                            flag = true;
                            //getchar();
                        }
                        else{
                            flag = false;
                            break;
                        }
                    }
                    if(flag) {
                        cout<<"flag break"<<endl;
                        break; // break while loop
                    }

                }

            }
            if(eary_terminate) {
                cout<<"eary_ break"<<endl;
                break; // break while loop
            }


        }
        //当队列内还有更优的元素时
        if(Queue.size()>0){
            top = Queue.top();
            Queue.pop();
            if (top.isvertex) {  //当前最优为某个路网上的兴趣点（顶点）
#ifdef TRACKLEAF
                cout<<"出队列，o="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
#endif
                set<int> term_idx;
                term_idx.clear();
                int poi = top.id;
                //access disk
                POI p = getPOIFromO2UOrgLeafData(poi);
                //对兴趣点所覆盖的关键词进行检查
                for(int i=0;i<qKey_multi.size();i++){
                    double _score = textRelevance(qKey_multi[i], p.keywords);
                    if(_score>0){
                        term_idx.insert(i);
                    }
                }

                for(int idx: term_idx) {
                    int term_id = qKey_multi[idx][0];
                    if(multi_Q[idx].size() < K){
                        int poi = top.id;

                        LCLEntry le(poi, top.score, 1);

                        multi_Q[idx].push(le);
                        //cout<<"simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;

                        if(multi_Q[idx].size() == K) {
                            score_k_multi[idx] = multi_Q[idx].top().score;
                            if(score_k_multi[idx]>score_max)
                                score_max = score_k_multi[idx];

                        }
                    }
                    else{  // update lcl[idx],  score_k[idx]
                        int poi = top.id;

                        LCLEntry le(poi, top.score, 1);
                        multi_Q[idx].push(le);
                        multi_Q[idx].pop();

                        score_k_multi[idx] = multi_Q[idx].top().score;
                        if(score_k_multi[idx]>score_max)
                            score_max = score_k_multi[idx];

                        //cout<<"更新score_k="<<score_k<<endl;
                    }

                }

            }
            else if (GTree[top.id].isleaf) {           //当前最优为某个叶节点
#ifdef  TRACKLEAF
                cout<<"出队列，leaf="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
#endif
                if (top.id == Nodes[locid].gtreepath[top.lca_pos]) {  //该object 叶节点为usr所在的叶节点
                    vector<int> cands;
                    vector<int> result;
                    set<int> object_set;
                    // 获得该叶节点下包含usr_leaf关键词的所有poi集合
                    for (int term: Ukeys) {
                        vector<int> olist = getObjectTermRelatedEntry(term,top.id);
                        if (olist.size() > 0) {
                            for (int poi: olist) {
                                object_set.insert(poi);
                            }
                        }
                    }
                    //获得user对这些poi的评分
                    vector<int> objects; vector<POI> poiData_list;
                    for (int poi : object_set) {
                        POI p = getPOIFromO2UOrgLeafData(poi);

                        //printPOIInfo(p);

                        poiData_list.push_back(p);
                        int Ni = p.Ni;
                        int Nj = p.Nj;
                        cands.push_back(Ni);
                        cands.push_back(Nj);
                        //objects.push_back(object_id);
                        objects.push_back(p.id);
                        objects.push_back(p.id);
                    }
                    result = dijkstra_candidate(locid, cands, Nodes);
                    //将叶节点内的各poi(顶点)插入
                    for (int i = 0; i < cands.size(); i += 2) {
                        int o_id = objects[i];
                        POI p = poiData_list[i/2];//getPOIFromO2UOrgLeafData(o_id);

                        double distance = result[i] / 1000;
                        distance = distance + u_dis + p.dis;
                        double distance2 = result[i + 1] / 1000;
                        distance2 = distance2 + u_dis + p.dist - p.dis;
                        double distance_min = min(distance, distance2);
                        //double score = getGSKScoreu2O_Upper(u_id, a, distance, node_uKey_Map[top.id]);

                        double score = getGSKScoreusrLeaf2o_Upper(a, alpha, distance_min,Ukeys,p.keywords);

                        if(score>0.0){
                            //cout<<"score="<<score<<endl;
                        } else{
                            continue;
                        }

                        if (score > score_bound) {
                            QEntry entry(o_id, true, top.lca_pos, result[i], score);
                            Queue.push(entry);
                            //cout<<"加入o"<<o_id<<"distance="<<distance_min<<",score="<<score<<endl;

                        }

                    }
                } else {  // 是叶节点，但与loc不在同一子图
                    set<int> object_set;
                    // 获得该叶节点下包含用户关键词的所有poi
                    for (int term: Ukeys) {
                        vector<int> olist = getObjectTermRelatedEntry(term,top.id);
                        if (olist.size() > 0) {
                            for (int poi: olist) {
                                object_set.insert(poi);
                            }
                        }
                    }
                    for (int id: object_set) {
                        int o_id = id;

                        POI p = getPOIFromO2UOrgLeafData(o_id);

                        int Ni = p.Ni;
                        int Nj = p.Nj;
                        double score = 0; int allmin = -1; double distance;
                        for (int i = 0; i < GTree[top.id].leafinvlist.size(); i++) {    //获取当前最优叶节点的顶点倒排表
                            int posa = GTree[top.id].leafinvlist[i];
                            int vertex = GTree[top.id].leafnodes[posa];
                            if (vertex == Ni || vertex == Nj) {
                                allmin = -1;
                                for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                    //计算loc 经过 叶节点中第 k个 border 到达 top.id叶节点内第i个顶点的距离
                                    dis = itm[top.id][k] + GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa];
                                    //cout<<"dis="<<dis<<",itm[top.id][k]="<<itm[top.id][k]<<",边界间距离="<<GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa]<<endl;
                                    if (allmin == -1) {
                                        allmin = dis;
                                    } else {
                                        if (dis < allmin) {
                                            allmin = dis;   // 提取loc 经过所有叶节点top.id中的boder, 而到达第i个顶点的最短距离
                                        }
                                    }

                                }
                                //把所有的allmin
                                distance = allmin / 1000;
                                distance = distance + u_dis + p.dis;

                                score = max(score, getGSKScoreusrLeaf2o_Upper(a, alpha,
                                                                              distance,Ukeys,p.keywords));

                            }

                        }

                        if (score > score_bound) {
                            QEntry entry(o_id, true, top.lca_pos, distance, score);
                            Queue.push(entry);
                            //cout<<"加入o"<<o_id<<"distance="<<distance<<",score="<<score<<endl;

                        }


                    }

                }

            }
            else {    //当前最优为某个中间节点
#ifdef  TRACKLEAF
                cout<<"出队列，n="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
#endif
                set<int> entry_set;
                // 获得该中间节点下包含用户关键词的所有child
                for (int term: Ukeys) {
                    vector<int> child_list = getObjectTermRelatedEntry(term,top.id);
                    if (child_list.size() > 0) {
                        for (int child: child_list) {
                            entry_set.insert(child);
                        }
                    }
                }
                //将各子条目加入优先队列
                for(int child: entry_set){

                    int son = Nodes[locid].gtreepath[top.lca_pos + 1];
                    TreeNode child_node = getGIMTreeNodeData(child, OnlyO);
                    // on gtreepath
                    if (child == son) {
                        double _min_Dist =0;
                        double score = getGSKScoreusrLeaf2O_Upper(a, alpha, _min_Dist, Ukeys, child_node.objectUKeySet);
                        //cout<<"child_score="<<score<<endl;


                        if(score > score_bound){
                            QEntry entry(child, false, top.lca_pos + 1, min_Dist, score);
                            Queue.push(entry);

                        }

                    }
                        // brothers
                    else if (GTree[child].father == GTree[son].father) {
                        itm[child].clear();
                        int allmin = -1;

                        for (int j = 0; j < GTree[child].borders.size(); j++) {   //把到当前最优节点的各个child的距离算出，把child加入
                            _min = -1;
                            posa = GTree[child].up_pos[j];
                            for (int k = 0; k < GTree[son].borders.size(); k++) {   //第k个border
                                posb = GTree[son].up_pos[k];
                                dis = itm[son][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                                double dist = dis;
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
                        double distance = allmin/WEIGHT_INFLATE_FACTOR;
                        double score = getGSKScoreusrLeaf2O_Upper(a, alpha, distance, Ukeys, child_node.objectUKeySet);

                        if(score > score_bound){
                            QEntry entry(child, false, top.lca_pos, allmin, score);
                            Queue.push(entry);

                        }
                    }
                        // downstream
                    else {
                        itm[child].clear();
                        int allmin = -1;

                        for (int j = 0; j < GTree[child].borders.size(); j++) {
                            _min = -1;
                            posa = GTree[child].up_pos[j];
                            for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                posb = GTree[top.id].current_pos[k];
                                dis = itm[top.id][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
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
                        double distance = allmin/WEIGHT_INFLATE_FACTOR;
                        double score = getGSKScoreusrLeaf2O_Upper(a, alpha, distance, Ukeys, child_node.objectUKeySet);

                        if(score > score_bound){
                            QEntry entry(child, false, top.lca_pos, allmin, score);
                            Queue.push(entry);
                        }

                    }

                }


            }


        }



    }//while

    for(int i=0;i<Ukeys.size();i++){
        LCLResult result;
        result.LCL = multi_Q[i];
        if(multi_Q[i].size()==K)
            result.topscore = multi_Q[i].top().score;
        else{
            result.topscore = 0;
            //cout<<"lcl for"<<Ukeys[i]<<"less than"<<K<<endl;
        }

        lcl_multi[i]= result;

    }

    return lcl_multi;
}


MultiLCLResult getLeafNode_LCL_bottom2up_multi(int leaf_id, vector<int> Ukeys,int K, float a, float alpha, double score_bound){
#ifndef DiskAccess
    return getLeafNode_LCL_bottom2up_multi_op_Memory(leaf_id, Ukeys,K, a,alpha, score_bound);
#else
    return getLeafNode_LCL_bottom2up_multi_op_Disk(leaf_id, Ukeys,K, a,alpha, score_bound);
#endif
}


//usr继承其所属的leaf点的LCL结果，并对自身的LCL进行求解
double get_usrLCL_light(int usr_id, LCLResult leaf_lcl, int Qk, float a, double alpha,double gsk){
    if(Users[usr_id].u_lcl !=-1){
        //cout<<"直接返回"<<endl;
        return Users[usr_id].u_lcl;
    }
    //cout<<"需要计算"<<endl;
    priority_queue<LCLEntry> LCL = leaf_lcl.LCL;
    //priority_queue<UCLEntry> UCL = node_filter.UCL;
    priority_queue<ResultCurrent> u_lcl;
    double current_score = 0.0;
    //cout<<LCL.size()<<endl;
    //getchar();
    while(LCL.size() > 0){
        LCLEntry le = LCL.top();
        LCL.pop();
        if(le.score< current_score) break;
        int o_id = le.id;
        double score = getGSKScoreo2u_Lower(a, alpha, o_id, usr_id);
        if(u_lcl.size()<Qk){
            u_lcl.push(ResultCurrent(o_id,score));
        }
        else{
            u_lcl.push(ResultCurrent(o_id,score));
            u_lcl.pop();
            ResultCurrent rc = u_lcl.top();
            current_score = rc.score;
            if(current_score > gsk){
                //cout<<"提前return:"<<current_score<<endl;
                return current_score;
            }

        }
    }
    //getchar();
    current_score = u_lcl.top().score;
    Users[usr_id].u_lcl = current_score;
    //cout<<"return:"<<current_score<<endl;
    //getchar();
    return current_score;
}


double get_usrLCL_light(int usr_id, LCLResult leaf_lcl, map<int,set<int>> node_poiCnt, int Qk, float a, double alpha,double gsk){
    if(Users[usr_id].u_lcl !=-1){
        //cout<<"直接返回"<<endl;
        return Users[usr_id].u_lcl;
    }
    //cout<<"需要计算"<<endl;
    priority_queue<LCLEntry> LCL = leaf_lcl.LCL;
    //priority_queue<UCLEntry> UCL = node_filter.UCL;
    priority_queue<ResultCurrent> u_lcl;
    double current_score = 0.0;
    //cout<<LCL.size()<<endl;
    while(LCL.size() > 0){
        LCLEntry le = LCL.top();
        LCL.pop();
        if(le.score< current_score) break;
        int o_id = le.id;
        double score = getGSKScoreo2u_Lower(a, alpha, o_id, usr_id);
        if(u_lcl.size()<Qk){
            u_lcl.push(ResultCurrent(o_id,score));
        }
        else{
            u_lcl.push(ResultCurrent(o_id,score));
            u_lcl.pop();
            ResultCurrent rc = u_lcl.top();
            current_score = rc.score;
            if(current_score > gsk){
                //cout<<"提前return:"<<current_score<<endl;
                return current_score;
            }

        }
    }
    //getchar();
    current_score = u_lcl.top().score;
    Users[usr_id].u_lcl = current_score;

    return current_score;
}

priority_queue<ResultCurrent> get_usrLCL_update(int usr_id, LCLResult leaf_lcl, map<int,set<int>> node_poiCnt, int Qk, float a, double alpha,double gsk){
    if(Users[usr_id].u_LCL.size()>0){
        //cout<<"直接返回"<<endl;
        return Users[usr_id].u_LCL;
    }
    //cout<<"需要计算"<<endl;
    priority_queue<LCLEntry> LCL = leaf_lcl.LCL;
    //priority_queue<UCLEntry> UCL = node_filter.UCL;
    priority_queue<ResultCurrent> u_lcl;
    double current_score = 0.0;
    //cout<<LCL.size()<<endl;
    while(LCL.size() > 0){
        LCLEntry le = LCL.top();
        LCL.pop();
        if(le.score< current_score) break;
        int o_id = le.id;
        double score = getGSKScoreo2u_Lower(a,alpha, o_id, usr_id);
        if(u_lcl.size()<Qk){
            u_lcl.push(ResultCurrent(o_id,score));
        }
        else{
            u_lcl.push(ResultCurrent(o_id,score));
            u_lcl.pop();
            ResultCurrent rc = u_lcl.top();
            current_score = rc.score;
            if(current_score > gsk){
                //cout<<"提前return:"<<current_score<<endl;
                return u_lcl;//current_score;
            }

        }
    }
    //getchar();
    current_score = u_lcl.top().score;
    Users[usr_id].u_lcl = current_score;
    //cout<<"return:"<<current_score<<endl;
    //getchar();
    return u_lcl;//current_score;
}


priority_queue<ResultCurrent> get_usrLCL_update_Plus(int usr_id, set<int>& update_o_set, int Qk, float a, double alpha,double gsk){
    if(Users[usr_id].u_LCL.size()>0){
        //cout<<"直接返回"<<endl;
        return Users[usr_id].u_LCL;
    }

    priority_queue<ResultCurrent> u_lcl;
    double current_score = 0.0;
    //cout<<LCL.size()<<endl;
    int _size;

    //cout<<"---------------u"<<usr_id<<"‘s keyword:";
    //printElements(Users[usr_id].keywords);
    //cout<<"update_o_set size="<<update_o_set.size()<<endl;
    for(int o_id:update_o_set){
        //cout<<"o"<<o_id<<":";
        //printElements(POIs[o_id].keywords);

        double score = getGSKScoreo2u_Lower(a,alpha, o_id, usr_id);
        //cout<<"score="<<score<<endl;
        if(score==0){
            //getchar();
            continue;
        }

        if(u_lcl.size()<Qk){
            u_lcl.push(ResultCurrent(o_id,score));
            _size = u_lcl.size();
        }
        else{
            u_lcl.push(ResultCurrent(o_id,score));
            u_lcl.pop();
            ResultCurrent rc = u_lcl.top();
            current_score = rc.score;


        }
    }
    //getchar();
    //current_score = u_lcl.top().score;

    Users[usr_id].u_lcl = current_score;
    //cout<<"lcl size="<<u_lcl.size();
    if(u_lcl.size()<Qk)
        current_score = 0;
    else
        current_score = u_lcl.top().score;
    //cout<<" ,lcl score:"<<current_score<<endl;

    return u_lcl;//current_score;
}







/*-------------------------------用GTree kNN算法求用户的top-k结果-----------------------------*/



vector<TopResult> TkGSKQ_process(int u_id, int K, float a, float alpha) {
    // init priority queue & result set
    //vector<Status_query> pq;
    //pq.clear();
    //vector<Topk_entry> Q;

    vector<TopResult> rstset;
    rstset.clear();
    priority_queue<QEntry> Queue;

    //获取user信息
    int locid = Users[u_id].Ni;
    double u_dis = Users[u_id].dis;
    vector<int> keywords = Users[u_id].keywords;
    //获取与用户查询关键词相关的兴趣点所在的leaf node
    set<int> usr_related_poiLeaf;
    set<int> term_Leaf;
    map<int,set<int>> occurence_list;
    map<int,set<int>> node_uKey_Map;
    for(int term_id:keywords){
        set<int> term_Leaf = poiTerm_leafSet[term_id];
        if(term_Leaf.size()>0){
            for(int leaf: term_Leaf){
                usr_related_poiLeaf.insert(leaf);
                node_uKey_Map[leaf].insert(term_id); // 把该关键词加入node_uKey_Map中
                for (int poi: leafPoiInv[term_id][leaf]) {
                    //构建 leaf 的 occurence_list
                    occurence_list[leaf].insert(poi);

                }
                if(true){  //（当与u.LCL对接时候， 记得变为false-重要！）
                    int father;
                    int current = leaf;
                    while (current != 0) {
                        father = GTree[current].father;
                        occurence_list[father].insert(current);
                        node_uKey_Map[father].insert(term_id); // 把该关键词加入node_uKey_Map中
                        current = father;
                    }
                }
            }
        }
    }
    //cout<<"usr_related_poiLeaf size="<<usr_related_poiLeaf.size()<<endl;

    /*-----------继续构建 occurence_list-----------*/
    //构建每个leaf的祖先节点的occurrence list  （当与u.LCL对接时候， 记得变为true-重要！）
    if(false){
        for(int leaf: usr_related_poiLeaf){
            //将leaf加入他祖先节点的 occurence_list
            int father;
            int current = leaf;
            while (current != 0) {
                father = GTree[current].father;
                occurence_list[father].insert(current);
                current = father;
            }
        }
    }
    //return rstset;


    /*-------------根据usr_related_poiLeaf构建 occurrence list----------------*/

    // 计算 loc 到 所有节点的border的距离
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

    //return rstset;

    // 自顶向下进行遍历
    int root_id = 0; double min_Dist = 0; int node_pos =0;
    double score_root = getGSKScoreu2O_Upper(u_id, a, alpha, min_Dist, node_uKey_Map[root_id]);
    //cout<<"score_root="<<score_root<<endl;

    //Status_query rootstatus = {root_id, false, node_pos, min_Dist, score_root};   //id; isvertex; lca_pos; dis; score;
    //pq.push_back(rootstatus);
    //make_heap(pq.begin(), pq.end(), Status_query_max_comp());
    QEntry entry(root_id,false,node_pos,min_Dist,score_root);
    Queue.push(entry);


    vector<int> cands;
    vector<int> result;
    int child, son, allmin, vertex;

    while (Queue.size() > 0 && rstset.size() < K) {  //&& rstset.size() < K
        //Status_query top = pq[0];
        //pop_heap(pq.begin(), pq.end(), Status_query_max_comp());
        //pq.pop_back();
        QEntry top = Queue.top();
        Queue.pop();

        if (top.isvertex) {  //当前最优为某个路网上的兴趣点（顶点）
            TopResult rs = {top.id, top.score};
            rstset.push_back(rs);
        } else {           //当前最优为某个叶节点
            if (GTree[top.id].isleaf) {
                // inner of leaf node, do dijkstra, 该叶节点为usr所在的叶节点
                if (top.id == Nodes[locid].gtreepath[top.lca_pos]) {
                    int node_id = top.id;
                    int u_leaf = Nodes[locid].gtreepath.back();
                    cands.clear();
                    //for (int i = 0; i < GTree[top.id].leafinvlist.size(); i++) {
                    //cands.push_back(GTree[top.id].leafnodes[GTree[top.id].leafinvlist[i]]);
                    //}
                    vector<int> objects;
                    for (int poi : occurence_list[top.id]) {

                        int Ni = POIs[poi].Ni;
                        int Nj = POIs[poi].Nj;
                        cands.push_back(Ni);
                        cands.push_back(Nj);
                        //objects.push_back(object_id);
                        objects.push_back(poi);
                        objects.push_back(poi);
                    }
                    result = dijkstra_candidate(locid, cands, Nodes);
                    //将叶节点内的各poi(顶点)插入
                    for (int i = 0; i < cands.size(); i+=2) {
                        int o_id = objects[i];
                        //if(o_id == 4128){
                        //getchar();
                        //}
                        double distance = result[i]/1000;
                        distance = distance + u_dis + POIs[o_id].dis;
                        double distance2 = result[i+1]/1000;
                        distance2 = distance + u_dis + POIs[o_id].dist-POIs[o_id].dis;
                        double distance_min = min(distance,distance2);
                        //double score = getGSKScoreu2O_Upper(u_id, a, distance, node_uKey_Map[top.id]);
                        double score = getGSKScoreu2o(u_id, o_id,a, alpha, distance_min);
                        //cout<<"o"<<o_id<<",(Ni)vertex"<<cands[i]<<" score = "<<score<<endl;
                        //Status_query status = {cands[i], true, top.lca_pos, result[i],score};
                        //Status_query status = {o_id, true, top.lca_pos, result[i],score};
                        //pq.push_back(status);
                        //push_heap(pq.begin(), pq.end(), Status_query_max_comp());
                        QEntry entry(o_id, true, top.lca_pos, result[i],score);
                        Queue.push(entry);
                        if(5721 == o_id){
                            cout<<"发现o5721在同一叶节点,distance"<<distance_min<<endl;
                        }
                        //else{cout<<"发现其他"<<endl;}

                    }


                }
                    // else do
                else {          //是叶节点，但与loc不在同一子图
                    for(int id :occurence_list[top.id]){
                        int o_id = id;

                        int Ni = POIs[o_id].Ni;
                        int Nj = POIs[o_id].Nj;
                        if(5721 == o_id){
                            cout<<"发现o5721"<<endl;
                        }
                        else{cout<<"发现其他"<<endl;}
                        double score = 0; double distance =0;
                        for (int i = 0; i < GTree[top.id].leafinvlist.size(); i++) {    //获取当前最优叶节点的顶点倒排表
                            posa = GTree[top.id].leafinvlist[i];
                            vertex = GTree[top.id].leafnodes[posa];
                            if(vertex == Ni || vertex == Nj){
                                allmin = -1;

                                for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                    //计算loc 经过 叶节点中第 k个 border 到达 top.id叶节点内第i个顶点的距离
                                    dis = itm[top.id][k] + GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa];
                                    cout<<"dis="<<dis<<",itm[top.id][k]="<<itm[top.id][k]<<",边界间距离="<<GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa]<<endl;
                                    if (allmin == -1) {
                                        allmin = dis;
                                    } else {
                                        if (dis < allmin) {
                                            allmin = dis;   // 提取loc 经过所有叶节点top.id中的boder, 而到达第i个顶点的最短距离
                                        }
                                    }

                                }
                                //把所有的allmin
                                distance = allmin /1000;
                                distance = distance + u_dis + POIs[o_id].dis;

                                double score_max = max(score,getGSKScoreu2o(u_id, o_id,a, alpha, distance));
                                score = score_max;

                                //double score = getGSKScoreu2O_Upper(u_id, a, allmin, node_uKey_Map[top.id]);
                                //cout<<"vertex_score="<<score<<endl;

                            }

                        }
                        if(score>0){
                            //Status_query status = {o_id, true, top.lca_pos, allmin, score};
                            //pq.push_back(status);
                            //push_heap(pq.begin(), pq.end(), Status_query_max_comp());
                            QEntry entry(o_id, true, top.lca_pos, allmin, score);
                            Queue.push(entry);
                            cout<<"o_"<<o_id<<"最终score = "<<score<<",distance = "<<distance<<endl;
                            /*if(1682==o_id){
                                cout<<"发现1682，score="<<entry.score<<endl;
                                cout<<"leaf"<<top.id<<",score="<<top.score<<endl;
                                getchar();
                            }*/
                        }


                    }
                }
            } else if(occurence_list[top.id].size()>0){    //为非叶节点,且occurence_list不为空
                //for (int i = 0; i < GTree[top.id].nonleafinvlist.size(); i++) {
                //child = GTree[top.id].nonleafinvlist[i];
                for(int child: occurence_list[top.id]){

                    son = Nodes[locid].gtreepath[top.lca_pos + 1];
                    // on gtreepath
                    if (child == son) {
                        double score = getGSKScoreu2O_Upper(u_id, a, alpha, min_Dist, node_uKey_Map[child]);
                        //cout<<"child_score="<<score<<endl;
                        double min_Dist =0;
                        //Status_query status = {child, false, top.lca_pos + 1, min_Dist,score};
                        //pq.push_back(status);
                        //push_heap(pq.begin(), pq.end(), Status_query_max_comp());
                        QEntry entry(child, false, top.lca_pos + 1, min_Dist,score);
                        Queue.push(entry);
                    }
                        // brothers
                    else if (GTree[child].father == GTree[son].father) {
                        itm[child].clear();
                        allmin = -1;

                        for (int j = 0; j < GTree[child].borders.size(); j++) {   //把到当前最优节点的各个child的距离算出，把child加入
                            _min = -1;
                            posa = GTree[child].up_pos[j];
                            for (int k = 0; k < GTree[son].borders.size(); k++) {   //第k个border
                                posb = GTree[son].up_pos[k];
                                dis = itm[son][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                                double dist = dis;
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
                        double distance = allmin/WEIGHT_INFLATE_FACTOR;
                        double score = getGSKScoreu2O_Upper(u_id, a, alpha, distance, node_uKey_Map[child]);
                        //cout<<"child_score(brother)="<<score<<endl;
                        //Status_query status = {child, false, top.lca_pos, allmin, score};
                        //pq.push_back(status);
                        //push_heap(pq.begin(), pq.end(), Status_query_max_comp());
                        QEntry entry(child, false, top.lca_pos, allmin, score);
                        Queue.push(entry);
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
                                dis = itm[top.id][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
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
                        double distance = allmin/WEIGHT_INFLATE_FACTOR;
                        double score = getGSKScoreu2O_Upper(u_id, a, alpha, distance, node_uKey_Map[child]);
                        //cout<<"child_score(downstream)="<<score<<endl;
                        //Status_query status = {child, false, top.lca_pos, allmin, score};
                        //pq.push_back(status);
                        //push_heap(pq.begin(), pq.end(), Status_query_max_comp());
                        QEntry entry(child, false, top.lca_pos, allmin, score);
                        Queue.push(entry);
                    }
                }
            }

        }
    }

    //vector<int> rst;
    for (int i = rstset.size()-1; i < rstset.size(); i-- ){
        cout<<"top"<<(i+1)<<":o"<<rstset[i].id<<",score="<<rstset[i].score<<endl;
    }
    if(rstset[K-1].score> Users[u_id].topkScore_current){
        Users[u_id].topkScore_current = rstset[K-1].score;  //记得cache user 的 topk结果
        Users[u_id].topkScore_Final = rstset[K-1].score;
    }

    return rstset;
}

vector<TopResult> TkGSKQ_top2down_verify(int u_id, int K, float a, float alpha, double score_bound) {
    // init priority queue & result set
    //vector<Status_query> pq;
    //pq.clear();
    //vector<Topk_entry> Q;

    vector<TopResult> rstset;
    rstset.clear();
    priority_queue<QEntry> Queue;
    int count =0;
    //获取user信息
    int locid = Users[u_id].Ni;
    double u_dis = Users[u_id].dis;
    vector<int> keywords = Users[u_id].keywords;
    //获取与用户查询关键词相关的兴趣点所在的leaf node
    set<int> usr_related_poiLeaf;
    set<int> term_Leaf;
    map<int,set<int>> occurence_list;
    map<int,set<int>> node_uKey_Map;
    for(int term_id:keywords){
        set<int> term_Leaf = poiTerm_leafSet[term_id];
        count++;
        if(term_Leaf.size()>0){
            for(int leaf: term_Leaf){
                usr_related_poiLeaf.insert(leaf);
                node_uKey_Map[leaf].insert(term_id); // 把该关键词加入node_uKey_Map中
                for (int poi: leafPoiInv[term_id][leaf]) {
                    //构建 leaf 的 occurence_list
                    occurence_list[leaf].insert(poi);
                    count++;
                    //block_num++;

                }
                if(true){  //（当与u.LCL对接时候， 记得变为false-重要！）
                    int father;
                    int current = leaf;
                    while (current != 0) {
                        father = GTree[current].father;
                        occurence_list[father].insert(current);
                        node_uKey_Map[father].insert(term_id); // 把该关键词加入node_uKey_Map中
                        current = father;
                        count++;
                    }
                }
            }
        }
    }
    //cout<<"usr_related_poiLeaf size="<<usr_related_poiLeaf.size()<<endl;



    /*-------------根据usr_related_poiLeaf构建 occurrence list----------------*/

    // 计算 loc 到 所有节点的border的距离
    unordered_map<int, vector<int> > itm; // intermediate answer, tree node -> array
    itm.clear();
    int tn, cid, posa, posb, _min, dis;
    for (int i = Nodes[locid].gtreepath.size() - 1; i > 0; i--) {   //自底向上
        tn = Nodes[locid].gtreepath[i];
        itm[tn].clear();
        count++;
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
                //count++;
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
                    //count++;
                }
                // update
                itm[tn].push_back(_min);      //自底向上记录 loc 到其所在 叶节点，以及叶节点的祖先节点的最小距离
            }
        }

    }

    //return rstset;

    // 自顶向下进行遍历
    int root_id = 0; double min_Dist = 0; int node_pos =0;
    double score_root = getGSKScoreu2O_Upper(u_id, a, alpha, min_Dist, node_uKey_Map[root_id]);
    //cout<<"score_root="<<score_root<<endl;

    //Status_query rootstatus = {root_id, false, node_pos, min_Dist, score_root};   //id; isvertex; lca_pos; dis; score;
    //pq.push_back(rootstatus);
    //make_heap(pq.begin(), pq.end(), Status_query_max_comp());
    QEntry entry(root_id,false,node_pos,min_Dist,score_root);
    Queue.push(entry);


    vector<int> cands;
    vector<int> result;
    int child, son, allmin, vertex;

    while (Queue.size() > 0 && rstset.size() < K) {  //rstset.size() < K
        //Status_query top = pq[0];
        //pop_heap(pq.begin(), pq.end(), Status_query_max_comp());
        //pq.pop_back();
        QEntry top = Queue.top();
        Queue.pop();

        if (top.isvertex) {  //当前最优为某个路网上的兴趣点（顶点）
            TopResult rs = {top.id, top.score};
            rstset.push_back(rs);
        } else {           //当前最优为某个叶节点
            if (GTree[top.id].isleaf) {
                // inner of leaf node, do dijkstra, 该叶节点为usr所在的叶节点
                if (top.id == Nodes[locid].gtreepath[top.lca_pos]) {
                    int node_id = top.id;
                    int u_leaf = Nodes[locid].gtreepath.back();
                    cands.clear();
                    //for (int i = 0; i < GTree[top.id].leafinvlist.size(); i++) {
                    //cands.push_back(GTree[top.id].leafnodes[GTree[top.id].leafinvlist[i]]);
                    //}
                    vector<int> objects;
                    for (int poi : occurence_list[top.id]) {

                        int Ni = POIs[poi].Ni;
                        int Nj = POIs[poi].Nj;
                        cands.push_back(Ni);
                        cands.push_back(Nj);
                        //objects.push_back(object_id);
                        objects.push_back(poi);
                        objects.push_back(poi);
                    }
                    result = dijkstra_candidate(locid, cands, Nodes);
                    //将叶节点内的各poi(顶点)插入
                    for (int i = 0; i < cands.size(); i+=2) {
                        int o_id = objects[i];
                        //if(o_id == 4128){
                        //getchar();
                        //}
                        double distance = result[i]/1000;
                        distance = distance + u_dis + POIs[o_id].dis;
                        double distance2 = result[i+1]/1000;
                        distance2 = distance + u_dis + POIs[o_id].dist-POIs[o_id].dis;
                        double distance_min = min(distance,distance2);
                        //double score = getGSKScoreu2O_Upper(u_id, a, distance, node_uKey_Map[top.id]);
                        double score = getGSKScoreu2o(u_id, o_id,a, alpha, distance_min);
                        //cout<<"o"<<o_id<<",(Ni)vertex"<<cands[i]<<" score = "<<score<<endl;
                        //Status_query status = {cands[i], true, top.lca_pos, result[i],score};
                        //Status_query status = {o_id, true, top.lca_pos, result[i],score};
                        //pq.push_back(status);
                        //push_heap(pq.begin(), pq.end(), Status_query_max_comp());
                        QEntry entry(o_id, true, top.lca_pos, result[i],score);
                        Queue.push(entry);
                        cout<<"加入o"<<o_id<<"distance="<<distance_min<<",score="<<score<<endl;

                    }


                }
                    // else do
                else {          //是叶节点，但与loc不在同一子图
                    for(int id :occurence_list[top.id]){
                        int o_id = id;

                        int Ni = POIs[o_id].Ni;
                        int Nj = POIs[o_id].Nj;
                        double score = 0;
                        for (int i = 0; i < GTree[top.id].leafinvlist.size(); i++) {    //获取当前最优叶节点的顶点倒排表
                            posa = GTree[top.id].leafinvlist[i];
                            vertex = GTree[top.id].leafnodes[posa];
                            if(vertex == Ni || vertex == Nj){
                                allmin = -1;

                                for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                    //计算loc 经过 叶节点中第 k个 border 到达 top.id叶节点内第i个顶点的距离
                                    dis = itm[top.id][k] + GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa];
                                    if (allmin == -1) {
                                        allmin = dis;
                                    } else {
                                        if (dis < allmin) {
                                            allmin = dis;   // 提取loc 经过所有叶节点top.id中的boder, 而到达第i个顶点的最短距离
                                        }
                                    }

                                }
                                //把所有的allmin
                                double distance = allmin /1000;
                                distance = distance + u_dis + POIs[o_id].dis;

                                score = max(score,getGSKScoreu2o(u_id, o_id,a, alpha, distance));
                                //double score = getGSKScoreu2O_Upper(u_id, a, allmin, node_uKey_Map[top.id]);
                                //cout<<"vertex_score="<<score<<endl;

                            }

                        }
                        if(score>0){
                            //Status_query status = {o_id, true, top.lca_pos, allmin, score};
                            //pq.push_back(status);
                            //push_heap(pq.begin(), pq.end(), Status_query_max_comp());
                            QEntry entry(o_id, true, top.lca_pos, allmin, score);
                            Queue.push(entry);
                            cout<<"加入o"<<o_id<<"allmin="<<allmin<<",score="<<score<<endl;
                            /*if(1682==o_id){
                                cout<<"发现1682，score="<<entry.score<<endl;
                                cout<<"leaf"<<top.id<<",score="<<top.score<<endl;
                                getchar();
                            }*/
                        }


                    }
                }
            } else if(occurence_list[top.id].size()>0){    //为非叶节点,且occurence_list不为空
                //for (int i = 0; i < GTree[top.id].nonleafinvlist.size(); i++) {
                //child = GTree[top.id].nonleafinvlist[i];
                for(int child: occurence_list[top.id]){

                    son = Nodes[locid].gtreepath[top.lca_pos + 1];
                    // on gtreepath
                    if (child == son) {
                        double score = getGSKScoreu2O_Upper(u_id, a, alpha, min_Dist, node_uKey_Map[child]);
                        //cout<<"child_score="<<score<<endl;
                        double min_Dist =0;

                        if(score > score_bound){
                            QEntry entry(child, false, top.lca_pos + 1, min_Dist,score);
                            Queue.push(entry);
                            //cout<<"加入o"<<o_id<<"distance="<<min_Dist<<",score="<<score<<endl;
                        }

                    }
                        // brothers
                    else if (GTree[child].father == GTree[son].father) {
                        itm[child].clear();
                        allmin = -1;

                        for (int j = 0; j < GTree[child].borders.size(); j++) {   //把到当前最优节点的各个child的距离算出，把child加入
                            _min = -1;
                            posa = GTree[child].up_pos[j];
                            for (int k = 0; k < GTree[son].borders.size(); k++) {   //第k个border
                                posb = GTree[son].up_pos[k];
                                dis = itm[son][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                                double dist = dis;
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
                        double distance = allmin/WEIGHT_INFLATE_FACTOR;
                        double score = getGSKScoreu2O_Upper(u_id, a, alpha,distance, node_uKey_Map[child]);
                        //cout<<"child_score(brother)="<<score<<endl;
                        //Status_query status = {child, false, top.lca_pos, allmin, score};
                        //pq.push_back(status);
                        //push_heap(pq.begin(), pq.end(), Status_query_max_comp());
                        if(score > score_bound){
                            QEntry entry(child, false, top.lca_pos, allmin, score);
                            Queue.push(entry);

                        }
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
                                dis = itm[top.id][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
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
                        double distance = allmin/WEIGHT_INFLATE_FACTOR;
                        double score = getGSKScoreu2O_Upper(u_id, a, alpha, distance, node_uKey_Map[child]);

                        if(score > score_bound){
                            QEntry entry(child, false, top.lca_pos, allmin, score);
                            Queue.push(entry);

                        }

                    }
                }
            }

        }
    }
    //block_num +=(count/1024);

    Users[u_id].topkScore_current = rstset[K-1].score;  //记得cache user 的 topk结果
    Users[u_id].topkScore_Final = rstset[K-1].score;
    return rstset;
}


/*
 * 基于GIM-tree 1. 求top-k结果， 2. 验证是否在top-k结果中
 */
//for update the range of subgrap in TkGSKQ_bottom2up_memory
void UpdateR_Memory(int u_id, int locid, int a, double alpha, double score_bound, int& node_highest,int& node_pos, vector<int> keywords, unordered_map<int, vector<int> >& itm, priority_queue<QEntry>& Queue, double& score_max, bool& eary_terminate){
    //表示当前范围内已全部遍历， 需扩大遍历子图范围
    int posa, posb, _min, dis;
    int root_id =0;
    int node_highest_pre = node_highest;
    int second_highest =0;
    //int node_highest = GTree[node_highest].father;
    int current = node_highest; int current_pre = current;
    set<int> entry_set;
    //cout<<"从"<<current_pre<<"向上寻找"<<endl;
    while(current != root_id){
        // 获得该中间节点下包含用户关键词的所有child
        current_pre = current;
        current = GTree[current].father;
        node_pos--;   //最高层节点层数往上一层
        for (int term: keywords) {
            if (GTree[current].inverted_list_o[term].size() > 0) {
                for (int child: GTree[current].inverted_list_o[term]) {
                    entry_set.insert(child);
                    //block_num++;
                }
            }
        }
        //父节点就有关键词
        if(entry_set.size()>0){
            node_highest = current;
            second_highest = current_pre;
            //cout<<"到"<<current<<",发现其有关键词！"<<endl;
            break;
        }
        else{

            //继续往上寻找关键词LCA
            //cout<<"到"<<current<<endl;
        }

    }
    //将扩大的查询范围子图加入
    for(int child: entry_set){
        if(child == current_pre) {
            //cout<<"跳过"<<child<<endl;
            //getchar();
            continue; //跳过已经遍历过的节点
        }
        //int son = Nodes[locid].gtreepath.back()(错误的！)
        //int node_pos = 0;
        int son = Nodes[locid].gtreepath[node_pos + 1];             //重要！容易出错
        // on gtreepath
        if (child == son) {
            double min_Dist =0;
            double score = getGSKScoreu2O_Upper(u_id, a, alpha, min_Dist, GTree[child].objectUKeySet);

            if(score > score_bound){
                QEntry entry(child, false, node_pos+1, min_Dist,score);
                Queue.push(entry);
                //cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
            }

        }
            // brothers
        else if (GTree[child].father == GTree[son].father) {
            itm[child].clear();
            int allmin = -1;

            for (int j = 0; j < GTree[child].borders.size(); j++) {   //把到当前最优节点的各个child的距离算出，把child加入
                _min = -1;
                posa = GTree[child].up_pos[j];
                for (int k = 0; k < GTree[son].borders.size(); k++) {   //第k个border
                    posb = GTree[son].up_pos[k];
                    dis = itm[son][k] + GTree[current].mind[posb * GTree[current].union_borders.size() + posa];
                    //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                    double dist = dis;
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
            double distance = allmin/WEIGHT_INFLATE_FACTOR;
            double score = getGSKScoreu2O_Upper(u_id, a, alpha, distance, GTree[child].objectUKeySet);

            if(score > score_bound){
                QEntry entry(child, false, node_pos, allmin, score);
                Queue.push(entry);
                //cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
            }
        }
            // downstream (也按brother的方式来)
        else {
            itm[child].clear();
            int allmin = -1;

            for (int j = 0; j < GTree[child].borders.size(); j++) {
                _min = -1;
                posa = GTree[child].up_pos[j];
                for (int k = 0; k < GTree[second_highest].borders.size(); k++) {  //重要！ 容易出错！
                    posb = GTree[second_highest].current_pos[k];
                    //int right_end = second_highest;
                    dis = itm[second_highest][k] + GTree[current].mind[posb * GTree[current].union_borders.size() + posa];
                    //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
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
            double distance = allmin/WEIGHT_INFLATE_FACTOR;
            double score = getGSKScoreu2O_Upper(u_id, a, alpha, distance, GTree[child].objectUKeySet);

            if(score > score_bound){
                QEntry entry(child, false, node_pos, allmin, score);
                Queue.push(entry);
                //cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
            }

        }

    }
    //cout<<"update后的 topscore="<<Queue.top().score<<"score_k="<<score_max<<endl;
    if(Queue.size()>0){
        if(score_max > Queue.top().score){
            //priority_queue<QEntry> emptyQ;
            //swap(emptyQ,Queue);     // 把当前队列清空!
            eary_terminate = true;
        }
    }


}


//1
priority_queue<Result> TkGSKQ_bottom2up_memory(int u_id, int K, float a, float alpha) {   //注意：这里的score_bound因为：lcl.rt,


    priority_queue<Result> resultFinal;
    priority_queue<QEntry> Queue;


    //获取user信息
    int locid = Users[u_id].Ni;
    double u_dis = Users[u_id].dis;
    vector<int> keywords = Users[u_id].keywords;
    int u_leaf = Nodes[Users[u_id].Ni].gtreepath.back();

    // 计算 loc 到 所有节点的border的距离
    unordered_map<int, vector<int> > itm; // intermediate answer, tree node -> array
    itm.clear();
    int cid, posa, posb, _min, dis;
    for (int i = Nodes[locid].gtreepath.size() - 1; i > 0; i--) {   //自底向上
        int tn = Nodes[locid].gtreepath[i];
        itm[tn].clear();
        //block_num++;

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


    //先把 u_leaf压入最大优先队列中
    int root_id = 0;
    int node_current = u_leaf;
    int node_highest = node_current;
    double score_max = -1;  double score_k = -1;
    double min_Dist = 0; int node_pos = Nodes[locid].gtreepath.size()-1;
    double score_leaf = getGSKScoreu2O_Upper(u_id, a, alpha, min_Dist, GTree[u_leaf].userUKeySet);
    QEntry entry(u_leaf,false,node_pos,0, score_leaf);
    Queue.push(entry);
    // 自底向上进行节点遍历迭代
    //cout<<"开始自底向上进行节点遍历迭代"<<endl;
    bool eary_terminate = false;
    while(Queue.size() > 0 || node_highest != root_id){//node_highest != root_id){  //Queue.size() > 0 ||

        //当前范围内已全部遍历， 需更新扩大路网图的遍历范围
        if(Queue.empty()){
            //cout<<"当前范围内已全部遍历， 需更新扩大路网图的遍历范围"<<endl;
            UpdateR_Memory(u_id, locid, a,  alpha, score_max, node_highest, node_pos, keywords, itm, Queue, score_max, eary_terminate);
        }
        //当前范围内遍历
        //cout<<"node_highest="<<node_highest<<endl;
        QEntry top = Queue.top();
        double topk_score = top.score;


        if(score_max > topk_score) {   // 当前队列中已无更优元素
            //priority_queue<QEntry> emptyQ;
            //swap(emptyQ,Queue);     // 把当前队列清空!
            if (node_highest!= root_id) { //如果还未到最高层，那就再拓展一波
                UpdateR_Memory(u_id, locid, a,  alpha, score_max, node_highest, node_pos, keywords, itm, Queue, score_max,
                               eary_terminate);
            }
            else {
                if(Queue.size()>0){
                    top = Queue.top();
                    if(score_k > top.score){   //提前终止条件！
                        if(top.isvertex){
                            //cout<<"提前终止条件，出队列，o="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                        }
                        else if(GTree[top.id].isleaf){
                            //cout<<"提前终止条件，出队列，leaf="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                        } else{
                            //cout<<"提前终止条件，出队列，n="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                        }
                        break;
                        //getchar();
                    }
                }

            }
            //UpdateR(u_id, locid, a, score_bound, node_highest, keywords, itm, Queue, score_max);
            //top = Queue.top();


        }
        //当队列内还有更优的元素时
        if(Queue.size()>0){
            top = Queue.top();
            Queue.pop();
            if (top.isvertex) {  //当前最优为某个路网上的兴趣点（顶点）
                //cout<<"出队列，o="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
                if(resultFinal.size()< K){
                    int poi = top.id;
                    Result tmpRlt(poi, POIs[poi].Ni, POIs[poi].Nj, POIs[poi].dis, POIs[poi].dist, top.score, POIs[poi].keywords);
                    resultFinal.push(tmpRlt);
                    if(resultFinal.size()== K) {
                        score_k = resultFinal.top().score;
                        score_max = score_k;
                    }
                }
                else{  // update result,  score_k
                    int poi = top.id;
                    Result tmpRlt(poi, POIs[poi].Ni, POIs[poi].Nj, POIs[poi].dis, POIs[poi].dist, top.score, POIs[poi].keywords);
                    resultFinal.push(tmpRlt);
                    resultFinal.pop();
                    score_k = resultFinal.top().score;
                    score_max = max(score_max, score_k);
                    //cout<<"更新score_k="<<score_k<<endl;
                }

            }
            else if (GTree[top.id].isleaf) {           //当前最优为某个叶节点
                //cout<<"出队列，leaf="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
                if (top.id == Nodes[locid].gtreepath[top.lca_pos]) {  //该叶节点为usr所在的叶节点
                    vector<int> cands;
                    vector<int> result;
                    set<int> object_set;
                    // 获得该叶节点下包含用户关键词的所有poi
                    for (int term: Users[u_id].keywords) {
                        if (GTree[top.id].inverted_list_o[term].size() > 0) {
                            for (int poi: GTree[top.id].inverted_list_o[term]) {
                                object_set.insert(poi);
                            }
                        }
                    }
                    //获得user对这些poi的评分
                    vector<int> objects;
                    for (int poi : object_set) {
                        int Ni = POIs[poi].Ni;
                        int Nj = POIs[poi].Nj;
                        cands.push_back(Ni);
                        cands.push_back(Nj);
                        //objects.push_back(object_id);
                        objects.push_back(poi);
                        objects.push_back(poi);
                    }
                    result = dijkstra_candidate(locid, cands, Nodes);
                    //将叶节点内的各poi(顶点)插入
                    for (int i = 0; i < cands.size(); i += 2) {
                        int o_id = objects[i];

                        double distance = result[i] / 1000;
                        distance = distance + u_dis + POIs[o_id].dis;
                        double distance2 = result[i + 1] / 1000;
                        distance2 = distance2 + u_dis + POIs[o_id].dist - POIs[o_id].dis;
                        double distance_min = min(distance, distance2);
                        //double score = getGSKScoreu2O_Upper(u_id, a, distance, node_uKey_Map[top.id]);
                        double score = getGSKScoreu2o(u_id, o_id, a, alpha, distance_min);
                        if (score > score_max) {
                            QEntry entry(o_id, true, top.lca_pos, result[i], score);
                            Queue.push(entry);
                            //cout<<"加入o"<<o_id<<"distance="<<distance_min<<",score="<<score<<endl;

                        }

                    }
                } else {  // 是叶节点，但与loc不在同一子图
                    set<int> object_set;
                    // 获得该叶节点下包含用户关键词的所有poi
                    for (int term: Users[u_id].keywords) {
                        if (GTree[top.id].inverted_list_o[term].size() > 0) {
                            for (int poi: GTree[top.id].inverted_list_o[term]) {
                                object_set.insert(poi);
                            }
                        }
                    }
                    for (int id: object_set) {
                        int o_id = id;

                        int Ni = POIs[o_id].Ni;
                        int Nj = POIs[o_id].Nj;
                        double score = 0; int allmin = -1; double distance;
                        for (int i = 0; i < GTree[top.id].leafinvlist.size(); i++) {    //获取当前最优叶节点的顶点倒排表
                            int posa = GTree[top.id].leafinvlist[i];
                            int vertex = GTree[top.id].leafnodes[posa];
                            if (vertex == Ni || vertex == Nj) {
                                allmin = -1;
                                for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                    //计算loc 经过 叶节点中第 k个 border 到达 top.id叶节点内第i个顶点的距离
                                    dis = itm[top.id][k] + GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa];
                                    //cout<<"dis="<<dis<<",itm[top.id][k]="<<itm[top.id][k]<<",边界间距离="<<GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa]<<endl;
                                    if (allmin == -1) {
                                        allmin = dis;
                                    } else {
                                        if (dis < allmin) {
                                            allmin = dis;   // 提取loc 经过所有叶节点top.id中的boder, 而到达第i个顶点的最短距离
                                        }
                                    }

                                }
                                //把所有的allmin
                                distance = allmin / 1000;
                                distance = distance + u_dis + POIs[o_id].dis;
                                score = max(score, getGSKScoreu2o(u_id, o_id, a, alpha, distance));

                            }

                        }

                        if (score > score_max) {
                            QEntry entry(o_id, true, top.lca_pos, distance, score);
                            Queue.push(entry);
                            //cout<<"加入o"<<o_id<<"distance="<<distance<<",score="<<score<<endl;

                        }

                    }

                }

            }
            else {    //当前最优为某个中间节点
                //cout<<"出队列，n="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
                set<int> entry_set;
                // 获得该中间节点下包含用户关键词的所有child
                for (int term: Users[u_id].keywords) {
                    if (GTree[top.id].inverted_list_o[term].size() > 0) {
                        for (int child: GTree[top.id].inverted_list_o[term]) {
                            entry_set.insert(child);
                        }
                    }
                }
                //将各子条目加入优先队列
                for(int child: entry_set){

                    int son = Nodes[locid].gtreepath[top.lca_pos + 1];
                    // on gtreepath
                    if (child == son) {
                        double score = getGSKScoreu2O_Upper(u_id, a, alpha, min_Dist, GTree[child].objectUKeySet);
                        //cout<<"child_score="<<score<<endl;
                        double min_Dist =0;

                        if(score > score_max){
                            QEntry entry(child, false, top.lca_pos + 1, min_Dist,score);
                            Queue.push(entry);

                        }

                    }
                        // brothers
                    else if (GTree[child].father == GTree[son].father) {
                        itm[child].clear();
                        int allmin = -1;

                        for (int j = 0; j < GTree[child].borders.size(); j++) {   //把到当前最优节点的各个child的距离算出，把child加入
                            _min = -1;
                            posa = GTree[child].up_pos[j];
                            for (int k = 0; k < GTree[son].borders.size(); k++) {   //第k个border
                                posb = GTree[son].up_pos[k];
                                dis = itm[son][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                                double dist = dis;
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
                        double distance = allmin/WEIGHT_INFLATE_FACTOR;
                        double score = getGSKScoreu2O_Upper(u_id, a, alpha, distance, GTree[child].objectUKeySet);

                        if(score > score_max){
                            QEntry entry(child, false, top.lca_pos, allmin, score);
                            Queue.push(entry);

                        }
                    }
                        // downstream
                    else {
                        itm[child].clear();
                        int allmin = -1;

                        for (int j = 0; j < GTree[child].borders.size(); j++) {
                            _min = -1;
                            posa = GTree[child].up_pos[j];
                            for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                posb = GTree[top.id].current_pos[k];
                                dis = itm[top.id][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
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
                        double distance = allmin/WEIGHT_INFLATE_FACTOR;
                        double score = getGSKScoreu2O_Upper(u_id, a, alpha, distance, GTree[child].objectUKeySet);

                        if(score > score_max){
                            QEntry entry(child, false, top.lca_pos, allmin, score);
                            Queue.push(entry);
                        }

                    }

                }


            }


        }




    }//while
    if(resultFinal.size()>0){    //倘若resultFinal为空，说明没有任何元素比score_bound更优
        double ss = resultFinal.top().score;
        if(ss> Users[u_id].topkScore_current){
            Users[u_id].topkScore_current = ss;  //记得cache user 的 topk结果
            Users[u_id].topkScore_Final = ss;

        }
    }
    return resultFinal;
}


//2
//验证在gsk评分下某个poi是否在u的topk结果中  当gsk未知时 设为无穷大， 如1000
/*提前终止条件：1. resultFinal.rk > gsk_score(u,p);
*剪枝条件： score_bound = lcl.r_t > gsk_score_upper(u,O), 当其未知时，设为0*/
priority_queue<Result> TkGSKQ_bottom2up_verify_memory(int u_id, int K, float a, float alpha, double gsk, double score_bound) {


    priority_queue<Result> resultFinal;
    priority_queue<QEntry> Queue;


    //获取user信息
    int locid = Users[u_id].Ni;
    double u_dis = Users[u_id].dis;
    vector<int> keywords = Users[u_id].keywords;
    int u_leaf = Nodes[Users[u_id].Ni].gtreepath.back();

    // 计算 loc 到 所有节点的border的距离
    unordered_map<int, vector<int> > itm; // intermediate answer, tree node -> array
    itm.clear();
    int cid, posa, posb, _min, dis;
    for (int i = Nodes[locid].gtreepath.size() - 1; i > 0; i--) {   //自底向上
        int tn = Nodes[locid].gtreepath[i];
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


    //先把 u_leaf压入最大优先队列中
    int root_id = 0;
    int node_current = u_leaf;
    int node_highest = node_current;
    double score_max = score_bound;  double score_k = score_bound;
    double min_Dist = 0; int node_pos = Nodes[locid].gtreepath.size()-1;
    double score_leaf = getGSKScoreu2O_Upper(u_id, a, alpha, min_Dist, GTree[u_leaf].userUKeySet);
    QEntry entry(u_leaf,false,node_pos,0, score_leaf);
    Queue.push(entry);
    // 自底向上进行节点遍历迭代
    //cout<<"开始自底向上进行节点遍历迭代"<<endl;
    bool eary_terminate = false;
    while(Queue.size() > 0 || node_highest != root_id){//node_highest != root_id){  //Queue.size() > 0 ||
        if(score_max>gsk)
            break; //提前终止条件：已发现k个对象比目标poi有更高的gsk评分

        //当前范围内已全部遍历， 需更新扩大路网图的遍历范围
        if(Queue.empty()){
            //cout<<"当前范围内已全部遍历， 需更新扩大路网图的遍历范围"<<endl;
            UpdateR_Memory(u_id, locid, a,  alpha, score_max, node_highest, node_pos, keywords, itm, Queue, score_max, eary_terminate);
        }
        //当前范围内遍历
        //cout<<"node_highest="<<node_highest<<endl;
        QEntry top = Queue.top();
        double topk_score = top.score;


        if(score_max > topk_score) {   // 当前队列中已无更优元素
            //priority_queue<QEntry> emptyQ;
            //swap(emptyQ,Queue);     // 把当前队列清空!
            if (node_highest!= root_id) { //如果还未到最高层，那就再拓展一波
                UpdateR_Memory(u_id, locid, a,  alpha, score_max, node_highest, node_pos, keywords, itm, Queue, score_max,
                               eary_terminate);
            }
            else {
                if(Queue.size()>0){
                    top = Queue.top();
                    if(score_k > top.score){   //提前终止条件！
                        if(top.isvertex){
                            //cout<<"提前终止条件，出队列，o="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                        }
                        else if(GTree[top.id].isleaf){
                            //cout<<"提前终止条件，出队列，leaf="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                        } else{
                            //cout<<"提前终止条件，出队列，n="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                        }

                        break;
                    }
                }

            }



        }
        //当队列内还有更优的元素时
        if(Queue.size()>0){
            top = Queue.top();
            Queue.pop();
            if (top.isvertex) {  //当前最优为某个路网上的兴趣点（顶点）
                //cout<<"出队列，o="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
                if(resultFinal.size()< K){
                    int poi = top.id;
                    Result tmpRlt(poi, POIs[poi].Ni, POIs[poi].Nj, POIs[poi].dis, POIs[poi].dist, top.score, POIs[poi].keywords);
                    resultFinal.push(tmpRlt);
                    if(resultFinal.size()== K) {
                        score_k = resultFinal.top().score;
                        score_max = score_k;
                    }
                }
                else{  // update result,  score_k
                    int poi = top.id;
                    Result tmpRlt(poi, POIs[poi].Ni, POIs[poi].Nj, POIs[poi].dis, POIs[poi].dist, top.score, POIs[poi].keywords);
                    resultFinal.push(tmpRlt);
                    resultFinal.pop();
                    score_k = resultFinal.top().score;
                    score_max = max(score_max, score_k);
                    //cout<<"更新score_k="<<score_k<<endl;
                }

            }
            else if (GTree[top.id].isleaf) {           //当前最优为某个叶节点
                //cout<<"出队列，leaf="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
                if (top.id == Nodes[locid].gtreepath[top.lca_pos]) {  //该叶节点为usr所在的叶节点
                    vector<int> cands;
                    vector<int> result;
                    set<int> object_set;
                    // 获得该叶节点下包含用户关键词的所有poi
                    for (int term: Users[u_id].keywords) {
                        if (GTree[top.id].inverted_list_o[term].size() > 0) {
                            for (int poi: GTree[top.id].inverted_list_o[term]) { //access disk (leaf_id, term)
                                object_set.insert(poi);
                            }
                        }
                    }
                    //获得user对这些poi的评分
                    vector<int> objects;
                    for (int poi : object_set) {
                        int Ni = POIs[poi].Ni;
                        int Nj = POIs[poi].Nj;
                        cands.push_back(Ni);
                        cands.push_back(Nj);
                        //objects.push_back(object_id);
                        objects.push_back(poi);
                        objects.push_back(poi);
                    }
                    result = dijkstra_candidate(locid, cands, Nodes);
                    //将叶节点内的各poi(顶点)插入
                    for (int i = 0; i < cands.size(); i += 2) {
                        int o_id = objects[i];

                        double distance = result[i] / 1000;
                        distance = distance + u_dis + POIs[o_id].dis;
                        double distance2 = result[i + 1] / 1000;
                        distance2 = distance2 + u_dis + POIs[o_id].dist - POIs[o_id].dis;
                        double distance_min = min(distance, distance2);
                        //double score = getGSKScoreu2O_Upper(u_id, a, distance, node_uKey_Map[top.id]);
                        double score = getGSKScoreu2o(u_id, o_id, a, alpha,distance_min);
                        if (score > score_max) {
                            QEntry entry(o_id, true, top.lca_pos, result[i], score);
                            Queue.push(entry);
                            //cout<<"加入o"<<o_id<<"distance="<<distance_min<<",score="<<score<<endl;

                        }

                    }
                } else {  // 是叶节点，但与loc不在同一子图
                    set<int> object_set;
                    // 获得该叶节点下包含用户关键词的所有poi
                    for (int term: Users[u_id].keywords) {
                        if (GTree[top.id].inverted_list_o[term].size() > 0) {
                            for (int poi: GTree[top.id].inverted_list_o[term]) {
                                object_set.insert(poi);
                            }
                        }
                    }
                    for (int id: object_set) {
                        int o_id = id;

                        //else cout<<"发现其他"<<endl;

                        int Ni = POIs[o_id].Ni;
                        int Nj = POIs[o_id].Nj;
                        double score = 0; int allmin = -1; double distance;
                        for (int i = 0; i < GTree[top.id].leafinvlist.size(); i++) {    //获取当前最优叶节点的顶点倒排表
                            int posa = GTree[top.id].leafinvlist[i];
                            int vertex = GTree[top.id].leafnodes[posa];
                            if (vertex == Ni || vertex == Nj) {
                                allmin = -1;
                                for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                    //计算loc 经过 叶节点中第 k个 border 到达 top.id叶节点内第i个顶点的距离
                                    dis = itm[top.id][k] + GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa];
                                    //cout<<"dis="<<dis<<",itm[top.id][k]="<<itm[top.id][k]<<",边界间距离="<<GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa]<<endl;
                                    if (allmin == -1) {
                                        allmin = dis;
                                    } else {
                                        if (dis < allmin) {
                                            allmin = dis;   // 提取loc 经过所有叶节点top.id中的boder, 而到达第i个顶点的最短距离
                                        }
                                    }

                                }
                                //把所有的allmin
                                distance = allmin / 1000;
                                distance = distance + u_dis + POIs[o_id].dis;
                                score = max(score, getGSKScoreu2o(u_id, o_id, a, alpha,distance));

                            }

                        }

                        if (score > score_max) {
                            QEntry entry(o_id, true, top.lca_pos, distance, score);
                            Queue.push(entry);
                            //cout<<"加入o"<<o_id<<"distance="<<distance<<",score="<<score<<endl;

                        }


                    }

                }

            }
            else {    //当前最优为某个中间节点
                //cout<<"出队列，n="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
                set<int> entry_set;
                // 获得该中间节点下包含用户关键词的所有child
                for (int term: Users[u_id].keywords) {
                    if (GTree[top.id].inverted_list_o[term].size() > 0) {
                        for (int child: GTree[top.id].inverted_list_o[term]) {
                            entry_set.insert(child);
                        }
                    }
                }
                //将各子条目加入优先队列
                for(int child: entry_set){

                    int son = Nodes[locid].gtreepath[top.lca_pos + 1];
                    // on gtreepath
                    if (child == son) {
                        double score = getGSKScoreu2O_Upper(u_id, a, alpha, min_Dist, GTree[child].objectUKeySet);
                        //cout<<"child_score="<<score<<endl;
                        double min_Dist =0;

                        if(score > score_max){
                            QEntry entry(child, false, top.lca_pos + 1, min_Dist,score);
                            Queue.push(entry);

                        }

                    }
                        // brothers
                    else if (GTree[child].father == GTree[son].father) {
                        itm[child].clear();
                        int allmin = -1;

                        for (int j = 0; j < GTree[child].borders.size(); j++) {   //把到当前最优节点的各个child的距离算出，把child加入
                            _min = -1;
                            posa = GTree[child].up_pos[j];
                            for (int k = 0; k < GTree[son].borders.size(); k++) {   //第k个border
                                posb = GTree[son].up_pos[k];
                                dis = itm[son][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                                double dist = dis;
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
                        double distance = allmin/WEIGHT_INFLATE_FACTOR;
                        double score = getGSKScoreu2O_Upper(u_id, a,  alpha,distance, GTree[child].objectUKeySet);

                        if(score > score_max){
                            QEntry entry(child, false, top.lca_pos, allmin, score);
                            Queue.push(entry);

                        }
                    }
                        // downstream
                    else {
                        itm[child].clear();
                        int allmin = -1;

                        for (int j = 0; j < GTree[child].borders.size(); j++) {
                            _min = -1;
                            posa = GTree[child].up_pos[j];
                            for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                posb = GTree[top.id].current_pos[k];
                                dis = itm[top.id][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
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
                        double distance = allmin/WEIGHT_INFLATE_FACTOR;
                        double score = getGSKScoreu2O_Upper(u_id, a, alpha,distance, GTree[child].objectUKeySet);

                        if(score > score_max){
                            QEntry entry(child, false, top.lca_pos, allmin, score);
                            Queue.push(entry);
                        }

                    }

                }


            }


        }




    }//while
    if(resultFinal.size()>0){    //倘若resultFinal为空，说明没有任何元素比score_bound更优
        double ss = resultFinal.top().score;
        if(ss> Users[u_id].topkScore_current){
            Users[u_id].topkScore_current = ss;  //记得cache user 的 topk结果
            //Users[u_id].topkScore_Final = ss;  //easy error

        }
    }
    return resultFinal;
}

int upper_jumper = 0;



//for update the range of subgrap in TkGSKQ_bottom2up_memory
void UpdateR_disk(int u_id, int locid, int a, double alpha, double score_bound, int& node_highest,int& node_pos, vector<int> keywords, unordered_map<int, vector<int> >& itm, priority_queue<QEntry>& Queue, double& score_max, bool& eary_terminate){
    //表示当前范围内已全部遍历， 需扩大遍历子图范围
    int posa, posb, _min, dis;
    int root_id =0;
    int node_highest_pre = node_highest;
    int second_highest =0;
    //int node_highest = GTree[node_highest].father;
    int current = node_highest; int current_pre = current;
    set<int> entry_set;
    //cout<<"从"<<current_pre<<"向上寻找"<<endl;
    while(current != root_id){
        // 获得该中间节点下包含用户关键词的所有child
        current_pre = current;
        current = GTree[current].father;
        upper_jumper++;
        node_pos--;   //最高层节点层数往上一层
        for (int term: keywords) {
            vector<int> o_entrylist = getObjectTermRelatedEntry(term,current);
            if (o_entrylist.size() > 0) {
                for (int child: o_entrylist) {
                    entry_set.insert(child);
                    //block_num++;
                }
            }
        }
        //父节点就有关键词
        if(entry_set.size()>0){
            node_highest = current;
            second_highest = current_pre;
            //cout<<"到"<<current<<",发现其有关键词！"<<endl;
            break;
        }
        else{

            //继续往上寻找关键词LCA
            //cout<<"到"<<current<<endl;
        }

    }
    //将扩大的查询范围子图加入
    for(int child: entry_set){
        if(child == current_pre) {
            //cout<<"跳过"<<child<<endl;
            //getchar();
            continue; //跳过已经遍历过的节点
        }
        //int son = Nodes[locid].gtreepath.back()(错误的！)
        //int node_pos = 0;
        int son = Nodes[locid].gtreepath[node_pos + 1];             //重要！容易出错
        TreeNode child_Node = getGIMTreeNodeData(child, OnlyO);  //gangc
        // on gtreepath
        if (child == son) {
            double min_Dist =0;
            double score = getGSKScoreu2O_Upper(u_id, a, alpha, min_Dist, child_Node.objectUKeySet);

            if(score > score_bound){
                QEntry entry(child, false, node_pos+1, min_Dist,score);
                Queue.push(entry);
                //cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
            }

        }
            // brothers
        else if (GTree[child].father == GTree[son].father) {
            itm[child].clear();
            int allmin = -1;

            for (int j = 0; j < GTree[child].borders.size(); j++) {   //把到当前最优节点的各个child的距离算出，把child加入
                _min = -1;
                posa = GTree[child].up_pos[j];
                for (int k = 0; k < GTree[son].borders.size(); k++) {   //第k个border
                    posb = GTree[son].up_pos[k];
                    dis = itm[son][k] + GTree[current].mind[posb * GTree[current].union_borders.size() + posa];
                    //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                    double dist = dis;
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
            double distance = allmin/WEIGHT_INFLATE_FACTOR;
            double score = getGSKScoreu2O_Upper(u_id, a, alpha, distance, child_Node.objectUKeySet);

            if(score > score_bound){
                QEntry entry(child, false, node_pos, allmin, score);
                Queue.push(entry);
                //cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
            }
        }
            // downstream (也按brother的方式来)
        else {
            itm[child].clear();
            int allmin = -1;

            for (int j = 0; j < GTree[child].borders.size(); j++) {
                _min = -1;
                posa = GTree[child].up_pos[j];
                for (int k = 0; k < GTree[second_highest].borders.size(); k++) {  //重要！ 容易出错！
                    posb = GTree[second_highest].current_pos[k];
                    //int right_end = second_highest;
                    dis = itm[second_highest][k] + GTree[current].mind[posb * GTree[current].union_borders.size() + posa];
                    //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
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
            double distance = allmin/WEIGHT_INFLATE_FACTOR;
            double score = getGSKScoreu2O_Upper(u_id, a, alpha, distance, child_Node.objectUKeySet);

            if(score > score_bound){
                QEntry entry(child, false, node_pos, allmin, score);
                Queue.push(entry);
                //cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
            }

        }

    }
    //cout<<"update后的 topscore="<<Queue.top().score<<"score_k="<<score_max<<endl;
    if(Queue.size()>0){
        if(score_max > Queue.top().score){
            //priority_queue<QEntry> emptyQ;
            //swap(emptyQ,Queue);     // 把当前队列清空!
            eary_terminate = true;
        }
    }


}


//for extending the valuation range of user in RkGSKQ_bottom2up_by_NVD
void UpdateR_disk_for_user(int u_id, int locid, int a, double alpha, double score_bound, int& node_highest,int& node_pos, vector<int> keywords, unordered_map<int, vector<int> >& itm, priority_queue<QEntry>& Queue, double& score_max, bool& eary_terminate){
    //表示当前范围内已全部遍历， 需扩大遍历子图范围
    int posa, posb, _min, dis;
    int root_id =0;
    int node_highest_pre = node_highest;
    int second_highest =0;
    //int node_highest = GTree[node_highest].father;
    int current = node_highest; int current_pre = current;
    set<int> entry_set;
    //cout<<"从"<<current_pre<<"向上寻找"<<endl;
    while(current != root_id){
        // 获得该中间节点下包含用户关键词的所有child
        current_pre = current;
        current = GTree[current].father;
        upper_jumper++;
        node_pos--;   //最高层节点层数往上一层
        for (int term: keywords) {
            vector<int> o_entrylist = getObjectTermRelatedEntry(term,current);
            if (o_entrylist.size() > 0) {
                for (int child: o_entrylist) {
                    entry_set.insert(child);
                    //block_num++;
                }
            }
        }
        //父节点就有关键词
        if(entry_set.size()>0){
            node_highest = current;
            second_highest = current_pre;
            //cout<<"到"<<current<<",发现其有关键词！"<<endl;
            break;
        }
        else{

            //继续往上寻找关键词LCA
            //cout<<"到"<<current<<endl;
        }

    }
    //将扩大的查询范围子图加入
    for(int child: entry_set){
        if(child == current_pre) {
            //cout<<"跳过"<<child<<endl;
            //getchar();
            continue; //跳过已经遍历过的节点
        }
        //int son = Nodes[locid].gtreepath.back()(错误的！)
        //int node_pos = 0;
        int son = Nodes[locid].gtreepath[node_pos + 1];             //重要！容易出错
        TreeNode child_Node = getGIMTreeNodeData(child, OnlyO);  //gangc
        // on gtreepath
        if (child == son) {
            double min_Dist =0;
            double score = getGSKScoreu2O_Upper(u_id, a, alpha, min_Dist, child_Node.objectUKeySet);

            if(score > score_bound){
                QEntry entry(child, false, node_pos+1, min_Dist,score);
                Queue.push(entry);
                //cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
            }

        }
            // brothers
        else if (GTree[child].father == GTree[son].father) {
            itm[child].clear();
            int allmin = -1;

            for (int j = 0; j < GTree[child].borders.size(); j++) {   //把到当前最优节点的各个child的距离算出，把child加入
                _min = -1;
                posa = GTree[child].up_pos[j];
                for (int k = 0; k < GTree[son].borders.size(); k++) {   //第k个border
                    posb = GTree[son].up_pos[k];
                    dis = itm[son][k] + GTree[current].mind[posb * GTree[current].union_borders.size() + posa];
                    //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                    double dist = dis;
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
            double distance = allmin/WEIGHT_INFLATE_FACTOR;
            double score = getGSKScoreu2O_Upper(u_id, a, alpha, distance, child_Node.objectUKeySet);

            if(score > score_bound){
                QEntry entry(child, false, node_pos, allmin, score);
                Queue.push(entry);
                //cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
            }
        }
            // downstream (也按brother的方式来)
        else {
            itm[child].clear();
            int allmin = -1;

            for (int j = 0; j < GTree[child].borders.size(); j++) {
                _min = -1;
                posa = GTree[child].up_pos[j];
                for (int k = 0; k < GTree[second_highest].borders.size(); k++) {  //重要！ 容易出错！
                    posb = GTree[second_highest].current_pos[k];
                    //int right_end = second_highest;
                    dis = itm[second_highest][k] + GTree[current].mind[posb * GTree[current].union_borders.size() + posa];
                    //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
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
            double distance = allmin/WEIGHT_INFLATE_FACTOR;
            double score = getGSKScoreu2O_Upper(u_id, a, alpha, distance, child_Node.objectUKeySet);

            if(score > score_bound){
                QEntry entry(child, false, node_pos, allmin, score);
                Queue.push(entry);
                //cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
            }

        }

    }
    //cout<<"update后的 topscore="<<Queue.top().score<<"score_k="<<score_max<<endl;
    if(Queue.size()>0){
        if(score_max > Queue.top().score){
            //priority_queue<QEntry> emptyQ;
            //swap(emptyQ,Queue);     // 把当前队列清空!
            eary_terminate = true;
        }
    }


}



priority_queue<Result> TkGSKQ_bottom2up_verify_disk(int u_id, int K, float a, float alpha, double gsk, double score_bound) {

    //double direction

    priority_queue<Result> resultFinal;
    priority_queue<QEntry> Queue;

    User u = getUserFromO2UOrgLeafData(u_id);
    User query_u = u;

    //获取user信息
    int locid = u.Ni; int locid2 = u.Nj;
    double u_dis = u.dis;
    vector<int> keywords = u.keywords;
    int u_leaf = Nodes[u.Ni].gtreepath.back();

    // 计算 loc 到 所有节点的border的距离
    unordered_map<int, vector<int> > itm; // intermediate answer, tree node -> array
    itm.clear();
    int cid, posa, posb;
    for (int i = Nodes[locid].gtreepath.size() - 1; i > 0; i--) {   //自底向上
        int tn = Nodes[locid].gtreepath[i];
        itm[tn].clear();

        if (GTree[tn].isleaf) {

            int pos_Ni = id2idx_vertexInLeaf[locid];
            int pos_Nj = id2idx_vertexInLeaf[locid2];

            for (int j = 0; j < GTree[tn].borders.size(); j++) {
                int dist_via_Ni = GTree[tn].mind[j * GTree[tn].leafnodes.size() + pos_Ni]+u.dis;
                int dist_via_Nj = GTree[tn].mind[j * GTree[tn].leafnodes.size() + pos_Nj]+(u.dist-u.dis);
                int dist_to_border = min(dist_via_Ni, dist_via_Nj);
                itm[tn].push_back(dist_to_border);
            }
        } else {
            cid = Nodes[locid].gtreepath[i + 1];
            for (int j = 0; j < GTree[tn].borders.size(); j++) {
                int _min = -1;
                posa = GTree[tn].current_pos[j];
                for (int k = 0; k < GTree[cid].borders.size(); k++) {
                    posb = GTree[cid].up_pos[k];
                    int dis = itm[cid][k] + GTree[tn].mind[posa * GTree[tn].union_borders.size() + posb];
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

#ifdef  TOPKTRACK
    cout<<"完成 loc 到 所有节点的border的距离计算"<<endl;
#endif


    //先把 u_leaf压入最大优先队列中
    int root_id = 0;
    int node_current = u_leaf;
    int node_highest = node_current;
    double score_max = score_bound;  double score_k = score_bound;
    double min_Dist = 0; int node_pos = Nodes[locid].gtreepath.size()-1;

    TreeNode u_leafNode = getGIMTreeNodeData(u_leaf,OnlyU);
    double score_leaf = getGSKScoreu2O_Upper(u_id, a, alpha, min_Dist, u_leafNode.userUKeySet);
    QEntry entry(u_leaf,false,node_pos,0, score_leaf);
    Queue.push(entry);
    // 自底向上进行节点遍历迭代
    //cout<<"开始自底向上进行节点遍历迭代"<<endl;
    bool eary_terminate = false;
    while(Queue.size() > 0 || node_highest != root_id){//node_highest != root_id){  //Queue.size() > 0 ||
        if(score_max>gsk)
            break; //提前终止条件：已发现k个对象比目标poi有更高的gsk评分

        //当前范围内已全部遍历， 需更新扩大路网图的遍历范围
        if(Queue.empty()){
#ifdef TOPKTRACK
            cout<<"当前范围内已全部遍历， 需更新扩大路网图的遍历范围"<<endl;
#endif
            UpdateR_disk(u_id, locid, a,  alpha, score_max, node_highest, node_pos, keywords, itm, Queue, score_max, eary_terminate);
        }
        //当前范围内遍历
        //cout<<"node_highest="<<node_highest<<endl;
        QEntry top = Queue.top();
        double topk_score = top.score;


        if(score_max > topk_score) {   // 当前队列中已无更优元素
            //priority_queue<QEntry> emptyQ;
            //swap(emptyQ,Queue);     // 把当前队列清空!
            if (node_highest!= root_id) { //如果还未到最高层，那就再拓展一波
                UpdateR_disk(u_id, locid, a,  alpha, score_max, node_highest, node_pos, keywords, itm, Queue, score_max,
                             eary_terminate);
            }
            else {
                if(Queue.size()>0){
                    top = Queue.top();
                    if(score_k > top.score){   //提前终止条件！
                        if(top.isvertex){
                            //cout<<"提前终止条件，出队列，o="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                        }
                        else if(GTree[top.id].isleaf){
                            //cout<<"提前终止条件，出队列，leaf="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                        } else{
                            //cout<<"提前终止条件，出队列，n="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                        }

                        break;
                    }
                }

            }
            //UpdateR(u_id, locid, a, score_bound, node_highest, keywords, itm, Queue, score_max);
            //top = Queue.top();


        }
        //当队列内还有更优的元素时
        if(Queue.size()>0){
            top = Queue.top();
            Queue.pop();
            if (top.isvertex) {  //当前最优为某个路网上的兴趣点（顶点）
#ifdef  TOPKTRACK
                cout<<"出队列，o="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
#endif
                int poi = top.id;


            }
            else if (GTree[top.id].isleaf) {           //当前最优为某个叶节点
                //cout<<"出队列，leaf="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
#ifdef  TOPKTRACK
                cout<<"叶节点出队列，leaf="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
#endif
                if (top.id == Nodes[locid].gtreepath[top.lca_pos]) {  //该叶节点为usr所在的叶节点
                    vector<int> cands;
                    vector<int> result;
                    vector<int> result2;
                    set<int> object_set;
                    // 获得该叶节点下包含用户关键词的所有poi
                    for (int term: query_u.keywords) {
                        vector<int> o_list = getObjectTermRelatedEntry(term,top.id);
                        if (o_list.size() > 0) {
                            for (int poi: o_list) { //access disk (leaf_id, term)
                                object_set.insert(poi);
                            }
                        }
                    }
                    //获得user对这些poi的评分
                    vector<POI> objects;
                    for (int poi : object_set) {
                        //cout<<"将检索p"<<poi<<endl;
                        POI _p = getPOIFromO2UOrgLeafData(poi);
                        //printPOIInfo(_p);
                        int Ni = _p.Ni;
                        int Nj = _p.Nj;
                        cands.push_back(Ni);
                        cands.push_back(Nj);
                        objects.push_back(_p);
                        objects.push_back(_p);
                    }
                    result = dijkstra_candidate(locid, cands, Nodes);
                    result2 = dijkstra_candidate(locid2, cands, Nodes);
                    //将叶节点内的各poi(顶点)插入
                    for (int i = 0; i < cands.size(); i += 2) {
                        POI  _p = objects[i];
                        int o_id = _p.id;


                        double distance_min;
                        if(_p.Ni==u.Ni&&_p.Nj==u.Nj){  //这里要注意兴趣点和user在同一条路段上的情况！
                            distance_min = fabs(u.dis - _p.dis);
                            //cout<<"在同一条边上！"<<endl;
                            //getchar();
                        }
                        else if(_p.Ni==u.Ni){
                            distance_min = u.dis + _p.dis;
                            //cout<<"在相邻（顺序）边上！"<<endl;
                        }
                        else if(_p.Nj==u.Ni){
                            distance_min = u.dis + (_p.dist - _p.dis);
                            //cout<<"在相邻（非顺序）边上！"<<endl;
                        }
                        else if(_p.Ni==u.Nj){
                            distance_min = (u.dist-u.dis) + _p.dis;
                            //cout<<"在相邻（非顺序）边上！"<<endl;
                        }

                        else{
                            //for pNi
                            double distance_uNi_pNi = result[i] / WEIGHT_INFLATE_FACTOR;
                            double distance_total_uNi = distance_uNi_pNi + u.dis + _p.dis;
                            double distance_uNj_pNi = result2[i] / WEIGHT_INFLATE_FACTOR;
                            double distance_total_uNj = distance_uNj_pNi + (u.dist-u.dis) + _p.dis;
                            double distance_total_via_pNi = min(distance_total_uNi,distance_total_uNj);
                            //cout<<"dist_total="<<distance_total<<", dist(u.Ni, p.Ni)="<<distance<<",u_dis="<<u_dis<<",p_dis="<<_p.dis<<endl;

                            //for pNj
                            double distance_uNi_pNj =  result[i + 1] / WEIGHT_INFLATE_FACTOR;
                            double distance_total_uNi2 = distance_uNi_pNj + u.dis + (_p.dist-_p.dis);
                            double distance_uNj_pNj =  result2[i + 1] / WEIGHT_INFLATE_FACTOR;
                            double distance_total_uNj2 = distance_uNj_pNj + (u.dist-u.dis) + (_p.dist-_p.dis);
                            double distance_total_via_pNj = min(distance_total_uNi2,distance_total_uNj2);

                            //overall
                            distance_min = min(distance_total_via_pNi, distance_total_via_pNj);
                        }

                        //double score = getGSKScoreu2O_Upper(u_id, a, distance, node_uKey_Map[top.id]);
                        double score = getGSKScoreu2o(u_id, o_id, a, alpha,distance_min);

                        if (score > score_max) {

                            if(resultFinal.size()< K){
                                int poi = top.id;
                                Result tmpRlt(_p.id, _p.Ni, _p.Nj, _p.dis,distance_min, score, _p.keywords);
                                resultFinal.push(tmpRlt);
                                if(resultFinal.size()== K) {
                                    score_k = resultFinal.top().score;
                                    score_max = score_k;
                                }
                            }
                            else{  // update result,  score_k
                                int poi = top.id;
                                Result tmpRlt(_p.id, _p.Ni, _p.Nj, _p.dis,distance_min, score, _p.keywords);
                                resultFinal.push(tmpRlt);
                                resultFinal.pop();
                                score_k = resultFinal.top().score;
                                score_max = max(score_max, score_k);
                                //cout<<"更新score_k="<<score_k<<endl;
                            }
#ifdef  TOPKTRACK
                            cout<<"加入o"<<o_id<<"distance="<<distance_min<<",score="<<score<<endl;
#endif


                        }

                    }
                } else {  // 是叶节点，但与loc不在同一子图
                    set<int> object_set;
                    // 获得该叶节点下包含用户关键词的所有poi
                    for (int term: query_u.keywords) {
                        vector<int> o_list = getObjectTermRelatedEntry(term,top.id);
                        if (o_list.size() > 0) {
                            for (int poi: o_list) { //access disk (leaf_id, term)
                                object_set.insert(poi);
                            }
                        }
                    }
                    for (int id: object_set) {
                        int o_id = id;
                        POI _p = getPOIFromO2UOrgLeafData(o_id);
#ifdef  TOPKTRACK
                        if(o_id==8594){
                            cout<<"find 8594!"<<endl;
                        }
#endif
                        //else cout<<"发现其他"<<endl;

                        int Ni = _p.Ni;
                        int Nj = _p.Nj;
                        double score = 0; int allmin = -1; double distance;
                        double distance_pNi=999999999999; double distance_pNj=999999999999;
                        for (int i = 0; i < GTree[top.id].leafinvlist.size(); i++) {    //获取当前最优叶节点的顶点倒排表
                            int posa = GTree[top.id].leafinvlist[i];
                            int vertex = GTree[top.id].leafnodes[posa];
                            if (vertex == Ni || vertex == Nj) {  //找见poi.Ni, Nj
                                int allmin = -1;
                                for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                    //计算loc 经过 叶节点中第 k个 border 到达 top.id叶节点内第i个顶点的距离
                                    int dis = itm[top.id][k] + GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa];
                                    //cout<<"dis="<<dis<<",itm[top.id][k]="<<itm[top.id][k]<<",边界间距离="<<GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa]<<endl;
                                    if (allmin == -1) {
                                        allmin = dis;
                                    } else {
                                        if (dis < allmin) {
                                            allmin = dis;   // 提取loc 经过所有叶节点top.id中的boder, 而到达第i个顶点的最短距离
                                        }
                                    }

                                }
                                //把所有的allmin
                                distance = allmin / WEIGHT_INFLATE_FACTOR;
                                if(vertex == Ni){
                                    distance_pNi = distance + _p.dis;
                                } else{
                                    distance_pNj = distance +(_p.dist - _p.dis);
                                }


                            }

                        }
                        if(distance_pNj==999999999999){
                            distance_pNj = usrToPOIDistance_phl(u,_p);
                        }

                        if(distance_pNi==999999999999){
                            distance_pNi = usrToPOIDistance_phl(u,_p);
                        }

                        double distance_final = min(distance_pNi,distance_pNj);
                        score = getGSKScoreu2o(u_id, o_id, a, alpha,distance_final);

                        if (score > score_max) {
                            if(resultFinal.size()< K){
                                //int poi = top.id;
                                Result tmpRlt(_p.id, _p.Ni, _p.Nj, _p.dis,distance_final, score, _p.keywords);
                                resultFinal.push(tmpRlt);
                                if(resultFinal.size()== K) {
                                    score_k = resultFinal.top().score;
                                    score_max = score_k;
                                }
                            }
                            else{  // update result,  score_k
                                int poi = top.id;
                                Result tmpRlt(_p.id, _p.Ni, _p.Nj, _p.dis,distance_final, score, _p.keywords);
                                resultFinal.push(tmpRlt);
                                resultFinal.pop();
                                score_k = resultFinal.top().score;
                                score_max = max(score_max, score_k);
                                //cout<<"更新score_k="<<score_k<<endl;
                            }
#ifdef  TOPKTRACK
                            cout<<"加入o"<<_p.id<<"distance="<<distance_final<<",score="<<score<<endl;
#endif
                        }

                    }

                }

            }
            else {    //当前最优为某个中间节点
#ifdef  TOPKTRACK
                cout<<"中间节点出队列，n="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
#endif
                //cout<<"出队列，n="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
                set<int> entry_set;
                // 获得该中间节点下包含用户关键词的所有child
                for (int term: query_u.keywords) {
                    vector<int> entry_list = getObjectTermRelatedEntry(term,top.id);
                    if (entry_list.size() > 0) {
                        for (int child: entry_list) { //access disk (leaf_id, term)
                            entry_set.insert(child);
                        }
                    }

                }
                //将各子条目加入优先队列
                for(int child: entry_set){

                    TreeNode child_Node = getGIMTreeNodeData(child, OnlyO);
                    int son = Nodes[locid].gtreepath[top.lca_pos + 1];
                    // on gtreepath
                    if (child == son) {
                        double score = getGSKScoreu2O_Upper(u_id, a, alpha, min_Dist, child_Node.objectUKeySet);
                        //cout<<"child_score="<<score<<endl;
                        double min_Dist =0;

                        if(score > score_max){
                            QEntry entry(child, false, top.lca_pos + 1, min_Dist,score);
                            Queue.push(entry);

                        }

                    }
                        // brothers
                    else if (GTree[child].father == GTree[son].father) {
                        itm[child].clear();
                        int allmin = -1;

                        for (int j = 0; j < GTree[child].borders.size(); j++) {   //把到当前最优节点的各个child的距离算出，把child加入
                            int _min = -1;
                            posa = GTree[child].up_pos[j];
                            for (int k = 0; k < GTree[son].borders.size(); k++) {   //第k个border
                                posb = GTree[son].up_pos[k];
                                int dis = itm[son][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                                double dist = dis;
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
                        double distance = allmin/WEIGHT_INFLATE_FACTOR;
                        double score = getGSKScoreu2O_Upper(u_id, a,  alpha,distance, child_Node.objectUKeySet);

                        if(score > score_max){
                            QEntry entry(child, false, top.lca_pos, allmin, score);
                            Queue.push(entry);

                        }
                    }
                        // downstream
                    else {
                        itm[child].clear();
                        int allmin = -1;

                        for (int j = 0; j < GTree[child].borders.size(); j++) {
                            int _min = -1;
                            posa = GTree[child].up_pos[j];
                            for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                posb = GTree[top.id].current_pos[k];
                                int dis = itm[top.id][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
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
                        double distance = allmin/WEIGHT_INFLATE_FACTOR;
                        double score = getGSKScoreu2O_Upper(u_id, a, alpha,distance, child_Node.objectUKeySet);

                        if(score > score_max){
                            QEntry entry(child, false, top.lca_pos, allmin, score);
                            Queue.push(entry);
                        }

                    }

                }


            }


        }




    }//while
    if(resultFinal.size()>0){    //倘若resultFinal为空，说明没有任何元素比score_bound更优
        double ss = resultFinal.top().score;
        if(ss> Users[u_id].topkScore_current){
            Users[u_id].topkScore_current = ss;  //记得cache user 的 topk结果
            //Users[u_id].topkScore_Final = ss;  //easy error

        }
    }
    return resultFinal;
}






void UpdateR_disk_group(vector<User> users, vector<int> b_list, int usrLeaf_id, int a, double alpha, double score_bound, int& node_highest,int& node_pos, set<int> KeysetDesire, unordered_map<int, vector<int> > itms[], priority_queue<QEntry>& Queue, double& score_max, bool& eary_terminate){
    //表示当前范围内已全部遍历， 需扩大遍历子图范围
    int posa, posb, _min, dis;
    int root_id =0;
    int node_highest_pre = node_highest;
    int second_highest =0;
    //int node_highest = GTree[node_highest].father;
    int current = node_highest; int current_pre = current;
    set<int> entry_set;
    //cout<<"从"<<current_pre<<"向上寻找"<<endl;
    while(current != root_id){
        // 获得该中间节点下包含用户关键词的所有child
        current_pre = current;
        current = GTree[current].father;
        node_pos--;   //最高层节点层数往上一层
        for (int term: KeysetDesire) {
            vector<int> o_entrylist = getObjectTermRelatedEntry(term,current);
            if (o_entrylist.size() > 0) {
                for (int child: o_entrylist) {
                    entry_set.insert(child);
                    //block_num++;
                }
            }
        }
        //父节点就有关键词
        if(entry_set.size()>0){
            node_highest = current;
            second_highest = current_pre;
            //cout<<"到"<<current<<",发现其有关键词！"<<endl;
            break;
        }
        else{

            //继续往上寻找关键词LCA
            //cout<<"到"<<current<<endl;
        }

    }
    //将扩大的查询范围子图加入
    for(int child: entry_set){
        if(child == current_pre) {
            //cout<<"跳过"<<child<<endl;
            //getchar();
            continue; //跳过已经遍历过的节点
        }
        //int son = Nodes[locid].gtreepath.back()(错误的！)
        //int node_pos = 0;
        int b = GTree[usrLeaf_id].borders[0];
        int son = Nodes[b].gtreepath[node_pos + 1];             //重要！容易出错
        TreeNode child_Node = getGIMTreeNodeData(child, OnlyO);  //gangc
        // on gtreepath
        if (child == son) {
            double min_Dist =0;
            double score = getGSKScoreU2O_Upper(users, a, alpha, min_Dist, child_Node.objectUKeySet);

            if(score > score_bound){
                QEntry entry(child, false, node_pos+1, min_Dist,score);
                Queue.push(entry);
                //cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
            }

        }
            // brothers
        else if (GTree[child].father == GTree[son].father) {
            double distance_overall = 99999999999999999;
            for(int b_th=0;b_th<b_list.size();b_th++){
                itms[b_th][child].clear();
                int allmin = -1;

                for (int j = 0; j < GTree[child].borders.size(); j++) {   //把到当前最优节点的各个child的距离算出，把child加入
                    _min = -1;
                    posa = GTree[child].up_pos[j];
                    for (int k = 0; k < GTree[son].borders.size(); k++) {   //第k个border
                        posb = GTree[son].up_pos[k];
                        dis = itms[b_th][son][k] + GTree[current].mind[posb * GTree[current].union_borders.size() + posa];
                        //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                        double dist = dis;
                        if (_min == -1) {
                            _min = dis;
                        } else {
                            if (dis < _min) {
                                _min = dis;
                            }
                        }
                    }
                    itms[b_th][child].push_back(_min);  //保存到child的距离
                    // update all min
                    if (allmin == -1) {
                        allmin = _min;
                    } else if (_min < allmin) {
                        allmin = _min;
                    }
                }
                double distance = allmin/WEIGHT_INFLATE_FACTOR;
                distance_overall = min(distance,distance_overall);

            }
            double score = getGSKScoreU2O_Upper(users, a, alpha, distance_overall, child_Node.objectUKeySet);
            if(score > score_bound){
                QEntry entry(child, false, node_pos, distance_overall, score);
                Queue.push(entry);
                //cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
            }

        }
            // downstream (也按brother的方式来)
        else {
            double distance_overall = 99999999999999999;
            for(int b_th=0; b_th<b_list.size();b_th++){
                itms[b_th][child].clear();
                int allmin = -1;

                for (int j = 0; j < GTree[child].borders.size(); j++) {
                    _min = -1;
                    posa = GTree[child].up_pos[j];
                    for (int k = 0; k < GTree[second_highest].borders.size(); k++) {  //重要！ 容易出错！
                        posb = GTree[second_highest].current_pos[k];
                        //int right_end = second_highest;
                        dis = itms[b_th][second_highest][k] + GTree[current].mind[posb * GTree[current].union_borders.size() + posa];
                        //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                        if (_min == -1) {
                            _min = dis;
                        } else {
                            if (dis < _min) {
                                _min = dis;
                            }
                        }
                    }
                    itms[b_th][child].push_back(_min);
                    // update all min
                    if (allmin == -1) {
                        allmin = _min;
                    } else if (_min < allmin) {
                        allmin = _min;
                    }
                }
                double distance = allmin/WEIGHT_INFLATE_FACTOR;
                distance_overall = min(distance,distance_overall);
            }

            double score = getGSKScoreU2O_Upper(users, a, alpha, distance_overall,child_Node.objectUKeySet);

            if(score > score_bound){
                QEntry entry(child, false, node_pos, distance_overall, score);
                Queue.push(entry);
                //cout<<"插入e, id="<<entry.id<<",score="<<entry.score<<"Queue.size="<<Queue.size()<<endl;
            }

        }

    }
    //cout<<"update后的 topscore="<<Queue.top().score<<"score_k="<<score_max<<endl;
    if(Queue.size()>0){
        if(score_max > Queue.top().score){
            //priority_queue<QEntry> emptyQ;
            //swap(emptyQ,Queue);     // 把当前队列清空!
            eary_terminate = true;
        }
    }


}

void print_groupTopkResults(vector<priority_queue<Result>> results, vector<int> u_list, int leaf){
    cout<<" usr_leaf" <<leaf<<"的groupTopk 查询处理结果"<<endl;
    if(u_list.size()!=results.size()) {
        cout<<"error!"<<endl;
        exit(-1);
    }
    for(int i=0;i<u_list.size();i++){
        cout<<"user"<<u_list[i]<<"的topk score="<<results[i].top().score<<"具体为："<<endl;
        printTopkResultsQueue(results[i]);
    }
}


vector<priority_queue<Result>> Group_TkGSKQ_bottom2up_verify_disk_slow(int u_leaf, vector<int> u_list, int K, float a, float alpha, double gsk , double score_bound ) {

    //叶节点内各个user的topk 结果集合
    vector<priority_queue<Result>> final_results;
    priority_queue<QEntry> Queue;

    vector<User> users; int u_size = u_list.size();


    //获取user信息
    int user_size = 0;  map<int, int> id_idx_map; map<int, int> idx_id_map;  set<int> userUkeySet; map<int,int> uKey_occuranceCntMap;
    map <int, set<int>> term_DesireUser;
    for (int u_id: u_list){
        User u = getUserFromO2UOrgLeafData(u_id);
        users.push_back(u);
        User query_u = u;
        priority_queue<Result> resultFinal;
        final_results.push_back(resultFinal);
        id_idx_map[u_id] = user_size;
        idx_id_map[user_size] = u_id;
        user_size ++;
        //获取group 用户组的查询关键词并集
        for(int term: u.keywords){
            userUkeySet.insert(term);
            //用户偏好（查询）关键词出现次数记录；
            if(uKey_occuranceCntMap.count(term))
                uKey_occuranceCntMap[term]++;
            else
                uKey_occuranceCntMap[term]=1;
            //所需关键词对应 用户集合列表
            term_DesireUser[term].insert(u_id);
        }


    }
    int border_size = GTree[u_leaf].borders.size();
    unordered_map<int, vector<int> > itm_list[border_size];
    int border_user_distanceMap[border_size][u_size];
    //计算各border到各user的距离
    for(int b_th=0;b_th<border_size;b_th++){
        int border = GTree[u_leaf].borders[b_th];
        vector<int> cands;
        vector<int> result;
        vector<User> tmp_users;
        for(User u: users){
            int Ni = u.Ni;
            int Nj = u.Nj;
            cands.push_back(Ni);
            cands.push_back(Nj);
            tmp_users.push_back(u);
            tmp_users.push_back(u);
        }
        //获取该border到各user的Ni, Nj的距离
        result = dijkstra_candidate(border, cands, Nodes);
        for(int i=0; i<cands.size();i+=2){
            double distance = result[i] / 1000;
            distance = distance + users[i/2].dis;
            double distance2 = result[i + 1] / 1000;
            distance2 = distance2  + users[i/2].dist - users[i].dis;
            double distance_min = min(distance, distance2);  //距离取最小值
            int u_th  = i/2;
            border_user_distanceMap[b_th][u_th] = distance_min;   //这里注意要除以2 -----------------------------------!!!!!!!!!!
        }
    }

    vector<int> border_list = GTree[u_leaf].borders;
    // 计算 叶节点的每个border 到 所有祖先节点border的距离
    for(int b_th = 0; b_th < border_size;b_th ++){
        //User u = users[u_th];
        int locid = border_list[b_th];

        // 计算 叶节点内 某user的loc 到Gtree所有节点border的距离
        unordered_map<int, vector<int> > itm; // intermediate answer, tree node -> array
        itm.clear();
        int cid, posa, posb, _min, dis;
        for (int i = Nodes[locid].gtreepath.size() - 1; i > 0; i--) {   //自底向上
            int tn = Nodes[locid].gtreepath[i];
            itm_list[b_th][tn].clear();

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
        itm_list[b_th] = itm;
    }


#ifdef  TOPKTRACK
    cout<<"完成 loc 到 所有节点的border的距离计算"<<endl;
#endif


    //先把 u_leaf压入最大优先队列中
    int root_id = 0;
    int node_current = u_leaf;
    int node_highest = node_current;
    double score_max = score_bound;  double score_k[u_size];
    for(int i=0;i<u_size;i++){
        score_k[u_size] = score_bound;  //每个user的当前scorek;
    }
    int locid = GTree[u_leaf].borders[0];
    double min_Dist = 0; int node_pos = Nodes[locid].gtreepath.size()-1;


    TreeNode oleaf_same = getGIMTreeNodeData(u_leaf,OnlyO);
    double score_leaf = getGSKScoreU2O_Upper(users, a, alpha, min_Dist, oleaf_same.objectUKeySet);
    QEntry entry(u_leaf,false,node_pos,0, score_leaf);
    Queue.push(entry);
    // 自底向上进行节点遍历迭代
    //cout<<"开始自底向上进行节点遍历迭代"<<endl;
    bool eary_terminate = false;
    set<int> keyset_desire = userUkeySet;
    map<int, vector<Result>> ossociate_ResultMap;//< o_id, 关联user 对其的Result tmpRlt(poi, p.Ni, p.Nj, p.dis, p.dist, top.score, p.keywords);
    map<int, vector<int>> ossociate_UserMap; //< o_id, 与之关联的user
    while(Queue.size() > 0 || node_highest != root_id){//node_highest != root_id){  //Queue.size() > 0 ||
        if(score_max> gsk)
            break; //提前终止条件：已发现k个对象比目标poi有更高的gsk评分

        //当前范围内已全部遍历， 需更新扩大路网图的遍历范围
        if(Queue.empty()){
            //cout<<"当前范围内已全部遍历， 需更新扩大路网图的遍历范围"<<endl; //mojito
            UpdateR_disk_group(users, border_list, u_leaf, a,  alpha, score_max, node_highest, node_pos, keyset_desire, itm_list, Queue, score_max, eary_terminate);
        }
        //当前范围内遍历
        //cout<<"node_highest="<<node_highest<<endl;
        QEntry top = Queue.top();
        double topk_score = top.score;


        if(score_max > topk_score) {   // 当前队列中已无更优元素
            //priority_queue<QEntry> emptyQ;
            //swap(emptyQ,Queue);     // 把当前队列清空!
            if (node_highest!= root_id) { //如果还未到最高层，那就再拓展一波
                UpdateR_disk_group(users, border_list, u_leaf, a,  alpha, score_max, node_highest, node_pos, keyset_desire, itm_list, Queue, score_max,
                                   eary_terminate);
            }
            else {
                if(Queue.size()>0){
                    top = Queue.top();
                    if(score_max > top.score){   //提前终止条件！
                        if(top.isvertex){
                            //cout<<"提前终止条件，出队列，o="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                        }
                        else if(GTree[top.id].isleaf){
                            //cout<<"提前终止条件，出队列，leaf="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                        } else{
                            //cout<<"提前终止条件，出队列，n="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                        }

                        break;
                    }
                }

            }
            //UpdateR(u_id, locid, a, score_bound, node_highest, keywords, itm, Queue, score_max);
            //top = Queue.top();


        }
        //当队列内还有更优的元素时
        if(Queue.size()>0){
            top = Queue.top();
            Queue.pop();
            if (GTree[top.id].isleaf) {           //当前最优为某个叶节点
                //cout<<"出队列，leaf="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
#ifdef  TOPKTRACK
                cout<<"叶节点出队列，leaf="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
#endif
                if (top.id == u_leaf) {  //该叶节点为usr所在的叶节点
                    set<int> object_set;

                    //求O与U集合的查询关键词交集
                    TreeNode oleaf_data = getGIMTreeNodeData(top.id,OnlyO);
                    set<int> interKeySet = obtain_itersection_jins(keyset_desire, oleaf_data.objectUKeySet);


                    // 获得该叶节点下包含用户所需关键词的所有poi
                    for (int term: interKeySet) {
                        vector<int> o_list = getObjectTermRelatedEntry(term,top.id);
                        if (o_list.size() > 0) {
                            for (int poi: o_list) { //access disk (leaf_id, term)
                                object_set.insert(poi);
                            }
                        }
                    }

                    //获取相关poi, 计算各关联用户user对该poi的评分， 并更新关联用户的topk结果表
                    for (int poi : object_set) {
                        POI _p = getPOIFromO2UOrgLeafData(poi);
                        int poi_loc = _p.Ni; int poi_loc2 = _p.Nj;
                        float pDist = _p.dist;
                        float p_dis = _p.dis;
                        //int Nj = _p.Nj;
                        set<int> interKeySet_poi = obtain_itersection_jins(keyset_desire, _p.keywordSet);
                        //获得与该兴趣点相关联的的user
                        set<int> related_user_set;
                        for(int term: interKeySet_poi){
                            for(int u_id: term_DesireUser[term])
                                related_user_set.insert(u_id);
                        }
                        vector<int> cands;  vector<User> related_users;
                        for(int u_id : related_user_set){
                            User u = users[id_idx_map[u_id]];
                            cands.push_back(u.Ni);
                            cands.push_back(u.Nj);
                            //objects.push_back(object_id);
                            related_users.push_back(u);
                            related_users.push_back(u);

                        }
                        //求poi到各个user的叶节点内距离
                        vector<int> result;  vector<int> result2;
                        result = dijkstra_candidate(poi_loc, cands, Nodes);
                        result2 = dijkstra_candidate(poi_loc2, cands, Nodes);

                        //将poi(顶点)更新各个关联user的topk结果集
                        double score_overall_min = 99999999999; bool flag = false;
                        for (int i = 0; i < cands.size(); i += 2) {
                            User  _u = related_users[i/2];
                            int u_id = _u.id;

                            double distanceII = result[i] / 1000;
                            distanceII = distanceII + _u.dis + p_dis;
                            double distanceIJ = result[i + 1] / 1000;
                            distanceIJ = distanceIJ + (_u.dist-_u.dis) + p_dis;
                            double distance_min = min(distanceII, distanceIJ);

                            double distanceJI = result2[i] / 1000;
                            distanceJI = distanceJI + _u.dis + (pDist-p_dis);
                            double distanceJJ = result[i + 1] / 1000;
                            distanceJJ = distanceJJ + (_u.dist-_u.dis) + (pDist-p_dis);
                            double distance_min2 = min(distanceJI, distanceJJ);

                            distance_min = min(distance_min, distance_min2);

                            //double score = getGSKScoreu2O_Upper(u_id, a, distance, node_uKey_Map[top.id]);
                            double score = getGSKScoreu2o(u_id, _p.id, a, alpha,distance_min);
                            //用该poi填充关联用户user的当前topk结果
                            int u_th = id_idx_map[u_id];
                            if(final_results[u_th].size()< K){
                                Result tmpRlt(poi, _p.Ni, _p.Nj, _p.dis, _p.dist, score, _p.keywords);
                                final_results[u_th].push(tmpRlt);
                                if(final_results[u_th].size()== K) {
                                    score_k[u_th] = final_results[u_th].top().score;

                                }
                            }
                            else{  // 更新update result,  score_k
                                Result tmpRlt(poi, _p.Ni, _p.Nj, _p.dis, _p.dist, score, _p.keywords);
                                final_results[u_th].push(tmpRlt);
                                final_results[u_th].pop();
                                score_k[u_th] = final_results[u_th].top().score;  //mojito
                                //cout<<"更新score_k="<<score_k<<endl;
                            }
                        }


                    }
                    for(int u_th=0;u_th<user_size;u_th++){
                        score_max = min(score_max, score_k[u_th]); //全局了
                    }


                } else {  // 是叶节点，但与loc不在同一子图
                    set<int> object_set;
                    // 获得该叶节点下包含用户关键词的所有poi
                    for (int term: userUkeySet) {
                        vector<int> o_list = getObjectTermRelatedEntry(term,top.id);
                        if (o_list.size() > 0) {
                            for (int poi: o_list) { //access disk (leaf_id, term)
                                object_set.insert(poi);
                            }
                        }
                    }
                    //获取相关poi, 计算各关联用户user对该poi的评分， 并更新关联用户的topk结果表
                    map<int, float> borderObject_DisMap[border_size];   //[b_th][poi_id] = distance
                    for (int o_id: object_set) {
                        POI _p = getPOIFromO2UOrgLeafData(o_id);
                        //else cout<<"发现其他"<<endl;
                        int Ni = _p.Ni;
                        int Nj = _p.Nj;
                        set<int> interKeySet_poi = obtain_itersection_jins(keyset_desire, _p.keywordSet);
                        //获得与该兴趣点相关联的的user
                        set<int> related_user_set;
                        for(int term: interKeySet_poi){
                            for(int u_id: term_DesireUser[term])
                                related_user_set.insert(u_id);
                        }
                        //获取该poi到所有user_leaf 中 border的距离
                        for(int b_th =0; b_th<border_size;b_th++){
                            int border = border_list[b_th];
                            double score = 0; int allmin = -1; double distanceI; double distanceJ;
                            for (int i = 0; i < GTree[top.id].leafinvlist.size(); i++) {    //获取当前最优叶节点（object所在）的顶点倒排表
                                int posa = GTree[top.id].leafinvlist[i];
                                int vertex = GTree[top.id].leafnodes[posa];
                                if (vertex == Ni) {
                                    allmin = -1;
                                    for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                        //计算loc 经过 叶节点中第 k个 border 到达 top.id叶节点内第i个顶点的距离
                                        double dis = itm_list[b_th][top.id][k] + GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa];
                                        //cout<<"dis="<<dis<<",itm[top.id][k]="<<itm[top.id][k]<<",边界间距离="<<GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa]<<endl;
                                        if (allmin == -1) {
                                            allmin = dis;
                                        } else {
                                            if (dis < allmin) {
                                                allmin = dis;   // 提取loc 经过所有叶节点top.id中的boder, 而到达第i个顶点的最短距离
                                            }
                                        }

                                    }
                                    //把所有的allmin
                                    distanceI = allmin / 1000;
                                    distanceI = distanceI + _p.dis;
                                    //score = max(score, getGSKScoreu2o(u_id, o_id, a, alpha,distance));

                                }
                                else if (vertex == Nj) {
                                    allmin = -1;
                                    for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                        //计算loc 经过 叶节点中第 k个 border 到达 top.id叶节点内第i个顶点的距离
                                        double dis = itm_list[b_th][top.id][k] + GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa];
                                        //cout<<"dis="<<dis<<",itm[top.id][k]="<<itm[top.id][k]<<",边界间距离="<<GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa]<<endl;
                                        if (allmin == -1) {
                                            allmin = dis;
                                        } else {
                                            if (dis < allmin) {
                                                allmin = dis;   // 提取loc 经过所有叶节点top.id中的boder, 而到达第i个顶点的最短距离
                                            }
                                        }

                                    }
                                    //把所有的allmin
                                    distanceJ = allmin / 1000;
                                    distanceJ = distanceJ + (_p.dist - _p.dis) ;

                                }
                                //poi到b_th个border的距离，取经Ni或经Nj的距离的较小值
                                borderObject_DisMap[b_th][_p.id] = min(distanceI,distanceJ);

                            }
                        }

                        //获取（相关user 到usr_leaf border距离）+ （user_leaf border到该poi距离）  的最小值
                        for(int related_u_id: related_user_set){
                            int u_th = id_idx_map[related_u_id];
                            float disu2p = 9999999999999;
                            //获取该u到p的最小距离
                            for(int b_th=0;b_th<border_size;b_th++){
                                float total_dist = border_user_distanceMap[b_th][u_th] + borderObject_DisMap[b_th][o_id];
                                disu2p = min(disu2p, total_dist);
                            }

                            double score = getGSKScoreu2o(related_u_id, _p.id, a, alpha,disu2p);
                            //将p对u的topk结果进行更新
                            if(final_results[u_th].size()< K){
                                Result tmpRlt(_p.id, _p.Ni, _p.Nj, _p.dis, _p.dist, score, _p.keywords);
                                final_results[u_th].push(tmpRlt);
                                if(final_results[u_th].size()== K) {
                                    score_k[u_th] = final_results[u_th].top().score;

                                }
                            }
                            else{  // 更新update result,  score_k
                                Result tmpRlt(_p.id, _p.Ni, _p.Nj, _p.dis, _p.dist, score, _p.keywords);
                                final_results[u_th].push(tmpRlt);
                                final_results[u_th].pop();
                                score_k[u_th] = final_results[u_th].top().score;  //mojito
                                //cout<<"更新score_k="<<score_k<<endl;
                            }

                        }

                    }
                    //所有object更新完毕后
                    for(int u_th=0;u_th<user_size;u_th++){
                        score_max = min(score_max, score_k[u_th]); //全局了
                    }


                }

            }
            else {    //当前最优为某个中间节点
#ifdef  TOPKTRACK
                cout<<"中间节点出队列，n="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
#endif
                //cout<<"出队列，n="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
                set<int> entry_set;
                // 获得该中间节点下包含用户关键词的所有child
                set<int> inter_keyset;
                TreeNode tn = getGIMTreeNodeData(top.id,OnlyO);
                inter_keyset = obtain_itersection_jins(userUkeySet, tn.objectUKeySet);
                for (int term: inter_keyset) {
                    vector<int> entry_list = getObjectTermRelatedEntry(term,top.id);
                    if (entry_list.size() > 0) {
                        for (int child: entry_list) { //access disk (leaf_id, term)
                            entry_set.insert(child);
                        }
                    }

                }
                //将各子条目加入优先队列
                for(int child: entry_set){

                    TreeNode child_Node = getGIMTreeNodeData(child, OnlyO);
                    int son = Nodes[locid].gtreepath[top.lca_pos + 1];
                    // on gtreepath
                    if (child == son) {

                        //cout<<"child_score="<<score<<endl;
                        double min_Dist =0;
                        double score = getGSKScoreU2O_Upper(users, a, alpha, min_Dist, child_Node.objectUKeySet);
                        if(score > score_max){
                            QEntry entry(child, false, top.lca_pos + 1, min_Dist,score);
                            Queue.push(entry);

                        }

                    }
                        // brothers
                    else if (GTree[child].father == GTree[son].father) {
                        int allmin = -1; double distance_overall= 999999999;
                        for(int b_th =0;b_th<border_size;b_th++){
                            itm_list[b_th][child].clear();
                            for (int j = 0; j < GTree[child].borders.size(); j++) {   //把到当前最优节点的各个child的距离算出，把child加入
                                int _min = -1;
                                int posa = GTree[child].up_pos[j];
                                for (int k = 0; k < GTree[son].borders.size(); k++) {   //第k个border
                                    int posb = GTree[son].up_pos[k];
                                    int dis = itm_list[b_th][son][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                    //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                                    double dist = dis;
                                    if (_min == -1) {
                                        _min = dis;
                                    } else {
                                        if (dis < _min) {
                                            _min = dis;
                                        }
                                    }
                                }
                                itm_list[b_th][child].push_back(_min);  //保存到child的距离
                                // update all min
                                if (allmin == -1) {
                                    allmin = _min;
                                } else if (_min < allmin) {
                                    allmin = _min;
                                }
                            }
                            double distance = allmin/WEIGHT_INFLATE_FACTOR;
                            distance_overall = min(distance_overall,distance);
                        }

                        double score = getGSKScoreU2O_Upper(users, a,  alpha, distance_overall, child_Node.objectUKeySet);

                        if(score > score_max){
                            QEntry entry(child, false, top.lca_pos, distance_overall, score);
                            Queue.push(entry);

                        }
                    }
                        // downstream
                    else {
                        int allmin = -1;
                        double distance_overall= 999999999;
                        for(int b_th =0;b_th<border_size;b_th++){
                            itm_list[b_th][child].clear();
                            for (int j = 0; j < GTree[child].borders.size(); j++) {
                                int _min = -1;
                                int posa = GTree[child].up_pos[j];
                                for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                    int posb = GTree[top.id].current_pos[k];
                                    int dis = itm_list[b_th][top.id][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                    //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                                    if (_min == -1) {
                                        _min = dis;
                                    } else {
                                        if (dis < _min) {
                                            _min = dis;
                                        }
                                    }
                                }
                                itm_list[b_th][child].push_back(_min);
                                // update all min
                                if (allmin == -1) {
                                    allmin = _min;
                                } else if (_min < allmin) {
                                    allmin = _min;
                                }
                            }
                            double distance = allmin/WEIGHT_INFLATE_FACTOR;
                            distance_overall = min(distance_overall,distance);
                        }

                        double score = getGSKScoreU2O_Upper(users, a, alpha,distance_overall, child_Node.objectUKeySet);

                        if(score > score_max){
                            QEntry entry(child, false, top.lca_pos, distance_overall, score);
                            Queue.push(entry);
                        }

                    }

                }


            }


        }




    }//while

    return final_results; //resultFinal;
}

vector<priority_queue<Result>> Group_TkGSKQ_bottom2up_verify_disk(int u_leaf, vector<int> u_list, int K, float a, float alpha, double gsk , double score_bound ) {

    //叶节点内各个user的topk 结果集合
    vector<priority_queue<Result>> final_results;
    priority_queue<QEntry> Queue;

    vector<User> users; int u_size = u_list.size();


    //获取user信息
    int user_size = 0;  map<int, int> id_idx_map; map<int, int> idx_id_map;  set<int> userUkeySet; map<int,int> uKey_occuranceCntMap;
    map <int, set<int>> term_DesireUser; map<int, float> user_BestPrefScore;  //<id, best_score>
    float group_BestPrefScore_max=0;
    for (int u_id: u_list){
        User u = getUserFromO2UOrgLeafData(u_id);
        users.push_back(u);
        User query_u = u;
        priority_queue<Result> resultFinal;
        final_results.push_back(resultFinal);
        id_idx_map[u_id] = user_size;
        idx_id_map[user_size] = u_id;
        user_size ++;
        //获取group 用户组的查询关键词并集
        for(int term: u.keywords){
            userUkeySet.insert(term);
            //用户偏好（查询）关键词出现次数记录；
            if(uKey_occuranceCntMap.count(term))
                uKey_occuranceCntMap[term]++;
            else
                uKey_occuranceCntMap[term]=1;
            //所需关键词对应 用户集合列表
            term_DesireUser[term].insert(u_id);
        }
        //计算user的最佳 textual-social 评分
        float textual_best = textRelevance(u.keywords,u.keywords);
        float social_best = 1+a;
        float textual_social_best = alpha* social_best + (1-alpha)*textual_best;
        user_BestPrefScore[u_id] = textual_social_best;
        group_BestPrefScore_max = max(textual_social_best, group_BestPrefScore_max);

    }
    int border_size = GTree[u_leaf].borders.size();
    cout<<"border_size="<<border_size<<endl;
    unordered_map<int, vector<int> > itm_list[border_size];
    int border_user_distanceMap[border_size][u_size];
    TEST_START
    //计算各border到各user的距离
    for(int b_th=0;b_th<border_size;b_th++){
        int border = GTree[u_leaf].borders[b_th];
        vector<int> cands;
        vector<int> result;
        vector<User> tmp_users;
        for(int u_th =0; u_th<users.size();u_th++){
            User u = users[u_th];
            int Ni = u.Ni;
            int Nj = u.Nj;
            cands.push_back(Ni);
            cands.push_back(Nj);
            tmp_users.push_back(u);
            tmp_users.push_back(u);
            /*int tn = u_leaf; int locid = border;
            int posaNi = lower_bound(GTree[tn].leafnodes.begin(), GTree[tn].leafnodes.end(), locid) -
                   GTree[tn].leafnodes.begin();
            int posaNj = lower_bound(GTree[tn].leafnodes.begin(), GTree[tn].leafnodes.end(), locid) -
                         GTree[tn].leafnodes.begin();
            int distanceNi   = GTree[tn].mind[b_th * GTree[tn].leafnodes.size() + posaNi] / WEIGHT_INFLATE_FACTOR + u.dis;
            int distanceNj   = GTree[tn].mind[b_th * GTree[tn].leafnodes.size() + posaNj] / WEIGHT_INFLATE_FACTOR + (u.dist-u.dis);
            int distance_u = min(distanceNi,distanceNj);
            border_user_distanceMap[b_th][u_th] = distance_u;*/
        }
        //获取该border到各user的Ni, Nj的距离
        result = dijkstra_candidate(border, cands, Nodes);
        for(int i=0; i<cands.size();i+=2){
            double distance = result[i] / 1000;
            distance = distance + users[i/2].dis;
            double distance2 = result[i + 1] / 1000;
            distance2 = distance2  + users[i/2].dist - users[i].dis;
            double distance_min = min(distance, distance2);  //距离取最小值
            int u_th  = i/2;
            border_user_distanceMap[b_th][u_th] = distance_min;   //这里注意要除以2 -----------------------------------!!!!!!!!!!
        }
        //这里注意要除以2 ---
    }
    TEST_END
    TEST_DURA_PRINT("计算各border到各user的距离")


    TEST_START
    vector<int> border_list = GTree[u_leaf].borders;
    // 计算 叶节点的每个border 到 所有祖先节点border的距离
    for(int b_th = 0; b_th < border_size;b_th ++){
        //User u = users[u_th];
        int locid = border_list[b_th];

        // 计算 叶节点内 某user的loc 到Gtree所有节点border的距离
        unordered_map<int, vector<int> > itm; // intermediate answer, tree node -> array
        itm.clear();
        int cid, posa, posb, _min, dis;
        for (int i = Nodes[locid].gtreepath.size() - 1; i > 0; i--) {   //自底向上
            int tn = Nodes[locid].gtreepath[i];
            itm_list[b_th][tn].clear();

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
        itm_list[b_th] = itm;
    }
    TEST_END
    TEST_DURA_PRINT("计算 叶节点的每个border 到 所有祖先节点border的距离")

#ifdef  TOPKTRACK
    cout<<"完成 loc 到 所有节点的border的距离计算"<<endl;
#endif


    //先把 u_leaf压入最大优先队列中
    int root_id = 0;
    int node_current = u_leaf;
    int node_highest = node_current;
    double score_max = score_bound;  double score_k[u_size];
    for(int i=0;i<u_size;i++){
        score_k[u_size] = score_bound;  //每个user的当前scorek;
    }
    int locid = GTree[u_leaf].borders[0];
    double min_Dist = 0; int node_pos = Nodes[locid].gtreepath.size()-1;


    TreeNode oleaf_same = getGIMTreeNodeData(u_leaf,OnlyO);
    double score_leaf = getGSKScoreU2O_Upper(users, a, alpha, min_Dist, oleaf_same.objectUKeySet);
    QEntry entry(u_leaf,false,node_pos,0, score_leaf);
    Queue.push(entry);
    // 自底向上进行节点遍历迭代
    //cout<<"开始自底向上进行节点遍历迭代"<<endl;
    bool eary_terminate = false;
    set<int> keyset_desire = userUkeySet;
    map<int, vector<Result>> ossociate_ResultMap;//< o_id, 关联user 对其的Result tmpRlt(poi, p.Ni, p.Nj, p.dis, p.dist, top.score, p.keywords);
    map<int, vector<int>> ossociate_UserMap; //< o_id, 与之关联的user
    TEST_START
    while(Queue.size() > 0 || node_highest != root_id){//node_highest != root_id){  //Queue.size() > 0 ||
        if(score_max> gsk)
            break; //提前终止条件：已发现k个对象比目标poi有更高的gsk评分

        //当前范围内已全部遍历， 需更新扩大路网图的遍历范围
        if(Queue.empty()){
            //cout<<"当前范围内已全部遍历， 需更新扩大路网图的遍历范围"<<endl; //mojito
            UpdateR_disk_group(users, border_list, u_leaf, a,  alpha, score_max, node_highest, node_pos, keyset_desire, itm_list, Queue, score_max, eary_terminate);
        }
        //当前范围内遍历
        //cout<<"node_highest="<<node_highest<<endl;
        QEntry top = Queue.top();
        double topk_score = top.score;


        if(score_max > topk_score) {   // 当前队列中已无更优元素
            //priority_queue<QEntry> emptyQ;
            //swap(emptyQ,Queue);     // 把当前队列清空!
            if (node_highest!= root_id) { //如果还未到最高层，那就再拓展一波
                UpdateR_disk_group(users, border_list, u_leaf, a,  alpha, score_max, node_highest, node_pos, keyset_desire, itm_list, Queue, score_max,
                                   eary_terminate);
            }
            else {
                if(Queue.size()>0){
                    top = Queue.top();
                    if(score_max > top.score){   //提前终止条件！
                        if(top.isvertex){
                            //cout<<"提前终止条件，出队列，o="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                        }
                        else if(GTree[top.id].isleaf){
                            //cout<<"提前终止条件，出队列，leaf="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                        } else{
                            //cout<<"提前终止条件，出队列，n="<<top.id<<",top-score="<<top.score<<",Queue.size="<<Queue.size()<<endl;
                        }

                        break;
                    }
                }

            }
            //UpdateR(u_id, locid, a, score_bound, node_highest, keywords, itm, Queue, score_max);
            //top = Queue.top();


        }
        //当队列内还有更优的元素时
        if(Queue.size()>0){
            top = Queue.top();
            Queue.pop();
            if (GTree[top.id].isleaf) {           //当前最优为某个叶节点
                //cout<<"出队列，leaf="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
#ifdef  TOPKTRACK
                cout<<"叶节点出队列，leaf="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
#endif
                if (top.id == u_leaf) {  //该叶节点为usr所在的叶节点
                    set<int> object_set;

                    //求O与U集合的查询关键词交集
                    TreeNode oleaf_data = getGIMTreeNodeData(top.id,OnlyO);
                    set<int> interKeySet = obtain_itersection_jins(keyset_desire, oleaf_data.objectUKeySet);


                    // 获得该叶节点下包含用户所需关键词的所有poi
                    for (int term: interKeySet) {
                        vector<int> o_list = getObjectTermRelatedEntry(term,top.id);
                        if (o_list.size() > 0) {
                            for (int poi: o_list) { //access disk (leaf_id, term)
                                object_set.insert(poi);
                            }
                        }
                    }

                    //获取相关poi, 计算各关联用户user对该poi的评分， 并更新关联用户的topk结果表
                    for (int poi : object_set) {
                        POI _p = getPOIFromO2UOrgLeafData(poi);
                        int poi_loc = _p.Ni; int poi_loc2 = _p.Nj;
                        float pDist = _p.dist;
                        float p_dis = _p.dis;
                        //int Nj = _p.Nj;
                        set<int> interKeySet_poi = obtain_itersection_jins(keyset_desire, _p.keywordSet);
                        //获得与该兴趣点相关联的的user
                        set<int> related_user_set;
                        for(int term: interKeySet_poi){
                            for(int u_id: term_DesireUser[term])
                                related_user_set.insert(u_id);
                        }
                        vector<int> cands;  vector<User> related_users;
                        for(int u_id : related_user_set){
                            User u = users[id_idx_map[u_id]];
                            cands.push_back(u.Ni);
                            cands.push_back(u.Nj);
                            //objects.push_back(object_id);
                            related_users.push_back(u);
                            related_users.push_back(u);

                        }
                        //求poi到各个user的叶节点内距离
                        vector<int> result;  vector<int> result2;
                        result = dijkstra_candidate(poi_loc, cands, Nodes);
                        result2 = dijkstra_candidate(poi_loc2, cands, Nodes);

                        //将poi(顶点)更新各个关联user的topk结果集
                        double score_overall_min = 99999999999; bool flag = false;
                        for (int i = 0; i < cands.size(); i += 2) {
                            User  _u = related_users[i/2];
                            int u_id = _u.id;

                            double distanceII = result[i] / 1000;
                            distanceII = distanceII + _u.dis + p_dis;
                            double distanceIJ = result[i + 1] / 1000;
                            distanceIJ = distanceIJ + (_u.dist-_u.dis) + p_dis;
                            double distance_min = min(distanceII, distanceIJ);

                            double distanceJI = result2[i] / 1000;
                            distanceJI = distanceJI + _u.dis + (pDist-p_dis);
                            double distanceJJ = result[i + 1] / 1000;
                            distanceJJ = distanceJJ + (_u.dist-_u.dis) + (pDist-p_dis);
                            double distance_min2 = min(distanceJI, distanceJJ);

                            distance_min = min(distance_min, distance_min2);

                            //double score = getGSKScoreu2O_Upper(u_id, a, distance, node_uKey_Map[top.id]);
                            double score = getGSKScoreu2o(u_id, _p.id, a, alpha,distance_min);
                            //用该poi填充关联用户user的当前topk结果
                            int u_th = id_idx_map[u_id];
                            if(final_results[u_th].size()< K){
                                Result tmpRlt(poi, _p.Ni, _p.Nj, _p.dis, _p.dist, score, _p.keywords);
                                final_results[u_th].push(tmpRlt);
                                if(final_results[u_th].size()== K) {
                                    score_k[u_th] = final_results[u_th].top().score;

                                }
                            }
                            else{  // 更新update result,  score_k
                                Result tmpRlt(poi, _p.Ni, _p.Nj, _p.dis, _p.dist, score, _p.keywords);
                                final_results[u_th].push(tmpRlt);
                                final_results[u_th].pop();
                                score_k[u_th] = final_results[u_th].top().score;  //mojito
                                //cout<<"更新score_k="<<score_k<<endl;
                            }
                        }


                    }
                    for(int u_th=0;u_th<user_size;u_th++){
                        score_max = min(score_max, score_k[u_th]); //全局了
                    }


                } else {  // 是叶节点，但与loc不在同一子图
                    set<int> object_set;
                    // 获得该叶节点下包含用户关键词的所有poi
                    for (int term: userUkeySet) {
                        vector<int> o_list = getObjectTermRelatedEntry(term,top.id);
                        if (o_list.size() > 0) {
                            for (int poi: o_list) { //access disk (leaf_id, term)
                                object_set.insert(poi);
                            }
                        }
                    }
                    //获取相关poi, 计算各关联用户user对该poi的评分， 并更新关联用户的topk结果表
                    map<int, float> borderObject_DisMap[border_size];   //[b_th][poi_id] = distance
                    for (int o_id: object_set) {
                        POI _p = getPOIFromO2UOrgLeafData(o_id);
                        //else cout<<"发现其他"<<endl;
                        int Ni = _p.Ni;
                        int Nj = _p.Nj;
                        set<int> interKeySet_poi = obtain_itersection_jins(keyset_desire, _p.keywordSet);
                        //获得与该兴趣点相关联的的user
                        set<int> related_user_set;
                        for(int term: interKeySet_poi){
                            for(int u_id: term_DesireUser[term])
                                related_user_set.insert(u_id);
                        }
                        //获取该poi到所有user_leaf 中 border的距离
                        double pToBordersDistance_min = 999999999999999999999999;
                        for(int b_th =0; b_th<border_size;b_th++){
                            int border = border_list[b_th];
                            double score = 0; int allmin = -1; double distanceI; double distanceJ;
                            for (int i = 0; i < GTree[top.id].leafinvlist.size(); i++) {    //获取当前最优叶节点（object所在）的顶点倒排表
                                int posa = GTree[top.id].leafinvlist[i];
                                int vertex = GTree[top.id].leafnodes[posa];
                                if (vertex == Ni) {
                                    allmin = -1;
                                    for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                        //计算loc 经过 叶节点中第 k个 border 到达 top.id叶节点内第i个顶点的距离
                                        double dis = itm_list[b_th][top.id][k] + GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa];
                                        //cout<<"dis="<<dis<<",itm[top.id][k]="<<itm[top.id][k]<<",边界间距离="<<GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa]<<endl;
                                        if (allmin == -1) {
                                            allmin = dis;
                                        } else {
                                            if (dis < allmin) {
                                                allmin = dis;   // 提取loc 经过所有叶节点top.id中的boder, 而到达第i个顶点的最短距离
                                            }
                                        }

                                    }
                                    //把所有的allmin
                                    distanceI = allmin / 1000;
                                    distanceI = distanceI + _p.dis;
                                    //score = max(score, getGSKScoreu2o(u_id, o_id, a, alpha,distance));

                                }
                                else if (vertex == Nj) {
                                    allmin = -1;
                                    for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                        //计算loc 经过 叶节点中第 k个 border 到达 top.id叶节点内第i个顶点的距离
                                        double dis = itm_list[b_th][top.id][k] + GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa];
                                        //cout<<"dis="<<dis<<",itm[top.id][k]="<<itm[top.id][k]<<",边界间距离="<<GTree[top.id].mind[k * GTree[top.id].leafnodes.size() + posa]<<endl;
                                        if (allmin == -1) {
                                            allmin = dis;
                                        } else {
                                            if (dis < allmin) {
                                                allmin = dis;   // 提取loc 经过所有叶节点top.id中的boder, 而到达第i个顶点的最短距离
                                            }
                                        }

                                    }
                                    //把所有的allmin
                                    distanceJ = allmin / 1000;
                                    distanceJ = distanceJ + (_p.dist - _p.dis) ;

                                }
                                //poi到b_th个border的距离，取经Ni或经Nj的距离的较小值
                                double _distance = min(distanceI,distanceJ);
                                borderObject_DisMap[b_th][_p.id] = _distance;
                                pToBordersDistance_min = min(pToBordersDistance_min, _distance);
                            }
                        }

                        float textual_score_group_best = 9999999;
                        if(group_BestPrefScore_max/(1.0+pToBordersDistance_min)<score_max&&score_max>0)
                            continue;


                        //获取（相关user 到usr_leaf border距离）+ （user_leaf border到该poi距离）  的最小值
                        for(int related_u_id: related_user_set){
                            int u_th = id_idx_map[related_u_id];
                            float disu2p = 9999999999999;
                            //获取该u到p的最小距离
                            for(int b_th=0;b_th<border_size;b_th++){
                                float total_dist = border_user_distanceMap[b_th][u_th] + borderObject_DisMap[b_th][o_id];
                                disu2p = min(disu2p, total_dist);
                            }

                            if(user_BestPrefScore[related_u_id]/(1.0+disu2p)<score_k[u_th]&&score_k[u_th]>0){  //评分不会更优
                                continue;
                            }
                            //POI p_check = getPOIFromO2UOrgLeafData(_p.id);
                            float textual_score = textRelevance(users[u_th].keywords,_p.keywords);
                            float _score_bound = alpha*(1+a)+( 1-alpha)*textual_score;
                            if(_score_bound/(1.0+disu2p)<score_k[u_th]&& score_k[u_th]>0) //评分不会更优
                                continue;


                            double score = getGSKScoreuData2oDataTextual(users[u_th], _p, a, alpha,textual_score, disu2p);
                            //将p对u的topk结果进行更新
                            if(final_results[u_th].size()< K){
                                Result tmpRlt(_p.id, _p.Ni, _p.Nj, _p.dis, _p.dist, score, _p.keywords);
                                final_results[u_th].push(tmpRlt);
                                if(final_results[u_th].size()== K) {
                                    score_k[u_th] = final_results[u_th].top().score;

                                }
                            }
                            else{  // 更新update result,  score_k
                                Result tmpRlt(_p.id, _p.Ni, _p.Nj, _p.dis, _p.dist, score, _p.keywords);
                                final_results[u_th].push(tmpRlt);
                                final_results[u_th].pop();
                                score_k[u_th] = final_results[u_th].top().score;  //mojito
                                //cout<<"更新score_k="<<score_k<<endl;
                            }

                        }

                    }
                    //所有object更新完毕后
                    for(int u_th=0;u_th<user_size;u_th++){
                        score_max = min(score_max, score_k[u_th]); //全局了
                    }


                }

            }
            else {    //当前最优为某个中间节点
#ifdef  TOPKTRACK
                cout<<"中间节点出队列，n="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
#endif
                //cout<<"出队列，n="<<top.id<<",top-score="<<top.score<<"Queue.size="<<Queue.size()<<endl;
                set<int> entry_set;
                // 获得该中间节点下包含用户关键词的所有child
                set<int> inter_keyset;
                TreeNode tn = getGIMTreeNodeData(top.id,OnlyO);
                inter_keyset = obtain_itersection_jins(userUkeySet, tn.objectUKeySet);
                for (int term: inter_keyset) {
                    vector<int> entry_list = getObjectTermRelatedEntry(term,top.id);
                    if (entry_list.size() > 0) {
                        for (int child: entry_list) { //access disk (leaf_id, term)
                            entry_set.insert(child);
                        }
                    }

                }
                //将各子条目加入优先队列
                for(int child: entry_set){

                    TreeNode child_Node = getGIMTreeNodeData(child, OnlyO);
                    int son = Nodes[locid].gtreepath[top.lca_pos + 1];
                    // on gtreepath
                    if (child == son) {

                        //cout<<"child_score="<<score<<endl;
                        double min_Dist =0;
                        double score = getGSKScoreU2O_Upper(users, a, alpha, min_Dist, child_Node.objectUKeySet);
                        if(score > score_max){
                            QEntry entry(child, false, top.lca_pos + 1, min_Dist,score);
                            Queue.push(entry);

                        }

                    }
                        // brothers
                    else if (GTree[child].father == GTree[son].father) {
                        int allmin = -1; double distance_overall= 999999999;
                        for(int b_th =0;b_th<border_size;b_th++){
                            itm_list[b_th][child].clear();
                            for (int j = 0; j < GTree[child].borders.size(); j++) {   //把到当前最优节点的各个child的距离算出，把child加入
                                int _min = -1;
                                int posa = GTree[child].up_pos[j];
                                for (int k = 0; k < GTree[son].borders.size(); k++) {   //第k个border
                                    int posb = GTree[son].up_pos[k];
                                    int dis = itm_list[b_th][son][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                    //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                                    double dist = dis;
                                    if (_min == -1) {
                                        _min = dis;
                                    } else {
                                        if (dis < _min) {
                                            _min = dis;
                                        }
                                    }
                                }
                                itm_list[b_th][child].push_back(_min);  //保存到child的距离
                                // update all min
                                if (allmin == -1) {
                                    allmin = _min;
                                } else if (_min < allmin) {
                                    allmin = _min;
                                }
                            }
                            double distance = allmin/WEIGHT_INFLATE_FACTOR;
                            distance_overall = min(distance_overall,distance);
                        }

                        double score = getGSKScoreU2O_Upper(users, a,  alpha, distance_overall, child_Node.objectUKeySet);

                        if(score > score_max){
                            QEntry entry(child, false, top.lca_pos, distance_overall, score);
                            Queue.push(entry);

                        }
                    }
                        // downstream
                    else {
                        int allmin = -1;
                        double distance_overall= 999999999;
                        for(int b_th =0;b_th<border_size;b_th++){
                            itm_list[b_th][child].clear();
                            for (int j = 0; j < GTree[child].borders.size(); j++) {
                                int _min = -1;
                                int posa = GTree[child].up_pos[j];
                                for (int k = 0; k < GTree[top.id].borders.size(); k++) {
                                    int posb = GTree[top.id].current_pos[k];
                                    int dis = itm_list[b_th][top.id][k] + GTree[top.id].mind[posb * GTree[top.id].union_borders.size() + posa];
                                    //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                                    if (_min == -1) {
                                        _min = dis;
                                    } else {
                                        if (dis < _min) {
                                            _min = dis;
                                        }
                                    }
                                }
                                itm_list[b_th][child].push_back(_min);
                                // update all min
                                if (allmin == -1) {
                                    allmin = _min;
                                } else if (_min < allmin) {
                                    allmin = _min;
                                }
                            }
                            double distance = allmin/WEIGHT_INFLATE_FACTOR;
                            distance_overall = min(distance_overall,distance);
                        }

                        double score = getGSKScoreU2O_Upper(users, a, alpha,distance_overall, child_Node.objectUKeySet);

                        if(score > score_max){
                            QEntry entry(child, false, top.lca_pos, distance_overall, score);
                            Queue.push(entry);
                        }

                    }

                }


            }


        }




    }//while
    TEST_END
    TEST_DURA_PRINT("while循环")

    return final_results; //resultFinal;
}



LCLResult getLeafNode_Topkbottom2up(int leaf_id, int term_id, int Qk,float a, float alpha,float score_bound) {

    int nodeCnt = 0;
    float Max_D = 25000;

    bool vis[VertexNum];
    bool rearch[VertexNum];
    bool poi_rearch[poi_num];
    float pointval[VertexNum];

    memset(vis, false, sizeof(vis));
    memset(rearch, false, sizeof(rearch));
    float maxDisValue = 9999999.9;
    //memset(pointval, maxDisValue, sizeof(float)*sizeof(pointval));
    //for (int i = 0; i < VertexNum; ++i)
    //pointval[i] = maxNum;

    priority_queue<nodeDistType> tmpPoint;
    priority_queue<LCLEntry> lcl;
    LCLResult result;
    float maxDisInsideLeaf = MaxDisBound[leaf_id][leaf_id];
    int leaf_border = GTree[leaf_id].borders[0];


    int qNi = leaf_border;
    float qDis = maxDisInsideLeaf;
    //cout<<maxDisInsideLeaf<<endl;
    vector<int> qKey;
    qKey.push_back(term_id);

    // expand the edge from Ni, Nj
    pointval[qNi] = qDis;

    nodeDistType tmp1(qNi, qDis);
    tmpPoint.push(tmp1);
    double  optimal_text_score = textRelevance(qKey, qKey);
    while (!tmpPoint.empty()) {
        nodeDistType tmp = tmpPoint.top();
        tmpPoint.pop();
        int tmpNi = tmp.ID;
        if (vis[tmpNi] == true)
            continue;
        nodeCnt++;
        vis[tmpNi] = true;

        vector<int> tmpAdj = adjList[tmpNi];
        for (size_t i = 0; i < tmpAdj.size(); i++) {  // expand each edge
            int tmpNj = tmpAdj[i];
            if (vis[tmpNj])
                continue;

            edge e = EdgeMap[getKey(tmpNi, tmpNj)];
            int Ni = e.Ni;
            float edgeDist = e.dist;
            vector<POI> poiList = e.pts;
            for (size_t j = 0; j < poiList.size(); j++) { // check each poi on the edge
                POI tmpPOI = poiList[j];

                //tmpPOI = POIs[60];

                if (tmpPOI.category == 0)
                    continue;

                int pId = tmpPOI.id;

                float pDis = tmpPOI.dis;
                vector<int> pKey = tmpPOI.keywords;
                //这里待加入社交

                float distance = -1;

                bool relevance = textCover(term_id, pKey);

                //必须文本相关
                if (!relevance) continue;

                double  optimal_text_score = textRelevance(qKey, qKey);

                // user and poi are on the same edge
                if ((tmpNi == qNi) || (tmpNj == qNi)) {
                    if (qDis > pDis)
                        distance = qDis - pDis;
                    else
                        distance = pDis - qDis;
                }
                    //user and poi are on the different edge
                else {
                    if (tmpNi == e.Ni)  //同一个点
                        distance = pointval[tmpNi] + pDis;
                    else //边的另一个点
                        distance = pointval[tmpNi] + edgeDist - pDis;

                }

                double simD = distance;

                float simT = 0;
                simT = tfIDF_term(term_id);

                //double pow(double x,double y) //求x的y次方
                float simS;
                //double social_score = getSocialScore(qId,pId);
                double social_score = 0;
                /*if(social_score == 0)
                    simS = 1.0;
                else
                    simS = 1.0 + pow(a,social_score);*/
                if(social_score > 0.0)
                    simS = 1.0 + pow(a,social_score);
                else
                    simS = 1.0;
                double social_textual_score = alpha*simS + (1-alpha)*simT;
                double score = social_textual_score / (1+simD);

                LCLEntry le(pId, score, 1);
                //cout<<"qId="<<qId<<"pId="<<pId<<endl;
                if (lcl.size() < Qk) {
                    lcl.push(le);
                    //cout<<"simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;
                }
                else {
                    LCLEntry _tmpNode = lcl.top();
                    float tmpScore = _tmpNode.score;

                    if (tmpScore < score) {
                        lcl.pop();
                        lcl.push(_tmpNode);
                        //cout<<"simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;
                    }
                    if(lcl.top().score > score_bound){//当前已有更优的k个结果，拓展可提前终止

                        result.LCL = lcl;
                        result.topscore = lcl.top().score;
                        return result;  //注意此时的score只是current score
                    }
                }
            }
            //首次拓展当前拓展顶点
            if(rearch[tmpNj]==false){
                pointval[tmpNj] = maxDisValue;
            }
            if (e.dist + tmp.dist < pointval[tmpNj]) {
                pointval[tmpNj] = e.dist + tmp.dist;
                nodeDistType tmp2(tmpNj, edgeDist + tmp.dist);
                //cout << tmp2.num << ' ' << tmp2.dist << endl;
                tmpPoint.push(tmp2);
                rearch[tmpNj]=true;
            }
        }
        //路网边拓展提前终止条件
        if (lcl.size() == Qk) {
            float tmpScore = lcl.top().score;
            float currentDistance = pointval[tmpNi];
            double optimal_social_textual_score = alpha*(a+1.0)+(1-alpha)*optimal_text_score;
            double bound = optimal_social_textual_score /tmpScore - 1;               //提前拓展终止条件出错
            if (currentDistance > bound){
                //cout<<"当前拓展终止"<<endl;
                break;
            }

        }
    }
    result.LCL = lcl;
    result.topscore = lcl.top().score;
    return result;
}


/*
 * 基于Dijastra， 1. 求top-k结果  2. 验证是否在top-k结果中
 */
/*************************1**************************/
//基于内存
TopkQueryResult topkSDijkstra_memory(query Query, int Qk,float a, float alpha) {

    int nodeCnt = 0;
    float Max_D = 25000;
    double score_bound = 1000000000000000;
    bool vis[VertexNum];
    bool rearch[VertexNum];
    bool poi_rearch[poi_num];
    float pointval[VertexNum];

    memset(vis, false, sizeof(vis));
    memset(rearch, false, sizeof(rearch));
    memset(poi_rearch, false, sizeof(poi_rearch));
    float maxDisValue = 9999999.9;
    //memset(pointval, maxDisValue, sizeof(float)*sizeof(pointval));
    //for (int i = 0; i < VertexNum; ++i)
    //pointval[i] = maxNum;
    //TopkQueryResult  topkCurrentResult2(0.5);
    //return  topkCurrentResult2;

    priority_queue<nodeDistType> tmpPoint;
    priority_queue<Result> resultFinal;

    int qId = Query.id;
    int qNi = Query.Ni;
    int qNj = Query.Nj;
    float qDist = Query.dist;
    float qDis = Query.dis;
    vector<int> qKey = Query.keywords;
    double  optimal_text_score = textRelevance(qKey, qKey);

    // expand the edge from Ni, Nj
    pointval[qNi] = qDis;
    pointval[qNj] = qDist - qDis;

    nodeDistType tmp1(qNi, qDis);
    tmpPoint.push(tmp1);

    nodeDistType tmp2(qNj, qDist - qDis);
    tmp2.dist = qDist - qDis;
    tmpPoint.push(tmp2);

    while (!tmpPoint.empty()) {
        nodeDistType tmp = tmpPoint.top();
        tmpPoint.pop();
        int tmpNi = tmp.ID;
        if (vis[tmpNi] == true)
            continue;
        nodeCnt++;
        vis[tmpNi] = true;

        vector<int> tmpAdj = adjList[tmpNi];
        for (size_t i = 0; i < tmpAdj.size(); i++) {  // expand each edge
            int tmpNj = tmpAdj[i];
            if (vis[tmpNj])
                continue;

            edge e = EdgeMap[getKey(tmpNi, tmpNj)];
            int Ni = e.Ni;
            float edgeDist = e.dist;
            vector<POI> poiList = e.pts;
            for (size_t j = 0; j < poiList.size(); j++) { // check each poi on the edge
                POI tmpPOI = poiList[j];
                //tmpPOI = POIs[60];
                if (tmpPOI.category == 0)
                    continue;

                int pId = tmpPOI.id;

                float pDis = tmpPOI.dis;
                vector<int> pKey = tmpPOI.keywords;
                //这里待加入社交

                float distance = -1;

                // user and poi are on the same edge
                if ((tmpNi == qNi && tmpNj == qNj) || (tmpNi == qNj && tmpNj == qNi)) {
                    if (qDis > pDis)
                        distance = qDis - pDis;
                    else
                        distance = pDis - qDis;

                }
                    //user and poi are on the different edge
                else {
                    if (tmpNi == e.Ni) { //同一个点 此时Ni的邻接点Nj是兴趣点所在边的Ni
                        distance = pointval[tmpNi] + pDis;   //edgeDist - pDis;//
                        //cout<<edgeDist - pDis;
                        //cout<<"同一个点"<<endl;
                    }
                    else{ //边的另一个点
                        distance = pointval[tmpNi] + edgeDist - pDis;
                        //cout<<(edgeDist - pDis);
                        //cout<<"边的另外一个点"<<endl;
                    }

                }

                double simD = distance;
                //记录拓展到兴趣点的当前最短距离
                if(poi_rearch[pId] == false){
                    poi_rearch[pId] = true;
                    poiDistance[pId] = distance;
                }else if(distance<poiDistance[pId])
                    poiDistance[pId] = distance;


                //float simT = tfIDF(qKey, pKey);
                double simT = textRelevance(qKey, pKey);
                if (simT==0) continue;
                //必须文本相关

                //double d2 = usrToPOIDistance(Users[qId],POIs[pId]); //getDistance(Users[qId],POIs[pId]);//
                //cout<<"o"<<pId<<"simD="<<poiDistance[pId]<<",u2o="<<d2<<endl;

                //double pow(double x,double y) //求x的y次方
                double simS;
                double social_score = calcSocialScore(qId,pId);

                if(social_score > 0)
                    simS = 1.0 + pow(a,social_score);
                else
                    simS = 1.0;

                double score = simT * simS / (1+simD);
                /*if(pId == 57){
                    cout<<"find o57! simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;
                }*/

                Result tmpRlt(pId, tmpNi, tmpNj, pDis, distance, score, pKey);
                //cout<<"qId="<<qId<<"pId="<<pId<<endl;
                if (resultFinal.size() < Qk) {
                    resultFinal.push(tmpRlt);
                    //
                }
                else {
                    Result _tmpNode = resultFinal.top();
                    double tmpScore = _tmpNode.score;

                    if (tmpScore < score) {
                        resultFinal.pop();
                        resultFinal.push(tmpRlt);
                        //cout<<"simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;
                    }
                    if(resultFinal.top().score > score_bound){//当前已有更优的k个结果，拓展可提前终止
                        if(resultFinal.top().score > Users[qId].topkScore_current)
                            Users[qId].topkScore_current = resultFinal.top().score;
                        TopkQueryResult  topkCurrentResult(resultFinal.top().score,resultFinal);
                        return topkCurrentResult;  //注意此时的score只是current score
                    }
                }
            }
            //首次拓展当前拓展顶点
            if(rearch[tmpNj]==false){
                pointval[tmpNj] = maxDisValue;
            }
            if (e.dist + tmp.dist < pointval[tmpNj]) {
                //cout<<e.dist<<endl;
                pointval[tmpNj] = e.dist + tmp.dist;
                nodeDistType tmp2(tmpNj, edgeDist + tmp.dist);
                //cout << tmp2.num << ' ' << tmp2.dist << endl;
                tmpPoint.push(tmp2);
                rearch[tmpNj]=true;
            }
        }
        //路网边拓展提前终止条件
        if (resultFinal.size() == Qk) {
            double tmpScore = resultFinal.top().score;
            double currentDistance = pointval[tmpNi];

            double optimal_social_textual_score = alpha*(a+1.0)+(1-alpha)*optimal_text_score;
            double bound = optimal_social_textual_score /tmpScore - 1;               //提前拓展终止条件出错
            if (currentDistance > bound){
                //cout<<"当前拓展终止"<<endl;
                break;
            }

        }
    }
    Result topkNode = resultFinal.top();
    double topkScore = topkNode.score;
    Users[qId].topkScore_current = max(Users[qId].topkScore_current,topkScore);
    Users[qId].topkScore_Final = topkScore;
    TopkQueryResult  topkCurrentResult(resultFinal.top().score,resultFinal);
    return topkCurrentResult;  //注意此时的score只是current score
}

//基于外存
TopkQueryResult topkSDijkstra(query Query, int Qk,float a, float alpha) {
    bool vis[VertexNum];
    bool rearch[VertexNum];
    bool poi_rearch[poi_num];
    float pointval[VertexNum];

    int nodeCnt = 0;
    float Max_D = 25000;
    double score_bound = 1000000000000000;
    memset(vis, false, sizeof(vis));
    memset(rearch, false, sizeof(rearch));
    memset(poi_rearch, false, sizeof(poi_rearch));
    float maxDisValue = 9999999.9;

    priority_queue<nodeDistType> tmpPoint;
    priority_queue<Result> resultFinal;

    int qId = Query.id;
    int qNi = Query.Ni;
    int qNj = Query.Nj;
    float qDist = Query.dist;
    float qDis = Query.dis;
    vector<int> qKey = Query.keywords;
    double  optimal_text_score = textRelevance(qKey, qKey);

    // expand the edge from Ni, Nj
    pointval[qNi] = qDis;
    pointval[qNj] = qDist - qDis;

    nodeDistType tmp1(qNi, qDis);
    tmpPoint.push(tmp1);

    nodeDistType tmp2(qNj, qDist - qDis);
    tmp2.dist = qDist - qDis;
    tmpPoint.push(tmp2);


    while (!tmpPoint.empty()) {
        nodeDistType tmp = tmpPoint.top();
        tmpPoint.pop();
        int tmpNi = tmp.ID;
        if (vis[tmpNi] == true)
            continue;
        nodeCnt++;
        vis[tmpNi] = true;

        int AdjGrpAddr,AdjListSize;
        vector<int> tmpAdj = adjList[tmpNi];
        /*1. 取得地址，   2. 根据地址从文件中取得数据*/
        //---------------取得顶点tmpNi的数据地址：AdjGrpAddr
        AdjGrpAddr=getAdjListGrpAddr(tmpNi);
        //cout<<"AdjGrpAddr="<<AdjGrpAddr<<endl;
        bool testanswer=true;
        //---------------取得顶点tmpNi的邻接点个数：AdjListSize
        getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
        //cout<<"开始N"<<tmpNi<<"的拓展，邻接点size= "<<AdjListSize<<endl;
        //for (size_t i = 0; i < tmpAdj.size(); i++) {  // expand each edge
        for (size_t i = 0; i < AdjListSize; i++) {  // expand each edge
            //int tmpNj = tmpAdj[i];
            //------------取得邻接顶点tmpNj的id：tmpNj
            int tmpNj;
            getVarE(ADJNODE_A,Ref(tmpNj),AdjGrpAddr,i);  //读到的Nj不对
            //cout<<"读到Nj"<<tmpNj<<endl;
            if (vis[tmpNj])
                continue;

            edge e = EdgeMap[getKey(tmpNi, tmpNj)];
            //cout<<"e.dist="<<e.dist<<endl;
            int Ni = e.Ni;
            //float edgeDist = e.dist;
            float edgeDist; getVarE(DIST_A,Ref(edgeDist),AdjGrpAddr,i);  //读取边的长度(成功)

            vector<POI> poiList = e.pts;

            int PoiGrpAddress =-1; getVarE(PTKEY_A, Ref(PoiGrpAddress),AdjGrpAddr,i);  //读取边上poi集的地址(成功）
            int poi_size = 0;
            //cout<<"poiList size= "<<poiList.size()<<endl;
            if(PoiGrpAddress !=-1){
                //读取poi索引位置
                int PoiIdxAddress =-1; getVarE(PIDX_A, Ref(PoiIdxAddress),AdjGrpAddr,i);

                poi_size = -3;
                int Ni=-1; int Nj=-1; float Ddist= 0;
                //读取该兴趣点信息，首先度边信息：Ni,Nj, Dist, poiGrpAddress--------------------------------
                getFixedF(NI_P,Ref(Ni),PoiGrpAddress);    //
                getFixedF(NJ_P,Ref(Nj),PoiGrpAddress);
                //cout <<" Ni="<<Ni<<",Nj="<<Nj<<endl;
                getFixedF(DIST_P,Ref(Ddist),PoiGrpAddress); //读取边上poi个数
                //cout<<"Ddist="<<Ddist<<endl;
                getFixedF(SIZE_P,Ref(poi_size),PoiGrpAddress); //读取边上poi个数
                //cout<<"poi_size="<<poi_size<<endl;

                for (size_t j = 0; j < poi_size; j++) { //------------取得边上poi集合的大小
                    //POI tmpPOI = poiList[j];   //------------取得第j个poi的地址
                    //tmpPOI = POIs[60];
                    //
                    //if (tmpPOI.category == 0)
                    //continue;
                    int pId = 0;    //------------取得第n个poi的id
                    int poi_indx_address = getPoiIdxAddr(PoiIdxAddress,j);  //索引表所在的地址
                    int poi_address;
                    getPOIAddr(&poi_address, poi_indx_address);  //获取poi数据所在地址
                    int poi_id = 0;
                    //getPOIID(&pId,poi_address);
                    POI p = getPOIDATA(poi_address);
                    /*if(p.id!=tmpPOI.id)
                        cout<<"读取p.id="<<p.id<<"真实id="<<tmpPOI.id<<endl;
                    if(p.keywords.size()!= tmpPOI.keywords.size()){
                        cout<<"keywords 读取情况与实际不符合"<<endl;
                        cout<<"p"<<p.id<<"keyword size="<<p.keywords.size()<<",tmpPOI.keywords.size="<<tmpPOI.keywords.size()<<endl;
                        //getchar();
                    }
                    if(p.check_ins.size()!=tmpPOI.check_ins.size()){
                        cout<<"check in读取情况与实际不符合"<<endl;
                        cout<<"p"<<p.id<<"check size="<<p.check_ins.size()<<",tmpPOI.check_ins.size="<<tmpPOI.check_ins.size()<<endl;
                        //printSeeds(p.check_ins);
                        //printSeeds(tmpPOI.check_ins);
                        //getchar();
                    }*/
                    pId = p.id;
                    float pDis = p.dis;   //------------取得第n个poi的距离Ni的距离
                    vector<int> pKey = p.keywords;//tmpPOI.keywords;
                    //这里待加入社交

                    float distance = -1;

                    // user and poi are on the same edge
                    if ((tmpNi == qNi && tmpNj == qNj) || (tmpNi == qNj && tmpNj == qNi)) {
                        if (qDis > pDis)
                            distance = qDis - pDis;
                        else
                            distance = pDis - qDis;

                    }
                        //user and poi are on the different edge
                    else {
                        if (tmpNi == Ni) { //同一个点 此时Ni的邻接点Nj是兴趣点所在边的Ni
                            distance = pointval[tmpNi] + pDis;   //edgeDist - pDis;//

                        }
                        else{ //边的另一个点
                            distance = pointval[tmpNi] + edgeDist - pDis;

                        }

                    }

                    double simD = distance;
                    //记录拓展到兴趣点的当前最短距离
                    if(poi_rearch[pId] == false){
                        poi_rearch[pId] = true;
                        poiDistance[pId] = distance;
                    }else if(distance<poiDistance[pId])
                        poiDistance[pId] = distance;

                    //float simT = tfIDF(qKey, pKey);
                    double simT = textRelevance(qKey, pKey);
                    if (simT==0) continue;
                    //必须文本相关

                    //double pow(double x,double y) //求x的y次方
                    double simS;
                    //double social_score = calcSocialScore(qId,pId);
                    double social_score = reCalcuSocialScore(qId,p.check_ins);

                    if(social_score > 0)
                        simS = 1.0 + pow(a,social_score);
                    else
                        simS = 1.0;

                    double score = simT * simS / (1+simD);
                    /*if(pId == 57){
                        cout<<"find o57! simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;
                    }*/

                    Result tmpRlt(pId, tmpNi, tmpNj, pDis, distance, score, pKey);
                    //cout<<"qId="<<qId<<"pId="<<pId<<endl;
                    if (resultFinal.size() < Qk) {
                        resultFinal.push(tmpRlt);
                        //
                    }
                    else {
                        Result _tmpNode = resultFinal.top();
                        double tmpScore = _tmpNode.score;

                        if (tmpScore < score) {
                            resultFinal.pop();
                            resultFinal.push(tmpRlt);
                            //cout<<"simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;
                        }
                        if(resultFinal.top().score > score_bound){//当前已有更优的k个结果，拓展可提前终止
                            if(resultFinal.top().score > Users[qId].topkScore_current)
                                Users[qId].topkScore_current = resultFinal.top().score;
                            TopkQueryResult  topkCurrentResult(resultFinal.top().score,resultFinal);
                            return topkCurrentResult;  //注意此时的score只是current score
                        }
                    }
                }

            }
            //cout<<poi_size<<","<<poiList.size()<<endl;
            //for (size_t j = 0; j < poi_size; j++) {

            //首次拓展当前拓展顶点
            if(rearch[tmpNj]==false){
                pointval[tmpNj] = maxDisValue;
            }
            if (e.dist + tmp.dist < pointval[tmpNj]) {
                //cout<<e.dist<<endl;
                pointval[tmpNj] = e.dist + tmp.dist;
                nodeDistType tmp2(tmpNj, edgeDist + tmp.dist);
                //cout << tmp2.num << ' ' << tmp2.dist << endl;
                tmpPoint.push(tmp2);
                rearch[tmpNj]=true;
            }
        }
        //路网边拓展提前终止条件
        if (resultFinal.size() == Qk) {
            double tmpScore = resultFinal.top().score;
            double currentDistance = pointval[tmpNi];

            double bound = optimal_text_score*(a+1.0)/tmpScore-1;               //提前拓展终止条件出错
            if (currentDistance > bound){
                //cout<<"当前拓展终止"<<endl;
                break;
            }

        }
    }
    Result topkNode = resultFinal.top();
    double topkScore = topkNode.score;
    Users[qId].topkScore_current = topkScore;
    Users[qId].topkScore_Final = topkScore;
    TopkQueryResult  topkCurrentResult(resultFinal.top().score,resultFinal);
    return topkCurrentResult;  //注意此时的score只是current score
}



/*************************2**************************/
//采用Dijastra算法验证某兴趣点是否为用户的top-k结果,返回当前求到的topk 评分

//---基于内存
double topkSDijkstra_verify_memory(query Query, int Qk,float a, float alpha,double score_bound) {

    bool vis[VertexNum];
    bool rearch[VertexNum];
    bool poi_rearch[poi_num];
    float pointval[VertexNum];

    int nodeCnt = 0;
    float Max_D = 25000;

    memset(vis, false, sizeof(vis));
    memset(rearch, false, sizeof(rearch));
    float maxDisValue = 9999999.9;
    //memset(pointval, maxDisValue, sizeof(float)*sizeof(pointval));
    //for (int i = 0; i < VertexNum; ++i)
    //pointval[i] = maxNum;

    priority_queue<nodeDistType> tmpPoint;
    priority_queue<Result> resultFinal;

    int qId = Query.id;
    int qNi = Query.Ni;
    int qNj = Query.Nj;
    float qDist = Query.dist;
    float qDis = Query.dis;
    vector<int> qKey = Query.keywords;

    // expand the edge from Ni, Nj
    pointval[qNi] = qDis;
    pointval[qNj] = qDist - qDis;

    nodeDistType tmp1(qNi, qDis);
    tmpPoint.push(tmp1);
    rearch[qNi] = true;

    nodeDistType tmp2(qNj, qDist - qDis);
    tmp2.dist = qDist - qDis;
    tmpPoint.push(tmp2);
    rearch[qNj] = true;


    double  optimal_text_score = textRelevance(qKey, qKey);

    while (!tmpPoint.empty()) {
        nodeDistType tmp = tmpPoint.top();
        tmpPoint.pop();
        int tmpNi = tmp.ID;
        if (vis[tmpNi] == true)
            continue;
        nodeCnt++;
        vis[tmpNi] = true;

        vector<int> tmpAdj = adjList[tmpNi];
        for (size_t i = 0; i < tmpAdj.size(); i++) {  // expand each edge
            int tmpNj = tmpAdj[i];
            if (vis[tmpNj])
                continue;

            edge e = EdgeMap[getKey(tmpNi, tmpNj)];
            int Ni = e.Ni;
            float edgeDist = e.dist;
            vector<POI> poiList = e.pts;
            for (size_t j = 0; j < poiList.size(); j++) { // check each poi on the edge
                POI tmpPOI = poiList[j];

                //tmpPOI = POIs[60];
                if (tmpPOI.category == 0)
                    continue;

                int pId = tmpPOI.id;

                float pDis = tmpPOI.dis;
                vector<int> pKey = tmpPOI.keywords;
                //这里待加入社交

                float distance = -1;

                // user and poi are on the same edge
                if ((tmpNi == qNi && tmpNj == qNj) || (tmpNi == qNj && tmpNj == qNi)) {
                    if (qDis > pDis)
                        distance = qDis - pDis;
                    else
                        distance = pDis - qDis;
                }
                    //user and poi are on the different edge
                else {
                    if (tmpNi == e.Ni)  //同一个点
                        distance = pointval[tmpNi] + pDis;
                    else //边的另一个点
                        distance = pointval[tmpNi] + edgeDist - pDis;

                }

                double simD = distance;

                //float simT = tfIDF(qKey, pKey);
                double simT = 0;
                simT = textRelevance(qKey, pKey);

                //必须文本相关
                if (simT==0) continue;

                //double pow(double x,double y) //求x的y次方
                double simS;
                //double social_score = getSocialScore(qId,pId);
                double social_score = calcSocialScore(qId,pId);

                if(social_score > 0.0)
                    simS = 1.0 + pow(a,social_score);
                else
                    simS = 1.0;

                double score = simT * simS / (1+simD);

                Result tmpRlt(pId, tmpNi, tmpNj, pDis, distance, score, pKey);
                //cout<<"qId="<<qId<<"pId="<<pId<<endl;
                if (resultFinal.size() < Qk) {
                    resultFinal.push(tmpRlt);
                    //cout<<"simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;
                }
                else {
                    Result _tmpNode = resultFinal.top();
                    double tmpScore = _tmpNode.score;

                    if (tmpScore < score) {
                        resultFinal.pop();
                        resultFinal.push(tmpRlt);
                        //cout<<"simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;
                    }
                    if(resultFinal.top().score > score_bound){//当前已有更优的k个结果，拓展可提前终止
                        if(resultFinal.top().score > Users[qId].topkScore_current)
                            Users[qId].topkScore_current = resultFinal.top().score;
                        return resultFinal.top().score;  //注意此时的score只是current score
                    }
                }
            }
            //首次拓展当前拓展顶点
            if(rearch[tmpNj]==false){
                pointval[tmpNj] = maxDisValue;
            }
            if (e.dist + tmp.dist < pointval[tmpNj]) {
                pointval[tmpNj] = e.dist + tmp.dist;
                nodeDistType tmp2(tmpNj, edgeDist + tmp.dist);
                //cout << tmp2.num << ' ' << tmp2.dist << endl;
                tmpPoint.push(tmp2);
                rearch[tmpNj]=true;
            }
        }
        //路网边拓展提前终止条件
        if (resultFinal.size() == Qk) {
            double tmpScore = resultFinal.top().score;
            float currentDistance = pointval[tmpNi];
            double optimal_social_textual_score = alpha*(a+1.0)+(1-alpha)*optimal_text_score;
            double bound = optimal_social_textual_score /tmpScore - 1;               //提前拓展终止条件出错
            if (currentDistance > bound){
                //cout<<"当前拓展终止"<<endl;
                break;
            }

        }
    }
    Result topkNode = resultFinal.top();
    double topkScore = topkNode.score;
    Users[qId].topkScore_current = topkScore;
    Users[qId].topkScore_Final = topkScore;
    return topkScore;
}

TopkQueryCurrentResult topkSDijkstra_verify_usr_memory(query Query, int Qk,float a, float alpha,float gsk_score, float usr_lcl_rt) {

    int nodeCnt = 0;
    float Max_D = 25000;
    bool vis[VertexNum];
    bool rearch[VertexNum];
    bool poi_rearch[poi_num];
    float pointval[VertexNum];

    memset(vis, false, sizeof(vis));
    memset(rearch, false, sizeof(rearch));
    float maxDisValue = 9999999.9;
    //memset(pointval, maxDisValue, sizeof(float)*sizeof(pointval));
    //for (int i = 0; i < VertexNum; ++i)
    //pointval[i] = maxNum;

    priority_queue<nodeDistType> tmpPoint;
    priority_queue<Result> resultFinal;


    int qId = Query.id;
    int qNi = Query.Ni;
    int qNj = Query.Nj;
    float qDist = Query.dist;
    float qDis = Query.dis;
    vector<int> qKey = Query.keywords;
    double  optimal_text_score = textRelevance(qKey, qKey);
    double optimal_social_textual_score = alpha*(a+1.0)+(1-alpha)*optimal_text_score;
    double score_rt = usr_lcl_rt;

    // expand the edge from Ni, Nj
    pointval[qNi] = qDis;
    pointval[qNj] = qDist - qDis;

    nodeDistType tmp1(qNi, qDis);
    tmpPoint.push(tmp1);

    nodeDistType tmp2(qNj, qDist - qDis);
    tmp2.dist = qDist - qDis;
    tmpPoint.push(tmp2);


    while (!tmpPoint.empty()) {
        nodeDistType tmp = tmpPoint.top();
        tmpPoint.pop();
        int tmpNi = tmp.ID;
        if (vis[tmpNi] == true)
            continue;
        nodeCnt++;
        vis[tmpNi] = true;

        vector<int> tmpAdj = adjList[tmpNi];    //取得拓展点的邻接表（地址）
        for (size_t i = 0; i < tmpAdj.size(); i++) {  // expand each edge
            int tmpNj = tmpAdj[i];
            if (vis[tmpNj])
                continue;

            edge e = EdgeMap[getKey(tmpNi, tmpNj)];   //取得邻接表中的edge（地址）
            int Ni = e.Ni;
            float edgeDist = e.dist;
            vector<POI> poiList = e.pts;              //取得edge上的poi集合 （地址）
            for (size_t j = 0; j < poiList.size(); j++) { // check each poi on the edge
                POI tmpPOI = poiList[j];           //获取poi （地址）

                //tmpPOI = POIs[60];

                if (tmpPOI.category == 0)
                    continue;

                int pId = tmpPOI.id;

                float pDis = tmpPOI.dis;
                vector<int> pKey = tmpPOI.keywords;
                //这里待加入社交

                float distance = -1;

                // user and poi are on the same edge
                if ((tmpNi == qNi && tmpNj == qNj) || (tmpNi == qNj && tmpNj == qNi)) {
                    if (qDis > pDis)
                        distance = qDis - pDis;
                    else
                        distance = pDis - qDis;
                }
                    //user and poi are on the different edge
                else {
                    if (tmpNi == e.Ni)  //同一个点
                        distance = pointval[tmpNi] + pDis;
                    else //边的另一个点
                        distance = pointval[tmpNi] + edgeDist - pDis;

                }

                double simD = distance;

                //float simT = tfIDF(qKey, pKey);
                double simT = textRelevance(qKey, pKey);

                //必须文本相关
                if (simT==0) continue;

                //double pow(double x,double y) //求x的y次方
                double simS;
                //double social_score = getSocialScore(qId,pId);
                double social_score = calcSocialScore(qId,pId);
                /*if(social_score == 0)
                    simS = 1.0;
                else
                    simS = 1.0 + pow(a,social_score);*/
                if(social_score > 0.0)
                    simS = 1.0 + pow(a,social_score);
                else
                    simS = 1.0;


                double social_textual_score = alpha*simS + (1-alpha)*simT;

                double score = social_textual_score / (1+simD);

                Result tmpRlt(pId, tmpNi, tmpNj, pDis, distance, score, pKey);
                //cout<<"qId="<<qId<<"pId="<<pId<<endl;
                if (resultFinal.size() < Qk) {
                    resultFinal.push(tmpRlt);
                    //cout<<"simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;
                }
                else {
                    Result _tmpNode = resultFinal.top();
                    double tmpScore = _tmpNode.score;

                    if (tmpScore < score) {
                        resultFinal.pop();
                        resultFinal.push(tmpRlt);
                        //cout<<"simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;
                    }
                    if(resultFinal.top().score>score_rt)  score_rt = resultFinal.top().score;


                    if(resultFinal.top().score > gsk_score){//当前已有更优的k个结果，拓展可提前终止
                        if(resultFinal.top().score > Users[qId].topkScore_current)
                            Users[qId].topkScore_current = resultFinal.top().score;
                        TopkQueryCurrentResult topkCurrentResult(resultFinal.top().score,resultFinal);
                        return topkCurrentResult;  //注意此时的score只是current score
                    }
                }
            }
            //首次拓展当前拓展顶点
            if(rearch[tmpNj]==false){
                pointval[tmpNj] = maxDisValue;
            }
            if (e.dist + tmp.dist < pointval[tmpNj]) {
                pointval[tmpNj] = e.dist + tmp.dist;
                nodeDistType tmp2(tmpNj, edgeDist + tmp.dist);
                //cout << tmp2.num << ' ' << tmp2.dist << endl;
                if(optimal_social_textual_score/(1.0+pointval[tmpNj])>score_rt){  //判断其是否还有拓展价值
                    tmpPoint.push(tmp2);
                    rearch[tmpNj]=true;
                }

            }
        }
        //路网边拓展提前终止条件
        if (resultFinal.size() == Qk) {
            double tmpScore = resultFinal.top().score;
            float currentDistance = pointval[tmpNi];
            double bound = optimal_social_textual_score /tmpScore - 1;
            if (currentDistance > bound){
                //cout<<"当前拓展终止"<<endl;
                break;
            }

        }
    }
    Result topkNode = resultFinal.top();
    double topkScore = topkNode.score;
    Users[qId].topkScore_current = topkScore;
    Users[qId].topkScore_Final = topkScore;
    TopkQueryCurrentResult topkCurrentResult(resultFinal.top().score,resultFinal);
    return topkCurrentResult;  //注意此时的score只是current score
}

//---基于外存
double topkSDijkstra_verify_disk(query Query, int Qk,float a, float alpha,double score_bound) {
    int nodeCnt = 0;
    float Max_D = 25000;
    //double score_bound = 1000000000000000;
    bool vis[VertexNum];
    bool rearch[VertexNum];
    bool poi_rearch[poi_num];
    float pointval[VertexNum];

    memset(vis, false, sizeof(vis));
    memset(rearch, false, sizeof(rearch));
    memset(poi_rearch, false, sizeof(poi_rearch));
    float maxDisValue = 9999999.9;

    priority_queue<nodeDistType> tmpPoint;
    priority_queue<Result> resultFinal;

    int qId = Query.id;
    int qNi = Query.Ni;
    int qNj = Query.Nj;
    float qDist = Query.dist;
    float qDis = Query.dis;
    vector<int> qKey = Query.keywords;
    double  optimal_text_score = textRelevance(qKey, qKey);

    // expand the edge from Ni, Nj
    pointval[qNi] = qDis;
    pointval[qNj] = qDist - qDis;

    nodeDistType tmp1(qNi, qDis);
    tmpPoint.push(tmp1);

    nodeDistType tmp2(qNj, qDist - qDis);
    tmp2.dist = qDist - qDis;
    tmpPoint.push(tmp2);

    while (!tmpPoint.empty()) {
        nodeDistType tmp = tmpPoint.top();
        tmpPoint.pop();
        int tmpNi = tmp.ID;
        if (vis[tmpNi] == true)
            continue;
        nodeCnt++;
        vis[tmpNi] = true;

        int AdjGrpAddr,AdjListSize;
        vector<int> tmpAdj = adjList[tmpNi];
        /*1. 取得地址，   2. 根据地址从文件中取得数据*/
        //---------------取得顶点tmpNi的数据地址：AdjGrpAddr
        AdjGrpAddr=getAdjListGrpAddr(tmpNi);
        //cout<<"AdjGrpAddr="<<AdjGrpAddr<<endl;
        bool testanswer=true;
        //---------------取得顶点tmpNi的邻接点个数：AdjListSize
        getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
        //cout<<"开始N"<<tmpNi<<"的拓展，邻接点size= "<<AdjListSize<<endl;
        //for (size_t i = 0; i < tmpAdj.size(); i++) {  // expand each edge
        for (size_t i = 0; i < AdjListSize; i++) {  // expand each edge
            //int tmpNj = tmpAdj[i];
            //------------取得邻接顶点tmpNj的id：tmpNj
            int tmpNj;
            getVarE(ADJNODE_A,Ref(tmpNj),AdjGrpAddr,i);  //读到的Nj不对
            //cout<<"读到Nj"<<tmpNj<<endl;
            if (vis[tmpNj])
                continue;

            edge e = EdgeMap[getKey(tmpNi, tmpNj)];
            //cout<<"e.dist="<<e.dist<<endl;
            int Ni = e.Ni;
            //float edgeDist = e.dist;
            float edgeDist; getVarE(DIST_A,Ref(edgeDist),AdjGrpAddr,i);  //读取边的长度(成功)

            vector<POI> poiList = e.pts;

            int PoiGrpAddress =-1; getVarE(PTKEY_A, Ref(PoiGrpAddress),AdjGrpAddr,i);  //读取边上poi集的地址(成功）
            int poi_size = 0;
            //cout<<"poiList size= "<<poiList.size()<<endl;
            if(PoiGrpAddress !=-1){
                //读取poi索引位置
                int PoiIdxAddress =-1; getVarE(PIDX_A, Ref(PoiIdxAddress),AdjGrpAddr,i);

                poi_size = -3;
                int Ni=-1; int Nj=-1; float Ddist= 0;
                //读取该兴趣点信息，首先度边信息：Ni,Nj, Dist, poiGrpAddress--------------------------------
                getFixedF(NI_P,Ref(Ni),PoiGrpAddress);    //
                getFixedF(NJ_P,Ref(Nj),PoiGrpAddress);
                //cout <<" Ni="<<Ni<<",Nj="<<Nj<<endl;
                getFixedF(DIST_P,Ref(Ddist),PoiGrpAddress); //读取边上poi个数
                //cout<<"Ddist="<<Ddist<<endl;
                getFixedF(SIZE_P,Ref(poi_size),PoiGrpAddress); //读取边上poi个数
                //cout<<"poi_size="<<poi_size<<endl;
                for (size_t j = 0; j < poi_size; j++) { //------------取得边上poi集合的大小
                    int pId = 0;    //------------取得第n个poi的id
                    int poi_indx_address = getPoiIdxAddr(PoiIdxAddress,j);  //索引表所在的地址
                    int poi_address;
                    getPOIAddr(&poi_address, poi_indx_address);  //获取poi数据所在地址
                    int poi_id = 0;
                    //getPOIID(&pId,poi_address);
                    POI p = getPOIDATA(poi_address);
                    pId = p.id;
                    float pDis = p.dis;   //------------取得第n个poi的距离Ni的距离
                    vector<int> pKey = p.keywords;//tmpPOI.keywords;
                    //这里待加入社交

                    float distance = -1;
                    // user and poi are on the same edge
                    if ((tmpNi == qNi && tmpNj == qNj) || (tmpNi == qNj && tmpNj == qNi)) {
                        if (qDis > pDis)
                            distance = qDis - pDis;
                        else
                            distance = pDis - qDis;

                    }
                        //user and poi are on the different edge
                    else {
                        if (tmpNi == Ni) { //同一个点 此时Ni的邻接点Nj是兴趣点所在边的Ni
                            distance = pointval[tmpNi] + pDis;   //edgeDist - pDis;//

                        }
                        else{ //边的另一个点
                            distance = pointval[tmpNi] + edgeDist - pDis;

                        }
                    }

                    double simD = distance;
                    //记录拓展到兴趣点的当前最短距离
                    if(poi_rearch[pId] == false){
                        poi_rearch[pId] = true;
                        poiDistance[pId] = distance;
                    }else if(distance<poiDistance[pId])
                        poiDistance[pId] = distance;

                    //float simT = tfIDF(qKey, pKey);
                    double simT = textRelevance(qKey, pKey);
                    if (simT==0) continue;  //必须文本相关

                    //double pow(double x,double y) //求x的y次方
                    double simS;
                    //double social_score = calcSocialScore(qId,pId);
                    double social_score = reCalcuSocialScore(qId,p.check_ins);

                    if(social_score > 0)
                        simS = 1.0 + pow(a,social_score);
                    else
                        simS = 1.0;

                    double score = simT * simS / (1+simD);

                    Result tmpRlt(pId, tmpNi, tmpNj, pDis, distance, score, pKey);
                    if (resultFinal.size() < Qk) {
                        resultFinal.push(tmpRlt);
                    }
                    else {
                        Result _tmpNode = resultFinal.top();
                        double tmpScore = _tmpNode.score;

                        if (tmpScore < score) {
                            resultFinal.pop();
                            resultFinal.push(tmpRlt);
                            //cout<<"simD="<<simD<<",simT="<<simT<<",score="<<score<<endl;
                        }
                        if(resultFinal.top().score > score_bound){//当前已有更优的k个结果，拓展可提前终止
                            if(resultFinal.top().score > Users[qId].topkScore_current)
                                Users[qId].topkScore_current = resultFinal.top().score;
                            //TopkQueryResult  topkCurrentResult(resultFinal.top().score,resultFinal);
                            return resultFinal.top().score;  //注意此时的score只是current score
                        }
                    }
                }

            }

            //首次拓展当前拓展顶点
            if(rearch[tmpNj]==false){
                pointval[tmpNj] = maxDisValue;
            }
            if (e.dist + tmp.dist < pointval[tmpNj]) {
                //cout<<e.dist<<endl;
                pointval[tmpNj] = e.dist + tmp.dist;
                nodeDistType tmp2(tmpNj, edgeDist + tmp.dist);
                //cout << tmp2.num << ' ' << tmp2.dist << endl;
                tmpPoint.push(tmp2);
                rearch[tmpNj]=true;
            }
        }
        //路网边拓展提前终止条件
        if (resultFinal.size() == Qk) {
            double tmpScore = resultFinal.top().score;
            double currentDistance = pointval[tmpNi];

            double bound = optimal_text_score*(a+1.0)/tmpScore-1;               //提前拓展终止条件出错
            if (currentDistance > bound){
                //cout<<"当前拓展终止"<<endl;
                break;
            }

        }
    }
    Result topkNode = resultFinal.top();
    double topkScore = topkNode.score;
    Users[qId].topkScore_current = topkScore;
    Users[qId].topkScore_Final = topkScore;
    return topkScore;

}



//原来的

float SPD_by_topkSDijkstra_disk(int s, int e) {


    int nodeCnt = 0;
    float Max_D = 25000;
    bool vis[VertexNum];
    bool rearch[VertexNum];
    float pointval[VertexNum];

    map<int,bool> EdgeVisFlag;    // <edge_id, pre_midresults>

    memset(vis, false, sizeof(vis));
    memset(rearch, false, sizeof(rearch));

    float maxDisValue = 9999999.9;

    priority_queue<nodeDistType> tmpPoint;
    priority_queue<Result> resultFinal;


    // expand the edge from Ni,
    pointval[s] = 0.0;
    nodeDistType tmp1(s, 0);
    tmpPoint.push(tmp1);
    rearch[e]=true;  //重要！



    while (!tmpPoint.empty()) {
        nodeDistType tmp = tmpPoint.top();
        tmpPoint.pop();
        int tmpNi = tmp.ID;

        if(tmpNi == e)
            return tmp.dist;


        nodeCnt++;
        vis[tmpNi] = true;

        int AdjGrpAddr,AdjListSize;
        vector<int> tmpAdj = adjList[tmpNi];
        /*1. 取得地址，   2. 根据地址从文件中取得数据*/
        //---------------取得顶点tmpNi的数据地址：AdjGrpAddr
        AdjGrpAddr=getAdjListGrpAddr(tmpNi);
        //cout<<"AdjGrpAddr="<<AdjGrpAddr<<endl;
        bool testanswer=true;
        //---------------取得顶点tmpNi的邻接点个数：AdjListSize
        getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
        //getFixedF_Self(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
        //cout<<"开始N"<<tmpNi<<"的拓展，邻接点size= "<<AdjListSize<<endl;
        for (size_t j = 0; j < AdjListSize; j++) {  // expand each edge
            //int tmpNj = tmpAdj[i];
            //------------取得邻接顶点tmpNj的id：tmpNj
            int tmpNj;
            getVarE(ADJNODE_A,Ref(tmpNj),AdjGrpAddr,j);  //读到的Nj不对
            //cout<<"for Ni"<<tmpNi<<"读到Nj"<<tmpNj<<endl;
            /* if(vis[tmpNj] ==true)  //
                 continue;*/
            edge e = EdgeMap[getKey(tmpNi, tmpNj)];
            int Ni = e.Ni;
            float edgeDist; getVarE(DIST_A,Ref(edgeDist),AdjGrpAddr,j);  //读取边的长度(成功)

            if(vis[tmpNj] ==true) continue; //之前已从其拓展过，则跳过不用拓展了
            //首次拓展当前拓展顶点
            if(rearch[tmpNj]==false){
                //pointval[tmpNj] = maxDisValue;
                pointval[tmpNj] = e.dist + tmp.dist;
                nodeDistType tmp2(tmpNj, edgeDist + tmp.dist);
                tmpPoint.push(tmp2);
                rearch[tmpNj]=true;
                //cout<<"tmpNj："<<tmpNj<<"入列1！"<<endl;
            }
            else if (e.dist + tmp.dist < pointval[tmpNj]) {
                //cout<<e.dist<<endl;
                pointval[tmpNj] = e.dist + tmp.dist;
                nodeDistType tmp2(tmpNj, edgeDist + tmp.dist);
                tmpPoint.push(tmp2);
                rearch[tmpNj]=true;
                //cout<<"tmpNj："<<tmpNj<<"入列2！"<<endl;
                //}

            }
        }

    }
    cout<<"v"<<s<<"到v"<<e<<",不可达!"<<endl;
    return 9999999999999999.9;  //注意此时的score只是current score

}

float SPD_by_topkSDijkstra_memory(int s, int e){

    int nodeCnt = 0;
    float Max_D = 25000;
    bool vis[VertexNum];
    bool rearch[VertexNum];
    float pointval[VertexNum];

    memset(vis, false, sizeof(vis));
    memset(rearch, false, sizeof(rearch));
    float maxDisValue = 9999999.9;
    //memset(pointval, maxDisValue, sizeof(float)*sizeof(pointval));
    //for (int i = 0; i < VertexNum; ++i)
    //pointval[i] = maxNum;

    priority_queue<nodeDistType> tmpPoint;
    priority_queue<Result> resultFinal;


    pointval[s] = 0;
    nodeDistType tmp1(s, 0);
    tmpPoint.push(tmp1);

    double p2p_dist = 999999999999999999999.9;
    while (!tmpPoint.empty()) {
        nodeDistType tmp = tmpPoint.top();
        tmpPoint.pop();
        int tmpNi = tmp.ID;
        if (vis[tmpNi] == true)
            continue;
        nodeCnt++;
        vis[tmpNi] = true;
        if(tmpNi==e) {//到达终点
            p2p_dist = tmp.dist;
            //return  tmp.dist;
        }


        vector<int> tmpAdj = adjList[tmpNi];    //取得拓展点的邻接表（地址）
        for (size_t i = 0; i < tmpAdj.size(); i++) {  // expand each edge
            int tmpNj = tmpAdj[i];
            if (vis[tmpNj])
                continue;

            edge e = EdgeMap[getKey(tmpNi, tmpNj)];   //取得邻接表中的edge（地址）
            int Ni = e.Ni;
            float edgeDist = e.dist;

            //首次拓展当前拓展顶点
            if(rearch[tmpNj]==false){
                pointval[tmpNj] = maxDisValue;
            }
            if (e.dist + tmp.dist < pointval[tmpNj]) {
                pointval[tmpNj] = e.dist + tmp.dist;
                nodeDistType tmp2(tmpNj, edgeDist + tmp.dist);
                //cout << tmp2.num << ' ' << tmp2.dist << endl;
                tmpPoint.push(tmp2);
                rearch[tmpNj]=true;

            }
        }
        //路网边拓展提前终止条件

    }
    //cout<<"不可达！"<<endl;
    return p2p_dist;  //注意此时的score只是current score
}




TopkQueryCurrentResult topkSDijkstra_verify_usr_disk(query Query, int Qk,float a, float alpha,float gsk_score, float score_bound) {


    int nodeCnt = 0;
    float Max_D = 25000;
    bool* vis = new bool[VertexNum];
    bool* rearch = new bool[VertexNum];
    float* pointval = new float[VertexNum];
    for(int i=0;i<VertexNum;i++){
        vis[i]=false;
        rearch[i]=false;
    }


    bool* poi_rearch = new bool[poi_num];
    bool* poi_evaluated = new bool [poi_num];
    for(int j=0;j<poi_num;j++){
        poi_rearch[j]=false;
        poi_evaluated[j]=false;
    }


    map<int,vector<Result>> resultCache;    // <edge_id, pre_midresults>
    map<int,bool> EdgeVisFlag;    // <edge_id, pre_midresults>

    float maxDisValue = 9999999.9;

    priority_queue<nodeDistType> tmpPoint;
    priority_queue<Result> resultFinal;
    double score_max =  score_bound;

    int qId = Query.id;
    int qNi = Query.Ni;
    int qNj = Query.Nj;
    float qDist = Query.dist;
    float qDis = Query.dis;
    vector<int> qKey = Query.keywords;
    set<int> qKeySet;
    for(int term:qKey)
        qKeySet.insert(term);
    double  optimal_text_score = textRelevance(qKey, qKey);
    //最优score
    double optimal_social_textual_score = alpha*(a+1.0)+(1-alpha)*optimal_text_score;

    //首先判断user所


    // expand the edge from Ni,
    pointval[qNi] = qDis;
    pointval[qNj] = qDist - qDis;

    nodeDistType tmp1(qNi, qDis);
    tmpPoint.push(tmp1);
    rearch[qNi]=true;  //重要！

    nodeDistType tmp2(qNj, qDist - qDis);
    tmp2.dist = qDist - qDis;
    tmpPoint.push(tmp2);
    rearch[qNj]=true;  //重要！

    while (!tmpPoint.empty()) {
        nodeDistType tmp = tmpPoint.top();
        tmpPoint.pop();
        int tmpNi = tmp.ID;

        //cout<<"当前队首顶点"<<tmpNi<<endl;
        /*if (vis[tmpNi] == true)
            continue;*/
        nodeCnt++;
        vis[tmpNi] = true;

        int AdjGrpAddr,AdjListSize;
        vector<int> tmpAdj = adjList[tmpNi];
        /*1. 取得地址，   2. 根据地址从文件中取得数据*/
        //---------------取得顶点tmpNi的数据地址：AdjGrpAddr
        AdjGrpAddr=getAdjListGrpAddr(tmpNi);
        //cout<<"AdjGrpAddr="<<AdjGrpAddr<<endl;
        bool testanswer=true;
        //---------------取得顶点tmpNi的邻接点个数：AdjListSize
        getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
        //getFixedF_Self(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
        //cout<<"开始N"<<tmpNi<<"的拓展，邻接点size= "<<AdjListSize<<endl;
        for (size_t j = 0; j < AdjListSize; j++) {  // expand each edge
            //int tmpNj = tmpAdj[i];
            //------------取得邻接顶点tmpNj的id：tmpNj
            int tmpNj;
            getVarE(ADJNODE_A,Ref(tmpNj),AdjGrpAddr,j);  //读到的Nj不对
            //cout<<"for Ni"<<tmpNi<<"读到Nj"<<tmpNj<<endl;
            /* if(vis[tmpNj] ==true)  //
                 continue;*/
            edge e = EdgeMap[getKey(tmpNi, tmpNj)];
            int Ni = e.Ni;
            float edgeDist; getVarE(DIST_A,Ref(edgeDist),AdjGrpAddr,j);  //读取边的长度(成功)
            if (vis[tmpNj] ==true){  //另外一侧的顶点也已拓展过， 将e(Ni, Nj)上相关的poi进行比较、进行最终更新，并加入
                for(Result rr: resultCache[getKey(tmpNi,tmpNj)]){
                    float dis2 = e.dist - rr.dis;
                    float distance2 = tmp.dist + dis2;
                    if(distance2<rr.dist){  //更新最终结果
                        double score2 = rr.score*(rr.dist/distance2);
                        Result tmpRlt(rr.id, tmpNi, tmpNj, dis2, distance2, score2, rr.key);
                        rr = tmpRlt;
                    }
                    if(poi_evaluated[rr.id]==true)
                        continue;
                    else
                        poi_evaluated[rr.id] = true;

                    if (resultFinal.size() < Qk) {
                        resultFinal.push(rr);
                    }
                    else {
                        Result _tmpNode = resultFinal.top();
                        double tmpScore = _tmpNode.score;

                        if (tmpScore < rr.score) {
                            resultFinal.pop();
                            resultFinal.push(rr);
                        }

                        if(resultFinal.top().score>score_max)  score_max = resultFinal.top().score;


                        if(resultFinal.top().score > gsk_score){//当前已有更优的k个结果，拓展可提前终止
                            if(resultFinal.top().score > Users[qId].topkScore_current)
                                Users[qId].topkScore_current = resultFinal.top().score;
                            TopkQueryCurrentResult  topkCurrentResult(resultFinal.top().score,resultFinal);

                            if(true){
                                delete []vis;
                                delete []rearch;
                                delete []pointval;
                                delete []poi_rearch ;
                                delete []poi_evaluated;
                            }

                            return  topkCurrentResult;  //注意此时的score只是current score
                        }
                    }

                }
            }
            else{  //另一侧顶点尚未拓展， 将Ni与相关poi关联起来

                //cout<<"edgeDist="<<edgeDist<<endl;

                //读取Okey psudo doc 的地址,并赋予okeyAddr


                set<int> OKeyset;
                OKeyset = getOKeyEdge(tmpNi,j);
                //cout<<"边上的关键词情况=";
                //printSetElements(OKeyset);
                //getchar();


                //检查边上的兴趣点的文本相关性，若不相关直接跳过该边
                set<int> iter = obtain_itersection_jins(qKeySet,OKeyset);  //ala
                if(iter.size() >0){ //

                    int poi_size = 0;
                    poi_size = getPOISizeOfEdge(tmpNi,j);
                    //cout<<"poi size = "<<poi_size<<endl;
                    if(poi_size >0){

                        for (size_t p_th = 0; p_th < poi_size; p_th++) { //------------取得边上poi集合数据

                            POI p = getPOIDataOnEdge(tmpNi,j,p_th);
                            int pId = p.id;

                            float pDis = p.dis;   //------------取得第n个poi的距离Ni的距离
                            vector<int> pKey = p.keywords;//tmpPOI.keywords;
                            //这里待加入社交
                            //float simT = tfIDF(qKey, pKey);
                            double simT = textRelevance(qKey, pKey);
                            //必须文本相关
                            if (simT==0) continue;

                            float distance = -1; float exp_dis;

                            // user and poi are on the same edge
                            if ((tmpNi == qNi && tmpNj == qNj) || (tmpNi == qNj && tmpNj == qNi)) {
                                if (qDis > pDis){
                                    distance = qDis - pDis;
                                    exp_dis = distance;
                                }

                                else{
                                    distance = pDis - qDis;
                                    exp_dis = distance;
                                }

                                double simD = distance;
                                //double pow(double x,double y) //求x的y次方
                                double simS;
                                //double social_score = calcSocialScore(qId,pId);
                                if(a==0){
                                    simS =1.0;
                                }
                                else{
                                    double social_score = reCalcuSocialScore(qId,p.check_ins);

                                    if(social_score > 0)
                                        simS = 1.0 + pow(a,social_score);
                                    else
                                        simS = 1.0;
                                }

                                double social_textual_score = alpha*simS + (1-alpha)*simT;

                                double score = social_textual_score / (1+simD);


                                Result tmpRlt(pId, tmpNi, tmpNj, pDis, distance, score, pKey);
                                resultCache[getKey(tmpNi,tmpNj)].push_back(tmpRlt);  //先缓存



                            }
                                //user and poi are on the different edge
                            else {
                                if (tmpNi == Ni) { //同一个点 此时Ni的邻接点Nj是兴趣点所在边的Ni
                                    distance = pointval[tmpNi] + pDis;   //edgeDist - pDis;//
                                    exp_dis = pDis;
                                }
                                else{ //边的另一个点
                                    distance = pointval[tmpNi] + (edgeDist - pDis);
                                    exp_dis = (edgeDist - pDis);
                                }
                                double simD = distance;
                                //double pow(double x,double y) //求x的y次方
                                double simS;
                                //double social_score = calcSocialScore(qId,pId);
                                if(a==0){
                                    simS =1.0;
                                }
                                else{
                                    double social_score = reCalcuSocialScore(qId,p.check_ins);

                                    if(social_score > 0)
                                        simS = 1.0 + pow(a,social_score);
                                    else
                                        simS = 1.0;
                                }

                                double social_textual_score = alpha*simS + (1-alpha)*simT;

                                double score = social_textual_score / (1+simD);


                                Result tmpRlt(pId, tmpNi, tmpNj, exp_dis, distance, score, pKey);
                                resultCache[getKey(tmpNi,tmpNj)].push_back(tmpRlt);  //先缓存


                            }


                            //记录拓展到兴趣点的当前最短距离
                            if(poi_rearch[pId] == false){
                                poi_rearch[pId] = true;
                                poiDistance[pId] = distance;
                            }else if(distance<poiDistance[pId])
                                poiDistance[pId] = distance;


                        }

                    }
                }
            }

            if(vis[tmpNj] ==true) continue; //之前已从其拓展过，则跳过不用拓展了
            //首次拓展当前拓展顶点
            if(rearch[tmpNj]==false){
                //pointval[tmpNj] = maxDisValue;
                pointval[tmpNj] = e.dist + tmp.dist;
                nodeDistType tmp2(tmpNj, edgeDist + tmp.dist);
                tmpPoint.push(tmp2);
                rearch[tmpNj]=true;
                //cout<<"tmpNj："<<tmpNj<<"入列1！"<<endl;
            }
            else if (e.dist + tmp.dist < pointval[tmpNj]) {
                //cout<<e.dist<<endl;
                pointval[tmpNj] = e.dist + tmp.dist;
                nodeDistType tmp2(tmpNj, edgeDist + tmp.dist);
                //cout << tmp2.num << ' ' << tmp2.dist << endl;

                //double score_later = optimal_social_textual_score/(1.0+pointval[tmpNj]);
                //if(score_later>score_max){  //当其仅当经Nj还能拓展到评分更优的poi时，则后续继续拓展Nj
                tmpPoint.push(tmp2);
                rearch[tmpNj]=true;
                //cout<<"tmpNj："<<tmpNj<<"入列2！"<<endl;
                //}

            }
        }
        //路网边拓展提前终止条件
        if (resultFinal.size() == Qk) {
            double tmpScore = resultFinal.top().score;
            double currentDistance = pointval[tmpNi];
            double bound = optimal_social_textual_score /tmpScore - 1;               //提前拓展终止条件出错
            if (currentDistance > bound){
                //cout<<"当前拓展终止"<<endl;
                break;
            }
        }
    }
    Result topkNode = resultFinal.top();
    double topkScore = topkNode.score;
    Users[qId].topkScore_current = topkScore;
    Users[qId].topkScore_Final = topkScore;
    TopkQueryCurrentResult  topkCurrentResult(resultFinal.top().score,resultFinal);
    if(true){
        delete []vis;
        delete []rearch;
        delete []pointval;
        delete []poi_rearch ;
        delete []poi_evaluated;
    }
    return topkCurrentResult;  //注意此时的score只是current score

}





TopkQueryCurrentResult topkBruteForce_disk(User user, int Qk,float a, float alpha,float gsk_score, float score_bound) {
    //TopkQueryCurrentResult  rr(0.08);
    //return  rr;

    //TopkQueryCurrentResult  topkCurrentResult;

    priority_queue<Result> resultFinal;

    for(int i=0;i<poi_num;i++){
        int poi_id = i;
        //poi_id = 348382;
        POI poi = getPOIFromO2UOrgLeafData(poi_id);
        //cout<<"p"<<poi_id<<"，Ni="<<poi.Ni<<", Nj="<<poi.Nj<<endl;
        //getchar();
        //getGSKScore_o2u(a, alpha, p, u);
        set<int> uKeyset; set<int> qKeyset;
        //for(int term : Users[usr_id].keywords) uKeyset.insert(term);
        //for(int term :qKeys)  qKeyset.insert(term);
        //double relevance = textRelevance(uKeyset,qKeyset);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
        vector<int> qKeys = poi.keywords;
        vector<int> check_usr = poiCheckInIDList[poi.id];
        double relevance = textRelevance(user.keywords,qKeys);
        if(relevance == 0){
            continue;
        }
        double influence;
        double social_Score = 0;
        if(a==0||alpha==0.0){
            influence = 1;
        }else{
            if(check_usr.size()>0)
                social_Score =  reCalcuSocialScore(user.id,check_usr);

            if(social_Score == 0)
                influence =1;
            else
                influence =1 + pow(a,social_Score);
        }


        double distance = usrToPOIDistance(user, poi);   //这里要注意区分 getDistance(vertex, vertex)
        double social_textual_score = alpha*influence + (1.0-alpha)*relevance;
        double gsk_score = social_textual_score / (1.0+distance);



        Result tmpRlt(poi.id, poi.Ni, poi.Nj, poi.dis, distance, gsk_score, poi.keywords);


        if (resultFinal.size() < Qk) {
            resultFinal.push(tmpRlt);
        }
        else {
            Result _tmpNode = resultFinal.top();
            double tmpScore = _tmpNode.score;

            if (tmpScore < tmpRlt.score) {
                resultFinal.pop();
                resultFinal.push(tmpRlt);
            }

        }


    }



    TopkQueryCurrentResult  topkCurrentResult(resultFinal.top().score,resultFinal);
    return topkCurrentResult;  //注意此时的score只是current score

}


TopkQueryCurrentResult topkBruteForce_disk_FromVertex(int vertex, int term, int Qk) {
    //TopkQueryCurrentResult  rr(0.08);
    //return  rr;

    vector<int> query_keywords; query_keywords.push_back(term);

    priority_queue<Result> resultFinal;

    for(int i=0;i<poi_num;i++){
        int poi_id = i;
        //poi_id = 348382;
        POI poi = getPOIFromO2UOrgLeafData(poi_id);
        //cout<<"p"<<poi_id<<"，Ni="<<poi.Ni<<", Nj="<<poi.Nj<<endl;
        //getchar();
        //getGSKScore_o2u(a, alpha, p, u);
        set<int> uKeyset; set<int> qKeyset;
        //for(int term : Users[usr_id].keywords) uKeyset.insert(term);
        //for(int term :qKeys)  qKeyset.insert(term);
        //double relevance = textRelevance(uKeyset,qKeyset);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
        vector<int> qKeys = poi.keywords;
        vector<int> check_usr = poiCheckInIDList[poi.id];
        double relevance = textRelevance(query_keywords,qKeys);
        if(relevance == 0){
            continue;
        }
        double influence =1;


        double dist1 = SPSP(vertex, poi.Ni)+poi.dis;
        double dist2 = SPSP(vertex, poi.Nj)+(poi.dist-poi.dis);

        double social_textual_score = relevance;
        double gsk_score = social_textual_score / (1+min(dist2,dist1));



        Result tmpRlt(poi.id, poi.Ni, poi.Nj, poi.dis, min(dist2,dist1), gsk_score, poi.keywords);


        if (resultFinal.size() < Qk) {
            resultFinal.push(tmpRlt);
        }
        else {
            Result _tmpNode = resultFinal.top();
            double tmpScore = _tmpNode.score;

            if (tmpScore < tmpRlt.score) {
                resultFinal.pop();
                resultFinal.push(tmpRlt);
            }

        }


    }

    TopkQueryCurrentResult  topkCurrentResult(resultFinal.top().score,resultFinal);
    return topkCurrentResult;  //注意此时的score只是current score

}




//generating network voronoi diagram


void NVD_generation_Keyword_aware(int target_term) {

    float term_weight = getTermIDFWeight(target_term);
    cout<<"----------------------term"<<target_term<<" idf_score="<<term_weight<<"----------------------"<<endl;
    vector<int> poi_ids = getTermOjectInvList(target_term);
    cout<<"object posting list size="<<poi_ids.size()<<"具体为："<<endl;
    printElements(poi_ids);

    int nodeCnt = 0;
    float Max_D = 25000;
    bool vis[VertexNum];
    bool rearch[VertexNum];
    bool poi_rearch[poi_num];
    float pointval[VertexNum];
    map<int,vector<Result>> resultCache;    // <edge_id, pre_midresults>
    map<int,bool> EdgeVisFlag;    // <edge_id, pre_midresults>

    memset(vis, false, sizeof(vis));
    memset(rearch, false, sizeof(rearch));
    memset(poi_rearch, false, sizeof(poi_rearch));
    float maxDisValue = 9999999.9;

    priority_queue<nodeDistType_NVD> Q;

    //hashmap: <vertex_id, poi_id>
    int vertex_point_pair[VertexNum];  //idx : vertex_id, value : poi_id;
    for(int i=0;i<VertexNum;i++)
        vertex_point_pair[i] = -100;
    //map<int, int> vertex_point_pairMap;
    //hashmap: <poi_id,  vertex_set> : road network partition by generating points (pois in the posting list for a specific keywords)
    map<int, set<int>> NVD_vertex_partition;
    //hashmap: <poi_id, adjacent_poi list>  : adjacent graph of NVD
    map<int, set<int>> NVD_adjacent_graph;


    //initialize for each generating points
    for(int p: poi_ids){
        POI poi = getPOIFromO2UOrgLeafData(p);
        //for poi.Ni
        int p_Ni = poi.Ni;
        float p_dis = poi.dis;
        pointval[p_Ni] = p_dis; //vertex_point_pairMap[p_Ni] = poi.id;
        nodeDistType_NVD element1(p_Ni, poi.id,poi.dis);
        rearch[p_Ni]=true;  //重要！
        Q.push(element1);

        //for poi.Nj
        int p_Nj = poi.Nj;
        pointval[p_Nj] = poi.dist-poi.dis; //vertex_point_pairMap[p_Ni] = poi.id;
        nodeDistType_NVD element2(p_Nj, poi.id,poi.dist-poi.dis);
        rearch[p_Nj]=true;  //重要！
        Q.push(element2);
    }



    while (!Q.empty()) {
        nodeDistType_NVD entry = Q.top();Q.pop();
        int cur_vertex = entry.vertex_id;
        int cur_poi = entry.poi_id;
        vertex_point_pair[cur_vertex] = cur_poi;  //当前vertex就处在当前poi的 Voronoi Cell中（因为距离最小）
        vis[cur_vertex] = true;
        //int AdjGrpAddr, AdjListSize;
        //获取当前顶点的邻接表信息
        vector<int> current_adj = adjList[cur_vertex];
        for(int i=0;i<current_adj.size();i++){
            int next_vertex =  current_adj[i];
            edge  e = EdgeMap[getKey(cur_vertex, next_vertex)];
            if(next_vertex==true) continue;
            if(rearch[next_vertex]==false){
                pointval[next_vertex] = entry.dist + e.dist;
                rearch[next_vertex] = true; //重要
                nodeDistType_NVD ele_next(next_vertex,entry.poi_id,entry.dist + e.dist);
                Q.push(ele_next);
            }
            else if(entry.dist + e.dist < pointval[next_vertex]){
                pointval[next_vertex] = entry.dist + e.dist;
                rearch[next_vertex] = true; //重要
                if(entry.poi_id==19324){
                    cout<<"find entry.poi_id==19324!"<<endl;
                }
                nodeDistType_NVD entry_next(next_vertex,entry.poi_id,entry.dist + e.dist);
                Q.push(entry_next);
            }
        }
    }
    //结束路网边遍历
    //打印出每个顶点对应的 Voronoi Cell 中, 将v加入到poi的关联顶点列表中
    map<int, set<int>> NVD_AdjGraph;
    map<int, set<int>> NVD_BorderMap;
    for(int i=0;i<VertexNum;i++){
        int v = i;
        int poi = vertex_point_pair[v];
        printf("v%d 对应的最近邻poi为o %d\n",v,poi);
        NVD_vertex_partition[poi].insert(v);
        vector<int> v_adj = adjList[v];
        for(int u : v_adj){
            int poi2 = vertex_point_pair[u];
            if(poi == poi2){
                continue;
            }
            else{  //建立generating point间的邻接关系
                NVD_AdjGraph[poi].insert(poi2);
                NVD_AdjGraph[poi2].insert(poi);
                //v为poi所在NV Cell的边界点
                NVD_BorderMap[poi].insert(v);
                //u为poi2所在NV Cell的边界点
                NVD_BorderMap[poi2].insert(u);
            }
        }
    }
    //打印NVD AdjGraph结果
    if(true){
        map<int, set<int>> leaf2POIMap;
        for(int poi:poi_ids){
            if(NVD_AdjGraph[poi].size()>0){
                printf("兴趣点o%d 的NVD结果是：border_size=%d, \n adjacent_poi size = %d, 毗邻的poi分别为：",poi, NVD_BorderMap[poi].size(), NVD_AdjGraph[poi].size());
                for(int id: NVD_AdjGraph[poi]){
                    cout<<"o"<<id<<", ";
                }
                cout<<endl;
                printf("partition_size(v) = %d,分别为：\n", NVD_vertex_partition[poi].size());
                for(int v: NVD_vertex_partition[poi]){
                    //cout<<"包括 v"<<v<<",";
                    //cout<<"验证 v"<<v<<"..."; int term = 2;
                    int leafNode = Nodes[v].gtreepath.back();
                    leaf2POIMap[leafNode].insert(poi);
                    //cout<<"位于Gtree 叶节点"<<leafNode<<";"<<endl;
                    //TopkQueryCurrentResult topk_r5 = topkBruteForce_disk_FromVertex(v,term,5);
                    //cout<<"具体为："<<endl;
                    //topk_r5.printResults(5);
                }
                //cout<<endl;

                getchar();
            }
            //cout<<endl;
        }
        //打印 leaf2POIMap中结果
        map<int,set<int>>:: iterator iter = leaf2POIMap.begin();
        int count = 0;
        while(iter != leaf2POIMap.end()){
            int leaf_id = iter->first; set<int> poi_Set = iter->second;
            cout<<"叶节点"<<leaf_id<<"对应poi size="<<poi_Set.size()<<"，具体为：";
            printSetElements(poi_Set);
            count ++;
            iter++;
        }
        cout<<"总共"<<count<<"个叶节点"<<endl;
    }

    //将NVD 内容输出到disk文件中
    //


}





//从文本文件中加载hash index 内容
/*
void Load_V2P_NVDG_AddressIdx(){
    clock_t startTime, endTime;
    startTime = clock();
    ifstream v2p_inputFile, nvd_inputFile;
    char pre_fileName[255];   char tmpFileName[255];
    sprintf(pre_fileName,"../../../data/%s/NVD/",dataset_name);
    sprintf(tmpFileName,"%sv2p.idx",pre_fileName);
    v2p_inputFile.open(tmpFileName);
    cout<<"v2p idx(数据地址索引)文件 open!"<<endl;

    sprintf(tmpFileName,"%sNVD.idx",pre_fileName);
    nvd_inputFile.open(tmpFileName);
    cout<<"nvd idx(数据地址索引)文件 open!"<<endl;


    string str;
    while(getline(v2p_inputFile, str)){  //先读到buffer上
        istringstream tt(str);
        int key = -1; int address_value = -1;
        tt>>key>>address_value;
        int v_id = key%VertexNum;  int term_id = (int)((key-v_id)/VertexNum);

        addressHashMap_v2p[key] = address_value;
        //cout<<"read pair from v2p.idx: <key="<<key<<",address="<<address_value<<">"<<endl;


    }
    v2p_inputFile.close();
    endTime = clock();
    cout<<"导入v2p idx(数据地址索引)文件 完毕!用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s" << endl;

    //读取nvd
    startTime = clock();
    while(getline(nvd_inputFile, str)){  //先读到buffer上
        istringstream tt(str);
        int key = -1; int address_value = -1;
        tt>>key>>address_value;
        int poi_id = key%poi_num;  int term_id = (int)((key-poi_id)/poi_num);

        addressHashMap_nvd[key] = address_value;


    }
    nvd_inputFile.close();
    endTime = clock();
    cout<<"导入nvd idx(数据地址索引)文件 完毕! 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s" << endl;


}
*/

//从序列化文件中加载hash index 内容
void Load_V2P_NVDG_AddressIdx_fast(){
    clock_t startTime, endTime;
    startTime = clock();
    string indexPath = getNVDHashIndexInputPath();
    V2PNVDHashIndex index = serialization::getIndexFromBinaryFile<V2PNVDHashIndex>(indexPath);

    index.initialIndexGiven(idHashMap_v2p,addressHashMap_nvd);
    endTime = clock();

    cout<<"快速 导入v2p and nvd hash idx(数据地址索引)文件 完毕! 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s" << endl;


}


void Load_Hybrid_NVDG_AddressIdx_fast(){
    clock_t startTime, endTime;
    startTime = clock();
    string indexPath = getHybridNVDHashIndexInputPath();
    HybridPOINNHashIndex index = serialization::getIndexFromBinaryFile<HybridPOINNHashIndex>(indexPath);

    index.initialIndexGiven(idHashMap_v2p_hybrid,idListHash_l2p_hybrid, addressHashMap_nvd_hybrid);
    endTime = clock();

    cout<<"快速 导入hybrid nvd hash idx(数据地址索引)文件 完毕! 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s" << endl;


}


void Load_Hybrid_NVDG_AddressIdx_fast_varyingO(float rate){
    clock_t startTime, endTime;
    startTime = clock();
    string indexPath = getHybridNVDHashIndexInputPath_varying(rate);
    HybridPOINNHashIndex index = serialization::getIndexFromBinaryFile<HybridPOINNHashIndex>(indexPath);

    index.initialIndexGiven(idHashMap_v2p_hybrid,idListHash_l2p_hybrid, addressHashMap_nvd_hybrid);
    endTime = clock();

    cout<<"快速 导入hybrid nvd hash idx(数据地址索引)文件 完毕! 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s" << endl;


}




void test_serializeThenLoad_V2P_NVDG_AddressIdx(){
    clock_t startTime, endTime;
    startTime = clock();
    V2PNVDHashIndex hashIndex;

    unordered_map<int,int> addressHashMap_v2p0;
    unordered_map<int,int> addressHashMap_nvd0;

    addressHashMap_v2p0[2]= 1333;
    addressHashMap_nvd0[1]= 1222;
    hashIndex.setNVDIndex(addressHashMap_nvd0);
    hashIndex.setV2PIndex(addressHashMap_v2p0);

    string indexOutputFilePath = getNVDHashIndexInputPath();
    serialization::outputIndexToBinaryFile<V2PNVDHashIndex>(hashIndex,indexOutputFilePath);
    endTime = clock();
    cout<<"序列化 输出nvd idx(数据地址索引)文件 完毕! 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s" << endl;

    getchar();
    startTime = clock();
    string indexPath = getNVDHashIndexInputPath();
    V2PNVDHashIndex index = serialization::getIndexFromBinaryFile<V2PNVDHashIndex>(indexPath);
    endTime = clock();
    cout<<"序列化文件 导入 完毕! 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s" << endl;
    index.initialIndexGiven(idHashMap_v2p,addressHashMap_nvd);
    unordered_map<int,int>::iterator iterator_v2p;
    cout<<"addressHashMap_v2p content:"<<endl;
    for(iterator_v2p = idHashMap_v2p.begin(); iterator_v2p!=idHashMap_v2p.end(); iterator_v2p++){
        int key = iterator_v2p->first;
        int value = iterator_v2p -> second;
        cout<<"<key="<<key<<",value="<<value<<">"<<endl;
    }

    unordered_map<int,int>::iterator iterator_nvd;
    cout<<"addressHashMap_nvd content:"<<endl;
    for(iterator_nvd = addressHashMap_nvd.begin(); iterator_nvd!=addressHashMap_nvd.end(); iterator_nvd++){
        int key = iterator_nvd->first;
        int value = iterator_nvd -> second;
        cout<<"<key="<<key<<",value="<<value<<">"<<endl;
    }



}




//注意一定要在Load_V2P_NVDG_AddressIdx之后处理

void serialize_V2PID_NVDG_AddressIdx_batch(int round){
    clock_t startTime, endTime;
    startTime = clock();
    string indexOutputPre = getNVDHashIndexInputPath();
    stringstream ss;
    ss<<indexOutputPre<<round;
    string indexOutputFilePath = ss.str();
    cout<<"*****************第"<<round<<"批, (new)序列化输出V2PNVDHashIndex...";

    V2PNVDHashIndex hashIndex;
    hashIndex.setNVDIndex(addressHashMap_nvd);
    hashIndex.setV2PIndex(idHashMap_v2p);


    serialization::outputIndexToBinaryFile<V2PNVDHashIndex>(hashIndex,indexOutputFilePath);
    endTime = clock();
    cout<<"输出完毕! 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s*****************" << endl;

}


void serialize_V2PID_NVDG_AddressIdx(){
    clock_t startTime, endTime;
    startTime = clock();
    string indexOutputPre = getNVDHashIndexInputPath();
    stringstream ss;
    ss<<indexOutputPre;
    string indexOutputFilePath = ss.str();
    cout<<"*****************序列化输出V2PNVDHashIndex...";

    V2PNVDHashIndex hashIndex;
    hashIndex.setNVDIndex(addressHashMap_nvd);
    hashIndex.setV2PIndex(idHashMap_v2p);


    serialization::outputIndexToBinaryFile<V2PNVDHashIndex>(hashIndex,indexOutputFilePath);

    endTime = clock();
    cout<<"输出完毕! 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s*****************" << endl;

    /*cout<<"*****************输出Leaf2PListHashIndex...";
    startTime = clock();
    ofstream  output;
    string outputPath = getLeaf2POIListIndexInputPath();
    output.open(outputPath.c_str());
    unordered_map<int,vector<int>>:: iterator iterator_l2ps;
    for(iterator_l2ps = idListHash_l2p_hybrid.begin();iterator_l2ps!=idListHash_l2p_hybrid.end();iterator_l2ps++){
        int key = iterator_l2ps->first;
        vector<int> id_list = iterator_l2ps->second;
        output<<key;
        for(int poi_id: id_list){
            output<<' '<<poi_id;
        }
        output<<endl;
    }

    output.close();
    endTime = clock();
    cout<<"输出完毕! 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s*****************" << endl;*/



}


void serialize_Hybrid_NVDG_AddressIdx(){
    clock_t startTime, endTime;
    startTime = clock();
    string indexOutputPre = getHybridNVDHashIndexInputPath();
    stringstream ss;
    ss<<indexOutputPre;
    string indexOutputFilePath = ss.str();
    cout<<"*****************序列化输出 HybridNVDHashIndex...";

    HybridPOINNHashIndex hashIndex;
    hashIndex.setNVDIndex(addressHashMap_nvd);
    hashIndex.setV2PIndex(idHashMap_v2p);
    hashIndex.setL2PIndex(idListHash_l2p_hybrid);


    serialization::outputIndexToBinaryFile<HybridPOINNHashIndex>(hashIndex,indexOutputFilePath);

    endTime = clock();
    cout<<"输出完毕! 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s*****************" << endl;



}

void serialize_Hybrid_NVDG_AddressIdx_varying(float ratio){
    clock_t startTime, endTime;
    startTime = clock();
    string indexOutputPre = getHybridNVDHashIndexInputPath_varying(ratio);
    stringstream ss;
    ss<<indexOutputPre;
    string indexOutputFilePath = ss.str();
    cout<<"*****************序列化输出 HybridNVDHashIndex...";

    HybridPOINNHashIndex hashIndex;
    hashIndex.setNVDIndex(addressHashMap_nvd);
    hashIndex.setV2PIndex(idHashMap_v2p);
    hashIndex.setL2PIndex(idListHash_l2p_hybrid);


    serialization::outputIndexToBinaryFile<HybridPOINNHashIndex>(hashIndex,indexOutputFilePath);

    endTime = clock();
    cout<<"输出完毕! 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s*****************" << endl;



}



//混合hash索引模式


void NVD_generation_AllKeyword_V2PHash_Enhenced_road_slow(){ // leaf node level,去掉了之前多余累赘的 address_v2p对应的索引与数据文件
    clock_t iter_startTime, iter_endTime;

    clock_t start_time, end_time;

    int shi = 0;
    cout<<"Begin NVD generation!"<<endl;

    char outfileprefix[200];
    sprintf(outfileprefix,"../../../data/%s/NVD/%s",dataset_name,road_data);
    cout<<"目录为:"<<outfileprefix<<endl;

    //打开NVDGraph文件
    char nvd_file_name[200];
    sprintf(nvd_file_name,"%s.nvd",outfileprefix);
    NVD_ADJGraph = fopen(nvd_file_name, "wb+");

    int key_nvd = 0;
    int key_v2p = 0;
    char idx_file_prefix[200];
    char nvd_idx_file_name[200];
    //char v2p_idx_file_name[200];
    sprintf(idx_file_prefix,"../../../data/%s/NVD",dataset_name);
    sprintf(nvd_idx_file_name,"%s/NVD.idx",idx_file_prefix);
    //sprintf(v2p_idx_file_name,"%s/v2p.idx",idx_file_prefix);



    int term_size = getTermSize();

    double dure_runtime = 0;

    int round = 0;
    for(int i=1;i<=term_size;i++){

        iter_startTime = clock();

        int target_term = i;
        float term_weight = getTermIDFWeight(target_term);
        int _size = getTermOjectInvListSize(target_term);
        cout<<"----------------------term"<<target_term<<", idf_score="<<term_weight<<", size="<<_size<<"----------------------"<<endl;
        if(_size<= posting_size_threshold) continue;  //注意！
        vector<int> poi_ids = getTermOjectInvList(target_term);
        cout<<"building NVD for objects in its posting list (size="<<poi_ids.size()<<")"<<endl;

        int nodeCnt = 0;
        float Max_D = 25000;

        start_time = clock();
        bool* vis  = new bool[VertexNum];
        bool* rearch  = new bool[VertexNum];
        float*  pointval = new float[VertexNum];
        for(int i=0;i<VertexNum;i++){
            vis[i] = false;
            rearch[i] = false;  //jins: 重要！
        }

        bool* poi_rearch = new bool [poi_num];
        for(int i=0;i<poi_num;i++){
            poi_rearch[i] = false;
        }

        //hashmap: <vertex_id, poi_id>
        int* vertex_point_pair = new int[VertexNum];//idx : vertex_id, value : poi_id;
        for(int i=0;i<VertexNum;i++)
            vertex_point_pair[i] = -100;  //jins: 重要！

        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"new bool [] 用时："<<dure_runtime<<"ms!"<<endl;


        map<int,vector<Result>> resultCache;    // <edge_id, pre_midresults>
        map<int,bool> EdgeVisFlag;    // <edge_id, pre_midresults>

        float maxDisValue = 9999999.9;

        priority_queue<Entry_NVD> Q;



        //map<int, int> vertex_point_pairMap;
        //hashmap: <poi_id,  vertex_set> : road network partition by generating points (pois in the posting list for a specific keywords)
        map<int, set<int>> NVD_vertex_partition;
        //hashmap: <poi_id, adjacent_poi list>  : adjacent graph of NVD
        map<int, set<int>> NVD_adjacent_graph;


        //initialize for each generating points
        start_time = clock();
        for(int p: poi_ids){
            //cout<<"需要检索p"<<p<<endl;
            /* if(p==13470){
                 cout<<"find p==13470!"<<endl;
             }*/
            POI poi = getPOIFromO2UOrgLeafData(p);
            //printPOIInfo(poi);
            //for poi.Ni
            int p_Ni = poi.Ni;
            float p_dis = poi.dis;
            int scale_dist = (int)(poi.dis*dist_scale_dj);
            pointval[p_Ni] = p_dis; //vertex_point_pairMap[p_Ni] = poi.id;
            Entry_NVD element1(p_Ni, poi.id,p_dis,scale_dist);
            rearch[p_Ni]=true;  //重要！
            Q.push(element1);

            //for poi.Nj
            int p_Nj = poi.Nj;
            float p_dis2 =  poi.dist-poi.dis;
            int scale_dis2 = (int)(p_dis2*dist_scale_dj);
            pointval[p_Nj] = poi.dist-poi.dis; //vertex_point_pairMap[p_Ni] = poi.id;
            Entry_NVD element2(p_Nj, poi.id,p_dis2,scale_dis2);
            rearch[p_Nj]=true;  //重要！
            Q.push(element2);
        }
        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"从inv list获得 相关 poi用时："<<dure_runtime <<"ms!"<<endl;


        start_time = clock();
        while (!Q.empty()) {
            Entry_NVD entry = Q.top();Q.pop();
            int cur_vertex = entry.vertex_id;
            int cur_poi = entry.poi_id;
            vertex_point_pair[cur_vertex] = cur_poi;  //当前vertex就处在当前poi的 Voronoi Cell中（因为距离最小）
            vis[cur_vertex] = true;
            //int AdjGrpAddr, AdjListSize;
            //获取当前顶点的邻接表信息
            vector<int> current_adj = adjList[cur_vertex];
            vector<float> tmpAdjW = adjWList[cur_vertex];
            for(int i=0;i<current_adj.size();i++){
                int next_vertex =  current_adj[i];
                float edge_dist = tmpAdjW[i];
                //edge  e = EdgeMap[getKey(cur_vertex, next_vertex)];
                if(vis[next_vertex]==true) continue;
                if(rearch[next_vertex]==false){
                    pointval[next_vertex] = entry.dist + edge_dist;
                    rearch[next_vertex] = true; //重要
                    int sca_dist = (int)(dist_scale_dj*pointval[next_vertex]);
                    Entry_NVD ele_next(next_vertex,entry.poi_id,pointval[next_vertex],sca_dist);
                    Q.push(ele_next);
                }
                else if(entry.dist + edge_dist < pointval[next_vertex]){
                    pointval[next_vertex] = entry.dist + edge_dist;
                    int sca_dist = (int)(dist_scale_dj*pointval[next_vertex]);
                    rearch[next_vertex] = true; //重要
                    Entry_NVD entry_next(next_vertex,entry.poi_id,pointval[next_vertex],sca_dist);
                    Q.push(entry_next);
                }
            }
        }
        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"dj下路网遍历用时："<<dure_runtime <<" ms!"<<endl;
        //结束路网边遍历
        //统计每个顶点对应的 Voronoi Cell 中, 将v加入到poi的关联顶点列表中

        start_time = clock();
        map<int, set<int>> NVD_AdjGraph;
        map<int, set<int>> NVD_BorderMap;
        for(int i=0;i<VertexNum;i++){
            int v = i;
            int poi = vertex_point_pair[v];
            //printf("v%d 对应的最近邻poi为o %d\n",v,poi);
            NVD_vertex_partition[poi].insert(v);
            vector<int> v_adj = adjList[v];
            for(int u : v_adj){
                int poi2 = vertex_point_pair[u];
                if(poi == poi2){
                    continue;
                }
                else{  //建立generating point间的邻接关系
                    NVD_AdjGraph[poi].insert(poi2);
                    NVD_AdjGraph[poi2].insert(poi);
                    //v为poi所在NV Cell的边界点
                    NVD_BorderMap[poi].insert(v);
                    //u为poi2所在NV Cell的边界点
                    NVD_BorderMap[poi2].insert(u);
                }
            }
        }

        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"统计NVD中cell的邻接情况："<<dure_runtime <<" ms!"<<endl;



        start_time = clock();
        //向-----------nvd文件中写入的数据内容------------
        //先向nvd btree 写入poi在当前term下将要在nvd文件中写入数据的地址
        //int key_poi_keyword = term*poi_num + poi.id;
        //nvd_btree.insert(key_poi_keyword, key_nvd);
        int term = target_term;
        for(int poi:poi_ids){
            int addr = key_nvd; //poi_Address[poi];
            //cout<<"addr="<<addr<<endl;
            //0. 先向nvd btree 写入poi在当前term下将要在nvd文件中写入数据的地址
            int key_poi_keyword = term*poi_num + poi;
            ////char key_poi_keyword_char[16] = { 0 };
            ////sprintf(key_poi_keyword_char,"%d",key_poi_keyword);
            ////nvd_btree.insert(key_poi_keyword_char, addr);
            addressHashMap_nvd[key_poi_keyword] = addr;
            //nvdIdx_outputFile <<key_poi_keyword<<' '<<addr<<endl;

            //1. 写入自身id
            fwrite(&poi, 1, sizeof(int), NVD_ADJGraph);
            key_nvd += sizeof(int);
            //2. 写入毗邻的cell的个数
            int _adjCellSize = NVD_AdjGraph[poi].size();
            fwrite(&_adjCellSize, 1, sizeof(int), NVD_ADJGraph);
            key_nvd += sizeof(int);
            if(NVD_AdjGraph[poi].size()>0){
                //printf("兴趣点o%d 的NVD结果是：border_size=%d, \n adjacent_poi size = %d, 毗邻的poi分别为：",poi, NVD_BorderMap[poi].size(), NVD_AdjGraph[poi].size());
                //3. 写入各个邻居的id
                for(int id: NVD_AdjGraph[poi]){
                    //cout<<"o"<<id<<", ";
                    fwrite(&id, 1, sizeof(int), NVD_ADJGraph);
                    key_nvd += sizeof(int);
                }
            }
            /*if(poi==0&&term==2){
                cout<<"poi==0&&term==2!"<<endl;
                getchar();
            }*/

        }
        //------------向v2p 文件中写入数据内容, 将相关key与<poi_id1,id2...>记录到-vertexHashMap或leafHashMap中---------------
        for(int node_id =0;node_id<GTree.size();node_id++){
            if(GTree[node_id].isleaf==false) continue;
            int current_leaf = node_id;

            ////1. 先判断当前leaf中是否有该term的相关poi
            set<int> poiNN_set;
            vector<int> pois_within = getObjectTermRelatedEntry(term,current_leaf);
            if(pois_within.size()>0){
                for(int poiNN: pois_within){
                    poiNN_set.insert(poiNN);
                }
            }
            ////3. 将当前leaf的border点对应的poiNN加入poiNN_set中
            for(int border: GTree[current_leaf].borders){
                int border_poiNN = vertex_point_pair[border];
                poiNN_set.insert(border_poiNN);
            }
            int within_size = poiNN_set.size();
            //cout<<"within_size= "<<within_size<<endl;
            ////4. 混合索引模式：取决于poiNN_set中所包含的不同poi的个数粒度，根据稀疏度判定是写leaf信息（poiNN_set），还是写vertex 信息（poiNN）
            if(poiNN_set.size()<= rau){  //关键词信息稀疏型节点，以leaf node整体为单位，统一写信息
                //lalala
                int key_node_keyword = term * GTree.size() + current_leaf;
                if(idListHash_l2p_hybrid.count(key_node_keyword)) continue;   //同节点中的其他顶点已经对该leaf的信息进行过写入了
                vector<int> poiNNList;
                for(int _poi: poiNN_set){
                    poiNNList.push_back(_poi);
                }
                idListHash_l2p_hybrid[key_node_keyword] = poiNNList;
            }
            else{  //关键词信息密集型节点，需对leaf中所有vertex的信息进行记录
                for(int vertex_id: GTree[current_leaf].vetexSet){
                    int addr = key_v2p;
                    //0. 向v2p btree 写入vertex poi pair数据地址
                    int key_vertex_keyword = term*VertexNum + vertex_id;
                    //1. 向idHashMap_v2p 表中记录 当前vertex 对应的最近邻 poi的id
                    int poi_id = vertex_point_pair[vertex_id];
                    idHashMap_v2p[key_vertex_keyword] = poi_id;
#ifdef  TRACKV2P
                    cout<<"顶点v"<<i<<"在关键词t"<<term<<"下的NN 为p"<<poi_id<<endl;
#endif
                    //fwrite(&poi,1, sizeof(int), VertexHashFile);
                    key_v2p += sizeof(int);

                }


            }



        }
        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"nvd 信息写入用时："<<dure_runtime <<"ms!"<<endl;
        iter_endTime = clock();

        cout<<"***finish nvd build for term"<<target_term<<" spend "<<(double)(iter_endTime - iter_startTime) / CLOCKS_PER_SEC * 1000<<"ms!"<<endl;

        delete []vis;
        delete []poi_rearch;
        delete []rearch;
        delete []pointval;
        delete []vertex_point_pair;

        //getchar();




    }

    //serialize_V2PID_NVDG_AddressIdx();
    serialize_Hybrid_NVDG_AddressIdx();

    //if(term_size % wordInter !=0)   //map中还有小部分内容要输出


    //关闭数据文件
    fclose(NVD_ADJGraph);
    //fclose(VertexHashFile);
    //ala


}


void NVD_generation_AllKeyword_V2PHash_Enhenced_road0(){ // leaf node level,去掉了之前多余累赘的 address_v2p对应的索引与数据文件
    clock_t iter_startTime, iter_endTime;

    clock_t start_time, end_time;

    int shi = 0;
    cout<<"Begin NVD generation!"<<endl;

    char outfileprefix[200];
    sprintf(outfileprefix,"../../../data/%s/NVD/%s",dataset_name,road_data);
    cout<<"目录为:"<<outfileprefix<<endl;

    //打开NVDGraph文件
    char nvd_file_name[200];
    sprintf(nvd_file_name,"%s.nvd",outfileprefix);
    NVD_ADJGraph = fopen(nvd_file_name, "wb+");

    int key_nvd = 0;
    int key_v2p = 0;
    char idx_file_prefix[200];
    char nvd_idx_file_name[200];
    //char v2p_idx_file_name[200];
    sprintf(idx_file_prefix,"../../../data/%s/NVD",dataset_name);
    sprintf(nvd_idx_file_name,"%s/NVD.idx",idx_file_prefix);
    //sprintf(v2p_idx_file_name,"%s/v2p.idx",idx_file_prefix);



    int term_size = getTermSize();

    double dure_runtime = 0;

    int round = 0;
    for(int i=1;i<=term_size;i++){

        iter_startTime = clock();

        int target_term = i;
        float term_weight = getTermIDFWeight(target_term);
        int _size = invListOfPOI[target_term].size();
        if(_size<= posting_size_threshold) continue;  //注意！
        cout<<"----------------------term"<<target_term<<", idf_score="<<term_weight<<", size="<<_size<<"----------------------"<<endl;
        vector<int> poi_ids = invListOfPOI[target_term];
        cout<<"building NVD for objects in its posting list (size="<<poi_ids.size()<<")"<<endl;

        int nodeCnt = 0;
        float Max_D = 25000;

        start_time = clock();
        bool* vis  = new bool[VertexNum];
        bool* rearch  = new bool[VertexNum];
        float*  pointval = new float[VertexNum];
        for(int i=0;i<VertexNum;i++){
            vis[i] = false;
            rearch[i] = false;  //jins: 重要！
        }

        bool* poi_rearch = new bool [poi_num];
        for(int i=0;i<poi_num;i++){
            poi_rearch[i] = false;
        }

        //hashmap: <vertex_id, poi_id>
        int* vertex_point_pair = new int[VertexNum];//idx : vertex_id, value : poi_id;
        for(int i=0;i<VertexNum;i++)
            vertex_point_pair[i] = -100;  //jins: 重要！

        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"new bool [] 用时："<<dure_runtime<<"ms!"<<endl;


        map<int,vector<Result>> resultCache;    // <edge_id, pre_midresults>
        map<int,bool> EdgeVisFlag;    // <edge_id, pre_midresults>

        float maxDisValue = 9999999.9;

        priority_queue<Entry_NVD> Q;



        //hashmap: <poi_id,  vertex_set> : road network partition by generating points (pois in the posting list for a specific keywords)
        map<int, set<int>> NVD_vertex_partition;
        //hashmap: <poi_id, adjacent_poi list>  : adjacent graph of NVD
        map<int, set<int>> NVD_adjacent_graph;


        //initialize for each generating points
        start_time = clock();
        for(int p: poi_ids){
            //cout<<"需要检索p"<<p<<endl;
            /* if(p==13470){
                 cout<<"find p==13470!"<<endl;
             }*/
            //POI poi = getPOIFromO2UOrgLeafData(p);
            POI poi = POIs[p];
            //printPOIInfo(poi);
            //for poi.Ni
            int p_Ni = poi.Ni;
            float p_dis = poi.dis;
            int scale_dist = (int)(poi.dis*dist_scale_dj);
            pointval[p_Ni] = p_dis; //vertex_point_pairMap[p_Ni] = poi.id;
            Entry_NVD element1(p_Ni, poi.id,p_dis,scale_dist);
            rearch[p_Ni]=true;  //重要！
            Q.push(element1);

            //for poi.Nj
            int p_Nj = poi.Nj;
            float p_dis2 =  poi.dist-poi.dis;
            int scale_dis2 = (int)(p_dis2*dist_scale_dj);
            pointval[p_Nj] = poi.dist-poi.dis; //vertex_point_pairMap[p_Ni] = poi.id;
            Entry_NVD element2(p_Nj, poi.id,p_dis2,scale_dis2);
            rearch[p_Nj]=true;  //重要！
            Q.push(element2);
        }
        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"从inv list获得 相关 poi用时："<<dure_runtime <<"ms!"<<endl;


        start_time = clock();
        while (!Q.empty()) {
            Entry_NVD entry = Q.top();Q.pop();
            int cur_vertex = entry.vertex_id;
            int cur_poi = entry.poi_id;
            vertex_point_pair[cur_vertex] = cur_poi;  //当前vertex就处在当前poi的 Voronoi Cell中（因为距离最小）
            vis[cur_vertex] = true;
            //int AdjGrpAddr, AdjListSize;
            //获取当前顶点的邻接表信息
            vector<int> current_adj = adjList[cur_vertex];
            vector<float> tmpAdjW = adjWList[cur_vertex];
            for(int i=0;i<current_adj.size();i++){
                int next_vertex =  current_adj[i];
                float edge_dist = tmpAdjW[i];
                //edge  e = EdgeMap[getKey(cur_vertex, next_vertex)];
                if(vis[next_vertex]==true) continue;
                if(rearch[next_vertex]==false){
                    pointval[next_vertex] = entry.dist + edge_dist;
                    rearch[next_vertex] = true; //重要
                    int sca_dist = (int)(dist_scale_dj*pointval[next_vertex]);
                    Entry_NVD ele_next(next_vertex,entry.poi_id,pointval[next_vertex],sca_dist);
                    Q.push(ele_next);
                }
                else if(entry.dist + edge_dist < pointval[next_vertex]){
                    pointval[next_vertex] = entry.dist + edge_dist;
                    int sca_dist = (int)(dist_scale_dj*pointval[next_vertex]);
                    rearch[next_vertex] = true; //重要
                    Entry_NVD entry_next(next_vertex,entry.poi_id,pointval[next_vertex],sca_dist);
                    Q.push(entry_next);
                }
            }
        }
        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"dj下路网遍历用时："<<dure_runtime <<" ms!"<<endl;
        //结束路网边遍历
        //统计每个顶点对应的 Voronoi Cell 中, 将v加入到poi的关联顶点列表中
        start_time = clock();
        map<int, set<int>> NVD_AdjGraph;
        map<int, set<int>> NVD_BorderMap;
        for(int i=0;i<VertexNum;i++){
            int v = i;
            int poi = vertex_point_pair[v];
            //printf("v%d 对应的最近邻poi为o %d\n",v,poi);
            NVD_vertex_partition[poi].insert(v);
            vector<int> v_adj = adjList[v];
            for(int u : v_adj){
                int poi2 = vertex_point_pair[u];
                if(poi == poi2){
                    continue;
                }
                else{  //建立generating point间的邻接关系
                    NVD_AdjGraph[poi].insert(poi2);
                    NVD_AdjGraph[poi2].insert(poi);
                    //v为poi所在NV Cell的边界点
                    NVD_BorderMap[poi].insert(v);
                    //u为poi2所在NV Cell的边界点
                    NVD_BorderMap[poi2].insert(u);
                }
            }
        }

        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"统计NVD中cell的邻接情况："<<dure_runtime <<" ms!"<<endl;



        start_time = clock();
        //向-----------nvd文件中写入的数据内容------------
        //int key_poi_keyword = term*poi_num + poi.id;
        //nvd_btree.insert(key_poi_keyword, key_nvd);
        int term = target_term;
        for(int poi:poi_ids){
            int addr = key_nvd; //poi_Address[poi];
            //cout<<"addr="<<addr<<endl;
            //0. 先向nvd btree 写入poi在当前term下将要在nvd文件中写入数据的地址
            int key_poi_keyword = term*poi_num + poi;
            ////char key_poi_keyword_char[16] = { 0 };
            ////sprintf(key_poi_keyword_char,"%d",key_poi_keyword);
            ////nvd_btree.insert(key_poi_keyword_char, addr);
            addressHashMap_nvd[key_poi_keyword] = addr;
            //nvdIdx_outputFile <<key_poi_keyword<<' '<<addr<<endl;

            //1. 写入自身id
            fwrite(&poi, 1, sizeof(int), NVD_ADJGraph);
            key_nvd += sizeof(int);
            //2. 写入毗邻的cell的个数
            int _adjCellSize = NVD_AdjGraph[poi].size();
            fwrite(&_adjCellSize, 1, sizeof(int), NVD_ADJGraph);
            key_nvd += sizeof(int);
            if(NVD_AdjGraph[poi].size()>0){
                //printf("兴趣点o%d 的NVD结果是：border_size=%d, \n adjacent_poi size = %d, 毗邻的poi分别为：",poi, NVD_BorderMap[poi].size(), NVD_AdjGraph[poi].size());
                //3. 写入各个邻居的id
                for(int id: NVD_AdjGraph[poi]){
                    //cout<<"o"<<id<<", ";
                    fwrite(&id, 1, sizeof(int), NVD_ADJGraph);
                    key_nvd += sizeof(int);
                }
            }


        }
        //------------向v2p 文件中写入数据内容, 将相关key与<poi_id1,id2...>记录到-vertexHashMap或leafHashMap中---------------
        for(int node_id =0;node_id<GTree.size();node_id++){
            if(GTree[node_id].isleaf==false) continue;
            int current_leaf = node_id;

            ////1. 先判断当前leaf中是否有该term的相关poi
            set<int> poiNN_set;
            //vector<int> pois_within = getObjectTermRelatedEntry(term,current_leaf);
            set<int> pois_within = GTree[current_leaf].inverted_list_o[term];
            if(pois_within.size()>0){
                for(int poiNN: pois_within){
                    poiNN_set.insert(poiNN);
                }
            }
            ////3. 将当前leaf的border点对应的poiNN加入poiNN_set中
            for(int border: GTree[current_leaf].borders){
                int border_poiNN = vertex_point_pair[border];
                poiNN_set.insert(border_poiNN);
            }
            int within_size = poiNN_set.size();
            //cout<<"within_size= "<<within_size<<endl;
            ////4. 混合索引模式：取决于poiNN_set中所包含的不同poi的个数粒度，根据稀疏度判定是写leaf信息（poiNN_set），还是写vertex 信息（poiNN）
            if(poiNN_set.size()<= rau){  //关键词信息稀疏型节点，以leaf node整体为单位，统一写信息
                //cout<<"关键词信息疏松型节点，within_size="<<within_size<<endl;
                int key_node_keyword = term * GTree.size() + current_leaf;
                if(idListHash_l2p_hybrid.count(key_node_keyword)) continue;   //同节点中的其他顶点已经对该leaf的信息进行过写入了
                vector<int> poiNNList;
                for(int _poi: poiNN_set){
                    poiNNList.push_back(_poi);
                }
                idListHash_l2p_hybrid[key_node_keyword] = poiNNList;
            }
            else{  //关键词信息密集型节点，需对leaf中所有vertex的信息进行记录
                //cout<<"关键词信息密集型节点,within_size="<<within_size<<endl;
                for(int vertex_id: GTree[current_leaf].vetexSet){
                    int addr = key_v2p;
                    //0. 向v2p btree 写入vertex poi pair数据地址
                    int key_vertex_keyword = term*VertexNum + vertex_id;
                    //1. 向idHashMap_v2p 表中记录 当前vertex 对应的最近邻 poi的id
                    int poi_id = vertex_point_pair[vertex_id];
                    idHashMap_v2p[key_vertex_keyword] = poi_id;
#ifdef  TRACKV2P
                    cout<<"顶点v"<<i<<"在关键词t"<<term<<"下的NN 为p"<<poi_id<<endl;
#endif
                    //fwrite(&poi,1, sizeof(int), VertexHashFile);
                    key_v2p += sizeof(int);

                }


            }



        }
        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"nvd 信息写入用时："<<dure_runtime <<"ms!"<<endl;
        iter_endTime = clock();

        cout<<"***finish nvd build for term"<<target_term<<" spend "<<(double)(iter_endTime - iter_startTime) / CLOCKS_PER_SEC * 1000<<"ms!"<<endl;

        delete []vis;
        delete []poi_rearch;
        delete []rearch;
        delete []pointval;
        delete []vertex_point_pair;

        //getchar();




    }

    //serialize_V2PID_NVDG_AddressIdx();
    serialize_Hybrid_NVDG_AddressIdx();

    //if(term_size % wordInter !=0)   //map中还有小部分内容要输出


    //关闭数据文件
    fclose(NVD_ADJGraph);
    //fclose(VertexHashFile);
    //ala


}

///IL-NVD构建函数
void NVD_generation_AllKeyword_V2PHash_Enhenced_road(){ // leaf node level,去掉了之前多余累赘的 address_v2p对应的索引与数据文件
    clock_t iter_startTime, iter_endTime;

    clock_t start_time, end_time;

    int shi = 0;
    cout<<"Begin NVD generation!"<<endl;

    char outfileprefix[200];
    sprintf(outfileprefix,"../../../data/%s/NVD/%s",dataset_name,road_data);
    cout<<"目录为:"<<outfileprefix<<endl;

    //打开NVDGraph文件
    char nvd_file_name[200];
    sprintf(nvd_file_name,"%s.nvd",outfileprefix);
    NVD_ADJGraph = fopen(nvd_file_name, "wb+");

    int key_nvd = 0;
    int key_v2p = 0;
    char idx_file_prefix[200];
    char nvd_idx_file_name[200];
    //char v2p_idx_file_name[200];
    sprintf(idx_file_prefix,"../../../data/%s/NVD",dataset_name);
    sprintf(nvd_idx_file_name,"%s/NVD.idx",idx_file_prefix);
    //sprintf(v2p_idx_file_name,"%s/v2p.idx",idx_file_prefix);



    int term_size = getTermSize();

    double dure_runtime = 0;

    int round = 0;
    for(int i=1;i<=term_size;i++){

        iter_startTime = clock();

        int target_term = i;
        float term_weight = getTermIDFWeight(target_term);
        int _size = invListOfPOI[target_term].size();
        if(_size<= posting_size_threshold) continue;  //注意！
        cout<<"----------------------term"<<target_term<<", idf_score="<<term_weight<<", size="<<_size<<"----------------------"<<endl;
        vector<int> poi_ids = invListOfPOI[target_term];
        cout<<"building NVD for objects in its posting list (size="<<poi_ids.size()<<")"<<endl;

        int nodeCnt = 0;
        float Max_D = 25000;

        start_time = clock();
        bool* vis  = new bool[VertexNum];
        bool* rearch  = new bool[VertexNum];
        float*  pointval = new float[VertexNum];
        for(int i=0;i<VertexNum;i++){
            vis[i] = false;
            rearch[i] = false;  //jins: 重要！
        }

        bool* poi_rearch = new bool [poi_num];
        for(int i=0;i<poi_num;i++){
            poi_rearch[i] = false;
        }

        //hashmap: <vertex_id, poi_id>
        int* vertex_point_pair = new int[VertexNum];//idx : vertex_id, value : poi_id;
        for(int i=0;i<VertexNum;i++)
            vertex_point_pair[i] = -100;  //jins: 重要！

        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"new bool [] 用时："<<dure_runtime<<"ms!"<<endl;


        map<int,vector<Result>> resultCache;    // <edge_id, pre_midresults>
        map<int,bool> EdgeVisFlag;    // <edge_id, pre_midresults>

        float maxDisValue = 9999999.9;

        priority_queue<Entry_NVD> Q;



        //hashmap: <poi_id,  vertex_set> : road network partition by generating points (pois in the posting list for a specific keywords)
        map<int, set<int>> NVD_vertex_partition;
        //hashmap: <poi_id, adjacent_poi list>  : adjacent graph of NVD
        map<int, set<int>> NVD_adjacent_graph;


        //initialize for each generating points
        start_time = clock();
        for(int p: poi_ids){
            //cout<<"需要检索p"<<p<<endl;

            //POI poi = getPOIFromO2UOrgLeafData(p);
            POI poi = POIs[p];
            //printPOIInfo(poi);
            int p_Ni = poi.Ni;
            float p_dis = poi.dis;
            int scale_dist = (int)(poi.dis*dist_scale_dj);
            pointval[p_Ni] = p_dis; //vertex_point_pairMap[p_Ni] = poi.id;
            Entry_NVD element1(p_Ni, poi.id,p_dis,scale_dist);
            rearch[p_Ni]=true;  //重要！
            Q.push(element1);

            //for poi.Nj
            int p_Nj = poi.Nj;
            float p_dis2 =  poi.dist-poi.dis;
            int scale_dis2 = (int)(p_dis2*dist_scale_dj);
            pointval[p_Nj] = poi.dist-poi.dis; //vertex_point_pairMap[p_Ni] = poi.id;
            Entry_NVD element2(p_Nj, poi.id,p_dis2,scale_dis2);
            rearch[p_Nj]=true;  //重要！
            Q.push(element2);
        }
        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"从inv list获得 相关 poi用时："<<dure_runtime <<"ms!"<<endl;


        start_time = clock();
        while (!Q.empty()) {
            Entry_NVD entry = Q.top();Q.pop();
            int cur_vertex = entry.vertex_id;
            int cur_poi = entry.poi_id;
            vertex_point_pair[cur_vertex] = cur_poi;  //当前vertex就处在当前poi的 Voronoi Cell中（因为距离最小）
            vis[cur_vertex] = true;
            //int AdjGrpAddr, AdjListSize;
            //获取当前顶点的邻接表信息
            vector<int> current_adj = adjList[cur_vertex];
            vector<float> tmpAdjW = adjWList[cur_vertex];
            for(int i=0;i<current_adj.size();i++){
                int next_vertex =  current_adj[i];
                float edge_dist = tmpAdjW[i];
                //edge  e = EdgeMap[getKey(cur_vertex, next_vertex)];
                if(vis[next_vertex]==true) continue;
                if(rearch[next_vertex]==false){
                    pointval[next_vertex] = entry.dist + edge_dist;
                    rearch[next_vertex] = true; //重要
                    int sca_dist = (int)(dist_scale_dj*pointval[next_vertex]);
                    Entry_NVD ele_next(next_vertex,entry.poi_id,pointval[next_vertex],sca_dist);
                    Q.push(ele_next);
                }
                else if(entry.dist + edge_dist < pointval[next_vertex]){
                    pointval[next_vertex] = entry.dist + edge_dist;
                    int sca_dist = (int)(dist_scale_dj*pointval[next_vertex]);
                    rearch[next_vertex] = true; //重要
                    Entry_NVD entry_next(next_vertex,entry.poi_id,pointval[next_vertex],sca_dist);
                    Q.push(entry_next);
                }
            }
        }
        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"dj下路网遍历用时："<<dure_runtime <<" ms!"<<endl;
        //结束路网边遍历
        //统计每个顶点对应的 Voronoi Cell 中, 将v加入到poi的关联顶点列表中
        start_time = clock();
        map<int, set<int>> NVD_AdjGraph;
        map<int, set<int>> NVD_BorderMap;
        for(int i=0;i<VertexNum;i++){
            int v = i;
            int poi = vertex_point_pair[v];
            //printf("v%d 对应的最近邻poi为o %d\n",v,poi);
            NVD_vertex_partition[poi].insert(v);
            vector<int> v_adj = adjList[v];
            for(int u : v_adj){
                int poi2 = vertex_point_pair[u];
                if(poi == poi2){
                    continue;
                }
                else{  //建立generating point间的邻接关系
                    NVD_AdjGraph[poi].insert(poi2);
                    NVD_AdjGraph[poi2].insert(poi);
                    //v为poi所在NV Cell的边界点
                    NVD_BorderMap[poi].insert(v);
                    //u为poi2所在NV Cell的边界点
                    NVD_BorderMap[poi2].insert(u);
                }
            }
        }
        //修正一部分孤立的poi,使他们在NVD上有邻接poi
        for(int p_id:poi_ids){
            if(p_id==936111){
                cout<<"find p==936111!"<<endl;
            }
            if(NVD_vertex_partition[p_id].size()==0){  ///该兴趣点 p 所在边的两个顶点，也未将p 作为 其 NN，
                POI poi = POIs[p_id];
                //printPOIInfo(poi);
                int p_Ni = poi.Ni;
                int p_Nj = poi.Nj;
                int poi_adj1 = vertex_point_pair[p_Ni];
                int poi_adj2 = vertex_point_pair[p_Nj];
                if(poi_adj1!=p_id){
                    NVD_AdjGraph[p_id].insert(poi_adj1);  //认为p 在 NVD grah上与 其所在边顶点 v(u)的NN:p',互为邻居
                    NVD_AdjGraph[poi_adj1].insert(p_id);
                }
                if(poi_adj2!=p_id){
                    NVD_AdjGraph[p_id].insert(poi_adj2);
                    NVD_AdjGraph[poi_adj2].insert(p_id);
                }

            }
        }

        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"统计NVD中cell的邻接情况："<<dure_runtime <<" ms!"<<endl;



        start_time = clock();
        //向-----------nvd文件中写入的数据内容------------
        //int key_poi_keyword = term*poi_num + poi.id;
        //nvd_btree.insert(key_poi_keyword, key_nvd);
        int term = target_term;
        for(int poi:poi_ids){
            if(poi==936111){
                cout<<"find p==936111!"<<endl;
            }
            int addr = key_nvd; //poi_Address[poi];
            //cout<<"addr="<<addr<<endl;
            //0. 先向nvd btree 写入poi在当前term下将要在nvd文件中写入数据的地址
            int key_poi_keyword = term*poi_num + poi;
            ////char key_poi_keyword_char[16] = { 0 };
            ////sprintf(key_poi_keyword_char,"%d",key_poi_keyword);
            ////nvd_btree.insert(key_poi_keyword_char, addr);
            addressHashMap_nvd[key_poi_keyword] = addr;
            //nvdIdx_outputFile <<key_poi_keyword<<' '<<addr<<endl;

            //1. 写入自身id
            fwrite(&poi, 1, sizeof(int), NVD_ADJGraph);
            key_nvd += sizeof(int);
            //2. 写入毗邻的cell的个数
            int _adjCellSize = NVD_AdjGraph[poi].size();
            fwrite(&_adjCellSize, 1, sizeof(int), NVD_ADJGraph);
            key_nvd += sizeof(int);
            if(NVD_AdjGraph[poi].size()>0){
                //3. 写入各个邻居的id
                for(int id: NVD_AdjGraph[poi]){
                    fwrite(&id, 1, sizeof(int), NVD_ADJGraph);
                    key_nvd += sizeof(int);
                    if(poi==936111){
                         printf("o%d\n",id);
                    }
                }
            }
           /*if(poi==936111){
                cout<<"find p==936111!"<<endl;
                getchar();
            }*/


        }
        //------------向v2p 文件中写入数据内容, 将相关key与<poi_id1,id2...>记录到-vertexHashMap或leafHashMap中---------------
        for(int node_id =0;node_id<GTree.size();node_id++){
            if(GTree[node_id].isleaf==false) continue;
            int current_leaf = node_id;

            ////1. 先判断当前leaf中是否有该term的相关poi
            set<int> poiNN_set;
            //vector<int> pois_within = getObjectTermRelatedEntry(term,current_leaf);
            set<int> pois_within = GTree[current_leaf].inverted_list_o[term];
            if(pois_within.size()>0){
                for(int poiNN: pois_within){
                    poiNN_set.insert(poiNN);
                }
            }
            ////3. 将当前leaf的border点对应的poiNN加入poiNN_set中
            for(int border: GTree[current_leaf].borders){
                int border_poiNN = vertex_point_pair[border];
                poiNN_set.insert(border_poiNN);
            }
            int within_size = poiNN_set.size();
            //cout<<"within_size= "<<within_size<<endl;
            ////4. 混合索引模式：取决于poiNN_set中所包含的不同poi的个数粒度，根据稀疏度判定是写leaf信息（poiNN_set），还是写vertex 信息（poiNN）
            if(poiNN_set.size()<= rau){  //关键词信息稀疏型节点，以leaf node整体为单位，统一写信息
                //cout<<"关键词信息疏松型节点，within_size="<<within_size<<endl;
                int key_node_keyword = term * GTree.size() + current_leaf;
                if(idListHash_l2p_hybrid.count(key_node_keyword)) continue;   //同节点中的其他顶点已经对该leaf的信息进行过写入了
                vector<int> poiNNList;
                for(int _poi: poiNN_set){
                    poiNNList.push_back(_poi);
                }
                idListHash_l2p_hybrid[key_node_keyword] = poiNNList;
            }
            else{  //关键词信息密集型节点，需对leaf中所有vertex的信息进行记录
                //cout<<"关键词信息密集型节点,within_size="<<within_size<<endl;
                for(int vertex_id: GTree[current_leaf].vetexSet){
                    int addr = key_v2p;
                    //0. 向v2p btree 写入vertex poi pair数据地址
                    int key_vertex_keyword = term*VertexNum + vertex_id;
                    //1. 向idHashMap_v2p 表中记录 当前vertex 对应的最近邻 poi的id
                    int poi_id = vertex_point_pair[vertex_id];
                    idHashMap_v2p[key_vertex_keyword] = poi_id;
#ifdef  TRACKV2P
                    ///cout<<"顶点v"<<vertex_id<<"在关键词t"<<term<<"下的NN 为p"<<poi_id<<endl;
#endif
                    //fwrite(&poi,1, sizeof(int), VertexHashFile);
                    key_v2p += sizeof(int);

                }


            }



        }
        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"nvd 信息写入用时："<<dure_runtime <<"ms!"<<endl;
        iter_endTime = clock();

        cout<<"***finish nvd build for term"<<target_term<<" spend "<<(double)(iter_endTime - iter_startTime) / CLOCKS_PER_SEC * 1000<<"ms!"<<endl;

        delete []vis;
        delete []poi_rearch;
        delete []rearch;
        delete []pointval;
        delete []vertex_point_pair;

        //getchar();




    }

    //serialize_V2PID_NVDG_AddressIdx();
    serialize_Hybrid_NVDG_AddressIdx();

    //if(term_size % wordInter !=0)   //map中还有小部分内容要输出


    //关闭数据文件
    fclose(NVD_ADJGraph);
    //fclose(VertexHashFile);
    //ala


}


void NVD_generation_AllKeyword_V2PHash_Enhenced_road_varyingO(float o_ratio){ // leaf node level,去掉了之前多余累赘的 address_v2p对应的索引与数据文件
    clock_t iter_startTime, iter_endTime;

    clock_t start_time, end_time;

    int shi = 0;
    cout<<"Begin NVD generation!"<<endl;

    char outfileprefix[200]; int _tmp_ratio = (int)(o_ratio*100);
    sprintf(outfileprefix,"../../../data/%s/NVD/VaryingO/%d%/%s",dataset_name,_tmp_ratio,road_data);
    cout<<"目录为:"<<outfileprefix<<endl;

    //打开NVDGraph文件
    char nvd_file_name[200];
    sprintf(nvd_file_name,"%s.nvd",outfileprefix);
    NVD_ADJGraph = fopen(nvd_file_name, "wb+");

    int key_nvd = 0;
    int key_v2p = 0;
    char idx_file_prefix[200];
    char nvd_idx_file_name[200];
    //char v2p_idx_file_name[200];
    sprintf(idx_file_prefix,"../../../data/%s/NVD",dataset_name);
    sprintf(nvd_idx_file_name,"%s/NVD.idx",idx_file_prefix);
    //sprintf(v2p_idx_file_name,"%s/v2p.idx",idx_file_prefix);



    int term_size = getTermSize();

    double dure_runtime = 0;

    int round = 0;
    for(int i=1;i<=term_size;i++){

        iter_startTime = clock();

        int target_term = i;
        float term_weight = getTermIDFWeight(target_term);
        int _size = invListOfPOI[target_term].size();
        if(_size<= posting_size_threshold) continue;  //注意！
        cout<<"----------------------term"<<target_term<<", idf_score="<<term_weight<<", size="<<_size<<"----------------------"<<endl;
        vector<int> poi_ids = invListOfPOI[target_term];
        cout<<"building NVD for objects in its posting list (size="<<poi_ids.size()<<")"<<endl;

        int nodeCnt = 0;
        float Max_D = 25000;

        start_time = clock();
        bool* vis  = new bool[VertexNum];
        bool* rearch  = new bool[VertexNum];
        float*  pointval = new float[VertexNum];
        for(int i=0;i<VertexNum;i++){
            vis[i] = false;
            rearch[i] = false;  //jins: 重要！
        }

        bool* poi_rearch = new bool [poi_num];
        for(int i=0;i<poi_num;i++){
            poi_rearch[i] = false;
        }

        //hashmap: <vertex_id, poi_id>
        int* vertex_point_pair = new int[VertexNum];//idx : vertex_id, value : poi_id;
        for(int i=0;i<VertexNum;i++)
            vertex_point_pair[i] = -100;  //jins: 重要！

        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"new bool [] 用时："<<dure_runtime<<"ms!"<<endl;


        map<int,vector<Result>> resultCache;    // <edge_id, pre_midresults>
        map<int,bool> EdgeVisFlag;    // <edge_id, pre_midresults>

        float maxDisValue = 9999999.9;

        priority_queue<Entry_NVD> Q;



        //hashmap: <poi_id,  vertex_set> : road network partition by generating points (pois in the posting list for a specific keywords)
        map<int, set<int>> NVD_vertex_partition;
        //hashmap: <poi_id, adjacent_poi list>  : adjacent graph of NVD
        map<int, set<int>> NVD_adjacent_graph;


        //initialize for each generating points
        start_time = clock();
        for(int p: poi_ids){

            POI poi = POIs[p];
            //printPOIInfo(poi);
            int p_Ni = poi.Ni;
            float p_dis = poi.dis;
            int scale_dist = (int)(poi.dis*dist_scale_dj);
            pointval[p_Ni] = p_dis; //vertex_point_pairMap[p_Ni] = poi.id;
            Entry_NVD element1(p_Ni, poi.id,p_dis,scale_dist);
            rearch[p_Ni]=true;  //重要！
            Q.push(element1);

            //for poi.Nj
            int p_Nj = poi.Nj;
            float p_dis2 =  poi.dist-poi.dis;
            int scale_dis2 = (int)(p_dis2*dist_scale_dj);
            pointval[p_Nj] = poi.dist-poi.dis; //vertex_point_pairMap[p_Ni] = poi.id;
            Entry_NVD element2(p_Nj, poi.id,p_dis2,scale_dis2);
            rearch[p_Nj]=true;  //重要！
            Q.push(element2);
        }
        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"从inv list获得 相关 poi用时："<<dure_runtime <<"ms!"<<endl;


        start_time = clock();
        while (!Q.empty()) {
            Entry_NVD entry = Q.top();Q.pop();
            int cur_vertex = entry.vertex_id;
            int cur_poi = entry.poi_id;
            vertex_point_pair[cur_vertex] = cur_poi;  //当前vertex就处在当前poi的 Voronoi Cell中（因为距离最小）
            vis[cur_vertex] = true;
            //int AdjGrpAddr, AdjListSize;
            //获取当前顶点的邻接表信息
            vector<int> current_adj = adjList[cur_vertex];
            vector<float> tmpAdjW = adjWList[cur_vertex];
            for(int i=0;i<current_adj.size();i++){
                int next_vertex =  current_adj[i];
                float edge_dist = tmpAdjW[i];
                //edge  e = EdgeMap[getKey(cur_vertex, next_vertex)];
                if(vis[next_vertex]==true) continue;
                if(rearch[next_vertex]==false){
                    pointval[next_vertex] = entry.dist + edge_dist;
                    rearch[next_vertex] = true; //重要
                    int sca_dist = (int)(dist_scale_dj*pointval[next_vertex]);
                    Entry_NVD ele_next(next_vertex,entry.poi_id,pointval[next_vertex],sca_dist);
                    Q.push(ele_next);
                }
                else if(entry.dist + edge_dist < pointval[next_vertex]){
                    pointval[next_vertex] = entry.dist + edge_dist;
                    int sca_dist = (int)(dist_scale_dj*pointval[next_vertex]);
                    rearch[next_vertex] = true; //重要
                    Entry_NVD entry_next(next_vertex,entry.poi_id,pointval[next_vertex],sca_dist);
                    Q.push(entry_next);
                }
            }
        }
        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"dj下路网遍历用时："<<dure_runtime <<" ms!"<<endl;
        //结束路网边遍历
        //统计每个顶点对应的 Voronoi Cell 中, 将v加入到poi的关联顶点列表中
        start_time = clock();
        map<int, set<int>> NVD_AdjGraph;
        map<int, set<int>> NVD_BorderMap;
        for(int i=0;i<VertexNum;i++){
            int v = i;
            int poi = vertex_point_pair[v];
            //printf("v%d 对应的最近邻poi为o %d\n",v,poi);
            NVD_vertex_partition[poi].insert(v);
            vector<int> v_adj = adjList[v];
            for(int u : v_adj){
                int poi2 = vertex_point_pair[u];
                if(poi == poi2){
                    continue;
                }
                else{  //建立generating point间的邻接关系
                    NVD_AdjGraph[poi].insert(poi2);
                    NVD_AdjGraph[poi2].insert(poi);
                    //v为poi所在NV Cell的边界点
                    NVD_BorderMap[poi].insert(v);
                    //u为poi2所在NV Cell的边界点
                    NVD_BorderMap[poi2].insert(u);
                }
            }
        }
        //修正一部分孤立的poi,使他们在NVD上有邻接poi
        for(int p_id:poi_ids){
            if(NVD_vertex_partition[p_id].size()==0){
                POI poi = POIs[p_id];
                //printPOIInfo(poi);
                int p_Ni = poi.Ni;
                int p_Nj = poi.Nj;
                int poi_adj1 = vertex_point_pair[p_Ni];
                int poi_adj2 = vertex_point_pair[p_Nj];
                if(poi_adj1!=p_id){
                    NVD_AdjGraph[p_id].insert(poi_adj1);
                    NVD_AdjGraph[poi_adj1].insert(p_id);
                }
                if(poi_adj2!=p_id){
                    NVD_AdjGraph[p_id].insert(poi_adj2);
                    NVD_AdjGraph[poi_adj2].insert(p_id);
                }

            }
        }

        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"统计NVD中cell的邻接情况："<<dure_runtime <<" ms!"<<endl;



        start_time = clock();
        //向-----------nvd文件中写入的数据内容------------
        //int key_poi_keyword = term*poi_num + poi.id;
        //nvd_btree.insert(key_poi_keyword, key_nvd);
        int term = target_term;
        for(int poi:poi_ids){
            int addr = key_nvd; //poi_Address[poi];
            //cout<<"addr="<<addr<<endl;
            //0. 先向nvd btree 写入poi在当前term下将要在nvd文件中写入数据的地址
            int key_poi_keyword = term*poi_num + poi;
            ////char key_poi_keyword_char[16] = { 0 };
            ////sprintf(key_poi_keyword_char,"%d",key_poi_keyword);
            ////nvd_btree.insert(key_poi_keyword_char, addr);
            addressHashMap_nvd[key_poi_keyword] = addr;
            //nvdIdx_outputFile <<key_poi_keyword<<' '<<addr<<endl;

            //1. 写入自身id
            fwrite(&poi, 1, sizeof(int), NVD_ADJGraph);
            key_nvd += sizeof(int);
            //2. 写入毗邻的cell的个数
            int _adjCellSize = NVD_AdjGraph[poi].size();
            fwrite(&_adjCellSize, 1, sizeof(int), NVD_ADJGraph);
            key_nvd += sizeof(int);
            if(NVD_AdjGraph[poi].size()>0){
                /*if(target_term == 14){
                    printf("兴趣点o%d 的NVD结果是：border_size=%d, adjacent_poi size = %d\n",poi, NVD_BorderMap[poi].size(), NVD_AdjGraph[poi].size());
                }*/
                //3. 写入各个邻居的id
                for(int id: NVD_AdjGraph[poi]){
                    //cout<<"o"<<id<<", ";
                    fwrite(&id, 1, sizeof(int), NVD_ADJGraph);
                    key_nvd += sizeof(int);
                }
            }


        }
        //------------向v2p 文件中写入数据内容, 将相关key与<poi_id1,id2...>记录到-vertexHashMap或leafHashMap中---------------
        for(int node_id =0;node_id<GTree.size();node_id++){
            if(GTree[node_id].isleaf==false) continue;
            int current_leaf = node_id;

            ////1. 先判断当前leaf中是否有该term的相关poi
            set<int> poiNN_set;
            //vector<int> pois_within = getObjectTermRelatedEntry(term,current_leaf);
            set<int> pois_within = GTree[current_leaf].inverted_list_o[term];
            if(pois_within.size()>0){
                for(int poiNN: pois_within){
                    poiNN_set.insert(poiNN);
                }
            }
            ////3. 将当前leaf的border点对应的poiNN加入poiNN_set中
            for(int border: GTree[current_leaf].borders){
                int border_poiNN = vertex_point_pair[border];
                poiNN_set.insert(border_poiNN);
            }
            int within_size = poiNN_set.size();
            //cout<<"within_size= "<<within_size<<endl;
            ////4. 混合索引模式：取决于poiNN_set中所包含的不同poi的个数粒度，根据稀疏度判定是写leaf信息（poiNN_set），还是写vertex 信息（poiNN）
            if(poiNN_set.size()<= rau){  //关键词信息稀疏型节点，以leaf node整体为单位，统一写信息
                //cout<<"关键词信息疏松型节点，within_size="<<within_size<<endl;
                int key_node_keyword = term * GTree.size() + current_leaf;
                if(idListHash_l2p_hybrid.count(key_node_keyword)) continue;   //同节点中的其他顶点已经对该leaf的信息进行过写入了
                vector<int> poiNNList;
                for(int _poi: poiNN_set){
                    poiNNList.push_back(_poi);
                }
                idListHash_l2p_hybrid[key_node_keyword] = poiNNList;
            }
            else{  //关键词信息密集型节点，需对leaf中所有vertex的信息进行记录
                //cout<<"关键词信息密集型节点,within_size="<<within_size<<endl;
                for(int vertex_id: GTree[current_leaf].vetexSet){
                    int addr = key_v2p;
                    //0. 向v2p btree 写入vertex poi pair数据地址
                    int key_vertex_keyword = term*VertexNum + vertex_id;
                    //1. 向idHashMap_v2p 表中记录 当前vertex 对应的最近邻 poi的id
                    int poi_id = vertex_point_pair[vertex_id];
                    idHashMap_v2p[key_vertex_keyword] = poi_id;
#ifdef  TRACKV2P
                    cout<<"顶点v"<<i<<"在关键词t"<<term<<"下的NN 为p"<<poi_id<<endl;
#endif
                    //fwrite(&poi,1, sizeof(int), VertexHashFile);
                    key_v2p += sizeof(int);

                }


            }



        }
        end_time = clock();
        dure_runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
        cout<<"nvd 信息写入用时："<<dure_runtime <<"ms!"<<endl;
        iter_endTime = clock();

        cout<<"***finish nvd build for term"<<target_term<<" spend "<<(double)(iter_endTime - iter_startTime) / CLOCKS_PER_SEC * 1000<<"ms!"<<endl;

        delete []vis;
        delete []poi_rearch;
        delete []rearch;
        delete []pointval;
        delete []vertex_point_pair;

        //getchar();




    }

    serialize_Hybrid_NVDG_AddressIdx_varying(o_ratio);

    //if(term_size % wordInter !=0)   //map中还有小部分内容要输出


    //关闭数据文件
    fclose(NVD_ADJGraph);
    //fclose(VertexHashFile);
    //ala


}




int getNVDAddress_by_Hash_vertexOnly(int poi_id, int term){
    //计算key的大小
    clock_t startTime, endTime;
    startTime = clock();
    int key;
    key = term*poi_num + poi_id;

    int addr_logic = -1;
    if(addressHashMap_nvd.count(key)){
        addr_logic = addressHashMap_nvd[key];
    }
    else{
        printf("find no such key=%d in NVD.idx\n");
        exit(-1);
    }
    endTime = clock();
#ifdef DEBUG
    cout << "getAddress in NVD.idx runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;
#endif
    return  addr_logic;
}




double dist_total_runtime = 0;
int dist_cmp_count = 0;




map<int, double> usr_SBCache;

//core
double getUserLCL_NVD(User user, int K, int a, double alpha){
    if(usr_SBCache.count(user.id))
        return usr_SBCache[user.id];


    priority_queue<nodeDistType_NVD> Queue;
    vector<POI> polled_object_lists;
    vector<float> polled_distance_lists;

    priority_queue<Result> resultFinal;
    int vi = user.Ni; int vj = user.Nj;
    map<int,bool> poiIsVist;

    map<int, vector<int>> posting_Map;
    map<int, int> posting_ListSize_Map;

    map<int, int> posting_ListCurrentPos_Map;

    clock_t startTime, endTime;
    clock_t idx_start, idx_end; double idx_time=0; int idx_access_count =0;
    clock_t nvd_start, nvd_end; double nvd_time=0;
    clock_t io_start, io_end; double io_time=0;
    clock_t queue_start, queue_end; double queue_time=0;
    clock_t dist_start, dist_end; double dist_time=0;

    startTime = clock();//计时开始
    bool hasFrequentKey = false;
                            ///printUsrInfo(user);
    for(int term: user.keywords){
        ///cout<<"t"<<term<<endl;
        int posting_size = getTermOjectInvListSize(term);
        //如果term为低频词汇
        if(posting_size<= posting_size_threshold){  //不大于 rau,即认为是低频词汇
            posting_ListSize_Map[term] = posting_size;
            vector<int> posting_list = getTermOjectInvList(term);
            posting_Map[term] = posting_list;
            io_start = clock();
            for(int poi_id: posting_list){
                //cout<<"access p"<<poi_id<<endl;
                POI poi = getPOIFromO2UOrgLeafData(poi_id);
#ifdef TRACK
                printPOIInfo(poi);
#endif

                if(poiIsVist.count(poi_id)) continue;
                //计算距离
                float dist = getDistanceUpper_Oracle(user,poi);  //getDistance(user,poi);
                nodeDistType_NVD tmp1(-1, poi_id, dist, term);
                Queue.push(tmp1);
                posting_ListCurrentPos_Map[term] = 1;
                polled_object_lists.push_back(poi);  //加入相关列表
                polled_distance_lists.push_back(dist);
                poiIsVist[poi_id] =true;

            }

        }
        else{ //否，则为高频词汇
            hasFrequentKey = true;
            clock_t s1,e1;
            int poi_id = getNNPOI_By_HybridVertex_Keyword(vi,term);  //core
            ///cout<<"poi_id="<<poi_id<<endl;
            if(poi_id!=-1){  //vi 在term下有对应的 NN poi信息
#ifdef TRACK
                cout<<"v"<<vi<<"在t"<<term<<"下对应的 NN poi为：p"<<poi_id<<endl;
#endif
                if(poiIsVist.count(poi_id)) continue;  //poiIsVist[poi_id]==true
                POI poi = getPOIFromO2UOrgLeafData(poi_id);
                float dist = getDistanceUpper_Oracle(user,poi); //SPSP(user.Ni,user.Nj)/1.0+user.dis+poi.dis;//getDistance(user,poi);
                dist_cmp_count++;
                nodeDistType_NVD tmp1(-1, poi_id, dist, term);
                Queue.push(tmp1);

                polled_object_lists.push_back(poi);  //加入相关列表
                polled_distance_lists.push_back(dist);
                poiIsVist[poi_id] =true;
            }
            else{ //从leaf hash中寻找！
                vector<int> pois;
                int leaf_node = Nodes[vi].gtreepath.back();
                pois = getNNPOIList_by_Hybrid_Hash(leaf_node,term);
#ifdef TRACK
                cout<<"v"<<vi<<"在t"<<term<<"下对应的 NN po ilist为:"<<endl;
                printElements(pois);
#endif
                for(int poi_id2: pois){  //将leaf node 对应的poiNNSet中元素加入
                    if(poiIsVist.count(poi_id2)) continue;  //poiIsVist[poi_id]==true
                    POI poi2 = getPOIFromO2UOrgLeafData(poi_id2);

                    float dist = getDistanceUpper_Oracle(user,poi2); //SPSP(user.Ni,user.Nj)/1.0+user.dis+poi.dis;//getDistance(user,poi);
                    dist_cmp_count++;
                    nodeDistType_NVD tmp1(-1, poi_id2, dist, term);
                    Queue.push(tmp1);
                    polled_object_lists.push_back(poi2);  //加入相关列表
                    polled_distance_lists.push_back(dist);
                    poiIsVist[poi_id2] =true;
                }
#ifdef TRACK
                cout<<"v"<<vi<<"在t"<<term<<"下对应的 NN poi组为："<<endl;
                printElements(pois);
#endif
            }

#else
            int poi_id = getNNPOI_By_Vertex_Keyword(vi,term);  //先得到距离vertex最近的带关键词term的poi

            if(poi_id ==-100) continue;   //此时说明 所有term相关的poi 相对user都不可达
            if(poiIsVist.count(poi_id)) continue;  //poiIsVist[poi_id]==true
            POI poi = getPOIFromO2UOrgLeafData(poi_id);
            float dist = getDistanceUpper_Oracle(user,poi); //SPSP(user.Ni,user.Nj)/1.0+user.dis+poi.dis;//getDistance(user,poi);
            dist_cmp_count++;
            nodeDistType_NVD tmp1(-1, poi_id, dist, term);
            Queue.push(tmp1);
            polled_object_lists.push_back(poi);  //加入相关列表
            polled_distance_lists.push_back(dist);
            poiIsVist[poi_id] =true;
#endif
        }



    }
    endTime = clock();//计时结束////


    //填充与排序
    if(Queue.size()==0) return -1;
    startTime = clock();//计时开始
    double score_k = 0;
    if(hasFrequentKey==false&&polled_object_lists.size()<K)  //若user包含的全是低频词汇，且与之关键词相关的object总体个数也小于K
        return 0;
    int aa_size = polled_object_lists.size();
    while(polled_object_lists.size()<K){
        nodeDistType_NVD ele = Queue.top(); Queue.pop();
        int _poi_id = ele.poi_id; int _term = ele.term_id;
#ifdef TRACK
        cout<<"p"<<_poi_id<<"出队列, Queue.size="<<Queue.size()<<endl;
        //如果term为低频词汇
#endif
        if(posting_ListSize_Map.count(_term)){
            //cout<<"低频词汇"<<endl;
            //加入下一个
            int _current_pos = posting_ListCurrentPos_Map[_term];
            vector<int> posting = posting_Map[_term];
            int _poi_id = posting[_current_pos];
            if(poiIsVist.count(_poi_id)) continue;
            POI poi = getPOIFromO2UOrgLeafData(_poi_id);
            //距离计算
            dist_start = clock();
            float dist = getDistanceUpper_Oracle(user,poi); //SPSP(user.Ni,user.Nj)/1.0+user.dis+poi.dis;//getDistance(user, poi);
            dist_end = clock();
            dist_total_runtime += (double)(dist_end - dist_start);
            dist_cmp_count++;

            nodeDistType_NVD poi_entry(-1,_poi_id,dist,_term);
            Queue.push(poi_entry);
            polled_object_lists.push_back(poi);  //加入相关列表
            polled_distance_lists.push_back(dist);
            //标记posting列表，当前队首下标
            posting_ListCurrentPos_Map[_term] = _current_pos++;
            poiIsVist[_poi_id] =true; //标记已被加入
            //cout<<"继续加入p"<<_poi_id<<endl;

        }
        else{  //为高频词汇
            //cout<<"load poi"<<_poi_id<<"ADJ..."<<endl;
            int _id = _poi_id;
            vector<int> poi_ADJ = getPOIAdj_NVD_By_Keyword(_id,_term);  //core
            //cout<<"获得poi"<<_poi_id<<"ADJ成功, size="<<poi_ADJ.size()<<endl;
            //NVD graph中相邻poi加入：
            for(int _poi_id : poi_ADJ){
                if(poiIsVist.count(_poi_id)) continue;
                if(_poi_id<0 ||_poi_id>poi_num) continue;
                //cout<<"load 兴趣点"<<_poi_id<<"数据..."<<endl;
                POI poi = getPOIFromO2UOrgLeafData(_poi_id); int id = poi.id;
                //cout<<"finish load!"<<endl;
                dist_start = clock();
                float dist = getDistanceUpper_Oracle(user,poi); //SPSP(user.Ni,user.Nj)/1.0+user.dis+poi.dis;//getDistance(user,poi);
                dist_end = clock();
                dist_total_runtime += (double)(dist_end - dist_start);
                dist_cmp_count++;

                //cout<<"dist computation runtime:"<< (double)(dist_end - dist_start) / CLOCKS_PER_SEC * 1000000 << "us" << endl;
                //getchar();

                //cout<<"dist="<<dist<<endl;
                //int adj_poi = getPOIAdj_NVD_By_Keyword(poi_id,term);
                nodeDistType_NVD poi_entry(-1, id, dist, _term);
                Queue.push(poi_entry);
                polled_object_lists.push_back(poi);  //加入相关列表
                polled_distance_lists.push_back(dist);
                poiIsVist[id] =true;
                //cout<<"继续加入p"<<_poi_id<<endl;
            }
        }


    }
    endTime = clock();//计时结束
    ////cout << "filling runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;

    startTime = clock();
    int _size = polled_object_lists.size();
    for(int j = 0; j<_size;j++){
        POI p = polled_object_lists[j];
        float dist = polled_distance_lists[j];
        double simT = textRelevance(user.keywords, p.keywords);
        double social_textual_score = alpha*1 + (1-alpha)*simT;

        double gsk_score = social_textual_score / (1.0+ dist);
        Result tmpRlt(p.id, p.Ni, p.Nj, p.dis, p.dist, gsk_score, p.keywords);
        resultFinal.push(tmpRlt);
    }
    endTime = clock();
    ////cout << "score computing and runking is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;


    if(resultFinal.size()==K){
        double lcl_rt = resultFinal.top().score;
        score_k = lcl_rt;
    }
    else if(resultFinal.size()>K){
        while(resultFinal.size()>K){
            //cout<<resultFinal.top().score<<endl;
            resultFinal.pop();
        }
        double lcl_rt = resultFinal.top().score;
        score_k = lcl_rt;
    }
    else
        score_k = 0.0;

    usr_SBCache[user.id] = score_k;
    return score_k;

}

double getUserLCL_NVDFaster(User user, int K, int a, double alpha){
    //return 0.0;

    if(usr_SBCache.count(user.id))
        return usr_SBCache[user.id];


    priority_queue<nodeDistType_NVD> Queue;
    vector<POI> polled_object_lists;
    vector<float> polled_distance_lists;

    priority_queue<Result> resultFinal;
    int vi = user.Ni; int vj = user.Nj;
    map<int,bool> poiIsVist;

    map<int, vector<int>> posting_Map;
    map<int, int> posting_ListSize_Map;

    map<int, int> posting_ListCurrentPos_Map;

    clock_t startTime, endTime;
    clock_t idx_start, idx_end; double idx_time=0; int idx_access_count =0;
    clock_t nvd_start, nvd_end; double nvd_time=0;
    clock_t io_start, io_end; double io_time=0;
    clock_t queue_start, queue_end; double queue_time=0;
    clock_t dist_start, dist_end; double dist_time=0;

    startTime = clock();//计时开始
    bool hasFrequentKey = false;
    set<int> evaluate_terms; int mostFrequentKey; int maxFrequency = 0;
    for(int term: user.keywords){
        int posting_size = getTermOjectInvListSize(term);
        if(posting_size>maxFrequency){
            mostFrequentKey = term;
            maxFrequency = posting_size;
        }
        if(posting_size>posting_size_threshold){
            evaluate_terms.insert(term);
        }
        if(evaluate_terms.size()==2) break;
    }
    if(evaluate_terms.size()==0) return 0.0;
    else
        evaluate_terms.insert(mostFrequentKey);


    for(int term: evaluate_terms){
        int posting_size = getTermOjectInvListSize(term);
        //如果term为低频词汇
        if(posting_size<= posting_size_threshold){  //不大于 rau,即认为是低频词汇
            posting_ListSize_Map[term] = posting_size;
            vector<int> posting_list = getTermOjectInvList(term);
            posting_Map[term] = posting_list;
            io_start = clock();
            for(int poi_id: posting_list){
                //cout<<"access p"<<poi_id<<endl;
                POI poi = getPOIFromO2UOrgLeafData(poi_id);
                //printPOIInfo(poi);

                if(poiIsVist.count(poi_id)) continue;
                //计算距离
                float dist = getDistanceUpper_Oracle(user,poi);  //getDistance(user,poi);
                nodeDistType_NVD tmp1(-1, poi_id, dist, term);
                Queue.push(tmp1);
                posting_ListCurrentPos_Map[term] = 1;
                polled_object_lists.push_back(poi);  //加入相关列表
                polled_distance_lists.push_back(dist);
                poiIsVist[poi_id] =true;
            }

        }
        else{ //否，则为高频词汇
            hasFrequentKey = true;
            clock_t s1,e1;
#ifdef HybridHash
             cout<<"v"<<vi<<",t"<<term<<endl;
            int poi_id = getNNPOI_By_HybridVertex_Keyword(vi,term);
            if(poi_id!=-1){  //vi 在term下有对应的 NN poi信息
#ifdef TRACK
                cout<<"v"<<vi<<"在t"<<term<<"下对应的 NN poi为：p"<<poi_id<<endl;
#endif
                if(poiIsVist.count(poi_id)) continue;  //poiIsVist[poi_id]==true
                POI poi = getPOIFromO2UOrgLeafData(poi_id);
                float dist = getDistanceUpper_Oracle(user,poi); //SPSP(user.Ni,user.Nj)/1.0+user.dis+poi.dis;//getDistance(user,poi);
                dist_cmp_count++;
                nodeDistType_NVD tmp1(-1, poi_id, dist, term);
                Queue.push(tmp1);
                polled_object_lists.push_back(poi);  //加入相关列表
                polled_distance_lists.push_back(dist);
                poiIsVist[poi_id] =true;
            }
            else{ //从leaf hash中寻找！
                vector<int> pois;
                int leaf_node = Nodes[vi].gtreepath.back();
                pois = getNNPOIList_by_Hybrid_Hash(leaf_node,term);
                for(int poi_id2: pois){  //将leaf node 对应的poiNNSet中元素加入
                    if(poiIsVist.count(poi_id2)) continue;  //poiIsVist[poi_id]==true
                    POI poi2 = getPOIFromO2UOrgLeafData(poi_id2);
                    float dist = getDistanceUpper_Oracle(user,poi2); //SPSP(user.Ni,user.Nj)/1.0+user.dis+poi.dis;//getDistance(user,poi);
                    dist_cmp_count++;
                    nodeDistType_NVD tmp1(-1, poi_id2, dist, term);
                    Queue.push(tmp1);
                    polled_object_lists.push_back(poi2);  //加入相关列表
                    polled_distance_lists.push_back(dist);
                    poiIsVist[poi_id2] =true;
                }
#ifdef TRACK
                cout<<"v"<<vi<<"在t"<<term<<"下对应的 NN poi组为："<<endl;
                printElements(pois);
#endif
            }

#else
            int poi_id = getNNPOI_By_Vertex_Keyword(vi,term);  //先得到距离vertex最近的带关键词term的poi
            if(poi_id ==-100) continue;   //此时说明 所有term相关的poi 相对user都不可达
            if(poiIsVist.count(poi_id)) continue;  //poiIsVist[poi_id]==true
            POI poi = getPOIFromO2UOrgLeafData(poi_id);
            float dist = getDistanceUpper_Oracle(user,poi); //SPSP(user.Ni,user.Nj)/1.0+user.dis+poi.dis;//getDistance(user,poi);
            dist_cmp_count++;
            nodeDistType_NVD tmp1(-1, poi_id, dist, term);
            Queue.push(tmp1);
            polled_object_lists.push_back(poi);  //加入相关列表
            polled_distance_lists.push_back(dist);
            poiIsVist[poi_id] =true;
#endif
        }



    }
    endTime = clock();//计时结束////


    //填充与排序
    if(Queue.size()==0) return -1;
    startTime = clock();//计时开始
    double score_k = 0;
    if(hasFrequentKey==false&&polled_object_lists.size()<K)  //若user包含的全是低频词汇，且与之关键词相关的object总体个数也小于K
        return 0;
    int aa_size = polled_object_lists.size();
    while(polled_object_lists.size()<K){
        nodeDistType_NVD ele = Queue.top(); Queue.pop();
        int _poi_id = ele.poi_id; int _term = ele.term_id;
        //cout<<"p"<<_poi_id<<"出队列"<<endl;
        //如果term为低频词汇
        if(posting_ListSize_Map.count(_term)){
            //cout<<"低频词汇"<<endl;
            //加入下一个
            int _current_pos = posting_ListCurrentPos_Map[_term];
            vector<int> posting = posting_Map[_term];
            int _poi_id = posting[_current_pos];
            if(poiIsVist.count(_poi_id)) continue;
            POI poi = getPOIFromO2UOrgLeafData(_poi_id);
            //距离计算
            dist_start = clock();
            float dist = getDistanceUpper_Oracle(user,poi); //SPSP(user.Ni,user.Nj)/1.0+user.dis+poi.dis;//getDistance(user, poi);
            dist_end = clock();
            dist_total_runtime += (double)(dist_end - dist_start);
            dist_cmp_count++;

            nodeDistType_NVD poi_entry(-1,_poi_id,dist,_term);
            Queue.push(poi_entry);
            polled_object_lists.push_back(poi);  //加入相关列表
            polled_distance_lists.push_back(dist);
            //标记posting列表，当前队首下标
            posting_ListCurrentPos_Map[_term] = _current_pos++;
            poiIsVist[_poi_id] =true; //标记已被加入
            //cout<<"继续加入p"<<_poi_id<<endl;

        }
        else{  //为高频词汇
            //cout<<"load poi"<<_poi_id<<"ADJ..."<<endl;
            vector<int> poi_ADJ = getPOIAdj_NVD_By_Keyword(_poi_id,_term);
            //cout<<"获得poi"<<_poi_id<<"ADJ成功, size="<<poi_ADJ.size()<<endl;
            //NVD graph中相邻poi加入：
            for(int _poi_id : poi_ADJ){
                if(poiIsVist.count(_poi_id)) continue;
                if(_poi_id<0 ||_poi_id>poi_num) continue;
                //cout<<"load 兴趣点"<<_poi_id<<"数据..."<<endl;
                POI poi = getPOIFromO2UOrgLeafData(_poi_id); int id = poi.id;
                //cout<<"finish load!"<<endl;
                dist_start = clock();
                float dist = getDistanceUpper_Oracle(user,poi); //SPSP(user.Ni,user.Nj)/1.0+user.dis+poi.dis;//getDistance(user,poi);
                dist_end = clock();
                dist_total_runtime += (double)(dist_end - dist_start);
                dist_cmp_count++;

                //cout<<"dist computation runtime:"<< (double)(dist_end - dist_start) / CLOCKS_PER_SEC * 1000000 << "us" << endl;
                //getchar();

                //cout<<"dist="<<dist<<endl;
                //int adj_poi = getPOIAdj_NVD_By_Keyword(poi_id,term);
                nodeDistType_NVD poi_entry(-1, id, dist, _term);
                Queue.push(poi_entry);
                polled_object_lists.push_back(poi);  //加入相关列表
                polled_distance_lists.push_back(dist);
                poiIsVist[id] =true;
                //cout<<"继续加入p"<<_poi_id<<endl;
            }
        }


    }
    endTime = clock();//计时结束
    ////cout << "filling runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;

    startTime = clock();
    int _size = polled_object_lists.size();
    for(int j = 0; j<_size;j++){
        POI p = polled_object_lists[j];
        float dist = polled_distance_lists[j];
        double simT = textRelevance(user.keywords, p.keywords);
        double social_textual_score = alpha*1 + (1-alpha)*simT;

        double gsk_score = social_textual_score / (1.0+ dist);
        Result tmpRlt(p.id, p.Ni, p.Nj, p.dis, p.dist, gsk_score, p.keywords);
        resultFinal.push(tmpRlt);
    }
    endTime = clock();
    ////cout << "score computing and runking is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;


    if(resultFinal.size()==K){
        double lcl_rt = resultFinal.top().score;
        score_k = lcl_rt;
    }
    else if(resultFinal.size()>K){
        while(resultFinal.size()>K){
            //cout<<resultFinal.top().score<<endl;
            resultFinal.pop();
        }
        double lcl_rt = resultFinal.top().score;
        score_k = lcl_rt;
    }
    else
        score_k = 0.0;

    usr_SBCache[user.id] = score_k;
    return score_k;

}



unordered_map<int, unordered_map<int, double>> border_KeyAwareCache;

set<int> getBorder_SB_NVD(CheckEntry border_entry,int K, int a, double alpha){

    int border_id = border_entry.id;
#ifdef TRACK_BorderSB_NVD
    cout<<"getBorder_SB_NVD for border: v"<<border_id<<endl;
#endif
    set<int> Keys = border_entry.keys_cover;
#ifdef TRACK_BorderSB_NVD
    cout<<"Keys:"; printSetElements(Keys);
#endif

    float current_dist = border_entry.dist;


    map<int, vector<int>> posting_Map;
    map<int, int> posting_ListSize_Map;

    map<int, int> posting_ListCurrentPos_Map;


    priority_queue<GlobalSBEntry>  rank_term;


    int idx = 0;

    //对各个term的 sb_list先进行计算；
    for(int term: Keys){  //对Keys中每个单词，进行评分下界队列（sb_list(term)）求解
#ifdef TRACK_BorderSB_NVD
        //cout<<"对t"<<term<<"的sb_list进行计算"<<endl;
#endif

        if(border_KeyAwareCache.count(border_id)){
            if(border_KeyAwareCache[border_id].count(term)){
                double sb_rk_t = border_KeyAwareCache[border_id][term];
                GlobalSBEntry termEntry(term,idx,sb_rk_t);
                rank_term.push(termEntry);
                //scoring_Bound_lists.push_back(scoring_Bound_list_t);
                idx++;
                continue;

            }


        }



        map<int,bool> poiIsVist;
        priority_queue<nodeDistType_NVD> H_t;
        priority_queue<SBEntry> scoring_Bound_list_t;

        int posting_size = getTermOjectInvListSize(term);
        //如果term为高频词汇
        if(posting_size > posting_size_threshold){

            vector<POI> polled_object_lists;
            vector<float> polled_distance_lists;

//#ifdef HybridHash  //在hybrid hash中读取
            cout<<"border: b"<<border_id<<",t"<<term<<endl;
            int nearest_poi_id = getNNPOI_By_HybridVertex_Keyword(border_id,term);  //先得到距离vertex最近的带关键词term的poi

            if(nearest_poi_id !=-1) {
                POI poi = getPOIFromO2UOrgLeafData(nearest_poi_id);
                float dist_upper = getUpperDistance_b2P_phl(border_id,poi);
                nodeDistType_NVD tmp1(-1, nearest_poi_id, dist_upper, term);
                H_t.push(tmp1);
            }
            else{
                int leaf_node = Nodes[border_id].gtreepath.back();
                vector<int> pois = getNNPOIList_by_Hybrid_Hash(leaf_node, term);
                for(int poi_id2:pois){
                    POI poi2 = getPOIFromO2UOrgLeafData(poi_id2);
                    float dist_upper = getUpperDistance_b2P_phl(border_id,poi2);
                    nodeDistType_NVD tmp1(-1, poi_id2, dist_upper, term);
                    H_t.push(tmp1);

                }

            }


//#endif


            while(scoring_Bound_list_t.size()<K){
                //H_t队首元素出列
                nodeDistType_NVD _tmp = H_t.top();
                H_t.pop();
                int poi_id = _tmp.poi_id;
                if(poiIsVist.count(poi_id)) continue;  //poiIsVist[poi_id]==true
                //计算scoring lower bound
                double simT_LB = 0; double simD_UB = 0; double simS_LB = 1;
                simT_LB = tfIDF_term(term);
                simD_UB = _tmp.dist;
                double social_textual_LB = alpha*simS_LB + (1-alpha)*simT_LB;
                double scoring_LB = social_textual_LB / (1+simD_UB);
                SBEntry sb_entry(poi_id,simD_UB,scoring_LB);
                scoring_Bound_list_t.push(sb_entry);
                //sb_entry.printRlt();
                poiIsVist[poi_id]= true;
                //将在NVD graph中与当前poi相邻的adj_poi加入
                vector<int> poi_ADJ = getPOIAdj_NVD_By_Keyword(poi_id,term);
                for(int adj_poi_id : poi_ADJ){
                    if(poiIsVist.count(adj_poi_id)) continue;
                    if(adj_poi_id<0 ||adj_poi_id>poi_num) continue;
                    POI adj_poi = getPOIFromO2UOrgLeafData(adj_poi_id);
                    float adj_dist_upper = getUpperDistance_b2P_phl(border_id,adj_poi);
                    nodeDistType_NVD tmp1(-1, adj_poi_id, adj_dist_upper, term);
                    H_t.push(tmp1);
                }
            }

        }
        else{ //否，则为低频词汇


            posting_ListSize_Map[term] = posting_size;
            vector<int> posting_list = getTermOjectInvList(term);
            posting_Map[term] = posting_list;
#ifdef TRACK_BorderSB_NVD
            cout<<"t"<<term<<"'s posting_list size="<<posting_list.size()<<endl;
#endif
            int _count =1;
            for(int p_id: posting_list){
                POI poi = getPOIFromO2UOrgLeafData(p_id);
                //cout<<_count<<"个poi: o"<<p_id<<endl;
                double simT_LB = 0; double simD_UB = 0; double simS_LB = 1;
                simT_LB = tfIDF_term(term);
                simD_UB = getUpperDistance_b2P_phl(border_id,poi);
                double social_textual_LB = alpha*simS_LB + (1-alpha)*simT_LB;
                double scoring_LB = social_textual_LB / (1+simD_UB);
                SBEntry sb_entry(p_id,simD_UB,scoring_LB);
                scoring_Bound_list_t.push(sb_entry);
                _count++;
            }
            while(scoring_Bound_list_t.size()>K)
                scoring_Bound_list_t.pop();


        }

        double sb_rk_t = scoring_Bound_list_t.top().score_bound;
        GlobalSBEntry termEntry(term,idx,sb_rk_t);
        rank_term.push(termEntry);
        //scoring_Bound_lists.push_back(scoring_Bound_list_t);
        idx++;
        border_KeyAwareCache[border_id][term] = sb_rk_t;


    }
    set<int> _Keys = Keys;

    //对按rank_term中关键词排序，根据scoring_Bound_lists中内容，对Key进行更新
#ifdef TRACK_BorderSB_NVD
    cout<<"对按rank_term中关键词排序，根据scoring_Bound_lists中内容，对Key进行更新"<<endl;
#endif
    while(rank_term.size()>0){
        GlobalSBEntry _termEntry = rank_term.top();
        rank_term.pop();
        int term_id = _termEntry.term_id;
        double sb_rk_t = _termEntry.sk_bound;

        //计算SB_upper(e, qo)
        double score_upper = -1;
        double max_relevance = tfIDF_termSet(_Keys);
        double min_distance = current_dist;

        double social_textual_maxscore = alpha*1 + (1.0-alpha)*max_relevance;
        score_upper = social_textual_maxscore / (1.0+ min_distance);
#ifdef TRACK_BorderSB_NVD
        cout<<"for t"<<term_id<<", sb_rk_t="<<sb_rk_t<<",score_upper="<<score_upper<<endl;
#endif
        if(sb_rk_t > score_upper){ //remove term from _Keys
#ifdef TRACK_BorderSB_NVD
            cout<<"term"<<term_id<<"可被去掉"<<endl;
            cout<<"此时_Keys:"; printSetElements(_Keys);
#endif
            set<int> :: iterator iter;
            for(iter= _Keys.begin();iter!=_Keys.end();){
                if(*iter == term_id){
                    _Keys.erase(iter);
                    break;
                }

                iter++;
            }
#ifdef TRACK_BorderSB_NVD
            cout<<"去掉后，_Keys:"; printSetElements(_Keys);
#endif
        }
        //cout<<", for term"<<term_id<<",sb_rk_t="<<sb_rk_t;
    }
#ifdef TRACK_BorderSB_NVD
    cout<<endl;
#endif

    return _Keys;

}

set<int> getBorder_SB_NVD_Cluster(CheckEntry border_entry,int K, int a, double alpha, vector<unordered_map<int,set<int>>>& nodes_uCLUSTERList){

    int border_id = border_entry.id;
    int node_id = border_entry.node_id;
#ifdef TRACK_BorderSB_NVD
    cout<<"getBorder_SB_NVD_Cluster for border: v"<<border_id<<endl;
#endif
    set<int> Keys = border_entry.keys_cover;
#ifdef TRACK_BorderSB_NVD
    cout<<"Keys:"; printSetElements(Keys);
#endif

    float current_dist = border_entry.dist;


    unordered_map<int, vector<int>> posting_Map;
    unordered_map<int, int> posting_ListSize_Map;

    unordered_map<int, int> posting_ListCurrentPos_Map;


    priority_queue<GlobalSBEntry>  rank_term;


    int idx = 0;

    //对各个term的 sb_list先进行计算；
    for(int term: Keys){  //对Keys中每个单词，进行评分下界队列（sb_list(term)）求解
#ifdef TRACK_BorderSB_NVD
        //cout<<"对t"<<term<<"的sb_list进行计算"<<endl;
#endif

        if(border_KeyAwareCache.count(border_id)){
            if(border_KeyAwareCache[border_id].count(term)){
                double sb_rk_t = border_KeyAwareCache[border_id][term];
                GlobalSBEntry termEntry(term,idx,sb_rk_t);
                rank_term.push(termEntry);
                //scoring_Bound_lists.push_back(scoring_Bound_list_t);
                idx++;
                continue;

            }


        }


        map<int,bool> poiIsVist;
        priority_queue<nodeDistType_NVD> H_t;
        priority_queue<SBEntry> scoring_Bound_list_t;

        int posting_size = getTermOjectInvListSize(term);
        //如果term为高频词汇
        if(posting_size > posting_size_threshold){

            vector<POI> polled_object_lists;
            vector<float> polled_distance_lists;

//#ifdef HybridHash  //在hybrid hash中读取
            int nearest_poi_id = getNNPOI_By_HybridVertex_Keyword(border_id,term);  //先得到距离vertex最近的带关键词term的poi

            if(nearest_poi_id !=-1) {
                POI poi = getPOIFromO2UOrgLeafData(nearest_poi_id);
                float dist_upper = getUpperDistance_b2P_phl(border_id,poi);
                nodeDistType_NVD tmp1(-1, nearest_poi_id, dist_upper, term);
                H_t.push(tmp1);
            }
            else{
                int leaf_node = Nodes[border_id].gtreepath.back();
                vector<int> pois = getNNPOIList_by_Hybrid_Hash(leaf_node, term);
                for(int poi_id2:pois){
                    POI poi2 = getPOIFromO2UOrgLeafData(poi_id2);
                    float dist_upper = getUpperDistance_b2P_phl(border_id,poi2);
                    nodeDistType_NVD tmp1(-1, poi_id2, dist_upper, term);
                    H_t.push(tmp1);

                }

            }


//#endif


            while(scoring_Bound_list_t.size()<K){
                //H_t队首元素出列
                nodeDistType_NVD _tmp = H_t.top();
                H_t.pop();
                int poi_id = _tmp.poi_id;
                if(poiIsVist.count(poi_id)) continue;  //poiIsVist[poi_id]==true
                //计算scoring lower bound
                double simT_LB = 0; double simD_UB = 0; double simS_LB = 1;
                simT_LB = tfIDF_term(term);
                simD_UB = _tmp.dist;
                double social_textual_LB = alpha*simS_LB + (1-alpha)*simT_LB;
                double scoring_LB = social_textual_LB / (1+simD_UB);
                SBEntry sb_entry(poi_id,simD_UB,scoring_LB);
                scoring_Bound_list_t.push(sb_entry);
                //sb_entry.printRlt();
                poiIsVist[poi_id]= true;
                //将在NVD graph中与当前poi相邻的adj_poi加入
                vector<int> poi_ADJ = getPOIAdj_NVD_By_Keyword(poi_id,term);  //jins
                for(int adj_poi_id : poi_ADJ){
                    if(poiIsVist.count(adj_poi_id)) continue;
                    if(adj_poi_id<0 ||adj_poi_id>poi_num) continue;
                    POI adj_poi = getPOIFromO2UOrgLeafData(adj_poi_id);
                    float adj_dist_upper = getUpperDistance_b2P_phl(border_id,adj_poi);
                    nodeDistType_NVD tmp1(-1, adj_poi_id, adj_dist_upper, term);
                    H_t.push(tmp1);
                }
            }

        }
        else{ //否，则为低频词汇


            posting_ListSize_Map[term] = posting_size;
            vector<int> posting_list = getTermOjectInvList(term);
            posting_Map[term] = posting_list;
#ifdef TRACK_BorderSB_NVD
            cout<<"t"<<term<<"'s posting_list size="<<posting_list.size()<<endl;
#endif
            int _count =1;
            for(int p_id: posting_list){
                POI poi = getPOIFromO2UOrgLeafData(p_id);
                //cout<<_count<<"个poi: o"<<p_id<<endl;
                double simT_LB = 0; double simD_UB = 0; double simS_LB = 1;
                simT_LB = tfIDF_term(term);
                simD_UB = getUpperDistance_b2P_phl(border_id,poi);
                double social_textual_LB = alpha*simS_LB + (1-alpha)*simT_LB;
                double scoring_LB = social_textual_LB / (1+simD_UB);
                SBEntry sb_entry(p_id,simD_UB,scoring_LB);
                scoring_Bound_list_t.push(sb_entry);
                _count++;
            }
            while(scoring_Bound_list_t.size()>K)
                scoring_Bound_list_t.pop();


        }

        double sb_rk_t = scoring_Bound_list_t.top().score_bound;
        GlobalSBEntry termEntry(term,idx,sb_rk_t);
        rank_term.push(termEntry);
        //scoring_Bound_lists.push_back(scoring_Bound_list_t);
        idx++;
        border_KeyAwareCache[border_id][term] = sb_rk_t;


    }
    set<int> _Keys = Keys;

    //对按rank_term中关键词排序，根据scoring_Bound_lists中内容，对Key进行更新
#ifdef TRACK_BorderSB_NVD
    cout<<"对按rank_term中关键词排序，根据scoring_Bound_lists中内容，对Key进行更新"<<endl;
#endif
    while(rank_term.size()>0){
        GlobalSBEntry _termEntry = rank_term.top();
        rank_term.pop();
        int term_id = _termEntry.term_id;
        double sb_rk_t = _termEntry.sk_bound;

        //计算SB_upper(e, qo)
        double score_upper = -1;  double score_upper2 = -1;
        set<int> cluster_Keys = nodes_uCLUSTERList[node_id][term_id];
        //printSetElements(cluster_Keys);
        //printSetElements(_Keys);
        double max_relevance = -1; double max_relevance2 = -1;
        set<int> interKeys = obtain_itersection_jins(cluster_Keys,_Keys);
        max_relevance = tfIDF_termSet(interKeys);
        //max_relevance2 = tfIDF_termSet(_Keys);

        double min_distance = current_dist;

        double social_textual_maxscore = alpha*1 + (1.0-alpha)*max_relevance;
        //double social_textual_maxscore2 = alpha*1 + (1.0-alpha)*max_relevance2;
        score_upper = social_textual_maxscore / (1.0+ min_distance);
        //score_upper2 = social_textual_maxscore2 / (1.0+ min_distance);;
#ifdef TRACK_BorderSB_NVD
        cout<<"for t"<<term_id<<", sb_rk_t="<<sb_rk_t<<",score_upper="<<score_upper<<endl;
        //cout<<"for t"<<term_id<<", sb_rk_t="<<sb_rk_t<<",score_upper="<<score_upper<<",score_upper_tight="<<score_upper2<<endl;
#endif
        if(sb_rk_t > score_upper){ //remove term from _Keys
#ifdef TRACK_BorderSB_NVD
            cout<<"term"<<term_id<<"可被去掉"<<endl;
            cout<<"此时_Keys:"; printSetElements(_Keys);
#endif
            set<int> :: iterator iter;
            for(iter= _Keys.begin();iter!=_Keys.end();){
                if(*iter == term_id){
                    _Keys.erase(iter);
                    break;
                }

                iter++;
            }
#ifdef TRACK_BorderSB_NVD
            cout<<"去掉后，_Keys:"; printSetElements(_Keys);
#endif
        }
        //cout<<", for term"<<term_id<<",sb_rk_t="<<sb_rk_t;
    }
#ifdef TRACK_BorderSB_NVD
    cout<<endl;
#endif

    return _Keys;

}

typedef vector<set<int>> UOCCURList;
typedef unordered_map<int, set<int>> UOCCURListMap;
typedef vector<unordered_map<int,set<int>>> UCLUSTERList;
typedef unordered_map<int, unordered_map<int,set<int>>> UCLUSTERListMap;


set<int> getBorder_SB_BatchNVD_Cluster(BatchCheckEntry border_entry,int K, int a, double alpha, unordered_map<int,UCLUSTERList>& pnodes_uCLUSTERListMap){  //key:poi_id, value:UCLUSTERList

    int border_id = border_entry.id;
    int node_id = border_entry.node_id;
    int poi_id = border_entry.p_id;
#ifdef TRACKBorderSB
    cout<<"getBorder_SB_BatchNVD for border: v"<<border_id<<endl;
#endif
    set<int> Keys = border_entry.keys_cover;
#ifdef TRACKBorderSB
    cout<<"Keys:"; printSetElements(Keys);
#endif

    float current_dist = border_entry.dist;


    map<int, vector<int>> posting_Map;
    map<int, int> posting_ListSize_Map;

    map<int, int> posting_ListCurrentPos_Map;


    priority_queue<GlobalSBEntry>  rank_term;


    int idx = 0;
#ifdef DEBUG
    if(460818==border_id){
        cout<<"find 460818!"<<endl;
    }
#endif

    //对各个term的 sb_list先进行计算；
    for(int term: Keys){  //对Keys中每个单词，进行评分下界队列（sb_list(term)）求解
#ifdef TRACKBorderSB
        //cout<<"对t"<<term<<"的sb_list进行计算"<<endl;
#endif

        map<int,bool> poiIsVist;
        priority_queue<nodeDistType_NVD> H_t;
        priority_queue<SBEntry> scoring_Bound_list_t;
        int _leaf_tmp; int _size;

        int posting_size = getTermOjectInvListSize(term);
        //如果term为高频词汇
        if(posting_size > posting_size_threshold){

            vector<POI> polled_object_lists;
            vector<float> polled_distance_lists;

#ifdef HybridHash   ////注意这里使用了混合哈希索引策略
            int nearest_poi_id = getNNPOI_By_HybridVertex_Keyword(border_id,term);  //先得到距离vertex最近的带关键词term的poi
            if(nearest_poi_id ==-1){
                vector<int> pois;
                int leaf_node = Nodes[border_id].gtreepath.back(); _leaf_tmp = leaf_node;
                pois = getNNPOIList_by_Hybrid_Hash(leaf_node,term); _size = pois.size();


                for(int poi_id2: pois){
                    POI poi = getPOIFromO2UOrgLeafData(poi_id2);
                    float dist_upper = getUpperDistance_b2P_phl(border_id,poi);
                    nodeDistType_NVD tmp1(-1, poi_id2, dist_upper, term);
                    H_t.push(tmp1);
#ifdef TRACK
                    cout<<"list加入o"<<poi_id2<<endl;
#endif
                }
            }
            else{
                POI poi = getPOIFromO2UOrgLeafData(nearest_poi_id);
                float dist_upper = getUpperDistance_b2P_phl(border_id,poi);
                nodeDistType_NVD tmp1(-1, nearest_poi_id, dist_upper, term);
                H_t.push(tmp1);

            }


#else
           /* int nearest_poi_id = getNNPOI_By_Vertex_Keyword(border_id,term);  //先得到距离vertex最近的带关键词term的poi
            if(nearest_poi_id ==-100) continue;   //此时说明 所有term相关的poi 相对user都不可达
            POI poi = getPOIFromO2UOrgLeafData(nearest_poi_id);
            float dist_upper = getUpperDistance_b2P_phl(border_id,poi);
            nodeDistType_NVD tmp1(-1, nearest_poi_id, dist_upper, term);
            H_t.push(tmp1);*/
#endif


            while(scoring_Bound_list_t.size()<K){
                //H_t队首元素出列
                nodeDistType_NVD _tmp = H_t.top();
                H_t.pop();
                int poi_id = _tmp.poi_id;
#ifdef TRACK
                cout<<"poi_id="<<poi_id<<", H_t size="<<H_t.size()<<", scoring_Bound_list_t.size()="<<scoring_Bound_list_t.size()<<endl;
                cout<<"pois.size="<<_size<<",leaf"<<_leaf_tmp<<",t="<<term<<endl;
                cout<<"border="<<border_id;
#endif

                if(poiIsVist.count(poi_id)) continue;  //poiIsVist[poi_id]==true
                //计算scoring lower bound
                double simT_LB = 0; double simD_UB = 0; double simS_LB = 1;
                simT_LB = tfIDF_term(term);
                simD_UB = _tmp.dist;
                double social_textual_LB = alpha*simS_LB + (1-alpha)*simT_LB;
                double scoring_LB = social_textual_LB / (1+simD_UB);
                SBEntry sb_entry(poi_id,simD_UB,scoring_LB);
                scoring_Bound_list_t.push(sb_entry);
                //sb_entry.printRlt();
                poiIsVist[poi_id]= true;
                //将在NVD graph中与当前poi相邻的adj_poi加入
                vector<int> poi_ADJ = getPOIAdj_NVD_By_Keyword(poi_id,term); //jins

#ifdef TRACK
                cout<<"getPOIAdj_NVD_By_Keyword(p"<<poi_id<<",t"<<term<<") size= "<<poi_ADJ.size()<<endl;
                if(26788==poi_id){
                    getchar();
                }
#endif
                for(int adj_poi_id : poi_ADJ){
                    ///cout<<"adj_poi_id: p"<<adj_poi_id<<endl;
                    if(poiIsVist.count(adj_poi_id)) continue;
                    if(adj_poi_id<0 ||adj_poi_id>poi_num) continue;
                    POI adj_poi = getPOIFromO2UOrgLeafData(adj_poi_id);
                    float adj_dist_upper = getUpperDistance_b2P_phl(border_id,adj_poi);
                    nodeDistType_NVD tmp1(-1, adj_poi_id, adj_dist_upper, term);
                    H_t.push(tmp1);
                }

            }

        }
        else{ //否，则为低频词汇

            posting_ListSize_Map[term] = posting_size;
            vector<int> posting_list = getTermOjectInvList(term);
            posting_Map[term] = posting_list;
            int _size = posting_list.size();


#ifdef TRACKBorderSB
            //cout<<"t"<<term<<"'s posting_list size="<<posting_list.size()<<endl;
#endif

            for(int p_id: posting_list){
                POI poi = getPOIFromO2UOrgLeafData(p_id);
                double simT_LB = 0; double simD_UB = 0; double simS_LB = 1;
                simT_LB = tfIDF_term(term);
                simD_UB = getUpperDistance_b2P_phl(border_id,poi);
                double social_textual_LB = alpha*simS_LB + (1-alpha)*simT_LB;
                double scoring_LB = social_textual_LB / (1+simD_UB);
                SBEntry sb_entry(p_id,simD_UB,scoring_LB);
                scoring_Bound_list_t.push(sb_entry);
            }
            while(scoring_Bound_list_t.size()>K)
                scoring_Bound_list_t.pop();


        }

        double sb_rk_t = scoring_Bound_list_t.top().score_bound;
        GlobalSBEntry termEntry(term,idx,sb_rk_t);
        rank_term.push(termEntry);
        //scoring_Bound_lists.push_back(scoring_Bound_list_t);
        idx++;


    }
    set<int> _Keys = Keys;

    //对按rank_term中关键词排序，根据scoring_Bound_lists中内容，对Key进行更新
#ifdef TRACKBorderSB
    cout<<"对按rank_term中关键词排序，根据scoring_Bound_lists中内容，对Key进行更新"<<endl;
#endif
    while(rank_term.size()>0){
        GlobalSBEntry _termEntry = rank_term.top();
        rank_term.pop();
        int term_id = _termEntry.term_id;
        double sb_rk_t = _termEntry.sk_bound;

        //计算SB_upper(e, qo)
        double score_upper = -1;
        set<int> cluster_Keys = pnodes_uCLUSTERListMap[poi_id][node_id][term_id];
        set<int> interKeys = obtain_itersection_jins(cluster_Keys,_Keys);
        //double max_relevance = tfIDF_termSet(_Keys);
        double max_relevance = tfIDF_termSet(interKeys);
        double min_distance = current_dist;

        double social_textual_maxscore = alpha*1 + (1.0-alpha)*max_relevance;
        score_upper = social_textual_maxscore / (1.0+ min_distance);
#ifdef TRACKBorderSB
        cout<<"for t"<<term_id<<", sb_rk_t="<<sb_rk_t<<",score_upper="<<score_upper<<endl;
#endif
        if(sb_rk_t > score_upper){ //remove term from _Keys
#ifdef TRACKBorderSB
            //cout<<"term"<<term_id<<"可被去掉"<<endl;
            //cout<<"此时_Keys:"; printSetElements(_Keys);
#endif
            set<int> :: iterator iter;
            for(iter= _Keys.begin();iter!=_Keys.end();){
                if(*iter == term_id){
                    _Keys.erase(iter);
                    break;
                }

                iter++;
            }
#ifdef TRACKBorderSB
            //cout<<"去掉后，_Keys:"; printSetElements(_Keys);
#endif
        }
        //cout<<", for term"<<term_id<<",sb_rk_t="<<sb_rk_t;
    }
#ifdef TRACKBorderSB
    cout<<endl;
#endif

    return _Keys;

}

set<int> getBorder_SB_BatchNVD_Cluster(BatchCheckEntry border_entry,int K, int a, double alpha, unordered_map<int,UCLUSTERListMap>& pnodes_uCLUSTERListMap){  //key:poi_id, value:UCLUSTERList

    int border_id = border_entry.id;
    int node_id = border_entry.node_id;
    int poi_id = border_entry.p_id;
#ifdef TRACKBorderSB
    cout<<"getBorder_SB_BatchNVD for border: v"<<border_id<<endl;
#endif
    set<int> Keys = border_entry.keys_cover;
#ifdef TRACKBorderSB
    cout<<"Keys:"; printSetElements(Keys);
#endif

    float current_dist = border_entry.dist;


    map<int, vector<int>> posting_Map;
    map<int, int> posting_ListSize_Map;

    map<int, int> posting_ListCurrentPos_Map;


    priority_queue<GlobalSBEntry>  rank_term;


    int idx = 0;

    //对各个term的 sb_list先进行计算；
    for(int term: Keys){  //对Keys中每个单词，进行评分下界队列（sb_list(term)）求解
#ifdef TRACKBorderSB
        //cout<<"对t"<<term<<"的sb_list进行计算"<<endl;
#endif

        map<int,bool> poiIsVist;
        priority_queue<nodeDistType_NVD> H_t;
        priority_queue<SBEntry> scoring_Bound_list_t;

        int posting_size = getTermOjectInvListSize(term);
        //如果term为高频词汇
        if(posting_size > posting_size_threshold){

            vector<POI> polled_object_lists;
            vector<float> polled_distance_lists;

#ifdef HybridHash   ////注意这里使用了混合哈希索引策略
            int nearest_poi_id = getNNPOI_By_HybridVertex_Keyword(border_id,term);  //先得到距离vertex最近的带关键词term的poi
            if(nearest_poi_id ==-1){
                vector<int> pois;
                int leaf_node = Nodes[border_id].gtreepath.back();
                pois = getNNPOIList_by_Hybrid_Hash(leaf_node,term);
                for(int poi_id2: pois){
                    POI poi = getPOIFromO2UOrgLeafData(poi_id2);
                    float dist_upper = getUpperDistance_b2P_phl(border_id,poi);
                    nodeDistType_NVD tmp1(-1, poi_id2, dist_upper, term);
                    H_t.push(tmp1);
#ifdef TRACK
                    cout<<"list加入o"<<poi_id2<<endl;
                    if(poi_id2==3144050){
                        cout<<"border_id="<<border_id<<"leaf="<<leaf_node<<"pois.size="<<pois.size()<<endl;
                        getchar();
                    }
#endif
                }
            }
            else{
                POI poi = getPOIFromO2UOrgLeafData(nearest_poi_id);
                float dist_upper = getUpperDistance_b2P_phl(border_id,poi);
                nodeDistType_NVD tmp1(-1, nearest_poi_id, dist_upper, term);
                H_t.push(tmp1);

            }


#else
            int nearest_poi_id = getNNPOI_By_Vertex_Keyword(border_id,term);  //先得到距离vertex最近的带关键词term的poi
            if(nearest_poi_id ==-100) continue;   //此时说明 所有term相关的poi 相对user都不可达
            POI poi = getPOIFromO2UOrgLeafData(nearest_poi_id);
            float dist_upper = getUpperDistance_b2P_phl(border_id,poi);
            nodeDistType_NVD tmp1(-1, nearest_poi_id, dist_upper, term);
            H_t.push(tmp1);
#endif


            while(scoring_Bound_list_t.size()<K){
                //H_t队首元素出列
                nodeDistType_NVD _tmp = H_t.top();
                H_t.pop();
                int poi_id = _tmp.poi_id;
                if(poiIsVist.count(poi_id)) continue;  //poiIsVist[poi_id]==true
                //计算scoring lower bound
                double simT_LB = 0; double simD_UB = 0; double simS_LB = 1;
                simT_LB = tfIDF_term(term);
                simD_UB = _tmp.dist;
                double social_textual_LB = alpha*simS_LB + (1-alpha)*simT_LB;
                double scoring_LB = social_textual_LB / (1+simD_UB);
                SBEntry sb_entry(poi_id,simD_UB,scoring_LB);
                scoring_Bound_list_t.push(sb_entry);
                //sb_entry.printRlt();
                poiIsVist[poi_id]= true;
                //将在NVD graph中与当前poi相邻的adj_poi加入
                vector<int> poi_ADJ = getPOIAdj_NVD_By_Keyword(poi_id,term);
                for(int adj_poi_id : poi_ADJ){
                    if(poiIsVist.count(adj_poi_id)) continue;
                    if(adj_poi_id<0 ||adj_poi_id>poi_num) continue;
                    POI adj_poi = getPOIFromO2UOrgLeafData(adj_poi_id);
                    float adj_dist_upper = getUpperDistance_b2P_phl(border_id,adj_poi);
                    nodeDistType_NVD tmp1(-1, adj_poi_id, adj_dist_upper, term);
                    H_t.push(tmp1);
                }
            }

        }
        else{ //否，则为低频词汇

            posting_ListSize_Map[term] = posting_size;
            vector<int> posting_list = getTermOjectInvList(term);
            posting_Map[term] = posting_list;
            int _size = posting_list.size();


#ifdef TRACKBorderSB
            //cout<<"t"<<term<<"'s posting_list size="<<posting_list.size()<<endl;
#endif

            for(int p_id: posting_list){
                POI poi = getPOIFromO2UOrgLeafData(p_id);
                double simT_LB = 0; double simD_UB = 0; double simS_LB = 1;
                simT_LB = tfIDF_term(term);
                simD_UB = getUpperDistance_b2P_phl(border_id,poi);
                double social_textual_LB = alpha*simS_LB + (1-alpha)*simT_LB;
                double scoring_LB = social_textual_LB / (1+simD_UB);
                SBEntry sb_entry(p_id,simD_UB,scoring_LB);
                scoring_Bound_list_t.push(sb_entry);
            }
            while(scoring_Bound_list_t.size()>K)
                scoring_Bound_list_t.pop();


        }

        double sb_rk_t = scoring_Bound_list_t.top().score_bound;
        GlobalSBEntry termEntry(term,idx,sb_rk_t);
        rank_term.push(termEntry);
        //scoring_Bound_lists.push_back(scoring_Bound_list_t);
        idx++;


    }
    set<int> _Keys = Keys;

    //对按rank_term中关键词排序，根据scoring_Bound_lists中内容，对Key进行更新
#ifdef TRACKBorderSB
    cout<<"对按rank_term中关键词排序，根据scoring_Bound_lists中内容，对Key进行更新"<<endl;
#endif
    while(rank_term.size()>0){
        GlobalSBEntry _termEntry = rank_term.top();
        rank_term.pop();
        int term_id = _termEntry.term_id;
        double sb_rk_t = _termEntry.sk_bound;

        //计算SB_upper(e, qo)
        double score_upper = -1;
        set<int> cluster_Keys = pnodes_uCLUSTERListMap[poi_id][node_id][term_id];
        set<int> interKeys = obtain_itersection_jins(cluster_Keys,_Keys);
        //double max_relevance = tfIDF_termSet(_Keys);
        double max_relevance = tfIDF_termSet(interKeys);
        double min_distance = current_dist;

        double social_textual_maxscore = alpha*1 + (1.0-alpha)*max_relevance;
        score_upper = social_textual_maxscore / (1.0+ min_distance);
#ifdef TRACKBorderSB
        cout<<"for t"<<term_id<<", sb_rk_t="<<sb_rk_t<<",score_upper="<<score_upper<<endl;
#endif
        if(sb_rk_t > score_upper){ //remove term from _Keys
#ifdef TRACKBorderSB
            //cout<<"term"<<term_id<<"可被去掉"<<endl;
            //cout<<"此时_Keys:"; printSetElements(_Keys);
#endif
            set<int> :: iterator iter;
            for(iter= _Keys.begin();iter!=_Keys.end();){
                if(*iter == term_id){
                    _Keys.erase(iter);
                    break;
                }

                iter++;
            }
#ifdef TRACKBorderSB
            //cout<<"去掉后，_Keys:"; printSetElements(_Keys);
#endif
        }
        //cout<<", for term"<<term_id<<",sb_rk_t="<<sb_rk_t;
    }
#ifdef TRACKBorderSB
    cout<<endl;
#endif

    return _Keys;

}



typedef unordered_map<int, vector<int>>  MaterilizedDistanceMap;


//自底向上，加local、global 减枝策略的filtering算法

//extend valuation range for users in the road network

void Extend_Range(int& node_highest,unordered_map<int, vector<int> >& itm, priority_queue<CheckEntry>& Queue,set<int>& nonRare_Keys, POI& poi, int K, float a, float alpha){
    //表示当前范围内已全部遍历， 需扩大遍历子图范围
    int posa, posb; float _min, dis;
    int root_id =0;
    int node_highest_pre = node_highest;
    int second_highest =0;
    vector<int> keywords;
    set<int> keyword_Set;
    for(int term: nonRare_Keys){
        keywords.push_back(term);
        keyword_Set.insert(term);
    }


    //寻找与qo有关键词交集且在Gtree中level最低的user 节点
    int current = node_highest; int current_pre = current;
    set<int> entry_set;
#ifdef TRACKEXTEND
    if(GTree[current_pre].isleaf)
        cout<<"从leaf"<<current_pre<<"向上寻找"<<endl;
    else
        cout<<"从node"<<current_pre<<"向上寻找"<<endl;
#endif

    while(current != root_id){
        // 获得该中间节点下包含用户关键词的所有child
        current_pre = current;
        current = GTree[current].father;
        upper_jumper++;
        //node_pos--;   //最高层节点层数往上一层

        //往上一层的border点的距离需要计算：
        for (int j = 0; j < GTree[current].borders.size(); j++) {  //从child node的每一个border 出发
            _min = -1;
            posa = GTree[current].current_pos[j];
            for (int k = 0; k < GTree[current_pre].borders.size(); k++) { //到parent node的每一个border
                posb = GTree[current_pre].up_pos[k];
                dis = itm[current_pre][k] + GTree[current].mind[posa * GTree[current].union_borders.size() + posb];
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
            itm[current].push_back(_min);      //自底向上记录 qo 到当前最大评估范围对应的node的各个边界点的距离
        }


        //判断节点上移后，其中用户是否与查询对象关键词相关
        for (int term: keywords) {
            vector<int> u_entrylist = getUsrTermRelatedEntry(term,current);//getObjectTermRelatedEntry(term,current);
            if (u_entrylist.size() > 0) {
                for (int child: u_entrylist) {
                    if(child ==current_pre)  //不能是同支路的孩子节点
                        continue;

                    entry_set.insert(child);
                    //block_num++;
                }
            }
        }
        //父节点就有关键词
        if(entry_set.size()>0){
            node_highest = current;
            second_highest = current_pre;
#ifdef TRACKEXTEND
            cout<<"到node"<<current<<",发现其有关键词！"<<endl;
#endif
            break;
        }
        else{

            //继续往上寻找关键词LCA
#ifdef TRACKEXTEND
            cout<<"到node"<<current<<endl;
#endif
        }

    }
    set<int> Key_upper;
#ifdef TRACKEXTEND
    //检查 拓展是否可以提前终止
    cout<<"当前最高节点为n"<<node_highest<<",对其进行evaluation!"<<endl;
#endif

    clock_t  evaluation_start, evaluation_end;
    evaluation_start = clock();

    int b_th = 0;
    for(int b_id: GTree[node_highest].borders){
#ifdef TRACKNVDBATCH
        cout<<"在n"<<node_highest<<"中检查b"<<b_id<<endl;
#endif
        float b_dist = getDistance_phl(b_id,poi);
        set<int> Key_b;
        //getchar();
        CheckEntry border_entry(b_id, node_highest, nonRare_Keys, b_dist,b_th);
        Key_b = getBorder_SB_NVD(border_entry,K,a,alpha);
#ifdef TRACKNVDBATCH
        cout<<"remained 关键词:"; printSetElements(Key_b);
#endif
        Key_upper.insert(Key_b.begin(),Key_b.end());
        b_th++;
    }
    evaluation_end = clock();
#ifdef TRACKEXTEND
    cout<<"最终，Key_upper："; printSetElements(Key_upper);
    cout<<"评估node"<<node_highest<<"用时："<<(double)(evaluation_end-evaluation_start)/CLOCKS_PER_SEC*1000<<" ms!"<<endl;

#endif
    if(Key_upper.size()==0&&node_highest!=0) {
#ifdef TRACKNVDBATCH
        cout<<"路网拓展可被终止!"<<endl;
#endif
        node_highest = -1;
    }


    //将扩大的查询范围子图加入
    evaluation_start = clock();

    for(int child: entry_set){  //每个child
        if(child == second_highest) {
            //cout<<"跳过"<<child<<endl;
            //getchar();
            continue; //跳过同支路节点
        }

        itm[child].clear();
        int allmin = -1;

        TreeNode child_node = getGIMTreeNodeData(child,OnlyU);

        set<int> inter_Key = obtain_itersection_jins(child_node.userUKeySet, keyword_Set);

        for (int j = 0; j < GTree[child].borders.size(); j++) {  //计算到child中的每个border的距离
            int border_id = GTree[child].borders[j];
            _min = -1;
#ifdef TRACK
            if(GTree[child].isleaf)
                cout<<"考虑加入leaf"<<child<<"的border"<<endl;
            else
                cout<<"考虑加入node"<<child<<"的border"<<endl;
#endif
            posa = GTree[child].up_pos[j];
            for (int k = 0; k < GTree[second_highest].borders.size(); k++) {  //重要！ 容易出错！
                posb = GTree[second_highest].current_pos[k];
                //int right_end = second_highest;
                dis = itm[second_highest][k] + GTree[current].mind[posb * GTree[current].union_borders.size() + posa];
                //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                if (_min == -1) {
                    _min = dis;
                } else {
                    if (dis < _min) {
                        _min = dis;
                    }
                }
            }
            itm[child].push_back(_min);
            float dist2 = _min; int border_th = j;
            float dist = getDistance_phl(border_id,poi);
            CheckEntry entry(border_id, child, inter_Key, dist,j);  //int check_id, int n_id, vector<int> _keys, float d
            Queue.push(entry);
#ifdef TRACK
            cout<<"在Extand中插入border entry:"; entry.printRlt();
#endif


        }


    }
    evaluation_end = clock();
#ifdef TRACKEXTEND
    cout<<"加入node"<<node_highest<<"的后继完毕，用时："<<(double)(evaluation_end-evaluation_start)/CLOCKS_PER_SEC*1000<<" ms!"<<endl;
#endif



}


void BatchExtend_Range(priority_queue<POIHighestNode>& T_max, unordered_map<int, MaterilizedDistanceMap>& disMapList, priority_queue<BatchCheckEntry>& Queue,unordered_map<int,set<int>>& nonRare_KeyMap, unordered_map<int, POI>& poisMap, int K, float a, float alpha){
    //表示当前范围内已全部遍历， 需扩大遍历子图范围

    priority_queue<POIHighestNode> T_max_new;
    while(T_max.size()>0){

        POIHighestNode pair = T_max.top();
        T_max.pop();
        int node_highest = pair.node_highest;
        if(node_highest==0)
            continue;

        int poi_id = pair.p_id;
        POI poi = poisMap[poi_id];

        int posa, posb; float _min, dis;
        int root_id =0;
        int node_highest_pre = node_highest;
        int second_highest =0;
        vector<int> keywords;
        set<int> keyword_Set;
        for(int term: nonRare_KeyMap[poi_id]){
            keywords.push_back(term);
            keyword_Set.insert(term);
        }


        //寻找与qo有关键词交集且在Gtree中level最低的user 节点
        int current = node_highest; int current_pre = current;
        set<int> entry_set;
#ifdef TRACKEXTEND
        if(GTree[current_pre].isleaf)
            cout<<"从leaf"<<current_pre<<"向上寻找"<<endl;
        else
            cout<<"从node"<<current_pre<<"向上寻找"<<endl;

        if(current_pre==1){
            cout<<"find node 1"<<endl;
        }
#endif

        while(current != root_id){
            // 获得该中间节点下包含用户关键词的所有child
            current_pre = current;
            current = GTree[current].father;
            upper_jumper++;
            //node_pos--;   //最高层节点层数往上一层

            //往上一层的border点的距离需要计算：
            for (int j = 0; j < GTree[current].borders.size(); j++) {  //从child node的每一个border 出发
                _min = -1;
                posa = GTree[current].current_pos[j];
                for (int k = 0; k < GTree[current_pre].borders.size(); k++) { //到parent node的每一个border
                    posb = GTree[current_pre].up_pos[k];
                    dis = disMapList[poi_id][current_pre][k] + GTree[current].mind[posa * GTree[current].union_borders.size() + posb];
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
                disMapList[poi_id][current].push_back(_min);      //自底向上记录 qo 到当前最大评估范围对应的node的各个边界点的距离
            }


            //判断节点上移后，其中用户是否与查询对象关键词相关
            for (int term: keywords) {
                vector<int> u_entrylist = getUsrTermRelatedEntry(term,current);//getObjectTermRelatedEntry(term,current);
                if (u_entrylist.size() > 0) {
                    for (int child: u_entrylist) {
                        if(child ==current_pre)  //不能是同支路的孩子节点
                            continue;

                        entry_set.insert(child);
                        //block_num++;
                    }
                }
            }
            //父节点就有关键词
            if(entry_set.size()>0){
                node_highest = current;
                second_highest = current_pre;
#ifdef TRACKEXTEND
                cout<<"到node"<<current<<",发现其有关键词！"<<endl;
#endif
                break;
            }
            else{

                //继续往上寻找关键词LCA
#ifdef TRACKEXTEND
                cout<<"到node"<<current<<endl;
#endif
            }

        }
        set<int> Key_upper;

        //检查 拓展是否可以提前终止
#ifdef TRACKEXTEND
        cout<<"当前最高节点为n"<<node_highest<<",其内部："<<endl;
#endif

        int b_th = 0;
        for(int b_id: GTree[node_highest].borders){
#ifdef TRACKEXTEND
            //cout<<"在n"<<node_highest<<"中检查b"<<b_id<<endl;
#endif
            float b_dist = getDistance_phl(b_id,poi);
            set<int> Key_b;
            //getchar();
            BatchCheckEntry border_entry(b_id, poi_id, node_highest, nonRare_KeyMap[poi_id], b_dist,b_th);
            Key_b = getBorder_SB_BatchNVD(border_entry,K,a,alpha);
#ifdef TRACKEXTEND
            //cout<<"remained 关键词:"; printSetElements(Key_b);
#endif
            Key_upper.insert(Key_b.begin(),Key_b.end());
            b_th++;
        }
#ifdef TRACKEXTEND
        cout<<"最终，Key_upper："; printSetElements(Key_upper);
#endif
        if(Key_upper.size()==0 && node_highest!=0) {
#ifdef TRACKEXTEND
            cout<<"p"<<poi_id<<"处的路网拓展可被终止!"<<endl;
#endif
            //continue;  //路网拓展可被终止
        }
        else{
            POIHighestNode pair_update(poi_id,node_highest);

            T_max_new.push(pair_update);
        }


        //将扩大的查询范围子图加入
        for(int child: entry_set){  //每个child
            if(child == second_highest) {

                continue; //跳过同支路节点
            }

            disMapList[poi_id][child].clear();
            int allmin = -1;

            TreeNode child_node = getGIMTreeNodeData(child,OnlyU);

            set<int> inter_Key = obtain_itersection_jins(child_node.userUKeySet, keyword_Set);
#ifdef TRACKEXTEND
            if(GTree[child].isleaf)
                cout<<"考虑加入leaf"<<child<<"的border"<<endl;
            else
                cout<<"考虑加入node"<<child<<"的border"<<endl;
#endif
            for (int j = 0; j < GTree[child].borders.size(); j++) {  //计算到child中的每个border的距离
                int border_id = GTree[child].borders[j];
                _min = -1;
                posa = GTree[child].up_pos[j];
                for (int k = 0; k < GTree[second_highest].borders.size(); k++) {  //重要！ 容易出错！
                    posb = GTree[second_highest].current_pos[k];
                    //int right_end = second_highest;
                    dis = disMapList[poi_id][second_highest][k] + GTree[current].mind[posb * GTree[current].union_borders.size() + posa];
                    //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                    if (_min == -1) {
                        _min = dis;
                    } else {
                        if (dis < _min) {
                            _min = dis;
                        }
                    }
                }
                disMapList[poi_id][child].push_back(_min);
                float dist2 = _min; int border_th = j;
                float dist = getDistance_phl(border_id,poi);
                BatchCheckEntry entry(border_id, poi_id, child, inter_Key, dist,j);  //int check_id, int n_id, vector<int> _keys, float d
                Queue.push(entry);
#ifdef TRACKEXTEND2
                cout<<"在BatchExtand中插入border entry:"; entry.printRlt();
                cout<<"Queue size="<<Queue.size()<<endl;
                if(border_id==22297&& 352==child){
                    cout<<"find border_id==22297&& 352==child"<<endl;
                    cout<<"Queue size="<<Queue.size()<<endl;

                }
#endif




            }


        }



    }//end while

    T_max = T_max_new;





}

//更为细分
void BatchExtend_Range(priority_queue<POIHighestNode>& T_max, unordered_map<int, MaterilizedDistanceMap>& disMapList, priority_queue<BatchCheckEntry>& Queue,unordered_map<int,set<int>>& nonRare_KeyMap, unordered_map<int, POI>& poisMap, unordered_map<int,UCLUSTERList>& pnodes_uCLUSTERListMap, int K, float a, float alpha){
    //表示当前范围内已全部遍历， 需扩大遍历子图范围

    priority_queue<POIHighestNode> T_max_new;
    while(T_max.size()>0){

        POIHighestNode pair = T_max.top();
        T_max.pop();
        int node_highest = pair.node_highest;
        if(node_highest==0)
            continue;

        int poi_id = pair.p_id;
        POI poi = poisMap[poi_id];

        int posa, posb; float _min, dis;
        int root_id =0;
        int node_highest_pre = node_highest;
        int second_highest =0;
        vector<int> keywords;
        set<int> keyword_Set;
        for(int term: nonRare_KeyMap[poi_id]){
            keywords.push_back(term);
            keyword_Set.insert(term);
        }


        //寻找与qo有关键词交集且在Gtree中level最低的user 节点
        int current = node_highest; int current_pre = current;
        set<int> entry_set;
#ifdef TRACKEXTEND
        if(GTree[current_pre].isleaf)
            cout<<"从leaf"<<current_pre<<"向上寻找"<<endl;
        else
            cout<<"从node"<<current_pre<<"向上寻找"<<endl;

        if(current_pre==1){
            cout<<"find node 1"<<endl;
        }
#endif

        while(current != root_id){
            // 获得该中间节点下包含用户关键词的所有child
            current_pre = current;
            current = GTree[current].father;
            upper_jumper++;
            //node_pos--;   //最高层节点层数往上一层

            //往上一层的border点的距离需要计算：
            for (int j = 0; j < GTree[current].borders.size(); j++) {  //从child node的每一个border 出发
                _min = -1;
                posa = GTree[current].current_pos[j];
                for (int k = 0; k < GTree[current_pre].borders.size(); k++) { //到parent node的每一个border
                    posb = GTree[current_pre].up_pos[k];
                    dis = disMapList[poi_id][current_pre][k] + GTree[current].mind[posa * GTree[current].union_borders.size() + posb];
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
                disMapList[poi_id][current].push_back(_min);      //自底向上记录 qo 到当前最大评估范围对应的node的各个边界点的距离
            }


            //判断节点上移后，其中用户是否与查询对象关键词相关
            for (int term: keywords) {
                vector<int> u_entrylist = getUsrTermRelatedEntry(term,current);//getObjectTermRelatedEntry(term,current);
                if (u_entrylist.size() > 0) {
                    for (int child: u_entrylist) {
                        if(child ==current_pre)  //不能是同支路的孩子节点
                            continue;

                        entry_set.insert(child);
                        //block_num++;
                    }
                }
            }
            //父节点就有关键词
            if(entry_set.size()>0){
                node_highest = current;
                second_highest = current_pre;
#ifdef TRACKEXTEND
                cout<<"到node"<<current<<",发现其有关键词！"<<endl;
#endif
                break;
            }
            else{

                //继续往上寻找关键词LCA
#ifdef TRACKEXTEND
                cout<<"到node"<<current<<endl;
#endif
            }

        }
        set<int> Key_upper;

        //检查 拓展是否可以提前终止
#ifdef TRACKEXTEND
        cout<<"当前最高节点为n"<<node_highest<<",其内部："<<endl;
#endif

        int b_th = 0;
        for(int b_id: GTree[node_highest].borders){
#ifdef TRACKEXTEND
            //cout<<"在n"<<node_highest<<"中检查b"<<b_id<<endl;
#endif
            float b_dist = getDistance_phl(b_id,poi);
            set<int> Key_b;
            //getchar();
            BatchCheckEntry border_entry(b_id, poi_id, node_highest, nonRare_KeyMap[poi_id], b_dist,b_th);
            //Key_b = getBorder_SB_BatchNVD(border_entry,K,a,alpha);
            Key_b = getBorder_SB_BatchNVD_Cluster(border_entry,K,a,alpha,pnodes_uCLUSTERListMap);
#ifdef TRACKEXTEND
            //cout<<"remained 关键词:"; printSetElements(Key_b);
#endif
            Key_upper.insert(Key_b.begin(),Key_b.end());
            b_th++;
        }
#ifdef TRACKEXTEND
        cout<<"最终，Key_upper："; printSetElements(Key_upper);
#endif
        if(Key_upper.size()==0 && node_highest!=0) {
#ifdef TRACKEXTEND
            cout<<"p"<<poi_id<<"处的路网拓展可被终止!"<<endl;
#endif
            //continue;  //路网拓展可被终止
        }
        else{
            POIHighestNode pair_update(poi_id,node_highest);

            T_max_new.push(pair_update);
        }


        //将扩大的查询范围子图加入
        for(int child: entry_set){  //每个child
            if(child == second_highest) {

                continue; //跳过同支路节点
            }

            disMapList[poi_id][child].clear();
            int allmin = -1;

            TreeNode child_node = getGIMTreeNodeData(child,OnlyU);

            set<int> inter_Key = obtain_itersection_jins(child_node.userUKeySet, keyword_Set);
#ifdef TRACKEXTEND
            if(GTree[child].isleaf)
                cout<<"考虑加入leaf"<<child<<"的border"<<endl;
            else
                cout<<"考虑加入node"<<child<<"的border"<<endl;
#endif
            for (int j = 0; j < GTree[child].borders.size(); j++) {  //计算到child中的每个border的距离
                int border_id = GTree[child].borders[j];
                _min = -1;
                posa = GTree[child].up_pos[j];
                for (int k = 0; k < GTree[second_highest].borders.size(); k++) {  //重要！ 容易出错！
                    posb = GTree[second_highest].current_pos[k];
                    //int right_end = second_highest;
                    dis = disMapList[poi_id][second_highest][k] + GTree[current].mind[posb * GTree[current].union_borders.size() + posa];
                    //GTree[top.id].mind[posa * GTree[top.id].union_borders.size() + posb];
                    if (_min == -1) {
                        _min = dis;
                    } else {
                        if (dis < _min) {
                            _min = dis;
                        }
                    }
                }
                disMapList[poi_id][child].push_back(_min);
                float dist2 = _min; int border_th = j;
                float dist = getDistance_phl(border_id,poi);
                BatchCheckEntry entry(border_id, poi_id, child, inter_Key, dist,j);  //int check_id, int n_id, vector<int> _keys, float d
                Queue.push(entry);
#ifdef TRACKEXTEND2
                cout<<"在BatchExtand中插入border entry:"; entry.printRlt();
                cout<<"Queue size="<<Queue.size()<<endl;
                if(border_id==22297&& 352==child){
                    cout<<"find border_id==22297&& 352==child"<<endl;
                    cout<<"Queue size="<<Queue.size()<<endl;

                }
#endif




            }


        }



    }//end while

    T_max = T_max_new;





}





/*-------------------------------全数据集top-k, 反topk结果输出-----------------------------*/


//将"../data/LV/LV.query"中每个user的top-k评分结果输出到../data/LV/LV.userScore文件中
int outputUserScore(int Qk, float a, float alpha) {
    vector<query> querySet;
    loadQuery(querySet);  //从/data/LV/LV.query"中读取每个user的信息
    ofstream outputScore;
    stringstream tt;
    tt<<FILE_USERSCORE<<"(Qk="<<Qk<<",a="<<a<<")";
    outputScore.open(tt.str());
    for (size_t i = 0; i < Users.size(); i++) {
        double topk = topkSDijkstra(transformToQ(Users[i]),Qk,a,alpha).topkScore;
        outputScore << i << ' ' << topk << endl;
    }
    outputScore.close();
    return 0;
}

//计算并输出兴趣点的反top-k结果
int outputPOIRkNNResults(int Qk, float a, float alpha) {
    cout<<"outputPOIRkNNResults(pre-procession work)...(Qk="<<Qk<<",a="<<a<<")"<<endl;
    //vector<query> querySet;
    //loadQuery(querySet);  //从/data/LV/LV.query"中读取每个user的信息
    ofstream outputScore;
    ofstream outputTopKResults;
    stringstream tt;
    stringstream tt2;
    tt<<FILE_USERSCORE<<"(Qk="<<Qk<<",a="<<a<<")";
    tt2<<FILE_USERTopKResults<<"(Qk="<<Qk<<",a="<<a<<")";

    outputScore.open(tt.str());
    outputTopKResults.open(tt2.str());

    //poi_ID, vector<u_id>, maintain a map for each poi and its reverse top-k results set
    map<int, vector<int>> reverseResults;
    for (User usr: Users) {
        //query Query = querySet[234];
        //printQuery(Query);
        query usr_q = transformToQ(usr);
        TopkQueryResult qr = topkSDijkstra(usr_q,Qk,a,alpha);
        //vector<Result> topkElements = qr.topkElements;
        float topkScore = qr.topkScore;
        //输出topk评分信息
        outputScore << usr.id  << ' ' << topkScore << endl;
        //cout << i << ' ' << topkScore << endl;
        //输出topk对象信息
        outputTopKResults<< "user"<< usr.id << "'s top-k elements:"<<endl;
        for(Result r : qr.topkElements){
            int o_id = r.id;
            float o_score = r.score;
            outputTopKResults<<"<p"<<o_id<<",score="<< o_score<<">,";
            /*if(reverseResults.count(r.id)==0){
                vector<int> reverse_usr;
                reverse_usr.push_back(i);
                reverseResults[r.id] = reverse_usr;
            }
            else
                reverseResults[r.id].push_back(i);*/   //add usr i into p[r.id]'s interested set list
            //topkElements.pop();
            reverseResults[o_id].push_back(usr.id);
        }

        outputTopKResults<<endl;


    }

    outputScore.close();
    outputTopKResults.close();

    //反top-k结果
    stringstream tt3;
    tt3<<FILE_POIReverseResults<<"(Qk="<<Qk<<",a="<<a<<")";
    ofstream outputResults;
    outputResults.open(tt3.str());

    stringstream tt4;
    tt4<<FILE_LargePOIReverse<<"(Qk="<<Qk<<",a="<<a<<")";
    ofstream outputLogs;
    outputLogs.open(tt4.str());

    map<int, vector<int>>:: iterator iter;
    int optimal =0; vector<int> optimalUserSet; int poi_optimal;

    //输出反top-k结果
    for(iter = reverseResults.begin(); iter!= reverseResults.end(); iter++){
        int p_id = iter ->first;
        vector<int> users = iter ->second;
        //cout<<"poi_"<<p_id<<"'s potential user :"<<users.size()<<endl<<"<";
        outputResults<<"poi_"<<p_id<<" has " << users.size() << "potential users:"<<endl<<"<";
        cout<<"<poi_"<<p_id<<" has" << users.size() << " potential users: ";
        for(int usr: users) cout<<"u"<<usr<<",";
        cout<<">"<<endl;
        //记录最多的反-top-k结果

        if(users.size()>20){
            cout<<"poi"<<p_id<<"has interseted user set of:"<< users.size()<<endl<<"<";
            outputLogs<<"poi"<<p_id<<"has interseted user set of:"<< users.size()<<endl;
            for(int i:users){
                //cout<<"u"<<i<<",";
                cout<<"u"<<i<<",";
            }
            cout<<">"<<endl;

        }
        if(users.size()>optimal){
            optimal = users.size();
            optimalUserSet = users;
            poi_optimal= p_id;
        }
        for(int i:users){
            //cout<<"u"<<i<<",";
            outputResults<<"u"<<i<<",";
        }

        //cout<<">"<<endl;
        outputResults<<">"<<endl;
        //getchar();
    }
    outputResults.close();

    cout<<"outputPOIRkNNResults COMPLETE !"<<endl;
    cout<<"poi"<<poi_optimal<<"has the largest interseted user set:"<< optimal<<endl<<"<";
    outputLogs<<"poi"<<poi_optimal<<"has the largest interseted user set:"<< optimal<<endl<<"<";
    outputLogs.close();
    for(int i:optimalUserSet){
        //cout<<"u"<<i<<",";
        cout<<"u"<<i<<",";
    }
    cout<<">"<<endl;
    //getchar();
    return 0;
}


//分析RkGSKQ 结果
SingleResults  resultsAnalysis_loc(vector<ResultDetail> candidates, vector<int> keywords, vector<int> checkin_usrs, int loc, int Qk, float a, float alpha){
    vector<ResultDetail> results;
    int verification = 0;
    map <int,vector<int>> invCandidateKeys;
    map <int,vector<int>> invCandidateSeeds;
    for(size_t j=0; j< candidates.size();j++){
        int usr_id = candidates[j].usr_id;
        //对每一个候选用户进行评估
        User usr = Users[usr_id];

        double score = candidates[j].score;

        query usr_q = transformToQ(usr);
        double Rk_usr = candidates[j].rk;
        //cout<<topkSDijkstra(usr_q,Qk,a, alpha).topkScore<<endl;
        verification ++;
        //Users[usr_id].isVisited = false; //表明该节点已验证过
        //cout<<"u"<<usr_id<<",score="<<score<<", Rk_usr="<<Rk_usr<<endl;
        //if (score > Rk_usr) {
        for(int t:usr.keywords)
            invCandidateKeys[t].push_back(usr.id);
        for(int u:friendshipMap[usr.id])
            invCandidateSeeds[u].push_back(usr.id);

        //sd.printResultDetail();
        //getchar();
        results.push_back(candidates[j]);
        //}
    }

    SingleResults filter(invCandidateKeys,invCandidateSeeds,results);
    return  filter;
}

SingleResults  resultsAnalysis_poi(vector<ResultDetail> candidates, POI poi, int Qk, float a, float alpha){
    vector<ResultDetail> results;
    int verification = 0;
    map <int,vector<int>> invCandidateKeys;
    map <int,vector<int>> invCandidateSeeds;
    for(size_t j=0; j< candidates.size();j++){
        int usr_id = candidates[j].usr_id;
        //对每一个候选用户进行评估
        User usr = Users[usr_id];

        double score = candidates[j].score;

        query usr_q = transformToQ(usr);
        double Rk_usr = candidates[j].rk;
        //cout<<topkSDijkstra(usr_q,Qk,a, alpha).topkScore<<endl;
        verification ++;
        //Users[usr_id].isVisited = false; //表明该节点已验证过
        //cout<<"u"<<usr_id<<",score="<<score<<", Rk_usr="<<Rk_usr<<endl;
        //if (score > Rk_usr) {
        for(int t:usr.keywords)
            invCandidateKeys[t].push_back(usr.id);
        for(int u:friendshipMap[usr.id])
            invCandidateSeeds[u].push_back(usr.id);

        //sd.printResultDetail();
        //getchar();
        results.push_back(candidates[j]);
        //}
    }

    SingleResults sr(invCandidateKeys,invCandidateSeeds,results);
    return  sr;
}


//用nvd进行topk求解
TopkQueryCurrentResult TkGSKQ_NVD(User user, int K, int a, double alpha){
    clock_t startTime, endTime, social_endTime;
    clock_t io_start, io_end; double io_time=0;
    clock_t textual_start, textual_end; double textual_time=0;

    priority_queue<GlobalEntry_Pseudo> Queue;
    priority_queue<Result> resultFinal;  //存放topk结果
    float score_rk = -1;
    vector<POI> polled_object_lists;
    vector<float> polled_distance_lists;

    int vi = user.Ni; int vj = user.Nj;
    //map<int,bool> poiIsVist;

    startTime = clock();//计时开始
    bool* poiMark = new bool[poi_num];
    bool* poiADJMark = new bool[poi_num];
    for(int i=0;i<poi_num;i++){
        poiMark[i]=false;
        poiADJMark[i] = false;
    }



    //先获得user的朋友
    int user_id = user.id;
    vector<int> friends = friendshipMap[user_id];
    //获得各个friends 留下check-in过的poi
    set<int> userRelated_poiSet;
    //vector<int> social_textualRelated_poiSet;

    int count = 0; int io_count = 0;
    //POI _poi = getPOIFromO2UOrgLeafData(1);
    if(true){
        for(int friend_id: friends){
            vector<int> checkin_objects = userCheckInIDList[friend_id];    //<poi_id...>
            for(int poi_id : checkin_objects){
                count++;

                if(poiMark[poi_id] == true) continue;
#ifdef TRACKNVDTOPK
                cout<<"social find related object: p"<<poi_id<<endl;
#endif

                POI poi = getPOIFromO2UOrgLeafData(poi_id);  //_poi ; //

                poiMark[poi_id] = true;
                //poiIsVist[poi_id] = true;
                //POIs[poi_id].isVisited = true;

                double relevance = 0;
                textual_start = clock();
                relevance = textRelevance(user.keywords, poi.keywords);
                //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
                if(relevance==0) continue;
                double distance =  usrToPOIDistance_phl(user,poi);//min(d1,d2);
                double social_Score = 0;
                social_Score =  reCalcuSocialScore(user.id,poiCheckInIDList[poi_id]);  //gazi
                double influence;
                if(social_Score == 0)
                    influence =1;
                else
                    influence =1 + pow(a,social_Score);
                double social_textual_score =  alpha*influence + (1-alpha)*relevance;
                double gsk_score = social_textual_score / (1+distance);

                userRelated_poiSet.insert(poi_id);
                //social_textualRelated_poiSet.push_back(poi_id);
                Result tmpRlt(poi_id,poi.Ni,poi.Nj,poi.dis,distance,gsk_score,poi.keywords);



                if(resultFinal.size()<K){
                    resultFinal.push(tmpRlt);
                }
                else{
                    Result _tmpNode = resultFinal.top();
                    double tmpScore = _tmpNode.score;
                    if(tmpScore<gsk_score){
                        resultFinal.pop();
                        resultFinal.push(tmpRlt);
                    }

                }

            }
        }
    }

    if(resultFinal.size()==K){
        score_rk = resultFinal.top().score;
    }
    social_endTime = clock();



    vector<int> qKey = user.keywords;
    set<int> qKeySet;
    for(int term:qKey)
        qKeySet.insert(term);
    double  optimal_text_score = textRelevance(qKey, qKey);
    //最优score
    double optimal_social_textual_score = alpha*(1.0)+(1-alpha)*optimal_text_score;


    ////基于K-SPIN策略进行topk结果查询
    //1.对每个user.keyword分配一个queue
    vector<priority_queue<nodeDistType_POI>>  H_list;
    int term_size = getTermSize();
    int termid2idx[term_size];
    bool termIsFrequent[term_size];
    memset(termIsFrequent, false, sizeof(termIsFrequent));
    for(int term_th = 0; term_th< user.keywords.size();term_th++){

#ifdef TRACKNVDTOPK
        cout<<"第"<<term_th<<"个单词！"<<endl;

#endif
        int term = user.keywords[term_th];

#ifdef TRACKNVDTOPK
        cout<<"该单词为t"<<term<<endl;

#endif
        termid2idx[term] = term_th;  //
        priority_queue<nodeDistType_POI> Hi;
        int posting_size = getTermOjectInvListSize(term);

#ifdef TRACKNVDTOPK
        cout<<"单词t"<<term<<"的posting size="<<posting_size<<endl;

#endif

#ifdef TRACKNVDTOPK
        cout<<"初始化填充H(t"<<term<<")的内容"<<endl;
        if(term==65){
            cout<<"find t65"<<endl;
        }
#endif

        //如果term为低频词汇
        if(posting_size<=posting_size_threshold){
            //termIsFrequent[term] = false;
            vector<int> posting_list = getTermOjectInvList(term);
            for(int o_id: posting_list){

                POI related_poi = getPOIFromO2UOrgLeafData(o_id);
                float dist = usrToPOIDistance_phl(user,related_poi);
                nodeDistType_POI tmp1(o_id, dist, related_poi.keywords);
                Hi.push(tmp1);
                poiADJMark[o_id] = true;
            }
        }
        else{ //高频词汇,则将最近邻POI加入
            termIsFrequent[term] = true;
            /*int poi_id = getNNPOI_By_Vertex_Keyword(user.Ni,term);  //先得到距离Ni最近的带关键词term的poi
            int poi_id2 = getNNPOI_By_Vertex_Keyword(user.Nj,term);*/

            int poi_id1 = getNNPOI_By_HybridVertex_Keyword(user.Ni,term);
            int poi_id2 = getNNPOI_By_HybridVertex_Keyword(user.Nj,term);

            if(poi_id1==-1){
                vector<int> pois_list1;
                int leaf1 = Nodes[user.Ni].gtreepath.back();
                pois_list1 = getNNPOIList_by_Hybrid_Hash(leaf1,term);
                for(int poi_within: pois_list1){
#ifdef TRACKNVDTOPK
                    cout<<"1.获得:o"<<poi_within<<endl;
#endif
                    POI poi = getPOIFromO2UOrgLeafData(poi_within);
                    float dist = usrToPOIDistance_phl(user,poi);
                    nodeDistType_POI tmp1(poi_within, dist, poi.keywords);
                    Hi.push(tmp1);
#ifdef TRACKNVDTOPK
                    cout<<"1.list中填充:o"<<poi_within<<endl;
#endif
                    poiADJMark[poi_within] = true;
#ifdef TRACKNVDTOPK
                    cout<<"1. 标记 list1中填充对象: o"<<poi_within<<endl;
#endif
                }
            }
            else{
#ifdef TRACKNVDTOPK
                cout<<"1.获得当个对象:o"<<poi_id1<<endl;
#endif
                POI poi = getPOIFromO2UOrgLeafData(poi_id1);
                float dist = usrToPOIDistance_phl(user,poi);
                nodeDistType_POI tmp1(poi_id1, dist, poi.keywords);
                Hi.push(tmp1);
#ifdef TRACKNVDTOPK
                cout<<"1填充单个对象：o"<<poi_id1<<endl;
#endif
                poiADJMark[poi_id1] = true;
#ifdef TRACKNVDTOPK
                cout<<"1标记填充对象：o"<<poi_id1<<endl;
#endif
            }

#ifdef TRACKNVDTOPK

            cout<<"将进行poi_id2的判断"<<endl;
#endif

            if(poi_id2==-1){
                vector<int> pois_list2;
                int leaf2 = Nodes[user.Nj].gtreepath.back();
#ifdef TRACKNVDTOPK
                cout<<"leaf2="<<leaf2<<endl;
#endif
                pois_list2 = getNNPOIList_by_Hybrid_Hash(leaf2,term);
#ifdef TRACKNVDTOPK
                cout<<"pois_list2 size="<<pois_list2.size()<<endl;
#endif
                for(int poi_within: pois_list2){
#ifdef TRACKNVDTOPK
                    cout<<"2.获得:o"<<poi_within<<endl;
#endif
                    if(poiADJMark[poi_within]==true) continue;
                    POI poi = getPOIFromO2UOrgLeafData(poi_within);
                    float dist = usrToPOIDistance_phl(user,poi);
                    nodeDistType_POI tmp1(poi_within, dist, poi.keywords);
                    Hi.push(tmp1);
#ifdef TRACKNVDTOPK
                    cout<<"2.list中填充对象: o"<<poi_within<<endl;
#endif
                    poiADJMark[poi_within] = true;
#ifdef TRACKNVDTOPK
                    cout<<"2.标记list中填充对象: o"<<poi_within<<endl;
#endif
                }
            }
            else{
                if(poiADJMark[poi_id2]!=true){
                    POI poi = getPOIFromO2UOrgLeafData(poi_id2);
                    float dist = usrToPOIDistance_phl(user,poi);
                    nodeDistType_POI tmp1(poi_id2, dist, poi.keywords);
                    Hi.push(tmp1);
#ifdef TRACKNVDTOPK
                    cout<<"2.填充：o"<<poi_id2<<endl;
#endif
                    poiADJMark[poi_id2] = true;
#ifdef TRACKNVDTOPK
                    cout<<"2.标记填充：o"<<poi_id2<<endl;
#endif
                }

            }

        }
        H_list.push_back(Hi);
        //priority_queue<nodeDistType_POI> _Hi = H_list[term_th];
        //int _size = _Hi.size();
    }

    //1. 初始化全局Queue，对各个Hi的优先级进行计算并排序
    for(int term_th =0; term_th< user.keywords.size(); term_th++){
        int term_i = user.keywords[term_th];
        //priority_queue<nodeDistType_POI> _Hi = H_list[term_th];
        //int _size = _Hi.size();
        float dist_i = H_list[term_th].top().dist;
        float textual_score = 0;
        ////cout<<"----t"<<term_i<<",dist_i="<<dist_i<<endl;

        float score = optimal_social_textual_score / (1.0+dist_i);//textual_score / dist_i;
        GlobalEntry_Pseudo entry_Q(term_i,dist_i,score);
        Queue.push(entry_Q);

    }
    float current_score_global = 10000;

    //2. 每次从全局Queue中选取当前最优的Hi中的元素
    while(Queue.size()>0 && Queue.top().score>score_rk){
        GlobalEntry_Pseudo entry = Queue.top();
        Queue.pop(); //global 出列
        int term_n = entry.term_id;
        int term_th = termid2idx[term_n];
        //cout<<"当前最佳Hi为t"<<term_n<<", score="<<entry.score<<endl;
        ////获得当前最佳Hi
        //priority_queue<nodeDistType_POI> Hn = H_list[term_th];
        ////Hi中队首元素出列，分析队首元素（1.将队首元素在NVD 中的邻居加入Hn；2.更新Hn当前的最佳评分；将Hn重新插回Queue）
        nodeDistType_POI entry_p  = H_list[term_th].top();
        int poi_current = entry_p.poi_id;

        vector<int> poi_current_keywords = entry_p.keywords;
        float dist_top_Hn = entry_p.dist;
        float score_Hn = entry.score;
        if(H_list[term_th].size()==0) continue;
        //cout<<"H"<<term_n<<"中队首元素为p"<<entry_p.poi_id<<"出列"<<endl;
        //Hn.pop();
        H_list[term_th].pop();

        if(termIsFrequent[term_n]){  //若为高频词汇
            ////1.将首元素（poi_current）的NVD_adj 加入Hn
            //cout<<"1.将首元素p"<<entry_p.poi_id<<"的NVD_adj 加入Hn"<<endl;
            //cout<<"尝试向H"<<term_n<<"加入新元素"<<endl;

#ifdef TRACKNVDTOPK
            if(entry_p.poi_id==2633){
                cout<<"while loop find o2633"<<endl;
            }

#endif
            vector<int> poi_ADJ = getPOIAdj_NVD_By_Keyword(entry_p.poi_id,term_n);
            //printElements(poi_ADJ);
            for(int _poi_id : poi_ADJ){
                //cout<<",nvd neighbor-p"<<_poi_id;
                if(poiADJMark[_poi_id] == true) continue;  //重要(Hn中不要重复加入)
                POI poi = getPOIFromO2UOrgLeafData(_poi_id); int id = poi.id;

                //float p_dist2 = usrToPOIDistance_phl(getUserFromO2UOrgLeafData(67),getPOIFromO2UOrgLeafData(3413));

                float p_dist = usrToPOIDistance_phl(user,poi);
                nodeDistType_POI tmp_adj(id, p_dist, poi.keywords);
                H_list[term_th].push(tmp_adj);
                poiADJMark[_poi_id] = true;
                //cout<<",确定加入(nvd_adj)p"<<id;

            }
        }

        //cout<<endl;
        ////2.更新Hn队列当前评分，重新压入全局队列
        //cout<<"2.更新Hn队列当前评分，重新压入全局队列"<<endl;
        //double pseudo_score_update = tfIDF_term(term_n);
        float dist_current = H_list[term_th].top().dist;
        //float score_updated = pseudo_score_update / dist;
        float score_updated = optimal_social_textual_score / (1+dist_current);
        GlobalEntry_Pseudo entry_Q(term_n,dist_current,score_updated);
        Queue.push(entry_Q);
        //cout<<"H"<<term_n<<"当前评分更新为"<<score_updated<<endl;
        current_score_global = Queue.top().score;

        ////3.若poi_current没有处理过，分析之（用于更新resultFinal）
        //cout<<"当前处理的poi为：p"<<entry_p.poi_id<<endl;
        if(poiMark[poi_current] == true) continue;
        poiMark[poi_current] = true;
        userRelated_poiSet.insert(poi_current);

        //尝试更新 resultFinal
        double simT = textRelevance(user.keywords, poi_current_keywords);
        double simD = entry_p.dist;
        double social_textual_score = alpha*1 + (1-alpha)*simT;

        double gsk_score = social_textual_score / (1+simD);
        Result tmpRlt(poi_current,-1,-1,-1,simD,gsk_score,poi_current_keywords);
        if(resultFinal.size()<K){
            resultFinal.push(tmpRlt);
            //cout<<"resultFinal中加入 p"<<poi_current<<endl;
        }
        else{  //更新current topk评分
            Result _tmpNode = resultFinal.top();
            double tmpScore = _tmpNode.score;
            if(tmpScore<gsk_score){
                resultFinal.pop();
                resultFinal.push(tmpRlt);
                float score_rk_update = resultFinal.top().score;
                score_rk = score_rk_update;
                //cout<<"resultFinal中加入 p"<<poi_current<<endl;
                //cout<<"score_rk更新为"<<score_rk<<endl;
            }

        }


        float best_score = Queue.top().score;
        //cout<<"best_score="<<best_score<<",score_rk="<<score_rk<<endl;
        //提前终止条件
        if(resultFinal.size()==K){
            double current_score_rk = resultFinal.top().score;
            int best_term = Queue.top().term_id;
            int term_th = termid2idx[best_term];
            float best_dist = H_list[term_th].top().dist;
            //double dist_best = resultFinal.top().dist;
            double score_bound = optimal_social_textual_score /(best_dist + 1);
            if(score_bound < current_score_rk) break;
        }


    }



    delete []poiMark; delete []poiADJMark;

    endTime = clock();//计时开始

    TopkQueryCurrentResult  topkCurrentResult(resultFinal.top().score,resultFinal);
    //topkCurrentResult.printResults(K);

    //cout << "runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000 << "ms，其中：" << endl;
    //cout << "social runtime is: " << (double)(social_endTime - startTime) / CLOCKS_PER_SEC * 1000 << "ms" << endl;
    /*cout << "io time is: " << io_time / CLOCKS_PER_SEC * 1000000 << "us" << endl;
    cout << "textual time is: " << textual_time / CLOCKS_PER_SEC * 1000000 << "us" << endl;
    cout<<"social_textualRelated_poiSet size="<<userRelated_poiSet.size()<<", count="<<count<<",io_count="<<io_count<<endl;
    cout<<"u"<<user.id<<"'s topk score="<<score_rk<<", resultFinal size = "<<resultFinal.size()<<",具体为："<<endl;*/

    return  topkCurrentResult;
    //printTopkResultsQueue(resultFinal);
}

/*-------------------------------verification 阶段功能函数-----------------------------*/

//用不同的方式验证候选用户是否在RkGSKQ结果中
enum verify_Method{Dijkstra_RAM, Dijkstra_IO, Gtree_top, Gtree_bottom_RAM, Gtree_bottom_IO,NVD_Verify};

double candidate_user_verify_single(int usr_id, int Qk,float a, float alpha,float gsk_score){
    double score_k = 0;
    int verify_way ; //Dijkstra_RAM; //Dijkstra_IO; //  // Gtree_bottom_RAM; // Gtree_bottom_IO;
#ifndef DiskAccess
    #ifdef LV
        verify_way = Dijkstra_RAM;

    #else
        verify_way = Gtree_bottom_RAM;

    #endif
#else

#ifdef LV
    verify_way = Dijkstra_IO;   //如果不是LV则为大图，edge expasion太慢，需基于GTree进行topk结果验证
#else
    verify_way = Gtree_bottom_IO;
#endif
#endif

    if(verify_way == Dijkstra_RAM){
        score_k = topkSDijkstra_verify_usr_memory(transformToQ(Users[usr_id]), Qk,a, alpha,gsk_score,0.0).topkScore;
    }
    else if (verify_way == Dijkstra_IO){

        score_k = topkSDijkstra_verify_usr_disk(transformToQ(getUserFromO2UOrgLeafData(usr_id)),
                                                Qk,a, alpha,gsk_score,0.0).topkScore;
    }
    else if (verify_way == Gtree_bottom_RAM){

        double score_bound = 0.0;
        priority_queue<Result> results = TkGSKQ_bottom2up_verify_memory(usr_id, Qk, a, alpha,gsk_score,score_bound);
        if(results.size()==Qk) {
            score_k = results.top().score;
        }
        else{ //倘若resultFinal为空，说明没有任何元素比score_bound更优
            score_k = 0;
        }

    }
    else if (verify_way == Gtree_bottom_IO){

        double score_bound = 0.0;
        priority_queue<Result> results = TkGSKQ_bottom2up_verify_disk(usr_id, Qk, a, alpha,gsk_score,score_bound);
        if(results.size()==Qk) {
            score_k = results.top().score;
        }
        else{ //倘若resultFinal为空，说明没有任何元素比score_bound更优
            score_k = 0;
        }

    }

    return  score_k;

}


double candidate_user_verify_for_batch(int usr_id, int Qk,float a, float alpha,float gsk_score, float score_bound){
    double score_k = 0;
    int verify_way =-1; //Dijkstra_RAM; //Dijkstra_IO; //  // Gtree_bottom_RAM; // Gtree_bottom_IO;
    User u;

#ifdef LV
    verify_way = Dijkstra_IO;
#else
    verify_way = Gtree_bottom_IO;
#endif


    //verify_way = Gtree_bottom_IO;
    verify_way = NVD_Verify;
    if(verify_way == NVD_Verify){
        User u = getUserFromO2UOrgLeafData(usr_id);
        TopkQueryCurrentResult topk_rs= TkGSKQ_NVD(u, Qk, a, alpha);
        score_k = topk_rs.topkScore;
    }
    else if (verify_way == Gtree_bottom_IO){

        priority_queue<Result> results = TkGSKQ_bottom2up_verify_disk(usr_id, Qk, a, alpha,gsk_score,score_bound);
        if(results.size()==Qk) {
            score_k = results.top().score;
        }
        else{ //倘若resultFinal为空，说明没有任何元素比score_bound更优
            score_k = 0;
        }

    }
    else if(verify_way == Dijkstra_IO){
        score_k = topkSDijkstra_verify_usr_disk(transformToQ(getUserFromO2UOrgLeafData(usr_id)),
                                                Qk,a, alpha,gsk_score, score_bound).topkScore;  //这里后续需完善
    }
    else if(verify_way == Dijkstra_RAM){
        score_k = topkSDijkstra_verify_usr_memory(transformToQ(Users[usr_id]), Qk,a, alpha,gsk_score, score_bound).topkScore;
    }

    else if (verify_way == Gtree_bottom_RAM){

        //priority_queue<Result> results = TkGSKQ_bottom2up_memory(usr_id, Qk, a, alpha);
        priority_queue<Result> results = TkGSKQ_bottom2up_verify_memory(usr_id, Qk, a, alpha,gsk_score,score_bound);
        if(results.size()==Qk) {
            score_k = results.top().score;
        }
        else{ //倘若resultFinal为空，说明没有任何元素比score_bound更优
            score_k = 0;
        }

    }

    return  score_k;

}



