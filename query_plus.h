//
// Created by jins on 10/31/19.
//

#ifndef MAXINFRKGSKQ_1_1_QUERY_PLUS_H

#include "query_basic.h"

#include <random>
#include "omp.h"
#include "SFMT/dSFMT/dSFMT.h"

struct timeval upperv;
long long upper_s, upper_e;
long upper_last=0;

struct timeval leafv;
long long leaf_s, leaf_e;
long leaf_last=0;


#define UPPER_START gettimeofday( &upperv, NULL ); upper_s = upperv.tv_sec * 1000000 + upperv.tv_usec;
#define UPPER_PAUSE gettimeofday( &upperv, NULL ); upper_e = upperv.tv_sec * 1000000 + upperv.tv_usec; upper_last += upper_e - upper_s;

#define UPPER_PRINT(T) printf("%s %lld (ms)\r\n", (#T), upper_last /1000 );
#define UPPER_CLEAR upper_last =0;


#define LEAF_START gettimeofday( &leafv, NULL ); leaf_s = leafv.tv_sec * 1000000 + leafv.tv_usec;
#define LEAF_PAUSE gettimeofday( &leafv, NULL ); leaf_e = leafv.tv_sec * 1000000 + leafv.tv_usec; leaf_last += leaf_e - leaf_s;

#define LEAF_PRINT(T) printf("%s %lld (ms)\r\n", (#T), leaf_last /1000 );
#define LEAF_CLEAR leaf_last =0;


//enum test_task{};  //this is for test
enum topk_method{NVD_topk, GTree_topk, DJ_topk, BruteForce_topk};  //this is for topk solution
enum reverse_method{Reverse_ALL, NVD_reverse, NVD_reverse_OP, Group_reverse, BruteForce_reverse};     //this is for revese knn solution chose


int getTopKMethodChose(std::string topk_chose){
    if(topk_chose == "dj")
        return  DJ_topk;
    else if (topk_chose== "gtree")
        return GTree_topk;
    else if(topk_chose == "nvd")
        return NVD_topk;
    else if(topk_chose == "bruteforce")
        return BruteForce_topk;

}




int getReverseMethod(std::string reverse_chose){
    if(reverse_chose == "group")
        return Group_reverse;
    else if (reverse_chose== "nvd")
        return NVD_reverse;
    else if(reverse_chose == "nvd_op")
        return NVD_reverse_OP;
    else if(reverse_chose == "bruteforce")
        return BruteForce_reverse;
    else if(reverse_chose =="all"){
        return Reverse_ALL;
    }

}


int Uc_size = 0;


/*------------------------- The algorithms below are for single RkGSKQ-------------------------*/


//brute force algorithm to retrieve RkGSKQ results (run topk for all users)  check: ok
int Naive_RkGSKQ(int poi_id, int Qk, float a, float alpha){
    TIME_TICK_START
#ifdef Brightkite
    poi_id = 83;
#endif
    int size = 0;
    POI p;
#ifdef DiskAccess
    p = getPOIFromO2UOrgLeafData(poi_id);
#else
    p = POIs[poi_id];
#endif
    int loc = p.Ni;
    vector<int> keywords = p.keywords;
#ifdef GIVENKEY_
    keywords = extractQueryKeyowrds(poi_id,DEFAULT_KEYWORDSIZE);
#endif
    cout<<"query keyword:"<<endl;
    printElements(keywords);

    vector<ResultDetail> results; int cnt =0;
    for (int i=0; i<UserID_MaxKey;i++) {
        User usr;
#ifdef DiskAccess
        usr = getUserFromO2UOrgLeafData(i);
#else
        usr= Users[i];
#endif

        double relevance = 0;
        relevance = textRelevance(usr.keywords, keywords);
        //relevance = 1;
        if (relevance == 0) continue;
        //cout<<"verify usr"<<usr.id<<endl; cnt++;

        //double distance = getDistance(POIs[poi_id].Nj, usr);//距离会变小
        double distance =  usrToPOIDistance(usr,p);//min(d1,d2);

        double social_Score = 0;
        social_Score =  reCalcuSocialScore(usr.id,poiCheckInIDList[poi_id]);  //gazi

        double influence;
        if(social_Score == 0)
            influence =1;
        else
            influence =1 + pow(a,social_Score);


        double social_textual_score =  alpha*influence + (1-alpha)*relevance;
        //double score = GSKScoreFunction(distance,usr.keywords,keywords,usr.id,0, poiCheckInIDList[0], a);
        double gsk_score = social_textual_score / (1+distance);
        //score = getGSKScore(a, POIs[poi_id].Nj, usr.id, usr.keywords, poiCheckInIDList[poi_id]);
        double simD = distance; double simT = relevance; double simS = influence;

        query usr_q = transformToQ(usr);
        double Rk_usr = -1;   double score_bound = -1.0;

#ifdef LV
        #ifdef  DiskAccess
        Rk_usr = topkSDijkstra_verify_usr_disk(usr_q,Qk,a,alpha,gsk_score,0.0).topkScore;
    #else
        Rk_usr = topkSDijkstra_verify_usr_memory(usr_q,Qk,a,alpha,gsk_score,0.0).topkScore;
    #endif

#else
#ifdef  DiskAccess


        priority_queue<Result> topkResults = TkGSKQ_bottom2up_verify_disk(usr.id,Qk,a,alpha,gsk_score,score_bound);
        if(topkResults.size()==Qk)
            Rk_usr = topkResults.top().score;
        else
            Rk_usr = 0;
#else
        priority_queue<Result> topkResults = TkGSKQ_bottom2up_verify_memory(usr.id,Qk,a,alpha,gsk_score,score_bound);
        if(topkResults.size()==Qk)
            Rk_usr = topkResults.top().score;
        else
            Rk_usr = 0;
#endif

#endif

        //if (score > Rk_usr && distance < 10000) {   // 10km范围内！(这个条件会减小结果)
        if (gsk_score > Rk_usr) {
            size++;
            ResultDetail sd(usr.id,influence,distance,relevance,gsk_score,Rk_usr);
            results.push_back(sd);
            //cout<<"u"<<usr.id<<",";
        }
        cnt++;
    }
    cout <<"p"<<p.id<<",RkGSKQ 用户结果个数=" << size << "，具体内容为："<<endl;
    printResults(results);
    cout<< "topk执行个数"<<cnt<<endl;
    TIME_TICK_END
    TIME_TICK_PRINT("Naive_RkGSKQ running time:")
    return 0;

}

//more fast
double getGSKScore_o2u_phl(int a, double alpha, POI poi, User user){   //poi与user的GSK评分

    set<int> uKeyset; set<int> qKeyset;
    //for(int term : Users[usr_id].keywords) uKeyset.insert(term);
    //for(int term :qKeys)  qKeyset.insert(term);
    //double relevance = textRelevance(uKeyset,qKeyset);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    int _id = poi.id;
    vector<int> qKeys = poi.keywords;
    vector<int> check_usr = poiCheckInIDList[poi.id];

    double relevance = textRelevance(user.keywords,qKeys);
    if(relevance == 0){
        return  0.0;
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


    double distance = usrToPOIDistance_phl(user, poi);   //这里要注意区分 getDistance(vertex, vertex)
    double social_textual_score = alpha*influence + (1.0-alpha)*relevance;
    double gsk_score = social_textual_score / (1+distance);
    //if(poi.id== 57&& user.id ==625 ) cout<<"getGSKScore_o2u，发现u57,p625,distance="<<distance<<",relevance="<<relevance<<"influence="<<influence<<",score="<<gsk_score<<endl;
    return gsk_score;

}


double getSKScore_o2u_phl(int a, double alpha, POI poi, User user){   //poi与user的GSK评分

    set<int> uKeyset; set<int> qKeyset;
    //for(int term : Users[usr_id].keywords) uKeyset.insert(term);
    //for(int term :qKeys)  qKeyset.insert(term);
    //double relevance = textRelevance(uKeyset,qKeyset);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    vector<int> qKeys = poi.keywords;
    vector<int> check_usr = poiCheckInIDList[poi.id];

    double relevance = textRelevance(user.keywords,qKeys);
    if(relevance == 0){
        return  0.0;
    }
    double influence = 1;


    double distance = usrToPOIDistance_phl(user, poi);   //这里要注意区分 getDistance(vertex, vertex)
    double social_textual_score = alpha*influence + (1.0-alpha)*relevance;
    double gsk_score = social_textual_score / (1+distance);
    //if(poi.id== 57&& user.id ==625 ) cout<<"getGSKScore_o2u，发现u57,p625,distance="<<distance<<",relevance="<<relevance<<"influence="<<influence<<",score="<<gsk_score<<endl;
    return gsk_score;

}




//给定某一路网位置的查询
SingleResults RkGSKQ_Single_loc_Plus(int loc, vector<int> keywords, vector<int> checkin_usr, int Qk, float a, float alpha){
    TIME_TICK_START
    cout<<"RkGSKQ_Single_Plus begin!"<<endl;

    int block_num_pre = block_num;

    vector<int> qKey = keywords;

    set<int> keyword_Set;
    for(int t: keywords)
        keyword_Set.insert(t);

    set<int> checkin_usr_Set;
    for(int u: checkin_usr)
        checkin_usr_Set.insert(u);
    //vector<int> candidateLeaf;
    // node by node
    vector<int> Leaf_node_User,Leaf_node_Poi;

    vector<ResultDetail> result_User;
    //提取每个查询关键词 在leafNodeInv表中对应信息

    //对每个有usr的叶节点进行上下界的评估
    int prune=0; int unprune=0; int user_access=0; int usr_prune =0;
    int leaf_prune=0; int leaf_unprune=0;
    int may_prune = 0;
    // 获取与查询关键词相关的用户叶子节点
    set<int> related_usrLeafSet;
    int t3=0; int root=0;

    //set<int> verification_User;  //凝练集
    set<VerifyEntry> verification_User;
    //map<int,map<int,set<VerifyEntry>>> verification_Map;  //凝练结果表   <term_id, leaf_id, <user_id....>>
    map<int,set<VerifyEntry>> verification_Map;    //凝练结果表   <term_id, <user_id....>>
    //对每个关键词做一次filter处理
    set<int> totalUsrSet;
    Filter_START
    for(int term_id:keywords){
        vector<int> related_usr_node;
        vector<int> related_poi_nodes;
        // 获取与查询关键词相关的****用户****所在子图对应的上层节点
        //cout<<"for w"<<term_id<<endl;
        set<int> nodes = usrNodeTermChild[root][term_id]; //get the child list of root whose text contains term_id
        for (int n: nodes) {
            //cout<<"node"<<n<<endl;
            related_usr_node.push_back(n);
            t3++;
        }
        //cout << "(usr)related_nodeSet size =" << related_usr_node.size() << ",usr node access=" << t3;
        // 获取与查询关键词相关的****兴趣点***** 所在子图对应的上层节点
        t3 = 0;
        set<int> nodes2 = poiNodeTermChild[root][term_id]; //get the child list of root whose text contains term_id
        for (int n: nodes2) {
            //cout<<"node"<<n<<endl;
            related_poi_nodes.push_back(n);
            t3++;
        }
        //cout << ",(poi)related_nodeSet size =" << related_poi_nodes.size() << ",poi node access=" << t3;

        // 自底向上： 获取与查询关键词相关的兴趣点所在各节点（首先是叶子节点）中 的非重复相关兴趣点（相关兴趣点计数）
        int i = 0;
        set<int> related_poi_leafSet;
        set<int> related_usr_leafSet;
        map<int, set<int>> node_poiCnt;
        map<int, set<int>> node_occurence_leaf;
        //统计相关兴趣点的信息
        set<int> term_OLeaf = poiTerm_leafSet[term_id];
        for (int leaf: term_OLeaf) {
            related_poi_leafSet.insert(leaf);  //提取相关poi_leaf
            for (int poi: leafPoiInv[term_id][leaf]) {
                //记录该leaf中的相关poi(非重复）
                node_poiCnt[leaf].insert(poi);
                //记录该leaf对应上层祖先节点中的相关poi(非重复）
                int current = leaf;
                int father;
                while (current != 0) {
                    father = GTree[current].father;
                    node_poiCnt[father].insert(poi);
                    node_occurence_leaf[father].insert(leaf);
                    current = father;
                }
            }
        }
        //统计相关用户的信息
        map<int, set<int>> node_relateUserCnt; //  < leaf_id, <usr_id...> >
        set<int> term_ULeaf = usrTerm_leafSet[term_id];
        map<int, set<int>> node_usr_occurence_leaf;
        int totalUsr =0;
        for (int leaf: term_ULeaf) {
            related_usr_leafSet.insert(leaf);  //提取相关usr_leaf
            for (int usr: leafUsrInv[term_id][leaf]) {
                //记录该leaf中的相关用户(非重复）
                node_relateUserCnt[leaf].insert(usr);
                totalUsrSet.insert(usr);
                //记录该leaf对应上层祖先节点中的相关poi(非重复）
                int current = leaf;
                int father;
                while (current != 0) {
                    father = GTree[current].father;
                    node_relateUserCnt[father].insert(usr);
                    node_usr_occurence_leaf[father].insert(leaf);
                    current = father;
                }
            }
        }

        //cout << " ,relate poi leaf size=" << related_poi_leafSet.size() << endl;
        //cout << " ,relate usr leaf size=" << related_usr_leafSet.size() << endl;

        set<int> related_children;
        set<int> related_usr;
        set<int> related_leaves;
        //verify the user node from upper to bottom
        int level = 1;  int cnt_ = 0;
        int nodeAccess = 0;
        while (related_usr_node.size() > 0) {
            related_children.clear();
            for (int node:related_usr_node) {    //
                nodeAccess++;
                //cout<<"iteration"<<i<<":"<<endl;
                i++;
                double max_Score = 0.0;

                if (GTree[node].isleaf) {  //若该节点为用户叶子节点
                    //cout<<"usr_leaf"<<node<<endl;
                    //计算 MaxFGSK(q, U)
                    double min_Distance = getMinMaxDistance(loc, node).min;  //这里可以被优化

                    double max_Score = getMaxScoreOfLeafOnce(loc, node, min_Distance, keyword_Set, checkin_usr_Set, a);

                    double Rt_U = 0;
                    LCLResult lclResult;
                    LCLResult lclResult2;

                    //有cache
                    if(GTree[node].termLCLDetail[term_id].topscore > 0){
                        lclResult = GTree[node].termLCLDetail[term_id];
                        //Rt_U = GTree[node].termLCLRt[term_id];
                    }
                        //无cache
                    else{

                        //lclResult = getLeafNode_LCL_object(node, keyword_Set, GTree[node].termParentLCLDetail[term_id].LCL , node_poiCnt, Qk, a, max_Score);
                        if(GTree[node].termParentLCLDetail[term_id].LCL.size()>0)
                            lclResult = getLeafNode_UCL_object(node, term_id, GTree[node].termParentLCLDetail[term_id].LCL , node_poiCnt, Qk, a, alpha, max_Score);
                        else{

                            lclResult = getLeafNode_TopkSDijkstra(node, term_id,Qk,a, alpha, max_Score);
                            //cout<<"getLeafNode_TopkSDijkstra"<<endl;

                        }
                        //lclResult = getLeafNode_TopkSDijkstra(node, term_id,Qk,a, alpha, max_Score);
                        //cout<<"max_Score="<<max_Score<<",lcl:"<<lclResult.topscore<<",dj: "<<lclResult2.topscore;
                        GTree[node].termLCLDetail[term_id] = lclResult;
                    }

                    Rt_U = lclResult.topscore;
                    //Rt_U = 0;
                    //用户叶子节点剪枝条件判断
                    //cout<<"max_Score="<<max_Score<<",Rt_U="<<Rt_U<<endl;
                    if (max_Score > Rt_U) {  // check if leaf can be pruned
                        unprune++;
                        leaf_unprune++;

                        //若该节点为用户叶子节点
                        if (leafUsrInv[term_id][node].size() != 0) {
                            for (int usr_id: leafUsrInv[term_id][node]) {
                                double gsk_score = getGSKScore_loc(a, alpha,loc, usr_id, qKey,checkin_usr);
                                //用Rt_U进行用户剪枝条件初步判断
                                if (gsk_score > Rt_U) {
                                    double u_lcl = 0;
                                    u_lcl = get_usrLCL_light(usr_id,lclResult,node_poiCnt,Qk,a, alpha, gsk_score);
                                    //u_lcl = get_usrLCL_heavy(usr_id,lclResult,node_poiCnt,Qk,a,gsk_score);
                                    //cout<<"u"<<usr_id<<", gsk_score="<<gsk_score<<",u_lcl="<<u_lcl<<endl;
                                    if(u_lcl > gsk_score){
                                        //cout<<"用户被"<<usr_id<<"被剪枝, gsk_score="<<gsk_score<<",u_lcl="<<u_lcl<<endl;
                                        usr_prune++;
                                        continue;   //该用户不是反top-k结果

                                    }
                                    VerifyEntry ve = VerifyEntry(usr_id,gsk_score,Rt_U);
                                    verification_User.insert(ve);
                                    verification_Map[term_id].insert(ve);

                                } else {
                                    //cout<<"usr"<<usr_id<<"prune"<<endl;
                                    //cout<<"gsk_Score="<<gsk_score<<endl;
                                    //cout<<"Rt_U="<<Rt_U<<endl;

                                }
                            }
                        }

                    }
                        //用户叶子节点整体被剪枝
                    else {
                        //cout<<"leaf "<<node<<"被剪枝！"<<endl;
                        prune++;
                        leaf_prune++;
                    }
                }
                    //该节点为用户上层节点
                else {
                    int usrNum = node_relateUserCnt[node].size();
                    //cout<<"该节点有相关用户："<< usrNum <<endl;
                    if(usrNum < 5){
                        for (int usr_id: node_relateUserCnt[node]){ //用户稀疏节点，直接相关用户进入凝练集
                            double gsk_score = getGSKScore_loc(a, alpha, loc, usr_id, qKey,checkin_usr);
                            VerifyEntry ve = VerifyEntry(usr_id,gsk_score,0.0);
                            verification_User.insert(ve);
                            verification_Map[term_id].insert(ve);
                        }
                        //用户稀疏节点，直接加入相关叶子节点(太慢)
                        /*for (int leaf_id: node_usr_occurence_leaf[node]){
                            related_children.insert(leaf_id);
                        }*/
                        prune++;
                        //cout<<"用户稀疏节点 usr_node"<<node<<",prune"<<endl;   //剪枝整棵子树

                    }
                    else{
                        //cout<<"用户密集型节点，需进行剪枝判断"<<endl;
                        max_Score = getMaxScoreOfNodeOnce_loc(loc, node, keyword_Set, checkin_usr_Set, a);
                        double Rt_U =0;
                        LCLResult lcl;
                        if(GTree[node].termLCLDetail[term_id].topscore>0){
                            LCLResult lcl = GTree[node].termLCLDetail[term_id];
                            Rt_U = lcl.topscore;
                        }

                        else {
                            lcl = getUpperNode_LCL(node, term_id, related_poi_nodes, node_poiCnt, Qk, a,  alpha, max_Score);
                            GTree[node].termLCLDetail[term_id] = lcl;
                            Rt_U = lcl.topscore;
                        }

                        //cout<<"Rt_U="<<Rt_U<<endl;
                        //cout<<"get the LCL of u leaf"<<leaf_id<<endl;
                        if (max_Score > Rt_U) {  // 该子树下所有用户不能被剪枝，加入列表后续进一步验证
                            unprune++;
                            if (usrNodeTermChild[node][term_id].size() != 0) {
                                for (int child: usrNodeTermChild[node][term_id]){
                                    //孩子节点继承父节点的lcl内容
                                    GTree[child].termParentLCLDetail[term_id] = lcl;
                                    //GTree[child].termParentFilter[term_id] = filter;
                                    related_children.insert(child);

                                }

                            }
                            //cout<<"usr_node"<<node<<",unprune"<<endl;
                            //cout<<"Rt_U="<<Rt_U<<endl;
                        }
                            //整棵用户子树被剪枝
                        else {
                            prune++;
                            cout<<"usr_node"<<node<<",prune"<<endl;
                            //cout<<"max_Score="<<max_Score<<endl;
                            //cout<<"Rt_U="<<Rt_U<<endl;

                        }
                    }
                }

            }
            // upfold the related poi nodes
            //set<int> related_poi_nodes_next;
            //set<int> related_poi_leaf_next;
            vector<int> related_poi_nodes_next;
            vector<int> related_poi_leaf_next;
            for (int pn: related_poi_nodes) {
                if (!GTree[pn].isleaf) {
                    if (poiNodeTermChild[pn][term_id].size() != 0) {
                        for (int child: poiNodeTermChild[pn][term_id]) {
                            //cout<<"insert node"<<child<<endl;
                            //related_poi_nodes_next.insert(child);
                            related_poi_nodes_next.push_back(child);
                        }
                    }

                } else {
                    //related_poi_leaf_next.insert(pn);
                    related_poi_leaf_next.push_back(pn);
                }
            }
            related_poi_nodes.clear();
            for (int cn:related_poi_nodes_next)
                related_poi_nodes.push_back(cn);
            for (int lf:related_poi_leaf_next)
                related_poi_nodes.push_back(lf);

            // upfold the related usr nodes
            related_usr_node.clear();

            for (int cn:related_children){
                related_usr_node.push_back(cn);
            }

            if(related_children.size()==1){
                cnt_++;
                if(cnt_==5){
                    //cout<<"直接跳到叶节点"<<endl;
                    //getchar();
                    related_usr_node.clear();
                    for(int node: related_children){
                        if(GTree[node].isleaf)
                            related_usr_node.push_back(node);
                        else{
                            for(int leaf: node_usr_occurence_leaf[node])
                                related_usr_node.push_back(leaf);
                        }
                    }
                    //cout<<"新加入"<<related_usr_node.size()<<endl;
                    //getchar();

                }

            }

            //cout<<"第"<<level<<"层用户节点剪枝验证结束！"<<endl;
            //cout<<"related_poi_node size in next iteration:"<<related_poi_nodes.size()<<endl;
            //cout<<"related_usr_node size in next iteration:"<<related_usr_node.size()<<endl;
            if(related_usr_node.size()==0) break;
            level++;
        }//end while
    }// end for each keywords
    Filter_END


    cout<<"verification_usr="<<verification_User.size()<<endl;
    user_access= 0;
    int topK_count = 0;
    //开始凝练过程
    Refine_START
    //for(int v: verification_User){
    for(VerifyEntry v: verification_User){
        user_access++;
        int usr_id = v.u_id;
        double Rk_u;
        double gsk_score = v.score;//getGSKScore(a,loc, usr_id, qKey,checkin_usr);

        //continue;
        //是否执行top-k下界更新机制
        bool flag = true;
        if(Users[usr_id].topkScore_current > -1){  //有cache

            Rk_u = Users[usr_id].topkScore_current;
            //if(usr_id==800){
            //cout<<"cache u800="<< Rk_u<<endl;
            //}
            if(gsk_score < Rk_u){
                //cout<<"必然不是反topk结果，prune u"<<usr_id<<", gsk_score="<<gsk_score<<",Rk_U="<<v.rk_current<<"Rk_u"<<Rk_u<<endl;
                continue;
            }  //必然不是topk结果

            else{   //已获得最终topk结果

                TopkQueryCurrentResult result = topkSDijkstra_verify_usr_memory(transformToQ(Users[usr_id]),Qk,a, alpha,gsk_score,0.0);
                topK_count++;
                Rk_u = result.topkScore;
                vector<Result> current_elements = result.topkElements;
                //updateTopk(a, usr_id, Qk, verification_Map, current_elements);


            }

        }
        else{   //无cache

            TopkQueryCurrentResult result = topkSDijkstra_verify_usr_memory(transformToQ(Users[usr_id]),Qk,a, alpha,gsk_score,0.0);
            Rk_u = result.topkScore;
            vector<Result> current_elements = result.topkElements;
            topK_count++;
            //updateTopk(a, usr_id, Qk, verification_Map, current_elements);

        }


        if(gsk_score < Rk_u){
            //prune;
            //cout<<"prune u"<<usr_id<<", gsk_score="<<gsk_score<<",Rk_U="<<v.rk_current<<"Rk_u"<<Rk_u<<endl;
        }
        else{
            //usr_id的topk结果当前已知,用该topk结果更新U.LCL
            ResultDetail rd(usr_id,0,0,0, gsk_score,Rk_u);
            result_User.push_back(rd);
        }

    }
    Refine_END

    //对每个候选usr进行分析并打印结果
    SingleResults results =  resultsAnalysis_loc(result_User,keywords,checkin_usr,loc, Qk,a,alpha);
    printResults(results.candidateUsr);
    cout<<"用户验证个数 =" << user_access <<endl;
    cout<<"与兴趣点关键词相关用户总数"<< totalUsrSet.size() <<endl;
    cout<<"topk查询执行个数 =" << topK_count <<endl;
    cout<<"用户结果个数= "<<results.candidateUsr.size()<<"，内容："<<endl;
    cout<<"g-tree node prune ="<<prune<<" ,unprune="<<unprune<<",usr_prune="<<usr_prune<<endl;
    double ratio = (1.0*prune) /(prune+unprune);
    double ratio2 = (1.0*leaf_prune) /(leaf_prune+leaf_unprune);
    cout<<"节点剪枝率= "<<ratio<<"，其中叶子节点剪枝率="<<ratio2<<endl;
    //cout<<"g-tree node access ="<<nodeAccess<<endl;
    TIME_TICK_END
    Filter_PRINT("filter lasts")
    Refine_PRINT("refinement lasts")
    TIME_TICK_PRINT("RkGSKQ_Single_Plus's ")

    return results;

}
//针对某一兴趣点
//内存算法

SingleResults RkGSKQ_Single_poi_Plus(POI poi, vector<int> keywords, vector<int> checkin_usr, int Qk, float a, float alpha){
    TIME_TICK_START

    cout<<"RkGSKQ_Single_poi_Plus !"<<poi.id<<endl;


    //提取该兴趣点的关键词与签到用户信息
    set<int> keyword_Set;
    for(int t: keywords)
        keyword_Set.insert(t);

    set<int> checkin_usr_Set;
    for(int u: checkin_usr)
        checkin_usr_Set.insert(u);
    vector<int> qKey = keywords;

    vector<int> Leaf_node_User,Leaf_node_Poi;

    vector<ResultDetail> result_User;
    //提取每个查询关键词 在leafNodeInv表中对应信息

    //对每个有usr的叶节点进行上下界的评估
    int prune=0; int unprune=0; int user_access=0; int usr_prune =0;
    int leaf_prune=0; int leaf_unprune=0;
    int may_prune = 0;
    // 获取与查询关键词相关的用户叶子节点
    set<int> related_usrLeafSet;
    int t3=0; int root=0;

    //set<int> verification_User;  //凝练集
    set<VerifyEntry> verification_User;
    //map<int,map<int,set<VerifyEntry>>> verification_Map;  //凝练结果表   <term_id, leaf_id, <user_id....>>
    map<int,set<VerifyEntry>> verification_Map;    //凝练结果表   <term_id, <user_id....>>
    //对每个关键词做一次filter处理
    set<int> totalUsrSet;
    int nodeAccess = 0;
    vector<int> related_usr_node;
    map<int,vector<int>> related_poi_nodes;

    for(int term:keyword_Set)
        for(int usr: invListOfUser[term])
            totalUsrSet.insert(usr);


    Filter_START

    //for(int term_id:keywords){
    // 获取与查询关键词相关的****用户***** 所在子图对应的上层节点

    vector<int> unodes = GTree[root].children; //get the child list of root whose text contains term_id
    for (int n : unodes) {
        //cout<<"node"<<n<<endl;
        related_usr_node.push_back(n);
    }
    // 获取与查询关键词相关的****兴趣点***** 所在子图对应的上层节点
    for(int keyword: keyword_Set){
        for(int child: GTree[root].inverted_list_o[keyword]){
            related_poi_nodes[keyword].push_back(child);
        }
    }



    set<int> related_children;
    set<int> related_usr;
    set<int> related_leaves;
    //verify the user node from upper to bottom
    int level = 1;  int cnt_ = 0;
    int term_id ; bool _flag = false;
    int lcl_leafCount = 0;
    while (related_usr_node.size() > 0) {

        related_children.clear();
        for (int node:related_usr_node) {    //


            //TEST_START
            set<int> inter_Key  = obtain_itersection_jins(GTree[node].userUKeySet, keyword_Set);
            //TEST_END
            //TEST_DURA_PRINT("求交集jins")
            //该节点为用户上层节点
            if (!GTree[node].isleaf) {

#ifdef TRACK
                cout<<"访问 node "<<node<<" ,covered query keywords:";
                //这里可以优化, 用户偏好关键词并集与查询关键词的交集
                printSetElements(inter_Key);
#endif

                int usrNum = 0; //node_relateUserCnt[node].size();
                //cout<<"该节点"<<"有相关用户："<< usrNum <<endl;
                if(true){
                    //cout<<"用户密集型节点，需进行剪枝判断"<<endl;
                    int loc = poi.Nj;
                    double max_Score =
                            //getMaxScoreOfNodeOnce_poi(poi, node, keyword_Set, checkin_usr_Set, alpha, a);
                            getMaxScoreOfNodeOnce_poi_InterKey(poi, node, inter_Key, checkin_usr_Set, alpha, a);
                    double Rt_U =100;

#ifdef TRACK
                    cout<<"max_Score ="<< max_Score<<endl;
#endif
                    for(int term: inter_Key){
                        LCLResult lcl;
                        map<int,LCLResult>::iterator iter;
                        iter = GTree[node].termLCLDetail.find(term);
                        if(iter != GTree[node].termLCLDetail.end() ){  //先访问cache
                            lcl = iter->second;
                        } else{ //先访问cache  若无再计算
                            lcl = getUpperNode_termLCL(node, term, related_poi_nodes, Qk, a,  alpha, max_Score);
                            GTree[node].termLCLDetail[term] = lcl;
#ifdef LCL_LOG
                            cout<<"完成 n"<<node<<" 的lcl计算（getUpperNode_termLCL）"<<endl;
                            //printLCLResultsInfo(lcl);

#endif

                        }

                        double current_score = lcl.topscore;
                        Rt_U = min(Rt_U,current_score);

                    }


                    if (max_Score > Rt_U) {  // 该子树下所有用户不能被剪枝，加入列表后续进一步验证
                        unprune++;
                        for(int term: inter_Key){
                            if (GTree[node].inverted_list_u[term].size() != 0) {
                                for (int child: GTree[node].inverted_list_u[term]){  //kobe
                                    //孩子节点继承父节点的lcl内容
                                    GTree[child].termParentLCLDetail[term] = GTree[node].termLCLDetail[term];;
                                    related_children.insert(child);

                                }

                            }
                        }
#ifdef TRACK
                        if(Rt_U*a > max_Score){
                            //cout<<"节点"<<node<<" unpruned, max_score="<<max_Score<<", Rt_U="<<Rt_U<<endl;
                            //getchar();
                        }

#endif

                    }
                        //整棵用户子树被剪枝
                    else {
                        prune++;
#ifdef PRUNELOG
                        cout<<"节点"<<node<<"下所有用户都被剪枝！"<<endl;
#endif

                    }
                }
            }
            else if (GTree[node].isleaf) {  //若该节点为用户叶子节点

#ifdef TRACK

                cout<<"访问 usr_leaf "<<node<<", covered query keywords:";
                printSetElements(inter_Key);
#endif
                //计算 MaxFGSK(q, U)

                int loc = poi.Nj;
                double min_Distance = getMinMaxDistance(loc, node).min;  //这里可以被优化


                double max_Score =
                        getMaxScoreOfLeaf_InterKey(loc, node, min_Distance, inter_Key, checkin_usr_Set, alpha,a);
#ifdef TRACK
                cout<<"max_Score="<<max_Score<<endl;
#endif

                //jordan
                MultiLCLResult lclResult_multi;
                vector<int> Ukeys;
                //这里可以优化
                for(int key:inter_Key)
                    Ukeys.push_back(key);

                //有cache
                if (GTree[node].nodeLCL_multi.size() > 0){
                    lclResult_multi = GTree[node].nodeLCL_multi;
                    //Rt_U = GTree[node].termLCLRt[term_id];
                }
                    //无cache
                else{

                    map<set<int>,MultiLCLResult>::iterator iter;
                    iter = GTree[node].cacheMultiLCLDetail.find(inter_Key);
                    if (iter!=GTree[node].cacheMultiLCLDetail.end()){
                        lclResult_multi = iter->second;
                    }
                    else{
#ifdef LV
                        lclResult_multi = getLeafNode_LCL_Dijkstra_multi(node, Ukeys, Qk, a, alpha, max_Score);
#endif

#ifndef LV
                        lclResult_multi = getLeafNode_LCL_bottom2up_multi(node, Ukeys, Qk, a, alpha, max_Score);
#endif

                        GTree[node].cacheMultiLCLDetail[inter_Key] = lclResult_multi;

                        lcl_leafCount++;
                    }


                }

                //获得 usr_lcl update_o_set
                set<int> update_o_set;
                for(LCLResult lcl: lclResult_multi){
                    priority_queue<LCLEntry> LCL = lcl.LCL;
                    while(LCL.size() > 0){
                        LCLEntry le = LCL.top();
                        LCL.pop();
                        int o = le.id;
                        update_o_set.insert(o);
                        //cout<<"o"<<o<<",";
                    }

                }
#ifdef TRACK
                cout<<", update_o_set size="<<update_o_set.size()<<endl;

#endif


                double Rt_U = -1.0;

                if(update_o_set.size()<Qk){
                    //cout<<"non-enough(empty) update o set"<<endl;
                    Rt_U = 0.0;
                }
                else{
                    for(LCLResult lclResult : lclResult_multi){
                        double _ss = lclResult.topscore;
                        if(-1.0==Rt_U)
                            Rt_U =_ss;
                        else
                            Rt_U = min(Rt_U,_ss);
                    }
                }
                //cout<<"Rt_U="<<Rt_U<<endl;

                //用户叶子节点剪枝条件判断
                //cout<<"用户叶子节点剪枝条件判断: max_Score="<<max_Score<<",Rt_U="<<Rt_U<<endl;
                if (max_Score > Rt_U) {  // check if leaf can be pruned
                    unprune++;
                    leaf_unprune++;
#ifdef TRACK
                    if(a*Rt_U > max_Score){
                        //cout<<"叶节点"<<node<<" unpruned, max_score="<<max_Score<<", Rt_U="<<Rt_U<<endl;

                    }
#endif

                    //若该节点为用户叶子节点
                    set<int> relate_usr_fromLeaf;
                    for(int term: Ukeys) {
                        for(int u : GTree[node].inverted_list_u[term]){
                            relate_usr_fromLeaf.insert(u);

                        }
                    }

                    for (int usr_id: relate_usr_fromLeaf) {

                        double gsk_score =-1;


#ifdef GIVENKEY
                        gsk_score = getGSKScore_o2u_givenKey(a,alpha, poi, Users[usr_id], qKey,checkin_usr);  //这里是给定关键词，而后结合该兴趣点的loc与check-in 进行的评分（givenKeyword）

#else
                        gsk_score = getGSKScore_o2u(a, alpha, poi, Users[usr_id]);  //这里注意是兴趣点的真实评分， 与给定关键词时的评分 不一样了
#endif



                        //用Rt_U进行用户剪枝条件初步判断
                        double u_lcl_rt =0;
                        if (gsk_score > Rt_U) {

                            if(update_o_set.size()>0) {
                                priority_queue<ResultCurrent> u_LCL = get_usrLCL_update_Plus(usr_id,
                                                                                             update_o_set, Qk, a, alpha, gsk_score);  //获取该用户的lcl

                                if (u_LCL.size() < Qk) u_lcl_rt = 0.0;
                                else if (u_LCL.size() == Qk) u_lcl_rt = u_LCL.top().score;
                            }

                            if(u_lcl_rt > gsk_score){
                                //cout<<"用户被"<<usr_id<<"被剪枝, gsk_score="<<gsk_score<<",u_lcl_rt="<<u_lcl_rt<<endl;
                                usr_prune++;
                                continue;   //该用户不是反top-k结果

                            }

                            Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current,u_lcl_rt);

                            VerifyEntry ve = VerifyEntry(usr_id,gsk_score,u_lcl_rt);
                            verification_User.insert(ve);
                            verification_Map[term_id].insert(ve);
#ifdef TRACK
                            cout<<"发现候选用户 u"<<usr_id<<", gsk_score="<<gsk_score<<",u_lcl_rt="<<u_lcl_rt<<endl;
#endif
                        } else {
#ifdef PRUNELOG
                            cout<<"用户"<<usr_id<<" is pruned"<<endl;
#endif
                            //cout<<"gsk_Score="<<gsk_score<<endl;
                            //cout<<"Rt_U="<<Rt_U<<endl;

                        }

                    }

                }
                    //用户叶子节点整体被剪枝
                else {
#ifdef PRUNELOG
                    cout<<"leaf "<< node<<" 被剪枝！"<<endl;
#endif
                    prune++;
                    leaf_prune++;
                }
            }

        }
        level++;
        // upfold the related poi nodes
        //set<int> related_poi_nodes_next;
        //set<int> related_poi_leaf_next;
        map<int,vector<int>> related_poi_nodes_next;
        map<int,vector<int>> related_poi_leaf_next;
        for(int term: keyword_Set){
            for (int pn: related_poi_nodes[term]) {
                if (!GTree[pn].isleaf) {
                    if (GTree[pn].inverted_list_o[term].size() != 0) {
                        for (int child: GTree[pn].inverted_list_o[term]) {
                            //cout<<"insert node"<<child<<endl;
                            //related_poi_nodes_next.insert(child);
                            related_poi_nodes_next[term].push_back(child);
                        }
                    }

                } else {
                    //related_poi_leaf_next.insert(pn);
                    related_poi_leaf_next[term].push_back(pn);
                }
            }
        }
        related_poi_nodes.clear();
        for(int term:keyword_Set){
            for(int cn: related_poi_nodes_next[term])
                related_poi_nodes[term].push_back(cn);
            for (int lf:related_poi_leaf_next[term])
                related_poi_nodes[term].push_back(lf);
        }


        // upfold the related usr nodes
        related_usr_node.clear();

        for (int cn:related_children){
            related_usr_node.push_back(cn);
        }

        if(related_children.size()==1){
            cnt_++;
            if(cnt_==3){
                //cout<<"直接跳到叶节点"<<endl;
                //getchar();
                related_usr_node.clear();
                for(int node: related_children){
                    if(GTree[node].isleaf)
                        related_usr_node.push_back(node);
                    else{
                        set<int> leaves;
                        for(int term: keyword_Set){
                            for(int leaf: GTree[node].term_usrLeaf_Map[term]){
                                leaves.insert(leaf);
                            }
                        }
                        for(int leaf: leaves)
                            related_usr_node.push_back(leaf);
                    }
                }
                //cout<<"新加入"<<related_usr_node.size()<<endl;
                //getchar();

            }

        }
#ifdef TRACK
        cout<<"第"<<level<<"层用户节点剪枝验证结束！"<<endl;
#endif
        //cout<<"related_poi_node size in next iteration:"<<related_poi_nodes.size()<<endl;
        //cout<<"related_usr_node size in next iteration:"<<related_usr_node.size()<<endl;
        if(related_usr_node.size()==0) break;
    }//end while
    //}// end for each keywords
    Filter_END


    cout<<"verification_usr="<<verification_User.size()<<endl;
    user_access= 0;
    int topK_count = 0;
    //开始凝练过程
    Refine_START
    //for(int v: verification_User){
    for(VerifyEntry v: verification_User){

        user_access++;
        int usr_id = v.u_id;
        double Rk_u=0;
        double gsk_score = v.score;//getGSKScore(a,loc, usr_id, qKey,checkin_usr);
        double u_lcl_rt = v.rk_current;
        //continue;
        //是否执行top-k下界更新机制
        bool flag = true;

        if(Users[usr_id].topkScore_current > -1){  //有cache

            Rk_u = Users[usr_id].topkScore_current;

            if(gsk_score < Rk_u){

#ifdef VERIFICATION_LOG
                cout<<"必然不是反topk结果，prune u"<<usr_id<<", gsk_score="<<gsk_score<<",Rk_U(cache)="<<Rk_u<<endl;

#endif
                continue;
            }  //必然不是topk结果

            else{   //已获得最终topk结果
                if(Users[usr_id].topkScore_Final > -1)
                    Rk_u = Users[usr_id].topkScore_Final;
                else{   //更新当前topk结果

                    Rk_u = candidate_user_verify_single(usr_id,Qk,a, alpha,gsk_score);

                    topK_count++;

                }
            }

        }
        else{   //无cache

            Rk_u = candidate_user_verify_single(usr_id,Qk,a, alpha,gsk_score);
            topK_count++;
            //updateTopk(a, usr_id, Qk, verification_Map, current_elements);

        }

#ifdef VERIFICATION_LOG
        cout<< "u"<<v.u_id<<":Rk_u="<<Rk_u<<"gsk(u,o)="<<gsk_score<<endl;
#endif

        if(gsk_score > Rk_u || gsk_score == Rk_u){
            //usr_id的topk结果当前已知,用该topk结果更新U.LCL
            ResultDetail rd(usr_id,0,0,0, gsk_score,Rk_u);
            result_User.push_back(rd);
        }

    }
    Refine_END

    //对每个候选usr进行分析并打印结果
    SingleResults results =  resultsAnalysis_poi(result_User,poi, Qk,a,alpha);

    cout<<"用户验证个数 =" << user_access <<endl;
    cout<<"与兴趣点关键词相关用户总数"<< totalUsrSet.size() <<endl;
    cout<<"topk查询执行个数 =" << topK_count <<endl;
    cout<<"lcl for leaf 执行个数"<<lcl_leafCount<<endl;
    cout<<"p"<<poi.id<<"用户结果个数= "<<results.candidateUsr.size()<<"，内容："<<endl;
    printResults(results.candidateUsr);
    /*cout<<"g-tree node prune ="<<prune<<" ,unprune="<<unprune<<",usr_prune="<<usr_prune<<endl;
    double ratio = (1.0*prune) /(prune+unprune);
    double ratio2 = (1.0*leaf_prune) /(leaf_prune+leaf_unprune);
    cout<<"节点剪枝率= "<<ratio<<"，其中叶子节点剪枝率="<<ratio2<<endl;
    //cout<<"g-tree node access ="<<nodeAccess<<endl;*/
    TIME_TICK_END
    Filter_PRINT("filter lasts")
    Refine_PRINT("refinement lasts")
    TIME_TICK_PRINT("RkGSKQ_Single_Plus end!")
    cout<<"磁盘IO: number of block accessed="<<block_num<<endl;

    return results;

}


//外存算法
SingleResults RkGSKQ_Single_poi_Disk(POI poi, vector<int> keywords, vector<int> checkin_usr, int Qk, float a, float alpha){
    TIME_TICK_START

    cout<<"---------Running RkGSKQ_by_Group_on_GIM-Tree!---------"<<endl;

    printPOIInfo(poi);

    //提取该兴趣点的关键词与签到用户信息
    set<int> keyword_Set;
    for(int t: keywords)
        keyword_Set.insert(t);

    set<int> checkin_usr_Set;
    for(int u: checkin_usr)
        checkin_usr_Set.insert(u);
    vector<int> qKey = keywords;

    vector<int> Leaf_node_User,Leaf_node_Poi;

    vector<ResultDetail> result_User;
    //提取每个查询关键词 在leafNodeInv表中对应信息

    //对每个有usr的叶节点进行上下界的评估
    int upper_prune; int upper_unprune=0;
    int leaf_prune=0; int leaf_unprune=0;
    int user_access=0; int usr_prune =0;

    int may_prune = 0;
    // 获取与查询关键词相关的用户叶子节点
    set<int> related_usrLeafSet;
    int t3=0; int root=0;

    //set<int> verification_User;  //凝练集
    set<VerifyEntry> verification_User;
    //map<int,map<int,set<VerifyEntry>>> verification_Map;  //凝练结果表   <term_id, leaf_id, <user_id....>>
    map<int,set<VerifyEntry>> verification_Map;    //凝练结果表   <term_id, <user_id....>>
    //对每个关键词做一次filter处理
    set<int> totalUsrSet;
    int nodeAccess = 0;
    vector<int> related_usr_node;
    map<int,vector<int>> related_poi_nodes;


    Filter_START

    // 获取与查询关键词相关的****用户***** 所在子图对应的上层节点

    //vector<int> unodes = GTree[root].children; //get the child list of root whose text contains term_id
    //printElements(unodes);
    vector<int> unodes = getGIMTreeNodeData(root,OnlyU).children;
    //printElements(unodes);


    for (int n : unodes) {
        //cout<<"node"<<n<<endl;
        related_usr_node.push_back(n);
    }
    // 获取与查询关键词相关的****兴趣点***** 所在子图对应的上层节点
    for(int keyword: keyword_Set){
        for(int child_entry: getObjectTermRelatedEntry(keyword, root)){
            related_poi_nodes[keyword].push_back(child_entry);
        }
        //printElements(related_poi_nodes[keyword]);
        //printElements(getObjectTermRelatedEntry(keyword, root));
    }

    int lcl_leaf_count =0;
    set<int> related_children;
    set<int> related_usr;
    set<int> related_leaves;
    //verify the user node from upper to bottom
    int level = 1;  int cnt_ = 0;
    int term_id ;
    while (related_usr_node.size() > 0) {

        related_children.clear();
        for (int node:related_usr_node) {    //


            //access GIM-Tree Node
#ifdef TRACK
            if(GTree[node].isleaf){
                cout<<"--------------访问 leaf "<<node<<endl;
            }
            else
                cout<<"--------------访问  node"<<node<<endl;
#endif
            TreeNode tn = getGIMTreeNodeData(node, OnlyU);
            //TEST_DURA_PRINT("求交集jins")
            set<int> inter_Key  = obtain_itersection_jins(tn.userUKeySet, keyword_Set);
#ifdef TRACK
            cout<<"读取 node"<<node<<"数据完毕 ,covered query keywords:";
            //这里可以优化, 用户偏好关键词并集与查询关键词的交集
            printSetElements(inter_Key);
#endif


            //该节点为用户上层节点
            if (!tn.isleaf) {

                UPPER_START

                int usrNum = 0; //node_relateUserCnt[node].size();
                //cout<<"该节点"<<"有相关用户："<< usrNum <<endl;
                if(true){
                    //cout<<"用户密集型节点，需进行剪枝判断"<<endl;
                    int loc = poi.Nj;
                    double max_Score = getMaxScoreOfNodeOnce_poi_InterKey(poi, node, inter_Key, checkin_usr_Set, alpha, a);
                    double Rt_U =100;

#ifdef  TRACK
                    cout<<"max_Score= "<<max_Score<<endl;
#endif

                    for(int term: inter_Key){
                        LCLResult lcl;
                        map<int,LCLResult>::iterator iter;
                        iter = GTree[node].termLCLDetail.find(term);
                        if(iter != GTree[node].termLCLDetail.end() ){  //先访问cache
                            lcl = iter->second;
                        } else{ //先访问cache  若无再计算

                            lcl = getUpperNode_termLCL_Disk(node, term, related_poi_nodes, Qk, a,  alpha, max_Score);
#ifdef TRACK
                            //cout<<"针对t"<<term<<", 完成 n"<<node<<" 的lcl计算（getUpperNode_termLCL_Disk）, lcl_rt="<<lcl.topscore<<endl;
                            //printLCLResultsInfo(lcl);

#endif

                            //GTree[node].termLCLDetail[term] = lcl;  这里取消cache

                        }

                        double current_score = lcl.topscore;
                        Rt_U = min(Rt_U,current_score);

                    }

#ifdef TRACK
                    cout<<"最终算得：max_Score="<<max_Score<<",Rt_U="<<Rt_U<<endl;
#endif
                    if (max_Score > Rt_U) {  // 该子树下所有用户不能被剪枝，加入列表后续进一步验证
                        upper_unprune++;
                        for(int term: inter_Key){
                            //int _size = GTree[node].inverted_list_u[term].size();
                            vector<int> child_entry = getUsrTermRelatedEntry(term,node);
                            //int _size2 = child_entry.size();
                            if (child_entry.size() != 0) {
                                for (int child: child_entry){  //kobe
                                    //孩子节点继承父节点的lcl内容
                                    GTree[child].termParentLCLDetail[term] = GTree[node].termLCLDetail[term];;
                                    related_children.insert(child);

                                }

                            }
                        }
#ifdef TRACK
                        if(Rt_U*a > max_Score){
                            //cout<<"节点"<<node<<" unpruned, max_score="<<max_Score<<", Rt_U="<<Rt_U<<endl;
                            //getchar();
                        }

#endif

                    }
                        //整棵用户子树被剪枝
                    else {
                        upper_prune++;
#ifdef TRACK
                        cout<<"节点"<<node<<"下所有用户都被剪枝！"<<endl;
#endif

                    }
                }

                UPPER_PAUSE
            }
            else if (GTree[node].isleaf) {  //若该节点为用户叶子节点

#ifdef TRACK
                cout<<"访问 usr_leaf "<<node<<", covered query keywords:";
                printSetElements(inter_Key);
#endif
                //计算 MaxFGSK(q, U)

                int loc = poi.Nj;
                double min_Distance = getMinMaxDistance(loc, node).min;  //这里可以被优化


                double max_Score =
                        getMaxScoreOfLeaf_InterKey(loc, node, min_Distance, inter_Key, checkin_usr_Set, alpha,a);

#ifdef  TRACK
                cout<<"max_Score= "<<max_Score<<endl;
#endif

                //jordan
                MultiLCLResult lclResult_multi;
                vector<int> Ukeys;
                //这里可以优化
                for(int key:inter_Key)
                    Ukeys.push_back(key);

                //有cache
                //map<set<int>,MultiLCLResult> :: iterator iter2;
                //iter2 = GTree[node].cacheMultiLCLDetail.find(inter_Key);
                if(GTree[node].cacheMultiLCLDetail.count(inter_Key)>0){  //先访问cache
                    lclResult_multi = GTree[node].cacheMultiLCLDetail[inter_Key];
                }
                    //无cache
                else{

                    map<set<int>,MultiLCLResult>::iterator iter;
                    iter = GTree[node].cacheMultiLCLDetail.find(inter_Key);
                    if (iter!=GTree[node].cacheMultiLCLDetail.end()){
                        lclResult_multi = iter->second;
                    }
                    else{
#ifdef LV
                        lclResult_multi = getLeafNode_LCL_Dijkstra_multi(node, Ukeys, Qk, a, alpha, max_Score);

#else
#ifdef LasVegas
                        lclResult_multi = getLeafNode_LCL_Dijkstra_multi(node, Ukeys, Qk, a, alpha, max_Score);
#else
                        lclResult_multi = getLeafNode_LCL_bottom2up_multi(node, Ukeys, Qk, a, alpha, max_Score);
#endif
#endif

                        //GTree[node].cacheMultiLCLDetail[inter_Key] = lclResult_multi;  取消cache

                        lcl_leaf_count++;
                    }


                }

                //获得 usr_lcl update_o_set
                set<int> update_o_set;
                for(LCLResult lcl: lclResult_multi){
                    priority_queue<LCLEntry> LCL = lcl.LCL;
                    while(LCL.size() > 0){
                        LCLEntry le = LCL.top();
                        LCL.pop();
                        int o = le.id;
                        update_o_set.insert(o);
                        //cout<<"o"<<o<<",";
                    }

                }
#ifdef TRACK
                cout<<", update_o_set size="<<update_o_set.size()<<endl;
#endif


                double Rt_U = -1.0;

                if(update_o_set.size()<Qk){
                    //cout<<"non-enough(empty) update o set"<<endl;
                    Rt_U = 0.0;
                }
                else{
                    for(LCLResult lclResult : lclResult_multi){
                        double _ss = lclResult.topscore;
                        if(-1.0==Rt_U)
                            Rt_U =_ss;
                        else
                            Rt_U = min(Rt_U,_ss);
                    }
                }
                //cout<<"Rt_U="<<Rt_U<<endl;

                //用户叶子节点剪枝条件判断
                //cout<<"用户叶子节点剪枝条件判断: max_Score="<<max_Score<<",Rt_U="<<Rt_U<<endl;
                if (max_Score > Rt_U) {  // check if leaf can be pruned  //
                    //unprune++;
                    leaf_unprune++;
#ifdef TRACK
                    //if(a*Rt_U > max_Score){
                        //cout<<"叶节点"<<node<<" unpruned, max_score="<<max_Score<<", Rt_U="<<Rt_U<<endl;

                    //}
#endif

                    //若该节点为用户叶子节点,取得该节点下关键词相关的用户
                    set<int> relate_usrID_set;
                    for(int term: Ukeys) {
                        vector<int> related_usrID_list = getUsrTermRelatedEntry(term,node);
                        for(int u: related_usrID_list)
                            relate_usrID_set.insert(u);
                    }
                    for(int usr_id: relate_usrID_set) {

                        //int id2 = users[i]; i++;
                        double gsk_score =-1;

                        User user = getUserFromO2UOrgLeafData(usr_id);


#ifdef GIVENKEY
                        //gsk_score = getGSKScore_o2u_givenKey(a,alpha, poi, Users[usr_id], qKey,checkin_usr);  //这里是给定关键词，而后结合该兴趣点的loc与check-in 进行的评分（givenKeyword）
                            gsk_score = getGSKScore_o2u_givenKey(a,alpha, poi, user, qKey,checkin_usr);


#else //joke phl
                        gsk_score = getGSKScore_o2u_phl(a, alpha, poi, user);  //这里注意是兴趣点的真实评分， 与给定关键词时的评分 不一样了
#endif

                        //cout<<"u"<<usr_id<<"gsk_score="<<gsk_score<<endl;
                        //用Rt_U进行用户剪枝条件初步判断
                        double u_lcl_rt =0;
                        if (true) {  //gsk_score > Rt_U

                            if(update_o_set.size()>0) {
                                priority_queue<ResultCurrent> u_LCL = get_usrLCL_update_Plus(usr_id,
                                                                                             update_o_set, Qk, a, alpha, gsk_score);  //获取该用户的lcl

                                if (u_LCL.size() < Qk) u_lcl_rt = 0.0;
                                else if (u_LCL.size() == Qk) u_lcl_rt = u_LCL.top().score;
                            }

                            if(u_lcl_rt > gsk_score){
                                //cout<<"用户被"<<usr_id<<"被剪枝, gsk_score="<<gsk_score<<",u_lcl_rt="<<u_lcl_rt<<endl;
                                usr_prune++;
                                continue;   //该用户不是反top-k结果

                            }

                            Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current,u_lcl_rt);

                            VerifyEntry ve = VerifyEntry(usr_id,gsk_score,u_lcl_rt);
                            verification_User.insert(ve);
                            verification_Map[term_id].insert(ve);
#ifdef TRACK
                            cout<<"发现候选用户 u"<<usr_id<<", gsk_score="<<gsk_score<<",u_lcl_rt="<<u_lcl_rt<<endl;
#endif
                        } else {
#ifdef PRUNELOG
                            cout<<"用户"<<usr_id<<" is pruned"<<endl;
#endif
                            //cout<<"gsk_Score="<<gsk_score<<endl;
                            //cout<<"Rt_U="<<Rt_U<<endl;

                        }
                    }

                }
                    //用户叶子节点整体被剪枝
                else {
#ifdef PRUNELOG
                    cout<<"leaf "<< node<<" 被剪枝！"<<endl;
#endif
                    //prune++;
                    leaf_prune++;
                }
            }

        }
        level++;
        // upfold the related poi nodes
#ifdef TRACK
        cout<<"upfold the related poi nodes"<<endl;
#endif
        map<int,vector<int>> related_poi_nodes_next;
        map<int,vector<int>> related_poi_leaf_next;
        for(int term: keyword_Set){
            //cout<<"for term "<<term<<endl;
            for (int pn: related_poi_nodes[term]) {
                if (!GTree[pn].isleaf) {
#ifdef TRACK
                    //cout<<"读取node"<<pn<<"在t"<<term<<" 下的后继："<<endl;
#endif
                    vector<int> child_entry = getObjectTermRelatedEntry(term,pn);

                    if (child_entry.size() != 0) {

                        for (int child: child_entry) {
                            related_poi_nodes_next[term].push_back(child);
                            //int child2 = child_entry[i]; i++;
                        }
                    }

                } else {
                    //related_poi_leaf_next.insert(pn);
                    related_poi_leaf_next[term].push_back(pn);
                }
            }
        }
#ifdef TRACK
        cout<<"更新related_poi_nodes"<<endl;
#endif
        related_poi_nodes.clear();

        for(int term:keyword_Set){
            //cout<<"for term"<<term<<endl;
            for(int cn: related_poi_nodes_next[term])
                related_poi_nodes[term].push_back(cn);
            for (int lf:related_poi_leaf_next[term])
                related_poi_nodes[term].push_back(lf);
        }
        //cout<<"更新完毕！"<<endl;

#ifdef TRACK
        // upfold the related usr nodes
        cout<<"upfold the related usr nodes"<<endl;
#endif
        related_usr_node.clear();

        for (int cn:related_children){
            related_usr_node.push_back(cn);

        }

        if(related_children.size()==1){
            cnt_++;

            if(cnt_==3){
                //cout<<"直接跳到叶节点"<<endl;
                //getchar();
                related_usr_node.clear();
                for(int node: related_children){
                    if(GTree[node].isleaf)
                        related_usr_node.push_back(node);
                    else{
                        set<int> leaves;
                        for(int term: keyword_Set){
                            vector<int> direct_uLeaves;
                            direct_uLeaves = getUsrTermRelatedLeafNode(term,node);
                            for(int leaf: direct_uLeaves){ //GTree[node].term_usrLeaf_Map[term]
                                leaves.insert(leaf);
                            }
                        }
                        for(int leaf: leaves)
                            related_usr_node.push_back(leaf);
                    }
                }

                //cout<<"新加入"<<related_usr_node.size()<<"具体为："<<endl;
                //printElements(related_usr_node);
                //getchar();

            }

        }
#ifdef TRACK
        cout<<"第"<<level<<"层用户节点剪枝验证结束！"<<endl;
#endif
        //cout<<"related_poi_node size in next iteration:"<<related_poi_nodes.size()<<endl;
        //cout<<"related_usr_node size in next iteration:"<<related_usr_node.size()<<endl;
        if(related_usr_node.size()==0) break;
    }//end while
    //}// end for each keywords
    Filter_END


    cout<<"verification_usr="<<verification_User.size()<<endl;
    user_access= 0;
    int topK_count = 0;
    //开始凝练过程yanz
    Refine_START
    //for(int v: verification_User){
    for(VerifyEntry v: verification_User){

        user_access++;
        int usr_id = v.u_id;
        double Rk_u=0;
        double gsk_score = v.score;//getGSKScore(a,loc, usr_id, qKey,checkin_usr);
        double u_lcl_rt = v.rk_current;
        //continue;
        //是否执行top-k下界更新机制
        bool flag = true;

        if(Users[usr_id].topkScore_current > -1){  //有cache

            Rk_u = Users[usr_id].topkScore_current;

            if(gsk_score < Rk_u){

#ifdef VERIFICATION_LOG
                cout<<"必然不是反topk结果，prune u"<<usr_id<<", gsk_score="<<gsk_score<<",Rk_U(cache)="<<Rk_u<<endl;

#endif
                continue;
            }  //必然不是topk结果

            else{   //已获得最终topk结果
                if(Users[usr_id].topkScore_Final > -1)
                    Rk_u = Users[usr_id].topkScore_Final;
                else{   //更新当前topk结果

                    Rk_u = candidate_user_verify_single(usr_id,Qk,a, alpha,gsk_score);

                    topK_count++;

                }
            }

        }
        else{   //无cache

            Rk_u = candidate_user_verify_single(usr_id,Qk,a, alpha,gsk_score);
            topK_count++;
            //updateTopk(a, usr_id, Qk, verification_Map, current_elements);

        }

#ifdef VERIFICATION_LOG
        cout<< "u"<<v.u_id<<":Rk_u="<<Rk_u<<"gsk(u,o)="<<gsk_score<<endl;
#endif

        if(gsk_score > Rk_u || gsk_score == Rk_u){
            //usr_id的topk结果当前已知,用该topk结果更新U.LCL
            ResultDetail rd(usr_id,0,0,0, gsk_score,Rk_u);
            result_User.push_back(rd);
        }

    }
    Refine_END

    //对每个候选usr进行分析并打印结果
    SingleResults results =  resultsAnalysis_poi(result_User,poi, Qk,a,alpha);
    TIME_TICK_END
    cout<<"用户验证个数 =" << user_access <<endl;
    cout<<"叶节点lcl computation个数："<<lcl_leaf_count<<endl;
    //cout<<"topk查询执行个数 =" << topK_count <<endl;
    cout<<"p"<<poi.id<<"用户结果个数= "<<results.candidateUsr.size()<<"，内容："<<endl;
    printResults(results.candidateUsr);
    //cout<<"上层节点:  prune="<<upper_prune<<" , unprune="<<upper_unprune<<endl;
    //cout<<"叶节点： prune="<<leaf_prune<<" , unprune="<<leaf_unprune<<endl;
    /*cout<<"g-tree node prune ="<<prune<<" ,unprune="<<unprune<<",usr_prune="<<usr_prune<<endl;
    double ratio = (1.0*prune) /(prune+unprune);
    double ratio2 = (1.0*leaf_prune) /(leaf_prune+leaf_unprune);
    cout<<"节点剪枝率= "<<ratio<<"，其中叶子节点剪枝率="<<ratio2<<endl;
    //cout<<"g-tree node access ="<<nodeAccess<<endl;*/

    //UPPER_PRINT("runtime for upper node:") UPPER_CLEAR
    Filter_PRINT("filter lasts")
    Refine_PRINT("refinement lasts")
    TIME_TICK_PRINT("RkGSKQ_Single_Plus end!")

    //cout<<"磁盘IO: number of block accessed="<<block_num<<endl;

    return results;

}





//外存算法，Filtering阶段基于NVD下的lcl, verification阶段基于GIM-Tree topk
void RkGSKQ_NVD_Naive_hash(int poi_id, int K, float a, float alpha) {

    //double direction

    string phl_fileName = getRoadInputPath(PHL);//"../exp/indexes/LV.phl";

    //char phl_idxFileName[255] ; sprintf(phl_idxFileName, "../exp/indexes/%s.phl", road_map);

    phl.LoadLabel(phl_fileName.c_str());


    priority_queue<Result> resultFinal;

    //获取poi信息；
    POI poi = getPOIFromO2UOrgLeafData(poi_id);
    printPOIInfo(poi);

    //获取与poi关键词相关的user的全集列表信息；
    set<int> poiRelated_usr;
    //候选验证用户
    vector<int> candidate_usr;
    vector<double> candidate_gsk;
    vector<double> candidate_rk_current;
    //最终结果： 反top-k最近邻用户
    vector<int> usr_results;
    for(int term: poi.keywords){
        vector<int> term_related_user = getTermUserInvList(term);
        for(int usr: term_related_user)
            poiRelated_usr.insert(usr);
    }
    printf("总共有poiRelated_usr： %d个...\n", poiRelated_usr.size());
    //对poi关键词相关的user进行一一验证
    clock_t startTime, endTime;

    clock_t filter_startTime, filter_endTime;
    clock_t verification_startTime, verification_endTime;

    filter_startTime = clock();
    for(int usr_id: poiRelated_usr){
        //usr_id = 346;
        User user = getUserFromO2UOrgLeafData(usr_id);
        //cout<<"验证";
        //printUsrInfo(user);
        startTime = clock();//计时开始
        //double dist = getDistanceLower_Oracle_PHL(user,poi);
        //double gsk_score = getGSKScore_o2u_distComputed(a,alpha,poi,user,dist);
        double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,user);
        endTime = clock();//计时结束
        ////cout << "getGSKScore_o2u runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;
        cout<<"gsk_score="<<gsk_score;
        startTime = clock();//计时开始
        double u_lcl_rk = getUserLCL_NVD(user,K,a,alpha);
        endTime = clock();//计时结束
        cout<<",u_lcl_rk="<<u_lcl_rk<<endl; ////
        ////cout << "getUserLCL_NVD runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;

        cout << "-------------------------------------------------------" << endl;


        //cout<<",u_lcl_rk="<<u_lcl_rk<<endl;
        if(gsk_score > u_lcl_rk){
            candidate_usr.push_back(usr_id);
            candidate_gsk.push_back(gsk_score);
            candidate_rk_current.push_back(u_lcl_rk);
        }
        //getchar();

    }

    filter_endTime = clock();
    cout<<"候选用户个数="<<candidate_usr.size()<<"，具体为："<<endl;
    printElements(candidate_usr);


    verification_startTime = clock();
    int candidate_size = candidate_usr.size();
    for(int i=0;i<candidate_size;i++){
        int usr_id = candidate_usr[i];
        //cout<<"验证候选用户"<<usr_id<<endl;
        double gsk_score = candidate_gsk[i];
        double score_bound = candidate_rk_current[i];
        double Rk_u = 0;

        priority_queue<Result> results = TkGSKQ_bottom2up_verify_disk(usr_id, K, a, alpha, 1000,0.0); //gsk_score, score_bound);
        if (results.size() == K) {
            Rk_u = results.top().score;
            //cout<<"gsk_score="<<gsk_score<<", Rk_u="<<Rk_u<<endl;
        }
        else{ //倘若resultFinal为空，说明没有任何元素比score_bound更优
            //cout<<"为空！"<<endl;
            Rk_u = score_bound;

        }

        //验证结果
        if(gsk_score >= Rk_u){
            usr_results.push_back(usr_id);
        }

        cout<<"验证候选用户"<<usr_id<<",gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;

    }
    verification_endTime = clock();
    cout<<"最终结果： 反top-k最近邻用户个数="<<usr_results.size()<<"，具体为："<<endl;
    printElements(usr_results);
    cout<<"Filter 阶段 用时："<< (double)(filter_endTime - filter_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;
    cout<<"Verification 阶段 用时："<< (double)(verification_endTime - verification_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;



}

TopkQueryCurrentResult TkGSKQ_NVD_vertexOnly(User user, int K, int a, double alpha, bool* poiMark, bool* poiADJMark){
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
                cout<<"social find related object: p<<"<<poi_id<<endl;
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
        int term = user.keywords[term_th];
        termid2idx[term] = term_th;  //
        priority_queue<nodeDistType_POI> Hi;
        int posting_size = getTermOjectInvListSize(term);
        //如果term为低频词汇
        if(posting_size<=posting_size_threshold){
            //termIsFrequent[term] = false;
            vector<int> posting_list = getTermOjectInvList(term);
            for(int o_id: posting_list){
                /*if(3413 == o_id){
                    cout<<"find poi is 3413!"<<endl;
                }*/
                POI related_poi = getPOIFromO2UOrgLeafData(o_id);
                float dist = usrToPOIDistance_phl(user,related_poi);
                nodeDistType_POI tmp1(o_id, dist, related_poi.keywords);
                Hi.push(tmp1);
                poiADJMark[o_id] = true;
            }
        }
        else{ //高频词汇,则将最近邻POI加入
            termIsFrequent[term] = true;
            int poi_id = getNNPOI_By_Vertex_Keyword(user.Ni,term);  //先得到距离Ni最近的带关键词term的poi
            int poi_id2 = getNNPOI_By_Vertex_Keyword(user.Nj,term);
#ifdef TRACKNVDTOPK
            if(poi_id==2633){
                cout<<"nvd initial find o2633"<<endl;
            }

#endif

            if(poi_id==poi_id2){
                POI poi = getPOIFromO2UOrgLeafData(poi_id);
                float dist = usrToPOIDistance_phl(user,poi);
                nodeDistType_POI tmp1(poi_id, dist, poi.keywords);
                Hi.push(tmp1);
                poiADJMark[poi_id] = true;
                //cout<<"H"<<term<<"初始化，加入p"<<poi_id<<endl;

            }else{

                if(true) {  //poiADJMark[poi_id] == false
                    POI poi1 = getPOIFromO2UOrgLeafData(poi_id);
                    float dist = usrToPOIDistance_phl(user,poi1);
                    nodeDistType_POI tmp1(poi_id, dist, poi1.keywords);
                    Hi.push(tmp1);
                    poiADJMark[poi_id] = true;
                    //cout<<"H"<<term<<"初始化，加入p"<<poi_id<<endl;
                }
                if(true) { //poiADJMark[poi_id2] == false
                    POI poi2 = getPOIFromO2UOrgLeafData(poi_id2);
                    float dist = usrToPOIDistance_phl(user,poi2);
                    nodeDistType_POI tmp2(poi_id2, dist, poi2.keywords);
                    Hi.push(tmp2);
                    poiADJMark[poi_id2] = true;
                    //cout<<"H"<<term<<"初始化，加入p"<<poi_id2<<endl;
                }

            }

        }
        H_list.push_back(Hi);
    }

    //1. 初始化全局Queue，对各个Hi的优先级进行计算并排序
    for(int term_th =0; term_th< user.keywords.size(); term_th++){
        int term_i = user.keywords[term_th];
        //priority_queue<nodeDistType_POI> Hi = ;
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
            //vector<int> poi_ADJ = getPOIAdj_NVD_By_Keyword(entry_p.poi_id,term_n);
            vector<int> poi_ADJ = getPOIAdj_NVD_By_Keyword_vertexOnly(entry_p.poi_id,term_n);

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



    for(int poi_id: userRelated_poiSet)
        poiMark[poi_id] = false;

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




//Filtering、verification阶段皆基于NVD,其中Filtering阶段直接对所有关键词相关user进行遍历，并求解各个user的lcl
vector<ResultDetail>RkGSKQ_NVD_Naive_optimized(int poi_id, int K, float a, float alpha) {

    cout<<"----------Running RkGSKQ_NVD_Naive_optimized---------------"<<endl;
    //double direction

    string phl_fileName = getRoadInputPath(PHL);//"../exp/indexes/LV.phl";

    //char phl_idxFileName[255] ; sprintf(phl_idxFileName, "../exp/indexes/%s.phl", road_map);

    phl.LoadLabel(phl_fileName.c_str());


    priority_queue<Result> resultFinal;

    clock_t algBeginTime, algEndTime;
    clock_t startTime, endTime;

    algBeginTime = clock();
    //获取poi信息；
    POI poi = getPOIFromO2UOrgLeafData(poi_id);
    printPOIInfo(poi);

    //获取与poi关键词相关的user的全集列表信息；
    set<int> poiRelated_usr;
    //候选验证用户
    vector<User> candidate_usr;
    vector<int> candidate_id;
    vector<double> candidate_gsk;
    vector<double> candidate_rk_current;
    //最终结果： 反top-k最近邻用户
    vector<ResultDetail> usr_results;
    vector<int> usrID_results;
    for(int term: poi.keywords){
        vector<int> term_related_user = getTermUserInvList(term);
        for(int usr: term_related_user)
            poiRelated_usr.insert(usr);
    }
    printf("总共有poiRelated_usr： %d个...\n", poiRelated_usr.size());
    //对poi关键词相关的user进行一一验证

    clock_t filter_startTime, filter_endTime;
    clock_t verification_startTime, verification_endTime;

    //startTime = clock();
    filter_startTime = clock();
    for(int usr_id: poiRelated_usr){
        //usr_id = 346;
        User user = getUserFromO2UOrgLeafData(usr_id);
#ifdef TRACK
        cout<<"评估poi文本相关用户：";
        printUsrInfo(user);
#endif

        startTime = clock();//计时开始
        double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,user);
        endTime = clock();//计时结束
        ////cout << "getGSKScore_o2u runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;
#ifdef TRACK
        cout<<"gsk_score="<<gsk_score;
#endif
        startTime = clock();//计时开始
        double u_lcl_rk = getUserLCL_NVD(user,K,a,alpha);
        endTime = clock();//计时结束
#ifdef TRACK
        cout<<",u_lcl_rk="<<u_lcl_rk<<endl; ////
        ////cout << "getUserLCL_NVD runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;
        cout << "-------------------------------------------------------" << endl;
        //cout<<",u_lcl_rk="<<u_lcl_rk<<endl;
#endif

        if(gsk_score > u_lcl_rk){
            candidate_usr.push_back(user);
            candidate_gsk.push_back(gsk_score);
            candidate_rk_current.push_back(u_lcl_rk);
            candidate_id.push_back(usr_id);
        }
        //getchar();

    }

    filter_endTime = clock();
    cout<<"候选用户个数="<<candidate_usr.size()<<"，具体为："<<endl;
    printElements(candidate_id);

    verification_startTime = clock();
    int candidate_size = candidate_usr.size();

    for(int i=0;i<candidate_size;i++){
        User u = candidate_usr[i];
        int u_id = u.id;
        double gsk_score = candidate_gsk[i];
        double score_bound = candidate_rk_current[i];
        double Rk_u = 0;


        TopkQueryCurrentResult topk_r6 = TkGSKQ_NVD(u, K, DEFAULT_A,alpha); //,poiMark,poiADJMark);
        Rk_u = topk_r6.topkScore;

        //验证结果
        if(gsk_score >= Rk_u){
            //usr_results.push_back(u_id);
            ResultDetail ur(u_id,-1,-1,-1,gsk_score,Rk_u);
            usr_results.push_back(ur);
            usrID_results.push_back(u_id);
        }
#ifdef TRACK
        cout<<"验证候选用户"<<u_id<<",gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
#endif


    }
    verification_endTime = clock();
    algEndTime = clock();


    cout<<"最终结果： 反top-k最近邻用户个数="<<usr_results.size()<<"，具体为："<<endl;
    //printElements(usrID_results);
    printPotentialUsers(usrID_results);
    cout<<"Filter 阶段 用时："<< (double)(filter_endTime - filter_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;
    cout<<"Verification 阶段 用时："<< (double)(verification_endTime - verification_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;
    cout << "---------------------RkGSKQ_NVD_userPure 用时：" << (double)(algEndTime - algBeginTime) / CLOCKS_PER_SEC * 1000 << "ms-----------" << endl;
    return  usr_results;


}



bool isCovered(int term, set<int>& keys){
    for(int _term: keys){
        if(term==_term)
            return true;
    }
    return false;
}


void RkGSKQ_bottom2up_by_NVD_noneReturn(int poi_id, int K, float a, float alpha) {

    cout<<"------------Running RkGSKQ_bottom2up_by_NVD!------------"<<endl;
    clock_t startTime, endTime;
    clock_t filter_startTime, while_startTime, filter_endTime;
    clock_t verification_startTime, verification_endTime;

    vector<int> RkGSKQ_results;


    //候选验证用户
    vector<User> candidate_usr;
    vector<int> candidate_id;
    vector<double> candidate_gsk;
    vector<double> candidate_rk_current;



    //user已被评估的标志向量
    bool* userMark = new bool[UserID_MaxKey];
    for(int i=0;i<UserID_MaxKey;i++)
        userMark[i] = false;

    int node_size = GTree.size();
    bool* nodeMark = new bool[node_size];
    for(int j=0;j<node_size;j++)
        nodeMark[j] = false;


    startTime = clock();
    filter_startTime = clock();

    //距离优先队列
    priority_queue<CheckEntry> Queue;
    //获取poi信息；
    POI poi = getPOIFromO2UOrgLeafData(poi_id);
    printPOIInfo(poi);
    int p_Ni = poi.Ni; int p_Nj = poi.Nj;
    vector<int> qo_keys = poi.keywords;
    int qo_leaf = Nodes[p_Ni].gtreepath.back();
    TreeNode tnData = getGIMTreeNodeData(qo_leaf,OnlyU);
#ifdef TRACK
    cout<<"o_leaf为leaf"<<qo_leaf<<endl;
    cout<<"o_leaf中user的关键词并集："<<endl;
    printSetElements(tnData.userUKeySet);
#endif


    //// 计算 poi 到 o_leaf上所有节点的border的距离
    unordered_map<int, vector<int>> itm; // intermediate answer, tree node -> array
    itm.clear();
    int cid, posa, posb, _min, dis;
    int o_leaf = Nodes[poi.Ni].gtreepath.back();
    itm[o_leaf].clear();
    int pos_Ni = id2idx_vertexInLeaf[p_Ni];
    int pos_Nj = id2idx_vertexInLeaf[p_Nj];
    int tn = o_leaf;
    for (int j = 0; j < GTree[tn].borders.size(); j++) {
        int dist_via_Ni = GTree[tn].mind[j * GTree[tn].leafnodes.size() + pos_Ni]+poi.dis;
        int dist_via_Nj = GTree[tn].mind[j * GTree[tn].leafnodes.size() + pos_Nj]+(poi.dist-poi.dis);
        int dist_to_border = min(dist_via_Ni, dist_via_Nj);
        itm[tn].push_back(dist_to_border);
    }

    ////获得与poi 社交文本相关的 user
    vector<int> checkin_users = poiCheckInIDList[poi_id];
    int user_textual_count = 0; int user_textual_useful_count = 0;
    for(int check_usr_id: checkin_users){
        vector<int> friends = friendshipMap[check_usr_id];
        for(int friend_id: friends){
            if(friend_id>=UserID_MaxKey) continue;
            if(userMark[friend_id] == true) continue;
            userMark[friend_id] = true;
            User u = getUserFromO2UOrgLeafData(friend_id);
            int usr_id = u.id;
            double relevance = 0;
            relevance = textRelevance(u.keywords, poi.keywords);
            //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
            if(relevance==0) continue;
            //该用户社交文本相关
            //cout<<"评估 社交文本相关用户：u"<<usr_id<<endl;
            user_textual_count ++;
            double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,u);
            double u_lcl_rk = getUserLCL_NVD(u,K,a,alpha);
            if(gsk_score > u_lcl_rk){
                candidate_usr.push_back(u);
                candidate_id.push_back(usr_id);
                candidate_gsk.push_back(gsk_score);
                candidate_rk_current.push_back(u_lcl_rk);
                user_textual_useful_count++;
            }

        }
    }
#ifdef TRACK
    cout<<"社交文本相关的 user总体个数="<<user_textual_count<<endl;
    cout<<"完成社交文本相关的 user的剪枝，此时 candidate user size="<<candidate_id.size()<<"，具体为："<<endl;
    printElements(candidate_id);
#endif
    //getchar();

    ////将q_o查询关键词分为 “（rare)稀有关键词”与“(none-rare)非稀有关键词” 两类,获得包含特稀疏keyword( |inv(term)|<K )的 user
    set<int> rareKeys;  set<int> nonRareKeys;
    set<int> rareKeywordRelated_usrSet;
    for(int term: qo_keys){
        int inv_olist_size = getTermOjectInvListSize(term);
        if(inv_olist_size>K){
            vector<int> usr_termRelated = getTermUserInvList(term);
            int inv_ulist_size = usr_termRelated.size();
            if(inv_ulist_size> uterm_SPARCITY_VALUE){
                nonRareKeys.insert(term);
                //printf("term %d 为用户高频词汇，u_frequency=%d\n", term,inv_ulist_size);
            }
            else{  //把一些带稀有关键词的user提前进行验证，以提升后续剪枝的效果
                for(int usr_id:usr_termRelated){
                    if(userMark[usr_id] == true) continue;
                    userMark[usr_id] = true;
                    User u = getUserFromO2UOrgLeafData(usr_id);
                    //int usr_id = u.id;
                    double relevance = 0;
                    relevance = textRelevance(u.keywords, poi.keywords);
                    if(relevance==0) continue;
                    double sk_score = getSKScore_o2u_phl(a,alpha,poi,u);
                    double u_lcl_rk = getUserLCL_NVD(u,K,a,alpha);
                    if(sk_score > u_lcl_rk){
                        candidate_usr.push_back(u);
                        candidate_id.push_back(usr_id);
                        candidate_gsk.push_back(sk_score);
                        candidate_rk_current.push_back(u_lcl_rk);
                        user_textual_useful_count++;
                    }
                }
            }

        }
        else{
            rareKeys.insert(term);
            vector<int> term_related_users = getTermUserInvList(term);
            for(int usr: term_related_users)
                rareKeywordRelated_usrSet.insert(usr);
        }

    }
    //getchar();
    ////对包含特稀疏关键词的user （这样的user通常较少）进行剪枝
    for(int usr_id:rareKeywordRelated_usrSet){
        if(userMark[usr_id] == true) continue;
        userMark[usr_id] = true;
        User u = getUserFromO2UOrgLeafData(usr_id);
        //int usr_id = u.id;
        double relevance = 0;
        relevance = textRelevance(u.keywords, poi.keywords);
        //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
        if(relevance==0) continue;
        //该用户社交文本相关
        //cout<<"评估 社交文本相关用户：u"<<usr_id<<endl;
        double sk_score = getSKScore_o2u_phl(a,alpha,poi,u);
        double u_lcl_rk = getUserLCL_NVD(u,K,a,alpha);
        if(sk_score > u_lcl_rk){
            candidate_usr.push_back(u);
            candidate_id.push_back(usr_id);
            candidate_gsk.push_back(sk_score);
            candidate_rk_current.push_back(u_lcl_rk);
            user_textual_useful_count++;
        }
    }

#ifdef TRACK
    cout<<"特稀疏关键词相关的 user 总体个数="<<rareKeywordRelated_usrSet.size()<<endl;
    cout<<"完成特稀疏关键词相关的 user的剪枝，此时 candidate user size="<<candidate_id.size()<<"，具体为："<<endl;
    printElements(candidate_id);
#endif

    set<int> poiRelated_usr;
    vector<int> usr_results;
    for(int term: nonRareKeys){
        vector<int> term_related_user = getTermUserInvList(term);
        for(int usr: term_related_user){
            if(userMark[usr] == true) continue;
            poiRelated_usr.insert(usr);
        }

    }
    clock_t occur_startTime, occur_endTime;
    printf("还余下与poi文本相关的用户： %d个， 在GTree中进行user occurant list构建...", poiRelated_usr.size());
    occur_startTime = clock();
    vector<set<int>> nodes_uOCCURList;
    vector<unordered_map<int,set<int>>> nodes_uCLUSTERList;  //vector_index= node_id,  term_id,
    for(int i=0;i<GTree.size();i++){
        set<int> node_occurList;
        nodes_uOCCURList.push_back(node_occurList);
        unordered_map<int,set<int>> clusters;
        nodes_uCLUSTERList.push_back(clusters);
    }
    for(int u_id: poiRelated_usr){
        User user = getUserFromO2UOrgLeafData(u_id);
        int u_Ni = user.Ni;
        int u_leaf = Nodes[u_Ni].gtreepath.back();
        nodes_uOCCURList[u_leaf].insert(u_id);
        vector<int> u_keywords = user.keywords;
        for(int term:u_keywords){
            if(isCovered(term,nonRareKeys)){
                for(int term2:u_keywords){
                    if(isCovered(term2,nonRareKeys)){
                        nodes_uCLUSTERList[u_leaf][term].insert(term2);
                    }

                }
            }

        }


        int father = GTree[u_leaf].father;
        while(father!=-1){
            nodes_uOCCURList[father].insert(u_id);
            for(int term:u_keywords){
                if(isCovered(term,nonRareKeys)){
                    for(int term2:u_keywords){
                        if(isCovered(term2,nonRareKeys)){
                            nodes_uCLUSTERList[father][term].insert(term2);
                        }

                    }

                }
            }

            father = GTree[father].father;
        }
    }
    occur_endTime = clock();
    cout<<"用时："<<(double)(occur_endTime-occur_startTime)/CLOCKS_PER_SEC*1000<<"ms!"<<endl;


    ////1. q_o 所在leaf 中的关键词相关 user先加入队列考察
    set<int> user_nearSet;
    for(int term: nonRareKeys){
        vector<int> u_list = getUsrTermRelatedEntry(term,qo_leaf);
        if(u_list.size()>0){
            for(int u_id: u_list){
                user_nearSet.insert(u_id);
            }
        }
    }

#ifdef TRACK
    if(user_nearSet.size() != 0){
        cout<<"q_o 所在leaf"<<qo_leaf<<"中的关键词相关 user先加入队列考察"<<endl;
    } else{
        cout<<"q_o 所在leaf "<<qo_leaf<<"中的无关键词相关user"<<endl;
    }
#endif

    for(int usr_id: user_nearSet){
        if(userMark[usr_id] == true) continue;
        User user = getUserFromO2UOrgLeafData(usr_id);
        double dist = usrToPOIDistance_phl(user, poi);
        set<int> u_keySet;
        /*for(int term: user.keywords)
            u_keySet.insert(term);*/
        CheckEntry entry(usr_id,qo_leaf,dist,user);
        Queue.push(entry);
    }
    //设定当前level最高的Gtree Node
    int node_highest = qo_leaf;
    int root = 0;

    ////2.从qo_leaf出发按距离优先,进行路网拓展（在Gtree中从qo_leaf开始，自顶向上，遍历各个关键词相关节点）含有相关关键词的 user,user_leaf, user_node

    //map<int,map<int,bool>> borderEntryVist;
    bool** borderEntryVist = new bool*[VertexNum]; //border_id, node_id
    for(int i=0;i<VertexNum;i++){
        int node_size = GTree.size();
        borderEntryVist[i] = new bool [node_size];
        for(int j=0;j<node_size;j++){
            borderEntryVist[i][j] = false;
        }
    }
    bool* userEntryVist = new bool [UserID_MaxKey];
    for(int i=0;i<UserID_MaxKey;i++){
        userEntryVist[i]= false;
    }

    while_startTime = clock();

    while(Queue.size()>0 ||node_highest!= root||node_highest!= -1){
        if(Queue.empty()){ //当前范围内已全部遍历， 需扩大路网遍历范围
#ifdef TRACK_FILTER
            cout<<"当前范围内已全部遍历,"<<"当前最高层节点node"<<node_highest<<"， 需更新扩大路网图的遍历范围"<<endl;
#endif
            //Extend_Range(node_highest,itm,Queue,nonRareKeys,poi, K);
            Extend_Range(node_highest,itm,Queue,nonRareKeys,poi, K,a,alpha);
            if(Queue.empty()) break;  //路网拓展可被终止

        }
        CheckEntry topEntry = Queue.top();
        Queue.pop();

        float dist_parentBorder = topEntry.dist;

        if(topEntry.isUser){ //// an entry for user
            User user = topEntry.user;
            int usr_id = user.id;
            //if(userEntryVist.count(usr_id)) continue;
            if(userMark[usr_id] == true) continue;
            userMark[usr_id] = true;

            double u_SB_rk = getUserLCL_NVD(user,K,a,alpha);

            ////cout << "getUserLCL_NVD runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;

            double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,user);
            ////cout << "getGSKScore_o2u runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;
#ifdef TRACK_FILTER
            cout<<"------DEQUEUE: user entry: u"<<usr_id<<"------"<<endl;

            cout<<"u_SB_rk="<<u_SB_rk<<",gsk_score="<<gsk_score<<endl;
#endif

            if(gsk_score > u_SB_rk){
                candidate_usr.push_back(user);
                candidate_id.push_back(usr_id);
                candidate_gsk.push_back(gsk_score);
                candidate_rk_current.push_back(u_SB_rk);
            }
        }
        else{  //// an entry for border

            int border_id = topEntry.id;
            int node_id =topEntry.node_id;
            //若该节点已经被aggresive prune了，则discard
            if(nodeMark[node_id]==true)
                continue;
            //判断是否为user稀疏节点， 若是则aggresive prune
#ifdef   userspacity
            if(nodes_uOCCURList[node_id].size()< NODE_SPARCITY_VALUE ){//&& dist_parentBorder>LONG_TERM
                for(int usr_id: nodes_uOCCURList[node_id]){
                    if(userMark[usr_id] != true){
                        User user = getUserFromO2UOrgLeafData(usr_id);
                        double dist = usrToPOIDistance_phl(user, poi);
                        set<int> u_keySet;
                        CheckEntry entry(usr_id,node_id,dist,user);
                        Queue.push(entry);
                    }
                }
                nodeMark[node_id] =true;
                continue;
            }
#endif


            if(borderEntryVist[border_id][node_id] == true) {
                continue;
            }
            borderEntryVist[border_id][node_id] = true;

#ifdef  TRACK_FILTER
            if(GTree[node_id].isleaf)
                cout<<"----DEQUEUE: border entry: b"<<border_id<<" in leaf"<<node_id<<" 出列！----"<<endl;
            else
                cout<<"----DEQUEUE: border entry: b"<<border_id<<" in node"<<node_id<<" 出列！----"<<endl;
#endif

            set<int> Keys;


            //Keys = getBorder_SB_NVD(topEntry,K,a, alpha);
            Keys = getBorder_SB_NVD_Cluster(topEntry,K,a, alpha,nodes_uCLUSTERList);

            if(Keys.size() ==0) {
#ifdef TRACK_FILTER
                cout<<"所有其后继被 ****prune！"<<endl;
#endif
                continue;
            }
#ifdef TRACK_FILTER
            cout<<"remained (user) Keywords"; printSetElements(Keys);
            cout<<"其后继**不能**被prune！"<<endl;
#endif

            //若为叶节点，加入叶节点下的各个user
            if(GTree[node_id].isleaf){
#ifdef TRACK
                cout<<"b"<<border_id<<", 位于leaf"<<node_id<<endl;
#endif
                set<int> user_within;
                for(int term: Keys){
                    vector<int> u_list = getUsrTermRelatedEntry(term,node_id);
                    if(u_list.size()>0){
                        for(int u_id: u_list){
                            user_within.insert(u_id);
                        }
                    }
                }
                //加入叶节点下的各个user
                for(int usr_id: user_within){
                    if(userMark[usr_id] == true) continue;
                    User user = getUserFromO2UOrgLeafData(usr_id);
                    double dist = usrToPOIDistance_phl(user, poi);
                    set<int> u_keySet;
                    CheckEntry entry(usr_id,node_id,dist,user);
                    Queue.push(entry);
#ifdef TRACK
                    cout<<"在main中加入user entry:"; entry.printRlt();
                    cout<<"Queue.size="<<Queue.size()<<endl;
#endif
                }
            }
            else{
                //获得父节点下所有关键词相关的孩子节点(line 22-25)
#ifdef TRACK
                cout<<"b"<<border_id<<", 位于node"<<node_id<<endl;
#endif

                set<int> child_set;
                for(int term:Keys){
                    vector<int> related_childs = getUsrTermRelatedEntry(term,node_id);
                    if(related_childs.size()>0){
                        for(int cn_id: related_childs){
                            child_set.insert(cn_id);
                        }
                    }
                }
#ifdef TRACK
                cout<<"关键词相关的后继："; printSetElements(child_set);
#endif

                int parent_border_th = topEntry.b_th;
                int posa = parent_border_th;
                for(int cn_id: child_set){
                    for(int b_th=0;b_th<GTree[cn_id].borders.size();b_th++){
                        int child_border_id = GTree[cn_id].borders[b_th];
                        int posb = GTree[cn_id].up_pos[b_th];
                        float b2b_dist = GTree[node_id].mind[posb * GTree[node_id].union_borders.size() + posa];

                        float dist_childBorder2 = topEntry.dist + b2b_dist;
                        float dist_childBorder = getDistance_phl(child_border_id,poi);  //jins

                        CheckEntry entry(child_border_id,cn_id,Keys,dist_childBorder,b_th);
                        Queue.push(entry);
#ifdef TRACK
                        cout<<"在main中加入border entry:"; entry.printRlt();
                        cout<<"Queue.size="<<Queue.size()<<endl;
#endif

                    }
                }

            }
        }

    }
    delete [] userMark; delete [] nodeMark;
    delete [] userEntryVist;
    for(int i=0;i<VertexNum;i++){
        delete [] borderEntryVist[i];
    }
    delete [] borderEntryVist;

    filter_endTime = clock();


    ////对不能被减枝的user, 进行最终TkGSKQ求解

#ifdef TRACK
    cout<<"候选用户个数="<<candidate_usr.size()<<"，具体为："<<endl;
    printElements(candidate_id);
#else
    cout<<"候选用户个数="<<candidate_usr.size()<<endl;
#endif

    //nvd_filtering
    verification_startTime = clock();
    int candidate_size = candidate_usr.size();
    for(int i=0;i<candidate_size;i++){
        User u = candidate_usr[i];
        int u_id = u.id;
        double gsk_score = candidate_gsk[i];
        double score_bound = candidate_rk_current[i];
        double Rk_u = 0;

        if(true){
            TopkQueryCurrentResult topk_r6 = TkGSKQ_NVD(u, K, DEFAULT_A,alpha);
            Rk_u = topk_r6.topkScore;

        }
        else{
            priority_queue<Result> results = TkGSKQ_bottom2up_verify_disk(u_id,K,DEFAULT_A,alpha,gsk_score,score_bound);  //1000,0.0);

            if(results.size()==K) {
                Rk_u = results.top().score;
            }
            else{ //倘若resultFinal为空，说明没有任何元素比score_bound更优
                Rk_u = score_bound;
            }

        }

        //验证结果
        if(gsk_score >= Rk_u){
            RkGSKQ_results.push_back(u_id);
        }
#ifdef TRACK
        cout<<"验证候选用户"<<u_id<<",gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
#endif

    }
    verification_endTime = clock();
    endTime = clock();
    cout<<"Filter 阶段 用时："<< (double)(filter_endTime - filter_startTime) / CLOCKS_PER_SEC*1000  << "ms, （其中while loop用时 "
        << (double)(filter_endTime - while_startTime) / CLOCKS_PER_SEC*1000 <<"ms)" << endl;
    cout<<"Verification 阶段 用时："<< (double)(verification_endTime - verification_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;
    cout<<"最终结果： 反top-k最近邻用户个数="<<RkGSKQ_results.size()<<"，具体为："<<endl;
    printPotentialUsers(RkGSKQ_results);

    //cout<<"user_textual_count="<<user_textual_count<<",user_textual_useful_count="<<user_textual_useful_count<<endl;
    cout<<"---------RkGSKQ_bottom2up_by_NVD *p"<<poi.id<<"* 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC*1000  << "ms-----------" << endl;
}


SingleResults RkGSKQ_bottom2up_by_NVD(POI& poi, int K, float a, float alpha) {  //gangcai

    cout<<"------------Running Single RkGSKQ bottom2up based on IG-tree&NVD!------------"<<endl;
    clock_t startTime, endTime;
    clock_t filter_startTime,filter_endTime;
    clock_t verification_startTime, verification_endTime;

    vector<int> RkGSKQ_results;


    //候选验证用户
    vector<User> candidate_usr;
    vector<int> candidate_id;
    vector<double> candidate_gsk;
    vector<double> candidate_rk_current;

    //user已被评估的标志向量
    bool userMark[UserID_MaxKey];
    memset(userMark,false, sizeof(userMark));


    startTime = clock();

    //距离优先队列
    priority_queue<CheckEntry> Queue;
    //获取poi信息；
    int poi_id = poi.id;
    printPOIInfo(poi);
    int p_Ni = poi.Ni; int p_Nj = poi.Nj;
    vector<int> qo_keys = poi.keywords;
    int qo_leaf = Nodes[p_Ni].gtreepath.back();
    TreeNode tnData = getGIMTreeNodeData(qo_leaf,OnlyU);
#ifdef TRACK
    cout<<"o_leaf为leaf"<<qo_leaf<<endl;
    cout<<"o_leaf中user的关键词并集："<<endl;
    printSetElements(tnData.userUKeySet);
#endif


    //// 计算 poi 到 o_leaf上所有节点的border的距离
    unordered_map<int, vector<int>> itm; // intermediate answer, tree node -> array
    itm.clear();
    int cid, posa, posb, _min, dis;
    int o_leaf = Nodes[poi.Ni].gtreepath.back();
    itm[o_leaf].clear();
    int pos_Ni = id2idx_vertexInLeaf[p_Ni];
    int pos_Nj = id2idx_vertexInLeaf[p_Nj];
    int tn = o_leaf;
    for (int j = 0; j < GTree[tn].borders.size(); j++) {
        int dist_via_Ni = GTree[tn].mind[j * GTree[tn].leafnodes.size() + pos_Ni]+poi.dis;
        int dist_via_Nj = GTree[tn].mind[j * GTree[tn].leafnodes.size() + pos_Nj]+(poi.dist-poi.dis);
        int dist_to_border = min(dist_via_Ni, dist_via_Nj);
        itm[tn].push_back(dist_to_border);
    }

    ////获得与poi 社交文本相关的 user
    vector<int> checkin_users = poiCheckInIDList[poi_id];
    int user_textual_count = 0; int user_textual_useful_count = 0;
    for(int check_usr_id: checkin_users){
        vector<int> friends = friendshipMap[check_usr_id];
        for(int friend_id: friends){
            if(friend_id>=UserID_MaxKey) continue;
            if(userMark[friend_id] == true) continue;
            userMark[friend_id] = true;
            User u = getUserFromO2UOrgLeafData(friend_id);
            int usr_id = u.id;
            double relevance = 0;
            relevance = textRelevance(u.keywords, poi.keywords);
            //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
            if(relevance==0) continue;
            //该用户社交文本相关
            //cout<<"评估 社交文本相关用户：u"<<usr_id<<endl;
            user_textual_count ++;
            double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,u);
            double u_lcl_rk = getUserLCL_NVD(u,K,a,alpha);
            if(gsk_score > u_lcl_rk){
                candidate_usr.push_back(u);
                candidate_id.push_back(usr_id);
                candidate_gsk.push_back(gsk_score);
                candidate_rk_current.push_back(u_lcl_rk);
                user_textual_useful_count++;
            }

        }
    }
#ifdef TRACK
    cout<<"社交文本相关的 user总体个数="<<user_textual_count<<endl;
    cout<<"完成社交文本相关的 user的剪枝，此时 candidate user size="<<candidate_id.size()<<"，具体为："<<endl;
    printElements(candidate_id);
#endif
    //getchar();

    ////将q_o查询关键词分为 “（rare)稀有关键词”与“(none-rare)非稀有关键词” 两类,获得包含特稀疏keyword( |inv(term)|<K )的 user
    set<int> rareKeys;  set<int> nonRareKeys;
    set<int> rareKeywordRelated_usrSet;
    for(int term: qo_keys){
        int inv_olist_size = getTermOjectInvListSize(term);
        if(inv_olist_size>K){
            nonRareKeys.insert(term);

        }
        else{
            rareKeys.insert(term);
            vector<int> term_related_users = getTermUserInvList(term);
            for(int usr: term_related_users)
                rareKeywordRelated_usrSet.insert(usr);
        }

    }
    ////对包含特稀疏关键词的user （这样的user通常较少）进行剪枝
    for(int usr_id:rareKeywordRelated_usrSet){
        if(userMark[usr_id] == true) continue;
        userMark[usr_id] = true;
        User u = getUserFromO2UOrgLeafData(usr_id);
        //int usr_id = u.id;
        double relevance = 0;
        relevance = textRelevance(u.keywords, poi.keywords);
        //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
        if(relevance==0) continue;
        //该用户社交文本相关
        //cout<<"评估 社交文本相关用户：u"<<usr_id<<endl;
        double sk_score = getSKScore_o2u_phl(a,alpha,poi,u);
        double u_lcl_rk = getUserLCL_NVD(u,K,a,alpha);
        if(sk_score > u_lcl_rk){
            candidate_usr.push_back(u);
            candidate_id.push_back(usr_id);
            candidate_gsk.push_back(sk_score);
            candidate_rk_current.push_back(u_lcl_rk);
            user_textual_useful_count++;
        }
    }

#ifdef TRACK
    cout<<"特稀疏关键词相关的 user 总体个数="<<rareKeywordRelated_usrSet.size()<<endl;
    cout<<"完成特稀疏关键词相关的 user的剪枝，此时 candidate user size="<<candidate_id.size()<<"，具体为："<<endl;
    printElements(candidate_id);
#endif
    //getchar();

    ////1. q_o 所在leaf 中的关键词相关 user先加入队列考察
    set<int> user_nearSet;
    for(int term: nonRareKeys){
        vector<int> u_list = getUsrTermRelatedEntry(term,qo_leaf);
        if(u_list.size()>0){
            for(int u_id: u_list){
                user_nearSet.insert(u_id);
            }
        }
    }

#ifdef TRACK
    if(user_nearSet.size() != 0){
        cout<<"q_o 所在leaf"<<qo_leaf<<"中的关键词相关 user先加入队列考察"<<endl;
    } else{
        cout<<"q_o 所在leaf "<<qo_leaf<<"中的无关键词相关user"<<endl;
    }
#endif

    for(int usr_id: user_nearSet){
        if(userMark[usr_id] == true) continue;
        User user = getUserFromO2UOrgLeafData(usr_id);
        double dist = usrToPOIDistance_phl(user, poi);
        set<int> u_keySet;
        /*for(int term: user.keywords)
            u_keySet.insert(term);*/
        CheckEntry entry(usr_id,qo_leaf,dist,user);
        Queue.push(entry);
    }
    //设定当前level最高的Gtree Node
    int node_highest = qo_leaf;
    int root = 0;

    ////2.从qo_leaf出发按距离优先,进行路网拓展（在Gtree中从qo_leaf开始，自顶向上，遍历各个关键词相关节点）含有相关关键词的 user,user_leaf, user_node

    map<int,map<int,bool>> borderEntryVist;  //border_id, node_id
    map<int,bool> userEntryVist;
    filter_startTime = clock();
    while(Queue.size()>0 ||node_highest!= root){
        if(Queue.empty()){ //当前范围内已全部遍历， 需扩大路网遍历范围
#ifdef TRACK
            cout<<"当前范围内已全部遍历,"<<"当前最高层节点node"<<node_highest<<"， 需更新扩大路网图的遍历范围"<<endl;
#endif
            //Extend_Range(node_highest,itm,Queue,nonRareKeys,poi, K);
            Extend_Range(node_highest,itm,Queue,nonRareKeys,poi, K,a,alpha);
            if(Queue.empty()) break;  //路网拓展可被终止

        }
        CheckEntry topEntry = Queue.top();
        Queue.pop();

        float dist_parentBorder = topEntry.dist;

        if(topEntry.isUser){ //// an entry for user
            User user = topEntry.user;
            int usr_id = user.id;
            //if(userEntryVist.count(usr_id)) continue;
            if(userMark[usr_id] == true) continue;
            userMark[usr_id] = true;

            double u_SB_rk = getUserLCL_NVD(user,K,a,alpha);

            ////cout << "getUserLCL_NVD runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;

            double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,user);
            ////cout << "getGSKScore_o2u runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;
#ifdef TRACK
            cout<<"------DEQUEUE: user entry: u"<<usr_id<<"------"<<endl;

            cout<<"u_SB_rk="<<u_SB_rk<<endl;
            cout<<"gsk_score="<<gsk_score<<endl;
#endif

            if(gsk_score > u_SB_rk){
                candidate_usr.push_back(user);
                candidate_id.push_back(usr_id);
                candidate_gsk.push_back(gsk_score);
                candidate_rk_current.push_back(u_SB_rk);
            }
        }
        else{  //// an entry for border
            int border_id = topEntry.id;
            int node_id =topEntry.node_id;
            if(borderEntryVist.count(border_id)) {
                if(borderEntryVist[border_id].count(node_id))
                    continue;
            }
            borderEntryVist[border_id][node_id] = true;

#ifdef  TRACK
            if(GTree[node_id].isleaf)
                cout<<"----DEQUEUE: border entry: b"<<border_id<<" in leaf"<<node_id<<" 出列！----"<<endl;
            else
                cout<<"----DEQUEUE: border entry: b"<<border_id<<" in node"<<node_id<<" 出列！----"<<endl;
#endif

            set<int> Keys;


            //CheckEntry border_entry,int K, int a, double alpha
            Keys = getBorder_SB_NVD(topEntry,K,a, alpha);
            //Keys = getBorder_SB_NVD(topEntry,poi,nonRareKeys,K,a, alpha);

            //cout<<"remained (user) Keywords"; printSetElements(Keys);
            if(Keys.size() ==0) {
#ifdef TRACK
                cout<<"所有其后继被 ****prune！"<<endl;
#endif
                continue;
            }
#ifdef TRACK
            cout<<"其后继**不能**被prune！"<<endl;
#endif

            //若为叶节点，加入叶节点下的各个user
            if(GTree[node_id].isleaf){
#ifdef TRACK
                cout<<"b"<<border_id<<", 位于leaf"<<node_id<<endl;
#endif
                set<int> user_within;
                for(int term: Keys){
                    vector<int> u_list = getUsrTermRelatedEntry(term,node_id);
                    if(u_list.size()>0){
                        for(int u_id: u_list){
                            user_within.insert(u_id);
                        }
                    }
                }
                //加入叶节点下的各个user
                for(int usr_id: user_within){
                    if(userMark[usr_id] == true) continue;
                    User user = getUserFromO2UOrgLeafData(usr_id);
                    double dist = usrToPOIDistance_phl(user, poi);
                    set<int> u_keySet;
                    CheckEntry entry(usr_id,node_id,dist,user);
                    Queue.push(entry);
#ifdef TRACK
                    cout<<"在main中加入user entry:"; entry.printRlt();
                    cout<<"Queue.size="<<Queue.size()<<endl;
#endif
                }
            }
            else{
                //获得父节点下所有关键词相关的孩子节点(line 22-25)
#ifdef TRACK
                cout<<"b"<<border_id<<", 位于node"<<node_id<<endl;
#endif
                set<int> child_set;
                for(int term:Keys){
                    vector<int> related_childs = getUsrTermRelatedEntry(term,node_id);
                    if(related_childs.size()>0){
                        for(int cn_id: related_childs){
                            child_set.insert(cn_id);
                        }
                    }
                }
#ifdef TRACK
                cout<<"关键词相关的后继："; printSetElements(child_set);
#endif

                int parent_border_th = topEntry.b_th;
                int posa = parent_border_th;
                for(int cn_id: child_set){
                    for(int b_th=0;b_th<GTree[cn_id].borders.size();b_th++){
                        int child_border_id = GTree[cn_id].borders[b_th];
                        int posb = GTree[cn_id].up_pos[b_th];
                        float b2b_dist = GTree[node_id].mind[posb * GTree[node_id].union_borders.size() + posa];

                        float dist_childBorder2 = topEntry.dist + b2b_dist;
                        float dist_childBorder = getDistance_phl(child_border_id,poi);  //jins

                        CheckEntry entry(child_border_id,cn_id,Keys,dist_childBorder,b_th);
                        Queue.push(entry);
#ifdef TRACK
                        cout<<"在main中加入border entry:"; entry.printRlt();
                        cout<<"Queue.size="<<Queue.size()<<endl;
#endif

                    }
                }

            }
        }

    }
    filter_endTime = clock();


    ////对不能被减枝的user, 进行最终TkGSKQ求解

#ifdef TRACK
    cout<<"候选用户个数="<<candidate_usr.size()<<"，具体为："<<endl;
    printElements(candidate_id);
#endif

    verification_startTime = clock();
    int candidate_size = candidate_usr.size();
    vector<ResultDetail> potential_customerDetails;

////yanzsingle
    for(int i=0;i<candidate_size;i++){
        User u = candidate_usr[i];
        int u_id = u.id;
        //cout<<"将验证候选用户"<<u_id<<endl;
        double gsk_score = candidate_gsk[i];
        //cout<<"gsk_score="<<gsk_score<<endl;
        double score_bound = candidate_rk_current[i];
        //cout<<"score_bound="<<score_bound<<endl;
        double Rk_u = 0;


        if(false){  //joke
            TopkQueryCurrentResult topk_r6 = TkGSKQ_NVD(u, K, DEFAULT_A,alpha);//,poiMark,poiADJMark);
            Rk_u = topk_r6.topkScore;
        }
        else{
            priority_queue<Result> results = TkGSKQ_bottom2up_verify_disk(u_id,K,DEFAULT_A,alpha,gsk_score,score_bound);//1000,0.0);

            if(results.size()==K){
                Rk_u= results.top().score;     //gangcai
            }
            else{
                Rk_u = score_bound;
            }
        }



        //验证结果
        if(gsk_score >= Rk_u){
            RkGSKQ_results.push_back(u_id);
            //加入detail细节
            ResultDetail rd(u_id,0,0,0, gsk_score,Rk_u);
            potential_customerDetails.push_back(rd);

        }
#ifdef TRACK
        cout<<"验证候选用户"<<u_id<<",gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
#endif
        //cout<<"验证候选用户"<<u_id<<",gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;

    }
    verification_endTime = clock();
    endTime = clock();
    SingleResults signle_results =  resultsAnalysis_poi(potential_customerDetails,poi, K,a,alpha);

    cout<<"最终结果： 反top-k最近邻用户个数="<<RkGSKQ_results.size()<<"，具体为："<<endl;
    printPotentialUsers(RkGSKQ_results);
    cout<<"Filter 阶段 用时："<< (double)(filter_endTime - filter_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;
    cout<<"Verification 阶段 用时："<< (double)(verification_endTime - verification_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;

    //cout<<"user_textual_count="<<user_textual_count<<",user_textual_useful_count="<<user_textual_useful_count<<endl;
    cout<<"---------RkGSKQ_bottom2up_by_NVD 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC*1000  << "ms-----------" << endl;
    return signle_results;
}


//基于NVD的Batch Reverse处理算法



BatchResults RkGSKQ_Batch_NVD_RealData(vector<int> stores, int K, float a, float alpha) {  //gangcai

    cout<<"------------Running Batch RkGSKQ using NVD&IG-Tree!------------"<<endl;
    //评估 社交文本相关用户：u22557 for poi13448, user信息：
    User _tmp_user = getUserFromO2UOrgLeafData(22557); printUsrInfo(_tmp_user);
    POI _tmp_poi = getPOIFromO2UOrgLeafData(13448); printPOIInfo(_tmp_poi);
    float relevance = textRelevance(_tmp_user.keywords, _tmp_poi.keywords);
    //getchar();

    clock_t startTime, endTime;
    clock_t initial_startTime, initial_endTime;
    clock_t filter_startTime,while_startTime, filter_endTime;
    clock_t verification_startTime, verification_endTime;

    //vector<int> RkGSKQ_results;
    BatchResults batch_results;



    unordered_map<int, POI> candidate_POIs;
    unordered_map<int, int> poi_Id2IdxMap;
    unordered_map<int, set<int>>  nonRareKeywordMap;


    //候选验证用户
    set<BatchVerifyEntry> candidate_User;  //candidate users
    map<int, priority_queue<ResultCurrent>> candidate_usr_related_store;   // user_id, {related poi1,2,....}



    //border已被评估的标志向量
    int node_size = GTree.size();
    int inner_size = node_size*VertexNum;  // func(bEntry_key)= v_id+v_num*node_id + v_num*node_num*poi_th
    int bEntry_MaxKey =  inner_size * stores.size();
    unordered_map<int,bool> borderMark;


    //user已被评估的标志向量表
    int upEntry_MaxKey =  UserID_MaxKey * stores.size();
    unordered_map<int,bool> upEntryMark;
    unordered_map<int,bool> userMark;
    //标记gtree node、 poi 对的标志向量
    int np_maxKey = GTree.size()*stores.size();
    bool nodepMark[np_maxKey];
    memset(nodepMark,false, sizeof(nodepMark));


    startTime = clock();
    filter_startTime = clock();
    initial_startTime = clock();

    //距离优先队列

    ////初始化gtree node的 uOCCURList以及 uCLUSTERList
    unordered_map<int,UOCCURList> pnodes_uOCCURListMap;
    unordered_map<int,UCLUSTERList> pnodes_uCLUSTERListMap;



    priority_queue<BatchCheckEntry> Queue;
    priority_queue<POIHighestNode> T_max;


    int root = 0;
    //// 1. 获取各个poi信息, 并计算各个 poi 到 o_leaf上所有节点的border的距离
    unordered_map<int, unordered_map<int, vector<int>>> global_itm;
    int poi_th = 0;
    set<int> midKeywordRelated_usrSet;
    for(int poi_id: stores){
        poi_Id2IdxMap[poi_id] = poi_th;

        if(poi_id==256){
            cout<<"find poi 256!"<<endl;
        }

        POI poi = getPOIFromO2UOrgLeafData(poi_id);

        candidate_POIs[poi_id] = poi;

        printPOIInfo(poi);


        //初始化并构建该poi对应的 nodes_uOCCURList 及 nodes_uCLUSTERList
        UOCCURList nodes_uOCCURList;
        UCLUSTERList nodes_uCLUSTERList;  //vector_index= node_id,  term_id,
        for(int i=0;i<GTree.size();i++){
            set<int> node_occurList;
            nodes_uOCCURList.push_back(node_occurList);
            unordered_map<int,set<int>> clusters;
            nodes_uCLUSTERList.push_back(clusters);
        }
        pnodes_uOCCURListMap[poi_id] = nodes_uOCCURList;
        pnodes_uCLUSTERListMap[poi_id] = nodes_uCLUSTERList;


        int p_Ni = poi.Ni; int p_Nj = poi.Nj;
        vector<int> qo_keys = poi.keywords;
        int qo_leaf = Nodes[p_Ni].gtreepath.back();
        TreeNode tnData = getGIMTreeNodeData(qo_leaf,OnlyU);
#ifdef TRACK
        cout<<"o_leaf为leaf"<<qo_leaf<<endl;
        cout<<"o_leaf中user的关键词并集："<<endl;
        printSetElements(tnData.userUKeySet);
#endif

        //// 计算 各个 poi 到 o_leaf上所有节点的border的距离
        MaterilizedDistanceMap itm; // intermediate answer, tree node -> array
        itm.clear();
        int cid, posa, posb, _min, dis;
        int o_leaf = Nodes[poi.Ni].gtreepath.back();
        itm[o_leaf].clear();
        int pos_Ni = id2idx_vertexInLeaf[p_Ni];
        int pos_Nj = id2idx_vertexInLeaf[p_Nj];
        int tn = o_leaf;
        for (int j = 0; j < GTree[tn].borders.size(); j++) {
            int dist_via_Ni = GTree[tn].mind[j * GTree[tn].leafnodes.size() + pos_Ni]+poi.dis;
            int dist_via_Nj = GTree[tn].mind[j * GTree[tn].leafnodes.size() + pos_Nj]+(poi.dist-poi.dis);
            int dist_to_border = min(dist_via_Ni, dist_via_Nj);
            itm[tn].push_back(dist_to_border);
        }
        global_itm[poi_id] = itm;

        ////获得与poi 社交文本相关的 user
        vector<int> checkin_users = poiCheckInIDList[poi_id];
        int user_textual_count = 0; int user_textual_useful_count = 0;
        for(int check_usr_id: checkin_users){
            vector<int> friends = friendshipMap[check_usr_id];
            for(int friend_id: friends){
                if(friend_id>=UserID_MaxKey) continue;
                int up_key = friend_id + UserID_MaxKey* poi_th;   //重要：每个poi对应（0,usrnum-1)

                if(upEntryMark[up_key] == true) continue;
                upEntryMark[up_key] = true;


                User u = getUserFromO2UOrgLeafData(friend_id);
                int usr_id = u.id;
                double relevance = 0;
                relevance = textRelevance(u.keywords, poi.keywords);
                //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
                if(relevance==0) continue;
                //该用户社交文本相关
#ifdef TRACKNVDBATCH
                cout<<"评估 社交文本相关用户1`：u"<<usr_id<<endl;
#endif
                user_textual_count ++;
                double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,u);
                double u_lcl_rk = getUserLCL_NVD(u,K,a,alpha);
                if(gsk_score > u_lcl_rk){

                    //更新 user-poi associating list
                    candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                    if(candidate_usr_related_store[usr_id].size()>K){
                        candidate_usr_related_store[usr_id].pop();
                    }
                    //加入candidate user set
                    BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_lcl_rk,u);
                    candidate_User.insert(ve);

                    user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                    cout<<"social_textual find candidate 用户： u"<<usr_id<<endl;
#endif
                }

            }
        }
#ifdef TRACKNVDBATCH
        cout<<"社交文本相关的 user总体个数="<<user_textual_count<<endl;
    //cout<<"完成社交文本相关的 user的剪枝，此时 candidate user size="<<candidate_id.size()<<"，具体为："<<endl;
    //printElements(candidate_id);
#endif
        //getchar();

        ////将q_o查询关键词分为 “（rare)稀有关键词”与“(none-rare)非稀有关键词” 两类,获得包含特稀疏keyword( |inv(term)|<K )的 user
        set<int> rareKeys;  set<int> nonRareKeys;
        set<int> rareKeywordRelated_usrSet;
        for(int term: qo_keys){
            int inv_olist_size = getTermOjectInvListSize(term);
            if(inv_olist_size>K){
                vector<int> usr_termRelated = getTermUserInvList(term);
                int inv_ulist_size = usr_termRelated.size();
                if(inv_ulist_size> uterm_SPARCITY_VALUE){
                    nonRareKeys.insert(term);
                    printf("term %d 为用户高频词汇，u_frequency=%d\n", term,inv_ulist_size);
                }
                else{ ////把一些带稀有关键词的user提前进行验证，以提升后续剪枝的效果
                    for(int user_id: usr_termRelated){
#ifdef TRACKNVDBATCH
                        cout<<"评估 社交文本相关用户2：u"<<user_id<<" for poi"<<poi.id<<", user信息：";
#endif
                        if(user_id>=UserID_MaxKey) continue;
                        int up_key = user_id + UserID_MaxKey* poi_th;   //重要：每个poi对应（0,usrnum-1)

                        if(upEntryMark[up_key] == true) continue;
                        upEntryMark[up_key] = true;
                        //midKeywordRelated_usrSet.insert(user_id);

                        User u = getUserFromO2UOrgLeafData(user_id);
#ifdef TRACKNVDBATCH
                        printUsrInfo(u);
#endif

                        int usr_id = u.id;
                        double relevance = 0;
                        relevance = textRelevance(u.keywords, poi.keywords);
                        //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
                        if(relevance==0) continue;
                        //该用户社交文本相关
#ifdef TRACKNVDBATCH
                        cout<<"relevance="<<relevance;
#endif

                        user_textual_count ++;
                        double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,u);
#ifdef TRACKNVDBATCH
                        printf("gsk_score(u%d, p%d)=%f, ", u.id,poi.id, gsk_score);
#endif
                        double u_lcl_rk = getUserLCL_NVD(u,K,a,alpha);
#ifdef TRACKNVDBATCH
                        printf(" u%d.lcl_rk=%f\n", u.id, u_lcl_rk);
#endif
                        if(gsk_score > u_lcl_rk){

                            //更新 user-poi associating list
                            candidate_usr_related_store[user_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                            if(candidate_usr_related_store[user_id].size()>K){
                                candidate_usr_related_store[user_id].pop();
                            }
                            //加入candidate user set
                            BatchVerifyEntry ve = BatchVerifyEntry(user_id,-1,u_lcl_rk,u);
                            candidate_User.insert(ve);
                            user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                            cout<<"social_textual find candidate 用户： u"<<usr_id<<endl;
#endif
                        }

                    }
                }

            }
            else{
                rareKeys.insert(term);
                vector<int> term_related_users = getTermUserInvList(term);
                for(int usr: term_related_users)
                    rareKeywordRelated_usrSet.insert(usr);
            }

        }

        nonRareKeywordMap[poi_id] = nonRareKeys;

        set<int> poiRelated_usr; poiRelated_usr.clear();

#ifdef TRACKNVDBATCH
        cout<<"把非低频关键词对应的用户提取出来，分析他们的特征分布"<<endl;
#endif
        ////把非低频关键词对应的用户提取出来，分析他们的特征分布
        for(int term: nonRareKeys){
            vector<int> term_related_user = getTermUserInvList(term);
            for(int user_id: term_related_user){
                if(user_id>=UserID_MaxKey) continue;
                int up_key = user_id + UserID_MaxKey* poi_th;   //重要：每个poi对应（0,usrnum-1)
                if(userMark[up_key] == true) continue; ////仅正对当前poi下
                poiRelated_usr.insert(user_id);
            }

        }
        for(int u_id: poiRelated_usr){
            User user = getUserFromO2UOrgLeafData(u_id);
            int u_Ni = user.Ni;
            int u_leaf = Nodes[u_Ni].gtreepath.back();
            pnodes_uOCCURListMap[poi_id][u_leaf].insert(u_id);  //更新叶节点的uOCCURList
            vector<int> u_keywords = user.keywords;
            for(int term:u_keywords){
                if(isCovered(term,nonRareKeys)){
                    for(int term2:u_keywords){
                        if(isCovered(term2,nonRareKeys)){
                            pnodes_uCLUSTERListMap[poi_id][u_leaf][term].insert(term2); //更新叶节点的uCLUSTERList
                        }

                    }
                }

            }

            //更新上层节点的uOCCURList、uCLUSTERList
            int father = GTree[u_leaf].father;
            while(father!=-1){
                pnodes_uOCCURListMap[poi_id][father].insert(u_id);
                for(int term:u_keywords){
                    if(isCovered(term,nonRareKeys)){
                        for(int term2:u_keywords){
                            if(isCovered(term2,nonRareKeys)){
                                pnodes_uCLUSTERListMap[poi_id][father][term].insert(term2);
                            }

                        }

                    }
                }

                father = GTree[father].father;
            }
        }


        ////对包含特稀疏关键词的user （这样的user通常较少）进行剪枝
        int _size = rareKeywordRelated_usrSet.size();
        //int _size2 = midKeywordRelated_usrSet.size();
        //cout<<"包含POI下低频分布关键词的user个数="<<_size<<endl;

        for(int usr_id:rareKeywordRelated_usrSet){

            int up_key = usr_id + UserID_MaxKey* poi_th;
            //cout<<"access u"<<usr_id<<endl;

            if(upEntryMark[up_key] == true) continue;
            upEntryMark[up_key] = true;
            User u = getUserFromO2UOrgLeafData(usr_id);
            //printUsrInfo(u);

            double relevance = 0;
            relevance = textRelevance(u.keywords, poi.keywords);
            //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
            if(relevance==0) continue;
            //该用户社交文本相关
            //cout<<"评估 社交文本相关用户：u"<<usr_id<<endl;
            double gsk_score = getSKScore_o2u_phl(a,alpha,poi,u);
            double u_lcl_rk = getUserLCL_NVD(u,K,a,alpha);
            if(gsk_score > u_lcl_rk){
                //更新 user-poi associating list
                candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                if(candidate_usr_related_store[usr_id].size()>K){
                    candidate_usr_related_store[usr_id].pop();
                }
                //加入candidate user set
                BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_lcl_rk,u);
                candidate_User.insert(ve);
                user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                cout<<"rare textual find candidate 用户： u"<<usr_id<<endl;
#endif
            }
        }

#ifdef TRACKNVDBATCH
        cout<<"特稀疏关键词相关的 user 总体个数="<<rareKeywordRelated_usrSet.size()<<endl;
    //cout<<"完成特稀疏关键词相关的 user的剪枝，此时 candidate user size="<<candidate_id.size()<<"，具体为："<<endl;
    //printElements(candidate_id);
#endif
        //getchar();

        ////1. q_o 所在leaf 中的关键词相关 user先加入队列考察
        set<int> user_nearSet;
        for(int term: nonRareKeys){
            vector<int> u_list = getUsrTermRelatedEntry(term,qo_leaf);
            if(u_list.size()>0){
                for(int u_id: u_list){
                    user_nearSet.insert(u_id);
                }
            }
        }



#ifdef TRACKNVDBATCH
        if(user_nearSet.size() != 0){
        cout<<"q_o 所在leaf"<<qo_leaf<<"中的关键词相关 user先加入队列考察"<<endl;
    } else{
        cout<<"q_o 所在leaf "<<qo_leaf<<"中的无关键词相关user"<<endl;
    }

        cout<<"插入user_nearSet前,Queue.size="<<Queue.size()<<endl;
#endif

        for(int usr_id: user_nearSet){
            int up_key = usr_id  + UserID_MaxKey * poi_th;
            if(upEntryMark[up_key] == true) continue;

            User user = getUserFromO2UOrgLeafData(usr_id);
            double dist = usrToPOIDistance_phl(user, poi);
            set<int> u_keySet;
            /*for(int term: user.keywords)
                u_keySet.insert(term);*/
            BatchCheckEntry entry(usr_id,poi_id, qo_leaf,dist,user);  //zhiq
            Queue.push(entry);
        }
#ifdef TRACKNVDBATCH
        cout<<"插入user_nearSet后,Queue.size="<<Queue.size()<<endl;
#endif
        //getchar();


        //设定当前level最高的Gtree Node
        int node_highest = qo_leaf;
        POIHighestNode pair(poi_id,node_highest);
        T_max.push(pair);

        poi_th ++;

    }

    initial_endTime = clock();

    //cout<<"总体上，包含用户下非高频分布关键词的user个数="<<midKeywordRelated_usrSet.size()<<endl;
    //cout<<"完成初始步骤，开始从各个qo_leaf出发进行拓展与验证"<<endl;
    double initial_time = (double) (initial_endTime-initial_startTime)/CLOCKS_PER_SEC*1000;

    //getchar();



    ////2.从各qo_leaf出发按距离优先,进行路网拓展（在Gtree中从qo_leaf开始，自顶向上，遍历各个关键词相关节点）含有相关关键词的 user,user_leaf, user_node

    map<int,map<int,bool>> borderEntryVist;  //border_id, node_id
    map<int,bool> userEntryVist;
    while_startTime = clock();
    while(Queue.size()>0 ||T_max.size()>0){
        if(Queue.empty()){ //当前范围内已全部遍历， 需扩大路网遍历范围
#ifdef TRACKNVDBATCH
            cout<<"当前所有poi的检查范围内的相关user皆已被全部遍历!"<<endl;
#endif

            BatchExtend_Range(T_max,global_itm,Queue,nonRareKeywordMap,candidate_POIs, K,a,alpha);
            if(Queue.empty()) break;  //路网拓展可被终止

        }
        BatchCheckEntry topEntry = Queue.top();
        Queue.pop();

        float dist_parentBorder = topEntry.dist;
        int poi_id = topEntry.p_id;
        POI poi = candidate_POIs[poi_id];
        int poi_th = poi_Id2IdxMap[poi_id];

        if(topEntry.isUser){ //// an entry for user
            User user = topEntry.user;
            int usr_id = user.id;


            int up_key = usr_id + UserID_MaxKey* poi_th;
            //cout<<"up_key="<<up_key<<endl;
            if(upEntryMark[up_key])
                continue;
            else
                upEntryMark[up_key] = true;


            double u_SB_rk = getUserLCL_NVD(user,K,a,alpha);

            ////cout << "getUserLCL_NVD runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;

            double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,user);
            ////cout << "getGSKScore_o2u runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;
#ifdef TRACK
            cout<<"------DEQUEUE: user entry: u"<<usr_id<<"for p"<<poi_id<<"------"<<endl;
            cout<<"u_SB_rk="<<u_SB_rk<<endl;
            cout<<"gsk_score="<<gsk_score<<endl;


#endif

            if(gsk_score > u_SB_rk){
                if(candidate_usr_related_store.count(usr_id)){
                    candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0)); //jinsg
                    if(candidate_usr_related_store[usr_id].size()>K){
                        candidate_usr_related_store[usr_id].pop();
                    }
                }
                else{
                    priority_queue<ResultCurrent> list;
                    list.push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                    candidate_usr_related_store[usr_id] = list;
                }

                //加入candidate user set
                BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_SB_rk,user);
                candidate_User.insert(ve);

            }
        }
        else{  //// an entry for border
            int border_id = topEntry.id;
            int poi_id = topEntry.p_id;
            int poi_th = poi_Id2IdxMap[poi_id];
            int node_id =topEntry.node_id;
            int bEntry_key =  border_id +  node_id*VertexNum  + poi_th*inner_size;

#ifdef  TRACK
            if(GTree[node_id].isleaf)
                cout<<"----DEQUEUE: border entry: b"<<border_id<<" in leaf"<<node_id<<"for p"<<poi_id<<"，出列！------"<<endl;
            else
                cout<<"----DEQUEUE: border entry: b"<<border_id<<" in node"<<node_id<<"for p"<<poi_id<<"，出列！------"<<endl;


            cout<<"border_id: v"<<border_id<<", poi_th="<<poi_th<<endl;
            cout<<"border_key:"<<bEntry_key<<endl;
#endif

            // func(bEntry_key)= v_id+v_num*node_id + v_num*node_num*poi_th



            if(borderMark[bEntry_key]==true) continue;
            borderMark[bEntry_key] = true;

#ifdef   userspacity
            if(pnodes_uOCCURListMap[poi_id][node_id].size()< NODE_SPARCITY_VALUE_BATCH || dist_parentBorder>LONG_TERM){//&& dist_parentBorder>LONG_TERM
                for(int usr_id: pnodes_uOCCURListMap[poi_id][node_id]){
                    int up_key = usr_id + UserID_MaxKey *poi_th;
                    if(upEntryMark[up_key] == true) continue;

                    User user = getUserFromO2UOrgLeafData(usr_id);
                    double dist = usrToPOIDistance_phl(user, poi);
                    set<int> u_keySet;
                    BatchCheckEntry entry(usr_id,poi_id, node_id,dist,user);
                    Queue.push(entry);
                }
                int np_key = node_id + GTree.size()*poi_th;
                nodepMark[np_key] =true;
                continue;
            }
#endif


            set<int> Keys;

            //CheckEntry border_entry,int K, int a, double alpha
            Keys = getBorder_SB_BatchNVD(topEntry,K,a, alpha);  //joke
            //Keys = getBorder_SB_BatchNVD_Cluster(topEntry,K,a, alpha,pnodes_uCLUSTERListMap);

            //cout<<"remained (user) Keywords"; printSetElements(Keys);
            if(Keys.size() ==0) {
#ifdef TRACK
                cout<<"所有其后继被 ****prune！"<<endl;
#endif
                continue;
            }
#ifdef TRACK
            cout<<"其后继**不能**被prune！"<<endl;
#endif

            //若为叶节点，加入叶节点下的各个user
            if(GTree[node_id].isleaf){
#ifdef TRACK
                cout<<"b"<<border_id<<", 位于leaf"<<node_id<<endl;
#endif
                set<int> user_within;
                for(int term: Keys){
                    vector<int> u_list = getUsrTermRelatedEntry(term,node_id);
                    if(u_list.size()>0){
                        for(int u_id: u_list){
                            user_within.insert(u_id);
                        }
                    }
                }
                //加入叶节点下的各个user
                for(int usr_id: user_within){
                    //if(userMarkMap[poi_id][usr_id] == true) continue;
                    int up_key = usr_id + UserID_MaxKey *poi_th;
                    if(upEntryMark[up_key] == true) continue;

                    User user = getUserFromO2UOrgLeafData(usr_id);
                    double dist = usrToPOIDistance_phl(user, poi);
                    set<int> u_keySet;
                    BatchCheckEntry entry(usr_id,poi_id, node_id,dist,user);
                    Queue.push(entry);
#ifdef TRACK
                    cout<<"在main中加入user entry:"; entry.printRlt();
                    cout<<"Queue.size="<<Queue.size()<<endl;
#endif
                }
            }
            else{
                //获得父节点下所有关键词相关的孩子节点(line 22-25)
#ifdef TRACK
                cout<<"b"<<border_id<<", 位于node"<<node_id<<endl;
#endif
                set<int> child_set;
                for(int term:Keys){
                    vector<int> related_childs = getUsrTermRelatedEntry(term,node_id);
                    if(related_childs.size()>0){
                        for(int cn_id: related_childs){
                            child_set.insert(cn_id);
                        }
                    }
                }
#ifdef TRACK
                cout<<"关键词相关的后继："; printSetElements(child_set);
#endif

                int parent_border_th = topEntry.b_th;
                int posa = parent_border_th;
                for(int cn_id: child_set){
                    for(int b_th=0;b_th<GTree[cn_id].borders.size();b_th++){
                        int child_border_id = GTree[cn_id].borders[b_th];
                        int posb = GTree[cn_id].up_pos[b_th];
                        float b2b_dist = GTree[node_id].mind[posb * GTree[node_id].union_borders.size() + posa];

                        float dist_childBorder2 = topEntry.dist + b2b_dist;
                        float dist_childBorder = getDistance_phl(child_border_id,poi);  //jins


                        int border_key =  child_border_id +  cn_id*VertexNum  + poi_th*inner_size;
                        if(borderMark[border_key] == true) continue;
#ifdef TRACK
                        cout<<"child_border_id: v"<<child_border_id<<", poi_th"<<poi_th<<endl;
                        cout<<"border_key:"<<border_key<<endl;
#endif
                        BatchCheckEntry entry(child_border_id,poi_id, cn_id,Keys,dist_childBorder,b_th);
                        Queue.push(entry);
#ifdef TRACK
                        cout<<"在main中加入border entry:"; entry.printRlt();
                        cout<<"Queue.size="<<Queue.size()<<endl;
#endif


                    }
                }

            }
        }

    }
    filter_endTime = clock();


    ////对不能被减枝的user, 进行最终TkGSKQ求解

//nvd_verification
    verification_startTime = clock();

    int candidate_size = candidate_User.size();
    vector<ResultDetail> potential_customerDetails;
    vector<BatchCardinalityFisrt> user_list;

    for(BatchVerifyEntry ve:candidate_User){
        int u_id = ve.u_id;
        double u_rk_score = ve.rk_current;
        User u = ve.u;

        priority_queue<ResultCurrent> relate_query_queue = candidate_usr_related_store[u_id];
        priority_queue<ResultLargerFisrt> Lu;
        while(!relate_query_queue.empty()){
            ResultCurrent rc = relate_query_queue.top();
            relate_query_queue.pop();
            //先算出真实gsk评分
            double rel = rc.relevance; double inf = rc.influence;
            double gsk_score = rc.score;
            double score_upper = rc.score_upper;
            double _score = rc.score;
            if(gsk_score < u_rk_score)  //jins
                continue;


            ResultLargerFisrt rlf(rc.o_id,gsk_score,rel, inf);
            Lu.push(rlf);

            //max_score = max(max_score,gsk_score);
        }
        if(Lu.size()>0){
            BatchCardinalityFisrt cf(u_id,u_rk_score, Lu,u);
            user_list.push_back(cf);
        }
    }
    //bool* poiMark = new bool[poi_num];
    //bool* poiADJMark = new bool[poi_num];

    for(int i=0;i< user_list.size();i++){

        BatchCardinalityFisrt cf = user_list[i];
        int usr_id = cf.usr_id;


#ifdef TRACKNVDBATCH_VERIFY
        cout<<"验证候选用户u"<<usr_id<<"的topk结果"<<endl;
        if(usr_id==1184962){
            cout<<"find u 1184962"<<endl;
        }
#endif

        double current_sk = cf.current_rk;
        double gsk_score_max = cf.max_score;


        User u = cf.u;
        //double usr_lcl_rt = cf.current_rk;
        priority_queue<ResultLargerFisrt> Lu = cf.Lu;
        double score_max = Lu.top().score;

        double Rk_u =0;
        if(true){

            TopkQueryCurrentResult topk_r6 = TkGSKQ_NVD(u, K, DEFAULT_A,alpha);//,poiMark,poiADJMark);
            Rk_u = topk_r6.topkScore;
        }
        else{
            //priority_queue<Result> results = TkGSKQ_bottom2up_verify_disk(usr_id,K,DEFAULT_A,alpha,gsk_score_max,current_sk);//1000,0.0);
            priority_queue<Result> results = TkGSKQ_bottom2up_verify_disk(usr_id,K,DEFAULT_A,alpha,1000,0.0);

            if(results.size()==K){
                Rk_u= results.top().score;     //gangcai
            }
            else{
                Rk_u = current_sk;
            }

            if(usr_id==1184962){
                priority_queue<Result> results2 = results;
                int rank = 1;
                while(results2.size()>0){
                    Result rr = results2.top();
                    results2.pop();
                    cout<<"**rank"<<rank<<" :  ";
                    rr.printResult();
                    rank++;
                }
                getchar();
            }
        }


        //验证结果
        if( gsk_score_max >= Rk_u){
#ifdef TRACKNVDBATCH_VERIFY
            cout<<"*************发现潜在用户u"<<usr_id<<",Rk_u="<<Rk_u<<endl;
#endif
            //对Lu中的各个query object进行评估
            while(!Lu.empty()){
                ResultLargerFisrt rlf = Lu.top();
                Lu.pop();
                double rel = 0; double inf = 0;
                double gsk_score = rlf.score; int query_object = rlf.o_id;
                if(gsk_score > Rk_u){
                    ResultDetail rd(usr_id, 0, 0, 0, gsk_score, Rk_u);
                    batch_results[query_object].push_back(rd);
#ifdef TRACK
                    cout<<"对于o"<<query_object<<"确实为有效用户，gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
#endif
                }
                else break;   //fanz
            }
        }




#ifdef TRACKNVDBATCH_VERIFY
        cout<<"验证候选用户"<<usr_id<<", gsk_score_max="<<gsk_score_max<<" ,Rk_u="<<Rk_u<<endl;
#endif

    }
    //delete []poiMark; delete []poiADJMark;

    verification_endTime = clock();
    endTime = clock();


    printBatchRkGSKQResults(stores, batch_results);
    cout<<"filter阶段总共发现"<<candidate_User.size()<<"个候选用户！"<<endl;
    Uc_size = candidate_User.size();
    cout<<"Filter 阶段 用时："<< (double)(filter_endTime - filter_startTime) / CLOCKS_PER_SEC*1000  <<
        "ms (其中 while loop用时:"<< (double)(filter_endTime - while_startTime) / CLOCKS_PER_SEC*1000<<"ms)"<<endl;
    cout<<"Verification 阶段 用时："<< (double)(verification_endTime - verification_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;

    //cout<<"user_textual_count="<<user_textual_count<<",user_textual_useful_count="<<user_textual_useful_count<<endl;
    cout<<"---------BatchRkGSKQ_by_NVD 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC*1000  << "ms-----------" << endl;


    return batch_results;
}



BatchResults RkGSKQ_Batch_NVD_SyntheticData(vector<int> stores, int K, float a, float alpha) {  //gangcai

    cout<<"------------Running Batch RkGSKQ using NVD&IG-Tree!------------"<<endl;
    clock_t startTime, endTime;
    clock_t filter_startTime,while_startTime, filter_endTime;
    clock_t verification_startTime, verification_endTime;

    //vector<int> RkGSKQ_results;
    BatchResults batch_results;



    unordered_map<int, POI> candidate_POIs;
    unordered_map<int, int> poi_Id2IdxMap;
    unordered_map<int, set<int>>  nonRareKeywordMap;


    //候选验证用户
    set<BatchVerifyEntry> candidate_User;  //candidate users
    unordered_map<int, priority_queue<ResultCurrent>> candidate_usr_related_store;   // user_id, {related poi1,2,....}


    //border已被评估的标志向量
    int node_size = GTree.size();
    int inner_size = node_size*VertexNum;  // func(bEntry_key)= v_id+v_num*node_id + v_num*node_num*poi_th
    long bEntry_MaxKey =  inner_size * stores.size();
    unordered_map<int,bool> borderMark;
    /*bool* borderMark = new bool[bEntry_MaxKey];
    for(int i=0;i< bEntry_MaxKey;i++){
        borderMark[i] = false;
    }*/


    //user已被评估的标志向量表
    int upEntry_MaxKey =  UserID_MaxKey * stores.size();
    //unordered_map<int,bool> upEntryMark;
    bool* upEntryMark = new bool[upEntry_MaxKey];
    for(int i=0;i<upEntry_MaxKey;i++){
        upEntryMark[i] = false;
    }
    //memset(upEntryMark,false,sizeof(upEntryMark));

    bool userMark[UserID_MaxKey];
    memset(userMark,false, sizeof(userMark));

    //标记gtree node、 poi 对的标志向量
    int np_maxKey = GTree.size()*stores.size();
    bool nodepMark[np_maxKey];
    memset(nodepMark,false, sizeof(nodepMark));


    startTime = clock();
    filter_startTime = clock();
    //距离优先队列

    ////初始化gtree node的 uOCCURList以及 uCLUSTERList
    unordered_map<int,UOCCURList> pnodes_uOCCURListMap;
    unordered_map<int,UCLUSTERList> pnodes_uCLUSTERListMap;



    priority_queue<BatchCheckEntry> Queue;
    priority_queue<POIHighestNode> T_max;


    int root = 0;
    //// 1. 获取各个poi信息, 并计算各个 poi 到 o_leaf上所有节点的border的距离
    unordered_map<int, unordered_map<int, vector<int>>> global_itm;
    int poi_th = 0;
    set<int> midKeywordRelated_usrSet;
    for(int poi_id: stores){
        poi_Id2IdxMap[poi_id] = poi_th;


        POI poi = getPOIFromO2UOrgLeafData(poi_id);

        candidate_POIs[poi_id] = poi;

        //printPOIInfo(poi);


        //初始化并构建该poi对应的 nodes_uOCCURList 及 nodes_uCLUSTERList
        UOCCURList nodes_uOCCURList;
        UCLUSTERList nodes_uCLUSTERList;  //vector_index= node_id,  term_id,
        for(int i=0;i<GTree.size();i++){
            set<int> node_occurList;
            nodes_uOCCURList.push_back(node_occurList);
            unordered_map<int,set<int>> clusters;
            nodes_uCLUSTERList.push_back(clusters);
        }
        pnodes_uOCCURListMap[poi_id] = nodes_uOCCURList;
        pnodes_uCLUSTERListMap[poi_id] = nodes_uCLUSTERList;


        int p_Ni = poi.Ni; int p_Nj = poi.Nj;
        vector<int> qo_keys = poi.keywords;
        int qo_leaf = Nodes[p_Ni].gtreepath.back();
        TreeNode tnData = getGIMTreeNodeData(qo_leaf,OnlyU);
#ifdef TRACK
        cout<<"POI信息："; printPOIInfo(poi);
        cout<<"o_leaf为leaf"<<qo_leaf<<endl;
        cout<<"o_leaf中user的关键词并集："<<endl;
        printSetElements(tnData.userUKeySet);
#endif

        //// 计算 各个 poi 到 o_leaf上所有节点的border的距离
        MaterilizedDistanceMap itm; // intermediate answer, tree node -> array
        itm.clear();
        int cid, posa, posb, _min, dis;
        int o_leaf = Nodes[poi.Ni].gtreepath.back();
        itm[o_leaf].clear();
        int pos_Ni = id2idx_vertexInLeaf[p_Ni];
        int pos_Nj = id2idx_vertexInLeaf[p_Nj];
        int tn = o_leaf;
        for (int j = 0; j < GTree[tn].borders.size(); j++) {
            int dist_via_Ni = GTree[tn].mind[j * GTree[tn].leafnodes.size() + pos_Ni]+poi.dis;
            int dist_via_Nj = GTree[tn].mind[j * GTree[tn].leafnodes.size() + pos_Nj]+(poi.dist-poi.dis);
            int dist_to_border = min(dist_via_Ni, dist_via_Nj);
            itm[tn].push_back(dist_to_border);
        }
        global_itm[poi_id] = itm;

        ////获得与poi 社交文本相关的 user
        vector<int> checkin_users = poiCheckInIDList[poi_id];
        int user_textual_count = 0; int user_textual_useful_count = 0;
        for(int check_usr_id: checkin_users){
            vector<int> friends = friendshipMap[check_usr_id];
            for(int friend_id: friends){
                if(friend_id>=UserID_MaxKey) continue;
                int up_key = friend_id + UserID_MaxKey* poi_th;   //重要：每个poi对应（0,usrnum-1)

                if(upEntryMark[up_key] == true) continue;
                upEntryMark[up_key] = true;


                User u = getUserFromO2UOrgLeafData(friend_id);
                int usr_id = u.id;
                double relevance = 0;
                relevance = textRelevance(u.keywords, poi.keywords);
                //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
                if(relevance==0) continue;
                //该用户社交文本相关
#ifdef TRACKNVDBATCH
                cout<<"评估 社交文本相关用户：u"<<usr_id<<endl;
#endif
                user_textual_count ++;
                double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,u);
                double u_lcl_rk = getUserLCL_NVD(u,K,a,alpha);
                if(gsk_score > u_lcl_rk){

                    //更新 user-poi associating list
                    candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                    if(candidate_usr_related_store[usr_id].size()>K){
                        candidate_usr_related_store[usr_id].pop();
                    }
                    //加入candidate user set
                    BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_lcl_rk,u);
                    candidate_User.insert(ve);

                    user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                    cout<<"social_textual find candidate 用户： u"<<usr_id<<endl;
#endif
                }

            }
        }
#ifdef TRACKNVDBATCH
        cout<<"社交文本相关的 user总体个数="<<user_textual_count<<endl;
    //cout<<"完成社交文本相关的 user的剪枝，此时 candidate user size="<<candidate_id.size()<<"，具体为："<<endl;
    //printElements(candidate_id);
#endif
        //getchar();

        ////将q_o查询关键词分为 “（rare)稀有关键词”与“(none-rare)非稀有关键词” 两类,获得包含特稀疏keyword( |inv(term)|<K )的 user
        set<int> rareKeys;  set<int> nonRareKeys;
        set<int> rareKeywordRelated_usrSet;
        for(int term: qo_keys){
            int inv_olist_size = getTermOjectInvListSize(term);
            if(inv_olist_size>K){
                vector<int> usr_termRelated = getTermUserInvList(term);
                int inv_ulist_size = usr_termRelated.size();
                if(inv_ulist_size> uterm_SPARCITY_VALUE){
                    nonRareKeys.insert(term);
                    //printf("term %d 为用户高频词汇，u_frequency=%d\n", term,inv_ulist_size);
                }
                else{ ////把一些带稀有关键词的user提前进行验证，以提升后续剪枝的效果
                    for(int user_id: usr_termRelated){
                        if(user_id>=UserID_MaxKey) continue;
                        int up_key = user_id + UserID_MaxKey* poi_th;   //重要：每个poi对应（0,usrnum-1)

                        if(upEntryMark[up_key] == true) continue;
                        upEntryMark[up_key] = true;
                        midKeywordRelated_usrSet.insert(user_id);

                        User u = getUserFromO2UOrgLeafData(user_id);
                        int usr_id = u.id;
                        double relevance = 0;
                        relevance = textRelevance(u.keywords, poi.keywords);
                        //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
                        if(relevance==0) continue;
                        //该用户社交文本相关
#ifdef TRACKNVDBATCH
                        cout<<"评估 社交文本相关用户：u"<<usr_id<<endl;
#endif
                        user_textual_count ++;
                        double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,u);
                        double u_lcl_rk = getUserLCL_NVD(u,K,a,alpha);
                        if(gsk_score > u_lcl_rk){

                            //更新 user-poi associating list
                            candidate_usr_related_store[user_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                            if(candidate_usr_related_store[user_id].size()>K){
                                candidate_usr_related_store[user_id].pop();
                            }
                            //加入candidate user set
                            BatchVerifyEntry ve = BatchVerifyEntry(user_id,-1,u_lcl_rk,u);
                            candidate_User.insert(ve);
                            user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                            cout<<"social_textual find candidate 用户： u"<<usr_id<<endl;
#endif
                        }

                    }
                }

            }
            else{
                rareKeys.insert(term);
                vector<int> term_related_users = getTermUserInvList(term);
                for(int usr: term_related_users)
                    rareKeywordRelated_usrSet.insert(usr);
            }

        }
        nonRareKeywordMap[poi_id] = nonRareKeys;

        set<int> poiRelated_usr; poiRelated_usr.clear();

        ////把非低频关键词对应的用户提取出来，分析他们的特征分布
        for(int term: nonRareKeys){
            vector<int> term_related_user = getTermUserInvList(term);
            for(int user_id: term_related_user){
                if(user_id>=UserID_MaxKey) continue;
                int up_key = user_id + UserID_MaxKey* poi_th;   //重要：每个poi对应（0,usrnum-1)
                if(upEntryMark[up_key] == true) continue; ////仅正对当前poi下
                poiRelated_usr.insert(user_id);
            }

        }
        for(int u_id: poiRelated_usr){
            User user = getUserFromO2UOrgLeafData(u_id);
            int u_Ni = user.Ni;
            int u_leaf = Nodes[u_Ni].gtreepath.back();
            pnodes_uOCCURListMap[poi_id][u_leaf].insert(u_id);  //更新叶节点的uOCCURList
            vector<int> u_keywords = user.keywords;
            for(int term:u_keywords){
                if(isCovered(term,nonRareKeys)){
                    for(int term2:u_keywords){
                        if(isCovered(term2,nonRareKeys)){
                            pnodes_uCLUSTERListMap[poi_id][u_leaf][term].insert(term2); //更新叶节点的uCLUSTERList
                        }

                    }
                }

            }

            //更新上层节点的uOCCURList、uCLUSTERList
            int father = GTree[u_leaf].father;
            while(father!=-1){
                pnodes_uOCCURListMap[poi_id][father].insert(u_id);
                for(int term:u_keywords){
                    if(isCovered(term,nonRareKeys)){
                        for(int term2:u_keywords){
                            if(isCovered(term2,nonRareKeys)){
                                pnodes_uCLUSTERListMap[poi_id][father][term].insert(term2);
                            }

                        }

                    }
                }

                father = GTree[father].father;
            }
        }


        ////对包含特稀疏关键词的user （这样的user通常较少）进行剪枝
        int _size = rareKeywordRelated_usrSet.size();
        //int _size2 = midKeywordRelated_usrSet.size();
        //cout<<"包含POI下低频分布关键词的user个数="<<_size<<endl;

        for(int usr_id:rareKeywordRelated_usrSet){

            int up_key = usr_id + UserID_MaxKey* poi_th;
            //cout<<"access u"<<usr_id<<endl;

            if(upEntryMark[up_key] == true) continue;
            upEntryMark[up_key] = true;
            User u = getUserFromO2UOrgLeafData(usr_id);
            //printUsrInfo(u);

            double relevance = 0;
            relevance = textRelevance(u.keywords, poi.keywords);
            //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
            if(relevance==0) continue;
            //该用户社交文本相关
            //cout<<"评估 社交文本相关用户：u"<<usr_id<<endl;
            double gsk_score = getSKScore_o2u_phl(a,alpha,poi,u);
            double u_lcl_rk = getUserLCL_NVD(u,K,a,alpha);
            if(gsk_score > u_lcl_rk){
                //更新 user-poi associating list
                candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                if(candidate_usr_related_store[usr_id].size()>K){
                    candidate_usr_related_store[usr_id].pop();
                }
                //加入candidate user set
                BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_lcl_rk,u);
                candidate_User.insert(ve);
                user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                cout<<"rare textual find candidate 用户： u"<<usr_id<<endl;
#endif
            }
        }

#ifdef TRACKNVDBATCH
        cout<<"特稀疏关键词相关的 user 总体个数="<<rareKeywordRelated_usrSet.size()<<endl;
    //cout<<"完成特稀疏关键词相关的 user的剪枝，此时 candidate user size="<<candidate_id.size()<<"，具体为："<<endl;
    //printElements(candidate_id);
#endif
        //getchar();

        ////1. q_o 所在leaf 中的关键词相关 user先加入队列考察
        set<int> user_nearSet;
        for(int term: nonRareKeys){
            vector<int> u_list = getUsrTermRelatedEntry(term,qo_leaf);
            if(u_list.size()>0){
                for(int u_id: u_list){
                    user_nearSet.insert(u_id);
                }
            }
        }



#ifdef TRACKNVDBATCH
        if(user_nearSet.size() != 0){
        cout<<"q_o 所在leaf"<<qo_leaf<<"中的关键词相关 user先加入队列考察"<<endl;
    } else{
        cout<<"q_o 所在leaf "<<qo_leaf<<"中的无关键词相关user"<<endl;
    }

        cout<<"插入user_nearSet前,Queue.size="<<Queue.size()<<endl;
#endif

        for(int usr_id: user_nearSet){
            int up_key = usr_id  + UserID_MaxKey * poi_th;  
#ifdef TRACKNVDBATCH
            cout<<"usr_id="<<usr_id<<",poi_th="<<poi_th<<",up_key="<<up_key<<endl;
#endif
            if(upEntryMark[up_key] == true) continue;

            User user = getUserFromO2UOrgLeafData(usr_id);

            double dist = usrToPOIDistance_phl(user, poi);
#ifdef TRACKNVDBATCH
            cout<<"dist(u"<<user.id<<",p"<<poi.id<<")="<<dist<<endl;
#endif
            set<int> u_keySet;
           
            BatchCheckEntry entry(usr_id,poi_id, qo_leaf,dist,user);  //zhiq
            Queue.push(entry);
        }
#ifdef TRACKNVDBATCH
        cout<<"插入user_nearSet后,Queue.size="<<Queue.size()<<endl;
#endif
        //getchar();


        //设定当前level最高的Gtree Node
        int node_highest = qo_leaf;
        POIHighestNode pair(poi_id,node_highest);
        T_max.push(pair);

        poi_th ++;

    }

    cout<<"总体上，包含用户下非高频分布关键词的user个数="<<midKeywordRelated_usrSet.size()<<endl;
    cout<<"完成初始步骤，开始从各个qo_leaf出发进行拓展与验证"<<endl;
    //getchar();



    ////2.从各qo_leaf出发按距离优先,进行路网拓展（在Gtree中从qo_leaf开始，自顶向上，遍历各个关键词相关节点）含有相关关键词的 user,user_leaf, user_node

    map<int,map<int,bool>> borderEntryVist;  //border_id, node_id
    map<int,bool> userEntryVist;
    while_startTime = clock();
    while(Queue.size()>0 ||T_max.size()>0){
        if(Queue.empty()){ //当前范围内已全部遍历， 需扩大路网遍历范围
#ifdef TRACKEXTEND
            cout<<"当前所有poi的检查范围内的相关user皆已被全部遍历!"<<endl;
#endif

            BatchExtend_Range(T_max,global_itm,Queue,nonRareKeywordMap,candidate_POIs, K,a,alpha);
            if(Queue.empty()) break;  //路网拓展可被终止

        }
        BatchCheckEntry topEntry = Queue.top();
        Queue.pop();

        float dist_parentBorder = topEntry.dist;
        int poi_id = topEntry.p_id;
        POI poi = candidate_POIs[poi_id];
        int poi_th = poi_Id2IdxMap[poi_id];

        if(topEntry.isUser){ //// an entry for user
            User user = topEntry.user;
            int usr_id = user.id;


            int up_key = usr_id + UserID_MaxKey* poi_th;
            //cout<<"up_key="<<up_key<<endl;
            if(upEntryMark[up_key])
                continue;
            else
                upEntryMark[up_key] = true;


            double u_SB_rk = getUserLCL_NVD(user,K,a,alpha);

            ////cout << "getUserLCL_NVD runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;

            double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,user);
            ////cout << "getGSKScore_o2u runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;
#ifdef WhileTRACK
            cout<<"------DEQUEUE: user entry: u"<<usr_id<<"for p"<<poi_id<<"------"<<endl;
            cout<<"u_SB_rk="<<u_SB_rk<<endl;
            cout<<"gsk_score="<<gsk_score<<endl;


#endif

            if(gsk_score > u_SB_rk){
                if(candidate_usr_related_store.count(usr_id)){
                    candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0)); //jinsg
                    if(candidate_usr_related_store[usr_id].size()>K){
                        candidate_usr_related_store[usr_id].pop();
                    }
                }
                else{
                    priority_queue<ResultCurrent> list;
                    list.push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                    candidate_usr_related_store[usr_id] = list;
                }

                //加入candidate user set
                BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_SB_rk,user);
                candidate_User.insert(ve);

            }
        }
        else{  //// an entry for border
            int border_id = topEntry.id;
            int poi_id = topEntry.p_id;
            int poi_th = poi_Id2IdxMap[poi_id];
            int node_id =topEntry.node_id;
            int bEntry_key =  border_id +  node_id*VertexNum  + poi_th*inner_size;

#ifdef  WhileTRACK
            if(GTree[node_id].isleaf)
                cout<<"----DEQUEUE: border entry: b"<<border_id<<" in leaf"<<node_id<<"for p"<<poi_id<<"，出列！------"<<endl;
            else
                cout<<"----DEQUEUE: border entry: b"<<border_id<<" in node"<<node_id<<"for p"<<poi_id<<"，出列！------"<<endl;

            if(node_id==1190){
                cout<<"node_id==1190！"<<endl;
            }

            cout<<"border_id: v"<<border_id<<", poi_th="<<poi_th<<endl;
            cout<<"border_key:"<<bEntry_key<<endl;
#endif

            // func(bEntry_key)= v_id+v_num*node_id + v_num*node_num*poi_th



            if(borderMark[bEntry_key]==true) continue;
            borderMark[bEntry_key] = true;

#ifdef   AGGRESSIVE_I
            int _size = pnodes_uOCCURListMap[poi_id][node_id].size();
            if(_size< NODE_SPARCITY_VALUE_BATCH){//&& dist_parentBorder>LONG_TERM
                for(int usr_id: pnodes_uOCCURListMap[poi_id][node_id]){
                    int up_key = usr_id + UserID_MaxKey *poi_th;
                    if(upEntryMark[up_key] == true) continue;
                    //cout<<"u"<<usr_id<<"被提前评估！"<<endl;
                    User user = getUserFromO2UOrgLeafData(usr_id);
                    double dist = usrToPOIDistance_phl(user, poi);
                    set<int> u_keySet;
                    BatchCheckEntry entry(usr_id,poi_id, node_id,dist,user);
                    Queue.push(entry);
                }
                int np_key = node_id + GTree.size()*poi_th;
                nodepMark[np_key] =true;
                continue;
            }
#endif


            set<int> Keys;

            //CheckEntry border_entry,int K, int a, double alpha
            //Keys = getBorder_SB_BatchNVD(topEntry,K,a, alpha);  //joke
            Keys = getBorder_SB_BatchNVD_Cluster(topEntry,K,a, alpha,pnodes_uCLUSTERListMap);

            //cout<<"remained (user) Keywords"; printSetElements(Keys);
            if(Keys.size() ==0) {
#ifdef TRACK
                cout<<"所有其后继被 ****prune！"<<endl;
#endif
                continue;
            }
#ifdef TRACK
            cout<<"其后继**不能**被prune！"<<endl;
#endif

            //若为叶节点，加入叶节点下的各个user
            if(GTree[node_id].isleaf){
#ifdef TRACK
                cout<<"b"<<border_id<<", 位于leaf"<<node_id<<endl;
#endif
                set<int> user_within;
                for(int term: Keys){
                    vector<int> u_list = getUsrTermRelatedEntry(term,node_id);
                    if(u_list.size()>0){
                        for(int u_id: u_list){
                            user_within.insert(u_id);
                        }
                    }
                }
                //加入叶节点下的各个user
                for(int usr_id: user_within){
                    //if(userMarkMap[poi_id][usr_id] == true) continue;
                    int up_key = usr_id + UserID_MaxKey *poi_th;
                    if(upEntryMark[up_key] == true) continue;

                    User user = getUserFromO2UOrgLeafData(usr_id);
                    double dist = usrToPOIDistance_phl(user, poi);
                    set<int> u_keySet;
                    BatchCheckEntry entry(usr_id,poi_id, node_id,dist,user);
                    Queue.push(entry);
#ifdef TRACK
                    cout<<"在main中加入user entry:"; entry.printRlt();
                    cout<<"Queue.size="<<Queue.size()<<endl;
#endif
                }
            }
            else{
                //获得父节点下所有关键词相关的孩子节点(line 22-25)
#ifdef WhileTRACK
                cout<<"b"<<border_id<<", 位于node"<<node_id<<endl;
#endif
                set<int> child_set;
                for(int term:Keys){
                    vector<int> related_childs = getUsrTermRelatedEntry(term,node_id);
                    if(related_childs.size()>0){
                        for(int cn_id: related_childs){
                            child_set.insert(cn_id);
                        }
                    }
                }
#ifdef TRACK
                cout<<"关键词相关的后继："; printSetElements(child_set);
#endif

                int parent_border_th = topEntry.b_th;
                int posa = parent_border_th;
                for(int cn_id: child_set){
                    for(int b_th=0;b_th<GTree[cn_id].borders.size();b_th++){
                        int child_border_id = GTree[cn_id].borders[b_th];
                        int posb = GTree[cn_id].up_pos[b_th];
                        float b2b_dist = GTree[node_id].mind[posb * GTree[node_id].union_borders.size() + posa];

                        float dist_childBorder2 = topEntry.dist + b2b_dist;
                        float dist_childBorder = getDistance_phl(child_border_id,poi);  //jins


                        int border_key =  child_border_id +  cn_id*VertexNum  + poi_th*inner_size;
                        if(borderMark[border_key] == true) continue;
#ifdef TRACK
                        cout<<"child_border_id: v"<<child_border_id<<", poi_th"<<poi_th<<endl;
                        cout<<"border_key:"<<border_key<<endl;
#endif
                        BatchCheckEntry entry(child_border_id,poi_id, cn_id,Keys,dist_childBorder,b_th);
                        Queue.push(entry);
#ifdef TRACK
                        cout<<"在main中加入border entry:"; entry.printRlt();
                        cout<<"Queue.size="<<Queue.size()<<endl;
#endif


                    }
                }

            }
        }

    }
    filter_endTime = clock();
    //delete []borderMark;
    delete [] upEntryMark;

    ////对不能被减枝的user, 进行最终TkGSKQ求解

//nvd_verification
    verification_startTime = clock();

    int candidate_size = candidate_User.size();
    vector<ResultDetail> potential_customerDetails;
    vector<BatchCardinalityFisrt> user_list;

    for(BatchVerifyEntry ve:candidate_User){
        int u_id = ve.u_id;
        double u_rk_score = ve.rk_current;
        User u = ve.u;

        priority_queue<ResultCurrent> relate_query_queue = candidate_usr_related_store[u_id];
        priority_queue<ResultLargerFisrt> Lu;
        while(!relate_query_queue.empty()){
            ResultCurrent rc = relate_query_queue.top();
            relate_query_queue.pop();
            //先算出真实gsk评分
            double rel = rc.relevance; double inf = rc.influence;
            double gsk_score = rc.score;
            double score_upper = rc.score_upper;
            double _score = rc.score;
            if(gsk_score < u_rk_score)  //jins
                continue;


            ResultLargerFisrt rlf(rc.o_id,gsk_score,rel, inf);
            Lu.push(rlf);

            //max_score = max(max_score,gsk_score);
        }
        if(Lu.size()>0){
            BatchCardinalityFisrt cf(u_id,u_rk_score, Lu,u);
            user_list.push_back(cf);
        }
    }
    //bool* poiMark = new bool[poi_num];
    //bool* poiADJMark = new bool[poi_num];

    for(int i=0;i< user_list.size();i++){

        BatchCardinalityFisrt cf = user_list[i];
        int usr_id = cf.usr_id;


#ifdef TRACKNVDBATCH_VERIFY
        cout<<"验证候选用户u"<<usr_id<<"的topk结果"<<endl;
        if(usr_id==1184962){
            cout<<"find u 1184962"<<endl;
        }
#endif

        double current_sk = cf.current_rk;
        double gsk_score_max = cf.max_score;


        User u = cf.u;
        //double usr_lcl_rt = cf.current_rk;
        priority_queue<ResultLargerFisrt> Lu = cf.Lu;
        double score_max = Lu.top().score;

        double Rk_u =0;
        if(true){

            TopkQueryCurrentResult topk_r6 = TkGSKQ_NVD(u, K, DEFAULT_A,alpha);//,poiMark,poiADJMark);
            Rk_u = topk_r6.topkScore;
        }
        else{
            //priority_queue<Result> results = TkGSKQ_bottom2up_verify_disk(usr_id,K,DEFAULT_A,alpha,gsk_score_max,current_sk);//1000,0.0);
            priority_queue<Result> results = TkGSKQ_bottom2up_verify_disk(usr_id,K,DEFAULT_A,alpha,1000,0.0);

            if(results.size()==K){
                Rk_u= results.top().score;     //gangcai
            }
            else{
                Rk_u = current_sk;
            }

        }


        //验证结果
        if( gsk_score_max >= Rk_u){
#ifdef TRACKNVDBATCH_VERIFY
            cout<<"*************发现潜在用户u"<<usr_id<<",Rk_u="<<Rk_u<<endl;
#endif
            //对Lu中的各个query object进行评估
            while(!Lu.empty()){
                ResultLargerFisrt rlf = Lu.top();
                Lu.pop();
                double rel = 0; double inf = 0;
                double gsk_score = rlf.score; int query_object = rlf.o_id;
                if(gsk_score > Rk_u){
                    ResultDetail rd(usr_id, 0, 0, 0, gsk_score, Rk_u);
                    batch_results[query_object].push_back(rd);
#ifdef TRACK
                    cout<<"对于o"<<query_object<<"确实为有效用户，gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
#endif
                }
                else break;   //fanz
            }
        }




#ifdef TRACKNVDBATCH_VERIFY
        cout<<"验证候选用户"<<usr_id<<", gsk_score_max="<<gsk_score_max<<" ,Rk_u="<<Rk_u<<endl;
#endif

    }
    //delete []poiMark; delete []poiADJMark;

    verification_endTime = clock();
    endTime = clock();


    printBatchRkGSKQResults(stores, batch_results);
    cout<<"filter阶段总共发现"<<candidate_User.size()<<"个候选用户！"<<endl;
    Uc_size = candidate_User.size();
    cout<<"Filter 阶段 用时："<< (double)(filter_endTime - filter_startTime) / CLOCKS_PER_SEC*1000  <<
        "ms (其中 while loop用时:"<< (double)(filter_endTime - while_startTime) / CLOCKS_PER_SEC*1000<<"ms)"<<endl;
    cout<<"Verification 阶段 用时："<< (double)(verification_endTime - verification_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;

    //cout<<"user_textual_count="<<user_textual_count<<",user_textual_useful_count="<<user_textual_useful_count<<endl;
    cout<<"---------BatchRkGSKQ_by_NVD 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC*1000  << "ms-----------" << endl;


    return batch_results;
}


//快速筛选各个poi的potential customers 的candidate


//Batch Reverse 快速过滤输出结果
typedef map<int, vector<ResultDetail>> FastBatch_filterResults;

typedef struct {
    set<BatchVerifyEntry> candidate_User;  //candidate users
    unordered_map<int, priority_queue<ResultCurrent>> candidate_usr_related_store; // user_id, {related poi1,2,....}
    FastBatch_filterResults batch_results;
    unordered_map<int, POI> candidate_POIs;

}BatchFastFilterResults;


BatchFastFilterResults Fast_Batch_NVDFilter0(vector<int> stores, int K, float a, float alpha) {  //gangcai

    cout<<"------------Running Fast_Batch_NVDFilter!------------"<<endl;
    clock_t startTime, endTime;
    clock_t filter_startTime,while_startTime, filter_endTime;
    clock_t verification_startTime, verification_endTime;
    clock_t initial_startTime, initial_endTime;

    //vector<int> RkGSKQ_results;
    FastBatch_filterResults batch_results;



    unordered_map<int, POI> candidate_POIs;
    unordered_map<int, int> poi_Id2IdxMap;
    unordered_map<int, set<int>>  nonRareKeywordMap;


    //候选验证用户
    set<BatchVerifyEntry> candidate_User;  //candidate users
    unordered_map<int, priority_queue<ResultCurrent>> candidate_usr_related_store;   // user_id, {related poi1,2,....}
    unordered_map<int, priority_queue<ResultCurrent>> candidate_usrLeaf_related_store;   // user_id, {related poi1,2,....}


    //border已被评估的标志向量
    int node_size = GTree.size();
    int inner_size = node_size*VertexNum;  // func(bEntry_key)= v_id+v_num*node_id + v_num*node_num*poi_th
    long bEntry_MaxKey =  inner_size * stores.size();
    unordered_map<int,bool> borderMark;



    //user已被评估的标志向量表
    int upEntry_MaxKey =  UserID_MaxKey * stores.size();
    //unordered_map<int,bool> upEntryMark;
    bool* upEntryMark = new bool[upEntry_MaxKey];
    for(int i=0;i<upEntry_MaxKey;i++){
        upEntryMark[i] = false;
    }

    bool userMark[UserID_MaxKey];
    memset(userMark,false, sizeof(userMark));

    //标记gtree node、 poi 对的标志向量
    int np_maxKey = GTree.size()*stores.size();
    bool nodepMark[np_maxKey];
    memset(nodepMark,false, sizeof(nodepMark));


    startTime = clock();
    filter_startTime = clock();
    //距离优先队列

    ////初始化gtree node的 uOCCURList以及 uCLUSTERList
    unordered_map<int,UOCCURList> pnodes_uOCCURListMap;
    unordered_map<int,UCLUSTERList> pnodes_uCLUSTERListMap;
    //unordered_map<int,UOCCURListMap> pnodes_uOCCURListMap;
    //unordered_map<int,UCLUSTERListMap> pnodes_uCLUSTERListMap;



    priority_queue<BatchCheckEntry> Queue;
    priority_queue<POIHighestNode> T_max;


    int root = 0;
    //// 1. 获取各个poi信息, 并计算各个 poi 到 o_leaf上所有节点的border的距离
    unordered_map<int, unordered_map<int, vector<int>>> global_itm;
    int poi_th = 0;
    set<int> midKeywordRelated_usrSet; double initial_time = -1;
    initial_startTime = clock();
    for(int poi_id: stores){
        poi_Id2IdxMap[poi_id] = poi_th;


        POI poi = getPOIFromO2UOrgLeafData(poi_id);

        candidate_POIs[poi_id] = poi;

        ////printPOIInfo(poi);


        //初始化并构建该poi对应的 nodes_uOCCURList 及 nodes_uCLUSTERList
        //UOCCURListMap nodes_uOCCURList;
        //UCLUSTERListMap nodes_uCLUSTERList;  //vector_index= node_id,  term_id,
        UOCCURList nodes_uOCCURList;
        UCLUSTERList nodes_uCLUSTERList;  //vector_index= node_id,  term_id,
        for(int i=0;i<GTree.size();i++){
            set<int> node_occurList;
            nodes_uOCCURList.push_back(node_occurList);
            unordered_map<int,set<int>> clusters;
            nodes_uCLUSTERList.push_back(clusters);
        }
        pnodes_uOCCURListMap[poi_id] = nodes_uOCCURList;
        pnodes_uCLUSTERListMap[poi_id] = nodes_uCLUSTERList;


        int p_Ni = poi.Ni; int p_Nj = poi.Nj;
        vector<int> qo_keys = poi.keywords;
        int qo_leaf = Nodes[p_Ni].gtreepath.back();
        TreeNode tnData = getGIMTreeNodeData(qo_leaf,OnlyU);
#ifdef TRACK
        cout<<"o_leaf为leaf"<<qo_leaf<<endl;
        cout<<"o_leaf中user的关键词并集："<<endl;
        printSetElements(tnData.userUKeySet);
#endif

        //// 计算 各个 poi 到 o_leaf上所有节点的border的距离
        MaterilizedDistanceMap itm; // intermediate answer, tree node -> array
        itm.clear();
        int cid, posa, posb, _min, dis;
        int o_leaf = Nodes[poi.Ni].gtreepath.back();
        itm[o_leaf].clear();
        int pos_Ni = id2idx_vertexInLeaf[p_Ni];
        int pos_Nj = id2idx_vertexInLeaf[p_Nj];
        int tn = o_leaf;
        for (int j = 0; j < GTree[tn].borders.size(); j++) {
            int dist_via_Ni = GTree[tn].mind[j * GTree[tn].leafnodes.size() + pos_Ni]+poi.dis;
            int dist_via_Nj = GTree[tn].mind[j * GTree[tn].leafnodes.size() + pos_Nj]+(poi.dist-poi.dis);
            int dist_to_border = min(dist_via_Ni, dist_via_Nj);
            itm[tn].push_back(dist_to_border);
        }
        global_itm[poi_id] = itm;

        ////获得与poi 社交文本相关的 user
        vector<int> checkin_users = poiCheckInIDList[poi_id];
        int user_textual_count = 0; int user_textual_useful_count = 0;
        for(int check_usr_id: checkin_users){
            vector<int> friends = friendshipMap[check_usr_id];
            for(int friend_id: friends){
                if(friend_id>=UserID_MaxKey) continue;
                int up_key = friend_id + UserID_MaxKey* poi_th;   //重要：每个poi对应（0,usrnum-1)

                if(upEntryMark[up_key] == true) continue;
                upEntryMark[up_key] = true;


                User u = getUserFromO2UOrgLeafData(friend_id);
                int usr_id = u.id;
                double relevance = 0;
                relevance = textRelevance(u.keywords, poi.keywords);
                //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
                if(relevance==0) continue;
                //该用户社交文本相关
#ifdef TRACKNVDBATCH
                cout<<"评估 社交文本相关用户：u"<<usr_id<<endl;
#endif
                user_textual_count ++;
                double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,u);
                double u_lcl_rk = getUserLCL_NVDFaster(u,K,a,alpha);;//getUserLCL_NVD(u,K,a,alpha);
                if(gsk_score > u_lcl_rk){

                    //更新 user-poi associating list
                    candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                    if(candidate_usr_related_store[usr_id].size()>K){
                        candidate_usr_related_store[usr_id].pop();
                    }
                    //加入candidate user set
                    BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_lcl_rk,u);
                    candidate_User.insert(ve);

                    user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                    cout<<"social_textual find candidate 用户： u"<<usr_id<<endl;
#endif
                }

            }
        }
#ifdef TRACKNVDBATCH
        cout<<"社交文本相关的 user总体个数="<<user_textual_count<<endl;
    //cout<<"完成社交文本相关的 user的剪枝，此时 candidate user size="<<candidate_id.size()<<"，具体为："<<endl;
    //printElements(candidate_id);
#endif
        //getchar();

        ////将q_o查询关键词分为 “（rare)稀有关键词”与“(none-rare)非稀有关键词” 两类,获得包含特稀疏keyword( |inv(term)|<K )的 user
        set<int> rareKeys;  set<int> nonRareKeys;
        set<int> rareKeywordRelated_usrSet;
        for(int term: qo_keys){
            int inv_olist_size = getTermOjectInvListSize(term);
            if(inv_olist_size>K){
                vector<int> usr_termRelated = getTermUserInvList(term);
                int inv_ulist_size = usr_termRelated.size();
                if(inv_ulist_size> uterm_SPARCITY_VALUE){
                    nonRareKeys.insert(term);
                    ////printf("term %d 为用户高频词汇，u_frequency=%d\n", term,inv_ulist_size);
                }
                else{ ////把一些带稀有关键词的user提前进行验证，以提升后续剪枝的效果
                    for(int user_id: usr_termRelated){
                        if(user_id>=UserID_MaxKey) continue;
                        int up_key = user_id + UserID_MaxKey* poi_th;   //重要：每个poi对应（0,usrnum-1)

                        if(upEntryMark[up_key] == true) continue;
                        upEntryMark[up_key] = true;
                        midKeywordRelated_usrSet.insert(user_id);

                        User u = getUserFromO2UOrgLeafData(user_id);
                        int usr_id = u.id;
                        double relevance = 0;
                        relevance = textRelevance(u.keywords, poi.keywords);
                        //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
                        if(relevance==0) continue;
                        //该用户社交文本相关
#ifdef TRACKNVDBATCH
                        cout<<"评估 社交文本相关用户：u"<<usr_id<<endl;
#endif
                        user_textual_count ++;
                        double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,u);
                        double u_lcl_rk = getUserLCL_NVDFaster(u,K,a,alpha);  //0
                        if(gsk_score > u_lcl_rk){

                            //更新 user-poi associating list
                            candidate_usr_related_store[user_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                            if(candidate_usr_related_store[user_id].size()>K){
                                candidate_usr_related_store[user_id].pop();
                            }
                            //加入candidate user set
                            BatchVerifyEntry ve = BatchVerifyEntry(user_id,-1,u_lcl_rk,u);
                            candidate_User.insert(ve);
                            user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                            cout<<"social_textual find candidate 用户： u"<<usr_id<<endl;
#endif
                        }

                    }
                }

            }
            else{
                rareKeys.insert(term);
                vector<int> term_related_users = getTermUserInvList(term);
                for(int usr: term_related_users)
                    rareKeywordRelated_usrSet.insert(usr);
            }

        }
        nonRareKeywordMap[poi_id] = nonRareKeys;

        set<int> poiRelated_usr; poiRelated_usr.clear();

        ////把非低频关键词对应的用户提取出来，分析他们的特征分布
        for(int term: nonRareKeys){
            vector<int> term_related_user = getTermUserInvList(term);
            for(int user_id: term_related_user){
                if(user_id>=UserID_MaxKey) continue;
                int up_key = user_id + UserID_MaxKey* poi_th;   //重要：每个poi对应（0,usrnum-1)
                if(upEntryMark[up_key] == true) continue; ////仅正对当前poi下
                poiRelated_usr.insert(user_id);
            }

        }
        for(int u_id: poiRelated_usr){
            User user = getUserFromO2UOrgLeafData(u_id);
            int u_Ni = user.Ni;
            int u_leaf = Nodes[u_Ni].gtreepath.back();
            pnodes_uOCCURListMap[poi_id][u_leaf].insert(u_id);  //更新叶节点的uOCCURList
            vector<int> u_keywords = user.keywords;
            for(int term:u_keywords){
                if(isCovered(term,nonRareKeys)){
                    for(int term2:u_keywords){
                        if(isCovered(term2,nonRareKeys)){
                            pnodes_uCLUSTERListMap[poi_id][u_leaf][term].insert(term2); //更新叶节点的uCLUSTERList
                        }

                    }
                }

            }

            //更新上层节点的uOCCURList、uCLUSTERList
            int father = GTree[u_leaf].father;
            while(father!=-1){
                pnodes_uOCCURListMap[poi_id][father].insert(u_id);
                for(int term:u_keywords){
                    if(isCovered(term,nonRareKeys)){
                        for(int term2:u_keywords){
                            if(isCovered(term2,nonRareKeys)){
                                pnodes_uCLUSTERListMap[poi_id][father][term].insert(term2);
                            }

                        }

                    }
                }

                father = GTree[father].father;
            }
        }


        ////对包含特稀疏关键词的user （这样的user通常较少）进行剪枝
        int _size = rareKeywordRelated_usrSet.size();
        //int _size2 = midKeywordRelated_usrSet.size();
        ////cout<<"包含POI下低频分布关键词的user个数="<<_size<<endl;

        for(int usr_id:rareKeywordRelated_usrSet){

            int up_key = usr_id + UserID_MaxKey* poi_th;
            //cout<<"access u"<<usr_id<<endl;

            if(upEntryMark[up_key] == true) continue;
            upEntryMark[up_key] = true;
            User u = getUserFromO2UOrgLeafData(usr_id);
            //printUsrInfo(u);

            double relevance = 0;
            relevance = textRelevance(u.keywords, poi.keywords);
            //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
            if(relevance==0) continue;
            //该用户社交文本相关
            //cout<<"评估 社交文本相关用户：u"<<usr_id<<endl;
            double gsk_score = getSKScore_o2u_phl(a,alpha,poi,u);
            double u_lcl_rk = getUserLCL_NVDFaster(u,K,a,alpha);
            if(gsk_score > u_lcl_rk){
                //更新 user-poi associating list
                candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                if(candidate_usr_related_store[usr_id].size()>K){
                    candidate_usr_related_store[usr_id].pop();
                }
                //加入candidate user set
                BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_lcl_rk,u);
                candidate_User.insert(ve);
                user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                cout<<"rare textual find candidate 用户： u"<<usr_id<<endl;
#endif
            }
        }

#ifdef TRACKNVDBATCH
        cout<<"特稀疏关键词相关的 user 总体个数="<<rareKeywordRelated_usrSet.size()<<endl;
    //cout<<"完成特稀疏关键词相关的 user的剪枝，此时 candidate user size="<<candidate_id.size()<<"，具体为："<<endl;
    //printElements(candidate_id);
#endif
        //getchar();

        ////1. q_o 所在leaf 中的关键词相关 user先加入队列考察
        set<int> user_nearSet;
        for(int term: nonRareKeys){
            vector<int> u_list = getUsrTermRelatedEntry(term,qo_leaf);
            if(u_list.size()>0){
                for(int u_id: u_list){
                    user_nearSet.insert(u_id);
                }
            }
        }



#ifdef TRACKNVDBATCH
        if(user_nearSet.size() != 0){
        cout<<"q_o 所在leaf"<<qo_leaf<<"中的关键词相关 user先加入队列考察"<<endl;
    } else{
        cout<<"q_o 所在leaf "<<qo_leaf<<"中的无关键词相关user"<<endl;
    }

        cout<<"插入user_nearSet前,Queue.size="<<Queue.size()<<endl;
#endif

        for(int usr_id: user_nearSet){
            int up_key = usr_id  + UserID_MaxKey * poi_th;
            if(upEntryMark[up_key] == true) continue;

            User user = getUserFromO2UOrgLeafData(usr_id);
            double dist = usrToPOIDistance_phl(user, poi);
            set<int> u_keySet;
            /*for(int term: user.keywords)
                u_keySet.insert(term);*/
            BatchCheckEntry entry(usr_id,poi_id, qo_leaf,dist,user);  //zhiq
            Queue.push(entry);
        }
#ifdef TRACKNVDBATCH
        cout<<"插入user_nearSet后,Queue.size="<<Queue.size()<<endl;
#endif
        //getchar();


        //设定当前level最高的Gtree Node
        int node_highest = qo_leaf;
        POIHighestNode pair(poi_id,node_highest);
        T_max.push(pair);

        poi_th ++;

    }
    initial_endTime = clock();
    initial_time = (double) (initial_endTime-initial_startTime)/CLOCKS_PER_SEC*1000;

    cout<<"总体上，包含用户下非高频分布关键词的user个数="<<midKeywordRelated_usrSet.size()<<endl;
    cout<<"完成初始步骤，开始从各个qo_leaf出发进行拓展与验证"<<endl;
    //getchar();


    ////2.从各qo_leaf出发按距离优先,进行路网拓展（在Gtree中从qo_leaf开始，自顶向上，遍历各个关键词相关节点）含有相关关键词的 user,user_leaf, user_node

    map<int,map<int,bool>> borderEntryVist;  //border_id, node_id
    map<int,bool> userEntryVist;
    while_startTime = clock();
    while(Queue.size()>0 ||T_max.size()>0){
        if(Queue.empty()){ //当前范围内已全部遍历， 需扩大路网遍历范围
#ifdef TRACKEXTEND
            cout<<"当前所有poi的检查范围内的相关user皆已被全部遍历!"<<endl;
#endif

            BatchExtend_Range(T_max,global_itm,Queue,nonRareKeywordMap,candidate_POIs, K,a,alpha);
            if(Queue.empty()) break;  //路网拓展可被终止

        }
        BatchCheckEntry topEntry = Queue.top();
        Queue.pop();

        float dist_parentBorder = topEntry.dist;
        int poi_id = topEntry.p_id;
        POI poi = candidate_POIs[poi_id];
        int poi_th = poi_Id2IdxMap[poi_id];

        if(topEntry.isUser){ //// an entry for user
            User user = topEntry.user;
            int usr_id = user.id;


            int up_key = usr_id + UserID_MaxKey* poi_th;
            //cout<<"up_key="<<up_key<<endl;
            if(upEntryMark[up_key])
                continue;
            else
                upEntryMark[up_key] = true;


            double u_SB_rk = getUserLCL_NVDFaster(user,K,a,alpha);

            ////cout << "getUserLCL_NVD runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;

            double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,user);
            ////cout << "getGSKScore_o2u runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;
#ifdef WhileTRACK
            cout<<"------DEQUEUE: user entry: u"<<usr_id<<"for p"<<poi_id<<"------"<<endl;
            cout<<"u_SB_rk="<<u_SB_rk<<endl;
            cout<<"gsk_score="<<gsk_score<<endl;


#endif

            if(gsk_score > u_SB_rk){
                if(candidate_usr_related_store.count(usr_id)){
                    candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0)); //jinsg
                    if(candidate_usr_related_store[usr_id].size()>K){
                        candidate_usr_related_store[usr_id].pop();
                    }
                }
                else{
                    priority_queue<ResultCurrent> list;
                    list.push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                    candidate_usr_related_store[usr_id] = list;
                }

                //加入candidate user set
                BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_SB_rk,user);
                candidate_User.insert(ve);

            }
        }
        else{  //// an entry for border
            int border_id = topEntry.id;
            int poi_id = topEntry.p_id;
            int poi_th = poi_Id2IdxMap[poi_id];
            int node_id =topEntry.node_id;
            int np_key = node_id + GTree.size()*poi_th;
            int bEntry_key =  border_id +  node_id*VertexNum  + poi_th*inner_size;


            if(borderMark[bEntry_key]==true) continue;
            borderMark[bEntry_key] = true;

            set<int> Keys = topEntry.keys_cover;
            //若为叶节点，直接快速评估叶节点下的各个user
            if(GTree[node_id].isleaf){
#ifdef TRACK
                cout<<"b"<<border_id<<", 位于leaf"<<node_id<<endl;
#endif
                if(nodepMark[np_key]!=true) {
                    nodepMark[np_key] =true;  //标记该node相对当前poi已被处理
                    //continue;
                    set<int> user_within;
                    for(int term: Keys){
                        vector<int> u_list = getUsrTermRelatedEntry(term,node_id);
                        if(u_list.size()>0){
                            for(int u_id: u_list){
                                user_within.insert(u_id);
                            }
                        }
                    }

                    for(int usr_id: user_within){
                        //if(userMarkMap[poi_id][usr_id] == true) continue;
                        int up_key = usr_id + UserID_MaxKey *poi_th;
                        if(upEntryMark[up_key] == true) continue;

                        User user = getUserFromO2UOrgLeafData(usr_id);
                        double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,user);
                        double u_SB_rk = getUserLCL_NVDFaster(user,K,a,alpha);
                        if(gsk_score > u_SB_rk){
                            if(candidate_usr_related_store.count(usr_id)){
                                candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0)); //jinsg
                                if(candidate_usr_related_store[usr_id].size()>K){
                                    candidate_usr_related_store[usr_id].pop();
                                }
                            }
                            else{
                                priority_queue<ResultCurrent> list;
                                list.push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                                candidate_usr_related_store[usr_id] = list;
                            }

                            //加入candidate user set
                            BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_SB_rk,user);
                            candidate_User.insert(ve);

                        }

                    }
                }


            }
            else{  //非叶子节点

#ifdef   AGGRESSIVE_I
                int _size = pnodes_uOCCURListMap[poi_id][node_id].size();
                if(_size < NODE_SPARCITY_VALUE_HEURIS){//&& dist_parentBorder>LONG_TERM
                    for(int usr_id: pnodes_uOCCURListMap[poi_id][node_id]){
                        int up_key = usr_id + UserID_MaxKey *poi_th;
                        if(upEntryMark[up_key] == true)
                            continue;
                        else
                            upEntryMark[up_key] = true;
                        //cout<<"u"<<usr_id<<"被提前评估！"<<endl;
                        User user = getUserFromO2UOrgLeafData(usr_id);
                        //double dist = usrToPOIDistance_phl(user, poi);
                        double u_SB_rk = getUserLCL_NVD(user,K,a,alpha);//0;

                        double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,user);

                        if(gsk_score > u_SB_rk){
                            if(candidate_usr_related_store.count(usr_id)){
                                candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0)); //jinsg
                                if(candidate_usr_related_store[usr_id].size()>K){
                                    candidate_usr_related_store[usr_id].pop();
                                }
                            }
                            else{
                                priority_queue<ResultCurrent> list;
                                list.push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                                candidate_usr_related_store[usr_id] = list;
                            }

                            //加入candidate user set
                            BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_SB_rk,user);
                            candidate_User.insert(ve);

                        }

                    }
                    nodepMark[np_key] =true;  //标记该node相对当前poi已被处理
                    continue;
                }
#endif


                set<int> Keys;

                Keys = getBorder_SB_BatchNVD_Cluster(topEntry,K,a, alpha,pnodes_uCLUSTERListMap);

                //cout<<"remained (user) Keywords"; printSetElements(Keys);
                if(Keys.size() ==0) {
#ifdef TRACK
                    cout<<"所有其后继被 ****prune！"<<endl;
#endif
                    continue;
                }
#ifdef TRACK
                cout<<"其后继**不能**被prune！"<<endl;
#endif


                //获得父节点下所有关键词相关的孩子节点(line 22-25)
#ifdef WhileTRACK
                cout<<"b"<<border_id<<", 位于node"<<node_id<<endl;
#endif
                set<int> child_set;
                for(int term:Keys){
                    vector<int> related_childs = getUsrTermRelatedEntry(term,node_id);
                    if(related_childs.size()>0){
                        for(int cn_id: related_childs){
                            child_set.insert(cn_id);
                        }
                    }
                }
#ifdef TRACK
                cout<<"关键词相关的后继："; printSetElements(child_set);
#endif

                int parent_border_th = topEntry.b_th;
                int posa = parent_border_th;
                for(int cn_id: child_set){
                    for(int b_th=0;b_th<GTree[cn_id].borders.size();b_th++){

                        int child_border_id = GTree[cn_id].borders[b_th];
                        float dist_childBorder = getDistance_phl(child_border_id,poi);  //jins

                        int border_key =  child_border_id +  cn_id*VertexNum  + poi_th*inner_size;
                        if(borderMark[border_key] == true) continue;
#ifdef TRACK
                        cout<<"child_border_id: v"<<child_border_id<<", poi_th"<<poi_th<<endl;
                        cout<<"border_key:"<<border_key<<endl;
#endif
                        BatchCheckEntry entry(child_border_id,poi_id, cn_id,Keys,dist_childBorder,b_th);
                        Queue.push(entry);
#ifdef TRACK
                        cout<<"在main中加入border entry:"; entry.printRlt();
                        cout<<"Queue.size="<<Queue.size()<<endl;
#endif


                    }
                }

            }
        }

    }
    filter_endTime = clock();
    endTime = clock();
    delete [] upEntryMark;

    BatchFastFilterResults bfOutput;
    bfOutput.candidate_User= candidate_User;
    bfOutput.candidate_usr_related_store = candidate_usr_related_store;
    bfOutput.batch_results = batch_results;
    bfOutput.candidate_POIs = candidate_POIs;
    printBatchRkGSKQResults(stores,batch_results);


    printf("初始化步骤耗时：%f ms !\n", initial_time);
    cout<<"fast filter阶段总共发现"<<candidate_User.size()<<"个候选用户！"<<endl;
    cout<<"Fast Filter 阶段 用时："<< (double)(filter_endTime - filter_startTime) / CLOCKS_PER_SEC*1000  <<
        "ms (其中 while loop用时:"<< (double)(filter_endTime - while_startTime) / CLOCKS_PER_SEC*1000<<"ms)"<<endl;
    //cout<<"Verification 阶段 用时："<< (double)(verification_endTime - verification_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;

    //cout<<"user_textual_count="<<user_textual_count<<",user_textual_useful_count="<<user_textual_useful_count<<endl;
    cout<<"---------Batch Fast Filter 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC*1000  << "ms-----------" << endl;


    return bfOutput;
}

BatchFastFilterResults Fast_Batch_NVDFilter(vector<int> stores, int K, float a, float alpha) {  //gangcai

    cout<<"------------Running Fast_Batch_NVDFilter!------------"<<endl;
    clock_t startTime, endTime;
    clock_t filter_startTime,while_startTime, filter_endTime;
    clock_t verification_startTime, verification_endTime;
    clock_t initial_startTime, initial_endTime;

    //vector<int> RkGSKQ_results;
    FastBatch_filterResults batch_results;



    unordered_map<int, POI> candidate_POIs;
    unordered_map<int, int> poi_Id2IdxMap;
    unordered_map<int, set<int>>  nonRareKeywordMap;


    //候选验证用户
    set<BatchVerifyEntry> candidate_User;  //candidate users
    unordered_map<int, priority_queue<ResultCurrent>> candidate_usr_related_store;   // user_id, {related poi1,2,....}
    unordered_map<int, priority_queue<ResultCurrent>> candidate_usrLeaf_related_store;   // user_id, {related poi1,2,....}


    //border已被评估的标志向量
    int node_size = GTree.size();
    int inner_size = node_size*VertexNum;  // func(bEntry_key)= v_id+v_num*node_id + v_num*node_num*poi_th
    long bEntry_MaxKey =  inner_size * stores.size();
    unordered_map<int,bool> borderMark;



    //user已被评估的标志向量表
    int upEntry_MaxKey =  UserID_MaxKey * stores.size();
    //unordered_map<int,bool> upEntryMark;
    bool* upEntryMark = new bool[upEntry_MaxKey];
    for(int i=0;i<upEntry_MaxKey;i++){
        upEntryMark[i] = false;
    }

    bool userMark[UserID_MaxKey];
    memset(userMark,false, sizeof(userMark));

    //标记gtree node、 poi 对的标志向量
    int np_maxKey = GTree.size()*stores.size();
    bool nodepMark[np_maxKey];
    memset(nodepMark,false, sizeof(nodepMark));


    startTime = clock();
    filter_startTime = clock();
    //距离优先队列

    ////初始化gtree node的 uOCCURList以及 uCLUSTERList
    unordered_map<int,UOCCURList> pnodes_uOCCURListMap;
    unordered_map<int,UCLUSTERList> pnodes_uCLUSTERListMap;
    //unordered_map<int,UOCCURListMap> pnodes_uOCCURListMap;
    //unordered_map<int,UCLUSTERListMap> pnodes_uCLUSTERListMap;



    priority_queue<BatchCheckEntry> Queue;
    priority_queue<POIHighestNode> T_max;


    int root = 0;
    //// 1. 获取各个poi信息, 并计算各个 poi 到 o_leaf上所有节点的border的距离
    unordered_map<int, unordered_map<int, vector<int>>> global_itm;
    int poi_th = 0;
    set<int> midKeywordRelated_usrSet; double initial_time = -1;
    initial_startTime = clock();
    for(int poi_id: stores){
        poi_Id2IdxMap[poi_id] = poi_th;


        POI poi = getPOIFromO2UOrgLeafData(poi_id);

        candidate_POIs[poi_id] = poi;

        ////printPOIInfo(poi);


        //初始化并构建该poi对应的 nodes_uOCCURList 及 nodes_uCLUSTERList
        //UOCCURListMap nodes_uOCCURList;
        //UCLUSTERListMap nodes_uCLUSTERList;  //vector_index= node_id,  term_id,
        UOCCURList nodes_uOCCURList;
        UCLUSTERList nodes_uCLUSTERList;  //vector_index= node_id,  term_id,
        for(int i=0;i<GTree.size();i++){
            set<int> node_occurList;
            nodes_uOCCURList.push_back(node_occurList);
            unordered_map<int,set<int>> clusters;
            nodes_uCLUSTERList.push_back(clusters);
        }
        pnodes_uOCCURListMap[poi_id] = nodes_uOCCURList;
        pnodes_uCLUSTERListMap[poi_id] = nodes_uCLUSTERList;


        int p_Ni = poi.Ni; int p_Nj = poi.Nj;
        vector<int> qo_keys = poi.keywords;
        int qo_leaf = Nodes[p_Ni].gtreepath.back();
        TreeNode tnData = getGIMTreeNodeData(qo_leaf,OnlyU);
#ifdef TRACK
        cout<<"o_leaf为leaf"<<qo_leaf<<endl;
        cout<<"o_leaf中user的关键词并集："<<endl;
        printSetElements(tnData.userUKeySet);
#endif

        //// 计算 各个 poi 到 o_leaf上所有节点的border的距离
        MaterilizedDistanceMap itm; // intermediate answer, tree node -> array
        itm.clear();
        int cid, posa, posb, _min, dis;
        int o_leaf = Nodes[poi.Ni].gtreepath.back();
        itm[o_leaf].clear();
        int pos_Ni = id2idx_vertexInLeaf[p_Ni];
        int pos_Nj = id2idx_vertexInLeaf[p_Nj];
        int tn = o_leaf;
        for (int j = 0; j < GTree[tn].borders.size(); j++) {
            int dist_via_Ni = GTree[tn].mind[j * GTree[tn].leafnodes.size() + pos_Ni]+poi.dis;
            int dist_via_Nj = GTree[tn].mind[j * GTree[tn].leafnodes.size() + pos_Nj]+(poi.dist-poi.dis);
            int dist_to_border = min(dist_via_Ni, dist_via_Nj);
            itm[tn].push_back(dist_to_border);
        }
        global_itm[poi_id] = itm;

        ////获得与poi 社交文本相关的 user
        vector<int> checkin_users = poiCheckInIDList[poi_id];
        int user_textual_count = 0; int user_textual_useful_count = 0;
        for(int check_usr_id: checkin_users){
            vector<int> friends = friendshipMap[check_usr_id];
            for(int friend_id: friends){
                if(friend_id>=UserID_MaxKey) continue;
                int up_key = friend_id + UserID_MaxKey* poi_th;   //重要：每个poi对应（0,usrnum-1)

                if(upEntryMark[up_key] == true) continue;
                upEntryMark[up_key] = true;


                User u = getUserFromO2UOrgLeafData(friend_id);
                int usr_id = u.id;
                double relevance = 0;
                relevance = textRelevance(u.keywords, poi.keywords);
                //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
                if(relevance==0) continue;
                //该用户社交文本相关
#ifdef TRACKNVDBATCH
                cout<<"评估 社交文本相关用户：u"<<usr_id<<endl;
#endif
                user_textual_count ++;
                double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,u);
                double u_lcl_rk = getUserLCL_NVDFaster(u,K,a,alpha);;//getUserLCL_NVD(u,K,a,alpha);
                if(gsk_score > u_lcl_rk){

                    //更新 user-poi associating list
                    candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                    if(candidate_usr_related_store[usr_id].size()>K){
                        candidate_usr_related_store[usr_id].pop();
                    }
                    //加入candidate user set
                    BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_lcl_rk,u);
                    candidate_User.insert(ve);

                    user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                    cout<<"social_textual find candidate 用户： u"<<usr_id<<endl;
#endif
                }

            }
        }
#ifdef TRACKNVDBATCH
        cout<<"社交文本相关的 user总体个数="<<user_textual_count<<endl;
    //cout<<"完成社交文本相关的 user的剪枝，此时 candidate user size="<<candidate_id.size()<<"，具体为："<<endl;
    //printElements(candidate_id);
#endif
        //getchar();

        ////将q_o查询关键词分为 “（rare)稀有关键词”与“(none-rare)非稀有关键词” 两类,获得包含特稀疏keyword( |inv(term)|<K )的 user
        set<int> rareKeys;  set<int> nonRareKeys;
        set<int> rareKeywordRelated_usrSet;
        for(int term: qo_keys){
            int inv_olist_size = getTermOjectInvListSize(term);
            if(inv_olist_size>K){
                vector<int> usr_termRelated = getTermUserInvList(term);
                int inv_ulist_size = usr_termRelated.size();
                if(inv_ulist_size> uterm_SPARCITY_VALUE){
                    nonRareKeys.insert(term);
                    ////printf("term %d 为用户高频词汇，u_frequency=%d\n", term,inv_ulist_size);
                }
                else{ ////把一些带稀有关键词的user提前进行验证，以提升后续剪枝的效果
                    for(int user_id: usr_termRelated){
                        if(user_id>=UserID_MaxKey) continue;
                        int up_key = user_id + UserID_MaxKey* poi_th;   //重要：每个poi对应（0,usrnum-1)

                        if(upEntryMark[up_key] == true) continue;
                        upEntryMark[up_key] = true;
                        midKeywordRelated_usrSet.insert(user_id);

                        User u = getUserFromO2UOrgLeafData(user_id);
                        int usr_id = u.id;
                        double relevance = 0;
                        relevance = textRelevance(u.keywords, poi.keywords);
                        //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
                        if(relevance==0) continue;
                        //该用户社交文本相关
#ifdef TRACKNVDBATCH
                        cout<<"评估 社交文本相关用户：u"<<usr_id<<endl;
#endif
                        user_textual_count ++;
                        double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,u);
                        double u_lcl_rk = getUserLCL_NVDFaster(u,K,a,alpha);  //0
                        if(gsk_score > u_lcl_rk){

                            //更新 user-poi associating list
                            candidate_usr_related_store[user_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                            if(candidate_usr_related_store[user_id].size()>K){
                                candidate_usr_related_store[user_id].pop();
                            }
                            //加入candidate user set
                            BatchVerifyEntry ve = BatchVerifyEntry(user_id,-1,u_lcl_rk,u);
                            candidate_User.insert(ve);
                            user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                            cout<<"social_textual find candidate 用户： u"<<usr_id<<endl;
#endif
                        }

                    }
                }

            }
            else{
                rareKeys.insert(term);
                vector<int> term_related_users = getTermUserInvList(term);
                for(int usr: term_related_users)
                    rareKeywordRelated_usrSet.insert(usr);
            }

        }
        nonRareKeywordMap[poi_id] = nonRareKeys;

        set<int> poiRelated_usr; poiRelated_usr.clear();

        ////把非低频关键词对应的用户提取出来，分析他们的特征分布
        for(int term: nonRareKeys){
            vector<int> term_related_user = getTermUserInvList(term);
            for(int user_id: term_related_user){
                if(user_id>=UserID_MaxKey) continue;
                int up_key = user_id + UserID_MaxKey* poi_th;   //重要：每个poi对应（0,usrnum-1)
                if(upEntryMark[up_key] == true) continue; ////仅正对当前poi下
                poiRelated_usr.insert(user_id);
            }

        }
        for(int u_id: poiRelated_usr){
            User user = getUserFromO2UOrgLeafData(u_id);
            int u_Ni = user.Ni;
            int u_leaf = Nodes[u_Ni].gtreepath.back();
            pnodes_uOCCURListMap[poi_id][u_leaf].insert(u_id);  //更新叶节点的uOCCURList
            vector<int> u_keywords = user.keywords;
            for(int term:u_keywords){
                if(isCovered(term,nonRareKeys)){
                    for(int term2:u_keywords){
                        if(isCovered(term2,nonRareKeys)){
                            pnodes_uCLUSTERListMap[poi_id][u_leaf][term].insert(term2); //更新叶节点的uCLUSTERList
                        }

                    }
                }

            }

            //更新上层节点的uOCCURList、uCLUSTERList
            int father = GTree[u_leaf].father;
            while(father!=-1){
                pnodes_uOCCURListMap[poi_id][father].insert(u_id);
                for(int term:u_keywords){
                    if(isCovered(term,nonRareKeys)){
                        for(int term2:u_keywords){
                            if(isCovered(term2,nonRareKeys)){
                                pnodes_uCLUSTERListMap[poi_id][father][term].insert(term2);
                            }

                        }

                    }
                }

                father = GTree[father].father;
            }
        }


        ////对包含特稀疏关键词的user （这样的user通常较少）进行剪枝
        int _size = rareKeywordRelated_usrSet.size();
        //int _size2 = midKeywordRelated_usrSet.size();
        ////cout<<"包含POI下低频分布关键词的user个数="<<_size<<endl;

        for(int usr_id:rareKeywordRelated_usrSet){

            int up_key = usr_id + UserID_MaxKey* poi_th;
            //cout<<"access u"<<usr_id<<endl;

            if(upEntryMark[up_key] == true) continue;
            upEntryMark[up_key] = true;
            User u = getUserFromO2UOrgLeafData(usr_id);
            //printUsrInfo(u);

            double relevance = 0;
            relevance = textRelevance(u.keywords, poi.keywords);
            //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
            if(relevance==0) continue;
            //该用户社交文本相关
            //cout<<"评估 社交文本相关用户：u"<<usr_id<<endl;
            double gsk_score = getSKScore_o2u_phl(a,alpha,poi,u);
            double u_lcl_rk = getUserLCL_NVDFaster(u,K,a,alpha);
            if(gsk_score > u_lcl_rk){
                //更新 user-poi associating list
                candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                if(candidate_usr_related_store[usr_id].size()>K){
                    candidate_usr_related_store[usr_id].pop();
                }
                //加入candidate user set
                BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_lcl_rk,u);
                candidate_User.insert(ve);
                user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                cout<<"rare textual find candidate 用户： u"<<usr_id<<endl;
#endif
            }
        }

#ifdef TRACKNVDBATCH
        cout<<"特稀疏关键词相关的 user 总体个数="<<rareKeywordRelated_usrSet.size()<<endl;
    //cout<<"完成特稀疏关键词相关的 user的剪枝，此时 candidate user size="<<candidate_id.size()<<"，具体为："<<endl;
    //printElements(candidate_id);
#endif
        //getchar();

        ////1. q_o 所在leaf 中的关键词相关 user先加入队列考察
        set<int> user_nearSet;
        for(int term: nonRareKeys){
            vector<int> u_list = getUsrTermRelatedEntry(term,qo_leaf);
            if(u_list.size()>0){
                for(int u_id: u_list){
                    user_nearSet.insert(u_id);
                }
            }
        }



#ifdef TRACKNVDBATCH
        if(user_nearSet.size() != 0){
        cout<<"q_o 所在leaf"<<qo_leaf<<"中的关键词相关 user先加入队列考察"<<endl;
    } else{
        cout<<"q_o 所在leaf "<<qo_leaf<<"中的无关键词相关user"<<endl;
    }

        cout<<"插入user_nearSet前,Queue.size="<<Queue.size()<<endl;
#endif

        for(int usr_id: user_nearSet){
            int up_key = usr_id  + UserID_MaxKey * poi_th;
            if(upEntryMark[up_key] == true) continue;

            User user = getUserFromO2UOrgLeafData(usr_id);
            double dist = usrToPOIDistance_phl(user, poi);
            set<int> u_keySet;
            /*for(int term: user.keywords)
                u_keySet.insert(term);*/
            BatchCheckEntry entry(usr_id,poi_id, qo_leaf,dist,user);  //zhiq
            Queue.push(entry);
        }
#ifdef TRACKNVDBATCH
        cout<<"插入user_nearSet后,Queue.size="<<Queue.size()<<endl;
#endif
        //getchar();


        //设定当前level最高的Gtree Node
        int node_highest = qo_leaf;
        POIHighestNode pair(poi_id,node_highest);
        T_max.push(pair);

        poi_th ++;

    }
    initial_endTime = clock();
    initial_time = (double) (initial_endTime-initial_startTime)/CLOCKS_PER_SEC*1000;

    cout<<"总体上，包含用户下非高频分布关键词的user个数="<<midKeywordRelated_usrSet.size()<<endl;
    cout<<"完成初始步骤，开始从各个qo_leaf出发进行拓展与验证"<<endl;
    //getchar();


    ////2.从各qo_leaf出发按距离优先,进行路网拓展（在Gtree中从qo_leaf开始，自顶向上，遍历各个关键词相关节点）含有相关关键词的 user,user_leaf, user_node

    map<int,map<int,bool>> borderEntryVist;  //border_id, node_id
    map<int,bool> userEntryVist;
    while_startTime = clock();
    while(Queue.size()>0 ||T_max.size()>0){
        if(Queue.empty()){ //当前范围内已全部遍历， 需扩大路网遍历范围
#ifdef TRACKEXTEND
            cout<<"当前所有poi的检查范围内的相关user皆已被全部遍历!"<<endl;
#endif

            BatchExtend_Range(T_max,global_itm,Queue,nonRareKeywordMap,candidate_POIs, K,a,alpha);
            if(Queue.empty()) break;  //路网拓展可被终止

        }
        BatchCheckEntry topEntry = Queue.top();
        Queue.pop();

        float dist_parentBorder = topEntry.dist;
        int poi_id = topEntry.p_id;
        POI poi = candidate_POIs[poi_id];
        int poi_th = poi_Id2IdxMap[poi_id];

        if(topEntry.isUser){ //// an entry for user
            User user = topEntry.user;
            int usr_id = user.id;


            int up_key = usr_id + UserID_MaxKey* poi_th;
            //cout<<"up_key="<<up_key<<endl;
            if(upEntryMark[up_key])
                continue;
            else
                upEntryMark[up_key] = true;


            double u_SB_rk = getUserLCL_NVDFaster(user,K,a,alpha);

            ////cout << "getUserLCL_NVD runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;

            double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,user);
            ////cout << "getGSKScore_o2u runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;
#ifdef WhileTRACK
            cout<<"------DEQUEUE: user entry: u"<<usr_id<<"for p"<<poi_id<<"------"<<endl;
            cout<<"u_SB_rk="<<u_SB_rk<<endl;
            cout<<"gsk_score="<<gsk_score<<endl;


#endif

            if(gsk_score > u_SB_rk){
                if(candidate_usr_related_store.count(usr_id)){
                    candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0)); //jinsg
                    if(candidate_usr_related_store[usr_id].size()>K){
                        candidate_usr_related_store[usr_id].pop();
                    }
                }
                else{
                    priority_queue<ResultCurrent> list;
                    list.push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                    candidate_usr_related_store[usr_id] = list;
                }

                //加入candidate user set
                BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_SB_rk,user);
                candidate_User.insert(ve);

            }
        }
        else{  //// an entry for border
            int border_id = topEntry.id;
            int poi_id = topEntry.p_id;
            int poi_th = poi_Id2IdxMap[poi_id];
            int node_id =topEntry.node_id;
            int np_key = node_id + GTree.size()*poi_th;
            int bEntry_key =  border_id +  node_id*VertexNum  + poi_th*inner_size;


            if(borderMark[bEntry_key]==true) continue;
            borderMark[bEntry_key] = true;

            set<int> Keys = topEntry.keys_cover;
            //若为叶节点，直接快速评估叶节点下的各个user
            if(GTree[node_id].isleaf){
#ifdef TRACK
                cout<<"b"<<border_id<<", 位于leaf"<<node_id<<endl;
#endif
                if(nodepMark[np_key]!=true) {
                    nodepMark[np_key] =true;  //标记该node相对当前poi已被处理
                    //continue;
                    set<int> user_within;
                    for(int term: Keys){
                        vector<int> u_list = getUsrTermRelatedEntry(term,node_id);
                        if(u_list.size()>0){
                            for(int u_id: u_list){
                                user_within.insert(u_id);
                            }
                        }
                    }

                    for(int usr_id: user_within){
                        //if(userMarkMap[poi_id][usr_id] == true) continue;
                        int up_key = usr_id + UserID_MaxKey *poi_th;
                        if(upEntryMark[up_key] == true) continue;

                        User user = getUserFromO2UOrgLeafData(usr_id);
                        double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,user);
                        double u_SB_rk = getUserLCL_NVDFaster(user,K,a,alpha);
                        if(gsk_score > u_SB_rk){
                            if(candidate_usr_related_store.count(usr_id)){
                                candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0)); //jinsg
                                if(candidate_usr_related_store[usr_id].size()>K){
                                    candidate_usr_related_store[usr_id].pop();
                                }
                            }
                            else{
                                priority_queue<ResultCurrent> list;
                                list.push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                                candidate_usr_related_store[usr_id] = list;
                            }

                            //加入candidate user set
                            BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_SB_rk,user);
                            candidate_User.insert(ve);

                        }

                    }
                }


            }
            else{  //非叶子节点

#ifdef   AGGRESSIVE_I
                int _size = pnodes_uOCCURListMap[poi_id][node_id].size();
                if(_size < NODE_SPARCITY_VALUE_HEURIS){//&& dist_parentBorder>LONG_TERM
                    for(int usr_id: pnodes_uOCCURListMap[poi_id][node_id]){
                        int up_key = usr_id + UserID_MaxKey *poi_th;
                        if(upEntryMark[up_key] == true)
                            continue;
                        else
                            upEntryMark[up_key] = true;
                        //cout<<"u"<<usr_id<<"被提前评估！"<<endl;
                        User user = getUserFromO2UOrgLeafData(usr_id);
                        //double dist = usrToPOIDistance_phl(user, poi);
                        double u_SB_rk = getUserLCL_NVD(user,K,a,alpha);//0;

                        double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,user);

                        if(gsk_score > u_SB_rk){
                            if(candidate_usr_related_store.count(usr_id)){
                                candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0)); //jinsg
                                if(candidate_usr_related_store[usr_id].size()>K){
                                    candidate_usr_related_store[usr_id].pop();
                                }
                            }
                            else{
                                priority_queue<ResultCurrent> list;
                                list.push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                                candidate_usr_related_store[usr_id] = list;
                            }

                            //加入candidate user set
                            BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_SB_rk,user);
                            candidate_User.insert(ve);

                        }

                    }
                    nodepMark[np_key] =true;  //标记该node相对当前poi已被处理
                    continue;
                }
#endif


                set<int> Keys;

                Keys = getBorder_SB_BatchNVD_Cluster(topEntry,K,a, alpha,pnodes_uCLUSTERListMap);

                //cout<<"remained (user) Keywords"; printSetElements(Keys);
                if(Keys.size() ==0) {
#ifdef TRACK
                    cout<<"所有其后继被 ****prune！"<<endl;
#endif
                    continue;
                }
#ifdef TRACK
                cout<<"其后继**不能**被prune！"<<endl;
#endif


                //获得父节点下所有关键词相关的孩子节点(line 22-25)
#ifdef WhileTRACK
                cout<<"b"<<border_id<<", 位于node"<<node_id<<endl;
#endif
                set<int> child_set;
                for(int term:Keys){
                    vector<int> related_childs = getUsrTermRelatedEntry(term,node_id);
                    if(related_childs.size()>0){
                        for(int cn_id: related_childs){
                            child_set.insert(cn_id);
                        }
                    }
                }
#ifdef TRACK
                cout<<"关键词相关的后继："; printSetElements(child_set);
#endif

                int parent_border_th = topEntry.b_th;
                int posa = parent_border_th;
                for(int cn_id: child_set){
                    for(int b_th=0;b_th<GTree[cn_id].borders.size();b_th++){

                        int child_border_id = GTree[cn_id].borders[b_th];
                        float dist_childBorder = getDistance_phl(child_border_id,poi);  //jins

                        int border_key =  child_border_id +  cn_id*VertexNum  + poi_th*inner_size;
                        if(borderMark[border_key] == true) continue;
#ifdef TRACK
                        cout<<"child_border_id: v"<<child_border_id<<", poi_th"<<poi_th<<endl;
                        cout<<"border_key:"<<border_key<<endl;
#endif
                        BatchCheckEntry entry(child_border_id,poi_id, cn_id,Keys,dist_childBorder,b_th);
                        Queue.push(entry);
#ifdef TRACK
                        cout<<"在main中加入border entry:"; entry.printRlt();
                        cout<<"Queue.size="<<Queue.size()<<endl;
#endif


                    }
                }

            }
        }

    }
    filter_endTime = clock();
    endTime = clock();
    delete [] upEntryMark;

    BatchFastFilterResults bfOutput;
    bfOutput.candidate_User= candidate_User;
    bfOutput.candidate_usr_related_store = candidate_usr_related_store;
    bfOutput.batch_results = batch_results;
    bfOutput.candidate_POIs = candidate_POIs;
    //printBatchRkGSKQResults(stores, batch_results);


    printf("初始化步骤耗时：%f ms !\n", initial_time);
    cout<<"fast filter阶段总共发现"<<candidate_User.size()<<"个候选用户！"<<endl;
    cout<<"Fast Filter 阶段 用时："<< (double)(filter_endTime - filter_startTime) / CLOCKS_PER_SEC*1000  <<
        "ms (其中 while loop用时:"<< (double)(filter_endTime - while_startTime) / CLOCKS_PER_SEC*1000<<"ms)"<<endl;
    //cout<<"Verification 阶段 用时："<< (double)(verification_endTime - verification_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;

    //cout<<"user_textual_count="<<user_textual_count<<",user_textual_useful_count="<<user_textual_useful_count<<endl;
    cout<<"---------Batch Fast Filter 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC*1000  << "ms-----------" << endl;


    return bfOutput;
}


// The algorithms below are for Multiple RkGSKQ

int Naive_MRkGSKQ_givenKeys(vector<int> stores, vector<int> keywords, int Qk, float a, float alpha){
    TIME_TICK_START

    vector<vector<ResultDetail>> MRkGSKQ_results;
    for(int poi_id: stores){
        cout<<"for p"<<poi_id<<endl;
        int size = 0;
        int loc = POIs[poi_id].Ni;
        //vector<int> keywords = keywords;

        vector<ResultDetail> results;
        for (User usr: Users) {
            int usr_id = usr.id;
            //usr_id = 42419;

            //double distance = getDistance(POIs[poi_id].Nj, usr);
            double distance =  usrToPOIDistance(usr,POIs[poi_id]);//min(d1,d2);
            double relevance = 0;
            relevance = textRelevance(usr.keywords, keywords);
            //relevance = 1;
            if (relevance == 0) continue;

            double social_Score = 0;

            social_Score = reCalcuSocialScore(usr.id,poiCheckInIDList[poi_id]);  //gazi
            double social_Score2 = calcSocialScore(usr.id,poi_id);


            double influence;
            if(social_Score == 0)
                influence =1;
            else
                influence =1 + pow(a,social_Score);


            double social_textual_score =  alpha*influence + (1-alpha)*relevance;
            //double score = GSKScoreFunction(distance,usr.keywords,keywords,usr.id,0, poiCheckInIDList[0], a);
            double score = social_textual_score / (1+distance);
            //score = getGSKScore(a, POIs[poi_id].Nj, usr.id, usr.keywords, poiCheckInIDList[poi_id]);
            double simD = distance; double simT = relevance; double simS = influence;

            query usr_q = transformToQ(usr);

            if(Users[usr.id].topkScore_current>score)
                continue;
#ifdef LV
            double Rk_usr = topkSDijkstra_verify_usr_memory(usr_q,Qk,a,alpha,score,0.0).topkScore;
            Users[usr.id].topkScore_current = max(Users[usr.id].topkScore_current,Rk_usr);
#endif
#ifndef LV
            double usr_lcl_rt = 0.0;
            double Rk_usr = TkGSKQ_bottom2up_verify_memory(usr.id,Qk,a,alpha,score,usr_lcl_rt).top().score;
            Users[usr.id].topkScore_current = max(Users[usr.id].topkScore_current,Rk_usr);
#endif

            //if (score > Rk_usr && distance < 10000) {   // 10km范围内！(这个条件会减小结果)
            if (score > Rk_usr) {
                size++;
                ResultDetail sd(usr.id,influence,distance,relevance,score,Rk_usr);
                results.push_back(sd);
                //cout<<"u"<<usr.id<<",";

            }
        }
        MRkGSKQ_results.push_back(results);
        cout << "RkGSKQ: size=" << size << endl;
        printResults(results);

    }

    //cout<<endl<<"for loc:"<<loc<<",";
    //cout<<"texual:";printQueryKeywords(keywords);
    for(int i=0; i<stores.size();i++){
        cout <<"store"<<stores[i]<< "'s RkGSKQ: size=" << MRkGSKQ_results[i].size() << endl;
        printResults(MRkGSKQ_results[i]);
        cout <<"---------------------------------------------------"<<endl;
    }



    TIME_TICK_END
    TIME_TICK_PRINT("Naive_MRkGSKQ running time:")
    return 0;

}


int Naive_MRkGSKQ(vector<int> stores, int Qk, float a, float alpha){
    cout<<"------------------------Naive_MRkGSKQ begin!....---------------------------"<<endl;
    TIME_TICK_START

    vector<vector<ResultDetail>> MRkGSKQ_results;
    for(int poi_id: stores){
        cout<<"for p"<<poi_id<<endl;
        int size = 0;
        POI p;
#ifdef  DiskAccess
        p = getPOIFromO2UOrgLeafData(poi_id);
#else
        p = POIs[poi_id];
#endif
        int loc = p.Ni;
        //vector<int> keywords = keywords;

        vector<ResultDetail> results;
        for (int i=0;i<Users.size();i++) {

#ifdef            NAVIE_DEBUG
            if(i==106&& poi_id==57){
                cout<<"find u"<<i<<",o"<<poi_id<<" evaluate!"<<endl;
                getchar();
            }
#endif

            User u;
#ifdef  DiskAccess
            u = getUserFromO2UOrgLeafData(i);
#else
            u = Users[i];
#endif
            double distance =  usrToPOIDistance(u,p);//min(d1,d2);
            double relevance = 0;
            relevance = textRelevance(u.keywords, p.keywords);
            //relevance = 1;
            if (relevance == 0) continue;

            double social_Score = 0;

            //social_Score = reCalcuSocialScore(usr.id,poiCheckInIDList[poi_id]);
            social_Score = calcSocialScore(u.id,poi_id);

            double influence;
            if(social_Score == 0)
                influence =1;
            else
                influence =1 + pow(a,social_Score);


            double social_textual_score =  alpha*influence + (1-alpha)*relevance;
            //double score = GSKScoreFunction(distance,usr.keywords,keywords,usr.id,0, poiCheckInIDList[0], a);
            double score = social_textual_score / (1+distance);
            //score = getGSKScore(a, POIs[poi_id].Nj, usr.id, usr.keywords, poiCheckInIDList[poi_id]);
            double simD = distance; double simT = relevance; double simS = influence;

            query usr_q = transformToQ(u);

            if(Users[i].topkScore_current>score)
                continue;
#ifdef LV
            #ifndef DiskAccess
            double Rk_usr = topkSDijkstra_verify_usr_memory(usr_q,Qk,a,alpha,score,0.0).topkScore;
    #else
            double Rk_usr = topkSDijkstra_verify_usr_disk(usr_q,Qk,a,alpha,score,0.0).topkScore;
    #endif
#else
#ifndef DiskAccess
            double usr_lcl_rt = 0.0;
            double Rk_usr = TkGSKQ_bottom2up_verify_memory(u.id,Qk,a,alpha,score,usr_lcl_rt).top().score;
#else
            double usr_lcl_rt = 0.0;
            double Rk_usr = TkGSKQ_bottom2up_verify_disk(u.id,Qk,a,alpha,score,usr_lcl_rt).top().score;
#endif
#endif
            Users[i].topkScore_current = max(Users[i].topkScore_current,Rk_usr);

            if (score > Rk_usr) {
                size++;
                ResultDetail sd(i,influence,distance,relevance,score,Rk_usr);
                results.push_back(sd);
                //cout<<"u"<<usr.id<<",";

            }
#ifdef            NAVIE_DEBUG
            if(i==106&& poi_id==57){
                cout<<"find u"<<i<<",o"<<poi_id<<"score"<<score<<",Rk_u="<<Rk_usr<<endl;
                getchar();
            }
#endif
        }
        MRkGSKQ_results.push_back(results);
        cout << "RkGSKQ: size=" << size << endl;
        printResults(results);

    }

    //cout<<endl<<"for loc:"<<loc<<",";
    //cout<<"texual:";printQueryKeywords(keywords);
    cout<<"------------------------Naive_MRkGSKQ Complete!....---------------------------"<<endl;

    for(int i=0; i<stores.size();i++){
        cout <<"store"<<stores[i]<< "'s RkGSKQ: size=" << MRkGSKQ_results[i].size() << endl;
        printResults(MRkGSKQ_results[i]);
    }



    TIME_TICK_END
    TIME_TICK_PRINT("Naive_MRkGSKQ running time:")
    return 0;

}



int topK_count = 0;


void updateTopkScore_Current(int a, double alpha, int Qk, TopkResults& current_elements, InvList& term_candidateUserInv){
    for(Result r: current_elements){
        int o_id = r.id;
        for(int term: POIs[o_id].keywords){ //对topk结果o的每个关键词
            if(term_candidateUserInv[term].size()>0){
                for(int usr_id:term_candidateUserInv[term]){ //对每个关键词相关的用户
                    double rel; double inf;
                    if(Users[usr_id].verified==true) continue; //已经验证过了
                    double score_upper = getGSKScore_o2u_Upper(a, alpha, POIs[o_id], Users[usr_id], rel, inf);
                    if(Users[usr_id].current_Results.size()<Qk){  //插入内容
                        double score = getGSKScore_o2u(a,  alpha, POIs[o_id], Users[usr_id], rel, inf);
                        Users[usr_id].current_Results.push(Top_result(o_id, score));
                    }
                    else{  //更新内容
                        if(score_upper > Users[usr_id].current_Results.top().score){
                            double score = getGSKScore_o2u(a,  alpha, POIs[o_id], Users[usr_id], rel, inf);
                            if(score > Users[usr_id].current_Results.top().score){
                                Users[usr_id].current_Results.push(Top_result(o_id, score));
                                Users[usr_id].current_Results.pop();
                                //不断更新候选用户的当前topk score
                                Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, Users[usr_id].current_Results.top().score);
                                //cout<<"有更新！"<<Users[usr_id].topkScore_current<<endl;
                            }
                        }
                    }
                }
            }
        }
    }
}





/*------------------------- The algorithms below are for Batch RkGSKQ-------------------------*/


//1. 过滤阶段
//filter-base

void  Group_Filter_Baseline_memory(map<int, map<int,bool>>& usr_store_map, int store_partition_id, vector<int> stores, vector<int> keywords, int Qk, float a, float alpha, set<int>& verification_User, map<int,set<VerifyEntry>>& verification_Map, map<int, priority_queue<ResultCurrent>>& candidate_usr_related_store, map<int, LCLResult>& candidate_usr_LCL) {
    vector<int> qKey = keywords;
    set<int> store_leafSet;
    set<int> store_nodes;
    map<int, vector<int>> leaf_store_list;//<leaf_id, <store_id....>>  leaf 中的store的记录
    map<int, set<int>> leaf_verify_list;//<leaf_id, <usr_id....>>  leaf 中所有store对应的凝练用户集合


    BatchResults result_Stores; //<store_id, <usr_id....>>
    int loc =0;

    //generate the query tree for candidate stores
    for(int store_id:stores){
        int leaf = Nodes[POIs[store_id].Nj].gtreepath.back();
        leaf_store_list[leaf].push_back(store_id);
        store_leafSet.insert(leaf);
        store_nodes.insert(leaf);
        int current = leaf;
        GTree[current].stores.insert(store_id);
        int father;
        int i = 0;
        while (current != 0) {

            father = GTree[current].father;
            GTree[father].child_hasStores.insert(current);
            GTree[father].stores.insert(store_id);
            store_nodes.insert(father);
            current = father;
            i++;
        }
        //cout<<"层数"<<i<<",leaf="<<leaf<<endl;
    }
    cout<<"候选商店分布在"<<store_leafSet.size()<<"个不同叶节点中,GTree的"<<store_nodes.size()<<"个不同节点中"<<endl;

    set<int> keyword_Set;
    for(int t: keywords)
        keyword_Set.insert(t);

    //整合partition中的用户签到情况
    set<int> check_in_partition;
    for(int s: stores){
        for(int u: POIs[s].check_ins){
            check_in_partition.insert(u);
        }
    }

    //统计所有关键词相关用户个数
    set<int> totalUsrSet;
    for(int term:keyword_Set){
        for(int usr: invListOfUser[term]){
            totalUsrSet.insert(usr);
        }

    }


    //提取每个查询关键词 在leafNodeInv表中对应信息

    int prune=0; int unprune=0; int user_access=0; int usr_prune =0;
    int leaf_prune=0; int leaf_unprune=0;
    int may_prune = 0;
    //int topK_count = 0; //topk执行次数
    // 获取与查询关键词相关的用户叶子节点
    int t3=0; int root=0;



    //cout << " ,relate poi leaf size=" << related_poi_leafSet.size() << endl;
    //cout << " ,relate usr leaf size=" << related_usr_leafSet.size() << endl;

    vector<int> related_usr_node;   //相关用户节点列表

    map<int,vector<int>> related_poi_nodes;  //相关兴趣点节点列表

    // 获取与查询关键词相关的****用户****所在子图对应的上层节点, 并与query partition 关联

    vector<int> unodes = GTree[root].children;
    for(int unode: unodes){
        related_usr_node.push_back(unode);
        GTree[unode].ossociate_queryNodes = store_leafSet; //先只是查询对象组叶节点
    }


    // 获取与查询关键词相关的****兴趣点***** 所在子图对应的上层节点
    for(int keyword: keyword_Set){
        for(int child: GTree[root].inverted_list_o[keyword]){
            related_poi_nodes[keyword].push_back(child);
        }
    }




    set<int> related_children;
    set<int> related_usr;
    set<int> related_leaves;
    //verify the user node from upper to bottom
    int level = 1;  int cnt_ = 0;

    while (true) {  //daxunhuanonce
        //related_children.clear();
        if (related_usr_node.capacity() == 0) break;
        set<int>().swap(related_children);

        for (int node:related_usr_node) {
            //cout << "node" << node << "出列！" << endl;

            //cout<<"iteration"<<i<<":"<<endl;

            set<int> inter_Key = obtain_itersection_jins(GTree[node].userUKeySet, keyword_Set);

            //该节点为用户上层节点
            if(!GTree[node].isleaf) {
                int usrNum = 0;
                int store_leaf = store_partition_id;
                if (false) {   //用户稀疏节点  //x
                    //cout << "用户稀疏节点" << node <<",usrNum="<<usrNum<<endl;

                    for(int term: keyword_Set){
                        for(int usr_leaf: GTree[node].term_usrLeaf_Map[term]){
                            related_children.insert(usr_leaf);
                            for(int store_leaf:GTree[node].ossociate_queryNodes){
                                GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                            }
                        }

                    }

                    prune++;

                    //cout<<"用户稀疏节点 usr_node"<<node<<",prune"<<endl;   //剪枝整棵子树

                } else {  //用户密集型节点，需进行剪枝判断
                    //计算 评分上界
                    double max_Score = -1;

                    //max_Score = getGSKScoreP2U_Upper_GivenKey(a, alpha, keywords, store_leaf, node);
                    max_Score = getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, check_in_partition, store_leaf, node); //注意：这里用 InterKey, 与联合check-in计算的评分上界

                    double Rt_U = 100;
                    //获取该用户节点LCL

                    for(int term: inter_Key){
                        LCLResult lcl;
                        map<int,LCLResult>::iterator iter;
                        iter = GTree[node].termLCLDetail.find(term);
                        if(iter != GTree[node].termLCLDetail.end()){ //先访问cache
                            lcl = iter->second;
                        }else{
                            lcl = getUpperNode_termLCL(node,term, related_poi_nodes,Qk,a,alpha,max_Score);
                            GTree[node].termLCLDetail[term] = lcl;
                        }
                        double current_score = lcl.topscore;
                        Rt_U = min(Rt_U, current_score);
                    }


                    //cout<<"GTree[node].related_queryNodes[term_id].size="<<GTree[node].related_queryNodes[term_id].size()<<endl;
                    //对每个相关联的查询节点

                    if (max_Score > Rt_U) {  // 该子树下所有用户不能被剪枝，加入列表后续进一步验证
                        unprune++;

                        for(int term: inter_Key){
                            if(GTree[node].inverted_list_u.count(term)>0){
                                for (int child: GTree[node].inverted_list_u[term]) {
                                    //for (int child: textual_child){
                                    //孩子节点继承父节点的lcl内容
                                    GTree[child].termParentLCLDetail[term] = GTree[node].termParentLCLDetail[term];
                                    //GTree[child].termParentFilter[term_id] = filter;
                                    related_children.insert(child);
                                    //cout<<"加入n"<<child<<endl;

                                }
                            }

                        }

                        //getchar();
                    }
                        //整棵用户子树被剪枝
                    else {
                        prune++;
#ifdef PRUNELOG
                        cout << "节点" << node << "与store_leaf" << store_leaf << "无关， 其下所有用户被剪枝" << endl;
#endif
                    }

                    //cout<<"Rt_U="<<Rt_U<<endl;
                    //cout<<"get the LCL of u leaf"<<leaf_id<<endl;


                }//else用户密集型
            }

            else if (GTree[node].isleaf) {  //若该节点为用户叶子节点   onceleafsuoz
#ifdef TRACK
                cout << "access usr_leaf" << node << endl;
#endif

                //获取每个store_leaf
                int store_leaf = store_partition_id;
                double max_Score =
                        getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, check_in_partition, store_leaf, node);



                MultiLCLResult lclResult_multi;
                vector<int> Ukeys;  //重要！
                for(int term_id: inter_Key)
                    Ukeys.push_back(term_id);

                map<set<int>,MultiLCLResult>::iterator iter;
                iter = GTree[node].cacheMultiLCLDetail.find(inter_Key);
                if (iter!=GTree[node].cacheMultiLCLDetail.end()){
                    lclResult_multi = iter->second;
                }

                    //无cache, 算之
                else {

                    //lclResult_multi = getLeafNode_LCL_Dijkstra_multi_disk(node, Ukeys, Qk, a, alpha, max_Score);
#ifdef LV
                    lclResult_multi = getLeafNode_LCL_Dijkstra_multi(node, Ukeys, Qk, a, alpha, max_Score);
#endif
#ifndef LV
                    lclResult_multi = getLeafNode_LCL_bottom2up_multi(node, Ukeys, Qk, a, alpha, max_Score);
#endif


                    GTree[node].cacheMultiLCLDetail[inter_Key] = lclResult_multi;
                }

                //获得 usr_lcl update_o_set
                set<int> update_o_set;
                for(LCLResult lcl: lclResult_multi){
                    priority_queue<LCLEntry> LCL = lcl.LCL;
                    while(LCL.size() > 0){
                        LCLEntry le = LCL.top();
                        LCL.pop();
                        int o = le.id;
                        update_o_set.insert(o);
                        //cout<<"o"<<o<<",";
                    }

                }

                double Rt_U = -1.0;
                if(update_o_set.size()<Qk){
                    //cout<<"non-enough(empty) update o set"<<endl;
                    Rt_U = 0.0;
                }
                else{
                    for(LCLResult lclResult : lclResult_multi){
                        double _ss = lclResult.topscore;
                        if(-1.0==Rt_U)
                            Rt_U =_ss;
                        else
                            Rt_U = min(Rt_U,_ss);
                    }
                }
                //cout<<"Rt_U="<<Rt_U<<endl;
                //对相关查询对象节点进行剪枝评估
                set<int> L;
                //用户叶子节点剪枝条件判断
                if (max_Score > Rt_U) {  // U2P的上界大于 U.LCL下界
                    unprune++;
                    leaf_unprune++;
                    L.insert(store_leaf);
                }
                    //用户叶子节点整体被剪枝
                else {
                    //cout<<"leaf "<<store_leaf<<"与n"<<node<<"无关！"<<endl;
                    prune++;
                    leaf_prune++;
                }
                //对该叶子节点下各个用户进行评估
                set<int> usr_withinLeaf;
                for (int term: Ukeys) {
                    for (int usr_id: GTree[node].inverted_list_u[term]) {   //获取叶节点下每个用户
                        usr_withinLeaf.insert(usr_id);
                    }
                }

                for (int usr_id: usr_withinLeaf) {   //对叶节点下每个用户进行评估
                    //double u_lcl = 0;
                    double u_lcl_rt =0;
                    priority_queue<ResultCurrent> u_LCL = get_usrLCL_update_Plus(usr_id,
                                                                                 update_o_set, Qk, a, alpha, 100);  //获取该用户的lcl
                    priority_queue<ResultCurrent> u_related_store;

                    if (u_LCL.size() < Qk) u_lcl_rt = 0.0;
                    else if (u_LCL.size() == Qk) u_lcl_rt = u_LCL.top().score;


                    //cout<<"u_lcl="<<u_lcl<<endl;
                    //记录
                    Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);
                    for (int store: stores) {  //获取每个store_leaf中的store
                        double rel = 0;
                        double inf = 0;
                        Score_Bound gsk_score_bound;
#ifdef  GIVENKEY
                        vector<int> givenKeys = keywords;
                       gsk_score_bound =
                               getGSKScore_o2u_Upper_and_Lower_givenKey(a, alpha, POIs[store], Users[usr_id], givenKeys,rel, inf);
#else
                        gsk_score_bound =
                                getGSKScore_o2u_Upper_and_Lower(a, alpha, POIs[store], Users[usr_id],rel, inf);
#endif
                        double gsk_score_Upper = gsk_score_bound.upper_score;
                        double gsk_score_Lower = gsk_score_bound.lower_score;
                        //TopkQueryCurrentResult result = topkSDijkstra_verify_usr(transformToQ(Users[usr_id]),Qk,a, alpha,gsk_score_Upper);
                        //u_lcl = max(u_lcl,result.topkScore);
                        //对用户与查询对象进行判断
                        if (u_lcl_rt < gsk_score_Upper) {
                            //update lcl
                            u_related_store.push(ResultCurrent(store, gsk_score_Lower, gsk_score_Upper, rel, inf));


                            //cout<<"u"<<usr_id<<", gsk_score_Lower="<<gsk_score_Lower<<",u_lcl="<<u_lcl_rt<<",加入store="<<store<<endl;
                        } else {
                            usr_prune++;
                            //continue;   //该用户不是反top-k结果
                        }
                    }

                    if (u_related_store.size() > 0) {
                        verification_User.insert(usr_id);

                        if (candidate_usr_related_store[usr_id].size() > 0) {  //之前已经有u_lcl的结果了，更新之

                            //将拥有更有评分的object与user关联
                            while (!u_related_store.empty()) {
                                ResultCurrent rc = u_related_store.top();
                                u_related_store.pop();

                                if(usr_store_map[usr_id][rc.o_id] ==true){

                                }else{
                                    candidate_usr_related_store[usr_id].push(rc);

                                    usr_store_map[usr_id][rc.o_id] = true;
                                }


                                while(candidate_usr_related_store[usr_id].size()>Qk){
                                    candidate_usr_related_store[usr_id].pop();
                                }

                            }

                        } else {  //若无，直接赋予之
                            //priority_queue<ResultCurrent> u_related_Query = u_related_store;
                            candidate_usr_related_store[usr_id] = u_related_store;
                        }


                    } else {
                        //cout << "u" << usr_id << "用户被减, as Lu为空,u_lcl=" << u_lcl << endl;
                    }
                    //getchar();
                }

                //cout<<"验证叶节点结束"<<node<<endl;
            }



        }// the first for loop
        level++;
        map<int, vector<int>> related_poi_nodes_next;
        map<int, vector<int>> related_poi_leaf_next;

        //upfold related poi nodes
        for(int term: keyword_Set){
            for(int pn: related_poi_nodes[term]){
                if(!GTree[pn].isleaf){
                    if(GTree[pn].inverted_list_o[term].size()!=0){
                        for(int child: GTree[pn].inverted_list_o[term]){
                            related_poi_nodes_next[term].push_back(child);
                        }
                    }
                }
                else{
                    related_poi_leaf_next[term].push_back(pn);
                }
            }
        }



        // upfold the related usr nodes
        vector<int>().swap(related_usr_node);
        related_usr_node.clear();

        for (int cn:related_children){       //jiachildyici

            related_usr_node.push_back(cn);
            //cout<<"新的usr node 加入:n"<<cn<<"加入"<<endl;
        }
        if(related_usr_node.size() ==0){
            cout<<"related_usr_node空了"<<endl;
            break;

        }
        //防止单路分支带来的重复验证
        if(related_children.size()==1){
            cnt_++;
            if(cnt_==3){
                //cout<<"直接跳到叶节点"<<endl;
                //getchar();
                related_usr_node.clear();
                for(int node: related_children){
                    if(GTree[node].isleaf)
                        //related_usr_node.push_back(node);
                        related_usr_node.push_back(node);
                    else{
                        set<int> leaves;
                        for(int term: keyword_Set){
                            if(!(GTree[node].term_usrLeaf_Map.count(term)>0))
                                continue;
                            for(int leaf: GTree[node].term_usrLeaf_Map[term]){
                                leaves.insert(leaf);
                            }
                        }
                        for(int usr_leaf: leaves){
                            related_usr_node.push_back(usr_leaf);
                        }

                    }
                }
                //cout<<"新加入"<<related_usr_node.size()<<endl;
                //getchar();

            }

        }
        //cout<<"第"<<level<<"层用户节点剪枝验证结束！"<<endl;
        //cout<<"related_poi_node size in next iteration:"<<related_poi_nodes.size()<<endl;
        //cout<<"related_usr_node size in next iteration:"<<related_usr_node.size()<<endl;
        if(related_usr_node.size()==0) break;

    }//end while
    //

}

void  Group_Filter_Baseline_disk(map<int, map<int,bool>>& usr_store_map, int store_partition_id, vector<int> stores, vector<int> keywords, int Qk, float a, float alpha, set<int>& verification_User, map<int,set<VerifyEntry>>& verification_Map, map<int, priority_queue<ResultCurrent>>& candidate_usr_related_store, map<int, LCLResult>& candidate_usr_LCL) {
    vector<int> qKey = keywords;
    set<int> store_leafSet;
    set<int> store_nodes;
    map<int, vector<int>> leaf_store_list;//<leaf_id, <store_id....>>  leaf 中的store的记录
    map<int, set<int>> leaf_verify_list;//<leaf_id, <usr_id....>>  leaf 中所有store对应的凝练用户集合


    BatchResults result_Stores; //<store_id, <usr_id....>>
    int loc =0;

    //generate the query tree for candidate stores
    for(int store_id:stores){
        int leaf = Nodes[POIs[store_id].Nj].gtreepath.back();
        leaf_store_list[leaf].push_back(store_id);
        store_leafSet.insert(leaf);
        store_nodes.insert(leaf);
        int current = leaf;
        GTree[current].stores.insert(store_id);
        int father;
        int i = 0;
        while (current != 0) {
#ifdef DiskAccess
            TreeNode tn = getGIMTreeNodeData(current,OnlyO);

#endif
            father = GTree[current].father;
            GTree[father].child_hasStores.insert(current);
            GTree[father].stores.insert(store_id);
            store_nodes.insert(father);
            current = father;
            i++;
        }
        //cout<<"层数"<<i<<",leaf="<<leaf<<endl;
    }
    cout<<"候选商店分布在"<<store_leafSet.size()<<"个不同叶节点中,GTree的"<<store_nodes.size()<<"个不同节点中"<<endl;

    set<int> keyword_Set;
    for(int t: keywords)
        keyword_Set.insert(t);

    //整合partition中的用户签到情况
    set<int> check_in_partition;
    for(int s: stores){
        for(int u: POIs[s].check_ins){
            check_in_partition.insert(u);
        }
    }

    //统计所有关键词相关用户个数
    set<int> totalUsrSet;
    for(int term:keyword_Set){
        for(int usr: invListOfUser[term]){
            totalUsrSet.insert(usr);
        }

    }


    //提取每个查询关键词 在leafNodeInv表中对应信息

    int prune=0; int unprune=0; int user_access=0; int usr_prune =0;
    int leaf_prune=0; int leaf_unprune=0;
    int may_prune = 0;
    //int topK_count = 0; //topk执行次数
    // 获取与查询关键词相关的用户叶子节点
    int t3=0; int root=0;



    //cout << " ,relate poi leaf size=" << related_poi_leafSet.size() << endl;
    //cout << " ,relate usr leaf size=" << related_usr_leafSet.size() << endl;

    vector<int> related_usr_node;   //相关用户节点列表

    map<int,vector<int>> related_poi_nodes;  //相关兴趣点节点列表

    // 获取与查询关键词相关的****用户****所在子图对应的上层节点, 并与query partition 关联

    vector<int> unodes = GTree[root].children;
    for(int unode: unodes){
        related_usr_node.push_back(unode);
        GTree[unode].ossociate_queryNodes = store_leafSet; //先只是查询对象组叶节点
    }


    // 获取与查询关键词相关的****兴趣点***** 所在子图对应的上层节点
    for(int keyword: keyword_Set){
        for(int child: getObjectTermRelatedEntry(keyword,root)){
            related_poi_nodes[keyword].push_back(child);
        }
    }




    set<int> related_children;
    set<int> related_usr;
    set<int> related_leaves;
    //verify the user node from upper to bottom
    int level = 1;  int cnt_ = 0;

    while (true) {  //daxunhuanonce
        //related_children.clear();
        if (related_usr_node.capacity() == 0) break;
        set<int>().swap(related_children);

        for (int node:related_usr_node) {
            //cout << "node" << node << "出列！" << endl;

            //cout<<"iteration"<<i<<":"<<endl;
            TreeNode un = getGIMTreeNodeData(node,OnlyU);

            set<int> inter_Key = obtain_itersection_jins(un.userUKeySet, keyword_Set);

            //该节点为用户上层节点
            if(!un.isleaf) {
                int usrNum = 0;
                int store_leaf = store_partition_id;
                if (false) {   //用户稀疏节点  //x
                    //cout << "用户稀疏节点" << node <<",usrNum="<<usrNum<<endl;

                    for(int term: keyword_Set){
                        for(int usr_leaf: GTree[node].term_usrLeaf_Map[term]){
                            related_children.insert(usr_leaf);
                            for(int store_leaf:GTree[node].ossociate_queryNodes){
                                GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                            }
                        }

                    }

                    prune++;

                    //cout<<"用户稀疏节点 usr_node"<<node<<",prune"<<endl;   //剪枝整棵子树

                } else {  //用户密集型节点，需进行剪枝判断
                    //计算 评分上界
                    double max_Score = -1;

                    //max_Score = getGSKScoreP2U_Upper_GivenKey(a, alpha, keywords, store_leaf, node);
                    max_Score = getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, check_in_partition, store_leaf, node); //注意：这里用 InterKey, 与联合check-in计算的评分上界

                    double Rt_U = 100;
                    //获取该用户节点LCL

                    for(int term: inter_Key){
                        LCLResult lcl;
                        map<int,LCLResult>::iterator iter;
                        iter = GTree[node].termLCLDetail.find(term);
                        if(iter != GTree[node].termLCLDetail.end()){ //先访问cache
                            lcl = iter->second;
                        }else{
                            lcl = getUpperNode_termLCL_Disk(node,term, related_poi_nodes,Qk,a,alpha,max_Score);
                            GTree[node].termLCLDetail[term] = lcl;
                        }
                        double current_score = lcl.topscore;
                        Rt_U = min(Rt_U, current_score);
                    }


                    //cout<<"GTree[node].related_queryNodes[term_id].size="<<GTree[node].related_queryNodes[term_id].size()<<endl;
                    //对每个相关联的查询节点

                    if (max_Score > Rt_U) {  // 该子树下所有用户不能被剪枝，加入列表后续进一步验证
                        unprune++;

                        for(int term: inter_Key){
                            vector<int> usr_leafList = getUsrTermRelatedEntry(term,node);
                            if(usr_leafList.size()>0){
                                for (int child: usr_leafList) {
                                    //for (int child: textual_child){
                                    //孩子节点继承父节点的lcl内容
                                    GTree[child].termParentLCLDetail[term] = GTree[node].termParentLCLDetail[term];
                                    //GTree[child].termParentFilter[term_id] = filter;
                                    related_children.insert(child);
                                    //cout<<"加入n"<<child<<endl;

                                }
                            }

                        }

                        //getchar();
                    }
                        //整棵用户子树被剪枝
                    else {
                        prune++;
#ifdef PRUNELOG
                        cout << "节点" << node << "与store_leaf" << store_leaf << "无关， 其下所有用户被剪枝" << endl;
#endif
                    }

                    //cout<<"Rt_U="<<Rt_U<<endl;
                    //cout<<"get the LCL of u leaf"<<leaf_id<<endl;


                }//else用户密集型
            }

            else if (un.isleaf) {  //若该节点为用户叶子节点   onceleafsuoz
#ifdef TRACK
                cout << "access usr_leaf" << node << endl;
#endif

                //获取每个store_leaf
                int store_leaf = store_partition_id;
                double max_Score =
                        getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, check_in_partition, store_leaf, node);



                MultiLCLResult lclResult_multi;
                vector<int> Ukeys;  //重要！
                for(int term_id: inter_Key)
                    Ukeys.push_back(term_id);

                map<set<int>,MultiLCLResult>::iterator iter;
                iter = GTree[node].cacheMultiLCLDetail.find(inter_Key);
                if (iter!=GTree[node].cacheMultiLCLDetail.end()){
                    lclResult_multi = iter->second;
                }

                    //无cache, 算之
                else {

                    //lclResult_multi = getLeafNode_LCL_Dijkstra_multi_disk(node, Ukeys, Qk, a, alpha, max_Score);
#ifdef LV
                    lclResult_multi = getLeafNode_LCL_Dijkstra_multi(node, Ukeys, Qk, a, alpha, max_Score);
#endif
#ifndef LV
                    lclResult_multi = getLeafNode_LCL_bottom2up_multi(node, Ukeys, Qk, a, alpha, max_Score);
#endif


                    GTree[node].cacheMultiLCLDetail[inter_Key] = lclResult_multi;
                }

                //获得 usr_lcl update_o_set
                set<int> update_o_set;
                for(LCLResult lcl: lclResult_multi){
                    priority_queue<LCLEntry> LCL = lcl.LCL;
                    while(LCL.size() > 0){
                        LCLEntry le = LCL.top();
                        LCL.pop();
                        int o = le.id;
                        update_o_set.insert(o);
                        //cout<<"o"<<o<<",";
                    }

                }

                double Rt_U = -1.0;
                if(update_o_set.size()<Qk){
                    //cout<<"non-enough(empty) update o set"<<endl;
                    Rt_U = 0.0;
                }
                else{
                    for(LCLResult lclResult : lclResult_multi){
                        double _ss = lclResult.topscore;
                        if(-1.0==Rt_U)
                            Rt_U =_ss;
                        else
                            Rt_U = min(Rt_U,_ss);
                    }
                }
                //cout<<"Rt_U="<<Rt_U<<endl;
                //对相关查询对象节点进行剪枝评估
                set<int> L;
                //用户叶子节点剪枝条件判断
                if (max_Score > Rt_U) {  // U2P的上界大于 U.LCL下界
                    unprune++;
                    leaf_unprune++;
                    L.insert(store_leaf);
                }
                    //用户叶子节点整体被剪枝
                else {
                    //cout<<"leaf "<<store_leaf<<"与n"<<node<<"无关！"<<endl;
                    prune++;
                    leaf_prune++;
                }
                //对该叶子节点下各个用户进行评估
                set<int> usr_withinLeaf;
                for (int term: Ukeys) {
                    vector<int> usr_list = getUsrTermRelatedEntry(term,node);
                    for (int usr_id: usr_list) {   //获取叶节点下每个用户
                        usr_withinLeaf.insert(usr_id);
                    }
                }

                for (int usr_id: usr_withinLeaf) {   //对叶节点下每个用户进行评估
                    //double u_lcl = 0;
                    double u_lcl_rt =0;
                    priority_queue<ResultCurrent> u_LCL = get_usrLCL_update_Plus(usr_id,
                                                                                 update_o_set, Qk, a, alpha, 100);  //获取该用户的lcl
                    priority_queue<ResultCurrent> u_related_store;

                    if (u_LCL.size() < Qk) u_lcl_rt = 0.0;
                    else if (u_LCL.size() == Qk) u_lcl_rt = u_LCL.top().score;


                    //cout<<"u_lcl="<<u_lcl<<endl;
                    //记录
                    Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);
                    for (int store: stores) {  //获取每个store_leaf中的store
                        double rel = 0;
                        double inf = 0;
                        Score_Bound gsk_score_bound;
                        POI p = getPOIFromO2UOrgLeafData(store);
                        User u = getUserFromO2UOrgLeafData(usr_id);

#ifdef  GIVENKEY
                        vector<int> givenKeys = keywords;
                        gsk_score_bound =
                                getGSKScore_o2u_Upper_and_Lower_givenKey(a, alpha, p, u, givenKeys,rel, inf);
#else
                        gsk_score_bound =
                                getGSKScore_o2u_Upper_and_Lower(a, alpha, p, u,rel, inf);
#endif
                        double gsk_score_Upper = gsk_score_bound.upper_score;
                        double gsk_score_Lower = gsk_score_bound.lower_score;
                        //TopkQueryCurrentResult result = topkSDijkstra_verify_usr(transformToQ(Users[usr_id]),Qk,a, alpha,gsk_score_Upper);
                        //u_lcl = max(u_lcl,result.topkScore);
                        //对用户与查询对象进行判断
                        if (u_lcl_rt < gsk_score_Upper) {
                            //update lcl
                            u_related_store.push(ResultCurrent(store, gsk_score_Lower, gsk_score_Upper, rel, inf));


                            //cout<<"u"<<usr_id<<", gsk_score_Lower="<<gsk_score_Lower<<",u_lcl="<<u_lcl_rt<<",加入store="<<store<<endl;
                        } else {
                            usr_prune++;
                            //continue;   //该用户不是反top-k结果
                        }
                    }

                    if (u_related_store.size() > 0) {
                        verification_User.insert(usr_id);

                        if (candidate_usr_related_store[usr_id].size() > 0) {  //之前已经有u_lcl的结果了，更新之

                            //将拥有更有评分的object与user关联
                            while (!u_related_store.empty()) {
                                ResultCurrent rc = u_related_store.top();
                                u_related_store.pop();

                                if(usr_store_map[usr_id][rc.o_id] ==true){

                                }else{
                                    candidate_usr_related_store[usr_id].push(rc);

                                    usr_store_map[usr_id][rc.o_id] = true;
                                }


                                while(candidate_usr_related_store[usr_id].size()>Qk){
                                    candidate_usr_related_store[usr_id].pop();
                                }

                            }

                        } else {  //若无，直接赋予之
                            //priority_queue<ResultCurrent> u_related_Query = u_related_store;
                            candidate_usr_related_store[usr_id] = u_related_store;
                        }


                    } else {
                        //cout << "u" << usr_id << "用户被减, as Lu为空,u_lcl=" << u_lcl << endl;
                    }
                    //getchar();
                }

                //cout<<"验证叶节点结束"<<node<<endl;
            }



        }// the first for loop
        level++;
        map<int, vector<int>> related_poi_nodes_next;
        map<int, vector<int>> related_poi_leaf_next;

        //upfold related poi nodes
        for(int term: keyword_Set){
            for(int pn: related_poi_nodes[term]){
                if(!GTree[pn].isleaf){
                    vector<int> oentry_list = getObjectTermRelatedEntry(term,pn);
                    if(oentry_list.size()!=0){
                        for(int child: oentry_list){
                            related_poi_nodes_next[term].push_back(child);
                        }
                    }
                }
                else{
                    related_poi_leaf_next[term].push_back(pn);
                }
            }
        }



        // upfold the related usr nodes
        vector<int>().swap(related_usr_node);
        related_usr_node.clear();

        for (int cn:related_children){       //jiachildyici

            related_usr_node.push_back(cn);
            //cout<<"新的usr node 加入:n"<<cn<<"加入"<<endl;
        }
        if(related_usr_node.size() ==0){
            cout<<"related_usr_node空了"<<endl;
            break;

        }
        //防止单路分支带来的重复验证
        if(related_children.size()==1){
            cnt_++;
            if(cnt_==3){
                //cout<<"直接跳到叶节点"<<endl;
                //getchar();
                related_usr_node.clear();
                for(int node: related_children){
                    if(GTree[node].isleaf)
                        //related_usr_node.push_back(node);
                        related_usr_node.push_back(node);
                    else{
                        set<int> leaves;  //必须是set, 否则会重复
                        for(int term: keyword_Set){

                            vector<int> uLeaf_list;
                            uLeaf_list = getUsrTermRelatedLeafNode(term,node);
                            for(int leaf: uLeaf_list){
                                leaves.insert(leaf);
                            }
                        }
                        for(int usr_leaf: leaves){
                            related_usr_node.push_back(usr_leaf);
                        }

                    }
                }
                //cout<<"新加入"<<related_usr_node.size()<<endl;
                //getchar();

            }

        }
        //cout<<"第"<<level<<"层用户节点剪枝验证结束！"<<endl;
        //cout<<"related_poi_node size in next iteration:"<<related_poi_nodes.size()<<endl;
        //cout<<"related_usr_node size in next iteration:"<<related_usr_node.size()<<endl;
        if(related_usr_node.size()==0) break;

    }//end while
    //

}


//2. 验证阶段
//内存操作
BRkGSKQResults Verification_BRkGSKQ_memory(set<int> verification_User, map<int, priority_queue<ResultCurrent>> candidate_usr_related_store, vector<int> qKey, int Qk, float a, float alpha){
    //存放RkGSKQ结果的
    BRkGSKQResults result_Stores;
    set<int> candidate_usr = verification_User;//iter->second; //
    //cout<<"candidate_usr.size="<<candidate_usr.size()<<endl;
    map<int, vector<int>> usr_group;
    //1. 按maxscore，将用户排序(OK)；2. 按关键词将候选用户分组
    vector<CardinalityFisrt> user_list;
    map<int, vector<int>> term_candidateUserInv;
    for(int usr_id : candidate_usr){
        int leaf_node = Nodes[Users[usr_id].Nj].gtreepath.back();
        usr_group[leaf_node].push_back(usr_id);
        priority_queue<ResultCurrent> relate_query_queue = candidate_usr_related_store[usr_id];
        int size= relate_query_queue.size();
        //cout<<"leaf_store_list[store_leaf].size="<<candidate_usr_related_store_leaf[usr_id].size()<<endl;
        double maximum_upper_score = 0; double maximum_lower_score = 0; double max_score = 0;
        priority_queue<ResultLargerFisrt> Lu;

        while(!relate_query_queue.empty()){
            ResultCurrent rc = relate_query_queue.top();
            relate_query_queue.pop();
            //先算出真实gsk评分
            double rel = rc.relevance; double inf = rc.influence;
            double gsk_score = getGSKScore_q2u(a,  alpha, POIs[rc.o_id], Users[usr_id], rel, inf);
            ResultLargerFisrt rlf(rc.o_id,gsk_score,rel, inf);
            Lu.push(rlf);

            //max_score = max(max_score,gsk_score);
        }
        CardinalityFisrt cf(usr_id,0.0,Lu);
        //user_queue.push(cf);
        user_list.push_back(cf);
        //构建候选用户的关键词倒排表
        /*for(int term: Users[usr_id].keywords){
            term_candidateUserInv[term].push_back(usr_id);
        }*/

    }
    //cout<<"候选用户被分成"<<usr_group.size()<<"组"<<endl;

    // verifiying the TkGSKQ results for each users in candidate users, report query objects which havie user as top-k
    int useful =0; int usr_count =0;

    for(int i=0; i<user_list.size();i++){
        //topkSDijkstra_verify_usr(transformToQ(Users[0]), Qk, a, alpha, 1);
        CardinalityFisrt cf = user_list[i];
        usr_count++;
        int usr_id = cf.usr_id;
        priority_queue<ResultLargerFisrt> Lu = cf.Lu;
        //cout<<"第"<<usr_count<<"个验证用户,关联"<<Lu.size()<<"家store,max_score="<<Lu.top().score<<endl;
        if(Lu.top().score > Users[usr_id].topkScore_current){


            double score_bound = Lu.top().score; //Users[usr_id].topkScore_current;
            double Rk_u = candidate_user_verify_single(usr_id,Qk,a, alpha,score_bound);

            //更新其他用户的当前topk score值
            //updateTopkScore_Current(a, Qk, current_elements,term_candidateUserInv);
            topK_count++;
            if(Lu.top().score > Rk_u){
                useful++;
                //对Lu中的各个query object进行评估
                while(!Lu.empty()){
                    ResultLargerFisrt rlf = Lu.top();
                    Lu.pop();
                    double rel = 0; double inf = 0;
                    double gsk_score = rlf.score; int query_object = rlf.o_id;
                    if(gsk_score > Rk_u){
                        ResultDetail rd(usr_id, 0, 0, 0, gsk_score, Rk_u);
                        result_Stores[query_object].push_back(rd);
                        //cout<<"u"<<usr_id<<"确实在o"<<query_object<<"的RkGSKQ结果集合中，gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
                    }
                    else break;
                }
            }
            else{
                //cout<<"topk score = "<<Rk_u<<",max_score="<<max_score<<endl;
            }
        }
        else{
            //cout<<"不用算！"<<endl;
        }
        //该候选用户已被验证
        Users[usr_id].verified = true;
    }
    //cout<<"useful ratio"<<useful*1.0/topK_count<<endl;
    //getchar();
    return  result_Stores;
}

//外存操作(ok)
BRkGSKQResults Verification_BRkGSKQ_disk(set<int> verification_User, map<int, priority_queue<ResultCurrent>> candidate_usr_related_store, vector<int> qKey, int Qk, float a, float alpha){
    //存放RkGSKQ结果的
    BRkGSKQResults result_Stores;
    set<int> candidate_usr = verification_User;//iter->second; //
    //cout<<"candidate_usr.size="<<candidate_usr.size()<<endl;
    map<int, vector<int>> usr_group;
    //1. 按maxscore，将用户排序(OK)；2. 按关键词将候选用户分组
    vector<CardinalityFisrt> user_list;
    map<int, vector<int>> term_candidateUserInv;
    for(int usr_id : candidate_usr){
        User u = getUserFromO2UOrgLeafData(usr_id);
        int leaf_node = Nodes[u.Nj].gtreepath.back();
        usr_group[leaf_node].push_back(usr_id);
        priority_queue<ResultCurrent> relate_query_queue = candidate_usr_related_store[usr_id];
        int size= relate_query_queue.size();

        double maximum_upper_score = 0; double maximum_lower_score = 0; double max_score = 0;
        priority_queue<ResultLargerFisrt> Lu;

        while(!relate_query_queue.empty()){
            ResultCurrent rc = relate_query_queue.top();
            relate_query_queue.pop();
            //先算出真实gsk评分
            double rel = rc.relevance; double inf = rc.influence;
            POI p = getPOIFromO2UOrgLeafData(rc.o_id);
            User u = getUserFromO2UOrgLeafData(usr_id);
            double gsk_score = getGSKScore_q2u(a,  alpha, p, u, rel, inf);
            ResultLargerFisrt rlf(rc.o_id,gsk_score,rel, inf);
            Lu.push(rlf);

            //max_score = max(max_score,gsk_score);
        }
        CardinalityFisrt cf(usr_id,0.0,Lu);
        //user_queue.push(cf);
        user_list.push_back(cf);
        //构建候选用户的关键词倒排表
        /*for(int term: Users[usr_id].keywords){
            term_candidateUserInv[term].push_back(usr_id);
        }*/

    }
    //cout<<"候选用户被分成"<<usr_group.size()<<"组"<<endl;

    // verifiying the TkGSKQ results for each users in candidate users, report query objects which havie user as top-k
    int useful =0; int usr_count =0;

    for(int i=0; i<user_list.size();i++){
        //topkSDijkstra_verify_usr(transformToQ(Users[0]), Qk, a, alpha, 1);
        CardinalityFisrt cf = user_list[i];
        usr_count++;
        int usr_id = cf.usr_id;
        priority_queue<ResultLargerFisrt> Lu = cf.Lu;
        //cout<<"第"<<usr_count<<"个验证用户,关联"<<Lu.size()<<"家store,max_score="<<Lu.top().score<<endl;
        if(Lu.top().score > Users[usr_id].topkScore_current){


            double score_bound = Lu.top().score; //Users[usr_id].topkScore_current;
            double Rk_u = candidate_user_verify_single(usr_id,Qk,a, alpha,score_bound);
            //cout<<"u"<<usr_id<<" Rk_u="<<Rk_u<<endl;
            //更新其他用户的当前topk score值
            //updateTopkScore_Current(a, Qk, current_elements,term_candidateUserInv);
            topK_count++;
            if(Lu.top().score > Rk_u){
                useful++;
                //对Lu中的各个query object进行评估
                while(!Lu.empty()){
                    ResultLargerFisrt rlf = Lu.top();
                    Lu.pop();
                    double rel = 0; double inf = 0;
                    double gsk_score = rlf.score; int query_object = rlf.o_id;
                    if(gsk_score > Rk_u){
                        ResultDetail rd(usr_id, 0, 0, 0, gsk_score, Rk_u);
                        result_Stores[query_object].push_back(rd);
                        //cout<<"u"<<usr_id<<"确实在o"<<query_object<<"的RkGSKQ结果集合中，gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
                    }
                    else break;
                }
            }
            else{
                //cout<<"topk score = "<<Rk_u<<",max_score="<<max_score<<endl;
            }
        }
        else{
            //cout<<"不用算！"<<endl;
        }
        //该候选用户已被验证
        Users[usr_id].verified = true;
    }
    //cout<<"useful ratio"<<useful*1.0/topK_count<<endl;
    //getchar();
    return  result_Stores;
}


//put above steps togather
BatchResults RkGSKQ_Batch_Baseline(vector<int> stores, vector<int> keywords, int Qk, float a, float alpha){
    cout<<"MRkGSKQ Baseline begin!"<<endl;
    int block_num_pre = block_num;

    vector<int> qKey = keywords;
    set<int> store_leafSet;
    set<int> store_nodes;
    map<int, vector<int>> leaf_store_list;//<leaf_id, <store_id....>>  leaf 中的store的记录
    map<int, set<int>> leaf_verify_list;//<leaf_id, <usr_id....>>  leaf 中所有store对应的凝练用户集合

    int loc =0;
    int node_visite = 0;


    //generate the query tree for candidate stores
    for(int store_id:stores){
        POI p;
#ifndef DiskAccess
        p = POIs[store_id];
#else
        p = getPOIFromO2UOrgLeafData(store_id);
#endif
        int leaf = Nodes[p.Nj].gtreepath.back();
        leaf_store_list[leaf].push_back(store_id);

        store_leafSet.insert(leaf);
        store_nodes.insert(leaf);
        int current = leaf;
        GTree[current].stores.insert(store_id);
        int father;
        int i = 0;
        while (current != 0) {
            father = GTree[current].father;
            GTree[father].child_hasStores.insert(current);
            GTree[father].stores.insert(store_id);
            store_nodes.insert(father);
            current = father;
            i++;
        }
        //cout<<"层数"<<i<<",leaf="<<leaf<<endl;
    }
    cout<<"候选商店分布在"<<store_leafSet.size()<<"个不同叶节点中,GTree的"<<store_nodes.size()<<"个不同节点中"<<endl;
    //getchar();

    set<int> keyword_Set;
    for(int t: keywords)
        keyword_Set.insert(t);

    set<int> checkin_usr_Set;


    vector<int> Leaf_node_User,Leaf_node_Poi;

    vector<ResultDetail> result_User;
    //提取每个查询关键词 在leafNodeInv表中对应信息

    //对每个有usr的叶节点进行上下界的评估
    int prune=0; int unprune=0; int user_access=0; int usr_prune =0;
    int leaf_prune=0; int leaf_unprune=0;
    int may_prune = 0;
    int topK_count = 0; //topk执行次数
    // 获取与查询关键词相关的用户叶子节点
    set<int> related_usrLeafSet;
    int t3=0; int root=0;

    set<int> totalUsrSet;
    Filter_START
    int nodeAccess = 0;

    /*----------分别对每个query cluster执行Group Filtering------------------*/
    set<int> verification_User;  //候选验证用户集
    map<int, set<VerifyEntry>> verification_Map;    //凝练结果表   <leaf_id, <user_id....>>
    map<int, priority_queue<ResultCurrent>> candidate_usr_related_store;   //候选用户与查询对象的关联表
    //map<int, vector<ResultDetail>> result_Stores; //<store_id, <usr_id....>>
    map<int, LCLResult> candidate_usr_LCL;     //候选用户的LCL表
    cout<<"共"<<store_leafSet.size()<<"个query partition"<<endl; //getchar();
    TIME_TICK_START
    Filter_START
    map<int, map<int,bool>> usr_store_map;   //避免 user与object重复关联，重要！
    int p_count = 0;
    for(int store_leaf: store_leafSet){   //对每个leaf中的store进行处理
        //cout<< store_leaf<<endl;
        p_count++;
        vector<int> stores_inpartition =  leaf_store_list[store_leaf];   //store_leaf
#ifndef DiskAccess
        Group_Filter_Baseline_memory(usr_store_map, store_leaf, stores_inpartition, keywords, Qk, a, alpha, verification_User, verification_Map, candidate_usr_related_store, candidate_usr_LCL);
#else
        Group_Filter_Baseline_disk(usr_store_map, store_leaf, stores_inpartition, keywords, Qk, a, alpha, verification_User, verification_Map, candidate_usr_related_store, candidate_usr_LCL);
#endif
        //PAUSE_START
        cout<<"query partition "<<store_leaf<<"的filter处理完毕,剩余"<<(store_leafSet.size()-p_count)<<"个partition"<<endl;
        cout<<"当前找到"<<verification_User.size()<<"个候选用户"<<endl;

    }
    Filter_END
    //对所有candidate
    cout<<"所有 query partition 的filter 步骤处理完毕"<<endl;
    cout<<"累计总共找到"<<verification_User.size()<<"个候选用户"<<endl;

    Refine_START

    //对所有query_partition对应汇总的candidate进行验证处理
#ifndef DiskAccess
    BatchResults batch_query_results = Verification_BRkGSKQ_memory(verification_User, candidate_usr_related_store, keywords, Qk, a, alpha);
#else
    BatchResults batch_query_results = Verification_BRkGSKQ_disk(verification_User, candidate_usr_related_store, keywords, Qk, a, alpha);
#endif
    cout<<"有关联用户的query个数"<<batch_query_results.size()<<endl;  //baseyan

    Refine_END


    //对每个候选usr进行分析并打印结果
    //FilterResults results =  resultsAnalysis(result_User,keywords,checkin_usr,loc, Qk,a,alpha);
    //printBatchRkGSKQResults(stores, batch_query_results);
    cout<<"累计总共找到"<<verification_User.size()<<"个候选用户"<<endl;
    //cout<<"与兴趣点关键词相关用户总数"<< totalUsrSet.size() <<endl;
    cout<<"topk查询执行个数 =" << topK_count <<endl;
    cout<<"g-tree node prune ="<<prune<<" ,unprune="<<unprune<<",usr_prune="<<usr_prune<<endl;
    double ratio = (1.0*prune) /(prune+unprune);
    double ratio2 = (1.0*leaf_prune) /(leaf_prune+leaf_unprune);
    cout<<"节点剪枝率= "<<ratio<<"，其中叶子节点剪枝率="<<ratio2<<endl;
    //cout<<"g-tree node access ="<<nodeAccess<<endl;
    TIME_TICK_END
    Filter_PRINT("filter lasts")
    Refine_PRINT("refinement lasts")
    int io_count = block_num - block_num_pre;
    cout<<"i/o cost = "<<io_count<<endl;
    TIME_TICK_PRINT("MRkGSKQ Baseline runtime is ")


    return batch_query_results;

}

vector<int> transforSet2Vector(set<int> _set){
    vector<int> _list;
    for(int element: _set)
        _list.push_back(element);
    return _list;
}

BatchResults RkGSKQ_Batch_Baseline_formal(vector<int> stores, int Qk, float a, float alpha){
    cout<<"MRkGSKQ Baseline begin!"<<endl;
    int block_num_pre = block_num;

    vector<int> keywords;
    set<int> store_leafSet;
    set<int> store_nodes;
    map<int, vector<int>> leaf_store_list;//<leaf_id, <store_id....>>  leaf 中的store的记录
    map<int, set<int>> leaf_unionKeyword;//<leaf_id, <store_id....>>  leaf 中的store的记录
    map<int, set<int>> leaf_verify_list;//<leaf_id, <usr_id....>>  leaf 中所有store对应的凝练用户集合

    int loc =0;
    int node_visite = 0;


    //generate the query tree for candidate stores
    for(int store_id:stores){
        POI p;
#ifndef DiskAccess
        p = POIs[store_id];
#else
        p = getPOIFromO2UOrgLeafData(store_id);
#endif
        int leaf = Nodes[p.Nj].gtreepath.back();
        leaf_store_list[leaf].push_back(store_id);
        //将poi的keyword加入所在叶节点的unionKeyword中
        for(int keyword: p.keywords)
            leaf_unionKeyword[leaf].insert(keyword);
        //将该leaf（cluster）加入到查询叶节点组中
        store_leafSet.insert(leaf);
        store_nodes.insert(leaf);
        int current = leaf;
        GTree[current].stores.insert(store_id);
        int father;
        int i = 0;
        while (current != 0) {
            father = GTree[current].father;
            GTree[father].child_hasStores.insert(current);
            GTree[father].stores.insert(store_id);
            store_nodes.insert(father);
            current = father;
            i++;
        }
        //cout<<"层数"<<i<<",leaf="<<leaf<<endl;
    }
    cout<<"候选商店分布在"<<store_leafSet.size()<<"个不同叶节点中,GTree的"<<store_nodes.size()<<"个不同节点中"<<endl;
    //getchar();

    set<int> keyword_Set;
    for(int t: keywords)
        keyword_Set.insert(t);

    set<int> checkin_usr_Set;


    vector<int> Leaf_node_User,Leaf_node_Poi;

    vector<ResultDetail> result_User;
    //提取每个查询关键词 在leafNodeInv表中对应信息

    //对每个有usr的叶节点进行上下界的评估
    int prune=0; int unprune=0; int user_access=0; int usr_prune =0;
    int leaf_prune=0; int leaf_unprune=0;
    int may_prune = 0;
    int topK_count = 0; //topk执行次数
    // 获取与查询关键词相关的用户叶子节点
    set<int> related_usrLeafSet;
    int t3=0; int root=0;

    set<int> totalUsrSet;
    Filter_START
    int nodeAccess = 0;

    /*----------分别对每个query cluster执行Group Filtering------------------*/
    set<int> verification_User;  //候选验证用户集
    map<int, set<VerifyEntry>> verification_Map;    //凝练结果表   <leaf_id, <user_id....>>
    map<int, priority_queue<ResultCurrent>> candidate_usr_related_store;   //候选用户与查询对象的关联表
    //map<int, vector<ResultDetail>> result_Stores; //<store_id, <usr_id....>>
    map<int, LCLResult> candidate_usr_LCL;     //候选用户的LCL表
    cout<<"共"<<store_leafSet.size()<<"个query partition"<<endl; //getchar();
    TIME_TICK_START
    Filter_START
    map<int, map<int,bool>> usr_store_map;   //避免 user与object重复关联，重要！
    int p_count = 0;
    for(int store_leaf: store_leafSet){   //对每个leaf中的store进行处理
        //cout<< store_leaf<<endl;
        p_count++;
        //获得partition的 unionKeyword
        keywords = transforSet2Vector(leaf_unionKeyword[store_leaf]);
        vector<int> stores_inpartition =  leaf_store_list[store_leaf];   //store_leaf
#ifndef DiskAccess
        Group_Filter_Baseline_memory(usr_store_map, store_leaf, stores_inpartition, keywords, Qk, a, alpha, verification_User, verification_Map, candidate_usr_related_store, candidate_usr_LCL);
#else
        Group_Filter_Baseline_disk(usr_store_map, store_leaf, stores_inpartition, keywords, Qk, a, alpha, verification_User, verification_Map, candidate_usr_related_store, candidate_usr_LCL);
#endif
        //PAUSE_START
        cout<<"query partition "<<store_leaf<<"的filter处理完毕,剩余"<<(store_leafSet.size()-p_count)<<"个partition"<<endl;
        cout<<"当前找到"<<verification_User.size()<<"个候选用户"<<endl;

    }
    Filter_END
    //对所有candidate
    cout<<"所有 query partition 的filter 步骤处理完毕"<<endl;
    cout<<"累计总共找到"<<verification_User.size()<<"个候选用户"<<endl;

    Refine_START

    //对所有query_partition对应汇总的candidate进行验证处理
#ifndef DiskAccess
    BatchResults batch_query_results = Verification_BRkGSKQ_memory(verification_User, candidate_usr_related_store, keywords, Qk, a, alpha);
#else
    BatchResults batch_query_results = Verification_BRkGSKQ_disk(verification_User, candidate_usr_related_store, keywords, Qk, a, alpha);
#endif
    cout<<"有关联用户的query个数"<<batch_query_results.size()<<endl;  //baseyan

    Refine_END


    //对每个候选usr进行分析并打印结果
    //FilterResults results =  resultsAnalysis(result_User,keywords,checkin_usr,loc, Qk,a,alpha);
    //printBatchRkGSKQResults(stores, batch_query_results);
    cout<<"累计总共找到"<<verification_User.size()<<"个候选用户"<<endl;
    //cout<<"与兴趣点关键词相关用户总数"<< totalUsrSet.size() <<endl;
    cout<<"topk查询执行个数 =" << topK_count <<endl;
    cout<<"g-tree node prune ="<<prune<<" ,unprune="<<unprune<<",usr_prune="<<usr_prune<<endl;
    double ratio = (1.0*prune) /(prune+unprune);
    double ratio2 = (1.0*leaf_prune) /(leaf_prune+leaf_unprune);
    cout<<"节点剪枝率= "<<ratio<<"，其中叶子节点剪枝率="<<ratio2<<endl;
    //cout<<"g-tree node access ="<<nodeAccess<<endl;
    TIME_TICK_END
    Filter_PRINT("filter lasts")
    Refine_PRINT("refinement lasts")
    int io_count = block_num - block_num_pre;
    cout<<"i/o cost = "<<io_count<<endl;
    TIME_TICK_PRINT("MRkGSKQ Baseline runtime is ")


    return batch_query_results;

}

//some optimizitios
//1. 优化的过滤阶段， Access GIM-Tree  Once

//内存操作

void  Group_Filter_Once_memory(vector<int> stores, vector<int> keywords, int Qk, float a, float alpha, set<VerifyEntry>& verification_User, map<int,set<VerifyEntry>>& verification_Map, map<int, priority_queue<ResultCurrent>>& candidate_usr_related_store, map<int, LCLResult>& candidate_usr_LCL) {
    vector<int> qKey = keywords;
    set<int> store_leafSet;
    map<int, set<int>> partition_checkinMap;

    set<int> store_nodes;
    map<int, vector<int>> leaf_store_list;//<leaf_id, <store_id....>>  leaf 中的store的记录
    map<int, set<int>> leaf_verify_list;//<leaf_id, <usr_id....>>  leaf 中所有store对应的凝练用户集合


    BatchResults result_Stores; //<store_id, <usr_id....>>
    int loc =0;

    //generate the query tree for candidate stores
    for(int store_id:stores){
        int leaf = Nodes[POIs[store_id].Nj].gtreepath.back();
        leaf_store_list[leaf].push_back(store_id);
        //加入兴趣点的叶节点列表
        store_leafSet.insert(leaf);
        //叶子节点层partition的check-in信息
        for(int u: POIs[store_id].check_ins){
            partition_checkinMap[leaf].insert(u);
        }

        store_nodes.insert(leaf);
        int current = leaf;
        GTree[current].stores.insert(store_id);
        int father;
        int i = 0;
        while (current != 0) {

            father = GTree[current].father;
            GTree[father].child_hasStores.insert(current);
            GTree[father].stores.insert(store_id);
            store_nodes.insert(father);
            current = father;
            i++;
        }
        //cout<<"层数"<<i<<",leaf="<<leaf<<endl;
    }
    cout<<"候选商店分布在"<<store_leafSet.size()<<"个不同叶节点中,GTree的"<<store_nodes.size()<<"个不同节点中"<<endl;

    TIME_TICK_START
    set<int> keyword_Set;
    for(int t: keywords)
        keyword_Set.insert(t);


    //提取每个查询关键词 在leafNodeInv表中对应信息

    int prune=0; int unprune=0; int user_access=0; int usr_prune =0;
    int leaf_prune=0; int leaf_unprune=0;
    int may_prune = 0;
    //int topK_count = 0; //topk执行次数
    // 获取与查询关键词相关的用户叶子节点
    int t3=0; int root=0;


    set<int> totalUsrSet;
    int nodeAccess = 0;
    //查询处理过程中，关键词相关 用户节点列表
    vector<int> related_usr_node;
    //查询处理过程中，关键词相关 兴趣点节点列表
    map<int,vector<int>> related_poi_nodes;

    for(int term:keyword_Set)
        for(int usr: invListOfUser[term])
            totalUsrSet.insert(usr);

    cout<<"与关键词相关的用户个数 num="<<totalUsrSet.size();


    Filter_START

    //for(int term_id:keywords){
    // 获取与查询关键词相关的****用户***** 所在子图对应的上层节点

    vector<int> unodes = GTree[root].children; //get the child list of root whose text contains term_id
    for (int n : unodes) {
        //cout<<"node"<<n<<endl;
        related_usr_node.push_back(n);
        GTree[n].ossociate_queryNodes = store_leafSet; //先只是查询对象组叶节点
    }
    // 获取与查询关键词相关的****兴趣点***** 所在子图对应的上层节点
    for(int keyword: keyword_Set){
        for(int child: GTree[root].inverted_list_o[keyword]){
            related_poi_nodes[keyword].push_back(child);
        }
    }


    set<int> related_children;
    set<int> related_usr;
    set<int> related_leaves;
    //verify the user node from upper to bottom
    TIME_TICK_END
    TIME_TICK_PRINT("初始化完成，runtime")


    //verify the user node from upper to bottom
    int level = 1;  int cnt_ = 0;


    while (true) {  //daxunhuanonce
        //related_children.clear();
        if (related_usr_node.capacity() == 0) break;
        set<int>().swap(related_children);

        for (int node:related_usr_node) {


            set<int> inter_Key  = obtain_itersection_jins(GTree[node].userUKeySet, keyword_Set);
            if(!inter_Key.size()>0) continue;

            //该节点为用户上层节点
            if(!GTree[node].isleaf){

#ifdef TRACK
                cout << "access usr_node" << node << endl;
#endif
                int usrNum = 0;
                double global_max_Score = -1;
                //cout<<"non-leaf node："<< node <<endl;
                if (false) {   //用户稀疏节点  //x

                    for(int term: keyword_Set){
                        for(int usr_leaf: GTree[node].term_usrLeaf_Map[term]){
                            related_children.insert(usr_leaf);
                            for(int store_leaf:GTree[node].ossociate_queryNodes){
                                GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                            }
                        }

                    }
                    prune++;
                    //cout << "用户稀疏节点" << node <<",usrNum="<<usrNum<<endl;

                } else {  //用户密集型节点，需进行剪枝判断
                    //global upper bound  computation
                    vector<double> max_score_list;
                    for (int store_leaf: GTree[node].ossociate_queryNodes) {
                        //upper bound computation
                        double  max_Score =
                                getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, partition_checkinMap[store_leaf], store_leaf, node); //注意：这里用 InterKey, 与联合check-in计算的评分上界

                        global_max_Score = max( max_Score, global_max_Score);

                        max_score_list.push_back(max_Score);

                    }
                    //obtaining lcl from cache or lcl computation on the fly
                    double Rt_U = 100;
                    for(int term: inter_Key){
                        LCLResult lcl;
                        map<int,LCLResult>::iterator iter;
                        iter = GTree[node].termLCLDetail.find(term);
                        if(iter != GTree[node].termLCLDetail.end()){ //先访问cache
                            lcl = iter->second;
                        }else{
                            lcl = getUpperNode_termLCL(node,term, related_poi_nodes,Qk,a,alpha,global_max_Score);
                            GTree[node].termLCLDetail[term] = lcl;
                        }
                        double current_score = lcl.topscore;
                        Rt_U = min(Rt_U, current_score);
                    }
                    // group pruning judgement
                    if (global_max_Score > Rt_U) {  // 该子树下所有用户不能被剪枝，加入列表后续进一步验证
                        unprune++;
                        int _idx = 0;
                        for (int store_leaf: GTree[node].ossociate_queryNodes) {
                            //upper bound computation
                            double max_score = max_score_list[_idx];
                            _idx++;

                            if(Rt_U>max_score)   //满足pruning 条件，跳过
                                continue;

                            for(int term: inter_Key){
                                if(GTree[node].inverted_list_u.count(term)>0){
                                    for (int child: GTree[node].inverted_list_u[term]) {

                                        //孩子节点继承父节点的lcl内容
                                        GTree[child].termParentLCLDetail[term] = GTree[node].termParentLCLDetail[term];
                                        //加入后续候选用户
                                        related_children.insert(child);
                                        //添加 用户组与 对象组的 关联
                                        GTree[child].ossociate_queryNodes.insert(store_leaf);

                                    }
                                }
                            }

                        }


                    }
                        //整棵用户子树被剪枝
                    else {
                        prune++;
#ifdef PRUNELOG
                        cout << "节点" << node << "与所有 store都不相关， 其下所有用户可被剪枝" << endl;
#endif

                    }


                }//else用户密集型
            }

            else if (GTree[node].isleaf) {  //若该节点为用户叶子节点
#ifdef TRACK
                cout << "access usr_leaf" << node << endl;
#endif


                //if(node==360){
                //cout<<"find leaf360!"<<endl;
                //}

                MultiLCLResult lclResult_multi;
                vector<int> Ukeys;
                for(int term: inter_Key)
                    Ukeys.push_back(term);

                //global upper bound  computation
                double global_max_Score = -1;
                vector<double> max_score_list;
                for (int store_leaf: GTree[node].ossociate_queryNodes) {
                    //upper bound computation
                    double  max_Score =
                            getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, partition_checkinMap[store_leaf], store_leaf, node); //注意：这里用 InterKey, 与联合check-in计算的评分上界

                    global_max_Score = max( max_Score, global_max_Score);

                    max_score_list.push_back(max_Score);

                }

                //lcl computation, obtaining the lcl-r_t

                map<set<int>,MultiLCLResult>::iterator iter;
                iter = GTree[node].cacheMultiLCLDetail.find(inter_Key);
                if (iter!=GTree[node].cacheMultiLCLDetail.end()){
                    lclResult_multi = iter->second;
                }

                else{  //无cache  算之


                    //lclResult = getLeafNode_TopkSDijkstra(node, term_id, Qk, a, alpha, max_Score);
#ifdef LV
                    lclResult_multi = getLeafNode_LCL_Dijkstra_multi(node, Ukeys, Qk, a, alpha, global_max_Score);
#endif
#ifndef LV
                    lclResult_multi = getLeafNode_LCL_bottom2up_multi(node, Ukeys, Qk, a, alpha, global_max_Score);
#endif


                    GTree[node].cacheMultiLCLDetail[inter_Key] = lclResult_multi;
                }

                //获得 usr_lcl update_o_set
                set<int> update_o_set;
                for(LCLResult lcl: lclResult_multi){
                    priority_queue<LCLEntry> LCL = lcl.LCL;
                    while(LCL.size() > 0){
                        LCLEntry le = LCL.top();
                        LCL.pop();
                        int o = le.id;
                        update_o_set.insert(o);
                        //cout<<"o"<<o<<",";
                    }

                }

                //获得lcl-r_t的数值
                double Rt_U = -1.0;

                if(update_o_set.size()<Qk){
                    //cout<<"non-enough(empty) update o set"<<endl;
                    Rt_U = 0.0;
                }
                else{
                    for(LCLResult lclResult : lclResult_multi){
                        double _ss = lclResult.topscore;
                        if(-1.0==Rt_U)
                            Rt_U =_ss;
                        else
                            Rt_U = min(Rt_U,_ss);
                    }
                }

                if(global_max_Score > Rt_U){
                    // evaluate each patition
                    int  _idx = 0; set<int> L;
                    for (int store_leaf: GTree[node].ossociate_queryNodes) {  //每个partition在max_score_list中都有对应的upper bound
                        //upper bound computation
                        double  max_Score = max_score_list[_idx];
                        _idx++;

                        if( max_Score > Rt_U)
                            L.insert(store_leaf);
                        else{
                            continue;   //该store_leaf 被剪枝
                        }

                    }

                    //获取该叶子节点下各个用户
                    set<int> usr_withinLeaf;
                    for (int term: inter_Key) {
                        for (int usr_id: GTree[node].inverted_list_u[term]) {   //获取叶节点下每个用户
                            usr_withinLeaf.insert(usr_id);
                        }
                    }

                    //对该叶子节点下各个用户进行评估
                    for (int usr_id: usr_withinLeaf) {   //获取每个用户

                        double u_lcl_rt = 0.0;


                        priority_queue<ResultCurrent> u_related_store;

                        if(update_o_set.size()>0) {
                            //obtain the lcl list of a user from leaf lcl
                            priority_queue<ResultCurrent> u_LCL = get_usrLCL_update_Plus(
                                    usr_id, update_o_set, Qk, a, alpha, global_max_Score);  //获取该用户的lcl

                            if (u_LCL.size() < Qk) u_lcl_rt = 0.0;
                            else if (u_LCL.size() == Qk) u_lcl_rt = u_LCL.top().score;

                        }

                        Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);
                        for (int store_leaf: L) {
                            for (int store: GTree[store_leaf].stores) {  //获取每个store_leaf中的store
                                double rel = 0;
                                double inf = 0;

                                Score_Bound gsk_score_bound = getGSKScore_o2u_Upper_and_Lower_givenKey(a, alpha, POIs[store],
                                                                                                       Users[usr_id], keywords, rel,
                                                                                                       inf);
                                double gsk_score_Upper = gsk_score_bound.upper_score;
                                double gsk_score_Lower = gsk_score_bound.lower_score;
                                //TopkQueryCurrentResult result = topkSDijkstra_verify_usr(transformToQ(Users[usr_id]),Qk,a, alpha,gsk_score_Upper);


                                if (u_lcl_rt < gsk_score_Upper) {
                                    //update lcl
                                    u_related_store.push(ResultCurrent(store, gsk_score_Lower, gsk_score_Upper, rel,
                                                                       inf));  //记录relevance与influence 避免重算
                                    if (u_related_store.size() > Qk) {
                                        u_related_store.pop();
                                        u_lcl_rt = max(u_lcl_rt, u_related_store.top().score);
                                    }
                                    //cout<<"u"<<usr_id<<", gsk_score_Lower="<<gsk_score_Lower<<",u_lcl="<<u_lcl<<",加入store="<<store<<endl;
                                    Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);

                                } else {
                                    usr_prune++;
                                    //continue;   //该用户不是反top-k结果
                                }

                            }
                        }
                        if (u_related_store.size() > 0) {
                            VerifyEntry ve = VerifyEntry(usr_id,-1,u_lcl_rt);
                            verification_User.insert(ve);
                            if (candidate_usr_related_store[usr_id].size() > 0) {  //之前已经有u_lcl的结果了，更新之
                                priority_queue<ResultCurrent> u_related_Query_pre = candidate_usr_related_store[usr_id];
                                while (!u_related_store.empty()) {
                                    ResultCurrent rc = u_related_store.top();
                                    u_related_store.pop();
                                    u_related_Query_pre.push(rc);
                                    if (u_related_Query_pre.size() > Qk)
                                        u_related_Query_pre.pop();
                                }
                            } else {  //若无，直接赋予之
                                priority_queue<ResultCurrent> u_related_Query = u_related_store;
                                candidate_usr_related_store[usr_id] = u_related_Query;
                            }
                            //Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current,
                            //candidate_usr_related_store[usr_id].top().score);

#ifdef PRINTINFO
                            cout << "发现候选用户 u" << usr_id << ", u_lcl=" << u_lcl_rt << endl;
#endif

                            //getchar();
                        } else {
                            //cout << "u" << usr_id << "用户被剪枝, as Lu为空,u_lcl=" << u_lcl << endl;
                        }
                        //getchar();
                    }


                }
                    //该用户叶子节点整体被剪枝
                else {
#ifdef PRUNELOG
                    cout<<"leaf "<<store_leaf<<"与n"<<node<<"无关！"<<endl;
#endif
                    prune++;
                    leaf_prune++;
                }

                //cout<<"验证叶节点结束"<<node<<endl;
            }


        }// the first for loop
        level++;
        map<int, vector<int>> related_poi_nodes_next;
        map<int, vector<int>> related_poi_leaf_next;

        //upfold related poi nodes
        for(int term: keyword_Set){
            for(int pn: related_poi_nodes[term]){
                if(!GTree[pn].isleaf){
                    if(GTree[pn].inverted_list_o[term].size()!=0){
                        for(int child: GTree[pn].inverted_list_o[term]){
                            related_poi_nodes_next[term].push_back(child);
                        }
                    }
                }
                else{
                    related_poi_leaf_next[term].push_back(pn);
                }
            }
        }

        // upfold the related usr nodes
        vector<int>().swap(related_usr_node);
        related_usr_node.clear();

        for (int cn:related_children){       //jiachildyici

            related_usr_node.push_back(cn);

        }
        if(related_usr_node.size() ==0){
            //cout<<"related_usr_node空了"<<endl;
            break;

        }
        //防止单路分支带来的重复验证
        if(related_children.size()==1){
            cnt_++;
            if(cnt_==3){
                //cout<<"直接跳到叶节点"<<endl;
                //getchar();
                related_usr_node.clear();
                for(int node: related_children){
                    if(GTree[node].isleaf){

                        related_usr_node.push_back(node);  //

                    }

                    else{
                        set<int> leaves;
                        for(int term: keyword_Set){
                            if(!(GTree[node].term_usrLeaf_Map.count(term)>0))
                                continue;
                            for(int leaf: GTree[node].term_usrLeaf_Map[term]){
                                leaves.insert(leaf);
                            }
                        }
                        for(int usr_leaf: leaves){
                            related_usr_node.push_back(usr_leaf);
                            //加入关联关系
                            for (int store_leaf: GTree[node].ossociate_queryNodes) {
                                GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                            }
                        }

                    }


                }
                //cout<<"新加入"<<related_usr_node.size()<<endl;
                //getchar();

            }

        }
        //cout<<"第"<<level<<"层用户节点剪枝验证结束！"<<endl;
        //cout<<"related_poi_node size in next iteration:"<<related_poi_nodes.size()<<endl;
        //cout<<"related_usr_node size in next iteration:"<<related_usr_node.size()<<endl;
        if(related_usr_node.size()==0) break;

    }//end while

    //getchar();
}

void  Group_Filter_Once_disk_origin(vector<int> stores, vector<int> keywords, int Qk, float a, float alpha, set<VerifyEntry>& verification_User, map<int,set<VerifyEntry>>& verification_Map, map<int, priority_queue<ResultCurrent>>& candidate_usr_related_store, map<int, LCLResult>& candidate_usr_LCL) {
    vector<int> qKey = keywords;
    set<int> store_leafSet;
    map<int, set<int>> partition_checkinMap;

    set<int> store_nodes;
    map<int, vector<int>> leaf_store_list;//<leaf_id, <store_id....>>  leaf 中的store的记录
    map<int, set<int>> leaf_verify_list;//<leaf_id, <usr_id....>>  leaf 中所有store对应的凝练用户集合


    BatchResults result_Stores; //<store_id, <usr_id....>>
    int loc =0;

    //generate the query tree for candidate stores
    for(int store_id:stores){
        POI p = getPOIFromO2UOrgLeafData(store_id);
        int leaf = Nodes[p.Nj].gtreepath.back();
        leaf_store_list[leaf].push_back(store_id);
        //加入兴趣点的叶节点列表
        store_leafSet.insert(leaf);
        //叶子节点层partition的check-in信息
        for(int u: p.check_ins){
            partition_checkinMap[leaf].insert(u);
        }

        store_nodes.insert(leaf);
        int current = leaf;
        GTree[current].stores.insert(store_id);
        int father;
        int i = 0;
        while (current != 0) {

            father = GTree[current].father;
            GTree[father].child_hasStores.insert(current);
            GTree[father].stores.insert(store_id);
            store_nodes.insert(father);
            current = father;
            i++;
        }
        //cout<<"层数"<<i<<",leaf="<<leaf<<endl;
    }
    cout<<"候选商店分布在"<<store_leafSet.size()<<"个不同叶节点中,GTree的"<<store_nodes.size()<<"个不同节点中"<<endl;

    TIME_TICK_START
    set<int> keyword_Set;
    for(int t: keywords)
        keyword_Set.insert(t);


    //提取每个查询关键词 在leafNodeInv表中对应信息

    int upper_prune=0; int upper_unprune=0;
    int leaf_prune=0; int leaf_unprune=0;
    int aggressive_prune =0;
    int user_access=0; int usr_prune =0;
    int may_prune = 0;
    //int topK_count = 0; //topk执行次数
    // 获取与查询关键词相关的用户叶子节点
    int t3=0; int root=0;


    set<int> totalUsrSet;
    int nodeAccess = 0;
    //查询处理过程中，关键词相关 用户节点列表
    vector<int> related_usr_node;
    //查询处理过程中，关键词相关 兴趣点节点列表
    map<int,vector<int>> related_poi_nodes;

    /*for(int term:keyword_Set)
        for(int usr: invListOfUser[term])
            totalUsrSet.insert(usr);

    cout<<"与关键词相关的用户个数 num="<<totalUsrSet.size();*/


    Filter_START


    // 获取与查询关键词相关的****用户***** 所在子图对应的上层节点
    vector<int> unodes = GTree[root].children; //get the child list of root whose text contains term_id
    for (int n : unodes) {
        //cout<<"node"<<n<<endl;
        related_usr_node.push_back(n);
        GTree[n].ossociate_queryNodes = store_leafSet; //先只是查询对象组叶节点
    }
    // 获取与查询关键词相关的****兴趣点***** 所在子图对应的上层节点
    for(int keyword: keyword_Set){
        vector<int> child_entry_o = getObjectTermRelatedEntry(keyword,root);
        for(int child: child_entry_o){
            related_poi_nodes[keyword].push_back(child);
        }
    }


    set<int> related_children;
    set<int> related_usr;
    set<int> related_leaves;
    //verify the user node from upper to bottom
    TIME_TICK_END
    TIME_TICK_PRINT("初始化完成，runtime")


    //verify the user node from upper to bottom
    int level = 1;  int cnt_ = 0;


    while (true) {  //daxunhuanonce
        //related_children.clear();
        if (related_usr_node.capacity() == 0) break;
        set<int>().swap(related_children);

        for (int node:related_usr_node) {

            TreeNode user_node = getGIMTreeNodeData(node,OnlyU);

            set<int> inter_Key  = obtain_itersection_jins(user_node.userUKeySet, keyword_Set);
            if(!inter_Key.size()>0) continue;

            //该节点为用户上层节点
            if(!user_node.isleaf){

                UPPER_START
                int usrNum = 0;
                double global_max_Score = -1;
#ifdef TRACK
                cout<<"non-leaf node："<< node <<endl;
#endif
                for(int cover_key: inter_Key){

                }

                if (false) {   //用户稀疏节点  //x

                    for(int term: keyword_Set){
                        vector<int> usr_leaf_descend = getUsrTermRelatedLeafNode(term,node);
                        for(int usr_leaf: usr_leaf_descend){
                            related_children.insert(usr_leaf);
                            for(int store_leaf:GTree[node].ossociate_queryNodes){
                                GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                            }
                        }

                    }
                    //upper_prune++;
                    //cout << "用户稀疏节点" << node <<",usrNum="<<usrNum<<endl;

                } else {  //用户密集型节点，需进行剪枝判断
                    //global upper bound  computation
                    vector<double> max_score_list;
                    for (int store_leaf: GTree[node].ossociate_queryNodes) {
                        //upper bound computation
                        double  max_Score =
                                getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, partition_checkinMap[store_leaf], store_leaf, node); //注意：这里用 InterKey, 与联合check-in计算的评分上界

                        global_max_Score = max( max_Score, global_max_Score);

                        max_score_list.push_back(max_Score);

                    }
                    //obtaining lcl from cache or lcl computation on the fly
                    double Rt_U = 100;
                    for(int term: inter_Key){
                        LCLResult lcl;
                        map<int,LCLResult>::iterator iter;
                        iter = GTree[node].termLCLDetail.find(term);
                        if(iter != GTree[node].termLCLDetail.end()){ //先访问cache
                            lcl = iter->second;
                        }else{
                            lcl = getUpperNode_termLCL_Disk(node,term, related_poi_nodes,Qk,a,alpha,global_max_Score);
                            GTree[node].termLCLDetail[term] = lcl;
                        }
                        double current_score = lcl.topscore;
                        Rt_U = min(Rt_U, current_score);
                    }
                    // group pruning judgement
                    if (global_max_Score > Rt_U) {  // 该子树下所有用户不能被剪枝，加入列表后续进一步验证
                        upper_unprune++;
                        int _idx = 0;
                        for (int store_leaf: GTree[node].ossociate_queryNodes) {
                            //upper bound computation
                            double max_score = max_score_list[_idx];
                            _idx++;

                            if(Rt_U>max_score)   //满足pruning 条件，跳过
                                continue;

                            for(int term: inter_Key){
                                vector<int> child_entry_u = getUsrTermRelatedEntry(term,node);
                                if(child_entry_u.size()>0){
                                    for (int child: child_entry_u) {

                                        //孩子节点继承父节点的lcl内容
                                        GTree[child].termParentLCLDetail[term] = GTree[node].termParentLCLDetail[term];
                                        //加入后续候选用户
                                        related_children.insert(child);
                                        //添加 用户组与 对象组的 关联
                                        GTree[child].ossociate_queryNodes.insert(store_leaf);

                                    }
                                }
                            }

                        }


                    }
                        //整棵用户子树被剪枝
                    else {
                        upper_prune++;
#ifdef PRUNELOG
                        cout << "节点" << node << "与所有 store都不相关， 其下所有用户可被剪枝" << endl;
#endif

                    }


                }//else用户密集型
                UPPER_PAUSE
            }

            else if (user_node.isleaf) {  //若该节点为用户叶子节点
#ifdef TRACK
                cout << "access usr_leaf" << node << endl;
#endif


                MultiLCLResult lclResult_multi;
                vector<int> Ukeys;
                for(int term: inter_Key)
                    Ukeys.push_back(term);

                //global upper bound  computation
                double global_max_Score = -1;
                vector<double> max_score_list;
                for (int store_leaf: GTree[node].ossociate_queryNodes) {
                    //upper bound computation
                    double  max_Score =
                            getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, partition_checkinMap[store_leaf], store_leaf, node); //注意：这里用 InterKey, 与联合check-in计算的评分上界

                    global_max_Score = max( max_Score, global_max_Score);

                    max_score_list.push_back(max_Score);

                }

                //lcl computation, obtaining the lcl-r_t

                map<set<int>,MultiLCLResult>::iterator iter;
                iter = GTree[node].cacheMultiLCLDetail.find(inter_Key);
                if (iter!=GTree[node].cacheMultiLCLDetail.end()){
                    lclResult_multi = iter->second;
                }

                else{  //无cache  算之


                    //lclResult = getLeafNode_TopkSDijkstra(node, term_id, Qk, a, alpha, max_Score);
#ifdef LV
                    lclResult_multi = getLeafNode_LCL_Dijkstra_multi(node, Ukeys, Qk, a, alpha, global_max_Score);
#endif
#ifndef LV
                    lclResult_multi = getLeafNode_LCL_bottom2up_multi(node, Ukeys, Qk, a, alpha, global_max_Score);
#endif


                    GTree[node].cacheMultiLCLDetail[inter_Key] = lclResult_multi;
                }

                //获得 usr_lcl update_o_set
                set<int> update_o_set;
                for(LCLResult lcl: lclResult_multi){
                    priority_queue<LCLEntry> LCL = lcl.LCL;
                    while(LCL.size() > 0){
                        LCLEntry le = LCL.top();
                        LCL.pop();
                        int o = le.id;
                        update_o_set.insert(o);
                        //cout<<"o"<<o<<",";
                    }

                }

                //获得lcl-r_t的数值
                double Rt_U = -1.0;

                if(update_o_set.size()<Qk){
                    //cout<<"non-enough(empty) update o set"<<endl;
                    Rt_U = 0.0;
                }
                else{
                    for(LCLResult lclResult : lclResult_multi){
                        double _ss = lclResult.topscore;
                        if(-1.0==Rt_U)
                            Rt_U =_ss;
                        else
                            Rt_U = min(Rt_U,_ss);
                    }
                }

                if(global_max_Score > Rt_U){
                    leaf_unprune;
                    // evaluate each patition
                    int  _idx = 0; set<int> L;
                    for (int store_leaf: GTree[node].ossociate_queryNodes) {  //每个partition在max_score_list中都有对应的upper bound
                        //upper bound computation
                        double  max_Score = max_score_list[_idx];
                        _idx++;

                        if( max_Score > Rt_U)
                            L.insert(store_leaf);
                        else{
                            continue;   //该store_leaf 被剪枝
                        }

                    }

                    //获取该叶子节点下各个用户
                    set<int> usr_withinLeaf;
                    for (int term: inter_Key) {
                        vector<int> usr_list = getUsrTermRelatedEntry(term,node);
                        for (int usr_id: usr_list) {   //获取叶节点下每个用户
                            usr_withinLeaf.insert(usr_id);
                        }
                    }

                    //对该叶子节点下各个用户进行评估
                    for (int usr_id: usr_withinLeaf) {   //获取每个用户

                        double u_lcl_rt = 0.0;


                        priority_queue<ResultCurrent> u_related_store;

                        if(update_o_set.size()>0) {
                            //obtain the lcl list of a user from leaf lcl
                            priority_queue<ResultCurrent> u_LCL = get_usrLCL_update_Plus(
                                    usr_id, update_o_set, Qk, a, alpha, global_max_Score);  //获取该用户的lcl

                            if (u_LCL.size() < Qk) u_lcl_rt = 0.0;
                            else if (u_LCL.size() == Qk) u_lcl_rt = u_LCL.top().score;

                        }

                        Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);
                        for (int store_leaf: L) {
                            for (int store: GTree[store_leaf].stores) {  //获取每个store_leaf中的store
                                double rel = 0;
                                double inf = 0;
                                POI p = getPOIFromO2UOrgLeafData(store);
                                User u = getUserFromO2UOrgLeafData(usr_id);
#ifdef GIVENKEY
                                Score_Bound gsk_score_bound = getGSKScore_o2u_Upper_and_Lower_givenKey(a, alpha, p,
                                                                                                       u, keywords, rel,
                                                                                                       inf);
#else
                                Score_Bound gsk_score_bound = getGSKScore_o2u_Upper_and_Lower(a, alpha, p,
                                                                                              u, rel,inf);
#endif
                                double gsk_score_Upper = gsk_score_bound.upper_score;
                                double gsk_score_Lower = gsk_score_bound.lower_score;
                                //TopkQueryCurrentResult result = topkSDijkstra_verify_usr(transformToQ(Users[usr_id]),Qk,a, alpha,gsk_score_Upper);


                                if (u_lcl_rt < gsk_score_Upper) {
                                    //update lcl
                                    u_related_store.push(ResultCurrent(store, gsk_score_Lower, gsk_score_Upper, rel,
                                                                       inf));  //记录relevance与influence 避免重算
                                    if (u_related_store.size() > Qk) {
                                        u_related_store.pop();
                                        u_lcl_rt = max(u_lcl_rt, u_related_store.top().score);
                                    }
                                    //cout<<"u"<<usr_id<<", gsk_score_Lower="<<gsk_score_Lower<<",u_lcl="<<u_lcl<<",加入store="<<store<<endl;
                                    Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);

                                } else {
                                    usr_prune++;
                                    //continue;   //该用户不是反top-k结果
                                }

                            }
                        }
                        if (u_related_store.size() > 0) {
                            VerifyEntry ve = VerifyEntry(usr_id,-1,u_lcl_rt);
                            verification_User.insert(ve);
                            if (candidate_usr_related_store[usr_id].size() > 0) {  //之前已经有u_lcl的结果了，更新之
                                priority_queue<ResultCurrent> u_related_Query_pre = candidate_usr_related_store[usr_id];
                                while (!u_related_store.empty()) {
                                    ResultCurrent rc = u_related_store.top();
                                    u_related_store.pop();
                                    u_related_Query_pre.push(rc);
                                    if (u_related_Query_pre.size() > Qk)
                                        u_related_Query_pre.pop();
                                }
                            } else {  //若无，直接赋予之
                                priority_queue<ResultCurrent> u_related_Query = u_related_store;
                                candidate_usr_related_store[usr_id] = u_related_Query;
                            }


#ifdef PRINTINFO
                            cout << "发现候选用户 u" << usr_id << ", u_lcl=" << u_lcl_rt << endl;
#endif

                            //getchar();
                        } else {
                            //cout << "u" << usr_id << "用户被剪枝, as Lu为空,u_lcl=" << u_lcl << endl;
                        }
                        //getchar();
                    }


                }
                    //该用户叶子节点整体被剪枝
                else {
#ifdef PRUNELOG
                    cout<<"leaf "<<store_leaf<<"与n"<<node<<"无关！"<<endl;
#endif
                    //prune++;
                    leaf_prune++;
                }

                //cout<<"验证叶节点结束"<<node<<endl;
            }


        }// the first for loop
        level++;
        map<int, vector<int>> related_poi_nodes_next;
        map<int, vector<int>> related_poi_leaf_next;

        //upfold related poi nodes
        for(int term: keyword_Set){
            for(int pn: related_poi_nodes[term]){
                if(!GTree[pn].isleaf){
                    vector<int> child_entry_o = getObjectTermRelatedEntry(term,pn);
                    if(child_entry_o.size()!=0){
                        for(int child: child_entry_o){
                            related_poi_nodes_next[term].push_back(child);
                        }
                    }
                }
                else{
                    related_poi_leaf_next[term].push_back(pn);
                }
            }
        }

        // upfold the related usr nodes
        vector<int>().swap(related_usr_node);
        related_usr_node.clear();

        for (int cn:related_children){       //jiachildyici

            related_usr_node.push_back(cn);

        }
        if(related_usr_node.size() ==0){
            //cout<<"related_usr_node空了"<<endl;
            break;

        }
        //防止单路分支带来的重复验证
        if(related_children.size()==1){
            cnt_++;
            if(cnt_==3){
                //cout<<"直接跳到叶节点"<<endl;
                //getchar();
                related_usr_node.clear();
                for(int node: related_children){
                    if(GTree[node].isleaf){
                        related_usr_node.push_back(node);  //

                    }

                    else{
                        set<int> leaves;
                        for(int term: keyword_Set){
                            vector<int> usr_leaf = getUsrTermRelatedEntry(term,node);
                            for(int leaf: usr_leaf){
                                leaves.insert(leaf);
                            }
                        }
                        for(int usr_leaf: leaves){
                            related_usr_node.push_back(usr_leaf);
                            //加入关联关系
                            for (int store_leaf: GTree[node].ossociate_queryNodes) {
                                GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                            }
                        }

                    }


                }
                //cout<<"新加入"<<related_usr_node.size()<<endl;
                //getchar();

            }

        }
        //cout<<"第"<<level<<"层用户节点剪枝验证结束！"<<endl;
        //cout<<"related_poi_node size in next iteration:"<<related_poi_nodes.size()<<endl;
        //cout<<"related_usr_node size in next iteration:"<<related_usr_node.size()<<endl;
        if(related_usr_node.size()==0) break;

    }//end while

    //getchar();
}

void  Group_Filter_Once_disk(vector<int> stores, vector<int> keywords, int Qk, float a, float alpha, set<VerifyEntry>& verification_User, map<int,set<VerifyEntry>>& verification_Map, map<int, priority_queue<ResultCurrent>>& candidate_usr_related_store, map<int, LCLResult>& candidate_usr_LCL) {
    vector<int> qKey = keywords;
    set<int> store_leafSet;
    map<int, set<int>> partition_checkinMap;

    set<int> store_nodes;
    map<int, vector<int>> leaf_store_list;//<leaf_id, <store_id....>>  leaf 中的store的记录
    map<int, set<int>> leaf_verify_list;//<leaf_id, <usr_id....>>  leaf 中所有store对应的凝练用户集合


    BatchResults result_Stores; //<store_id, <usr_id....>>
    int loc =0;

    //generate the query tree for candidate stores
    for(int store_id:stores){
        POI p = getPOIFromO2UOrgLeafData(store_id);
        int leaf = Nodes[p.Nj].gtreepath.back();
        leaf_store_list[leaf].push_back(store_id);
        //加入兴趣点的叶节点列表
        store_leafSet.insert(leaf);
        //叶子节点层partition的check-in信息
        for(int u: p.check_ins){
            partition_checkinMap[leaf].insert(u);
        }

        store_nodes.insert(leaf);
        int current = leaf;
        GTree[current].stores.insert(store_id);
        int father;
        int i = 0;
        while (current != 0) {

            father = GTree[current].father;
            GTree[father].child_hasStores.insert(current);
            GTree[father].stores.insert(store_id);
            store_nodes.insert(father);
            current = father;
            i++;
        }
        //cout<<"层数"<<i<<",leaf="<<leaf<<endl;
    }
    cout<<"候选商店分布在"<<store_leafSet.size()<<"个不同叶节点中,GTree的"<<store_nodes.size()<<"个不同节点中"<<endl;

    TIME_TICK_START
    set<int> keyword_Set;
    for(int t: keywords){  //
        keyword_Set.insert(t);
    }



    //提取每个查询关键词 在leafNodeInv表中对应信息

    int upper_prune=0; int upper_unprune=0;
    int leaf_prune=0; int leaf_unprune=0;
    int aggressive_prune =0;
    int user_access=0; int usr_prune =0;
    int may_prune = 0;
    //int topK_count = 0; //topk执行次数
    // 获取与查询关键词相关的用户叶子节点
    int t3=0; int root=0;


    set<int> totalUsrSet;
    int nodeAccess = 0;
    //查询处理过程中，关键词相关 用户节点列表
    vector<int> related_usr_node;
    //查询处理过程中，关键词相关 兴趣点节点列表
    map<int,vector<int>> related_poi_nodes;

    /*for(int term:keyword_Set)
        for(int usr: invListOfUser[term])
            totalUsrSet.insert(usr);

    cout<<"与关键词相关的用户个数 num="<<totalUsrSet.size();*/


    Filter_START


    // 获取与查询关键词相关的****用户***** 所在子图对应的上层节点
    vector<int> unodes = GTree[root].children; //get the child list of root whose text contains term_id
    for (int n : unodes) {
        related_usr_node.push_back(n);
        GTree[n].ossociate_queryNodes = store_leafSet; //先只是查询对象组叶节点
        /*printf("GTree[%d].ossociate_queryNodes\n:",n);
        printSetElements(GTree[n].ossociate_queryNodes);*/
    }
    // 获取与查询关键词相关的****兴趣点***** 所在子图对应的上层节点
    for(int keyword: keyword_Set){
        vector<int> child_entry_o = getObjectTermRelatedEntry(keyword,root);
        for(int child: child_entry_o){
            related_poi_nodes[keyword].push_back(child);
        }
    }


    set<int> related_children;
    set<int> related_usr;
    set<int> related_leaves;
    //verify the user node from upper to bottom
    TIME_TICK_END
    TIME_TICK_PRINT("初始化完成，runtime")


    //verify the user node from upper to bottom
    int level = 1;  int cnt_ = 0;
    int leafLcl_cnt = 0;

    while (true) {  //daxunhuanonce
        //related_children.clear();
        if (related_usr_node.capacity() == 0) break;
        set<int>().swap(related_children);

        for (int node:related_usr_node) {

            TreeNode user_node = getGIMTreeNodeData(node,OnlyU);

            set<int> inter_Key  = obtain_itersection_jins(user_node.userUKeySet, keyword_Set);
            if(!inter_Key.size()>0) continue;

            //该节点为用户上层节点
            if(!user_node.isleaf){

                UPPER_START
                int usrNum = 0;
                double global_max_Score = -1;
#ifdef TRACK
                cout<<"non-leaf node："<< node <<endl;
#endif
                int user_capacity = -1;
                for(int keyword: inter_Key){
                    int u_cnt = user_node.term_usr_countMap[keyword];
                    user_capacity = max(user_capacity,u_cnt);
                }
                bool cardinality_flag = false;
                bool aggressive_flag =false;
                set<int> usr_leaf_direct;
#ifdef  AGGRESSIVE_I
                cardinality_flag = (user_capacity<2*SPARCITY_VALUE);
                if(cardinality_flag){
                    for(int term: inter_Key){
                        vector<int> usr_leaf_descend = getUsrTermRelatedLeafNode(term,node);
                        for(int usr_leaf: usr_leaf_descend){
                            usr_leaf_direct.insert(usr_leaf);
                        }
                    }
                    if(usr_leaf_direct.size()<SPARCITY_VALUE)
                        aggressive_flag = true;
                }

#endif

                if (aggressive_flag) {   //用户稀疏节点  //aggressive prune //user_capacity<4

                    for(int usr_leaf: usr_leaf_direct){
                        related_children.insert(usr_leaf);
                        for(int store_leaf:GTree[node].ossociate_queryNodes){
                            GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                        }

                    }
                    aggressive_prune++;
                    //cout << "用户稀疏节点" << node <<",usrNum="<<usrNum<<endl;
#ifdef TRACK
                    cout << "用户稀疏节点" << node <<",usrNum="<<usrNum<<endl;
#endif

                } else {  //用户密集型节点，需进行剪枝判断
#ifdef TRACK
                    cout << "用户密集型节点，需进行剪枝判断"<<endl;
#endif
                    //global upper bound  computation
                    vector<double> max_score_list;

                    /*printf("GTree[node]: n%d 's ossociate_queryNodes:\n",node);
                    printSetElements(GTree[node].ossociate_queryNodes);
                    getchar();*/

                    for (int store_leaf: GTree[node].ossociate_queryNodes) {
                        //upper bound computation
                        double  max_Score =
                                getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, partition_checkinMap[store_leaf], store_leaf, node); //注意：这里用 InterKey, 与联合check-in计算的评分上界
                        //cout<<"max_Score="<<max_Score<<endl;

                        global_max_Score = max( max_Score, global_max_Score);

                        max_score_list.push_back(max_Score);

                    }
                    //obtaining lcl from cache or lcl computation on the fly
                    double Rt_U = 100;
                    for(int term: inter_Key){
                        LCLResult lcl;
                        map<int,LCLResult>::iterator iter;
                        iter = GTree[node].termLCLDetail.find(term);
                        if(iter != GTree[node].termLCLDetail.end()){ //先访问cache
                            lcl = iter->second;
                        }else{
                            lcl = getUpperNode_termLCL_Disk(node,term, related_poi_nodes,Qk,a,alpha,global_max_Score);
                            GTree[node].termLCLDetail[term] = lcl;
                        }
                        double current_score = lcl.topscore;
                        Rt_U = min(Rt_U, current_score);
                    }
                    //cout<<"Rt_U="<<Rt_U<<endl;
                    // group pruning judgement
                    if (global_max_Score > Rt_U) {  // 该子树下所有用户不能被剪枝，加入列表后续进一步验证
                        upper_unprune++;

                        int _idx = 0;
                        for (int store_leaf: GTree[node].ossociate_queryNodes) {
                            //upper bound computation
                            double max_score = max_score_list[_idx];
                            _idx++;

                            if(Rt_U>max_score)   //满足pruning 条件，跳过
                                continue;

                            for(int term: inter_Key){
                                vector<int> child_entry_u = getUsrTermRelatedEntry(term,node);
                                if(child_entry_u.size()>0){
                                    for (int child: child_entry_u) {

                                        //孩子节点继承父节点的lcl内容
                                        GTree[child].termParentLCLDetail[term] = GTree[node].termParentLCLDetail[term];
                                        //加入后续候选用户
                                        related_children.insert(child);
                                        //添加 用户组与 对象组的 关联
                                        GTree[child].ossociate_queryNodes.insert(store_leaf);

                                    }
                                }
                            }

                        }


                    }
                        //整棵用户子树被剪枝
                    else {
                        upper_prune++;
#ifdef PRUNELOG
                        cout << "节点" << node << "与所有 store都不相关， 其下所有用户可被剪枝" << endl;
#endif

                    }


                }//else用户密集型
                UPPER_PAUSE
            }

            else if (user_node.isleaf) {  //若该节点为用户叶子节点
                LEAF_START
#ifdef TRACK
                cout << "access usr_leaf" << node << endl;
#endif


                MultiLCLResult lclResult_multi;
                vector<int> Ukeys;
                for(int term: inter_Key)
                    Ukeys.push_back(term);

                //global upper bound  computation
                double global_max_Score = -1;
                vector<double> max_score_list;
                for (int store_leaf: GTree[node].ossociate_queryNodes) {
                    //upper bound computation
                    double  max_Score =
                            getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, partition_checkinMap[store_leaf], store_leaf, node); //注意：这里用 InterKey, 与联合check-in计算的评分上界

                    global_max_Score = max( max_Score, global_max_Score);

                    max_score_list.push_back(max_Score);

                }

                //lcl computation, obtaining the lcl-r_t

                map<set<int>,MultiLCLResult>::iterator iter;
                iter = GTree[node].cacheMultiLCLDetail.find(inter_Key);
                if (iter!=GTree[node].cacheMultiLCLDetail.end()){
                    lclResult_multi = iter->second;
                }

                else{  //无cache  算之


                    lclResult_multi = getLeafNode_LCL_Dijkstra_multi(node, Ukeys, Qk, a, alpha, global_max_Score);
/*
#ifdef LV
                    lclResult_multi = getLeafNode_LCL_Dijkstra_multi(node, Ukeys, Qk, a, alpha, global_max_Score);
#endif
#ifndef LV
        #ifdef LasVegas
                    lclResult_multi = getLeafNode_LCL_Dijkstra_multi(node, Ukeys, Qk, a, alpha, global_max_Score);
        #else
                    lclResult_multi = getLeafNode_LCL_bottom2up_multi(node, Ukeys, Qk, a, alpha, global_max_Score);
        #endif

#endif
*/


                    GTree[node].cacheMultiLCLDetail[inter_Key] = lclResult_multi;

                    leafLcl_cnt++;
                }

                //获得 usr_lcl update_o_set
                set<int> update_o_set;
                for(LCLResult lcl: lclResult_multi){
                    priority_queue<LCLEntry> LCL = lcl.LCL;
                    while(LCL.size() > 0){
                        LCLEntry le = LCL.top();
                        LCL.pop();
                        int o = le.id;
                        update_o_set.insert(o);
                        //cout<<"o"<<o<<",";
                    }

                }

                //获得lcl-r_t的数值
                double Rt_U = -1.0;

                if(update_o_set.size()<Qk){
                    //cout<<"non-enough(empty) update o set"<<endl;
                    Rt_U = 0.0;
                }
                else{
                    for(LCLResult lclResult : lclResult_multi){
                        double _ss = lclResult.topscore;
                        if(-1.0==Rt_U)
                            Rt_U =_ss;
                        else
                            Rt_U = min(Rt_U,_ss);
                    }
                }

                if(global_max_Score > Rt_U){
                    leaf_unprune++;
                    // evaluate each patition
                    int  _idx = 0; set<int> L;
#ifdef LEAF_UNPRUNELOG_
                    if(true){ //global_max_Score<2*Rt_U
                        cout<<"usr_leaf"<<node<<",global_max_Score="<<global_max_Score<<", Rt_u="<<Rt_U<<endl;
                    }
#endif

                    for (int store_leaf: GTree[node].ossociate_queryNodes) {  //每个partition在max_score_list中都有对应的upper bound
                        //upper bound computation
                        double  max_Score = max_score_list[_idx];
                        _idx++;

                        if( max_Score > Rt_U)
                            L.insert(store_leaf);
                        else{
                            continue;   //该store_leaf 被剪枝
                        }

                    }

                    //获取该叶子节点下各个用户
                    set<int> usr_withinLeaf;
                    for (int term: inter_Key) {
                        vector<int> usr_list = getUsrTermRelatedEntry(term,node);
                        for (int usr_id: usr_list) {   //获取叶节点下每个用户
                            usr_withinLeaf.insert(usr_id);
                        }
                    }

                    //对该叶子节点下各个用户进行评估
                    for (int usr_id: usr_withinLeaf) {   //获取每个用户

                        double u_lcl_rt = 0.0;


                        priority_queue<ResultCurrent> u_related_store;

                        if(update_o_set.size()>0) {
                            //obtain the lcl list of a user from leaf lcl
                            priority_queue<ResultCurrent> u_LCL = get_usrLCL_update_Plus(
                                    usr_id, update_o_set, Qk, a, alpha, global_max_Score);  //获取该用户的lcl

                            if (u_LCL.size() < Qk) u_lcl_rt = 0.0;
                            else if (u_LCL.size() == Qk) u_lcl_rt = u_LCL.top().score;

                        }

                        Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);
                        for (int store_leaf: L) {
                            for (int store: GTree[store_leaf].stores) {  //获取每个store_leaf中的store
                                double rel = 0;
                                double inf = 0;
                                POI p = getPOIFromO2UOrgLeafData(store);
                                User u = getUserFromO2UOrgLeafData(usr_id);
#ifdef GIVENKEY
                                Score_Bound gsk_score_bound = getGSKScore_o2u_Upper_and_Lower_givenKey(a, alpha, p,
                                                                                                       u, keywords, rel,
                                                                                                       inf);
#else
                                Score_Bound gsk_score_bound = getGSKScore_o2u_Upper_and_Lower(a, alpha, p,
                                                                                              u, rel,inf);
#endif
                                double gsk_score_Upper = gsk_score_bound.upper_score;
                                double gsk_score_Lower = gsk_score_bound.lower_score;
                                //TopkQueryCurrentResult result = topkSDijkstra_verify_usr(transformToQ(Users[usr_id]),Qk,a, alpha,gsk_score_Upper);


                                if (u_lcl_rt < gsk_score_Upper) {
                                    //update lcl
                                    u_related_store.push(ResultCurrent(store, gsk_score_Lower, gsk_score_Upper, rel,
                                                                       inf));  //记录relevance与influence 避免重算
                                    if (u_related_store.size() > Qk) {
                                        u_related_store.pop();
                                        u_lcl_rt = max(u_lcl_rt, u_related_store.top().score);
                                    }
                                    //cout<<"u"<<usr_id<<", gsk_score_Lower="<<gsk_score_Lower<<",u_lcl="<<u_lcl<<",加入store="<<store<<endl;
                                    Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);

                                } else {
                                    usr_prune++;
                                    //continue;   //该用户不是反top-k结果
                                }

                            }
                        }
                        if (u_related_store.size() > 0) {
                            VerifyEntry ve = VerifyEntry(usr_id,-1,u_lcl_rt);
                            verification_User.insert(ve);
                            if (candidate_usr_related_store[usr_id].size() > 0) {  //之前已经有u_lcl的结果了，更新之
                                priority_queue<ResultCurrent> u_related_Query_pre = candidate_usr_related_store[usr_id];
                                while (!u_related_store.empty()) {
                                    ResultCurrent rc = u_related_store.top();
                                    u_related_store.pop();
                                    u_related_Query_pre.push(rc);
                                    if (u_related_Query_pre.size() > Qk)
                                        u_related_Query_pre.pop();
                                }
                            } else {  //若无，直接赋予之
                                priority_queue<ResultCurrent> u_related_Query = u_related_store;
                                candidate_usr_related_store[usr_id] = u_related_Query;
                            }


#ifdef PRINTINFO
                            cout << "发现候选用户 u" << usr_id << ", u_lcl=" << u_lcl_rt << endl;
#endif

                            //getchar();
                        } else {
                            //cout << "u" << usr_id << "用户被剪枝, as Lu为空,u_lcl=" << u_lcl << endl;
                        }
                        //getchar();
                    }


                }
                    //该用户叶子节点整体被剪枝
                else {
#ifdef LEAF_PRUNELOG_
                    cout<<"uleaf "<<node<<"与所有查询对象都无关！"<<endl;
#endif
                    //prune++;
                    leaf_prune++;
                }

                //cout<<"验证叶节点结束"<<node<<endl;
                LEAF_PAUSE
            }


        }// the first for loop
        level++;
        map<int, vector<int>> related_poi_nodes_next;
        map<int, vector<int>> related_poi_leaf_next;

        //upfold related poi nodes
        for(int term: keyword_Set){
            for(int pn: related_poi_nodes[term]){
                if(!GTree[pn].isleaf){
                    vector<int> child_entry_o = getObjectTermRelatedEntry(term,pn);
                    if(child_entry_o.size()!=0){
                        for(int child: child_entry_o){
                            related_poi_nodes_next[term].push_back(child);
                        }
                    }
                }
                else{
                    related_poi_leaf_next[term].push_back(pn);
                }
            }
        }

        // upfold the related usr nodes
        vector<int>().swap(related_usr_node);
        related_usr_node.clear();

        for (int cn:related_children){       //jiachildyici

            related_usr_node.push_back(cn);

        }
        if(related_usr_node.size() ==0){
            //cout<<"related_usr_node空了"<<endl;
            break;

        }
        //防止单路分支带来的重复验证
        if(related_children.size()==1){
            cnt_++;
            if(cnt_==3){
                //cout<<"直接跳到叶节点"<<endl;
                //getchar();
                related_usr_node.clear();
                for(int node: related_children){
                    if(GTree[node].isleaf){
                        related_usr_node.push_back(node);  //

                    }

                    else{
                        set<int> leaves;
                        for(int term: keyword_Set){
                            vector<int> usr_leaf = getUsrTermRelatedEntry(term,node);
                            for(int leaf: usr_leaf){
                                leaves.insert(leaf);
                            }
                        }
                        for(int usr_leaf: leaves){
                            related_usr_node.push_back(usr_leaf);
                            //加入关联关系
                            for (int store_leaf: GTree[node].ossociate_queryNodes) {
                                GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                            }
                        }

                    }


                }
                //cout<<"新加入"<<related_usr_node.size()<<endl;
                //getchar();

            }

        }
        //cout<<"第"<<level<<"层用户节点剪枝验证结束！"<<endl;
        //cout<<"related_poi_node size in next iteration:"<<related_poi_nodes.size()<<endl;
        //cout<<"related_usr_node size in next iteration:"<<related_usr_node.size()<<endl;
        if(related_usr_node.size()==0) break;

    }//end while
    UPPER_PRINT("upper nodes cost")
    UPPER_CLEAR
    LEAF_PRINT("leaf nodes cost")
    LEAF_CLEAR
    cout<<"upper prune="<<upper_prune<<",upper unprune="<<upper_unprune<<endl;
    cout<<"leaf prune="<<leaf_prune<<",leaf unprune="<<leaf_unprune<<endl;
    cout<<"leaf lcl computation cnt:"<<leafLcl_cnt<<endl;
    cout<<"aggressive_prune="<<aggressive_prune<<endl;
    //getchar();
}


void  Group_Filter_Once_disk_error(vector<int> stores, vector<int> keywords, int Qk, float a, float alpha, set<VerifyEntry>& verification_User, map<int,set<VerifyEntry>>& verification_Map, map<int, priority_queue<ResultCurrent>>& candidate_usr_related_store, map<int, LCLResult>& candidate_usr_LCL) {
    vector<int> qKey = keywords;
    set<int> store_leafSet;
    map<int, set<int>> partition_checkinMap;

    set<int> store_nodes;
    map<int, vector<int>> leaf_store_list;//<leaf_id, <store_id....>>  leaf 中的store的记录
    map<int, set<int>> leaf_verify_list;//<leaf_id, <usr_id....>>  leaf 中所有store对应的凝练用户集合
    unordered_map<int,unordered_map<int,bool>> upMask;

    BatchResults result_Stores; //<store_id, <usr_id....>>
    int loc =0;

    //generate the query tree for candidate stores
    for(int store_id:stores){
        POI p = getPOIFromO2UOrgLeafData(store_id);
        vector<int> p_keywords = p.keywords;
        for(int term: p_keywords){
            vector<int> users_related = getTermUserInvList(term);
            if(users_related.size()<=uterm_SPARCITY_VALUE){//先把用户稀疏型关键词对应的用户筛选出来，优先验证
                for(int usr_id:users_related){
                    double rel = 0;
                    double inf = 0;
                    //if(upMask[usr_id][store_id]==true) continue;
                    if(upMask.count(usr_id)){
                        if(upMask[usr_id].count(store_id)){
                            continue;
                        }
                    }
                    User u = getUserFromO2UOrgLeafData(usr_id);
                    Score_Bound gsk_score_bound = getGSKScore_o2u_Upper_and_Lower(a, alpha, p,
                                                                                  u, rel,inf);
                    double gsk_score_Upper = gsk_score_bound.upper_score;
                    double gsk_score_Lower = gsk_score_bound.lower_score;


                    if (candidate_usr_related_store.count(usr_id)) {
                        candidate_usr_related_store[usr_id].push(ResultCurrent(store_id, gsk_score_Lower, gsk_score_Upper, rel,
                                                                               inf));  //记录relevance与influence 避免重算
                        if (candidate_usr_related_store[usr_id].size() > Qk) {
                            candidate_usr_related_store[usr_id].pop();
                        }

                    }
                    else{
                        priority_queue<ResultCurrent> u_related_store;
                        u_related_store.push(ResultCurrent(store_id, gsk_score_Lower, gsk_score_Upper, rel,
                                                           inf));  //记录relevance与influence 避免重算
                        candidate_usr_related_store[usr_id] = u_related_store;

                    }

                    upMask[usr_id][store_id] = true;
                    VerifyEntry ve = VerifyEntry(usr_id,-1,0);
                    verification_User.insert(ve);

                    if(usr_id==7209){
                        cout<<"find usr_id==7209!"<<endl;//getchar();
                    }

                }
            }
        }


        int leaf = Nodes[p.Nj].gtreepath.back();
        leaf_store_list[leaf].push_back(store_id);
        //加入兴趣点的叶节点列表
        store_leafSet.insert(leaf);
        //叶子节点层partition的check-in信息
        for(int u: p.check_ins){
            partition_checkinMap[leaf].insert(u);
        }

        store_nodes.insert(leaf);
        int current = leaf;
        GTree[current].stores.insert(store_id);
        int father;
        int i = 0;
        while (current != 0) {

            father = GTree[current].father;
            GTree[father].child_hasStores.insert(current);
            GTree[father].stores.insert(store_id);
            store_nodes.insert(father);
            current = father;
            i++;
        }
        //cout<<"层数"<<i<<",leaf="<<leaf<<endl;
    }
    cout<<"候选商店分布在"<<store_leafSet.size()<<"个不同叶节点中,GTree的"<<store_nodes.size()<<"个不同节点中"<<endl;

    TIME_TICK_START
    set<int> keyword_Set;
    for(int t: keywords){  //只剩下用户密集型关键词
        vector<int> users_related = getTermUserInvList(t);
        if(users_related.size()>uterm_SPARCITY_VALUE){
            keyword_Set.insert(t);
        }
    }
    cout<<"剩下的需要filter验证的keywords:"; printSetElements(keyword_Set);



    //提取每个查询关键词 在leafNodeInv表中对应信息

    int upper_prune=0; int upper_unprune=0;
    int leaf_prune=0; int leaf_unprune=0;
    int aggressive_prune =0;
    int user_access=0; int usr_prune =0;
    int may_prune = 0;
    //int topK_count = 0; //topk执行次数
    // 获取与查询关键词相关的用户叶子节点
    int t3=0; int root=0;


    set<int> totalUsrSet;
    int nodeAccess = 0;
    //查询处理过程中，关键词相关 用户节点列表
    vector<int> related_usr_node;
    //查询处理过程中，关键词相关 兴趣点节点列表
    map<int,vector<int>> related_poi_nodes;

    /*for(int term:keyword_Set)
        for(int usr: invListOfUser[term])
            totalUsrSet.insert(usr);

    cout<<"与关键词相关的用户个数 num="<<totalUsrSet.size();*/


    Filter_START


    // 获取与查询关键词相关的****用户***** 所在子图对应的上层节点
    vector<int> unodes = GTree[root].children; //get the child list of root whose text contains term_id
    for (int n : unodes) {
        related_usr_node.push_back(n);
        GTree[n].ossociate_queryNodes = store_leafSet; //先只是查询对象组叶节点
        /*printf("GTree[%d].ossociate_queryNodes\n:",n);
        printSetElements(GTree[n].ossociate_queryNodes);*/
    }
    // 获取与查询关键词相关的****兴趣点***** 所在子图对应的上层节点
    for(int keyword: keyword_Set){
        vector<int> child_entry_o = getObjectTermRelatedEntry(keyword,root);
        for(int child: child_entry_o){
            related_poi_nodes[keyword].push_back(child);
        }
    }


    set<int> related_children;
    set<int> related_usr;
    set<int> related_leaves;
    //verify the user node from upper to bottom
    TIME_TICK_END
    TIME_TICK_PRINT("初始化完成，runtime")


    //verify the user node from upper to bottom
    int level = 1;  int cnt_ = 0;
    int leafLcl_cnt = 0;

    while (true) {  //daxunhuanonce
        //related_children.clear();
        if (related_usr_node.capacity() == 0) break;
        set<int>().swap(related_children);

        for (int node:related_usr_node) {

            TreeNode user_node = getGIMTreeNodeData(node,OnlyU);

            set<int> inter_Key  = obtain_itersection_jins(user_node.userUKeySet, keyword_Set);
            if(!inter_Key.size()>0) continue;

            //该节点为用户上层节点
            if(!user_node.isleaf){

                UPPER_START
                int usrNum = 0;
                double global_max_Score = -1;
#ifdef TRACK
                cout<<"non-leaf node："<< node <<endl;
#endif
                int user_capacity = -1;
                for(int keyword: inter_Key){
                    int u_cnt = user_node.term_usr_countMap[keyword];
                    user_capacity = max(user_capacity,u_cnt);
                }
                bool cardinality_flag = false;
                bool aggressive_flag =false;
                set<int> usr_leaf_direct;
#ifdef  AGGRESSIVE_I
                cardinality_flag = (user_capacity<2*SPARCITY_VALUE);
                if(cardinality_flag){
                    for(int term: inter_Key){
                        vector<int> usr_leaf_descend = getUsrTermRelatedLeafNode(term,node);
                        for(int usr_leaf: usr_leaf_descend){
                            usr_leaf_direct.insert(usr_leaf);
                        }
                    }
                    if(usr_leaf_direct.size()<SPARCITY_VALUE)
                        aggressive_flag = true;
                }

#endif

                if (aggressive_flag) {   //用户稀疏节点  //aggressive prune //user_capacity<4

                    for(int usr_leaf: usr_leaf_direct){
                        related_children.insert(usr_leaf);
                        for(int store_leaf:GTree[node].ossociate_queryNodes){
                            GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                        }

                    }
                    aggressive_prune++;
                    //cout << "用户稀疏节点" << node <<",usrNum="<<usrNum<<endl;
#ifdef TRACK
                    cout << "用户稀疏节点" << node <<",usrNum="<<usrNum<<endl;
#endif

                } else {  //用户密集型节点，需进行剪枝判断
#ifdef TRACK
                    cout << "用户密集型节点，需进行剪枝判断"<<endl;
#endif
                    //global upper bound  computation
                    vector<double> max_score_list;

                    /*printf("GTree[node]: n%d 's ossociate_queryNodes:\n",node);
                    printSetElements(GTree[node].ossociate_queryNodes);
                    getchar();*/

                    for (int store_leaf: GTree[node].ossociate_queryNodes) {
                        //upper bound computation
                        double  max_Score =
                                getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, partition_checkinMap[store_leaf], store_leaf, node); //注意：这里用 InterKey, 与联合check-in计算的评分上界
                        //cout<<"max_Score="<<max_Score<<endl;

                        global_max_Score = max( max_Score, global_max_Score);

                        max_score_list.push_back(max_Score);

                    }
                    //obtaining lcl from cache or lcl computation on the fly
                    double Rt_U = 100;
                    for(int term: inter_Key){
                        LCLResult lcl;
                        map<int,LCLResult>::iterator iter;
                        iter = GTree[node].termLCLDetail.find(term);
                        if(iter != GTree[node].termLCLDetail.end()){ //先访问cache
                            lcl = iter->second;
                        }else{
                            lcl = getUpperNode_termLCL_Disk(node,term, related_poi_nodes,Qk,a,alpha,global_max_Score);
                            GTree[node].termLCLDetail[term] = lcl;
                        }
                        double current_score = lcl.topscore;
                        Rt_U = min(Rt_U, current_score);
                    }
                    //cout<<"Rt_U="<<Rt_U<<endl;
                    // group pruning judgement
                    if (global_max_Score > Rt_U) {  // 该子树下所有用户不能被剪枝，加入列表后续进一步验证
                        upper_unprune++;

                        int _idx = 0;
                        for (int store_leaf: GTree[node].ossociate_queryNodes) {
                            //upper bound computation
                            double max_score = max_score_list[_idx];
                            _idx++;

                            if(Rt_U>max_score)   //满足pruning 条件，跳过
                                continue;

                            for(int term: inter_Key){
                                vector<int> child_entry_u = getUsrTermRelatedEntry(term,node);
                                if(child_entry_u.size()>0){
                                    for (int child: child_entry_u) {

                                        //孩子节点继承父节点的lcl内容
                                        GTree[child].termParentLCLDetail[term] = GTree[node].termParentLCLDetail[term];
                                        //加入后续候选用户
                                        related_children.insert(child);
                                        //添加 用户组与 对象组的 关联
                                        GTree[child].ossociate_queryNodes.insert(store_leaf);

                                    }
                                }
                            }

                        }


                    }
                        //整棵用户子树被剪枝
                    else {
                        upper_prune++;
#ifdef PRUNELOG
                        cout << "节点" << node << "与所有 store都不相关， 其下所有用户可被剪枝" << endl;
#endif

                    }


                }//else用户密集型
                UPPER_PAUSE
            }

            else if (user_node.isleaf) {  //若该节点为用户叶子节点
                LEAF_START
#ifdef TRACK
                cout << "access usr_leaf" << node << endl;
#endif


                MultiLCLResult lclResult_multi;
                vector<int> Ukeys;
                for(int term: inter_Key)
                    Ukeys.push_back(term);

                //global upper bound  computation
                double global_max_Score = -1;
                vector<double> max_score_list;
                for (int store_leaf: GTree[node].ossociate_queryNodes) {
                    //upper bound computation
                    double  max_Score =
                            getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, partition_checkinMap[store_leaf], store_leaf, node); //注意：这里用 InterKey, 与联合check-in计算的评分上界

                    global_max_Score = max( max_Score, global_max_Score);

                    max_score_list.push_back(max_Score);

                }

                //lcl computation, obtaining the lcl-r_t

                map<set<int>,MultiLCLResult>::iterator iter;
                iter = GTree[node].cacheMultiLCLDetail.find(inter_Key);
                if (iter!=GTree[node].cacheMultiLCLDetail.end()){
                    lclResult_multi = iter->second;
                }

                else{  //无cache  算之


                    //lclResult = getLeafNode_TopkSDijkstra(node, term_id, Qk, a, alpha, max_Score);

/*#ifdef LasVegas
                    lclResult_multi = getLeafNode_LCL_Dijkstra_multi(node, Ukeys, Qk, a, alpha, global_max_Score);
#else
                    lclResult_multi = getLeafNode_LCL_bottom2up_multi(node, Ukeys, Qk, a, alpha, global_max_Score);
#endif*/
                    lclResult_multi = getLeafNode_LCL_Dijkstra_multi(node, Ukeys, Qk, a, alpha, global_max_Score);



                    GTree[node].cacheMultiLCLDetail[inter_Key] = lclResult_multi;

                    leafLcl_cnt++;
                }

                //获得 usr_lcl update_o_set
                set<int> update_o_set;
                for(LCLResult lcl: lclResult_multi){
                    priority_queue<LCLEntry> LCL = lcl.LCL;
                    while(LCL.size() > 0){
                        LCLEntry le = LCL.top();
                        LCL.pop();
                        int o = le.id;
                        update_o_set.insert(o);
                        //cout<<"o"<<o<<",";
                    }

                }

                //获得lcl-r_t的数值
                double Rt_U = -1.0;

                if(update_o_set.size()<Qk){
                    //cout<<"non-enough(empty) update o set"<<endl;
                    Rt_U = 0.0;
                }
                else{
                    for(LCLResult lclResult : lclResult_multi){
                        double _ss = lclResult.topscore;
                        if(-1.0==Rt_U)
                            Rt_U =_ss;
                        else
                            Rt_U = min(Rt_U,_ss);
                    }
                }

                if(global_max_Score > Rt_U){
                    leaf_unprune++;
                    // evaluate each patition
                    int  _idx = 0; set<int> L;
#ifdef LEAF_UNPRUNELOG_
                    if(true){ //global_max_Score<2*Rt_U
                        cout<<"usr_leaf"<<node<<",global_max_Score="<<global_max_Score<<", Rt_u="<<Rt_U<<endl;
                    }
#endif

                    for (int store_leaf: GTree[node].ossociate_queryNodes) {  //每个partition在max_score_list中都有对应的upper bound
                        //upper bound computation
                        double  max_Score = max_score_list[_idx];
                        _idx++;

                        if( max_Score > Rt_U)
                            L.insert(store_leaf);
                        else{
                            continue;   //该store_leaf 被剪枝
                        }

                    }

                    //获取该叶子节点下各个用户
                    set<int> usr_withinLeaf;
                    for (int term: inter_Key) {
                        vector<int> usr_list = getUsrTermRelatedEntry(term,node);
                        for (int usr_id: usr_list) {   //获取叶节点下每个用户
                            usr_withinLeaf.insert(usr_id);
                        }
                    }

                    //对该叶子节点下各个用户进行评估
                    for (int usr_id: usr_withinLeaf) {   //获取每个用户

                        double u_lcl_rt = 0.0;


                        priority_queue<ResultCurrent> u_related_store;

                        if(update_o_set.size()>0) {
                            //obtain the lcl list of a user from leaf lcl
                            priority_queue<ResultCurrent> u_LCL = get_usrLCL_update_Plus(
                                    usr_id, update_o_set, Qk, a, alpha, global_max_Score);  //获取该用户的lcl

                            if (u_LCL.size() < Qk) u_lcl_rt = 0.0;
                            else if (u_LCL.size() == Qk) u_lcl_rt = u_LCL.top().score;

                        }

                        Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);
                        for (int store_leaf: L) {
                            for (int store: GTree[store_leaf].stores) {  //获取每个store_leaf中的store
                                double rel = 0;
                                double inf = 0;
                                POI p = getPOIFromO2UOrgLeafData(store);
                                User u = getUserFromO2UOrgLeafData(usr_id);

                                if(upMask.count(usr_id)){
                                    if(upMask[usr_id].count(store)){
                                        continue;
                                    }
                                }
                                upMask[usr_id][store] = true;
#ifdef GIVENKEY
                                Score_Bound gsk_score_bound = getGSKScore_o2u_Upper_and_Lower_givenKey(a, alpha, p,
                                                                                                       u, keywords, rel,
                                                                                                       inf);
#else
                                Score_Bound gsk_score_bound = getGSKScore_o2u_Upper_and_Lower(a, alpha, p,
                                                                                              u, rel,inf);
#endif
                                double gsk_score_Upper = gsk_score_bound.upper_score;
                                double gsk_score_Lower = gsk_score_bound.lower_score;
                                //TopkQueryCurrentResult result = topkSDijkstra_verify_usr(transformToQ(Users[usr_id]),Qk,a, alpha,gsk_score_Upper);


                                if (u_lcl_rt < gsk_score_Upper) {
                                    //update lcl
                                    u_related_store.push(ResultCurrent(store, gsk_score_Lower, gsk_score_Upper, rel,
                                                                       inf));  //记录relevance与influence 避免重算
                                    if (u_related_store.size() > Qk) {
                                        u_related_store.pop();
                                        u_lcl_rt = max(u_lcl_rt, u_related_store.top().score);
                                    }
                                    //cout<<"u"<<usr_id<<", gsk_score_Lower="<<gsk_score_Lower<<",u_lcl="<<u_lcl<<",加入store="<<store<<endl;
                                    Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);

                                } else {
                                    usr_prune++;
                                    //continue;   //该用户不是反top-k结果
                                }

                            }
                        }
                        if (u_related_store.size() > 0) {
                            VerifyEntry ve = VerifyEntry(usr_id,-1,u_lcl_rt);
                            verification_User.insert(ve);
                            if (candidate_usr_related_store[usr_id].size() > 0) {  //之前已经有u_lcl的结果了，更新之
                                priority_queue<ResultCurrent> u_related_Query_pre = candidate_usr_related_store[usr_id];
                                while (!u_related_store.empty()) {
                                    ResultCurrent rc = u_related_store.top();
                                    u_related_store.pop();
                                    u_related_Query_pre.push(rc);
                                    if (u_related_Query_pre.size() > Qk)
                                        u_related_Query_pre.pop();
                                }
                            } else {  //若无，直接赋予之
                                priority_queue<ResultCurrent> u_related_Query = u_related_store;
                                candidate_usr_related_store[usr_id] = u_related_Query;
                            }


#ifdef PRINTINFO
                            cout << "发现候选用户 u" << usr_id << ", u_lcl=" << u_lcl_rt << endl;
#endif

                            //getchar();
                        } else {
                            //cout << "u" << usr_id << "用户被剪枝, as Lu为空,u_lcl=" << u_lcl << endl;
                        }
                        //getchar();
                    }


                }
                    //该用户叶子节点整体被剪枝
                else {
#ifdef LEAF_PRUNELOG_
                    cout<<"uleaf "<<node<<"与所有查询对象都无关！"<<endl;
#endif
                    //prune++;
                    leaf_prune++;
                }

                //cout<<"验证叶节点结束"<<node<<endl;
                LEAF_PAUSE
            }


        }// the first for loop
        level++;
        map<int, vector<int>> related_poi_nodes_next;
        map<int, vector<int>> related_poi_leaf_next;

        //upfold related poi nodes
        for(int term: keyword_Set){
            for(int pn: related_poi_nodes[term]){
                if(!GTree[pn].isleaf){
                    vector<int> child_entry_o = getObjectTermRelatedEntry(term,pn);
                    if(child_entry_o.size()!=0){
                        for(int child: child_entry_o){
                            related_poi_nodes_next[term].push_back(child);
                        }
                    }
                }
                else{
                    related_poi_leaf_next[term].push_back(pn);
                }
            }
        }

        // upfold the related usr nodes
        vector<int>().swap(related_usr_node);
        related_usr_node.clear();

        for (int cn:related_children){       //jiachildyici

            related_usr_node.push_back(cn);

        }
        if(related_usr_node.size() ==0){
            //cout<<"related_usr_node空了"<<endl;
            break;

        }
        //防止单路分支带来的重复验证
        if(related_children.size()==1){
            cnt_++;
            if(cnt_==3){
                //cout<<"直接跳到叶节点"<<endl;
                //getchar();
                related_usr_node.clear();
                for(int node: related_children){
                    if(GTree[node].isleaf){
                        related_usr_node.push_back(node);  //

                    }

                    else{
                        set<int> leaves;
                        for(int term: keyword_Set){
                            vector<int> usr_leaf = getUsrTermRelatedEntry(term,node);
                            for(int leaf: usr_leaf){
                                leaves.insert(leaf);
                            }
                        }
                        for(int usr_leaf: leaves){
                            related_usr_node.push_back(usr_leaf);
                            //加入关联关系
                            for (int store_leaf: GTree[node].ossociate_queryNodes) {
                                GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                            }
                        }

                    }


                }
                //cout<<"新加入"<<related_usr_node.size()<<endl;
                //getchar();

            }

        }
        //cout<<"第"<<level<<"层用户节点剪枝验证结束！"<<endl;
        //cout<<"related_poi_node size in next iteration:"<<related_poi_nodes.size()<<endl;
        //cout<<"related_usr_node size in next iteration:"<<related_usr_node.size()<<endl;
        if(related_usr_node.size()==0) break;

    }//end while
    UPPER_PRINT("upper nodes cost")
    UPPER_CLEAR
    LEAF_PRINT("leaf nodes cost")
    LEAF_CLEAR
    cout<<"upper prune="<<upper_prune<<",upper unprune="<<upper_unprune<<endl;
    cout<<"leaf prune="<<leaf_prune<<",leaf unprune="<<leaf_unprune<<endl;
    cout<<"leaf lcl computation cnt:"<<leafLcl_cnt<<endl;
    cout<<"aggressive_prune="<<aggressive_prune<<endl;
    //getchar();
}


//2. 优化的验证阶段 OPT
//内存操作
BRkGSKQResults Group_Verification_Optimized_memory(set<VerifyEntry> verification_User, map<int, priority_queue<ResultCurrent>> candidate_usr_related_store, vector<int> qKey, int Qk, float a, float alpha){
    //存放RkGSKQ结果的
    BRkGSKQResults result_Stores;
    set<VerifyEntry> candidate_usr = verification_User;//iter->second; //
    //cout<<"candidate_usr.size="<<candidate_usr.size()<<endl;
    map<int, vector<int>> usr_group;
    //1. 按maxscore，将用户排序(OK)；2. 按关键词将候选用户分组
    vector<CardinalityFisrt> user_list;
    map<int, vector<int>> term_candidateUserInv;
    for(VerifyEntry ve: candidate_usr){
        int usr_id = ve.u_id;
        double u_lcl_rt = ve.rk_current;
        int leaf_node = Nodes[Users[usr_id].Nj].gtreepath.back();
        usr_group[leaf_node].push_back(usr_id);
        priority_queue<ResultCurrent> relate_query_queue = candidate_usr_related_store[usr_id];
        int size= relate_query_queue.size();
        //cout<<"leaf_store_list[store_leaf].size="<<candidate_usr_related_store_leaf[usr_id].size()<<endl;
        double maximum_upper_score = 0; double maximum_lower_score = 0; double max_score = 0;
        priority_queue<ResultLargerFisrt> Lu;


        double rk_score = 0;
        if(Users[usr_id].topkScore_current>0)
            rk_score = Users[usr_id].topkScore_current;
        while(!relate_query_queue.empty()){
            ResultCurrent rc = relate_query_queue.top();
            relate_query_queue.pop();
            //先算出真实gsk评分
            double rel = rc.relevance; double inf = rc.influence; //shabi
            double gsk_score = getGSKScore_q2u(a,  alpha, POIs[rc.o_id], Users[usr_id], rel, inf);

            double score_upper = rc.score_upper;
            double _score = rc.score;
            if(gsk_score < rk_score)  //jins
                continue;
            //if(score_upper<rk_score)
            //break;

            ResultLargerFisrt rlf(rc.o_id,gsk_score,rel, inf);
            Lu.push(rlf);

            //max_score = max(max_score,gsk_score);
        }

        if(Lu.size()>0){
            CardinalityFisrt cf(usr_id,u_lcl_rt, Lu);
            //user_queue.push(cf);
            user_list.push_back(cf);
        }

        //构建候选用户的关键词倒排表
        for(int term: Users[usr_id].keywords){
            term_candidateUserInv[term].push_back(usr_id);
        }

    }

    // verifiying the TkGSKQ results for each users in candidate users, report query objects which havie user as top-k
    int useful =0; int usr_count =0;

    //可并行化

#ifdef Multicore
    omp_set_num_threads(4);
#pragma omp parallel for
#endif
    for(int i=0; i<user_list.size();i++){
        //topkSDijkstra_verify_usr(transformToQ(Users[0]), Qk, a, alpha, 1);
        CardinalityFisrt cf = user_list[i];
        usr_count++;
        int usr_id = cf.usr_id;
        double usr_lcl_rt = cf.current_rk;
        priority_queue<ResultLargerFisrt> Lu = cf.Lu;
        //cout<<"第"<<usr_count<<"个验证用户,关联"<<Lu.size()<<"家store,max_score="<<Lu.top().score<<endl;
        double Rk_u=0.0;

#ifdef TRACK
        cout<<"验证u"<<usr_id<<endl;

#endif

        if(Lu.top().score > Users[usr_id].topkScore_current){

//jordan
            double score_bound = usr_lcl_rt;
            //cout<<"u"<<usr_id<<", rt="<<usr_lcl_rt<<endl;


            Rk_u = candidate_user_verify_for_batch(usr_id,Qk,a, alpha,Lu.top().score,score_bound);

            //更新其他用户的当前topk score值
            //updateTopkScore_Current(a, Qk, current_elements,term_candidateUserInv);
            topK_count++;
            if(Lu.top().score > Rk_u){
                useful++;
                //对Lu中的各个query object进行评估
                while(!Lu.empty()){
                    ResultLargerFisrt rlf = Lu.top();
                    Lu.pop();
                    double rel = 0; double inf = 0;
                    double gsk_score = rlf.score; int query_object = rlf.o_id;
                    if(gsk_score > Rk_u){
                        ResultDetail rd(usr_id, 0, 0, 0, gsk_score, Rk_u);
                        result_Stores[query_object].push_back(rd);
                        //cout<<"确实为有效用户，gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
                    }
                    else break;
                }
            }
            else{
                //cout<<"topk score = "<<Rk_u<<",max_score="<<max_score<<endl;
            }
        }
        else{
            //cout<<"不用算！"<<endl;
        }

        //该候选用户已被验证
        Users[usr_id].verified = true;
    }
    //cout<<"useful ratio"<<useful*1.0/topK_count<<endl;
    //getchar();
    return  result_Stores;
}

//外存操作
BRkGSKQResults Group_Verification_Optimized_disk(set<VerifyEntry> verification_User, map<int, priority_queue<ResultCurrent>> candidate_usr_related_store, vector<int> qKey, int Qk, float a, float alpha){
    //存放RkGSKQ结果的
    BRkGSKQResults result_Stores;
    set<VerifyEntry> candidate_usr = verification_User;//iter->second; //
    //cout<<"candidate_usr.size="<<candidate_usr.size()<<endl;
    map<int, vector<int>> usr_group;
    //1. 按maxscore，将用户排序(OK)；2. 按关键词将候选用户分组
    vector<CardinalityFisrt> user_list;
    map<int, vector<int>> term_candidateUserInv;
    for(VerifyEntry ve: candidate_usr){
        int usr_id = ve.u_id;
        /* if(7209==usr_id){
             cout<<"verify find u7209！"<<endl;
         }*/
        double u_lcl_rt = ve.rk_current;
        User u = getUserFromO2UOrgLeafData(usr_id);
        int leaf_node = Nodes[u.Nj].gtreepath.back();
        usr_group[leaf_node].push_back(usr_id);
        priority_queue<ResultCurrent> relate_query_queue = candidate_usr_related_store[usr_id];
        int size= relate_query_queue.size();
        //cout<<"leaf_store_list[store_leaf].size="<<candidate_usr_related_store_leaf[usr_id].size()<<endl;
        double maximum_upper_score = 0; double maximum_lower_score = 0; double max_score = 0;
        priority_queue<ResultLargerFisrt> Lu;


        double rk_score = 0;
        if(Users[usr_id].topkScore_current>0)
            rk_score = Users[usr_id].topkScore_current;
        while(!relate_query_queue.empty()){
            ResultCurrent rc = relate_query_queue.top();
            relate_query_queue.pop();
            //先算出真实gsk评分
            double rel = rc.relevance; double inf = rc.influence;
            POI p = getPOIFromO2UOrgLeafData(rc.o_id);
            double gsk_score = getGSKScore_q2u(a,  alpha, p, u, rel, inf);

            double score_upper = rc.score_upper;
            double _score = rc.score;
            if(gsk_score < rk_score)  //jins
                continue;
            //if(score_upper<rk_score)
            //break;

            ResultLargerFisrt rlf(rc.o_id,gsk_score,rel, inf);
            Lu.push(rlf);

            //max_score = max(max_score,gsk_score);
        }

        if(Lu.size()>0){
            CardinalityFisrt cf(usr_id,u_lcl_rt, Lu);
            //user_queue.push(cf);
            user_list.push_back(cf);
        }

        //构建候选用户的关键词倒排表
        for(int term: u.keywords){
            term_candidateUserInv[term].push_back(usr_id);
        }

    }

    // verifiying the TkGSKQ results for each users in candidate users, report query objects which havie user as top-k
    int usr_count =0; int potential =0;

#ifdef  Verification_DEBUG
    cout<<"总共有"<<user_list.size()<<"个candidate user!"<<endl;
#endif

    //可并行化

#ifdef Multicore
    // omp_set_num_threads(4);
//#pragma omp parallel for
#endif
    for(int i=0; i<user_list.size();i++){
        //topkSDijkstra_verify_usr(transformToQ(Users[0]), Qk, a, alpha, 1);
        CardinalityFisrt cf = user_list[i];
        usr_count++;
        int usr_id = cf.usr_id;
        double usr_lcl_rt = cf.current_rk;
        priority_queue<ResultLargerFisrt> Lu = cf.Lu;
        //cout<<"第"<<usr_count<<"个验证用户,关联"<<Lu.size()<<"家store,max_score="<<Lu.top().score<<endl;
        double Rk_u=0.0;
#ifdef  Verification_DEBUG
        cout<<"尝试验证用户：u"<<usr_id<<endl;
#endif
/*#ifdef TRACK
        cout<<"发现用户：u"<<usr_id<<endl;
#endif*/
        if(Lu.top().score > Users[usr_id].topkScore_current){

//jordan
#ifdef TRACK
            cout<<"验证用户：u"<<usr_id<<endl;
#endif
            double score_bound = usr_lcl_rt;
            //cout<<"u"<<usr_id<<", rt="<<usr_lcl_rt<<endl;


            Rk_u = candidate_user_verify_for_batch(usr_id,Qk,a, alpha,Lu.top().score,score_bound);
            //Rk_u = candidate_user_verify_for_OGPM(usr_id,Qk,a, alpha,10000,0);


            //更新其他用户的当前topk score值
            //updateTopkScore_Current(a, Qk, current_elements,term_candidateUserInv);
            topK_count++;
            if(Lu.top().score > Rk_u){
#ifdef TRACK
                cout<<"发现潜在用户u"<<usr_id<<",Rk_u="<<Rk_u<<endl;

#endif

                potential++;
                //对Lu中的各个query object进行评估
                while(!Lu.empty()){
                    ResultLargerFisrt rlf = Lu.top();
                    Lu.pop();
                    double rel = 0; double inf = 0;
                    double gsk_score = rlf.score; int query_object = rlf.o_id;
                    if(gsk_score > Rk_u){
                        ResultDetail rd(usr_id, 0, 0, 0, gsk_score, Rk_u);
                        result_Stores[query_object].push_back(rd);
#ifdef TRACK
                        cout<<",对于o"<<query_object<<"确实为有效用户，gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
#endif
                    }
                    else break;
                }
            }
            else{
                //cout<<"topk score = "<<Rk_u<<",max_score="<<max_score<<endl;
            }
        }
        else{
            //cout<<"不用算！"<<endl;
        }

        //该候选用户已被验证
        Users[usr_id].verified = true;

#ifdef TRACK
        //cout<<"验证候选user"<<usr_id<<",Rk_u="<<Rk_u<<endl;

#endif
    }
    cout<<"real potential user cnt="<<potential<<", ratio="<<potential*1.0/topK_count<<endl;
    //getchar();
    return  result_Stores;
}




//put above two steps togather
BatchResults RkGSKQ_Batch_Group(vector<int> stores, vector<int> keywords, int Qk, float a, float alpha){
#ifndef DiskAccess
    cout<<"RkGSKQ_Batch_Group (memory) begin!"<<endl;
#else
    cout<<"RkGSKQ_Batch_Group (disk) begin, search potential custorms for "<<stores.size()<<" chain stores"<<endl;
#endif
    int block_num_pre = block_num;

    clock_t batch_start, batch_end;
    Filter_START;
    batch_start = clock();

    /*----------开始执行Group Filtering phase------------------*/
    set<VerifyEntry> verification_User;  //候选验证用户集
    map<int,set<VerifyEntry>> verification_Map;    //凝练结果表   <leaf_id, <user_id....>>
    map<int, priority_queue<ResultCurrent>> candidate_usr_related_store;   //候选用户与查询对象的关联表
    map<int, LCLResult> candidate_usr_LCL;     //候选用户的LCL表
    BatchResults batch_query_results;
    //组合过滤阶段
    Filter_START
#ifndef DiskAccess
    Group_Filter_Once_memory(stores, keywords, Qk, a, alpha, verification_User, verification_Map, candidate_usr_related_store, candidate_usr_LCL);
#else
    Group_Filter_Once_disk(stores, keywords, Qk, a, alpha, verification_User, verification_Map, candidate_usr_related_store, candidate_usr_LCL);

#endif
    Filter_END

    //开始凝练过程
    Refine_START
#ifndef DiskAccess
    batch_query_results = Group_Verification_Optimized_memory(verification_User, candidate_usr_related_store, keywords, Qk,a,alpha);
#else
    batch_query_results = Group_Verification_Optimized_disk(verification_User, candidate_usr_related_store, keywords, Qk,a,alpha);
#endif
    batch_end = clock();
    Refine_END
    //对每个候选usr进行分析并打印结果
    //FilterResults results =  resultsAnalysis(result_User,keywords,checkin_usr,loc, Qk,a,alpha);
    printBatchRkGSKQResults(stores, batch_query_results);
    cout<<"filter阶段总共找到"<<verification_User.size()<<"个候选用户"<<endl;
    Uc_size = verification_User.size();
    cout<<"verification user size: "<<verification_User.size()<<endl;
    int io_count = block_num - block_num_pre;
    //cout<<"i/o cost is "<<io_count<<endl;
    Filter_PRINT("filter lasts")
    Refine_PRINT("refinement lasts")

    cout<<"RkGSKQ_Batch_ by GIM-Tree runtime:"<<(double)(batch_end - batch_start)/CLOCKS_PER_SEC*1000<<" ms!"<<endl;

    return batch_query_results;

}





#define MAXINFRKGSKQ_1_1_QUERY_PLUS_H

#endif //MAXINFRKGSKQ_1_1_QUERY_PLUS_H