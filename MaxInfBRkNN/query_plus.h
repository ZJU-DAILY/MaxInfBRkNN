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

        //double distance = getDistance(POIs[poi_id].Nj, usr);//???????????????
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

        //if (score > Rk_usr && distance < 10000) {   // 10km????????????(???????????????????????????)
        if (gsk_score > Rk_usr) {
            size++;
            ResultDetail sd(usr.id,influence,distance,relevance,gsk_score,Rk_usr);
            results.push_back(sd);
            //cout<<"u"<<usr.id<<",";
        }
        cnt++;
    }
    cout <<"p"<<p.id<<",RkGSKQ ??????????????????=" << size << "?????????????????????"<<endl;
    printResults(results);
    cout<< "topk????????????"<<cnt<<endl;
    TIME_TICK_END
    TIME_TICK_PRINT("Naive_RkGSKQ running time:")
    return 0;

}

//more fast
double getGSKScore_o2u_phl(int a, double alpha, POI poi, User user){   //poi???user???GSK??????

    set<int> uKeyset; set<int> qKeyset;
    //for(int term : Users[usr_id].keywords) uKeyset.insert(term);
    //for(int term :qKeys)  qKeyset.insert(term);
    //double relevance = textRelevance(uKeyset,qKeyset);  //?????? u.keywords?????????????????????qKey)??????????????????????????????
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


    double distance = usrToPOIDistance_phl(user, poi);   //????????????????????? getDistance(vertex, vertex)
    double social_textual_score = alpha*influence + (1.0-alpha)*relevance;
    double gsk_score = social_textual_score / (1+distance);
    //if(poi.id== 57&& user.id ==625 ) cout<<"getGSKScore_o2u?????????u57,p625,distance="<<distance<<",relevance="<<relevance<<"influence="<<influence<<",score="<<gsk_score<<endl;
    return gsk_score;

}


double getSKScore_o2u_phl(int a, double alpha, POI poi, User user){   //poi???user???GSK??????

    set<int> uKeyset; set<int> qKeyset;
    //for(int term : Users[usr_id].keywords) uKeyset.insert(term);
    //for(int term :qKeys)  qKeyset.insert(term);
    //double relevance = textRelevance(uKeyset,qKeyset);  //?????? u.keywords?????????????????????qKey)??????????????????????????????
    vector<int> qKeys = poi.keywords;
    vector<int> check_usr = poiCheckInIDList[poi.id];

    double relevance = textRelevance(user.keywords,qKeys);
    if(relevance == 0){
        return  0.0;
    }
    double influence = 1;


    double distance = usrToPOIDistance_phl(user, poi);   //????????????????????? getDistance(vertex, vertex)
    double social_textual_score = alpha*influence + (1.0-alpha)*relevance;
    double gsk_score = social_textual_score / (1+distance);
    //if(poi.id== 57&& user.id ==625 ) cout<<"getGSKScore_o2u?????????u57,p625,distance="<<distance<<",relevance="<<relevance<<"influence="<<influence<<",score="<<gsk_score<<endl;
    return gsk_score;

}




//?????????????????????????????????
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
    //??????????????????????????? ???leafNodeInv??????????????????

    //????????????usr????????????????????????????????????
    int prune=0; int unprune=0; int user_access=0; int usr_prune =0;
    int leaf_prune=0; int leaf_unprune=0;
    int may_prune = 0;
    // ???????????????????????????????????????????????????
    set<int> related_usrLeafSet;
    int t3=0; int root=0;

    //set<int> verification_User;  //?????????
    set<VerifyEntry> verification_User;
    //map<int,map<int,set<VerifyEntry>>> verification_Map;  //???????????????   <term_id, leaf_id, <user_id....>>
    map<int,set<VerifyEntry>> verification_Map;    //???????????????   <term_id, <user_id....>>
    //???????????????????????????filter??????
    set<int> totalUsrSet;
    Filter_START
    for(int term_id:keywords){
        vector<int> related_usr_node;
        vector<int> related_poi_nodes;
        // ?????????????????????????????????****??????****?????????????????????????????????
        //cout<<"for w"<<term_id<<endl;
        set<int> nodes = usrNodeTermChild[root][term_id]; //get the child list of root whose text contains term_id
        for (int n: nodes) {
            //cout<<"node"<<n<<endl;
            related_usr_node.push_back(n);
            t3++;
        }
        //cout << "(usr)related_nodeSet size =" << related_usr_node.size() << ",usr node access=" << t3;
        // ?????????????????????????????????****?????????***** ?????????????????????????????????
        t3 = 0;
        set<int> nodes2 = poiNodeTermChild[root][term_id]; //get the child list of root whose text contains term_id
        for (int n: nodes2) {
            //cout<<"node"<<n<<endl;
            related_poi_nodes.push_back(n);
            t3++;
        }
        //cout << ",(poi)related_nodeSet size =" << related_poi_nodes.size() << ",poi node access=" << t3;

        // ??????????????? ??????????????????????????????????????????????????????????????????????????????????????? ??????????????????????????????????????????????????????
        int i = 0;
        set<int> related_poi_leafSet;
        set<int> related_usr_leafSet;
        map<int, set<int>> node_poiCnt;
        map<int, set<int>> node_occurence_leaf;
        //??????????????????????????????
        set<int> term_OLeaf = poiTerm_leafSet[term_id];
        for (int leaf: term_OLeaf) {
            related_poi_leafSet.insert(leaf);  //????????????poi_leaf
            for (int poi: leafPoiInv[term_id][leaf]) {
                //?????????leaf????????????poi(????????????
                node_poiCnt[leaf].insert(poi);
                //?????????leaf????????????????????????????????????poi(????????????
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
        //???????????????????????????
        map<int, set<int>> node_relateUserCnt; //  < leaf_id, <usr_id...> >
        set<int> term_ULeaf = usrTerm_leafSet[term_id];
        map<int, set<int>> node_usr_occurence_leaf;
        int totalUsr =0;
        for (int leaf: term_ULeaf) {
            related_usr_leafSet.insert(leaf);  //????????????usr_leaf
            for (int usr: leafUsrInv[term_id][leaf]) {
                //?????????leaf??????????????????(????????????
                node_relateUserCnt[leaf].insert(usr);
                totalUsrSet.insert(usr);
                //?????????leaf????????????????????????????????????poi(????????????
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

                if (GTree[node].isleaf) {  //?????????????????????????????????
                    //cout<<"usr_leaf"<<node<<endl;
                    //?????? MaxFGSK(q, U)
                    double min_Distance = getMinMaxDistance(loc, node).min;  //?????????????????????

                    double max_Score = getMaxScoreOfLeafOnce(loc, node, min_Distance, keyword_Set, checkin_usr_Set, a);

                    double Rt_U = 0;
                    LCLResult lclResult;
                    LCLResult lclResult2;

                    //???cache
                    if(GTree[node].termLCLDetail[term_id].topscore > 0){
                        lclResult = GTree[node].termLCLDetail[term_id];
                        //Rt_U = GTree[node].termLCLRt[term_id];
                    }
                        //???cache
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
                    //????????????????????????????????????
                    //cout<<"max_Score="<<max_Score<<",Rt_U="<<Rt_U<<endl;
                    if (max_Score > Rt_U) {  // check if leaf can be pruned
                        unprune++;
                        leaf_unprune++;

                        //?????????????????????????????????
                        if (leafUsrInv[term_id][node].size() != 0) {
                            for (int usr_id: leafUsrInv[term_id][node]) {
                                double gsk_score = getGSKScore_loc(a, alpha,loc, usr_id, qKey,checkin_usr);
                                //???Rt_U????????????????????????????????????
                                if (gsk_score > Rt_U) {
                                    double u_lcl = 0;
                                    u_lcl = get_usrLCL_light(usr_id,lclResult,node_poiCnt,Qk,a, alpha, gsk_score);
                                    //u_lcl = get_usrLCL_heavy(usr_id,lclResult,node_poiCnt,Qk,a,gsk_score);
                                    //cout<<"u"<<usr_id<<", gsk_score="<<gsk_score<<",u_lcl="<<u_lcl<<endl;
                                    if(u_lcl > gsk_score){
                                        //cout<<"?????????"<<usr_id<<"?????????, gsk_score="<<gsk_score<<",u_lcl="<<u_lcl<<endl;
                                        usr_prune++;
                                        continue;   //??????????????????top-k??????

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
                        //?????????????????????????????????
                    else {
                        //cout<<"leaf "<<node<<"????????????"<<endl;
                        prune++;
                        leaf_prune++;
                    }
                }
                    //??????????????????????????????
                else {
                    int usrNum = node_relateUserCnt[node].size();
                    //cout<<"???????????????????????????"<< usrNum <<endl;
                    if(usrNum < 5){
                        for (int usr_id: node_relateUserCnt[node]){ //??????????????????????????????????????????????????????
                            double gsk_score = getGSKScore_loc(a, alpha, loc, usr_id, qKey,checkin_usr);
                            VerifyEntry ve = VerifyEntry(usr_id,gsk_score,0.0);
                            verification_User.insert(ve);
                            verification_Map[term_id].insert(ve);
                        }
                        //???????????????????????????????????????????????????(??????)
                        /*for (int leaf_id: node_usr_occurence_leaf[node]){
                            related_children.insert(leaf_id);
                        }*/
                        prune++;
                        //cout<<"?????????????????? usr_node"<<node<<",prune"<<endl;   //??????????????????

                    }
                    else{
                        //cout<<"?????????????????????????????????????????????"<<endl;
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
                        if (max_Score > Rt_U) {  // ???????????????????????????????????????????????????????????????????????????
                            unprune++;
                            if (usrNodeTermChild[node][term_id].size() != 0) {
                                for (int child: usrNodeTermChild[node][term_id]){
                                    //??????????????????????????????lcl??????
                                    GTree[child].termParentLCLDetail[term_id] = lcl;
                                    //GTree[child].termParentFilter[term_id] = filter;
                                    related_children.insert(child);

                                }

                            }
                            //cout<<"usr_node"<<node<<",unprune"<<endl;
                            //cout<<"Rt_U="<<Rt_U<<endl;
                        }
                            //???????????????????????????
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
                    //cout<<"?????????????????????"<<endl;
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
                    //cout<<"?????????"<<related_usr_node.size()<<endl;
                    //getchar();

                }

            }

            //cout<<"???"<<level<<"????????????????????????????????????"<<endl;
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
    //??????????????????
    Refine_START
    //for(int v: verification_User){
    for(VerifyEntry v: verification_User){
        user_access++;
        int usr_id = v.u_id;
        double Rk_u;
        double gsk_score = v.score;//getGSKScore(a,loc, usr_id, qKey,checkin_usr);

        //continue;
        //????????????top-k??????????????????
        bool flag = true;
        if(Users[usr_id].topkScore_current > -1){  //???cache

            Rk_u = Users[usr_id].topkScore_current;
            //if(usr_id==800){
            //cout<<"cache u800="<< Rk_u<<endl;
            //}
            if(gsk_score < Rk_u){
                //cout<<"???????????????topk?????????prune u"<<usr_id<<", gsk_score="<<gsk_score<<",Rk_U="<<v.rk_current<<"Rk_u"<<Rk_u<<endl;
                continue;
            }  //????????????topk??????

            else{   //???????????????topk??????

                TopkQueryCurrentResult result = topkSDijkstra_verify_usr_memory(transformToQ(Users[usr_id]),Qk,a, alpha,gsk_score,0.0);
                topK_count++;
                Rk_u = result.topkScore;
                vector<Result> current_elements = result.topkElements;
                //updateTopk(a, usr_id, Qk, verification_Map, current_elements);


            }

        }
        else{   //???cache

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
            //usr_id???topk??????????????????,??????topk????????????U.LCL
            ResultDetail rd(usr_id,0,0,0, gsk_score,Rk_u);
            result_User.push_back(rd);
        }

    }
    Refine_END

    //???????????????usr???????????????????????????
    SingleResults results =  resultsAnalysis_loc(result_User,keywords,checkin_usr,loc, Qk,a,alpha);
    printResults(results.candidateUsr);
    cout<<"?????????????????? =" << user_access <<endl;
    cout<<"???????????????????????????????????????"<< totalUsrSet.size() <<endl;
    cout<<"topk?????????????????? =" << topK_count <<endl;
    cout<<"??????????????????= "<<results.candidateUsr.size()<<"????????????"<<endl;
    cout<<"g-tree node prune ="<<prune<<" ,unprune="<<unprune<<",usr_prune="<<usr_prune<<endl;
    double ratio = (1.0*prune) /(prune+unprune);
    double ratio2 = (1.0*leaf_prune) /(leaf_prune+leaf_unprune);
    cout<<"???????????????= "<<ratio<<"??????????????????????????????="<<ratio2<<endl;
    //cout<<"g-tree node access ="<<nodeAccess<<endl;
    TIME_TICK_END
    Filter_PRINT("filter lasts")
    Refine_PRINT("refinement lasts")
    TIME_TICK_PRINT("RkGSKQ_Single_Plus's ")

    return results;

}
//?????????????????????
//????????????

SingleResults RkGSKQ_Single_poi_Plus(POI poi, vector<int> keywords, vector<int> checkin_usr, int Qk, float a, float alpha){
    TIME_TICK_START

    cout<<"RkGSKQ_Single_poi_Plus !"<<poi.id<<endl;


    //???????????????????????????????????????????????????
    set<int> keyword_Set;
    for(int t: keywords)
        keyword_Set.insert(t);

    set<int> checkin_usr_Set;
    for(int u: checkin_usr)
        checkin_usr_Set.insert(u);
    vector<int> qKey = keywords;

    vector<int> Leaf_node_User,Leaf_node_Poi;

    vector<ResultDetail> result_User;
    //??????????????????????????? ???leafNodeInv??????????????????

    //????????????usr????????????????????????????????????
    int prune=0; int unprune=0; int user_access=0; int usr_prune =0;
    int leaf_prune=0; int leaf_unprune=0;
    int may_prune = 0;
    // ???????????????????????????????????????????????????
    set<int> related_usrLeafSet;
    int t3=0; int root=0;

    //set<int> verification_User;  //?????????
    set<VerifyEntry> verification_User;
    //map<int,map<int,set<VerifyEntry>>> verification_Map;  //???????????????   <term_id, leaf_id, <user_id....>>
    map<int,set<VerifyEntry>> verification_Map;    //???????????????   <term_id, <user_id....>>
    //???????????????????????????filter??????
    set<int> totalUsrSet;
    int nodeAccess = 0;
    vector<int> related_usr_node;
    map<int,vector<int>> related_poi_nodes;

    for(int term:keyword_Set)
        for(int usr: invListOfUser[term])
            totalUsrSet.insert(usr);


    Filter_START

    //for(int term_id:keywords){
    // ?????????????????????????????????****??????***** ?????????????????????????????????

    vector<int> unodes = GTree[root].children; //get the child list of root whose text contains term_id
    for (int n : unodes) {
        //cout<<"node"<<n<<endl;
        related_usr_node.push_back(n);
    }
    // ?????????????????????????????????****?????????***** ?????????????????????????????????
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
            //TEST_DURA_PRINT("?????????jins")
            //??????????????????????????????
            if (!GTree[node].isleaf) {

#ifdef TRACK
                cout<<"?????? node "<<node<<" ,covered query keywords:";
                //??????????????????, ??????????????????????????????????????????????????????
                printSetElements(inter_Key);
#endif

                int usrNum = 0; //node_relateUserCnt[node].size();
                //cout<<"?????????"<<"??????????????????"<< usrNum <<endl;
                if(true){
                    //cout<<"?????????????????????????????????????????????"<<endl;
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
                        if(iter != GTree[node].termLCLDetail.end() ){  //?????????cache
                            lcl = iter->second;
                        } else{ //?????????cache  ???????????????
                            lcl = getUpperNode_termLCL(node, term, related_poi_nodes, Qk, a,  alpha, max_Score);
                            GTree[node].termLCLDetail[term] = lcl;
#ifdef LCL_LOG
                            cout<<"?????? n"<<node<<" ???lcl?????????getUpperNode_termLCL???"<<endl;
                            //printLCLResultsInfo(lcl);

#endif

                        }

                        double current_score = lcl.topscore;
                        Rt_U = min(Rt_U,current_score);

                    }


                    if (max_Score > Rt_U) {  // ???????????????????????????????????????????????????????????????????????????
                        unprune++;
                        for(int term: inter_Key){
                            if (GTree[node].inverted_list_u[term].size() != 0) {
                                for (int child: GTree[node].inverted_list_u[term]){  //kobe
                                    //??????????????????????????????lcl??????
                                    GTree[child].termParentLCLDetail[term] = GTree[node].termLCLDetail[term];;
                                    related_children.insert(child);

                                }

                            }
                        }
#ifdef TRACK
                        if(Rt_U*a > max_Score){
                            //cout<<"??????"<<node<<" unpruned, max_score="<<max_Score<<", Rt_U="<<Rt_U<<endl;
                            //getchar();
                        }

#endif

                    }
                        //???????????????????????????
                    else {
                        prune++;
#ifdef PRUNELOG
                        cout<<"??????"<<node<<"??????????????????????????????"<<endl;
#endif

                    }
                }
            }
            else if (GTree[node].isleaf) {  //?????????????????????????????????

#ifdef TRACK

                cout<<"?????? usr_leaf "<<node<<", covered query keywords:";
                printSetElements(inter_Key);
#endif
                //?????? MaxFGSK(q, U)

                int loc = poi.Nj;
                double min_Distance = getMinMaxDistance(loc, node).min;  //?????????????????????


                double max_Score =
                        getMaxScoreOfLeaf_InterKey(loc, node, min_Distance, inter_Key, checkin_usr_Set, alpha,a);
#ifdef TRACK
                cout<<"max_Score="<<max_Score<<endl;
#endif

                //jordan
                MultiLCLResult lclResult_multi;
                vector<int> Ukeys;
                //??????????????????
                for(int key:inter_Key)
                    Ukeys.push_back(key);

                //???cache
                if (GTree[node].nodeLCL_multi.size() > 0){
                    lclResult_multi = GTree[node].nodeLCL_multi;
                    //Rt_U = GTree[node].termLCLRt[term_id];
                }
                    //???cache
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

                //?????? usr_lcl update_o_set
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

                //????????????????????????????????????
                //cout<<"????????????????????????????????????: max_Score="<<max_Score<<",Rt_U="<<Rt_U<<endl;
                if (max_Score > Rt_U) {  // check if leaf can be pruned
                    unprune++;
                    leaf_unprune++;
#ifdef TRACK
                    if(a*Rt_U > max_Score){
                        //cout<<"?????????"<<node<<" unpruned, max_score="<<max_Score<<", Rt_U="<<Rt_U<<endl;

                    }
#endif

                    //?????????????????????????????????
                    set<int> relate_usr_fromLeaf;
                    for(int term: Ukeys) {
                        for(int u : GTree[node].inverted_list_u[term]){
                            relate_usr_fromLeaf.insert(u);

                        }
                    }

                    for (int usr_id: relate_usr_fromLeaf) {

                        double gsk_score =-1;


#ifdef GIVENKEY
                        gsk_score = getGSKScore_o2u_givenKey(a,alpha, poi, Users[usr_id], qKey,checkin_usr);  //??????????????????????????????????????????????????????loc???check-in ??????????????????givenKeyword???

#else
                        gsk_score = getGSKScore_o2u(a, alpha, poi, Users[usr_id]);  //?????????????????????????????????????????? ?????????????????????????????? ????????????
#endif



                        //???Rt_U????????????????????????????????????
                        double u_lcl_rt =0;
                        if (gsk_score > Rt_U) {

                            if(update_o_set.size()>0) {
                                priority_queue<ResultCurrent> u_LCL = get_usrLCL_update_Plus(usr_id,
                                                                                             update_o_set, Qk, a, alpha, gsk_score);  //??????????????????lcl

                                if (u_LCL.size() < Qk) u_lcl_rt = 0.0;
                                else if (u_LCL.size() == Qk) u_lcl_rt = u_LCL.top().score;
                            }

                            if(u_lcl_rt > gsk_score){
                                //cout<<"?????????"<<usr_id<<"?????????, gsk_score="<<gsk_score<<",u_lcl_rt="<<u_lcl_rt<<endl;
                                usr_prune++;
                                continue;   //??????????????????top-k??????

                            }

                            Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current,u_lcl_rt);

                            VerifyEntry ve = VerifyEntry(usr_id,gsk_score,u_lcl_rt);
                            verification_User.insert(ve);
                            verification_Map[term_id].insert(ve);
#ifdef TRACK
                            cout<<"?????????????????? u"<<usr_id<<", gsk_score="<<gsk_score<<",u_lcl_rt="<<u_lcl_rt<<endl;
#endif
                        } else {
#ifdef PRUNELOG
                            cout<<"??????"<<usr_id<<" is pruned"<<endl;
#endif
                            //cout<<"gsk_Score="<<gsk_score<<endl;
                            //cout<<"Rt_U="<<Rt_U<<endl;

                        }

                    }

                }
                    //?????????????????????????????????
                else {
#ifdef PRUNELOG
                    cout<<"leaf "<< node<<" ????????????"<<endl;
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
                //cout<<"?????????????????????"<<endl;
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
                //cout<<"?????????"<<related_usr_node.size()<<endl;
                //getchar();

            }

        }
#ifdef TRACK
        cout<<"???"<<level<<"????????????????????????????????????"<<endl;
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
    //??????????????????
    Refine_START
    //for(int v: verification_User){
    for(VerifyEntry v: verification_User){

        user_access++;
        int usr_id = v.u_id;
        double Rk_u=0;
        double gsk_score = v.score;//getGSKScore(a,loc, usr_id, qKey,checkin_usr);
        double u_lcl_rt = v.rk_current;
        //continue;
        //????????????top-k??????????????????
        bool flag = true;

        if(Users[usr_id].topkScore_current > -1){  //???cache

            Rk_u = Users[usr_id].topkScore_current;

            if(gsk_score < Rk_u){

#ifdef VERIFICATION_LOG
                cout<<"???????????????topk?????????prune u"<<usr_id<<", gsk_score="<<gsk_score<<",Rk_U(cache)="<<Rk_u<<endl;

#endif
                continue;
            }  //????????????topk??????

            else{   //???????????????topk??????
                if(Users[usr_id].topkScore_Final > -1)
                    Rk_u = Users[usr_id].topkScore_Final;
                else{   //????????????topk??????

                    Rk_u = candidate_user_verify_single(usr_id,Qk,a, alpha,gsk_score);

                    topK_count++;

                }
            }

        }
        else{   //???cache

            Rk_u = candidate_user_verify_single(usr_id,Qk,a, alpha,gsk_score);
            topK_count++;
            //updateTopk(a, usr_id, Qk, verification_Map, current_elements);

        }

#ifdef VERIFICATION_LOG
        cout<< "u"<<v.u_id<<":Rk_u="<<Rk_u<<"gsk(u,o)="<<gsk_score<<endl;
#endif

        if(gsk_score > Rk_u || gsk_score == Rk_u){
            //usr_id???topk??????????????????,??????topk????????????U.LCL
            ResultDetail rd(usr_id,0,0,0, gsk_score,Rk_u);
            result_User.push_back(rd);
        }

    }
    Refine_END

    //???????????????usr???????????????????????????
    SingleResults results =  resultsAnalysis_poi(result_User,poi, Qk,a,alpha);

    cout<<"?????????????????? =" << user_access <<endl;
    cout<<"???????????????????????????????????????"<< totalUsrSet.size() <<endl;
    cout<<"topk?????????????????? =" << topK_count <<endl;
    cout<<"lcl for leaf ????????????"<<lcl_leafCount<<endl;
    cout<<"p"<<poi.id<<"??????????????????= "<<results.candidateUsr.size()<<"????????????"<<endl;
    printResults(results.candidateUsr);
    /*cout<<"g-tree node prune ="<<prune<<" ,unprune="<<unprune<<",usr_prune="<<usr_prune<<endl;
    double ratio = (1.0*prune) /(prune+unprune);
    double ratio2 = (1.0*leaf_prune) /(leaf_prune+leaf_unprune);
    cout<<"???????????????= "<<ratio<<"??????????????????????????????="<<ratio2<<endl;
    //cout<<"g-tree node access ="<<nodeAccess<<endl;*/
    TIME_TICK_END
    Filter_PRINT("filter lasts")
    Refine_PRINT("refinement lasts")
    TIME_TICK_PRINT("RkGSKQ_Single_Plus end!")
    cout<<"??????IO: number of block accessed="<<block_num<<endl;

    return results;

}


//????????????
SingleResults RkGSKQ_Single_poi_Disk(POI poi, vector<int> keywords, vector<int> checkin_usr, int Qk, float a, float alpha){
    TIME_TICK_START

    cout<<"---------Running RkGSKQ_by_Group_on_GIM-Tree!---------"<<endl;

    printPOIInfo(poi);

    //???????????????????????????????????????????????????
    set<int> keyword_Set;
    for(int t: keywords)
        keyword_Set.insert(t);

    set<int> checkin_usr_Set;
    for(int u: checkin_usr)
        checkin_usr_Set.insert(u);
    vector<int> qKey = keywords;

    vector<int> Leaf_node_User,Leaf_node_Poi;

    vector<ResultDetail> result_User;
    //??????????????????????????? ???leafNodeInv??????????????????

    //????????????usr????????????????????????????????????
    int upper_prune; int upper_unprune=0;
    int leaf_prune=0; int leaf_unprune=0;
    int user_access=0; int usr_prune =0;

    int may_prune = 0;
    // ???????????????????????????????????????????????????
    set<int> related_usrLeafSet;
    int t3=0; int root=0;

    //set<int> verification_User;  //?????????
    set<VerifyEntry> verification_User;
    //map<int,map<int,set<VerifyEntry>>> verification_Map;  //???????????????   <term_id, leaf_id, <user_id....>>
    map<int,set<VerifyEntry>> verification_Map;    //???????????????   <term_id, <user_id....>>
    //???????????????????????????filter??????
    set<int> totalUsrSet;
    int nodeAccess = 0;
    vector<int> related_usr_node;
    map<int,vector<int>> related_poi_nodes;


    Filter_START

    // ?????????????????????????????????****??????***** ?????????????????????????????????

    //vector<int> unodes = GTree[root].children; //get the child list of root whose text contains term_id
    //printElements(unodes);
    vector<int> unodes = getGIMTreeNodeData(root,OnlyU).children;
    //printElements(unodes);


    for (int n : unodes) {
        //cout<<"node"<<n<<endl;
        related_usr_node.push_back(n);
    }
    // ?????????????????????????????????****?????????***** ?????????????????????????????????
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
                cout<<"--------------?????? leaf "<<node<<endl;
            }
            else
                cout<<"--------------??????  node"<<node<<endl;
#endif
            TreeNode tn = getGIMTreeNodeData(node, OnlyU);
            //TEST_DURA_PRINT("?????????jins")
            set<int> inter_Key  = obtain_itersection_jins(tn.userUKeySet, keyword_Set);
#ifdef TRACK
            cout<<"?????? node"<<node<<"???????????? ,covered query keywords:";
            //??????????????????, ??????????????????????????????????????????????????????
            printSetElements(inter_Key);
#endif


            //??????????????????????????????
            if (!tn.isleaf) {

                UPPER_START

                int usrNum = 0; //node_relateUserCnt[node].size();
                //cout<<"?????????"<<"??????????????????"<< usrNum <<endl;
                if(true){
                    //cout<<"?????????????????????????????????????????????"<<endl;
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
                        if(iter != GTree[node].termLCLDetail.end() ){  //?????????cache
                            lcl = iter->second;
                        } else{ //?????????cache  ???????????????

                            lcl = getUpperNode_termLCL_Disk(node, term, related_poi_nodes, Qk, a,  alpha, max_Score);
#ifdef TRACK
                            //cout<<"??????t"<<term<<", ?????? n"<<node<<" ???lcl?????????getUpperNode_termLCL_Disk???, lcl_rt="<<lcl.topscore<<endl;
                            //printLCLResultsInfo(lcl);

#endif

                            //GTree[node].termLCLDetail[term] = lcl;  ????????????cache

                        }

                        double current_score = lcl.topscore;
                        Rt_U = min(Rt_U,current_score);

                    }

#ifdef TRACK
                    cout<<"???????????????max_Score="<<max_Score<<",Rt_U="<<Rt_U<<endl;
#endif
                    if (max_Score > Rt_U) {  // ???????????????????????????????????????????????????????????????????????????
                        upper_unprune++;
                        for(int term: inter_Key){
                            //int _size = GTree[node].inverted_list_u[term].size();
                            vector<int> child_entry = getUsrTermRelatedEntry(term,node);
                            //int _size2 = child_entry.size();
                            if (child_entry.size() != 0) {
                                for (int child: child_entry){  //kobe
                                    //??????????????????????????????lcl??????
                                    GTree[child].termParentLCLDetail[term] = GTree[node].termLCLDetail[term];;
                                    related_children.insert(child);

                                }

                            }
                        }
#ifdef TRACK
                        if(Rt_U*a > max_Score){
                            //cout<<"??????"<<node<<" unpruned, max_score="<<max_Score<<", Rt_U="<<Rt_U<<endl;
                            //getchar();
                        }

#endif

                    }
                        //???????????????????????????
                    else {
                        upper_prune++;
#ifdef TRACK
                        cout<<"??????"<<node<<"??????????????????????????????"<<endl;
#endif

                    }
                }

                UPPER_PAUSE
            }
            else if (GTree[node].isleaf) {  //?????????????????????????????????

#ifdef TRACK
                cout<<"?????? usr_leaf "<<node<<", covered query keywords:";
                printSetElements(inter_Key);
#endif
                //?????? MaxFGSK(q, U)

                int loc = poi.Nj;
                double min_Distance = getMinMaxDistance(loc, node).min;  //?????????????????????


                double max_Score =
                        getMaxScoreOfLeaf_InterKey(loc, node, min_Distance, inter_Key, checkin_usr_Set, alpha,a);

#ifdef  TRACK
                cout<<"max_Score= "<<max_Score<<endl;
#endif

                //jordan
                MultiLCLResult lclResult_multi;
                vector<int> Ukeys;
                //??????????????????
                for(int key:inter_Key)
                    Ukeys.push_back(key);

                //???cache
                //map<set<int>,MultiLCLResult> :: iterator iter2;
                //iter2 = GTree[node].cacheMultiLCLDetail.find(inter_Key);
                if(GTree[node].cacheMultiLCLDetail.count(inter_Key)>0){  //?????????cache
                    lclResult_multi = GTree[node].cacheMultiLCLDetail[inter_Key];
                }
                    //???cache
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

                        //GTree[node].cacheMultiLCLDetail[inter_Key] = lclResult_multi;  ??????cache

                        lcl_leaf_count++;
                    }


                }

                //?????? usr_lcl update_o_set
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

                //????????????????????????????????????
                //cout<<"????????????????????????????????????: max_Score="<<max_Score<<",Rt_U="<<Rt_U<<endl;
                if (max_Score > Rt_U) {  // check if leaf can be pruned  //
                    //unprune++;
                    leaf_unprune++;
#ifdef TRACK
                    //if(a*Rt_U > max_Score){
                        //cout<<"?????????"<<node<<" unpruned, max_score="<<max_Score<<", Rt_U="<<Rt_U<<endl;

                    //}
#endif

                    //?????????????????????????????????,??????????????????????????????????????????
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
                        //gsk_score = getGSKScore_o2u_givenKey(a,alpha, poi, Users[usr_id], qKey,checkin_usr);  //??????????????????????????????????????????????????????loc???check-in ??????????????????givenKeyword???
                            gsk_score = getGSKScore_o2u_givenKey(a,alpha, poi, user, qKey,checkin_usr);


#else //joke phl
                        gsk_score = getGSKScore_o2u_phl(a, alpha, poi, user);  //?????????????????????????????????????????? ?????????????????????????????? ????????????
#endif

                        //cout<<"u"<<usr_id<<"gsk_score="<<gsk_score<<endl;
                        //???Rt_U????????????????????????????????????
                        double u_lcl_rt =0;
                        if (true) {  //gsk_score > Rt_U

                            if(update_o_set.size()>0) {
                                priority_queue<ResultCurrent> u_LCL = get_usrLCL_update_Plus(usr_id,
                                                                                             update_o_set, Qk, a, alpha, gsk_score);  //??????????????????lcl

                                if (u_LCL.size() < Qk) u_lcl_rt = 0.0;
                                else if (u_LCL.size() == Qk) u_lcl_rt = u_LCL.top().score;
                            }

                            if(u_lcl_rt > gsk_score){
                                //cout<<"?????????"<<usr_id<<"?????????, gsk_score="<<gsk_score<<",u_lcl_rt="<<u_lcl_rt<<endl;
                                usr_prune++;
                                continue;   //??????????????????top-k??????

                            }

                            Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current,u_lcl_rt);

                            VerifyEntry ve = VerifyEntry(usr_id,gsk_score,u_lcl_rt);
                            verification_User.insert(ve);
                            verification_Map[term_id].insert(ve);
#ifdef TRACK
                            cout<<"?????????????????? u"<<usr_id<<", gsk_score="<<gsk_score<<",u_lcl_rt="<<u_lcl_rt<<endl;
#endif
                        } else {
#ifdef PRUNELOG
                            cout<<"??????"<<usr_id<<" is pruned"<<endl;
#endif
                            //cout<<"gsk_Score="<<gsk_score<<endl;
                            //cout<<"Rt_U="<<Rt_U<<endl;

                        }
                    }

                }
                    //?????????????????????????????????
                else {
#ifdef PRUNELOG
                    cout<<"leaf "<< node<<" ????????????"<<endl;
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
                    //cout<<"??????node"<<pn<<"???t"<<term<<" ???????????????"<<endl;
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
        cout<<"??????related_poi_nodes"<<endl;
#endif
        related_poi_nodes.clear();

        for(int term:keyword_Set){
            //cout<<"for term"<<term<<endl;
            for(int cn: related_poi_nodes_next[term])
                related_poi_nodes[term].push_back(cn);
            for (int lf:related_poi_leaf_next[term])
                related_poi_nodes[term].push_back(lf);
        }
        //cout<<"???????????????"<<endl;

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
                //cout<<"?????????????????????"<<endl;
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

                //cout<<"?????????"<<related_usr_node.size()<<"????????????"<<endl;
                //printElements(related_usr_node);
                //getchar();

            }

        }
#ifdef TRACK
        cout<<"???"<<level<<"????????????????????????????????????"<<endl;
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
    //??????????????????yanz
    Refine_START
    //for(int v: verification_User){
    for(VerifyEntry v: verification_User){

        user_access++;
        int usr_id = v.u_id;
        double Rk_u=0;
        double gsk_score = v.score;//getGSKScore(a,loc, usr_id, qKey,checkin_usr);
        double u_lcl_rt = v.rk_current;
        //continue;
        //????????????top-k??????????????????
        bool flag = true;

        if(Users[usr_id].topkScore_current > -1){  //???cache

            Rk_u = Users[usr_id].topkScore_current;

            if(gsk_score < Rk_u){

#ifdef VERIFICATION_LOG
                cout<<"???????????????topk?????????prune u"<<usr_id<<", gsk_score="<<gsk_score<<",Rk_U(cache)="<<Rk_u<<endl;

#endif
                continue;
            }  //????????????topk??????

            else{   //???????????????topk??????
                if(Users[usr_id].topkScore_Final > -1)
                    Rk_u = Users[usr_id].topkScore_Final;
                else{   //????????????topk??????

                    Rk_u = candidate_user_verify_single(usr_id,Qk,a, alpha,gsk_score);

                    topK_count++;

                }
            }

        }
        else{   //???cache

            Rk_u = candidate_user_verify_single(usr_id,Qk,a, alpha,gsk_score);
            topK_count++;
            //updateTopk(a, usr_id, Qk, verification_Map, current_elements);

        }

#ifdef VERIFICATION_LOG
        cout<< "u"<<v.u_id<<":Rk_u="<<Rk_u<<"gsk(u,o)="<<gsk_score<<endl;
#endif

        if(gsk_score > Rk_u || gsk_score == Rk_u){
            //usr_id???topk??????????????????,??????topk????????????U.LCL
            ResultDetail rd(usr_id,0,0,0, gsk_score,Rk_u);
            result_User.push_back(rd);
        }

    }
    Refine_END

    //???????????????usr???????????????????????????
    SingleResults results =  resultsAnalysis_poi(result_User,poi, Qk,a,alpha);
    TIME_TICK_END
    cout<<"?????????????????? =" << user_access <<endl;
    cout<<"?????????lcl computation?????????"<<lcl_leaf_count<<endl;
    //cout<<"topk?????????????????? =" << topK_count <<endl;
    cout<<"p"<<poi.id<<"??????????????????= "<<results.candidateUsr.size()<<"????????????"<<endl;
    printResults(results.candidateUsr);
    //cout<<"????????????:  prune="<<upper_prune<<" , unprune="<<upper_unprune<<endl;
    //cout<<"???????????? prune="<<leaf_prune<<" , unprune="<<leaf_unprune<<endl;
    /*cout<<"g-tree node prune ="<<prune<<" ,unprune="<<unprune<<",usr_prune="<<usr_prune<<endl;
    double ratio = (1.0*prune) /(prune+unprune);
    double ratio2 = (1.0*leaf_prune) /(leaf_prune+leaf_unprune);
    cout<<"???????????????= "<<ratio<<"??????????????????????????????="<<ratio2<<endl;
    //cout<<"g-tree node access ="<<nodeAccess<<endl;*/

    //UPPER_PRINT("runtime for upper node:") UPPER_CLEAR
    Filter_PRINT("filter lasts")
    Refine_PRINT("refinement lasts")
    TIME_TICK_PRINT("RkGSKQ_Single_Plus end!")

    //cout<<"??????IO: number of block accessed="<<block_num<<endl;

    return results;

}





//???????????????Filtering????????????NVD??????lcl, verification????????????GIM-Tree topk
void RkGSKQ_NVD_Naive_hash(int poi_id, int K, float a, float alpha) {

    //double direction

    string phl_fileName = getRoadInputPath(PHL);//"../exp/indexes/LV.phl";

    //char phl_idxFileName[255] ; sprintf(phl_idxFileName, "../exp/indexes/%s.phl", road_map);

    phl.LoadLabel(phl_fileName.c_str());


    priority_queue<Result> resultFinal;

    //??????poi?????????
    POI poi = getPOIFromO2UOrgLeafData(poi_id);
    printPOIInfo(poi);

    //?????????poi??????????????????user????????????????????????
    set<int> poiRelated_usr;
    //??????????????????
    vector<int> candidate_usr;
    vector<double> candidate_gsk;
    vector<double> candidate_rk_current;
    //??????????????? ???top-k???????????????
    vector<int> usr_results;
    for(int term: poi.keywords){
        vector<int> term_related_user = getTermUserInvList(term);
        for(int usr: term_related_user)
            poiRelated_usr.insert(usr);
    }
    printf("?????????poiRelated_usr??? %d???...\n", poiRelated_usr.size());
    //???poi??????????????????user??????????????????
    clock_t startTime, endTime;

    clock_t filter_startTime, filter_endTime;
    clock_t verification_startTime, verification_endTime;

    filter_startTime = clock();
    for(int usr_id: poiRelated_usr){
        //usr_id = 346;
        User user = getUserFromO2UOrgLeafData(usr_id);
        //cout<<"??????";
        //printUsrInfo(user);
        startTime = clock();//????????????
        //double dist = getDistanceLower_Oracle_PHL(user,poi);
        //double gsk_score = getGSKScore_o2u_distComputed(a,alpha,poi,user,dist);
        double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,user);
        endTime = clock();//????????????
        ////cout << "getGSKScore_o2u runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;
        cout<<"gsk_score="<<gsk_score;
        startTime = clock();//????????????
        double u_lcl_rk = getUserLCL_NVD(user,K,a,alpha);
        endTime = clock();//????????????
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
    cout<<"??????????????????="<<candidate_usr.size()<<"???????????????"<<endl;
    printElements(candidate_usr);


    verification_startTime = clock();
    int candidate_size = candidate_usr.size();
    for(int i=0;i<candidate_size;i++){
        int usr_id = candidate_usr[i];
        //cout<<"??????????????????"<<usr_id<<endl;
        double gsk_score = candidate_gsk[i];
        double score_bound = candidate_rk_current[i];
        double Rk_u = 0;

        priority_queue<Result> results = TkGSKQ_bottom2up_verify_disk(usr_id, K, a, alpha, 1000,0.0); //gsk_score, score_bound);
        if (results.size() == K) {
            Rk_u = results.top().score;
            //cout<<"gsk_score="<<gsk_score<<", Rk_u="<<Rk_u<<endl;
        }
        else{ //??????resultFinal????????????????????????????????????score_bound??????
            //cout<<"?????????"<<endl;
            Rk_u = score_bound;

        }

        //????????????
        if(gsk_score >= Rk_u){
            usr_results.push_back(usr_id);
        }

        cout<<"??????????????????"<<usr_id<<",gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;

    }
    verification_endTime = clock();
    cout<<"??????????????? ???top-k?????????????????????="<<usr_results.size()<<"???????????????"<<endl;
    printElements(usr_results);
    cout<<"Filter ?????? ?????????"<< (double)(filter_endTime - filter_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;
    cout<<"Verification ?????? ?????????"<< (double)(verification_endTime - verification_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;



}

TopkQueryCurrentResult TkGSKQ_NVD_vertexOnly(User user, int K, int a, double alpha, bool* poiMark, bool* poiADJMark){
    clock_t startTime, endTime, social_endTime;
    clock_t io_start, io_end; double io_time=0;
    clock_t textual_start, textual_end; double textual_time=0;

    priority_queue<GlobalEntry_Pseudo> Queue;
    priority_queue<Result> resultFinal;  //??????topk??????
    float score_rk = -1;
    vector<POI> polled_object_lists;
    vector<float> polled_distance_lists;

    int vi = user.Ni; int vj = user.Nj;
    //map<int,bool> poiIsVist;

    startTime = clock();//????????????
    //?????????user?????????
    int user_id = user.id;
    vector<int> friends = friendshipMap[user_id];
    //????????????friends ??????check-in??????poi
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
    //??????score
    double optimal_social_textual_score = alpha*(1.0)+(1-alpha)*optimal_text_score;


    ////??????K-SPIN????????????topk????????????
    //1.?????????user.keyword????????????queue
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
        //??????term???????????????
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
        else{ //????????????,???????????????POI??????
            termIsFrequent[term] = true;
            int poi_id = getNNPOI_By_Vertex_Keyword(user.Ni,term);  //???????????????Ni?????????????????????term???poi
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
                //cout<<"H"<<term<<"??????????????????p"<<poi_id<<endl;

            }else{

                if(true) {  //poiADJMark[poi_id] == false
                    POI poi1 = getPOIFromO2UOrgLeafData(poi_id);
                    float dist = usrToPOIDistance_phl(user,poi1);
                    nodeDistType_POI tmp1(poi_id, dist, poi1.keywords);
                    Hi.push(tmp1);
                    poiADJMark[poi_id] = true;
                    //cout<<"H"<<term<<"??????????????????p"<<poi_id<<endl;
                }
                if(true) { //poiADJMark[poi_id2] == false
                    POI poi2 = getPOIFromO2UOrgLeafData(poi_id2);
                    float dist = usrToPOIDistance_phl(user,poi2);
                    nodeDistType_POI tmp2(poi_id2, dist, poi2.keywords);
                    Hi.push(tmp2);
                    poiADJMark[poi_id2] = true;
                    //cout<<"H"<<term<<"??????????????????p"<<poi_id2<<endl;
                }

            }

        }
        H_list.push_back(Hi);
    }

    //1. ???????????????Queue????????????Hi?????????????????????????????????
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

    //2. ???????????????Queue????????????????????????Hi????????????
    while(Queue.size()>0 && Queue.top().score>score_rk){
        GlobalEntry_Pseudo entry = Queue.top();
        Queue.pop(); //global ??????
        int term_n = entry.term_id;
        int term_th = termid2idx[term_n];
        //cout<<"????????????Hi???t"<<term_n<<", score="<<entry.score<<endl;
        ////??????????????????Hi
        //priority_queue<nodeDistType_POI> Hn = H_list[term_th];
        ////Hi?????????????????????????????????????????????1.??????????????????NVD ??????????????????Hn???2.??????Hn???????????????????????????Hn????????????Queue???
        nodeDistType_POI entry_p  = H_list[term_th].top();
        int poi_current = entry_p.poi_id;

        vector<int> poi_current_keywords = entry_p.keywords;
        float dist_top_Hn = entry_p.dist;
        float score_Hn = entry.score;
        if(H_list[term_th].size()==0) continue;
        //cout<<"H"<<term_n<<"??????????????????p"<<entry_p.poi_id<<"??????"<<endl;
        //Hn.pop();
        H_list[term_th].pop();

        if(termIsFrequent[term_n]){  //??????????????????
            ////1.???????????????poi_current??????NVD_adj ??????Hn
            //cout<<"1.????????????p"<<entry_p.poi_id<<"???NVD_adj ??????Hn"<<endl;
            //cout<<"?????????H"<<term_n<<"???????????????"<<endl;

#ifdef TRACKNVDTOPK
            if(entry_p.poi_id==2633){
                cout<<"while loop find o2633"<<endl;
            }

#endif
            //vector<int> poi_ADJ = getPOIAdj_NVD_By_Keyword(entry_p.poi_id,term_n);
            vector<int> poi_ADJ = getPOIAdj_NVD_By_Keyword_vertexOnly(entry_p.poi_id,term_n);

            for(int _poi_id : poi_ADJ){
                //cout<<",nvd neighbor-p"<<_poi_id;
                if(poiADJMark[_poi_id] == true) continue;  //??????(Hn?????????????????????)
                POI poi = getPOIFromO2UOrgLeafData(_poi_id); int id = poi.id;

                //float p_dist2 = usrToPOIDistance_phl(getUserFromO2UOrgLeafData(67),getPOIFromO2UOrgLeafData(3413));

                float p_dist = usrToPOIDistance_phl(user,poi);
                nodeDistType_POI tmp_adj(id, p_dist, poi.keywords);
                H_list[term_th].push(tmp_adj);
                poiADJMark[_poi_id] = true;
                //cout<<",????????????(nvd_adj)p"<<id;

            }
        }

        //cout<<endl;
        ////2.??????Hn?????????????????????????????????????????????
        //cout<<"2.??????Hn?????????????????????????????????????????????"<<endl;
        //double pseudo_score_update = tfIDF_term(term_n);
        float dist_current = H_list[term_th].top().dist;
        //float score_updated = pseudo_score_update / dist;
        float score_updated = optimal_social_textual_score / (1+dist_current);
        GlobalEntry_Pseudo entry_Q(term_n,dist_current,score_updated);
        Queue.push(entry_Q);
        //cout<<"H"<<term_n<<"?????????????????????"<<score_updated<<endl;
        current_score_global = Queue.top().score;

        ////3.???poi_current??????????????????????????????????????????resultFinal???
        //cout<<"???????????????poi??????p"<<entry_p.poi_id<<endl;
        if(poiMark[poi_current] == true) continue;
        poiMark[poi_current] = true;
        userRelated_poiSet.insert(poi_current);

        //???????????? resultFinal
        double simT = textRelevance(user.keywords, poi_current_keywords);
        double simD = entry_p.dist;
        double social_textual_score = alpha*1 + (1-alpha)*simT;

        double gsk_score = social_textual_score / (1+simD);
        Result tmpRlt(poi_current,-1,-1,-1,simD,gsk_score,poi_current_keywords);
        if(resultFinal.size()<K){
            resultFinal.push(tmpRlt);
            //cout<<"resultFinal????????? p"<<poi_current<<endl;
        }
        else{  //??????current topk??????
            Result _tmpNode = resultFinal.top();
            double tmpScore = _tmpNode.score;
            if(tmpScore<gsk_score){
                resultFinal.pop();
                resultFinal.push(tmpRlt);
                float score_rk_update = resultFinal.top().score;
                score_rk = score_rk_update;
                //cout<<"resultFinal????????? p"<<poi_current<<endl;
                //cout<<"score_rk?????????"<<score_rk<<endl;
            }

        }


        float best_score = Queue.top().score;
        //cout<<"best_score="<<best_score<<",score_rk="<<score_rk<<endl;
        //??????????????????
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

    endTime = clock();//????????????


    TopkQueryCurrentResult  topkCurrentResult(resultFinal.top().score,resultFinal);
    //topkCurrentResult.printResults(K);

    //cout << "runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000 << "ms????????????" << endl;
    //cout << "social runtime is: " << (double)(social_endTime - startTime) / CLOCKS_PER_SEC * 1000 << "ms" << endl;
    /*cout << "io time is: " << io_time / CLOCKS_PER_SEC * 1000000 << "us" << endl;
    cout << "textual time is: " << textual_time / CLOCKS_PER_SEC * 1000000 << "us" << endl;
    cout<<"social_textualRelated_poiSet size="<<userRelated_poiSet.size()<<", count="<<count<<",io_count="<<io_count<<endl;
    cout<<"u"<<user.id<<"'s topk score="<<score_rk<<", resultFinal size = "<<resultFinal.size()<<",????????????"<<endl;*/

    return  topkCurrentResult;
    //printTopkResultsQueue(resultFinal);
}




//Filtering???verification???????????????NVD,??????Filtering????????????????????????????????????user??????????????????????????????user???lcl
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
    //??????poi?????????
    POI poi = getPOIFromO2UOrgLeafData(poi_id);
    printPOIInfo(poi);

    //?????????poi??????????????????user????????????????????????
    set<int> poiRelated_usr;
    //??????????????????
    vector<User> candidate_usr;
    vector<int> candidate_id;
    vector<double> candidate_gsk;
    vector<double> candidate_rk_current;
    //??????????????? ???top-k???????????????
    vector<ResultDetail> usr_results;
    vector<int> usrID_results;
    for(int term: poi.keywords){
        vector<int> term_related_user = getTermUserInvList(term);
        for(int usr: term_related_user)
            poiRelated_usr.insert(usr);
    }
    printf("?????????poiRelated_usr??? %d???...\n", poiRelated_usr.size());
    //???poi??????????????????user??????????????????

    clock_t filter_startTime, filter_endTime;
    clock_t verification_startTime, verification_endTime;

    //startTime = clock();
    filter_startTime = clock();
    for(int usr_id: poiRelated_usr){
        //usr_id = 346;
        User user = getUserFromO2UOrgLeafData(usr_id);
#ifdef TRACK
        cout<<"??????poi?????????????????????";
        printUsrInfo(user);
#endif

        startTime = clock();//????????????
        double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,user);
        endTime = clock();//????????????
        ////cout << "getGSKScore_o2u runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;
#ifdef TRACK
        cout<<"gsk_score="<<gsk_score;
#endif
        startTime = clock();//????????????
        double u_lcl_rk = getUserLCL_NVD(user,K,a,alpha);
        endTime = clock();//????????????
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
    cout<<"??????????????????="<<candidate_usr.size()<<"???????????????"<<endl;
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

        //????????????
        if(gsk_score >= Rk_u){
            //usr_results.push_back(u_id);
            ResultDetail ur(u_id,-1,-1,-1,gsk_score,Rk_u);
            usr_results.push_back(ur);
            usrID_results.push_back(u_id);
        }
#ifdef TRACK
        cout<<"??????????????????"<<u_id<<",gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
#endif


    }
    verification_endTime = clock();
    algEndTime = clock();


    cout<<"??????????????? ???top-k?????????????????????="<<usr_results.size()<<"???????????????"<<endl;
    //printElements(usrID_results);
    printPotentialUsers(usrID_results);
    cout<<"Filter ?????? ?????????"<< (double)(filter_endTime - filter_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;
    cout<<"Verification ?????? ?????????"<< (double)(verification_endTime - verification_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;
    cout << "---------------------RkGSKQ_NVD_userPure ?????????" << (double)(algEndTime - algBeginTime) / CLOCKS_PER_SEC * 1000 << "ms-----------" << endl;
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


    //??????????????????
    vector<User> candidate_usr;
    vector<int> candidate_id;
    vector<double> candidate_gsk;
    vector<double> candidate_rk_current;



    //user???????????????????????????
    bool* userMark = new bool[UserID_MaxKey];
    for(int i=0;i<UserID_MaxKey;i++)
        userMark[i] = false;

    int node_size = GTree.size();
    bool* nodeMark = new bool[node_size];
    for(int j=0;j<node_size;j++)
        nodeMark[j] = false;


    startTime = clock();
    filter_startTime = clock();

    //??????????????????
    priority_queue<CheckEntry> Queue;
    //??????poi?????????
    POI poi = getPOIFromO2UOrgLeafData(poi_id);
    printPOIInfo(poi);
    int p_Ni = poi.Ni; int p_Nj = poi.Nj;
    vector<int> qo_keys = poi.keywords;
    int qo_leaf = Nodes[p_Ni].gtreepath.back();
    TreeNode tnData = getGIMTreeNodeData(qo_leaf,OnlyU);
#ifdef TRACK
    cout<<"o_leaf???leaf"<<qo_leaf<<endl;
    cout<<"o_leaf???user?????????????????????"<<endl;
    printSetElements(tnData.userUKeySet);
#endif


    //// ?????? poi ??? o_leaf??????????????????border?????????
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

    ////?????????poi ????????????????????? user
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
            //???????????????????????????
            //cout<<"?????? ???????????????????????????u"<<usr_id<<endl;
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
    cout<<"????????????????????? user????????????="<<user_textual_count<<endl;
    cout<<"??????????????????????????? user?????????????????? candidate user size="<<candidate_id.size()<<"???????????????"<<endl;
    printElements(candidate_id);
#endif
    //getchar();

    ////???q_o????????????????????? ??????rare)????????????????????????(none-rare)????????????????????? ??????,?????????????????????keyword( |inv(term)|<K )??? user
    set<int> rareKeys;  set<int> nonRareKeys;
    set<int> rareKeywordRelated_usrSet;
    for(int term: qo_keys){
        int inv_olist_size = getTermOjectInvListSize(term);
        if(inv_olist_size>K){
            vector<int> usr_termRelated = getTermUserInvList(term);
            int inv_ulist_size = usr_termRelated.size();
            if(inv_ulist_size> uterm_SPARCITY_VALUE){
                nonRareKeys.insert(term);
                //printf("term %d ????????????????????????u_frequency=%d\n", term,inv_ulist_size);
            }
            else{  //??????????????????????????????user???????????????????????????????????????????????????
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
    ////??????????????????????????????user ????????????user???????????????????????????
    for(int usr_id:rareKeywordRelated_usrSet){
        if(userMark[usr_id] == true) continue;
        userMark[usr_id] = true;
        User u = getUserFromO2UOrgLeafData(usr_id);
        //int usr_id = u.id;
        double relevance = 0;
        relevance = textRelevance(u.keywords, poi.keywords);
        //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
        if(relevance==0) continue;
        //???????????????????????????
        //cout<<"?????? ???????????????????????????u"<<usr_id<<endl;
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
    cout<<"??????????????????????????? user ????????????="<<rareKeywordRelated_usrSet.size()<<endl;
    cout<<"????????????????????????????????? user?????????????????? candidate user size="<<candidate_id.size()<<"???????????????"<<endl;
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
    printf("????????????poi???????????????????????? %d?????? ???GTree?????????user occurant list??????...", poiRelated_usr.size());
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
    cout<<"?????????"<<(double)(occur_endTime-occur_startTime)/CLOCKS_PER_SEC*1000<<"ms!"<<endl;


    ////1. q_o ??????leaf ????????????????????? user?????????????????????
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
        cout<<"q_o ??????leaf"<<qo_leaf<<"????????????????????? user?????????????????????"<<endl;
    } else{
        cout<<"q_o ??????leaf "<<qo_leaf<<"????????????????????????user"<<endl;
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
    //????????????level?????????Gtree Node
    int node_highest = qo_leaf;
    int root = 0;

    ////2.???qo_leaf?????????????????????,????????????????????????Gtree??????qo_leaf???????????????????????????????????????????????????????????????????????????????????? user,user_leaf, user_node

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
        if(Queue.empty()){ //????????????????????????????????? ???????????????????????????
#ifdef TRACK_FILTER
            cout<<"??????????????????????????????,"<<"?????????????????????node"<<node_highest<<"??? ???????????????????????????????????????"<<endl;
#endif
            //Extend_Range(node_highest,itm,Queue,nonRareKeys,poi, K);
            Extend_Range(node_highest,itm,Queue,nonRareKeys,poi, K,a,alpha);
            if(Queue.empty()) break;  //????????????????????????

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
            //?????????????????????aggresive prune?????????discard
            if(nodeMark[node_id]==true)
                continue;
            //???????????????user??????????????? ?????????aggresive prune
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
                cout<<"----DEQUEUE: border entry: b"<<border_id<<" in leaf"<<node_id<<" ?????????----"<<endl;
            else
                cout<<"----DEQUEUE: border entry: b"<<border_id<<" in node"<<node_id<<" ?????????----"<<endl;
#endif

            set<int> Keys;


            //Keys = getBorder_SB_NVD(topEntry,K,a, alpha);
            Keys = getBorder_SB_NVD_Cluster(topEntry,K,a, alpha,nodes_uCLUSTERList);

            if(Keys.size() ==0) {
#ifdef TRACK_FILTER
                cout<<"?????????????????? ****prune???"<<endl;
#endif
                continue;
            }
#ifdef TRACK_FILTER
            cout<<"remained (user) Keywords"; printSetElements(Keys);
            cout<<"?????????**??????**???prune???"<<endl;
#endif

            //?????????????????????????????????????????????user
            if(GTree[node_id].isleaf){
#ifdef TRACK
                cout<<"b"<<border_id<<", ??????leaf"<<node_id<<endl;
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
                //???????????????????????????user
                for(int usr_id: user_within){
                    if(userMark[usr_id] == true) continue;
                    User user = getUserFromO2UOrgLeafData(usr_id);
                    double dist = usrToPOIDistance_phl(user, poi);
                    set<int> u_keySet;
                    CheckEntry entry(usr_id,node_id,dist,user);
                    Queue.push(entry);
#ifdef TRACK
                    cout<<"???main?????????user entry:"; entry.printRlt();
                    cout<<"Queue.size="<<Queue.size()<<endl;
#endif
                }
            }
            else{
                //??????????????????????????????????????????????????????(line 22-25)
#ifdef TRACK
                cout<<"b"<<border_id<<", ??????node"<<node_id<<endl;
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
                cout<<"???????????????????????????"; printSetElements(child_set);
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
                        cout<<"???main?????????border entry:"; entry.printRlt();
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


    ////?????????????????????user, ????????????TkGSKQ??????

#ifdef TRACK
    cout<<"??????????????????="<<candidate_usr.size()<<"???????????????"<<endl;
    printElements(candidate_id);
#else
    cout<<"??????????????????="<<candidate_usr.size()<<endl;
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
            else{ //??????resultFinal????????????????????????????????????score_bound??????
                Rk_u = score_bound;
            }

        }

        //????????????
        if(gsk_score >= Rk_u){
            RkGSKQ_results.push_back(u_id);
        }
#ifdef TRACK
        cout<<"??????????????????"<<u_id<<",gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
#endif

    }
    verification_endTime = clock();
    endTime = clock();
    cout<<"Filter ?????? ?????????"<< (double)(filter_endTime - filter_startTime) / CLOCKS_PER_SEC*1000  << "ms, ?????????while loop?????? "
        << (double)(filter_endTime - while_startTime) / CLOCKS_PER_SEC*1000 <<"ms)" << endl;
    cout<<"Verification ?????? ?????????"<< (double)(verification_endTime - verification_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;
    cout<<"??????????????? ???top-k?????????????????????="<<RkGSKQ_results.size()<<"???????????????"<<endl;
    printPotentialUsers(RkGSKQ_results);

    //cout<<"user_textual_count="<<user_textual_count<<",user_textual_useful_count="<<user_textual_useful_count<<endl;
    cout<<"---------RkGSKQ_bottom2up_by_NVD *p"<<poi.id<<"* ?????????"<< (double)(endTime - startTime) / CLOCKS_PER_SEC*1000  << "ms-----------" << endl;
}


SingleResults RkGSKQ_bottom2up_by_NVD(POI& poi, int K, float a, float alpha) {  //gangcai

    cout<<"------------Running Single RkGSKQ bottom2up based on IG-tree&NVD!------------"<<endl;
    clock_t startTime, endTime;
    clock_t filter_startTime,filter_endTime;
    clock_t verification_startTime, verification_endTime;

    vector<int> RkGSKQ_results;


    //??????????????????
    vector<User> candidate_usr;
    vector<int> candidate_id;
    vector<double> candidate_gsk;
    vector<double> candidate_rk_current;

    //user???????????????????????????
    bool userMark[UserID_MaxKey];
    memset(userMark,false, sizeof(userMark));


    startTime = clock();

    //??????????????????
    priority_queue<CheckEntry> Queue;
    //??????poi?????????
    int poi_id = poi.id;
    printPOIInfo(poi);
    int p_Ni = poi.Ni; int p_Nj = poi.Nj;
    vector<int> qo_keys = poi.keywords;
    int qo_leaf = Nodes[p_Ni].gtreepath.back();
    TreeNode tnData = getGIMTreeNodeData(qo_leaf,OnlyU);
#ifdef TRACK
    cout<<"o_leaf???leaf"<<qo_leaf<<endl;
    cout<<"o_leaf???user?????????????????????"<<endl;
    printSetElements(tnData.userUKeySet);
#endif


    //// ?????? poi ??? o_leaf??????????????????border?????????
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

    ////?????????poi ????????????????????? user
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
            //???????????????????????????
            //cout<<"?????? ???????????????????????????u"<<usr_id<<endl;
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
    cout<<"????????????????????? user????????????="<<user_textual_count<<endl;
    cout<<"??????????????????????????? user?????????????????? candidate user size="<<candidate_id.size()<<"???????????????"<<endl;
    printElements(candidate_id);
#endif
    //getchar();

    ////???q_o????????????????????? ??????rare)????????????????????????(none-rare)????????????????????? ??????,?????????????????????keyword( |inv(term)|<K )??? user
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
    ////??????????????????????????????user ????????????user???????????????????????????
    for(int usr_id:rareKeywordRelated_usrSet){
        if(userMark[usr_id] == true) continue;
        userMark[usr_id] = true;
        User u = getUserFromO2UOrgLeafData(usr_id);
        //int usr_id = u.id;
        double relevance = 0;
        relevance = textRelevance(u.keywords, poi.keywords);
        //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
        if(relevance==0) continue;
        //???????????????????????????
        //cout<<"?????? ???????????????????????????u"<<usr_id<<endl;
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
    cout<<"??????????????????????????? user ????????????="<<rareKeywordRelated_usrSet.size()<<endl;
    cout<<"????????????????????????????????? user?????????????????? candidate user size="<<candidate_id.size()<<"???????????????"<<endl;
    printElements(candidate_id);
#endif
    //getchar();

    ////1. q_o ??????leaf ????????????????????? user?????????????????????
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
        cout<<"q_o ??????leaf"<<qo_leaf<<"????????????????????? user?????????????????????"<<endl;
    } else{
        cout<<"q_o ??????leaf "<<qo_leaf<<"????????????????????????user"<<endl;
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
    //????????????level?????????Gtree Node
    int node_highest = qo_leaf;
    int root = 0;

    ////2.???qo_leaf?????????????????????,????????????????????????Gtree??????qo_leaf???????????????????????????????????????????????????????????????????????????????????? user,user_leaf, user_node

    map<int,map<int,bool>> borderEntryVist;  //border_id, node_id
    map<int,bool> userEntryVist;
    filter_startTime = clock();
    while(Queue.size()>0 ||node_highest!= root){
        if(Queue.empty()){ //????????????????????????????????? ???????????????????????????
#ifdef TRACK
            cout<<"??????????????????????????????,"<<"?????????????????????node"<<node_highest<<"??? ???????????????????????????????????????"<<endl;
#endif
            //Extend_Range(node_highest,itm,Queue,nonRareKeys,poi, K);
            Extend_Range(node_highest,itm,Queue,nonRareKeys,poi, K,a,alpha);
            if(Queue.empty()) break;  //????????????????????????

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
                cout<<"----DEQUEUE: border entry: b"<<border_id<<" in leaf"<<node_id<<" ?????????----"<<endl;
            else
                cout<<"----DEQUEUE: border entry: b"<<border_id<<" in node"<<node_id<<" ?????????----"<<endl;
#endif

            set<int> Keys;


            //CheckEntry border_entry,int K, int a, double alpha
            Keys = getBorder_SB_NVD(topEntry,K,a, alpha);
            //Keys = getBorder_SB_NVD(topEntry,poi,nonRareKeys,K,a, alpha);

            //cout<<"remained (user) Keywords"; printSetElements(Keys);
            if(Keys.size() ==0) {
#ifdef TRACK
                cout<<"?????????????????? ****prune???"<<endl;
#endif
                continue;
            }
#ifdef TRACK
            cout<<"?????????**??????**???prune???"<<endl;
#endif

            //?????????????????????????????????????????????user
            if(GTree[node_id].isleaf){
#ifdef TRACK
                cout<<"b"<<border_id<<", ??????leaf"<<node_id<<endl;
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
                //???????????????????????????user
                for(int usr_id: user_within){
                    if(userMark[usr_id] == true) continue;
                    User user = getUserFromO2UOrgLeafData(usr_id);
                    double dist = usrToPOIDistance_phl(user, poi);
                    set<int> u_keySet;
                    CheckEntry entry(usr_id,node_id,dist,user);
                    Queue.push(entry);
#ifdef TRACK
                    cout<<"???main?????????user entry:"; entry.printRlt();
                    cout<<"Queue.size="<<Queue.size()<<endl;
#endif
                }
            }
            else{
                //??????????????????????????????????????????????????????(line 22-25)
#ifdef TRACK
                cout<<"b"<<border_id<<", ??????node"<<node_id<<endl;
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
                cout<<"???????????????????????????"; printSetElements(child_set);
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
                        cout<<"???main?????????border entry:"; entry.printRlt();
                        cout<<"Queue.size="<<Queue.size()<<endl;
#endif

                    }
                }

            }
        }

    }
    filter_endTime = clock();


    ////?????????????????????user, ????????????TkGSKQ??????

#ifdef TRACK
    cout<<"??????????????????="<<candidate_usr.size()<<"???????????????"<<endl;
    printElements(candidate_id);
#endif

    verification_startTime = clock();
    int candidate_size = candidate_usr.size();
    vector<ResultDetail> potential_customerDetails;

////yanzsingle
    for(int i=0;i<candidate_size;i++){
        User u = candidate_usr[i];
        int u_id = u.id;
        //cout<<"?????????????????????"<<u_id<<endl;
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



        //????????????
        if(gsk_score >= Rk_u){
            RkGSKQ_results.push_back(u_id);
            //??????detail??????
            ResultDetail rd(u_id,0,0,0, gsk_score,Rk_u);
            potential_customerDetails.push_back(rd);

        }
#ifdef TRACK
        cout<<"??????????????????"<<u_id<<",gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
#endif
        //cout<<"??????????????????"<<u_id<<",gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;

    }
    verification_endTime = clock();
    endTime = clock();
    SingleResults signle_results =  resultsAnalysis_poi(potential_customerDetails,poi, K,a,alpha);

    cout<<"??????????????? ???top-k?????????????????????="<<RkGSKQ_results.size()<<"???????????????"<<endl;
    printPotentialUsers(RkGSKQ_results);
    cout<<"Filter ?????? ?????????"<< (double)(filter_endTime - filter_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;
    cout<<"Verification ?????? ?????????"<< (double)(verification_endTime - verification_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;

    //cout<<"user_textual_count="<<user_textual_count<<",user_textual_useful_count="<<user_textual_useful_count<<endl;
    cout<<"---------RkGSKQ_bottom2up_by_NVD ?????????"<< (double)(endTime - startTime) / CLOCKS_PER_SEC*1000  << "ms-----------" << endl;
    return signle_results;
}


//??????NVD???Batch Reverse????????????



BatchResults RkGSKQ_Batch_NVD_RealData(vector<int> stores, int K, float a, float alpha) {  //gangcai

    cout<<"------------Running Batch RkGSKQ using NVD&IG-Tree!------------"<<endl;
    //?????? ???????????????????????????u22557 for poi13448, user?????????
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


    //??????????????????
    set<BatchVerifyEntry> candidate_User;  //candidate users
    map<int, priority_queue<ResultCurrent>> candidate_usr_related_store;   // user_id, {related poi1,2,....}



    //border???????????????????????????
    int node_size = GTree.size();
    int inner_size = node_size*VertexNum;  // func(bEntry_key)= v_id+v_num*node_id + v_num*node_num*poi_th
    int bEntry_MaxKey =  inner_size * stores.size();
    unordered_map<int,bool> borderMark;


    //user??????????????????????????????
    int upEntry_MaxKey =  UserID_MaxKey * stores.size();
    unordered_map<int,bool> upEntryMark;
    unordered_map<int,bool> userMark;
    //??????gtree node??? poi ??????????????????
    int np_maxKey = GTree.size()*stores.size();
    bool nodepMark[np_maxKey];
    memset(nodepMark,false, sizeof(nodepMark));


    startTime = clock();
    filter_startTime = clock();
    initial_startTime = clock();

    //??????????????????

    ////?????????gtree node??? uOCCURList?????? uCLUSTERList
    unordered_map<int,UOCCURList> pnodes_uOCCURListMap;
    unordered_map<int,UCLUSTERList> pnodes_uCLUSTERListMap;



    priority_queue<BatchCheckEntry> Queue;
    priority_queue<POIHighestNode> T_max;


    int root = 0;
    //// 1. ????????????poi??????, ??????????????? poi ??? o_leaf??????????????????border?????????
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


        //?????????????????????poi????????? nodes_uOCCURList ??? nodes_uCLUSTERList
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
        cout<<"o_leaf???leaf"<<qo_leaf<<endl;
        cout<<"o_leaf???user?????????????????????"<<endl;
        printSetElements(tnData.userUKeySet);
#endif

        //// ?????? ?????? poi ??? o_leaf??????????????????border?????????
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

        ////?????????poi ????????????????????? user
        vector<int> checkin_users = poiCheckInIDList[poi_id];
        int user_textual_count = 0; int user_textual_useful_count = 0;
        for(int check_usr_id: checkin_users){
            vector<int> friends = friendshipMap[check_usr_id];
            for(int friend_id: friends){
                if(friend_id>=UserID_MaxKey) continue;
                int up_key = friend_id + UserID_MaxKey* poi_th;   //???????????????poi?????????0,usrnum-1)

                if(upEntryMark[up_key] == true) continue;
                upEntryMark[up_key] = true;


                User u = getUserFromO2UOrgLeafData(friend_id);
                int usr_id = u.id;
                double relevance = 0;
                relevance = textRelevance(u.keywords, poi.keywords);
                //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
                if(relevance==0) continue;
                //???????????????????????????
#ifdef TRACKNVDBATCH
                cout<<"?????? ????????????????????????1`???u"<<usr_id<<endl;
#endif
                user_textual_count ++;
                double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,u);
                double u_lcl_rk = getUserLCL_NVD(u,K,a,alpha);
                if(gsk_score > u_lcl_rk){

                    //?????? user-poi associating list
                    candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                    if(candidate_usr_related_store[usr_id].size()>K){
                        candidate_usr_related_store[usr_id].pop();
                    }
                    //??????candidate user set
                    BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_lcl_rk,u);
                    candidate_User.insert(ve);

                    user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                    cout<<"social_textual find candidate ????????? u"<<usr_id<<endl;
#endif
                }

            }
        }
#ifdef TRACKNVDBATCH
        cout<<"????????????????????? user????????????="<<user_textual_count<<endl;
    //cout<<"??????????????????????????? user?????????????????? candidate user size="<<candidate_id.size()<<"???????????????"<<endl;
    //printElements(candidate_id);
#endif
        //getchar();

        ////???q_o????????????????????? ??????rare)????????????????????????(none-rare)????????????????????? ??????,?????????????????????keyword( |inv(term)|<K )??? user
        set<int> rareKeys;  set<int> nonRareKeys;
        set<int> rareKeywordRelated_usrSet;
        for(int term: qo_keys){
            int inv_olist_size = getTermOjectInvListSize(term);
            if(inv_olist_size>K){
                vector<int> usr_termRelated = getTermUserInvList(term);
                int inv_ulist_size = usr_termRelated.size();
                if(inv_ulist_size> uterm_SPARCITY_VALUE){
                    nonRareKeys.insert(term);
                    printf("term %d ????????????????????????u_frequency=%d\n", term,inv_ulist_size);
                }
                else{ ////??????????????????????????????user???????????????????????????????????????????????????
                    for(int user_id: usr_termRelated){
#ifdef TRACKNVDBATCH
                        cout<<"?????? ????????????????????????2???u"<<user_id<<" for poi"<<poi.id<<", user?????????";
#endif
                        if(user_id>=UserID_MaxKey) continue;
                        int up_key = user_id + UserID_MaxKey* poi_th;   //???????????????poi?????????0,usrnum-1)

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
                        //???????????????????????????
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

                            //?????? user-poi associating list
                            candidate_usr_related_store[user_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                            if(candidate_usr_related_store[user_id].size()>K){
                                candidate_usr_related_store[user_id].pop();
                            }
                            //??????candidate user set
                            BatchVerifyEntry ve = BatchVerifyEntry(user_id,-1,u_lcl_rk,u);
                            candidate_User.insert(ve);
                            user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                            cout<<"social_textual find candidate ????????? u"<<usr_id<<endl;
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
        cout<<"??????????????????????????????????????????????????????????????????????????????"<<endl;
#endif
        ////??????????????????????????????????????????????????????????????????????????????
        for(int term: nonRareKeys){
            vector<int> term_related_user = getTermUserInvList(term);
            for(int user_id: term_related_user){
                if(user_id>=UserID_MaxKey) continue;
                int up_key = user_id + UserID_MaxKey* poi_th;   //???????????????poi?????????0,usrnum-1)
                if(userMark[up_key] == true) continue; ////???????????????poi???
                poiRelated_usr.insert(user_id);
            }

        }
        for(int u_id: poiRelated_usr){
            User user = getUserFromO2UOrgLeafData(u_id);
            int u_Ni = user.Ni;
            int u_leaf = Nodes[u_Ni].gtreepath.back();
            pnodes_uOCCURListMap[poi_id][u_leaf].insert(u_id);  //??????????????????uOCCURList
            vector<int> u_keywords = user.keywords;
            for(int term:u_keywords){
                if(isCovered(term,nonRareKeys)){
                    for(int term2:u_keywords){
                        if(isCovered(term2,nonRareKeys)){
                            pnodes_uCLUSTERListMap[poi_id][u_leaf][term].insert(term2); //??????????????????uCLUSTERList
                        }

                    }
                }

            }

            //?????????????????????uOCCURList???uCLUSTERList
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


        ////??????????????????????????????user ????????????user???????????????????????????
        int _size = rareKeywordRelated_usrSet.size();
        //int _size2 = midKeywordRelated_usrSet.size();
        //cout<<"??????POI???????????????????????????user??????="<<_size<<endl;

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
            //???????????????????????????
            //cout<<"?????? ???????????????????????????u"<<usr_id<<endl;
            double gsk_score = getSKScore_o2u_phl(a,alpha,poi,u);
            double u_lcl_rk = getUserLCL_NVD(u,K,a,alpha);
            if(gsk_score > u_lcl_rk){
                //?????? user-poi associating list
                candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                if(candidate_usr_related_store[usr_id].size()>K){
                    candidate_usr_related_store[usr_id].pop();
                }
                //??????candidate user set
                BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_lcl_rk,u);
                candidate_User.insert(ve);
                user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                cout<<"rare textual find candidate ????????? u"<<usr_id<<endl;
#endif
            }
        }

#ifdef TRACKNVDBATCH
        cout<<"??????????????????????????? user ????????????="<<rareKeywordRelated_usrSet.size()<<endl;
    //cout<<"????????????????????????????????? user?????????????????? candidate user size="<<candidate_id.size()<<"???????????????"<<endl;
    //printElements(candidate_id);
#endif
        //getchar();

        ////1. q_o ??????leaf ????????????????????? user?????????????????????
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
        cout<<"q_o ??????leaf"<<qo_leaf<<"????????????????????? user?????????????????????"<<endl;
    } else{
        cout<<"q_o ??????leaf "<<qo_leaf<<"????????????????????????user"<<endl;
    }

        cout<<"??????user_nearSet???,Queue.size="<<Queue.size()<<endl;
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
        cout<<"??????user_nearSet???,Queue.size="<<Queue.size()<<endl;
#endif
        //getchar();


        //????????????level?????????Gtree Node
        int node_highest = qo_leaf;
        POIHighestNode pair(poi_id,node_highest);
        T_max.push(pair);

        poi_th ++;

    }

    initial_endTime = clock();

    //cout<<"??????????????????????????????????????????????????????user??????="<<midKeywordRelated_usrSet.size()<<endl;
    //cout<<"????????????????????????????????????qo_leaf???????????????????????????"<<endl;
    double initial_time = (double) (initial_endTime-initial_startTime)/CLOCKS_PER_SEC*1000;

    //getchar();



    ////2.??????qo_leaf?????????????????????,????????????????????????Gtree??????qo_leaf???????????????????????????????????????????????????????????????????????????????????? user,user_leaf, user_node

    map<int,map<int,bool>> borderEntryVist;  //border_id, node_id
    map<int,bool> userEntryVist;
    while_startTime = clock();
    while(Queue.size()>0 ||T_max.size()>0){
        if(Queue.empty()){ //????????????????????????????????? ???????????????????????????
#ifdef TRACKNVDBATCH
            cout<<"????????????poi???????????????????????????user?????????????????????!"<<endl;
#endif

            BatchExtend_Range(T_max,global_itm,Queue,nonRareKeywordMap,candidate_POIs, K,a,alpha);
            if(Queue.empty()) break;  //????????????????????????

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

                //??????candidate user set
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
                cout<<"----DEQUEUE: border entry: b"<<border_id<<" in leaf"<<node_id<<"for p"<<poi_id<<"????????????------"<<endl;
            else
                cout<<"----DEQUEUE: border entry: b"<<border_id<<" in node"<<node_id<<"for p"<<poi_id<<"????????????------"<<endl;


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
                cout<<"?????????????????? ****prune???"<<endl;
#endif
                continue;
            }
#ifdef TRACK
            cout<<"?????????**??????**???prune???"<<endl;
#endif

            //?????????????????????????????????????????????user
            if(GTree[node_id].isleaf){
#ifdef TRACK
                cout<<"b"<<border_id<<", ??????leaf"<<node_id<<endl;
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
                //???????????????????????????user
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
                    cout<<"???main?????????user entry:"; entry.printRlt();
                    cout<<"Queue.size="<<Queue.size()<<endl;
#endif
                }
            }
            else{
                //??????????????????????????????????????????????????????(line 22-25)
#ifdef TRACK
                cout<<"b"<<border_id<<", ??????node"<<node_id<<endl;
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
                cout<<"???????????????????????????"; printSetElements(child_set);
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
                        cout<<"???main?????????border entry:"; entry.printRlt();
                        cout<<"Queue.size="<<Queue.size()<<endl;
#endif


                    }
                }

            }
        }

    }
    filter_endTime = clock();


    ////?????????????????????user, ????????????TkGSKQ??????

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
            //???????????????gsk??????
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
        cout<<"??????????????????u"<<usr_id<<"???topk??????"<<endl;
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


        //????????????
        if( gsk_score_max >= Rk_u){
#ifdef TRACKNVDBATCH_VERIFY
            cout<<"*************??????????????????u"<<usr_id<<",Rk_u="<<Rk_u<<endl;
#endif
            //???Lu????????????query object????????????
            while(!Lu.empty()){
                ResultLargerFisrt rlf = Lu.top();
                Lu.pop();
                double rel = 0; double inf = 0;
                double gsk_score = rlf.score; int query_object = rlf.o_id;
                if(gsk_score > Rk_u){
                    ResultDetail rd(usr_id, 0, 0, 0, gsk_score, Rk_u);
                    batch_results[query_object].push_back(rd);
#ifdef TRACK
                    cout<<"??????o"<<query_object<<"????????????????????????gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
#endif
                }
                else break;   //fanz
            }
        }




#ifdef TRACKNVDBATCH_VERIFY
        cout<<"??????????????????"<<usr_id<<", gsk_score_max="<<gsk_score_max<<" ,Rk_u="<<Rk_u<<endl;
#endif

    }
    //delete []poiMark; delete []poiADJMark;

    verification_endTime = clock();
    endTime = clock();


    printBatchRkGSKQResults(stores, batch_results);
    cout<<"filter??????????????????"<<candidate_User.size()<<"??????????????????"<<endl;
    Uc_size = candidate_User.size();
    cout<<"Filter ?????? ?????????"<< (double)(filter_endTime - filter_startTime) / CLOCKS_PER_SEC*1000  <<
        "ms (?????? while loop??????:"<< (double)(filter_endTime - while_startTime) / CLOCKS_PER_SEC*1000<<"ms)"<<endl;
    cout<<"Verification ?????? ?????????"<< (double)(verification_endTime - verification_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;

    //cout<<"user_textual_count="<<user_textual_count<<",user_textual_useful_count="<<user_textual_useful_count<<endl;
    cout<<"---------BatchRkGSKQ_by_NVD ?????????"<< (double)(endTime - startTime) / CLOCKS_PER_SEC*1000  << "ms-----------" << endl;


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


    //??????????????????
    set<BatchVerifyEntry> candidate_User;  //candidate users
    unordered_map<int, priority_queue<ResultCurrent>> candidate_usr_related_store;   // user_id, {related poi1,2,....}


    //border???????????????????????????
    int node_size = GTree.size();
    int inner_size = node_size*VertexNum;  // func(bEntry_key)= v_id+v_num*node_id + v_num*node_num*poi_th
    long bEntry_MaxKey =  inner_size * stores.size();
    unordered_map<int,bool> borderMark;
    /*bool* borderMark = new bool[bEntry_MaxKey];
    for(int i=0;i< bEntry_MaxKey;i++){
        borderMark[i] = false;
    }*/


    //user??????????????????????????????
    int upEntry_MaxKey =  UserID_MaxKey * stores.size();
    //unordered_map<int,bool> upEntryMark;
    bool* upEntryMark = new bool[upEntry_MaxKey];
    for(int i=0;i<upEntry_MaxKey;i++){
        upEntryMark[i] = false;
    }
    //memset(upEntryMark,false,sizeof(upEntryMark));

    bool userMark[UserID_MaxKey];
    memset(userMark,false, sizeof(userMark));

    //??????gtree node??? poi ??????????????????
    int np_maxKey = GTree.size()*stores.size();
    bool nodepMark[np_maxKey];
    memset(nodepMark,false, sizeof(nodepMark));


    startTime = clock();
    filter_startTime = clock();
    //??????????????????

    ////?????????gtree node??? uOCCURList?????? uCLUSTERList
    unordered_map<int,UOCCURList> pnodes_uOCCURListMap;
    unordered_map<int,UCLUSTERList> pnodes_uCLUSTERListMap;



    priority_queue<BatchCheckEntry> Queue;
    priority_queue<POIHighestNode> T_max;


    int root = 0;
    //// 1. ????????????poi??????, ??????????????? poi ??? o_leaf??????????????????border?????????
    unordered_map<int, unordered_map<int, vector<int>>> global_itm;
    int poi_th = 0;
    set<int> midKeywordRelated_usrSet;
    for(int poi_id: stores){
        poi_Id2IdxMap[poi_id] = poi_th;


        POI poi = getPOIFromO2UOrgLeafData(poi_id);

        candidate_POIs[poi_id] = poi;

        //printPOIInfo(poi);


        //?????????????????????poi????????? nodes_uOCCURList ??? nodes_uCLUSTERList
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
        cout<<"POI?????????"; printPOIInfo(poi);
        cout<<"o_leaf???leaf"<<qo_leaf<<endl;
        cout<<"o_leaf???user?????????????????????"<<endl;
        printSetElements(tnData.userUKeySet);
#endif

        //// ?????? ?????? poi ??? o_leaf??????????????????border?????????
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

        ////?????????poi ????????????????????? user
        vector<int> checkin_users = poiCheckInIDList[poi_id];
        int user_textual_count = 0; int user_textual_useful_count = 0;
        for(int check_usr_id: checkin_users){
            vector<int> friends = friendshipMap[check_usr_id];
            for(int friend_id: friends){
                if(friend_id>=UserID_MaxKey) continue;
                int up_key = friend_id + UserID_MaxKey* poi_th;   //???????????????poi?????????0,usrnum-1)

                if(upEntryMark[up_key] == true) continue;
                upEntryMark[up_key] = true;


                User u = getUserFromO2UOrgLeafData(friend_id);
                int usr_id = u.id;
                double relevance = 0;
                relevance = textRelevance(u.keywords, poi.keywords);
                //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
                if(relevance==0) continue;
                //???????????????????????????
#ifdef TRACKNVDBATCH
                cout<<"?????? ???????????????????????????u"<<usr_id<<endl;
#endif
                user_textual_count ++;
                double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,u);
                double u_lcl_rk = getUserLCL_NVD(u,K,a,alpha);
                if(gsk_score > u_lcl_rk){

                    //?????? user-poi associating list
                    candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                    if(candidate_usr_related_store[usr_id].size()>K){
                        candidate_usr_related_store[usr_id].pop();
                    }
                    //??????candidate user set
                    BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_lcl_rk,u);
                    candidate_User.insert(ve);

                    user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                    cout<<"social_textual find candidate ????????? u"<<usr_id<<endl;
#endif
                }

            }
        }
#ifdef TRACKNVDBATCH
        cout<<"????????????????????? user????????????="<<user_textual_count<<endl;
    //cout<<"??????????????????????????? user?????????????????? candidate user size="<<candidate_id.size()<<"???????????????"<<endl;
    //printElements(candidate_id);
#endif
        //getchar();

        ////???q_o????????????????????? ??????rare)????????????????????????(none-rare)????????????????????? ??????,?????????????????????keyword( |inv(term)|<K )??? user
        set<int> rareKeys;  set<int> nonRareKeys;
        set<int> rareKeywordRelated_usrSet;
        for(int term: qo_keys){
            int inv_olist_size = getTermOjectInvListSize(term);
            if(inv_olist_size>K){
                vector<int> usr_termRelated = getTermUserInvList(term);
                int inv_ulist_size = usr_termRelated.size();
                if(inv_ulist_size> uterm_SPARCITY_VALUE){
                    nonRareKeys.insert(term);
                    //printf("term %d ????????????????????????u_frequency=%d\n", term,inv_ulist_size);
                }
                else{ ////??????????????????????????????user???????????????????????????????????????????????????
                    for(int user_id: usr_termRelated){
                        if(user_id>=UserID_MaxKey) continue;
                        int up_key = user_id + UserID_MaxKey* poi_th;   //???????????????poi?????????0,usrnum-1)

                        if(upEntryMark[up_key] == true) continue;
                        upEntryMark[up_key] = true;
                        midKeywordRelated_usrSet.insert(user_id);

                        User u = getUserFromO2UOrgLeafData(user_id);
                        int usr_id = u.id;
                        double relevance = 0;
                        relevance = textRelevance(u.keywords, poi.keywords);
                        //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
                        if(relevance==0) continue;
                        //???????????????????????????
#ifdef TRACKNVDBATCH
                        cout<<"?????? ???????????????????????????u"<<usr_id<<endl;
#endif
                        user_textual_count ++;
                        double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,u);
                        double u_lcl_rk = getUserLCL_NVD(u,K,a,alpha);
                        if(gsk_score > u_lcl_rk){

                            //?????? user-poi associating list
                            candidate_usr_related_store[user_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                            if(candidate_usr_related_store[user_id].size()>K){
                                candidate_usr_related_store[user_id].pop();
                            }
                            //??????candidate user set
                            BatchVerifyEntry ve = BatchVerifyEntry(user_id,-1,u_lcl_rk,u);
                            candidate_User.insert(ve);
                            user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                            cout<<"social_textual find candidate ????????? u"<<usr_id<<endl;
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

        ////??????????????????????????????????????????????????????????????????????????????
        for(int term: nonRareKeys){
            vector<int> term_related_user = getTermUserInvList(term);
            for(int user_id: term_related_user){
                if(user_id>=UserID_MaxKey) continue;
                int up_key = user_id + UserID_MaxKey* poi_th;   //???????????????poi?????????0,usrnum-1)
                if(upEntryMark[up_key] == true) continue; ////???????????????poi???
                poiRelated_usr.insert(user_id);
            }

        }
        for(int u_id: poiRelated_usr){
            User user = getUserFromO2UOrgLeafData(u_id);
            int u_Ni = user.Ni;
            int u_leaf = Nodes[u_Ni].gtreepath.back();
            pnodes_uOCCURListMap[poi_id][u_leaf].insert(u_id);  //??????????????????uOCCURList
            vector<int> u_keywords = user.keywords;
            for(int term:u_keywords){
                if(isCovered(term,nonRareKeys)){
                    for(int term2:u_keywords){
                        if(isCovered(term2,nonRareKeys)){
                            pnodes_uCLUSTERListMap[poi_id][u_leaf][term].insert(term2); //??????????????????uCLUSTERList
                        }

                    }
                }

            }

            //?????????????????????uOCCURList???uCLUSTERList
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


        ////??????????????????????????????user ????????????user???????????????????????????
        int _size = rareKeywordRelated_usrSet.size();
        //int _size2 = midKeywordRelated_usrSet.size();
        //cout<<"??????POI???????????????????????????user??????="<<_size<<endl;

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
            //???????????????????????????
            //cout<<"?????? ???????????????????????????u"<<usr_id<<endl;
            double gsk_score = getSKScore_o2u_phl(a,alpha,poi,u);
            double u_lcl_rk = getUserLCL_NVD(u,K,a,alpha);
            if(gsk_score > u_lcl_rk){
                //?????? user-poi associating list
                candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                if(candidate_usr_related_store[usr_id].size()>K){
                    candidate_usr_related_store[usr_id].pop();
                }
                //??????candidate user set
                BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_lcl_rk,u);
                candidate_User.insert(ve);
                user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                cout<<"rare textual find candidate ????????? u"<<usr_id<<endl;
#endif
            }
        }

#ifdef TRACKNVDBATCH
        cout<<"??????????????????????????? user ????????????="<<rareKeywordRelated_usrSet.size()<<endl;
    //cout<<"????????????????????????????????? user?????????????????? candidate user size="<<candidate_id.size()<<"???????????????"<<endl;
    //printElements(candidate_id);
#endif
        //getchar();

        ////1. q_o ??????leaf ????????????????????? user?????????????????????
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
        cout<<"q_o ??????leaf"<<qo_leaf<<"????????????????????? user?????????????????????"<<endl;
    } else{
        cout<<"q_o ??????leaf "<<qo_leaf<<"????????????????????????user"<<endl;
    }

        cout<<"??????user_nearSet???,Queue.size="<<Queue.size()<<endl;
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
        cout<<"??????user_nearSet???,Queue.size="<<Queue.size()<<endl;
#endif
        //getchar();


        //????????????level?????????Gtree Node
        int node_highest = qo_leaf;
        POIHighestNode pair(poi_id,node_highest);
        T_max.push(pair);

        poi_th ++;

    }

    cout<<"??????????????????????????????????????????????????????user??????="<<midKeywordRelated_usrSet.size()<<endl;
    cout<<"????????????????????????????????????qo_leaf???????????????????????????"<<endl;
    //getchar();



    ////2.??????qo_leaf?????????????????????,????????????????????????Gtree??????qo_leaf???????????????????????????????????????????????????????????????????????????????????? user,user_leaf, user_node

    map<int,map<int,bool>> borderEntryVist;  //border_id, node_id
    map<int,bool> userEntryVist;
    while_startTime = clock();
    while(Queue.size()>0 ||T_max.size()>0){
        if(Queue.empty()){ //????????????????????????????????? ???????????????????????????
#ifdef TRACKEXTEND
            cout<<"????????????poi???????????????????????????user?????????????????????!"<<endl;
#endif

            BatchExtend_Range(T_max,global_itm,Queue,nonRareKeywordMap,candidate_POIs, K,a,alpha);
            if(Queue.empty()) break;  //????????????????????????

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

                //??????candidate user set
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
                cout<<"----DEQUEUE: border entry: b"<<border_id<<" in leaf"<<node_id<<"for p"<<poi_id<<"????????????------"<<endl;
            else
                cout<<"----DEQUEUE: border entry: b"<<border_id<<" in node"<<node_id<<"for p"<<poi_id<<"????????????------"<<endl;

            if(node_id==1190){
                cout<<"node_id==1190???"<<endl;
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
                    //cout<<"u"<<usr_id<<"??????????????????"<<endl;
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
                cout<<"?????????????????? ****prune???"<<endl;
#endif
                continue;
            }
#ifdef TRACK
            cout<<"?????????**??????**???prune???"<<endl;
#endif

            //?????????????????????????????????????????????user
            if(GTree[node_id].isleaf){
#ifdef TRACK
                cout<<"b"<<border_id<<", ??????leaf"<<node_id<<endl;
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
                //???????????????????????????user
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
                    cout<<"???main?????????user entry:"; entry.printRlt();
                    cout<<"Queue.size="<<Queue.size()<<endl;
#endif
                }
            }
            else{
                //??????????????????????????????????????????????????????(line 22-25)
#ifdef WhileTRACK
                cout<<"b"<<border_id<<", ??????node"<<node_id<<endl;
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
                cout<<"???????????????????????????"; printSetElements(child_set);
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
                        cout<<"???main?????????border entry:"; entry.printRlt();
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

    ////?????????????????????user, ????????????TkGSKQ??????

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
            //???????????????gsk??????
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
        cout<<"??????????????????u"<<usr_id<<"???topk??????"<<endl;
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


        //????????????
        if( gsk_score_max >= Rk_u){
#ifdef TRACKNVDBATCH_VERIFY
            cout<<"*************??????????????????u"<<usr_id<<",Rk_u="<<Rk_u<<endl;
#endif
            //???Lu????????????query object????????????
            while(!Lu.empty()){
                ResultLargerFisrt rlf = Lu.top();
                Lu.pop();
                double rel = 0; double inf = 0;
                double gsk_score = rlf.score; int query_object = rlf.o_id;
                if(gsk_score > Rk_u){
                    ResultDetail rd(usr_id, 0, 0, 0, gsk_score, Rk_u);
                    batch_results[query_object].push_back(rd);
#ifdef TRACK
                    cout<<"??????o"<<query_object<<"????????????????????????gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
#endif
                }
                else break;   //fanz
            }
        }




#ifdef TRACKNVDBATCH_VERIFY
        cout<<"??????????????????"<<usr_id<<", gsk_score_max="<<gsk_score_max<<" ,Rk_u="<<Rk_u<<endl;
#endif

    }
    //delete []poiMark; delete []poiADJMark;

    verification_endTime = clock();
    endTime = clock();


    printBatchRkGSKQResults(stores, batch_results);
    cout<<"filter??????????????????"<<candidate_User.size()<<"??????????????????"<<endl;
    Uc_size = candidate_User.size();
    cout<<"Filter ?????? ?????????"<< (double)(filter_endTime - filter_startTime) / CLOCKS_PER_SEC*1000  <<
        "ms (?????? while loop??????:"<< (double)(filter_endTime - while_startTime) / CLOCKS_PER_SEC*1000<<"ms)"<<endl;
    cout<<"Verification ?????? ?????????"<< (double)(verification_endTime - verification_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;

    //cout<<"user_textual_count="<<user_textual_count<<",user_textual_useful_count="<<user_textual_useful_count<<endl;
    cout<<"---------BatchRkGSKQ_by_NVD ?????????"<< (double)(endTime - startTime) / CLOCKS_PER_SEC*1000  << "ms-----------" << endl;


    return batch_results;
}


//??????????????????poi???potential customers ???candidate


//Batch Reverse ????????????????????????
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


    //??????????????????
    set<BatchVerifyEntry> candidate_User;  //candidate users
    unordered_map<int, priority_queue<ResultCurrent>> candidate_usr_related_store;   // user_id, {related poi1,2,....}
    unordered_map<int, priority_queue<ResultCurrent>> candidate_usrLeaf_related_store;   // user_id, {related poi1,2,....}


    //border???????????????????????????
    int node_size = GTree.size();
    int inner_size = node_size*VertexNum;  // func(bEntry_key)= v_id+v_num*node_id + v_num*node_num*poi_th
    long bEntry_MaxKey =  inner_size * stores.size();
    unordered_map<int,bool> borderMark;



    //user??????????????????????????????
    int upEntry_MaxKey =  UserID_MaxKey * stores.size();
    //unordered_map<int,bool> upEntryMark;
    bool* upEntryMark = new bool[upEntry_MaxKey];
    for(int i=0;i<upEntry_MaxKey;i++){
        upEntryMark[i] = false;
    }

    bool userMark[UserID_MaxKey];
    memset(userMark,false, sizeof(userMark));

    //??????gtree node??? poi ??????????????????
    int np_maxKey = GTree.size()*stores.size();
    bool nodepMark[np_maxKey];
    memset(nodepMark,false, sizeof(nodepMark));


    startTime = clock();
    filter_startTime = clock();
    //??????????????????

    ////?????????gtree node??? uOCCURList?????? uCLUSTERList
    unordered_map<int,UOCCURList> pnodes_uOCCURListMap;
    unordered_map<int,UCLUSTERList> pnodes_uCLUSTERListMap;
    //unordered_map<int,UOCCURListMap> pnodes_uOCCURListMap;
    //unordered_map<int,UCLUSTERListMap> pnodes_uCLUSTERListMap;



    priority_queue<BatchCheckEntry> Queue;
    priority_queue<POIHighestNode> T_max;


    int root = 0;
    //// 1. ????????????poi??????, ??????????????? poi ??? o_leaf??????????????????border?????????
    unordered_map<int, unordered_map<int, vector<int>>> global_itm;
    int poi_th = 0;
    set<int> midKeywordRelated_usrSet; double initial_time = -1;
    initial_startTime = clock();
    for(int poi_id: stores){
        poi_Id2IdxMap[poi_id] = poi_th;


        POI poi = getPOIFromO2UOrgLeafData(poi_id);

        candidate_POIs[poi_id] = poi;

        ////printPOIInfo(poi);


        //?????????????????????poi????????? nodes_uOCCURList ??? nodes_uCLUSTERList
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
        cout<<"o_leaf???leaf"<<qo_leaf<<endl;
        cout<<"o_leaf???user?????????????????????"<<endl;
        printSetElements(tnData.userUKeySet);
#endif

        //// ?????? ?????? poi ??? o_leaf??????????????????border?????????
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

        ////?????????poi ????????????????????? user
        vector<int> checkin_users = poiCheckInIDList[poi_id];
        int user_textual_count = 0; int user_textual_useful_count = 0;
        for(int check_usr_id: checkin_users){
            vector<int> friends = friendshipMap[check_usr_id];
            for(int friend_id: friends){
                if(friend_id>=UserID_MaxKey) continue;
                int up_key = friend_id + UserID_MaxKey* poi_th;   //???????????????poi?????????0,usrnum-1)

                if(upEntryMark[up_key] == true) continue;
                upEntryMark[up_key] = true;


                User u = getUserFromO2UOrgLeafData(friend_id);
                int usr_id = u.id;
                double relevance = 0;
                relevance = textRelevance(u.keywords, poi.keywords);
                //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
                if(relevance==0) continue;
                //???????????????????????????
#ifdef TRACKNVDBATCH
                cout<<"?????? ???????????????????????????u"<<usr_id<<endl;
#endif
                user_textual_count ++;
                double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,u);
                double u_lcl_rk = getUserLCL_NVDFaster(u,K,a,alpha);;//getUserLCL_NVD(u,K,a,alpha);
                if(gsk_score > u_lcl_rk){

                    //?????? user-poi associating list
                    candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                    if(candidate_usr_related_store[usr_id].size()>K){
                        candidate_usr_related_store[usr_id].pop();
                    }
                    //??????candidate user set
                    BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_lcl_rk,u);
                    candidate_User.insert(ve);

                    user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                    cout<<"social_textual find candidate ????????? u"<<usr_id<<endl;
#endif
                }

            }
        }
#ifdef TRACKNVDBATCH
        cout<<"????????????????????? user????????????="<<user_textual_count<<endl;
    //cout<<"??????????????????????????? user?????????????????? candidate user size="<<candidate_id.size()<<"???????????????"<<endl;
    //printElements(candidate_id);
#endif
        //getchar();

        ////???q_o????????????????????? ??????rare)????????????????????????(none-rare)????????????????????? ??????,?????????????????????keyword( |inv(term)|<K )??? user
        set<int> rareKeys;  set<int> nonRareKeys;
        set<int> rareKeywordRelated_usrSet;
        for(int term: qo_keys){
            int inv_olist_size = getTermOjectInvListSize(term);
            if(inv_olist_size>K){
                vector<int> usr_termRelated = getTermUserInvList(term);
                int inv_ulist_size = usr_termRelated.size();
                if(inv_ulist_size> uterm_SPARCITY_VALUE){
                    nonRareKeys.insert(term);
                    ////printf("term %d ????????????????????????u_frequency=%d\n", term,inv_ulist_size);
                }
                else{ ////??????????????????????????????user???????????????????????????????????????????????????
                    for(int user_id: usr_termRelated){
                        if(user_id>=UserID_MaxKey) continue;
                        int up_key = user_id + UserID_MaxKey* poi_th;   //???????????????poi?????????0,usrnum-1)

                        if(upEntryMark[up_key] == true) continue;
                        upEntryMark[up_key] = true;
                        midKeywordRelated_usrSet.insert(user_id);

                        User u = getUserFromO2UOrgLeafData(user_id);
                        int usr_id = u.id;
                        double relevance = 0;
                        relevance = textRelevance(u.keywords, poi.keywords);
                        //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
                        if(relevance==0) continue;
                        //???????????????????????????
#ifdef TRACKNVDBATCH
                        cout<<"?????? ???????????????????????????u"<<usr_id<<endl;
#endif
                        user_textual_count ++;
                        double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,u);
                        double u_lcl_rk = getUserLCL_NVDFaster(u,K,a,alpha);  //0
                        if(gsk_score > u_lcl_rk){

                            //?????? user-poi associating list
                            candidate_usr_related_store[user_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                            if(candidate_usr_related_store[user_id].size()>K){
                                candidate_usr_related_store[user_id].pop();
                            }
                            //??????candidate user set
                            BatchVerifyEntry ve = BatchVerifyEntry(user_id,-1,u_lcl_rk,u);
                            candidate_User.insert(ve);
                            user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                            cout<<"social_textual find candidate ????????? u"<<usr_id<<endl;
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

        ////??????????????????????????????????????????????????????????????????????????????
        for(int term: nonRareKeys){
            vector<int> term_related_user = getTermUserInvList(term);
            for(int user_id: term_related_user){
                if(user_id>=UserID_MaxKey) continue;
                int up_key = user_id + UserID_MaxKey* poi_th;   //???????????????poi?????????0,usrnum-1)
                if(upEntryMark[up_key] == true) continue; ////???????????????poi???
                poiRelated_usr.insert(user_id);
            }

        }
        for(int u_id: poiRelated_usr){
            User user = getUserFromO2UOrgLeafData(u_id);
            int u_Ni = user.Ni;
            int u_leaf = Nodes[u_Ni].gtreepath.back();
            pnodes_uOCCURListMap[poi_id][u_leaf].insert(u_id);  //??????????????????uOCCURList
            vector<int> u_keywords = user.keywords;
            for(int term:u_keywords){
                if(isCovered(term,nonRareKeys)){
                    for(int term2:u_keywords){
                        if(isCovered(term2,nonRareKeys)){
                            pnodes_uCLUSTERListMap[poi_id][u_leaf][term].insert(term2); //??????????????????uCLUSTERList
                        }

                    }
                }

            }

            //?????????????????????uOCCURList???uCLUSTERList
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


        ////??????????????????????????????user ????????????user???????????????????????????
        int _size = rareKeywordRelated_usrSet.size();
        //int _size2 = midKeywordRelated_usrSet.size();
        ////cout<<"??????POI???????????????????????????user??????="<<_size<<endl;

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
            //???????????????????????????
            //cout<<"?????? ???????????????????????????u"<<usr_id<<endl;
            double gsk_score = getSKScore_o2u_phl(a,alpha,poi,u);
            double u_lcl_rk = getUserLCL_NVDFaster(u,K,a,alpha);
            if(gsk_score > u_lcl_rk){
                //?????? user-poi associating list
                candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                if(candidate_usr_related_store[usr_id].size()>K){
                    candidate_usr_related_store[usr_id].pop();
                }
                //??????candidate user set
                BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_lcl_rk,u);
                candidate_User.insert(ve);
                user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                cout<<"rare textual find candidate ????????? u"<<usr_id<<endl;
#endif
            }
        }

#ifdef TRACKNVDBATCH
        cout<<"??????????????????????????? user ????????????="<<rareKeywordRelated_usrSet.size()<<endl;
    //cout<<"????????????????????????????????? user?????????????????? candidate user size="<<candidate_id.size()<<"???????????????"<<endl;
    //printElements(candidate_id);
#endif
        //getchar();

        ////1. q_o ??????leaf ????????????????????? user?????????????????????
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
        cout<<"q_o ??????leaf"<<qo_leaf<<"????????????????????? user?????????????????????"<<endl;
    } else{
        cout<<"q_o ??????leaf "<<qo_leaf<<"????????????????????????user"<<endl;
    }

        cout<<"??????user_nearSet???,Queue.size="<<Queue.size()<<endl;
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
        cout<<"??????user_nearSet???,Queue.size="<<Queue.size()<<endl;
#endif
        //getchar();


        //????????????level?????????Gtree Node
        int node_highest = qo_leaf;
        POIHighestNode pair(poi_id,node_highest);
        T_max.push(pair);

        poi_th ++;

    }
    initial_endTime = clock();
    initial_time = (double) (initial_endTime-initial_startTime)/CLOCKS_PER_SEC*1000;

    cout<<"??????????????????????????????????????????????????????user??????="<<midKeywordRelated_usrSet.size()<<endl;
    cout<<"????????????????????????????????????qo_leaf???????????????????????????"<<endl;
    //getchar();


    ////2.??????qo_leaf?????????????????????,????????????????????????Gtree??????qo_leaf???????????????????????????????????????????????????????????????????????????????????? user,user_leaf, user_node

    map<int,map<int,bool>> borderEntryVist;  //border_id, node_id
    map<int,bool> userEntryVist;
    while_startTime = clock();
    while(Queue.size()>0 ||T_max.size()>0){
        if(Queue.empty()){ //????????????????????????????????? ???????????????????????????
#ifdef TRACKEXTEND
            cout<<"????????????poi???????????????????????????user?????????????????????!"<<endl;
#endif

            BatchExtend_Range(T_max,global_itm,Queue,nonRareKeywordMap,candidate_POIs, K,a,alpha);
            if(Queue.empty()) break;  //????????????????????????

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

                //??????candidate user set
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
            //?????????????????????????????????????????????????????????user
            if(GTree[node_id].isleaf){
#ifdef TRACK
                cout<<"b"<<border_id<<", ??????leaf"<<node_id<<endl;
#endif
                if(nodepMark[np_key]!=true) {
                    nodepMark[np_key] =true;  //?????????node????????????poi????????????
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

                            //??????candidate user set
                            BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_SB_rk,user);
                            candidate_User.insert(ve);

                        }

                    }
                }


            }
            else{  //???????????????

#ifdef   AGGRESSIVE_I
                int _size = pnodes_uOCCURListMap[poi_id][node_id].size();
                if(_size < NODE_SPARCITY_VALUE_HEURIS){//&& dist_parentBorder>LONG_TERM
                    for(int usr_id: pnodes_uOCCURListMap[poi_id][node_id]){
                        int up_key = usr_id + UserID_MaxKey *poi_th;
                        if(upEntryMark[up_key] == true)
                            continue;
                        else
                            upEntryMark[up_key] = true;
                        //cout<<"u"<<usr_id<<"??????????????????"<<endl;
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

                            //??????candidate user set
                            BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_SB_rk,user);
                            candidate_User.insert(ve);

                        }

                    }
                    nodepMark[np_key] =true;  //?????????node????????????poi????????????
                    continue;
                }
#endif


                set<int> Keys;

                Keys = getBorder_SB_BatchNVD_Cluster(topEntry,K,a, alpha,pnodes_uCLUSTERListMap);

                //cout<<"remained (user) Keywords"; printSetElements(Keys);
                if(Keys.size() ==0) {
#ifdef TRACK
                    cout<<"?????????????????? ****prune???"<<endl;
#endif
                    continue;
                }
#ifdef TRACK
                cout<<"?????????**??????**???prune???"<<endl;
#endif


                //??????????????????????????????????????????????????????(line 22-25)
#ifdef WhileTRACK
                cout<<"b"<<border_id<<", ??????node"<<node_id<<endl;
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
                cout<<"???????????????????????????"; printSetElements(child_set);
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
                        cout<<"???main?????????border entry:"; entry.printRlt();
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


    printf("????????????????????????%f ms !\n", initial_time);
    cout<<"fast filter??????????????????"<<candidate_User.size()<<"??????????????????"<<endl;
    cout<<"Fast Filter ?????? ?????????"<< (double)(filter_endTime - filter_startTime) / CLOCKS_PER_SEC*1000  <<
        "ms (?????? while loop??????:"<< (double)(filter_endTime - while_startTime) / CLOCKS_PER_SEC*1000<<"ms)"<<endl;
    //cout<<"Verification ?????? ?????????"<< (double)(verification_endTime - verification_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;

    //cout<<"user_textual_count="<<user_textual_count<<",user_textual_useful_count="<<user_textual_useful_count<<endl;
    cout<<"---------Batch Fast Filter ?????????"<< (double)(endTime - startTime) / CLOCKS_PER_SEC*1000  << "ms-----------" << endl;


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


    //??????????????????
    set<BatchVerifyEntry> candidate_User;  //candidate users
    unordered_map<int, priority_queue<ResultCurrent>> candidate_usr_related_store;   // user_id, {related poi1,2,....}
    unordered_map<int, priority_queue<ResultCurrent>> candidate_usrLeaf_related_store;   // user_id, {related poi1,2,....}


    //border???????????????????????????
    int node_size = GTree.size();
    int inner_size = node_size*VertexNum;  // func(bEntry_key)= v_id+v_num*node_id + v_num*node_num*poi_th
    long bEntry_MaxKey =  inner_size * stores.size();
    unordered_map<int,bool> borderMark;



    //user??????????????????????????????
    int upEntry_MaxKey =  UserID_MaxKey * stores.size();
    //unordered_map<int,bool> upEntryMark;
    bool* upEntryMark = new bool[upEntry_MaxKey];
    for(int i=0;i<upEntry_MaxKey;i++){
        upEntryMark[i] = false;
    }

    bool userMark[UserID_MaxKey];
    memset(userMark,false, sizeof(userMark));

    //??????gtree node??? poi ??????????????????
    int np_maxKey = GTree.size()*stores.size();
    bool nodepMark[np_maxKey];
    memset(nodepMark,false, sizeof(nodepMark));


    startTime = clock();
    filter_startTime = clock();
    //??????????????????

    ////?????????gtree node??? uOCCURList?????? uCLUSTERList
    unordered_map<int,UOCCURList> pnodes_uOCCURListMap;
    unordered_map<int,UCLUSTERList> pnodes_uCLUSTERListMap;
    //unordered_map<int,UOCCURListMap> pnodes_uOCCURListMap;
    //unordered_map<int,UCLUSTERListMap> pnodes_uCLUSTERListMap;



    priority_queue<BatchCheckEntry> Queue;
    priority_queue<POIHighestNode> T_max;


    int root = 0;
    //// 1. ????????????poi??????, ??????????????? poi ??? o_leaf??????????????????border?????????
    unordered_map<int, unordered_map<int, vector<int>>> global_itm;
    int poi_th = 0;
    set<int> midKeywordRelated_usrSet; double initial_time = -1;
    initial_startTime = clock();
    for(int poi_id: stores){
        poi_Id2IdxMap[poi_id] = poi_th;


        POI poi = getPOIFromO2UOrgLeafData(poi_id);

        candidate_POIs[poi_id] = poi;

        ////printPOIInfo(poi);


        //?????????????????????poi????????? nodes_uOCCURList ??? nodes_uCLUSTERList
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
        cout<<"o_leaf???leaf"<<qo_leaf<<endl;
        cout<<"o_leaf???user?????????????????????"<<endl;
        printSetElements(tnData.userUKeySet);
#endif

        //// ?????? ?????? poi ??? o_leaf??????????????????border?????????
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

        ////?????????poi ????????????????????? user
        vector<int> checkin_users = poiCheckInIDList[poi_id];
        int user_textual_count = 0; int user_textual_useful_count = 0;
        for(int check_usr_id: checkin_users){
            vector<int> friends = friendshipMap[check_usr_id];
            for(int friend_id: friends){
                if(friend_id>=UserID_MaxKey) continue;
                int up_key = friend_id + UserID_MaxKey* poi_th;   //???????????????poi?????????0,usrnum-1)

                if(upEntryMark[up_key] == true) continue;
                upEntryMark[up_key] = true;


                User u = getUserFromO2UOrgLeafData(friend_id);
                int usr_id = u.id;
                double relevance = 0;
                relevance = textRelevance(u.keywords, poi.keywords);
                //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
                if(relevance==0) continue;
                //???????????????????????????
#ifdef TRACKNVDBATCH
                cout<<"?????? ???????????????????????????u"<<usr_id<<endl;
#endif
                user_textual_count ++;
                double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,u);
                double u_lcl_rk = getUserLCL_NVDFaster(u,K,a,alpha);;//getUserLCL_NVD(u,K,a,alpha);
                if(gsk_score > u_lcl_rk){

                    //?????? user-poi associating list
                    candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                    if(candidate_usr_related_store[usr_id].size()>K){
                        candidate_usr_related_store[usr_id].pop();
                    }
                    //??????candidate user set
                    BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_lcl_rk,u);
                    candidate_User.insert(ve);

                    user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                    cout<<"social_textual find candidate ????????? u"<<usr_id<<endl;
#endif
                }

            }
        }
#ifdef TRACKNVDBATCH
        cout<<"????????????????????? user????????????="<<user_textual_count<<endl;
    //cout<<"??????????????????????????? user?????????????????? candidate user size="<<candidate_id.size()<<"???????????????"<<endl;
    //printElements(candidate_id);
#endif
        //getchar();

        ////???q_o????????????????????? ??????rare)????????????????????????(none-rare)????????????????????? ??????,?????????????????????keyword( |inv(term)|<K )??? user
        set<int> rareKeys;  set<int> nonRareKeys;
        set<int> rareKeywordRelated_usrSet;
        for(int term: qo_keys){
            int inv_olist_size = getTermOjectInvListSize(term);
            if(inv_olist_size>K){
                vector<int> usr_termRelated = getTermUserInvList(term);
                int inv_ulist_size = usr_termRelated.size();
                if(inv_ulist_size> uterm_SPARCITY_VALUE){
                    nonRareKeys.insert(term);
                    ////printf("term %d ????????????????????????u_frequency=%d\n", term,inv_ulist_size);
                }
                else{ ////??????????????????????????????user???????????????????????????????????????????????????
                    for(int user_id: usr_termRelated){
                        if(user_id>=UserID_MaxKey) continue;
                        int up_key = user_id + UserID_MaxKey* poi_th;   //???????????????poi?????????0,usrnum-1)

                        if(upEntryMark[up_key] == true) continue;
                        upEntryMark[up_key] = true;
                        midKeywordRelated_usrSet.insert(user_id);

                        User u = getUserFromO2UOrgLeafData(user_id);
                        int usr_id = u.id;
                        double relevance = 0;
                        relevance = textRelevance(u.keywords, poi.keywords);
                        //textual_end = clock(); textual_time += (double) (textual_end - textual_start);
                        if(relevance==0) continue;
                        //???????????????????????????
#ifdef TRACKNVDBATCH
                        cout<<"?????? ???????????????????????????u"<<usr_id<<endl;
#endif
                        user_textual_count ++;
                        double gsk_score = getGSKScore_o2u_phl(a,alpha,poi,u);
                        double u_lcl_rk = getUserLCL_NVDFaster(u,K,a,alpha);  //0
                        if(gsk_score > u_lcl_rk){

                            //?????? user-poi associating list
                            candidate_usr_related_store[user_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                            if(candidate_usr_related_store[user_id].size()>K){
                                candidate_usr_related_store[user_id].pop();
                            }
                            //??????candidate user set
                            BatchVerifyEntry ve = BatchVerifyEntry(user_id,-1,u_lcl_rk,u);
                            candidate_User.insert(ve);
                            user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                            cout<<"social_textual find candidate ????????? u"<<usr_id<<endl;
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

        ////??????????????????????????????????????????????????????????????????????????????
        for(int term: nonRareKeys){
            vector<int> term_related_user = getTermUserInvList(term);
            for(int user_id: term_related_user){
                if(user_id>=UserID_MaxKey) continue;
                int up_key = user_id + UserID_MaxKey* poi_th;   //???????????????poi?????????0,usrnum-1)
                if(upEntryMark[up_key] == true) continue; ////???????????????poi???
                poiRelated_usr.insert(user_id);
            }

        }
        for(int u_id: poiRelated_usr){
            User user = getUserFromO2UOrgLeafData(u_id);
            int u_Ni = user.Ni;
            int u_leaf = Nodes[u_Ni].gtreepath.back();
            pnodes_uOCCURListMap[poi_id][u_leaf].insert(u_id);  //??????????????????uOCCURList
            vector<int> u_keywords = user.keywords;
            for(int term:u_keywords){
                if(isCovered(term,nonRareKeys)){
                    for(int term2:u_keywords){
                        if(isCovered(term2,nonRareKeys)){
                            pnodes_uCLUSTERListMap[poi_id][u_leaf][term].insert(term2); //??????????????????uCLUSTERList
                        }

                    }
                }

            }

            //?????????????????????uOCCURList???uCLUSTERList
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


        ////??????????????????????????????user ????????????user???????????????????????????
        int _size = rareKeywordRelated_usrSet.size();
        //int _size2 = midKeywordRelated_usrSet.size();
        ////cout<<"??????POI???????????????????????????user??????="<<_size<<endl;

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
            //???????????????????????????
            //cout<<"?????? ???????????????????????????u"<<usr_id<<endl;
            double gsk_score = getSKScore_o2u_phl(a,alpha,poi,u);
            double u_lcl_rk = getUserLCL_NVDFaster(u,K,a,alpha);
            if(gsk_score > u_lcl_rk){
                //?????? user-poi associating list
                candidate_usr_related_store[usr_id].push(ResultCurrent(poi_id, gsk_score, gsk_score, 0, 0));
                if(candidate_usr_related_store[usr_id].size()>K){
                    candidate_usr_related_store[usr_id].pop();
                }
                //??????candidate user set
                BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_lcl_rk,u);
                candidate_User.insert(ve);
                user_textual_useful_count++;
#ifdef TRACKNVDBATCH
                cout<<"rare textual find candidate ????????? u"<<usr_id<<endl;
#endif
            }
        }

#ifdef TRACKNVDBATCH
        cout<<"??????????????????????????? user ????????????="<<rareKeywordRelated_usrSet.size()<<endl;
    //cout<<"????????????????????????????????? user?????????????????? candidate user size="<<candidate_id.size()<<"???????????????"<<endl;
    //printElements(candidate_id);
#endif
        //getchar();

        ////1. q_o ??????leaf ????????????????????? user?????????????????????
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
        cout<<"q_o ??????leaf"<<qo_leaf<<"????????????????????? user?????????????????????"<<endl;
    } else{
        cout<<"q_o ??????leaf "<<qo_leaf<<"????????????????????????user"<<endl;
    }

        cout<<"??????user_nearSet???,Queue.size="<<Queue.size()<<endl;
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
        cout<<"??????user_nearSet???,Queue.size="<<Queue.size()<<endl;
#endif
        //getchar();


        //????????????level?????????Gtree Node
        int node_highest = qo_leaf;
        POIHighestNode pair(poi_id,node_highest);
        T_max.push(pair);

        poi_th ++;

    }
    initial_endTime = clock();
    initial_time = (double) (initial_endTime-initial_startTime)/CLOCKS_PER_SEC*1000;

    cout<<"??????????????????????????????????????????????????????user??????="<<midKeywordRelated_usrSet.size()<<endl;
    cout<<"????????????????????????????????????qo_leaf???????????????????????????"<<endl;
    //getchar();


    ////2.??????qo_leaf?????????????????????,????????????????????????Gtree??????qo_leaf???????????????????????????????????????????????????????????????????????????????????? user,user_leaf, user_node

    map<int,map<int,bool>> borderEntryVist;  //border_id, node_id
    map<int,bool> userEntryVist;
    while_startTime = clock();
    while(Queue.size()>0 ||T_max.size()>0){
        if(Queue.empty()){ //????????????????????????????????? ???????????????????????????
#ifdef TRACKEXTEND
            cout<<"????????????poi???????????????????????????user?????????????????????!"<<endl;
#endif

            BatchExtend_Range(T_max,global_itm,Queue,nonRareKeywordMap,candidate_POIs, K,a,alpha);
            if(Queue.empty()) break;  //????????????????????????

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

                //??????candidate user set
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
            //?????????????????????????????????????????????????????????user
            if(GTree[node_id].isleaf){
#ifdef TRACK
                cout<<"b"<<border_id<<", ??????leaf"<<node_id<<endl;
#endif
                if(nodepMark[np_key]!=true) {
                    nodepMark[np_key] =true;  //?????????node????????????poi????????????
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

                            //??????candidate user set
                            BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_SB_rk,user);
                            candidate_User.insert(ve);

                        }

                    }
                }


            }
            else{  //???????????????

#ifdef   AGGRESSIVE_I
                int _size = pnodes_uOCCURListMap[poi_id][node_id].size();
                if(_size < NODE_SPARCITY_VALUE_HEURIS){//&& dist_parentBorder>LONG_TERM
                    for(int usr_id: pnodes_uOCCURListMap[poi_id][node_id]){
                        int up_key = usr_id + UserID_MaxKey *poi_th;
                        if(upEntryMark[up_key] == true)
                            continue;
                        else
                            upEntryMark[up_key] = true;
                        //cout<<"u"<<usr_id<<"??????????????????"<<endl;
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

                            //??????candidate user set
                            BatchVerifyEntry ve = BatchVerifyEntry(usr_id,-1,u_SB_rk,user);
                            candidate_User.insert(ve);

                        }

                    }
                    nodepMark[np_key] =true;  //?????????node????????????poi????????????
                    continue;
                }
#endif


                set<int> Keys;

                Keys = getBorder_SB_BatchNVD_Cluster(topEntry,K,a, alpha,pnodes_uCLUSTERListMap);

                //cout<<"remained (user) Keywords"; printSetElements(Keys);
                if(Keys.size() ==0) {
#ifdef TRACK
                    cout<<"?????????????????? ****prune???"<<endl;
#endif
                    continue;
                }
#ifdef TRACK
                cout<<"?????????**??????**???prune???"<<endl;
#endif


                //??????????????????????????????????????????????????????(line 22-25)
#ifdef WhileTRACK
                cout<<"b"<<border_id<<", ??????node"<<node_id<<endl;
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
                cout<<"???????????????????????????"; printSetElements(child_set);
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
                        cout<<"???main?????????border entry:"; entry.printRlt();
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


    printf("????????????????????????%f ms !\n", initial_time);
    cout<<"fast filter??????????????????"<<candidate_User.size()<<"??????????????????"<<endl;
    cout<<"Fast Filter ?????? ?????????"<< (double)(filter_endTime - filter_startTime) / CLOCKS_PER_SEC*1000  <<
        "ms (?????? while loop??????:"<< (double)(filter_endTime - while_startTime) / CLOCKS_PER_SEC*1000<<"ms)"<<endl;
    //cout<<"Verification ?????? ?????????"<< (double)(verification_endTime - verification_startTime) / CLOCKS_PER_SEC*1000  << "ms" << endl;

    //cout<<"user_textual_count="<<user_textual_count<<",user_textual_useful_count="<<user_textual_useful_count<<endl;
    cout<<"---------Batch Fast Filter ?????????"<< (double)(endTime - startTime) / CLOCKS_PER_SEC*1000  << "ms-----------" << endl;


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

            //if (score > Rk_usr && distance < 10000) {   // 10km????????????(???????????????????????????)
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
        for(int term: POIs[o_id].keywords){ //???topk??????o??????????????????
            if(term_candidateUserInv[term].size()>0){
                for(int usr_id:term_candidateUserInv[term]){ //?????????????????????????????????
                    double rel; double inf;
                    if(Users[usr_id].verified==true) continue; //??????????????????
                    double score_upper = getGSKScore_o2u_Upper(a, alpha, POIs[o_id], Users[usr_id], rel, inf);
                    if(Users[usr_id].current_Results.size()<Qk){  //????????????
                        double score = getGSKScore_o2u(a,  alpha, POIs[o_id], Users[usr_id], rel, inf);
                        Users[usr_id].current_Results.push(Top_result(o_id, score));
                    }
                    else{  //????????????
                        if(score_upper > Users[usr_id].current_Results.top().score){
                            double score = getGSKScore_o2u(a,  alpha, POIs[o_id], Users[usr_id], rel, inf);
                            if(score > Users[usr_id].current_Results.top().score){
                                Users[usr_id].current_Results.push(Top_result(o_id, score));
                                Users[usr_id].current_Results.pop();
                                //?????????????????????????????????topk score
                                Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, Users[usr_id].current_Results.top().score);
                                //cout<<"????????????"<<Users[usr_id].topkScore_current<<endl;
                            }
                        }
                    }
                }
            }
        }
    }
}





/*------------------------- The algorithms below are for Batch RkGSKQ-------------------------*/


//1. ????????????
//filter-base

void  Group_Filter_Baseline_memory(map<int, map<int,bool>>& usr_store_map, int store_partition_id, vector<int> stores, vector<int> keywords, int Qk, float a, float alpha, set<int>& verification_User, map<int,set<VerifyEntry>>& verification_Map, map<int, priority_queue<ResultCurrent>>& candidate_usr_related_store, map<int, LCLResult>& candidate_usr_LCL) {
    vector<int> qKey = keywords;
    set<int> store_leafSet;
    set<int> store_nodes;
    map<int, vector<int>> leaf_store_list;//<leaf_id, <store_id....>>  leaf ??????store?????????
    map<int, set<int>> leaf_verify_list;//<leaf_id, <usr_id....>>  leaf ?????????store???????????????????????????


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
        //cout<<"??????"<<i<<",leaf="<<leaf<<endl;
    }
    cout<<"?????????????????????"<<store_leafSet.size()<<"?????????????????????,GTree???"<<store_nodes.size()<<"??????????????????"<<endl;

    set<int> keyword_Set;
    for(int t: keywords)
        keyword_Set.insert(t);

    //??????partition????????????????????????
    set<int> check_in_partition;
    for(int s: stores){
        for(int u: POIs[s].check_ins){
            check_in_partition.insert(u);
        }
    }

    //???????????????????????????????????????
    set<int> totalUsrSet;
    for(int term:keyword_Set){
        for(int usr: invListOfUser[term]){
            totalUsrSet.insert(usr);
        }

    }


    //??????????????????????????? ???leafNodeInv??????????????????

    int prune=0; int unprune=0; int user_access=0; int usr_prune =0;
    int leaf_prune=0; int leaf_unprune=0;
    int may_prune = 0;
    //int topK_count = 0; //topk????????????
    // ???????????????????????????????????????????????????
    int t3=0; int root=0;



    //cout << " ,relate poi leaf size=" << related_poi_leafSet.size() << endl;
    //cout << " ,relate usr leaf size=" << related_usr_leafSet.size() << endl;

    vector<int> related_usr_node;   //????????????????????????

    map<int,vector<int>> related_poi_nodes;  //???????????????????????????

    // ?????????????????????????????????****??????****?????????????????????????????????, ??????query partition ??????

    vector<int> unodes = GTree[root].children;
    for(int unode: unodes){
        related_usr_node.push_back(unode);
        GTree[unode].ossociate_queryNodes = store_leafSet; //?????????????????????????????????
    }


    // ?????????????????????????????????****?????????***** ?????????????????????????????????
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
            //cout << "node" << node << "?????????" << endl;

            //cout<<"iteration"<<i<<":"<<endl;

            set<int> inter_Key = obtain_itersection_jins(GTree[node].userUKeySet, keyword_Set);

            //??????????????????????????????
            if(!GTree[node].isleaf) {
                int usrNum = 0;
                int store_leaf = store_partition_id;
                if (false) {   //??????????????????  //x
                    //cout << "??????????????????" << node <<",usrNum="<<usrNum<<endl;

                    for(int term: keyword_Set){
                        for(int usr_leaf: GTree[node].term_usrLeaf_Map[term]){
                            related_children.insert(usr_leaf);
                            for(int store_leaf:GTree[node].ossociate_queryNodes){
                                GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                            }
                        }

                    }

                    prune++;

                    //cout<<"?????????????????? usr_node"<<node<<",prune"<<endl;   //??????????????????

                } else {  //?????????????????????????????????????????????
                    //?????? ????????????
                    double max_Score = -1;

                    //max_Score = getGSKScoreP2U_Upper_GivenKey(a, alpha, keywords, store_leaf, node);
                    max_Score = getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, check_in_partition, store_leaf, node); //?????????????????? InterKey, ?????????check-in?????????????????????

                    double Rt_U = 100;
                    //?????????????????????LCL

                    for(int term: inter_Key){
                        LCLResult lcl;
                        map<int,LCLResult>::iterator iter;
                        iter = GTree[node].termLCLDetail.find(term);
                        if(iter != GTree[node].termLCLDetail.end()){ //?????????cache
                            lcl = iter->second;
                        }else{
                            lcl = getUpperNode_termLCL(node,term, related_poi_nodes,Qk,a,alpha,max_Score);
                            GTree[node].termLCLDetail[term] = lcl;
                        }
                        double current_score = lcl.topscore;
                        Rt_U = min(Rt_U, current_score);
                    }


                    //cout<<"GTree[node].related_queryNodes[term_id].size="<<GTree[node].related_queryNodes[term_id].size()<<endl;
                    //?????????????????????????????????

                    if (max_Score > Rt_U) {  // ???????????????????????????????????????????????????????????????????????????
                        unprune++;

                        for(int term: inter_Key){
                            if(GTree[node].inverted_list_u.count(term)>0){
                                for (int child: GTree[node].inverted_list_u[term]) {
                                    //for (int child: textual_child){
                                    //??????????????????????????????lcl??????
                                    GTree[child].termParentLCLDetail[term] = GTree[node].termParentLCLDetail[term];
                                    //GTree[child].termParentFilter[term_id] = filter;
                                    related_children.insert(child);
                                    //cout<<"??????n"<<child<<endl;

                                }
                            }

                        }

                        //getchar();
                    }
                        //???????????????????????????
                    else {
                        prune++;
#ifdef PRUNELOG
                        cout << "??????" << node << "???store_leaf" << store_leaf << "????????? ???????????????????????????" << endl;
#endif
                    }

                    //cout<<"Rt_U="<<Rt_U<<endl;
                    //cout<<"get the LCL of u leaf"<<leaf_id<<endl;


                }//else???????????????
            }

            else if (GTree[node].isleaf) {  //?????????????????????????????????   onceleafsuoz
#ifdef TRACK
                cout << "access usr_leaf" << node << endl;
#endif

                //????????????store_leaf
                int store_leaf = store_partition_id;
                double max_Score =
                        getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, check_in_partition, store_leaf, node);



                MultiLCLResult lclResult_multi;
                vector<int> Ukeys;  //?????????
                for(int term_id: inter_Key)
                    Ukeys.push_back(term_id);

                map<set<int>,MultiLCLResult>::iterator iter;
                iter = GTree[node].cacheMultiLCLDetail.find(inter_Key);
                if (iter!=GTree[node].cacheMultiLCLDetail.end()){
                    lclResult_multi = iter->second;
                }

                    //???cache, ??????
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

                //?????? usr_lcl update_o_set
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
                //?????????????????????????????????????????????
                set<int> L;
                //????????????????????????????????????
                if (max_Score > Rt_U) {  // U2P??????????????? U.LCL??????
                    unprune++;
                    leaf_unprune++;
                    L.insert(store_leaf);
                }
                    //?????????????????????????????????
                else {
                    //cout<<"leaf "<<store_leaf<<"???n"<<node<<"?????????"<<endl;
                    prune++;
                    leaf_prune++;
                }
                //?????????????????????????????????????????????
                set<int> usr_withinLeaf;
                for (int term: Ukeys) {
                    for (int usr_id: GTree[node].inverted_list_u[term]) {   //??????????????????????????????
                        usr_withinLeaf.insert(usr_id);
                    }
                }

                for (int usr_id: usr_withinLeaf) {   //???????????????????????????????????????
                    //double u_lcl = 0;
                    double u_lcl_rt =0;
                    priority_queue<ResultCurrent> u_LCL = get_usrLCL_update_Plus(usr_id,
                                                                                 update_o_set, Qk, a, alpha, 100);  //??????????????????lcl
                    priority_queue<ResultCurrent> u_related_store;

                    if (u_LCL.size() < Qk) u_lcl_rt = 0.0;
                    else if (u_LCL.size() == Qk) u_lcl_rt = u_LCL.top().score;


                    //cout<<"u_lcl="<<u_lcl<<endl;
                    //??????
                    Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);
                    for (int store: stores) {  //????????????store_leaf??????store
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
                        //????????????????????????????????????
                        if (u_lcl_rt < gsk_score_Upper) {
                            //update lcl
                            u_related_store.push(ResultCurrent(store, gsk_score_Lower, gsk_score_Upper, rel, inf));


                            //cout<<"u"<<usr_id<<", gsk_score_Lower="<<gsk_score_Lower<<",u_lcl="<<u_lcl_rt<<",??????store="<<store<<endl;
                        } else {
                            usr_prune++;
                            //continue;   //??????????????????top-k??????
                        }
                    }

                    if (u_related_store.size() > 0) {
                        verification_User.insert(usr_id);

                        if (candidate_usr_related_store[usr_id].size() > 0) {  //???????????????u_lcl????????????????????????

                            //????????????????????????object???user??????
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

                        } else {  //????????????????????????
                            //priority_queue<ResultCurrent> u_related_Query = u_related_store;
                            candidate_usr_related_store[usr_id] = u_related_store;
                        }


                    } else {
                        //cout << "u" << usr_id << "????????????, as Lu??????,u_lcl=" << u_lcl << endl;
                    }
                    //getchar();
                }

                //cout<<"?????????????????????"<<node<<endl;
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
            //cout<<"??????usr node ??????:n"<<cn<<"??????"<<endl;
        }
        if(related_usr_node.size() ==0){
            cout<<"related_usr_node??????"<<endl;
            break;

        }
        //???????????????????????????????????????
        if(related_children.size()==1){
            cnt_++;
            if(cnt_==3){
                //cout<<"?????????????????????"<<endl;
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
                //cout<<"?????????"<<related_usr_node.size()<<endl;
                //getchar();

            }

        }
        //cout<<"???"<<level<<"????????????????????????????????????"<<endl;
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
    map<int, vector<int>> leaf_store_list;//<leaf_id, <store_id....>>  leaf ??????store?????????
    map<int, set<int>> leaf_verify_list;//<leaf_id, <usr_id....>>  leaf ?????????store???????????????????????????


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
        //cout<<"??????"<<i<<",leaf="<<leaf<<endl;
    }
    cout<<"?????????????????????"<<store_leafSet.size()<<"?????????????????????,GTree???"<<store_nodes.size()<<"??????????????????"<<endl;

    set<int> keyword_Set;
    for(int t: keywords)
        keyword_Set.insert(t);

    //??????partition????????????????????????
    set<int> check_in_partition;
    for(int s: stores){
        for(int u: POIs[s].check_ins){
            check_in_partition.insert(u);
        }
    }

    //???????????????????????????????????????
    set<int> totalUsrSet;
    for(int term:keyword_Set){
        for(int usr: invListOfUser[term]){
            totalUsrSet.insert(usr);
        }

    }


    //??????????????????????????? ???leafNodeInv??????????????????

    int prune=0; int unprune=0; int user_access=0; int usr_prune =0;
    int leaf_prune=0; int leaf_unprune=0;
    int may_prune = 0;
    //int topK_count = 0; //topk????????????
    // ???????????????????????????????????????????????????
    int t3=0; int root=0;



    //cout << " ,relate poi leaf size=" << related_poi_leafSet.size() << endl;
    //cout << " ,relate usr leaf size=" << related_usr_leafSet.size() << endl;

    vector<int> related_usr_node;   //????????????????????????

    map<int,vector<int>> related_poi_nodes;  //???????????????????????????

    // ?????????????????????????????????****??????****?????????????????????????????????, ??????query partition ??????

    vector<int> unodes = GTree[root].children;
    for(int unode: unodes){
        related_usr_node.push_back(unode);
        GTree[unode].ossociate_queryNodes = store_leafSet; //?????????????????????????????????
    }


    // ?????????????????????????????????****?????????***** ?????????????????????????????????
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
            //cout << "node" << node << "?????????" << endl;

            //cout<<"iteration"<<i<<":"<<endl;
            TreeNode un = getGIMTreeNodeData(node,OnlyU);

            set<int> inter_Key = obtain_itersection_jins(un.userUKeySet, keyword_Set);

            //??????????????????????????????
            if(!un.isleaf) {
                int usrNum = 0;
                int store_leaf = store_partition_id;
                if (false) {   //??????????????????  //x
                    //cout << "??????????????????" << node <<",usrNum="<<usrNum<<endl;

                    for(int term: keyword_Set){
                        for(int usr_leaf: GTree[node].term_usrLeaf_Map[term]){
                            related_children.insert(usr_leaf);
                            for(int store_leaf:GTree[node].ossociate_queryNodes){
                                GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                            }
                        }

                    }

                    prune++;

                    //cout<<"?????????????????? usr_node"<<node<<",prune"<<endl;   //??????????????????

                } else {  //?????????????????????????????????????????????
                    //?????? ????????????
                    double max_Score = -1;

                    //max_Score = getGSKScoreP2U_Upper_GivenKey(a, alpha, keywords, store_leaf, node);
                    max_Score = getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, check_in_partition, store_leaf, node); //?????????????????? InterKey, ?????????check-in?????????????????????

                    double Rt_U = 100;
                    //?????????????????????LCL

                    for(int term: inter_Key){
                        LCLResult lcl;
                        map<int,LCLResult>::iterator iter;
                        iter = GTree[node].termLCLDetail.find(term);
                        if(iter != GTree[node].termLCLDetail.end()){ //?????????cache
                            lcl = iter->second;
                        }else{
                            lcl = getUpperNode_termLCL_Disk(node,term, related_poi_nodes,Qk,a,alpha,max_Score);
                            GTree[node].termLCLDetail[term] = lcl;
                        }
                        double current_score = lcl.topscore;
                        Rt_U = min(Rt_U, current_score);
                    }


                    //cout<<"GTree[node].related_queryNodes[term_id].size="<<GTree[node].related_queryNodes[term_id].size()<<endl;
                    //?????????????????????????????????

                    if (max_Score > Rt_U) {  // ???????????????????????????????????????????????????????????????????????????
                        unprune++;

                        for(int term: inter_Key){
                            vector<int> usr_leafList = getUsrTermRelatedEntry(term,node);
                            if(usr_leafList.size()>0){
                                for (int child: usr_leafList) {
                                    //for (int child: textual_child){
                                    //??????????????????????????????lcl??????
                                    GTree[child].termParentLCLDetail[term] = GTree[node].termParentLCLDetail[term];
                                    //GTree[child].termParentFilter[term_id] = filter;
                                    related_children.insert(child);
                                    //cout<<"??????n"<<child<<endl;

                                }
                            }

                        }

                        //getchar();
                    }
                        //???????????????????????????
                    else {
                        prune++;
#ifdef PRUNELOG
                        cout << "??????" << node << "???store_leaf" << store_leaf << "????????? ???????????????????????????" << endl;
#endif
                    }

                    //cout<<"Rt_U="<<Rt_U<<endl;
                    //cout<<"get the LCL of u leaf"<<leaf_id<<endl;


                }//else???????????????
            }

            else if (un.isleaf) {  //?????????????????????????????????   onceleafsuoz
#ifdef TRACK
                cout << "access usr_leaf" << node << endl;
#endif

                //????????????store_leaf
                int store_leaf = store_partition_id;
                double max_Score =
                        getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, check_in_partition, store_leaf, node);



                MultiLCLResult lclResult_multi;
                vector<int> Ukeys;  //?????????
                for(int term_id: inter_Key)
                    Ukeys.push_back(term_id);

                map<set<int>,MultiLCLResult>::iterator iter;
                iter = GTree[node].cacheMultiLCLDetail.find(inter_Key);
                if (iter!=GTree[node].cacheMultiLCLDetail.end()){
                    lclResult_multi = iter->second;
                }

                    //???cache, ??????
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

                //?????? usr_lcl update_o_set
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
                //?????????????????????????????????????????????
                set<int> L;
                //????????????????????????????????????
                if (max_Score > Rt_U) {  // U2P??????????????? U.LCL??????
                    unprune++;
                    leaf_unprune++;
                    L.insert(store_leaf);
                }
                    //?????????????????????????????????
                else {
                    //cout<<"leaf "<<store_leaf<<"???n"<<node<<"?????????"<<endl;
                    prune++;
                    leaf_prune++;
                }
                //?????????????????????????????????????????????
                set<int> usr_withinLeaf;
                for (int term: Ukeys) {
                    vector<int> usr_list = getUsrTermRelatedEntry(term,node);
                    for (int usr_id: usr_list) {   //??????????????????????????????
                        usr_withinLeaf.insert(usr_id);
                    }
                }

                for (int usr_id: usr_withinLeaf) {   //???????????????????????????????????????
                    //double u_lcl = 0;
                    double u_lcl_rt =0;
                    priority_queue<ResultCurrent> u_LCL = get_usrLCL_update_Plus(usr_id,
                                                                                 update_o_set, Qk, a, alpha, 100);  //??????????????????lcl
                    priority_queue<ResultCurrent> u_related_store;

                    if (u_LCL.size() < Qk) u_lcl_rt = 0.0;
                    else if (u_LCL.size() == Qk) u_lcl_rt = u_LCL.top().score;


                    //cout<<"u_lcl="<<u_lcl<<endl;
                    //??????
                    Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);
                    for (int store: stores) {  //????????????store_leaf??????store
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
                        //????????????????????????????????????
                        if (u_lcl_rt < gsk_score_Upper) {
                            //update lcl
                            u_related_store.push(ResultCurrent(store, gsk_score_Lower, gsk_score_Upper, rel, inf));


                            //cout<<"u"<<usr_id<<", gsk_score_Lower="<<gsk_score_Lower<<",u_lcl="<<u_lcl_rt<<",??????store="<<store<<endl;
                        } else {
                            usr_prune++;
                            //continue;   //??????????????????top-k??????
                        }
                    }

                    if (u_related_store.size() > 0) {
                        verification_User.insert(usr_id);

                        if (candidate_usr_related_store[usr_id].size() > 0) {  //???????????????u_lcl????????????????????????

                            //????????????????????????object???user??????
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

                        } else {  //????????????????????????
                            //priority_queue<ResultCurrent> u_related_Query = u_related_store;
                            candidate_usr_related_store[usr_id] = u_related_store;
                        }


                    } else {
                        //cout << "u" << usr_id << "????????????, as Lu??????,u_lcl=" << u_lcl << endl;
                    }
                    //getchar();
                }

                //cout<<"?????????????????????"<<node<<endl;
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
            //cout<<"??????usr node ??????:n"<<cn<<"??????"<<endl;
        }
        if(related_usr_node.size() ==0){
            cout<<"related_usr_node??????"<<endl;
            break;

        }
        //???????????????????????????????????????
        if(related_children.size()==1){
            cnt_++;
            if(cnt_==3){
                //cout<<"?????????????????????"<<endl;
                //getchar();
                related_usr_node.clear();
                for(int node: related_children){
                    if(GTree[node].isleaf)
                        //related_usr_node.push_back(node);
                        related_usr_node.push_back(node);
                    else{
                        set<int> leaves;  //?????????set, ???????????????
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
                //cout<<"?????????"<<related_usr_node.size()<<endl;
                //getchar();

            }

        }
        //cout<<"???"<<level<<"????????????????????????????????????"<<endl;
        //cout<<"related_poi_node size in next iteration:"<<related_poi_nodes.size()<<endl;
        //cout<<"related_usr_node size in next iteration:"<<related_usr_node.size()<<endl;
        if(related_usr_node.size()==0) break;

    }//end while
    //

}


//2. ????????????
//????????????
BRkGSKQResults Verification_BRkGSKQ_memory(set<int> verification_User, map<int, priority_queue<ResultCurrent>> candidate_usr_related_store, vector<int> qKey, int Qk, float a, float alpha){
    //??????RkGSKQ?????????
    BRkGSKQResults result_Stores;
    set<int> candidate_usr = verification_User;//iter->second; //
    //cout<<"candidate_usr.size="<<candidate_usr.size()<<endl;
    map<int, vector<int>> usr_group;
    //1. ???maxscore??????????????????(OK)???2. ?????????????????????????????????
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
            //???????????????gsk??????
            double rel = rc.relevance; double inf = rc.influence;
            double gsk_score = getGSKScore_q2u(a,  alpha, POIs[rc.o_id], Users[usr_id], rel, inf);
            ResultLargerFisrt rlf(rc.o_id,gsk_score,rel, inf);
            Lu.push(rlf);

            //max_score = max(max_score,gsk_score);
        }
        CardinalityFisrt cf(usr_id,0.0,Lu);
        //user_queue.push(cf);
        user_list.push_back(cf);
        //???????????????????????????????????????
        /*for(int term: Users[usr_id].keywords){
            term_candidateUserInv[term].push_back(usr_id);
        }*/

    }
    //cout<<"?????????????????????"<<usr_group.size()<<"???"<<endl;

    // verifiying the TkGSKQ results for each users in candidate users, report query objects which havie user as top-k
    int useful =0; int usr_count =0;

    for(int i=0; i<user_list.size();i++){
        //topkSDijkstra_verify_usr(transformToQ(Users[0]), Qk, a, alpha, 1);
        CardinalityFisrt cf = user_list[i];
        usr_count++;
        int usr_id = cf.usr_id;
        priority_queue<ResultLargerFisrt> Lu = cf.Lu;
        //cout<<"???"<<usr_count<<"???????????????,??????"<<Lu.size()<<"???store,max_score="<<Lu.top().score<<endl;
        if(Lu.top().score > Users[usr_id].topkScore_current){


            double score_bound = Lu.top().score; //Users[usr_id].topkScore_current;
            double Rk_u = candidate_user_verify_single(usr_id,Qk,a, alpha,score_bound);

            //???????????????????????????topk score???
            //updateTopkScore_Current(a, Qk, current_elements,term_candidateUserInv);
            topK_count++;
            if(Lu.top().score > Rk_u){
                useful++;
                //???Lu????????????query object????????????
                while(!Lu.empty()){
                    ResultLargerFisrt rlf = Lu.top();
                    Lu.pop();
                    double rel = 0; double inf = 0;
                    double gsk_score = rlf.score; int query_object = rlf.o_id;
                    if(gsk_score > Rk_u){
                        ResultDetail rd(usr_id, 0, 0, 0, gsk_score, Rk_u);
                        result_Stores[query_object].push_back(rd);
                        //cout<<"u"<<usr_id<<"?????????o"<<query_object<<"???RkGSKQ??????????????????gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
                    }
                    else break;
                }
            }
            else{
                //cout<<"topk score = "<<Rk_u<<",max_score="<<max_score<<endl;
            }
        }
        else{
            //cout<<"????????????"<<endl;
        }
        //???????????????????????????
        Users[usr_id].verified = true;
    }
    //cout<<"useful ratio"<<useful*1.0/topK_count<<endl;
    //getchar();
    return  result_Stores;
}

//????????????(ok)
BRkGSKQResults Verification_BRkGSKQ_disk(set<int> verification_User, map<int, priority_queue<ResultCurrent>> candidate_usr_related_store, vector<int> qKey, int Qk, float a, float alpha){
    //??????RkGSKQ?????????
    BRkGSKQResults result_Stores;
    set<int> candidate_usr = verification_User;//iter->second; //
    //cout<<"candidate_usr.size="<<candidate_usr.size()<<endl;
    map<int, vector<int>> usr_group;
    //1. ???maxscore??????????????????(OK)???2. ?????????????????????????????????
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
            //???????????????gsk??????
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
        //???????????????????????????????????????
        /*for(int term: Users[usr_id].keywords){
            term_candidateUserInv[term].push_back(usr_id);
        }*/

    }
    //cout<<"?????????????????????"<<usr_group.size()<<"???"<<endl;

    // verifiying the TkGSKQ results for each users in candidate users, report query objects which havie user as top-k
    int useful =0; int usr_count =0;

    for(int i=0; i<user_list.size();i++){
        //topkSDijkstra_verify_usr(transformToQ(Users[0]), Qk, a, alpha, 1);
        CardinalityFisrt cf = user_list[i];
        usr_count++;
        int usr_id = cf.usr_id;
        priority_queue<ResultLargerFisrt> Lu = cf.Lu;
        //cout<<"???"<<usr_count<<"???????????????,??????"<<Lu.size()<<"???store,max_score="<<Lu.top().score<<endl;
        if(Lu.top().score > Users[usr_id].topkScore_current){


            double score_bound = Lu.top().score; //Users[usr_id].topkScore_current;
            double Rk_u = candidate_user_verify_single(usr_id,Qk,a, alpha,score_bound);
            //cout<<"u"<<usr_id<<" Rk_u="<<Rk_u<<endl;
            //???????????????????????????topk score???
            //updateTopkScore_Current(a, Qk, current_elements,term_candidateUserInv);
            topK_count++;
            if(Lu.top().score > Rk_u){
                useful++;
                //???Lu????????????query object????????????
                while(!Lu.empty()){
                    ResultLargerFisrt rlf = Lu.top();
                    Lu.pop();
                    double rel = 0; double inf = 0;
                    double gsk_score = rlf.score; int query_object = rlf.o_id;
                    if(gsk_score > Rk_u){
                        ResultDetail rd(usr_id, 0, 0, 0, gsk_score, Rk_u);
                        result_Stores[query_object].push_back(rd);
                        //cout<<"u"<<usr_id<<"?????????o"<<query_object<<"???RkGSKQ??????????????????gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
                    }
                    else break;
                }
            }
            else{
                //cout<<"topk score = "<<Rk_u<<",max_score="<<max_score<<endl;
            }
        }
        else{
            //cout<<"????????????"<<endl;
        }
        //???????????????????????????
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
    map<int, vector<int>> leaf_store_list;//<leaf_id, <store_id....>>  leaf ??????store?????????
    map<int, set<int>> leaf_verify_list;//<leaf_id, <usr_id....>>  leaf ?????????store???????????????????????????

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
        //cout<<"??????"<<i<<",leaf="<<leaf<<endl;
    }
    cout<<"?????????????????????"<<store_leafSet.size()<<"?????????????????????,GTree???"<<store_nodes.size()<<"??????????????????"<<endl;
    //getchar();

    set<int> keyword_Set;
    for(int t: keywords)
        keyword_Set.insert(t);

    set<int> checkin_usr_Set;


    vector<int> Leaf_node_User,Leaf_node_Poi;

    vector<ResultDetail> result_User;
    //??????????????????????????? ???leafNodeInv??????????????????

    //????????????usr????????????????????????????????????
    int prune=0; int unprune=0; int user_access=0; int usr_prune =0;
    int leaf_prune=0; int leaf_unprune=0;
    int may_prune = 0;
    int topK_count = 0; //topk????????????
    // ???????????????????????????????????????????????????
    set<int> related_usrLeafSet;
    int t3=0; int root=0;

    set<int> totalUsrSet;
    Filter_START
    int nodeAccess = 0;

    /*----------???????????????query cluster??????Group Filtering------------------*/
    set<int> verification_User;  //?????????????????????
    map<int, set<VerifyEntry>> verification_Map;    //???????????????   <leaf_id, <user_id....>>
    map<int, priority_queue<ResultCurrent>> candidate_usr_related_store;   //???????????????????????????????????????
    //map<int, vector<ResultDetail>> result_Stores; //<store_id, <usr_id....>>
    map<int, LCLResult> candidate_usr_LCL;     //???????????????LCL???
    cout<<"???"<<store_leafSet.size()<<"???query partition"<<endl; //getchar();
    TIME_TICK_START
    Filter_START
    map<int, map<int,bool>> usr_store_map;   //?????? user???object????????????????????????
    int p_count = 0;
    for(int store_leaf: store_leafSet){   //?????????leaf??????store????????????
        //cout<< store_leaf<<endl;
        p_count++;
        vector<int> stores_inpartition =  leaf_store_list[store_leaf];   //store_leaf
#ifndef DiskAccess
        Group_Filter_Baseline_memory(usr_store_map, store_leaf, stores_inpartition, keywords, Qk, a, alpha, verification_User, verification_Map, candidate_usr_related_store, candidate_usr_LCL);
#else
        Group_Filter_Baseline_disk(usr_store_map, store_leaf, stores_inpartition, keywords, Qk, a, alpha, verification_User, verification_Map, candidate_usr_related_store, candidate_usr_LCL);
#endif
        //PAUSE_START
        cout<<"query partition "<<store_leaf<<"???filter????????????,??????"<<(store_leafSet.size()-p_count)<<"???partition"<<endl;
        cout<<"????????????"<<verification_User.size()<<"???????????????"<<endl;

    }
    Filter_END
    //?????????candidate
    cout<<"?????? query partition ???filter ??????????????????"<<endl;
    cout<<"??????????????????"<<verification_User.size()<<"???????????????"<<endl;

    Refine_START

    //?????????query_partition???????????????candidate??????????????????
#ifndef DiskAccess
    BatchResults batch_query_results = Verification_BRkGSKQ_memory(verification_User, candidate_usr_related_store, keywords, Qk, a, alpha);
#else
    BatchResults batch_query_results = Verification_BRkGSKQ_disk(verification_User, candidate_usr_related_store, keywords, Qk, a, alpha);
#endif
    cout<<"??????????????????query??????"<<batch_query_results.size()<<endl;  //baseyan

    Refine_END


    //???????????????usr???????????????????????????
    //FilterResults results =  resultsAnalysis(result_User,keywords,checkin_usr,loc, Qk,a,alpha);
    //printBatchRkGSKQResults(stores, batch_query_results);
    cout<<"??????????????????"<<verification_User.size()<<"???????????????"<<endl;
    //cout<<"???????????????????????????????????????"<< totalUsrSet.size() <<endl;
    cout<<"topk?????????????????? =" << topK_count <<endl;
    cout<<"g-tree node prune ="<<prune<<" ,unprune="<<unprune<<",usr_prune="<<usr_prune<<endl;
    double ratio = (1.0*prune) /(prune+unprune);
    double ratio2 = (1.0*leaf_prune) /(leaf_prune+leaf_unprune);
    cout<<"???????????????= "<<ratio<<"??????????????????????????????="<<ratio2<<endl;
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
    map<int, vector<int>> leaf_store_list;//<leaf_id, <store_id....>>  leaf ??????store?????????
    map<int, set<int>> leaf_unionKeyword;//<leaf_id, <store_id....>>  leaf ??????store?????????
    map<int, set<int>> leaf_verify_list;//<leaf_id, <usr_id....>>  leaf ?????????store???????????????????????????

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
        //???poi???keyword????????????????????????unionKeyword???
        for(int keyword: p.keywords)
            leaf_unionKeyword[leaf].insert(keyword);
        //??????leaf???cluster?????????????????????????????????
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
        //cout<<"??????"<<i<<",leaf="<<leaf<<endl;
    }
    cout<<"?????????????????????"<<store_leafSet.size()<<"?????????????????????,GTree???"<<store_nodes.size()<<"??????????????????"<<endl;
    //getchar();

    set<int> keyword_Set;
    for(int t: keywords)
        keyword_Set.insert(t);

    set<int> checkin_usr_Set;


    vector<int> Leaf_node_User,Leaf_node_Poi;

    vector<ResultDetail> result_User;
    //??????????????????????????? ???leafNodeInv??????????????????

    //????????????usr????????????????????????????????????
    int prune=0; int unprune=0; int user_access=0; int usr_prune =0;
    int leaf_prune=0; int leaf_unprune=0;
    int may_prune = 0;
    int topK_count = 0; //topk????????????
    // ???????????????????????????????????????????????????
    set<int> related_usrLeafSet;
    int t3=0; int root=0;

    set<int> totalUsrSet;
    Filter_START
    int nodeAccess = 0;

    /*----------???????????????query cluster??????Group Filtering------------------*/
    set<int> verification_User;  //?????????????????????
    map<int, set<VerifyEntry>> verification_Map;    //???????????????   <leaf_id, <user_id....>>
    map<int, priority_queue<ResultCurrent>> candidate_usr_related_store;   //???????????????????????????????????????
    //map<int, vector<ResultDetail>> result_Stores; //<store_id, <usr_id....>>
    map<int, LCLResult> candidate_usr_LCL;     //???????????????LCL???
    cout<<"???"<<store_leafSet.size()<<"???query partition"<<endl; //getchar();
    TIME_TICK_START
    Filter_START
    map<int, map<int,bool>> usr_store_map;   //?????? user???object????????????????????????
    int p_count = 0;
    for(int store_leaf: store_leafSet){   //?????????leaf??????store????????????
        //cout<< store_leaf<<endl;
        p_count++;
        //??????partition??? unionKeyword
        keywords = transforSet2Vector(leaf_unionKeyword[store_leaf]);
        vector<int> stores_inpartition =  leaf_store_list[store_leaf];   //store_leaf
#ifndef DiskAccess
        Group_Filter_Baseline_memory(usr_store_map, store_leaf, stores_inpartition, keywords, Qk, a, alpha, verification_User, verification_Map, candidate_usr_related_store, candidate_usr_LCL);
#else
        Group_Filter_Baseline_disk(usr_store_map, store_leaf, stores_inpartition, keywords, Qk, a, alpha, verification_User, verification_Map, candidate_usr_related_store, candidate_usr_LCL);
#endif
        //PAUSE_START
        cout<<"query partition "<<store_leaf<<"???filter????????????,??????"<<(store_leafSet.size()-p_count)<<"???partition"<<endl;
        cout<<"????????????"<<verification_User.size()<<"???????????????"<<endl;

    }
    Filter_END
    //?????????candidate
    cout<<"?????? query partition ???filter ??????????????????"<<endl;
    cout<<"??????????????????"<<verification_User.size()<<"???????????????"<<endl;

    Refine_START

    //?????????query_partition???????????????candidate??????????????????
#ifndef DiskAccess
    BatchResults batch_query_results = Verification_BRkGSKQ_memory(verification_User, candidate_usr_related_store, keywords, Qk, a, alpha);
#else
    BatchResults batch_query_results = Verification_BRkGSKQ_disk(verification_User, candidate_usr_related_store, keywords, Qk, a, alpha);
#endif
    cout<<"??????????????????query??????"<<batch_query_results.size()<<endl;  //baseyan

    Refine_END


    //???????????????usr???????????????????????????
    //FilterResults results =  resultsAnalysis(result_User,keywords,checkin_usr,loc, Qk,a,alpha);
    //printBatchRkGSKQResults(stores, batch_query_results);
    cout<<"??????????????????"<<verification_User.size()<<"???????????????"<<endl;
    //cout<<"???????????????????????????????????????"<< totalUsrSet.size() <<endl;
    cout<<"topk?????????????????? =" << topK_count <<endl;
    cout<<"g-tree node prune ="<<prune<<" ,unprune="<<unprune<<",usr_prune="<<usr_prune<<endl;
    double ratio = (1.0*prune) /(prune+unprune);
    double ratio2 = (1.0*leaf_prune) /(leaf_prune+leaf_unprune);
    cout<<"???????????????= "<<ratio<<"??????????????????????????????="<<ratio2<<endl;
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
//1. ???????????????????????? Access GIM-Tree  Once

//????????????

void  Group_Filter_Once_memory(vector<int> stores, vector<int> keywords, int Qk, float a, float alpha, set<VerifyEntry>& verification_User, map<int,set<VerifyEntry>>& verification_Map, map<int, priority_queue<ResultCurrent>>& candidate_usr_related_store, map<int, LCLResult>& candidate_usr_LCL) {
    vector<int> qKey = keywords;
    set<int> store_leafSet;
    map<int, set<int>> partition_checkinMap;

    set<int> store_nodes;
    map<int, vector<int>> leaf_store_list;//<leaf_id, <store_id....>>  leaf ??????store?????????
    map<int, set<int>> leaf_verify_list;//<leaf_id, <usr_id....>>  leaf ?????????store???????????????????????????


    BatchResults result_Stores; //<store_id, <usr_id....>>
    int loc =0;

    //generate the query tree for candidate stores
    for(int store_id:stores){
        int leaf = Nodes[POIs[store_id].Nj].gtreepath.back();
        leaf_store_list[leaf].push_back(store_id);
        //?????????????????????????????????
        store_leafSet.insert(leaf);
        //???????????????partition???check-in??????
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
        //cout<<"??????"<<i<<",leaf="<<leaf<<endl;
    }
    cout<<"?????????????????????"<<store_leafSet.size()<<"?????????????????????,GTree???"<<store_nodes.size()<<"??????????????????"<<endl;

    TIME_TICK_START
    set<int> keyword_Set;
    for(int t: keywords)
        keyword_Set.insert(t);


    //??????????????????????????? ???leafNodeInv??????????????????

    int prune=0; int unprune=0; int user_access=0; int usr_prune =0;
    int leaf_prune=0; int leaf_unprune=0;
    int may_prune = 0;
    //int topK_count = 0; //topk????????????
    // ???????????????????????????????????????????????????
    int t3=0; int root=0;


    set<int> totalUsrSet;
    int nodeAccess = 0;
    //??????????????????????????????????????? ??????????????????
    vector<int> related_usr_node;
    //??????????????????????????????????????? ?????????????????????
    map<int,vector<int>> related_poi_nodes;

    for(int term:keyword_Set)
        for(int usr: invListOfUser[term])
            totalUsrSet.insert(usr);

    cout<<"????????????????????????????????? num="<<totalUsrSet.size();


    Filter_START

    //for(int term_id:keywords){
    // ?????????????????????????????????****??????***** ?????????????????????????????????

    vector<int> unodes = GTree[root].children; //get the child list of root whose text contains term_id
    for (int n : unodes) {
        //cout<<"node"<<n<<endl;
        related_usr_node.push_back(n);
        GTree[n].ossociate_queryNodes = store_leafSet; //?????????????????????????????????
    }
    // ?????????????????????????????????****?????????***** ?????????????????????????????????
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
    TIME_TICK_PRINT("??????????????????runtime")


    //verify the user node from upper to bottom
    int level = 1;  int cnt_ = 0;


    while (true) {  //daxunhuanonce
        //related_children.clear();
        if (related_usr_node.capacity() == 0) break;
        set<int>().swap(related_children);

        for (int node:related_usr_node) {


            set<int> inter_Key  = obtain_itersection_jins(GTree[node].userUKeySet, keyword_Set);
            if(!inter_Key.size()>0) continue;

            //??????????????????????????????
            if(!GTree[node].isleaf){

#ifdef TRACK
                cout << "access usr_node" << node << endl;
#endif
                int usrNum = 0;
                double global_max_Score = -1;
                //cout<<"non-leaf node???"<< node <<endl;
                if (false) {   //??????????????????  //x

                    for(int term: keyword_Set){
                        for(int usr_leaf: GTree[node].term_usrLeaf_Map[term]){
                            related_children.insert(usr_leaf);
                            for(int store_leaf:GTree[node].ossociate_queryNodes){
                                GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                            }
                        }

                    }
                    prune++;
                    //cout << "??????????????????" << node <<",usrNum="<<usrNum<<endl;

                } else {  //?????????????????????????????????????????????
                    //global upper bound  computation
                    vector<double> max_score_list;
                    for (int store_leaf: GTree[node].ossociate_queryNodes) {
                        //upper bound computation
                        double  max_Score =
                                getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, partition_checkinMap[store_leaf], store_leaf, node); //?????????????????? InterKey, ?????????check-in?????????????????????

                        global_max_Score = max( max_Score, global_max_Score);

                        max_score_list.push_back(max_Score);

                    }
                    //obtaining lcl from cache or lcl computation on the fly
                    double Rt_U = 100;
                    for(int term: inter_Key){
                        LCLResult lcl;
                        map<int,LCLResult>::iterator iter;
                        iter = GTree[node].termLCLDetail.find(term);
                        if(iter != GTree[node].termLCLDetail.end()){ //?????????cache
                            lcl = iter->second;
                        }else{
                            lcl = getUpperNode_termLCL(node,term, related_poi_nodes,Qk,a,alpha,global_max_Score);
                            GTree[node].termLCLDetail[term] = lcl;
                        }
                        double current_score = lcl.topscore;
                        Rt_U = min(Rt_U, current_score);
                    }
                    // group pruning judgement
                    if (global_max_Score > Rt_U) {  // ???????????????????????????????????????????????????????????????????????????
                        unprune++;
                        int _idx = 0;
                        for (int store_leaf: GTree[node].ossociate_queryNodes) {
                            //upper bound computation
                            double max_score = max_score_list[_idx];
                            _idx++;

                            if(Rt_U>max_score)   //??????pruning ???????????????
                                continue;

                            for(int term: inter_Key){
                                if(GTree[node].inverted_list_u.count(term)>0){
                                    for (int child: GTree[node].inverted_list_u[term]) {

                                        //??????????????????????????????lcl??????
                                        GTree[child].termParentLCLDetail[term] = GTree[node].termParentLCLDetail[term];
                                        //????????????????????????
                                        related_children.insert(child);
                                        //?????? ???????????? ???????????? ??????
                                        GTree[child].ossociate_queryNodes.insert(store_leaf);

                                    }
                                }
                            }

                        }


                    }
                        //???????????????????????????
                    else {
                        prune++;
#ifdef PRUNELOG
                        cout << "??????" << node << "????????? store??????????????? ??????????????????????????????" << endl;
#endif

                    }


                }//else???????????????
            }

            else if (GTree[node].isleaf) {  //?????????????????????????????????
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
                            getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, partition_checkinMap[store_leaf], store_leaf, node); //?????????????????? InterKey, ?????????check-in?????????????????????

                    global_max_Score = max( max_Score, global_max_Score);

                    max_score_list.push_back(max_Score);

                }

                //lcl computation, obtaining the lcl-r_t

                map<set<int>,MultiLCLResult>::iterator iter;
                iter = GTree[node].cacheMultiLCLDetail.find(inter_Key);
                if (iter!=GTree[node].cacheMultiLCLDetail.end()){
                    lclResult_multi = iter->second;
                }

                else{  //???cache  ??????


                    //lclResult = getLeafNode_TopkSDijkstra(node, term_id, Qk, a, alpha, max_Score);
#ifdef LV
                    lclResult_multi = getLeafNode_LCL_Dijkstra_multi(node, Ukeys, Qk, a, alpha, global_max_Score);
#endif
#ifndef LV
                    lclResult_multi = getLeafNode_LCL_bottom2up_multi(node, Ukeys, Qk, a, alpha, global_max_Score);
#endif


                    GTree[node].cacheMultiLCLDetail[inter_Key] = lclResult_multi;
                }

                //?????? usr_lcl update_o_set
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

                //??????lcl-r_t?????????
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
                    for (int store_leaf: GTree[node].ossociate_queryNodes) {  //??????partition???max_score_list??????????????????upper bound
                        //upper bound computation
                        double  max_Score = max_score_list[_idx];
                        _idx++;

                        if( max_Score > Rt_U)
                            L.insert(store_leaf);
                        else{
                            continue;   //???store_leaf ?????????
                        }

                    }

                    //????????????????????????????????????
                    set<int> usr_withinLeaf;
                    for (int term: inter_Key) {
                        for (int usr_id: GTree[node].inverted_list_u[term]) {   //??????????????????????????????
                            usr_withinLeaf.insert(usr_id);
                        }
                    }

                    //?????????????????????????????????????????????
                    for (int usr_id: usr_withinLeaf) {   //??????????????????

                        double u_lcl_rt = 0.0;


                        priority_queue<ResultCurrent> u_related_store;

                        if(update_o_set.size()>0) {
                            //obtain the lcl list of a user from leaf lcl
                            priority_queue<ResultCurrent> u_LCL = get_usrLCL_update_Plus(
                                    usr_id, update_o_set, Qk, a, alpha, global_max_Score);  //??????????????????lcl

                            if (u_LCL.size() < Qk) u_lcl_rt = 0.0;
                            else if (u_LCL.size() == Qk) u_lcl_rt = u_LCL.top().score;

                        }

                        Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);
                        for (int store_leaf: L) {
                            for (int store: GTree[store_leaf].stores) {  //????????????store_leaf??????store
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
                                                                       inf));  //??????relevance???influence ????????????
                                    if (u_related_store.size() > Qk) {
                                        u_related_store.pop();
                                        u_lcl_rt = max(u_lcl_rt, u_related_store.top().score);
                                    }
                                    //cout<<"u"<<usr_id<<", gsk_score_Lower="<<gsk_score_Lower<<",u_lcl="<<u_lcl<<",??????store="<<store<<endl;
                                    Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);

                                } else {
                                    usr_prune++;
                                    //continue;   //??????????????????top-k??????
                                }

                            }
                        }
                        if (u_related_store.size() > 0) {
                            VerifyEntry ve = VerifyEntry(usr_id,-1,u_lcl_rt);
                            verification_User.insert(ve);
                            if (candidate_usr_related_store[usr_id].size() > 0) {  //???????????????u_lcl????????????????????????
                                priority_queue<ResultCurrent> u_related_Query_pre = candidate_usr_related_store[usr_id];
                                while (!u_related_store.empty()) {
                                    ResultCurrent rc = u_related_store.top();
                                    u_related_store.pop();
                                    u_related_Query_pre.push(rc);
                                    if (u_related_Query_pre.size() > Qk)
                                        u_related_Query_pre.pop();
                                }
                            } else {  //????????????????????????
                                priority_queue<ResultCurrent> u_related_Query = u_related_store;
                                candidate_usr_related_store[usr_id] = u_related_Query;
                            }
                            //Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current,
                            //candidate_usr_related_store[usr_id].top().score);

#ifdef PRINTINFO
                            cout << "?????????????????? u" << usr_id << ", u_lcl=" << u_lcl_rt << endl;
#endif

                            //getchar();
                        } else {
                            //cout << "u" << usr_id << "???????????????, as Lu??????,u_lcl=" << u_lcl << endl;
                        }
                        //getchar();
                    }


                }
                    //????????????????????????????????????
                else {
#ifdef PRUNELOG
                    cout<<"leaf "<<store_leaf<<"???n"<<node<<"?????????"<<endl;
#endif
                    prune++;
                    leaf_prune++;
                }

                //cout<<"?????????????????????"<<node<<endl;
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
            //cout<<"related_usr_node??????"<<endl;
            break;

        }
        //???????????????????????????????????????
        if(related_children.size()==1){
            cnt_++;
            if(cnt_==3){
                //cout<<"?????????????????????"<<endl;
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
                            //??????????????????
                            for (int store_leaf: GTree[node].ossociate_queryNodes) {
                                GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                            }
                        }

                    }


                }
                //cout<<"?????????"<<related_usr_node.size()<<endl;
                //getchar();

            }

        }
        //cout<<"???"<<level<<"????????????????????????????????????"<<endl;
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
    map<int, vector<int>> leaf_store_list;//<leaf_id, <store_id....>>  leaf ??????store?????????
    map<int, set<int>> leaf_verify_list;//<leaf_id, <usr_id....>>  leaf ?????????store???????????????????????????


    BatchResults result_Stores; //<store_id, <usr_id....>>
    int loc =0;

    //generate the query tree for candidate stores
    for(int store_id:stores){
        POI p = getPOIFromO2UOrgLeafData(store_id);
        int leaf = Nodes[p.Nj].gtreepath.back();
        leaf_store_list[leaf].push_back(store_id);
        //?????????????????????????????????
        store_leafSet.insert(leaf);
        //???????????????partition???check-in??????
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
        //cout<<"??????"<<i<<",leaf="<<leaf<<endl;
    }
    cout<<"?????????????????????"<<store_leafSet.size()<<"?????????????????????,GTree???"<<store_nodes.size()<<"??????????????????"<<endl;

    TIME_TICK_START
    set<int> keyword_Set;
    for(int t: keywords)
        keyword_Set.insert(t);


    //??????????????????????????? ???leafNodeInv??????????????????

    int upper_prune=0; int upper_unprune=0;
    int leaf_prune=0; int leaf_unprune=0;
    int aggressive_prune =0;
    int user_access=0; int usr_prune =0;
    int may_prune = 0;
    //int topK_count = 0; //topk????????????
    // ???????????????????????????????????????????????????
    int t3=0; int root=0;


    set<int> totalUsrSet;
    int nodeAccess = 0;
    //??????????????????????????????????????? ??????????????????
    vector<int> related_usr_node;
    //??????????????????????????????????????? ?????????????????????
    map<int,vector<int>> related_poi_nodes;

    /*for(int term:keyword_Set)
        for(int usr: invListOfUser[term])
            totalUsrSet.insert(usr);

    cout<<"????????????????????????????????? num="<<totalUsrSet.size();*/


    Filter_START


    // ?????????????????????????????????****??????***** ?????????????????????????????????
    vector<int> unodes = GTree[root].children; //get the child list of root whose text contains term_id
    for (int n : unodes) {
        //cout<<"node"<<n<<endl;
        related_usr_node.push_back(n);
        GTree[n].ossociate_queryNodes = store_leafSet; //?????????????????????????????????
    }
    // ?????????????????????????????????****?????????***** ?????????????????????????????????
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
    TIME_TICK_PRINT("??????????????????runtime")


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

            //??????????????????????????????
            if(!user_node.isleaf){

                UPPER_START
                int usrNum = 0;
                double global_max_Score = -1;
#ifdef TRACK
                cout<<"non-leaf node???"<< node <<endl;
#endif
                for(int cover_key: inter_Key){

                }

                if (false) {   //??????????????????  //x

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
                    //cout << "??????????????????" << node <<",usrNum="<<usrNum<<endl;

                } else {  //?????????????????????????????????????????????
                    //global upper bound  computation
                    vector<double> max_score_list;
                    for (int store_leaf: GTree[node].ossociate_queryNodes) {
                        //upper bound computation
                        double  max_Score =
                                getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, partition_checkinMap[store_leaf], store_leaf, node); //?????????????????? InterKey, ?????????check-in?????????????????????

                        global_max_Score = max( max_Score, global_max_Score);

                        max_score_list.push_back(max_Score);

                    }
                    //obtaining lcl from cache or lcl computation on the fly
                    double Rt_U = 100;
                    for(int term: inter_Key){
                        LCLResult lcl;
                        map<int,LCLResult>::iterator iter;
                        iter = GTree[node].termLCLDetail.find(term);
                        if(iter != GTree[node].termLCLDetail.end()){ //?????????cache
                            lcl = iter->second;
                        }else{
                            lcl = getUpperNode_termLCL_Disk(node,term, related_poi_nodes,Qk,a,alpha,global_max_Score);
                            GTree[node].termLCLDetail[term] = lcl;
                        }
                        double current_score = lcl.topscore;
                        Rt_U = min(Rt_U, current_score);
                    }
                    // group pruning judgement
                    if (global_max_Score > Rt_U) {  // ???????????????????????????????????????????????????????????????????????????
                        upper_unprune++;
                        int _idx = 0;
                        for (int store_leaf: GTree[node].ossociate_queryNodes) {
                            //upper bound computation
                            double max_score = max_score_list[_idx];
                            _idx++;

                            if(Rt_U>max_score)   //??????pruning ???????????????
                                continue;

                            for(int term: inter_Key){
                                vector<int> child_entry_u = getUsrTermRelatedEntry(term,node);
                                if(child_entry_u.size()>0){
                                    for (int child: child_entry_u) {

                                        //??????????????????????????????lcl??????
                                        GTree[child].termParentLCLDetail[term] = GTree[node].termParentLCLDetail[term];
                                        //????????????????????????
                                        related_children.insert(child);
                                        //?????? ???????????? ???????????? ??????
                                        GTree[child].ossociate_queryNodes.insert(store_leaf);

                                    }
                                }
                            }

                        }


                    }
                        //???????????????????????????
                    else {
                        upper_prune++;
#ifdef PRUNELOG
                        cout << "??????" << node << "????????? store??????????????? ??????????????????????????????" << endl;
#endif

                    }


                }//else???????????????
                UPPER_PAUSE
            }

            else if (user_node.isleaf) {  //?????????????????????????????????
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
                            getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, partition_checkinMap[store_leaf], store_leaf, node); //?????????????????? InterKey, ?????????check-in?????????????????????

                    global_max_Score = max( max_Score, global_max_Score);

                    max_score_list.push_back(max_Score);

                }

                //lcl computation, obtaining the lcl-r_t

                map<set<int>,MultiLCLResult>::iterator iter;
                iter = GTree[node].cacheMultiLCLDetail.find(inter_Key);
                if (iter!=GTree[node].cacheMultiLCLDetail.end()){
                    lclResult_multi = iter->second;
                }

                else{  //???cache  ??????


                    //lclResult = getLeafNode_TopkSDijkstra(node, term_id, Qk, a, alpha, max_Score);
#ifdef LV
                    lclResult_multi = getLeafNode_LCL_Dijkstra_multi(node, Ukeys, Qk, a, alpha, global_max_Score);
#endif
#ifndef LV
                    lclResult_multi = getLeafNode_LCL_bottom2up_multi(node, Ukeys, Qk, a, alpha, global_max_Score);
#endif


                    GTree[node].cacheMultiLCLDetail[inter_Key] = lclResult_multi;
                }

                //?????? usr_lcl update_o_set
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

                //??????lcl-r_t?????????
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
                    for (int store_leaf: GTree[node].ossociate_queryNodes) {  //??????partition???max_score_list??????????????????upper bound
                        //upper bound computation
                        double  max_Score = max_score_list[_idx];
                        _idx++;

                        if( max_Score > Rt_U)
                            L.insert(store_leaf);
                        else{
                            continue;   //???store_leaf ?????????
                        }

                    }

                    //????????????????????????????????????
                    set<int> usr_withinLeaf;
                    for (int term: inter_Key) {
                        vector<int> usr_list = getUsrTermRelatedEntry(term,node);
                        for (int usr_id: usr_list) {   //??????????????????????????????
                            usr_withinLeaf.insert(usr_id);
                        }
                    }

                    //?????????????????????????????????????????????
                    for (int usr_id: usr_withinLeaf) {   //??????????????????

                        double u_lcl_rt = 0.0;


                        priority_queue<ResultCurrent> u_related_store;

                        if(update_o_set.size()>0) {
                            //obtain the lcl list of a user from leaf lcl
                            priority_queue<ResultCurrent> u_LCL = get_usrLCL_update_Plus(
                                    usr_id, update_o_set, Qk, a, alpha, global_max_Score);  //??????????????????lcl

                            if (u_LCL.size() < Qk) u_lcl_rt = 0.0;
                            else if (u_LCL.size() == Qk) u_lcl_rt = u_LCL.top().score;

                        }

                        Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);
                        for (int store_leaf: L) {
                            for (int store: GTree[store_leaf].stores) {  //????????????store_leaf??????store
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
                                                                       inf));  //??????relevance???influence ????????????
                                    if (u_related_store.size() > Qk) {
                                        u_related_store.pop();
                                        u_lcl_rt = max(u_lcl_rt, u_related_store.top().score);
                                    }
                                    //cout<<"u"<<usr_id<<", gsk_score_Lower="<<gsk_score_Lower<<",u_lcl="<<u_lcl<<",??????store="<<store<<endl;
                                    Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);

                                } else {
                                    usr_prune++;
                                    //continue;   //??????????????????top-k??????
                                }

                            }
                        }
                        if (u_related_store.size() > 0) {
                            VerifyEntry ve = VerifyEntry(usr_id,-1,u_lcl_rt);
                            verification_User.insert(ve);
                            if (candidate_usr_related_store[usr_id].size() > 0) {  //???????????????u_lcl????????????????????????
                                priority_queue<ResultCurrent> u_related_Query_pre = candidate_usr_related_store[usr_id];
                                while (!u_related_store.empty()) {
                                    ResultCurrent rc = u_related_store.top();
                                    u_related_store.pop();
                                    u_related_Query_pre.push(rc);
                                    if (u_related_Query_pre.size() > Qk)
                                        u_related_Query_pre.pop();
                                }
                            } else {  //????????????????????????
                                priority_queue<ResultCurrent> u_related_Query = u_related_store;
                                candidate_usr_related_store[usr_id] = u_related_Query;
                            }


#ifdef PRINTINFO
                            cout << "?????????????????? u" << usr_id << ", u_lcl=" << u_lcl_rt << endl;
#endif

                            //getchar();
                        } else {
                            //cout << "u" << usr_id << "???????????????, as Lu??????,u_lcl=" << u_lcl << endl;
                        }
                        //getchar();
                    }


                }
                    //????????????????????????????????????
                else {
#ifdef PRUNELOG
                    cout<<"leaf "<<store_leaf<<"???n"<<node<<"?????????"<<endl;
#endif
                    //prune++;
                    leaf_prune++;
                }

                //cout<<"?????????????????????"<<node<<endl;
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
            //cout<<"related_usr_node??????"<<endl;
            break;

        }
        //???????????????????????????????????????
        if(related_children.size()==1){
            cnt_++;
            if(cnt_==3){
                //cout<<"?????????????????????"<<endl;
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
                            //??????????????????
                            for (int store_leaf: GTree[node].ossociate_queryNodes) {
                                GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                            }
                        }

                    }


                }
                //cout<<"?????????"<<related_usr_node.size()<<endl;
                //getchar();

            }

        }
        //cout<<"???"<<level<<"????????????????????????????????????"<<endl;
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
    map<int, vector<int>> leaf_store_list;//<leaf_id, <store_id....>>  leaf ??????store?????????
    map<int, set<int>> leaf_verify_list;//<leaf_id, <usr_id....>>  leaf ?????????store???????????????????????????


    BatchResults result_Stores; //<store_id, <usr_id....>>
    int loc =0;

    //generate the query tree for candidate stores
    for(int store_id:stores){
        POI p = getPOIFromO2UOrgLeafData(store_id);
        int leaf = Nodes[p.Nj].gtreepath.back();
        leaf_store_list[leaf].push_back(store_id);
        //?????????????????????????????????
        store_leafSet.insert(leaf);
        //???????????????partition???check-in??????
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
        //cout<<"??????"<<i<<",leaf="<<leaf<<endl;
    }
    cout<<"?????????????????????"<<store_leafSet.size()<<"?????????????????????,GTree???"<<store_nodes.size()<<"??????????????????"<<endl;

    TIME_TICK_START
    set<int> keyword_Set;
    for(int t: keywords){  //
        keyword_Set.insert(t);
    }



    //??????????????????????????? ???leafNodeInv??????????????????

    int upper_prune=0; int upper_unprune=0;
    int leaf_prune=0; int leaf_unprune=0;
    int aggressive_prune =0;
    int user_access=0; int usr_prune =0;
    int may_prune = 0;
    //int topK_count = 0; //topk????????????
    // ???????????????????????????????????????????????????
    int t3=0; int root=0;


    set<int> totalUsrSet;
    int nodeAccess = 0;
    //??????????????????????????????????????? ??????????????????
    vector<int> related_usr_node;
    //??????????????????????????????????????? ?????????????????????
    map<int,vector<int>> related_poi_nodes;

    /*for(int term:keyword_Set)
        for(int usr: invListOfUser[term])
            totalUsrSet.insert(usr);

    cout<<"????????????????????????????????? num="<<totalUsrSet.size();*/


    Filter_START


    // ?????????????????????????????????****??????***** ?????????????????????????????????
    vector<int> unodes = GTree[root].children; //get the child list of root whose text contains term_id
    for (int n : unodes) {
        related_usr_node.push_back(n);
        GTree[n].ossociate_queryNodes = store_leafSet; //?????????????????????????????????
        /*printf("GTree[%d].ossociate_queryNodes\n:",n);
        printSetElements(GTree[n].ossociate_queryNodes);*/
    }
    // ?????????????????????????????????****?????????***** ?????????????????????????????????
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
    TIME_TICK_PRINT("??????????????????runtime")


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

            //??????????????????????????????
            if(!user_node.isleaf){

                UPPER_START
                int usrNum = 0;
                double global_max_Score = -1;
#ifdef TRACK
                cout<<"non-leaf node???"<< node <<endl;
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

                if (aggressive_flag) {   //??????????????????  //aggressive prune //user_capacity<4

                    for(int usr_leaf: usr_leaf_direct){
                        related_children.insert(usr_leaf);
                        for(int store_leaf:GTree[node].ossociate_queryNodes){
                            GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                        }

                    }
                    aggressive_prune++;
                    //cout << "??????????????????" << node <<",usrNum="<<usrNum<<endl;
#ifdef TRACK
                    cout << "??????????????????" << node <<",usrNum="<<usrNum<<endl;
#endif

                } else {  //?????????????????????????????????????????????
#ifdef TRACK
                    cout << "?????????????????????????????????????????????"<<endl;
#endif
                    //global upper bound  computation
                    vector<double> max_score_list;

                    /*printf("GTree[node]: n%d 's ossociate_queryNodes:\n",node);
                    printSetElements(GTree[node].ossociate_queryNodes);
                    getchar();*/

                    for (int store_leaf: GTree[node].ossociate_queryNodes) {
                        //upper bound computation
                        double  max_Score =
                                getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, partition_checkinMap[store_leaf], store_leaf, node); //?????????????????? InterKey, ?????????check-in?????????????????????
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
                        if(iter != GTree[node].termLCLDetail.end()){ //?????????cache
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
                    if (global_max_Score > Rt_U) {  // ???????????????????????????????????????????????????????????????????????????
                        upper_unprune++;

                        int _idx = 0;
                        for (int store_leaf: GTree[node].ossociate_queryNodes) {
                            //upper bound computation
                            double max_score = max_score_list[_idx];
                            _idx++;

                            if(Rt_U>max_score)   //??????pruning ???????????????
                                continue;

                            for(int term: inter_Key){
                                vector<int> child_entry_u = getUsrTermRelatedEntry(term,node);
                                if(child_entry_u.size()>0){
                                    for (int child: child_entry_u) {

                                        //??????????????????????????????lcl??????
                                        GTree[child].termParentLCLDetail[term] = GTree[node].termParentLCLDetail[term];
                                        //????????????????????????
                                        related_children.insert(child);
                                        //?????? ???????????? ???????????? ??????
                                        GTree[child].ossociate_queryNodes.insert(store_leaf);

                                    }
                                }
                            }

                        }


                    }
                        //???????????????????????????
                    else {
                        upper_prune++;
#ifdef PRUNELOG
                        cout << "??????" << node << "????????? store??????????????? ??????????????????????????????" << endl;
#endif

                    }


                }//else???????????????
                UPPER_PAUSE
            }

            else if (user_node.isleaf) {  //?????????????????????????????????
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
                            getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, partition_checkinMap[store_leaf], store_leaf, node); //?????????????????? InterKey, ?????????check-in?????????????????????

                    global_max_Score = max( max_Score, global_max_Score);

                    max_score_list.push_back(max_Score);

                }

                //lcl computation, obtaining the lcl-r_t

                map<set<int>,MultiLCLResult>::iterator iter;
                iter = GTree[node].cacheMultiLCLDetail.find(inter_Key);
                if (iter!=GTree[node].cacheMultiLCLDetail.end()){
                    lclResult_multi = iter->second;
                }

                else{  //???cache  ??????


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

                //?????? usr_lcl update_o_set
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

                //??????lcl-r_t?????????
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

                    for (int store_leaf: GTree[node].ossociate_queryNodes) {  //??????partition???max_score_list??????????????????upper bound
                        //upper bound computation
                        double  max_Score = max_score_list[_idx];
                        _idx++;

                        if( max_Score > Rt_U)
                            L.insert(store_leaf);
                        else{
                            continue;   //???store_leaf ?????????
                        }

                    }

                    //????????????????????????????????????
                    set<int> usr_withinLeaf;
                    for (int term: inter_Key) {
                        vector<int> usr_list = getUsrTermRelatedEntry(term,node);
                        for (int usr_id: usr_list) {   //??????????????????????????????
                            usr_withinLeaf.insert(usr_id);
                        }
                    }

                    //?????????????????????????????????????????????
                    for (int usr_id: usr_withinLeaf) {   //??????????????????

                        double u_lcl_rt = 0.0;


                        priority_queue<ResultCurrent> u_related_store;

                        if(update_o_set.size()>0) {
                            //obtain the lcl list of a user from leaf lcl
                            priority_queue<ResultCurrent> u_LCL = get_usrLCL_update_Plus(
                                    usr_id, update_o_set, Qk, a, alpha, global_max_Score);  //??????????????????lcl

                            if (u_LCL.size() < Qk) u_lcl_rt = 0.0;
                            else if (u_LCL.size() == Qk) u_lcl_rt = u_LCL.top().score;

                        }

                        Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);
                        for (int store_leaf: L) {
                            for (int store: GTree[store_leaf].stores) {  //????????????store_leaf??????store
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
                                                                       inf));  //??????relevance???influence ????????????
                                    if (u_related_store.size() > Qk) {
                                        u_related_store.pop();
                                        u_lcl_rt = max(u_lcl_rt, u_related_store.top().score);
                                    }
                                    //cout<<"u"<<usr_id<<", gsk_score_Lower="<<gsk_score_Lower<<",u_lcl="<<u_lcl<<",??????store="<<store<<endl;
                                    Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);

                                } else {
                                    usr_prune++;
                                    //continue;   //??????????????????top-k??????
                                }

                            }
                        }
                        if (u_related_store.size() > 0) {
                            VerifyEntry ve = VerifyEntry(usr_id,-1,u_lcl_rt);
                            verification_User.insert(ve);
                            if (candidate_usr_related_store[usr_id].size() > 0) {  //???????????????u_lcl????????????????????????
                                priority_queue<ResultCurrent> u_related_Query_pre = candidate_usr_related_store[usr_id];
                                while (!u_related_store.empty()) {
                                    ResultCurrent rc = u_related_store.top();
                                    u_related_store.pop();
                                    u_related_Query_pre.push(rc);
                                    if (u_related_Query_pre.size() > Qk)
                                        u_related_Query_pre.pop();
                                }
                            } else {  //????????????????????????
                                priority_queue<ResultCurrent> u_related_Query = u_related_store;
                                candidate_usr_related_store[usr_id] = u_related_Query;
                            }


#ifdef PRINTINFO
                            cout << "?????????????????? u" << usr_id << ", u_lcl=" << u_lcl_rt << endl;
#endif

                            //getchar();
                        } else {
                            //cout << "u" << usr_id << "???????????????, as Lu??????,u_lcl=" << u_lcl << endl;
                        }
                        //getchar();
                    }


                }
                    //????????????????????????????????????
                else {
#ifdef LEAF_PRUNELOG_
                    cout<<"uleaf "<<node<<"?????????????????????????????????"<<endl;
#endif
                    //prune++;
                    leaf_prune++;
                }

                //cout<<"?????????????????????"<<node<<endl;
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
            //cout<<"related_usr_node??????"<<endl;
            break;

        }
        //???????????????????????????????????????
        if(related_children.size()==1){
            cnt_++;
            if(cnt_==3){
                //cout<<"?????????????????????"<<endl;
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
                            //??????????????????
                            for (int store_leaf: GTree[node].ossociate_queryNodes) {
                                GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                            }
                        }

                    }


                }
                //cout<<"?????????"<<related_usr_node.size()<<endl;
                //getchar();

            }

        }
        //cout<<"???"<<level<<"????????????????????????????????????"<<endl;
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
    map<int, vector<int>> leaf_store_list;//<leaf_id, <store_id....>>  leaf ??????store?????????
    map<int, set<int>> leaf_verify_list;//<leaf_id, <usr_id....>>  leaf ?????????store???????????????????????????
    unordered_map<int,unordered_map<int,bool>> upMask;

    BatchResults result_Stores; //<store_id, <usr_id....>>
    int loc =0;

    //generate the query tree for candidate stores
    for(int store_id:stores){
        POI p = getPOIFromO2UOrgLeafData(store_id);
        vector<int> p_keywords = p.keywords;
        for(int term: p_keywords){
            vector<int> users_related = getTermUserInvList(term);
            if(users_related.size()<=uterm_SPARCITY_VALUE){//????????????????????????????????????????????????????????????????????????
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
                                                                               inf));  //??????relevance???influence ????????????
                        if (candidate_usr_related_store[usr_id].size() > Qk) {
                            candidate_usr_related_store[usr_id].pop();
                        }

                    }
                    else{
                        priority_queue<ResultCurrent> u_related_store;
                        u_related_store.push(ResultCurrent(store_id, gsk_score_Lower, gsk_score_Upper, rel,
                                                           inf));  //??????relevance???influence ????????????
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
        //?????????????????????????????????
        store_leafSet.insert(leaf);
        //???????????????partition???check-in??????
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
        //cout<<"??????"<<i<<",leaf="<<leaf<<endl;
    }
    cout<<"?????????????????????"<<store_leafSet.size()<<"?????????????????????,GTree???"<<store_nodes.size()<<"??????????????????"<<endl;

    TIME_TICK_START
    set<int> keyword_Set;
    for(int t: keywords){  //?????????????????????????????????
        vector<int> users_related = getTermUserInvList(t);
        if(users_related.size()>uterm_SPARCITY_VALUE){
            keyword_Set.insert(t);
        }
    }
    cout<<"???????????????filter?????????keywords:"; printSetElements(keyword_Set);



    //??????????????????????????? ???leafNodeInv??????????????????

    int upper_prune=0; int upper_unprune=0;
    int leaf_prune=0; int leaf_unprune=0;
    int aggressive_prune =0;
    int user_access=0; int usr_prune =0;
    int may_prune = 0;
    //int topK_count = 0; //topk????????????
    // ???????????????????????????????????????????????????
    int t3=0; int root=0;


    set<int> totalUsrSet;
    int nodeAccess = 0;
    //??????????????????????????????????????? ??????????????????
    vector<int> related_usr_node;
    //??????????????????????????????????????? ?????????????????????
    map<int,vector<int>> related_poi_nodes;

    /*for(int term:keyword_Set)
        for(int usr: invListOfUser[term])
            totalUsrSet.insert(usr);

    cout<<"????????????????????????????????? num="<<totalUsrSet.size();*/


    Filter_START


    // ?????????????????????????????????****??????***** ?????????????????????????????????
    vector<int> unodes = GTree[root].children; //get the child list of root whose text contains term_id
    for (int n : unodes) {
        related_usr_node.push_back(n);
        GTree[n].ossociate_queryNodes = store_leafSet; //?????????????????????????????????
        /*printf("GTree[%d].ossociate_queryNodes\n:",n);
        printSetElements(GTree[n].ossociate_queryNodes);*/
    }
    // ?????????????????????????????????****?????????***** ?????????????????????????????????
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
    TIME_TICK_PRINT("??????????????????runtime")


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

            //??????????????????????????????
            if(!user_node.isleaf){

                UPPER_START
                int usrNum = 0;
                double global_max_Score = -1;
#ifdef TRACK
                cout<<"non-leaf node???"<< node <<endl;
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

                if (aggressive_flag) {   //??????????????????  //aggressive prune //user_capacity<4

                    for(int usr_leaf: usr_leaf_direct){
                        related_children.insert(usr_leaf);
                        for(int store_leaf:GTree[node].ossociate_queryNodes){
                            GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                        }

                    }
                    aggressive_prune++;
                    //cout << "??????????????????" << node <<",usrNum="<<usrNum<<endl;
#ifdef TRACK
                    cout << "??????????????????" << node <<",usrNum="<<usrNum<<endl;
#endif

                } else {  //?????????????????????????????????????????????
#ifdef TRACK
                    cout << "?????????????????????????????????????????????"<<endl;
#endif
                    //global upper bound  computation
                    vector<double> max_score_list;

                    /*printf("GTree[node]: n%d 's ossociate_queryNodes:\n",node);
                    printSetElements(GTree[node].ossociate_queryNodes);
                    getchar();*/

                    for (int store_leaf: GTree[node].ossociate_queryNodes) {
                        //upper bound computation
                        double  max_Score =
                                getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, partition_checkinMap[store_leaf], store_leaf, node); //?????????????????? InterKey, ?????????check-in?????????????????????
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
                        if(iter != GTree[node].termLCLDetail.end()){ //?????????cache
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
                    if (global_max_Score > Rt_U) {  // ???????????????????????????????????????????????????????????????????????????
                        upper_unprune++;

                        int _idx = 0;
                        for (int store_leaf: GTree[node].ossociate_queryNodes) {
                            //upper bound computation
                            double max_score = max_score_list[_idx];
                            _idx++;

                            if(Rt_U>max_score)   //??????pruning ???????????????
                                continue;

                            for(int term: inter_Key){
                                vector<int> child_entry_u = getUsrTermRelatedEntry(term,node);
                                if(child_entry_u.size()>0){
                                    for (int child: child_entry_u) {

                                        //??????????????????????????????lcl??????
                                        GTree[child].termParentLCLDetail[term] = GTree[node].termParentLCLDetail[term];
                                        //????????????????????????
                                        related_children.insert(child);
                                        //?????? ???????????? ???????????? ??????
                                        GTree[child].ossociate_queryNodes.insert(store_leaf);

                                    }
                                }
                            }

                        }


                    }
                        //???????????????????????????
                    else {
                        upper_prune++;
#ifdef PRUNELOG
                        cout << "??????" << node << "????????? store??????????????? ??????????????????????????????" << endl;
#endif

                    }


                }//else???????????????
                UPPER_PAUSE
            }

            else if (user_node.isleaf) {  //?????????????????????????????????
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
                            getGSKScoreP2U_Upper_InterKey(a, alpha, inter_Key, partition_checkinMap[store_leaf], store_leaf, node); //?????????????????? InterKey, ?????????check-in?????????????????????

                    global_max_Score = max( max_Score, global_max_Score);

                    max_score_list.push_back(max_Score);

                }

                //lcl computation, obtaining the lcl-r_t

                map<set<int>,MultiLCLResult>::iterator iter;
                iter = GTree[node].cacheMultiLCLDetail.find(inter_Key);
                if (iter!=GTree[node].cacheMultiLCLDetail.end()){
                    lclResult_multi = iter->second;
                }

                else{  //???cache  ??????


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

                //?????? usr_lcl update_o_set
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

                //??????lcl-r_t?????????
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

                    for (int store_leaf: GTree[node].ossociate_queryNodes) {  //??????partition???max_score_list??????????????????upper bound
                        //upper bound computation
                        double  max_Score = max_score_list[_idx];
                        _idx++;

                        if( max_Score > Rt_U)
                            L.insert(store_leaf);
                        else{
                            continue;   //???store_leaf ?????????
                        }

                    }

                    //????????????????????????????????????
                    set<int> usr_withinLeaf;
                    for (int term: inter_Key) {
                        vector<int> usr_list = getUsrTermRelatedEntry(term,node);
                        for (int usr_id: usr_list) {   //??????????????????????????????
                            usr_withinLeaf.insert(usr_id);
                        }
                    }

                    //?????????????????????????????????????????????
                    for (int usr_id: usr_withinLeaf) {   //??????????????????

                        double u_lcl_rt = 0.0;


                        priority_queue<ResultCurrent> u_related_store;

                        if(update_o_set.size()>0) {
                            //obtain the lcl list of a user from leaf lcl
                            priority_queue<ResultCurrent> u_LCL = get_usrLCL_update_Plus(
                                    usr_id, update_o_set, Qk, a, alpha, global_max_Score);  //??????????????????lcl

                            if (u_LCL.size() < Qk) u_lcl_rt = 0.0;
                            else if (u_LCL.size() == Qk) u_lcl_rt = u_LCL.top().score;

                        }

                        Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);
                        for (int store_leaf: L) {
                            for (int store: GTree[store_leaf].stores) {  //????????????store_leaf??????store
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
                                                                       inf));  //??????relevance???influence ????????????
                                    if (u_related_store.size() > Qk) {
                                        u_related_store.pop();
                                        u_lcl_rt = max(u_lcl_rt, u_related_store.top().score);
                                    }
                                    //cout<<"u"<<usr_id<<", gsk_score_Lower="<<gsk_score_Lower<<",u_lcl="<<u_lcl<<",??????store="<<store<<endl;
                                    Users[usr_id].topkScore_current = max(Users[usr_id].topkScore_current, u_lcl_rt);

                                } else {
                                    usr_prune++;
                                    //continue;   //??????????????????top-k??????
                                }

                            }
                        }
                        if (u_related_store.size() > 0) {
                            VerifyEntry ve = VerifyEntry(usr_id,-1,u_lcl_rt);
                            verification_User.insert(ve);
                            if (candidate_usr_related_store[usr_id].size() > 0) {  //???????????????u_lcl????????????????????????
                                priority_queue<ResultCurrent> u_related_Query_pre = candidate_usr_related_store[usr_id];
                                while (!u_related_store.empty()) {
                                    ResultCurrent rc = u_related_store.top();
                                    u_related_store.pop();
                                    u_related_Query_pre.push(rc);
                                    if (u_related_Query_pre.size() > Qk)
                                        u_related_Query_pre.pop();
                                }
                            } else {  //????????????????????????
                                priority_queue<ResultCurrent> u_related_Query = u_related_store;
                                candidate_usr_related_store[usr_id] = u_related_Query;
                            }


#ifdef PRINTINFO
                            cout << "?????????????????? u" << usr_id << ", u_lcl=" << u_lcl_rt << endl;
#endif

                            //getchar();
                        } else {
                            //cout << "u" << usr_id << "???????????????, as Lu??????,u_lcl=" << u_lcl << endl;
                        }
                        //getchar();
                    }


                }
                    //????????????????????????????????????
                else {
#ifdef LEAF_PRUNELOG_
                    cout<<"uleaf "<<node<<"?????????????????????????????????"<<endl;
#endif
                    //prune++;
                    leaf_prune++;
                }

                //cout<<"?????????????????????"<<node<<endl;
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
            //cout<<"related_usr_node??????"<<endl;
            break;

        }
        //???????????????????????????????????????
        if(related_children.size()==1){
            cnt_++;
            if(cnt_==3){
                //cout<<"?????????????????????"<<endl;
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
                            //??????????????????
                            for (int store_leaf: GTree[node].ossociate_queryNodes) {
                                GTree[usr_leaf].ossociate_queryNodes.insert(store_leaf);
                            }
                        }

                    }


                }
                //cout<<"?????????"<<related_usr_node.size()<<endl;
                //getchar();

            }

        }
        //cout<<"???"<<level<<"????????????????????????????????????"<<endl;
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


//2. ????????????????????? OPT
//????????????
BRkGSKQResults Group_Verification_Optimized_memory(set<VerifyEntry> verification_User, map<int, priority_queue<ResultCurrent>> candidate_usr_related_store, vector<int> qKey, int Qk, float a, float alpha){
    //??????RkGSKQ?????????
    BRkGSKQResults result_Stores;
    set<VerifyEntry> candidate_usr = verification_User;//iter->second; //
    //cout<<"candidate_usr.size="<<candidate_usr.size()<<endl;
    map<int, vector<int>> usr_group;
    //1. ???maxscore??????????????????(OK)???2. ?????????????????????????????????
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
            //???????????????gsk??????
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

        //???????????????????????????????????????
        for(int term: Users[usr_id].keywords){
            term_candidateUserInv[term].push_back(usr_id);
        }

    }

    // verifiying the TkGSKQ results for each users in candidate users, report query objects which havie user as top-k
    int useful =0; int usr_count =0;

    //????????????

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
        //cout<<"???"<<usr_count<<"???????????????,??????"<<Lu.size()<<"???store,max_score="<<Lu.top().score<<endl;
        double Rk_u=0.0;

#ifdef TRACK
        cout<<"??????u"<<usr_id<<endl;

#endif

        if(Lu.top().score > Users[usr_id].topkScore_current){

//jordan
            double score_bound = usr_lcl_rt;
            //cout<<"u"<<usr_id<<", rt="<<usr_lcl_rt<<endl;


            Rk_u = candidate_user_verify_for_batch(usr_id,Qk,a, alpha,Lu.top().score,score_bound);

            //???????????????????????????topk score???
            //updateTopkScore_Current(a, Qk, current_elements,term_candidateUserInv);
            topK_count++;
            if(Lu.top().score > Rk_u){
                useful++;
                //???Lu????????????query object????????????
                while(!Lu.empty()){
                    ResultLargerFisrt rlf = Lu.top();
                    Lu.pop();
                    double rel = 0; double inf = 0;
                    double gsk_score = rlf.score; int query_object = rlf.o_id;
                    if(gsk_score > Rk_u){
                        ResultDetail rd(usr_id, 0, 0, 0, gsk_score, Rk_u);
                        result_Stores[query_object].push_back(rd);
                        //cout<<"????????????????????????gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
                    }
                    else break;
                }
            }
            else{
                //cout<<"topk score = "<<Rk_u<<",max_score="<<max_score<<endl;
            }
        }
        else{
            //cout<<"????????????"<<endl;
        }

        //???????????????????????????
        Users[usr_id].verified = true;
    }
    //cout<<"useful ratio"<<useful*1.0/topK_count<<endl;
    //getchar();
    return  result_Stores;
}

//????????????
BRkGSKQResults Group_Verification_Optimized_disk(set<VerifyEntry> verification_User, map<int, priority_queue<ResultCurrent>> candidate_usr_related_store, vector<int> qKey, int Qk, float a, float alpha){
    //??????RkGSKQ?????????
    BRkGSKQResults result_Stores;
    set<VerifyEntry> candidate_usr = verification_User;//iter->second; //
    //cout<<"candidate_usr.size="<<candidate_usr.size()<<endl;
    map<int, vector<int>> usr_group;
    //1. ???maxscore??????????????????(OK)???2. ?????????????????????????????????
    vector<CardinalityFisrt> user_list;
    map<int, vector<int>> term_candidateUserInv;
    for(VerifyEntry ve: candidate_usr){
        int usr_id = ve.u_id;
        /* if(7209==usr_id){
             cout<<"verify find u7209???"<<endl;
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
            //???????????????gsk??????
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

        //???????????????????????????????????????
        for(int term: u.keywords){
            term_candidateUserInv[term].push_back(usr_id);
        }

    }

    // verifiying the TkGSKQ results for each users in candidate users, report query objects which havie user as top-k
    int usr_count =0; int potential =0;

#ifdef  Verification_DEBUG
    cout<<"?????????"<<user_list.size()<<"???candidate user!"<<endl;
#endif

    //????????????

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
        //cout<<"???"<<usr_count<<"???????????????,??????"<<Lu.size()<<"???store,max_score="<<Lu.top().score<<endl;
        double Rk_u=0.0;
#ifdef  Verification_DEBUG
        cout<<"?????????????????????u"<<usr_id<<endl;
#endif
/*#ifdef TRACK
        cout<<"???????????????u"<<usr_id<<endl;
#endif*/
        if(Lu.top().score > Users[usr_id].topkScore_current){

//jordan
#ifdef TRACK
            cout<<"???????????????u"<<usr_id<<endl;
#endif
            double score_bound = usr_lcl_rt;
            //cout<<"u"<<usr_id<<", rt="<<usr_lcl_rt<<endl;


            Rk_u = candidate_user_verify_for_batch(usr_id,Qk,a, alpha,Lu.top().score,score_bound);
            //Rk_u = candidate_user_verify_for_OGPM(usr_id,Qk,a, alpha,10000,0);


            //???????????????????????????topk score???
            //updateTopkScore_Current(a, Qk, current_elements,term_candidateUserInv);
            topK_count++;
            if(Lu.top().score > Rk_u){
#ifdef TRACK
                cout<<"??????????????????u"<<usr_id<<",Rk_u="<<Rk_u<<endl;

#endif

                potential++;
                //???Lu????????????query object????????????
                while(!Lu.empty()){
                    ResultLargerFisrt rlf = Lu.top();
                    Lu.pop();
                    double rel = 0; double inf = 0;
                    double gsk_score = rlf.score; int query_object = rlf.o_id;
                    if(gsk_score > Rk_u){
                        ResultDetail rd(usr_id, 0, 0, 0, gsk_score, Rk_u);
                        result_Stores[query_object].push_back(rd);
#ifdef TRACK
                        cout<<",??????o"<<query_object<<"????????????????????????gsk_score="<<gsk_score<<",Rk_u="<<Rk_u<<endl;
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
            //cout<<"????????????"<<endl;
        }

        //???????????????????????????
        Users[usr_id].verified = true;

#ifdef TRACK
        //cout<<"????????????user"<<usr_id<<",Rk_u="<<Rk_u<<endl;

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

    /*----------????????????Group Filtering phase------------------*/
    set<VerifyEntry> verification_User;  //?????????????????????
    map<int,set<VerifyEntry>> verification_Map;    //???????????????   <leaf_id, <user_id....>>
    map<int, priority_queue<ResultCurrent>> candidate_usr_related_store;   //???????????????????????????????????????
    map<int, LCLResult> candidate_usr_LCL;     //???????????????LCL???
    BatchResults batch_query_results;
    //??????????????????
    Filter_START
#ifndef DiskAccess
    Group_Filter_Once_memory(stores, keywords, Qk, a, alpha, verification_User, verification_Map, candidate_usr_related_store, candidate_usr_LCL);
#else
    Group_Filter_Once_disk(stores, keywords, Qk, a, alpha, verification_User, verification_Map, candidate_usr_related_store, candidate_usr_LCL);

#endif
    Filter_END

    //??????????????????
    Refine_START
#ifndef DiskAccess
    batch_query_results = Group_Verification_Optimized_memory(verification_User, candidate_usr_related_store, keywords, Qk,a,alpha);
#else
    batch_query_results = Group_Verification_Optimized_disk(verification_User, candidate_usr_related_store, keywords, Qk,a,alpha);
#endif
    batch_end = clock();
    Refine_END
    //???????????????usr???????????????????????????
    //FilterResults results =  resultsAnalysis(result_User,keywords,checkin_usr,loc, Qk,a,alpha);
    printBatchRkGSKQResults(stores, batch_query_results);
    cout<<"filter??????????????????"<<verification_User.size()<<"???????????????"<<endl;
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