//
// Created by jins on 10/31/19.
//

#ifndef MAXINFRKGSKQ_1_1_FUNCTIONAL_H

#include <stdio.h>
#include <vector>
#include <algorithm>
#include <queue>
#include <iostream>
#include <map>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <string.h>
#include <cmath>
#include "config.h"
#include "bichromatic.h"
#include "road_network.h"
#include "gim_tree.h"

using namespace std;




class UsrScore{
public:
    int u_id;
    double score;
public:
    UsrScore(int id, double s):u_id(id),score(s){

    }

};


//查询构造函数
vector<int> extractQueryKeyowrds(int poi_id,  int capacity){
    vector<int> keypool;
    cout<<"query keyword pool:";
#ifdef	 DiskAccess
    //int address = getPOIAddrByIdxFile(poi_id);
    POI poi = getPOIFromO2UOrgLeafData(poi_id);
    keypool =  poi.keywords;
    printElements(keypool);
#else
    keypool = POIs[poi_id].keywords;
    printElements(keypool);
#endif

    set<int> _set; vector<int> target;
    for(int term: keypool)
        _set.insert(term);
    for(int term2: _set){
        target.push_back(term2);
        if(target.size()==capacity) break;
    }
    cout<<"exact query keyword:";
    printElements(target);

    return 	target;
}


vector<int> extractQueryGroupKeywords(vector<int>& query_objects){
    cout<<"extractQueryGroupKeyowrds:";
    POI poi; set<int> QG_KeywordSet;
    for(int qo:query_objects){
        vector<int> curr_keywords;
#ifdef	 DiskAccess
        //int address = getPOIAddrByIdxFile(poi_id);
        poi = getPOIFromO2UOrgLeafData(qo);
        curr_keywords =  poi.keywords;
#else
        int poi_id = qo;
        curr_keywords = POIs[poi_id].keywords;
#endif
        //printElements(curr_keywords);

        for(int term: curr_keywords)
            QG_KeywordSet.insert(term);
        /*for(int term2: _set){
            target.push_back(term2);
            if(target.size()==capacity) break;
        }*/

    }

    vector<int> QG_Keywords;
    for(int term: QG_KeywordSet)
        QG_Keywords.push_back(term);
    cout<<"exact Query Group KeywordSet:";
    printElements(QG_Keywords);


    return 	QG_Keywords;
}

vector<int> extractQueryGroupKeywords_refine(vector<int>& query_objects){
    cout<<"extractQueryGroupKeyowrds:";
    POI poi; set<int> QG_KeywordSet;
    for(int qo:query_objects){
        vector<int> curr_keywords;
#ifdef	 DiskAccess
        //int address = getPOIAddrByIdxFile(poi_id);
        poi = getPOIFromO2UOrgLeafData(qo);
        curr_keywords =  poi.keywords;
#else
        int poi_id = qo;
        curr_keywords = POIs[poi_id].keywords;
#endif
        //printElements(curr_keywords);

        for(int term: curr_keywords){
            if(term>200) continue;
            QG_KeywordSet.insert(term);
        }



    }

    vector<int> QG_Keywords;
    for(int term: QG_KeywordSet)
        QG_Keywords.push_back(term);
    cout<<"exact Query Group KeywordSet:";
    printElements(QG_Keywords);


    return 	QG_Keywords;
}

vector<int> extractQueryGroupKeywords_refine_KQConstraint(vector<int>& query_objects, int KQ_num){
    cout<<"extractQueryGroupKeyowrds, need keyowrd_size="<<KQ_num<<endl;
    POI poi; set<int> QG_KeywordSet;
    for(int qo:query_objects){
        vector<int> curr_keywords;
#ifdef	 DiskAccess
        //int address = getPOIAddrByIdxFile(poi_id);
        poi = getPOIFromO2UOrgLeafData(qo);
        curr_keywords =  poi.keywords;
#else
        int poi_id = qo;
        curr_keywords = POIs[poi_id].keywords;
#endif
        //printElements(curr_keywords);

        for(int term: curr_keywords){
            if(term>200) continue;
            QG_KeywordSet.insert(term);
        }

        if(QG_KeywordSet.size() > KQ_num) break;

    }

    vector<int> QG_Keywords;
    for(int term: QG_KeywordSet)
        QG_Keywords.push_back(term);
    cout<<"exact Query Group KeywordSet:";
    printElements(QG_Keywords);


    return 	QG_Keywords;
}


vector<int> extractAllQueryKeywords(vector<int>& query_objects){

    POI poi; set<int> QG_KeywordSet;
    for(int qo:query_objects){
        //vector<int> curr_keywords;

        poi = getPOIFromO2UOrgLeafData(qo);
        //printElements(poi.keywords);
        for(int term: poi.keywords){
            QG_KeywordSet.insert(term);
        }


    }

    vector<int> QG_Keywords;
    for(int term: QG_KeywordSet)
        QG_Keywords.push_back(term);
    cout<<"The keywords in all query objects:";
    printElements(QG_Keywords);

    return 	QG_Keywords;
}



vector<int> extractCheckins(int poi_id){
    vector<int> check_ins;
    //cout<<"check_in records:";
#ifdef	 DiskAccess

    POI poi = getPOIFromO2UOrgLeafData(poi_id);
    check_ins =  poi.check_ins;
    //printElements(check_ins);

#else
    //check_ins = poiCheckInList[poi_id];//POIs[poi_id].check_ins;
    //printElements(check_ins);

    check_ins = POIs[poi_id].check_ins;
    printElements(check_ins);
#endif
    return  check_ins;

}

int extractLocation(int poi_id){
    int loc;
    cout<<"loc, Nj:";
#ifdef	 DiskAccess
    POI poi = getPOIFromO2UOrgLeafData(poi_id);
    loc = poi.Nj;
    cout<<loc<<endl;

#else
    loc = POIs[poi_id].Nj;
    cout<<loc;
#endif
    return  loc;

}
//查询构造函数





// 以下为 Geo-social textual 评分计算函数

double calcSocialScore(int u_id, int p_id){
    vector<int> friends = friendshipMap[u_id];

    int optimal_Checkin = 0;
    int optimal_poi = 0;
    double friend_on_Pi = 0.0;

    //求u朋友中在 p_id 处签到的人数
    for(int usr:friends){
        if(userCheckInMap[usr].count(p_id)){
            friend_on_Pi += 1;
        }
    }

    if(friend_on_Pi == 0.0)
        return 0;


    if(usrMaxSS.count(u_id)){
        //if(Users[u_id].maxPopular !=-1.0){    //判断当键u_id下，是否存在值 （重要）
        optimal_Checkin = usrMaxSS[u_id];
        //optimal_Checkin = Users[u_id].maxPopular;

        //optimal_poi = usrOptimalPOI[u_id];

        double influence = friend_on_Pi / optimal_Checkin;
        /*if(influence == 1){
            cout<<"optimal_Checkin="<<optimal_Checkin<<", optimal_poi="<<optimal_poi<<", ";
            cout<<"getSocialScore("<<u_id<<","<<p_id<<")="<<influence<<endl;
        }*/
        return influence;

    }


    map<int,int> FriendPOIMap;  //(poi_id, numbers)

    //分析用户朋友的兴趣点check-in情况
    for(int usr:friends){
        if(userCheckInIDList.count(usr) == 0) continue;
        vector<int> user_Checkin = userCheckInIDList[usr];
        for(int poi: user_Checkin){
            if(userCheckInMap[usr].count(poi)){   //判断当键为poi时，是否存在值 （重要）
                if(FriendPOIMap.count(poi)){
                    int sum = FriendPOIMap[poi];
                    FriendPOIMap[poi] = sum + 1;

                } else
                    FriendPOIMap[poi] = 1;

                if(optimal_Checkin < FriendPOIMap[poi]){
                    optimal_Checkin = FriendPOIMap[poi];
                    optimal_poi = poi;
                }
            }
        }
    }

    usrMaxSS[u_id] = optimal_Checkin;


    double influence = friend_on_Pi / optimal_Checkin;

    return influence;



}


double reCalcuSocialScore(int u_id, vector<int> checkin_usrs){
    vector<int> friends = friendshipMap[u_id];
    map<int,int> FriendPOIMap;
    int optimal_Checkin = 0;
    int optimal_poi = 0;
    double friend_on_Pi = 0.0;
    //求poi的签到用户中为user朋友的总人数
    for(int usr:friends){
        for(int u2: checkin_usrs){
            if(usr==u2)
                //friend_on_Pi += 2000;
                friend_on_Pi += 1;
        }
    }
    if(friend_on_Pi==0.0) return 0;

    //求所有poi的签到用户中包含user的朋友的最大人数
    if(usrMaxSS.count(u_id)){    //利用之前计算的结果
        //if(Users[u_id].maxPopular != -1.0){
        optimal_Checkin = usrMaxSS[u_id];

        double influence = friend_on_Pi / optimal_Checkin;

        return influence;

    }
    //之前未计算过，则计算
    for(int usr:friends){
        if(userCheckInIDList.count(usr) == 0) continue;
        vector<int> user_Checkin = userCheckInIDList[usr];
        for(int poi: user_Checkin){
            if(userCheckInMap[usr].count(poi)){
                int checkin_Count = 1;
                if(FriendPOIMap.count(poi) == 0){
                    int sum = FriendPOIMap[poi];
                    FriendPOIMap[poi] = sum + checkin_Count;

                } else
                    FriendPOIMap[poi] = checkin_Count;

                if(optimal_Checkin < FriendPOIMap[poi]){
                    optimal_Checkin = FriendPOIMap[poi];
                    optimal_poi = poi;
                }
            }
        }
    }
    usrMaxSS[u_id] = optimal_Checkin;
    //Users[u_id].maxPopular = optimal_Checkin;

    //usrOptimalPOI[u_id] = optimal_poi;

    double influence = friend_on_Pi / optimal_Checkin;


    return influence;

}


double getGSKScore_loc(int a, double alpha, int loc, int usr_id, vector<int> qKeys, vector<int> check_usr){
    set<int> uKeyset; set<int> qKeyset;
    //for(int term : Users[usr_id].keywords) uKeyset.insert(term);
    //for(int term :qKeys)  qKeyset.insert(term);
    //double relevance = textRelevance(uKeyset,qKeyset);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    double relevance = textRelevance(Users[usr_id].keywords,qKeys);
    if(relevance == 0){
        return  0.0;
    }
    double influence;
    double social_Score = 0;
    if(a==0||alpha==0.0){
        influence = 1;
    }
    else{
        if(check_usr.size()>0)
            social_Score =  reCalcuSocialScore(usr_id,check_usr);

        if(social_Score == 0)
            influence =1;
        else
            influence =1 + pow(a,social_Score);


    }
    double distance = getDistance(loc,Users[usr_id]);  //距离有问题！

    double social_textual_score = (alpha*influence) + (1.0-alpha)*relevance;
    double gsk_score = social_textual_score / (1+distance);
    return gsk_score;

}


double getGSKScore_o2u_givenKey(int a, double alpha, POI poi, User user, vector<int> qKeys, vector<int> check_usr){
    set<int> uKeyset; set<int> qKeyset;
    //for(int term : Users[usr_id].keywords) uKeyset.insert(term);
    //for(int term :qKeys)  qKeyset.insert(term);
    //double relevance = textRelevance(uKeyset,qKeyset);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    int usr_id = user.id;
    double relevance = textRelevance(user.keywords,qKeys);
    if(relevance == 0){
        return  0.0;
    }
    double influence;
    double social_Score = 0;
    if(a==0||alpha==0.0){
        influence = 1;
    }
    else{
        if(check_usr.size()>0)
            social_Score =  reCalcuSocialScore(usr_id,check_usr);

        if(social_Score == 0)
            influence =1;
        else
            influence =1 + pow(a,social_Score);


    }
    //double distance = getDistance(loc,Users[usr_id]);  //距离有问题！
    double distance = usrToPOIDistance(user, poi);



    double social_textual_score = (alpha*influence) + (1.0-alpha)*relevance;
    double gsk_score = social_textual_score / (1+distance);
    return gsk_score;

}



double getGSKScore_o2u(int a, double alpha, POI poi, User user){   //poi与user的GSK评分

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
    double gsk_score = social_textual_score / (1+distance);
    //if(poi.id== 57&& user.id ==625 ) cout<<"getGSKScore_o2u，发现u57,p625,distance="<<distance<<",relevance="<<relevance<<"influence="<<influence<<",score="<<gsk_score<<endl;
    return gsk_score;

}

double getGSKScore_o2u_distComputed(int a, double alpha, POI poi, User user,double distance){   //poi与user的GSK评分

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


    //double distance = usrToPOIDistance(user, poi);   //这里要注意区分 getDistance(vertex, vertex)
    double social_textual_score = alpha*influence + (1.0-alpha)*relevance;
    double gsk_score = social_textual_score / (1+distance);
    //if(poi.id== 57&& user.id ==625 ) cout<<"getGSKScore_o2u，发现u57,p625,distance="<<distance<<",relevance="<<relevance<<"influence="<<influence<<",score="<<gsk_score<<endl;
    return gsk_score;

}




double getGSKScore_o2u(int a, double alpha, POI poi, User user, double relevance, double influence){   //poi与user的GSK评分

    double distance = usrToPOIDistance(user, poi);   //这里要注意区分 getDistance(vertex, vertex)

    double social_textual_score = alpha*influence + (1.0-alpha)*relevance;
    double gsk_score = social_textual_score / (1+distance);
    //if(poi.id== 57&& user.id ==625 ) cout<<"getGSKScore_o2u，发现u57,p625,distance="<<distance<<",relevance="<<relevance<<"influence="<<influence<<",score="<<gsk_score<<endl;
    return gsk_score;

}


double getGSKScore_o2u_givenkey(int a, double alpha, POI poi, User user, vector<int> keys){   //poi与user的GSK评分

    set<int> uKeyset; set<int> qKeyset;
    //for(int term : Users[usr_id].keywords) uKeyset.insert(term);
    //for(int term :qKeys)  qKeyset.insert(term);
    //double relevance = textRelevance(uKeyset,qKeyset);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    vector<int> qKeys =  keys;
    vector<int> check_usr = poiCheckInIDList[poi.id];
    double gsk_score = 0.0;

    double relevance = textRelevance(user.keywords,qKeys);
    if(relevance == 0.0)
        return 0.0;

    double influence;
    double social_Score = 0;
    if(a==0||alpha==0.0){
        influence = 1;
    }else{
        if(check_usr.size()>0)
            social_Score =  reCalcuSocialScore(user.id,check_usr);  //gazi
            //social_Score = calcSocialScore(usr.id,poi_id);

        if(social_Score == 0)
            influence =1;
        else
            influence =1 + pow(a,social_Score);
    }


    double distance = usrToPOIDistance(user, poi);   //这里要注意区分 getDistance(vertex, vertex)
    double social_textual_score = alpha*influence + (1.0-alpha)*relevance;
    gsk_score = social_textual_score / (1+distance);
    //if(poi.id== 57&& user.id ==625 ) cout<<"getGSKScore_o2u，发现u57,p625,distance="<<distance<<",relevance="<<relevance<<"influence="<<influence<<",score="<<gsk_score<<endl;
    return gsk_score;

}




double getGSKScoreO2U_Lower(int a, double alpha, int o_node, int u_node, int term){


    double relevance = 0;

    relevance = tfIDF_term(term) ;  //考虑 u.keywords与 p.qKey 至少在关键词term上相关

    if(relevance > 0){

        //cout<<"read MaxDisBound["<<o_node<<"]["<< u_node <<"]=";
        double max_distance = MaxDisBound[o_node][u_node];//100000;//getMinMaxDistanceBetweenNode(leaf, usr_node).max;
#ifdef  LCL_LOG
        cout<<"max_dist[n"<<o_node<<",n"<<u_node<<"]="<<max_distance<<endl;
#endif

        double min_relevance;
        //double relevance = textRelevance(usr.keywords,qKey);
        min_relevance = relevance;    //这里可以再提高
        //cout<<"UkeySizeMax="<<UkeySizeMax<<endl;

        double min_influence = 1.0;

        double min_score = 0;
        double social_textual_minscore = alpha*min_influence + (1.0-alpha)*min_relevance;
        min_score = social_textual_minscore / (1.0+max_distance);
        //score =0;
        return min_score;

    }
    else
        return  0;

}


double getGSKScoreO2U_Upper(int a, double alpha, set<int> o_node_keySet, int o_node, int u_node){

    //提取usr_leaf 最大 keyword size

    //double max_relevance = 1.0 / minKeySize;
    double max_relevance = 0 ;
#ifndef DiskAccess
    max_relevance = textRelevanceSet(GTree[u_node].userUKeySet, o_node_keySet);
#else
    TreeNode u_Node = getGIMTreeNodeData(u_node,OnlyU);
    vector<int> o_node_keys;
    for(int oterm: o_node_keySet)
        o_node_keys.push_back(oterm);
    max_relevance = textRelevance(u_Node.userUKey, o_node_keys);

#endif
    if(max_relevance>0){
        double max_influence = 1.0 + pow(a,1);
        if(a==0) max_influence = 1.0;

        double min_distance = MinDisBound[o_node][u_node];//100000;//getMinMaxDistanceBetweenNode(leaf, usr_node).max;
        double social_textual_maxscore = alpha*max_influence + (1.0-alpha)*max_relevance;
        double max_score = social_textual_maxscore / (1.0+ min_distance);
        return max_score;
    }
    else return  0;
    //score =0;

}



double getGSKScoreP2U_Upper_Plus(int a, double alpha, vector<int> Ukeys, set<int> o_node_keySet, int o_node, int u_node){


    double max_relevance = 0 ;

    set<int> _keyset;
    for(int term:Ukeys)
        _keyset.insert(term);

    //max_relevance = tfIDF(_keyset, o_node_keySet);

#ifndef DiskAccess
    max_relevance = textRelevanceSet(_keyset, o_node_keySet);
#else
    TreeNode u_Node = getGIMTreeNodeData(u_node,OnlyU);
    set<int> UkeySet;
    for(int term: u_Node.userUKey)
        UkeySet.insert(term);
    max_relevance = textRelevanceSet(UkeySet, o_node_keySet);

#endif



    if(max_relevance>0){
        double max_influence = 1.0 + pow(a,1); //这里可以优化jordan
        if(a==0)  max_influence = 1.0;

        double min_distance = MinDisBound[o_node][u_node];//100000;//getMinMaxDistanceBetweenNode(leaf, usr_node).max;

        double social_textual_maxscore = alpha*max_influence + (1.0-alpha)*max_relevance;

        double max_score = social_textual_maxscore / (1.0+ min_distance);
        //score =0;
        return max_score;
    }
    return  0;

}



double getGSKScoreP2U_Upper(int a, double alpha, set<int> o_node_keySet, int o_node, int u_node){


    double max_relevance = 0 ;

#ifndef DiskAccess
    max_relevance = tfIDF(GTree[u_node].userUKeySet, o_node_keySet);

#else
    TreeNode u_Node = getGIMTreeNodeData(u_node,OnlyU);
    set<int> UkeySet;
    for(int term: u_Node.userUKey)
        UkeySet.insert(term);
    max_relevance = textRelevanceSet(UkeySet, o_node_keySet);

#endif


    if(max_relevance>0){
        double max_influence = 1.0 + pow(a,1); //这里可以优化jordan
        if(a==0)  max_influence = 1.0;

        double min_distance = MinDisBound[o_node][u_node];//100000;//getMinMaxDistanceBetweenNode(leaf, usr_node).max;

        double social_textual_maxscore = alpha*max_influence + (1.0-alpha)*max_relevance;

        double max_score = social_textual_maxscore / (1.0+ min_distance);
        //score =0;
        return max_score;
    }
    return  0;

}

//ok
double getGSKScoreP2U_Upper_GivenKey(int a, double alpha, vector<int> given_keys, int o_node, int u_node){


    double max_relevance = 0 ;

    set<int> givenkeySet;
    for(int term_id:given_keys)
        givenkeySet.insert(term_id);

    set<int> UkeySet;

#ifndef DiskAccess
    UkeySet = GTree[u_node].userUKeySet;

#else
    TreeNode u_Node = getGIMTreeNodeData(u_node,OnlyU);

    for(int term: u_Node.userUKey)
        UkeySet.insert(term);
#endif

    max_relevance = textRelevanceSet(UkeySet, givenkeySet);


    if(max_relevance>0){
        double max_influence = 1.0 + pow(a,1); //这里可以优化jordan
        if(a==0)  max_influence = 1.0;

        double min_distance = MinDisBound[o_node][u_node];//100000;//getMinMaxDistanceBetweenNode(leaf, usr_node).max;

        double social_textual_maxscore = alpha*max_influence + (1.0-alpha)*max_relevance;

        double max_score = social_textual_maxscore / (1.0+ min_distance);
        //score =0;
        return max_score;
    }
    return  0;

}

//ok
double getGSKScoreP2U_Upper_InterKey(int a, double alpha, set<int> inter_keys, set<int> check_ins, int o_node, int u_node){


    double max_relevance = 0.0 ;

    set<int> givenkeySet;
    for(int term_id:inter_keys)
        givenkeySet.insert(term_id);


    max_relevance = tfIDF_termSet(inter_keys);
            //tfIDF(GTree[u_node].userUKeySet, givenkeySet);

    if(max_relevance==0)
        return 0.0;

    set<int> unionFriends ;
#ifndef DiskAccess
    unionFriends = GTree[u_node].unionFriends; //这里可以优化（optimization）
#else
    unionFriends = GTree[u_node].unionFriends_load;
#endif
    set<int> inter;
    inter = obtain_itersection_jins(unionFriends,check_ins);
    int _size = inter.size();
    int  min_popular = GTree[u_node].minSocialCount;
    float max_ratio =1.0;
    if(_size==0)
        max_ratio = 0;
    else if(min_popular == INT32_MAX)
        max_ratio = 1.0;
    else{
        double rate = _size*1.0 / min_popular;
        max_ratio = min(1.0,rate);
    }
    double max_influence = -1;
    if(max_ratio ==0)
        max_influence = 1.0;
    else
        max_influence = 1.0 + pow(a,max_ratio); //这里可以优化jordan


    double min_distance = -1;
    min_distance = MinDisBound[o_node][u_node];//100000;//getMinMaxDistanceBetweenNode(leaf, usr_node).max;
    //printf("MinDisBound[n%d][n%d]=%f\n",o_node,u_node,min_distance);


    double social_textual_maxscore = alpha*max_influence + (1.0-alpha)*max_relevance;

    double max_score = social_textual_maxscore / (1.0+ min_distance);
        //score =0;
    return max_score;


}



double getGSKScoreu2P_Upper(int a, double alpha, int store_leaf, int usr_id, vector<int> qKeys, vector<int> check_usr){
    double relevance = 0;
    User u;
#ifndef DiskAccess
    u = Users[usr_id];
#else
    u = getUserFromO2UOrgLeafData(usr_id);
#endif
    relevance = textRelevance(u.keywords,qKeys);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    if(relevance >0 ){
        double influence;
        double social_Score = 0.0;
        if(a==0){
            influence =1.0;
        }else{
            if(check_usr.size()>0)
                social_Score =  reCalcuSocialScore(usr_id,check_usr);

            if(social_Score == 0)
                influence =1.0;
            else
                influence =1.0 + pow(a,social_Score);

        }

        int usr_leaf = Nodes[u.Nj].gtreepath.back();

        double distance_min = MinDisBound[usr_leaf][store_leaf];

        double social_textual_score = alpha*influence + (1.0-alpha)*relevance;
        double gsk_score = social_textual_score / (1.0 +distance_min);
        return gsk_score;

    }
    else
        return  0;
}

//ok
double getGSKScoreu2O_Upper(int usr_id, int a, double alpha, double min_Dist, set<int> QKeySet){ //vector<int> UCheck
    double relevance = 0;
    vector<int> QKeys;
    //把set 转 vector
    for(int term: QKeySet)
        QKeys.push_back(term);
    User u;
#ifndef DiskAccess
    u = Users[usr_id];
#else
    u = getUserFromO2UOrgLeafData(usr_id);
#endif
    relevance = textRelevance(u.keywords,QKeys);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    if(relevance >0 ){
        double influence = 1.0 + a ;
        double social_textual_score = alpha*influence + (1.0-alpha)*relevance;
        double gsk_score = social_textual_score / (1.0 + min_Dist);
        return gsk_score;

    }
    else
        return  0;
}

double getGSKScoreU2O_Upper(vector<User> u_list, int a, double alpha, double min_Dist, set<int> OUnionKeys){ //vector<int> UCheck
    double relevance_max = 0;

    for(User u: u_list){
        set<int> uKeySet;
        for(int term: u.keywords)
            uKeySet.insert(term);
        double relevance = textRelevanceSet(uKeySet,OUnionKeys);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
        relevance_max = max(relevance_max, relevance);
    }


    return  relevance_max;
}




double getGSKScoreusrLeaf2O_Lower(int a, double alpha, double min_Dist, vector<int> Ukeys, set<int> OKeySet){ //vector<int> UCheck
    double relevance = 0;
    vector<int> OKeys;
    //把set 转 vector
    for(int term: OKeySet)
        OKeys.push_back(term);
    relevance = textRelevance(Ukeys,OKeys);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    if(relevance >0 ){
        double influence = 1.0; //可以优化
        double social_textual_score = alpha*influence + (1.0-alpha)*relevance;
        double gsk_score = social_textual_score / (1.0 + min_Dist);
        return gsk_score;

    }
    else
        return  0;
}


double getGSKScoreusrLeaf2O_Upper(int a, double alpha, double min_Dist, vector<int> Ukeys, set<int> OKeySet){ //vector<int> UCheck


    double relevance = 0;
    vector<int> OKeys;
    //把set 转 vector
    for(int term: OKeySet)
        OKeys.push_back(term);
    relevance = textRelevance(Ukeys,OKeys);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    if(relevance >0 ){
        double influence = 1.0 + a ; //可以优化
        double social_textual_score = alpha*influence + (1.0-alpha)*relevance;
        double gsk_score = social_textual_score / (1.0 + min_Dist);
        return gsk_score;

    }
    else
        return  0.0;
}





double getGSKScoreusrLeaf2o_Upper(int a, double alpha, double min_Dist, vector<int> Ukeys, vector<int> OKeys){ //vector<int> UCheck
    double relevance = 0;
    /*vector<int> OKeys;
    //把set 转 vector
    for(int term: OKeySet)
        OKeys.push_back(term);*/
    relevance = textRelevance(Ukeys,OKeys);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    if(relevance >0 ){
        double influence = 1.0 +a ; //可以优化
        double social_textual_score = alpha*influence + (1.0-alpha)*relevance;
        double gsk_score = social_textual_score / (1.0 + min_Dist);
        return gsk_score;

    }
    else
        return  0.0;
}



double getGSKScoreo2U_Lower(int a, double alpha, int o_id, int u_leaf, int term_id){
    POI p;
#ifndef DiskAccess
    p = POIs[o_id];
#else
    p = getPOIFromO2UOrgLeafData(o_id);
#endif
    int loc = p.Ni;

    //double relevance = 5.0 / UkeySizeMax;  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    double relevance = 0;
    relevance = tfIDF_term(term_id);
    if(relevance >0){
        double influence;
        double social_Score = 0;
        //social_Score =  reCalcuSocialScore(usr_id,poiCheckInList[o_id]); //calcSocialScore(usr_id,o_id);
        if(social_Score == 0)
            influence =1.0;
        else
            influence =1.0 + pow(a,social_Score);
        int o_leaf = Nodes[p.Ni].gtreepath.back();
        //int u_leaf = Nodes[u_leaf].gtreepath.back();
        double  distance = MaxDisBound[o_leaf][u_leaf];
        //double distance = getDistance(loc,Users[usr_id]);
        double score =0;
        double social_textual_score = alpha*influence + (1.0-alpha)*relevance;
        score = social_textual_score  / (1.0+distance);
        //score =0;
        return score;

    }
    else
        return  0.0;

}

double getGSKScoreo2U_Upper(int a, double alpha, int o_id, int u_leaf){
    int loc = POIs[o_id].Ni;

    double relevance = 0.0;
    POI p; TreeNode u_leafNode;
#ifndef DiskAccess
    relevance = tfIDF(POIs[o_id].keywordSet, GTree[u_leaf].userUKeySet);   //考虑兴趣点包含多个用户关键词（超级用户的关键词）
#else
    p = getPOIFromO2UOrgLeafData(o_id);
    u_leafNode = getGIMTreeNodeData(u_leaf,OnlyU);
    relevance = textRelevance(p.keywords, u_leafNode.userUKey);
#endif
    if(relevance >0){
        double influence;
        double social_Score = 1.0;  //社交最优情况
        //social_Score =  reCalcuSocialScore(usr_id,poiCheckInList[o_id]); //calcSocialScore(usr_id,o_id);
        if(social_Score == 0)
            influence =1.0;
        else
            influence =1.0 + pow(a,social_Score);
        int o_leaf = Nodes[POIs[o_id].Ni].gtreepath.back();
        //int u_leaf = Nodes[u_leaf].gtreepath.back();
        double  distance = MinDisBound[o_leaf][u_leaf];
        double score_upper =0.0;
        double social_textual_score = alpha*influence + (1.0-alpha)*relevance;
        score_upper = social_textual_score / (1.0+distance);
        //score =0;
        return score_upper;

    }
    else
        return  0;

}

//disk ok
double getGSKScoreo2u_Lower(int a, double alpha, int o_id, int usr_id){
    POI p; User u;
#ifndef DiskAccess
    p = POIs[o_id];
    u = Users[usr_id];
#else
    p = getPOIFromO2UOrgLeafData(o_id);
    u = getUserFromO2UOrgLeafData(usr_id);
#endif

    int loc = p.Ni;
    double relevance = textRelevance(u.keywords,p.keywords);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    //cout<<"relevance ="<<relevance<<endl;
    if(relevance >0){
        double influence;
        double social_Score = 0.0;
        if(alpha>0 && a>0){
            social_Score =  reCalcuSocialScore(usr_id,poiCheckInIDList[o_id]); //calcSocialScore(usr_id,o_id);
        }

        if(social_Score == 0)
            influence =1.0;
        else
            influence =1.0 + pow(a,social_Score);
        int o_leaf = Nodes[p.Ni].gtreepath.back();
        int u_leaf = Nodes[u.Ni].gtreepath.back();
        double  distance = MaxDisBound[o_leaf][u_leaf];
        //double distance = usrToPOIDistance(Users[usr_id],POIs[o_id]);//getDistance(loc,Users[usr_id]);
        double score =0.0;
        double social_textual_score = alpha*influence + (1.0-alpha)*relevance;
        score = social_textual_score / (1.0+distance);
        //score =0;
        return score;

    }
    else
        return  0;

}



double getGSKScoreuData2oData(User u, POI p, int a, double alpha, double distance){

    int loc = p.Ni;

    double relevance = textRelevance(u.keywords,p.keywords);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    if(relevance >0){
        double influence;
        double social_Score = 0.0;
        if(alpha>0 &&a >0){
            social_Score =  reCalcuSocialScore(u.id,poiCheckInIDList[p.id]); //calcSocialScore(usr_id,o_id);
        }

        if(social_Score == 0)
            influence =1.0;
        else
            influence =1.0 + pow(a,social_Score);

        int o_leaf = Nodes[p.Ni].gtreepath.back();
        int u_leaf = Nodes[u.Ni].gtreepath.back();
        //double  distance = MaxDisBound[o_leaf][u_leaf];
        //double distance = usrToPOIDistance(Users[usr_id],POIs[o_id]);//getDistance(loc,Users[usr_id]);
        double social_textual_score = alpha*influence + (1.0-alpha)*relevance;

        double score =0.0;
        score = social_textual_score / (1.0+distance);
        //score =0;
        return score;

    }
    else
        return  0;

}

double getGSKScoreuData2oDataTextual(User u, POI p, int a, double alpha, double textual_relevance, double distance){

    int loc = p.Ni;

    double relevance = textual_relevance;  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    if(relevance >0){
        double influence;
        double social_Score = 0.0;
        if(alpha>0 &&a >0){
            social_Score =  reCalcuSocialScore(u.id,poiCheckInIDList[p.id]); //calcSocialScore(usr_id,o_id);
        }

        if(social_Score == 0)
            influence =1.0;
        else
            influence =1.0 + pow(a,social_Score);

        int o_leaf = Nodes[p.Ni].gtreepath.back();
        int u_leaf = Nodes[u.Ni].gtreepath.back();
        //double  distance = MaxDisBound[o_leaf][u_leaf];
        //double distance = usrToPOIDistance(Users[usr_id],POIs[o_id]);//getDistance(loc,Users[usr_id]);
        double social_textual_score = alpha*influence + (1.0-alpha)*relevance;

        double score =0.0;
        score = social_textual_score / (1.0+distance);
        //score =0;
        return score;

    }
    else
        return  0;

}


double getGSKScoreu2o(int usr_id, int o_id, int a, double alpha, double distance){
    POI p; User u;
#ifndef DiskAccess
    p = POIs[o_id];
    u = Users[usr_id];
#else
    p = getPOIFromO2UOrgLeafData(o_id);
    u = getUserFromO2UOrgLeafData(usr_id);
#endif

    int loc = p.Ni;

    double relevance = textRelevance(u.keywords,p.keywords);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    if(relevance >0){
        double influence;
        double social_Score = 0.0;
        if(alpha>0 &&a >0){
            social_Score =  reCalcuSocialScore(usr_id,poiCheckInIDList[o_id]); //calcSocialScore(usr_id,o_id);
        }

        if(social_Score == 0)
            influence =1.0;
        else
            influence =1.0 + pow(a,social_Score);

        int o_leaf = Nodes[p.Ni].gtreepath.back();
        int u_leaf = Nodes[u.Ni].gtreepath.back();
        //double  distance = MaxDisBound[o_leaf][u_leaf];
        //double distance = usrToPOIDistance(Users[usr_id],POIs[o_id]);//getDistance(loc,Users[usr_id]);
        //ts = relevance; ss = influence;
        double social_textual_score = alpha*influence + (1.0-alpha)*relevance;
        double score =0.0;
        score = social_textual_score / (1.0+distance);
        //score =0;
        return score;

    }
    else
        return  0;

}


double getGSKScore_o2u_Upper(int a, double alpha, POI poi, User user, double& rel, double& inf){   //poi与user的GSK评分
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
    double influence;
    double social_Score = 0;
    if(a==0||alpha==0){
        influence =1.0;
    }else{
        if(check_usr.size()>0)
            social_Score =  reCalcuSocialScore(user.id,check_usr);

        if(social_Score == 0)
            influence =1;
        else
            influence =1 + pow(a,social_Score);
    }

    //double distance = usrToPOIDistance(user, poi);   //这里要注意区分 getDistance(vertex, vertex)

    int usr_leaf = Nodes[user.Nj].gtreepath.back();

    int poi_leaf = Nodes[poi.Nj].gtreepath.back();

    float distance_min = MinDisBound[usr_leaf][poi_leaf];

    double social_textual_score = alpha*influence + (1.0-alpha)*relevance;

    double gsk_score_upper = social_textual_score / (1.0+distance_min);
    rel = relevance; inf = influence;
    //if(poi.id== 57&& user.id ==625 ) cout<<"getGSKScore_o2u，发现u57,p625,distance="<<distance<<",relevance="<<relevance<<"influence="<<influence<<",score="<<gsk_score<<endl;
    return gsk_score_upper;

}


double getGSKScore_q2u_Lower(int a, double alpha, vector<int>qKey, POI poi, User user, double& rel, double& inf){   //poi与user的GSK评分
    set<int> uKeyset; set<int> qKeyset;
    //for(int term : Users[usr_id].keywords) uKeyset.insert(term);
    //for(int term :qKeys)  qKeyset.insert(term);
    //double relevance = textRelevance(uKeyset,qKeyset);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    vector<int> qKeys = qKey; //poi.keywords;
    vector<int> check_usr = poiCheckInIDList[poi.id];
    double relevance = textRelevance(user.keywords,qKeys);
    if(relevance == 0){
        return  0.0;
    }
    double influence;
    double social_Score = 0;
    if(a==0||alpha==0){
        influence =1.0;
    }
    else{
        if(check_usr.size()>0)
            social_Score =  reCalcuSocialScore(user.id,check_usr);

        if(social_Score == 0)
            influence =1.0;
        else
            influence =1.0+ pow(a,social_Score);
    }

    //double distance = usrToPOIDistance(user, poi);   //这里要注意区分 getDistance(vertex, vertex)

    int usr_leaf = Nodes[user.Nj].gtreepath.back();

    int poi_leaf = Nodes[poi.Nj].gtreepath.back();

    //float distance_min = MinDisBound[usr_leaf][poi_leaf];

    float distance_max = MaxDisBound[usr_leaf][poi_leaf];

    double social_textual_score = alpha*influence + (1.0-alpha)*relevance;

    double gsk_score_lower = social_textual_score / (1.0+distance_max);
    rel = relevance; inf = influence;
    //if(poi.id== 57&& user.id ==625 ) cout<<"getGSKScore_o2u，发现u57,p625,distance="<<distance<<",relevance="<<relevance<<"influence="<<influence<<",score="<<gsk_score<<endl;
    return gsk_score_lower;

}


double getGSKScore_q2u_Upper(int a, double alpha, vector<int>qKey, POI poi, User user, double& rel, double& inf){   //poi与user的GSK评分
    set<int> uKeyset; set<int> qKeyset;
    //for(int term : Users[usr_id].keywords) uKeyset.insert(term);
    //for(int term :qKeys)  qKeyset.insert(term);
    //double relevance = textRelevance(uKeyset,qKeyset);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    vector<int> qKeys = qKey; //poi.keywords;
    vector<int> check_usr = poiCheckInIDList[poi.id];
    double relevance = textRelevance(user.keywords,qKeys);
    if(relevance == 0){
        return  0.0;
    }
    double influence;
    double social_Score = 0;
    if(a==0||alpha==0){
        influence =1.0;
    }
    else{
        if(check_usr.size()>0)
            social_Score =  reCalcuSocialScore(user.id,check_usr);

        if(social_Score == 0)
            influence =1.0;
        else
            influence =1.0 + pow(a,social_Score);
    }

    //double distance = usrToPOIDistance(user, poi);   //这里要注意区分 getDistance(vertex, vertex)

    int usr_leaf = Nodes[user.Nj].gtreepath.back();

    int poi_leaf = Nodes[poi.Nj].gtreepath.back();

    float distance_min = MinDisBound[usr_leaf][poi_leaf];

    double social_textual_score = alpha*influence + (1.0-alpha)*relevance;

    double gsk_score_upper = social_textual_score / (1.0+distance_min);
    rel = relevance; inf = influence;
    //if(poi.id== 57&& user.id ==625 ) cout<<"getGSKScore_o2u，发现u57,p625,distance="<<distance<<",relevance="<<relevance<<"influence="<<influence<<",score="<<gsk_score<<endl;
    return gsk_score_upper;

}


double getGSKScore_q2u(int a, double alpha, POI poi, User user, double relevance, double influence){   //poi与user的GSK评分

    double distance = usrToPOIDistance(user, poi);   //这里要注意区分 getDistance(vertex, vertex)

    double social_textual_score = alpha*influence + (1.0-alpha)*relevance;

    double gsk_score = social_textual_score / (1+distance);
    //if(poi.id== 57&& user.id ==625 ) cout<<"getGSKScore_o2u，发现u57,p625,distance="<<distance<<",relevance="<<relevance<<"influence="<<influence<<",score="<<gsk_score<<endl;
    return gsk_score;

}




typedef  struct {
    double upper_score;
    double lower_score;
}Score_Bound;

Score_Bound getGSKScore_q2u_Upper_and_Lower(int a, double alpha, vector<int>qKey, POI poi, User user, double& rel, double& inf){   //poi与user的GSK评分
    set<int> uKeyset; set<int> qKeyset;
    //for(int term : Users[usr_id].keywords) uKeyset.insert(term);
    //for(int term :qKeys)  qKeyset.insert(term);
    //double relevance = textRelevance(uKeyset,qKeyset);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    vector<int> qKeys = qKey; //poi.keywords;
    vector<int> check_usr = poiCheckInIDList[poi.id];
    double relevance = textRelevance(user.keywords,qKeys);

    if(relevance==0){
        Score_Bound score_bound;
        score_bound.upper_score = 0;
        score_bound.lower_score = 0;
        return score_bound;
    }
    
    double influence;
    double social_Score = 0;
    if(a==0 || alpha==0){
        influence =1.0;
    }
    else{
        if(check_usr.size()>0)
            social_Score =  reCalcuSocialScore(user.id,check_usr);

        if(social_Score == 0)
            influence =1;
        else
            influence =1 + pow(a,social_Score);
    }

    //double distance = usrToPOIDistance(user, poi);   //这里要注意区分 getDistance(vertex, vertex)

    int usr_leaf = Nodes[user.Nj].gtreepath.back();

    int poi_leaf = Nodes[poi.Nj].gtreepath.back();

    float distance_min = MinDisBound[usr_leaf][poi_leaf];
    float distance_max = MaxDisBound[usr_leaf][poi_leaf];
    //double distance = usrToPOIDistance(user, poi);   //这里要注意区分 getDistance(vertex, vertex)

    double social_textual_score = alpha*influence + (1-alpha)*relevance;

    double gsk_score_upper = social_textual_score / (1.0+distance_min);
    double gsk_score_lower = social_textual_score / (1.0+distance_max);
    rel = relevance; inf = influence;
    Score_Bound score_bound;
    score_bound.upper_score = gsk_score_upper;
    score_bound.lower_score = gsk_score_lower;
    //if(poi.id== 57&& user.id ==625 ) cout<<"getGSKScore_o2u，发现u57,p625,distance="<<distance<<",relevance="<<relevance<<"influence="<<influence<<",score="<<gsk_score<<endl;
    return score_bound;

}

#define Score_Debug1

Score_Bound getGSKScore_o2u_Upper_and_Lower(int a,  double alpha, POI poi, User user, double& rel, double& inf){   //poi与user的GSK评分
    set<int> uKeyset; set<int> qKeyset;
#ifdef Score_Debug
    if(poi.id==57&&user.id==106){
        cout<<"poi.id==57&&user.id==106"<<endl;
    }
#endif
    //for(int term : Users[usr_id].keywords) uKeyset.insert(term);
    //for(int term :qKeys)  qKeyset.insert(term);
    //double relevance = textRelevance(uKeyset,qKeyset);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    vector<int> qKeys = poi.keywords;
    vector<int> check_usr = poi.check_ins;//poiCheckInList[poi.id];
    double relevance = textRelevance(user.keywords,qKeys);
    if(relevance==0){
        Score_Bound score_bound;
        score_bound.upper_score = 0;
        score_bound.lower_score = 0;
        return score_bound;
    }
        
    double influence;
    double social_Score = 0;
    if(a==0||alpha==0){
        influence =1.0;
    }
    else{
        if(check_usr.size()>0)
            social_Score =  reCalcuSocialScore(user.id,check_usr);

        if(social_Score == 0)
            influence =1;
        else
            influence =1 + pow(a,social_Score);
    }


    double distance = usrToPOIDistance(user, poi);   //这里要注意区分 getDistance(vertex, vertex)

    int usr_leaf = Nodes[user.Nj].gtreepath.back();

    int poi_leaf = Nodes[poi.Nj].gtreepath.back();

    double distance_min = distance;
    double  distance_max = distance;
    //float distance = usrToPOIDistance(user, poi);   //这里要注意区分 getDistance(vertex, vertex)
    //distance_max = distance;
    //cout<<"distance="<<distance<<",distance_min="<<distance_min<<",distance_max="<<distance_max<<endl;
    //getchar();

    double social_textual_score = alpha*influence + (1-alpha)*relevance;

    double gsk_score_upper = social_textual_score / (1.0+distance_min);
    double gsk_score_lower = social_textual_score / (1.0+distance_max);
    rel = relevance; inf = influence;
    Score_Bound score_bound;
    score_bound.upper_score = gsk_score_upper;
    score_bound.lower_score = gsk_score_lower;
    //if(poi.id== 57&& user.id ==625 ) cout<<"getGSKScore_o2u，发现u57,p625,distance="<<distance<<",relevance="<<relevance<<"influence="<<influence<<",score="<<gsk_score<<endl;
    return score_bound;

}


Score_Bound getGSKScore_o2u_Upper_and_Lower_givenKey(int a,  double alpha, POI poi, User user, vector<int> qKeys, double& rel, double& inf){   //poi与user的GSK评分
    set<int> uKeyset; set<int> qKeyset;
    //for(int term : Users[usr_id].keywords) uKeyset.insert(term);
    //for(int term :qKeys)  qKeyset.insert(term);
    //double relevance = textRelevance(uKeyset,qKeyset);  //考虑 u.keywords与候选关键词（qKey)的交集对应相关的对象
    vector<int> check_usr = poi.check_ins; //poiCheckInList[poi.id];
    double relevance = textRelevance(user.keywords,qKeys);
    double influence;
    double social_Score = 0;
    if(a==0||alpha==0){
        influence =1.0;
    }
    else{
        if(check_usr.size()>0)
            social_Score =  reCalcuSocialScore(user.id,check_usr);

        if(social_Score == 0)
            influence =1;
        else
            influence =1 + pow(a,social_Score);
    }


    double distance = usrToPOIDistance(user, poi);   //这里要注意区分 getDistance(vertex, vertex)

    int usr_leaf = Nodes[user.Nj].gtreepath.back();

    int poi_leaf = Nodes[poi.Nj].gtreepath.back();

    double distance_min = distance;
    double  distance_max = distance;
    //float distance = usrToPOIDistance(user, poi);   //这里要注意区分 getDistance(vertex, vertex)
    //distance_max = distance;
    //cout<<"distance="<<distance<<",distance_min="<<distance_min<<",distance_max="<<distance_max<<endl;
    //getchar();

    double social_textual_score = alpha*influence + (1-alpha)*relevance;

    double gsk_score_upper = social_textual_score / (1.0+distance_min);
    double gsk_score_lower = social_textual_score / (1.0+distance_max);
    rel = relevance; inf = influence;
    Score_Bound score_bound;
    score_bound.upper_score = gsk_score_upper;
    score_bound.lower_score = gsk_score_lower;
    //if(poi.id== 57&& user.id ==625 ) cout<<"getGSKScore_o2u，发现u57,p625,distance="<<distance<<",relevance="<<relevance<<"influence="<<influence<<",score="<<gsk_score<<endl;
    return score_bound;

}





#define MAXINFRKGSKQ_1_1_FUNCTIONAL_H

#endif //MAXINFRKGSKQ_1_1_FUNCTIONAL_H
