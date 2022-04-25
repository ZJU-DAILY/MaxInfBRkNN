//
// Created by jins on 12/7/19.
//
#include <iostream>
#include "query_plus.h"
#include "gim_tree.h"
#include "InfluenceModels.h"
//#include "process.h"

#ifndef TEST_H

namespace constants{

    ////for batch reverse
    std::string const GRP_ALG = "group";  //对应第二个算法
    std::string const NVD_ALG = "nvd"; //基于nvd的RkGSKQ算法

    //// for max-inf
    std::string const BASE_ALG = "base";  //基于group + OPIMC 的baseline算法
    std::string const APPRO_ALG = "appro"; // 基于nvd + hybrid 下的近似算法
    std::string const HEURISTIC_ALG = "heuris"; //基于预采样下的启发式算法


    //for experiment selections



}

enum Query_Type{Single, Multi};
enum Processing_Pattern{Naive, ONE_BY_ONE, BATCH};
enum Single_Method{GROUP_SINGLE, NVD_SINGLE};
enum Batch_Method{CLUSTER, SEPRATE, GROUP_BATCH, NVD_BATCH};

enum Keyword_Popularity{Popular, Mid, Rare};
enum Solution{NAIVE, BASELINE, APPROXIMATE, HEURISTIC, ONEHOP, TWOHOP, RANDOM, RELEVANCE, CARDINALITY, INFLUENCER};  //  in BASE, we retrieve the RkGSKQ for each object one by one
//产生一组查询对象
#define DEFAULT_ALG    BASELINE


void test_load_V2P_NVDG_AddressIdx(){
    clock_t startTime, endTime;

    startTime = clock();
    string indexPath = getNVDHashIndexInputPath();
    V2PNVDHashIndex index = serialization::getIndexFromBinaryFile<V2PNVDHashIndex>(indexPath);
    endTime = clock();
    cout<<"序列化文件 导入 完毕! 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s" << endl;
    index.initialIndexGiven(idHashMap_v2p,addressHashMap_nvd);
    unordered_map<int,int>::iterator iterator_v2p;
    cout<<"idHashMap_v2p content:"<<endl;
    int i = 0;
    for(iterator_v2p = idHashMap_v2p.begin(); iterator_v2p!=idHashMap_v2p.end(); iterator_v2p++){
        int key = iterator_v2p->first;
        int value = iterator_v2p -> second;

        if (i%100==0){
            //cout<<"<key="<<key<<",value="<<value<<">"<<endl;
        }
        i++;
    }

    unordered_map<int,int>::iterator iterator_nvd;
    cout<<"addressHashMap_nvd content:"<<endl;
    i =0;
    for(iterator_nvd = addressHashMap_nvd.begin(); iterator_nvd!=addressHashMap_nvd.end(); iterator_nvd++){
        int key = iterator_nvd->first;
        int value = iterator_nvd -> second;
        //cout<<"<key="<<key<<",value="<<value<<">"<<endl;
        if (i%50==0){
            cout<<"<nvd adj key="<<key<<",value="<<value<<">"<<endl;
        }
        i++;
    }



}

void test_write_HybirdPOINN_IDX(){
    clock_t startTime, endTime;
    startTime = clock();
    string indexOutputPre = getHybridNVDHashIndexInputPath();
    stringstream ss;
    ss<<indexOutputPre;
    string indexOutputFilePath = ss.str();
    cout<<"*****************序列化输出V2PNVDHashIndex...";

    addressHashMap_nvd[1] = 12;
    addressHashMap_nvd[2] = 12;

    idHashMap_v2p[10] = 88;
    idHashMap_v2p[20] = 118;

    vector<int> _list;
    for(int i=0;i<8;i++){
        _list.push_back(i);
    }
    idListHash_l2p_hybrid[30] = _list;



    HybridPOINNHashIndex hashIndex;
    hashIndex.setNVDIndex(addressHashMap_nvd);
    hashIndex.setV2PIndex(idHashMap_v2p);
    hashIndex.setL2PIndex(idListHash_l2p_hybrid);


    serialization::outputIndexToBinaryFile<HybridPOINNHashIndex>(hashIndex,indexOutputFilePath);
    endTime = clock();
    cout<<"输出完毕! 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s*****************" << endl;

}

void test_load_HybirdPOINN_AddressIdx(){
    clock_t startTime, endTime;

    startTime = clock();
    string indexPath = getHybridNVDHashIndexInputPath();
    HybridPOINNHashIndex index = serialization::getIndexFromBinaryFile<HybridPOINNHashIndex>(indexPath);
    endTime = clock();
    cout<<"序列化文件 导入 hybrid poiNN index 完毕! 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s" << endl;
    index.initialIndexGiven(idHashMap_v2p,idListHash_l2p_hybrid, addressHashMap_nvd);
    unordered_map<int,int>::iterator iterator_v2p;
    cout<<"idHashMap_v2p content:"<<endl;
    int i = 0;
    for(iterator_v2p = idHashMap_v2p.begin(); iterator_v2p!=idHashMap_v2p.end(); iterator_v2p++){
        int key = iterator_v2p->first;
        int value = iterator_v2p -> second;

        //if (i%100==0){
            cout<<"<key="<<key<<",poi_id="<<value<<">"<<endl;
        //}
        i++;
    }

    unordered_map<int,vector<int>>:: iterator iterator_l2ps;
    cout<<"idListHashMap_l2p 中内容:"<<endl;
    for(iterator_l2ps = idListHash_l2p_hybrid.begin();iterator_l2ps!=idListHash_l2p_hybrid.end();iterator_l2ps++){
        int key = iterator_l2ps->first;
        vector<int> pois_list = iterator_l2ps->second;
        cout<<"key="<<key; cout<<", id_list:";
        printElements(pois_list);

    }


    unordered_map<int,int>::iterator iterator_nvd;
    cout<<"addressHashMap_nvd content:"<<endl;
    i =0;
    for(iterator_nvd = addressHashMap_nvd.begin(); iterator_nvd!=addressHashMap_nvd.end(); iterator_nvd++){
        int key = iterator_nvd->first;
        int value = iterator_nvd -> second;
        //cout<<"<key="<<key<<",value="<<value<<">"<<endl;
        //if (i%50==0){
            cout<<"<nvd adj key="<<key<<",value="<<value<<">"<<endl;
        //}
        i++;
    }



}



void test_getBorder_SB_NVD(int K, int a, double alpha){
    POI poi = getPOIFromO2UOrgLeafData(57);
    set<int> qo_nonRarekeys;
    for(int term: poi.keywords){
        int _size = getTermOjectInvListSize(term);
        if(_size>K)
            qo_nonRarekeys.insert(term);
    }

    CheckEntry border_entry(1337,80,qo_nonRarekeys,986,0);
    //getBorder_SB_NVD(border_entry, poi, qo_nonRarekeys, K, a,alpha);
    getBorder_SB_NVD(border_entry, K, a,alpha);
}




vector<int> generate_query(){
    vector <int > chainStores;
    int test =0;
    if(test){
        string str = "57 60 155 900 80 39";
        /*stringstream ss = stringstream(str);
        int store_id;
        while(ss>>store_id)
            chainStores.push_back(store_id);*/
    }
    else{
        int query_num = DEFAULT_NUMBEROFSTORE;   //50, 100, 200, 400, 800  RkGSKQ查询的数量
        while(query_num){
            chainStores.push_back(query_num);
            query_num--;
        }
    }
    return  chainStores;
}


vector<int> generate_query_popular(int poi_id, int candidate_poi_size){
    vector <int > chainStores;

    int query_num = candidate_poi_size;   //50, 100, 200, 400, 800  RkGSKQ查询的数量
    cout<<"*******************************************BRkNN serarch for:"<<endl;
    for(int i=0;i<query_num;i++){
        chainStores.push_back(poi_id+i);
        cout<<"p"<<poi_id+i<<",";
    }
    cout<<endl;
    return  chainStores;
}

vector<int> generate_query_popular(int poi_id){
    vector <int > chainStores;

    int query_num = DEFAULT_NUMBEROFSTORE;   //50, 100, 200, 400, 800  RkGSKQ查询的数量

    for(int i=0;i<query_num;i++){
        chainStores.push_back(poi_id+i);
    }

    return  chainStores;
}


vector<int> generate_query_PopularDominant(int query_type, int Qn){
    int pivot_keyword;  bool isPopular = false;
    set<int> poi_poolset; vector<int> poi_poolList;
    vector<int> pois;
    cout<<"创建"<<Qn<<"个RkGSKQ查询,";  //Qn: 25, 50, 100, 200, 400;   |KQ|: 20, 40, 60, 80, 80, 100
    switch(query_type){
        case 1: //group 1 (popular POI with dominant keywords

            pivot_keyword = 1; isPopular = true;

            cout<<"关键词类型：dominant, 兴趣点类型： popular!"<<endl;
            break;
        case 2:

            pivot_keyword = 2; isPopular = false;

            cout<<"关键词类型：dominant, 兴趣点类型： non-popular!"<<endl;
            break;
        case 3:
            pivot_keyword = 20;  isPopular = true;
            cout<<"关键词类型：non-dominant, 兴趣点类型： popular!"<<endl;
            break;
        case 4:
            pivot_keyword = 25;  isPopular = false;
            cout<<"关键词类型：non-dominant, 兴趣点类型： non-popular!"<<endl;
            break;

        default:
            pivot_keyword = 20;  isPopular = true;
            cout<<"关键词类型：non-dominant, 兴趣点类型： popular!"<<endl;
            break;

    }
    //pivot_keyword2 = pivot_keyword+10;
    int KQ_num = Qn* 0.4;


    for(int target_term = pivot_keyword; target_term< pivot_keyword+5;target_term++){
        cout<<"加入t"<<(target_term)<<"的相关用户"<<endl;
        vector<int> current_usersList = getTermOjectInvList(target_term);
        for(int poi_id: current_usersList){
            poi_poolset.insert(poi_id);
        }
    }


    set<int> extract_pois; int i_th = 0;
    map<int, POI> poi_cache;
    int _count =0;
    for(int poi_id: poi_poolset){
        _count++;
        POI poi = getPOIFromO2UOrgLeafData(poi_id);
        poi_cache[poi_id]= poi;
        //cout<<"line"<<_count<<", p"<<poi.id<<".id="<<poi.id<<", checkin_count="<<poi.check_ins.size()<<endl;
        if(query_type ==2 || query_type ==4){  // check-in 稀疏
            if(poi.check_ins.size()==0){
                extract_pois.insert(poi_id);
            }

        }
        else{   //check-in 密集 (必须筛选出有check-in的poi)
#ifdef LasVegas
            if(poi.check_ins.size()>=0){  // check-in 密集
                extract_pois.insert(poi_id);
                //getchar();

            }
#else
            if(poi.check_ins.size()>=0){  // check-in 密集 (>=1)
                extract_pois.insert(poi_id);
                //getchar();
            }
#endif


        }
        if(extract_pois.size()==200) break;
    }

    //getchar();
    //从中抽取 |P_c|个POI作为candidate
    vector<int> query_objects;
    for(int poi_id: extract_pois){
        //if(curr_keywords.size()>keyword_perpoi) continue;
        query_objects.push_back(poi_id);
        cout<<"加入了兴趣点为：";
        printPOIInfo(poi_cache[poi_id]);
        if(query_objects.size()==Qn) break;
    }
    printf("*********总共加入了: %d个 兴趣点!*****************\n",query_objects.size());
    //getchar();

    return query_objects;

}


vector<int> generate_query_PopularDominantLasVegas(int query_type, int Qn){
    int pivot_keyword;  bool isPopular = false;
    set<int> poi_poolset;
    set<int> poi_poolset2; vector<int> poi_poolList2;
    vector<int> poi_poolList;
    vector<int> pois;
    cout<<"创建"<<Qn<<"个RkGSKQ查询,";  //Qn: 25, 50, 100, 200, 400;   |KQ|: 20, 40, 60, 80, 80, 100
    switch(query_type){
        case 1: //group 1 (popular POI with dominant keywords

            pivot_keyword = 1; isPopular = true;

            cout<<"关键词类型：dominant, 兴趣点类型： popular!"<<endl;
            break;
        case 2:

            pivot_keyword = 2; isPopular = false;

            cout<<"关键词类型：dominant, 兴趣点类型： non-popular!"<<endl;
            break;
        case 3:
            pivot_keyword = 20;  isPopular = true;
            cout<<"关键词类型：non-dominant, 兴趣点类型： popular!"<<endl;
            break;
        case 4:
            pivot_keyword = 25;  isPopular = false;
            cout<<"关键词类型：non-dominant, 兴趣点类型： non-popular!"<<endl;
            break;

        default:
            pivot_keyword = 20;  isPopular = true;
            cout<<"关键词类型：non-dominant, 兴趣点类型： popular!"<<endl;
            break;

    }
    //pivot_keyword2 = pivot_keyword+10;
    int KQ_num = Qn* 0.4;

    ifstream fin;
    //stringstream tt;
    //tt<<FILE_DISBOUND;
    stringstream ss;
    ss<<"../../../data/LasVegas/poi_candidate/Q"<<query_type;
    string fileName = ss.str();

    //第一部分poi
    fin.open(fileName);
    string str;
    while(getline(fin, str)){
        istringstream tt(str);
        int query_id = -1;
        tt>>query_id;
        poi_poolset.insert(query_id);

    }
    fin.close();

    set<int> extract_pois; int i_th = 0;
    map<int, POI> poi_cache;
    int _count =0;
    for(int poi_id: poi_poolset){
        _count++;
        POI poi = getPOIFromO2UOrgLeafData(poi_id);
        poi_cache[poi_id]= poi;
        //cout<<"line"<<_count<<", p"<<poi.id<<".id="<<poi.id<<", checkin_count="<<poi.check_ins.size()<<endl;
        if(query_type ==2 || query_type ==4){  // check-in 稀疏
            if(poi.check_ins.size()==0){
                extract_pois.insert(poi_id);
            }

        }
        else{   //check-in 密集 (必须筛选出有check-in的poi)
#ifdef LasVegas
            if(poi.check_ins.size()>=0){  // check-in 密集
                extract_pois.insert(poi_id);
                //getchar();

            }
#else
            if(poi.check_ins.size()>=0){  // check-in 密集 (>=1)
                extract_pois.insert(poi_id);
                //getchar();
            }
#endif


        }
        if(extract_pois.size()==200) break;
    }
    ////第二部分poi
    for(int target_term = pivot_keyword; target_term< pivot_keyword+5;target_term++){
        cout<<"加入t"<<(target_term)<<"的相关用户"<<endl;
        vector<int> current_usersList = getTermOjectInvList(target_term);
        for(int poi_id: current_usersList){
            poi_poolset2.insert(poi_id);
        }
    }
    for(int poi_id:poi_poolset2)
        poi_poolList2.push_back(poi_id);

    //getchar();
    //从中抽取 |P_c|个POI作为candidate
    vector<int> query_objects; set<int> query_objectSet;
    for(int poi_id: extract_pois){
        //if(curr_keywords.size()>keyword_perpoi) continue;
        query_objectSet.insert(poi_id);
       /* if(query_objectSet.size()>=5){
            int _random_th = random()%poi_poolList2.size();
            query_objectSet.insert(poi_poolList2[_random_th]);
        }*/

        if(query_objectSet.size()==Qn) break;
    }
    printf("*********总共加入了: %d个 兴趣点!*****************\n",query_objects.size());
    //getchar();

    return query_objects;

}


vector<int> generate_query_PopularDominantNonEmpty(int query_type, int Qn){
    int pivot_keyword;  bool isPopular = false;
    set<int> poi_poolset;
    set<int> poi_poolset2; vector<int> poi_poolList2;
    vector<int> poi_poolList;
    vector<int> pois;
    cout<<"创建"<<Qn<<"个RkGSKQ查询,";  //Qn: 25, 50, 100, 200, 400;   |KQ|: 20, 40, 60, 80, 80, 100
    switch(query_type){
        case 1: //group 1 (popular POI with dominant keywords

            pivot_keyword = 1; isPopular = true;

            cout<<"关键词类型：dominant, 兴趣点类型： popular!"<<endl;
            break;
        case 2:

            pivot_keyword = 2; isPopular = false;

            cout<<"关键词类型：dominant, 兴趣点类型： non-popular!"<<endl;
            break;
        case 3:
            pivot_keyword = 20;  isPopular = true;
            cout<<"关键词类型：non-dominant, 兴趣点类型： popular!"<<endl;
            break;
        case 4:
            pivot_keyword = 25;  isPopular = false;
            cout<<"关键词类型：non-dominant, 兴趣点类型： non-popular!"<<endl;
            break;

        default:
            pivot_keyword = 20;  isPopular = true;
            cout<<"关键词类型：non-dominant, 兴趣点类型： popular!"<<endl;
            break;

    }
    //pivot_keyword2 = pivot_keyword+10;
    int KQ_num = Qn* 0.4;

    ifstream fin;
    //stringstream tt;
    //tt<<FILE_DISBOUND;
    stringstream ss;
#ifdef LasVegas
    ss<<"../../../data/LasVegas/poi_candidate/";
#endif
#ifdef Brightkite
    ss<<"../../../data/Brightkite/poi_candidate/";
#endif
#ifdef Gowalla
    ss<<"../../../data/Gowalla/poi_candidate/";
#endif
    ss<<"Q"<<query_type;
    string fileName = ss.str();

    //第一部分poi
    fin.open(fileName);
    string str;
    while(getline(fin, str)){
        istringstream tt(str);
        int query_id = -1;
        tt>>query_id;
        poi_poolset.insert(query_id);
        if(poi_poolset.size()==200) break;

    }
    fin.close();


    //从中抽取 |P_c|个POI作为candidate
    vector<int> query_objects; set<int> query_objectSet;
    for(int poi_id: poi_poolset){
        //if(curr_keywords.size()>keyword_perpoi) continue;
        query_objectSet.insert(poi_id);

        if(query_objectSet.size()==Qn) break;
    }
    for(int poi_id: query_objectSet)
        query_objects.push_back(poi_id);
    printf("*********总共加入了: %d个 兴趣点!*****************\n",query_objects.size());
    //getchar();

    return query_objects;

}


////获得包含 关键词 term的 poi
vector<int> get_term_related_objectID(int term){
    vector<int> objects;
#ifndef DiskAccess
    objects = invListOfPOI[term];
#else
    set<int> poi_poolSet;
    vector<int> related_oleaf = getObjectTermRelatedLeaves(term);
    for(int o_leaf: related_oleaf){
        vector<int> objects = getObjectTermRelatedEntry(term,o_leaf);
        for(int o: objects)
            poi_poolSet.insert(o);
    }
    for(int o: poi_poolSet)
        objects.push_back(o);
#endif
    return objects;
}

void combine_set_collection(set<int>& set, vector<int>& list){
    for(int e: list)
        set.insert(e);
}

void poiSetToList(set<int>& set, vector<int>& list){
    for(int e: set)
        list.push_back(e);
}

vector<int> generate_query_According2Popular_lv(int keyword_type, int Qn, int KQ){
    int pivot_keyword; int pivot_keyword2;
    set<int> poi_poolset;
    vector<int> pois;
    cout<<"创建"<<Qn<<"个RkGSKQ查询,";  //Qn: 25, 50, 100, 200, 400;   |KQ|: 20, 40, 60, 80, 80, 100
    switch(keyword_type){
        case Popular:
            pivot_keyword = 5;   // 1-20
#ifdef LV
            pivot_keyword = 1;
#endif
            cout<<"关键词类型：popular"<<endl;
            break;
        case Mid:
            pivot_keyword = 50;   // 50-100
#ifdef LV
            pivot_keyword = 20;
#endif
            cout<<"关键词类型：mid"<<endl;
            break;
        case Rare:
            pivot_keyword = 100;  //100以外
#ifdef LV
            pivot_keyword = 50;
#endif
            cout<<"关键词类型：rare"<<endl;
            break;
        default:
            pivot_keyword = 50;
            cout<<"关键词(默认)类型：mid"<<endl;
            break;

    }
    //pivot_keyword2 = pivot_keyword+10;
    int KQ_num = Qn* 0.4;
    pois = get_term_related_objectID(pivot_keyword);
    combine_set_collection(poi_poolset,pois);
    cout<<"加入t"<<(pivot_keyword)<<"的相关用户"<<endl;
    int i =1;
    //当候选兴趣点个数不够大时
    while(poi_poolset.size()<2*Qn){
        pois = get_term_related_objectID(pivot_keyword+i*5);
        combine_set_collection(poi_poolset,pois);
        cout<<"加入t"<<(pivot_keyword+i*5)<<"的相关用户"<<endl;
        i++;
    }
    vector<int> poi_pool;
    poiSetToList(poi_poolset,poi_pool);


    set<int> extract_pois;
    while(extract_pois.size()<Qn){
        //随机从poi_pool选一个object
        int o_random = poi_pool[rand()%poi_pool.size()];
        extract_pois.insert(o_random);
    }
    vector<int> query_objects;
    //int keyword_perpoi = 7;
    //if(Qn<30) keyword_perpoi = 4;
    //if(Qn>=100) keyword_perpoi = 10;
    for(int poi_id: extract_pois){
            POI poi;
#ifdef DiskAccess
            poi = getPOIFromO2UOrgLeafData(poi_id);
#else
            poi = POIs[poi_id];
#endif

            vector<int> curr_keywords =  poi.keywords;
            //if(curr_keywords.size()>keyword_perpoi) continue;
            query_objects.push_back(poi_id);
    }

    return query_objects;

}





vector<int> generate_query_According2Popular_big(int keyword_type, int Qn, int KQ){
    int pivot_keyword;
    int pivot_keyword2;
    vector<int> pivot_keywords;
    //int pivot_keywords[5];
    vector<set<int>> poi_pools;
    cout<<"创建"<<Qn<<"个RkGSKQ查询,";  //Qn: 25, 50, 100, 200, 400;   |KQ|: 40, 60, 80, 100, 120
    switch(keyword_type){
        case Popular:
            pivot_keyword = 5;   // 1-20
#ifdef LV
            pivot_keyword = 1;
#endif
            cout<<"关键词类型：popular"<<endl;
            break;
        case Mid:
            pivot_keyword = 50;   // 50-100
#ifdef LV
            pivot_keyword = 25;
#endif
            cout<<"关键词类型：mid"<<endl;
            break;
        case Rare:
            pivot_keyword = 100;  //100以外
#ifdef LV
            pivot_keyword = 50;
#endif
            cout<<"关键词类型：rare"<<endl;
            break;
        default:
            pivot_keyword = 50;
            cout<<"关键词(默认)类型：mid"<<endl;
            break;

    }
    pivot_keywords.push_back(pivot_keyword);
    pivot_keyword2 = pivot_keyword+5;
    pivot_keywords.push_back(pivot_keyword2);
    pivot_keywords.push_back(pivot_keyword+10);

    int KQ_num = Qn* 0.4;
    for(int term: pivot_keywords){

        set<int> poi_pool;

#ifndef DiskAccess
        vector<int> objects = invListOfPOI[term];
        for(int o: objects)
            poi_pool.insert(o);
#else
        //set<int> poi_poolSet;
        vector<int> related_oleaf = getObjectTermRelatedLeaves(term);
        for(int o_leaf: related_oleaf){
            vector<int> objects = getObjectTermRelatedEntry(term,o_leaf);
            for(int o: objects)
                poi_pool.insert(o);
        }
#endif
        poi_pools.push_back(poi_pool);

    }

    set<int> inter_object;
    set_intersection(poi_pools[0].begin(),poi_pools[0].end(),poi_pools[1].begin(), poi_pools[1].end(),inserter(inter_object,inter_object.begin()));
    cout<<"t"<<pivot_keywords[0]<<"与t"<<pivot_keywords[1]<<"相关用户个数="<<inter_object.size()<<" 具体为："<<endl;
    printSetElements(inter_object);

    set<int> inter_object2;
    set_intersection(poi_pools[2].begin(), poi_pools[2].end(),inter_object.begin(),inter_object.end(),inserter(inter_object2,inter_object2.begin()));
    cout<<"进一步与t"<<pivot_keywords[2]<<"相交后，相关用户个数="<<inter_object2.size()<<" 具体为："<<endl;
    printSetElements(inter_object2);




    set<int> group_keywords;
    int cnt = 0;
    for(int poi_id: inter_object2){
        POI poi;
#ifdef DiskAccess
        poi = getPOIFromO2UOrgLeafData(poi_id);
#else
        poi = POIs[poi_id];
#endif

        vector<int> curr_keywords =  poi.keywords;
        if(curr_keywords.size()>5) continue;
        cnt++;
        for(int term: curr_keywords)
            group_keywords.insert(term);
    }
    cout<<"只留下关键词个数较小的对象, 剩余="<<cnt;
    cout<<"他们关键词的并集,size="<<group_keywords.size()<<endl;
    printSetElements(group_keywords);



    vector<int> query_objects;
    /*
    set<int> extract_pois;
    while(extract_pois.size()<Qn){
        //随机从poi_pool选一个object
        int o_random = poi_pool[rand()%poi_pool.size()];
        extract_pois.insert(o_random);
    }

    for(int o: extract_pois)
        query_objects.push_back(o);*/

    return query_objects;

}



vector<int> obtain_query_keyword(int poi){
    vector<int> candidateKeywords;
    set<int> _keyset;
    POI p;
#ifdef DiskAccess
    p = getPOIFromO2UOrgLeafData(poi);
#else
    p = POIs[poi];
#endif
    for(int term:p.keywords)
        _keyset.insert(term);
    for(int term2:_keyset){
        candidateKeywords.push_back(term2);
        if(candidateKeywords.size()==DEFAULT_KEYWORDSIZE) break;
    }

    printElements(candidateKeywords);
}





vector<int> generate_keywords_pool(int KQ_size, set<int>& poi_collection, int frequent_thethod, int Q_size){
    vector <int > chainStores;
    int test =0;

    map<int, vector<int>>::iterator iter;
    iter = invListOfPOI.begin();
    vector<int> key_pool;
    //cout<<"KQ_size="<<KQ_size<<endl;
    int size2;

    while(iter != invListOfPOI.end()){
        int term = iter->first;
        int size =  iter->second.size();
        //cout<<"t"<<term<<",size="<<size<<endl;
        if(size > DEFAULT_THRETHEHOLD && size < DEFAULT_THRETHEHOLD_MAX){
            key_pool.push_back(iter->first);
            //cout<<"key_pool_size="<<key_pool.size()<<endl;

            vector<int> poi_list = invListOfPOI[term];
            for(int poi: poi_list){
                poi_collection.insert(poi);
            }

        }
        size2 = key_pool.size();
        if(key_pool.size()> 10*KQ_size && poi_collection.size()> 10*Q_size){
            cout<<"key_pool_size="<<key_pool.size()<<endl;
            cout<<"break"<<endl;
            break;
        }

        iter++;

    }
    //getchar();
    if(key_pool.size()>KQ_size && poi_collection.size()>Q_size){
        return key_pool;
    }
    else{
        if(key_pool.size()<KQ_size)
            cout<<"关键词个数不够, 请适当减小frequent_thethod！"<<endl;
        else
            cout<<"候选兴趣点不够，请适当减小frequent_thethod！"<<endl;
        exit(1);
    }

}


vector<int> generate_keywords(vector<int> keywords_pool, int KQ_size){
    if(keywords_pool.size()<KQ_size){
        cout<<"关键词候选集不够！"<<endl;
    }
    vector<int> keywords;
    for(int i=0;i<KQ_size;i++){
        keywords.push_back(keywords_pool[i]);
    }
    return  keywords;
}

vector<int> generate_location(set<int>& poi_pool, vector<int>& checkIn_usr, int Q_size){

    vector <int > chainStores;


    //extract |Q| query objects 构建查询对象位置
    if(poi_pool.size()> Q_size){
        for(int o: poi_pool){
            chainStores.push_back(o);
            if(chainStores.size()==Q_size)
                break;
        }
    }
    //retain some check-in records
    for(int o: chainStores){
        if(poiCheckInIDList.count(o)==1){
            vector<int> usrs = poiCheckInIDList[o];
            for(int u: usrs){  //这里可以适当控制check in 个数
                checkIn_usr.push_back(u);
                if(checkIn_usr.size()>100)
                    break;
            }

        }
    }

    return  chainStores;
}


vector<int> query_generation(vector<int>& keyword_pool, vector<int>& KQ, vector<int>& checkIn_usr, int frequent_thethod, int KQ_size, int Q_size){
    set<int> poi_pool;
    vector <int > chainStores;
    //获取高频词汇
    keyword_pool = generate_keywords_pool(KQ_size, poi_pool, frequent_thethod, Q_size);


    //构建查询对象位置
    chainStores = generate_location(poi_pool, checkIn_usr, Q_size);

    //构建查询关键词
    KQ = generate_keywords(keyword_pool, KQ_size);

    return  chainStores;
}



double testTopkUser(int usr_id,int Qk, float a, float alpha){
    query user0;
    user0 = transformToQ(Users[usr_id]);
    //TopkQueryResult qr= topkSDijkstra(user0, Qk, a, alpha);
    TopkQueryResult qr= topkSDijkstra_memory(user0, Qk, a, alpha);
    //cout<<qr.topkScore<<endl;
    //cout<<"usr_id"<<user0.id<<endl;
    //qr.printResults(usr_id, Qk);
    //cout<<topkSDijkstra_verify_usr(user0,Qk,a,alpha,10000000000000.0).topkScore<<endl;
    return qr.topkScore;

}

void testAllTopk(int Qk, float a, float alpha){
    int usr =Users.size();
    //设置线程数，一般设置的线程数不超过CPU核心数，这里开4个线程执行并行代码段
    omp_set_num_threads(4);
    //#pragma omp parallel for
    for(int id=0;id<usr;id++){
        cout<<"usr"<<id<<",当前处理线程id:"<<omp_get_thread_num()<<endl;
        testTopkUser(id, Qk,a, alpha);
    }
}



void testObjectTermRelatedEntry_all_coverage(){



    set<int> queryKeys;
    //queryKeys.insert(1);
    queryKeys.insert(5);
    //test access node
    for(int i=0;i<GTree.size();i++){
        //i=382;
        set<int> coverdKeys; set<int> usr_coverdKeys;
        int node_id = i; //20
        if(GTree[i].isleaf)
            cout<<"---------------------for leaf"<<i<<"-----------------------"<<endl;
        else
            cout<<"---------------------for node"<<i<<"-----------------------"<<endl;

        TreeNode node;
        //测试object
        cout<<"测试object node"<<endl;

        node = getGIMTreeNodeData(node_id,OnlyO);


        //测试object
        for(int keyword: node.objectUKeySet){
            cout<<"for oterm"<<keyword<<endl;
            cout<<"term_weight="<<getTermIDFWeight(keyword)<<endl;
            vector<int> related = getObjectTermRelatedEntry(keyword,node_id);
            cout<<"O node "<<node_id<<"关键词"<<keyword<<"的后继:";
            printElements(related);

        }

        if(GTree[i].isleaf)
            cout<<"---------------------finish leaf"<<i<<"-----------------------"<<endl;
        else
            cout<<"---------------------finish node"<<i<<"-----------------------"<<endl;
        //getchar();


    }



}


void testObjectTermRelatedEntry_all_self(){


    //test access node
    for(int i=0;i<GTree.size();i++){
        //i=382;
        set<int> coverdKeys; set<int> usr_coverdKeys;
        int node_id = i; //20
                string ss;
        if(GTree[i].isleaf){
            cout<<"---------------------for leaf"<<i<<"-----------------------"<<endl;
            ss="O leaf";

        }
        else{
            cout<<"---------------------for node"<<i<<"-----------------------"<<endl;
            ss="O node";
        }

        TreeNode node;
        //测试object
        cout<<"测试object node"<<endl;

        node = getGIMTreeNodeData(node_id,OnlyO);
        cout<<"node object keyword size:"<<node.objectUKeySet.size()<<endl;

        //测试object
        for(int keyword: node.objectUKeySet){
            cout<<"for oterm"<<keyword<<endl;
            cout<<"term_weight="<<getTermIDFWeight(keyword)<<endl;
#ifdef DEBUG
            if(i==0 && keyword==2){
                cout<<"find term 2 in node 0 !"<<endl;
            }
#endif
            vector<int> related = getObjectTermRelatedEntry(keyword,node_id);
            cout<<ss<<node_id<<" 关键词"<<keyword<<"的后继:";
            printElements(related);
            if(related.size()>50){
                getchar();
            }

        }

        if(GTree[i].isleaf)
            cout<<"---------------------finish leaf"<<i<<"-----------------------"<<endl;
        else
            cout<<"---------------------finish node"<<i<<"-----------------------"<<endl;

        if(i%100==0){
            cout<<"node"<<i<<",Enter..."<<endl;
            getchar();
        }


    }

    cout<<"testObjectTermRelatedEntry_all_self! Enter..."<<endl;
    //getchar();



}



void testObjectTermRelatedEntry_specific(int term, int node){

#ifdef DiskAccess
    TEST_START
    vector<int> related;
    related= getObjectTermRelatedEntry(term,node);
    cout<<"o node"<<node<<" 在t"<<term<<"下的后继：";
    printElements(related);
    TEST_END
    TEST_DURA_PRINT(getObjectTermRelatedEntry)

#else
    set<int> relatedSet;
    relatedSet= GTree[node].inverted_list_o[term];
    printSetElements(relatedSet);
#endif

    cout<<"finnish test term"<<term<<"object node"<<node<<endl;

}

void testUserTermRelatedEntry_specific(int term, int node){


    vector<int> related = getUsrTermRelatedEntry(term,node);
    printElements(related);



}


void testGetGIMTreeNode(int node_id, int pattern){
    TreeNode n = getGIMTreeNodeData(node_id,pattern);
}

void testAccessPOI(){
    for(int id=0;id<poi_num;id++){
        cout<<"检索兴趣点p"<<id<<endl;
        POI p = getPOIFromO2UOrgLeafData(id);
        printPOIInfo(p);
        if(id%100==0){
            getchar();
        }
    }
}


void testUserTermRelatedEntry_all(){


    set<int> queryKeys;
    //queryKeys.insert(1);
    queryKeys.insert(1);
    //test access node
    for(int i=0;i<GTree.size();i++){

        set<int> coverdKeys; set<int> usr_coverdKeys;
        int node_id = i; //20
        string ss;
        if(GTree[i].isleaf){
            cout<<"---------------------for leaf"<<i<<"-----------------------"<<endl;
            ss="U leaf";

        }
        else{
            cout<<"---------------------for node"<<i<<"-----------------------"<<endl;
            ss="U node";
        }
        TreeNode node;

        node = getGIMTreeNodeData(node_id,OnlyU);
        if(GTree[i].isleaf)
            cout<<"user leaf's union keywords:";
        else
            cout<<"user node's union keywords:";
         printSetElements(node.userUKeySet);

        //测试usr
        for(int keyword: node.userUKeySet){
            cout<<"for uterm"<<keyword<<endl;
            vector<int> related = getUsrTermRelatedEntry(keyword,node_id);
            cout<<ss<<node_id<<"关键词"<<keyword<<"的后继条目:";
            printElements(related);
            cout<<"具体为："<<endl;

            vector<int> related_leaves = getUsrTermRelatedLeafNode(keyword,node_id);
            cout<<ss<<node_id<<"关键词"<<keyword<<"的后继叶节点:";
            printElements(related_leaves);
            cout<<"具体为："<<endl;

        }


        if(GTree[i].isleaf)
            cout<<"---------------------finish leaf"<<i<<"-----------------------"<<endl;
        else
            cout<<"---------------------finish node"<<i<<"-----------------------"<<endl;

        if(i%100==0){
            cout<<"node"<<i<<" ,Enter..."<<endl;
            getchar();
        }

    }

    cout<<"testUserTermRelatedEntry_all! Enter..."<<endl;
    getchar();



}



void testTackNodePathForUser(int usr){
    User u = Users[usr];
    int loc = u.Ni;
    int usr_leaf = Nodes[loc].gtreepath.back();
    cout<<"leaf "<<usr_leaf<<endl;
    int current = GTree[usr_leaf].father;
    while(current!=0){
        cout<<"node "<<current<<endl;
        current = GTree[current].father;
    }
}


void testTackNodePathForLeaf(int leaf_node){

    cout<<"leaf "<<leaf_node<<endl;
    int current = GTree[leaf_node].father;
    while(current!=0){
        cout<<"node "<<current<<endl;
        current = GTree[current].father;
    }
}


//测试GIMtree节点地址获取函数
void testGetGTreeNodeAddr(){
    for(int i=0;i<GTree.size();i++){
        int addr2 = getGTreeNodeAddr(i);
        if(true){ //
            cout<<"for node"<<i;
            cout<<",addr2="<<addr2<<endl;
        }

    }
    cout<<"finish testGetGTreeNodeAddr"<<endl;
    //getchar();
}

//测试GIMtree节点数据内容获取函数
void testGetGIMTreeNodeData(){
    for(int i=0;i<GTree.size();i++){

        TreeNode od = getGIMTreeNodeData(i,OnlyO);
        printONodeInfo(od,i);

    }
    cout<<"object: finish testGetGIMTreeNode data, Enter... "<<endl;
    getchar();

    for(int i=0;i<GTree.size();i++){

        TreeNode ud = getGIMTreeNodeData(i,OnlyU);
        printUNodeInfo(ud,i);

    }
    cout<<"user: finish testGetGIMTreeNode data "<<endl;

    cout<<"finish all procedures! Enter..."<<endl;
    getchar();
}




//测试邻边的关键词获取
void testFunctionForEdge(){
    for(int Ni=0;Ni<VertexNum;Ni++){
        int Nk = adjList[Ni][0];
        cout<<"-------------对edge("<<Ni<<","<<Nk<<") 上数据获取的功能函数进行测试-------------"<<endl;
        //测试读取边的关键词信息
        set<int> OKeyset;
        set<int> keys1 = EdgeMap[getKey(Ni,Nk)].OUnionKeys;
        cout<<"edge("<<Ni<<","<<Nk<<") 's okey=";
        printSetElements(keys1);
        OKeyset = getOKeyEdge(Ni,0);
        cout<<"read edge("<<Ni<<","<<Nk<<") 's okey=";
        printSetElements(OKeyset);
        cout<<"edge keyword content test over!"<<endl;

        if(!OKeyset.size()>0) continue;

        //测试读取边上兴趣点个数
        int poi_size = EdgeMap[getKey(Ni,Nk)].pts.size();
        cout<<"边上实际有"<<poi_size<<"个兴趣点"<<endl;
        int read_poi_size =  getPOISizeOfEdge(Ni, 0);
        cout<<"读取边上有"<<read_poi_size<<"个兴趣点"<<endl;

        //测试边上兴趣点数据
        for(int i=0;i<poi_size;i++){
            //int idx_baseAddr = getPOIIdxAddrOfEdge(Ni,0);
            POI p = getPOIDataOnEdge(Ni,0, i);
            printPOIInfo(p);
        }

        cout<<"-----------------------对edge("<<Ni<<","<<Nk<<") 的测试结束-----------------------"<<endl;
        getchar();
    }
}


void testTermWeight(){
    for(int i=1;i<15;i++){
        int term = i;
        float term_weight = getTermIDFWeight(term);
        cout<<"term"<<term<<" idf_score="<<term_weight<<endl;
    }
}

void testTermInvInfo(){
    int term_size = getTermSize();
    for(int i=2;i<term_size;i++){
        int term = i;
        float term_weight = getTermIDFWeight(term);
        cout<<"----------------------term"<<term<<" idf_score="<<term_weight<<"----------------------"<<endl;
        vector<int> posting_list = getTermOjectInvList(term);
        cout<<"object posting list size="<<posting_list.size()<<endl;
        printElements(posting_list);

        vector<int> user_posting_list = getTermUserInvList(term);
        cout<<"user posting list size="<<user_posting_list.size()<<endl;
        //printElements(user_posting_list);
        return;
    }
}

void testSPSP(){
    while(true){
        TIME_TICK_START
        int v = random()%VertexNum;
        int u = random()%VertexNum;
        int dist = SPSP(v,u);
        cout<<"dist(v"<<v<<",v"<<u<<")="<<dist<<endl;

        TIME_TICK_END
        TIME_TICK_PRINT_us("SPSP runtime: ")
        getchar();
    }
}


void test_topk_original(int k, float alpha){
    cout<<"进行test_topk()测试"<<endl;
    for(int i=1;i<UserID_MaxKey;i++){
        cout<<"run top-k for u"<<i<<endl;
        TIME_TICK_START
        User u;
        double score_topk = -1;
#ifndef DiskAccess
        u = Users[i];
    #ifdef LV
        score_topk = topkSDijkstra_verify_usr_memory(transformToQ(u),k,DEFAULT_A,alpha,1000,0.0).topkScore;
    #else
        priority_queue<Result> results = TkGSKQ_bottom2up_verify_memory(i,k,DEFAULT_A,alpha,1000,0.0);
        score_topk = results.top().score;
    #endif
#else
        u = getUserFromO2UOrgLeafData(i);
        printUsrInfo(u);

#ifdef LV
        score_topk = topkSDijkstra_verify_usr_disk(transformToQ(u),k,DEFAULT_A,alpha,1000,0.0).topkScore;
#else
        priority_queue<Result> results = TkGSKQ_bottom2up_verify_disk(i,k,DEFAULT_A,alpha,1000,0.0);
        score_topk = results.top().score;
#endif
        //cout<<"u"<<i<<": score_k2="<<results.top().score;
#endif

        cout<<"u"<<i<<": score_topk="<<score_topk;
        TIME_TICK_END
        TIME_TICK_PRINT(" runtime: ")

        cout<<endl;
        getchar();
    }
    cout<<"test_topk()测试结束"<<endl;
    //getchar();
}

void test_u2o(int u_id, int o_id){
    User u = getUserFromO2UOrgLeafData(u_id);
    POI p = getPOIFromO2UOrgLeafData(o_id);
    cout<<"u.Ni="<<u.Ni<<", p.Ni="<<p.Ni<<endl;
    cout<<"u.Nj="<<u.Nj<<", p.Nj="<<p.Nj<<endl;
    //u.Ni, p.Ni
    double dist_v2v = SPSP(u.Ni,p.Ni)/1000.0; double u_dis = u.dis; double p_dis = p.dis;
    double total_dist = dist_v2v + u_dis + p_dis;
    cout<<"dist_total="<<total_dist<<", dist(u.Ni, p.Ni)="<<dist_v2v<<",u_dis="<<u_dis<<",p_dis="<<p_dis<<endl;

    //u.Ni, p.Nj
    dist_v2v = SPSP(u.Ni,p.Nj)/1000.0;
    p_dis = p.dist - p.dis;
    total_dist = dist_v2v + u_dis + p_dis; double total_tmp2 = total_dist;
    cout<<"dist_total="<<total_dist<<", dist(u.Ni, p.Nj)="<<dist_v2v<<",u_dis="<<u_dis<<",p_dis="<<p_dis<<endl;

    cout<<"------------------min from u.Ni="<<min(total_tmp2, total_dist)<<"-----------------------"<<endl;




    //u.Nj, p.Ni
    dist_v2v = SPSP(u.Nj,p.Ni)/1000.0;
    u_dis = u.dist - u.dis; p_dis = p.dis;
    total_dist = dist_v2v + u_dis + p_dis;  double total_tmp3 = total_dist;
    cout<<"dist_total="<<total_dist<<", dist(u.Nj, p.Ni)="<<dist_v2v<<",u_dis="<<u_dis<<",p_dis="<<p_dis<<endl;

    //u.Nj, p.Nj
    dist_v2v = SPSP(u.Nj,p.Nj)/1000.0;
    u_dis = u.dist - u.dis; p_dis = p.dist - p.dis;
    total_dist = dist_v2v + u_dis + p_dis;  double total_tmp4 = total_dist;
    cout<<"dist_total="<<total_dist<<", dist(u.Nj, p.Nj)="<<dist_v2v<<",u_dis="<<u_dis<<",p_dis="<<p_dis<<endl;
    cout<<"------------------min from u.Nj="<<min(total_tmp3, total_tmp4)<<"-----------------------"<<endl;
}


void test_single_topk(int u_id, int k, float alpha, int method){
    cout<<"进行test_single_topk()测试"<<endl;
    clock_t  startTime, endTime;
        cout<<"run top-k for u"<<u_id<<endl;
        User u;
        double score_topk = -1;


        u = getUserFromO2UOrgLeafData(u_id);
        printUsrInfo(u);


        //string method = "gtree";
        if(method ==  GTree_topk ){
            TIME_TICK_START
            priority_queue<Result> results3 = TkGSKQ_bottom2up_verify_disk(u_id,k,DEFAULT_A,alpha,1000,0.0);
            score_topk = results3.top().score;
            cout<<"u"<<u_id<<": score_topk(g-tree-Ni&Nj)="<<score_topk<<"具体为："<<endl;
            //cout<<"(ni.leaf="<<Nodes[u.Ni].gtreepath.back()<<",nj.leaf="<<Nodes[u.Nj].gtreepath.back()<<")"<<endl;
            printTopkResultsQueue(results3);
            //cout<<"upper_jumper="<<upper_jumper<<endl;
            upper_jumper =0;
            TIME_TICK_END
            TIME_TICK_PRINT(" runtime(gtree-Ni&Nj) runtime")
            cout<<endl;
        }
        else if(method == DJ_topk){
            TIME_TICK_START
            TopkQueryCurrentResult topk_r4 = topkSDijkstra_verify_usr_disk(transformToQ(u),k,DEFAULT_A,alpha,1000,0.0); //topkSDijkstra_verify_usr_disk
            score_topk = topk_r4.topkScore;
            cout<<"u"<<u_id<<": score_topk(dj-Ni&Nj)="<<score_topk<<"具体为："<<endl;
            cout<<"(ni.leaf="<<Nodes[u.Ni].gtreepath.back()<<",nj.leaf="<<Nodes[u.Nj].gtreepath.back()<<")"<<endl;
            topk_r4.printResults(k);
            TIME_TICK_END
            TIME_TICK_PRINT(" topkdjstra(dj-Ni&Nj) runtime: ")
            cout<<endl;
        }
        else if(method == BruteForce_topk){
            TIME_TICK_START
            TopkQueryCurrentResult topk_r5 = topkBruteForce_disk(u,k,DEFAULT_A,alpha,1000,0.0);
            score_topk = topk_r5.topkScore;
            cout<<"u"<<u_id<<": score_topk(bruteForce)="<<score_topk<<"具体为："<<endl;
            cout<<"(ni.leaf="<<Nodes[u.Ni].gtreepath.back()<<",nj.leaf="<<Nodes[u.Nj].gtreepath.back()<<")"<<endl;
            topk_r5.printResults(k);
            TIME_TICK_END
            TIME_TICK_PRINT(" topkBruteForce runtime: ")
            //cout<<endl;
        }
        else if(method == NVD_topk){

            Load_Hybrid_NVDG_AddressIdx_fast();
            //TIME_TICK_START
            startTime = clock();

            TopkQueryCurrentResult topk_r6 = TkGSKQ_NVD(u, k, DEFAULT_A,alpha);//,poiMark,poiADJMark);

            endTime = clock();
            score_topk = topk_r6.topkScore;
            cout<<"u"<<u_id<<": score_topk(nvd)="<<score_topk<<"具体为："<<endl;
            cout<<"(ni.leaf="<<Nodes[u.Ni].gtreepath.back()<<",nj.leaf="<<Nodes[u.Nj].gtreepath.back()<<")"<<endl;
            topk_r6.printResults(k);
            cout<<"topkNVD runtime: "<<(double)(endTime - startTime)/CLOCKS_PER_SEC * 1000  << "ms" << endl;

        }


    cout<<"------------------------------------------------------------------------"<<endl;
    cout<<"test_topk()测试结束"<<endl;
    //getchar();
}

void test_checkUserWhere(int u_id){
    User u = getUserFromO2UOrgLeafData(u_id);
    int u_Ni = u.Ni;
    int u_leaf = Nodes[u_Ni].gtreepath.back();
    cout<<"u"<<u.id<<" in leaf"<<u_leaf<<endl;
    int u_node = u_leaf; int root = 0;
    while(u_node!=root){
        u_node = GTree[u_node].father;
        cout<<"in n"<<u_node<<endl;
    }
    cout<<"结束！"<<endl;

}

void test_topk(int k, float alpha, int method){
    // = "dj";
    cout<<"进行test_topk()测试"<<endl;

    clock_t  startTime, endTime; double total_Time = 0;

    for(int i=65;i<100;i++){
        //int i= u_id;
        cout<<"run top-k for u"<<i<<endl;
        User u;
        double score_topk = -1;
        u = getUserFromO2UOrgLeafData(i);

        bool* poiMark = new bool[poi_num];
        bool* poiADJMark = new bool[poi_num];
        for(int i=0;i<poi_num;i++){
            poiMark[i]=false;
            poiADJMark[i] = false;
        }


        startTime = clock();
        if(method == BruteForce_topk){
            //TIME_TICK_START
            //TIME_TICK_START
            TopkQueryCurrentResult topk_r5 = topkBruteForce_disk(u,k,DEFAULT_A,alpha,1000,0.0);
            score_topk = topk_r5.topkScore;
            cout<<"u"<<i<<": score_topk(bruteForce)="<<score_topk<<"具体为："<<endl;
            cout<<"(ni.leaf="<<Nodes[u.Ni].gtreepath.back()<<",nj.leaf="<<Nodes[u.Nj].gtreepath.back()<<")"<<endl;
            topk_r5.printResults(k);
            //TIME_TICK_END
            //TIME_TICK_PRINT(" topkBruteForce runtime: ")
            cout<<endl;
        }
        else if(method == DJ_topk){
            //TIME_TICK_START
            TopkQueryCurrentResult topk_r4 = topkSDijkstra_verify_usr_disk(transformToQ(u),k,DEFAULT_A,alpha,1000,0.0); //topkSDijkstra_verify_usr_disk
            score_topk = topk_r4.topkScore;
            cout<<"u"<<i<<": score_topk(dj-Ni&Nj)="<<score_topk<<"具体为："<<endl;
            cout<<"(ni.leaf="<<Nodes[u.Ni].gtreepath.back()<<",nj.leaf="<<Nodes[u.Nj].gtreepath.back()<<")"<<endl;
            topk_r4.printResults(k);
            //TIME_TICK_END
            //TIME_TICK_PRINT(" topkdjstra(dj-Ni&Nj) runtime: ")
            cout<<endl;
        }
        else if(method == GTree_topk){

            priority_queue<Result> results3 = TkGSKQ_bottom2up_verify_disk(i,k,DEFAULT_A,alpha,1000,0.0);

            score_topk = results3.top().score;
            cout<<"u"<<i<<": score_topk(g-tree-Ni&Nj)="<<score_topk<<"具体为："<<endl;
            //cout<<"(ni.leaf="<<Nodes[u.Ni].gtreepath.back()<<",nj.leaf="<<Nodes[u.Nj].gtreepath.back()<<")"<<endl;
            printTopkResultsQueue(results3);
            //cout<<"upper_jumper="<<upper_jumper<<endl;
            upper_jumper =0;
            //TIME_TICK_END
            //TIME_TICK_PRINT(" runtime(gtree-Ni&Nj) runtime")
            cout<<endl;
        }
        else if(method == NVD_topk){
            //TIME_TICK_START
#ifdef HybridHash
            TopkQueryCurrentResult topk_r6 = TkGSKQ_NVD(u, k, DEFAULT_A,alpha);//,poiMark,poiADJMark);
#else
            TopkQueryCurrentResult topk_r6 = TkGSKQ_NVD_vertexOnly(u, k, DEFAULT_A,alpha,poiMark,poiADJMark);
#endif
            score_topk = topk_r6.topkScore;
            cout<<"u"<<i<<": score_topk(nvd)="<<score_topk<<"具体为："<<endl;
            cout<<"(ni.leaf="<<Nodes[u.Ni].gtreepath.back()<<",nj.leaf="<<Nodes[u.Nj].gtreepath.back()<<")"<<endl;
            topk_r6.printResults(k);

            //TIME_TICK_END
            //TIME_TICK_PRINT(" topkNVD runtime: ")
        }
        endTime = clock();
        double runTime = (double)(endTime - startTime)/CLOCKS_PER_SEC * 1000;
        total_Time += runTime;
        cout<<"topk-"<<method<<" runtime: "<<runTime<< "ms" << endl;

        delete []poiMark; delete []poiADJMark;
        cout<<"------------------------------------------------------------------------"<<endl;
    }



    cout<<"*****************test_topk()测试结束, average runTime="<<total_Time/(100.0-65)<<"ms**************************"<<endl;
    //getchar();
}

void test_ComparationMethods_topk(int k, float alpha){
    string method = "gtree";
    cout<<"进行test_topk()测试"<<endl;

    clock_t  startTime, endTime; double total_Time = 0;
    Load_Hybrid_NVDG_AddressIdx_fast();

    for(int i=24;i<100;i++){

        cout<<"run top-k for u"<<i<<endl;

        bool* poiMark = new bool[poi_num];
        bool* poiADJMark = new bool[poi_num];
        for(int i=0;i<poi_num;i++){
            poiMark[poi_num]=false;
            poiADJMark[poi_num] = false;
        }

        User u;
        double score_topk = -1;
        u = getUserFromO2UOrgLeafData(i);


        startTime = clock();

        /*TIME_TICK_START
        TopkQueryCurrentResult topk_r5 = topkBruteForce_disk(u,k,DEFAULT_A,alpha,1000,0.0);
        score_topk = topk_r5.topkScore;
        cout<<"u"<<i<<": score_topk(bruteForce)="<<score_topk<<"具体为："<<endl;
        cout<<"(ni.leaf="<<Nodes[u.Ni].gtreepath.back()<<",nj.leaf="<<Nodes[u.Nj].gtreepath.back()<<")"<<endl;
        topk_r5.printResults(k);
        TIME_TICK_END
        TIME_TICK_PRINT(" topkBruteForce runtime: ")*/

        TIME_TICK_START
        TopkQueryCurrentResult topk_r4 = topkSDijkstra_verify_usr_disk(transformToQ(u),k,DEFAULT_A,alpha,1000,0.0); //topkSDijkstra_verify_usr_disk
        score_topk = topk_r4.topkScore;
        cout<<"----------------u"<<i<<": score_topk(dj-Ni&Nj)="<<score_topk<<"具体为："<<endl;
        cout<<"(ni.leaf="<<Nodes[u.Ni].gtreepath.back()<<",nj.leaf="<<Nodes[u.Nj].gtreepath.back()<<")"<<endl;
        topk_r4.printResults(k);
        TIME_TICK_END
        TIME_TICK_PRINT(" topkdjstra(dj-Ni&Nj) runtime: ")

        TIME_TICK_START
        priority_queue<Result> results3 = TkGSKQ_bottom2up_verify_disk(i,k,DEFAULT_A,alpha,1000,0.0);
        score_topk = results3.top().score;
        cout<<"----------------u"<<i<<": score_topk(g-tree-Ni&Nj)="<<score_topk<<"具体为："<<endl;
        printTopkResultsQueue(results3);
        TIME_TICK_END
        TIME_TICK_PRINT(" TkGSKQ_bottom2up_verify_disk runtime: ")



        TIME_TICK_START
        TopkQueryCurrentResult topk_r6 = TkGSKQ_NVD(u, k, DEFAULT_A,alpha);//,poiMark,poiADJMark);
        score_topk = topk_r6.topkScore;
        cout<<"----------------u"<<i<<": score_topk(nvd)="<<score_topk<<"具体为："<<endl;
        cout<<"(ni.leaf="<<Nodes[u.Ni].gtreepath.back()<<",nj.leaf="<<Nodes[u.Nj].gtreepath.back()<<")"<<endl;
        topk_r6.printResults(k);
        TIME_TICK_END
        TIME_TICK_PRINT(" TkGSKQ_NVD runtime: ")


        endTime = clock();
        double runTime = (double)(endTime - startTime)/CLOCKS_PER_SEC * 1000;
        total_Time += runTime;
        cout<<"topk-"<<method<<" runtime: "<<runTime<< "ms" << endl;

        delete  []poiADJMark; delete [] poiMark;
        cout<<"------------------------------------------------------------------------"<<endl;
        cout<<"*****************test_topk()测试结束, average runTime="<<total_Time/(100.0-65)<<"ms**************************"<<endl;
        getchar();

    }




}

void poi_selection_heuristic(int k, float alpha){

    Load_Hybrid_NVDG_AddressIdx_fast();
    map<int, set<int>> poi_PotentialUserMap;
    for(int i=0;i<reviewed_userNum;i++){  //reviewed_userNum

        cout<<"run top-k for u"<<i<<endl;
        if(i==25){
            cout<<"25!"<<endl;
        }
        User u;
        double score_topk = -1;
        u = getUserFromO2UOrgLeafData(i);
        TopkQueryCurrentResult topk_r6 = TkGSKQ_NVD(u, k, DEFAULT_A,alpha);//,poiMark,poiADJMark);
        score_topk = topk_r6.topkScore;
        topk_r6.topkElements;
        cout<<"----------------u"<<i<<": score_topk(nvd)="<<score_topk<<"具体为："<<endl;

        for(int j=0;j<topk_r6.topkElements.size();j++){
            Result rr = topk_r6.topkElements[j];
            int p_id = rr.id; int u_id = i;
            poi_PotentialUserMap[p_id].insert(i);
        }

    }
    int _size0 = poi_PotentialUserMap[526].size();
    cout<<"poi_PotentialUserMap[526].size="<<_size0<<",具体为"<<endl;
    for(int _id:poi_PotentialUserMap[526])
        cout<<"u"<<_id<<endl;
    //getchar();

    ////top-k计算完毕， 开始为 q(1~4)选兴趣点
    int pivot_keyword  =-1;
    int lower_size = -1; int upper_size = -1;
    for(int query_type=1; query_type<=4;query_type++){
        set<int> poi_poolset;
        switch(query_type){
            case 1: //group 1 (popular POI with dominant keywords
                pivot_keyword = 1; lower_size = 15; upper_size = 18;
                cout<<"关键词类型：dominant, 兴趣点类型： popular!"<<endl;
                break;
            case 2:
                pivot_keyword = 2; lower_size = 12; upper_size = 15;
                cout<<"关键词类型：dominant, 兴趣点类型： non-popular!"<<endl;
                break;
            case 3:
                pivot_keyword = 20; lower_size = 10; upper_size = 12;
                cout<<"关键词类型：non-dominant, 兴趣点类型： popular!"<<endl;
                break;
            case 4:
                pivot_keyword = 25; lower_size = 8; upper_size = 10;
                cout<<"关键词类型：non-dominant, 兴趣点类型： non-popular!"<<endl;
                break;

            default:
                pivot_keyword = 20;
                cout<<"关键词类型：non-dominant, 兴趣点类型： popular!"<<endl;
                break;

        }
        bool flag = false;
        for(int target_term = pivot_keyword; target_term< pivot_keyword+5;target_term++){
            cout<<"加入t"<<(target_term)<<"的相关用户"<<endl;
            vector<int> current_usersList = getTermOjectInvList(target_term);
            for(int poi_id: current_usersList){
                int _size = poi_PotentialUserMap[poi_id].size();
                if(_size>0) //if(_size>lower_size&& _size<=upper_size)
                    poi_poolset.insert(poi_id);
                if(poi_poolset.size()==200)
                    flag = true;
            }
            if(flag==true)
                break;
        }
        stringstream ss;
#ifdef LasVegas
        ss<<"../../../data/LasVegas/poi_candidate/";
#endif
#ifdef Brightkite
        ss<<"../../../data/Brightkite/poi_candidate/";
#endif
#ifdef Gowalla
        ss<<"../../../data/Gowalla/poi_candidate/";
#endif


        ss<<"Q"<<query_type;

        string fileName = ss.str();
        ofstream fout;
        fout.open(fileName.c_str());
        for(int selected_id:poi_poolset){
            fout<<selected_id<<endl;
        }
        fout.close();

    }




}


void test_ComparationMethods_topk_single(int u_id, int k, float alpha){
    string method = "gtree";
    cout<<"进行test_topk()测试"<<endl;

    clock_t  startTime, endTime; double total_Time = 0;

    //for(int i=65;i<100;i++){
        int i = u_id;
        cout<<"run top-k for u"<<i<<endl;

        bool* poiMark = new bool[poi_num];
        bool* poiADJMark = new bool[poi_num];
        for(int i=0;i<poi_num;i++){
            poiMark[poi_num]=false;
            poiADJMark[poi_num] = false;
        }

        User u;
        double score_topk = -1;
        u = getUserFromO2UOrgLeafData(i);


        startTime = clock();

        /*TIME_TICK_START
        TopkQueryCurrentResult topk_r5 = topkBruteForce_disk(u,k,DEFAULT_A,alpha,1000,0.0);
        score_topk = topk_r5.topkScore;
        cout<<"u"<<i<<": score_topk(bruteForce)="<<score_topk<<"具体为："<<endl;
        cout<<"(ni.leaf="<<Nodes[u.Ni].gtreepath.back()<<",nj.leaf="<<Nodes[u.Nj].gtreepath.back()<<")"<<endl;
        topk_r5.printResults(k);
        TIME_TICK_END
        TIME_TICK_PRINT(" topkBruteForce runtime: ")

        TIME_TICK_START
        TopkQueryCurrentResult topk_r4 = topkSDijkstra_verify_usr_disk(transformToQ(u),k,DEFAULT_A,alpha,1000,0.0); //topkSDijkstra_verify_usr_disk
        score_topk = topk_r4.topkScore;
        cout<<"u"<<i<<": score_topk(dj-Ni&Nj)="<<score_topk<<"具体为："<<endl;
        cout<<"(ni.leaf="<<Nodes[u.Ni].gtreepath.back()<<",nj.leaf="<<Nodes[u.Nj].gtreepath.back()<<")"<<endl;
        topk_r4.printResults(k);
        TIME_TICK_END
        TIME_TICK_PRINT(" topkdjstra(dj-Ni&Nj) runtime: ")*/

        /*TIME_TICK_START
        priority_queue<Result> results3 = TkGSKQ_bottom2up_verify_disk(i,k,DEFAULT_A,alpha,1000,0.0);
        score_topk = results3.top().score;
        cout<<"u"<<i<<": score_topk(g-tree-Ni&Nj)="<<score_topk<<"具体为："<<endl;
        printTopkResultsQueue(results3);
        TIME_TICK_END
        TIME_TICK_PRINT(" TkGSKQ_bottom2up_verify_disk runtime: ")*/


        TIME_TICK_START
        Load_Hybrid_NVDG_AddressIdx_fast();
#ifdef HybridHash
        TopkQueryCurrentResult topk_r6 = TkGSKQ_NVD(u, k, DEFAULT_A,alpha);//,poiMark,poiADJMark);
#else
        TopkQueryCurrentResult topk_r6 = TkGSKQ_NVD_vertexOnly(u, k, DEFAULT_A,alpha,poiMark,poiADJMark);
#endif

        score_topk = topk_r6.topkScore;
        cout<<"u"<<i<<": score_topk(nvd)="<<score_topk<<"具体为："<<endl;
        cout<<"(ni.leaf="<<Nodes[u.Ni].gtreepath.back()<<",nj.leaf="<<Nodes[u.Nj].gtreepath.back()<<")"<<endl;
        topk_r6.printResults(k);
        TIME_TICK_END
        TIME_TICK_PRINT(" TkGSKQ_NVD runtime: ")


        endTime = clock();
        double runTime = (double)(endTime - startTime)/CLOCKS_PER_SEC * 1000;
        total_Time += runTime;
        cout<<"topk-"<<method<<" runtime: "<<runTime<< "ms" << endl;

        delete  []poiADJMark; delete [] poiMark;
        cout<<"------------------------------------------------------------------------"<<endl;
        cout<<"*****************test_topk()测试结束, average runTime="<<total_Time/(100.0-65)<<"ms**************************"<<endl;
        getchar();

   // }




}


void test_GSKScoreAndTopk(int u_id, int p_id,int A, float alpha, int k){

    Load_Hybrid_NVDG_AddressIdx_fast();

    User user = getUserFromO2UOrgLeafData(u_id);
    printUsrInfo(user);
    POI poi = getPOIFromO2UOrgLeafData(p_id);
    printPOIInfo(poi);

    cout<<"usrToPOIDistance(user, poi)="<<usrToPOIDistance(user, poi)<<endl;
    cout<<"usrToPOIDistance_phl(user, poi)="<<usrToPOIDistance_phl(user, poi)<<endl;



    double gsk_score =  getGSKScore_o2u(DEFAULT_A, alpha, poi, user);
    double gsk_score_phl =  getGSKScore_o2u_phl(DEFAULT_A, alpha, poi, user);
    double score_topk_base = -1;

#ifdef LV
    TopkQueryCurrentResult topk_r4 = topkSDijkstra_verify_usr_disk(transformToQ(user),k,A,alpha,1000,0.0); //topkSDijkstra_verify_usr_disk
    score_topk_base = topk_r4.topkScore;
    cout<<"dj下结果:"; topk_r4.printResults(k);
#else
    priority_queue<Result> results3 = TkGSKQ_bottom2up_verify_disk(u_id,k,A,alpha,1000,0.0);
    score_topk_base= results3.top().score;
    cout<<"gimtree下结果:"; printTopkResultsQueue(results3);
#endif


    User u = getUserFromO2UOrgLeafData(u_id);
    TopkQueryCurrentResult topk_r6 = TkGSKQ_NVD(user, k, A,alpha);//,poiMark,poiADJMark);
    double score_topk_nvd = topk_r6.topkScore;
    cout<<"nvd下结果：：";topk_r6.printResults(k);

    cout<<"gsk_score="<<gsk_score<<",gsk_score_phl ="<<gsk_score_phl<<endl;
    cout<<"sk(baseline)="<<score_topk_base<<", sk(nvd)="<<score_topk_nvd<<endl;

}


void test_topk_naive(int k, float alpha){
    cout<<"进行test_topk()测试"<<endl;
    for(int i=64;i<100;i++){
        cout<<"run top-k for u"<<i<<endl;
        User u;
        double score_topk = -1;
#ifndef DiskAccess
        u = Users[i];
        TIME_TICK_START
        TopkQueryCurrentResult topk_r2 = topkSDijkstra_verify_usr_memory(transformToQ(u),k,DEFAULT_A,alpha,1000,0.0); //topkSDijkstra_verify_usr_disk
        score_topk = topk_r2.topkScore;
        cout<<"u"<<i<<": score_topk memeory(dj-现在)="<<score_topk<<"具体为："<<endl;
        topk_r2.printResults(k);
        TIME_TICK_END
        TIME_TICK_PRINT(" topkdjstra-memory 现在 runtime: ")


#else
        u = getUserFromO2UOrgLeafData(i);



        TIME_TICK_START
        TopkQueryCurrentResult topk_r5 = topkBruteForce_disk(u,k,DEFAULT_A,alpha,1000,0.0);
        score_topk = topk_r5.topkScore;
        cout<<"u"<<i<<": score_topk(bruteForce)="<<score_topk<<"具体为："<<endl;
        cout<<"(ni.leaf="<<Nodes[u.Ni].gtreepath.back()<<",nj.leaf="<<Nodes[u.Nj].gtreepath.back()<<")"<<endl;
        topk_r5.printResults(k);
        TIME_TICK_END
        TIME_TICK_PRINT(" topkBruteForce runtime: ")
        cout<<endl;
#endif



        cout<<endl;
        //if(i%2==0)
        //getchar();
    }
    cout<<"test_topk()测试结束"<<endl;
    //getchar();
}


void test_topk_gtree(int k, float alpha){
    cout<<"进行test_topk()测试"<<endl;
    for(int i=65;i<100;i++){
        //i = 65;
        cout<<"run top-k for u"<<i<<endl;
        User u;
        double score_topk = -1;

        u = getUserFromO2UOrgLeafData(i);


        TIME_TICK_START
        priority_queue<Result> results1 = TkGSKQ_bottom2up_verify_disk(i,k,DEFAULT_A,alpha,1000,0.0);
        score_topk = results1.top().score;
        cout<<"u"<<i<<": score_topk(g-tree-Ni&Nj)="<<score_topk<<"具体为："<<endl;
        cout<<"(ni.leaf="<<Nodes[u.Ni].gtreepath.back()<<",nj.leaf="<<Nodes[u.Nj].gtreepath.back()<<")"<<endl;
        printTopkResultsQueue(results1);
        //cout<<"upper_jumper="<<upper_jumper<<endl;
        upper_jumper =0;
        TIME_TICK_END
        TIME_TICK_PRINT(" runtime(gtree-Ni&Nj) runtime")
        cout<<endl;

        //if(i%2==0)
        //getchar();
    }
    cout<<"test_topk()测试结束"<<endl;
    //getchar();
}

void test_topk_dijstra(int k, float alpha){
    cout<<"进行test_topk()测试"<<endl;
    for(int i=65;i<100;i++){
        //i = 65;
        cout<<"run top-k for u"<<i<<endl;
        User u;
        double score_topk = -1;

#ifndef DiskAccess
        u = Users[i];
        TIME_TICK_START
        TopkQueryCurrentResult topk_r2 = topkSDijkstra_verify_usr_memory(transformToQ(u),k,DEFAULT_A,alpha,1000,0.0); //topkSDijkstra_verify_usr_disk
        score_topk = topk_r2.topkScore;
        cout<<"u"<<i<<": score_topk memeory(dj-双向)="<<score_topk<<"具体为："<<endl;
        topk_r2.printResults(k);
        TIME_TICK_END
        TIME_TICK_PRINT(" topkdjstra-memory 双向 runtime: ")
        cout<<endl;


#else

        u = getUserFromO2UOrgLeafData(i);

        TIME_TICK_START
        TopkQueryCurrentResult topk_r4 = topkSDijkstra_verify_usr_disk(transformToQ(u),k,DEFAULT_A,alpha,1000,0.0); //topkSDijkstra_verify_usr_disk
        score_topk = topk_r4.topkScore;
        cout<<"u"<<i<<": score_topk(dj-Ni&Nj)="<<score_topk<<"具体为："<<endl;
        cout<<"(ni.leaf="<<Nodes[u.Ni].gtreepath.back()<<",nj.leaf="<<Nodes[u.Nj].gtreepath.back()<<")"<<endl;
        topk_r4.printResults(k);
        TIME_TICK_END
        TIME_TICK_PRINT(" topkdjstra-双向 runtime: ")
        cout<<endl;

#endif

        //if(i%2==0)
        //getchar();
    }
    cout<<"test_topk()测试结束"<<endl;
    //getchar();
}

float SPD_by_topkSDijkstra_enhance(int s, int e){

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

    //priority_queue<nodeDistType> tmpPoint;
    priority_queue<Result> resultFinal;
    BinaryMinHeap<EdgeWeight,NodeID> pqueue;


    pointval[s] = 0;
    pqueue.insert(s,0);

    double p2p_dist = 999999999999999999999.9;
    while (pqueue.size()>0) {
        int entry_dist = pqueue.getMinKey();

        int tmpNi = pqueue.extractMinElement();


        if (vis[tmpNi] == true)
            continue;
        nodeCnt++;
        vis[tmpNi] = true;
        if(tmpNi==e) {//到达终点
            p2p_dist = entry_dist;
            //return  tmp.dist;
        }


        vector<int> tmpAdj = adjList[tmpNi];    //取得拓展点的邻接表（地址）
        vector<float> tmpAdjW = adjWList[tmpNi];
        for (size_t i = 0; i < tmpAdj.size(); i++) {  // expand each edge
            int tmpNj = tmpAdj[i];
            if (vis[tmpNj])
                continue;

            float edge_dist = tmpAdjW[i];   //取得邻接表中的edge（地址）
            //int Ni = e.Ni;

            //首次拓展当前拓展顶点
            if(rearch[tmpNj]==false){
                pointval[tmpNj] = maxDisValue;
            }
            if (edge_dist + entry_dist < pointval[tmpNj]) {
                pointval[tmpNj] = edge_dist + entry_dist;
                int newDistance = edge_dist + entry_dist;
                pqueue.insert(tmpNj,newDistance);

                rearch[tmpNj]=true;

            }
        }
        //路网边拓展提前终止条件

    }
    //cout<<"不可达！"<<endl;
    return p2p_dist;  //注意此时的score只是current score
}


void test_P2P(){

    // 加载 Gtree索引
    initGtree();
    // 加载 edge information
    /// loadEdgeMap(); jins记得补上

    //加载 PHL 索引
    PrunedHighwayLabeling phl;
    //stringstream phlFilePath;
    //phlFilePath << DATASET_HOME <<"road/"<<dataset_name<<".phl";
    phl.LoadLabel(getRoadInputPath(PHL).c_str());


    clock_t start_time, end_time;
    double runtime;
    double accumulation1 = 0; double accumulation2 = 0; double accumulation3 = 0;
    int count = 0;
    /*int dist1  = SPSP(7220,2412);
    int dist2  = phl.Query(7220,2412);
    cout<<"dist1="<<dist1<<",dist2="<<dist2<<endl;
    dist1  = SPSP(7220,9390);
    dist2  = phl.Query(7220,9390);
    cout<<"dist1="<<dist1<<",dist2="<<dist2<<endl;

    exit(1);*/

    while(true){
        int u  = rand()% 100000;  //65; //
        int v = rand()% 100000;   //7813; //VertexNum

        u = 909415; v = 1958394;
        /*int u = rand()%VertexNum;
        int v = rand()%VertexNum;*/
        //TIME_TICK_START
        start_time = clock();
        int dist_Gtree = SPSP(u, v);
        end_time = clock();
        double dur_time =  (double)(end_time - start_time);
        runtime =dur_time /CLOCKS_PER_SEC * 1000000.0;
        cout<<"SPDist(n"<<u<<",n"<<v<<") by Gtree ="<<dist_Gtree<<", runtime="<<runtime<<"us"<<endl;
        accumulation1+= runtime;
        //TIME_TICK_END
        //TIME_TICK_PRINT_us(" runtime(Gtree_SPSP):")
        //TIME_TICK_START
        start_time = clock();
        int dist_phl = phl.Query(u, v);
        end_time = clock();
        double dur_time2 =  (double)(end_time - start_time);
        runtime = dur_time2 /CLOCKS_PER_SEC * 1000000;
        cout<<"SPDist(n"<<u<<",n"<<v<<") by PHL ="<<dist_phl<<", runtime="<<runtime<<"us"<<endl;
        accumulation2+= runtime;
        // TIME_TICK_END
        //TIME_TICK_PRINT_us(" runtime(PHL):")
        count++;

        if(false){
            start_time = clock();
            float dist_dj = SPD_by_topkSDijkstra_memory(u, v);
            end_time = clock();
            double dur_time3 =  (double)(end_time - start_time);
            runtime = dur_time3 /CLOCKS_PER_SEC * 1000;
            cout<<"SPDist(n"<<u<<",n"<<v<<") by -dj ="<<dist_dj<<", runtime="<<runtime<<"ms"<<endl;
            //cout<<"u.Ni="<<u.Ni<<",p.Nj="<<p.Nj<<endl;
            accumulation3+= runtime;

            start_time = clock();
            dist_dj = SPD_by_topkSDijkstra_enhance(u, v);
            end_time = clock();
            double dur_time4 =  (double)(end_time - start_time);
            runtime = dur_time4 /CLOCKS_PER_SEC * 1000;
            cout<<"SPDist_enhance(n"<<u<<",n"<<v<<") by -dj ="<<dist_dj<<", runtime="<<runtime<<"ms"<<endl;
            //cout<<"u.Ni="<<u.Ni<<",p.Nj="<<p.Nj<<endl;
            accumulation3+= runtime;
        }





        count++;
        //getchar();

       /* if(dist_phl!=dist_Gtree){
            cout<<"SPDist("<<u.Ni<<","<<p.Nj<<"),dist_Gtree="<<dist_Gtree<<",dist_phl="<<dist_phl<<endl;
            getchar();
        }*/
        if(true){ //count%5==0
            cout<<"---------------------average Gtree runtime="<<accumulation1/count<<"us, average PHL runtime="<<accumulation2/count<<"us,average dijsra runtime="<<accumulation3/count<<"ms--------------------"<<endl;
            getchar();
        }



        /*TIME_TICK_START
        int leaf1 = Nodes[u.Ni].gtreepath.back();
        int leaf2 = Nodes[p.Ni].gtreepath.back();

        float dist2 ;
        if(leaf1==leaf2){
            //dist2 = SPSP(u.Ni, p.Nj)/1000.0;
            cout<<"算过了！"<<endl;
        }
        else{
            float dist_bordersMax = MaxDisBound[leaf1][leaf2];
            float dist_bordersMin = MinDisBound[leaf1][leaf2];
            cout<<"SPDist_by_Matrix("<<u.id<<","<<p.id<<")="<<dist_bordersMin<<endl;
        }

        TIME_TICK_END
        TIME_TICK_PRINT_us(" runtime(new):")*/

        //getchar();
        /*if(dist<500)
             getchar();*/
    }

}


void test_u2o(){
    stringstream phlFilePath;
    phlFilePath << DATASET_HOME <<"road/"<<dataset_name<<".phl";
    string phl_fileName = phlFilePath.str();

    PrunedHighwayLabeling phl;
    phl.LoadLabel(phl_fileName.c_str());
    clock_t start_time, end_time;
    double runtime;
    double accumulation1 = 0; double accumulation2 = 0; double accumulation3 = 0;
    int count = 0;
    while(true){
        int u_id = rand()% UserID_MaxKey;  //65; //
        int p_id = rand()% poi_num;   //7813; //
        User u = getUserFromO2UOrgLeafData(u_id);
        POI  p = getPOIFromO2UOrgLeafData(p_id);
        /*int u = rand()%VertexNum;
        int v = rand()%VertexNum;*/
        //TIME_TICK_START
        start_time = clock();
        float dist_Gtree = SPSP(u.Ni, p.Nj)/1000.0;
        end_time = clock();
        double dur_time =  (double)(end_time - start_time);
        runtime =dur_time /CLOCKS_PER_SEC * 1000000.0;
        cout<<"SPDist(u"<<u.id<<",p"<<p.id<<") by Gtree_SPSP ="<<dist_Gtree<<", runtime="<<runtime<<"us"<<endl;
        accumulation1+= runtime;
        //TIME_TICK_END
        //TIME_TICK_PRINT_us(" runtime(Gtree_SPSP):")

        //TIME_TICK_START
        start_time = clock();
        float dist_phl = phl.Query(u.Ni, p.Nj)/1000.0;
        end_time = clock();
        double dur_time2 =  (double)(end_time - start_time);
        runtime = dur_time2 /CLOCKS_PER_SEC * 1000000;
        cout<<"SPDist(u"<<u.id<<",p"<<p.id<<") by PHL ="<<dist_phl<<", runtime="<<runtime<<"us"<<endl;
        accumulation2+= runtime;
        // TIME_TICK_END
        //TIME_TICK_PRINT_us(" runtime(PHL):")
        count++;

        start_time = clock();
        float dist_dj = SPD_by_topkSDijkstra_memory(u.Ni, p.Nj);
        end_time = clock();
        double dur_time3 =  (double)(end_time - start_time);
        runtime = dur_time3 /CLOCKS_PER_SEC * 1000000;
        cout<<"SPDist(u"<<u.id<<",p"<<p.id<<") by dj ="<<dist_dj<<", runtime="<<runtime<<"us"<<endl;
        //cout<<"u.Ni="<<u.Ni<<",p.Nj="<<p.Nj<<endl;
        accumulation3+= runtime;
        // TIME_TICK_END
        //TIME_TICK_PRINT_us(" runtime(PHL):")
        count++;
        //getchar();

        /* if(dist_phl!=dist_Gtree){
             cout<<"SPDist("<<u.Ni<<","<<p.Nj<<"),dist_Gtree="<<dist_Gtree<<",dist_phl="<<dist_phl<<endl;
             getchar();
         }*/
        if(count%5==0){
            cout<<"---------------------average Gtree runtime="<<accumulation1/count<<"us, average PHL runtime="<<accumulation2/count<<"us,average dijsra runtime="<<accumulation3/count<<"--------------------"<<endl;
            getchar();
        }



        /*TIME_TICK_START
        int leaf1 = Nodes[u.Ni].gtreepath.back();
        int leaf2 = Nodes[p.Ni].gtreepath.back();

        float dist2 ;
        if(leaf1==leaf2){
            //dist2 = SPSP(u.Ni, p.Nj)/1000.0;
            cout<<"算过了！"<<endl;
        }
        else{
            float dist_bordersMax = MaxDisBound[leaf1][leaf2];
            float dist_bordersMin = MinDisBound[leaf1][leaf2];
            cout<<"SPDist_by_Matrix("<<u.id<<","<<p.id<<")="<<dist_bordersMin<<endl;
        }

        TIME_TICK_END
        TIME_TICK_PRINT_us(" runtime(new):")*/

        //getchar();
        /*if(dist<500)
             getchar();*/
    }

}

void test_CheckSocialConnection(){
    //读取社交网络图信息
    loadFriendShipData();

    //读取用户在兴趣点的check-in
    loadCheckinData();

    //加载user信息
    int User_id, User_Ni, User_Nj;
    float User_dist, User_dis;


    int usrKeyCnt = 0;

    //空间初始化
    Users.clear();
    for(int i=0;i<UserID_MaxKey;i++){
        User u;
        Users.push_back(u);
    }

    ifstream finUsr;
    //finUsr.open(FILE_USER);
    string userPath;

    userPath = getObjectNewInputPath(USER);

    finUsr.open(userPath.c_str());

    //计算User信息
    string str;
    while (getline(finUsr, str)) {
        istringstream tt(str);
        tt >> User_id >> User_Ni >> User_Nj >> User_dist >> User_dis;
        float exchangeNum = (int) (User_dist * WEIGHT_INFLATE_FACTOR);
        User_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (User_dis * WEIGHT_INFLATE_FACTOR);
        User_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;


        //更新节点的关键词属性
        vector<int> uKey;
        int keyTmp;

        while (tt >> keyTmp) {
            if (-1 == keyTmp) {
                break;
            }
            usrKeyCnt++;
            uKey.push_back(keyTmp);

        }

        User tmpNode;
        tmpNode.id = User_id;
        tmpNode.Ni = min(User_Ni,User_Nj);
        tmpNode.Nj = max(User_Ni,User_Nj);
        //tmpNode.category = 0;
        tmpNode.dis = User_dis;
        tmpNode.dist = User_dist;
        tmpNode.keywords = uKey;
        tmpNode.isVisited = false;
        Users[tmpNode.id] =tmpNode;
        usrCnt++;


    }

    finUsr.close();
    cout << " DONE! USER #:" << usrCnt;
    cout << " , USER KEY #: " << usrKeyCnt << endl;
    LOAD_END
    LOAD_PRINT(" load User ")

    int count = 0; int nonISO_num=0;
    for(int i=0;i<Users.size();i++){
        if(Users[i].keywords.size()>0){
            int u_id = Users[i].id;
            //User u = getUserFromO2UOrgLeafData(u_id);
            if(friendshipMap.count(u_id)){
                cout<<"用户u"<<u_id<<"的朋友为:";
                printElements(friendshipMap[u_id]);
                nonISO_num++;
            }
            else{
                cout<<"用户"<<u_id<<"为孤立点"<<endl;
            }
            /*if(count%20==0){
                getchar();
            }*/
            count++;
            if(count>reviewed_userNum) break;

        }
        else continue;
    }
    cout<<reviewed_userNum<<"个有评论的用户中，有社交连接的用户个数："<<nonISO_num<<endl;
}


void test_u2PDist(){ //int poi_id, int usr_id
    while(true){
        int usr_id = rand()% UserID_MaxKey;  //65; //
        int poi_id = rand()% poi_num;   //7813; //
        cout<<"test Dist(u"<<usr_id<<",p)"<<poi_id<<endl;
        POI p = getPOIFromO2UOrgLeafData(poi_id); User u = getUserFromO2UOrgLeafData(usr_id);
        cout<<"by gtree:"<<usrToPOIDistance(u,p);
        cout<<",by phl:"<<usrToPOIDistance_phl(u,p)<<endl;
    }

}

void test_u2PDist_single(int poi_id, int usr_id){ //
    //while(true){

        cout<<"test Dist(u"<<usr_id<<",p"<<poi_id<<endl;
        POI p = getPOIFromO2UOrgLeafData(poi_id); User u = getUserFromO2UOrgLeafData(usr_id);
        cout<<usrToPOIDistance(u,p)<<endl;
        cout<<usrToPOIDistance_phl(u,p)<<endl;
    //}

}

void test_enumerating_invertedList(){
    int poi_id = 0; vector<int> adj_poi;
    int term_size = getTermSize();
    for(int i= 1;i<=term_size;i++){
         int term = i;
         vector<int> posting_list = getTermOjectInvList(term);
         for(int poi_id: posting_list){
                    vector<int> adj_poi = getPOIAdj_NVD_By_Keyword(poi_id,term);
                    cout<<"p"<<poi_id<<"在NVD中的neighbor为：";
                    printElements(adj_poi);
                    //getchar();
          }
          cout<<"t"<<term<<" 's posting list size="<<posting_list.size()<<endl;
        getchar();
    }

}


//NVD测试

void test_v2p_NVDAdj(){
    int term_size = getTermSize();
    while(true){
        int vertex = random() % VertexNum;
        int term = random()% term_size+1;
        term = min(term,term_size-1);
        int _size = getTermOjectInvListSize(term);
        if(_size<=20) continue;
        int poi_id = getNNPOI_By_Vertex_Keyword(vertex,term);
        cout<<"v"<<vertex<<"在路网中，在关键词"<<term<<"下的最近邻兴趣点，为："<<endl;
        cout<<"o"<<poi_id<<endl;
        getchar();

        if(poi_id<0){
            cout<<"兴趣点id: "<<poi_id<<", 内容错误！"<<endl;
            exit(-1);
        }
        poi_id = 9619;term=1;
        cout<<"o"<<poi_id<<"在关键词"<<term<<"下在NVD中的neighbor为："<<endl;
        vector<int> adj_poi = getPOIAdj_NVD_By_Keyword(poi_id,term);
        printElements(adj_poi);


        poi_id = 9619;term=20;
        cout<<"o"<<poi_id<<"在关键词"<<term<<"下在NVD中的neighbor为："<<endl;
        adj_poi = getPOIAdj_NVD_By_Keyword(poi_id,term);
        printElements(adj_poi);
        getchar();
    }
}


void test_v2p_Then_NVDAdj(int vertex, int term){
    //int term_size = getTermSize();
    //while(true){

    cout<<"v"<<vertex<<"在路网中，在关键词"<<term<<"下的最近邻兴趣点，为："<<endl;
    int poi_id_first = getNNPOI_By_Vertex_Keyword(vertex,term);
    cout<<"o"<<poi_id_first<<endl;
    //getchar();

}


void test_poiNN_by_Hash(int vertex, int term){
    //int term_size = getTermSize();
    //while(true){
    cout<<"*****在vertex only hash索引中检索:"<<endl;
    cout<<"v"<<vertex<<"在路网中，在关键词"<<term<<"下的最近邻兴趣点，为："<<endl;
    int poi_id = getNNPOI_By_Vertex_Keyword(vertex,term);
    cout<<"最近邻poi： o"<<poi_id<<endl;
    //getchar();
    cout<<"******在hybrid hash索引中检索:"<<endl;
    vector<int> pois;
    int poi_id2 = getNNPOI_By_HybridVertex_Keyword(vertex,term);
    if(poi_id2==-1){
        cout<<"最近邻poi找不到！";
        int leaf_node = Nodes[vertex].gtreepath.back();
        pois = getNNPOIList_by_Hybrid_Hash(leaf_node, term);
        cout<<",转而找最近邻poi group list, 具体为：";
        printElements(pois);
    }

    else
        cout<<"最近邻poi： o"<<poi_id2<<endl;


}


void test_poiNN_by_Hash_all_round(){
    while(true){
        int vertex = (int)random()% VertexNum;
        int term = (int)random()%vocabularySize;
        term = max(1,term);term = 1;
        int size = getTermOjectInvListSize(term);
        if(size<=posting_size_threshold) continue;
        cout<<"-------------for vertex"<<vertex<<" ,term"<<term<<",size="<<size<<"---------------"<<endl;
        test_poiNN_by_Hash(vertex,term);
        getchar();
    }

}


void test_u2p_Then_NVDAdj(int u_id, int term){
    //int term_size = getTermSize();
    //while(true){
        User user = getUserFromO2UOrgLeafData(u_id);
        int vertex = user.Ni;
        int _size = getTermOjectInvListSize(term);
        if(_size<=posting_size_threshold) {
            cout<<"低频词汇"<<endl;
            return;
        }
        cout<<"v"<<vertex<<"在路网中，在关键词"<<term<<"下的最近邻兴趣点，为："<<endl;
        int poi_id_first = getNNPOI_By_Vertex_Keyword(vertex,term);
        cout<<"o"<<poi_id_first<<endl;
        //getchar();

        if(poi_id_first<0){
            cout<<"兴趣点id: "<<poi_id_first<<", 内容错误！"<<endl;
            exit(-1);
        }
        //poi_id = 9619;term=1;
        bool poiADJMask[poi_num];
        memset(poiADJMask,false, sizeof(poiADJMask));
        cout<<"遍历起始poi, o"<<poi_id_first<<"when 关键词=t"<<term<<endl;
        poiADJMask[poi_id_first] = true;
        vector<int> adj_poi = getPOIAdj_NVD_By_Keyword(poi_id_first,term);
        //printElements(adj_poi);
        priority_queue<int> access_list;
        for(int poi: adj_poi){
            if(poiADJMask[poi] == true) continue;
            access_list.push(poi);
            poiADJMask[poi] = true;
        }

        while(access_list.size()>0){
            int poi_list_first = access_list.top();
            access_list.pop();
            cout<<"遍历o"<<poi_list_first<<endl;
            vector<int> adj_poi2 = getPOIAdj_NVD_By_Keyword(poi_list_first,term);
            for(int poi: adj_poi2){
                if(poiADJMask[poi] == true)
                    continue;
                access_list.push(poi);
                poiADJMask[poi] = true;
            }
        }
        //从poi_id_first开始在term对应的NVD中遍历


    //}
}

void test_NVDAdj_GivenKeyword(int term){


    //获得term对应的posting list中的poi

    vector<int> poi_related = getTermOjectInvList(term);
    cout<<"关键词t"<<term<<"的posting list size="<<poi_related.size()<<",其中的poi为：";
    printElements(poi_related);

    for(int i=0;i<poi_related.size();i++){
        int poi_id = poi_related[i];
        vector<int> adj_poi = getPOIAdj_NVD_By_Keyword(poi_id,term);
        if(poi_id>=936111){
            cout<<i<<": p"<<poi_id<<"在NVD Graph("<<term<<"中的邻居为：" ;
            printElements(adj_poi);
            getchar();
        }
        //cout<<endl;
    }



    //}
}




void test_NVD(){
    int term = 104;

    if(false){
        TIME_TICK_START
        NVD_generation_Keyword_aware(term);
        TIME_TICK_END
        TIME_TICK_PRINT_us(" NVD_generation_Keyword_aware")
    }
    else{
        test_v2p_NVDAdj();
    }


}


int test_Btree(int endNum)
{
    int start = 0;
    int end = endNum;
    char* file_name = "../test.bt"; //../test.db

    clock_t startTime,endTime;
    bool btree_new = false;
    bpt::bplus_tree database(file_name, btree_new);
    if(btree_new){
        for (int i = start; i <= end; i++) {
            if (i % 1000 == 0)
                printf("%d\n", i);
            char key_char[16] = { 0 };
            sprintf(key_char, "%d", i);
            startTime = clock();//计时开始
            database.insert(key_char, i);
            printf("插入内容：key= %s, value= %d \n", key_char, i);
            endTime = clock();//计时结束
            cout << "insert runtime is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC *1000 << "ms" << endl;
        }
        //printf("%d\n", end);
        printf("btree完成插入%d条记录\n",(end-start));
    }
    else{
        cout<<"加载btree success!"<<endl;
        for(int i=0;i<200;i++){
            char key_char[16] = { 0 };
            sprintf(key_char, "%d", i);
            startTime = clock();//计时开始
            int value = -1;
            if (database.search(key_char, &value) != 0)
                printf("----------------------------------Key %d not found\n-----------------------------------", i);
            else
                printf("find key=%s, value=%d\n", key_char,value);
            endTime = clock();//计时结束
            cout << "Search run time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC *1000 << "ms" << endl;
        }
    }


    return 0;
}




void test_offline_compuation(int node){
    int v_idx1 = lower_bound(GTree[17].leafnodes.begin(), GTree[17].leafnodes.end(), 23550) -
                 GTree[17].leafnodes.begin();
    int v_idx2 = id2idx_vertexInLeaf[23550];

    int b_th =1; int border1 = GTree[17].borders[b_th];
    int distance_offline = GTree[17].mind[b_th * GTree[17].leafnodes.size() + v_idx1];
    vector<int> vcands; vcands.push_back(23550);
    vector<int> results = dijkstra_candidate(border1, vcands, Nodes);
    int tmp = results[0];
    double distance1 =  tmp /1000.0;
    double distance2 =  distance_offline /1000.0;
    cout<<"distance1="<<distance1<<", distance2="<<distance2<<endl;

}



void test_edgeInfo(int tmpNi) {
    int AdjGrpAddr, AdjListSize;
    vector<int> tmpAdj = adjList[tmpNi];
    /*1. 取得地址，   2. 根据地址从文件中取得数据*/
    //---------------取得顶点tmpNi的数据地址：AdjGrpAddr
    AdjGrpAddr = getAdjListGrpAddr(tmpNi);
    //cout<<"AdjGrpAddr="<<AdjGrpAddr<<endl;
    bool testanswer = true;
    //---------------取得顶点tmpNi的邻接点个数：AdjListSize
    getFixedF(SIZE_A, Ref(AdjListSize), AdjGrpAddr);
    //getFixedF_Self(SIZE_A,Ref(AdjListSize),AdjGrpAddr);
    //cout<<"开始N"<<tmpNi<<"的拓展，邻接点size= "<<AdjListSize<<endl;
    for (size_t j = 0; j < AdjListSize; j++) {  // expand each edge
        //int tmpNj = tmpAdj[i];
        //------------取得邻接顶点tmpNj的id：tmpNj
        int tmpNj;
        getVarE(ADJNODE_A, Ref(tmpNj), AdjGrpAddr, j);  //读到的Nj不对
        //cout<<"for Ni"<<tmpNi<<"读到Nj"<<tmpNj<<endl;

        edge e = EdgeMap[getKey(tmpNi, tmpNj)];
        int Ni = e.Ni;




        //cout<<"Ni"<<Ni<<endl;
        //float edgeDist = e.dist;
        float edgeDist;
        getVarE(DIST_A, Ref(edgeDist), AdjGrpAddr, j);  //读取边的长度(成功)
        //cout<<"edgeDist="<<edgeDist<<endl;

        //读取Okey psudo doc 的地址,并赋予okeyAddr


        set<int> OKeyset;
        //j=1;
        OKeyset = getOKeyEdge(tmpNi, j);

        if (tmpNi == 96697) {
            cout << "find tmpNi==96697,此时 Nj=" << tmpNj << endl;
            OKeyset = getOKeyEdge(tmpNi, j);
            cout << "Okeyset:";
            printSetElements(OKeyset);
            //getchar();
        }

    }
    cout<<"test is over!"<<endl;
}


void test_groupTopk(int k, float alpha){
    cout<<"进行test_topk()测试"<<endl;
    TreeNode tn; int leaf;
    for(int i=1;i<GTree.size();i++){
        tn = GTree[i]; leaf = i;
        if(tn.isleaf)
            break;
    }
    cout<<"抽取叶节点"<<leaf<<"中所有用户..."<<endl;
    set<int> users_withinLeaf;
    tn = getGIMTreeNodeData(leaf,OnlyU);
    set<int> UnionUkeywords = tn.userUKeySet;
    for(int term: UnionUkeywords){
        vector<int> related_users = getUsrTermRelatedEntry(term,leaf);
        for(int u: related_users)
            users_withinLeaf.insert(u);
    }
    int _size = users_withinLeaf.size();
    cout<<"叶节点"<<leaf<<"中共有用户"<<_size<<"个"<<endl;
    vector<query> queries; float score_bounds[_size];
    vector<int> u_list;
    int i = 0;
    for(int u: users_withinLeaf){
        if(u%2==0) continue;
        query Query = transformToQ(getUserFromO2UOrgLeafData(u));
        queries.push_back(Query);
        u_list.push_back(u);
        score_bounds[i] = 100;
    }
    cout<<"提取叶节点"<<leaf<<"中"<<queries.size()<<"个用户"<<endl;
    /*cout<<"-----------------------------------------------topk one by one -----------------------------------------------"<<endl;
    TIME_TICK_START
    for(int u: u_list){
        priority_queue<Result> single_results = TkGSKQ_bottom2up_verify_disk(u,k,DEFAULT_A,alpha,1000,0.0);
        double score_topk = single_results.top().score;
        cout<<"u"<<u<<": score_topk="<<score_topk<<"具体为："<<endl;
        printTopkResultsQueue(single_results);
    }
    TIME_TICK_END
    TIME_TICK_PRINT("topk one by one")*/

    cout<<"-------------------------------------------------group topk ---------------------------------------------------"<<endl;
    //group_topkSDijkstra_verify_usr_disk(leaf, queries, k,DEFAULT_A, alpha, score_bounds);
    TIME_TICK_START
    vector<priority_queue<Result>> group_results = Group_TkGSKQ_bottom2up_verify_disk(leaf,u_list,k,DEFAULT_A, alpha, 100, 0);
    print_groupTopkResults(group_results,u_list,leaf);
    TIME_TICK_END
    TIME_TICK_PRINT("group topk")


}


/*
void test_single_topk(int id, int k, float alpha){

        int u_id = id;
        cout<<"run top-k for u"<<u_id<<endl;
        TIME_TICK_START
        User u;
        double score_topk = -1;
#ifndef DiskAccess
        u = Users[u_id];
    #ifdef LV
        TopkQueryCurrentResult results = topkSDijkstra_verify_usr_memory(transformToQ(u),k,DEFAULT_A,alpha,1000,0.0).topkScore;
        score_topk = results.topkScore;
    #else
        priority_queue<Result> results = TkGSKQ_bottom2up_verify_memory(id,k,DEFAULT_A,alpha,1000,0.0);
        score_topk = results.top().score;
    #endif
#else
        u = getUserFromO2UOrgLeafData(u_id);
        printUsrInfo(u);

#ifdef LV
        TopkQueryCurrentResult results = topkSDijkstra_verify_usr_disk(transformToQ(u),k,DEFAULT_A,alpha,1000,0.0).topkScore;
        score_topk = results.topkScore;
#else
        priority_queue<Result> results = TkGSKQ_bottom2up_verify_disk(u_id,k,DEFAULT_A,alpha,1000,0.0);
        score_topk = results.top().score;
#endif
        //cout<<"u"<<i<<": score_k2="<<results.top().score;
#endif

        cout<<"u"<<u_id<<": score_topk="<<score_topk;
#ifdef LV
        printTopkResultsVector(results);
#else
        printTopkResultsQueue(results);
#endif
        TIME_TICK_END
        TIME_TICK_PRINT(" runtime: ")

        cout<<endl;


    cout<<"test_topk()测试结束"<<endl;
    //getchar();
}
*/

/*
void test_allReverseTopk(int k, float alpha){
    cout<<"进行test_topk()测试"<<endl; map<int,set<int>> all_reverse;
    for(int i=1;i<UserID_MaxKey;i++){
        cout<<"run top-k for u"<<i<<endl;
        TIME_TICK_START
        User u; int u_id = i;
        double score_topk = -1;
#ifndef DiskAccess
        u = Users[i];
    #ifdef LV
        score_topk = topkSDijkstra_verify_usr_memory(transformToQ(u),k,DEFAULT_A,alpha,1000,0.0).topkScore;
    #else
        priority_queue<Result> results = TkGSKQ_bottom2up_verify_memory(i,k,DEFAULT_A,alpha,1000,0.0);
        score_topk = results.top().score;
    #endif
#else
        u = getUserFromO2UOrgLeafData(i);
        printUsrInfo(u);

#ifdef LV
        TopkQueryCurrentResult results =  topkSDijkstra_verify_usr_disk(transformToQ(u),k,DEFAULT_A,alpha,1000,0.0).topkScore;
        score_topk = results.topkScore;
#else
        priority_queue<Result> results = TkGSKQ_bottom2up_verify_disk(i,k,DEFAULT_A,alpha,1000,0.0);
        //score_topk = results.top().score;
#endif
        //cout<<"u"<<i<<": score_k2="<<results.top().score;
#endif

        //cout<<"u"<<i<<": score_topk="<<score_topk;
        TIME_TICK_END
        TIME_TICK_PRINT(" runtime: ")
#ifdef LV
        for(Result rr: results.topkElements){
            int o_id = rr.id;
            all_reverse[o_id].insert(u_id);
        }

#else
        while(results.size()>0){
            Result ele = results.top();
            results.pop();
            int o_id = ele.id;
            all_reverse[o_id].insert(u_id);
        }

#endif

    }
    cout<<"test_allReverseTopk 测试结束, 结果为"<<endl;
    int current_max = -1; int current_max_id=-1;
    map<int,set<int>>::iterator it; it = all_reverse.begin();
    while(it!= all_reverse.end()){
        int o_id = it->first;
        set<int> reverse_u = it->second;
        int total_inf =0;
        for(int u: reverse_u){
            for(int v : friendshipMap[u])
            total_inf += 1.0/followedMap[v].size();
        }
        if(total_inf>current_max){
            current_max = total_inf;
            current_max_id = o_id;
        }
        //cout<<"o"<<o_id<<" has"<<reverse_u.size()<<" potential users:"<<endl;
        //printSetElements(reverse_u);
        it++;
    }
    cout<<"o"<<current_max<<" has"<<all_reverse[current_max].size()<<" potential users:"<<endl;
    printSetElements(all_reverse[current_max]);


    //getchar();
}
*/


void test_makeAdjListFiles(const char* outfileprefix) {

    //printf("将邻接表信息转化为二进制文件！\n");
    char idxFileName[200];
    sprintf(idxFileName,"%s.adj",outfileprefix);	remove(idxFileName); // remove existing file
    FILE * alFile=fopen(idxFileName,"w+");
    int current_point = 0; int addr =0;

    int key=0,size,PtSize;
    //写入图顶点的个数
    int num = VertexNum;
    fwrite(&num,1,sizeof(int),alFile);
    current_point += sizeof(int);
    addr += sizeof(int);
    addr= addr +sizeof(int)*VertexNum; //跳过第一个路网顶点数目，与随后的VertexNum行顶点数据指针

    //写入图中各个顶点数据块的基地址（注意上行语句含义）
    for (int Ni=0;Ni < VertexNum;Ni++) {
        fwrite(&addr,1,sizeof(int),alFile);  //把Ni数据块基地址写入
        current_point += sizeof(int);
        //统计edge size
        int _size = adjList[Ni].size();
        addr += sizeof(int);
        for (int j=0;j< adjList[Ni].size();j++) {
            int Nk = adjList[Ni][j];
            if(Ni==96697&& Nk==96712){
                cout<<"统计邻接表地址信息：find 96697&& Nk==96712"<<endl;
                cout<<"e("<<Ni<<","<<Nk<<")"<<endl;
                cout<<"addr="<<addr<<endl;
                getchar();
            }

            //统计edge data space
            int addr_delta = (4*sizeof(int)+sizeof(float));  // 4个int分别为：Nj, POIDataAddr, POIIdxAddr, OKeyAddr ， float:dist
            addr += addr_delta;  //原来是*(3*sizeof(int)+sizeof(float))  //下一个路网顶点的地址更新（1.当前顶点的邻接点数， size个邻边内容（Nj_ID, Dist, poigrp_address, poiIdxList_address））

        }

    }

    cout<<"Okey信息基地址：addr="<<addr<<endl;

    //统计 e.OUnionKeys 的地址内容
    for (int Ni=0;Ni < VertexNum;Ni++) {
        for (int j=0;j< adjList[Ni].size();j++) {
            int Nj = adjList[Ni][j];
            edge& e = EdgeMap[getKey(Ni,Nj)];
            if(Ni==96697&& Nj==96712){
                cout<<"统计邻接边关键词内容信息：find 96697&& Nk==96712"<<endl;
                cout<<"e.okeyAddress="<<addr<<endl;
                cout<<"_okeyword_onEdge_size="<<e.OUnionKeys.size()<<endl;
                //getchar();
            }
            if(Ni>Nj) continue;  //反方向重复边，跳过！

            int okeyword_onEdge_size = e.OUnionKeys.size();
            addr += sizeof(int);
            for(int term: e.OUnionKeys){
                int _term = term;
                addr += sizeof(int);
            }


        }
    }

    //写完文件头后，逐个写入顶点邻接信息
    float distsum=0;
    for (int Ni=0;Ni<VertexNum;Ni++) {
        size=adjList[Ni].size();
        fwrite(&(size),1,sizeof(int),alFile);               //----首先写入邻接点个数
        current_point += sizeof(int);
        for (int j=0;j < adjList[Ni].size();j++) {
            int Nj = adjList[Ni][j];	// Nk can be smaller or greater than Ni !!!
            edge e = EdgeMap[getKey(Ni,Nj)];

            if(Ni==96697&& Nj==96712){
                cout<<"write 邻接信息：Ni==96697&& Nk==96712"<<endl;
                cout<<"write e.okeyAddress="<<e.okeyAddress<<endl;
                cout<<"current_point="<<current_point<<endl;
            }

            PtSize=e.pts.size();

            fwrite(&Nj,1,sizeof(int),alFile);               //----写入邻接点ID
            fwrite(&(e.dist),1,sizeof(float),alFile);      //----写入该邻边的长度

            fwrite(&e.poiAddress,1,sizeof(int),alFile);                //----写入邻接边对应的POI data group的地址
            //cout<<"e.poiAddress="<<e.poiAddress;
            fwrite(&e.poiIdxAddress,1,sizeof(int),alFile);   //----写入邻接边上POI 索引表的地址


            int _address = e.okeyAddress;
            fwrite(&e.okeyAddress,1,sizeof(int),alFile);  //----写入邻接边 text psudo doc的地址
            current_point = current_point + (4*sizeof(int)+sizeof(float));


        }
    }
    cout<<"写完邻接信息后current address="<<current_point<<endl;
    //getchar();
    //写完顶点邻接表信息后，逐个写入邻边上的关键词信息
    for (int Ni=0;Ni<VertexNum;Ni++) {
        for (int j=0;j < adjList[Ni].size();j++) {

            int Nj = adjList[Ni][j];	// Nk can be smaller or greater than Ni !!!
            // 考虑到无向图上 (Ni, Nk) 与  (Nk, Ni)是同一条边，此时注意 j>i的情况（同一条边）
            if(Nj < Ni) continue; //careful（和上个careful对应）


            edge& e = EdgeMap[getKey(Ni,Nj)];

            //写入边上所有poi keyword的 vocabulary信息
            set<int> okeyword_onEdge = e.OUnionKeys;
            int _size = okeyword_onEdge.size();
            //写入边上的oKeySize_Edge
            if(Ni==96697&& Nj==96712){
                cout<<"write 边上关键词信息：Ni==96697 && Nk==96712"<<endl;
                cout<<"_size="<<_size<<endl;
                cout<<"边上关键词："; printSetElements(okeyword_onEdge);
                cout<<"current_point="<<current_point<<endl;
                getchar();
            }

            fwrite(&_size, 1, sizeof(int), alFile);
            current_point += sizeof(int);

            //写入edge 上的各keywords
            for(int term: okeyword_onEdge){
                fwrite(&term, 1, sizeof(int), alFile);
                current_point += sizeof(int);
            }

        }
    }


    distsum=distsum/2;
    fclose(alFile);
    printf("making AdjList success!\n");
    //getchar();
    //fout.close();
}






void test_lclLeaf(){
    vector<int> Ukeys;
    Ukeys.push_back(1);  Ukeys.push_back(3);  Ukeys.push_back(22);
    MultiLCLResult lclResult_multi = getLeafNode_LCL_Dijkstra_multi_disk(328, Ukeys, 20,DEFAULT_A, DEFAULT_ALPHA,0.005141);
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
    cout<<", update_o_set size="<<update_o_set.size()<<endl;
}



void test_getLeafNode_LCL_bottom2up_multi(){
    int node = 1765;
    vector<int> Ukeys;//1,76,
    Ukeys.push_back(1);
    Ukeys.push_back(76);
    //Ukeys.push_back(203);
    double score_bound = 0.00583063;
    MultiLCLResult lclResult_multi = getLeafNode_LCL_bottom2up_multi(node, Ukeys, DEFAULT_K, DEFAULT_A, DEFAULT_ALPHA, score_bound);
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
    cout<<"update_o_set size="<<update_o_set.size()<<endl;
    //printLCLResultsInfo()
}




//叶节点路网距离上下界计算测试函数
void test_DisBoundWithinLeaf(){

    //加载G-Tree索引
    ConfigType setting;
    AddConfigFromFile(setting,"../map_exp");
    //AddConfigFromCmdLine(setting,argc,argv);
    initGtree();

    // 加载 edge information
    loadEdgeMap();
    MinMaxD mm;


    cout<<"test_outputNodeDisBound map:"<<dataset_name<<endl;

    vector<int> leafNodes;
    vector<int> UpperNodes;

    map<int,map<int,MinMaxD>> DisBoundMap;
    map<int,map<int,bool>> DisBOOLMap;
    //map<int,map<int,bool >> DisBoundMap
    //取出叶子节点
    for(int i=1; i<GTree.size(); i++){
        if(GTree[i].isleaf){
            leafNodes.push_back(i);
            //将叶节点i， 加入其祖先节点的 leafnodes列表中
            int current = i;
            while(current!=0){
                int father =  GTree[current].father;
                GTree[father].leafnodes.push_back(i);
                current = father;
            }
        }
        else{
            UpperNodes.push_back(i);
        }
    }
    //先对各叶子节点间的距离上下界进行计算
    cout<<"对各叶节点间距离上下界进行计算"<<endl;
    for(int idx=0; idx<leafNodes.size();idx++){
        //omp_set_num_threads(4);
//#pragma omp parallel for
        int leaf_id = leafNodes[idx];
        TIME_TICK_START
        mm = getMinMaxDistanceWithinLeaf_slow(leaf_id);
        cout<<"leaf"<<leaf_id<<" mm.max="<<mm.max;
        TIME_TICK_END
        TIME_TICK_PRINT("  ");

        //getchar();

    }

    //对叶节点与上层节点间的距离上下界进行计算
    //cout<<"对叶节点与上层节点间距离上下界进行计算"<<endl;


    TIME_TICK_END
    TIME_TICK_PRINT("test_outputNodeDisBound time:");
    exit(1);


}


void test_DisBoundBetweenLeaf(){

    //加载G-Tree索引
    ConfigType setting;
    AddConfigFromFile(setting,"../map_exp");
    //AddConfigFromCmdLine(setting,argc,argv);
    initGtree();

    // 加载 edge information
    loadEdgeMap();
    MinMaxD mm;


    cout<<"test_outputNodeDisBound map:"<<dataset_name<<endl;

    vector<int> leafNodes;
    vector<int> UpperNodes;

    map<int,map<int,MinMaxD>> DisBoundMap;
    map<int,map<int,bool>> DisBOOLMap;
    //map<int,map<int,bool >> DisBoundMap
    //取出叶子节点
    for(int i=1; i<GTree.size(); i++){
        if(GTree[i].isleaf){
            leafNodes.push_back(i);
            //将叶节点i， 加入其祖先节点的 leafnodes列表中
            int current = i;
            while(current!=0){
                int father =  GTree[current].father;
                GTree[father].leafnodes.push_back(i);
                current = father;
            }
        }
        else{
            UpperNodes.push_back(i);
        }
    }
    //先对各叶子节点间的距离上下界进行计算
    cout<<"对各叶节点间距离上下界进行计算"<<endl;
    for(int idx=0; idx<leafNodes.size();idx++){
        int leaf_id1 = leafNodes[idx];
        for(int idx2=0; idx2<leafNodes.size();idx2++){
            int leaf_id2 = leafNodes[idx2];
            if(idx==idx2) continue;
            TIME_TICK_START
            mm = getMinMaxDistanceBetweenLeaves_slow(leaf_id1, leaf_id2);
            cout<<"between leaf"<<leaf_id1<<" and "<<"leaf"<<leaf_id2<<"by new, mm.max="<<mm.max<<", mm.min"<<mm.min<<endl;
            TIME_TICK_END
            TIME_TICK_PRINT(" getMinMaxDistanceBetweenLeaves ");
            //getchar();
        }


    }

    //对叶节点与上层节点间的距离上下界进行计算
    //cout<<"对叶节点与上层节点间距离上下界进行计算"<<endl;


    TIME_TICK_END
    TIME_TICK_PRINT("test_outputNodeDisBound time:");
    exit(1);


}



void test_outputNodeDisBound(){
    //test_DisBoundWithinLeaf();
    test_DisBoundBetweenLeaf();

}


void test_MinMaxDistanceWithinLeaf(){
    TIME_TICK_START
    cout<<"test_MinMaxDistanceWithinLeaf:"<<dataset_name<<endl;

    vector<int> leafNodes;
    vector<int> UpperNodes;

    int node_size = GTree.size();

    //map<int,map<int,MinMaxD>> DisBoundMap;
    //map<int,map<int,bool>> DisBOOLMap;
    MinMaxD ** DisBoundMap = new MinMaxD*[node_size];
    for(int i=0;i<node_size;i++){
        DisBoundMap[i] = new MinMaxD[node_size];
    }

    bool ** DisBOOLMap = new bool*[node_size];
    for(int i=0;i<node_size;i++){
        DisBOOLMap[i] = new bool[node_size];
    }

    //取出叶子节点
    for(int i=1; i<GTree.size(); i++){
        if(GTree[i].isleaf){
            leafNodes.push_back(i);
            //将叶节点i， 加入其祖先节点的 leafnodes列表中
            int current = i;
            while(current!=0){
                int father =  GTree[current].father;
                GTree[father].leafnodes.push_back(i);
                current = father;
            }
        }
        else{
            UpperNodes.push_back(i);
        }
    }

    cout<<"对各叶节点间距离上下界进行计算"<<endl;

    for(int leaf:leafNodes){ //对各叶节点内部距离上下界计算
        MinMaxD mm;
        TIME_TICK_START
        mm = getMinMaxDistanceWithinLeaf_slow(leaf);
        DisBoundMap[leaf][leaf] = mm;
        DisBOOLMap[leaf][leaf] = true;
        TIME_TICK_END
        cout<<"(leaf"<<leaf<<",leaf"<<leaf<<")=(min:"<<mm.min<<",max:"<<mm.max<<")"<<endl;
        TIME_TICK_PRINT("getMinMaxDistanceWithinLeaf:")

    }




}



//bound map 测试函数
void test_BoundMap(int unode, int onode){
    float _dismin = MinDisBound[unode][onode];
    float _dismax = MaxDisBound[unode][onode];
    cout<<" _dismin="<< _dismin<<endl;
    cout<<" _dismax="<< _dismax<<endl;
}

void transferWholeSocialNetwork4OPC(unordered_map<int,vector<int>>& wholeSocialLinkMap, int type, OPIMC* _opp){
    //  对wholeSocialLinkMap进行遍历
    // 采用weight cascade 对边上概率赋值并将边信息重新转入AM
    cout<<"In function: transferWholeSocialNetwork..."<<endl;
    //map<int,set<int>>::iterator it;
    auto it = wholeSocialLinkMap.begin();
    int edgeNum =0;
    while(it != wholeSocialLinkMap.end()) {
        int u2 = it->first;
        vector<int> infriendsOfU2 = it->second;
        int size =0; size= infriendsOfU2.size();
        if(size>0){
            double prob = 1.0/ infriendsOfU2.size();
            for(int u1: infriendsOfU2){
                edgeNum++;
                _opp->addSocialLinkIMM(u1,u2,prob);

                if (edgeNum%100000==0) {
                    cout << "(u1, u2, degree, weight,  AM size till now, edges till now, mem) = " << u1 << ", " << u2 <<","<<size << ", " << prob <<endl;
                }
            }
            it ++;
        }else{
            it ++;
            continue;
        }

    }
    cout<<"social graph 导入完毕！"<<endl;
}



void test_OPIMC(unordered_map<int,vector<int>>& wholeGraph, vector<int>& users, float accuracy, int type){
    OPIMC* _opp;
    _opp = new OPIMC(wholeGraph, accuracy,  users.size(), -1, IM);

    transferWholeSocialNetwork4OPC(wholeGraph,IM,_opp);

    _opp -> InfluenceMaximize_by_opimc(5);
    delete _opp;
}

void test_InfluenceModel(map<int,set<int>>& wholeGraph, set<int>& users, int type){



}



//functions for exprimetal evaluation
float count_time = 0;
int repeat = 5;
float repeat_baseline = 1.0;
int repeat_random = 5;


enum test_target{QK, PC, Budget,KQ, ALPHA, QUERYSET, URatio, ORatio};

enum exp_task{Pre_Process, SPD_Test, DisBound, TextModify, VaryingObject, SocialPROCESS,GIMTree_Build, GIMTreeU_Build, GIMTreeO_Build, BuildO_NVD, Build_NVD, Pre_MC, TestCode, Topk_Test, Reverse_Single, Reverse_Batch, MaxInf,
    EXP, Effectiveness, Performance,PerformanceU,PerformanceO, SelectPOI};

enum varying_object{vary_u, vary_o};


string getPreprocessTask(std::string task_chose){
    if(task_chose == "clean")
        return "clean";
    else if (task_chose == "gtree")
        return "gtree";
    else if (task_chose == "phl")
        return "phl";
    else if (task_chose == "modify")
        return "modify";
}

int getExpChose(std::string exp_chose){
    if(exp_chose == "pre_process")
        return Pre_Process;
    else if (exp_chose== "spdist")
        return SPD_Test;
    else if(exp_chose == "disbound")
        return DisBound;
    else if(exp_chose == "text")
        return TextModify;
    else if(exp_chose == "varying") //增加双色体对象个数
        return VaryingObject;
    else if(exp_chose == "social")
        return SocialPROCESS;
    else if (exp_chose =="gimtree")  //直接对每个相关user进行nvd剪枝
        return GIMTree_Build;
    else if (exp_chose =="gimtreeU")  //直接对每个相关user进行nvd剪枝
        return GIMTreeU_Build;
    else if (exp_chose =="gimtreeO")  //直接对每个相关user进行nvd剪枝
        return GIMTreeO_Build;
    /*else if (exp_chose =="gimtreeU")  //直接对每个相关user进行nvd剪枝
        return GIMTree_Build;*/
    else if (exp_chose =="nvd")  //直接对每个相关user进行nvd剪枝
        return Build_NVD;
    else if (exp_chose =="nvdO")  //直接对每个相关user进行nvd剪枝
        return BuildO_NVD;

    else if (exp_chose =="premc")  //对整副社交网络图进行蒙特卡洛预采样
        return Pre_MC;

    else if(exp_chose =="test")
        return TestCode;

        //进行topk测试
    else if (exp_chose =="topk")  //直接对每个相关user进行nvd剪枝
        return Topk_Test;

    else if (exp_chose =="reverse")  //直接对每个相关user进行nvd剪枝
        return Reverse_Single;

    else if(exp_chose == "batch")
        return Reverse_Batch;
    else if(exp_chose == "maxinf")
        return MaxInf;

    else if(exp_chose == "exp")
        return EXP;

    else if(exp_chose =="effective")
        return Effectiveness;
    else if(exp_chose =="performance")
        return Performance;
    else if(exp_chose =="performanceU")
        return PerformanceU;
    else if(exp_chose =="performanceO")
        return PerformanceO;
    else if(exp_chose =="select")
        return SelectPOI;
}


string getExpResultsOutputPath(int exp_taskType){
    stringstream ss;
    if(exp_taskType==MaxInf)
        ss<<CODE_HOME<<dataset_name<<"/maxinf.exp";
    else if(exp_taskType==Reverse_Batch)
        ss<<CODE_HOME<<dataset_name<<"/batch.exp";
    else if(exp_taskType==Reverse_Single)
        ss<<CODE_HOME<<dataset_name<<"/single.exp";
    string path = ss.str();
    return path;
}



int getTargetParameterSetting(std::string testing_chose){
    if(testing_chose == "k")
        return  QK;
    else if (testing_chose== "pc")
        return PC;
    else if(testing_chose == "b")
        return Budget;
    else if(testing_chose == "g")
        return QUERYSET;
    else if(testing_chose == "alpha")
        return ALPHA;
    else if(testing_chose == "u")
        return URatio;
    else if(testing_chose == "o")
        return ORatio;

}






void test_MaxInfBRGSTkNN_effectiveness(int Qk, int KQ_size, int Q_size, float alpha, int alg, char flag){
    cout<<"执行effectiveness 测试"<<endl;
    vector<int> keyword_pool;
    vector<int> chainStores;
    vector<int> query_keywords; vector<int> check_in;
    chainStores = query_generation(keyword_pool,query_keywords,check_in, 200, KQ_size, Q_size);
    float a = DEFAULT_A;
    string file_pre;
    int tmp; float tmp2;
    stringstream ss2;
    file_pre = "effectiveness";
    ss2 << Q_size;


    ofstream exp;
    stringstream tt;
    tt<<FILE_EXP<<file_pre<<"_"<<dataset_name;
    exp.open(tt.str(), ios::app);
    exp<<file_pre<<"\t\t\t";
    //exp<<tmp<<"\t\t\t"<<endl;
    string ss ;  ss2>>ss;
    exp<<ss<<"\t\t\t"<<endl;


    exp.close();
}



void test_MaxInfBRGSTkNN_performance(int Qk, int KQ_size, int Q_size, float alpha, int alg, char flag){
    // 多查询处理
    cout<<"执行performance test!"<<endl;
    vector<int> keyword_pool;
    vector<int> chainStores;
    vector<int> query_keywords; vector<int> check_in;
    chainStores = query_generation(keyword_pool,query_keywords,check_in, 200, KQ_size, Q_size);
    float a = DEFAULT_A;
    string file_pre;
    int tmp; float tmp2;
    stringstream ss2;
    switch(flag) {
        case 'k':
            file_pre = "k";
            //tmp = Qk;
            ss2 << Qk;
            break;
        case 'c':
            file_pre = "QC";
            //tmp = Q_size;
            ss2 << Q_size;
            break;
        case 'v':
            file_pre = "KQ";
            //tmp = KQ_size;
            ss2 << KQ_size;
            break;
        case 'p':
            file_pre = "alpha";
            ss2 << alpha;
            break;
        default :
            break;

    }

    ofstream exp;
    stringstream tt;
    tt<<FILE_EXP<<file_pre<<"_"<<dataset_name;
    exp.open(tt.str(), ios::app);
    exp<<file_pre<<"\t\t\t";
    //exp<<tmp<<"\t\t\t"<<endl;
    string ss ;  ss2>>ss;
    exp<<ss<<"\t\t\t"<<endl;



    exp.close();

}



void test_scalability(){

}




#define TEST_H

#endif /
