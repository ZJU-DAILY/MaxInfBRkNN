//
// Created by jins on 10/31/19.
//

#include "cmath"
#include "config.h"
#include "bichromatic.h"
#include "road_network.h"
#include "gim_tree.h"
#include "utility.cc"
#include <iostream>
#include "IOcontroller.h"
#include <queue>
#include "diskbased.h"
#include "btree.cc"
#include "omp.h"
#include <boost/archive/binary_oarchive.hpp>
#include "Tenindra/utility/serialization.h"
#include "Tenindra/processing/V2PNVD.h"
#include "BTree/predefined.h"
#include "BTree/bpt.cc"
#include "Tenindra/processing/pruned_highway_labeling.h"


PrunedHighwayLabeling phl;


double getDistanceUpper_PHL(User user, POI poi){
    int v1 = user.Ni;
    int v2 = poi.Ni;
    double dis1 = user.dis; double dis2 = poi.dis;
    int spdist = phl.Query(v1,v2);
    //int spdist2 = SPSP(v1,v2);
    //cout<<"("<<v1<<","<<v2<<").spdist="<<spdist<<", spdist2="<<spdist2<<endl;
    double distance = spdist + dis1 + dis2;
    return distance;
    //(user.Ni,poi.Nj)/1.0+user.dis+poi.dis;
}


double getDistanceUpper_Oracle(User user, POI poi){

    //return getDistance_upper_Gtree(user, poi);

    return getDistanceUpper_PHL(user, poi);

    //return

}

double getDistanceLower_Oracle_PHL(User user, POI poi){

    //return getDistance_upper_Gtree(user, poi);
    double dis1 = phl.Query(user.Ni,poi.Ni)/WEIGHT_INFLATE_FACTOR+user.dis+poi.dis;
    double dis2 = phl.Query(user.Ni,poi.Nj)/WEIGHT_INFLATE_FACTOR+user.dis+(poi.dist-poi.dis);

    return min(dis1,dis2);

    //return

}


double getDistanceu2o_phl_old(User user, POI poi){

    //return getDistance_upper_Gtree(user, poi);
    double dis1 = phl.Query(user.Ni,poi.Ni)/WEIGHT_INFLATE_FACTOR+user.dis+poi.dis;
    double dis2 = phl.Query(user.Ni,poi.Nj)/WEIGHT_INFLATE_FACTOR+user.dis+(poi.dist-poi.dis);
    double dd1 = min(dis1,dis2);

    double dis3 = phl.Query(user.Nj,poi.Ni)/WEIGHT_INFLATE_FACTOR+(user.dist-user.dis)+poi.dis;
    double dis4 = phl.Query(user.Nj,poi.Nj)/WEIGHT_INFLATE_FACTOR+(user.dis-user.dis)+(poi.dist-poi.dis);
    double dd2 = min(dis3,dis4);

    double dist = min(dd1,dd2);
    return dist;

}


double getDistance_phl(int loc, User usr){
    int Ni = min(usr.Ni,usr.Nj);
    int Nj = max(usr.Ni,usr.Nj);
    double distance;
    double d1 = phl.Query(loc,Ni)/WEIGHT_INFLATE_FACTOR +usr.dis;
    double d2 = phl.Query(loc,Nj)/WEIGHT_INFLATE_FACTOR + (usr.dist - usr.dis);

    if(loc == Ni)
        return usr.dis;
    if(loc == Nj)
        return usr.dist-usr.dis;

    return min(d1,d2);


}

double getDistance_phl(int loc, POI poi){
    int Ni = min(poi.Ni,poi.Nj);
    int Nj = max(poi.Ni,poi.Nj);
    double distance;
    double d1 = phl.Query(loc,Ni)/WEIGHT_INFLATE_FACTOR +poi.dis;
    double d2 = phl.Query(loc,Nj)/WEIGHT_INFLATE_FACTOR + (poi.dist - poi.dis);

    if(loc == Ni)
        return poi.dis;
    if(loc == Nj)
        return poi.dist-poi.dis;

    return min(d1,d2);


}





#ifndef LOADINFO_H

//// ??????????????????unordered_map<int, int> addressHashMap_v2p;  //<vertex_term_key, address>
unordered_map<int, int> idHashMap_v2p_hybrid;      //<vertex_term_key, poi_id>
unordered_map<int, int> addressHashMap_nvd_hybrid;  //<poi_term_key, address>
unordered_map<int, vector<int>> idListHash_l2p_hybrid;  //???leaf node???????????? NN poi?????????

unordered_map<int, int> idHashMap_v2p;
unordered_map<int, int> addressHashMap_nvd;

bpt::bplus_tree openV2PBtree(){
    char btree_file_prefix[200];
    char nvd_btree_file_name[200];
    char v2p_btree_file_name[200];
    sprintf(btree_file_prefix,"../../../data/%s/NVD",dataset_name);
    sprintf(nvd_btree_file_name,"%s/NVD.bt",btree_file_prefix);
    sprintf(v2p_btree_file_name,"%s/v2p.bt",btree_file_prefix);

    bpt::bplus_tree v2p_btree(v2p_btree_file_name, false);
    return  v2p_btree;
}

bpt::bplus_tree openNVDADJBtree(){
    char btree_file_prefix[200];
    char nvd_btree_file_name[200];
    char v2p_btree_file_name[200];
    sprintf(btree_file_prefix,"../../../data/%s/NVD",dataset_name);
    sprintf(nvd_btree_file_name,"%s/NVD.bt",btree_file_prefix);
    //sprintf(v2p_btree_file_name,"%s/v2p.bt",btree_file_prefix);
    bpt::bplus_tree nvd_btree(nvd_btree_file_name, false);
    return nvd_btree;
}




//access NVD data on disk
//below is for NVD


int getV2PID_by_Hash(int vertex, int term){
    //clock_t startTime, endTime;

    //??????key?????????
    int key;
    key = term*VertexNum + vertex;

    int poi_id = -1 ;
    if(idHashMap_v2p.count(key)){
        poi_id = idHashMap_v2p[key];
    }
    return poi_id;
}

int getV2PID_by_HybridHash(int vertex, int term){
    //clock_t startTime, endTime;

    //??????key?????????
    int key;
    key = term*VertexNum + vertex;

    int poi_id = -1 ;
    if(idHashMap_v2p_hybrid.count(key)){
        poi_id = idHashMap_v2p_hybrid[key];
    }
    return poi_id;
}


vector<int> getNNPOIList_by_Hybrid_Hash(int leaf, int term ){
    //??????key?????????
    int leaf_key= -1;
    leaf_key = term*GTree.size() + leaf;
    //cout<<"Hybrid Hash leaf_key="<<leaf_key<<endl;
    vector<int> pois;
    if(idListHash_l2p_hybrid.count(leaf_key)){
        pois = idListHash_l2p_hybrid[leaf_key];
        //cout<<"find the entry whose key="<<leaf_key<<endl;
    }

    return pois;

}



int getNNPOI_By_Vertex_Keyword(int vertex, int term){

#ifdef TRACKV2P
    cout<<"????????????v"<<vertex<<"???????????????????????????t"<<term<<"???poi"<<endl;
#endif
    int nearest_poi_id = getV2PID_by_Hash(vertex, term);
#ifdef TRACKV2P
    cout<<"???poi???p"<<nearest_poi_id<<endl;
#endif

    return nearest_poi_id;
}

int getNNPOI_By_HybridVertex_Keyword(int vertex, int term){  //core

#ifdef TRACKV2P
    //cout<<"????????????v"<<vertex<<"???????????????????????????t"<<term<<"???poi"<<endl;
#endif
    clock_t startTime, endTime;
    startTime = clock();
    int nearest_poi_id = getV2PID_by_HybridHash(vertex, term);
    endTime = clock();
    double _time = (double)(endTime-startTime)/CLOCKS_PER_SEC*1000000;
#ifdef TIME_NVD
    cout<<"getNNPOI_By_HybridVertex_Keyword for v"<<vertex<<",t"<<term<<"time="<<_time<<" us!"<<endl;
#endif

#ifdef TRACKV2P
    cout<<"???poi???p"<<nearest_poi_id<<endl;
    cout<<"time="<<_time<<" us!"<<endl;
    ///getchar();
    if(nearest_poi_id==3144050){
        cout<<"find nearest_poi_id is 3144050!"<<endl;
    }
#endif



    return nearest_poi_id;
}



/// below is for IL-NVD


int getNVDAddress_by_Hash(int poi_id, int term){
    //??????key?????????
    clock_t startTime, endTime;
    startTime = clock();

    int key;

    key = term*poi_num + poi_id;

    int addr_logic = -1;
    if(addressHashMap_nvd_hybrid.count(key)){
        addr_logic = addressHashMap_nvd_hybrid[key];
    }
    else{
        printf("find no such key=%ld in NVD.idx\n",key);
        return -10;
    }
    endTime = clock();
#ifdef TIME_NVD
    cout << "getAddress in NVD.idx runtime is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC * 1000000 << "us" << endl;
#endif
    return  addr_logic;
}


vector<int> getPOIAdj_NVD_By_Keyword_vertexOnly(int poi_id, int term){


    vector<int> adjList_OfPOI;

    int addr_logic = -1;
    //addr_logic = getNVDAddress_by_Btree(poi_id,term);

    //int addr_logic2 = -1;
    addr_logic = getNVDAddress_by_Hash(poi_id,term);
    if(addr_logic==-10) return adjList_OfPOI;
    //printf("addr_logic=%d, addr_logic2=%d\n",addr_logic,addr_logic2);


    //read id(?????????????????????  joke_jins: (??????????????????)
    int read_id = -1;
    getIdxBlockByjins(NVD_ADJGraph,addr_logic, sizeof(int),Ref(read_id));
    addr_logic += sizeof(int);
    //read size
    int adj_size =-1;
    getIdxBlockByjins(NVD_ADJGraph,addr_logic, sizeof(int),Ref(adj_size));
    addr_logic += sizeof(int);
    //read each adjacent poi
    for(int i=0;i<adj_size;i++){
        int adj_poi_id = -1;
        getIdxBlockByjins(NVD_ADJGraph,addr_logic, sizeof(int),Ref(adj_poi_id));
        adjList_OfPOI.push_back(adj_poi_id);
        addr_logic += sizeof(int);
    }


    return adjList_OfPOI;
}


vector<int> getPOIAdj_NVD_By_Keyword(int poi_id, int term){


    vector<int> adjList_OfPOI;

    int addr_logic = -1;
    //addr_logic = getNVDAddress_by_Btree(poi_id,term);

    //int addr_logic2 = -1;
    addr_logic = getNVDAddress_by_Hash(poi_id,term);
    if(addr_logic==-10) return adjList_OfPOI;
    //printf("addr_logic=%d, addr_logic2=%d\n",addr_logic,addr_logic2);


    //read id(?????????????????????  joke_jins: (??????????????????)
    int read_id = -1;
    getIdxBlockByjins(NVD_ADJGraph,addr_logic, sizeof(int),Ref(read_id));
    addr_logic += sizeof(int);
    //read size
    int adj_size =-1;
    getIdxBlockByjins(NVD_ADJGraph,addr_logic, sizeof(int),Ref(adj_size));
    addr_logic += sizeof(int);
    //read each adjacent poi
    for(int i=0;i<adj_size;i++){
        int adj_poi_id = -1;
        getIdxBlockByjins(NVD_ADJGraph,addr_logic, sizeof(int),Ref(adj_poi_id));
        adjList_OfPOI.push_back(adj_poi_id);
        addr_logic += sizeof(int);
    }


    return adjList_OfPOI;
}



double getUpperDistance_b2P_phl(int border_id, POI poi){
    int Ni = min(poi.Ni,poi.Nj);
    int Nj = max(poi.Ni,poi.Nj);
    double distance_upper = -1;
    distance_upper = phl.Query(border_id,Ni)/WEIGHT_INFLATE_FACTOR +poi.dis;
    //double d2 = phl.Query(border_id,Nj)/WEIGHT_INFLATE_FACTOR + (poi.dist - poi.dis);

    return distance_upper;


}



double usrToPOIDistance_phl(User usr,POI poi){
    float distance =0;
    if(poi.Ni == usr.Ni && poi.Nj == usr.Nj){
        if (usr.dis > poi.dis)
            distance = usr.dis - poi.dis;
        else
            distance = poi.dis - usr.dis;
        //cout<<"??????"<<endl;
        return distance;
    }
    if(poi.Ni == usr.Ni && poi.Nj != usr.Nj){
        distance = usr.dis+poi.dis;
        //cout<<"??????"<<endl;
        return distance;
    }

    double d1 = getDistance_phl(poi.Ni,usr) + poi.dis;
    double d2 = getDistance_phl(poi.Nj,usr)+poi.dist-poi.dis;

    return min(d1,d2);
}



set<int> getBorder_SB_BatchNVD(BatchCheckEntry border_entry,int K, int a, double alpha){

    int border_id = border_entry.id;
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

    //?????????term??? sb_list??????????????????
    for(int term: Keys){  //???Keys?????????????????????????????????????????????sb_list(term)?????????
#ifdef TRACKBorderSB
        //cout<<"???t"<<term<<"???sb_list????????????"<<endl;
#endif

        map<int,bool> poiIsVist;
        priority_queue<nodeDistType_NVD> H_t;
        priority_queue<SBEntry> scoring_Bound_list_t;

        int posting_size = getTermOjectInvListSize(term);
        //??????term???????????????
        if(posting_size > posting_size_threshold){

            vector<POI> polled_object_lists;
            vector<float> polled_distance_lists;

#ifdef HybridHash   ////?????????????????????????????????????????????
            int nearest_poi_id = getNNPOI_By_HybridVertex_Keyword(border_id,term);  //???????????????vertex?????????????????????term???poi
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
                    cout<<"list??????o"<<poi_id2<<endl;
                    if(poi_id2==3144050){
                        cout<<leaf_node<<",t"<<term<<"3144050, size="<<pois.size()<<endl;
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
            int nearest_poi_id = getNNPOI_By_Vertex_Keyword(border_id,term);  //???????????????vertex?????????????????????term???poi
            if(nearest_poi_id ==-100) continue;   //???????????? ??????term?????????poi ??????user????????????
            POI poi = getPOIFromO2UOrgLeafData(nearest_poi_id);
            float dist_upper = getUpperDistance_b2P_phl(border_id,poi);
            nodeDistType_NVD tmp1(-1, nearest_poi_id, dist_upper, term);
            H_t.push(tmp1);
#endif


            while(scoring_Bound_list_t.size()<K){
                //H_t??????????????????
                nodeDistType_NVD _tmp = H_t.top();
                H_t.pop();
                int poi_id = _tmp.poi_id;
                if(poiIsVist.count(poi_id)) continue;  //poiIsVist[poi_id]==true
                //??????scoring lower bound
                double simT_LB = 0; double simD_UB = 0; double simS_LB = 1;
                simT_LB = tfIDF_term(term);
                simD_UB = _tmp.dist;
                double social_textual_LB = alpha*simS_LB + (1-alpha)*simT_LB;
                double scoring_LB = social_textual_LB / (1+simD_UB);
                SBEntry sb_entry(poi_id,simD_UB,scoring_LB);
                scoring_Bound_list_t.push(sb_entry);
                //sb_entry.printRlt();
                poiIsVist[poi_id]= true;
                //??????NVD graph????????????poi?????????adj_poi??????
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
        else{ //????????????????????????

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

    //??????rank_term???????????????????????????scoring_Bound_lists???????????????Key????????????
#ifdef TRACKBorderSB
    cout<<"??????rank_term???????????????????????????scoring_Bound_lists???????????????Key????????????"<<endl;
#endif
    while(rank_term.size()>0){
        GlobalSBEntry _termEntry = rank_term.top();
        rank_term.pop();
        int term_id = _termEntry.term_id;
        double sb_rk_t = _termEntry.sk_bound;

        //??????SB_upper(e, qo)
        double score_upper = -1;
        double max_relevance = tfIDF_termSet(_Keys);
        double min_distance = current_dist;

        double social_textual_maxscore = alpha*1 + (1.0-alpha)*max_relevance;
        score_upper = social_textual_maxscore / (1.0+ min_distance);
#ifdef TRACKBorderSB
        cout<<"for t"<<term_id<<", sb_rk_t="<<sb_rk_t<<",score_upper="<<score_upper<<endl;
#endif
        if(sb_rk_t > score_upper){ //remove term from _Keys
#ifdef TRACKBorderSB
            //cout<<"term"<<term_id<<"????????????"<<endl;
            //cout<<"??????_Keys:"; printSetElements(_Keys);
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
            //cout<<"????????????_Keys:"; printSetElements(_Keys);
#endif
        }
        //cout<<", for term"<<term_id<<",sb_rk_t="<<sb_rk_t;
    }
#ifdef TRACKBorderSB
    cout<<endl;
#endif

    return _Keys;

}




//?????????????????????
//map<int,set<int>> wholeSocialLinkMap;
//unordered_map<int,set<int>> wholeSocialLinkMap;
unordered_map<int,vector<int>> wholeSocialLinkMap;
//set<int> online_users;
vector<int> online_users;

GIM_Tree gim_tree;

//--------????????????
int num_D;
float STEP_SIZE;
float MAX_SEG_DIST;
float AVG_DEG;
#define MAXLEVEL 10

// to be initialized
char **cur_block;
int *cur_block_offset, *cur_node_maxkey, *cur_node_entries;
char *block;
int root_block_id,num_written_blocks,top_level;
int PtMaxKey=0;
bool PtMaxUsing=false;	// for solving the bug that was to find
float FACTOR;
bool IS_NODE_SIZE_GIVEN=false;


BTree* initialize(char *treename) {   //btree initialize
    BTree *bt;

    cur_block = new char*[MAXLEVEL]; //ptr to current node at each level
    cur_block_offset = new int[MAXLEVEL];
    cur_node_maxkey = new int[MAXLEVEL];
    cur_node_entries = new int[MAXLEVEL];
    top_level = 0;

    for (int i=0;i<MAXLEVEL; i++) cur_block[i] = NULL;
    block = new char[BlkLen];

    // number of cached file entries: 128
    bt = new BTree(treename, BlkLen, 128);
    i_capacity = (BlkLen - sizeof(char) - sizeof(int))/(sizeof(int)+sizeof(int));
    printf("i_capacity=%d\n", i_capacity);
    root_block_id = bt->root_ptr->block;
    num_written_blocks = 0;
    return  bt;
}

void addentry(BTree *bt,int *top_level,int capacity,int level,int key,int *block_id,int RawAddr=0) {
    if (cur_block[level] == NULL) { //new node to be created
        if ((*top_level) < level) //new root
            *top_level = level;
        cur_block[level] = new char[BlkLen];

        char l = (char)level;
        memcpy(cur_block[level], &l, sizeof(char));

        cur_block_offset[level]=sizeof(char)+sizeof(int);
        cur_node_entries[level] = 0;
    }
    cur_node_maxkey[level]= key;
    if ((level==1)&&PtMaxUsing) cur_node_maxkey[level]=max(PtMaxKey,key);	// new change !!!

    //copy key as new current entry and also the pointer to lower node
    memcpy(cur_block[level]+cur_block_offset[level], &key, sizeof(int));
    cur_block_offset[level]+=sizeof(int);

    //********* (Xblock_id for raw sequential page !)
    int Xblock_id=(level==1)?(RawAddr):(*block_id);
    memcpy(cur_block[level]+cur_block_offset[level], &Xblock_id, sizeof(int));

    cur_block_offset[level]+=sizeof(int);
    cur_node_entries[level]++;

    if (cur_node_entries[level] == capacity) { //node is full
        //copy capacity information
        memcpy(cur_block[level]+sizeof(char), &capacity, sizeof(int));
        //add maxkey of this node to upper level node
        bt->file->append_block(cur_block[level]);
        (*block_id)++;
        bt->num_of_inodes++;
        addentry(bt, top_level, capacity, level+1,
                 cur_node_maxkey[level], block_id);
        delete [] cur_block[level];
        cur_block[level] = NULL;
    }
}

void finalize(BTree* bt) {
    //flush non-empty blocks
    for (int level=1; level<= top_level; level++) {
        if (cur_block[level] != NULL) {
            //copy capacity information
            memcpy(cur_block[level]+sizeof(char), &cur_node_entries[level], sizeof(int));
            //add mbr of this node to upper level node
            if (level == top_level) {
                //root
                bt->file->write_block(cur_block[level], root_block_id);
                bt->num_of_inodes++;
                bt->root_ptr = NULL;
                bt->load_root();
                //printf("root written, id=%d\n", root_block_id);
            } else {
                bt->file->append_block(cur_block[level]);
                num_written_blocks++;
                bt->num_of_inodes++;
                addentry(bt, &top_level, i_capacity, level+1,
                         cur_node_maxkey[level], &num_written_blocks);
            }
            delete [] cur_block[level];
            cur_block[level] = NULL;
        }
    }
    delete [] block;
}



void BuildBinaryStorage(const char* fileprefix) {
    BlkLen=getBlockLength();
    char idxFileName[255];
    //FILE *ptFile,*edgeFile;
    sprintf(idxFileName,"%s.p_d",fileprefix);	remove(idxFileName); // remove existing file ?????????????????????
    /*ptFile=fopen(idxFileName,"w+");
    if(ptFile==NULL) cout<<"faile to open";
    sprintf(idxFileName,"%s.p_bt",fileprefix);	remove(idxFileName); // remove existing file ?????????????????????
    makePtFiles(ptFile,idxFileName);*/
    //cout <<"00000";
    //printf("????????????ad?????????????????????????????????\n");

    //fclose(ptFile);
}


//??????????????????Gtree initial???
int loadEdgeMap() {

    LOAD_START

    int id=0, Ni=0, Nj=0;
    int dist_scale=-1;float dist =-1;
    for (int i = 0; i < VertexNum; i++) {
        vector<int> emptyVector;
        vector<float> emtyVector2;
        adjList.push_back(emptyVector);
        adjWList.push_back(emtyVector2);
    }
    cout << "LOADING EDGEMAP...";


    //stringstream ss;
    //ss<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name <<".cedge";
    std::string edge_path = getRoadInputPath(EDGE);

    int edgeCnt = 0;
    ifstream finEdge;
    finEdge.open(edge_path.c_str());
    string str;
    while (getline(finEdge, str)) {
        istringstream tt(str);
        tt >> id >> Ni >> Nj >> dist;

        if (Ni < VertexNum && Nj < VertexNum) {    // ignore edges outside the range
            //printf( "%d %d %d %f\n", id, Ni, Nj, dist);
            edge e;
            e.dist = dist;
            e.Ni = _MIN(Ni, Nj);
            e.Nj = _MAX(Ni, Nj);    // enforce the constriant Ni<Nj
            adjList[Ni].push_back(Nj);
            adjList[Nj].push_back(Ni);
            //???????????????????????????
            adjWList[Ni].push_back(dist);
            adjWList[Nj].push_back(dist);

            EdgeMap[getKey(Ni, Nj)] = e;    // should be ok
            edgeCnt++;
        }
        //????????????????????????????????????(??????)??????
        int leafNode = Nodes[Ni].gtreepath.back();// Ni belong to the leaf : leafNode
        if(GTree[leafNode].maxEdgeDis==NULL)
            GTree[leafNode].maxEdgeDis = dist;
        else if(GTree[leafNode].maxEdgeDis<dist)
            GTree[leafNode].maxEdgeDis = dist;
        //???Ni?????????????????????????????????????????????
        GTree[leafNode].vetexSet.insert(Ni);
        int current = leafNode; int father;
        while(current!=0){
            father =  GTree[current].father;
            GTree[father].vetexSet.insert(Ni);
            current = father;
        }

        int leafNode2 = Nodes[Nj].gtreepath.back();// Nj belong to the leaf : leafNode
        if(GTree[leafNode2].maxEdgeDis==NULL)
            GTree[leafNode2].maxEdgeDis = dist;
        else if(GTree[leafNode2].maxEdgeDis<dist)
            GTree[leafNode2].maxEdgeDis = dist;
        //???Nj?????????????????????????????????????????????
        GTree[leafNode2].vetexSet.insert(Ni);
        current = leafNode2;
        while(current!=0){
            father =  GTree[current].father;
            GTree[father].vetexSet.insert(Nj);
            current = father;
        }
    }

    finEdge.close();
    cout << " DONE! EDGE #: " << edgeCnt;
    LOAD_END
    LOAD_PRINT()

    return 0;
}


int loadEdgeMap_light() {

    LOAD_START
    initGtree();
    int id=0, Ni=0, Nj=0;
    int dist_scale=-1;float dist =-1;
    for (int i = 0; i < VertexNum; i++) {
        vector<int> emptyVector;
        vector<float> emtyVector2;
        adjList.push_back(emptyVector);
        adjWList.push_back(emtyVector2);
    }
    memset(IsBorder,false, sizeof(IsBorder));
    cout << "LOADING EDGEMAP...";


    //stringstream ss;
    //ss<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name <<".cedge";
    std::string edge_path = getRoadInputPath(EDGE);

    int edgeCnt = 0;
    ifstream finEdge;
    finEdge.open(edge_path.c_str());
    string str;
    while (getline(finEdge, str)) {
        istringstream tt(str);
        tt >> id >> Ni >> Nj >> dist;

        if (Ni < VertexNum && Nj < VertexNum) {    // ignore edges outside the range
            //printf( "%d %d %d %f\n", id, Ni, Nj, dist);
            edge e;
            e.dist = dist;
            e.Ni = _MIN(Ni, Nj);
            e.Nj = _MAX(Ni, Nj);    // enforce the constriant Ni<Nj
            adjList[Ni].push_back(Nj);
            adjList[Nj].push_back(Ni);
            adjWList[Ni].push_back(dist);
            adjWList[Nj].push_back(dist);

            int Ni_leaf = Nodes[Ni].gtreepath.back();
            int Nj_leaf = Nodes[Nj].gtreepath.back();
            if(Ni_leaf!=Nj_leaf){
                IsBorder[Ni] = true;
                IsBorder[Nj] = true;
            }

            EdgeMap[getKey(Ni, Nj)] = e;    // should be ok
            int leafNode_Ni = Nodes[Ni].gtreepath.back();// Ni belong to the leaf : leafNode
            GTree[leafNode_Ni].vetexSet.insert(Ni);
            edgeCnt++;
        }
    }

    finEdge.close();
    cout << " DONE! EDGE #: " << edgeCnt;
    LOAD_END
    LOAD_PRINT()

    return 0;
}


//??????????????????????????????

int loadFriendShipData(){

    printf("LOADING FRIENDSHIP...");

    LOAD_START

    string linkPath = getSocialInputPath(LINK);
    FILE *fin;
    fin = fopen(linkPath.c_str(), "r");

    string str;
    //?????????IM??????????????????????????????
    // Graph vecGRev(graph_vertex_num);  //????????????, ????????????????????????
    //std::vector<size_t> vecInDeg(graph_vertex_num);  //???????????????
    //???????????????
    wholeSocialLinkMap.clear();online_users.clear();
    friendshipMap.clear(); followedMap.clear();
    friendShipTable.clear();

    set<int> socialNodeSet;  map<int, set<int>> link_tmp;

    int ui,uj;
    while (fscanf(fin,"%d %d ", &ui,&uj)==2) {

        if(400793==ui||uj==400793){
            cout<<endl;
        }
#ifdef LV
        if(ui<social_UserID_MaxKey ){
            socialNodeSet.insert(ui);
        }
        if(ui<social_UserID_MaxKey && uj<social_UserID_MaxKey){
            socialNodeSet.insert(uj);
            //if(ui<(UserID_MaxKey-1)&&uj<(UserID_MaxKey-1)){
            friendshipMap[ui].push_back(uj);
            followedMap[uj].push_back(ui);
            friendShipTable[ui][uj] =1;
            friendShipTable[uj][ui] =1;
            link_tmp[uj].insert(ui);
        }

#else
        socialNodeSet.insert(ui);
        socialNodeSet.insert(uj);
        friendshipMap[ui].push_back(uj);
        followedMap[uj].push_back(ui);
        friendShipTable[ui][uj] =1;
        friendShipTable[uj][ui] =1;
        link_tmp[uj].insert(ui);
        link_tmp[ui].insert(uj);
#endif

        //???IM???graph????????????
        float weight = 0.0;
        //vecGRev[uj].push_back(Edge(ui, weight));
    }

    fclose(fin);

    //??????set???list
    for(int u_id: socialNodeSet){
        online_users.push_back(u_id);
    }
    for(auto it=link_tmp.begin();it!=link_tmp.end();it++){
        int u = it->first;
        set<int> neighbors= it->second;
        for(int v:neighbors){
            wholeSocialLinkMap[u].push_back(v);
        }
    }

    printf("COMPLETE! ");

    LOAD_END
    LOAD_PRINT()
    return  0;
}

void loadFriendShipBinaryData(){
    printf("load SocialLink binary file...");
    clock_t startTime, endTime;
    startTime = clock();
    //???????????????
    wholeSocialLinkMap.clear();online_users.clear();
    friendshipMap.clear(); followedMap.clear();
    friendShipTable.clear();
    //???????????????????????????
    string dataPath = getSocialLinkBinaryInputPath();
    SocialGraphBinary Gs = serialization::getIndexFromBinaryFile<SocialGraphBinary>(dataPath);
    online_users = Gs.user_nodes; wholeSocialLinkMap = Gs._wholeSocialGraph;
    friendshipMap = Gs._friendshipMap; followedMap = Gs._followedMap;
    endTime = clock();



    //startTime = clock();
    for(auto it=wholeSocialLinkMap.begin();it!=wholeSocialLinkMap.end();it++){
        int u = it->first;
        vector<int> neighbors = it->second;
        for(int v: neighbors){
            friendShipTable[u][v] =1;
        }
    }

    endTime = clock();
    cout<<"?????? ??????????????????????????????????????????????????? ??????! ?????????"<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s" << endl;

    //cout<<"?????? friendShipTable ??????! ?????????"<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s" << endl;


}



//???????????????????????????checkin ??????(???????????????)
int loadCheckinData(){
    LOAD_START

    printf("LOADING CHECKIN...");

    FILE *fin;

    string checkinPath = getSocialInputPath(CHECKIN);
//#endif
    fin = fopen(checkinPath.c_str(),"r");
    //finFri.open(FILE_POI);

    string str;
    int poi,u,count;
    while (fscanf(fin,"%d %d %d", &poi,&u,&count)==3) {
        //istringstream tt(str);
        //tt>>poi>>u>>count;
        if(poi< (poi_num-1) && u<(UserID_MaxKey-1)){
            userCheckInIDList[u].push_back(poi);
            userCheckInCountList[u].push_back(count);


            poiCheckInIDList[poi].push_back(u);
            poiCheckInCountList[poi].push_back(count);

            pair<int, int> pcheck_info(make_pair(poi, count));
            userCheckInfoList[u].push_back(pcheck_info);
            pair<int, int> ucheck_info(make_pair(u, count));
            poiCheckInfoList[u].push_back(ucheck_info);

            userCheckInMap[u][poi] = count;
            poiCheckInMap[poi][u] = count;
            //POIs[poi].check_ins.push_back(u);
            //POIs[poi].check_ins.push_back(u);  ////????????????????????????????????????
        }
    }
    //finCk.close();
    fclose(fin);


    printf("COMPLETE!\t");
    LOAD_END
    LOAD_PRINT()
    //cout<<"p57's check-in:"<<endl;
    //printElements(poiCheckInList[57]);
    return  0;

}



void loadCheckinBinaryData(){
    printf("load Check-in binary file...");
    clock_t startTime, endTime;
    startTime = clock();
    //???????????????
    userCheckInfoList.clear(); poiCheckInfoList.clear();
    userCheckInMap.clear(); poiCheckInMap.clear();
    //???????????????????????????
    string dataPath = getCheckInBinaryInputPath();
    CheckInBinary checkIns = serialization::getIndexFromBinaryFile<CheckInBinary>(dataPath);
    userCheckInfoList = checkIns._userCheckInfoList;
    poiCheckInfoList = checkIns._poiCheckInfoList;


    for(auto it=userCheckInfoList.begin();it!=userCheckInfoList.end();it++){
        int u = it->first;
        vector<pair<int,int>> ucheckInfo = it->second;
        for(int i=0;i< ucheckInfo.size();i++){
            pair<int,int> ele = ucheckInfo[i];
            int p= ele.first;
            int count = ele.second;
            userCheckInMap[u][p] = count;
            poiCheckInMap[p][u] = count;
        }
    }

    endTime = clock();

    cout<<"?????? ??????????????????????????????????????? ??????! ?????????"<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s" << endl;



}

void loadCheckinBinaryDataAndRecord(){
    printf("load Check-in binary file...");
    clock_t startTime, endTime;
    startTime = clock();
    //???????????????
    userCheckInfoList.clear(); poiCheckInfoList.clear();
    userCheckInMap.clear(); poiCheckInMap.clear();
    //???????????????????????????
    string dataPath = getCheckInBinaryInputPath();
    CheckInBinary checkIns = serialization::getIndexFromBinaryFile<CheckInBinary>(dataPath);
    userCheckInfoList = checkIns._userCheckInfoList;
    poiCheckInfoList = checkIns._poiCheckInfoList;


    for(auto it=userCheckInfoList.begin();it!=userCheckInfoList.end();it++){
        int u = it->first;
        vector<pair<int,int>> ucheckInfo = it->second;
        for(int i=0;i< ucheckInfo.size();i++){
            pair<int,int> ele = ucheckInfo[i];
            int p= ele.first;
            int count = ele.second;
            userCheckInMap[u][p] = count;
            poiCheckInMap[p][u] = count;
            POIs[p].check_ins.push_back(u);  ////?????????
        }
    }


    endTime = clock();

    cout<<"?????? ??????????????????????????????????????? ??????! ?????????"<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s" << endl;



}



/*---------------------------------------1206-------------------------------*/
int popularPOI_Count(int u){
    //??????????????????????????????
    vector<int> friends = friendshipMap[u];
    map<int,int> FriendPOIMap;
    int optimal_Checkin = 0;
    for(int usr:friends){
        vector<int> user_Checkin = userCheckInIDList[usr];
        for(int poi: user_Checkin){
            if(userCheckInMap[usr][poi]!= 0){
                int checkin_Count = 1;
                if(FriendPOIMap[poi]!= 0){
                    int sum = FriendPOIMap[poi];
                    FriendPOIMap[poi] = sum + checkin_Count;

                } else
                    FriendPOIMap[poi] = checkin_Count;

                if(optimal_Checkin < FriendPOIMap[poi]){
                    optimal_Checkin = FriendPOIMap[poi];
                    //optimal_poi = poi;
                }
            }
        }
    }
    //??????????????????????????????
    usrMaxSS[u] = optimal_Checkin;
    return optimal_Checkin;
}


void outputUsrMaxSocial(){
    cout<<"Output UsrMaxSocial..."<<endl;
    cout<<"map:"<<dataset_name<<endl;
    ofstream output;
    stringstream tt;
    tt<<FILE_MAXSOCIAL;
    output.open(tt.str());
    int popular_poi_count =0;


    for(auto& a: friendshipMap){
        int u_id = a.first;
        popular_poi_count = popularPOI_Count(u_id);
        usrMaxSS[u_id] = popular_poi_count;
        //usr_FriendPOIMaxCheckin.push_back(popular_poi_count);
        output<< u_id << ' ' << popular_poi_count << endl;
    }

    output.close();
    cout<<"Output UsrMaxSocial COMPLETE!"<<endl;

}

void outputUsrMaxSocial(const char* fileprefix){
    cout<<"Output UsrMaxSocial..."<<endl;
    cout<<"map:"<<dataset_name<<endl;
    char smFile_Name[255];
    sprintf(smFile_Name,"%s.MAXSOCIAL",fileprefix);
    remove(smFile_Name); // remove existing file ?????????????????????
    ofstream output;
    stringstream tt;
    string smfile_name = smFile_Name;
    tt<<smfile_name;
    output.open(tt.str());



    int popular_poi_count =0;


    for(auto& a: friendshipMap){
        int u_id = a.first;
        popular_poi_count = popularPOI_Count(u_id);
        usrMaxSS[u_id] = popular_poi_count;
        //usr_FriendPOIMaxCheckin.push_back(popular_poi_count);
        output<< u_id << ' ' << popular_poi_count << endl;
    }

    output.close();
    cout<<"Output UsrMaxSocial COMPLETE!"<<endl;

}




void loadUsrMaxSocial(){
    LOAD_START
    cout<<"LOADING UsrMaxSocial...";
    ifstream load;
    stringstream tt;
    tt<<FILE_MAXSOCIAL;
    string path = tt.str();
    load.open(path);
    string str;
    int cnt = 0;
    while (getline(load, str)) {
        istringstream tt(str);
        int u_id, popular_poi_count;

        tt >> u_id >> popular_poi_count;
        usrMaxSS[u_id]= popular_poi_count;


        //getchar();
    }
    load.close();
    cout<<"COMPLETE! ";
    LOAD_END
    LOAD_PRINT()
}

unordered_map<int, set<int>> term_ONodeMap;
unordered_map<int, set<int>> term_UNodeMap;

////??????poi??? usr????????????????????????????????????????????????G-Tree???????????????
int loadUPMap() {

    LOAD_START
    // read UserMap and insert User to edges
    printf("LOADING USERMAP...");

    int User_id, User_Ni, User_Nj;
    float User_dist, User_dis;


    int usrKeyCnt = 0;

    ifstream finUsr;
    //finUsr.open(FILE_USER);
    string userPath;
    #ifdef LV
        userPath = getObjectInputPath(USER);
    #else
        userPath = getObjectNewInputPath(USER);
    #endif
    finUsr.open(userPath.c_str());

    //??????User??????
    string str;
    while (getline(finUsr, str)) {
        istringstream tt(str);
        tt >> User_id >> User_Ni >> User_Nj >> User_dist >> User_dis;
        float exchangeNum = (int) (User_dist * WEIGHT_INFLATE_FACTOR);
        User_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (User_dis * WEIGHT_INFLATE_FACTOR);
        User_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;

        //????????????leaf???user?????????????????????????????????
        int leafNode;
        leafNode = Nodes[User_Ni].gtreepath.back();// User_Ni belong to the leaf : leafNode
        if(GTree[leafNode].minDis == maximum)
            GTree[leafNode].minDis = User_dist;
        else if(GTree[leafNode].minDis > User_dist)
            GTree[leafNode].minDis = User_dist;
        // ???????????????????????????
        int leafNode2 = Nodes[User_Nj].gtreepath.back();// User_Ni belong to the leaf : leafNode
        if(GTree[leafNode2].minDis== maximum)
            GTree[leafNode2].minDis = User_dist;
        else if(GTree[leafNode2].minDis > User_dist)
            GTree[leafNode2].minDis = User_dist;

        //???????????????????????????????????????
        int popular_poi_count = 0;
            //?????????usr???????????????????????????????????????????????????
        if (usrMaxSS[User_id]>0){
            popular_poi_count = usrMaxSS[User_id];
        }
        else
            popular_poi_count = popularPOI_Count(User_id);
        for(int usr_node: Nodes[User_Ni].gtreepath){
            GTree[usr_node].maxSocialCount= max( GTree[leafNode].maxSocialCount, popular_poi_count);
            GTree[usr_node].minSocialCount= min( GTree[leafNode].minSocialCount, popular_poi_count);

        }


        //??????????????????????????????
        vector<int> uKey;
        int keyTmp;

        while (tt >> keyTmp) {
            if (-1 == keyTmp) {
                break;
            }
            usrKeyCnt++;
            uKey.push_back(keyTmp);
            //??????usr????????????
            invListOfUser[keyTmp].push_back(User_id);

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
        EdgeMap[getKey(User_Ni, User_Nj)].usrs.push_back(tmpNode);
        Users.push_back(tmpNode);

        usrCnt++;

        // create U I key for Gtree
        for (size_t i = 0; i < uKey.size(); i++){

            GTree[leafNode].userUKey.push_back(uKey[i]); // build U key;
            GTree[leafNode].userUKeySet.insert(uKey[i]);
            leafUsrInv[uKey[i]][leafNode].push_back(User_id);

            // ?????? ????????????????????????????????????

            usrTerm_leafSet[uKey[i]].insert(leafNode);  // ?????? (usr)??????????????????????????????


            for(int t:tmpNode.keywords)  // ??????uKey[i]???leafNode?????????????????????????????? <(term_id, leaf_id), {usr....}>
                term_leaf_KeySet[uKey[i]][leafNode].insert(t);


            //??????entry
            GTree[leafNode].inverted_list_u[uKey[i]].insert(User_id);
            term_UNodeMap[uKey[i]].insert(leafNode);
            //??????usr_cnt
            GTree[leafNode].term_usr_Map[uKey[i]].insert(User_id);
            GTree[leafNode].term_usrCnt[uKey[i]] = GTree[leafNode].term_usr_Map[uKey[i]].size();

            //??????leaf node???????????????????????????????????????????????????????????? //(???????????????hierarchy prune?????????)
            int current = leafNode; int father=-1;
            while(current!=0){
                father =  GTree[current].father;
                //???????????????????????????????????????
                GTree[father].userUKeySet.insert(uKey[i]);
                usrNodeTermChild[father][uKey[i]].insert(current); // father????????? uKey[i] ??????????????????????????????

                //??????entry
                GTree[father].inverted_list_u[uKey[i]].insert(current);
                term_UNodeMap[uKey[i]].insert(father);
                //??????usr_cnt
                GTree[father].term_usr_Map[uKey[i]].insert(User_id);
                GTree[father].term_usrCnt[uKey[i]] = GTree[father].term_usr_Map[uKey[i]].size();
                //??????term relevant usr_leaf
                set<int> _set = GTree[father].term_usrLeaf_Map[uKey[i]];
                if(_set.size()>0){
                    _set.insert(leafNode);
                    GTree[father].term_usrLeaf_Map[uKey[i]]= _set;
                } else{
                    _set.insert(leafNode);
                    GTree[father].term_usrLeaf_Map[uKey[i]]= _set;
                }

                current = father;
            }

        }
        GTree[leafNode].userSet.push_back(User_id);
        // pair userId, leafNode
        userLeafMap[User_id] = leafNode;
        // build user map
        userMap[User_id] = uKey;

        //?????????????????????????????????????????????
        vector<int> friends = friendshipMap[User_id];
        Users[User_id].friends = friends;

        for(int f:friends){
            GTree[leafNode].unionFriends.insert(f);
            //??????leaf node????????????????????????????????????
            int current = leafNode; int father;
            while(current!=0) {
                father = GTree[current].father;
                GTree[father].unionFriends.insert(f);
                current = father;
            }
        }


    }
    finUsr.close();
    cout << " DONE! USER #:" << usrCnt;
    cout << " , USER KEY #: " << usrKeyCnt << endl;
    LOAD_END
    LOAD_PRINT(" load User ")



    LOAD_START
    // load PoiMap
    printf("LOADING POIMAP...");
    ifstream finPoi;
    string poiPath;
#ifdef LV
    poiPath = getObjectInputPath(BUSINESS);
#else
    poiPath = getObjectNewInputPath(BUSINESS);
#endif

    //finPoi.open(FILE_POI);
    finPoi.open(poiPath.c_str());
    while (getline(finPoi, str)) {
        istringstream tt(str);
        int POI_id, POI_Ni, POI_Nj;
        float POI_dist, POI_dis;
        tt >> POI_id >> POI_Ni >> POI_Nj >> POI_dist >> POI_dis;
        float exchangeNum = (int) (POI_dist * WEIGHT_INFLATE_FACTOR);
        POI_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (POI_dis * WEIGHT_INFLATE_FACTOR);
        POI_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;

        //????????????leaf???poi?????????????????????????????????
        // ?????????????????????
        int leafNode;
        leafNode = Nodes[POI_Ni].gtreepath.back();// User_Ni belong to the leaf : leafNode

        if(GTree[leafNode].minDis == maximum)
            GTree[leafNode].minDis = POI_dist;
        else if(GTree[leafNode].minDis > POI_dist)
            GTree[leafNode].minDis = POI_dist;     //???????????????????????????????????????????????????
        // ??????????????????
        int leafNode2 = Nodes[POI_Nj].gtreepath.back();// User_Ni belong to the leaf : leafNode
        if(GTree[leafNode2].minDis == maximum)
            GTree[leafNode2].minDis = POI_dist;
        else if(GTree[leafNode2].minDis > POI_dist)
            GTree[leafNode2].minDis = POI_dist;

        vector<int> uKey;
        set<int> keywordSet;
        int keyTmp;
        while (tt >> keyTmp) {
            if (-1 == keyTmp) {
                break;
            }

            uKey.push_back(keyTmp);
            keywordSet.insert(keyTmp);
            //??????poi????????????

        }
        for(int term_id:keywordSet)
            invListOfPOI[term_id].push_back(POI_id);

        POI tmpPoi;
        tmpPoi.id = POI_id;
        tmpPoi.Ni = min(POI_Ni,POI_Nj);
        tmpPoi.Nj = max(POI_Ni,POI_Nj);
        tmpPoi.dist = POI_dist;
        tmpPoi.dis = POI_dis;
        tmpPoi.category = 1;
        tmpPoi.keywords = uKey;
        tmpPoi.keywordSet = keywordSet;
        if(poiCheckInIDList[POI_id].size()>0){
            tmpPoi.check_ins = poiCheckInIDList[POI_id];

        }

        EdgeMap[getKey(POI_Ni, POI_Nj)].pts.push_back(tmpPoi);
        POIs.push_back(tmpPoi);
        leafNode = Nodes[tmpPoi.Ni].gtreepath.back();

        //omp_set_num_threads(4);
//#pragma omp parallel for
        for (size_t i = 0; i < uKey.size(); i++){
            //???????????????????????????
            EdgeMap[getKey(POI_Ni, POI_Nj)].OUnionKeys.insert(uKey[i]);
            //??????GIM-tree???????????????????????????
            GTree[leafNode].inverted_list_o[uKey[i]].insert(POI_id);
            term_ONodeMap[uKey[i]].insert(leafNode);
            //??????term_object_map?????????
            GTree[leafNode].term_object_Map[uKey[i]].insert(POI_id);

            GTree[leafNode].term_poiCnt[uKey[i]] =
                    GTree[leafNode].term_object_Map[uKey[i]].size();

            int _tmp = GTree[leafNode].term_poiCnt[uKey[i]];


            //???????????????????????????
            GTree[leafNode].objectUKeySet.insert(uKey[i]);


            leafPoiInv[uKey[i]][leafNode].push_back(POI_id); // ?????? ?????????POI?????????????????????
            poiTerm_leafSet[uKey[i]].insert(leafNode);  // ?????? (poi)??????????????????????????????
            int current = leafNode; int father;
            while(current!=0){
                father =  GTree[current].father;
                GTree[father].objectUKeySet.insert(uKey[i]);
                //??????????????????????????????????????????
                if(true){
                        poiNodeTermChild[father][uKey[i]].insert(current); // father????????? uKey[i] ??????????????????????????????
                        //??????GIM-tree???????????????????????????
                        GTree[father].inverted_list_o[uKey[i]].insert(current);
                        term_ONodeMap[uKey[i]].insert(father);
                        //??????GIM-tree?????????term_object???cnt??????
                        GTree[father].term_object_Map[uKey[i]].insert(POI_id);
                        GTree[father].term_poiCnt[uKey[i]] =
                            GTree[father].term_object_Map[uKey[i]].size();


                }


                current = father;
            }

        }

        GTree[leafNode].poiSet.push_back(tmpPoi.id);
        poiCnt++;
    }
    finPoi.close();
    cout << " DONE! POI #: " << poiCnt << endl;
    LOAD_END
    LOAD_PRINT(" load POI ")

    return 0;
}



int loadUPMap_varying_object(ObjectType type, float rate) {

    LOAD_START


    ifstream finUsr;
    string userPath;
    ifstream finPoi;
    string poiPath;



    if(type==USER){
        userPath = getObjectVaryingInputPath(USER,rate);
        poiPath = getObjectNewInputPath(BUSINESS);
    }
    else if(type==BUSINESS){
        userPath = getUserTermUpdatingInputPath(rate);
        poiPath = getObjectVaryingInputPath(BUSINESS,rate);
    }


    finUsr.open(userPath.c_str());

    cout<<"userPath="<<userPath<<endl;
    cout<<"poiPath="<<poiPath<<endl;

    // read UserMap and insert User to edges
    cout<<"LOADING USERMAP..."<<endl;

    int User_id, User_Ni, User_Nj;
    float User_dist, User_dis;


    int usrKeyCnt = 0;


    Users.clear();
    //??????User??????
    string str;
    while (getline(finUsr, str)) {
        istringstream tt(str);
        tt >> User_id >> User_Ni >> User_Nj >> User_dist >> User_dis;
        float exchangeNum = (int) (User_dist * WEIGHT_INFLATE_FACTOR);
        User_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (User_dis * WEIGHT_INFLATE_FACTOR);
        User_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;

        //????????????leaf???user?????????????????????????????????
        int leafNode;
        leafNode = Nodes[User_Ni].gtreepath.back();// User_Ni belong to the leaf : leafNode
        if(GTree[leafNode].minDis == maximum)
            GTree[leafNode].minDis = User_dist;
        else if(GTree[leafNode].minDis > User_dist)
            GTree[leafNode].minDis = User_dist;
        // ???????????????????????????
        int leafNode2 = Nodes[User_Nj].gtreepath.back();// User_Ni belong to the leaf : leafNode
        if(GTree[leafNode2].minDis== maximum)
            GTree[leafNode2].minDis = User_dist;
        else if(GTree[leafNode2].minDis > User_dist)
            GTree[leafNode2].minDis = User_dist;

        //???????????????????????????????????????
        int popular_poi_count = 0;
        //?????????usr???????????????????????????????????????????????????
        if (usrMaxSS[User_id]>0){
            popular_poi_count = usrMaxSS[User_id];
        }
        else
            popular_poi_count = popularPOI_Count(User_id);
        for(int usr_node: Nodes[User_Ni].gtreepath){
            GTree[usr_node].maxSocialCount= max( GTree[leafNode].maxSocialCount, popular_poi_count);
            GTree[usr_node].minSocialCount= min( GTree[leafNode].minSocialCount, popular_poi_count);

        }


        //??????????????????????????????
        vector<int> uKey;
        int keyTmp;

        while (tt >> keyTmp) {
            if (-1 == keyTmp) {
                break;
            }
            usrKeyCnt++;
            uKey.push_back(keyTmp);
            //??????usr????????????
            invListOfUser[keyTmp].push_back(User_id);

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
        EdgeMap[getKey(User_Ni, User_Nj)].usrs.push_back(tmpNode);
        Users.push_back(tmpNode);

        usrCnt++;

        // create U I key for Gtree
        for (size_t i = 0; i < uKey.size(); i++){

            GTree[leafNode].userUKey.push_back(uKey[i]); // build U key;
            GTree[leafNode].userUKeySet.insert(uKey[i]);
            leafUsrInv[uKey[i]][leafNode].push_back(User_id);
            // ?????? ????????????????????????????????????

            usrTerm_leafSet[uKey[i]].insert(leafNode);  // ?????? (usr)??????????????????????????????


            for(int t:tmpNode.keywords)  // ??????uKey[i]???leafNode?????????????????????????????? <(term_id, leaf_id), {usr....}>
                term_leaf_KeySet[uKey[i]][leafNode].insert(t);


            //??????entry
            GTree[leafNode].inverted_list_u[uKey[i]].insert(User_id);
            term_UNodeMap[uKey[i]].insert(leafNode);
            //??????usr_cnt
            GTree[leafNode].term_usr_Map[uKey[i]].insert(User_id);
            GTree[leafNode].term_usrCnt[uKey[i]] = GTree[leafNode].term_usr_Map[uKey[i]].size();

            //??????leaf node???????????????????????????????????????????????????????????? //(???????????????hierarchy prune?????????)
            int current = leafNode; int father=-1;
            while(current!=0){
                father =  GTree[current].father;
                //???????????????????????????????????????
                GTree[father].userUKeySet.insert(uKey[i]);
                usrNodeTermChild[father][uKey[i]].insert(current); // father????????? uKey[i] ??????????????????????????????

                //??????entry
                GTree[father].inverted_list_u[uKey[i]].insert(current);
                term_UNodeMap[uKey[i]].insert(father);
                //??????usr_cnt
                GTree[father].term_usr_Map[uKey[i]].insert(User_id);
                GTree[father].term_usrCnt[uKey[i]] = GTree[father].term_usr_Map[uKey[i]].size();
                //??????term relevant usr_leaf
                set<int> _set = GTree[father].term_usrLeaf_Map[uKey[i]];
                if(_set.size()>0){
                    _set.insert(leafNode);
                    GTree[father].term_usrLeaf_Map[uKey[i]]= _set;
                } else{
                    _set.insert(leafNode);
                    GTree[father].term_usrLeaf_Map[uKey[i]]= _set;
                }

                current = father;
            }

        }
        GTree[leafNode].userSet.push_back(User_id);
        // pair userId, leafNode
        userLeafMap[User_id] = leafNode;
        // build user map
        userMap[User_id] = uKey;

        //?????????????????????????????????????????????
        vector<int> friends = friendshipMap[User_id];
        Users[User_id].friends = friends;

        for(int f:friends){
            GTree[leafNode].unionFriends.insert(f);
            //??????leaf node????????????????????????????????????
            int current = leafNode; int father;
            while(current!=0) {
                father = GTree[current].father;
                GTree[father].unionFriends.insert(f);
                current = father;
            }
        }


    }
    finUsr.close();
    cout << " DONE! USER SIZE #:" << usrCnt;
    cout << " , USER KEY #: " << usrKeyCnt << endl;
    LOAD_END
    LOAD_PRINT(" load User ")



    LOAD_START
    // load PoiMap
    cout<<"LOADING POIMAP..."<<endl;


    POIs.clear();
    finPoi.open(poiPath.c_str());
    while (getline(finPoi, str)) {
        istringstream tt(str);
        int POI_id, POI_Ni, POI_Nj;
        float POI_dist, POI_dis;
        tt >> POI_id >> POI_Ni >> POI_Nj >> POI_dist >> POI_dis;
        float exchangeNum = (int) (POI_dist * WEIGHT_INFLATE_FACTOR);
        POI_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (POI_dis * WEIGHT_INFLATE_FACTOR);
        POI_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;

        //????????????leaf???poi?????????????????????????????????
        // ?????????????????????
        int leafNode;
        leafNode = Nodes[POI_Ni].gtreepath.back();// User_Ni belong to the leaf : leafNode

        if(GTree[leafNode].minDis == maximum)
            GTree[leafNode].minDis = POI_dist;
        else if(GTree[leafNode].minDis > POI_dist)
            GTree[leafNode].minDis = POI_dist;     //???????????????????????????????????????????????????
        // ??????????????????
        int leafNode2 = Nodes[POI_Nj].gtreepath.back();// User_Ni belong to the leaf : leafNode
        if(GTree[leafNode2].minDis == maximum)
            GTree[leafNode2].minDis = POI_dist;
        else if(GTree[leafNode2].minDis > POI_dist)
            GTree[leafNode2].minDis = POI_dist;

        vector<int> uKey;
        set<int> keywordSet;
        int keyTmp;
        while (tt >> keyTmp) {
            if (-1 == keyTmp) {
                break;
            }

            uKey.push_back(keyTmp);
            keywordSet.insert(keyTmp);
            //??????poi????????????

        }
        for(int term_id:keywordSet)
            invListOfPOI[term_id].push_back(POI_id);

        POI tmpPoi;
        tmpPoi.id = POI_id;
        tmpPoi.Ni = min(POI_Ni,POI_Nj);
        tmpPoi.Nj = max(POI_Ni,POI_Nj);
        tmpPoi.dist = POI_dist;
        tmpPoi.dis = POI_dis;
        tmpPoi.category = 1;
        tmpPoi.keywords = uKey;
        tmpPoi.keywordSet = keywordSet;
        if(poiCheckInIDList[POI_id].size()>0){
            tmpPoi.check_ins = poiCheckInIDList[POI_id];

        }

        EdgeMap[getKey(POI_Ni, POI_Nj)].pts.push_back(tmpPoi);
        POIs.push_back(tmpPoi);
        leafNode = Nodes[tmpPoi.Ni].gtreepath.back();

        //omp_set_num_threads(4);
//#pragma omp parallel for
        for (size_t i = 0; i < uKey.size(); i++){
            //???????????????????????????
            EdgeMap[getKey(POI_Ni, POI_Nj)].OUnionKeys.insert(uKey[i]);
            //??????GIM-tree???????????????????????????
            GTree[leafNode].inverted_list_o[uKey[i]].insert(POI_id);
            term_ONodeMap[uKey[i]].insert(leafNode);
            //??????term_object_map?????????
            GTree[leafNode].term_object_Map[uKey[i]].insert(POI_id);

            GTree[leafNode].term_poiCnt[uKey[i]] =
                    GTree[leafNode].term_object_Map[uKey[i]].size();

            int _tmp = GTree[leafNode].term_poiCnt[uKey[i]];


            //???????????????????????????
            GTree[leafNode].objectUKeySet.insert(uKey[i]);


            leafPoiInv[uKey[i]][leafNode].push_back(POI_id); // ?????? ?????????POI?????????????????????
            poiTerm_leafSet[uKey[i]].insert(leafNode);  // ?????? (poi)??????????????????????????????
            int current = leafNode; int father;
            while(current!=0){
                father =  GTree[current].father;
                GTree[father].objectUKeySet.insert(uKey[i]);
                //??????????????????????????????????????????
                if(true){
                    poiNodeTermChild[father][uKey[i]].insert(current); // father????????? uKey[i] ??????????????????????????????
                    //??????GIM-tree???????????????????????????
                    GTree[father].inverted_list_o[uKey[i]].insert(current);
                    term_ONodeMap[uKey[i]].insert(father);
                    //??????GIM-tree?????????term_object???cnt??????
                    GTree[father].term_object_Map[uKey[i]].insert(POI_id);
                    GTree[father].term_poiCnt[uKey[i]] =
                            GTree[father].term_object_Map[uKey[i]].size();


                }


                current = father;
            }

        }

        GTree[leafNode].poiSet.push_back(tmpPoi.id);
        poiCnt++;
    }
    finPoi.close();
    cout << " DONE! POI SIZE =" << poiCnt << endl;
    LOAD_END
    LOAD_PRINT(" load POI ")

    return 0;
}



int loadUPMap_LasVegas() {

    LOAD_START
    // read UserMap and insert User to edges
    printf("LOADING USERMAP for LasVegas...");

    int User_id, User_Ni, User_Nj;
    float User_dist, User_dis;


    int usrKeyCnt = 0;

    //???????????????
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

    //??????User??????
    string str;
    while (getline(finUsr, str)) {
        istringstream tt(str);
        tt >> User_id >> User_Ni >> User_Nj >> User_dist >> User_dis;
        float exchangeNum = (int) (User_dist * WEIGHT_INFLATE_FACTOR);
        User_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (User_dis * WEIGHT_INFLATE_FACTOR);
        User_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;

        //????????????leaf???user?????????????????????????????????
        int leafNode;
        leafNode = Nodes[User_Ni].gtreepath.back();// User_Ni belong to the leaf : leafNode
        if(GTree[leafNode].minDis == maximum)
            GTree[leafNode].minDis = User_dist;
        else if(GTree[leafNode].minDis > User_dist)
            GTree[leafNode].minDis = User_dist;
        // ???????????????????????????
        int leafNode2 = Nodes[User_Nj].gtreepath.back();// User_Ni belong to the leaf : leafNode
        if(GTree[leafNode2].minDis== maximum)
            GTree[leafNode2].minDis = User_dist;
        else if(GTree[leafNode2].minDis > User_dist)
            GTree[leafNode2].minDis = User_dist;

        //???????????????????????????????????????
        int popular_poi_count = 0;
        //?????????usr???????????????????????????????????????????????????
        if (usrMaxSS[User_id]>0){
            popular_poi_count = usrMaxSS[User_id];
        }
        else
            popular_poi_count = popularPOI_Count(User_id);
        for(int usr_node: Nodes[User_Ni].gtreepath){
            GTree[usr_node].maxSocialCount= max( GTree[leafNode].maxSocialCount, popular_poi_count);
            GTree[usr_node].minSocialCount= min( GTree[leafNode].minSocialCount, popular_poi_count);

        }


        //??????????????????????????????
        vector<int> uKey;
        int keyTmp;

        while (tt >> keyTmp) {
            if (-1 == keyTmp) {
                break;
            }
            usrKeyCnt++;
            uKey.push_back(keyTmp);
            //??????usr????????????
            invListOfUser[keyTmp].push_back(User_id);

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
        EdgeMap[getKey(User_Ni, User_Nj)].usrs.push_back(tmpNode);
        Users[tmpNode.id] =tmpNode;

        //printUsrInfo(tmpNode);

        usrCnt++;

        // create U I key for Gtree
        for (size_t i = 0; i < uKey.size(); i++){

            GTree[leafNode].userUKey.push_back(uKey[i]); // build U key;
            GTree[leafNode].userUKeySet.insert(uKey[i]);
            leafUsrInv[uKey[i]][leafNode].push_back(User_id);
            // ?????? ????????????????????????????????????

            usrTerm_leafSet[uKey[i]].insert(leafNode);  // ?????? (usr)??????????????????????????????


            for(int t:tmpNode.keywords)  // ??????uKey[i]???leafNode?????????????????????????????? <(term_id, leaf_id), {usr....}>
                term_leaf_KeySet[uKey[i]][leafNode].insert(t);


            //??????entry
            GTree[leafNode].inverted_list_u[uKey[i]].insert(User_id);
            term_UNodeMap[uKey[i]].insert(leafNode);
            //??????usr_cnt
            GTree[leafNode].term_usr_Map[uKey[i]].insert(User_id);
            GTree[leafNode].term_usrCnt[uKey[i]] = GTree[leafNode].term_usr_Map[uKey[i]].size();

            //??????leaf node???????????????????????????????????????????????????????????? //(???????????????hierarchy prune?????????)
            int current = leafNode; int father=-1;
            while(current!=0){
                father =  GTree[current].father;
                //???????????????????????????????????????
                GTree[father].userUKeySet.insert(uKey[i]);
                usrNodeTermChild[father][uKey[i]].insert(current); // father????????? uKey[i] ??????????????????????????????

                //??????entry
                GTree[father].inverted_list_u[uKey[i]].insert(current);
                term_UNodeMap[uKey[i]].insert(father);
                //??????usr_cnt
                GTree[father].term_usr_Map[uKey[i]].insert(User_id);
                GTree[father].term_usrCnt[uKey[i]] = GTree[father].term_usr_Map[uKey[i]].size();
                //??????term relevant usr_leaf
                set<int> _set = GTree[father].term_usrLeaf_Map[uKey[i]];
                if(_set.size()>0){
                    _set.insert(leafNode);
                    GTree[father].term_usrLeaf_Map[uKey[i]]= _set;
                } else{
                    _set.insert(leafNode);
                    GTree[father].term_usrLeaf_Map[uKey[i]]= _set;
                }

                current = father;
            }

        }
        GTree[leafNode].userSet.push_back(User_id);
        // pair userId, leafNode
        userLeafMap[User_id] = leafNode;
        // build user map
        userMap[User_id] = uKey;

        //?????????????????????????????????????????????
        vector<int> friends = friendshipMap[User_id];
        Users[User_id].friends = friends;

        for(int f:friends){
            GTree[leafNode].unionFriends.insert(f);
            //??????leaf node????????????????????????????????????
            int current = leafNode; int father;
            while(current!=0) {
                father = GTree[current].father;
                GTree[father].unionFriends.insert(f);
                current = father;
            }
        }


    }
    finUsr.close();
    cout << " DONE! USER #:" << usrCnt;
    cout << " , USER KEY #: " << usrKeyCnt << endl;
    LOAD_END
    LOAD_PRINT(" load User ")



    LOAD_START
    // load PoiMap
    printf("LOADING POIMAP...");
    ifstream finPoi;
    string poiPath;
#ifdef LV
    poiPath = getObjectInputPath(BUSINESS);
#else
    poiPath = getObjectNewInputPath(BUSINESS);
#endif

    //finPoi.open(FILE_POI);
    POIs.clear();
    finPoi.open(poiPath.c_str());
    while (getline(finPoi, str)) {
        istringstream tt(str);
        int POI_id, POI_Ni, POI_Nj;
        float POI_dist, POI_dis;
        tt >> POI_id >> POI_Ni >> POI_Nj >> POI_dist >> POI_dis;
        float exchangeNum = (int) (POI_dist * WEIGHT_INFLATE_FACTOR);
        POI_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (POI_dis * WEIGHT_INFLATE_FACTOR);
        POI_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;

        //????????????leaf???poi?????????????????????????????????
        // ?????????????????????
        int leafNode;
        leafNode = Nodes[POI_Ni].gtreepath.back();// User_Ni belong to the leaf : leafNode

        if(GTree[leafNode].minDis == maximum)
            GTree[leafNode].minDis = POI_dist;
        else if(GTree[leafNode].minDis > POI_dist)
            GTree[leafNode].minDis = POI_dist;     //???????????????????????????????????????????????????
        // ??????????????????
        int leafNode2 = Nodes[POI_Nj].gtreepath.back();// User_Ni belong to the leaf : leafNode
        if(GTree[leafNode2].minDis == maximum)
            GTree[leafNode2].minDis = POI_dist;
        else if(GTree[leafNode2].minDis > POI_dist)
            GTree[leafNode2].minDis = POI_dist;

        vector<int> uKey;
        set<int> keywordSet;
        int keyTmp;
        while (tt >> keyTmp) {
            if (-1 == keyTmp) {
                break;
            }

            uKey.push_back(keyTmp);
            keywordSet.insert(keyTmp);
            //??????poi????????????

        }
        for(int term_id:keywordSet)
            invListOfPOI[term_id].push_back(POI_id);

        POI tmpPoi;
        tmpPoi.id = POI_id;
        tmpPoi.Ni = min(POI_Ni,POI_Nj);
        tmpPoi.Nj = max(POI_Ni,POI_Nj);
        tmpPoi.dist = POI_dist;
        tmpPoi.dis = POI_dis;
        tmpPoi.category = 1;
        tmpPoi.keywords = uKey;
        tmpPoi.keywordSet = keywordSet;
        if(poiCheckInIDList[POI_id].size()>0){
            tmpPoi.check_ins = poiCheckInIDList[POI_id];

        }

        EdgeMap[getKey(POI_Ni, POI_Nj)].pts.push_back(tmpPoi);
        POIs.push_back(tmpPoi);
        leafNode = Nodes[tmpPoi.Ni].gtreepath.back();

        //omp_set_num_threads(4);
//#pragma omp parallel for
        for (size_t i = 0; i < uKey.size(); i++){
            //???????????????????????????
            EdgeMap[getKey(POI_Ni, POI_Nj)].OUnionKeys.insert(uKey[i]);
            //??????GIM-tree???????????????????????????
            GTree[leafNode].inverted_list_o[uKey[i]].insert(POI_id);
            term_ONodeMap[uKey[i]].insert(leafNode);
            //??????term_object_map?????????
            GTree[leafNode].term_object_Map[uKey[i]].insert(POI_id);

            GTree[leafNode].term_poiCnt[uKey[i]] =
                    GTree[leafNode].term_object_Map[uKey[i]].size();

            int _tmp = GTree[leafNode].term_poiCnt[uKey[i]];


            //???????????????????????????
            GTree[leafNode].objectUKeySet.insert(uKey[i]);


            leafPoiInv[uKey[i]][leafNode].push_back(POI_id); // ?????? ?????????POI?????????????????????
            poiTerm_leafSet[uKey[i]].insert(leafNode);  // ?????? (poi)??????????????????????????????
            int current = leafNode; int father;
            while(current!=0){
                father =  GTree[current].father;
                GTree[father].objectUKeySet.insert(uKey[i]);
                //??????????????????????????????????????????
                if(true){
                    poiNodeTermChild[father][uKey[i]].insert(current); // father????????? uKey[i] ??????????????????????????????
                    //??????GIM-tree???????????????????????????
                    GTree[father].inverted_list_o[uKey[i]].insert(current);
                    term_ONodeMap[uKey[i]].insert(father);
                    //??????GIM-tree?????????term_object???cnt??????
                    GTree[father].term_object_Map[uKey[i]].insert(POI_id);
                    GTree[father].term_poiCnt[uKey[i]] =
                            GTree[father].term_object_Map[uKey[i]].size();


                }


                current = father;
            }

        }

        GTree[leafNode].poiSet.push_back(tmpPoi.id);
        poiCnt++;
    }
    finPoi.close();
    cout << " DONE! POI #: " << poiCnt << endl;
    LOAD_END
    LOAD_PRINT(" load POI ")

    return 0;
}

int loadUPMap_Twitter() {

    LOAD_START
    // read UserMap and insert User to edges
#ifdef Twitter
    printf("LOADING USERMAP for Twitter...");
#endif

    int User_id, User_Ni, User_Nj;
    float User_dist, User_dis;


    int usrKeyCnt = 0;

    //???????????????
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

    //??????User??????
    string str;
    while (getline(finUsr, str)) {
        istringstream tt(str);
        tt >> User_id >> User_Ni >> User_Nj >> User_dist >> User_dis;
        float exchangeNum = (int) (User_dist * WEIGHT_INFLATE_FACTOR);
        User_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (User_dis * WEIGHT_INFLATE_FACTOR);
        User_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;

        //????????????leaf???user?????????????????????????????????
        int leafNode;
        leafNode = Nodes[User_Ni].gtreepath.back();// User_Ni belong to the leaf : leafNode
        if(GTree[leafNode].minDis == maximum)
            GTree[leafNode].minDis = User_dist;
        else if(GTree[leafNode].minDis > User_dist)
            GTree[leafNode].minDis = User_dist;
        // ???????????????????????????
        int leafNode2 = Nodes[User_Nj].gtreepath.back();// User_Ni belong to the leaf : leafNode
        if(GTree[leafNode2].minDis== maximum)
            GTree[leafNode2].minDis = User_dist;
        else if(GTree[leafNode2].minDis > User_dist)
            GTree[leafNode2].minDis = User_dist;

        //???????????????????????????????????????
        int popular_poi_count = 0;
        //?????????usr???????????????????????????????????????????????????
        if (usrMaxSS[User_id]>0){
            popular_poi_count = usrMaxSS[User_id];
        }
        else
            popular_poi_count = popularPOI_Count(User_id);
        for(int usr_node: Nodes[User_Ni].gtreepath){
            GTree[usr_node].maxSocialCount= max( GTree[leafNode].maxSocialCount, popular_poi_count);
            GTree[usr_node].minSocialCount= min( GTree[leafNode].minSocialCount, popular_poi_count);

        }


        //??????????????????????????????
        vector<int> uKey;
        int keyTmp;

        while (tt >> keyTmp) {
            if (-1 == keyTmp) {
                break;
            }
            usrKeyCnt++;
            uKey.push_back(keyTmp);
            //??????usr????????????
            invListOfUser[keyTmp].push_back(User_id);

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
        EdgeMap[getKey(User_Ni, User_Nj)].usrs.push_back(tmpNode);
        Users[tmpNode.id] =tmpNode;

        if(tmpNode.id%1000==0)
            printUsrInfo(tmpNode);

        usrCnt++;

        // create U I key for Gtree
        for (size_t i = 0; i < uKey.size(); i++){

            GTree[leafNode].userUKey.push_back(uKey[i]); // build U key;
            GTree[leafNode].userUKeySet.insert(uKey[i]);
            leafUsrInv[uKey[i]][leafNode].push_back(User_id);
            // ?????? ????????????????????????????????????

            usrTerm_leafSet[uKey[i]].insert(leafNode);  // ?????? (usr)??????????????????????????????


            for(int t:tmpNode.keywords)  // ??????uKey[i]???leafNode?????????????????????????????? <(term_id, leaf_id), {usr....}>
                term_leaf_KeySet[uKey[i]][leafNode].insert(t);


            //??????entry
            GTree[leafNode].inverted_list_u[uKey[i]].insert(User_id);
            term_UNodeMap[uKey[i]].insert(leafNode);
            //??????usr_cnt
            GTree[leafNode].term_usr_Map[uKey[i]].insert(User_id);
            GTree[leafNode].term_usrCnt[uKey[i]] = GTree[leafNode].term_usr_Map[uKey[i]].size();

            //??????leaf node???????????????????????????????????????????????????????????? //(???????????????hierarchy prune?????????)
            int current = leafNode; int father=-1;
            while(current!=0){
                father =  GTree[current].father;
                //???????????????????????????????????????
                GTree[father].userUKeySet.insert(uKey[i]);
                usrNodeTermChild[father][uKey[i]].insert(current); // father????????? uKey[i] ??????????????????????????????

                //??????entry
                GTree[father].inverted_list_u[uKey[i]].insert(current);
                term_UNodeMap[uKey[i]].insert(father);
                //??????usr_cnt
                GTree[father].term_usr_Map[uKey[i]].insert(User_id);
                GTree[father].term_usrCnt[uKey[i]] = GTree[father].term_usr_Map[uKey[i]].size();
                //??????term relevant usr_leaf
                set<int> _set = GTree[father].term_usrLeaf_Map[uKey[i]];
                if(_set.size()>0){
                    _set.insert(leafNode);
                    GTree[father].term_usrLeaf_Map[uKey[i]]= _set;
                } else{
                    _set.insert(leafNode);
                    GTree[father].term_usrLeaf_Map[uKey[i]]= _set;
                }

                current = father;
            }

        }
        GTree[leafNode].userSet.push_back(User_id);
        // pair userId, leafNode
        userLeafMap[User_id] = leafNode;
        // build user map
        userMap[User_id] = uKey;

        //?????????????????????????????????????????????
        vector<int> friends = friendshipMap[User_id];
        Users[User_id].friends = friends;

        for(int f:friends){
            GTree[leafNode].unionFriends.insert(f);
            //??????leaf node????????????????????????????????????
            int current = leafNode; int father;
            while(current!=0) {
                father = GTree[current].father;
                GTree[father].unionFriends.insert(f);
                current = father;
            }
        }


    }
    finUsr.close();
    cout << " DONE! USER #:" << usrCnt;
    cout << " , USER KEY #: " << usrKeyCnt << endl;
    LOAD_END
    LOAD_PRINT(" load User ")



    LOAD_START
    // load PoiMap
    printf("LOADING POIMAP...");
    ifstream finPoi;
    string poiPath;
#ifdef LV
    poiPath = getObjectInputPath(BUSINESS);
#else
    poiPath = getObjectNewInputPath(BUSINESS);
#endif

    //finPoi.open(FILE_POI);
    POIs.clear();
    finPoi.open(poiPath.c_str());
    while (getline(finPoi, str)) {
        istringstream tt(str);
        int POI_id, POI_Ni, POI_Nj;
        float POI_dist, POI_dis;
        tt >> POI_id >> POI_Ni >> POI_Nj >> POI_dist >> POI_dis;
        float exchangeNum = (int) (POI_dist * WEIGHT_INFLATE_FACTOR);
        POI_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (POI_dis * WEIGHT_INFLATE_FACTOR);
        POI_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;

        //????????????leaf???poi?????????????????????????????????
        // ?????????????????????
        int leafNode;
        leafNode = Nodes[POI_Ni].gtreepath.back();// User_Ni belong to the leaf : leafNode

        if(GTree[leafNode].minDis == maximum)
            GTree[leafNode].minDis = POI_dist;
        else if(GTree[leafNode].minDis > POI_dist)
            GTree[leafNode].minDis = POI_dist;     //???????????????????????????????????????????????????
        // ??????????????????
        int leafNode2 = Nodes[POI_Nj].gtreepath.back();// User_Ni belong to the leaf : leafNode
        if(GTree[leafNode2].minDis == maximum)
            GTree[leafNode2].minDis = POI_dist;
        else if(GTree[leafNode2].minDis > POI_dist)
            GTree[leafNode2].minDis = POI_dist;

        vector<int> uKey;
        set<int> keywordSet;
        int keyTmp;
        while (tt >> keyTmp) {
            if (-1 == keyTmp) {
                break;
            }

            uKey.push_back(keyTmp);
            keywordSet.insert(keyTmp);
            //??????poi????????????

        }
        for(int term_id:keywordSet)
            invListOfPOI[term_id].push_back(POI_id);

        POI tmpPoi;
        tmpPoi.id = POI_id;
        tmpPoi.Ni = min(POI_Ni,POI_Nj);
        tmpPoi.Nj = max(POI_Ni,POI_Nj);
        tmpPoi.dist = POI_dist;
        tmpPoi.dis = POI_dis;
        tmpPoi.category = 1;
        tmpPoi.keywords = uKey;
        tmpPoi.keywordSet = keywordSet;
        if(poiCheckInIDList[POI_id].size()>0){
            tmpPoi.check_ins = poiCheckInIDList[POI_id];

        }

        EdgeMap[getKey(POI_Ni, POI_Nj)].pts.push_back(tmpPoi);
        POIs.push_back(tmpPoi);
        leafNode = Nodes[tmpPoi.Ni].gtreepath.back();

        //omp_set_num_threads(4);
//#pragma omp parallel for
        for (size_t i = 0; i < uKey.size(); i++){
            //???????????????????????????
            EdgeMap[getKey(POI_Ni, POI_Nj)].OUnionKeys.insert(uKey[i]);
            //??????GIM-tree???????????????????????????
            GTree[leafNode].inverted_list_o[uKey[i]].insert(POI_id);
            term_ONodeMap[uKey[i]].insert(leafNode);
            //??????term_object_map?????????
            GTree[leafNode].term_object_Map[uKey[i]].insert(POI_id);

            GTree[leafNode].term_poiCnt[uKey[i]] =
                    GTree[leafNode].term_object_Map[uKey[i]].size();

            int _tmp = GTree[leafNode].term_poiCnt[uKey[i]];


            //???????????????????????????
            GTree[leafNode].objectUKeySet.insert(uKey[i]);


            leafPoiInv[uKey[i]][leafNode].push_back(POI_id); // ?????? ?????????POI?????????????????????
            poiTerm_leafSet[uKey[i]].insert(leafNode);  // ?????? (poi)??????????????????????????????
            int current = leafNode; int father;
            while(current!=0){
                father =  GTree[current].father;
                GTree[father].objectUKeySet.insert(uKey[i]);
                //??????????????????????????????????????????
                if(true){
                    poiNodeTermChild[father][uKey[i]].insert(current); // father????????? uKey[i] ??????????????????????????????
                    //??????GIM-tree???????????????????????????
                    GTree[father].inverted_list_o[uKey[i]].insert(current);
                    term_ONodeMap[uKey[i]].insert(father);
                    //??????GIM-tree?????????term_object???cnt??????
                    GTree[father].term_object_Map[uKey[i]].insert(POI_id);
                    GTree[father].term_poiCnt[uKey[i]] =
                            GTree[father].term_object_Map[uKey[i]].size();


                }


                current = father;
            }

        }

        GTree[leafNode].poiSet.push_back(tmpPoi.id);
        poiCnt++;
    }
    finPoi.close();
    cout << " DONE! POI #: " << poiCnt << endl;
    LOAD_END
    LOAD_PRINT(" load POI ")

    return 0;
}


int loadOnlyPMap() {

    LOAD_START
    // load PoiMap
    printf("LOADING POIMAP...");
    ifstream finPoi;
    string poiPath;
#ifdef LV
    poiPath = getObjectInputPath(BUSINESS);
#else
    poiPath = getObjectNewInputPath(BUSINESS);
#endif

    string str;
    finPoi.open(poiPath.c_str());
    while (getline(finPoi, str)) {
        istringstream tt(str);
        int POI_id, POI_Ni, POI_Nj;
        float POI_dist, POI_dis;
        tt >> POI_id >> POI_Ni >> POI_Nj >> POI_dist >> POI_dis;
        float exchangeNum = (int) (POI_dist * WEIGHT_INFLATE_FACTOR);
        POI_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (POI_dis * WEIGHT_INFLATE_FACTOR);
        POI_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;

        //????????????leaf???poi?????????????????????????????????
        // ?????????????????????
        int leafNode;
        leafNode = Nodes[POI_Ni].gtreepath.back();// User_Ni belong to the leaf : leafNode

        if(GTree[leafNode].minDis == maximum)
            GTree[leafNode].minDis = POI_dist;
        else if(GTree[leafNode].minDis > POI_dist)
            GTree[leafNode].minDis = POI_dist;     //???????????????????????????????????????????????????


        vector<int> uKey;
        set<int> keywordSet;
        int keyTmp;
        while (tt >> keyTmp) {
            if (-1 == keyTmp) {
                break;
            }

            uKey.push_back(keyTmp);
            keywordSet.insert(keyTmp);
            //??????poi????????????

        }
        for(int term_id:keywordSet){
            invListOfPOI[term_id].push_back(POI_id);
            GTree[leafNode].inverted_list_o[term_id].insert(POI_id);
        }


        POI tmpPoi;
        tmpPoi.id = POI_id;
        tmpPoi.Ni = min(POI_Ni,POI_Nj);
        tmpPoi.Nj = max(POI_Ni,POI_Nj);
        tmpPoi.dist = POI_dist;
        tmpPoi.dis = POI_dis;
        tmpPoi.category = 1;
        tmpPoi.keywords = uKey;
        tmpPoi.keywordSet = keywordSet;
        if(poiCheckInIDList[POI_id].size()>0){
            tmpPoi.check_ins = poiCheckInIDList[POI_id];

        }

        EdgeMap[getKey(POI_Ni, POI_Nj)].pts.push_back(tmpPoi);
        POIs.push_back(tmpPoi);
        leafNode = Nodes[tmpPoi.Ni].gtreepath.back();

        //omp_set_num_threads(4);
//#pragma omp parallel for
        poiCnt++;
    }
    finPoi.close();
    cout << " DONE! POI #: " << poiCnt << endl;
    LOAD_END
    LOAD_PRINT(" load POI ")

    return 0;
}

int loadOnlyPMap(float o_ratio) {

    LOAD_START
    // load PoiMap
    printf("LOADING POIMAP...");
    ifstream finPoi;
    string poiPath;
#ifdef LV
    poiPath = getObjectInputPath(BUSINESS);
#else
    poiPath = getObjectVaryingInputPath(BUSINESS,o_ratio);
#endif

    string str;
    finPoi.open(poiPath.c_str());
    while (getline(finPoi, str)) {
        istringstream tt(str);
        int POI_id, POI_Ni, POI_Nj;
        float POI_dist, POI_dis;
        tt >> POI_id >> POI_Ni >> POI_Nj >> POI_dist >> POI_dis;
        float exchangeNum = (int) (POI_dist * WEIGHT_INFLATE_FACTOR);
        POI_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (POI_dis * WEIGHT_INFLATE_FACTOR);
        POI_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;

        //????????????leaf???poi?????????????????????????????????
        // ?????????????????????
        int leafNode;
        leafNode = Nodes[POI_Ni].gtreepath.back();// User_Ni belong to the leaf : leafNode

        if(GTree[leafNode].minDis == maximum)
            GTree[leafNode].minDis = POI_dist;
        else if(GTree[leafNode].minDis > POI_dist)
            GTree[leafNode].minDis = POI_dist;     //???????????????????????????????????????????????????


        vector<int> uKey;
        set<int> keywordSet;
        int keyTmp;
        while (tt >> keyTmp) {
            if (-1 == keyTmp) {
                break;
            }

            uKey.push_back(keyTmp);
            keywordSet.insert(keyTmp);
            //??????poi????????????

        }
        for(int term_id:keywordSet){
            invListOfPOI[term_id].push_back(POI_id);
            GTree[leafNode].inverted_list_o[term_id].insert(POI_id);
        }


        POI tmpPoi;
        tmpPoi.id = POI_id;
        tmpPoi.Ni = min(POI_Ni,POI_Nj);
        tmpPoi.Nj = max(POI_Ni,POI_Nj);
        tmpPoi.dist = POI_dist;
        tmpPoi.dis = POI_dis;
        tmpPoi.category = 1;
        tmpPoi.keywords = uKey;
        tmpPoi.keywordSet = keywordSet;
        if(poiCheckInIDList[POI_id].size()>0){
            tmpPoi.check_ins = poiCheckInIDList[POI_id];

        }

        EdgeMap[getKey(POI_Ni, POI_Nj)].pts.push_back(tmpPoi);
        POIs.push_back(tmpPoi);
        leafNode = Nodes[tmpPoi.Ni].gtreepath.back();

        //omp_set_num_threads(4);
//#pragma omp parallel for
        poiCnt++;
    }
    finPoi.close();
    cout << " DONE! POI #: " << poiCnt << endl;
    LOAD_END
    LOAD_PRINT(" load POI ")

    return 0;
}



int modifyUPMap() {

    // ?????? edge information
    //ConfigType setting;
    //AddConfigFromFile(setting,"../map_exp");
    // ?????? Gtree??????
    loadEdgeMap_light();
    stringstream ss;
    ofstream  output;  string str;

    LOAD_START
    // read UserMap and insert User to edges
    printf("LOADING USERMAP...");

    int usrKeyCnt = 0;

    //user???????????????
    ss.str("");
    ss<<DATASET_HOME<<dataset_name<<"/source_dataset/""user_poi/"<<dataset_name <<".user";
    std::string user_path = ss.str();

    //??????????????????user????????????
    ifstream finUsr;
    ss.str("");
    ss<<DATASET_HOME<<dataset_name<<"/objects/raw/"<<dataset_name <<".user";
    std::string user_path_new = ss.str();

    finUsr.open(user_path.c_str());
    output.open(user_path_new.c_str());

    //??????User??????
    while (getline(finUsr, str)) {
        int User_id, User_Ni, User_Nj;
        float User_dist, User_dis;
        istringstream tt(str);
        tt >> User_id >> User_Ni >> User_Nj >> User_dist >> User_dis;

        //??????user?????????????????????????????? edge_dist
        edge e = EdgeMap[getKey(User_Ni,User_Nj)];
        float e_dist = e.dist;

        if(User_Ni >=VertexNum || User_Nj >=VertexNum){
            cout<<"?????????Edge("<<User_Ni<<","<<User_Nj<<")!"<<endl;
            User_Ni = User_Ni % VertexNum;
            vector<int> adj_vertex = adjList[User_Ni];
            int adj_size = adj_vertex.size();
            int idx = random() % adj_size;
            User_Nj = adj_vertex[idx];
            User_dist = EdgeMap[getKey(User_Ni,User_Nj)].dist;
            cout<<"??????edge("<<User_Ni<<","<<User_Nj<<")="<<User_dist<<endl;
        }
        else{
            User_dist = e_dist;
            User_dis = min(e_dist,User_dis);
        }

        //??????????????????????????????
        vector<int> uKey;
        int keyTmp;

        while (tt >> keyTmp) {
            if (-1 == keyTmp) {
                break;
            }
            usrKeyCnt++;
            uKey.push_back(keyTmp);

        }

        usrCnt++;
        output << User_id <<' '<< User_Ni<<' '<< User_Nj<<' '<< User_dist<<' '<< User_dis;
        for(int term: uKey)
            output <<' '<< term;
        output<<endl;

    }
    finUsr.close();
    output.close();

    cout << " DONE! Modify USER #:" << usrCnt;
    cout << " , USER KEY #: " << usrKeyCnt << endl;
    LOAD_END
    LOAD_PRINT(" load User ")


    //??????poi??????
    LOAD_START
    // load PoiMap
    printf("LOADING POIMAP...");
    ifstream finPoi;

    //poi???????????????
    ss.str("");
    ss<<DATASET_HOME<<dataset_name<<"/source_dataset/""user_poi/"<<dataset_name <<".poi";
    std::string poi_path = ss.str();

    //??????????????????poi????????????
    ss.str("");
    ss<<DATASET_HOME<<dataset_name<<"/objects/raw/"<<dataset_name <<".poi";
    std::string poi_path_new = ss.str();

    finPoi.open(poi_path.c_str());
    output.open(poi_path_new.c_str());


    while (getline(finPoi, str)) {
        istringstream tt(str);
        int POI_id, POI_Ni, POI_Nj;
        float POI_dist, POI_dis;
        tt >> POI_id >> POI_Ni >> POI_Nj >> POI_dist >> POI_dis;


        edge e = EdgeMap[getKey(POI_Ni,POI_Nj)];
        float e_dist = e.dist;

        if(POI_Ni >=VertexNum ||POI_Nj >=VertexNum){
            POI_Ni = POI_Ni % VertexNum;
            vector<int> adj_vertex = adjList[POI_Ni];
            int adj_size = adj_vertex.size();
            int idx = random() % adj_size;
            POI_Nj = adj_vertex[idx];
            POI_dist = EdgeMap[getKey(POI_Ni,POI_Nj)].dist;
        }
        else{
            POI_dist = e_dist;
            POI_dis = min(e_dist,POI_dis);
        }

        set<int> keywordSet;
        int keyTmp;
        while (tt >> keyTmp) {
            if (-1 == keyTmp) {
                break;
            }
            keywordSet.insert(keyTmp);
            //??????poi????????????

        }

        output << POI_id <<' '<< POI_Ni<<' '<< POI_Nj<<' '<< POI_dist<<' '<< POI_dis;

        for(int term_id:keywordSet){
            output <<' '<< term_id;
        }
        output<<endl;

        poiCnt++;
    }
    finPoi.close();
    output.close();
    cout << " DONE! Modify POI #: " << poiCnt << endl;
    LOAD_END
    LOAD_PRINT(" load POI ")



    return 0;
}

int modifyLasVegasUserPOI() {

    // ?????? edge information
    //ConfigType setting;
    //AddConfigFromFile(setting,"../map_exp");
    initGtree();
    loadEdgeMap();
    stringstream ss;
    ofstream  output;  string str;

    LOAD_START
    // load PoiMap
    printf("LOADING POIMAP...");
    ifstream finPoi;

    //poi???????????????
    ss.str("");
    ss<<DATASET_HOME<<dataset_name<<"/source_dataset/""user_poi/"<<dataset_name <<".poi";
    std::string poi_path = ss.str();

    //??????????????????poi????????????
    ss.str("");
    ss<<DATASET_HOME<<dataset_name<<"/objects/"<<dataset_name <<".poi";
    std::string poi_path_new = ss.str();

    finPoi.open(poi_path.c_str());
    output.open(poi_path_new.c_str());


    while (getline(finPoi, str)) {
        istringstream tt(str);
        int POI_id, POI_Ni, POI_Nj;
        float POI_dist, POI_dis;
        tt >> POI_id >> POI_Ni >> POI_Nj >> POI_dist >> POI_dis;

        edge e = EdgeMap[getKey(POI_Ni,POI_Nj)];
        float e_dist = e.dist;
        POI_dist = e_dist;
        POI_dis = min(e_dist,POI_dis);

        set<int> keywordSet;
        int keyTmp;
        while (tt >> keyTmp) {
            if (-1 == keyTmp) {
                break;
            }
            keywordSet.insert(keyTmp);
            //??????poi????????????

        }

        output << POI_id <<' '<< POI_Ni<<' '<< POI_Nj<<' '<< POI_dist<<' '<< POI_dis;

        for(int term_id:keywordSet){
            output <<' '<< term_id;
        }
        output<<endl;

        poiCnt++;
    }
    finPoi.close();
    cout << " DONE! Modify POI #: " << poiCnt << endl;
    LOAD_END
    LOAD_PRINT(" load POI ")



    LOAD_START
    // read UserMap and insert User to edges
    printf("LOADING USERMAP...");

    int usrKeyCnt = 0;

    //user???????????????

    ss<<DATASET_HOME<<dataset_name<<"/source_dataset/""user_poi/"<<dataset_name <<".user";
    std::string user_path = ss.str();

    //??????????????????user????????????
    ifstream finUsr;
    ss.str("");
    ss<<DATASET_HOME<<dataset_name<<"/objects/"<<dataset_name <<".user";
    std::string user_path_new = ss.str();

    finUsr.open(user_path.c_str());
    output.open(user_path_new.c_str());

    //??????User??????
    while (getline(finUsr, str)) {
        int User_id, User_Ni, User_Nj;
        float User_dist, User_dis;
        istringstream tt(str);
        tt >> User_id >> User_Ni >> User_Nj >> User_dist >> User_dis;
        //??????user?????????????????????????????? edge_dist
        edge e = EdgeMap[getKey(User_Ni,User_Nj)];
        float e_dist = e.dist;
        User_dist = e_dist;
        User_dis = min(e_dist,User_dis);

        //??????????????????????????????
        vector<int> uKey;
        int keyTmp;

        while (tt >> keyTmp) {
            if (-1 == keyTmp) {
                break;
            }
            usrKeyCnt++;
            uKey.push_back(keyTmp);

        }

        usrCnt++;
        output << User_id <<' '<< User_Ni<<' '<< User_Nj<<' '<< User_dist<<' '<< User_dis;
        for(int term: uKey)
            output <<' '<< term;
        output<<endl;

    }
    finUsr.close();
    output.close();
    cout << " DONE! Modify USER #:" << usrCnt;
    cout << " , USER KEY #: " << usrKeyCnt << endl;
    LOAD_END
    LOAD_PRINT(" load User ")





    return 0;
}


int initialO2UData(){
    for(int i=0;i<UserID_MaxKey;i++){
        User u;
        Users.push_back(u);
    }
    int _sizeu = Users.size();
    cout<<"initial Users list cache space!"<<endl;
    for(int j=0;j<poi_num;j++){
        POI p;
        POIs.push_back(p);
    }
    int _size0 = poi_num;
    int _size = POIs.size();
    cout<<"initial POIs list cache space!"<<endl;
}


int initialUserData(){
    for(int i=0;i<UserID_MaxKey;i++){
        User u;
        Users.push_back(u);
    }
    cout<<"initial Users list cache space!"<<endl;
}


/*---------------------------------------1206-------------------------------*/

//?????????????????????????????????????????????????????????????????????

//?????????????????????
void makePtFiles(char *fileprefix) {
    //???????????????????????????????????????????????????????????????.p_d???????????????.p_bt???
    bool PtMaxUsing = true;
    //????????????????????????????????????
    char poiDataFile_Name[255]; char poiIndxFile_Name[255];
    sprintf(poiDataFile_Name,"%s.poi_d",fileprefix);
    //sprintf(poiIndxFile_Name,"%s.poi_bt",fileprefix);
    remove(poiDataFile_Name); // remove existing file ?????????????????????
    //remove(poiIndxFile_Name); // remove existing file ?????????????????????
    //??????????????????
    FILE* ptFile=fopen(poiDataFile_Name,"w+");
    if(ptFile==NULL) cout<<"faile to open";
    //??????????????????
    //BTree* bt=initialize(poiIndxFile_Name);             //btree initialize
    printf("making PtFiles\n");
    int RawAddr=0,size;	// changed
    EdgeMapType::iterator iter=EdgeMap.begin();
    //?????????????????????poi????????????
    edge tmpedg;
    int key=0;
    int Ni =0; int Nk = adjList[Ni][0];
    while (iter!=EdgeMap.end()) {
        int _mapKey = iter->first;
        tmpedg = iter->second;
        edge* e= &(iter->second);
        int poi_size = 0;
        poi_size = tmpedg.pts.size();
        if (poi_size > 0) {	// do not index empty groups
            //sort(e->pts.begin(),e->pts.end());

            RawAddr=ftell(ptFile);	// set addr to correct amt.
            int poi_size=e->pts.size();
            int key_0 = key;
            e->poiAddress = key;
            //cout<<"_mapKey="<<_mapKey<<",poiAddress="<<e->poiAddress<<endl;


            //?????????????????????????????????????????????????????????
            fwrite(&(e->Ni),1,sizeof(int),ptFile);      //1.??????Ni, NI_P
            //cout<<"??????Ni="<<e->Ni<<",key="<<key<<endl;
            fwrite(&(e->Nj),1,sizeof(int),ptFile);      //2.??????Nj, NJ_P
            fwrite(&(e->dist),1,sizeof(float),ptFile);  //3.??????e.dist, DIST_P
            fwrite(&(poi_size),1,sizeof(int),ptFile);       //4.??????e.pois.size, SIZE_P
            //????????????(???????????????????????????????????????key)
            key+=sizeof(int)*3+sizeof(float);
            //???????????????????????????????????????????????????????????????????????????{<poi_id, dis, keyword{num1....}, checkin_usr{num2...}>}
            int poi_key = key; int poi_idx[poi_size];
            for(int i=0;i<poi_size; i++){
                POI p= (e->pts)[i];
                poi_idx[i] = key; //???????????????????????????????????????address
                //cout<<"count poi_idx["<<i<<"]="<<poi_idx[i]<<endl;

                fwrite(&(p.id),1,sizeof(int),ptFile);   //??????poi id
                key += sizeof(int);
                fwrite(&(p.dis),1,sizeof(float),ptFile);  //??????poi???????????????Ni?????????
                key += sizeof(float); //1???????????????dis??? (?????????)
                //?????????????????????
                int num1 = p.keywords.size();
                fwrite(&(num1),1,sizeof(int),ptFile);  //?????????????????????
                for(int term: p.keywords){
                    fwrite(&(term),1,sizeof(int),ptFile); //??????????????????
                }
                key+=sizeof(int)*num1+sizeof(int); //??????1????????????num1?????????num1????????????????????????
                //????????????????????????
                int num2 =0;
                if(p.check_ins.size()>0){
                    num2 = p.check_ins.size();
                    //cout<<"p.id="<<p.id<<"check size (num2)="<<num2<<endl;
                    fwrite(&(num2),1,sizeof(int),ptFile);  //????????????????????????
                    for(int checkin_usr: p.check_ins){
                        //if(p.id==12156) cout<<checkin_usr<<",";
                        fwrite(&(checkin_usr),1,sizeof(int),ptFile); //?????????????????????
                    }

                }
                else{

                    fwrite(&(num2),1,sizeof(int),ptFile);  //????????????????????????
                     //?????????????????????????????????num1????????????*/
                }
                key+=sizeof(int)*num2+sizeof(int); //?????????????????????????????????num1????????????*/

            }
            e->poiIdxAddress = key;   //?????????????????????????????????????????????key??????????????????????????????????????? (??????)
            //cout<<e->poiIdxAddress<<","<<getFileSize(ptFile)<<endl;
            for(int i=0;i<poi_size; i++){
                int addr = poi_idx[i];
                //cout<<"write poi_idx["<<i<<"]="<<addr<<",file_size="<<getFileSize(ptFile)<<endl;
                fwrite(&(addr),1,sizeof(int),ptFile);   //??????poi id
                key += sizeof(int);

            }
            //getchar();
            PtMaxKey=key+e->pts.size()-1;	// useful for our special ordering !
            //printf("(key,value)=(%d,%d)\n",key,RawAddr);

            //key+=e->pts.size()*sizeof(InvNode); // changed, moved to after addentry
        } else if(poi_size==0){
            e->poiAddress =-1;		// also later used by AdjFile
            e->poiIdxAddress =-1;
        }
        iter++;
    }
    //finalize(bt);
    //bt->UserField = num_D;
    //delete bt;
    fclose(ptFile);
    PtMaxUsing=false;
    printf("making PtFiles success!\n");
    //getchar();

}




// ?????????????????????????????????????????????????????????????????????????????????makePtFiles?????????,????????????????????????????????????poiAddress???
void makeAdjListFiles_withoutKeyword(const char* outfileprefix) {

    //printf("?????????????????????????????????????????????\n");
    char idxFileName[200];
    sprintf(idxFileName,"%s.adj",outfileprefix);	remove(idxFileName); // remove existing file
    FILE * alFile=fopen(idxFileName,"w+");

    int key=0,size,PtSize;
    //ofstream fout;                  //????????????????????
    //fout.open("??/makeAdj.txt");
    //????????????????????????
    int num = VertexNum;
    fwrite(&num,1,sizeof(int),alFile);

    // slotted header info.
    int addr=sizeof(int)+sizeof(int)*VertexNum; //????????????????????????????????????????????????VertexNum?????????????????????
    //????????????????????????????????????????????????????????????????????????????????????
    for (int Ni=0;Ni < VertexNum;Ni++) {
        fwrite(&addr,1,sizeof(int),alFile);  //???????????????
        int addr_delta = sizeof(int)+adjList[Ni].size()*(3*sizeof(int)+sizeof(float));  //4+n*(12+4)
        addr += addr_delta;  //?????????*(3*sizeof(int)+sizeof(float))  //???????????????????????????????????????1.?????????????????????????????? size??????????????????Nj_ID, Dist, poigrp_address, poiIdxList_address??????
        //cout<<"addr="<<addr<<",delta="<<addr_delta<<endl;
        //}
    }
    //???????????????????????????????????????????????????
    float distsum=0;
    for (int Ni=0;Ni<VertexNum;Ni++) {
        size=adjList[Ni].size();
        fwrite(&(size),1,sizeof(int),alFile);               //----???????????????????????????
        for (int k=0;k < adjList[Ni].size();k++) {
            int Nk = adjList[Ni][k];	// Nk can be smaller or greater than Ni !!!
            edge e = EdgeMap[getKey(Ni,Nk)];
            PtSize=e.pts.size();
            fwrite(&Nk,1,sizeof(int),alFile);               //----???????????????ID
            fwrite(&(e.dist),1,sizeof(float),alFile);      //----????????????????????????

            fwrite(&e.poiAddress,1,sizeof(int),alFile);                //----????????????????????????POI group ??????
            //cout<<"e.poiAddress="<<e.poiAddress;
            fwrite(&e.poiIdxAddress,1,sizeof(int),alFile);
            //cout<<",  e.poiIdxAddress="<<e.poiIdxAddress<<endl;
            //fwrite(&e.poiIdxAddress,1,sizeof(int),alFile);
            //fwrite(&(e->FirstRow),1,sizeof(int),alFile); // use FirstRow for other purpose ...//----????????????????????????POI group ??????
            //printf("(Ni,Nj,dataAddr)=(%d,%d,%d)\n",Ni,Nk,e->FirstRow);
            distsum+=e.dist;
            //fout<<"N"<<Nk<<"\t\tdist="<<e->dist<<"\t\tFirstRow="<<e->FirstRow<<endl;  //????????????????????
        }
        key=Ni;
    }
    distsum=distsum/2;
    fclose(alFile);
    printf("making AdjList success!\n");
    //getchar();
    //fout.close();
}


//???????????????
void makeAdjListFiles2(const char* outfileprefix) {

    //printf("?????????????????????????????????????????????\n");
    char idxFileName[200];
    sprintf(idxFileName,"%s.adj",outfileprefix);	remove(idxFileName); // remove existing file
    FILE * alFile=fopen(idxFileName,"w+");
    int current_point = 0; int addr =0;

    int key=0,size,PtSize;
    //????????????????????????
    int num = VertexNum;
    fwrite(&num,1,sizeof(int),alFile);
    current_point += sizeof(int);
    addr += sizeof(int);
    addr= addr +sizeof(int)*VertexNum; //????????????????????????????????????????????????VertexNum?????????????????????

     //???????????????????????????????????????????????????????????????????????????
    for (int Ni=0;Ni < VertexNum;Ni++) {
        fwrite(&addr,1,sizeof(int),alFile);  //???Ni????????????????????????
        current_point += sizeof(int);
        //??????edge size
        int _size = adjList[Ni].size();
        addr += sizeof(int);
        for (int j=0;j< adjList[Ni].size();j++) {
            int Nk = adjList[Ni][j];
            if(Ni==96697&& Nk==96712){
                cout<<"??????????????????????????????find 96697&& Nk==96712"<<endl;
                cout<<"e("<<Ni<<","<<Nk<<")"<<endl;
                cout<<"addr="<<addr<<endl;
                getchar();
            }

            //??????edge data space
            int addr_delta = (4*sizeof(int)+sizeof(float));  // 4???int????????????Nj, POIDataAddr, POIIdxAddr, OKeyAddr ??? float:dist
            addr += addr_delta;  //?????????*(3*sizeof(int)+sizeof(float))  //???????????????????????????????????????1.?????????????????????????????? size??????????????????Nj_ID, Dist, poigrp_address, poiIdxList_address??????

        }

    }

    cout<<"Okey??????????????????addr="<<addr<<endl;

    //?????? e.OUnionKeys ???????????????
    for (int Ni=0;Ni < VertexNum;Ni++) {
        for (int j=0;j< adjList[Ni].size();j++) {
            int Nk = adjList[Ni][j];
            edge& e = EdgeMap[getKey(Ni,Nk)];
            if(Ni==96697&& Nk==96712){
                cout<<"???????????????????????????????????????find 96697&& Nk==96712"<<endl;
                cout<<"e.okeyAddress="<<addr<<endl;
                cout<<"_okeyword_onEdge_size="<<e.OUnionKeys.size()<<endl;
                getchar();
            }


            int okeyword_onEdge_size = e.OUnionKeys.size();
            if(okeyword_onEdge_size>0 && e.okeyAddress ==-1){  // careful: ????????????????????? (Ni, Nk) ???  (Nk, Ni)???????????????????????????????????????e.okeyAddress?????????????????????
                //?????? edge text psudo doc?????????
                e.okeyAddress = addr;
                //??????????????????
                int _okeyword_onEdge_size = e.OUnionKeys.size();
                //???????????? keyword size ??????0???????????????????????????????????????
                if(_okeyword_onEdge_size>0){
                    addr += sizeof(int);
                    for(int term: e.OUnionKeys){
                        int _term = term;
                        addr += sizeof(int);
                    }
                }


            }

        }
    }

    //???????????????????????????????????????????????????
    float distsum=0;
    for (int Ni=0;Ni<VertexNum;Ni++) {
        size=adjList[Ni].size();
        fwrite(&(size),1,sizeof(int),alFile);               //----???????????????????????????
        current_point += sizeof(int);
        for (int j=0;j < adjList[Ni].size();j++) {
            int Nk = adjList[Ni][j];	// Nk can be smaller or greater than Ni !!!
            edge e = EdgeMap[getKey(Ni,Nk)];

            if(Nk < Ni) continue;  //????????????????????????????????????
            if(Ni==96697&& Nk==96712){
                cout<<"write ???????????????Ni==96697&& Nk==96712"<<endl;
                cout<<"write e.okeyAddress="<<e.okeyAddress<<endl;
                cout<<"current_point="<<current_point<<endl;
            }

            PtSize=e.pts.size();

            fwrite(&Nk,1,sizeof(int),alFile);               //----???????????????ID
            fwrite(&(e.dist),1,sizeof(float),alFile);      //----????????????????????????

            fwrite(&e.poiAddress,1,sizeof(int),alFile);                //----????????????????????????POI data group?????????
            //cout<<"e.poiAddress="<<e.poiAddress;
            fwrite(&e.poiIdxAddress,1,sizeof(int),alFile);   //----??????????????????POI ??????????????????


            int _address = e.okeyAddress;
            fwrite(&e.okeyAddress,1,sizeof(int),alFile);  //----??????????????? text psudo doc?????????
            current_point = current_point + (4*sizeof(int)+sizeof(float));


        }
    }
    cout<<"?????????????????????current address="<<current_point<<endl;
    //getchar();
    //????????????????????????????????????????????????????????????????????????
    for (int Ni=0;Ni<VertexNum;Ni++) {
        for (int j=0;j < adjList[Ni].size();j++) {

            int Nk = adjList[Ni][j];	// Nk can be smaller or greater than Ni !!!
            // ????????????????????? (Ni, Nk) ???  (Nk, Ni)?????????????????????????????? j>i???????????????????????????
            if(Nk < Ni) continue; //careful????????????careful?????????


            edge& e = EdgeMap[getKey(Ni,Nk)];

            //??????????????????poi keyword??? vocabulary??????
            set<int> okeyword_onEdge = e.OUnionKeys;
            int _size = okeyword_onEdge.size();
            if(_size>0){   //careful????????????careful?????????: ??????????????????????????????????????????????????????
                //???????????????oKeySize_Edge
                if(Ni==96697&& Nk==96712){
                    cout<<"write ????????????????????????Ni==96697 && Nk==96712"<<endl;
                    cout<<"_size="<<_size<<endl;
                    cout<<"??????????????????"; printSetElements(okeyword_onEdge);
                    cout<<"current_point="<<current_point<<endl;
                    getchar();
                }


                fwrite(&_size, 1, sizeof(int), alFile);
                current_point += sizeof(int);

                //??????edge ?????????keywords
                for(int term: okeyword_onEdge){
                    fwrite(&term, 1, sizeof(int), alFile);
                    current_point += sizeof(int);
                }

            }

        }
    }


    distsum=distsum/2;
    fclose(alFile);
    printf("making AdjList success!\n");
    //getchar();
    //fout.close();
}


void makeAdjListFiles(const char* outfileprefix) {

    //printf("?????????????????????????????????????????????\n");
    char idxFileName[200];
    sprintf(idxFileName,"%s.adj",outfileprefix);	remove(idxFileName); // remove existing file
    FILE * alFile=fopen(idxFileName,"w+");
    int current_point = 0; int addr =0;

    int key=0,size,PtSize;
    //????????????????????????
    int num = VertexNum;
    fwrite(&num,1,sizeof(int),alFile);
    current_point += sizeof(int);
    addr += sizeof(int);
    addr= addr +sizeof(int)*VertexNum; //????????????????????????????????????????????????VertexNum?????????????????????

    //???????????????????????????????????????????????????????????????????????????
    for (int Ni=0;Ni < VertexNum;Ni++) {
        fwrite(&addr,1,sizeof(int),alFile);  //???Ni????????????????????????
        current_point += sizeof(int);
        //??????edge size
        int _size = adjList[Ni].size();
        addr += sizeof(int);
        for (int j=0;j< adjList[Ni].size();j++) {
            int Nk = adjList[Ni][j];
            if(Ni==96697&& Nk==96712){
                cout<<"??????????????????????????????find 96697&& Nk==96712"<<endl;
                cout<<"e("<<Ni<<","<<Nk<<")"<<endl;
                cout<<"addr="<<addr<<endl;
                //getchar();
            }

            //??????edge data space
            int addr_delta = (4*sizeof(int)+sizeof(float));  // 4???int????????????Nj, POIDataAddr, POIIdxAddr, OKeyAddr ??? float:dist
            addr += addr_delta;  //?????????*(3*sizeof(int)+sizeof(float))  //???????????????????????????????????????1.?????????????????????????????? size??????????????????Nj_ID, Dist, poigrp_address, poiIdxList_address??????

        }

    }

    cout<<"Okey??????????????????addr="<<addr<<endl<<endl;

    //?????? e.OUnionKeys ???????????????
    for (int Ni=0;Ni < VertexNum;Ni++) {
        for (int j=0;j< adjList[Ni].size();j++) {
            int Nj = adjList[Ni][j];
            edge& e = EdgeMap[getKey(Ni,Nj)];
           /* if(Ni==96697&& Nj==96712){
                cout<<"???????????????????????????????????????find 96697&& Nk==96712"<<endl;
                cout<<"e.okeyAddress="<<addr<<endl;
                cout<<"_okeyword_onEdge_size="<<e.OUnionKeys.size()<<endl;
                //getchar();
            }*/
            if(Ni>Nj) continue;  //??????????????????????????????

            //?????? edge text psudo doc?????????
            e.okeyAddress = addr;
            //???????????????????????????????????????
            int okeyword_onEdge_size = e.OUnionKeys.size();
            addr += sizeof(int);
            //??????????????????????????????????????????
            for(int term: e.OUnionKeys){
                int _term = term;
                addr += sizeof(int);
            }


        }
    }

    //???????????????????????????????????????????????????
    float distsum=0;
    for (int Ni=0;Ni<VertexNum;Ni++) {
        size=adjList[Ni].size();
        fwrite(&(size),1,sizeof(int),alFile);               //----???????????????????????????
        current_point += sizeof(int);
        for (int j=0;j < adjList[Ni].size();j++) {
            int Nj = adjList[Ni][j];	// Nk can be smaller or greater than Ni !!!
            edge e = EdgeMap[getKey(Ni,Nj)];

           /* if(Ni==96697&& Nj==96712){
                cout<<"write ???????????????Ni==96697&& Nk==96712"<<endl;
                cout<<"write e.okeyAddress="<<e.okeyAddress<<endl;
                cout<<"current_point="<<current_point<<endl;
            }*/

            PtSize=e.pts.size();

            fwrite(&Nj,1,sizeof(int),alFile);               //----???????????????ID
            fwrite(&(e.dist),1,sizeof(float),alFile);      //----????????????????????????

            fwrite(&e.poiAddress,1,sizeof(int),alFile);                //----????????????????????????POI data group?????????
            //cout<<"e.poiAddress="<<e.poiAddress;
            fwrite(&e.poiIdxAddress,1,sizeof(int),alFile);   //----??????????????????POI ??????????????????


            int _address = e.okeyAddress;
            fwrite(&e.okeyAddress,1,sizeof(int),alFile);  //----??????????????? text psudo doc?????????
            current_point = current_point + (4*sizeof(int)+sizeof(float));


        }
    }
    cout<<"?????????????????????current address="<<current_point<<endl;
    //getchar();
    //????????????????????????????????????????????????????????????????????????
    for (int Ni=0;Ni<VertexNum;Ni++) {
        for (int j=0;j < adjList[Ni].size();j++) {

            int Nj = adjList[Ni][j];	// Nk can be smaller or greater than Ni !!!
            // ????????????????????? (Ni, Nk) ???  (Nk, Ni)?????????????????????????????? j>i???????????????????????????
            if(Nj < Ni) continue; //careful????????????careful?????????


            edge& e = EdgeMap[getKey(Ni,Nj)];

            //??????????????????poi keyword??? vocabulary??????
            set<int> okeyword_onEdge = e.OUnionKeys;
            int _size = okeyword_onEdge.size();
            //???????????????oKeySize_Edge

            fwrite(&_size, 1, sizeof(int), alFile);
            current_point += sizeof(int);

            //??????edge ?????????keywords
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






//???????????????????????????????????????G-tree???????????????????????????????????????????????????.pl_d???????????????.pl_bt???
void makePtFileByLeaf(const char* fileprefix){
    //???????????????????????????????????????leaf node????????????????????????????????????.p_d_leaf
    bool PtMaxUsing = true;
    //????????????????????????????????????
    char poiDataFile_Name[255]; char poiIndxFile_Name[255];
    sprintf(poiDataFile_Name,"%s.poi_ld",fileprefix);
    sprintf(poiIndxFile_Name,"%s.poi_idx",fileprefix);
    remove(poiDataFile_Name); // remove existing file ?????????????????????
    remove(poiIndxFile_Name); // remove existing file ?????????????????????
    //??????????????????
    FILE* ptFile=fopen(poiDataFile_Name,"w+");
    if(ptFile==NULL) cout<<"faile to open";
    //??????????????????
    //BTree* bt=initialize(poiIndxFile_Name);             //btree initialize
    FILE* idxFile=fopen(poiIndxFile_Name,"w+");
    if(idxFile==NULL) cout<<"faile to open";


    printf("making PtLeafDataFiles, PtLeafIdxFiles\n");
    int RawAddr=0,size;	// changed
    int key=0;
    for(int i=0; i<GTree.size(); i++){
        TreeNode& node = GTree[i];
        if(node.isleaf){
            if(node.poiSet.size()>0){
                //???????????????????????????
                for(int poi_id: node.poiSet){
                    RawAddr = ftell(ptFile);

                    //addentry(bt,&top_level,i_capacity,1,key,&num_written_blocks,RawAddr);

                    POI& p = POIs[poi_id];
                    POIs[poi_id].addressDataOrgByLeaf = key; //??????????????????????????????????????????
                    //?????????????????????
                    fwrite(&(p.id),1,sizeof(int),ptFile);   //??????poi id
                    key += sizeof(int);

                    int _Ni = p.Ni;
                    fwrite(&(_Ni),1,sizeof(int),ptFile);   //??????poi id
                    key += sizeof(int);

                    int _Nj = p.Nj;
                    fwrite(&(_Nj),1,sizeof(int),ptFile);   //??????poi id
                    key += sizeof(int);

                    fwrite(&(p.dist),1,sizeof(float),ptFile);   //??????poi id
                    key += sizeof(float);

                    fwrite(&(p.dis),1,sizeof(float),ptFile);  //??????poi???????????????Ni?????????
                    key += sizeof(float); //1???????????????dis??? (?????????)
                    //?????????????????????
                    int keyNum  = 0;
                    keyNum = p.keywords.size();
                    fwrite(&(keyNum),1,sizeof(int),ptFile);  //?????????????????????
                    for(int term: p.keywords){
                        fwrite(&(term),1,sizeof(int),ptFile); //??????????????????
                    }
                    key+=sizeof(int)*keyNum + sizeof(int); //??????1????????????num1?????????num1????????????????????????
                    //????????????????????????
                    int checkInNum =0;
                    if(p.check_ins.size()>0){
                        checkInNum = p.check_ins.size();
                        //cout<<"p.id="<<p.id<<"check size (num2)="<<num2<<endl;
                        fwrite(&(checkInNum),1,sizeof(int),ptFile);  //????????????????????????
                        for(int checkin_usr: p.check_ins){
                            //if(p.id==12156) cout<<checkin_usr<<",";
                            fwrite(&(checkin_usr),1,sizeof(int),ptFile); //?????????????????????
                        }

                    }
                    else{

                        fwrite(&(checkInNum),1,sizeof(int),ptFile);  //????????????????????????
                        //?????????????????????????????????num1????????????*/
                    }
                    key+=sizeof(int)*checkInNum+sizeof(int); //?????????????????????????????????num1????????????*/

                }
            }


        }
    }
    //???poi???????????????????????????
    int poi_Size = POIs.size();
    fwrite(&poi_Size, 1, sizeof(int),idxFile);
    for(int i=0;i<POIs.size();i++){
        int address = POIs[i].addressDataOrgByLeaf;
        fwrite(&address, 1, sizeof(int),idxFile);
    }


    fclose(ptFile);
    fclose(idxFile);
    printf("making PtFiles and IdxFile (Organized By Leaf) success !\n");



}

void makeO2UFileByLeaf(const char* fileprefix){
    //???????????????????????????????????????leaf node????????????????????????????????????.p_d_leaf
    bool PtMaxUsing = true;
    //????????????????????????????????????
    char o2uDataFile_Name[255]; char o2uIndxFile_Name[255];
    sprintf(o2uDataFile_Name,"%s.o2u_ld",fileprefix);
    sprintf(o2uIndxFile_Name,"%s.o2u_idx",fileprefix);
    remove(o2uDataFile_Name); // remove existing file ?????????????????????
    remove(o2uIndxFile_Name); // remove existing file ?????????????????????
    //??????????????????
    FILE* o2uFile=fopen(o2uDataFile_Name,"w+");
    if(o2uFile==NULL) cout<<"faile to open";
    //??????????????????
    //BTree* bt=initialize(o2uIndxFile_Name);             //btree initialize
    FILE* idxFile=fopen(o2uIndxFile_Name,"w+");
    if(idxFile==NULL) cout<<"faile to open";


    printf("making O2ULeafDataFiles, O2ULeafIdxFiles\n");
    int RawAddr=0,size;	// changed
    int key=0;
    //???leaf??????????????????????????????????????????????????????
    for(int i=0; i<GTree.size(); i++){
        TreeNode& node = GTree[i];
        if(node.isleaf){
            //??????????????? poi ?????????
            if(node.poiSet.size()>0){
                //???????????????????????????
                for(int poi_id: node.poiSet){
                    RawAddr = ftell(o2uFile);

                    //addentry(bt,&top_level,i_capacity,1,key,&num_written_blocks,RawAddr);

                    POI& p = POIs[poi_id];
                    POIs[poi_id].addressDataOrgByLeaf = key; //??????????????????????????????????????????
                    //?????????????????????
                    fwrite(&(p.id),1,sizeof(int),o2uFile);   //??????poi id
                    key += sizeof(int);

                    int _Ni = p.Ni;
                    fwrite(&(_Ni),1,sizeof(int),o2uFile);   //??????poi id
                    key += sizeof(int);

                    int _Nj = p.Nj;
                    fwrite(&(_Nj),1,sizeof(int),o2uFile);   //??????poi id
                    key += sizeof(int);

                    fwrite(&(p.dist),1,sizeof(float),o2uFile);   //??????poi id
                    key += sizeof(float);

                    fwrite(&(p.dis),1,sizeof(float),o2uFile);  //??????poi???????????????Ni?????????
                    key += sizeof(float); //1???????????????dis??? (?????????)
                    //?????????????????????
                    int keyNum  = 0;
                    keyNum = p.keywords.size();
                    fwrite(&(keyNum),1,sizeof(int),o2uFile);  //?????????????????????
                    for(int term: p.keywords){
                        fwrite(&(term),1,sizeof(int),o2uFile); //??????????????????
                    }
                    key+=sizeof(int)*keyNum + sizeof(int); //??????1????????????num1?????????num1????????????????????????
                    //????????????????????????
                    int checkInNum =0;
                    checkInNum = p.check_ins.size();
                    fwrite(&(checkInNum),1,sizeof(int),o2uFile);  //????????????????????????
                    for(int checkin_usr: p.check_ins){
                        //if(p.id==12156) cout<<checkin_usr<<",";
                        fwrite(&(checkin_usr),1,sizeof(int),o2uFile); //?????????????????????
                    }

                    key+=sizeof(int)*checkInNum+sizeof(int); //?????????????????????????????????num1????????????*/

                }
            }
            //??????????????? user ?????????
            if(node.userSet.size()>0){
                //???????????????????????????
                for(int usr_id: node.userSet){
                    RawAddr = ftell(o2uFile);

                    User& u = Users[usr_id];
                    Users[usr_id].addressDataOrgByLeaf = key; //??????????????????????????????????????????
                    //?????????????????????
                    int u_id = u.id;
                    fwrite(&(u_id),1,sizeof(int),o2uFile);   //??????u id
                    key += sizeof(int);

                    int _Ni = u.Ni;
                    fwrite(&(_Ni),1,sizeof(int),o2uFile);   //??????u Ni
                    key += sizeof(int);

                    int _Nj = u.Nj;
                    fwrite(&(_Nj),1,sizeof(int),o2uFile);   //??????u Nj
                    key += sizeof(int);

                    float _dist = u.dist;
                    fwrite(&(_dist),1,sizeof(float),o2uFile);   //??????dist
                    key += sizeof(float);

                    fwrite(&(u.dis),1,sizeof(float),o2uFile);  //??????u ???????????????Ni?????????
                    key += sizeof(float); //1???????????????dis??? (?????????)
                    //?????????????????????
                    int keyNum  = 0;
                    keyNum = u.keywords.size();
                    fwrite(&(keyNum),1,sizeof(int),o2uFile);  //?????????????????????
                    for(int term: u.keywords){
                        fwrite(&(term),1,sizeof(int),o2uFile); //??????????????????
                    }
                    key+=sizeof(int)*keyNum + sizeof(int); //??????1????????????num1?????????num1????????????????????????
                    //???????????????????????????id??????
                    int friendNum =0;
                    friendNum = u.friends.size();
                    fwrite(&(friendNum),1,sizeof(int),o2uFile);  //????????????????????????
                    for(int friend_usr: u.friends){
                        //if(p.id==12156) cout<<checkin_usr<<",";
                        fwrite(&(friend_usr),1,sizeof(int),o2uFile); //?????????????????????
                    }
                    key+=sizeof(int)*friendNum+sizeof(int); //?????????????????????????????????num1????????????*/

                }
            }


        }
    }//?????????????????????????????????

    /*------------------------???Address Idx------------------------*/
    //???o2U??????????????????????????????  poi_size; poi_base; usr_size; usr_base;
    int poi_Size = POIs.size();
    fwrite(&poi_Size, 1, sizeof(int),idxFile); //??????poi_size;

    int poi_base = 4* sizeof(int);
    fwrite(&poi_base, 1, sizeof(int),idxFile); //??????poi_base;

    int user_Size = Users.size();
    fwrite(&user_Size, 1, sizeof(int),idxFile); //??????usr_size

    int usr_base = 4* sizeof(int)+ sizeof(int)*poi_Size;  // ???????????????poi_Size*in
    fwrite(&usr_base, 1, sizeof(int),idxFile); //usr_base


    //???o2u?????????????????????poi??????
    for(int i=0;i<POIs.size();i++){
        int address = POIs[i].addressDataOrgByLeaf;
        fwrite(&address, 1, sizeof(int),idxFile);
    }

    //???o2u?????????????????????user??????
    for(int i=0;i<Users.size();i++){

        int address = Users[i].addressDataOrgByLeaf;
        fwrite(&address, 1, sizeof(int),idxFile);
    }


    fclose(o2uFile);
    fclose(idxFile);
    printf("making o2UFiles and IdxFile (Organized By Leaf) success !\n");



}



/*-----------------------------?????????????????????????????????------------------------------*/
//1. ???gimtree???poi???????????????inverted_list???????????????
//?????????????????????????????????


void makeInvIndex_poi2Usr(const char* fileprefix){
    //??????Gtree?????????poi???keyword ???????????????????????????????????????????????????.p_d_leaf
    bool PtMaxUsing = true;
    //????????????????????????????????????
    char poiInvFile_Name[255]; char poiInvFile_Name2[255];
    sprintf(poiInvFile_Name,"%s.o2u_invList",fileprefix);
    //sprintf(poiInvIdxFile_Name,"%s.poi_invList_idx",fileprefix);
    remove(poiInvFile_Name); // remove existing file ?????????????????????
    //remove(poiInvIdxFile_Name); // remove existing file ?????????????????????
    //??????????????????????????????
    FILE* invFile=fopen(poiInvFile_Name,"w+");
    //??????????????????
    /*FILE* idxFile= fopen(poiInvIdxFile_Name,"w+");
    if(invFile==NULL||idxFile==NULL) cout<<"faile to open";*/

    printf("making o2u invertedList Files\n");

    int RawAddr=0,size;	// changed
    int key = 0;
    int term_size= invListOfPOI.size();  int idx_address = 0;

    key += sizeof(int); //term_size
    key = key + term_size*(sizeof(float)+ 5 * sizeof(int));  //term_idx, every row: <tfIDF-score,
                                                            // baseAddress,
                                                            // addrObase, addrOleaf,
                                                            // addrUbase, addrUleaf >


    //????????????
    map<int,int> termBaseAddress;
    map<int,int> termOBaseAddress; map<int,int> termOLeafAddress;
    map<int,int> termUBaseAddress; map<int,int> termULeafAddress;

    //???term?????????
    map<int, vector<int>>::iterator iter;
    int number = 1;
    //??????object???user???global_vocabulary
    set<int> global_vocabulary;
    set<int> object_vocabulary;
    set<int> user_vocabulary;
    for(iter=invListOfPOI.begin();iter!=invListOfPOI.end();iter++){
        int term_id = iter->first;
        global_vocabulary.insert(term_id);
        object_vocabulary.insert(term_id);
    }
    for(iter=invListOfUser.begin();iter!=invListOfUser.end();iter++){
        int term_id = iter->first;
        global_vocabulary.insert(term_id);
        user_vocabulary.insert(term_id);
    }

    cout<<"object_vocabulary (size="<<object_vocabulary.size()<<"):"<<endl; printSetElements(object_vocabulary);
    cout<<"user_vocabulary (size="<<user_vocabulary.size()<<"):"<<endl; printSetElements(user_vocabulary);
    //getchar();

    bool flag = false;

    //???global_vocabulary
    for(int term_id : global_vocabulary){
        int innerkey = 0;
        int current = key;

        /*if(term_id==178){
            cout<<"t178 termBaseAddress="<<key<<endl;
            getchar();
        }*/
        termBaseAddress[term_id]=key;

        //1. ?????? global_vocabulary?????????keywords????????? object ??? entry ??????
        termOBaseAddress[term_id]=key;
        int _temp = key;
        for(int i=0; i<GTree.size(); i++){

            //?????? object???entry??????
            if(GTree[i].inverted_list_o.count(term_id) >0){

                int postList_size =
                        GTree[i].inverted_list_o[term_id].size(); //write
                //??????term_id???post list?????? ??????????????????
                int _innner = innerkey;

                GTree[i].invertedAddress_list_o[term_id].push_back(innerkey);


                if(!GTree[i].isleaf){ //?????????????????? entry ?????? chile node
                    //cout<<"<post list of t"<<term_id<<" in node"<<i<<": ";
                    for(int node: GTree[i].inverted_list_o[term_id]){
                        //cout<<"n"<<node<<",";  //write
                    }
                    //cout<<">"<<endl;
                }
                else{
                    //cout<<"<post list of t"<<term_id<<" in leaf"<<i<<": ";
                    for(int poi: GTree[i].inverted_list_o[term_id]){
                        //cout<<"p"<<poi<<",";   //write
                    }
                    //cout<<">"<<endl;
                }
                //????????????
                key += sizeof(int)*(1+postList_size);  // size + size* childrens'id
                innerkey += sizeof(int)*(1+postList_size);

            }



        }



        //2. for user: ?????? global_vocabulary?????????keywords????????? usr ??? entry ??????
        termUBaseAddress[term_id] = key;
        int _key_ = key;
        innerkey = 0;
        for(int i=0; i<GTree.size(); i++){


            //?????? user ???entry??????
            if(GTree[i].inverted_list_u.count(term_id) >0){

                //??????term_id???post list?????? ??????????????????
                int _innner = innerkey;

                GTree[i].invertedAddress_list_u[term_id].push_back(innerkey);



                int u_postList_size =
                        GTree[i].inverted_list_u[term_id].size(); //write


                int uLeaf_size = GTree[i].term_usrLeaf_Map[term_id].size();


                if(!GTree[i].isleaf){ //?????????????????? entry ?????? chile node
                    //cout<<"<post list of t"<<term_id<<" in node"<<i<<": ";
                    for(int unode: GTree[i].inverted_list_u[term_id]){
                        //cout<<"n"<<node<<",";  //write
                    }
                    //cout<<">"<<endl;
                }
                else{
                    //cout<<"<post list of t"<<term_id<<" in leaf"<<i<<": ";
                    for(int usr: GTree[i].inverted_list_u[term_id]){
                        //cout<<"p"<<poi<<",";   //write
                    }
                    //cout<<">"<<endl;
                }
                //????????????
                key += sizeof(int)*(1+u_postList_size+1+uLeaf_size);  // chile_size + chile_size* childrens'id + leaf_size+ chile_size*leaves_id
                innerkey += sizeof(int)*(1+u_postList_size+1+uLeaf_size);

            }

        }

        //3. ?????? term_id ????????? object leaf?????????
        termOLeafAddress[term_id]= key;
        int oleaf_size = 0;
        if(poiTerm_leafSet.count(term_id)>0){
            oleaf_size = poiTerm_leafSet[term_id].size();// write
            for(int leaf_id: poiTerm_leafSet[term_id]){
                //write;
            }
            key+= sizeof(int)*(1+oleaf_size);  //????????????
        }


        //4. ?????? term_id ????????? user leaf?????????
        termULeafAddress[term_id] = key;
        int uleaf_size = 0;
        if(usrTerm_leafSet.count(term_id)>0){
            uleaf_size = usrTerm_leafSet[term_id].size();// write
            for(int leaf_id: usrTerm_leafSet[term_id]){
                //write;
            }
            //????????????
            key+= sizeof(int)*(1+uleaf_size);
        }

    }

    //1. ????????????:term_size
    fwrite(&term_size,1, sizeof(int),invFile);
    //2. ????????????-??????????????????every row: <tfIDF-score, addrbase, obase, oleaf,ubase, uleaf>
    int _size = global_vocabulary.size();
    for(int keyword: global_vocabulary){
        //???????????????weight
        int poi_size = POIs.size();
        int IDF_size = invListOfPOI[keyword].size();
        float idfWeight = log2(poi_size*1.0/IDF_size);
        fwrite(&idfWeight,1, sizeof(float),invFile);

        int address_base = -1;
        int address_obase = -1;
        int address_oleaf = -1;
        if(termBaseAddress.count(keyword)>0){
            address_base = termBaseAddress[keyword];
            address_obase = termOBaseAddress[keyword];
            address_oleaf = termOLeafAddress[keyword];
        }
        int address_ubase = -1;
        int address_uleaf = -1;
        if(termUBaseAddress.count(keyword)>0){
            address_ubase = termUBaseAddress[keyword];
            address_uleaf = termULeafAddress[keyword];
        }


        //????????????????????????
        fwrite(&address_base,1, sizeof(int),invFile);
        //??????object?????????
        fwrite(&address_obase,1, sizeof(int),invFile);
        //??????object???????????????
        fwrite(&address_oleaf,1, sizeof(int),invFile);
        //??????user?????????
        fwrite(&address_ubase,1, sizeof(int),invFile);
        //??????user???????????????
        fwrite(&address_uleaf,1, sizeof(int),invFile);
    }

    int key2 = 0+ sizeof(int); //term_size
    key2 += (term_size*5* sizeof(int)); //???????????????????????????term_base, o_base, oleafAddr, u_base, oleafAddr???


    //3. ???????????????????????????????????????
    //???term?????????
    map<int, vector<int>>::iterator iter2;
    for(int term_id : global_vocabulary){

        int innerkey = 0;

        //1. ??????term ?????????object entry ??????
        for(int i=0; i<GTree.size(); i++){


            //1. ??????term ?????????object entry ??????
            if(GTree[i].inverted_list_o.count(term_id)>0){
                int entryList_size = GTree[i].inverted_list_o[term_id].size(); //write
                //??????term_id???post list?????? ??????????????????
                fwrite(&entryList_size,1, sizeof(int),invFile);

                if(!GTree[i].isleaf){ //?????????????????? entry ?????? chile node
                    for(int node: GTree[i].inverted_list_o[term_id]){
                        int entry= node;
                        fwrite(&entry,1, sizeof(int),invFile);
                    }
                }
                else{
                    for(int poi: GTree[i].inverted_list_o[term_id]){
                        //cout<<"p"<<poi<<",";   //write
                        int entry= poi;
                        fwrite(&entry,1, sizeof(int),invFile);
                    }
                }
                //????????????
                key2 = key2 + sizeof(int)*(1+entryList_size);  // size + size* childrens'id
                innerkey = innerkey + sizeof(int)*(1+entryList_size);

            }

        }
        //2. ??????term ?????????user entry ?????? user leaf ??????
        innerkey = 0;
        for(int i=0; i<GTree.size(); i++){
            if(GTree[i].inverted_list_u.count(term_id)>0){
                //??????term_id???child entry ???size
                int entryList_size = GTree[i].inverted_list_u[term_id].size(); //write
                fwrite(&entryList_size,1, sizeof(int),invFile);
                //??????term_id????????????child entry???id
                if(!GTree[i].isleaf){ //?????????????????? entry ?????? chile node
                    for(int unode: GTree[i].inverted_list_u[term_id]){
                        int entry= unode;
                        fwrite(&entry,1, sizeof(int),invFile);
                    }
                }
                else{
                    for(int usr: GTree[i].inverted_list_u[term_id]){
                        //cout<<"p"<<poi<<",";   //write
                        int entry= usr;
                        fwrite(&entry,1, sizeof(int),invFile);
                    }
                }

                //??????term_id????????? usr_leaf???size
                int uLeaf_size = GTree[i].term_usrLeaf_Map[term_id].size();
                fwrite(&uLeaf_size,1, sizeof(int),invFile);
                //??????term_id???????????????usr_leaf id
                for(int uLeaf_id: GTree[i].term_usrLeaf_Map[term_id]){
                    fwrite(&uLeaf_id,1, sizeof(int),invFile);
                }

                //????????????
                key2 = key2 + sizeof(int)*(1+entryList_size+1+uLeaf_size);  // chile_size + chile_size* childrens'id + leaf_size+ chile_size*leaves_id
                innerkey = innerkey + sizeof(int)*(1+entryList_size+1+uLeaf_size);

            }

        }

        //3. ????????????term???????????????poi leaf???id
        if(poiTerm_leafSet.count(term_id)>0){// write
            int leaf_size = poiTerm_leafSet[term_id].size();// write
            fwrite(&leaf_size,1, sizeof(int),invFile);
            for(int leaf_id: poiTerm_leafSet[term_id]){
                //write;
                fwrite(&leaf_id,1, sizeof(int),invFile);
            }
            //????????????
            key2 = key2+ sizeof(int)*(1+leaf_size);

        }

        //4. ????????????term???????????????user leaf???id
        if(usrTerm_leafSet.count(term_id)>0){// write
            int leaf_size = usrTerm_leafSet[term_id].size();// write
            fwrite(&leaf_size,1, sizeof(int),invFile);
            for(int leaf_id: usrTerm_leafSet[term_id]){
                //write;
                fwrite(&leaf_id,1, sizeof(int),invFile);
            }
            //????????????
            key2 = key2+ sizeof(int)*(1+leaf_size);
        }



    }

    fclose(invFile);
    cout<<"makeInvIndex_poi2Usr_new"<<endl;

}


void makeInvIndex_poi2Usr_new(const char* fileprefix){
    //??????Gtree?????????poi???keyword ???????????????????????????????????????????????????.p_d_leaf
    bool PtMaxUsing = true;
    //????????????????????????????????????
    char poiInvFile_Name[255]; char poiInvFile_Name2[255];
    sprintf(poiInvFile_Name,"%s.o2u_invList",fileprefix);
    //sprintf(poiInvIdxFile_Name,"%s.poi_invList_idx",fileprefix);
    remove(poiInvFile_Name); // remove existing file ?????????????????????
    //remove(poiInvIdxFile_Name); // remove existing file ?????????????????????
    //??????????????????????????????
    FILE* invFile=fopen(poiInvFile_Name,"w+");
    //??????????????????
    /*FILE* idxFile= fopen(poiInvIdxFile_Name,"w+");
    if(invFile==NULL||idxFile==NULL) cout<<"faile to open";*/

    printf("making o2u invertedList Files\n");

    int RawAddr=0,size;	// changed
    int key = 0;
    int term_size= invListOfPOI.size();  int idx_address = 0;

    key += sizeof(int); //term_size
    int index_record_size = 7; //address record, every row: <
    // baseAddress
    // addrObase, addrUbase
    // addrOleaf, addrUleaf
    // addrObject,addruser// >: ??? 7 ???
    key = key + term_size*(sizeof(float)+ index_record_size * sizeof(int));  //sizeof(term_weight)+ address_record_size


    //????????????
    unordered_map<int,int> termBaseAddress;
    unordered_map<int,int> termOBaseAddress; unordered_map<int,int> termOLeafAddress; unordered_map<int,int> termObjectAddress;
    unordered_map<int,int> termUBaseAddress; unordered_map<int,int> termULeafAddress; unordered_map<int,int> termUserAddress;

    //???term?????????
    map<int, vector<int>>::iterator iter;
    int number = 1;
    //??????object???user???global_vocabulary
    set<int> global_vocabulary;
    set<int> object_vocabulary;
    set<int> user_vocabulary;
    for(iter=invListOfPOI.begin();iter!=invListOfPOI.end();iter++){
        int term_id = iter->first;
        global_vocabulary.insert(term_id);
        object_vocabulary.insert(term_id);
    }
    for(iter=invListOfUser.begin();iter!=invListOfUser.end();iter++){
        int term_id = iter->first;
        global_vocabulary.insert(term_id);
        user_vocabulary.insert(term_id);
    }

    cout<<"object_vocabulary (size="<<object_vocabulary.size()<<"):"<<endl; printSetElements(object_vocabulary);
    cout<<"user_vocabulary (size="<<user_vocabulary.size()<<"):"<<endl; printSetElements(user_vocabulary);
    //getchar();

    bool flag = false;

    //???global_vocabulary
    for(int term_id : global_vocabulary){
        int innerkey = 0;
        int current = key;
        if(term_id%10==1)
            cout<<"????????????t"<<term_id<<"????????????"<<endl;
        termBaseAddress[term_id]=key;

        //1. ?????? global_vocabulary?????????keywords????????? object ??? entry ??????
        termOBaseAddress[term_id]=key;
        int _temp = key;
        set<int> _ONodes = term_ONodeMap[term_id];
        //for(int i=0; i<GTree.size(); i++){
        for(int i:_ONodes){

            //?????? object???entry??????
            if(GTree[i].inverted_list_o.count(term_id) >0){

                int postList_size =
                        GTree[i].inverted_list_o[term_id].size(); //write
                //??????term_id???post list?????? ??????????????????
                int _innner = innerkey;

                GTree[i].invertedAddress_list_o[term_id].push_back(innerkey);


                if(!GTree[i].isleaf){ //?????????????????? entry ?????? chile node
                    //cout<<"<post list of t"<<term_id<<" in node"<<i<<": ";
                    for(int node: GTree[i].inverted_list_o[term_id]){
                        //cout<<"n"<<node<<",";  //write
                    }
                    //cout<<">"<<endl;
                }
                else{
                    //cout<<"<post list of t"<<term_id<<" in leaf"<<i<<": ";
                    for(int poi: GTree[i].inverted_list_o[term_id]){
                        //cout<<"p"<<poi<<",";   //write
                    }
                    //cout<<">"<<endl;
                }
                //????????????
                key += sizeof(int)*(1+postList_size);  // size + size* children's id
                innerkey += sizeof(int)*(1+postList_size);

            }


        }



        //2. for user: ?????? global_vocabulary?????????keywords????????? usr ??? entry ??????
        termUBaseAddress[term_id] = key;
        int _key_ = key;
        innerkey = 0;
        set<int> _UNodes = term_UNodeMap[term_id];
        //for(int i=0; i<GTree.size(); i++){
        for(int i:_UNodes){

            //?????? user ???entry??????
            if(GTree[i].inverted_list_u.count(term_id) >0){

                //??????term_id???post list?????? ??????????????????
                int _innner = innerkey;

                GTree[i].invertedAddress_list_u[term_id].push_back(innerkey);



                int u_postList_size =
                        GTree[i].inverted_list_u[term_id].size(); //write


                int uLeaf_size = GTree[i].term_usrLeaf_Map[term_id].size();


                if(!GTree[i].isleaf){ //?????????????????? entry ?????? chile node
                    //cout<<"<post list of t"<<term_id<<" in node"<<i<<": ";
                    for(int unode: GTree[i].inverted_list_u[term_id]){
                        //cout<<"n"<<node<<",";  //write
                    }
                    //cout<<">"<<endl;
                }
                else{
                    //cout<<"<post list of t"<<term_id<<" in leaf"<<i<<": ";
                    for(int usr: GTree[i].inverted_list_u[term_id]){
                        //cout<<"p"<<poi<<",";   //write
                    }
                    //cout<<">"<<endl;
                }
                //????????????
                key += sizeof(int)*(1+u_postList_size+1+uLeaf_size);  // chile_size + chile_size* childrens'id + leaf_size+ chile_size*leaves_id
                innerkey += sizeof(int)*(1+u_postList_size+1+uLeaf_size);

            }

        }

        //3. ?????? term_id ????????? object leaf?????????
        termOLeafAddress[term_id]= key;
        int oleaf_size = 0;
        if(poiTerm_leafSet.count(term_id)>0){
            oleaf_size = poiTerm_leafSet[term_id].size();// write
            for(int leaf_id: poiTerm_leafSet[term_id]){
                //write;
            }
            key+= sizeof(int)*(1+oleaf_size);  //????????????
        }


        //4. ?????? term_id ????????? user leaf?????????
        termULeafAddress[term_id] = key;
        int uleaf_size = 0;
        if(usrTerm_leafSet.count(term_id)>0){
            uleaf_size = usrTerm_leafSet[term_id].size();// write
            for(int leaf_id: usrTerm_leafSet[term_id]){
                //write;
            }
            //????????????
            key+= sizeof(int)*(1+uleaf_size);
        }


        //5. ?????? term_id ????????? object ?????????????????????
        termObjectAddress[term_id] = key;
        int object_size = 0;
        if(invListOfPOI[term_id].size()>0){
            object_size = invListOfPOI[term_id].size();

            for(int object_id: invListOfPOI[term_id]){
                //write
            }
            key+= sizeof(int)*(1+object_size);
        }

        //6. ?????? term_id ????????? user ?????????????????????
        termUserAddress[term_id] = key;
        int user_size = 0;
        if(invListOfUser[term_id].size()>0){
            user_size = invListOfUser[term_id].size();
            for(int user_id: invListOfUser[term_id]){
                //write
            }
            key+= sizeof(int)*(1+user_size);
        }

    }

    //1. ????????????:term_size
    fwrite(&term_size,1, sizeof(int),invFile);
    //2. ????????????-??????????????????every row: <tfIDF-score, addrbase, obase, oleaf,ubase, uleaf>
    int _size = global_vocabulary.size();
    for(int keyword: global_vocabulary){
        if(keyword%10==0)
            cout<<"???inverted file??? ??????t"<<keyword<<"????????????????????????"<<endl;
        //???????????????weight
        int poi_size = POIs.size();
        int IDF_size = invListOfPOI[keyword].size();
        float idfWeight = log2(poi_size*1.0/IDF_size);
        if(keyword==1){
            cout<<"find term1! poi_size="<<poi_size<<",IDF_size="<<IDF_size<<",idWeight="<<idfWeight<<endl;
            //getchar();
        }

        fwrite(&idfWeight,1, sizeof(float),invFile);

        int address_base = -1;
        int address_obase = -1;
        int address_oleaf = -1;
        int address_object_inv = -1;
        if(termBaseAddress.count(keyword)>0){
            address_base = termBaseAddress[keyword];
            address_obase = termOBaseAddress[keyword];
            address_oleaf = termOLeafAddress[keyword];
            address_object_inv =termObjectAddress[keyword];
        }
        int address_ubase = -1;
        int address_uleaf = -1;
        int address_user_inv = -1;
        if(termUBaseAddress.count(keyword)>0){
            address_ubase = termUBaseAddress[keyword];
            address_uleaf = termULeafAddress[keyword];
            address_user_inv = termUserAddress[keyword];
        }



        //????????????????????????
        fwrite(&address_base,1, sizeof(int),invFile);

        //??????object?????????
        fwrite(&address_obase,1, sizeof(int),invFile);
        //??????object???????????????
        fwrite(&address_oleaf,1, sizeof(int),invFile);

        //??????user?????????
        fwrite(&address_ubase,1, sizeof(int),invFile);
        //??????user???????????????
        fwrite(&address_uleaf,1, sizeof(int),invFile);

        //??????object?????????????????????????????????
        fwrite(&address_object_inv,1,sizeof(int),invFile);
        //??????user?????????????????????????????????
        fwrite(&address_user_inv,1,sizeof(int),invFile);


    }

    int key2 = 0+ sizeof(int); //term_size
    key2 += term_size*(sizeof(float)+7* sizeof(int)); //????????????record ???????????????term_weight, term_base, o_base, oleafAddr, u_base, oleafAddr, object_addr, user_addr???


    //3. ???????????????????????????????????????
    //???term?????????
    map<int, vector<int>>::iterator iter2;
    for(int term_id : global_vocabulary){

        int innerkey = 0;
        if(term_id%10==0)
            cout<<"??? inverted file?????????t"<<term_id<<"??????????????????????????????"<<endl;

        //1. ??????term ?????????object entry ??????
        set<int> _ONodes = term_ONodeMap[term_id];
        //for(int i=0; i<GTree.size(); i++){
        for(int i:_ONodes){

            //1. ??????term ?????????object entry ??????
            if(GTree[i].inverted_list_o.count(term_id)>0){
                int entryList_size = GTree[i].inverted_list_o[term_id].size(); //write
                //??????term_id???post list?????? ??????????????????
                fwrite(&entryList_size,1, sizeof(int),invFile);

                if(!GTree[i].isleaf){ //?????????????????? entry ?????? chile node
                    for(int node: GTree[i].inverted_list_o[term_id]){
                        int entry= node;
                        fwrite(&entry,1, sizeof(int),invFile);
                    }
                }
                else{
                    for(int poi: GTree[i].inverted_list_o[term_id]){
                        //cout<<"p"<<poi<<",";   //write
                        int entry= poi;
                        fwrite(&entry,1, sizeof(int),invFile);
                    }
                }
                //????????????
                key2 = key2 + sizeof(int)*(1+entryList_size);  // size + size* childrens'id
                innerkey = innerkey + sizeof(int)*(1+entryList_size);

            }

        }
        //2. ??????term ?????????user entry ?????? user leaf ??????
        innerkey = 0;
        set<int> _UNodes = term_UNodeMap[term_id];
        //for(int i=0; i<GTree.size(); i++){
        for(int i:_UNodes){
            if(GTree[i].inverted_list_u.count(term_id)>0){
                //??????term_id???child entry ???size
                int entryList_size = GTree[i].inverted_list_u[term_id].size(); //write
                fwrite(&entryList_size,1, sizeof(int),invFile);
                //??????term_id????????????child entry???id
                if(!GTree[i].isleaf){ //?????????????????? entry ?????? chile node
                    for(int unode: GTree[i].inverted_list_u[term_id]){
                        int entry= unode;
                        fwrite(&entry,1, sizeof(int),invFile);
                    }
                }
                else{
                    for(int usr: GTree[i].inverted_list_u[term_id]){
                        //cout<<"p"<<poi<<",";   //write
                        int entry= usr;
                        fwrite(&entry,1, sizeof(int),invFile);
                    }
                }

                //??????term_id????????? usr_leaf???size
                int uLeaf_size = GTree[i].term_usrLeaf_Map[term_id].size();
                fwrite(&uLeaf_size,1, sizeof(int),invFile);
                //??????term_id???????????????usr_leaf id
                for(int uLeaf_id: GTree[i].term_usrLeaf_Map[term_id]){
                    fwrite(&uLeaf_id,1, sizeof(int),invFile);
                }

                //????????????
                key2 = key2 + sizeof(int)*(1+entryList_size+1+uLeaf_size);  // chile_size + chile_size* childrens'id + leaf_size+ chile_size*leaves_id
                innerkey = innerkey + sizeof(int)*(1+entryList_size+1+uLeaf_size);

            }

        }

        //3. ????????????term???????????????poi leaf???id
        if(poiTerm_leafSet.count(term_id)>0){// write
            int leaf_size = poiTerm_leafSet[term_id].size();// write
            fwrite(&leaf_size,1, sizeof(int),invFile);
            for(int leaf_id: poiTerm_leafSet[term_id]){
                //write;
                fwrite(&leaf_id,1, sizeof(int),invFile);
            }
            //????????????
            key2 = key2+ sizeof(int)*(1+leaf_size);

        }

        //4. ????????????term???????????????user leaf???id
        if(usrTerm_leafSet.count(term_id)>0){// write
            int leaf_size = usrTerm_leafSet[term_id].size();// write
            fwrite(&leaf_size,1, sizeof(int),invFile);
            for(int leaf_id: usrTerm_leafSet[term_id]){
                //write;
                fwrite(&leaf_id,1, sizeof(int),invFile);
            }
            //????????????
            key2 = key2+ sizeof(int)*(1+leaf_size);
        }

        //if(term_id==1){
            //cout<<"going to write posting list for term"<<term_id<<endl;
        //}

        //5. ????????????term???????????????geo-social object???id
        if(invListOfPOI[term_id].size()>0){
            int object_size = invListOfPOI[term_id].size();
            //cout<<"the size of object posting list of term"<<term_id<<"="<<object_size<<endl;
            fwrite(&object_size,1, sizeof(int),invFile);
            key2 += sizeof(int); int _tmp = key2;
            for(int object_id: invListOfPOI[term_id]){
                fwrite(&object_id,1, sizeof(int),invFile);
                int _current = key2;
                key2 += sizeof(int);
            }
        }

        //6.
        if(invListOfUser[term_id].size()>0){
            int user_size = invListOfUser[term_id].size();
            fwrite(&user_size,1, sizeof(int),invFile);
            for(int user_id: invListOfUser[term_id]){
                fwrite(&user_id,1, sizeof(int),invFile);
            }
            key2 = key2 + sizeof(int)*(1+user_size);
        }


    }

    fclose(invFile);

}



//2. ???gimtree?????????????????????????????????


void gimtree_write2disk() {
    // FILE_GTREE_disk
    char disk_file_name[200];
    sprintf(disk_file_name,"%s-disk",FILE_GTREE);
    FILE *fin = fopen(disk_file_name, "wb+");
    int tree_nodeAddress = 0;
    int first_nodeAddress = 0+sizeof(int)*(1+GTree.size()); // node number + node number*(address)
    tree_nodeAddress = first_nodeAddress;
    vector<int> nodeAddressList;
    int point = 0;
    //?????????GIMtree???????????????
    int node_size = GTree.size();
    fwrite(&node_size,1, sizeof(int),fin);
    point+= sizeof(int);
    //?????????????????????????????????????????????????????????????????????
    for(int i=0;i<node_size;i++){

        int current_addr = tree_nodeAddress;
        nodeAddressList.push_back(tree_nodeAddress);
        //????????????i?????????
        fwrite(&tree_nodeAddress,1, sizeof(int),fin);
        point+= sizeof(int);
        int delta_address =0;

        TreeNode& tn = GTree[i];

        //???????????????
        // borders
        int border_size = tn.borders.size();
        delta_address += (1+border_size)* sizeof(int);
        // children
        int children_size = tn.children.size();
        delta_address += (1+children_size)* sizeof(int);
        // isleaf
        int _isleaf = tn.isleaf;
        delta_address += sizeof(bool);
        // leafnodes
        int leafnode_size = tn.leafnodes.size();
        delta_address += (1+leafnode_size)* sizeof(int);
        // father
        int father =  tn.father;
        delta_address += sizeof(int);

        //GTree[node_id].minSocialCount;
        int _tmp = tn.minSocialCount;
        delta_address += sizeof(int);


        // poi vocabulary: oterm_size
        // each row: <term_id, innerAddressInpostListing, object_size>
        int o_vocabulary_size = tn.objectUKeySet.size();
        delta_address += sizeof(int)*(1+3*o_vocabulary_size); // vocabulary_size, {<term_id, oinneraddress, count>...}

        // usr vocabulary: uterm_size
        // (term_id, innerAddressInpostListing, user_size)
        int u_vocabulary_size = tn.userUKeySet.size();
        delta_address += sizeof(int)*(1+3*u_vocabulary_size); // vocabulary_size, {<term_id, uinneraddress, count>...}

        //update address
        tree_nodeAddress += delta_address;
    }


    //??????Gtree??????????????????????????????
    for(int i=0;i<GTree.size();i++){
        ///cout<<"?????????"<<i<<endl;

        TreeNode& tn = GTree[i];
        //record node address
        int _current = point;
        // borders
        int border_size = tn.borders.size();
        fwrite(&border_size,1, sizeof(int),fin);
        point+= sizeof(int);
        for (int i = 0; i < border_size; i++) {
            int border_id = tn.borders[i];
            fwrite(&border_id,1, sizeof(int),fin);
            point+= sizeof(int);
        }

        // children
        int children_size = tn.children.size();
        fwrite(&children_size,1, sizeof(int),fin);
        point+= sizeof(int);
        for (int i = 0; i < children_size; i++) {
            int childNode_id = tn.children[i];
            fwrite(&childNode_id,1,sizeof(int),fin);
            point+= sizeof(int);

        }

        // isleaf
        bool isleaf = tn.isleaf;
        fwrite(&isleaf,1, sizeof(bool),fin);
        point+= sizeof(bool);

        // leafnodes
        int leafnode_size = tn.leafnodes.size();
        fwrite(&leafnode_size,1, sizeof(int),fin);
        for (int i = 0; i < leafnode_size; i++) {
            int leaf_id = tn.leafnodes[i];
            fwrite(&leaf_id,1, sizeof(int),fin);
        }
        point+= sizeof(int)*(1+leafnode_size);
        // father
        int father_id = tn.father;
        fwrite(&father_id,1, sizeof(int),fin);
        point+= sizeof(int);


        //GTree[node_id].minSocialCount;
        int _check_count = tn.minSocialCount;
        fwrite(&_check_count,1,sizeof(int),fin);
        point+= sizeof(int);


        //poi keyword vocabulary
        int keyword_size = tn.objectUKeySet.size();
        fwrite(&keyword_size,1, sizeof(int),fin);
        point+= sizeof(int);
        ///cout<<"objectUKeySet:"<<endl;
        for(int keyword: tn.objectUKeySet){
            //okeyword
            ///cout<<"t"<<keyword<<endl;
            fwrite(&keyword,1, sizeof(int),fin);
            //innerAddress
            int innerAddr = tn.invertedAddress_list_o[keyword][0];
            fwrite(&innerAddr,1, sizeof(int),fin);


            //number of objects having this keyword
            int object_size =  GTree[i].term_object_Map[keyword].size();
            fwrite(&object_size,1, sizeof(int),fin);



            point+= sizeof(int)*3;
        }


        //usr keyword vocabulary
        int keyword_size2 = tn.userUKeySet.size();
        fwrite(&keyword_size2,1, sizeof(int),fin);
        point+= sizeof(int);
        ///cout<<"userUKeySet"<<endl;
        for(int keyword: tn.userUKeySet){
            //ukeyword
            ///cout<<"t"<<keyword<<endl;
            fwrite(&keyword,1, sizeof(int),fin);
            //innerAddress
            int innerAddr = tn.invertedAddress_list_u[keyword][0];
            fwrite(&innerAddr,1, sizeof(int),fin);
            //number of users having this keyword
            int user_size = GTree[i].term_usr_Map[keyword].size();
            fwrite(&user_size,1, sizeof(int),fin);

            //number of usr_leaf and ids of these leaves
            int usrLeaf_size = GTree[i].term_usrLeaf_Map[keyword].size();

            point+= sizeof(int)*3;
        }




    }


    fclose(fin);

    printf("saving gimtree-disk file success!\n");


}

void gimtree_write2disk(const char* fileprefix) {
    // FILE_GTREE_disk
    char disk_file_name[200];
    sprintf(disk_file_name,"%s.gtree-disk",fileprefix);
    FILE *fin = fopen(disk_file_name, "wb+");
    int tree_nodeAddress = 0;
    int first_nodeAddress = 0+sizeof(int)*(1+GTree.size()); // node number + node number*(address)
    tree_nodeAddress = first_nodeAddress;
    vector<int> nodeAddressList;
    int point = 0;
    //?????????GIMtree???????????????
    int node_size = GTree.size();
    fwrite(&node_size,1, sizeof(int),fin);
    point+= sizeof(int);
    //?????????????????????????????????????????????????????????????????????
    for(int i=0;i<node_size;i++){

        int current_addr = tree_nodeAddress;
        nodeAddressList.push_back(tree_nodeAddress);
        //????????????i?????????
        fwrite(&tree_nodeAddress,1, sizeof(int),fin);
        point+= sizeof(int);
        int delta_address =0;

        TreeNode& tn = GTree[i];

        //???????????????
        // borders
        int border_size = tn.borders.size();
        delta_address += (1+border_size)* sizeof(int);
        // children
        int children_size = tn.children.size();
        delta_address += (1+children_size)* sizeof(int);
        // isleaf
        int _isleaf = tn.isleaf;
        delta_address += sizeof(bool);
        // leafnodes
        int leafnode_size = tn.leafnodes.size();
        delta_address += (1+leafnode_size)* sizeof(int);
        // father
        int father =  tn.father;
        delta_address += sizeof(int);

        //GTree[node_id].minSocialCount;
        int _tmp = tn.minSocialCount;
        delta_address += sizeof(int);


        // poi vocabulary: oterm_size
        // each row: <term_id, innerAddressInpostListing, object_size>
        int o_vocabulary_size = tn.objectUKeySet.size();
        delta_address += sizeof(int)*(1+3*o_vocabulary_size); // vocabulary_size, {<term_id, oinneraddress, count>...}

        // usr vocabulary: uterm_size
        // (term_id, innerAddressInpostListing, user_size)
        int u_vocabulary_size = tn.userUKeySet.size();
        delta_address += sizeof(int)*(1+3*u_vocabulary_size); // vocabulary_size, {<term_id, uinneraddress, count>...}

        //update address
        tree_nodeAddress += delta_address;
    }


    //??????Gtree??????????????????????????????
    for(int i=0;i<GTree.size();i++){

        TreeNode& tn = GTree[i];
        //record node address
        int _current = point;
        // borders
        int border_size = tn.borders.size();
        fwrite(&border_size,1, sizeof(int),fin);
        point+= sizeof(int);
        for (int i = 0; i < border_size; i++) {
            int border_id = tn.borders[i];
            fwrite(&border_id,1, sizeof(int),fin);
            point+= sizeof(int);
        }

        // children
        int children_size = tn.children.size();
        fwrite(&children_size,1, sizeof(int),fin);
        point+= sizeof(int);
        for (int i = 0; i < children_size; i++) {
            int childNode_id = tn.children[i];
            fwrite(&childNode_id,1,sizeof(int),fin);
            point+= sizeof(int);

        }

        // isleaf
        bool isleaf = tn.isleaf;
        fwrite(&isleaf,1, sizeof(bool),fin);
        point+= sizeof(bool);

        // leafnodes
        int leafnode_size = tn.leafnodes.size();
        fwrite(&leafnode_size,1, sizeof(int),fin);
        for (int i = 0; i < leafnode_size; i++) {
            int leaf_id = tn.leafnodes[i];
            fwrite(&leaf_id,1, sizeof(int),fin);
        }
        point+= sizeof(int)*(1+leafnode_size);
        // father
        int father_id = tn.father;
        fwrite(&father_id,1, sizeof(int),fin);
        point+= sizeof(int);


        //GTree[node_id].minSocialCount;
        int _check_count = tn.minSocialCount;
        fwrite(&_check_count,1,sizeof(int),fin);
        point+= sizeof(int);


        //poi keyword vocabulary
        int keyword_size = tn.objectUKeySet.size();
        fwrite(&keyword_size,1, sizeof(int),fin);
        point+= sizeof(int);
        for(int keyword: tn.objectUKeySet){
            //okeyword
            fwrite(&keyword,1, sizeof(int),fin);
            //innerAddress
            int innerAddr = tn.invertedAddress_list_o[keyword][0];
            fwrite(&innerAddr,1, sizeof(int),fin);


            //number of objects having this keyword
            int object_size =  GTree[i].term_object_Map[keyword].size();
            fwrite(&object_size,1, sizeof(int),fin);



            point+= sizeof(int)*3;
        }


        //usr keyword vocabulary
        int keyword_size2 = tn.userUKeySet.size();
        fwrite(&keyword_size2,1, sizeof(int),fin);
        point+= sizeof(int);
        for(int keyword: tn.userUKeySet){
            //ukeyword
            fwrite(&keyword,1, sizeof(int),fin);
            //innerAddress
            int innerAddr = tn.invertedAddress_list_u[keyword][0];
            fwrite(&innerAddr,1, sizeof(int),fin);
            //number of users having this keyword
            int user_size = GTree[i].term_usr_Map[keyword].size();
            fwrite(&user_size,1, sizeof(int),fin);

            //number of usr_leaf and ids of these leaves
            int usrLeaf_size = GTree[i].term_usrLeaf_Map[keyword].size();

            point+= sizeof(int)*3;
        }




    }


    fclose(fin);

    printf("saving gimtree-disk file success!\n");


}



void save_socialMatrix(){

    char fileprefix[200];
    sprintf(fileprefix,"../../../data/%s/%s",dataset_name,road_data);
    char smFile_Name[255];
    sprintf(smFile_Name,"%s.sm",fileprefix);
    remove(smFile_Name); // remove existing file ?????????????????????

    FILE* smFile=fopen(smFile_Name,"w+");

    for(int i=0;i<GTree.size();i++){

        //minDis
        float _minDis = GTree[i].minDis;
        fwrite(&_minDis,1, sizeof(float),smFile);

        //unionFriend
        int unionFriend_Size = 0;
        unionFriend_Size = GTree[i].unionFriends.size();
        fwrite(&unionFriend_Size, 1, sizeof(int),smFile);
        for(int f: GTree[i].unionFriends){
            fwrite(&f, 1, sizeof(int),smFile);
        }

        //maxSocialCount
        int _maxSocialCount = GTree[i].maxSocialCount;
        fwrite(&_maxSocialCount,1, sizeof(int),smFile);

        //minSocialCount
        int _minSocialCount = GTree[i].minSocialCount;
        fwrite(&_minSocialCount,1, sizeof(int),smFile);
    }
    fclose(smFile);
    printf("saving save_socialMatrix success!\n");
}


void save_socialMatrix(const char* fileprefix){

    char smFile_Name[255];
    sprintf(smFile_Name,"%s.sm",fileprefix);
    remove(smFile_Name); // remove existing file ?????????????????????

    FILE* smFile=fopen(smFile_Name,"w+");

    for(int i=0;i<GTree.size();i++){

        //minDis
        float _minDis = GTree[i].minDis;
        fwrite(&_minDis,1, sizeof(float),smFile);

        //unionFriend
        int unionFriend_Size = 0;
        unionFriend_Size = GTree[i].unionFriends.size();
        fwrite(&unionFriend_Size, 1, sizeof(int),smFile);
        for(int f: GTree[i].unionFriends){
            fwrite(&f, 1, sizeof(int),smFile);
        }

        //maxSocialCount
        int _maxSocialCount = GTree[i].maxSocialCount;
        fwrite(&_maxSocialCount,1, sizeof(int),smFile);

        //minSocialCount
        int _minSocialCount = GTree[i].minSocialCount;
        fwrite(&_minSocialCount,1, sizeof(int),smFile);
    }
    fclose(smFile);
    printf("saving save_socialMatrix success!\n");
}



void load_socialMatrix(){

    char fileprefix[200];
    sprintf(fileprefix,"../../../data/%s/%s",dataset_name,road_data);

    char smFile_Name[255];
    sprintf(smFile_Name,"%s.sm",fileprefix);
    //remove(smFile_Name); // remove existing file ?????????????????????


    FILE *fin = fopen(smFile_Name, "rb");
    //int *buf = new int[5];
    int friend_Size = -1; float _minDis = -1.0;
    int node_id = 0;

    while (fread(&_minDis, sizeof(float), 1, fin)) {  //minDis

        GTree[node_id].minDis = _minDis;

        //union friend size
        fread(&friend_Size, sizeof(int), 1, fin);

        int size = friend_Size;
        int *buf = new int[size];
        //union friend
        fread(buf, sizeof(int), size, fin);
        for (int i = 0; i < size; i++) {
            int f = buf[i];
            GTree[node_id].unionFriends_load.insert(f);
        }

        //maxSocialCount
        int _maxSocialCount = -1;
        fread(&_maxSocialCount, sizeof(int), 1, fin);
        GTree[node_id].maxSocialCount_load = _maxSocialCount;

        //minSocialCount
        int _minSocialCount = INT32_MAX;
        fread(&_minSocialCount, sizeof(int), 1, fin);
        GTree[node_id].minSocialCount_load = _minSocialCount;


        node_id++;
        if(node_id==GTree.size()) break;
    }
    fclose(fin);
    printf("load_socialMatrix success!\n");

}


/*

void testBtree(BTree* test_tree) {
    //???????????????????????????????????????????????????????????????.p_d???????????????.p_bt???

    int PtNum = PtTree->UserField;                         // UserField ???????????????
    cout<<"PtNum="<<PtNum<<endl;
    //getchar();
    int  treekey = -1; int data= -1;
    for(int i= 0; i<100; i++){
        data = pointQuery(test_tree,i, treekey);
        cout<<"data = "<<data<<", treekey="<<treekey<<endl;
    }


    //getchar();

}

*/



#define LOADINFO_H

#endif //LOADINFO_H
