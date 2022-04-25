#ifndef CONFIG_H


//**************************map， 对应数据集的选择
#define  LasVegas1
#define  Brightkite1
#define  Gowalla1
#define  Twitter

#define LOAD_START gettimeofday( &lv, NULL ); ls = lv.tv_sec * 1000000 + lv.tv_usec;
#define LOAD_END gettimeofday( &lv, NULL ); le = lv.tv_sec * 1000000 + lv.tv_usec;
#define LOAD_PRINT(T) printf("%s runtime: %lld (ms)\r\n", (#T), (le - ls)/1000 );

#define  TRACK1
#define  NAVIE_DEBUG1

int inv_idx_record_size = 7;
#define  frequency_threshold 500
#define  dist_scale_dj 10000
#define  textual_capacity 12   //每个兴趣点的关键词最大个数

#define Global_Inf1
#define MC_LOG


//磁盘访问标志
#define  DiskAccess
#define  NONGIVENKEY

#define  small1


#define HybridHash
#define wordInter 20
//#define urau  5

#define PRINTINFO1
#define DEFAULT_KEYWORDSIZE  6

#define DEFAULT_NUMBEROFSTORE     50  //100  the size of candidate list 25
#define DEFAULT_RATIO 0.5
//the largest expansion threthold in leaf lcl computing
#define EXPANSION_MAX   5000 //(km)

#define AGGRESSIVE_I
#define AGGRESSIVE_II

//实验测试参数
#define DEFAULT_BUDGET  5



#define  Multicore_  // this is the flag for optimized running using multiple threads
#define  scratch
#define  DJ
#define  B2U_   //(B2U方法不可行，弃之)

//打印中间结果，辅助调试

#define  TRACKLEAF1
#define  LCL_LOG_
#define  LCLLeaf_LOG
#define  VERIFICATION_LOG1
#define  PRUNELOG1
#define  LEAF_UNPRUNELOG

#define  RIS_TRACK1
#define  Selection_TRACK1
#define  IO_TRACK_

//topk相关
#define  TOPKTRACK1
#define  TRACKDETAIL1
#define   TRACKNVDTOPK1

#define   WhileTRACK1
#define  TRACKEXTEND

//reverse相关
#define   TIME_NVD1
#define   LEAF_PRUNELOG1  //group
#define   TRACKV2P11
#define   DEBUG11
#define   TRACKNVDBATCH1
#define   TRACKNVDBATCH_VERIFY1
#define   TRACKBorderSB11
#define   TRACK_FILTER1      //filter的while循环调试flag
#define   TRACK_BorderSB_NVD1   //border based sb list computation 调试flag
#define   Verification_DEBUG1
#define   userspacity
#define   NODE_SPARCITY_VALUE     1
#define   LONG_TERM   5000 //(km)


#define  WriteDATA_TRACK
#define  TESTCODE
#define  TESTGIMTree
#define  TESTNVD_
#define  TEST_REVERSE
#define  MC_ESTIMATION




#define  datanew

#define MAXTEXTUAL 500

#define DATASET_HOME  "../../../data/"
#define CODE_HOME  "../../exp/results/"
#define WEIGHT_INFLATE_FACTOR 1



#ifdef LasVegas
#define dataset_name "LasVegas"
#define road_data "LasVegas_data"
#define VertexNum 43401
#define EdgeNum 358204
#define UserID_MaxKey 1512459 //(1512460为社交网络节点个数，实际在LasVegas有评论记录的用户为：28866个）
#define reviewed_userNum 28866
#define poi_num 28851
#define vocabularySize 4237
#define posting_size_threshold 30
#define rau 5
#define   NODE_SPARCITY_VALUE_BATCH     2
#define   NODE_SPARCITY_VALUE_HEURIS    3
#define   uterm_SPARCITY_VALUE     800
#define   SPARCITY_VALUE     2

size_t  graph_vertex_num = 58228;

#define  social_UserID_MaxKey 1512460
#define UserID_MaxKey2 28877
#define DATA_HOME
#define DEFAULT_THRETHEHOLD   5

#define DEFAULT_KQ_SIZE   40
#define DEFAULT_Q_SIZE  60

#define DEFAULT_K   20
#define DEFAULT_ALPHA    0.6
#define DEFAULT_A   6


#define DEFAULT_THRETHEHOLD_MAX   5000


// config for gtree

#define FILE_NODE  "../../../data/LasVegas/road/LasVegas.cnode"
#define FILE_EDGE  "../../../data/LasVegas/road/LasVegas.cedge"
// set all edge weight to 1(unweighted graph)
#define ADJWEIGHT_SET_TO_ALL_ONE true
// we assume edge weight is integer, thus (input edge) * WEIGHT_INFLATE_FACTOR = (our edge weight)
// gtree fanout
#define PARTITION_PART 4
// gtree leaf node capacity = tau(in paper)
#define LEAF_CAP 128
// gtree index disk storage
#define FILE_NODES_GTREE_PATH  "../../../data/LasVegas/LasVegas.paths"
#define FILE_GTREE  "../../../data/LasVegas/LasVegas.gtree"
#define FILE_SM  "../../../data/LasVegas/LasVegas.sm"
#define FILE_ONTREE_MIND  "../../../data/LasVegas/LasVegas.minds"
// input
#define FILE_OBJECT  "../../../data/LasVegas/LasVegas.object"
#define FILE_USER  "../../../data/LasVegas/LasVegas.user"
#define FILE_POI  "../../../data/LasVegas/LasVegas.poi"
#define FILE_FRIENDSHIP  "../../../data/LasVegas/LasVegas.friendship"
#define FILE_CHECKIN  "../../../data/LasVegas/LasVegas.checkin"


#define FILE_LOCATION "../../../data/LasVegas/LasVegas.location200"

#define FILE_DISBOUND "../../../data/LasVegas/LasVegas.DISBOUND"
#define FILE_BorderPair "../../../data/LasVegas/LasVegas.borderPair"
#define FILE_MAXSOCIAL "../../../data/LasVegas/LasVegas.MAXSOCIAL"
#define FILE_TEXTUALBOUND "../../../data/LasVegas/LasVegas.TEXTUALBOUND"
#define FILE_QUERY "../../../data/LasVegas/LasVegas.query"


#define FILE_GIMTREE "../../../data/LasVegas/LasVegas.gimtree"
#define FILE_USERSCORE "../../../data/LasVegas/pre-process/LasVegas.userScore"
#define FILE_USERTopKResults "../../../data/LasVegas/pre-process/LasVegas.userTokResults"
#define FILE_POIReverseResults "../../../data/LasVegas/pre-process/LasVegas.poiReverseResults"
#define FILE_LargePOIReverse "../../../data/LasVegas/pre-process/LasVegas.largeReverseResults"
#define FILE_INVLISTOFPOI   "../../../data/LasVegas/pre-process/LasVegas.poiInv"
#define FILE_INVLISTOFUSR   "../../../data/LasVegas/pre-process/LasVegas.usrInv"


#define FILE_EXP "../exp/results/LasVegas/"
#define FILE_LOG "../exp/results/LasVegas/"

#endif






#ifdef Brightkite
#define dataset_name "Brightkite"
#define road_data "Brightkite_data"
#define VertexNum 175813
#define EdgeNum 358204
#define UserID_MaxKey 58228
#define reviewed_userNum 58228
#define poi_num 364059
#define vocabularySize 31112
#define posting_size_threshold 40
#define rau 10
#define   NODE_SPARCITY_VALUE_BATCH     1
#define   NODE_SPARCITY_VALUE_HEURIS     4
#define   uterm_SPARCITY_VALUE     30
#define   SPARCITY_VALUE     100   // this is for group

size_t  graph_vertex_num = 58228;
#define DEFAULT_THRETHEHOLD   300


#define DEFAULT_KQ_SIZE   15
#define DEFAULT_Q_SIZE  10
#define DEFAULT_K   20
#define DEFAULT_ALPHA    0.6
#define DEFAULT_ALG    1
#define DEFAULT_A   6


#define DEFAULT_THRETHEHOLD_MAX   4000



//数据集源文件所在文件夹


//数据集（未处理过的）源文件
#define FILE_SOURCE  "../../../data/Brightkite/source_dataset/"

//路网部分数据集 (主要包括 路网点集 + 路网边集 + gtree索引 + phl索引)
#define FILE_INPUT_ROAD  "../../../data/Brightkite/road_network"
#define FILE_NODE  "../../../data/Brightkite/Brightkite.cnode"
#define FILE_EDGE  "../../../data/Brightkite/Brightkite.cedge"
#define FILE_NODES_GTREE_PATH  "../../../data/Brightkite/Brightkite.paths"
#define FILE_GTREE  "../../../data/Brightkite/Brightkite.gtree"
#define FILE_ONTREE_MIND  "../../../data/Brightkite/Brightkite.minds"

// 双色体部分数据集 (主要包括 LBSN users + poi of interests)
#define FILE_BICHROMATIC  "../../../data/Brightkite/bichromatic"
#define FILE_OBJECT  "../../../data/Brightkite/Brightkite.object"
#define FILE_USER  "../../../data/Brightkite/Brightkite.user"
#define FILE_POI  "../../../data/Brightkite/Brightkite.poi"

//社交部分数据集(主要包括 social network + check-in records )
#define FILE_SOCIAL  "../../../data/Brightkite/social"
#define FILE_FRIENDSHIP  "../../../data/Brightkite/Brightkite.friendship"
#define FILE_CHECKIN  "../../../data/Brightkite/Brightkite.checkin"

// set all edge weight to 1(unweighted graph)
#define ADJWEIGHT_SET_TO_ALL_ONE true
// we assume edge weight is integer, thus (input edge) * WEIGHT_INFLATE_FACTOR = (our edge weight)
// gtree fanout
#define PARTITION_PART 4
// gtree leaf node capacity = tau(in paper)
#define LEAF_CAP 128
// gtree index disk storage



//附加预处理信息输出文件所在
#define FILE_LOCATION "../../../data/LV/LV.location200"
#define FILE_DISBOUND "../../../data/Brightkite/Brightkite.DISBOUND"
#define FILE_BorderPair "../../../data/Brightkite/Brightkite.borderPair"
#define FILE_MAXSOCIAL "../../../data/Brightkite/Brightkite.MAXSOCIAL"
#define FILE_TEXTUALBOUND "../../../data/LV/LV.TEXTUALBOUND"
#define FILE_QUERY "../../../data/LV/LV.query"

#define FILE_GIMTREE "../../../data/Brightkite/Brightkite.gimtree"
#define FILE_USERSCORE "../../../data/Brightkite/pre-process/Brightkite.userScore"
#define FILE_USERTopKResults "../../../data/Brightkite/pre-process/Brightkite.userTokResults"
#define FILE_POIReverseResults "../../../data/Brightkite/pre-process/Brightkite.poiReverseResults"
#define FILE_LargePOIReverse "../../../data/Brightkite/pre-process/Brightkite.largeReverseResults"
#define FILE_INVLISTOFPOI   "../../../data/Brightkite/pre-process/Brightkite.poiInv"
#define FILE_INVLISTOFUSR   "../../../data/Brightkite/pre-process/Brightkite.usrInv"

#define FILE_EXP "../exp/results/Brightkite/"
#define FILE_LOG "../exp/results/Brightkite/"
#endif



#ifdef Gowalla
#define dataset_name "Gowalla"
#define road_data "Gowalla_data"
#define VertexNum 1070376
#define UserID_MaxKey 196591
#define reviewed_userNum 196591
#define poi_num 1280970
#define vocabularySize 115618
#define posting_size_threshold 60
#define rau 15
#define   NODE_SPARCITY_VALUE_BATCH     1
#define   NODE_SPARCITY_VALUE_HEURIS     3
#define   uterm_SPARCITY_VALUE     40
 #define  SPARCITY_VALUE    200

size_t  graph_vertex_num = 196591;
#define DEFAULT_THRETHEHOLD   300


#define DEFAULT_KQ_SIZE   15
#define DEFAULT_Q_SIZE  200
#define DEFAULT_K   20
#define DEFAULT_ALPHA    0.6
#define DEFAULT_ALG    1
#define DEFAULT_A   6

#define DEFAULT_THRETHEHOLD_MAX   10000

//数据集源文件所在文件夹

//数据集（未处理过的）源文件
#define FILE_SOURCE  "../../../data/Gowalla/source_dataset/"

//路网部分数据集 (主要包括 路网点集 + 路网边集 + gtree索引 + phl索引)
#define FILE_INPUT_ROAD  "../../../data/Gowalla/road_network"
#define FILE_NODE  "../../../data/Gowalla/Gowalla.cnode"
#define FILE_EDGE  "../../../data/Gowalla/Gowalla.cedge"
#define FILE_NODES_GTREE_PATH  "../../../data/Gowalla/Gowalla.paths"
#define FILE_GTREE  "../../../data/Gowalla/Gowalla.gtree"
#define FILE_ONTREE_MIND  "../../../data/Gowalla/Gowalla.minds"

// 双色体部分数据集 (主要包括 LBSN users + poi of interests)
#define FILE_BICHROMATIC  "../../../data/Gowalla/bichromatic"
#define FILE_OBJECT  "../../../data/Gowalla/Gowalla.object"
#define FILE_USER  "../../../data/Gowalla/Gowalla.user"
#define FILE_POI  "../../../data/Gowalla/Gowalla.poi"

//社交部分数据集(主要包括 social network + check-in records )
#define FILE_SOCIAL  "../../../data/Gowalla/social"
#define FILE_FRIENDSHIP  "../../../data/Gowalla/Gowalla.friendship"
#define FILE_CHECKIN  "../../../data/Gowalla/Gowalla.checkin"



// set all edge weight to 1(unweighted graph)
#define ADJWEIGHT_SET_TO_ALL_ONE true
// we assume edge weight is integer, thus (input edge) * WEIGHT_INFLATE_FACTOR = (our edge weight)
// gtree fanout
#define PARTITION_PART 4
// gtree leaf node capacity = tau(in paper)
#define LEAF_CAP 256 //512
// gtree index disk storage



//附加预处理信息输出文件所在
#define FILE_LOCATION "../../../data/LV/LV.location200"
#define FILE_DISBOUND "../../../data/Gowalla/Gowalla.DISBOUND"
#define FILE_BorderPair "../../../data/Gowalla/Gowalla.borderPair"
#define FILE_MAXSOCIAL "../../../data/Gowalla/Gowalla.MAXSOCIAL"
#define FILE_TEXTUALBOUND "../../../data/LV/LV.TEXTUALBOUND"
#define FILE_QUERY "../../../data/LV/LV.query"

#define FILE_GIMTREE "../../../data/Gowalla/Gowalla.gimtree"
#define FILE_USERSCORE "../../../data/Gowalla/pre-process/Gowalla.userScore"
#define FILE_USERTopKResults "../../../data/Gowalla/pre-process/Gowalla.userTokResults"
#define FILE_POIReverseResults "../../../data/Gowalla/pre-process/Gowalla.poiReverseResults"
#define FILE_LargePOIReverse "../../../data/Gowalla/pre-process/Gowalla.largeReverseResults"
#define FILE_INVLISTOFPOI   "../../../data/Gowalla/pre-process/Gowalla.poiInv"
#define FILE_INVLISTOFUSR   "../../../data/Gowalla/pre-process/Gowalla.usrInv"

#define FILE_EXP "../exp/results/Gowalla/"
#define FILE_LOG "../exp/results/Gowalla/"
#endif




#ifdef Twitter
#define dataset_name "Twitter"
#define road_data "Twitter_data"
#define VertexNum 1890815
#define UserID_MaxKey 554372 //(554372为社交网络节点个数，实际在.poi文件中有评论记录的用户同样为：554372个）
#define social_UserID_MaxKey 554372
#define reviewed_userNum 554372
#define poi_num 3266619
#define vocabularySize 83618
#define posting_size_threshold 200  //jins
#define rau 15
#define   NODE_SPARCITY_VALUE_BATCH     1
#define   NODE_SPARCITY_VALUE_HEURIS     3
#define   uterm_SPARCITY_VALUE     40
#define  SPARCITY_VALUE    200

size_t  graph_vertex_num = 1768149;
#define DEFAULT_THRETHEHOLD   300


#define DEFAULT_KQ_SIZE   15
#define DEFAULT_Q_SIZE  200
#define DEFAULT_K   20
#define DEFAULT_ALPHA    0.6
#define DEFAULT_ALG    1
#define DEFAULT_A   6

#define DEFAULT_THRETHEHOLD_MAX   10000

//数据集源文件所在文件夹

//数据集（未处理过的）源文件
#define FILE_SOURCE  "../../../data/Twitter/source_dataset/"

//路网部分数据集 (主要包括 路网点集 + 路网边集 + gtree索引 + phl索引)
#define FILE_INPUT_ROAD  "../../../data/Twitter/road_network"
#define FILE_NODE  "../../../data/Twitter/Twitter.cnode"
#define FILE_EDGE  "../../../data/Twitter/Twitter.cedge"
#define FILE_NODES_GTREE_PATH  "../../../data/Twitter/Twitter.paths"
#define FILE_GTREE  "../../../data/Twitter/Twitter.gtree"
#define FILE_ONTREE_MIND  "../../../data/Twitter/Twitter.minds"

// 双色体部分数据集 (主要包括 LBSN users + poi of interests)
#define FILE_BICHROMATIC  "../../../data/Twitter/bichromatic"
#define FILE_OBJECT  "../../../data/Twitter/Twitter.object"
#define FILE_USER  "../../../data/Twitter/Twitter.user"
#define FILE_POI  "../../../data/Twitter/Twitter.poi"

//社交部分数据集(主要包括 social network + check-in records )
#define FILE_SOCIAL  "../../../data/Twitter/social"
#define FILE_FRIENDSHIP  "../../../data/Twitter/Twitter.friendship"
#define FILE_CHECKIN  "../../../data/Twitter/Twitter.checkin"



// set all edge weight to 1(unweighted graph)
#define ADJWEIGHT_SET_TO_ALL_ONE true
// we assume edge weight is integer, thus (input edge) * WEIGHT_INFLATE_FACTOR = (our edge weight)
// gtree fanout
#define PARTITION_PART 4
// gtree leaf node capacity = tau(in paper)
#define LEAF_CAP 256 //512
// gtree index disk storage



//附加预处理信息输出文件所在
#define FILE_LOCATION "../../../data/Twitter/Twitter.location200"
#define FILE_DISBOUND "../../../data/Twitter/Twitter.DISBOUND"
#define FILE_BorderPair "../../../data/Twitter/Twitter.borderPair"
#define FILE_MAXSOCIAL "../../../data/Twitter/Twitter.MAXSOCIAL"
#define FILE_TEXTUALBOUND "../../../data/LV/LV.TEXTUALBOUND"
#define FILE_QUERY "../../../data/LV/LV.query"

#define FILE_GIMTREE "../../../data/Twitter/Twitter.gimtree"
#define FILE_USERSCORE "../../../data/Twitter/pre-process/Twitter.userScore"
#define FILE_USERTopKResults "../../../data/Twitter/pre-process/Twitter.userTokResults"
#define FILE_POIReverseResults "../../../data/Twitter/pre-process/Twitter.poiReverseResults"
#define FILE_LargePOIReverse "../../../data/Twitter/pre-process/Twitter.largeReverseResults"
#define FILE_INVLISTOFPOI   "../../../data/Twitter/pre-process/Twitter.poiInv"
#define FILE_INVLISTOFUSR   "../../../data/Twitter/pre-process/Twitter.usrInv"

#define FILE_EXP "../exp/results/Twitter/"
#define FILE_LOG "../exp/results/Twitter/"
#endif





#ifdef LV
#define dataset_name "LV"
#define road_data "LV_data"
#define VertexNum 25745
#define vocabularySize 664

#ifdef small
#define UserID_MaxKey 1000
#endif

#ifndef small
#endif

#define  social_UserID_MaxKey 44292
size_t  graph_vertex_num = 1000;
#define poi_num 13468
#define UserID_MaxKey 1001
#define UserID_MaxKey2 28877
#define DATA_HOME
#define DEFAULT_THRETHEHOLD   5

#define DEFAULT_KQ_SIZE   40
#define DEFAULT_Q_SIZE  60

#define DEFAULT_K   20
#define DEFAULT_ALPHA    0.6
#define DEFAULT_A   6
#define SPARCITY_VALUE     5


#define DEFAULT_THRETHEHOLD_MAX   5000
#define posting_size_threshold 20
#define rau 5


// config for gtree

#define FILE_NODE  "../../../data/LV/road/LV.cnode"
#define FILE_EDGE  "../../../data/LV/road/LV.cedge"
// set all edge weight to 1(unweighted graph)
#define ADJWEIGHT_SET_TO_ALL_ONE true
// we assume edge weight is integer, thus (input edge) * WEIGHT_INFLATE_FACTOR = (our edge weight)
// gtree fanout
#define PARTITION_PART 4
// gtree leaf node capacity = tau(in paper)
#define LEAF_CAP 128
// gtree index disk storage
#define FILE_NODES_GTREE_PATH  "../../../data/LV/LV.paths"
#define FILE_GTREE  "../../../data/LV/LV.gtree"
#define FILE_SM  "../../../data/LV/LV.sm"
#define FILE_ONTREE_MIND  "../../../data/LV/LV.minds"
// input
#define FILE_OBJECT  "../../../data/LV/LV.object"
#define FILE_USER  "../../../data/LV/LV.user"
#define FILE_POI  "../../../data/LV/LV.poi"
#define FILE_FRIENDSHIP  "../../../data/LV/LV.friendship"
#define FILE_CHECKIN  "../../../data/LV/LV.checkin"


#define FILE_LOCATION "../../../data/LV/LV.location200"

#define FILE_DISBOUND "../../../data/LV/LV.DISBOUND"
#define FILE_BorderPair "../../../data/LV/LV.borderPair"
#define FILE_MAXSOCIAL "../../../data/LV/LV.MAXSOCIAL"
#define FILE_TEXTUALBOUND "../../../data/LV/LV.TEXTUALBOUND"
#define FILE_QUERY "../../../data/LV/LV.query"


#define FILE_GIMTREE "../../../data/LV/LV.gimtree"
#define FILE_USERSCORE "../../../data/LV/pre-process/LV.userScore"
#define FILE_USERTopKResults "../../../data/LV/pre-process/LV.userTokResults"
#define FILE_POIReverseResults "../../../data/LV/pre-process/LV.poiReverseResults"
#define FILE_LargePOIReverse "../../../data/LV/pre-process/LV.largeReverseResults"
#define FILE_INVLISTOFPOI   "../../../data/LV/pre-process/LV.poiInv"
#define FILE_INVLISTOFUSR   "../../../data/LV/pre-process/LV.usrInv"


/*
 //for NVD
#define FILE_NVD   "../../../data/LV/LV.nvd"
#define FILE_V2P   "../../../data/LV/LV.v2p"
*/


#define FILE_EXP "../exp/results/LV/"
#define FILE_LOG "../exp/results/LV/"
#endif



enum RoadType{HOME,NODE,EDGE,GTREE,MINDS,GPATHS,OBJECTS, PHL,BIN,DISBOUND};
string getRoadInputPath(RoadType roadType){
    stringstream ss;
    string type;
    if(roadType==HOME)
        ss<<DATASET_HOME<<dataset_name<<"/road/";
    else if(roadType==NODE)
        ss<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name <<".cnode";
    else if(roadType==EDGE)
        ss<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name <<".cedge";
    else if(roadType==GTREE)
        ss<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name <<".gtree";
    else if(roadType==MINDS)
        ss<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name <<".minds";
    else if(roadType==GPATHS)
        ss<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name <<".paths";
    else if(roadType==PHL)
        ss<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name <<".phl";
    else if(roadType==BIN)
        ss<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name <<".bin";
    else if(roadType==DISBOUND)
        ss<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name <<".DISBOUND";
    else if(roadType==OBJECTS)
        ss<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name <<".object";
    string path = ss.str();
    return path;

}

enum ObjectType{USER,BUSINESS};

string getObjectInputPath(ObjectType objectType){
    stringstream ss;
    string type;
    if(objectType ==USER)
        type = ".user";
    else if(objectType==BUSINESS)
        type = ".poi";
    ss<<DATASET_HOME<<dataset_name<<"/objects/raw/"<<dataset_name <<type;
    string path = ss.str();
    return path;

}

string getObjectNewInputPath(ObjectType objectType){
    stringstream ss;
    string type;
    if(objectType ==USER)
        type = ".user";
    else if(objectType==BUSINESS)
        type = ".poi";
    ss<<DATASET_HOME<<dataset_name<<"/objects/"<<dataset_name<<type;
    string path = ss.str();
    return path;

}


string getObjectVaryingInputPath(ObjectType objectType,float rate){
    stringstream ss;
    string type;
    if(objectType ==USER)
        type = ".user";
    else if(objectType==BUSINESS)
        type = ".poi";
    ss<<DATASET_HOME<<dataset_name<<"/objects/Varying/"<<dataset_name<<"-"<<(int)(rate*100)<<"%"<<type;
    string path = ss.str();
    return path;

}


string getUserTermUpdatingInputPath(float rate){
    stringstream ss;

    ss<<DATASET_HOME<<dataset_name<<"/objects/Varying/"<<dataset_name<<"-"<<(int)(rate*100)<<"%poi.user";
    string path = ss.str();
    return path;

}




enum SocialType{LINK,CHECKIN};
string getSocialInputPath(SocialType objectType){
    stringstream ss;
    string type;
    if(objectType == LINK)
        type = ".friendship";
    else if(objectType==CHECKIN)
        type = ".checkin";
    ss<<DATASET_HOME<<dataset_name<<"/social/"<<dataset_name <<type;
    string path = ss.str();
    return path;

}

string getSocialRawInputPath(SocialType objectType){
    stringstream ss;
    string type;
    if(objectType == LINK)
        type = ".friendship";
    else if(objectType==CHECKIN)
        type = ".checkin";
    ss<<DATASET_HOME<<dataset_name<<"/source_dataset/social/"<<dataset_name <<type;
    string path = ss.str();
    return path;

}

string getNVDHashIndexInputPath(){
    stringstream ss;
    ss<<DATASET_HOME<<dataset_name<<"/NVD/"<<"V2PNVDHash.idx";
    string path = ss.str();
    return path;
}

string getLeaf2POIListIndexInputPath(){
    stringstream ss;
    ss<<DATASET_HOME<<dataset_name<<"/NVD/"<<"L2PListHash.idx";
    string path = ss.str();
    return path;
}

string getHybridNVDHashIndexInputPath(){
    stringstream ss;
    ss<<DATASET_HOME<<dataset_name<<"/NVD/"<<"HYPNVDHash.idx";
    string path = ss.str();
    return path;
}

string getHybridNVDHashIndexInputPath_varying(float ratio){
    stringstream ss;
    ss<<DATASET_HOME<<dataset_name<<"/NVD/"<<"VaryingO/"<<(int)(ratio*100)<<"%/HYPNVDHash.idx";
    string path = ss.str();
    return path;
}



string getSocialLinkBinaryInputPath(){
    stringstream ss;
    ss<<DATASET_HOME<<dataset_name<<"/social/"<<"SocialLink.bin";
    string path = ss.str();
    return path;
}


string getCheckInBinaryInputPath(){
    stringstream ss;
    ss<<DATASET_HOME<<dataset_name<<"/social/"<<"CheckIn.bin";
    string path = ss.str();
    return path;
}




#define CONFIG_H

#endif //BATCH_RKGSKQ_CONFIG_H