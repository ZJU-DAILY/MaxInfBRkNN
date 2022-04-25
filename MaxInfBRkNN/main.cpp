//
// Created by jins on 5/10/20.
//
#include "process.h"
#include "gim_tree.h"
#include "RkGSKQ.h"
#include "diskbased.h"
#include "Tenindra/command/TransformInputCommand.h"
#include "Tenindra/command/GraphCommand.h"
#include "Tenindra/command/ExperimentsCommand.h"
#include "query_plus.h"

#define NeedDisbound_  1
#define NeedPrepocess_  1



//exp<<average_overall_runtime<<"  "<<average_perquery_runtime<<"  "<<IO<<"\t\t\t";
ofstream exp_results;



//-c ../config_test.txt -s nvd_op  -p mid -n 100 reverse
int main_prepocessing(string& _task) { //_Reverse

    int process_source = 1;

    clock_t buildIndex_startTime,buildIndex_endTime;
    double gtree_time; double phl_time;
    if(process_source){
        ////1.transform (clean data)
        //-m osm -e ../exp/data/osm/USA-road-d.LV.gr -n ../exp/data/osm/USA-road-d.LV.co -o ../exp/data
        if(_task=="clean"){ ///命令行参数 -e pre_process -j clean
            TransformInputCommand* command0 = new TransformInputCommand();
            command0 -> cleanRoadNetworkSourceData();

            ////2. generate binary file (.bin)
            //-e ../exp/data/LV-d.graph  -n ../exp/data/LV-d.coordinates  -o ../exp/indexes/LV.bin  -s ../exp/stats/index_stats.txt
            GraphCommand* command1 = new GraphCommand();   // GraphCommand*
            command1 -> buildBinaryGraph();
        }

        if(_task=="gtree"){ ///命令行参数 -e pre_process -t gtree
            ////3.建立GTree索引
            buildIndex_startTime = clock();
            buildGTreeIndex();
            buildIndex_endTime = clock();
            gtree_time = (double)(buildIndex_endTime-buildIndex_startTime)/CLOCKS_PER_SEC;
            cout<<"**************gtree build time:"<<gtree_time<<"sec!"<<endl;
        }

        if(_task=="phl"){ ///命令行参数 -e pre_process -t phl
            ////4. 建PHL索引：(-e indexes -g ../exp/indexes/LV.bin -p gtree=0,road=0,silc=0,phl=1,ch=0,tnr=0,alt=1,gtree_fanout=4,gtree_maxleafsize=64,road_fanout=4,road_levels=7,silc_maxquadtreeleafsize=1000,tnr_gridsize=128,alt_numlandmarks=16,alt_landmarktype=random
            buildIndex_startTime = clock();
            ExperimentsCommand* command2 = new ExperimentsCommand();
            command2 -> buildingPHL();
            buildIndex_endTime = clock();
            phl_time = (double)(buildIndex_endTime-buildIndex_startTime)/CLOCKS_PER_SEC;
            cout<<"**************phl build time:"<<phl_time<<"sec!"<<endl;


        }


        if(_task=="modify"){
            ////5.进行User poi dataset modify
            clock_t _startTime = clock();
            modifyUPMap();
            clock_t _endTime = clock();
            double modify_time = (double)(_endTime-_startTime)/CLOCKS_PER_SEC;
            cout<<"**************modifyUPMap time:"<<modify_time<<"sec!"<<endl;
        }


    }


    return 1;


}

//-c ../config_test.txt -s group -p mid -n 100 reverse
////-e maxinf -r group -t nvd  -a 0.6 -c ../config_test.txt -p mid -n 60 -v 0.05 -s appro
////-e effective -t g -c ../config_test.txt -k 20 -a 0.6  -p mid -n 60 -b 5 -v 0.05 -s appro -g 1
///　几种基本操作（路网距离，ｒeverse knn )测试的主函数
int main_test(int argc, char* argv[]) { //the main funtion for all the steps in our framework  _coding



    int k=DEFAULT_K; int KQ =DEFAULT_KQ_SIZE; int QN =DEFAULT_Q_SIZE; int b = DEFAULT_BUDGET;
    float alpha = DEFAULT_ALPHA;
    float accuracy = 0.1;

    int alg = DEFAULT_ALG; int popularity = Popular;
    int exp_task = -1;
    int topk_method; int reverse_method;

    string preprocess_task="";
    parameter_settings(argc,argv,k, KQ, QN, b, accuracy,alpha, popularity, alg, exp_task,topk_method,reverse_method,preprocess_task);


//// -e nvd -c ../config_test.txt -s group -p mid -n 100 reverse
    if(exp_task==Pre_Process){ ////1. -e pre_process 对源数据进行清洗处理后输出，并对路网建立索引:

        main_prepocessing(preprocess_task);
        cout<<"main_prepocessing 结束!"<<endl;
        exit(1);
    }
    else if(exp_task==SPD_Test){  ////****2. -e spdist: 进行最短路径距离结果测试:

        test_P2P();
        exit(1);
    }

    else if(exp_task==DisBound){  ////3. -e disbound: 进行node间distance 上下界预计算(为加快baselin):

        TIME_TICK_START
        ConfigType setting;
        AddConfigFromFile(setting,"../map_exp");
        //AddConfigFromCmdLine(setting,argc,argv);
        initGtree();
        // 加载 edge information
        loadEdgeMap();

        initialPHL();

        outputNodeDisBound_Enhance(); //计算距离边界
        //outputNodeDisBound_gtree(); //计算距离边界
        TIME_TICK_END  TIME_TICK_PRINT("outputNodeDisBound by gtree runtime: ")
        exit(1);
    }

    else if(exp_task == TextModify){  ////4. -e text: 对合成数据集进行Zipfan定律下的关键词分布修正:
#ifndef LV_Source
#ifdef LasVegas_Source
        CleanLasVegasTextInfo();
#else
        Zipfan_UserTextInfo();
        Zipfan_POITextInfo();
#endif
#endif
        exit(1);
    }
    else if(exp_task == SocialPROCESS){ //// 5 -e social
#ifdef LasVegas_Source
        ModifyFriendShip();
        ModifyLasVegasUserIDValue();
        //exit(1);
#endif
#ifdef Twitter_Source
        ModifyFriendShip();
        ModifyLasVegasUserIDValue();
        //exit(1);
#endif

        //序列化输出社交网络图信息
        wholeSocialLinkMap.clear(); friendshipMap.clear();followedMap.clear(); friendShipTable.clear();
        serialize_SocialLink();
        //序列化输出check-in信息
        userCheckInIDList.clear(); userCheckInCountList.clear(); userCheckInMap.clear();
        poiCheckInIDList.clear(); poiCheckInCountList.clear(); poiCheckInMap.clear();
        serialize_CheckIn();

    }

    else if(exp_task==GIMTree_Build){  ////5. -e gimtree: 根据当前双色体数据集user, poi+路网划分情况，构建GIM-Tree

        TIME_TICK_START

        processData();

        TIME_TICK_END   TIME_TICK_PRINT(" processData(GIM-Tree) runtime: ")
        exit(1);
    }

    else if(exp_task==Build_NVD){  ////6. -e nvd: 构建 Network Voronoi Diagram 相关的索引信息文件（输出 NVD 文件 以及 v2P 文件）
        //(-e nvd -c ../config_test.txt -s group -p mid -n 100 reverse)
        int need_nvd_cache = 0;
        loadRoadNetworkIndexData();
        loadOnlyPMap();
        initialIOCache(need_nvd_cache);  //初始化Cache
        TIME_TICK_START


#ifdef HybridHash
        //NVD_generation_AllKeyword_V2PHash_Enhenced_gtree();
        NVD_generation_AllKeyword_V2PHash_Enhenced_road();
        TIME_TICK_END  TIME_TICK_PRINT("building (hybrid) NVD runtime: ")
#else
        NVD_generation_AllKeyword_Hash_serial();
        TIME_TICK_END  TIME_TICK_PRINT("building NVD runtime: ")
#endif

        //test_load_V2P_NVDG_AddressIdx();

        exit(1);

    }

    else if(exp_task == TestCode){ //ceshi   ////7. 进行 NVD、V2P部分的测试：-e test

        // *********************** below is for test *************************

        /*test_write_HybirdPOINN_IDX();
        test_load_HybirdPOINN_AddressIdx();*/
        //loadFriendShipBinaryData();
        //loadCheckinBinaryData();

        test_CheckSocialConnection();
        exit(1);



        //loadData();
        int need_nvd_cache = 1;
        loadData();
        initialIOCache(need_nvd_cache);  //初始化Cache
        Load_Hybrid_NVDG_AddressIdx_fast();


#ifdef TESTGIMTree
        //test_lclLeaf();
        //testFunctionForEdge();
        //testTermWeight();
        //testTermInvInfo();
        //testSPSP();
        //getTermSize();
        //testUserTermRelatedEntry_all();
        testObjectTermRelatedEntry_all_self(); ////测试GIM Tree 必须用
        //testAccessPOI();
#endif


        ////测试NVD必须用!
        //test_poiNN_by_Hash_all_round();
        test_NVDAdj_GivenKeyword(1);
        ////测试NVD必须用!
        exit(1);
        test_u2PDist();
        k = 20;
        test_getBorder_SB_NVD(k, DEFAULT_A, alpha);
        //test_u2PDist_single(3413,67);   //o_id, u_id
        int u_id = 53946; int p_id = 161656; int term = 40;
        //test_GSKScoreAndTopk (u_id, p_id,DEFAULT_A,alpha,k);
        int vertex_id = 802; int term_id = 80;
        //test_enumerating_invertedList();
        //test_NVD();
        u_id = 45044;
        test_checkUserWhere(u_id);

        //cout<<SPSP(2510,350)<<endl;
        //test_u2PDist(7813,65);
        //test_u2PDist(64,65);




    }

////
    else {

        termWeight_Map.clear();
        for(int i=0;i<vocabularySize;i++){  //term weight 的cache表
            termWeight_Map.push_back(-1);
        }

        if(exp_task==Topk_Test){  ////8.进行topk的测试：  -e topk -t gtree -c ../config_test.txt -p mid -n 100 reverse
            //-e topk -t gtree|dj|nvd|bruteforce
            int method = topk_method;

            int need_nvd_cache = 1;
            TIME_TICK_START
            loadData();
            TIME_TICK_END   TIME_TICK_PRINT(" loadData total runtime: ")
            initialIOCache(need_nvd_cache);  //初始化Cache


            //test_single_topk(1184962,k, alpha,topk_method); ////先single下测试
            //test_single_topk(625,k, alpha,topk_method); ////先single下测试
            //test_single_topk(970,k, alpha,topk_method); ////先single下测试
            //test_topk(k,alpha,method); ////后 multiple 下测试

            //test_ComparationMethods_topk_single(24,k,alpha);
            //test_GSKScoreAndTopk(1184962, 264, DEFAULT_A,alpha,k);
            //test_GSKScoreAndTopk(625, 57, DEFAULT_A,alpha,k); ////当发现BRkNN结果不对时，进行测试
            //test_topk(k,alpha);
            test_ComparationMethods_topk(k,alpha);
            exit(1);
        }



        else if(exp_task==Reverse_Single){  ////9. 进行reverse 的测试： -e reverse -s group -c ../config_test.txt -p mid -n 100 reverse
            /////主程序运行段落                                    ////(-e reverse -r nvd_op -c ../config_test.txt -p mid -n 100 )
            int method = reverse_method;
            int need_nvd_cache = 0;
            clock_t  load_startTime,load_endTime;
            load_startTime = clock();
            loadData();  //jiazai
            load_socialMatrix();
            load_endTime = clock();
            cout<<"loadData total running time:"<<(double)(load_endTime-load_startTime)/CLOCKS_PER_SEC<<"sec!"<<endl;


            if(method== Group_reverse){  // -s group
                //加载距离边界
                load_startTime = clock();
                loadDisBound();   //这步可以序列化
                load_endTime = clock();
                cout<<"load disBound running time:"<<(double)(load_endTime-load_startTime)/CLOCKS_PER_SEC<<"sec!"<<endl;
            }
            else if(method== NVD_reverse_OP||method== NVD_reverse){  //若是nvd下算法，要加载 v2p_id, nvd_adj  索引文件 //-s nvd or nvd_op

                need_nvd_cache = 1;
                Load_Hybrid_NVDG_AddressIdx_fast();

            }
            else if(method == Reverse_ALL){
                need_nvd_cache = 1;
                TIME_TICK_START
                cout<<"开始load DISBOUND..."<<endl;
                loadDisBound();   //这步可以序列化
                cout<<"DONE! ";
                TIME_TICK_END   TIME_TICK_PRINT()
                Load_Hybrid_NVDG_AddressIdx_fast();
            }

            initialIOCache(need_nvd_cache);  //初始化Cache


#ifdef LV
            int poi = 57;//57;
#endif
#ifdef LasVegas
            int poi = 40; //265;  //265(7)  266(10) 269(19)
#endif
#ifdef Brightkite
            int poi = 161656;//83;
#endif
#ifdef Gowalla
            int poi = 10;
#endif
#ifdef Twitter
            int poi = 10;
#endif

            if(method== Group_reverse){
                //TIME_TICK_START
                RkGSKQ_Processing_single_test(poi,k, DEFAULT_A, alpha, Single,-1,-1); //int Qk, float a, float alpha, int mode, int alg
                //TIME_TICK_END
                //TIME_TICK_PRINT("---------------------RkGSKQ GIM-tree running")
                clearDisBoundArray();
            }
            else if(method== NVD_reverse){
                //TIME_TICK_START
                //RkGSKQ_NVD_Naive_hash(poi,k,DEFAULT_A, alpha);
                RkGSKQ_NVD_Naive_optimized(poi,k,DEFAULT_A, alpha);
                //cout<<"dist computation average runtime: "<< (dist_total_runtime) / dist_cmp_count / CLOCKS_PER_SEC  * 1000000 << "us" << endl;
                //TIME_TICK_END
                //TIME_TICK_PRINT("---------------------RkGSKQ_NVD_Naive_hash running")

            }
            else if(method== NVD_reverse_OP){
                //TIME_TICK_START
                RkGSKQ_bottom2up_by_NVD_noneReturn(poi, k,DEFAULT_A,alpha);
                //cout<<"dist computation average runtime: "<< (dist_total_runtime) / dist_cmp_count / CLOCKS_PER_SEC  * 1000000 << "us" << endl;
                //TIME_TICK_END
                //TIME_TICK_PRINT("---------------------RkGSKQ_NVD_Naive_optimized running")
            }

            else if(method == BruteForce_reverse){
                TIME_TICK_START
                Naive_RkGSKQ(poi, k, DEFAULT_A, alpha);
                TIME_TICK_END
                TIME_TICK_PRINT("RkGSKQ(Naive) running")
            }
            else if(method == Reverse_ALL){
#ifdef TEST_REVERSE
                for(int qo=poi;qo<poi+50;qo++){

                    RkGSKQ_Processing_single_test(qo,k, DEFAULT_A, alpha, Single,-1,-1); //int Qk, float a, float alpha, int mode, int alg

                    //RkGSKQ_bottom2up_by_NVD_noneReturn(qo,k,DEFAULT_A, alpha);
                    cout<<"*********************** query object:p"<<qo<<"处理完毕！*************************"<<endl;
                    //getchar();
                }

#endif
            }


            CloseDiskComm();
            exit(1);

        }

        else if(exp_task== Reverse_Batch  || exp_task == MaxInf){  //pichuli   -e batch -s group -c ../config_test.txt -p mid -n 100

            /////主程序运行段落                        ////(-e batch -s nvd -c ../config_test.txt -p mid -n 100 )

            int need_nvd_cache = 1;
            clock_t  load_startTime,load_endTime;
            load_startTime = clock();
            loadData();  //jiazai
            load_socialMatrix();
            load_endTime = clock();
            cout<<"loadData total running time:"<<(double)(load_endTime-load_startTime)/CLOCKS_PER_SEC<<"sec!"<<endl;

            initialIOCache(need_nvd_cache);  //初始化Cache
            Load_Hybrid_NVDG_AddressIdx_fast();

            // 设定随机数的随机种子
            dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));


            int poi;

#ifdef LasVegas
            poi = 265;
#endif
#ifdef LV
            poi = 57;  //57; ////57(14)  //60( );//155( ) ; //79

#endif

#ifdef Brightkite   //for Brightkite map
            poi = 161656; //161662; //161656; //323200 ;//

#endif

#ifdef Gowalla
            poi =  10;
#endif
#ifdef Twitter
            poi =  10;
#endif
            vector<int> candidateKeywords ;
            vector<int> checkIn_usr;

            vector<int> chainStores;
            //产生查询点的空间位置
            chainStores = generate_query_popular(poi,QN);//generate_query();
            //产生查询点的关键词信息
            candidateKeywords = extractQueryGroupKeywords(chainStores);//extractQueryKeyowrds(poi,DEFAULT_KEYWORDSIZE); //POIs[poi].keywords;//
            stringstream tt;



            if(exp_task == Reverse_Batch){
                tt<<FILE_LOG<<"BatchRkGSKQ_result_"<<dataset_name;
                exp_results.open(tt.str(), ios::app);
                if(alg ==GROUP_BATCH){  // group processing  //// -e base
                    //加载距离边界
                    load_startTime = clock();
                    loadDisBound();   //这步可以序列化
                    load_endTime = clock();
                    cout<<"load disbound running time:"<<(double)(load_endTime-load_startTime)/CLOCKS_PER_SEC<<" sec!"<<endl;

                    MQueryResults results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, chainStores, candidateKeywords, BATCH, GROUP_BATCH, exp_results); //CLUSTER
                    clearDisBoundArray();
                }

                else if (alg == NVD_BATCH){  // efficient processing by nvd  //// -e nvd

                    BatchResults results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, chainStores, candidateKeywords, BATCH, NVD_BATCH, exp_results); //CLUSTER

                }
                else if(alg == SEPRATE){

                    BatchResults results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, chainStores, candidateKeywords, ONE_BY_ONE, NVD_SINGLE, exp_results);
                }

            }

            else if (exp_task==MaxInf){

                tt<<FILE_LOG<<"MaxInfBRGSTkNN_result_"<<dataset_name;
                exp_results.open(tt.str(), ios::app);
                clock_t maxinf_startTime, maxinf_endTime;
                OVERALL_START
                if(alg == BASELINE){ //true){

                    MQueryResults results;
                    if(false){
                        maxinf_startTime = clock();
                        results = sythetic_Multi_BRGSKTkNNResults(chainStores);

                    }
                    else{
                        //加载距离边界
                        load_startTime = clock();
                        loadDisBound();   //这步可以序列化
                        load_endTime = clock();
                        cout<<"load disbound running time:"<<(double)(load_endTime-load_startTime)/CLOCKS_PER_SEC<<" sec!"<<endl;
                        //1.baseline方法
                        maxinf_startTime = clock();
                        results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, chainStores, candidateKeywords, BATCH, GROUP_BATCH, exp_results);
                        //results =  Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, chainStores, candidateKeywords, BATCH, NVD_BATCH, exp_results);
                        //results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, chainStores, candidateKeywords, ONE_BY_ONE, GROUP_SINGLE, exp_results); //CLUSTER
                    }

                    //BatchResults results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, chainStores, candidateKeywords, BATCH, NVD_BATCH, exp_results);
                    //用于社交传播影响力分析
                    InfluenceModels *  L = new InfluenceModels();
                    MaxInfBRGSTkNNQ_Batch_OPIMC(argc, argv, accuracy, online_users, wholeSocialLinkMap,chainStores, results, b, L,exp_results);
                    maxinf_endTime = clock();
                    cout<<"*********************MaxInf by Baseline(Group+OPIMC), runtime="<<(double)(maxinf_endTime-maxinf_startTime)/CLOCKS_PER_SEC*1000  << "ms!" << endl;


#ifdef MC_ESTIMATION
                    TIME_TICK_START
                    L -> influence_esitimation_ByMC(10000,OPC);
                    TIME_TICK_END
                    TIME_TICK_PRINT("MC 10000")
#endif
                    if(L!=NULL) {
                        cout<<"删除L模块..."<<endl;
                        delete L; cout<<"删除L成功！"<<endl;
                    }



                }
                if(alg == APPROXIMATE){ //true){ // -s appro-p pop -n 25 -c ../config_test.txt

                    //用于社交传播影响力分析
                    InfluenceModels *  L = new InfluenceModels();
                    L->doInitial(argc, argv, accuracy, UserID_MaxKey,online_users, wholeSocialLinkMap, chainStores, OPC);

                    maxinf_startTime = clock();
                    BatchResults results;
                    if(false){
                        results = sythetic_Multi_BRGSKTkNNResults(chainStores);

                    }
                    else{
                        results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, chainStores, candidateKeywords, BATCH, NVD_BATCH, exp_results);
                    }


                    MaxInfBRGSTkNNQ_Batch_HYBRID(argc, argv, accuracy, online_users, wholeSocialLinkMap,chainStores, results, b, L,exp_results);
                    maxinf_endTime = clock();
                    cout<<"*********************MaxInf by Appro(NVD+Hybrid), runtime="<<(double)(maxinf_endTime-maxinf_startTime)/CLOCKS_PER_SEC*1000  << "ms!" << endl;


#ifdef MC_ESTIMATION
                    TIME_TICK_START
                    L -> influence_esitimation_ByMC(10000,OPC);
                    TIME_TICK_END
                    TIME_TICK_PRINT("MC 10000")
#endif
                    if(L!=NULL) {
                        cout<<"删除L模块..."<<endl;
                        delete L; cout<<"删除L成功！"<<endl;
                    }


                }

                if(alg == HEURISTIC){

                    InfluenceModels *  L = new InfluenceModels();
                    maxinf_startTime = clock();
                    BatchFastFilterResults fastFilter_output= Fast_Batch_NVDFilter(chainStores, k,DEFAULT_A,alpha);
                    //printBatchRkGSKQResults(chainStores, fastFilter_output.batch_results);

                    if(true){
                        float quality_tradoff = 0.5;
                        vector<int> poiSeeds = MaxInfBRGSTkNNQ_POISelect_Heuristic(argc, argv, k,alpha, quality_tradoff, online_users, wholeSocialLinkMap, chainStores, fastFilter_output, b, L,exp_results);
                        maxinf_endTime = clock();
                        cout<<"*********************MaxInf by heuristic, runtime="<<(double)(maxinf_endTime-maxinf_startTime)/CLOCKS_PER_SEC*1000  << "ms!" << endl;
                        ////验证poiSeeds的影响力结果质量
                        set<int> userSeedSet;
                        vector<int> poiSeedsTextUnion = extractQueryGroupKeywords(poiSeeds);//extractQueryKeyowrds(poi,DEFAULT_KEYWORDSIZE); //POIs[poi].keywords;//
                        BatchResults results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, poiSeeds, poiSeedsTextUnion, BATCH, NVD_BATCH, exp_results);
                        for(int poi:poiSeeds){
                            for(ResultDetail r: results[poi]){
                                int u_id = r.usr_id;
                                userSeedSet.insert(r.usr_id);
                            }
                        }////获得poi seed中所有兴趣点的潜在用户并集
                        ////蒙特卡洛采样以评估给定用户集的social influence
                        double poi_inf = L-> opp-> MCSimulation4GivenSeeds(10000, userSeedSet);
                        printf("poi seed的 influence set为：\n");
                        printSeedSet(userSeedSet);
                    }

                }


                if(alg == ONEHOP){ //true){ // -s appro-p pop -n 25 -c ../config_test.txt

                    //用于社交传播影响力分析
                    InfluenceModels *  L = new InfluenceModels();

                    maxinf_startTime = clock();
                    BatchResults results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, chainStores, candidateKeywords, BATCH, NVD_BATCH, exp_results);

                    MaxInfBRGSTkNNQ_Batch_HopHeuristic(argc, argv, accuracy, online_users, wholeSocialLinkMap,chainStores, results, b, L,exp_results);

                    maxinf_endTime = clock();
                    cout<<"*********************MaxInf by HOPHeuris(NVD+OneHop), runtime="<<(double)(maxinf_endTime-maxinf_startTime)/CLOCKS_PER_SEC*1000  << "ms!" << endl;


#ifdef MC_ESTIMATION
                    TIME_TICK_START
                    L -> influence_esitimation_ByMC(10000,OPC);
                    TIME_TICK_END
                    TIME_TICK_PRINT("MC 10000")
#endif
                    if(L!=NULL) {
                        cout<<"删除L模块..."<<endl;
                        delete L; cout<<"删除L成功！"<<endl;
                    }


                }




                //cout<<"have "<<chainStores.size()<<" chain stores for BRGSTkNN Quering"<<endl;



            }

            CloseDiskComm();


        }

        else if(exp_task== Pre_MC){  ////9. 进行 蒙特卡洛预采样：-e premc
            /////主程序运行段落                                    ////(-e reverse -r nvd_op -c ../config_test.txt -p mid -n 100 )

            // 设定随机数的随机种子
            dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
            TIME_TICK_START
            //读取社交网络图信息
            loadFriendShipData();  //jiazai
            TIME_TICK_END   TIME_TICK_PRINT(" loadFriendShip runtime: ")
            InfluenceModels *  L = new InfluenceModels();
            L-> doInitial_PreSampling_MC(argc, argv, accuracy,online_users, wholeSocialLinkMap);
            int sampleSize = 50;
            L-> preSampling_MC(sampleSize);

        }

        exp_results.close();


    }


}


////几种简单的poi selection 方法
vector<int> poiSelectionRkNNCardinalityFirst(vector<int>& stores,int k, int b,float alpha,ofstream& exp_results){
    //统计每个兴趣点与之文本相关的用户个数，并进行排序
    priority_queue<CardinarlityEntry> CardinalityFirst_Q;
    vector<int> poiSTextUnion = extractQueryGroupKeywords(stores);//extractQueryKeyowrds(poi,DEFAULT_KEYWORDSIZE); //POIs[poi].keywords;//
    BatchResults results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, stores, poiSTextUnion, BATCH, NVD_BATCH, exp_results);
    for(int store_id: stores){
        int _size = results[store_id].size();
        CardinarlityEntry _entry(store_id,_size);
        CardinalityFirst_Q.push(_entry);
    }
    vector<int> _poi_select;
    while(_poi_select.size()<b){
        CardinarlityEntry _entry = CardinalityFirst_Q.top();
        CardinalityFirst_Q.pop();
        int _id = _entry.poi_id;
        int _score = _entry.cardinality;
        _poi_select.push_back(_id);
        cout<<"选择p"<<_id<<", score="<<_score<<endl;
    }
    return _poi_select;
}


#define MillionSecond_ 1

////-e performance  -t k  -a 0.6 -c ../config_test.txt -p mid -n 60 -v 0.05 -s appro
//进行实验测试时的命令！
int main(int argc, char* argv[]) { // the main function for exp testing

    int k=DEFAULT_K; int KQ =DEFAULT_KQ_SIZE; int QN =DEFAULT_Q_SIZE; int b = DEFAULT_BUDGET;
    float alpha = DEFAULT_ALPHA;
    float accuracy = 0.1;

    int alg = DEFAULT_ALG; int popularity = Popular; float o_ratio = 1.0;
    int exp_task = -1; int target_parameter = -1; int group_id = 1; string preprocess_task="";

        //int argc, char* argv[], int& k, int& KQ, int& QN, int& b, float& accuracy, float& alpha, int& popularity, int& alg,int& exp_task, int& test_target
    performanceParameter_settings(argc,argv,k, KQ, QN, b, accuracy,alpha, o_ratio, popularity, alg, exp_task,target_parameter,group_id,preprocess_task);

    if(exp_task==SelectPOI){  // -e select 为数据集选择 poi 候选点
        int need_nvd_cache = 1;
        TIME_TICK_START
        loadData();
        TIME_TICK_END   TIME_TICK_PRINT(" loadData total runtime: ")
        initialIOCache(need_nvd_cache);  //初始化Cache
        poi_selection_heuristic(k,alpha);
        return  0;

    }

//// -e nvd -c ../config_test.txt -s group -p mfnvd_oid -n 100 reverse

    if(exp_task== Performance || exp_task== Effectiveness) {
        /////主程序运行前的初始程序                        ////(-e batch -s nvd -c ../config_test.txt -p mid -n 100 )
        termWeight_Map.clear();
        for(int i=0;i<vocabularySize;i++){  //term weight 的cache表
            termWeight_Map.push_back(-1);
        }

        int need_nvd_cache = 1;
        clock_t  load_startTime,load_endTime;
        load_startTime = clock();
        loadData_exp();  //jiazai
        load_socialMatrix();
        load_endTime = clock();
        cout<<"loadData total running time:"<<(double)(load_endTime-load_startTime)/CLOCKS_PER_SEC<<"sec!"<<endl;

        initialIOCache(need_nvd_cache);  //初始化Cache
        Load_Hybrid_NVDG_AddressIdx_fast();

        // 设定随机数的随机种子
        dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));

        int poi;

#ifdef LasVegas
        poi = 265;
#endif
#ifdef LV
        poi = 57;  //57; ////57(14)  //60( );//155( ) ; //79

#endif

#ifdef Brightkite   //for Brightkite map
        poi = 161656; //161662; //161656; //323200 ;//

#endif

#ifdef Gowalla
        poi =  10;
#endif
#ifdef Twitter
        poi =  10;
#endif
        clock_t maxinf_startTime, maxinf_endTime;
        double accumulation_time = 0;
        double accumulation_inf = 0; double accumulation_rrset = 0;
        double batch_time = 0; double POI_selectTime = 0;
        stringstream tt;  stringstream tt2; string outputPathPre;
        string target_info;
        vector<int> candidateKeywords;
        vector<int> checkIn_usr;
        vector<int> chainStores;

        if(exp_task== Performance) {  //pichuli   -e batch -s group -c ../config_test.txt -p mid -n 100


            if(false){
                //产生查询点的空间位置
                chainStores = generate_query_popular(poi, QN);//generate_query();
                //产生查询点的关键词信息
                candidateKeywords = extractQueryGroupKeywords(
                        chainStores);//extractQueryKeyowrds(poi,DEFAULT_KEYWORDSIZE); //POIs[poi].keywords;//
            }
            else{
#ifdef LasVegas
                chainStores = generate_query_PopularDominantNonEmpty(group_id, QN);
                //chainStores = generate_query_PopularDominant(group_id, QN);
#else
                chainStores = generate_query_PopularDominant(group_id, QN);
#endif
                //产生查询点的关键词信息
                candidateKeywords = extractQueryGroupKeywords(
                        chainStores);
            }


            if (target_parameter == QK) {
                tt2 << "k=" << k;
                target_info = tt2.str();
                tt << FILE_LOG << "MaxInf_" << dataset_name<<".k";
            }
            else if (target_parameter == PC) {
                tt2 << "|Pc|=" << QN;
                target_info = tt2.str();
                tt << FILE_LOG << "MaxInf_" << dataset_name<<".|Pc|";
            }
            else if (target_parameter == Budget) {
                tt2 << "budget=" << b;
                target_info = tt2.str();
                tt << FILE_LOG << "MaxInf_" << dataset_name<<".budget";
            }
            else if(target_parameter == ALPHA){
                tt2 << "alpha =" << alpha;
                target_info = tt2.str();
                tt << FILE_LOG << "MaxInf_" << dataset_name<<".alpha";
            }
            else if(target_parameter == URatio){
                tt2 << "u_ratio =" << o_ratio;
                target_info = tt2.str();
                tt << FILE_LOG << "MaxInf_" << dataset_name<<".u_raio";
            }

            else if(target_parameter == ORatio){
                tt2 << "o_ratio =" << o_ratio;
                target_info = tt2.str();
                tt << FILE_LOG << "MaxInf_" << dataset_name<<".o_raio";
            }


        }
        else if(exp_task== Effectiveness){  ////-e effective -g 1   -a 0.6 -c ../config_test.txt -p mid -n 60 -v 0.1 -s appro
#ifdef LasVegas
            chainStores = generate_query_PopularDominantNonEmpty(group_id, QN);
            //chainStores = generate_query_PopularDominantLasVegas(group_id, QN);
            //chainStores = generate_query_PopularDominant(group_id, QN);
#else
            chainStores = generate_query_PopularDominant(group_id, QN);
#endif

            //产生查询点的关键词信息
            candidateKeywords = extractQueryGroupKeywords(
                    chainStores);

            if (target_parameter == QUERYSET) {
                tt2 << "query set =" << group_id<<", |Ps|="<<b;
                target_info = tt2.str();
                tt << FILE_LOG << "MaxInf_" << dataset_name<<".effectiveness";
            }

        }

        outputPathPre = tt.str();
        tt.str("");
        OVERALL_START
        if(alg == BASELINE){ //true){
            tt<<outputPathPre<<"-baseline";
            exp_results.open(tt.str(), ios::app);
            MQueryResults results;
            //加载距离边界
            load_startTime = clock();
            loadDisBound();   //这步可以序列化
            load_endTime = clock();
            cout<<"load disbound running time:"<<(double)(load_endTime-load_startTime)/CLOCKS_PER_SEC<<" sec!"<<endl;
            //用于社交传播影响力分析
            InfluenceModels *  L = new InfluenceModels();
            L->doInitial(argc, argv, accuracy, UserID_MaxKey,online_users, wholeSocialLinkMap, chainStores, OPC);    //初始化社交影响力评估模块
            batch_time = 0; POI_selectTime = 0;
            for(int i =0;i<repeat_baseline;i++){
                //1.baseline方法
                exp_results <<"--------"<<i<<"th testing: "<<target_info<<" of  baseline solution begins!--------"<< endl;
                maxinf_startTime = clock();
                clock_t batch_startTime = clock();
                results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, chainStores, candidateKeywords, BATCH, GROUP_BATCH, exp_results);
                clock_t batch_endTime = clock();
                batch_time += (double)(batch_endTime-batch_startTime)/CLOCKS_PER_SEC*1000;

                clock_t POI_startTime = clock();
                int rrset_size = MaxInfBRGSTkNNQ_Batch_OPIMC(argc, argv, accuracy, online_users, wholeSocialLinkMap,chainStores, results, b, L, exp_results);
                clock_t POI_endTime = clock();
                accumulation_rrset +=rrset_size;
                POI_selectTime+= (double)(POI_endTime-POI_startTime)/CLOCKS_PER_SEC*1000;
                maxinf_endTime = clock();
#ifdef MillionSecond
                double current_runtime = (double)(maxinf_endTime-maxinf_startTime)/CLOCKS_PER_SEC*1000;
                cout<<"*********************MaxInf by Baseline(Group+OPIMC), runtime="<< current_runtime << "ms!" << endl;
                exp_results <<"--------"<<i<<"th testing: "<<target_info<<"runtime of  baseline solution:"<< current_runtime << " ms!";
#else
                double current_runtime = (double)(maxinf_endTime-maxinf_startTime)/CLOCKS_PER_SEC;
                cout<<"*********************MaxInf by Baseline(Group+OPIMC), runtime="<< current_runtime << "sec!" << endl;
                exp_results <<"--------"<<i<<"th testing: "<<target_info<<", runtime of  baseline solution:"<< current_runtime << " sec!--------";
#endif
                accumulation_time += current_runtime;



#ifdef MC_ESTIMATION
                TIME_TICK_START
                double current_inf = L -> influence_esitimation_ByMC(10000,OPC);
                accumulation_inf += current_inf;
                exp_results<<"current_inf:"<< current_inf <<endl;
                TIME_TICK_END
                TIME_TICK_PRINT("MC 10000")
#endif
            }
            if(L!=NULL) {
                cout<<"删除L模块..."<<endl;
                delete L; cout<<"删除L成功！"<<endl;
            }

#ifdef MillionSecond
            exp_results <<"***********"<<target_info<<", average runtime baseline solution:"<< accumulation_time/repeat_baseline <<"ms, average rrset size="<<accumulation_rrset/repeat_baseline<<",average inf="<<(accumulation_inf/repeat_baseline) << "*********************"<<endl;
#else
            exp_results <<"***********"<<target_info<<", average total runtime of baseline:"<< accumulation_time/repeat_baseline <<"sec, average rrset size="<<accumulation_rrset/repeat_baseline<<",average inf="<<(accumulation_inf/repeat_baseline) << "*********************"<<endl;
            exp_results <<"average batch time:"<< batch_time /repeat_baseline <<"ms, average POI selection time:"<<POI_selectTime /repeat_baseline<<"ms!"<<"*********************"<<endl;
            batch_time = 0;POI_selectTime=0;
#endif


        }
        else if(alg == APPROXIMATE){ //true){ // -s appro-p pop -n 25 -c ../config_test.txt

            tt<<outputPathPre<<"-appro";
            exp_results.open(tt.str(), ios::app);
            //用于社交传播影响力分析
            InfluenceModels *  L = new InfluenceModels();
            L->doInitial(argc, argv, accuracy, UserID_MaxKey,online_users, wholeSocialLinkMap, chainStores, OPC);    //初始化社交影响力评估模块

            for(int i =0;i<repeat;i++){
                exp_results <<"--------"<<i<<"th testing :"<<target_info<<", runtime of  approximate solution begins!--------"<<endl;
                maxinf_startTime = clock();
                clock_t batch_startTime = clock();
                BatchResults results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, chainStores, candidateKeywords, BATCH, NVD_BATCH, exp_results);
                clock_t batch_endTime = clock();
                batch_time += (double)(batch_endTime-batch_startTime)/CLOCKS_PER_SEC*1000;

                clock_t POI_startTime = clock();
                int rrset_size = MaxInfBRGSTkNNQ_Batch_HYBRID(argc, argv, accuracy, online_users, wholeSocialLinkMap,chainStores, results, b, L, exp_results);
                clock_t POI_endTime = clock();
                POI_selectTime+= (double)(POI_endTime-POI_startTime)/CLOCKS_PER_SEC*1000;
                maxinf_endTime = clock();

                double current_runtime = (double)(maxinf_endTime-maxinf_startTime)/CLOCKS_PER_SEC*1000;
                accumulation_time += current_runtime;
                accumulation_rrset+= rrset_size;

                cout<<"*********************MaxInf by Appro(NVD+Hybrid), runtime="<<current_runtime << "ms!*********************" << endl;

#ifdef MC_ESTIMATION
                TIME_TICK_START
                //L -> influence_esitimation_ByMC(10000,OPC);
                double current_inf = L -> influence_esitimation_ByMC(10000,OPC);
                accumulation_inf += current_inf; //exp_results<<"current_inf:"<< current_inf <<endl;
                TIME_TICK_END
                TIME_TICK_PRINT("MC 10000")
#endif
                exp_results <<"--------"<<i<<"th testing :"<<target_info<<", runtime of  approximate solution:"<< current_runtime << "ms!"<<"current_inf:"<< current_inf <<"--------"<<endl;

            }

            if(L!=NULL) {
                cout<<"删除L模块..."<<endl;
                delete L; cout<<"删除L成功！"<<endl;
            }

            exp_results <<"***********"<<target_info<<",average runtime of approximate solution:"<< accumulation_time /repeat << "ms, average rrset size"<<accumulation_rrset/repeat<<", average inf="<<(accumulation_inf/repeat) << "*********************"<<endl;
            exp_results <<"average batch time:"<< batch_time /repeat <<"ms, average POI selection time:"<<POI_selectTime /repeat<<"ms!"<<"*********************"<<endl;
            batch_time = 0;POI_selectTime=0;


        }

        else if(alg == HEURISTIC){
            tt<<outputPathPre<<"-heuristic";
            string _tmp = tt.str();
            exp_results.open(tt.str(), ios::app);
            InfluenceModels *  L = new InfluenceModels();
            //printBatchRkGSKQResults(chainStores, fastFilter_output.batch_results);
            ////执行heuristic solution，输出poi seeds
            vector<int> poiSeeds; accumulation_inf=0;
            for(int i =0;i<repeat;i++){
                float quality_tradoff = 0.5;
                maxinf_startTime = clock();
                BatchFastFilterResults fastFilter_output= Fast_Batch_NVDFilter(chainStores, k,DEFAULT_A,alpha);
                poiSeeds = MaxInfBRGSTkNNQ_POISelect_Heuristic(argc, argv, k,alpha, quality_tradoff, online_users, wholeSocialLinkMap, chainStores, fastFilter_output, b, L,exp_results);
                maxinf_endTime = clock();
                double current_runtime = (double)(maxinf_endTime-maxinf_startTime)/CLOCKS_PER_SEC*1000;
                accumulation_time += current_runtime;
                cout<<"*********************MaxInf by heuristic, runtime="<< current_runtime  << "ms!" << endl;
                exp_results <<"-------- "<<i<<"th testing ("<<target_info<<"):"<<"runtime of  heuristic solution: "<< current_runtime << "ms" << endl;
                cout<<"heuristic solution 所获得的poi seeds:"<<endl;
                for(int poi:poiSeeds){
                    cout<<"p"<<poi<<",";
                    exp_results<<"p"<<poi<<",";
                }
                cout<<endl; exp_results<<endl;
                ////验证poiSeeds的结果质量,体现在其影响力的大小上
                set<int> userSeedSet;
                vector<int> poiSeedsTextUnion = extractQueryGroupKeywords(poiSeeds);//extractQueryKeyowrds(poi,DEFAULT_KEYWORDSIZE); //POIs[poi].keywords;//
                BatchResults results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, poiSeeds, poiSeedsTextUnion, BATCH, NVD_BATCH, exp_results);
                for(int poi:poiSeeds){
                    for(ResultDetail r: results[poi]){
                        int u_id = r.usr_id;
                        userSeedSet.insert(r.usr_id);
                    }
                }////获得poi seed中所有兴趣点的潜在用户并集
                ////蒙特卡洛采样以评估给定用户集的social influence
                exp_results<<"heuristic 获得的 poi set 的 influence set size="<<userSeedSet.size()<<endl;
                double poi_inf = L-> opp-> MCSimulation4GivenSeeds(10000, userSeedSet);
                printf("heuristic 获得的 poi seed的 influence set为：\n");
                printSeedSet(userSeedSet);
                accumulation_inf += poi_inf;
            }

            exp_results <<target_info<<"***********average runtime of heuristic solution:"<< accumulation_time/repeat << "ms，" <<",heuristic solution 所获得的 poi seed 的 influence:"<< (accumulation_inf/repeat) << endl;

        }

        else if(alg == RELEVANCE){
            tt<<outputPathPre<<"-relevance";
            string _tmp = tt.str();
            exp_results.open(tt.str(), ios::app);
            InfluenceModels *  L = new InfluenceModels();
            L->doInitial(argc, argv, accuracy, UserID_MaxKey,online_users, wholeSocialLinkMap, chainStores, OPC);    //初始化社交影响力评估模块
            //printBatchRkGSKQResults(chainStores, fastFilter_output.batch_results);
            ////执行relevance solution，输出poi seeds
            vector<int> poiSeeds;
            for(int i =0;i<repeat;i++){

                poiSeeds = MaxInfBRGSTkNNQ_POISelect_Relevance(argc, argv, chainStores, b, exp_results);
                maxinf_endTime = clock();
                double current_runtime = (double)(maxinf_endTime-maxinf_startTime)/CLOCKS_PER_SEC*1000;
                accumulation_time += current_runtime;
                cout<<"*********************MaxInf by heuristic, runtime="<< current_runtime  << "ms!" << endl;
                exp_results <<"-------- "<<i<<"th testing ("<<target_info<<"):"<<"runtime of  relevance solution: "<< current_runtime << "ms" << endl;
                cout<<"relevance solution 所获得的poi seeds:"<<endl;
                for(int poi:poiSeeds){
                    cout<<"p"<<poi<<",";
                    exp_results<<"p"<<poi<<",";
                }
                cout<<endl; exp_results<<endl;
                ////验证poiSeeds的结果质量,体现在其影响力的大小上
                set<int> userSeedSet;
                vector<int> poiSeedsTextUnion = extractQueryGroupKeywords(poiSeeds);//extractQueryKeyowrds(poi,DEFAULT_KEYWORDSIZE); //POIs[poi].keywords;//
                BatchResults results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, poiSeeds, poiSeedsTextUnion, BATCH, NVD_BATCH, exp_results);
                for(int poi:poiSeeds){
                    for(ResultDetail r: results[poi]){
                        int u_id = r.usr_id;
                        userSeedSet.insert(r.usr_id);
                    }
                }////获得poi seed中所有兴趣点的潜在用户并集
                ////蒙特卡洛采样以评估给定用户集的social influence
                exp_results<<"Relevance 获得的 poi set 的 influence set size="<<userSeedSet.size()<<endl;
                double poi_inf = L-> opp-> MCSimulation4GivenSeeds(10000, userSeedSet);
                printf("Relevance 获得的 poi seed的 influence set为：\n");
                printSeedSet(userSeedSet);
                accumulation_inf += poi_inf;
            }

            exp_results <<target_info<<"***********average runtime of relevance solution:"<< (accumulation_time/repeat) << "ms，" <<",relevance solution 所获得的 poi seed 的 influence:"<< (accumulation_inf/repeat) << endl;

        }

        else if(alg == CARDINALITY){
            tt<<outputPathPre<<"-cardinality";
            string _tmp = tt.str();
            exp_results.open(tt.str(), ios::app);
            InfluenceModels *  L = new InfluenceModels();
            L->doInitial(argc, argv, accuracy, UserID_MaxKey,online_users, wholeSocialLinkMap, chainStores, OPC);    //初始化社交影响力评估模块
            //printBatchRkGSKQResults(chainStores, fastFilter_output.batch_results);
            ////执行relevance solution，输出poi seeds
            vector<int> poiSeeds;
            for(int i =0;i<repeat;i++){

                maxinf_startTime = clock();
                poiSeeds = poiSelectionRkNNCardinalityFirst(chainStores,k, b, alpha, exp_results);
                maxinf_endTime = clock();
                double current_runtime = (double)(maxinf_endTime-maxinf_startTime)/CLOCKS_PER_SEC*1000;
                accumulation_time += current_runtime;
                cout<<"*********************MaxInf by cardinality, runtime="<< current_runtime  << "ms!" << endl;
                exp_results <<"-------- "<<i<<"th testing ("<<target_info<<"):"<<"runtime of  heuristic solution: "<< current_runtime << "ms" << endl;
                cout<<"relevance solution 所获得的poi seeds:"<<endl;
                for(int poi:poiSeeds){
                    cout<<"p"<<poi<<",";
                    exp_results<<"p"<<poi<<",";
                }
                cout<<endl; exp_results<<endl;
                ////验证poiSeeds的结果质量,体现在其影响力的大小上
                set<int> userSeedSet;
                vector<int> poiSeedsTextUnion = extractQueryGroupKeywords(poiSeeds);//extractQueryKeyowrds(poi,DEFAULT_KEYWORDSIZE); //POIs[poi].keywords;//
                BatchResults results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, poiSeeds, poiSeedsTextUnion, BATCH, NVD_BATCH, exp_results);
                for(int poi:poiSeeds){
                    for(ResultDetail r: results[poi]){
                        int u_id = r.usr_id;
                        userSeedSet.insert(r.usr_id);
                    }
                }////获得poi seed中所有兴趣点的潜在用户并集
                ////蒙特卡洛采样以评估给定用户集的social influence
                exp_results<<"Cardinality获得的poi 获得的 poi set 的 influence set size="<<userSeedSet.size()<<endl;
                double poi_inf = L-> opp-> MCSimulation4GivenSeeds(10000, userSeedSet);
                printf("cardinality获得的poi seed的 influence set为：\n");
                printSeedSet(userSeedSet);
                accumulation_inf += poi_inf;
            }

            exp_results <<target_info<<"***********average runtime of cardinality solution:"<< (accumulation_time/repeat) << "ms，" <<",cardinality solution 所获得的 poi seed 的 influence:"<< (accumulation_inf/repeat) << endl;

        }

        else if(alg == INFLUENCER){
            tt<<outputPathPre<<"-influencer";
            string _tmp = tt.str();
            exp_results.open(tt.str(), ios::app);
            InfluenceModels *  L = new InfluenceModels();
            //printBatchRkGSKQResults(chainStores, fastFilter_output.batch_results);
            ////执行heuristic solution，输出poi seeds
            vector<int> poiSeeds; accumulation_inf=0;
            for(int i =0;i<repeat;i++){
                float quality_tradoff = 0.5;
                maxinf_startTime = clock();
                BatchFastFilterResults fastFilter_output= Fast_Batch_NVDFilter(chainStores, k,DEFAULT_A,alpha);
                poiSeeds = MaxInfBRGSTkNNQ_POISelect_Influencer(argc, argv, k,alpha, quality_tradoff, online_users, wholeSocialLinkMap, chainStores, fastFilter_output, b, L,exp_results);
                maxinf_endTime = clock();
                double current_runtime = (double)(maxinf_endTime-maxinf_startTime)/CLOCKS_PER_SEC*1000;
                accumulation_time += current_runtime;
                cout<<"*********************MaxInf by Influencer first , runtime="<< current_runtime  << "ms!" << endl;
                exp_results <<"-------- "<<i<<"th testing ("<<target_info<<"):"<<"runtime of  influencer solution: "<< current_runtime << "ms" << endl;
                cout<<"influencer first solution 所获得的poi seeds:"<<endl;
                for(int poi:poiSeeds){
                    cout<<"p"<<poi<<",";
                    exp_results<<"p"<<poi<<",";
                }
                cout<<endl; exp_results<<endl;
                ////验证poiSeeds的结果质量,体现在其影响力的大小上
                set<int> userSeedSet;
                vector<int> poiSeedsTextUnion = extractQueryGroupKeywords(poiSeeds);//extractQueryKeyowrds(poi,DEFAULT_KEYWORDSIZE); //POIs[poi].keywords;//
                BatchResults results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, poiSeeds, poiSeedsTextUnion, BATCH, NVD_BATCH, exp_results);
                for(int poi:poiSeeds){
                    for(ResultDetail r: results[poi]){
                        int u_id = r.usr_id;
                        userSeedSet.insert(r.usr_id);
                    }
                }////获得poi seed中所有兴趣点的潜在用户并集
                ////蒙特卡洛采样以评估给定用户集的social influence
                exp_results<<"Influencer 获得的 poi set 的 influence set size="<<userSeedSet.size()<<endl;double poi_inf = L-> opp-> MCSimulation4GivenSeeds(10000, userSeedSet);
                printf("Influencer 获得的poi seed的 influence set为：\n");
                printSeedSet(userSeedSet);
                accumulation_inf += poi_inf;
            }

            exp_results <<target_info<<"***********average runtime of influencer solution:"<< accumulation_time/repeat << "ms，" <<",influencer solution 所获得的 poi seed 的 influence:"<< (accumulation_inf/repeat) << endl;

        }

        else if(alg == RANDOM){
            tt<<outputPathPre<<"-random";
            string _tmp = tt.str();
            exp_results.open(tt.str(), ios::app);
            InfluenceModels *  L = new InfluenceModels();
            L->doInitial(argc, argv, accuracy, UserID_MaxKey,online_users, wholeSocialLinkMap, chainStores, OPC);    //初始化社交影响力评估模块
            ////执行heuristic solution，输出poi seeds
            vector<int> poiSeeds; accumulation_inf=0;
            for(int i =0;i<repeat_random;i++){
                float quality_tradoff = 0.5;
                maxinf_startTime = clock();
                poiSeeds = MaxInfBRGSTkNNQ_POISelect_Random(argc, argv, k,alpha, quality_tradoff, online_users, wholeSocialLinkMap, chainStores, b, L,exp_results);
                maxinf_endTime = clock();
                double current_runtime = (double)(maxinf_endTime-maxinf_startTime)/CLOCKS_PER_SEC*1000;
                accumulation_time += current_runtime;
                cout<<"*********************MaxInf by Random, runtime="<< current_runtime  << "ms!" << endl;
                exp_results <<"-------- "<<i<<"th testing ("<<target_info<<"):"<<"runtime of  random solution: "<< current_runtime << "ms" << endl;
                cout<<"Random selection所获得的poi seeds:"<<endl;
                for(int poi:poiSeeds){
                    cout<<"p"<<poi<<",";
                    exp_results<<"p"<<poi<<",";
                }
                cout<<endl; exp_results<<endl;
                ////验证poiSeeds的结果质量,体现在其影响力的大小上
                set<int> userSeedSet;
                vector<int> poiSeedsTextUnion = extractQueryGroupKeywords(poiSeeds);//extractQueryKeyowrds(poi,DEFAULT_KEYWORDSIZE); //POIs[poi].keywords;//
                BatchResults results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, poiSeeds, poiSeedsTextUnion, BATCH, NVD_BATCH, exp_results);
                for(int poi:poiSeeds){
                    for(ResultDetail r: results[poi]){
                        int u_id = r.usr_id;
                        userSeedSet.insert(r.usr_id);
                    }
                }////获得poi seed中所有兴趣点的潜在用户并集
                ////蒙特卡洛采样以评估给定用户集的social influence
                exp_results<<"RANDOM 获得的 poi set 的 influence set size="<<userSeedSet.size()<<endl;
                double poi_inf = L-> opp-> MCSimulation4GivenSeeds(10000, userSeedSet);
                printf("RANDOM 获得的poi seed的 influence set为：\n");
                printSeedSet(userSeedSet);
                accumulation_inf += poi_inf;
            }

            exp_results <<target_info<<"***********average runtime of Random select solution:"<< accumulation_time/repeat_random << "ms，" <<",random solution 所获得的 poi seed 的 influence:"<< (accumulation_inf/repeat_random) << endl;

        }

        CloseDiskComm();
        exp_results.close();
        return 1;

    }


    else if(exp_task==Pre_Process){ ////1. -e pre_process 对源数据进行清洗处理后输出，并对路网建立索引:

        string _task="";
        preprocessParameter_settings(argc,argv,_task);
        main_prepocessing(_task);
        cout<<"main_prepocessing 结束!"<<endl;
        exit(1);
    }
    else if(exp_task==SPD_Test){  ////****2. -e spdist: 进行最短路径距离结果测试:

        test_P2P();
        exit(1);
    }

    else if(exp_task==DisBound){  ////3. -e disbound: 进行node间distance 上下界预计算(为加快baselin):

        TIME_TICK_START
        ConfigType setting;
        AddConfigFromFile(setting,"../map_exp");
        //AddConfigFromCmdLine(setting,argc,argv);
        initGtree();
        // 加载 edge information
        loadEdgeMap(); initialPHL();

        outputNodeDisBound_Enhance(); //计算距离边界
        //outputNodeDisBound_gtree(); //计算距离边界
        TIME_TICK_END  TIME_TICK_PRINT("outputNodeDisBound by gtree runtime: ")
        exit(1);
    }

    else if(exp_task == TextModify){  ////4. -e text: 对合成数据集进行Zipfan定律下的关键词分布修正:
#ifndef LV_Source
#ifdef LasVegas_Source
        CleanLasVegasTextInfo();
#else
        Zipfan_UserTextInfo();
        Zipfan_POITextInfo();
#endif
#endif
        exit(1);
    }

    else if(exp_task == VaryingObject){  ////4. -e varying : 改变双色体对象（P or U）集合的个数:

        Varying_POIInfo();
        Varying_UserInfo();
        Varying_buildFold();
        exit(1);
    }

    else if(exp_task == SocialPROCESS){ ////5. -e social
#ifdef LasVegas_Source
        ModifyFriendShip();
        ModifyLasVegasUserIDValue();  //进一步更新social link 以及 check-in
        //exit(1);
#endif
        //序列化输出社交网络图信息
        wholeSocialLinkMap.clear(); friendshipMap.clear();followedMap.clear(); friendShipTable.clear();
        serialize_SocialLink();
        //序列化输出check-in信息
        userCheckInIDList.clear(); userCheckInCountList.clear(); userCheckInMap.clear();
        poiCheckInIDList.clear(); poiCheckInCountList.clear(); poiCheckInMap.clear();
        serialize_CheckIn();

    }

    else if(exp_task==GIMTree_Build||exp_task==GIMTreeU_Build ||exp_task==GIMTreeO_Build){  ////5. -e gimtree: 根据当前双色体数据集user, poi+路网划分情况，构建GIM-Tree

        TIME_TICK_START
        cout<<"o_ratio="<<o_ratio<<endl;
        //getchar();

        if(GIMTree_Build==exp_task)
            processData();
        else if(GIMTreeU_Build==exp_task)
            processData_varyingUP(USER,o_ratio);
        else if(GIMTreeO_Build==exp_task)
            processData_varyingUP(BUSINESS,o_ratio);

        TIME_TICK_END   TIME_TICK_PRINT(" processData(GIM-Tree) runtime: ")
        exit(1);
    }

    else if(exp_task== PerformanceU || exp_task== PerformanceO) {  //pu
        /////主程序运行前的初始程序                        ////(-e batch -s nvd -c ../config_test.txt -p mid -n 100 )
        termWeight_Map.clear();
        for(int i=0;i<vocabularySize;i++){  //term weight 的cache表
            termWeight_Map.push_back(-1);
        }


        clock_t  load_startTime,load_endTime;
        load_startTime = clock();
        loadData();  //jiazai
        load_socialMatrix();
        load_endTime = clock();
        cout<<"loadData total running time:"<<(double)(load_endTime-load_startTime)/CLOCKS_PER_SEC<<"sec!"<<endl;

        int objectType = -1; int need_nvd_cache = 1;
        if(exp_task==PerformanceU) {
            objectType = USER;
            initialIOCache_varyingO(need_nvd_cache, objectType, o_ratio);  //初始化Cache
            Load_Hybrid_NVDG_AddressIdx_fast();
        }
        else if(exp_task==PerformanceO){
            objectType = BUSINESS;
            initialIOCache_varyingO(need_nvd_cache, objectType, o_ratio);  //初始化Cache
            Load_Hybrid_NVDG_AddressIdx_fast_varyingO(o_ratio);
        }


        // 设定随机数的随机种子
        dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));

        int poi;

#ifdef LasVegas
        poi = 265;
#endif
#ifdef LV
        poi = 57;  //57; ////57(14)  //60( );//155( ) ; //79

#endif

#ifdef Brightkite   //for Brightkite map
        poi = 161656; //161662; //161656; //323200 ;//

#endif

#ifdef Gowalla
        poi =  10;
#endif
        clock_t maxinf_startTime, maxinf_endTime;
        double accumulation_time = 0;
        double accumulation_inf = 0; double accumulation_rrset = 0;
        stringstream tt;  stringstream tt2; string outputPathPre;
        string target_info;
        vector<int> candidateKeywords;
        vector<int> checkIn_usr;
        vector<int> chainStores;

        if(true) {  //pichuli   -e batch -s group -c ../config_test.txt -p mid -n 100

            //产生查询点的空间位置

            //chainStores = generate_query_popular(poi, QN);//generate_query();
#ifdef LasVegas
            chainStores = generate_query_PopularDominantNonEmpty(group_id, QN);
                //chainStores = generate_query_PopularDominant(group_id, QN);
#else
            chainStores = generate_query_PopularDominant(group_id, QN);
#endif

            //产生查询点的关键词信息
            candidateKeywords = extractQueryGroupKeywords(
                    chainStores);//extractQueryKeyowrds(poi,DEFAULT_KEYWORDSIZE); //POIs[poi].keywords;//

            if(exp_task== PerformanceU){
                tt2 << "u_ratio =" << o_ratio;
                target_info = tt2.str();
                tt << FILE_LOG << "MaxInf_" << dataset_name<<".u_raio";
            }

            else if(exp_task== PerformanceO){
                tt2 << "o_ratio =" << o_ratio;
                target_info = tt2.str();
                tt << FILE_LOG << "MaxInf_" << dataset_name<<".o_raio";
            }


        }

        outputPathPre = tt.str();
        tt.str("");
        OVERALL_START
        if(alg == BASELINE){ //true){
            tt<<outputPathPre<<"-baseline";
            exp_results.open(tt.str(), ios::app);
            MQueryResults results;
            //加载距离边界
            load_startTime = clock();
            loadDisBound();   //这步可以序列化
            load_endTime = clock();
            cout<<"load disbound running time:"<<(double)(load_endTime-load_startTime)/CLOCKS_PER_SEC<<" sec!"<<endl;
            //用于社交传播影响力分析
            InfluenceModels *  L = new InfluenceModels();
            L->doInitial(argc, argv, accuracy, UserID_MaxKey,online_users, wholeSocialLinkMap, chainStores, OPC);    //初始化社交影响力评估模块

            for(int i =0;i<repeat_baseline;i++){
                //1.baseline方法
                maxinf_startTime = clock();
                results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, chainStores, candidateKeywords, BATCH, GROUP_BATCH, exp_results);

                int rrset_size = MaxInfBRGSTkNNQ_Batch_OPIMC(argc, argv, accuracy, online_users, wholeSocialLinkMap,chainStores, results, b, L, exp_results);
                accumulation_rrset +=rrset_size;

                maxinf_endTime = clock();
#ifdef MillionSecond
                double current_runtime = (double)(maxinf_endTime-maxinf_startTime)/CLOCKS_PER_SEC*1000;
                cout<<"*********************MaxInf by Baseline(Group+OPIMC), runtime="<< current_runtime << "ms!" << endl;
                exp_results <<"--------"<<i<<"th testing: "<<target_info<<"runtime of  baseline solution:"<< current_runtime << " ms!";
#else
                double current_runtime = (double)(maxinf_endTime-maxinf_startTime)/CLOCKS_PER_SEC;
                cout<<"*********************MaxInf by Baseline(Group+OPIMC), runtime="<< current_runtime << "sec!" << endl;
                exp_results <<"--------"<<i<<"th testing: "<<target_info<<", runtime of  baseline solution:"<< current_runtime << " sec!";
#endif
                accumulation_time += current_runtime;



#ifdef MC_ESTIMATION
                TIME_TICK_START
                double current_inf = L -> influence_esitimation_ByMC(10000,OPC);
                accumulation_inf += current_inf;
                exp_results<<"current_inf:"<< current_inf <<endl;
                TIME_TICK_END
                TIME_TICK_PRINT("MC 10000")
#endif
            }
            if(L!=NULL) {
                cout<<"删除L模块..."<<endl;
                delete L; cout<<"删除L成功！"<<endl;
            }

#ifdef MillionSecond
            exp_results <<"***********"<<target_info<<", average runtime baseline solution:"<< accumulation_time/repeat_baseline <<"ms, average rrset size="<<accumulation_rrset/repeat_baseline<<",average inf="<<(accumulation_inf/repeat_baseline) << "*********************"<<endl;
#else
            exp_results <<"***********"<<target_info<<", average runtime baseline solution:"<< accumulation_time/repeat <<"sec, average rrset size="<<accumulation_rrset/repeat<<",average inf="<<(accumulation_inf/repeat) << "*********************"<<endl;
#endif


        }
        else if(alg == APPROXIMATE){ //true){ // -s appro-p pop -n 25 -c ../config_test.txt

            tt<<outputPathPre<<"-appro";
            exp_results.open(tt.str(), ios::app);
            //用于社交传播影响力分析
            InfluenceModels *  L = new InfluenceModels();
            L->doInitial(argc, argv, accuracy, UserID_MaxKey,online_users, wholeSocialLinkMap, chainStores, OPC);    //初始化社交影响力评估模块

            for(int i =0;i<repeat;i++){
                maxinf_startTime = clock();
                BatchResults results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, chainStores, candidateKeywords, BATCH, NVD_BATCH, exp_results);

                int rrset_size = MaxInfBRGSTkNNQ_Batch_HYBRID(argc, argv, accuracy, online_users, wholeSocialLinkMap,chainStores, results, b, L, exp_results);
                maxinf_endTime = clock();
                double current_runtime = (double)(maxinf_endTime-maxinf_startTime)/CLOCKS_PER_SEC*1000;
                accumulation_time += current_runtime;
                accumulation_rrset+= rrset_size;
                cout<<"*********************MaxInf by Appro(NVD+Hybrid), runtime="<<current_runtime << "ms!*********************" << endl;

#ifdef MC_ESTIMATION
                TIME_TICK_START
                //L -> influence_esitimation_ByMC(10000,OPC);
                double current_inf = L -> influence_esitimation_ByMC(10000,OPC);
                accumulation_inf += current_inf; //exp_results<<"current_inf:"<< current_inf <<endl;
                TIME_TICK_END
                TIME_TICK_PRINT("MC 10000")
#endif
                exp_results <<"--------"<<i<<"th testing :"<<target_info<<", runtime of  approximate solution:"<< current_runtime << "ms!"<<"current_inf:"<< current_inf <<endl;

            }

            if(L!=NULL) {
                cout<<"删除L模块..."<<endl;
                delete L; cout<<"删除L成功！"<<endl;
            }

            exp_results <<"***********"<<target_info<<",average runtime of approximate solution:"<< accumulation_time /repeat << "ms, average rrset size="<<accumulation_rrset/repeat<<", average inf="<<(accumulation_inf/repeat) << "*********************"<<endl;



        }

        else if(alg == HEURISTIC){
            tt<<outputPathPre<<"-heuristic";
            string _tmp = tt.str();
            exp_results.open(tt.str(), ios::app);
            InfluenceModels *  L = new InfluenceModels();
            //printBatchRkGSKQResults(chainStores, fastFilter_output.batch_results);
            ////执行heuristic solution，输出poi seeds
            vector<int> poiSeeds;
            for(int i =0;i<repeat;i++){
                float quality_tradoff = 0.5;
                maxinf_startTime = clock();
                BatchFastFilterResults fastFilter_output= Fast_Batch_NVDFilter(chainStores, k,DEFAULT_A,alpha);
                poiSeeds = MaxInfBRGSTkNNQ_POISelect_Heuristic(argc, argv, k,alpha, quality_tradoff, online_users, wholeSocialLinkMap, chainStores, fastFilter_output, b, L,exp_results);
                maxinf_endTime = clock();
                double current_runtime = (double)(maxinf_endTime-maxinf_startTime)/CLOCKS_PER_SEC*1000;
                accumulation_time += current_runtime;
                cout<<"*********************MaxInf by heuristic, runtime="<<(double)(maxinf_endTime-maxinf_startTime)  << "ms!" << endl;
                exp_results <<"-------- "<<i<<"th testing ("<<target_info<<"):"<<"runtime of  heuristic solution: "<< current_runtime << "ms" << endl;
                cout<<"heuristic solution 所获得的poi seeds:"<<endl;
                for(int poi:poiSeeds){
                    cout<<"p"<<poi<<",";
                    exp_results<<"p"<<poi<<",";
                }
                cout<<endl; exp_results<<endl;
            }
            ////验证poiSeeds的结果质量,体现在其影响力的大小上
            set<int> userSeedSet;
            vector<int> poiSeedsTextUnion = extractQueryGroupKeywords(poiSeeds);//extractQueryKeyowrds(poi,DEFAULT_KEYWORDSIZE); //POIs[poi].keywords;//
            BatchResults results = Answering_Multi_BRGSKTkNN(k, DEFAULT_A, alpha, poiSeeds, poiSeedsTextUnion, BATCH, NVD_BATCH, exp_results);
            for(int poi:poiSeeds){
                for(ResultDetail r: results[poi]){
                    int u_id = r.usr_id;
                    userSeedSet.insert(r.usr_id);
                }
            }////获得poi seed中所有兴趣点的潜在用户并集
            ////蒙特卡洛采样以评估给定用户集的social influence
            double poi_inf = L-> opp-> MCSimulation4GivenSeeds(10000, userSeedSet);
            printf("poi seed的 influence set为：\n");
            printSeedSet(userSeedSet);
            exp_results <<target_info<<",heuristic solution 所获得的 poi seed 的 influence:"<< poi_inf << endl;
            exp_results <<"***********average runtime of heuristic solution:"<< accumulation_time/repeat << "ms！***********" << endl;

        }


        CloseDiskComm();
        exp_results.close();

    }

    else if(exp_task==Build_NVD){  ////6. -e nvd: 构建 Network Voronoi Diagram 相关的索引信息文件（输出 NVD 文件 以及 v2P 文件）
        //(-e nvd -c ../config_test.txt -s group -p mid -n 100 reverse)
        int need_nvd_cache = 0;
        loadRoadNetworkIndexData();
        loadOnlyPMap();
        initialIOCache(need_nvd_cache);  //初始化Cache
        TIME_TICK_START


#ifdef HybridHash
        //NVD_generation_AllKeyword_V2PHash_Enhenced_gtree();
        NVD_generation_AllKeyword_V2PHash_Enhenced_road();
        TIME_TICK_END  TIME_TICK_PRINT("building (hybrid) NVD runtime: ")
#else
        NVD_generation_AllKeyword_Hash_serial();
        TIME_TICK_END  TIME_TICK_PRINT("building NVD runtime: ")
#endif

        //test_load_V2P_NVDG_AddressIdx();

        exit(1);

    }

    else if(exp_task==BuildO_NVD){  ////6. -e nvdO: 在不同大小的P集合下构建 Network Voronoi Diagram 相关的索引信息文件（输出 NVD 文件 以及 v2P 文件）
        //(-e nvd -c ../config_test.txt -s group -p mid -n 100 reverse)
        int need_nvd_cache = 0;
        loadRoadNetworkIndexData();
        loadOnlyPMap(o_ratio);
        initialIOCache_varyingO(need_nvd_cache,BUSINESS,o_ratio);  //初始化Cache
        TIME_TICK_START


#ifdef HybridHash
        //NVD_generation_AllKeyword_V2PHash_Enhenced_gtree();
        NVD_generation_AllKeyword_V2PHash_Enhenced_road_varyingO(o_ratio);
        TIME_TICK_END  TIME_TICK_PRINT("building(varying |O|) (hybrid) NVD runtime: ")
#else
        NVD_generation_AllKeyword_Hash_serial();
        TIME_TICK_END  TIME_TICK_PRINT("building NVD runtime: ")
#endif

        //test_load_V2P_NVDG_AddressIdx();

        exit(1);

    }

    else if(exp_task == TestCode){ //ceshi   ////7. 进行 NVD、V2P部分的测试：-e test

        // *********************** below is for test *************************

        /*test_write_HybirdPOINN_IDX();
        test_load_HybirdPOINN_AddressIdx();*/
        //loadFriendShipBinaryData();
        //loadCheckinBinaryData();

        /*test_CheckSocialConnection();
        exit(1);*/



        //loadData();
        int need_nvd_cache = 1;
        loadData_test();
        initialIOCache(need_nvd_cache);  //初始化Cache
        Load_Hybrid_NVDG_AddressIdx_fast();


#ifdef TESTGIMTree
        //test_lclLeaf();
        //testFunctionForEdge();
        //testTermWeight();
        //testTermInvInfo();
        //testSPSP();
        //getTermSize();
        //testUserTermRelatedEntry_all();
        //testObjectTermRelatedEntry_all_self(); ////测试GIM Tree 必须用
        //testAccessPOI();
#endif



      /*  test_u2PDist();*/
        k = 20;
        ///test_getBorder_SB_NVD(k, DEFAULT_A, alpha);
        //test_u2PDist_single(3413,67);   //o_id, u_id
        int u_id = 53946; int p_id = 161656; int term = 40;
        //test_GSKScoreAndTopk (u_id, p_id,DEFAULT_A,alpha,k);
        int vertex_id = 802; int term_id = 80;
        u_id = 45044;
        //test_checkUserWhere(u_id);
        ////测试NVD必须用!
        ///test_poiNN_by_Hash_all_round();
        test_NVDAdj_GivenKeyword(1);
        ////测试NVD必须用!
        exit(1);
        //cout<<SPSP(2510,350)<<endl;
        //test_u2PDist(7813,65);
        //test_u2PDist(64,65);




    }


}
