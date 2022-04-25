//
// Created by jins on 4/16/20.
//

#ifndef MAXINFBRGSTKNN_RKGSKQ_H

#include "process.h"

//多查询处理
void RkGSKQ_Multi(int Qk, float a, float alpha, vector<int> query_objects, vector<int> candidateKeywords, int processing_pattern, int algorithm){



    if(processing_pattern == BATCH){
        //for LV map
        //cout<<"待处理"<<POIs.size()<<"个poi的RkGSKQ"<<endl;getchar();
        Batch_START

        //loc = POIs[poi].Nj;
        cout<<"-----------------------------Batch RkGSKQ Retrieving!---------------------------------"<<endl;
        //cout<<"POI:"<<poi<<endl;
        cout<<"keywords, size ="<<candidateKeywords.size()<<"  elements:"; printQueryKeywords(candidateKeywords);


        // xuanze
        if(algorithm == CLUSTER){
#ifdef NONGIVENKEY
            RkGSKQ_Batch_Baseline_formal(query_objects, Qk,a,alpha);
#else
            RkGSKQ_Batch_Baseline(query_objects, candidateKeywords,Qk,a,alpha);
#endif


        }
        else if(algorithm == GROUP_BATCH){
            BatchResults  batch_results = RkGSKQ_Batch_Group(query_objects, candidateKeywords,Qk,a,alpha);
            printBatchRkGSKQResults(query_objects, batch_results);

        }
        else {


        }

        Batch_END;
        cout<<"total process:"<<query_objects.size()<<",";
        Batch_PRINT("Batch RkGSKQ running time:");
        cout<<"磁盘IO: number of block accessed="<<block_num<<endl;
        cout<<"by"<<algorithm<<endl;
        cout<<"--------------------------------------------------------------------------------------"<<endl;
    }

    else if(processing_pattern ==Naive){

        //Naive_MRkGSKQ_givenKeys(query_objects, candidateKeywords,Qk, DEFAULT_A, alpha);
        Naive_MRkGSKQ(query_objects, Qk, DEFAULT_A, alpha);


    }

    else if(processing_pattern == ONE_BY_ONE){  // process RkGSKQs one by one
        Group_START
        cout<<"-----------------------------Seprate_RkGSKQ Processing!---------------------------------"<<endl;
        for(int store :query_objects){
            //store = 85;
            cout<<"******************开始 POI:"<<store<<"的处理******************"<<endl;
            cout<<"keywords, size ="<<candidateKeywords.size()<<"  elements:";
            POI p;
#ifndef DiskAccess
            p = POIs[store];
#else
            p = getPOIFromO2UOrgLeafData(store);
#endif

#ifdef NONGIVENKEY
            candidateKeywords = p.keywords;
#endif
            printQueryKeywords(candidateKeywords);
            vector<int> checkIn_usr = poiCheckInIDList[store];
#ifndef DiskAccess
            RkGSKQ_Single_poi_Plus(p, candidateKeywords, checkIn_usr, Qk,a,alpha);
#else
            RkGSKQ_Single_poi_Disk(p, candidateKeywords, checkIn_usr, Qk,a,alpha);
#endif
            //cout<<"**************完成 POI:"<<store<<"的处理**************"<<endl;

        }
        Group_END;
        cout<<"total process:"<< query_objects.size()<<",";
        Group_PRINT("SepRkGSK running time:");
        cout<<"by single RkGSKQ"<<endl;
        cout<<"--------------------------------------------------------------------------------------"<<endl;
    }

}

MQueryResults sythetic_Multi_BRGSKTkNNResults(vector<int> query_objects){
    int Candidate_num = query_objects.size();
    map<int, vector<ResultDetail>> batch_Results;
#ifdef LasVegas
    loadLasVegasUserInfo();
#endif
    for(int poi_id: query_objects){
        int potential_userNum = random() % 20; //随机设定对poi感兴趣的潜在用户个数

        for(int i=0;i<potential_userNum;i++){
#ifdef LasVegas
            int u_th = random() % Users.size(); //随机选择感兴趣用户
            int u_id = Users[u_th].id;
#else
            int u_id = random() % UserID_MaxKey;
#endif
            ResultDetail rd(u_id, 0, 0, 0, -1, -1);
            batch_Results[poi_id].push_back(rd);
        }
    }
    printBatchRkGSKQResults(query_objects, batch_Results);
    return batch_Results;
}



MQueryResults Answering_Multi_BRGSKTkNN(int Qk, float a, float alpha, vector<int> query_objects, vector<int> candidateKeywords, int processing_pattern, int algorithm, ofstream& exp){

    MQueryResults mr;

    if(processing_pattern == BATCH){
        //for LV map
        //cout<<"待处理"<<POIs.size()<<"个poi的RkGSKQ"<<endl;getchar();


        BatchResults br;
        string info_pre;


        //cout<<"POI:"<<poi<<endl;
        cout<<"keywords, size ="<<candidateKeywords.size()<<"  elements:"; printQueryKeywords(candidateKeywords);

        Batch_START
        // xuanze
        if(algorithm == CLUSTER){
            cout<<"-----------------------------Batch RkGSKQ Retrieving (by cluster scheme)!---------------------------------"<<endl;
#ifdef  NONGIVENKEY
            br = RkGSKQ_Batch_Baseline_formal(query_objects,Qk,a,alpha);
#else
            br = RkGSKQ_Batch_Baseline(query_objects, candidateKeywords,Qk,a,alpha);
            info_pre ="cluster";
#endif


        }
        else if(algorithm == GROUP_BATCH){
            cout<<"-----------------------------Batch RkGSKQ Retrieving (by group scheme)!---------------------------------"<<endl;
            br = RkGSKQ_Batch_Group(query_objects, candidateKeywords,Qk,a,alpha);
            info_pre ="group";
        }
        else if(algorithm == NVD_BATCH){
            cout<<"-----------------------------Batch RkGSKQ Retrieving (by nvd scheme)!---------------------------------"<<endl;
#ifdef LasVegas
            br = RkGSKQ_Batch_NVD_RealData(query_objects, Qk,a,alpha);
#else
            br = RkGSKQ_Batch_NVD_SyntheticData(query_objects, Qk,a,alpha);
            info_pre ="nvd";
#endif

        }

        mr = br;

        Batch_END;
        //cout<<"total process:"<<query_objects.size()<<",";
        Batch_PRINT("Batch RkGSKQ running time:");
        //cout<<"磁盘IO: number of block accessed="<<block_num<<endl;
        //exp<<"*"<<info_pre<<" batch runtime="<<(be-bs)/1000000.0<<"sec, io cost ="<<block_num<<endl;
        exp<<"**"<<info_pre<<" batch runtime="<<(be-bs)/1000000.0<<"sec, |Uc| ="<<Uc_size<<endl; Uc_size = 0;

        cout<<"--------------------------------------------------------------------------------------"<<endl;
    }



    else if(processing_pattern == ONE_BY_ONE){  // process RkGSKQs one by one
        Group_START
        cout<<"-----------------------------Seprate_RkGSKQ Processing!---------------------------------"<<endl;
        for(int store :query_objects){
            //store = 85;
            //SingleResults* sr;
            cout<<"******************开始 POI:"<<store<<"的处理******************"<<endl;
            cout<<"keywords, size ="<<candidateKeywords.size()<<"  elements:";
            POI p;
#ifndef DiskAccess
            p = POIs[store];
#else
            p = getPOIFromO2UOrgLeafData(store);
#endif

#ifdef NONGIVENKEY
            candidateKeywords = p.keywords;
#endif
            printQueryKeywords(candidateKeywords);
            vector<int> checkIn_usr = poiCheckInIDList[store];

            if(algorithm == GROUP_SINGLE){

                SingleResults sr = RkGSKQ_Single_poi_Disk(p, candidateKeywords, checkIn_usr, Qk,a,alpha);
                mr[store] = sr.candidateUsr;
            }

            else if(algorithm == NVD_SINGLE){

                //SingleResults sr = RkGSKQ_bottom2up_by_NVD(p, Qk,a,alpha);
                //mr[store] = sr.candidateUsr;
                vector<ResultDetail> candidate_usr = RkGSKQ_NVD_Naive_optimized(p.id, Qk,a,alpha);
                mr[store] = candidate_usr;
            }


            //cout<<"**************完成 POI:"<<store<<"的处理**************"<<endl;

        }
        Group_END;
        printBatchRkGSKQResults(query_objects, mr);
        cout<<"total process:"<< query_objects.size()<<",";
        Group_PRINT("SepRkGSK running time:");
        cout<<"by single RkGSKQ"<<endl;
        cout<<"--------------------------------------------------------------------------------------"<<endl;
    }

    //printBatchRkGSKQResults(query_objects, mr);

    return mr;

}



/*
MQueryResults Faster_SepFilter_BRGSKTkNN(int Qk, float a, float alpha, vector<int> query_objects,ofstream& exp){



    //Group_START
    //cout<<"-----------------------------Seprate_RkGSKQ Processing!---------------------------------"<<endl;
    for(int store :query_objects){
        //store = 85;
        //SingleResults* sr;
        cout<<"******************开始 POI:"<<store<<"的处理******************"<<endl;
        //cout<<"keywords, size ="<<candidateKeywords.size()<<"  elements:";
        POI p;

        p = getPOIFromO2UOrgLeafData(store);

        vector<int> checkIn_usr = poiCheckInIDList[store];

        //SingleResults sr = RkGSKQ_bottom2up_by_NVD(p, Qk,a,alpha);
        vector<int> fasterFilter_candidate = Faster_SingleFilter_by_NVD(p, Qk,a,alpha);
        cout<<"fasterFilter_candidate size="<<fasterFilter_candidate.size()<<endl;

        //mr[store] = sr.candidateUsr;


        //cout<<"**************完成 POI:"<<store<<"的处理**************"<<endl;

    }

    //Group_END;
    cout<<"total process:"<< query_objects.size()<<",";
    //Group_PRINT("SepRkGSK running time:");
    cout<<"by single RkGSKQ"<<endl;
    cout<<"--------------------------------------------------------------------------------------"<<endl;


    //printBatchRkGSKQResults(query_objects, mr);

    return mr;

}
*/


MQueryResults Faster_BatchFilter_BRGSKTkNN(int Qk, float a, float alpha, vector<int> query_objects, ofstream& exp){

    MQueryResults mr; //Faster_BatchFilter_NVD(query_objects, Qk,a,alpha);;

    return mr;

}



//RkGSKQ查询处理模块，支持单查询与多查询
void RkGSKQ_Processing(int Qk, float a, float alpha, int task, int mode, int alg){

    int poi;
    int loc;
    vector<int> candidateKeywords ;
    vector<int> checkIn_usr;
    //程序选择分支
    bool sing_query = 0; bool multiple_queries = 1;
    bool SepRkGSK = 1; bool GrpRkGSK = 0;
    bool test = 0;

#ifdef LV
    poi = 57; ////57(14)  //60( );//155( ) ; //60
    string str = "57 60 155 900 80 39";  // user28 info
    //str = "57";
    stringstream ss = stringstream(str);
#endif

#ifdef Brightkite   //for Brightkite map
    poi = 83; //212897 ;//  //83(18)  60(15) 1277(23) 212897(22)
    string str = "83 60 1277 212897 770895";  // user28 info
    stringstream ss = stringstream(str);
#endif

#ifdef Gowalla   //for Brightkite map
    poi = 10; //
    //string str = "83 60 1277 212897 770895";  // user28 info
    //stringstream ss = stringstream(str);
#endif

    if(task == Multi){  // 多查询处理
        //int KQ_size,vector<int>& KQ, vector<int>& checkIn_usr, int frequent_thethod, int Q_size
        vector<int> candidateKeywords, checkIn_usr;
        vector<int> chainStores;
        if(true){
            //产生查询点的空间位置
            chainStores = generate_query_popular(poi);//generate_query();
            //产生查询点的关键词信息

            candidateKeywords = extractQueryKeyowrds(poi,DEFAULT_KEYWORDSIZE); //POIs[poi].keywords;//
            //checkIn_usr = POIs[poi].check_ins;
            if(poiCheckInIDList.count(poi)==1)
                checkIn_usr = poiCheckInIDList[poi];

        }

        RkGSKQ_Multi(Qk, a, alpha, chainStores, candidateKeywords, mode, alg);



    }
    else if(task == Single){   // 执行单个RkGSKQ

        ofstream log;
        stringstream tt;
        tt<<FILE_LOG<<"RkGSKQ_result_"<<dataset_name;
        log.open(tt.str(), ios::app);
        POI p;
#ifdef DiskAccess
        p = getPOIFromO2UOrgLeafData(poi);
#else
        p = POIs[poi]; //POIs[poi];
#endif

        candidateKeywords = p.keywords;
        //printElements(candidateKeywords);
#ifdef GIVENKEY_
        candidateKeywords = extractQueryKeyowrds(poi,DEFAULT_KEYWORDSIZE);
#endif
        checkIn_usr = extractCheckins(poi);

#ifdef DiskAccess
        SingleResults results = RkGSKQ_Single_poi_Disk(p, candidateKeywords, checkIn_usr, Qk,a,alpha);
        cout<<"RkGSKQ_Single_poi_Disk 获得POI:"<<poi<<"的RkGSKQ结果，size="<<results.candidateUsr.size()<<endl;
        log<<"RkGSKQ_Single_poi_Disk 获得POI:"<<poi<<"的RkGSKQ结果，size="<<results.candidateUsr.size()<<endl;
#else
        SingleResults results = RkGSKQ_Single_poi_Plus(p, candidateKeywords, checkIn_usr, Qk,a,alpha);
        cout<<"RkGSKQ_Single_poi_Plus 获得POI:"<<poi<<"的RkGSKQ结果，size="<<results.candidateUsr.size()<<endl;
        log<<"RkGSKQ_Single_poi_Plus 获得POI:"<<poi<<"的RkGSKQ结果，size="<<results.candidateUsr.size()<<endl;
#endif

        log.close();

    }

    exit(0);
}


void RkGSKQ_Processing_single_test(int poi, int Qk, float a, float alpha, int task, int mode, int alg){


    clock_t startTime, endTime;

    startTime = clock();
    int loc;
    vector<int> candidateKeywords ;
    vector<int> checkIn_usr;
    //程序选择分支
    bool sing_query = 0; bool multiple_queries = 1;
    bool SepRkGSK = 1; bool GrpRkGSK = 0;
    bool test = 0;


    if(true){   // 执行单个RkGSKQ

        ofstream log;
        stringstream tt;
        tt<<FILE_LOG<<"SingleRkGSKQ_result_"<<dataset_name;
        log.open(tt.str(), ios::app);
        POI p;
#ifdef DiskAccess
        p = getPOIFromO2UOrgLeafData(poi);
#else
        p = POIs[poi]; //POIs[poi];
#endif

        candidateKeywords = p.keywords;
        //printElements(candidateKeywords);
#ifdef GIVENKEY_
        candidateKeywords = extractQueryKeyowrds(poi,DEFAULT_KEYWORDSIZE);
#endif
        checkIn_usr = extractCheckins(poi);

#ifdef DiskAccess
        SingleResults results = RkGSKQ_Single_poi_Disk(p, candidateKeywords, checkIn_usr, Qk,a,alpha);
        //cout<<"RkGSKQ_Single_poi_Disk 获得POI:"<<poi<<" RkGSKQ结果，size="<<results.candidateUsr.size()<<"， 具体为："<<endl;
        log<<"RkGSKQ_Single_poi_Disk 获得POI:"<<poi<<"的RkGSKQ结果，size="<<results.candidateUsr.size()<<endl;
        //printResults(results.candidateUsr);

#else
        SingleResults results = RkGSKQ_Single_poi_Plus(p, candidateKeywords, checkIn_usr, Qk,a,alpha);
        cout<<"RkGSKQ_Single_poi_Plus 获得POI:"<<poi<<"的RkGSKQ结果，size="<<results.candidateUsr.size()<<"， 具体为："<<endl;endl;
        log<<"RkGSKQ_Single_poi_Plus 获得POI:"<<poi<<"的RkGSKQ结果，size="<<results.candidateUsr.size()<<endl;
        printResults(results.candidateUsr);

#endif

        log.close();

    }
    endTime = clock();


    cout<<"------RkGSKQ_Processing_single for *p"<<poi<<"* 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC*1000  << "ms------" << endl;

}


void test_query_generation(int Qn, int Popularity){


    vector<int> candidateKeywords, checkIn_usr;
    vector<int> chainStores;
        //产生查询点的空间位置
    int KQ = 60;
    chainStores = generate_query_According2Popular_lv(Popularity,Qn,60);//generate_query();
        //提取查询点的关键词信息
    candidateKeywords = extractQueryGroupKeywords(chainStores);





}





#define MAXINFBRGSTKNN_RKGSKQ_H

#endif //MAXINFBRGSTKNN_RKGSKQ_H
