//
// Created by jins on 12/7/19.
//
//
// Created by jins on 12/7/19.
//
#include "dijkstra.h"
#include "gim_tree.h"
#include <math.h>
#include "config.h"
#include  <metis.h>

#ifndef BATCH_RKGSKQ_BUILDIDX_H


// MACRO for timing
struct timeval tv2;
long long ts2, te2;
#define TIME_TICK_START gettimeofday( &tv2, NULL ); ts = tv2.tv_sec * 1000 + tv2.tv_usec / 1000;
#define TIME_TICK_END gettimeofday( &tv2, NULL ); te = tv2.tv_sec * 1000 + tv2.tv_usec / 1000;
#define TIME_TICK_PRINT(T) printf("%s RESULT: %lld (MS)\r\n", (#T), te2 - ts2 );
// ----------



int randomUsrTextSize(){
    return random()%4+1;

}

void EntheticUsr(){
    ifstream finUsr;
    finUsr.open(FILE_USER);
    ofstream foutUsr;
    foutUsr.open(FILE_USER2);
    ifstream finDoc;
    finDoc.open(FILE_DOC2);
    string str;string str2;
    int User_id, User_Ni, User_Nj;
    float User_dist, User_dis;

    while (getline(finUsr, str)) {
        getline(finDoc, str2);
        istringstream tt(str);
        istringstream tt2(str2);
        tt >> User_id >> User_Ni >> User_Nj >> User_dist >> User_dis;
        foutUsr << User_id << " " << User_Ni << " " << User_Nj << " " << User_dist << " " << User_dis << " ";
        //cout << User_id <<" "<< User_Ni <<" "<< User_Nj <<" "<< User_dist <<" "<< User_dis<<"\n";
        vector<int> uKey;
        int keyTmp;
        int popular_poi_count = 0;
        int id;
        tt2 >> id;
        while (tt2 >> keyTmp) {
            if (-1 == keyTmp) {
                break;
            }
            uKey.push_back(keyTmp);

        }
        for (int i = 0; i < randomUsrTextSize(); i++) {
            if (i == uKey.size()) break;
            int temp = uKey[i];
            if (temp > 300)
                temp = temp % 300;
            foutUsr << temp << " ";
        }
        foutUsr << endl;
    }

    finUsr.close();
    finDoc.close();
    foutUsr.close();

}

void EntheticPoi(){
    ifstream finPoi;
    finPoi.open(FILE_POI);
    ifstream finDoc;
    finDoc.open(FILE_DOC);
    ofstream foutPoi;
    foutPoi.open(FILE_POI2);

    string str; string str2;
    int POI_id, POI_Ni, POI_Nj;
    float POI_dist, POI_dis;

    while (getline(finPoi, str)) {
        getline(finDoc, str2);
        istringstream tt(str);
        istringstream tt2(str2);
        tt >> POI_id >> POI_Ni >> POI_Nj >> POI_dist >> POI_dis;
        //foutPoi << POI_id <<" "<< POI_Ni <<" "<< POI_Nj <<" "<< POI_dist <<" "<< POI_dis<<"\n";
        foutPoi << POI_id <<" "<< POI_Ni <<" "<< POI_Nj <<" "<< POI_dist <<" "<< POI_dis<<" ";
        vector<int> uKey;
        int keyTmp;
        int id; tt2 >> id;
        while (tt2 >> keyTmp) {
            if (-1 == keyTmp) {
                break;
            }
            if(keyTmp>500)
                keyTmp = keyTmp % 500;
            foutPoi<< keyTmp<<" ";
        }
        //cout<<endl;
        foutPoi<<endl;
        //break;
        //getchar();
    }
    finPoi.close();
    finDoc.close();
    foutPoi.close();

}


void buildIndex(){
    TIME_TICK_START
    init();
    TIME_TICK_END
    TIME_TICK_PRINT("INIT")

    // gtree_build
    TIME_TICK_START
    build();
    TIME_TICK_END
    TIME_TICK_PRINT("BUILD")

    // dump gtree
    gtree_save();

    // calculate distance matrix
    TIME_TICK_START
    hierarchy_shortest_path_calculation();
    TIME_TICK_END
    TIME_TICK_PRINT("MIND")

    // dump distance matrix
    hierarchy_shortest_path_save();

    TIME_TICK_START
    generateObject();
    TIME_TICK_END
    TIME_TICK_PRINT("OBJECT")

    int leafCnt = 0;

    for(size_t i = 0;i<GTree.size();i++)
    {
        if(GTree[i].isleaf == true)
        {
            leafCnt++;
        }
    }

    cout<<" leaf cnt is "<<leafCnt<<endl;

}

void loadIndex(){
    init();
    // load gtree index
    gtree_load();

    // load distance matrix
    hierarchy_shortest_path_load();

    //test SPSP
    //cout<<SPSP(10,12)<<endl;
}


int randomTextSize(){
    return random()%12+0;

}

double ZipfMaxVal;

void InitZipfMaxVal(int maxnum,double theta) {
    ZipfMaxVal = 0.0;
    for(int i=1;i<=maxnum;i++)
        ZipfMaxVal += 1.0/pow( (double)i, theta);
}

int zipf(double theta) {
    double r= drand48()*ZipfMaxVal, sum= 1.0;
    int i=1;
    while( sum<r){
        i++;
        sum += 1.0/pow( (double)i, theta);
    }
    return (i-1);
}




#define BATCH_RKGSKQ_BUILDIDX_H

#endif //BATCH_RKGSKQ_BUILDIDX_H
