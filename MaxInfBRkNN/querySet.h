//
// Created by jins on 10/31/19.
//
#ifndef MAXINFRKGSKQ_1_1_QUERYSET2_H
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <queue>
#include <sstream>
#include "config.h"

using namespace std;

struct query {
    int Ni, Nj, id;
    vector<int> keywords;
    float dist, dis;
};

// ************************ load user query log file *****************************
int loadQuery(vector<query> &querySet) {
    ifstream fin;
    fin.open(FILE_QUERY);
    string str;
    while (getline(fin, str)) {
        istringstream tt(str);
        query tmp_Q;
        tt >> tmp_Q.id >> tmp_Q.Ni >> tmp_Q.Nj >> tmp_Q.dist >> tmp_Q.dis;
        int key;
        vector<int> keyTmp;
        while (tt >> key) {
            if (-1 == key) {
                break;
            }
            keyTmp.push_back(key);
        }

        tmp_Q.keywords = keyTmp;

        querySet.push_back(tmp_Q);
        //if(querySet.size() == QueryNum)
        //   break;
    }
    fin.close();
    cout<<"LOADING QUERY... DONE! "<<endl;
    return 0;
}

// print query file
int printQuery(query queryItem) {
    cout << queryItem.id << " " << queryItem.Ni << " " << queryItem.Nj << " " << queryItem.dist << " " << queryItem.dis;
    for (size_t j = 0; j < queryItem.keywords.size(); j++) {
        cout << " " << queryItem.keywords[j];
    }
    cout << endl;
    return 0;
}

int testQuery() {
    vector<query> querySet;
    loadQuery(querySet);
    cout<<querySet.size()<<endl;
    for (size_t i = 0; i < 100; i++) {
        printQuery(querySet[i]);
    }
    return 0;
}


#define MAXINFRKGSKQ_1_1_QUERYSET2_H

#endif //MAXINFRKGSKQ_1_1_QUERYSET2_H
