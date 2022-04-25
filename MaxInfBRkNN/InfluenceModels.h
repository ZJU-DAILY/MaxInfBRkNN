#ifndef InfluenceModels_H
#define InfluenceModels_H

//class core;
#include "./SocialInfluenceModule/CELF_PLUS/common.h"
#include "./SocialInfluenceModule/CELF_PLUS/anyoption.h"
#include "./SocialInfluenceModule/CELF_PLUS/MC.h"
#include "./SocialInfluenceModule/PMC/pmc.hpp"
#include "imm.h"
#include "OPIMC.h"
#include "unordered_map"


//#include "utility.h"


using namespace _MC;



class InfluenceModels {

    int pi, debug;
    unsigned int maxTuples;
public :
    bool isInitial = false;

    MC* mc = NULL;

    Imm* imm= NULL;

    OPIMC* opp = NULL;

    InfluenceMaximizer_PMC* pmc = NULL;

    AnyOption* opt = NULL;

    // Data Structures
//	bool AM_graph[users][users];
//	unsigned int IM[users][actions] ;

    time_t start, wastedTime;
    const char* outdir;
    string training_dir;

    int actionsProcessTillNow;

    int dataset_size;

    // phase
    // 1 : training only
    // 2 : both training and testing
    // 3 : testing only
    int phase;

public:
    InfluenceModels();
    ~InfluenceModels();

    void doAll(int argc, char* argv[]);



    void doInitial(int argc, char* argv[], int online_MaxKey, vector<int>& online_users, unordered_map<int,vector<int>>&  wholeGraph, vector<int>& pois, int alg )
    {
        cout<<"InfluenceModels:: doInitial"<<endl;
        // Lets harcode oneAction option for now
        actionsProcessTillNow = 0;

        time (&start);
        wastedTime = 0;


        //cout<<"enter readOptions"<<endl;
        opt = readOptions(argc, argv);



        //cout << "PID of the program: " << getpid() << endl;

        string command = string("mkdir -p ") + outdir ;
        system(command.c_str());

        command = string("mkdir -p temp") ;
        system(command.c_str());


        //MC* mc = new MC(opt);
        if(alg == CELF){
            mc = new MC(opt);

        }
        else if (alg == PMC){
            pmc = new InfluenceMaximizer_PMC(opt);
        }


        transferWholeSocialNetwork(wholeGraph,alg);
        //mc->doAll();
        this -> isInitial = true;
    }


    void doInitial(int argc, char* argv[], float accuracy, int online_MaxKey, vector<int>& online_users, unordered_map<int,vector<int>>&  wholeGraph, vector<int>& pois, int alg )
    {
        cout<<"InfluenceModels:: doInitial"<<endl;
        // Lets harcode oneAction option for now
        actionsProcessTillNow = 0;

        time (&start);
        wastedTime = 0;


        //cout<<"enter readOptions"<<endl;
        opt = readOptions(argc, argv);



        //cout << "PID of the program: " << getpid() << endl;

        string command = string("mkdir -p ") + outdir ;
        system(command.c_str());

        command = string("mkdir -p temp") ;
        system(command.c_str());


        //MC* mc = new MC(opt);
        if(alg == CELF){
            mc = new MC(opt);

        }
        else if (alg == PMC){
            pmc = new InfluenceMaximizer_PMC(opt);
        }
        else{
            int u_num = -1;
            int _size = online_users.size();
            cout<<_size<<","<<online_MaxKey<<endl;
            u_num = online_MaxKey;

            if(alg == IMM){

                int p_num = pois.size();
                imm = new Imm(wholeGraph, u_num, p_num, MaxInfBRGST); // IM
            }
            else if(alg==OPC){
                int p_num = pois.size();
                opp = new OPIMC(wholeGraph, accuracy, u_num, p_num, MaxInfBRGST);  //IM
            }

        }



        transferWholeSocialNetwork(wholeGraph,alg);

        //mc->doAll();
    }

    void doModuleClear()
    {

        //mc->doAll();
    }


    void doInitial4Heuris(float accuracy, int userID_MaxKey, vector<int>& online_users, unordered_map<int,vector<int>>&  wholeGraph, vector<int>& pois, int alg )
    {

        int p_num = pois.size();
        opp = new OPIMC(wholeGraph, accuracy, userID_MaxKey, p_num, MaxInfBRGST);  //IM


        transferWholeSocialNetwork(wholeGraph,alg);
        //mc->doAll();
    }

    void doInitial_PreSampling_MC(int argc, char* argv[], float ratio, vector<int>& online_users, AdjListSocial&  wholeGraph)
    {
        cout<<"InfluenceModels:: doInitial for MC simulation"<<endl;
        int u_num = online_users.size();
        opp = new OPIMC(wholeGraph, ratio,u_num);  //IM

        transferWholeSocialNetwork(wholeGraph,OPC);
        //mc->doAll();
    }



    void transferWholeSocialNetwork(unordered_map<int,vector<int>>& wholeSocialLinkMap, int type){
        //  对wholeSocialLinkMap进行遍历
        // 采用weight cascade 对边上概率赋值并将边信息重新转入AM
        cout<<"transferWholeSocialNetwork...";
        //unordered_map<int,set<int>>::iterator it;
        auto it = wholeSocialLinkMap.begin();
        int edgeNum =0;
        while(it != wholeSocialLinkMap.end()) {
            int u2 = it->first;

            vector<int> infriendsOfU2 = it->second;
            int size =0; size= infriendsOfU2.size();
            if(size>0){
                double prob = 1.0/ infriendsOfU2.size();
                for(int u1: infriendsOfU2){

                    if(type == CELF){
                        //addSocialLinkMC(u1,u2,prob);
                        mc->addSocialLink(u1,u2,prob);
                    }
                    else if (type == PMC){
                        pmc->addSocialLinkPMC(u1,u2,prob);
                    }
                    else if(type == IMM){
                        imm->addSocialLinkIMM(u1,u2,prob);
                    }
                    else if(type == OPC){
                        opp->addSocialLinkIMM(u1,u2,prob);
                    }

                   /* if (edgeNum%20000==0) {
                        cout<<"插入第"<<edgeNum<<"条 social link:(u"<<u1<<",u"<<u2<<")"<<endl;
                    }*/
                    edgeNum++;
                }
                it ++;
            }else{
                it ++;
                continue;
            }


        }
        cout<<"COMPLETE！"<<endl;
    }



    void seperateThenPMC (vector<int>& stores, BatchResults& results, int b, int R);


    void batchThenCelf (vector<int>& stores, BatchResults& results, int b);


    void batchThenIMM (vector<int>& stores, BatchResults& results, int b);


    int batchThenOPIMC (vector<int>& stores, BatchResults& results, int b,ofstream& exp_resultOutput);

    int batchThenOPIMC_Hybrid (vector<int>& stores, BatchResults& results, int b,ofstream& exp_resultOutput);

    void batchThenHopbased (vector<int>& stores, BatchResults& results, int b,ofstream& exp_resultOutput);

    void test_batchThenOPIMC (vector<int>& stores, BatchResults& results, int b);


    double influence_esitimation_ByMC(int round, int method){
        if(method==OPC){
            return  opp -> MCSimulation4Seeds(round);
        }

        else if(method ==PMC){
            return  pmc -> MCSimulation4Seeds(round);
        }
    }


    void preSampling_MC(int R){
        opp -> Pre_Sampling_MC(R);
    }

    int getCurrentSize (int alg);



    void writeSummary(size_t tcount);

    time_t getStartTime() const {
        return start;
    }

    time_t getWastedTime() const {
        return wastedTime;
    }

    void setWastedTime(time_t t) {
        wastedTime = t;
    }
//	float getCurrentMemoryUsageNew();

private:
    AnyOption* readOptions(int argc, char* argv[]);

};
#endif
