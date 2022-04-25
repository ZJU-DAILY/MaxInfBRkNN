//
// Created by jin on 20-3-2.
//

#ifndef MAXINFBRGSTQ_OPIMC_H
#define MAXINFBRGSTQ_OPIMC_H

#endif //MAXINFBRGSTQ_OPIMC_H

#include "imm.h"
#include "infgraph.h"
#include <iostream>
//#include "config.h"
using namespace std;

#define  POI_SELECTION_TRACK

class OPIMC {   //: public Imm{
public:
    InfGraph* infGraph  = NULL;
    SocialGraph* social_graph = NULL;
    Argument arg;
    double sqr(double v){
        return v*v;
    }

    OPIMC(AdjListSocial& socialGraph, float accuracy, int u_num, int p_num, Problem p) {  //jordan
        social_graph = new  SocialGraph(u_num);
        if(p== MaxInfBRGST){
            infGraph = new InfGraph(u_num, p_num,social_graph);
        }
        else{
            infGraph = new InfGraph(u_num,social_graph);
        }

        arg.epsilon = accuracy;
        //arg.k = 5;
        arg.model = string("IC");
    }

    OPIMC(AdjListSocial& socialGraph, float accuracy, int u_num) {  //jordan
        social_graph = new  SocialGraph(u_num);
        infGraph = new InfGraph(u_num,social_graph);
        arg.epsilon = accuracy;
        //arg.k = 5;
        arg.model = string("IC");
    }

    ~OPIMC(){
        printf("begin ~OPIMC()...");
        if(social_graph!=NULL)
            delete social_graph;
        if(infGraph!=NULL)
            delete infGraph;
        printf("complete!\n");
    }

    void addSocialLinkIMM(int u, int v, double prob){
        infGraph->add_edge(u,v,prob);
        //infGraphVldt->add_edge(u,v,prob);
    }

    double InfluenceMaximize_by_opimc(int b)
    {
        cout<<"InfluenceMaximize_by_opimc..."<<endl;
        //arg.epsilon = 0.05;
        //Timer t(2, "step2");
        //ASSERT(OPT_prime > 0);
        int targetSize = b; arg.k = b;
        InfGraph &g = *infGraph;
        int n = g.social_graph->n;
        //InfGraph &g2 = *infGraphVldt;

        double e = exp(1);
        double approx = 1 - 1.0/e;

        double delta = 1.0 / (2*n);
        double alpha = sqrt(log(6.0/delta));
        double beta = sqrt((1-1/e) * (Math::logcnk(n, arg.k) + log(6.0/delta)));


        const auto numRbase = size_t(2.0 * Math::pow2((1 - 1 / e) * alpha + beta));
        const auto maxNumR = size_t(2.0 * n * Math::pow2((1 - 1 / e) * alpha + beta) / targetSize / Math::pow2(arg.epsilon)) + 1;
        auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
        const double a1 = log(numIter * 3.0 / delta);
        const double a2 = log(numIter * 3.0 / delta);

        g.init_hyper_graph();
        g.init_hyper_graph_vld();
        //numIter = 15;
        for(int round = 0; round<numIter;round++){
            cout<<"--------------------------round"<<round<<"-------------------------"<<endl;
            const auto numR = numRbase << round;
            cout<<"need "<<numR<<" RR sets"<<endl;

            g.build_hyper_graph_r(numR, arg); // R1
#ifdef Selection_TRACK
            cout<<"R1 build ("<<numR<<"RR set)"<<endl;
#endif


            //dsfmt_gv_init_gen_rand(idx + 1000000);
            g.build_hyper_graph_vld_r(numR, arg); // R2
#ifdef Selection_TRACK
            cout<<"R2 build ("<<numR<<"RR set)"<<endl;
#endif
            if(true){
                int numRRsets = g.get_RR_sets_size();

                double infSelf_ratio = g.build_seedset(arg.k);
#ifdef Selection_TRACK
                cout<<"infSelf_ratio="<<infSelf_ratio<<endl;
#endif
                vector<int> cur_seeds = g.seedSet;

                //const auto infVldt = 0 ;//__hyperGVldt.self_inf_cal(__vecSeed); //2. obtain verified Inf(R2: hyperGVldt)
                const auto coverVldt = g.check_RR_Cover(cur_seeds); //obtain Cove(R2)
                //cout<<"coverVldt="<<coverVldt<<endl;


                auto upperBound = infSelf_ratio / approx;  // use 1;
                //if (true) upperBound = upperBound;  // how(default): jins
                const auto upperCover = upperBound * numRRsets;
                const auto lowerbound_S = Math::pow2(sqrt(coverVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
                const auto upperbound_OPT = Math:: pow2(sqrt(upperCover + a2 / 2.0) + sqrt(a2 / 2.0));
                const auto approxOPIMC = lowerbound_S / upperbound_OPT;

                cout<<"lowerbound_S = "<<lowerbound_S<<", upperbound_OPT ="<<upperbound_OPT<<",rate="<<approxOPIMC<<" rr set:"<<numRRsets<<endl;

                if(approxOPIMC>approx-arg.epsilon){
                    cout<<"find satisable seeds selection:"<<endl;
                    for(int u:cur_seeds){
                        cout<<"u"<<u<<endl;
                    }
                    cout<<"number of round:"<<round<<endl;
                    cout<<"number of rr set:"<<numR<<endl;
                    cout<<"expected inf= "<<infSelf_ratio*n<<endl;
                    return 1.0;
                }
            }


        }

        return 0.0;
    }



    void checkReachable(vector<int> pois, BatchResults& results){
        //将所有兴趣的点，与之关联的反近邻用户  关联
        social_graph->rearchInfo_clear();
        for(int p: pois){
            vector<ResultDetail> rrs = results[p];
            if(rrs.size()>0){
                for(ResultDetail r: rrs){
                    int u = r.usr_id;
                    infGraph->social_graph->setRearchable(u,p);
                }
            }
        }
        infGraph->social_graph->traFromSource();
        //查看 source user 的个数
        int rearch_size = infGraph->social_graph->getNs();
    }

    void checkReachable_Upper(vector<int> pois, map<int,map<int,double>>& Lo_Upper_Map){
        //将所有兴趣的点，与之关联的反近邻用户  关联
        social_graph->rearchInfoUpper_clear();

        //set<int> aa;
        map<int,map<int,double>>:: iterator iter;
        for(iter=Lo_Upper_Map.begin();iter != Lo_Upper_Map.end();iter++){
            int p = iter -> first;
            map<int,double> Lo_Upper = iter -> second;
            map<int,double>:: iterator _iter;
            for(_iter= Lo_Upper.begin(); _iter != Lo_Upper.end();_iter++){
                int u = _iter -> first;
                social_graph->setRearchable_Upper(u,p);
                //aa.insert(u);
            }
        }
        //cout<<"aa-size="<<aa.size()<<endl;
        //getchar();

        social_graph->traFromSource_Upper();
        //查看 source user 的个数
        social_graph->ns = social_graph->getNs_upper();
    }


    void setPUConnection(vector<int> pois, map<int,vector<int>>& results){
        //将所有兴趣的点，与之关联的反近邻用户  关联
        social_graph->rearchInfo_clear();
        for(int p: pois){
            vector<int> rrs = results[p];
            if(rrs.size()>0){
                for(int u: rrs){
                    infGraph->social_graph->setRearchable(u,p);
                }
            }
        }
        social_graph->initialPOISelectionInfo();


    }

    void setPUConnection_Upper(vector<int> pois, map<int,vector<int>>& results, int& remained_th,
                               vector<int>& remained_candidate_usersID, unordered_map<int,set<int>>& associated_poiIDMap){
        //将当前已验证的所有兴趣的点与潜在用户间建立 链接关系
        social_graph->rearchInfoUpper_clear();
        for(int p: pois){
            vector<int> rrs = results[p];
            if(rrs.size()>0){
                for(int u: rrs){
                    //int u = r.usr_id;
                    social_graph->setRearchable_Upper(u,p);
                }
            }
        }
        //将"剩余" 所有候选用户与其可能相关的兴趣点间建立 （上界下的）链接关系
        for(int j =remained_th;j<remained_candidate_usersID.size();j++){
            int uj = remained_candidate_usersID[j];
            set<int> _usr_associated_poiIDSet = associated_poiIDMap[uj];
            for(int o: _usr_associated_poiIDSet){
                social_graph->setRearchable_Upper(uj,o);
            }
        }
        social_graph->initialPOISelectionInfo_Upper();

    }


   // void setPUSubAndSuperConnection(vector<int> pois, map<int, vector<int>> results_map, map<int, vector<int>> results_map)

    typedef struct {
        int iteration;
        int rrset_num;
        vector<int> poiSeeds;
    }RISInfo;

    RISInfo minePOIsOPIMC(int b, int& iteration, int& rr_num, vector<int>& pois)
    {
        //arg.epsilon = 0.05;
        //Timer t(2, "step2");
        //ASSERT(OPT_prime > 0);
        int targetSize = b; //arg.k;
        arg.k =b;
        InfGraph &g = *infGraph;
        g.poiSet.clear();
        g.seedSet.clear();

        double e = exp(1);
        double approx = 1 - 1.0/e;
        int nc = social_graph->nc;
        int np = social_graph->np;
        double delta = 1.0 / (2*nc);
        double alpha = sqrt(log(6.0/delta));
        double beta = sqrt((1-1/e) * (Math::logcnk(np, arg.k) + log(6.0/delta)));


        const auto numRbase = size_t(2.0 * Math::pow2((1 - 1 / e) * alpha + beta));
        const auto maxNumR = size_t(2.0 * nc * Math::pow2((1 - 1 / e) * alpha + beta) / targetSize / Math::pow2(arg.epsilon)) + 1;
        auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
        const double a1 = log(numIter * 3.0 / delta);
        const double a2 = log(numIter * 3.0 / delta);
        cout<<"最大需要"<<maxNumR<<"个RR set, 需要"<<numIter<<"轮迭代验证"<<endl;

        g.init_Two_hyper_graph_poi(np); // intial R1, R2;

        RISInfo info;
        for(int round = 0; round<numIter;round++){
#ifdef POI_SELECTION_TRACK
            cout<<"-------------------------round :"<<round<<" ------------------------"<<endl;
#endif
            const auto numR = numRbase << round;

            g.build_hyper_graph_r_poi(numR, arg); // R1
            //dsfmt_gv_init_gen_rand(idx + 1000000);
            g.build_hyper_graph_vld_r_poi(numR, arg); // R2


            int numRRsets = g.get_RR_POI_sets_size();

            //double infSelf_ratio = g.do_poiSelection(arg.k);
            double infSelf_ratio = g.do_poiSelection_UB(arg.k); // using R1 to obtain S and upper bound of I(S)
#ifdef POI_SELECTION_TRACK
            cout<<"infSelf_ratio="<<infSelf_ratio<<endl;
#endif
            vector<int> cur_pois = g.poiSet;

            const auto coverVldt = g.check_RR_POI_Cover(cur_pois); //obtain Cove(R2)

            auto upperBound = infSelf_ratio / approx;  // use 1;
            //if (true) upperBound = upperBound;  // how(default): jins
            const auto upperCover_loose = upperBound * numRRsets;
            //cout<<"jins:upperCover="<<upperCover<<endl;
            const auto upperCover_tight = g.influence_upper;
            //cout<<"upperCover_tight="<<upperCover_tight<<endl;


            const auto lowerbound_S = Math::pow2(sqrt(coverVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
            const auto upperbound_S_loose = Math:: pow2(sqrt(upperCover_loose + a2 / 2.0) + sqrt(a2 / 2.0));
            const auto upperbound_S_tight = Math:: pow2(sqrt(upperCover_tight + a2 / 2.0) + sqrt(a2 / 2.0));

            const auto approxOPIMC_loose = lowerbound_S / upperbound_S_loose;   //upperbound_OPT
            const auto approxOPIMC_tight = lowerbound_S / upperbound_S_tight;
#ifdef POI_SELECTION_TRACK
            cout<<"ratio_tight="<<approxOPIMC_tight<<endl;
            cout<<"coverVldt= "<<coverVldt<<", upperCover_loose ="<<upperCover_loose<<", upperCover_tight="<<upperCover_tight<<endl;
            cout<<"lowerbound_S="<<lowerbound_S<<", upperbound_S_loose="<<upperbound_S_loose<<", upperbound_S_tight ="<<upperbound_S_tight<<endl;
            cout<<"needs rr set:"<<numRRsets<<endl;
#endif

            if(approxOPIMC_tight > approx-arg.epsilon){
                cout<<"accuracy guarantees ="<<(approx-arg.epsilon)<<endl;
                cout<<"find satisable seeds selection:"<<endl;

                for(int p:cur_pois){
                    cout<<"poi"<<p<<", ";
                }
                cout<<"number of round:"<<round<<endl;
                cout<<"number of rr set:"<<numR<<endl;
                cout<<"expected inf= "<<infSelf_ratio*nc<<endl;
                info.iteration =round; iteration = round;
                info.rrset_num = numR; rr_num = numR;
                info.poiSeeds = cur_pois; pois = cur_pois;
                return info;
            }

        }

        info.iteration = numIter; iteration = numIter;
        info.rrset_num = maxNumR; rr_num = maxNumR;
        info.poiSeeds = g.poiSet; pois = g.poiSet;

        return info;
    }



    double minePOIsHopBased(int b)
    {

        //arg.epsilon = 0.05;
        //Timer t(2, "step2");
        //ASSERT(OPT_prime > 0);
        int targetSize = b; //arg.k;
        arg.k =b;
        InfGraph &g = *infGraph;
        g.poiSet.clear();
        g.seedSet.clear();

        double e = exp(1);
        double approx = 1 - 1.0/e;
        int nc = social_graph->nc;
        int np = social_graph->np;
        double delta = 1.0 / (2*nc);
        double alpha = sqrt(log(6.0/delta));
        double beta = sqrt((1-1/e) * (Math::logcnk(np, arg.k) + log(6.0/delta)));


        const auto numRbase = size_t(2.0 * Math::pow2((1 - 1 / e) * alpha + beta));
        const auto maxNumR = size_t(2.0 * nc * Math::pow2((1 - 1 / e) * alpha + beta) / targetSize / Math::pow2(arg.epsilon)) + 1;
        auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
        const double a1 = log(numIter * 3.0 / delta);
        const double a2 = log(numIter * 3.0 / delta);

        double inf = g.do_poiSelection_UB_ByOneHop(arg.k);
        cout<<"POI seeds(by HopBased):"<<endl;
        g.printPOISelectionInfo();
        cout<<"inf(by HopBased)="<<inf<<endl;



        return inf;
    }

    double minePOIsHopBased4Heuristic(int b)
    {


        int targetSize = b; //arg.k;
        arg.k =b;
        InfGraph &g = *infGraph;
        g.poiSet.clear();
        g.seedSet.clear();

        double inf = g.do_poiSelection_Heuristic(arg.k);
        cout<<"POI seeds(by HopBased Heuristic):"<<endl;
        g.printPOISelectionInfo();
        cout<<"inf(by HopBased)="<<inf<<endl;


        return inf;
    }

    double minePOIsHopBased4Heuristic_upper(int b)
    {

        //arg.epsilon = 0.05;
        //Timer t(2, "step2");
        //ASSERT(OPT_prime > 0);
        int targetSize = b; //arg.k;
        arg.k =b;
        InfGraph &g = *infGraph;
        g.poiSet2.clear();
        g.seedSet2.clear();

        double inf_upper = g.do_poiSelection_Heuristic_upper(arg.k);
        cout<<"POI super-seeds(by HopBased Heuristic):"<<endl;
        g.printPOISelectionInfo();
        cout<<"inf_upper(by HopBased)="<<inf_upper<<endl;


        return inf_upper;
    }

//fix
    RISInfo minePOIsbyHybrid(int b,int& iteration, int& rr_num, vector<int>& pois)
    {
        //arg.epsilon = 0.05;
        //Timer t(2, "step2");
        //ASSERT(OPT_prime > 0);
        printf("****begin minePOIsbyHybrid!\n");
        int targetSize = b; //arg.k;
        arg.k =b;
        InfGraph &g = *infGraph;
        g.poiSet.clear();
        g.seedSet.clear();

        double e = exp(1);
        double approx = 1 - 1.0/e;
        int nc = social_graph->nc;
        int np = social_graph->np;
        double delta = 1.0 / (2*nc);
        double alpha = sqrt(log(6.0/delta));
        double beta = sqrt((1-1/e) * (Math::logcnk(np, arg.k) + log(6.0/delta)));


        const auto numRbase = size_t(2.0 * Math::pow2((1 - 1 / e) * alpha + beta));
        const auto maxNumR = size_t(2.0 * nc * Math::pow2((1 - 1 / e) * alpha + beta) / 1 / Math::pow2(arg.epsilon)) + 1;//size_t(2.0 * nc * Math::pow2((1 - 1 / e) * alpha + beta) / targetSize / Math::pow2(arg.epsilon)) + 1;
        auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
        const double a1 = log(numIter * 3.0 / delta);
        const double a2 = log(numIter * 3.0 / delta);
        //1. 获得最佳
        //double local_influence = g.do_poiSelection_UB_ByOneHop(arg.k);

        g.init_Two_hyper_graph_poi_out(np); // intial R1, R2;

        //numIter = (size_t)11;
        g.do_poiSelection_UB_ByOneHop(arg.k);
        //g.do_poiSelection_UB(arg.k);
        double local_inf_upper = g.local_influence_upper;
        RISInfo info;
        for(int round = 0; round<numIter;round++){
#ifdef POI_SELECTION_TRACK
            cout<<"-------------------------round :"<<round<<" ------------------------"<<endl;
#endif

            const auto numR = numRbase << round;
#ifdef POI_SELECTION_TRACK
            cout<<"Hybrid RR num="<<numR<<endl;
#endif

            g.build_hyper_graph_r_poi_outer(numR, arg); // R1
            //cout<<"buld R1"<<endl;
            g.build_hyper_graph_vld_r_poi_outer(numR, arg); // R2
            //cout<<"buld R2"<<endl;

            int numRRsets = g.get_RR_POI_sets_size();

            ////得到 在 R1上的seeds情况，并得到 Upper (local_real + outer_estimated)
            Inf_Pair infSelf_overall = g.do_poiSelection_UB_ByHybrid(arg.k); // using R1 to obtain S and upper bound of I(S)
            //g.do_poiSelection_UB_OUT(arg.k);

            double inf_Outer_R1_cover = g.overall_influence_upperDetails.outer_inf;//outer_influence_only_upper; //obtain Cove(R2)
            double inf_local_R1 = g.overall_influence_upperDetails.local_inf;
            double inf_Outer_upper = (nc / numR)*Math:: pow2(sqrt(inf_Outer_R1_cover + a2 / 2.0) + sqrt(a2 / 2.0));
            double upperbound_overall_tight = inf_local_R1 + inf_Outer_upper ;



            vector<int> cur_pois = g.poiSet;
#ifdef POI_SELECTION_TRACK
            cout<<"此时poi seed为："<<endl;
            g.printPOISelectionInfo();
#endif

            ////在R2下对得到的seeds进行验证，并得到 lower (local_R2(same to loca_R1) , outer_estimated(注意这里应是 coverage_o(R2,Sp)))
            //cout<<"check_RR_POI_Cover_OUT..."<<endl;
            double inf_Outer_R2_cover = g.check_RR_POI_Cover_OUT(cur_pois); //obtain Cove(R2),即outer_estimated
            //cout<<"finish check_RR_POI_Cover_OUT!"<<endl;
            double inf_local_R2 = infSelf_overall.local_inf;
            double inf_Outer_lower = (nc/numR)*(Math::pow2(sqrt(inf_Outer_R2_cover + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0);
            double lowerbound_overall = inf_local_R2 + inf_Outer_lower ;

            const auto approxOPIMC_tight = lowerbound_overall / upperbound_overall_tight;


#ifdef POI_SELECTION_TRACK

            cout<<"inf_Outer_R2_cover(vld)= "<<inf_Outer_R2_cover<<", inf_local_R2 ="<<inf_local_R2<<", lowerbound_overall="<<lowerbound_overall<<",其中:exact_local="<<inf_local_R2<<", outer_lower="<<inf_Outer_lower<<endl;
            cout<<"inf_Outer_R1_cover(select)="<<inf_Outer_R1_cover<<", inf_local_R1="<<inf_local_R1<<", upperbound_overall_tight ="<<upperbound_overall_tight<<",其中:exact_local="<<inf_local_R1<<", outer_upper="<<inf_Outer_upper<<endl;
            cout<<"approxOPIMC_tight:"<<approxOPIMC_tight<<endl;
#endif

            if(approxOPIMC_tight > approx-arg.epsilon){
                cout<<"accuracy guarantees ="<<(approx-arg.epsilon)<<endl;
                cout<<"find satisable seeds selection:"<<endl;
                for(int p:cur_pois){
                    cout<<"poi"<<p<<", ";
                }
                cout<<endl;
                cout<<"number of round:"<<round<<endl;
                cout<<"number of rr set:"<<numR<<endl;
                cout<<"expected inf= "<<(infSelf_overall.outer_inf*(nc/numR)+infSelf_overall.local_inf)<<endl;

                info.iteration =round; iteration = round;
                info.rrset_num = numR; rr_num = numR;
                info.poiSeeds = cur_pois; pois = cur_pois;

                return info;
            }

        }

        info.iteration = numIter; iteration = numIter;
        info.rrset_num = maxNumR; rr_num = maxNumR;
        info.poiSeeds = g.poiSet; pois = g.poiSet;

        return info;
    }




   double MCSimulation4Seeds(int r){
#ifdef  SEED
        infGraph->MC_Evaluation(r);
#else
       return infGraph->MC_Evaluation_POI(r);
#endif
    }

    double MCSimulation4GivenSeeds(int r, set<int>& seedSet){
        return  infGraph->MC_Evaluation_GivenSeeds(r, seedSet);
    }


    void Pre_Sampling_MC(int r){
        infGraph->pre_sampling_FIS(r);
    }




};