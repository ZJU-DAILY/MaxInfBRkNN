//
// Created by jin on 20-2-27.
//

#ifndef MAXINFBRGSTQ_IMM_H

#define MAXINFBRGSTQ_IMM_H

#define INFO(...);


#include "math.h"
#include "infgraph.h"
#include "SocialInfluenceModule/SocialSimulator.h"

class Math{
public:
    static double log2(int n){
        return log(n) / log(2);
    }
    static double logcnk(int n, int k) {
        double ans = 0;
        for (int i = n - k + 1; i <= n; i++)
        {
            ans += log(i);
        }
        for (int i = 1; i <= k; i++)
        {
            ans -= log(i);
        }
        return ans;
    }
    static double pow2(const double x){
        return x*x;
    }
};


enum Problem{ IM, MaxInfBRGST} ;

class Imm{
public:
    InfGraph* infGraph  = NULL;
    SocialGraph* social_graph = NULL;
    Argument arg;
    double sqr(double v){
        return v*v;
    }

    Imm(AdjListSocial& socialGraph, int u_num, int p_num, Problem p){
        //int i = num;
        social_graph = new  SocialGraph(u_num);
        if(p== MaxInfBRGST){
            infGraph = new InfGraph(u_num, p_num,social_graph);
        }
        else{
            infGraph = new InfGraph(u_num,social_graph);
        }

        arg.epsilon = 0.1;
        arg.k = 5;
        arg.model = string("IC");

    }
    ~Imm(){
        if(social_graph!=NULL)
            delete social_graph;
        if(infGraph!=NULL)
            delete infGraph;
    }

    void addSocialLinkIMM(int u, int v, double prob){
        infGraph->add_edge(u,v,prob);
    }

    //mark the reachable user nodes from poi
    void checkReachable(vector<int> pois, BatchResults& results){
        //将所有兴趣的点，与之关联的反近邻用户  关联
        for(int p: pois){
            vector<ResultDetail> rrs = results[p];
            if(rrs.size()>0){
                for(ResultDetail r: rrs){
                    int u = r.usr_id;
                    infGraph->setRearchable(u,p);
                }
            }
        }
        social_graph->traFromSource();
        //查看 source user 的个数
        int rearch_size = social_graph->getNs();
    }


    double step1()
    {
        InfGraph &g = *infGraph;
        double epsilon_prime = arg.epsilon * sqrt(2);
        //Timer t(1, "step1");
        for (int x = 1; ; x++)
        {
            int n = social_graph->n;
            int ci = (2+2/3 * epsilon_prime)* ( log(n) + Math::logcnk(n, arg.k) + log(Math::log2(n))) * pow(2.0, x) / (epsilon_prime* epsilon_prime);
            g.build_hyper_graph_r(ci, arg);

            double ept = g.build_seedset(arg.k);
            //double estimate_influence = ept * g.n;
            //INFO(x, estimate_influence);

            double  bb = (1+epsilon_prime)/pow(2.0,x);

            cout<<"ept="<<ept<<",bb="<<bb<<",ci="<<ci<<",x="<<x<<endl;
            //if (ept > 1 / pow(2.0, x))  //(1+epsilon_prime)/pow(2.0,x)
            if(ept > bb)
            {
                double OPT_prime = ept * g.social_graph->n / (1+epsilon_prime);
                //INFO("step1", OPT_prime);
                //INFO("step1", OPT_prime * (1+epsilon_prime));
                return OPT_prime;
            }
        }
        //ASSERT(false);
        return -1;
    }

    double step2(double OPT_prime)
    {
        //Timer t(2, "step2");
        //ASSERT(OPT_prime > 0);
        InfGraph &g = *infGraph;
        double e = exp(1);
        int n = g.social_graph->n;
        double alpha = sqrt(log(n) + log(2));
        double beta = sqrt((1-1/e) * (Math::logcnk(n, arg.k) + log(n) + log(2)));

        int R = 2.0 * n *  sqr((1-1/e) * alpha + beta) /  OPT_prime / arg.epsilon / arg.epsilon ;
        cout<<"total need "<<R<<" RR set for IM"<<endl;
        g.build_hyper_graph_r(R, arg);

        double opt = g.build_seedset(arg.k)*n;

        return opt;
    }

    void InfluenceMaximize()
    {
        InfGraph &g = *infGraph;
        //Timer t(100, "InfluenceMaximize(Total Time)");

        //step 1 diyi
        printf("########## Step1 begin ##########\n");
        // debugging mode lalala
        double OPT_prime;
        OPT_prime = step1(); //estimate OPT  guji
        double opt_lower_bound = OPT_prime;
        cout<<"opt_lower_bound="<<opt_lower_bound<<endl;
        printf("########## Step1 finish ##########\n");

        //step 2 dier
        printf("########## Step2 begin ##########\n");  //xuanseed
        double influence_est = step2(OPT_prime);
        printSeedSet();
        cout<<"estimated_influence = "<<influence_est<<endl;
        printf("########## Step2 finish ##########\n");

    }

    void printSeedSet(){
        cout<<"seeds users:"<<endl;
        for(int u_id: infGraph -> seedSet){
            cout<<"u"<<u_id<<endl;
        }
    }

    void printPOISet(){
        cout<<"poi selection:"<<endl;
        for(int p_id: infGraph -> poiSet){
            cout<<"p"<<p_id<<endl;
        }
    }

    void MCSimulation4Seeds(int r){
        infGraph->MC_Evaluation(r);
    }

    void MCSimulation4POIs(int r){
        infGraph->MC_Evaluation_POI(r);
    }


    double sampling()
    {
        InfGraph &g = *infGraph;
        double epsilon_prime = arg.epsilon * sqrt(2);
        //Timer t(1, "step1");
        for (int x = 1; ; x++)
        {
            //ci = almta / x(=n/2i)
            int ci = (2+2/3 * epsilon_prime)* ( log(g.social_graph->nc) + Math::logcnk(g.social_graph->np, arg.k) + log(Math::log2(g.social_graph->nc))) * pow(2.0, x) / (epsilon_prime* epsilon_prime);
            g.build_hyper_graph_r_poi(ci, arg);

            double ept = g.do_poiSelection(arg.k);
            //double estimate_influence = ept * g.n;
            //INFO(x, estimate_influence);

            double  bb = (1+epsilon_prime)/pow(2.0,x);

            cout<<"ept="<<ept<<",bb="<<bb<<",ci="<<ci<<",x="<<x<<endl;
            //if (ept > 1 / pow(2.0, x))  //(1+epsilon_prime)/pow(2.0,x)
            if(ept > bb)
            {
                double OPT_prime = ept * g.social_graph->nc / (1+epsilon_prime);
                //INFO("step1", OPT_prime);
                //INFO("step1", OPT_prime * (1+epsilon_prime));
                return OPT_prime;
            }
        }
        //ASSERT(false);
        return -1;
    }

    double poiSelection(double OPT_prime)
    {
        //Timer t(2, "step2");
        //ASSERT(OPT_prime > 0);
        InfGraph &g = *infGraph;
        double e = exp(1);
        int v_count = g.social_graph->getNc();
        int nc =  g.social_graph->nc;

        double alpha = sqrt(log(nc) + log(2));
        double beta = sqrt((1-1/e) * (Math::logcnk(nc, arg.k) + log(nc) + log(2)));

        //int R = 2.0 * g.n *  sqr((1-1/e) * alpha + beta) /  OPT_prime / arg.epsilon / arg.epsilon ;

        int R = 2.0 * nc *  sqr((1-1/e) * alpha + beta) /  OPT_prime / arg.epsilon / arg.epsilon ;

        cout<<"total need "<<R<<" RR set for IM"<<endl;
        g.build_hyper_graph_r_poi(R, arg);

        double opt = g.do_poiSelection(arg.k)*nc;

        return opt;
    }

    void minePOIsIMM(){
        InfGraph &g = *infGraph;
        //Timer t(100, "InfluenceMaximize(Total Time)");

        //step 1 diyi
        printf("########## Step1: Sampling begin ##########\n");
        // debugging mode lalala
        double OPT_estimate;
        OPT_estimate = sampling(); //estimate OPT  guji
        double opt_lower_bound = OPT_estimate;
        cout<<"opt_lower_bound(poi)="<<opt_lower_bound<<endl;
        printf("########## Step1: Sampling finish ##########\n");

        //step 2 dier
        printf("########## Step2: POI Selecting begin ##########\n");  //xuanseed
        double influence_est = poiSelection(opt_lower_bound);
        printPOISet();
        cout<<"estimated_influence = "<<influence_est<<endl;
        printf("########## Step2: POI Selecting finish ##########\n");

    }


};



#endif //MAXINFBRGSTQ_IMM_H
