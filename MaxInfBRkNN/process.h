#ifndef   MAXINFRKGSKQ_PORCESS


//#pragma GCC optimize(2)
//#pragma G++ optimize(2)

#include <iostream>
#include<sys/time.h>
#include "dijisktra.h"
#include "query_plus.h"
#include "functional.h"
#include "diskbased.h"
#include "gim_tree.h"
#include "test.h"
#include "unistd.h"
#include <getopt.h>
#include <cstdlib>
#include "InfluenceModels.h"



//influence evaluation
struct timeval iv;
long long is, ie;
#define SOCIAL_START gettimeofday( &iv, NULL ); is = iv.tv_sec * 1000000 + iv.tv_usec;
#define SOCIAL_END gettimeofday( &iv, NULL ); ie = iv.tv_sec * 1000000 + iv.tv_usec;
#define SOCIAL_PRINT(T) printf("%s runtime: %lld (ms)\r\n", (#T), (ie - is)/1000 );


using namespace std;

//将G-Tree中各节点对应子图的距离上下界，作为预处理结果提前计算

void initialPHL(){

    LOAD_START
    //stringstream phlFilePath;
    //phlFilePath << DATASET_HOME <<"road/"<<dataset_name<<".phl";
    cout<<"load PHL indexing... ";
    string phl_fileName = getRoadInputPath(PHL);
    phl.LoadLabel(phl_fileName.c_str());

    LOAD_END
    LOAD_PRINT()
}



MinMaxD getMinMaxDistanceBetweenLeaves_phl(int leaf1, int leaf2){
    MinMaxD mm;
    vector<int> bs1 = GTree[leaf1].borders;
    vector<int> bs2 = GTree[leaf2].borders;
    int b1 = -1; int b2 = -1;
    float maxDis = 0;
    float minDis = 9999999999;

    for(int _b1:bs1){
        for(int _b2: bs2){
           ///cout<<"border:v"<<_b1<<",v"<<_b2<<"=";
           float distance = phl.Query(_b1, _b2)/ WEIGHT_INFLATE_FACTOR;
           ///cout<<distance<<endl;
            if(distance < minDis) {
                minDis = distance; b1 = _b1; b2 = _b2;
            }
        }
    }

    float maxDis1 = 0;
    vector<int> distance_results;
    vector<int> cands1;

    for(int v :GTree[leaf1].vetexSet){
        ///cout<<"dist:b"<<b1<<",v"<<v<<"=";
        float distance_normalized1 = phl.Query(b1, v)/ WEIGHT_INFLATE_FACTOR;
        ///cout<<distance_normalized1<<endl;
        maxDis1 = max(distance_normalized1,maxDis1);

    }

    float maxDis2 = 0;

    for(int  u: GTree[leaf2].vetexSet){
        ///cout<<"dist2:u"<<u<<",b"<<b2<<"=";
        float distance_normalized2 = phl.Query(b2, u)/ WEIGHT_INFLATE_FACTOR;
        ///cout<<distance_normalized2<<endl;
        maxDis2 = max(distance_normalized2,maxDis2);
        ///cout<<"maxDis2="<<maxDis2<<endl;

    }

    mm.min = minDis;
    mm.max = minDis+maxDis1+maxDis2;
    return  mm;

}

MinMaxD  getMinMaxDistanceWithinLeaf_phl(int leaf){
    MinMaxD mm;
    float maxDis = 0;
    float minDis = 10000000;

    vector<int> vertexes;
    for(int v :GTree[leaf].vetexSet)
        vertexes.push_back(v);
    //omp_set_num_threads(4);
//#pragma omp parallel for
    for(int i=0; i<vertexes.size();i++){
        int s = vertexes[i];
        for(int j=0; j<vertexes.size();j++){
            int d = vertexes[j];
            float distance_normalized = phl.Query(s,d);
            maxDis = max(maxDis,distance_normalized);
        }
    }

    mm.min = 0;
    mm.max = maxDis;
    return  mm;
}



//Pair **LeafBorderPair;

void outputNodeDisBound_Enhance_other(){  //经过phl索引加速
    TIME_TICK_START
    cout<<"output disbound, map:"<<dataset_name<<endl;
    ofstream outputBound;
    //stringstream tt;
    //tt<<FILE_DISBOUND;
    outputBound.open(getRoadInputPath(DISBOUND));
    vector<int> leafNodes;
    vector<int> UpperNodes;

    //map<int,map<int,MinMaxD>> DisBoundMap;
    int node_size = GTree.size();
    MinMaxD** DisBoundMap = new MinMaxD*[node_size];
    for(int i=0;i<node_size;i++){
        DisBoundMap[i] = new MinMaxD[node_size];
    }


    //map<int,map<int,bool>> DisBOOLMap;
    //map<int,map<int,bool >> DisBoundMap;
    //取出叶子节点
    for(int i=1; i<GTree.size(); i++){
        if(GTree[i].isleaf){
            leafNodes.push_back(i);
            cout<<"加入叶节点"<<i<<endl;
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
    //getchar();
    //先对各叶子节点间的距离上下界进行计算
    cout<<"对各叶节点间距离上下界进行计算"<<endl;
    int count = 1;
    for(int idx=0; idx<leafNodes.size();idx++){
        //omp_set_num_threads(2);
//#pragma omp parallel for
        //if(leafNodes[idx]<3178) continue;

        for(int idx2=0; idx2<leafNodes.size();idx2++){
            //cout<<"for "<<i<<","<<j<<endl;
            MinMaxD mm; MinMaxD mm2;
            int i= leafNodes[idx];
            int j = leafNodes[idx2];
            if(i==17&&j==57){
                cout<<"i==17&&j==57"<<endl;
            }
            cout<<"计算(leaf"<<i<<",leaf"<<j<<")=";
            if(i>j){
                cout<<"(leaf"<<i<<",leaf"<<j<<")=";
                cout<<"已计算过"<<endl;
                continue;
            }
            else{


                if(i==j){  //内部距离上下界计算
                    mm = getMinMaxDistanceWithinLeaf_phl(i);//getMinMaxDistanceWithinLeaf_slow(i);
                    DisBoundMap[i][i] = mm;
                    //DisBOOLMap[i][i] = true;
                }
                else{  //跨节点距离上下界计算
                    mm = getMinMaxDistanceBetweenLeaves_phl(i,j); //getMinMaxDistanceBetweenLeaves_slow(i,j);
                    DisBoundMap[i][j] = mm;
                    //DisBOOLMap[i][j] = true;
                    DisBoundMap[j][i] = mm;
                    //DisBOOLMap[j][i] = true;
                }
                outputBound << i << ' ' <<j<< ' ' << mm.max <<' '<< mm.min<< endl;
                cout<<"gtree:(max:"<<mm.max<<",min:"<<mm.min<<")"<<endl;
                //cout<<"_phl_:(max:"<<mm2.max<<",min:"<<mm2.min<<")"<<endl;
                if(count%5==0){
                    count++;
                    getchar();
                }

            }



        }
        cout<<"---------------------------------------------------"<<endl;
    }
    cout<<"完成所有叶节点间距离bound的计算"<<endl;

    //对叶节点与上层节点间的距离上下界进行计算
    cout<<"对叶节点与上层节点间距离上下界进行计算"<<endl;
    for(int l: leafNodes){
        //omp_set_num_threads(4);
//#pragma omp parallel for
        for(int idx = 0; idx<UpperNodes.size();idx++){
            int  n = UpperNodes[idx];
            float min= 9999999999; float max=0;
            for(int nl: GTree[n].leafnodes){   //上层节点的所有叶子节点
                MinMaxD mm = DisBoundMap[l][nl];
                if(mm.max > max) max = mm.max;
                if(mm.min < min) min = mm.min;
            }
            MinMaxD mm;
            mm.max = max; mm.min = min;
            DisBoundMap[l][n] = mm;

            outputBound << l << ' ' <<n<< ' ' << mm.max <<' '<< mm.min<< endl;
            cout<<"(leaf"<<l<<",node"<<n<<")=(max:"<<mm.max<<",min:"<<mm.min<<")"<<endl;
        }

    }
    //对上层节点间的距离上下界进行计算
    cout<< "对上层节点间距离上下界进行计算"<<endl;
    for(int n1: UpperNodes){
        for(int n2: UpperNodes){
            MinMaxD mm;
            //已计算
            if(n1>n2){
                mm = DisBoundMap[n1][n2];
                cout<<"(node"<<n1<<",node"<<n2<<")=已计算"<<endl;
                continue;
            }
                //未被计算
            else{
                float min= 9999999; float max=0;
                for(int nl1:GTree[n1].leafnodes){
                    //for(int nl2:GTree[n2].leafnodes){
                    MinMaxD mm_temp = DisBoundMap[nl1][n2];
                    if(mm_temp.max > max) max = mm_temp.max;
                    if(mm_temp.min < min) min = mm_temp.min;
                }
                mm.max = max; mm.min = min;
                DisBoundMap[n1][n2] = mm;
                DisBoundMap[n2][n1] = mm;
                outputBound << n1 << ' ' <<n2<< ' ' << mm.max <<' '<< mm.min<< endl;
                cout<<"(node"<<n1<<",node"<<n2<<")=(max:"<<mm.max<<",min:"<<mm.min<<")"<<endl;
            }

        }
    }




    cout<<"释放内存"<<endl;
    //释放堆内存
    for(int i=0;i<node_size;i++){
        delete [] DisBoundMap[i];
    }
    delete [] DisBoundMap;

    outputBound.close();
    TIME_TICK_END
    TIME_TICK_PRINT("outputNodeDisBound_Enhance new time:");

}



void outputNodeDisBound_Enhance_Gowalla(){
    TIME_TICK_START
    cout<<"output disbound, map:"<<dataset_name<<endl;
    ofstream outputBound;
    stringstream tt;
    //tt<<FILE_DISBOUND;
    //outputBound.open(tt.str());

    vector<int> leafNodes;
    vector<int> UpperNodes;

   /* cout<<"getMinMaxDistanceBetweenLeaves_phl(25,10090)..."<<endl;
    MinMaxD mm = getMinMaxDistanceBetweenLeaves_phl(25,10090);
    cout<<"(leaf"<<25<<",leaf"<<10090<<")=(max:"<<mm.max<<",min:"<<mm.min<<")"<<endl;
    getchar();*/

    //map<int,map<int,MinMaxD>> DisBoundMap;
    int node_size = GTree.size();
    MinMaxD** DisBoundMap = new MinMaxD*[node_size];
    for(int i=0;i<node_size;i++){
        DisBoundMap[i] = new MinMaxD[node_size];
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
    if(true){
        outputBound.open(getRoadInputPath(DISBOUND));
        //先对各叶子节点间的距离上下界进行计算
        cout<<"对各叶节点间距离上下界进行计算"<<endl;
        for(int idx=0; idx<leafNodes.size();idx++){
            //omp_set_num_threads(2);
//#pragma omp parallel for
            //if(leafNodes[idx]<3178) continue;

            for(int idx2=idx; idx2<leafNodes.size();idx2++){
                MinMaxD mm;
                int i= leafNodes[idx];
                int j = leafNodes[idx2];
                cout<<"for "<<i<<","<<j<<endl;
                if(i>j){
                    cout<<"(leaf"<<i<<",leaf"<<j<<")=";
                    cout<<"已计算过"<<endl;
                    continue;
                }
                else{

                    if(i==j){  //内部距离上下界计算
                        mm = getMinMaxDistanceWithinLeaf_phl(i);
                        DisBoundMap[i][i] = mm;
                        //DisBOOLMap[i][i] = true;
                    }
                    else{  //跨节点距离上下界计算
                        mm = getMinMaxDistanceBetweenLeaves_phl(i,j);
                        DisBoundMap[i][j] = mm;
                        //DisBOOLMap[i][j] = true;
                        DisBoundMap[j][i] = mm;
                        //DisBOOLMap[j][i] = true;
                    }
                    outputBound << i << ' ' <<j<< ' ' << mm.max <<' '<< mm.min<< endl;
                    cout<<"(leaf"<<i<<",leaf"<<j<<")=(max:"<<mm.max<<",min:"<<mm.min<<")"<<endl;
                }



            }
        }
        cout<<"完成所有叶节点间距离bound的计算"<<endl;
    }
    //先load 各叶子节点间的距离上下界,在计算上层节点bound
    if(true){
        //string dis_path = outputBound.open();
        //输出目录
        stringstream ott;
        ott<<getRoadInputPath(DISBOUND)<<"-leafupper";
        outputBound.open(ott.str());
        //for(int i=1;i<2;i++){
            cout<<"LOADING DISBOUND from ";
            ifstream inDisBoundFile;  ////叶节点间bound文件目录
            string  in_Path;
            in_Path = getRoadInputPath(DISBOUND);//<<i;
            cout<<in_Path<<"..."<<endl;
            inDisBoundFile.open(in_Path);
            string str;
            while(getline(inDisBoundFile,str)){
                istringstream tt2(str);
                int leaf1,leaf2;
                float maxDis, minDis;
                tt2>>leaf1>>leaf2>>maxDis>>minDis;
                DisBoundMap[leaf1][leaf2].max = maxDis;
                DisBoundMap[leaf1][leaf2].min = minDis;
                DisBoundMap[leaf2][leaf1].max = maxDis;
                DisBoundMap[leaf2][leaf1].min = minDis;
            }
            inDisBoundFile.close();
            cout<<"part0"<<" LOADING COMPLETE!"<<endl;
        //}

        //上层节点bound计算
        //对叶节点与上层节点间的距离上下界进行计算
        cout<<"对叶节点与上层节点间距离上下界进行计算"<<endl;

        for(int l: leafNodes){
            //omp_set_num_threads(4);
//#pragma omp parallel for
            for(int idx = 0; idx<UpperNodes.size();idx++){
                int  n = UpperNodes[idx];
                float min= 9999999999; float max=0;
                for(int nl: GTree[n].leafnodes){   //上层节点的所有叶子节点
                    MinMaxD mm = DisBoundMap[l][nl];
                    if(mm.max > max) max = mm.max;
                    if(mm.min < min) min = mm.min;
                }
                MinMaxD mm;
                mm.max = max; mm.min = min;
                DisBoundMap[l][n] = mm;

                outputBound << l << ' ' <<n<< ' ' << mm.max <<' '<< mm.min<< endl;
                cout<<"(leaf"<<l<<",node"<<n<<")=(max:"<<mm.max<<",min:"<<mm.min<<")"<<endl;
            }


        }
        outputBound.close();


        //对上层节点间的距离上下界进行计算
        stringstream tt;
        tt<<getRoadInputPath(DISBOUND)<<"-upper";
        outputBound.open(tt.str());
        cout<< "对上层节点间距离上下界进行计算"<<endl;
        for(int n1: UpperNodes){
            for(int n2: UpperNodes){
                MinMaxD mm;
                //已计算
                if(n1>n2){
                    mm = DisBoundMap[n1][n2];
                    cout<<"(node"<<n1<<",node"<<n2<<")=已计算"<<endl;
                    continue;
                }
                    //未被计算
                else{
                    float min= 9999999; float max=0;
                    for(int nl1:GTree[n1].leafnodes){
                        //for(int nl2:GTree[n2].leafnodes){
                        MinMaxD mm_temp = DisBoundMap[nl1][n2];
                        if(mm_temp.max > max) max = mm_temp.max;
                        if(mm_temp.min < min) min = mm_temp.min;
                    }
                    mm.max = max; mm.min = min;
                    DisBoundMap[n1][n2] = mm;
                    DisBoundMap[n2][n1] = mm;
                    outputBound << n1 << ' ' <<n2<< ' ' << mm.max <<' '<< mm.min<< endl;
                    cout<<"(node"<<n1<<",node"<<n2<<")=(max:"<<mm.max<<",min:"<<mm.min<<")"<<endl;
                }

            }
        }

    }
    outputBound.close();

    cout<<"释放内存"<<endl;
    //释放堆内存
    for(int i=0;i<node_size;i++){
        delete [] DisBoundMap[i];
    }
    delete [] DisBoundMap;

    TIME_TICK_END
    TIME_TICK_PRINT("outputNodeDisBound_Enhance new time:");

}


void outputNodeDisBound_Enhance_Twitter(){
    TIME_TICK_START
    cout<<"output disbound, map:"<<dataset_name<<endl;
    ofstream outputBound;
    stringstream tt;
    //tt<<FILE_DISBOUND;
    //outputBound.open(tt.str());

    vector<int> leafNodes;
    vector<int> UpperNodes;


    //map<int,map<int,MinMaxD>> DisBoundMap;
    int node_size = GTree.size();
    MinMaxD** DisBoundMap = new MinMaxD*[node_size];
    for(int i=0;i<node_size;i++){
        DisBoundMap[i] = new MinMaxD[node_size];
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
    if(true){
        outputBound.open(getRoadInputPath(DISBOUND),ios::app);
        //先对各叶子节点间的距离上下界进行计算
        cout<<"对各叶节点间距离上下界进行计算"<<endl;
        for(int idx=0; idx<leafNodes.size();idx++){
            //omp_set_num_threads(2);
//#pragma omp parallel for
            //if(leafNodes[idx]<3178) continue;

            for(int idx2=idx; idx2<leafNodes.size();idx2++){
                MinMaxD mm;
                int i= leafNodes[idx];
                int j = leafNodes[idx2];
                cout<<"for "<<i<<","<<j<<endl;
                if(i>j){
                    cout<<"(leaf"<<i<<",leaf"<<j<<")=";
                    cout<<"已计算过"<<endl;
                    continue;
                }
                else{

                    if(i==j){  //内部距离上下界计算
                        mm = getMinMaxDistanceWithinLeaf_phl(i);
                        DisBoundMap[i][i] = mm;
                        //DisBOOLMap[i][i] = true;
                    }
                    else{  //跨节点距离上下界计算
                        mm = getMinMaxDistanceBetweenLeaves_phl(i,j);
                        DisBoundMap[i][j] = mm;
                        //DisBOOLMap[i][j] = true;
                        DisBoundMap[j][i] = mm;
                        //DisBOOLMap[j][i] = true;
                    }
                    outputBound << i << ' ' <<j<< ' ' << mm.max <<' '<< mm.min<< endl;
                    cout<<"(leaf"<<i<<",leaf"<<j<<")=(max:"<<mm.max<<",min:"<<mm.min<<")"<<endl;
                }



            }
        }
        cout<<"完成所有叶节点间距离bound的计算"<<endl;
    }
    //先load 各叶子节点间的距离上下界,在计算上层节点bound
    if(true){
        //string dis_path = outputBound.open();
        //输出目录
        stringstream ott;
        ott<<getRoadInputPath(DISBOUND)<<"-leafupper";
        outputBound.open(ott.str());
        //for(int i=1;i<2;i++){
        cout<<"LOADING DISBOUND from ";
        ifstream inDisBoundFile;  ////叶节点间bound文件目录
        string  in_Path;
        in_Path = getRoadInputPath(DISBOUND);//<<i;
        cout<<in_Path<<"..."<<endl;
        inDisBoundFile.open(in_Path);
        string str;
        while(getline(inDisBoundFile,str)){
            istringstream tt2(str);
            int leaf1,leaf2;
            float maxDis, minDis;
            tt2>>leaf1>>leaf2>>maxDis>>minDis;
            DisBoundMap[leaf1][leaf2].max = maxDis;
            DisBoundMap[leaf1][leaf2].min = minDis;
            DisBoundMap[leaf2][leaf1].max = maxDis;
            DisBoundMap[leaf2][leaf1].min = minDis;
        }
        inDisBoundFile.close();
        cout<<"part0"<<" LOADING COMPLETE!"<<endl;
        //}

        //上层节点bound计算
        //对叶节点与上层节点间的距离上下界进行计算
        cout<<"对叶节点与上层节点间距离上下界进行计算"<<endl;

        for(int l: leafNodes){
            //omp_set_num_threads(4);
//#pragma omp parallel for
            for(int idx = 0; idx<UpperNodes.size();idx++){
                int  n = UpperNodes[idx];
                float min= 9999999999; float max=0;
                for(int nl: GTree[n].leafnodes){   //上层节点的所有叶子节点
                    MinMaxD mm = DisBoundMap[l][nl];
                    if(mm.max > max) max = mm.max;
                    if(mm.min < min) min = mm.min;
                }
                MinMaxD mm;
                mm.max = max; mm.min = min;
                DisBoundMap[l][n] = mm;

                outputBound << l << ' ' <<n<< ' ' << mm.max <<' '<< mm.min<< endl;
                cout<<"(leaf"<<l<<",node"<<n<<")=(max:"<<mm.max<<",min:"<<mm.min<<")"<<endl;
            }


        }
        outputBound.close();


        //对上层节点间的距离上下界进行计算
        stringstream tt;
        tt<<getRoadInputPath(DISBOUND)<<"-upper";
        outputBound.open(tt.str());
        cout<< "对上层节点间距离上下界进行计算"<<endl;
        for(int n1: UpperNodes){
            for(int n2: UpperNodes){
                MinMaxD mm;
                //已计算
                if(n1>n2){
                    mm = DisBoundMap[n1][n2];
                    cout<<"(node"<<n1<<",node"<<n2<<")=已计算"<<endl;
                    continue;
                }
                    //未被计算
                else{
                    float min= 9999999; float max=0;
                    for(int nl1:GTree[n1].leafnodes){
                        //for(int nl2:GTree[n2].leafnodes){
                        MinMaxD mm_temp = DisBoundMap[nl1][n2];
                        if(mm_temp.max > max) max = mm_temp.max;
                        if(mm_temp.min < min) min = mm_temp.min;
                    }
                    mm.max = max; mm.min = min;
                    DisBoundMap[n1][n2] = mm;
                    DisBoundMap[n2][n1] = mm;
                    outputBound << n1 << ' ' <<n2<< ' ' << mm.max <<' '<< mm.min<< endl;
                    cout<<"(node"<<n1<<",node"<<n2<<")=(max:"<<mm.max<<",min:"<<mm.min<<")"<<endl;
                }

            }
        }

    }
    outputBound.close();

    cout<<"释放内存"<<endl;
    //释放堆内存
    for(int i=0;i<node_size;i++){
        delete [] DisBoundMap[i];
    }
    delete [] DisBoundMap;

    TIME_TICK_END
    TIME_TICK_PRINT("outputNodeDisBound_Enhance new time:");

}


void outputNodeDisBound_Enhance(){
#ifdef Gowalla
    outputNodeDisBound_Enhance_Gowalla();
#else
    #ifdef Twitter
        outputNodeDisBound_Enhance_Twitter();
    #else
        outputNodeDisBound_Enhance_other();
    #endif
#endif
}






int do_some_preprocess(int Qk, int a, int Qk_values[], float a_values[], float alpha){
    //0.预处理(输出所有user的top-k结果)
    if(1){
        outputUserScore(Qk, a, alpha);
        /*for(int Qk_s: Qk_values)
            for(float a_s: a_values)
                outputPOIRkNNResults(Qk_s,a_s, alpha);
        outputPOIRkNNResults(Qk,a, alpha);*/
        return 0;
    }
    if(0){
        outputNodeDisBound_Enhance();
        return 0;
    }
}


void initialIOCache(int need_nvd_cache){
    //ListConfig(setting);
    //设定输入数据的文件目录
    ConfigType setting;
    AddConfigFromFile(setting,"../map_exp");
    string filename="../../../data/";
    filename=filename+dataset_name+"/"+road_data;
    const char* fileprefix=filename.c_str();
    //获取参数信息中的cache页数
    int _cachepages=getConfigInt("cachepages",setting,false,DEFAULT_CACHESIZE);


    srand((unsigned) time(NULL));
    //打开磁盘上数据文件
    if(need_nvd_cache ==0)
        OpenDiskComm(fileprefix,_cachepages);	// open the files
    else
        OpenDiskComm_Plus_NVD(fileprefix,_cachepages);	// open the files
    cout<<"磁盘文件打开！"<<endl;
    //为term weight table 开辟缓存
    for(int i=1;i<=vocabularySize;i++){
        termWeight_Map.push_back(-1);
    }
    for(int i=0;i<GTree.size();i++){
        TreeNode nodeCache;
        GTreeCache.push_back(nodeCache);
    }



}

void initialIOCache_varyingO(int need_nvd_cache,int objectType, float o_ratio){
    //ListConfig(setting);
    //设定输入数据的文件目录
    ConfigType setting;
    AddConfigFromFile(setting,"../map_exp");
    string filename="../../../data/";
    if(objectType==BUSINESS){
        filename= filename + dataset_name +"/Varying_IGtree/VaryingO/";
    }
    else if(objectType==USER){
        filename= filename + dataset_name +"/Varying_IGtree/VaryingU/";
    }

    stringstream ss; ss<<(int)(o_ratio*100)<<"%/";

    filename=filename+ss.str()+road_data;
    const char* fileprefix=filename.c_str();
    //获取参数信息中的cache页数
    int _cachepages=getConfigInt("cachepages",setting,false,DEFAULT_CACHESIZE);


    srand((unsigned) time(NULL));
    //打开磁盘上数据文件
    if(need_nvd_cache ==0)
        OpenDiskComm(fileprefix,_cachepages);	// open the files
    else{
        OpenDiskComm_Plus_NVD_varyingO(fileprefix,_cachepages,objectType,o_ratio);	// open the files
    }

    cout<<"磁盘文件打开！"<<endl;
    //为term weight table 开辟缓存
    for(int i=1;i<=vocabularySize;i++){
        termWeight_Map.push_back(-1);
    }
    for(int i=0;i<GTree.size();i++){
        TreeNode nodeCache;
        GTreeCache.push_back(nodeCache);
    }



}



/*---------------------------------------将poi、路网邻接表索引数据写到磁盘-------------------------------*/
void processData(){
    //加载G-Tree索引
    ConfigType setting;
    AddConfigFromFile(setting,"../map_exp");
    //AddConfigFromCmdLine(setting,argc,argv);
    initGtree();

    initialPHL();

    // 加载 edge information
    loadEdgeMap();

    //读取社交网络图信息
    //loadFriendShipData();
    loadFriendShipBinaryData();

    //读取用户在兴趣点的check-in
    //loadCheckinData();
    loadCheckinBinaryData();


    //读取用户与兴趣点集合的路网位置与文本信息,必须最后
    TIME_TICK_START
#ifdef LasVegas
    loadUPMap_LasVegas();
#else
    #ifdef Twitter
        loadUPMap_Twitter();
    #else
        loadUPMap();
    #endif
#endif
    TIME_TICK_END
    TIME_TICK_PRINT(" loadUPMap runtime: ")
    //getchar();

    //下面开始将数据写入磁盘
    BlkLen = getBlockLength();
    cout<<"BlkLen="<<BlkLen<<endl;

    char outfileprefix[200];
    sprintf(outfileprefix,"../../../data/%s/%s",dataset_name,road_data);


    //将路网图边上的兴趣点数据存入二进制数据文件.p_d和索引文件.p_bt中
    makePtFiles(outfileprefix);

    //路网图信息存入adj二进制文件中（或者gim-tree对象）(注意：这一步一定要在makePtFiles之后)
    makeAdjListFiles(outfileprefix);

    //将路网上的双色体对象对应数据以GIM-tree“叶节点”为单位存入二进制数据文件.p_d和索引文件.p_bt中
    makeO2UFileByLeaf(outfileprefix);

    //gim-tree上双色体对象的keyword的倒排索引信息存入.p_invList二进制文件中
    //makeInvIndex_poi2Usr(outfileprefix);
    makeInvIndex_poi2Usr_new(outfileprefix);

    //将gtree节点信息存入.gtree二进制文件中
    gimtree_write2disk();

    //输出联合社交关系SM
    save_socialMatrix();

    //输出用户朋友的最大兴趣点签到数
    outputUsrMaxSocial();   //记得加


}

void processData_varyingUP(ObjectType type, float rate){
    //加载G-Tree索引
    ConfigType setting;
    AddConfigFromFile(setting,"../map_exp");
    //AddConfigFromCmdLine(setting,argc,argv);

    initGtree();

    initialPHL();

    // 加载 edge information
    loadEdgeMap();

    //读取社交网络图信息
    //loadFriendShipData();
    loadFriendShipBinaryData();

    //读取用户在兴趣点的check-in
    //loadCheckinData();
    loadCheckinBinaryData();


    //读取用户与兴趣点集合的路网位置与文本信息,必须最后
    TIME_TICK_START

    loadUPMap_varying_object(type, rate);


    TIME_TICK_END
    TIME_TICK_PRINT(" loadUPMap runtime: ")
    //getchar();

    //下面开始将数据写入磁盘
    BlkLen = getBlockLength();
    cout<<"BlkLen="<<BlkLen<<endl;

    char outfileprefix[200];
    if(type==USER){
        sprintf(outfileprefix,"../../../data/%s/Varying_IGtree/VaryingU/%d%/%s",dataset_name,(int)(rate*100),road_data);
    }
    else if(type==BUSINESS){
        sprintf(outfileprefix,"../../../data/%s/Varying_IGtree/VaryingO/%d%/%s",dataset_name,(int)(rate*100),road_data);
    }

    printf("目录为：，%s\n", outfileprefix);

    //将路网图边上的兴趣点数据存入二进制数据文件.p_d和索引文件.p_bt中
    makePtFiles(outfileprefix);

    //路网图信息存入adj二进制文件中（或者gim-tree对象）(注意：这一步一定要在makePtFiles之后)
    makeAdjListFiles(outfileprefix);

    //将路网上的双色体对象对应数据以GIM-tree“叶节点”为单位存入二进制数据文件.p_d和索引文件.p_bt中
    makeO2UFileByLeaf(outfileprefix);

    //gim-tree上双色体对象的keyword的倒排索引信息存入.p_invList二进制文件中
    //makeInvIndex_poi2Usr(outfileprefix);
    makeInvIndex_poi2Usr_new(outfileprefix);

    //将gtree节点信息存入.gtree二进制文件中
    gimtree_write2disk(outfileprefix);

    //输出联合社交关系SM
    save_socialMatrix(outfileprefix);

    //输出用户朋友的最大兴趣点签到数
    outputUsrMaxSocial(outfileprefix);   //记得加


}




int computingTerFreZipfan(int rank, int C_value){
    int frequency = C_value/rank;
    return frequency;
}

void Zipfan_POITextInfo(){ // leaf node level,去掉了之前多余累赘的 address_v2p对应的索引与数据文件
    int shi = 0;
    cout<<"Begin Zipfan_POITextInfo!"<<endl;
    LOAD_START
    // load PoiMap
    printf("LOADING POI INFO...");
    ifstream finPoi;
    string poiPath = getObjectInputPath(BUSINESS);
    string str;
    finPoi.open(poiPath.c_str());


    map<int,map<int,bool>> termPOI_Mark;
    map<int,map<int,bool>> poiTerm_Mark;
    map<int,set<int>> posting_list;
    set<int> current_words;
    while (getline(finPoi, str)) {
        istringstream tt(str);
        int POI_id, POI_Ni, POI_Nj;
        float POI_dist, POI_dis;
        tt >> POI_id >> POI_Ni >> POI_Nj >> POI_dist >> POI_dis;
        float exchangeNum = (int) (POI_dist * WEIGHT_INFLATE_FACTOR);
        POI_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (POI_dis * WEIGHT_INFLATE_FACTOR);
        POI_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;


        vector<int> uKey;
        set<int> keywordSet;
        int keyTmp;
        //cout<<"keyword:";
        int term = -1;
        while (tt >> term) {
            if (0 > term) {
                break;
            }
            if(term==0) continue;  //关键词id从1开始
            uKey.push_back(term); keywordSet.insert(term);
            current_words.insert(term);
            posting_list[term].insert(POI_id);
            //cout<<",w"<<term;
            poiTerm_Mark[POI_id][term] = true;
            termPOI_Mark[term][POI_id] = true;

        }
       // cout<<endl;

        POI tmpPoi;
        tmpPoi.id = POI_id;
        tmpPoi.Ni = min(POI_Ni,POI_Nj);
        tmpPoi.Nj = max(POI_Ni,POI_Nj);
        tmpPoi.dist = POI_dist;
        tmpPoi.dis = POI_dis;
        tmpPoi.category = 1;
        tmpPoi.keywords = uKey;
        tmpPoi.keywordSet = keywordSet;
        POIs.push_back(tmpPoi);

        poiCnt++;
    }
    finPoi.close();
    cout << " DONE! POI #: " << poiCnt << endl;
    LOAD_END
    LOAD_PRINT(" load POI ")



    int term_size = current_words.size();


    int poi_size = POIs.size();
    cout<<"poi_size="<<poi_size<<endl;
    cout<<"keyword_size="<<term_size<<endl;
    ///getchar();

    float C_ratio = 0.15;

    int C_value = C_ratio * poi_size;

    ////先对已有的keyword进行处理：
    for(int i=1;i<=term_size;i++){

        //if(i>15) break;
        int target_term = i;
        //float term_weight = getTermIDFWeight(target_term);
        int current_frequency = posting_list[target_term].size();
        cout<<"----------------------term"<<target_term<<"----------------------"<<endl;

        int rank = i;
        int supposed_frequency = computingTerFreZipfan(rank,C_value);
        cout<<"supposed_frequency ="<<supposed_frequency;
        cout<<",current_frequency ="<<current_frequency<<endl;
        if(supposed_frequency>current_frequency ){  //需要加入

            int add_num = supposed_frequency - current_frequency;
            int added_count = 0;
            while(added_count<add_num){
                int random_id = (int)random() % poi_size;  //任意取一个poi
                int _size = POIs[random_id].keywordSet.size();
                if(_size > 20) continue;
                if(poiTerm_Mark[random_id].count(target_term)){ //若这poi已经有该keyword了， 则discard
                    continue;
                }
                else{
                    POIs[random_id].keywordSet.insert(target_term);  //否则向poi中加入该keyword
                    poiTerm_Mark[random_id][target_term] = true;
                    termPOI_Mark[target_term][random_id] = true;
                    added_count++;
                }
            }

        }
        else if(supposed_frequency<current_frequency){  //需要删除
            set<int> _list = posting_list[target_term];
            vector<int> poi_ids;
            for(int id: _list)
                poi_ids.push_back(id);

            int reduced_num = current_frequency-supposed_frequency;
            int reduced_count = 0;

            for(int i=0;i<reduced_num;i++){
                int target_id = poi_ids[i];
                int _size2 = POIs[target_id].keywordSet.size();
                if(_size2<2) {
                    POIs[target_id].keywordSet.erase(target_term) ;
                    poiTerm_Mark[target_id][target_term] = false;
                    termPOI_Mark[target_term][target_id]= false;
                    int  _term = (int) random() % 55;
                    int _coin = (int)random() % 5;
                    if(_coin!=1){  //0.8的概率是前25
                        int _th = random()% 25;
                        int _coin2 = random()%2;
                        if(_coin2==0){  //1/2的概率前10
                            _term = random()% 10;
                            int  _coin3 = (int)random()%2;
                            if(_coin3==0){ //1/2概率前5
                                _term = (int)random()%5;
                            }
                        }
                        else{
                            _term = _th;
                        }

                    }

                    _term = max(1, _term);

                    POIs[target_id].keywordSet.insert(_term);
                    poiTerm_Mark[target_id][_term] = true;
                    termPOI_Mark[_term][target_id]= true;
                }
                else{
                    POIs[target_id].keywordSet.erase(target_term) ;
                    poiTerm_Mark[target_id][target_term] = false;
                    termPOI_Mark[target_term][target_id]= false;
                }
            }



        }



    }

    for(int t = term_size+1;t<=vocabularySize;t++){
        int rank = t;
        int target_term = t;
        int supposed_frequency = computingTerFreZipfan(rank,C_value);
        //cout<<"----------------------增加term"<<t<<", 需要词频="<<supposed_frequency<<"----------------------"<<endl;
        int add_num = supposed_frequency;
        int added_count = 0;
        while(added_count<add_num){
            int random_id = (int)random() % poi_size;  //任意取一个poi
            int _size = POIs[random_id].keywordSet.size();
            if(_size > 20) continue;
            if(poiTerm_Mark[random_id].count(target_term)){ //若这poi已经有该keyword了， 则discard
                continue;
            }
            else{
                POIs[random_id].keywordSet.insert(target_term);  //否则向poi中加入该keyword
                poiTerm_Mark[random_id][target_term] = true;
                termPOI_Mark[target_term][random_id] = true;
                added_count++;
            }
        }


    }

    map<int, set<int>> inverted_list;
    for(POI p: POIs){
        for(int term: p.keywordSet){
            inverted_list[term].insert(p.id);
        }
    }

    for(int rank=1; rank<=vocabularySize;rank++){
        int term = rank;
        int frequency = inverted_list[term].size();
        ///cout<<"************************rank"<<rank<<"************************"<<endl;
        if(rank%1000==0){
            cout<<"term="<<term<<", 当前frequency="<<frequency<<",估计frequency="<<computingTerFreZipfan(rank,C_value)<<endl;
        }
        //getchar();
    }


    cout<<"End Zipfan_POITextInfo!"<<endl;
    cout<<"Output new POIInfor...";
    ofstream foutPOI;
    string newPOIPath = getObjectNewInputPath(BUSINESS);
    foutPOI.open(newPOIPath);

    for(int i=0;i<POIs.size();i++){
        POI p = POIs[i];
        foutPOI << p.id << " " << p.Ni << " " << p.Nj << " " << p.dist << " " << p.dis << " ";
        for(int term: p.keywordSet){
            foutPOI << term << " ";
        }
        foutPOI << endl;
    }

    foutPOI.close();
    cout<<"DONE!"<<endl;




}

void Zipfan_UserTextInfo(){ // leaf node level,去掉了之前多余累赘的 address_v2p对应的索引与数据文件
    int shi = 0;
    cout<<"Begin Zipfan_UserTextInfo!"<<endl;
    LOAD_START
    // load PoiMap
    printf("LOADING User INFO...");
    ifstream finUsr;
    string poiPath = getObjectInputPath(USER);
    string str;
    finUsr.open(poiPath.c_str());


    map<int,map<int,bool>> termUsr_Mark;
    map<int,map<int,bool>> usrTerm_Mark;
    map<int,set<int>> posting_list;
    set<int> current_words;
    while (getline(finUsr, str)) {
        istringstream tt(str);
        int usr_id, usr_Ni, usr_Nj;
        float usr_dist, usr_dis;
        tt >> usr_id >> usr_Ni >> usr_Nj >> usr_dist >> usr_dis;
        float exchangeNum = (int) (usr_dist * WEIGHT_INFLATE_FACTOR);
        usr_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (usr_dis * WEIGHT_INFLATE_FACTOR);
        usr_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;


        vector<int> uKey;
        set<int> keywordSet;
        int keyTmp;
        //cout<<"keyword:";
        int term = -1;
        while (tt >> term) {
            if (0 > term) {
                break;
            }
            if(term==0) continue;  //关键词id从1开始
            uKey.push_back(term); keywordSet.insert(term);
            current_words.insert(term);
            posting_list[term].insert(usr_id);
            //cout<<",w"<<term;
            usrTerm_Mark[usr_id][term] = true;
            termUsr_Mark[term][usr_id] = true;

        }
        // cout<<endl;

        User tmpUsr;
        tmpUsr.id = usr_id;
        tmpUsr.Ni = min(usr_Ni,usr_Nj);
        tmpUsr.Nj = max(usr_Ni,usr_Nj);
        tmpUsr.dist = usr_dist;
        tmpUsr.dis = usr_dis;
        tmpUsr.keywords = uKey;
        tmpUsr.keywordSet = keywordSet;
        Users.push_back(tmpUsr);

        usrCnt++;
    }
    finUsr.close();
    cout << " DONE! USER #: " << poiCnt << endl;
    LOAD_END
    LOAD_PRINT(" load User ")



    int term_size = current_words.size();
    cout<<"current term_size="<<term_size<<endl;


    int usr_size = Users.size();
    cout<<"user_size="<<usr_size<<endl;
    ///getchar();

    float C_ratio = 0.13;

    int C_value = C_ratio * usr_size;

    int user_vocabulary = (int)(vocabularySize*0.8);

    ////先对已有的keyword进行处理：
    for(int i=1;i<=term_size;i++){

        //if(i>15) break;
        int target_term = i;
        //float term_weight = getTermIDFWeight(target_term);
        int current_frequency = posting_list[target_term].size();
        cout<<"----------------------term"<<target_term<<"----------------------"<<endl;

        int rank = i;
        int supposed_frequency = computingTerFreZipfan(rank,C_value);
        cout<<"supposed_frequency ="<<supposed_frequency;
        cout<<",current_frequency ="<<current_frequency<<endl;
        if(supposed_frequency>current_frequency ){  //需要加入

            int add_num = supposed_frequency - current_frequency;
            int added_count = 0;
            while(added_count<add_num){
                int random_id = (int)random() % usr_size;  //任意取一个usr
                int _size = Users[random_id].keywordSet.size();
                if(_size > 8) continue;
                if(usrTerm_Mark[random_id].count(target_term)){ //若这user已经有该keyword了， 则discard
                    continue;
                }
                else{
                    Users[random_id].keywordSet.insert(target_term);  //否则向user中加入该keyword
                    usrTerm_Mark[random_id][target_term] = true;
                    termUsr_Mark[target_term][random_id] = true;
                    added_count++;
                }
            }

        }
        else if(supposed_frequency<current_frequency){  //需要删除
            set<int> _list = posting_list[target_term];
            vector<int> users_ids;
            for(int id: _list)
                users_ids.push_back(id);

            int reduced_num = current_frequency-supposed_frequency;
            int reduced_count = 0;

            for(int i=0;i<reduced_num;i++){
                int target_id = users_ids[i];
                int _size2 = Users[target_id].keywordSet.size();
                if(_size2<2) {
                    Users[target_id].keywordSet.erase(target_term) ;
                    usrTerm_Mark[target_id][target_term] = false;
                    termUsr_Mark[target_term][target_id]= false;
                    int  _term = (int) random() % user_vocabulary;
                    int _coin = (int)random() % 4;
                    if(_coin!=1){  //0.75的概率是前500
                        int _th = random()% 500;
                        int _coin2 = random()%2;
                        if(_coin2==0){  //1/2的概率前100
                            _term = random()% 100;
                            int  _coin3 = (int)random()%2;
                            if(_coin3==0){ //1/2概率前50
                                _term = (int)random()%20;
                            }
                        }
                        else{
                            _term = _th;
                        }

                    }

                    _term = max(1, _term);

                    Users[target_id].keywordSet.insert(_term);
                    usrTerm_Mark[target_id][_term] = true;
                    termUsr_Mark[_term][target_id]= true;
                }
                else{
                    Users[target_id].keywordSet.erase(target_term) ;
                    usrTerm_Mark[target_id][target_term] = false;
                    termUsr_Mark[target_term][target_id]= false;
                }
            }



        }



    }



    for(int t = term_size+1;t<=user_vocabulary;t++){
        int rank = t;
        int target_term = t;
        int supposed_frequency = max(1,computingTerFreZipfan(rank,C_value));
        cout<<"----------------------增加term"<<t<<", 需要词频="<<supposed_frequency<<"----------------------"<<endl;
        int add_num = supposed_frequency;
        int added_count = 0;
        while(added_count<add_num){
            int random_id = (int)random() % usr_size;  //任意取一个user
            int _size = Users[random_id].keywordSet.size();
            if(_size > 20) continue;
            if(usrTerm_Mark[random_id].count(target_term)){ //若这user已经有该keyword了， 则discard
                continue;
            }
            else{
                Users[random_id].keywordSet.insert(target_term);  //否则向poi中加入该keyword
                usrTerm_Mark[random_id][target_term] = true;
                termUsr_Mark[target_term][random_id] = true;
                added_count++;
            }
        }


    }

    map<int, set<int>> inverted_list;
    for(User u: Users){
        for(int term: u.keywordSet){
            inverted_list[term].insert(u.id);
        }
    }

    for(int rank=1; rank<=user_vocabulary;rank++){
        int term = rank;
        int frequency = inverted_list[term].size();
        cout<<"************************rank"<<rank<<"************************"<<endl;
        cout<<"term="<<term<<", 当前frequency="<<frequency<<",估计frequency="<<computingTerFreZipfan(rank,C_value)<<endl;
        //getchar();
    }


    cout<<"End Zipfan_UserTextInfo!"<<endl;
    cout<<"Output new UserInfor...";
    ofstream foutUsr;
    string newUserPath = getObjectNewInputPath(USER);
    foutUsr.open(newUserPath);

    for(int i=0;i<Users.size();i++){
        User u = Users[i];
        foutUsr << u.id << " " << u.Ni << " " << u.Nj << " " << u.dist << " " << u.dis << " ";
        for(int term: u.keywordSet){
            foutUsr << term << " ";
        }
        foutUsr << endl;
    }

    foutUsr.close();
    cout<<"DONE!"<<endl;




}

void Varying_POIInfo_old(){ // leaf node level,去掉了之前多余累赘的 address_v2p对应的索引与数据文件
    int shi = 0;
    cout<<"Begin Adding_POIInfo!"<<endl;
    LOAD_START
    // load PoiMap
    printf("LOADING POI INFO...");
    ifstream finPoi;
    string poiPath = getObjectNewInputPath(BUSINESS); //getObjectInputPath(BUSINESS);
    string str;
    finPoi.open(poiPath.c_str());


    map<int,map<int,bool>> termPOI_Mark;
    map<int,map<int,bool>> poiTerm_Mark;
    map<int,set<int>> posting_list;
    set<int> current_words; POIs.clear();
    while (getline(finPoi, str)) {
        istringstream tt(str);
        int POI_id, POI_Ni, POI_Nj;
        float POI_dist, POI_dis;
        tt >> POI_id >> POI_Ni >> POI_Nj >> POI_dist >> POI_dis;
        float exchangeNum = (int) (POI_dist * WEIGHT_INFLATE_FACTOR);
        POI_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (POI_dis * WEIGHT_INFLATE_FACTOR);
        POI_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;

        vector<int> uKey;
        set<int> keywordSet;
        int keyTmp;
        //cout<<"keyword:";
        int term = -1;
        while (tt >> term) {
            if (0 > term) {
                break;
            }
            uKey.push_back(term); keywordSet.insert(term);
            current_words.insert(term);
            posting_list[term].insert(POI_id);
            //cout<<",w"<<term;
            poiTerm_Mark[POI_id][term] = true;
            termPOI_Mark[term][POI_id] = true;

        }
        // cout<<endl;

        POI tmpPoi;
        tmpPoi.id = POI_id;
        tmpPoi.Ni = min(POI_Ni,POI_Nj);
        tmpPoi.Nj = max(POI_Ni,POI_Nj);
        tmpPoi.dist = POI_dist;
        tmpPoi.dis = POI_dis;
        tmpPoi.category = 1;
        tmpPoi.keywords = uKey;
        tmpPoi.keywordSet = keywordSet;
        POIs.push_back(tmpPoi);

        poiCnt++;
    }
    finPoi.close();
    cout << " DONE! POI #: " << poiCnt << endl;
    LOAD_END
    LOAD_PRINT(" load POI ")
    cout <<"Add new POI info based on different adding rate!"<<endl;
    float rate_0 = 1.0; float delta_rate = 0.2;
    float rate = 1.0;
    for(int i=1;i<=5;i++){
        rate =  i*(delta_rate);
        map<int, set<int>> posting_list;
        map<int,map<int,bool>> termPOI_Mark;
        map<int,map<int,bool>> poiTerm_Mark;
        int _size = rate * POIs.size();
        vector<POI> POIs_extract;
        while(POIs_extract.size()<_size){
            POI new_poi;
            new_poi.id = POIs_extract.size() ;
            int random_id = (int)random() % POIs.size();  //任意取一个poi，提取空间信息
            new_poi.Ni=  POIs[random_id].Ni;
            new_poi.Nj=  POIs[random_id].Nj;
            new_poi.dist = POIs[random_id].dist;
            new_poi.dis = POIs[random_id].dis;
            int random_id2 = (int)random() % POIs.size();  //再任意取一个poi， 提取文本信息
            new_poi.keywordSet = POIs[random_id2].keywordSet;
            POIs_extract.push_back(new_poi);
            for(int term: new_poi.keywordSet){
                posting_list[term].insert(new_poi.id);
            }
        }
        //对关键词rank进行重新排序，更新term_id
        map<int,set<int>>:: iterator iter;
        priority_queue<CardinalityEntry> Q;
        for(iter=posting_list.begin();iter!=posting_list.end();iter++){
            int term_id = iter->first;
            int term_posting_size = iter->second.size();
            CardinalityEntry entry(term_id,term_posting_size);
            Q.push(entry);
        }
        int rank=1; map<int,int> termNewIDMap;
        while(!Q.empty()){
            CardinalityEntry entry = Q.top();
            Q.pop();
            int current_term = entry.id;
            int current_size = entry.size;
            termNewIDMap[current_term]=rank;
            rank++;
        }


        cout<<"End extracting part of POIs Info!"<<endl;
        cout<<"Generate "<<POIs_extract.size()<<" new POIs!"<<endl;
        cout<<"Output extracted POIInfor...";
        ofstream foutPOI;
        string updatedPOIPath = getObjectVaryingInputPath(BUSINESS,rate);
        foutPOI.open(updatedPOIPath);




        for(int j=0;j<POIs_extract.size();j++){
            POI p = POIs_extract[j];
            foutPOI << p.id << " " << p.Ni << " " << p.Nj << " " << p.dist << " " << p.dis << " ";
            for(int term: p.keywordSet){
                int term_new = termNewIDMap[term];
                foutPOI << term_new << " ";
            }
            foutPOI << endl;
        }

        foutPOI.close();
        cout<<"POI DONE!"<<endl;


        //-----------------------------更新user中的关键词情况
        printf("LOADING User INFO..."); Users.clear();
        ifstream finUsr;
        string userPath = getObjectInputPath(USER);
        string str;
        finUsr.open(userPath.c_str());

        set<int> current_words;
        while (getline(finUsr, str)) {
            istringstream tt(str);
            int usr_id, usr_Ni, usr_Nj;
            float usr_dist, usr_dis;
            tt >> usr_id >> usr_Ni >> usr_Nj >> usr_dist >> usr_dis;
            float exchangeNum = (int) (usr_dist * WEIGHT_INFLATE_FACTOR);
            usr_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
            exchangeNum = (int) (usr_dis * WEIGHT_INFLATE_FACTOR);
            usr_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;


            vector<int> uKey;
            set<int> keywordSet;
            int keyTmp;
            //cout<<"keyword:";
            int term = -1;
            while (tt >> term) {
                if (0 > term) {
                    break;
                }
                uKey.push_back(term); keywordSet.insert(term);
                current_words.insert(term);
            }
            // cout<<endl;

            User tmpUsr;
            tmpUsr.id = usr_id;
            tmpUsr.Ni = min(usr_Ni,usr_Nj);
            tmpUsr.Nj = max(usr_Ni,usr_Nj);
            tmpUsr.dist = usr_dist;
            tmpUsr.dis = usr_dis;
            tmpUsr.keywords = uKey;
            tmpUsr.keywordSet = keywordSet;
            Users.push_back(tmpUsr);

            usrCnt++;
        }
        finUsr.close();
        LOAD_END
        LOAD_PRINT(" load User ")
        string updatedUserPath = getUserTermUpdatingInputPath(rate);
        cout<<"Updating term id in User Info...";
        ofstream foutUsr;
        foutUsr.open(updatedUserPath);

        for(int j=0;j<Users.size();j++){
            User u = Users[j];
            foutUsr << u.id << " " << u.Ni << " " << u.Nj << " " << u.dist << " " << u.dis << " ";
            for(int term: u.keywordSet){
                if(termNewIDMap.count(term)==true){
                    int termNew = termNewIDMap[term];
                    foutUsr << termNew << " ";
                }

            }
            foutUsr << endl;
        }

        foutUsr.close();
        cout<<"DONE!"<<endl;




    }



}


void Varying_POIInfo(){ // leaf node level,去掉了之前多余累赘的 address_v2p对应的索引与数据文件
    int shi = 0;
    cout<<"Begin Adding_POIInfo!"<<endl;
    LOAD_START
    // load PoiMap
    printf("LOADING POI INFO...");
    ifstream finPoi;
    string poiPath = getObjectNewInputPath(BUSINESS); //getObjectInputPath(BUSINESS);
    string str;
    finPoi.open(poiPath.c_str());


    map<int,map<int,bool>> termPOI_Mark;
    map<int,map<int,bool>> poiTerm_Mark;
    map<int,set<int>> posting_list;
    set<int> current_words; POIs.clear();
    while (getline(finPoi, str)) {
        istringstream tt(str);
        int POI_id, POI_Ni, POI_Nj;
        float POI_dist, POI_dis;
        tt >> POI_id >> POI_Ni >> POI_Nj >> POI_dist >> POI_dis;
        float exchangeNum = (int) (POI_dist * WEIGHT_INFLATE_FACTOR);
        POI_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (POI_dis * WEIGHT_INFLATE_FACTOR);
        POI_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;

        vector<int> uKey;
        set<int> keywordSet;
        int keyTmp;
        //cout<<"keyword:";
        int term = -1;
        while (tt >> term) {
            if (0 > term) {
                break;
            }
            uKey.push_back(term); keywordSet.insert(term);
            current_words.insert(term);
            posting_list[term].insert(POI_id);
            //cout<<",w"<<term;
            poiTerm_Mark[POI_id][term] = true;
            termPOI_Mark[term][POI_id] = true;

        }
        // cout<<endl;

        POI tmpPoi;
        tmpPoi.id = POI_id;
        tmpPoi.Ni = min(POI_Ni,POI_Nj);
        tmpPoi.Nj = max(POI_Ni,POI_Nj);
        tmpPoi.dist = POI_dist;
        tmpPoi.dis = POI_dis;
        tmpPoi.category = 1;
        tmpPoi.keywords = uKey;
        tmpPoi.keywordSet = keywordSet;
        POIs.push_back(tmpPoi);

        poiCnt++;
    }
    finPoi.close();
    cout << " DONE! POI #: " << poiCnt << endl;
    LOAD_END
    LOAD_PRINT(" load POI ")
    cout <<"Add new POI info based on different adding rate!"<<endl;
    float rate_0 = 1.0; float delta_rate = 0.1;//0.2;
    float rate = 1.0;
    for(int i=1;i<=5;i++){
        //rate =  i*(delta_rate);
        rate =  0.6+ (i-1)*delta_rate;
        map<int, set<int>> posting_list;
        map<int,map<int,bool>> termPOI_Mark;
        map<int,map<int,bool>> poiTerm_Mark;
        int _size = rate * POIs.size();
        vector<POI> POIs_extract;
        if(i==5){
            cout<<"no need to extract POI!"<<endl;
            cout<<"Output directly POIIn ...";
            ofstream foutPOI;
            string updatedPOIPath = getObjectVaryingInputPath(BUSINESS,rate);
            foutPOI.open(updatedPOIPath);

            for(int j=0;j<POIs.size();j++){
                POI p = POIs[j];
                foutPOI << p.id << " " << p.Ni << " " << p.Nj << " " << p.dist << " " << p.dis << " ";
                for(int term: p.keywordSet){
                    foutPOI << term << " ";
                }
                foutPOI << endl;
            }

            foutPOI.close();
            cout<<"POI DONE!"<<endl;


        }
        else{
            map<int, bool> poi_flag;
            while(POIs_extract.size()<_size){
                POI new_poi;
                new_poi.id = POIs_extract.size() ;
                int random_id = (int)random() % POIs.size();  //任意取一个poi，提取空间信息
                if(poi_flag.count(random_id)==true) continue;
                else{
                    new_poi.Ni=  POIs[random_id].Ni;
                    new_poi.Nj=  POIs[random_id].Nj;
                    new_poi.dist = POIs[random_id].dist;
                    new_poi.dis = POIs[random_id].dis;
                    new_poi.keywordSet = POIs[random_id].keywordSet;
                    POIs_extract.push_back(new_poi);
                    for(int term: new_poi.keywordSet){
                        posting_list[term].insert(new_poi.id);
                    }
                    poi_flag[random_id]= true;
                }

            }
            //对关键词rank进行重新排序，更新term_id
            map<int,set<int>>:: iterator iter;
            priority_queue<CardinalityEntry> Q;
            for(iter=posting_list.begin();iter!=posting_list.end();iter++){
                int term_id = iter->first;
                int term_posting_size = iter->second.size();
                CardinalityEntry entry(term_id,term_posting_size);
                Q.push(entry);
            }
            int rank=1; map<int,int> termNewIDMap;
            while(!Q.empty()){
                CardinalityEntry entry = Q.top();
                Q.pop();
                int current_term = entry.id;
                int current_size = entry.size;
                termNewIDMap[current_term]=rank;
                rank++;
            }


            cout<<"End extracting part of POIs Info!"<<endl;
            cout<<"extract "<<POIs_extract.size()<<" new POIs!"<<endl;
            cout<<"Output extracted POIInfor...";
            ofstream foutPOI;
            string updatedPOIPath = getObjectVaryingInputPath(BUSINESS,rate);
            foutPOI.open(updatedPOIPath);

            for(int j=0;j<POIs_extract.size();j++){
                POI p = POIs_extract[j];
                foutPOI << p.id << " " << p.Ni << " " << p.Nj << " " << p.dist << " " << p.dis << " ";
                for(int term: p.keywordSet){
                    int term_new = termNewIDMap[term];
                    foutPOI << term_new << " ";
                }
                foutPOI << endl;
            }

            foutPOI.close();
            cout<<"POI DONE!"<<endl;

            //-----------------------------更新user中的关键词情况
            printf("LOADING User INFO..."); Users.clear();
            ifstream finUsr;
            string userPath = getObjectNewInputPath(USER);
            string str;
            finUsr.open(userPath.c_str());

            set<int> current_words;
            while (getline(finUsr, str)) {
                istringstream tt(str);
                int usr_id, usr_Ni, usr_Nj;
                float usr_dist, usr_dis;
                tt >> usr_id >> usr_Ni >> usr_Nj >> usr_dist >> usr_dis;
                float exchangeNum = (int) (usr_dist * WEIGHT_INFLATE_FACTOR);
                usr_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
                exchangeNum = (int) (usr_dis * WEIGHT_INFLATE_FACTOR);
                usr_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;


                vector<int> uKey;
                set<int> keywordSet;
                int keyTmp;
                //cout<<"keyword:";
                int term = -1;
                while (tt >> term) {
                    if (0 > term) {
                        break;
                    }
                    uKey.push_back(term); keywordSet.insert(term);
                    current_words.insert(term);
                }
                // cout<<endl;

                User tmpUsr;
                tmpUsr.id = usr_id;
                tmpUsr.Ni = min(usr_Ni,usr_Nj);
                tmpUsr.Nj = max(usr_Ni,usr_Nj);
                tmpUsr.dist = usr_dist;
                tmpUsr.dis = usr_dis;
                tmpUsr.keywords = uKey;
                tmpUsr.keywordSet = keywordSet;
                Users.push_back(tmpUsr);

                usrCnt++;
            }
            finUsr.close();
            LOAD_END
            LOAD_PRINT(" load User ")
            string updatedUserPath = getUserTermUpdatingInputPath(rate);
            cout<<"Updating term id in User Info...";
            ofstream foutUsr;
            foutUsr.open(updatedUserPath);

            for(int j=0;j<Users.size();j++){
                User u = Users[j];
                foutUsr << u.id << " " << u.Ni << " " << u.Nj << " " << u.dist << " " << u.dis << " ";
                for(int term: u.keywordSet){
                    if(termNewIDMap.count(term)==true){
                        int termNew = termNewIDMap[term];
                        foutUsr << termNew << " ";
                    }

                }
                foutUsr << endl;
            }

            foutUsr.close();
            cout<<"DONE!"<<endl;

        }





    }



}




void Varying_UserInfo(){ // leaf node level,去掉了之前多余累赘的 address_v2p对应的索引与数据文件
    int shi = 0;
    cout<<"Begin Varying User Info!"<<endl;
    LOAD_START
    // load PoiMap
    printf("LOADING User INFO...");
    ifstream finUsr;
    string poiPath = getObjectNewInputPath(USER);
    string str;
    finUsr.open(poiPath.c_str());


    map<int,map<int,bool>> termUsr_Mark;
    map<int,map<int,bool>> usrTerm_Mark;
    map<int,set<int>> posting_list;
    set<int> current_words;
    Users.clear();
    while (getline(finUsr, str)) {
        istringstream tt(str);
        int usr_id, usr_Ni, usr_Nj;
        float usr_dist, usr_dis;
        tt >> usr_id >> usr_Ni >> usr_Nj >> usr_dist >> usr_dis;
        float exchangeNum = (int) (usr_dist * WEIGHT_INFLATE_FACTOR);
        usr_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (usr_dis * WEIGHT_INFLATE_FACTOR);
        usr_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;


        vector<int> uKey;
        set<int> keywordSet;
        int keyTmp;
        //cout<<"keyword:";
        int term = -1;
        while (tt >> term) {
            if (0 > term) {
                break;
            }
            uKey.push_back(term); keywordSet.insert(term);
            current_words.insert(term);
            posting_list[term].insert(usr_id);
            //cout<<",w"<<term;
            usrTerm_Mark[usr_id][term] = true;
            termUsr_Mark[term][usr_id] = true;

        }
        // cout<<endl;

        User tmpUsr;
        tmpUsr.id = usr_id;
        tmpUsr.Ni = min(usr_Ni,usr_Nj);
        tmpUsr.Nj = max(usr_Ni,usr_Nj);
        tmpUsr.dist = usr_dist;
        tmpUsr.dis = usr_dis;
        tmpUsr.keywords = uKey;
        tmpUsr.keywordSet = keywordSet;
        Users.push_back(tmpUsr);

        usrCnt++;
    }
    finUsr.close();
    cout << " DONE! USER #: " << poiCnt << endl;
    LOAD_END
    LOAD_PRINT(" load User ")
    cout<<"End adding new User Info!"<<endl;

    float rate_0 = 1.0; float delta_rate = 0.1;
    float rate = 1.0;
    for(int i=1;i<=5;i++){
        rate =  0.6+(i-1)*delta_rate;
        int _size = rate *Users.size();
        vector<User> Users_varying;

        if(i==5){
            cout<<"Directly output UserInfor...";
            ofstream foutUser;
            string updatedUserPath = getObjectVaryingInputPath(USER,rate);
            ofstream foutUsr;
            foutUsr.open(updatedUserPath);

            for(int j=0;j<Users.size();j++){
                User u = Users[j];
                foutUsr << u.id << " " << u.Ni << " " << u.Nj << " " << u.dist << " " << u.dis << " ";
                for(int term: u.keywordSet){
                    foutUsr << term << " ";
                }
                foutUsr << endl;
            }
            foutUsr.close();
            cout<<"DONE!"<<endl;
        }
        else{
            map<int, bool> user_flag;
            while(Users_varying.size()<=_size){
                User new_user;
                new_user.id = Users_varying.size();
                int random_id = (int)random() % Users.size();  //任意取一个user，提取空间信息
                if(user_flag.count(random_id)==true) continue;

                new_user.Ni=  Users[random_id].Ni;
                new_user.Nj=  Users[random_id].Nj;
                new_user.dist = Users[random_id].dist;
                new_user.dis = Users[random_id].dis;
                int random_id2 = random_id;//(保留他的文本信息) //(int)random() % Users.size();  //再任意取一个user， 提取文本信息
                new_user.keywordSet = Users[random_id2].keywordSet;
                Users_varying.push_back(new_user);
                user_flag[random_id] = true;
            }

            cout<<"End extract partial User Info!"<<endl;
            cout<<"extract "<<Users_varying.size()<<" Users!"<<endl;
            cout<<"Output updated UserInfor...";
            ofstream foutUser;
            string updatedUserPath = getObjectVaryingInputPath(USER,rate);
            cout<<"Output new UserInfor...";
            ofstream foutUsr;
            foutUsr.open(updatedUserPath);


            for(int j=0;j<Users_varying.size();j++){
                User u = Users_varying[j];
                foutUsr << u.id << " " << u.Ni << " " << u.Nj << " " << u.dist << " " << u.dis << " ";
                for(int term: u.keywordSet){
                    foutUsr << term << " ";
                }
                foutUsr << endl;
            }

            foutUsr.close();
            cout<<"DONE!"<<endl;
        }


    }


}

void Varying_buildFold(){
    /*----------------IGtree 对应 varying user 以及 poi----------------------------------------*/
    ////新建 Varying_IGtree文件夹
    string folderPath ;
    stringstream ss;
    ss<<DATASET_HOME<<dataset_name<<"/Varying_IGtree/";
    folderPath = ss.str();
    cout<<folderPath<<endl;
    string command;
    command = string("mkdir -p ") + folderPath;
    system(command.c_str());
    ////在Varying_IGtree(以及NVD)文件夹下建立 VaryingO , VaryingU
    ss.str("");
    ss<<DATASET_HOME<<dataset_name<<"/Varying_IGtree/VaryingO";
    folderPath = ss.str();
    string o_basepath = folderPath;
    cout<<o_basepath<<endl;
    command = string("mkdir -p ") + folderPath;
    system(command.c_str());
    ss.str("");
    ss<<DATASET_HOME<<dataset_name<<"/NVD/VaryingO";
    folderPath = ss.str();
    string nvdo_basepath = folderPath;
    cout<<nvdo_basepath<<endl;
    command = string("mkdir -p ") + nvdo_basepath;
    system(command.c_str());


    ss.str("");
    ss<<DATASET_HOME<<dataset_name<<"/Varying_IGtree/VaryingU";
    folderPath = ss.str();
    string u_basepath = folderPath;
    cout<<u_basepath<<endl;
    command = string("mkdir -p ") + folderPath;
    system(command.c_str());
    ////在Varying_IGtree/VaryingO (VaryingU) 下建立各个百分比对应的数据集子文件夹
    for(int i=1;i<=5;i++){
        float rate =  0.6+(i-1)*0.1;
        char command_str_u[200];
        sprintf(command_str_u,"mkdir -p %s/%d%/",u_basepath.c_str(),(int)(rate*100));
        system(command_str_u);
        char command_str_o[200];
        sprintf(command_str_o,"mkdir -p %s/%d%/",o_basepath.c_str(),(int)(rate*100));
        system(command_str_o);
        char command_str_nvdo[200];
        sprintf(command_str_nvdo,"mkdir -p %s/%d%/",nvdo_basepath.c_str(),(int)(rate*100));
        system(command_str_nvdo);
    }

}


int ModifyFriendShip(){  // 从source_dataset中读取 源数据内容， 处理后输入到social 文件夹下
    printf("ModifyFriendShip...");

    string linkPath = getSocialRawInputPath(LINK);
    string linkNewPath = getSocialInputPath(LINK);
    FILE *fin;
    fin = fopen(linkPath.c_str(), "r");
    FILE *fout;
    fout = fopen(linkNewPath.c_str(), "w");


    string str;
    //专门为IM问题建立的社交网络图
    // Graph vecGRev(graph_vertex_num);  //邻接边表, 注意图顶点的个数
    //std::vector<size_t> vecInDeg(graph_vertex_num);  //顶点入度表
    //读取图数据
    wholeSocialLinkMap.clear();

    int ui,uj;
    while (fscanf(fin,"%d %d ", &ui,&uj)==2) {


#ifdef LasVegas
        if(ui<social_UserID_MaxKey && uj<social_UserID_MaxKey){
            fprintf(fout,"%d %d\n",ui,uj);
        }
#endif

#ifdef Twitter
        if(ui<social_UserID_MaxKey && uj<social_UserID_MaxKey){
            fprintf(fout,"%d %d\n",ui,uj);
        }
#endif

        //为IM的graph加载数据
        float weight = 0.0;
        //vecGRev[uj].push_back(Edge(ui, weight));
    }

    fclose(fin);
    fclose(fout);



    printf("COMPLETE!\n");
    return  0;
}

typedef struct {
    int u;
    int v;
} SocialLink;

typedef struct{
    int p;
    int u;
    int count;
}CheckInRecord;

void loadLasVegasUserInfo(){
    printf("LOADING User INFO...");
    ifstream finUsr; string str; Users.clear();
    string userPath = getObjectNewInputPath(USER);
    finUsr.open(userPath.c_str());
    while (getline(finUsr, str)) {
        istringstream tt(str);
        int usr_id, usr_Ni, usr_Nj;
        float usr_dist, usr_dis;
        tt >> usr_id >> usr_Ni >> usr_Nj >> usr_dist >> usr_dis;

        float exchangeNum = (int) (usr_dist * WEIGHT_INFLATE_FACTOR);
        usr_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (usr_dis * WEIGHT_INFLATE_FACTOR);
        usr_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;

        vector<int> uKey;
        set<int> keywordSet;
        int keyTmp;
        //cout<<"keyword:";
        int term = -1;
        while (tt >> term) {
            if (0 > term) {
                break;
            }
            keywordSet.insert(term);
            uKey.push_back(term);

        }
        User tmpUsr;
        tmpUsr.id = usr_id;
        tmpUsr.Ni = min(usr_Ni,usr_Nj);
        tmpUsr.Nj = max(usr_Ni,usr_Nj);
        tmpUsr.dist = usr_dist;
        tmpUsr.dis = usr_dis;
        tmpUsr.keywords = uKey;
        tmpUsr.keywordSet = keywordSet;
        Users.push_back(tmpUsr);
        usrCnt++;

    }
    finUsr.close();
    cout << " DONE! USER #: " << usrCnt << endl;
}

int ModifyLasVegasUserIDValue(){
    online_users.clear(); Users.clear();
    set<int> review_userSet; review_userSet.clear();
    set<int> online_userSet; online_userSet.clear();
    ////1.加载.user中的用户信息
    printf("LOADING User INFO...");
    ifstream finUsr; string str;
    string userPath = getObjectNewInputPath(USER);
    finUsr.open(userPath.c_str());
    unordered_map<int,bool> idMaskMap;
    while (getline(finUsr, str)) {
        istringstream tt(str);
        int usr_id, usr_Ni, usr_Nj;
        float usr_dist, usr_dis;
        tt >> usr_id >> usr_Ni >> usr_Nj >> usr_dist >> usr_dis;
        ////将有评论的用户加入review_userSet
        review_userSet.insert(usr_id);
        idMaskMap[usr_id]=true;

        float exchangeNum = (int) (usr_dist * WEIGHT_INFLATE_FACTOR);
        usr_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (usr_dis * WEIGHT_INFLATE_FACTOR);
        usr_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;

        vector<int> uKey;
        set<int> keywordSet;
        int keyTmp;
        //cout<<"keyword:";
        int term = -1;
        while (tt >> term) {
            if (0 > term) {
                break;
            }
            keywordSet.insert(term);
            uKey.push_back(term);

        }
        User tmpUsr;
        tmpUsr.id = usr_id;
        tmpUsr.Ni = min(usr_Ni,usr_Nj);
        tmpUsr.Nj = max(usr_Ni,usr_Nj);
        tmpUsr.dist = usr_dist;
        tmpUsr.dis = usr_dis;
        tmpUsr.keywords = uKey;
        tmpUsr.keywordSet = keywordSet;
        Users.push_back(tmpUsr);

        usrCnt++;

    }
    finUsr.close();
    cout << " DONE! USER #: " << usrCnt << endl;


    ////2.加载FriendShip中的用户ID信息
    printf("loadFriendShip...");
    string linkPath = getSocialInputPath(LINK);
    FILE *fin;
    fin = fopen(linkPath.c_str(), "r");

    int ui=-1; int uj=-1; vector<SocialLink> slinkCollection;
    unordered_map<int,set<int>> user_neigbors;
    while (fscanf(fin,"%d %d ", &ui,&uj)==2) {

        int u = ui; int v= uj;
        ////若该用户没有评论，则加入online_userSet
        if(idMaskMap.count(u)==false)
            online_userSet.insert(u);
        if(idMaskMap.count(v)==false)
            online_userSet.insert(v);
        SocialLink link;
        link.u =ui;link.v=uj;
        slinkCollection.push_back(link);
        user_neigbors[ui].insert(uj);
        //为IM的graph加载数据
        float weight = 0.0;
        //vecGRev[uj].push_back(Edge(ui, weight));
    }

    fclose(fin);
    printf("COMPLETE!\n");

    ////3.加载Check-in中的用户ID信息
    int poi,u,count; vector<CheckInRecord> ckCollection;
    printf("LOADING CHECKIN...");
#ifdef LasVegas
    string checkinPath = getSocialRawInputPath(CHECKIN);
#else
    #ifdef Twitter
        string checkinPath = getSocialRawInputPath(CHECKIN);
    #else
        string checkinPath = getSocialInputPath(CHECKIN);
    #endif
#endif
    fin = fopen(checkinPath.c_str(),"r");
    while (fscanf(fin,"%d %d %d", &poi,&u,&count)==3) {
        //istringstream tt(str);
        //tt>>poi>>u>>count;
        if(poi< (poi_num-1) && u<(UserID_MaxKey-1)){
            if(user_neigbors.count(u)){ //只留下有社交关系的用户（筛除掉孤立点）
                if(idMaskMap.count(u)==false)
                    online_userSet.insert(u);
                CheckInRecord checkINR;
                checkINR.p = poi; checkINR.u=u; checkINR.count = count;
                ckCollection.push_back(checkINR);
            }


        }
    }
    fclose(fin); printf("COMPLETE!\t");


    ////对所有review_userSet中所有user的id与idx建立映射关系
    unordered_map<int,int> online_user_id2IdxMap;
    //unordered_map<int,int> online_user_idx2IdMap;
    int idx = 0;
    for(int u_id: review_userSet){
        online_user_id2IdxMap[u_id] = idx;
        //online_user_idx2IdMap[idx] = u_id;
        idx++;
    }
    idx = 0;
    for(int u_id: online_userSet){
        online_user_id2IdxMap[u_id] = review_userSet.size()+idx;
        idx++;
    }
    cout<<"总共有"<<(review_userSet.size()+online_userSet.size())<<"个不同用户的ID"<<endl;


    ////4. 输出更新用户ID后的相关文件
    cout<<"输出更新用户ID后的相关文件..."<<endl;
    //1. .user文件更新
    ofstream foutUser;
    string newUserPath = getObjectNewInputPath(USER);
    foutUser.open(newUserPath);
    for(int i=0;i<Users.size();i++){
        User u = Users[i];
        int old_ID = u.id;
        int newID = online_user_id2IdxMap[old_ID];
        foutUser << newID << " " << u.Ni << " " << u.Nj << " " << u.dist << " " << u.dis << " ";
        int count = 0;
        for(int term: u.keywordSet){
            foutUser << term << " ";

        }
        foutUser << endl;
    }
    foutUser.close();
    cout<<"Update .user 文件！"<<endl;

    //2. social link 文件更新

    unordered_map<int,set<int>> wholeSocialLink;
    for(int i=0;i<slinkCollection.size();i++){
        SocialLink link = slinkCollection[i];
        int u_newID = online_user_id2IdxMap[link.u];
        int v_newID = online_user_id2IdxMap[link.v];
        cout<<"原来是：u"<<link.u<<", u"<<link.v<<endl;
        cout<<"现在是：u"<<u_newID<<", u"<<v_newID<<endl;
        //fprintf(fout,"%d %d\n",u_newID,v_newID);
        wholeSocialLink[u_newID].insert(v_newID);
        wholeSocialLink[v_newID].insert(u_newID);
    }
    string linkNewPath = getSocialInputPath(LINK);
    FILE *fout;
    fout = fopen(linkNewPath.c_str(), "w");
    unordered_map<int,set<int>>::iterator iter ;iter=wholeSocialLink.begin();
    while(iter!=wholeSocialLink.end()){
        int u = iter->first;
        set<int> neighbors = iter->second;
        for(int v: neighbors){
            fprintf(fout,"%d\t%d\n",u,v);
        }
        iter++;
    }
    fclose(fout);
    cout<<"Update .friendship 文件！"<<endl;

    //3. check in 文件更新
    checkinPath = getSocialInputPath(CHECKIN);
    fout = fopen(checkinPath.c_str(), "w");
    for(int i=0;i<ckCollection.size();i++){
        CheckInRecord checkIn = ckCollection[i];
        int check_user_newID = online_user_id2IdxMap[checkIn.u];
        fprintf(fout,"%d %d %d\n",checkIn.p,check_user_newID, checkIn.count);
    }
    fclose(fout);
    cout<<"Update .checkin 文件！"<<endl;

    return  0;
}



void serialize_SocialLink(){
    clock_t startTime, endTime;
    //从（文本）文件中读取社交网络数据内容
    loadFriendShipData();

    startTime = clock();
    string indexOutputPre = getSocialLinkBinaryInputPath();
    stringstream ss;
    ss<<indexOutputPre;
    string indexOutputFilePath = ss.str();
    cout<<"*****************序列化输出 SocialLink 数据...";

    SocialGraphBinary socialBin;
    socialBin.setGraphData(wholeSocialLinkMap, online_users);
    //socialBin.setSocialLink(friendshipMap,followedMap,friendShipTable);
    socialBin.setSocialLink(friendshipMap,followedMap);


    serialization::outputIndexToBinaryFile<SocialGraphBinary>(socialBin,indexOutputFilePath);

    endTime = clock();
    cout<<"输出完毕! 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s*****************" << endl;



}


void serialize_CheckIn(){
    clock_t startTime, endTime;
    //从（文本）文件中读取社交网络数据内容
    loadCheckinData();
    startTime = clock();
    string indexOutputPre = getCheckInBinaryInputPath();
    stringstream ss;
    ss<<indexOutputPre;
    string dataOutputPath = ss.str();
    cout<<"*****************序列化输出 Check in 数据...";

    CheckInBinary checkInBin;
    checkInBin.setCheckInData(userCheckInfoList,poiCheckInfoList);


    serialization::outputIndexToBinaryFile<CheckInBinary>(checkInBin,dataOutputPath);

    endTime = clock();
    cout<<"输出完毕! 用时："<< (double)(endTime - startTime) / CLOCKS_PER_SEC  << "s*****************" << endl;
}



void CleanLasVegasTextInfo(){ // leaf node level,去掉了之前多余累赘的 address_v2p对应的索引与数据文件
    cout<<"Begin Modify LasVegasPOITextInfo!"<<endl;
    LOAD_START
    // load PoiMap
    printf("LOADING POI INFO...");
    ifstream finPoi;
    string poiPath = getObjectInputPath(BUSINESS);  //getObjectNewInputPath(BUSINESS);   //
    string str;
    finPoi.open(poiPath.c_str());

    map<int,set<int>> inv_list;
    set<int> current_words;
    while (getline(finPoi, str)) {
        istringstream tt(str);
        int POI_id, POI_Ni, POI_Nj;
        float POI_dist, POI_dis;
        tt >> POI_id >> POI_Ni >> POI_Nj >> POI_dist >> POI_dis;
        float exchangeNum = (int) (POI_dist * WEIGHT_INFLATE_FACTOR);
        POI_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (POI_dis * WEIGHT_INFLATE_FACTOR);
        POI_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;


        vector<int> uKey;
        set<int> keywordSet;
        int keyTmp;
        //cout<<"keyword:";
        int term = -1;
        while (tt >> term) {
            if (0 > term) {
                break;
            }
            uKey.push_back(term); keywordSet.insert(term);
            current_words.insert(term);
            inv_list[term].insert(POI_id);
            //cout<<",w"<<term;

        }
        // cout<<endl;

        POI tmpPoi;
        tmpPoi.id = POI_id;
        tmpPoi.Ni = min(POI_Ni,POI_Nj);
        tmpPoi.Nj = max(POI_Ni,POI_Nj);
        tmpPoi.dist = POI_dist;
        tmpPoi.dis = POI_dis;
        tmpPoi.category = 1;
        tmpPoi.keywords = uKey;
        tmpPoi.keywordSet = keywordSet;
        POIs.push_back(tmpPoi);

        poiCnt++;
    }
    finPoi.close();
    cout << " DONE! POI #: " << poiCnt << endl;
    LOAD_END
    LOAD_PRINT(" load POI ")

    map<int,set<int>> ::iterator iter;
    priority_queue<CardinalityEntry> Q;
    for(iter= inv_list.begin();iter!= inv_list.end();iter++){
        int term_id = iter->first;
        set<int> posting_list = iter->second;
        int term_posting_size = posting_list.size();
        CardinalityEntry entry(term_id,term_posting_size);
        Q.push(entry);
    }

    ////第一步排序
    int i = 1; map<int, int> termNewIDMap;
    while(!Q.empty()){
        CardinalityEntry entry = Q.top();
        Q.pop();
        int current_term = entry.id;
        int current_size = entry.size;
        termNewIDMap[current_term] = i;
        i ++;
        //getchar();
    }

    cout<<"End Modify LasVegas POI TextInfo!"<<endl;

    //清理poi中部分冗余的单词
    for(int i=0;i<POIs.size();i++){
        POI p = POIs[i];
        set<int> new_keys;
        for(int term: p.keywordSet){
            int newKey = termNewIDMap[term];
            new_keys.insert(newKey);

        }
        int current_textual_size = (int)random()%textual_capacity;
        current_textual_size = max(1,current_textual_size);
        set<int> refine_keys;refine_keys.clear();
        for(int term:new_keys){
            if(refine_keys.size()>= current_textual_size) continue;
            refine_keys.insert(term);
        }
        POIs[i].keywordSet = refine_keys;

    }

    ////2.重新排序
    map<int,set<int>> inv_listRefined;
    map<int, int> refinedNewIDMap;
    for(int i=0;i<POIs.size();i++){
        POI p = POIs[i];
        for(int term:p.keywordSet){
            inv_listRefined[term].insert(p.id);
        }
    }
    priority_queue<CardinalityEntry> Q_r;
    for(iter= inv_listRefined.begin();iter!= inv_listRefined.end();iter++){
        int term_id = iter->first;
        set<int> posting_list = iter->second;
        int term_posting_size = posting_list.size();
        CardinalityEntry entry(term_id,term_posting_size);
        Q_r.push(entry);
    }
    i = 1;
    while(!Q_r.empty()){
        CardinalityEntry entry = Q_r.top();
        Q_r.pop();
        int current_term = entry.id;
        int current_size = entry.size;
        if(i%50==0){
            //getchar();
            cout<<"rank "<<i<<": t"<<current_term<<", size="<<current_size<<endl;
        }
        refinedNewIDMap[current_term] = i;
        i ++;
    }


    cout<<"Output new POIInfo ...";

    ofstream foutPOI;
    string newPOIPath = getObjectNewInputPath(BUSINESS);
    foutPOI.open(newPOIPath);
    int total_keySize =0;
    for(int i=0;i<POIs.size();i++){
        POI p = POIs[i];
        foutPOI << p.id << " " << p.Ni << " " << p.Nj << " " << p.dist << " " << p.dis << " ";
        set<int> new_keys;  set<int> refine_keys;
        for(int term: p.keywordSet){
            int term_refinement = refinedNewIDMap[term];
            if(!term_refinement>0) continue;
            refine_keys.insert(term_refinement);
        }

        for(int term: refine_keys){
            foutPOI << term << " ";
        }
        total_keySize += refine_keys.size();

        foutPOI << endl;
    }

    foutPOI.close();
    float average_keySize = total_keySize / (1.0*POIs.size());
    printf("DONE! %f keywords per POI \n",average_keySize);

    printf("LOADING User INFO...");
    ifstream finUsr;
    string userPath = getObjectInputPath(USER);
    finUsr.open(userPath.c_str());
    while (getline(finUsr, str)) {
        istringstream tt(str);
        int usr_id, usr_Ni, usr_Nj;
        float usr_dist, usr_dis;
        tt >> usr_id >> usr_Ni >> usr_Nj >> usr_dist >> usr_dis;
        float exchangeNum = (int) (usr_dist * WEIGHT_INFLATE_FACTOR);
        usr_dist = exchangeNum / WEIGHT_INFLATE_FACTOR;
        exchangeNum = (int) (usr_dis * WEIGHT_INFLATE_FACTOR);
        usr_dis = exchangeNum / WEIGHT_INFLATE_FACTOR;

        vector<int> uKey;
        set<int> keywordSet;
        int keyTmp;
        //cout<<"keyword:";
        int term = -1;
        while (tt >> term) {
            if (0 > term) {
                break;
            }
            int term_new = termNewIDMap[term];
            int term_refine = refinedNewIDMap[term_new];
            if(term_refine>0){
                uKey.push_back(term_refine); keywordSet.insert(term_refine);
            }

        }
        if(!uKey.size()>0){
            int add_term = random()% 5 +1;
            uKey.push_back(add_term);
            keywordSet.insert(add_term);
        }
        User tmpUsr;
        tmpUsr.id = usr_id;
        tmpUsr.Ni = min(usr_Ni,usr_Nj);
        tmpUsr.Nj = max(usr_Ni,usr_Nj);
        tmpUsr.dist = usr_dist;
        tmpUsr.dis = usr_dis;
        tmpUsr.keywords = uKey;
        tmpUsr.keywordSet = keywordSet;
        Users.push_back(tmpUsr);

        usrCnt++;

    }
    finUsr.close();
    cout << " DONE! USER #: " << usrCnt << endl;
    LOAD_END
    LOAD_PRINT(" load User ")

    cout<<"Output new UserInfo ...";

    ofstream foutUser;
    string newUserPath = getObjectNewInputPath(USER);
    foutUser.open(newUserPath);
    total_keySize = 0;
    for(int i=0;i<Users.size();i++){
        User u = Users[i];
        foutUser << u.id << " " << u.Ni << " " << u.Nj << " " << u.dist << " " << u.dis << " ";
        int count = 0;
        for(int term: u.keywordSet){
            if(count>6) continue;
            if(term>=20000) continue;
            foutUser << term << " ";
            count++;
            total_keySize++;

        }
        foutUser << endl;
    }

    foutUser.close();
    average_keySize = total_keySize / (1.0*Users.size());
    printf("DONE! %f keywords per user \n",average_keySize);





}



/*---------------------------------------控制台程序参数获取-------------------------------*/

int getAlgChose(std::string alg_chose){  //reverse knn processing method
    if(alg_chose == constants::BASE_ALG)
        return BASELINE;
    else if(alg_chose== constants::APPRO_ALG)
        return APPROXIMATE;
    else if(alg_chose== constants::HEURISTIC_ALG)
        return HEURISTIC;

    else if(alg_chose=="onehop")
        return ONEHOP;
    else if(alg_chose=="twohop")
        return TWOHOP;
    else if(alg_chose=="rele")
        return RELEVANCE;
    else if(alg_chose=="cardi")
        return CARDINALITY;
    else if(alg_chose=="influ")
        return INFLUENCER;
    else if(alg_chose=="random")
        return RANDOM;

    else if(alg_chose=="sep")
        return SEPRATE;
    else if (alg_chose== constants::GRP_ALG)
        return GROUP_BATCH;

    else if (alg_chose==constants::NVD_ALG)  //直接对每个相关user进行nvd剪枝
        return NVD_BATCH;
}


float getObjectSizeRatioSettings(std::string ratio){
    if(ratio == "0.2")
        return 0.2;
    else if(ratio == "0.4")
        return 0.4;
    else if(ratio == "0.6")
        return 0.6;
    else if(ratio == "0.7")
        return 0.7;
    else if(ratio == "0.8")
        return 0.8;
    else if(ratio == "0.9")
        return 0.9;
    else if(ratio == "1.0")
        return 1.0;
}



float getAlphaSettings(std::string alpha_value){
    if(alpha_value == "0.0")
        return 0.0;
    else if(alpha_value == "0.2")
        return 0.2;
    else if(alpha_value == "0.4")
        return 0.4;
    else if(alpha_value == "0.6")
        return 0.6;
    else if(alpha_value == "0.8")
        return 0.8;
}

int getBudgetSettings(std::string qn_value){
    if(qn_value == "1")
        return 1;
    else if(qn_value == "3")
        return 3;
    else if(qn_value == "5")
        return 5;
    else if(qn_value == "7")
        return 7;
    else if(qn_value == "9")
        return 9;
    else
        return DEFAULT_BUDGET;
}

int getQNSettings(std::string qn_value){
    if(qn_value == "25")
        return 25;
    else if(qn_value == "50")
        return 50;
    else if(qn_value == "100")
        return 100;
    else if(qn_value == "200")
        return 200;
    else if(qn_value == "400")
        return 400;
    else
        return DEFAULT_NUMBEROFSTORE;
}

int getOcSettings(std::string qn_value){
    if(qn_value == "10")
        return  10;
    else if(qn_value == "20")
        return 20;
    else if(qn_value == "40")
        return 40;
    else if(qn_value == "60")
        return 60;
    else if(qn_value == "80")
        return 80;
    else if(qn_value == "100")
        return 100;
    else if(qn_value == "160")
        return 160;
    else if(qn_value == "200")
        return 200;
    else if(qn_value == "1")
        return 1;
    else if(qn_value == "2")
        return 2;
    else if(qn_value == "3")
        return 3;
    else if(qn_value == "4")
        return 4;
    else if(qn_value == "5")
        return 5;
    else if(qn_value == "6")
        return 6;
    else if(qn_value == "7")
        return 7;
    else if(qn_value == "8")
        return 8;
    else if(qn_value == "9")
        return 9;
    else
        return DEFAULT_NUMBEROFSTORE;
}

int getKQSettings(std::string KQ_value){
    if(KQ_value == "10")
        return  10;
    else if(KQ_value == "20")
        return 20;
    else if(KQ_value == "40")
        return 40;
    else if(KQ_value == "60")
        return 60;
    else if(KQ_value == "80")
        return 80;
    else
        return DEFAULT_NUMBEROFSTORE;
}

float getAccuracySettings(std::string ratio){
    if(ratio == "0.02")
        return  0.02;
    else if(ratio == "0.05")
        return  0.05;
    else if(ratio == "0.1")
        return 0.1;
    else if(ratio == "0.15")
        return 0.15;
    else if(ratio == "0.2")
        return 0.2;
    else if(ratio == "0.3")
        return 0.3;
    else
        return 0.1;
}




int getOkSettings(std::string k_value){
    if(k_value == "10")
        return  10;
    else if(k_value == "15")
        return 15;
    else if(k_value == "20")
        return 20;
    else if(k_value == "25")
        return 25;
    else if(k_value == "30")
        return 30;
    else
        return DEFAULT_NUMBEROFSTORE;
}


int getQueryGroupChose(std::string group_id){
    if(group_id == "1")
        return  1;
    else if(group_id == "2")
        return 2;
    else if(group_id == "3")
        return 3;
    else if(group_id == "4")
        return 4;

    else
        return 1;
}

/*int getPopularityChose(std::string popularity_chose){
    if(popularity_chose == "g1")
        return 1;
    else if(popularity_chose== "g2")
        return 2;
    else if(popularity_chose== "g3")
        return 3;
    else if(popularity_chose== "g4")
        return 4;
    else
        return 1;
}*/

int getPopularityChose(std::string popularity_chose){
    if(popularity_chose == "pop")
        return Popular;
    else if(popularity_chose== "mid")
        return Mid;
    else if(popularity_chose== "rare")
        return Rare;
}


void loadRoadNetworkIndexData(){


    ConfigType setting;
    AddConfigFromFile(setting,"../map_exp");

    //导入gtree
    initGtree();

    initialPHL();

    // 加载 edge information
    loadEdgeMap_light();



}



void loadData(){


    ConfigType setting;
    AddConfigFromFile(setting,"../map_exp");

    //导入gtree
    initGtree();

    cout<<"gtree node size="<<GTree.size()<<endl;
    //getchar();

    //导入phl索引
    initialPHL();
    
    // 加载 edge information
    loadEdgeMap();

#ifndef DiskAccess
    loadUPMap();
#else
    initialO2UData();
#endif


    //读取社交网络图信息
    //loadFriendShipData();
    loadFriendShipBinaryData();

    //读取用户在兴趣点的check-in
    //loadCheckinRawData();
    //loadCheckinBinaryData();
    loadCheckinBinaryDataAndRecord();

    //读取用户朋友的最大签到
    loadUsrMaxSocial();




}


void loadData_light(){


    ConfigType setting;
    AddConfigFromFile(setting,"../map_exp");

    //导入gtree
    initGtree();

    cout<<"gtree node size="<<GTree.size()<<endl;
    //getchar();

    //导入phl索引
    initialPHL();

    // 加载 edge information
    loadEdgeMap_light();

    //读取社交网络图信息
    loadFriendShipData();

    //读取用户在兴趣点的check-in
       ////loadCheckinData();

    //读取用户朋友的最大签到
    loadUsrMaxSocial();

#ifndef DiskAccess
    loadUPMap();
#else
    initialO2UData();
#endif




}


void loadData_test(){


    ConfigType setting;
    AddConfigFromFile(setting,"../map_exp");

    //导入gtree
    initGtree();

    cout<<"gtree node size="<<GTree.size()<<endl;
    //getchar();

    //导入phl索引
    initialPHL();

    // 加载 edge information
    ///loadEdgeMap_light();

    //读取社交网络图信息
    ///loadFriendShipData();

    //读取用户在兴趣点的check-in
    ////loadCheckinData();

    //读取用户朋友的最大签到
    loadUsrMaxSocial();

#ifndef DiskAccess
    loadUPMap();
#else
    initialO2UData();
#endif




}

void loadData_exp(){


    ConfigType setting;
    AddConfigFromFile(setting,"../map_exp");

    //导入gtree
    initGtree();

    cout<<"gtree node size="<<GTree.size()<<endl;
    //getchar();

    //导入phl索引
    initialPHL();

    // 加载 edge information
    loadEdgeMap_light();

#ifndef DiskAccess
    loadUPMap();
#else
    initialO2UData();
#endif


    //读取社交网络图信息
    //loadFriendShipData();
    loadFriendShipBinaryData();

    //读取用户在兴趣点的check-in
    //loadCheckinRawData();
    //loadCheckinBinaryData();
    loadCheckinBinaryDataAndRecord();

    //读取用户朋友的最大签到
    loadUsrMaxSocial();




}

void loadDisBound_Other() {  //float**& MaxDisBound, float**& MinDisBound
    cout<<"LOADING DISBOUND..."<<endl;
    ifstream inDisBoundFile;
    //stringstream tt;
    //tt<<FILE_DISBOUND;
    string disBoundFile = getRoadInputPath(DISBOUND);
    //inDisBoundFile.open(tt.str());
    inDisBoundFile.open(disBoundFile);

    string str;
    int cnt = 0;
    unordered_map<int, float> userScoreMap;

    //MaxDisBound = *new map<int,map<int,float>>();
    //MinDisBound = *new map<int,map<int,float>>();
    //新建 距离上下接的 二维数组
    int node_size = GTree.size();
    MaxDisBound = new float*[node_size];
    MinDisBound = new float*[node_size];
    for(int i=0;i<node_size;i++){
        MaxDisBound[i] = new float[node_size];
        MinDisBound[i] = new float[node_size];
    }


    while (getline(inDisBoundFile, str)) {
        istringstream tt(str);
        int node1,node2;
        float minDis,maxDis;
        tt >> node1 >> node2 >> maxDis >> minDis;
        MaxDisBound[node1][node2] = maxDis;
        MinDisBound[node1][node2] = minDis;
        MaxDisBound[node2][node1] = maxDis;
        MinDisBound[node2][node1] = minDis;
        //cout<<MaxDisBound[node1][node2]<<endl;
        //getchar();
    }
    inDisBoundFile.close();
    cout<<"LOADING DISBOUND COMPLETE!"<<endl;

}



void loadDisBound_Gowalla() {  //float**& MaxDisBound, float**& MinDisBound
    cout<<"LOADING DISBOUND";
    int cnt = 0;
    unordered_map<int, float> userScoreMap;

    //新建 距离上下接的 二维数组
    int node_size = GTree.size();
    MaxDisBound = new float*[node_size];
    MinDisBound = new float*[node_size];
    for(int i=0;i<node_size;i++){
        MaxDisBound[i] = new float[node_size];
        MinDisBound[i] = new float[node_size];
    }
    for(int i=0;i<node_size;i++){
        for(int j=0;j<node_size;j++){
            MaxDisBound[i][j] = 999999999.9;
            MinDisBound[i][j] = 0.0;
        }
    }
    //从各个文件中load bound

    ////1. 读取 leaf bound
    cout<<" LOADING ";
    ifstream inDisBoundFile;

    string disBoundFile = getRoadInputPath(DISBOUND);
    cout<<"from "<<disBoundFile<<" ...";
    inDisBoundFile.open(disBoundFile);

    string str;
    while (getline(inDisBoundFile, str)) {
        istringstream tt2(str);
        int node1,node2;
        float minDis,maxDis;
        tt2 >> node1 >> node2 >> maxDis >> minDis;
        MaxDisBound[node1][node2] = maxDis;
        MinDisBound[node1][node2] = minDis;
        MaxDisBound[node2][node1] = maxDis;
        MinDisBound[node2][node1] = minDis;
        //cout<<MaxDisBound[node1][node2]<<endl;
        //getchar();
    }
    inDisBoundFile.close();
    cout<<"COMPLETE!"<<endl;


    ////2. 读取 leaf upper bound
    stringstream tt;
    tt.str("");
    tt << getRoadInputPath(DISBOUND)<<"-leafupper";
    disBoundFile = tt.str();
    cout<<"LOADING DISBOUND from "<<disBoundFile<<" ...";
    inDisBoundFile.open(disBoundFile);

    while (getline(inDisBoundFile, str)) {
        istringstream tt2(str);
        int node1,node2;
        float minDis,maxDis;
        tt2 >> node1 >> node2 >> maxDis >> minDis;
        MaxDisBound[node1][node2] = maxDis;
        MinDisBound[node1][node2] = minDis;
        MaxDisBound[node2][node1] = maxDis;
        MinDisBound[node2][node1] = minDis;
    }
    inDisBoundFile.close();
    cout<<"COMPLETE!"<<endl;


    ////3. 读取 upper bound
    tt.str("");
    tt << getRoadInputPath(DISBOUND)<<"-upper";
    disBoundFile = tt.str();
    cout<<"LOADING DISBOUND from "<<disBoundFile<<" ...";
    inDisBoundFile.open(disBoundFile);

    while (getline(inDisBoundFile, str)) {
        istringstream tt2(str);
        int node1,node2;
        float minDis,maxDis;
        tt2 >> node1 >> node2 >> maxDis >> minDis;
        MaxDisBound[node1][node2] = maxDis;
        MinDisBound[node1][node2] = minDis;
        MaxDisBound[node2][node1] = maxDis;
        MinDisBound[node2][node1] = minDis;
    }
    inDisBoundFile.close();
    cout<<"COMPLETE!"<<endl;


}



void loadDisBound() {  //float**& MaxDisBound, float**& MinDisBound
#ifdef  Gowalla
    loadDisBound_Gowalla();
#else
    loadDisBound_Other();
#endif

}


void clearDisBoundArray(){
    cout<<"clearDisBoundArray..."<<endl;
    for(int i = 0; i < GTree.size(); i++){
        delete []MaxDisBound[i];
        delete []MinDisBound[i];
    }

    delete []MaxDisBound;
    delete []MinDisBound;
    cout<<"clearDisBoundArray COMPLETE!"<<endl;
}


void loadDisBound_Transfer2Binary(){
    loadDisBound();
    char binaryFile_Name[255];
    sprintf(binaryFile_Name,"%s._binary",FILE_DISBOUND);

    remove(binaryFile_Name); // remove existing file
    FILE* binaryFile = fopen(binaryFile_Name,"w+");

    int node_size = GTree.size();
    printf("begin transfering... \n");
    for(int i=1;i<node_size;i++){
        for(int j=1;j<node_size;j++){
            if(i>j) continue;
            float _max = MaxDisBound[i][j];
            float _min = MinDisBound[i][j];
            fwrite(&i, 1, sizeof(int),binaryFile);
            fwrite(&j, 1, sizeof(int),binaryFile);
            fwrite(&_max, 1, sizeof(float),binaryFile);
            fwrite(&_min, 1, sizeof(float),binaryFile);

        }
    }
    fclose(binaryFile);
    printf("Transfer to Binary file Success!\n");

}

//有问题！
void loadDisBound_FromBinary(){
    cout<<"LOADING DISBOUND From BinaryFile..."<<endl;

    //MaxDisBound = *new map<int,map<int,float>>();
    //MinDisBound = *new map<int,map<int,float>>();
    //新建 距离上下接的 二维数组
    int node_size = GTree.size();
    MaxDisBound = new float*[node_size];
    MinDisBound = new float*[node_size];
    for(int i=0;i<node_size;i++){
        MaxDisBound[i] = new float[node_size];
        MinDisBound[i] = new float[node_size];
    }

    char binaryFileName[255];
    sprintf(binaryFileName, "%s.binary", FILE_DISBOUND);
    FILE* binaryFile = fopen(binaryFileName,"rb");
    printf("disbound binaryFile open!\n");
    int node1; int node2; float _max; float _min;
    while (fread(&node1, sizeof(int), 1, binaryFile)) {
        fread(&node2, sizeof(int), 1, binaryFile);
        fread(&_max, sizeof(float), 1, binaryFile);
        fread(&_min, sizeof(float), 1, binaryFile);
        if(node1!=node2){
            MaxDisBound[node1][node2] = _max;
            MinDisBound[node1][node2] = _min;
            MaxDisBound[node2][node1] = _max;
            MinDisBound[node2][node1] = _min;
        } else{
            MaxDisBound[node1][node2] = _max;
            MinDisBound[node1][node2] = _min;
        }

    }
    fclose(binaryFile);
    cout<<"LOADING DISBOUND From BinaryFile Complete!"<<endl;
    cout<<"1242 1277(max)="<<MaxDisBound[1242][1277]<<endl;
    getchar();

}





/*---------------------------MaxInfBRGSTkNN query processing---------------------------*/

void MaxInfBRGSTkNNQ_Base_GELF(int argc, char* argv[], vector<int>& users, unordered_map<int,vector<int>>& wholeSocialLinkMap, vector<int>& stores, BatchResults& results, int b, InfluenceModels* L){
    cout<<"In function: MaxInfBRGSTkNNQ_Batch_GELF..."<<endl;

    L->doInitial(argc, argv, UserID_MaxKey, users, wholeSocialLinkMap, stores, CELF);    //初始化社交影响力评估模块


    L->batchThenCelf(stores, results, b);
    //cout<<endl;

}

void MaxInfBRGSTkNNQ_Base_PMC(int argc, char* argv[], vector<int>& users, unordered_map<int,vector<int>>& wholeSocialLinkMap, vector<int>& stores, BatchResults& results, int b, InfluenceModels* L){
    cout<<"In function: MaxInfBRGSTkNNQ_Base_PMC..."<<endl;

    SOCIAL_START

    L->doInitial(argc, argv, UserID_MaxKey, users, wholeSocialLinkMap, stores, PMC);    //初始化社交影响力评估模块

    int R = 200;
    L->seperateThenPMC(stores, results, b, R);

    SOCIAL_END
    SOCIAL_PRINT("Influence evaluation by PMC:");
    //cout<<endl;

}


void MaxInfBRGSTkNNQ_Batch_IMM(int argc, char* argv[], vector<int>& users, unordered_map<int,vector<int>>& wholeSocialLinkMap, vector<int>& stores, BatchResults& results, int b, InfluenceModels* L){
    cout<<"In function: MaxInfBRGSTkNNQ_Batch_IMM..."<<endl;

    TIME_TICK_START
    L->doInitial(argc, argv, UserID_MaxKey, users, wholeSocialLinkMap, stores, IMM);    //初始化社交影响力评估模块
    L->batchThenIMM(stores, results, b);
    TIME_TICK_END
    TIME_TICK_PRINT("IMM's runtime")

}


int MaxInfBRGSTkNNQ_Batch_OPIMC(int argc, char* argv[], float accuracy, vector<int>& users, unordered_map<int,vector<int>>& wholeSocialLinkMap, vector<int>& stores, BatchResults& results, int b, InfluenceModels* L,ofstream& exp_results){
    cout<<"In function: MaxInfBRGSTkNNQ_Batch_OPIMC..."<<endl;

    SOCIAL_START
    //L->doInitial(argc, argv, accuracy, UserID_MaxKey, users, wholeSocialLinkMap, stores, OPC);    //初始化社交影响力评估模块
    if(L->isInitial==false)
        L->doInitial(argc, argv, accuracy, UserID_MaxKey,users, wholeSocialLinkMap, stores, OPC);

    int rrSet_num = L->batchThenOPIMC(stores, results, b,exp_results);


    SOCIAL_END
    SOCIAL_PRINT("Influence evaluation by OPIMC runtime:");

    return rrSet_num;

}



int MaxInfBRGSTkNNQ_Batch_HYBRID(int argc, char* argv[], float accuracy, vector<int>& users, unordered_map<int,vector<int>>& wholeSocialLinkMap, vector<int>& stores, BatchResults& results, int b, InfluenceModels* L,ofstream& exp_results){
    cout<<"In function: MaxInfBRGSTkNNQ_Batch_Hybrid..."<<endl;

    SOCIAL_START
    if(L->isInitial==false)
        L->doInitial(argc, argv, accuracy, UserID_MaxKey,users, wholeSocialLinkMap, stores, OPC);
    int rrSet_num = L->batchThenOPIMC_Hybrid(stores, results, b,exp_results);

    SOCIAL_END
    SOCIAL_PRINT("Influence evaluation by Hybrid runtime:");
    return rrSet_num;

}


class BatchVerifyInfEntry {
public:
    int u_id;
    float u_inf;
    double rk_current;
    vector<int> key;
    User u;

    BatchVerifyInfEntry(int id, float inf, double r, User user) : u_id(id), u_inf(inf){
        rk_current = r; u = user;
    }
    ~BatchVerifyInfEntry(){
        //cout<<"Result被释放"<<endl;
    }
    void printResult() {
        printf("id: %d , inf= %f\n", u_id, u_inf);
        printf("key: ");

    }

    bool operator<(const BatchVerifyInfEntry &a) const {
        return a.u_inf > u_inf;
    }
};

class RelevanceDencityEntry {
public:
    int poi_id;
    float density;


    RelevanceDencityEntry(int id, float _density) : poi_id(id), density(_density){
    }
    ~RelevanceDencityEntry(){
        //cout<<"Result被释放"<<endl;
    }
    void printResult() {
        printf("id: %d , inf= %f\n", poi_id, density);

    }

    bool operator<(const RelevanceDencityEntry &a) const {
        return a.density > density;
    }
};



#define HeuristicSelection 1
vector<int> batchFilterThenHeuristic (int K, float alpha, vector<int>& stores, float error_parameter, vector<int>& users, unordered_map<int,vector<int>>& wholeSocialLinkMap,BatchFastFilterResults& filterOutput, int b,InfluenceModels* L, ofstream& exp_results){

    set<BatchVerifyEntry> candidate_User;  //candidate users
    unordered_map<int, priority_queue<ResultCurrent>> candidate_usr_related_store;
    map<int, vector<int>> incremental_results;   //<o_id, {u_id...}
    unordered_map<int, POI> candidate_POIs;


    candidate_User = filterOutput.candidate_User;
    candidate_usr_related_store = filterOutput.candidate_usr_related_store;
    candidate_POIs = filterOutput.candidate_POIs;
    //incremental_results = filterOutput.batch_results;

    ////opp ->checkReachable(stores, results);

    double runtime = 0; double opimc_time; double hop_time;
    clock_t start_time; clock_t end_time; clock_t transfer_endTime;
    start_time = clock();
    //先初始话社交评估模块
    L->doInitial4Heuris(error_parameter, UserID_MaxKey, users, wholeSocialLinkMap, stores, OPC);    //初始化社交影响力评估模块
    transfer_endTime = clock();
    double _time = (double)(transfer_endTime-start_time)/CLOCKS_PER_SEC*1000;
    printf("transfer runtime: %f ms!\n", _time);

    ////opp ->minePOIsHopBased(b);

    int candidate_size = candidate_User.size();
    vector<ResultDetail> potential_customerDetails;
    vector<int> user_list;
    unordered_map<int,priority_queue<ResultLargerFisrt>> associated_poisMap;
    unordered_map<int,set<int>> associated_poiIDMap;
    priority_queue<BatchVerifyInfEntry> Q_inf;
    float total_inf_upper =0; ////出错了！
    float total_inf_base =0;
    map<int, set<int>> o_brupper;
    for(BatchVerifyEntry ve:candidate_User){
        int u_id = ve.u_id;
        double u_rk_score = ve.rk_current;
        User u = ve.u;

        priority_queue<ResultCurrent> relate_query_queue = candidate_usr_related_store[u_id];
        priority_queue<ResultLargerFisrt> user_associating_list;
        set<int> associated_poiIDSet;
        while(!relate_query_queue.empty()){
            ResultCurrent rc = relate_query_queue.top();
            relate_query_queue.pop();
            //先算出真实gsk评分
            double rel = rc.relevance; double inf = rc.influence;
            double gsk_score = rc.score;
            double score_upper = rc.score_upper;
            double _score = rc.score;
            if(gsk_score < u_rk_score)  //jins
                continue;
            int o_id = rc.o_id;

            ResultLargerFisrt rlf(rc.o_id,gsk_score,rel, inf);
            user_associating_list.push(rlf);
            associated_poiIDSet.insert(o_id);
            o_brupper[o_id].insert(u_id);
            //max_score = max(max_score,gsk_score);
        }
        if(user_associating_list.size()>0){
            float u_inf = 1.0; //wholeSocialLinkMap[usr_id].size();
            for(int v: wholeSocialLinkMap[u_id]){   //0507
                //cout<<"followedMap["<<v<<"].size="<<followedMap[v].size()<<endl;  //1-hop neighbor
                float activate_v = 1.0/ followedMap[v].size();
                //cout<<"activate_v1="<<activate_v<<endl;
                u_inf += activate_v;
                for(int v_neighbor: wholeSocialLinkMap[v]){  //2-hop neighbor
                    //cout<<"followedMap["<<v_neighbor<<"].size="<<followedMap[v_neighbor].size()<<endl;
                    float activate_vn = activate_v * (1.0/ followedMap[v_neighbor].size());
                    //cout<<"activate_vn="<<activate_vn<<endl;
                    u_inf += activate_vn;
                    //cout<<"u_inf="<<u_inf<<endl;
                }
                //cout<<"u_inf="<<u_inf<<endl;
            }
            user_list.push_back(u_id);

            BatchVerifyInfEntry inf_entry(u_id,u_inf,u_rk_score,u);
            Q_inf.push(inf_entry);
            total_inf_upper += u_inf;
            associated_poisMap[u_id] =user_associating_list;
            associated_poiIDMap[u_id] = associated_poiIDSet;
        }
    }

    /*for(int store: stores){
        cout<<"store"<<store<<"的RkGSKQ, 结果个数="<<o_brupper[store].size()<<",内容："<<endl;
        printSetElements(o_brupper[store]);
    }*/



    //按照用户社交影响力优先级，对候选用户进行最后的验证评估
    vector<BatchVerifyInfEntry> remained_candidate_usersInfo;
    vector<int> remained_candidate_usersID;
    vector<User> remained_users;   vector<float> remained_userInf;
    while(Q_inf.size()>0){
        //for(int i=0;i< user_list.size();i++){
        BatchVerifyInfEntry inf_entry = Q_inf.top();
        Q_inf.pop();

        int usr_id = inf_entry.u_id;
        User _user = inf_entry.u;
        float _inf = inf_entry.u_inf;
        remained_candidate_usersInfo.push_back(inf_entry);
        remained_candidate_usersID.push_back(usr_id);
        remained_users.push_back(_user);
        remained_userInf.push_back(_inf);


    }
    double e = exp(1);
    double miu = 0.2; double relative_ratio = (1-1.0/e)-0.2; //1-1.0/e; 0.8; //1-1.0/e; 0.5
    double delta_everyroud = (1.0-miu)/20.0;
    int check_count = 0;
    for(int i=0;i<remained_candidate_usersInfo.size();i++){
        BatchVerifyInfEntry inf_entry = remained_candidate_usersInfo[i];
        User  user = remained_users[i];
        float u_inf = remained_userInf[i];
        //printf("将检测第%d个候选用户， u_inf=%f!",(1+i),u_inf);
        int usr_id = user.id;

        priority_queue<ResultLargerFisrt> associated_list = associated_poisMap[usr_id];
        float Rk_u = -1;
        TopkQueryCurrentResult topk_r = TkGSKQ_NVD(user, K, DEFAULT_A,alpha);//,poiMark,poiADJMark);
        Rk_u = topk_r.topkScore;
        //对Lu中的各个query object进行评估
        float gsk_score_max =  associated_list.top().score;
        if(gsk_score_max >Rk_u){
            total_inf_base += u_inf;
            //printf("用户为潜在用户！\n");
        }
        else{
            total_inf_upper = total_inf_upper - u_inf;
            //printf("用户不是潜在用户！\n");
        }
        while(!associated_list.empty()){
            //for(int o_id: associated_list){
            ResultLargerFisrt associated_entry = associated_list.top();
            associated_list.pop();
            int o_id = associated_entry.o_id;
            POI poi = candidate_POIs[o_id];
            float gsk_score = associated_entry.score;
            if(gsk_score > Rk_u){
                incremental_results[o_id].push_back(usr_id);

            }
            else break;

        }
        float rate = total_inf_base/total_inf_upper;
        //printf("total_inf_base=%f,total_inf_upper=%f, rate=%f\n", total_inf_base,total_inf_upper,rate);
        if(rate>= miu+check_count*delta_everyroud){  //进行 poi seed的选择  //total_inf_base/total_inf_upper>=miu  //remained_candidate_usersInfo.size()-1==i
#ifdef HeuristicSelection
            printf("当前检测了%d个候选用户， u_inf=%f,",(1+i),u_inf);
            printf("total_inf_base=%f,total_inf_upper=%f, rate=%f\n", total_inf_base,total_inf_upper,rate);
#endif
            check_count++;
            clock_t  select_startTime, select_endTime;
            select_startTime = clock();


            L->opp->setPUConnection(stores, incremental_results);

            ////更新 PUConnection_upper 关系 (理想情况下: 已有结果 + 全部剩余结果)
            int next = i+1;
            L->opp->setPUConnection_Upper(stores,incremental_results,next,remained_candidate_usersID,associated_poiIDMap);

            float inf_base = L->opp ->minePOIsHopBased4Heuristic(b);

            float inf_upper = L->opp ->minePOIsHopBased4Heuristic_upper(b);


            double current_ratio = inf_base / inf_upper;


            select_endTime = clock();
            double runTime = (double)(select_endTime-select_startTime)/CLOCKS_PER_SEC*1000;
#ifdef HeuristicSelection
            printf("ratio:%f \n", current_ratio);
            printf("当前下界情况下所选的seeds为:\n");
            printElements(L->opp->infGraph->poiSet);
            printf("checkReachable rumtime: %f ms!, checked user num: %d\n", runTime,i);
#endif

            if(current_ratio> relative_ratio){
                printf("满足相对精度范围（%f>%f）， 验证了%d 个候选用户，验证循环可提前终止!\n",current_ratio,relative_ratio,(i+1));
                exp_results <<"满足相对精度范围,此时验证了"<<(i+1)<<"个候选用户，验证循环可提前终止!"<<endl;
                break;
            } else{
                //getchar();
            }

        }
    }

    end_time = clock();
    runtime = (double)(end_time- start_time);
    //printBatchRkGSKQResults(stores, batch_results);
    cout << "batchFilter Then Heuristic seletion runtime : " << runtime / CLOCKS_PER_SEC * 1000 << "ms" << endl;
    vector<int> pois_selected = L->opp->infGraph->poiSet;
#ifdef HeuristicSelection
    printf("heuristic 最终所选seeds为:\n");
    printElements(pois_selected);
#endif
    exp_results <<"heuristic 最终所选seeds为:"<<endl;

    return pois_selected;

}

class CardinarlityEntry {
public:
    int poi_id;
    float cardinality;


    CardinarlityEntry(int id, float _car) : poi_id(id), cardinality(_car){

    }
    ~CardinarlityEntry(){
        //cout<<"Result被释放"<<endl;
    }
    void printResult() {
        printf("id: %d , cardinality= %f\n", poi_id, cardinality);

    }

    bool operator<(const CardinarlityEntry &a) const {
        return a.cardinality > cardinality;
    }
};

vector<int> poiSelectionRelevanceFirst (vector<int>& stores,int b){
    //统计每个兴趣点与之文本相关的用户个数，并进行排序
    priority_queue<CardinarlityEntry> RelevanceFirst_Q;
    for(int poi_id: stores){
        POI poi = getPOIFromO2UOrgLeafData(poi_id);
        set<int> relevance_userSet;
        for(int term: poi.keywords){
            vector<int> termRelevant_users = get_term_related_objectID(term);
            int _size = termRelevant_users.size();
            for(int u_id:termRelevant_users){
                User usr = getUserFromO2UOrgLeafData(u_id);
                float distance = getDistance_phl(poi.Ni,usr);
                //float distance = phl.Query(usr.Ni,poi.Ni);
                if(distance < 10000){
                    relevance_userSet.insert(u_id);
                }

            }

        }
        int _count = relevance_userSet.size();
        CardinarlityEntry entry(poi_id,_count);
        RelevanceFirst_Q.push(entry);

    }
    vector<int> _poi_select;
    while(_poi_select.size()<b){
        CardinarlityEntry _entry = RelevanceFirst_Q.top();
        RelevanceFirst_Q.pop();
        int _id = _entry.poi_id;
        int _score = _entry.cardinality;
        _poi_select.push_back(_id);
        cout<<"选择p"<<_id<<", score="<<_score<<endl;
    }
    return _poi_select;
}

vector<int> batchFilterThenInfluencer (int K, float alpha, vector<int>& stores, float error_parameter, vector<int>& users, unordered_map<int,vector<int>>& wholeSocialLinkMap,BatchFastFilterResults& filterOutput, int b,InfluenceModels* L, ofstream& exp_results){

    set<BatchVerifyEntry> candidate_User;  //candidate users
    unordered_map<int, priority_queue<ResultCurrent>> candidate_usr_related_store;
    map<int, vector<int>> incremental_results;   //<o_id, {u_id...}
    unordered_map<int, POI> candidate_POIs;


    candidate_User = filterOutput.candidate_User;
    candidate_usr_related_store = filterOutput.candidate_usr_related_store;
    candidate_POIs = filterOutput.candidate_POIs;
    //incremental_results = filterOutput.batch_results;

    ////opp ->checkReachable(stores, results);

    double runtime = 0; double opimc_time; double hop_time;
    clock_t start_time; clock_t end_time; clock_t transfer_endTime;
    start_time = clock();
    //先初始话社交评估模块
    L->doInitial4Heuris(error_parameter, UserID_MaxKey, users, wholeSocialLinkMap, stores, OPC);    //初始化社交影响力评估模块
    transfer_endTime = clock();
    double _time = (double)(transfer_endTime-start_time)/CLOCKS_PER_SEC*1000;
    printf("transfer runtime: %f ms!\n", _time);

    ////opp ->minePOIsHopBased(b);

    int candidate_size = candidate_User.size();
    vector<ResultDetail> potential_customerDetails;
    vector<int> user_list;
    unordered_map<int,priority_queue<ResultLargerFisrt>> associated_poisMap;
    unordered_map<int,set<int>> associated_poiIDMap;
    priority_queue<BatchVerifyInfEntry> Q_inf;
    float total_inf_upper =0; ////出错了！
    float total_inf_base =0;
    map<int, set<int>> o_brupper;
    for(BatchVerifyEntry ve:candidate_User){
        int u_id = ve.u_id;
        double u_rk_score = ve.rk_current;
        User u = ve.u;

        priority_queue<ResultCurrent> relate_query_queue = candidate_usr_related_store[u_id];
        priority_queue<ResultLargerFisrt> user_associating_list;
        set<int> associated_poiIDSet;
        while(!relate_query_queue.empty()){
            ResultCurrent rc = relate_query_queue.top();
            relate_query_queue.pop();
            //先算出真实gsk评分
            double rel = rc.relevance; double inf = rc.influence;
            double gsk_score = rc.score;
            double score_upper = rc.score_upper;
            double _score = rc.score;
            if(gsk_score < u_rk_score)  //jins
                continue;
            int o_id = rc.o_id;

            ResultLargerFisrt rlf(rc.o_id,gsk_score,rel, inf);
            user_associating_list.push(rlf);
            associated_poiIDSet.insert(o_id);
            o_brupper[o_id].insert(u_id);
            //max_score = max(max_score,gsk_score);
        }
        if(user_associating_list.size()>0){
            float u_inf = 1.0; //wholeSocialLinkMap[usr_id].size();
            for(int v: wholeSocialLinkMap[u_id]){   //0507
                //cout<<"followedMap["<<v<<"].size="<<followedMap[v].size()<<endl;  //1-hop neighbor
                float activate_v = 1.0/ followedMap[v].size();
                //cout<<"activate_v1="<<activate_v<<endl;
                u_inf += activate_v;
                for(int v_neighbor: wholeSocialLinkMap[v]){  //2-hop neighbor
                    //cout<<"followedMap["<<v_neighbor<<"].size="<<followedMap[v_neighbor].size()<<endl;
                    float activate_vn = activate_v * (1.0/ followedMap[v_neighbor].size());
                    //cout<<"activate_vn="<<activate_vn<<endl;
                    u_inf += activate_vn;
                    //cout<<"u_inf="<<u_inf<<endl;
                }
                //cout<<"u_inf="<<u_inf<<endl;
            }
            user_list.push_back(u_id);

            BatchVerifyInfEntry inf_entry(u_id,u_inf,u_rk_score,u);
            Q_inf.push(inf_entry);
            total_inf_upper += u_inf;
            associated_poisMap[u_id] =user_associating_list;
            associated_poiIDMap[u_id] = associated_poiIDSet;
        }
    }



    //按照用户社交影响力优先级，对候选用户进行最后的验证评估
    vector<BatchVerifyInfEntry> remained_candidate_usersInfo;
    vector<int> remained_candidate_usersID;
    vector<User> remained_users;   vector<float> remained_userInf;
    int topk_influencer = 200;
    while(remained_candidate_usersInfo.size()<=topk_influencer){
        //for(int i=0;i< user_list.size();i++){
        BatchVerifyInfEntry inf_entry = Q_inf.top();
        Q_inf.pop();
        int usr_id = inf_entry.u_id;
        User _user = inf_entry.u;
        float _inf = inf_entry.u_inf;
        remained_candidate_usersInfo.push_back(inf_entry);
        remained_candidate_usersID.push_back(usr_id);
        remained_users.push_back(_user);
        remained_userInf.push_back(_inf);


    }
    double e = exp(1);
    double miu = 0.2; double relative_ratio = (1-1.0/e)-0.2; //1-1.0/e; 0.8; //1-1.0/e; 0.5
    double delta_everyroud = (1.0-miu)/20.0;
    int check_count = 0;
    for(int i=0;i<remained_candidate_usersInfo.size();i++){
        BatchVerifyInfEntry inf_entry = remained_candidate_usersInfo[i];
        User  user = remained_users[i];
        float u_inf = remained_userInf[i];
        //printf("将检测第%d个候选用户， u_inf=%f!",(1+i),u_inf);
        int usr_id = user.id;

        priority_queue<ResultLargerFisrt> associated_list = associated_poisMap[usr_id];
        float Rk_u = -1;
        TopkQueryCurrentResult topk_r = TkGSKQ_NVD(user, K, DEFAULT_A,alpha);//,poiMark,poiADJMark);
        Rk_u = topk_r.topkScore;
        //对Lu中的各个query object进行评估
        float gsk_score_max =  associated_list.top().score;

        while(!associated_list.empty()){
            //for(int o_id: associated_list){
            ResultLargerFisrt associated_entry = associated_list.top();
            associated_list.pop();
            int o_id = associated_entry.o_id;
            POI poi = candidate_POIs[o_id];
            float gsk_score = associated_entry.score;
            if(gsk_score > Rk_u){
                incremental_results[o_id].push_back(usr_id);

            }
            else break;
        }


    }
    priority_queue<CardinarlityEntry>  Rank_Q;
    map<int, vector<int>>:: iterator iter = incremental_results.begin();
    while(iter!=incremental_results.end()){
        int poi_id = iter->first;
        vector<int> influencer_list = iter->second;
        CardinarlityEntry _entry(poi_id, influencer_list.size());
        Rank_Q.push(_entry);
        iter++;
    }
    vector<int> pois_selected;
    while(pois_selected.size()<b){
        CardinarlityEntry _entry = Rank_Q.top();
        Rank_Q.pop();
        int _storeID = _entry.poi_id;
        int _score = _entry.cardinality;
        pois_selected.push_back(_storeID);
    }



    end_time = clock();
    runtime = (double)(end_time- start_time);
    //printBatchRkGSKQResults(stores, batch_results);
    cout << "influencer first based seletion runtime : " << runtime / CLOCKS_PER_SEC * 1000 << "ms" << endl;
#ifdef HeuristicSelection
    printf("influencer first 最终所选seeds为:\n");
    printElements(pois_selected);
#endif
    exp_results <<"influencer 最终所选seeds为:"<<endl;

    return pois_selected;

}





vector<int> MaxInfBRGSTkNNQ_POISelect_Relevance(int argc, char* argv[], vector<int>& stores, int b, ofstream& exp_results){
    cout<<"In function: MaxInfBRGSTkNNQ_POISelect_Relevance..."<<endl;

    return poiSelectionRelevanceFirst(stores,b);


}



vector<int> MaxInfBRGSTkNNQ_POISelect_Heuristic(int argc, char* argv[], int K, float alpha,float quality_tradoff, vector<int>& users, unordered_map<int,vector<int>>& wholeSocialLinkMap,vector<int>& stores, BatchFastFilterResults& filterOutput, int b, InfluenceModels* L,ofstream& exp_results){
    cout<<"In function: MaxInfBRGSTkNNQ_POISelect_Heuristic..."<<endl;

    return batchFilterThenHeuristic(K, alpha, stores, quality_tradoff, users, wholeSocialLinkMap,filterOutput, b, L,exp_results);


}

vector<int> MaxInfBRGSTkNNQ_POISelect_Influencer(int argc, char* argv[], int K, float alpha,float quality_tradoff, vector<int>& users, unordered_map<int,vector<int>>& wholeSocialLinkMap,vector<int>& stores, BatchFastFilterResults& filterOutput, int b, InfluenceModels* L,ofstream& exp_results){
    cout<<"In function: MaxInfBRGSTkNNQ_POISelect_Influencer..."<<endl;

    return batchFilterThenInfluencer(K, alpha, stores, quality_tradoff, users, wholeSocialLinkMap,filterOutput, b, L,exp_results);


}

vector<int> MaxInfBRGSTkNNQ_POISelect_Random(int argc, char* argv[], int K, float alpha,float quality_tradoff, vector<int>& users, unordered_map<int,vector<int>>& wholeSocialLinkMap,vector<int>& stores, int b, InfluenceModels* L,ofstream& exp_results){
    cout<<"In function: MaxInfBRGSTkNNQ_POISelect_Random..."<<endl;

    vector<int> random_select;
    int _candidate_size = stores.size();
    map<int, bool> selected_flag;
    while(random_select.size()<b){
        int _selected_th = random() % _candidate_size;
        if(selected_flag.count(_selected_th)){
            ;
        }
        else{
            int _selected_ID = stores[_selected_th];
            random_select.push_back(_selected_ID);
        }
    }
    return random_select;

}


void MaxInfBRGSTkNNQ_Batch_HopHeuristic(int argc, char* argv[], float accuracy, vector<int>& users, unordered_map<int,vector<int>>& wholeSocialLinkMap, vector<int>& stores, BatchResults& results, int b, InfluenceModels* L,ofstream& exp_results){
    cout<<"In function: MaxInfBRGSTkNNQ_Batch_HopHeuristic..."<<endl;

    SOCIAL_START
    L->doInitial(argc, argv, accuracy, UserID_MaxKey,users, wholeSocialLinkMap, stores, OPC);    //初始化社交影响力评估模块

    L->batchThenHopbased(stores, results, b,exp_results);

    SOCIAL_END
    SOCIAL_PRINT("Influence evaluation by Hybrid runtime:");

}




typedef priority_queue<pair<int, double>,vector<pair<int, double>>, CompareBySecondDouble> Priority;



void preprocessParameter_settings(int argc, char* argv[], string& _task){
    stringstream ss1;
    int ch;   char flag; int tmp = 0;

    //根据控制台参数进行相应设定
    while ((ch = getopt(argc, argv, "t:j:"))!=-1){
        switch(ch) {
            case 't': //experimant: 选择任务
                flag = ch;
                _task = getPreprocessTask(optarg);
                ss1.clear();
                break;
            case 'j': //experimant: 选择任务
                flag = ch;
                break;

            default :
                break;

        }
    }
}


void parameter_settings(int argc, char* argv[], int& k, int& KQ, int& QN, int& b, float& accuracy, float& alpha, int& popularity, int& alg,int& exp_task, int& topk_method,int& reverse_method,string& _task){
    stringstream ss1;
    int ch;   char flag; int tmp = 0;

    //根据控制台参数进行相应设定
    while ((ch = getopt(argc, argv, "c:k:n:w:a:s:p:e:t:r:v:j:"))!=-1){
        switch(ch) {
            case 'e': //experimant: 选择任务
                flag = ch;
                exp_task = getExpChose(optarg);
                ss1.clear();
                break;

            case 'k':   //top-k参数
                k = getOkSettings(optarg);
                flag = ch;
                break;
            case 'n':  //number of query 代表查询对象个数
                QN = getOcSettings(optarg) ;//getQNSettings(optarg);
                flag = ch;
                break;
            case 'w':  //word count:代表查询关键词个数
                KQ = getKQSettings(optarg);
                flag = ch;
                break;
            case 'a':  //alpha ranking function 偏好参数：
                alpha = getAlphaSettings(optarg);
                flag = ch;

                break;
            case 's':  //solution: 代表使用的处理算法
                alg = getAlgChose(optarg);
                ss1.clear();
                break;
            case 'p': //popularity: 代表查询对象所带关键词的热度
                ss1<< optarg;
                popularity = getPopularityChose(optarg);
                flag = ch;
                ss1.clear();
                break;

            case 't': //选择top-k处理方法
                flag = ch;
                topk_method = getTopKMethodChose(optarg);
                ss1.clear();
                break;

            case 'r': //选择reverse top-k处理方法
                flag = ch;
                reverse_method = getReverseMethod(optarg);
                ss1.clear();
                break;

            case 'v': //选择reverse top-k处理方法
                flag = ch;
                accuracy = getAccuracySettings(optarg);
                ss1.clear();
                break;
            case 'j': //选择reverse top-k处理方法
                flag = ch;
                _task = getPreprocessTask(optarg);
                ss1.clear();
                break;

            default :
                break;

        }
    }
}


void performanceParameter_settings(int argc, char* argv[], int& k, int& KQ, int& QN, int& b, float& accuracy, float& alpha, float& o_ratio, int& popularity, int& alg,int& exp_task, int& test_target, int& group_type,string& _task){
    stringstream ss1;
    int ch;   char flag; int tmp = 0;

    //根据控制台参数进行相应设定
    while ((ch = getopt(argc, argv, "b:c:k:n:w:a:s:p:e:g:o:t:r:v:j:"))!=-1){
        switch(ch) {
            case 'e': //experimant: 选择任务
                flag = ch;
                exp_task = getExpChose(optarg);
                ss1.clear();
                break;
            case 'g': //选择reverse top-k处理方法

                group_type = getQueryGroupChose(optarg);
                flag = ch;
                break;

            case 't': //选择top-k处理方法
                flag = ch;
                test_target = getTargetParameterSetting(optarg);
                ss1.clear();
                break;

            case 'k':   //top-k参数
                k = getOkSettings(optarg);
                flag = ch;
                break;
            case 'n':  //number of query 代表查询对象个数
                QN = getOcSettings(optarg) ;//getQNSettings(optarg);
                flag = ch;
                break;
            case 'w':  //word count:代表查询关键词个数
                KQ = getKQSettings(optarg);
                flag = ch;
                break;
            case 'a':  //alpha ranking function 偏好参数：
                alpha = getAlphaSettings(optarg);
                flag = ch;
                break;

            case 'b':  //alpha ranking function 偏好参数：
                b = getBudgetSettings(optarg);
                flag = ch;
                break;
            case 's':  //solution: 代表使用的处理算法
                alg = getAlgChose(optarg);
                ss1.clear();
                break;
            case 'p': //popularity: 代表查询对象所带关键词的热度
                ss1<< optarg;
                popularity = getPopularityChose(optarg);
                flag = ch;
                ss1.clear();
                break;

            case 'o': //选择reverse top-k处理方法
                o_ratio = getObjectSizeRatioSettings(optarg);
                ss1.clear();
                break;


            case 'v': //选择reverse top-k处理方法
                flag = ch;
                accuracy = getAccuracySettings(optarg);
                ss1.clear();
                break;
            case 'j': //选择reverse top-k处理方法
                flag = ch;
                _task = getPreprocessTask(optarg);
                ss1.clear();
                break;

            default :
                break;

        }
    }
}





void runAllExperiment(string exp, int Qk, int alg, int KQ, int QC, float alpha, char flag){

    if(exp == "effectiveness"){
        test_MaxInfBRGSTkNN_effectiveness(Qk, KQ, QC, alpha, alg, flag);
    }
    else{
        //性能评估
        test_MaxInfBRGSTkNN_performance(Qk, KQ, QC, alpha, alg, flag);
    }

    //getchar();

}





//below this for test

#define  MAXINFRKGSKQ_PORCESS
#endif