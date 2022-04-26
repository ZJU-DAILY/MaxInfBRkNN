#ifndef MC_H
#define MC_H

#include "common.h"
#include "HashTreeMap.cc"
#include <vector>
#include <cmath>
#include <ctime>
#include <map>
#include <queue>
#include "../../results.h"
#include "../SocialSimulator.h"
#include "anyoption.h"

namespace _MC {

typedef HashTreeMap<UID, FriendsMap*> HashTreeCube;
typedef HashTreeMap<UID, multimap<float, UID>*> HashTreeCube2;
typedef map<int, unsigned int> RadiusCovMap;

struct MGStruct {
	UID nodeID; // user ID
	float gain; // MG(u|S)
	UID v_best; // user that gives best MG till now
	
	// marginal gain of user u w.r.t the curSeedSet + u_best. i.e., MG(u|S+v_best)
	float gain_next;
	int flag; // size of the seed set in current run
};

// to make a list sorted by MGs (or coverage)
typedef multimap<float, MGStruct*> Gains;  //按float由大到小排序

struct NodeParams {
	float threshold;
	float inWeight;
	bool active;
};



class MC : public SocialSimulator{

public:
	AnyOption* opt;
	
	UserList curSeedSet;
    POIList curPOISet;
    UserList users;
	UserList covBestNode;	 // newly made for binary prob. case 
	//UserList covSeedSet;
    unsigned int numEdges;

	HashTreeCube *AM = NULL;
    //HashTreeCube* AM_in;
    HashTreeCube2 *revAM;
    FriendsMap seedSetNeighbors;
	Gains mgs;

	string outdir;
	const char* probGraphFile;
	string m;
	PropModels model; 
	GraphType graphType;

	string problem;
	ofstream outFile;
	int startIt;
	int binaryProb; // 1: yes, binary probability. 0: no, real probability.


	int countIterations;
	float totalCov; 

	// parameters to monitor the progress of experiments
	time_t startTime;
	time_t stime_mintime;

	MC (AnyOption* opt);
	~MC ();

	//for MaxInfBRGSTQ processing

    float mineSeedSetPlus();

    float minePOI_CLEFPlus(vector<int>& stores, BatchResults& results, int b);

	void readInputData(float alpha=0);
	PropModels getModel();

	float LTCov(UserList& list);
	float computeLTCov(UserList& list);
	float ICCov(UserList& list);
	

	void setAM(HashTreeCube* AM1) {AM = AM1;}
    HashTreeCube *getAM() {return this->AM;}
    //HashTreeCube *getAMIn() {return this->AM_in;}
    HashTreeCube2 *getRevAM() {return this->revAM;}
    UserList *getUsers() {return &this->users;}
    unsigned int getNumEdges() {return this->numEdges;}

	void addSocialLink(int u, int v, double prob);
    int getAMSize(){
		return AM->size();
    }

    double MCSimulation4Seeds(int r);

    void printPOISelections(){
    	for(int poi: curPOISet){
    		cout<<"p"<<poi<<endl;
    	}
    }

	void printRealSeeds(){
		cout<<"{ " ;
		for(int u: curSeedSet){
			cout<<"u"<<u<<",";
		}
		cout<<"};"<<endl;
	}


private:
	// functions called from constructor
	void setModel(AnyOption* opt);

	// functions called from doAll
	float mineSeedSet(int t_ub=0);
	void mintime();
	void genMintimeTable();
	void clear();
	void computeCov();

	/****** For CELF++ on Monte Carlo Simulation!  *****/
	bool ICCovPlus(MGStruct *pMG, MGStruct *pBestMG); // IC model + CELF++
    bool LTCovPlus(MGStruct *pMG, MGStruct *pBestMG, UserList &); // LT model + CELF++	

    /****** For CELF++ on Monte Carlo Simulation to evaluate the influence of a POI!  *****/
    bool ICPOICovPlus(MGStruct* pMG,  BatchResults& results, MGStruct* pBestMG);



	// other private functions
	float getTime() const;
	float getTime_cur() const;
	void writeInFile(UID v, float cov, float marginal_gain, int curTimeStep, float actualCov, float actualMG, int countUsers);

    void writePOIChooseToFile(int p, float cov, float marginal_gain, int curTimeStep, float actualCov, float actualMG, int countUsers);



    void openOutputFiles();
	void openExpResutlsFiles();
	void printVector(vector<UID>& vec, float pp);



};
}
#endif
