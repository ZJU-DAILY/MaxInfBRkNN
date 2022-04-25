#ifndef __NET_SHARED
#define __NET_SHARED

#include "utility.h"
#include <queue>

#define MAX_DIST 9999999
#define REALTOTALKEY 64


typedef FastList<int> IDset;

int clusNum;
//float** NodeDist;

void printCluster(IDset &clus) {
	//printf("Size: %d\n",clus.size());
	/*sort(clus.begin(),clus.end());	// sorted if print
	printf("{ ");
	//for (int r=0;r<clus.size();r++) printf("%d ",clus[r]);
	int lb=0,s=0,r=0;

	int checkSum=0;
	//printf("\n");
	while (r<clus.size()) {
		if (clus[r]==clus[lb]+s) {	// adjacent
			s++;
		} else {
			checkSum+=s;
			printf("[%d,%d] ",clus[lb],clus[lb]+s-1);
			lb=r;
			s=1;
		}
		r++;
	}
	checkSum+=s;
	printf("[%d,%d] ",clus[lb],clus[lb]+s-1);

	printf("} with size %d\n",checkSum);*/
}

//-------
// Step
//-------
struct StepEvent {
	float dist;
	int node,ClusID;	// ClusID for multiple expansion
	float accDist;		// posDist for gendata.cc only
	int Ni,Nj;	// for gendata use only
	bool isLeft; // for gendata use only
};

struct StepComparison {
	bool operator () (StepEvent left,StepEvent right) const
    { return left.dist > right.dist; }
};

typedef	priority_queue<StepEvent,vector<StepEvent>,StepComparison> StepQueue;

void printStepEvent(StepEvent& event)  {
	printf("(%d,%0.3f)\n",event.node,event.dist);
}

#endif // __NET_SHARED


