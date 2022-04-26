#ifndef NODE_H 
#define NODE_H

#include "common.h"
#include <string>

class Node {
public:
	
	int out;
	int in; // number of actions v followed u
	float out_prob; // number of actions u followed v
	float in_prob;
	float out_time;
	float in_time;

	Node() : out(0), in(0), out_prob(0), in_prob(0), out_time(0), in_time(0) {};
	~Node() {};
};

#endif
