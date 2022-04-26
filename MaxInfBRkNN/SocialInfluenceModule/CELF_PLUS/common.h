#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <utility>
#include <set>
#include <stack>
#include <sys/resource.h>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>
#include <climits>
#include "anyoption.h"
#include <stdint.h>
#include "unistd.h"
#include "../imbasic.h"

using namespace std;
using std::ifstream;

typedef int64_t int_64;

typedef int UID; // userid
typedef uint32_t AID; // actionid 
typedef uint32_t TS; // timestamp
#define TSMAX UINT_MAX

typedef set<UID> UserList;
typedef set<int> POIList;
typedef map<UID, int> UserCounts;
typedef map<UID, int> UsersCounts;
typedef map<UID, UserList> Matrix; // matrix
typedef map<AID, Matrix> Cube; // Cube of <action, user, list_of_users>
typedef map<AID, unsigned int> InfCounts;
typedef map<UID, float> FriendsMap;



const int MAXLINE=256;
const int wordsize = 32;
const unsigned int AVGDELAY = 7862400; // 13 weeks






enum GraphType {
	DIRECTED, UNDIRECTED
};

inline unsigned int strToInt(string s) {
	unsigned int i;
	istringstream myStream(s);

	if (myStream>>i) {
		return i;
	} else {
		cout << "String " << s << " is not a number." << endl;
		return atoi(s.c_str());
	}
	return i;
}

inline unsigned int strToInt(const char* s) {
	unsigned int i;
	istringstream myStream(s);

	if (myStream>>i) {
		return i;
	} else {
		cout << "String " << s << " is not a number." << endl;
		return atoi(s);
	}
	return i;
}


inline float strToFloat(const char* s) {
	return atof(s);
	float i;
	istringstream myStream(s);

	if (myStream>>i) {
		return i;
	} else {
		cout << "String " << s << " is not a float." << endl;
		return atof(s);
	}
	return i;
}

inline float strToFloat(string s) {
	return atof(s.c_str());
	float i;
	istringstream myStream(s);

	if (myStream>>i) {
		return i;
	} else {
		cout << "String " << s << " is not a float." << endl;
		return atof(s.c_str());
	}
	return i;
}

inline double strToDouble(string s) {
    return atof(s.c_str());
    double i;
    istringstream myStream(s);

    if (myStream>>i) {
        return i;
    } else {
        cout << "String " << s << " is not a float." << endl;
		return atof(s.c_str());
    }
    return i;
}

inline string floatToStr(float f) {
	stringstream ss;
	ss << f;
	return ss.str();
}

inline int_64 strToInt64(string s) {
	int_64 i;
	istringstream myStream(s);

	if (myStream>>i) {
		return i;
	} else {
		cout << "String " << s << " is not a number." << endl;
		exit(1);
	}
	return i;
}

inline string intToStr(int i) {
	stringstream ss;
	ss << i;
	return ss.str();
}



float getCurrentMemoryUsage();
#endif
