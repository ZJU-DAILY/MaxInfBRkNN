//
// Created by jin on 20-2-26.
//

#ifndef MAXINFBRGSTQ_ARGUMENT_H



#define MAXINFBRGSTQ_ARGUMENT_H

#include <string>
using namespace std;

class Argument{
public:
    int k;
    string dataset;
    float epsilon=0.1;
    string model;
    double T;
};

#endif //MAXINFBRGSTQ_ARGUMENT_H
