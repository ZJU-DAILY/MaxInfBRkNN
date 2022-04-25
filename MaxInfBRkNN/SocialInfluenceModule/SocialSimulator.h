//
// Created by jin on 20-2-26.
//

#ifndef MAXINFBRGSTQ_SOCIALSIMULATOR_H

enum Influce_Simulator{CELF, PMC, IMM, OPC, KIM};

#define  SEED1


class SocialSimulator{
public: int simulator_type;
        void setSimulatorType(int tmp){
            simulator_type = (tmp);
        }
};

#define MAXINFBRGSTQ_SOCIALSIMULATOR_H

#endif //MAXINFBRGSTQ_SOCIALSIMULATOR_H
