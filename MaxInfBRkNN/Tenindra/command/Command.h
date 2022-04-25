

#ifndef _COMMAND_H
#define _COMMAND_H

#include <string>
#include <unistd.h>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include "../../source_config.h"

class Command {

    public:
        virtual void execute(int argc, char* argv[]) = 0;
        virtual void showCommandUsage(std::string programName) = 0;
        virtual ~Command() {};
        
    protected:
        void outputCommandStats(std::string statsFilePath, std::string tupleString) {
            std::ofstream statsFS(statsFilePath, std::ios::out | std::ios::app);
            if (statsFS.is_open()) {
                statsFS << tupleString << std::endl;
            } else {
                std::cerr << "Cannot open stats file " << statsFilePath << std::endl;
            }
        }
};

#endif // _COMMAND_H