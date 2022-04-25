/* Copyright (C) 2015 Tenindra Abeywickrama
 *
 * This file is part of Road Network kNN Experimental Evaluation.
 *
 * Road Network kNN Experimental Evaluation is free software; you can
 * redistribute it and/or modify it under the terms of the GNU Affero 
 * General Public License as published by the Free Software Foundation; 
 * either version 3 of the License, or (at your option) any later version.
 *
 * Road Network kNN Experimental Evaluation is distributed in the hope 
 * that it will be useful, but WITHOUT ANY WARRANTY; without even the 
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE. See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public 
 * License along with Road Network kNN Experimental Evaluation; see 
 * LICENSE.txt; if not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _GRAPHCOMMAND_H
#define _GRAPHCOMMAND_H

#include "Command.h"

#include "../processing/Graph.h"
#include "../tuple/IndexTuple.h"
#include "../common.h"
#include "../utility/StopWatch.h"
#include "../utility/serialization.h"
#include "../processing/Junction.h"

#include <fstream>
#include <boost/archive/binary_oarchive.hpp>


class GraphCommand: public Command {

    public:
        //void execute(int argc, char* argv[]);

        //void showCommandUsage(std::string programName);

        void execute(int argc, char* argv[]) {
        std::string coFilePath = "";
        std::string grFilePath = "";
        std::string indexOutputFilePath = "";
        std::string statsOutputFilePath = "";

        /*
         * Process Command Line Arguments
         */
        /*if (argc < 9) {
            // Arguments: -e <graph file> -n <coordinates file> -o <output file> -s <stats file>
            std::cerr << "Too few arguments!\n\n";
            this->showCommandUsage(argv[0]);
            return;
        }*/

        int opt;
        while ((opt = getopt (argc, argv, "e:n:o:s:")) != -1) {
            switch (opt) {
                case 'e':
                    grFilePath = optarg;
                    break;
                case 'n':
                    coFilePath = optarg;
                    break;
                case 'o':
                    indexOutputFilePath = optarg;
                    break;
                case 's':
                    statsOutputFilePath = optarg;
                    break;
                default:
                    std::cerr << "Unknown option(s) provided!\n\n";
                    showCommandUsage(argv[0]);
                    exit(1);
            }
        }

        if (grFilePath == "" || indexOutputFilePath == ""
            || coFilePath == "" || statsOutputFilePath == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showCommandUsage(argv[0]);
            exit(1);
        }

        /*
         * Parse Graph from Text File and Serialize to Binary File
         */
        Graph_RoadNetwork graph;
        StopWatch sw;

        sw.start();
        graph.parseGraphFile(grFilePath,coFilePath);
        sw.stop();
        double processingTimeMs = sw.getTimeMs();

        // Print Additional Graph Statistics
        graph.printGraphStats();

        /*
         * Serialize to Binary File
         */
        serialization::outputIndexToBinaryFile<Graph_RoadNetwork>(graph,indexOutputFilePath);

        /*
         * Collect Stats and Output
         */
        int nodes = graph.getNumNodes();
        int edges = graph.getNumEdges();
//  double memoryUsage = graph.computeIndexSize();
        double memoryUsage = graph.computeINEIndexSize(); // For experiments we report INE size

        IndexTuple stats(graph.getNetworkName(),nodes,edges,constants::IDX_GRAPH_CMD,processingTimeMs,memoryUsage);

        std::cout << stats.getMultilineTupleString();

        this->outputCommandStats(statsOutputFilePath,stats.getTupleString());

        std::cout << "Binary graph index successfully created!" << std::endl;

        /*
         * Supplementary Junction Index
         */

        sw.reset();
        sw.start();
        Junction junc(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges());
        junc.buildJunction(graph);
        sw.stop();

        // Print Additional Graph Statistics
        //junc.printGraphStats();

        /*
         * Serialize to Binary File
         */
        int extIndex = indexOutputFilePath.find_last_of(".");
        std::string idxFileName = indexOutputFilePath.substr(0,extIndex)+"."+constants::IDX_JUNC_CMD;
        serialization::outputIndexToBinaryFile<Junction>(junc,idxFileName);

        /*
         * Collect Stats and Output
         */
        processingTimeMs = sw.getTimeMs();
        memoryUsage = junc.computeIndexSize();

        IndexTuple stats2(junc.getNetworkName(),junc.getNumNodes(),junc.getNumEdges(),constants::IDX_JUNC_CMD,processingTimeMs,memoryUsage);

        //this->outputCommandStats(statsOutputFilePath,stats2.getTupleString()); // No need to output junction stats

        std::cout << "Junction index successfully created!" << std::endl;

        std::cout << stats2.getMultilineTupleString();
    }


        void buildBinaryGraph() {


            std::string method = "synthetic";
#ifdef LV
            method = "real";
#endif

#ifdef LasVegas
            method = "real";
#endif
            //对应-o下参数
            stringstream binaryOutput;
            binaryOutput<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name<<".bin";
            std::string outputFilePath = binaryOutput.str();
            std::string regionName = dataset_name;
            std::string subRegionName = "";

            //对应-e, -n 下参数
            stringstream nodePath; stringstream edgePath;
            nodePath<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name<<".cnode";
            std::string coFilePath = nodePath.str();
            edgePath<<DATASET_HOME<<dataset_name<<"/road/"<<dataset_name<<".cedge";
            std::string grFilePath = edgePath.str();

            //对应-s 下参数
            std::string statsOutputFilePath = "";
            stringstream statsPath;
            statsPath<<DATASET_HOME<<dataset_name<<"/road/"<<"stats/index_stats.txt";
            statsOutputFilePath = statsPath.str();

            /*
            * Parse Graph from Text File and Serialize to Binary File
            */
            Graph_RoadNetwork graph;
            StopWatch sw;

            sw.start();
            cout<<"加载路网图数据！"<<endl;
            graph.parseGraphFile(grFilePath,coFilePath);
            sw.stop();
            double processingTimeMs = sw.getTimeMs();
            // Print Additional Graph Statistics
            graph.printGraphStats();

            /*
            * Serialize to Binary File
            */
            serialization::outputIndexToBinaryFile<Graph_RoadNetwork>(graph,outputFilePath);

            /*
            * Collect Stats and Output
            */
            int nodes = graph.getNumNodes();
            int edges = graph.getNumEdges();
        //  double memoryUsage = graph.computeIndexSize();
            double memoryUsage = graph.computeINEIndexSize(); // For experiments we report INE size

            IndexTuple stats(dataset_name,nodes,edges,constants::IDX_GRAPH_CMD,processingTimeMs,memoryUsage);

            std::cout << stats.getMultilineTupleString();

            this->outputCommandStats(statsOutputFilePath,stats.getTupleString());

            std::cout << "Binary graph index successfully created!" << std::endl;


        }

        void showCommandUsage(std::string programName) {
        std::cerr << "Usage: " << programName << " -c " + constants::IDX_GRAPH_CMD + " -e <text graph file>"
                  << " -n <coordinates file>\n-o <binary graph output file> -s <stats file>\n";
        }

};

#endif // _GRAPHCOMMAND_H