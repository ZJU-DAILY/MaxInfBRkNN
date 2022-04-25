
#ifndef _TRANSFORMINPUTCOMMAND_H
#define _TRANSFORMINPUTCOMMAND_H

#include "Command.h"
#include "../common.h"

#include <unordered_set>
#include <cctype>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <limits>
#include <regex>
#include <unordered_map>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <vector>
#include <iostream>

#include "../processing/MortonList.h"
#include "../queue/BinaryMinHeap.h"
#include "../utility/StopWatch.h"
#include "../utility/utility.h"
#include "../utility/geometry.h"




using namespace std;

class TransformInputCommand: public Command {

    public:

        void execute(int argc, char* argv[]){

        }


        void cleanRoadNetworkSourceData(){
            //对应-m, -n 下参数

            std::string dataset = "synthetic";

            #ifdef LV_Source
                dataset = "LV";
            #endif

            #ifdef LasVegas_Source
                dataset = "LasVegas";
            #endif

            #ifdef Brightkite_Source
                dataset = "BRI";
            #endif

            #ifdef Gowalla_Source
                dataset = "GOW";
            #endif

            #ifdef Twitter_Source
                dataset = "EUS";
            #endif

            //对应-o下参数
            stringstream transformPath;
            transformPath<<DATASET_HOME<<dataset_name<<"/road/";
            std::string outputFilePrefix = transformPath.str();
            std::string regionName = dataset_name;
            std::string subRegionName = "";

            //对应-e, -n 下参数
            stringstream tt; stringstream tt2;
            tt<<FILE_SOURCE<<"road/"<<"road-d."<<dataset_name<<".co";
            std::string coFilePath = tt.str();
            tt2<<FILE_SOURCE<<"road/"<<"road-d."<<dataset_name<<".eg";
            std::string grFilePath = tt2.str();

            //std::cout<<coFilePath<<std::endl;
            //std::cout<<grFilePath<<std::endl;

            int edgeWeightInflateFactor = 0, coordinateInflateFactor = 0;


            if (dataset == "LasVegas" || dataset == "BRI" || dataset == "LV" ) {  //|| dataset == "LV"

                StopWatch sw;

                sw.start();
                //this->transformOSMFiles(grFilePath,coFilePath,outputFilePrefix);
                this->cleanRawRoadNetworkFiles(grFilePath,coFilePath,outputFilePrefix);
                sw.stop();

                double processingTimeMs = sw.getTimeMs();

                std::cout<<endl<<dataset_name<< "Dataset source processe SUCCESS!, runtime=" << processingTimeMs << "ms" << std::endl;

            }
            else if (dataset == "GOW"|| dataset =="EUS") {

                StopWatch sw;

                sw.start();
                this->cleanRawRoadNetworkFiles(grFilePath,coFilePath,outputFilePrefix);  //这里是big的

                sw.stop();

                double processingTimeMs = sw.getTimeMs();

                std::cout << "\nSynthetic Dataset source("<< dataset_name<<") processed in " << processingTimeMs << "ms" << std::endl;
            }


            else {
                std::cerr << "Invalid query method!\n\n";
                exit(1);
            }


    }


        void showCommandUsage(std::string programName);
        void showMethodUsage(std::string method, std::string programName);
        
        // Input Types for Transform Command
        std::string const TRANSFORM_DIMACS = "dimacs";
        std::string const TRANSFORM_OSM = "osm";
        std::string const TRANSFORM_TPQ = "tpq";
        std::string const TRANSFORM_OSM_DIMACS = "osm_dimacs";
    
    private:
        struct Point {
            int x;
            int y;
            Point(int _x, int _y): x(_x), y(_y) {}
        };
        struct Edge {
            NodeID target;
            EdgeWeight weight;
            Edge(NodeID _target, EdgeWeight _weight): target(_target), weight(_weight) {}
        };

    void outputStandardFormatFiles(std::string outputFilePrefix){
        std::string newCoFilePath = outputFilePrefix + dataset_name + "." + constants::NODE_EXT;

        std::ofstream outputCoordFile(newCoFilePath, std::ios::out | std::ios::trunc);
        if (outputCoordFile.is_open()) {
            double scaled_times = scale_coordinate;
            //outputCoordFile << this->networkName << " " << this->edgeType << " " << this->numNodes << "\n";
            for (std::size_t i = 0; i < this->nodes.size(); ++i) {  //注意这里的坐标数值是放大了10^6倍的
                outputCoordFile << this->nodes[i]<< " " << this->coordinates[i].first/scaled_times << " " << this->coordinates[i].second/scaled_times << "\n";
            }
        } else {
            std::cerr << "Cannot open coordinate output file " << newCoFilePath << std::endl;
        }

        std::string newGrFilePath = outputFilePrefix + "/" + dataset_name + "." + constants::EDGE_EXT;

        std::ofstream outputGrFile(newGrFilePath, std::ios::out | std::ios::trunc);
        if (outputGrFile.is_open()) {
            //outputGrFile << this->networkName << " " << this->edgeType << " " << this->numNodes << " " << this->numEdges << "\n";
            int count = 0;
            for (std::size_t i = 0; i < this->nodes.size(); ++i) {
                for (std::size_t j = 0; j < this->neighbours[i].size(); ++j) {
                    //float dist_original = this->neighbours[i][j].second*1.0/scale_edge;
                    //float dist =  scaled_dist *1.0 / scale_edge;
                    outputGrFile << count<<"\t"<<this->nodes[i] << "\t" << this->neighbours[i][j].first << "\t" << this->neighbours[i][j].second << "\n";
                    count++;
                }
            }
        } else {
            std::cerr << "Cannot open graph output file " << newGrFilePath << std::endl;
        }
    }



    void cleanRawRoadNetworkFiles(std::string grFilePath, std::string coFilePath, std::string outputFilePrefix){  //路网数据为NA shared by feifeili in 2005 sstd
        std::cout << "\nProcessing Coordinates File for " << dataset_name << std::endl;

        std::string line = "", lineType = "", dummyField;
        NodeID sourceID, targetID;
        EdgeWeight weight;
        int x, y, rawWeight; set<int> vertex_set;
        std::string coordEdgeType;
        std::ifstream coordFile(coFilePath, std::ios::in);

        if (coordFile.is_open()) {
            // Note: Boost regular expression has to match whole input string
            // partial matches will not occur
            std::regex expression(".*\\-([dt])\\.([^.]+)\\.co$");
            std::smatch matches;
            if (std::regex_match(coFilePath,matches,expression)) {
                coordEdgeType = matches[1];
                this->networkName = matches[2];
            } else {
                std::cout << "Input coordinate file not named in expected OSM format!" << std::endl;
                exit(1);
            }


            int vertex_id=0; float x_value,y_value;

            while (std::getline(coordFile, line))
            {
                //std::stringstream ss(line);
                istringstream tt(line);
                tt >> vertex_id >> x_value >> y_value ;
                Coordinate x = 0;
                Coordinate y = 0;
                x =   x_value*scale_coordinate ;
                y =   y_value*scale_coordinate ;
                // id should be in order
                this->addNode(vertex_id,x,y);
                vertex_set.insert(vertex_id);
            }
            this->numNodes = vertex_set.size();
            if (this->numNodes == 0) {
                std::cerr << "File header information not found or does not contain total number of nodes" << std::endl;
                exit(1);
            }

            if (this->numNodes != this->nodes.size() || this->numNodes != this->coordinates.size()) {
                std::cerr << "Number of nodes read from file does not match number of nodes given in header" << std::endl;
                exit(1);
            }

        } else {
            std::cerr << "Cannot open coordinates file " << coFilePath << std::endl;
            exit(1);
        }

        unsigned int edgesFileNumNodes, edgesAdded = 0;
        this->neighbours.resize(this->numNodes);

        map<NodeID,map<NodeID,bool>> edge_accessed;

        std::ifstream grFile(grFilePath, std::ios::in);
        if (grFile.is_open()) {
            std::regex expression(".*\\-([dt])\\.([^.]+)\\.gr$");
            std::smatch matches;

            int edge_id; int Ni; int Nj; double dist; int count = 0;

            vector<vector<int>> adjList; vector<vector<int>> adjDist_list;
            adjList.resize(this->numNodes); adjDist_list.resize(this->numNodes);
            while (std::getline(grFile, line))
            {
                istringstream tt(line);
#ifdef Twitter_Source
                tt >> Ni >> Nj >> dist ; edge_id = count; count++;
#else
                tt >> edge_id >> Ni >> Nj >> dist ;
#endif

                if(Ni==Nj) continue;

                // Validate Edges
               /* if(Ni<2){
                    cout<<"Ni="<< Ni<<",Nj="<<Nj<<",dist="<<dist<<endl;
                    getchar();
                }*/

                ////去除重复边（两顶点间有多重边：Parallel edge detected）
                if(edge_accessed[Ni].count(Nj)){
                    ///cout<<"有重复边！"<<endl;
                    continue;
                }
                edge_accessed[Ni][Nj] = true;
                edge_accessed[Nj][Ni] = true;

                double euclidDist_cmp =
                        geometry::getEuclideanDist_BRI(this->coordinates[Ni].first,this->coordinates[Ni].second,
                                                                  this->coordinates[Nj].first,this->coordinates[Nj].second) / scale_coordinate;


                double raw_weight = max(dist, euclidDist_cmp);

#ifdef LasVegas_Source
                raw_weight = dist;
#endif

#ifdef LV_Source
                raw_weight = dist;
#endif


#ifdef Brightkite_Source
                raw_weight = 100* dist;  //max(100*dist, euclidDist_cmp);
#endif

#ifdef Gowalla_Source
                raw_weight = 100* dist;  //max(100*dist, euclidDist_cmp);
#endif

#ifdef Twitter_Source
                raw_weight = 100* dist;  //max(100*dist, euclidDist_cmp);
#endif



                //int scaled_weight = (int)(raw_weight * scale_edge);

                int edge_dist = ceil(raw_weight);
                adjList[Ni].push_back(Nj);
                adjList[Nj].push_back(Ni);
                adjDist_list[Ni].push_back(edge_dist);
                adjDist_list[Nj].push_back(edge_dist);


            }
            //将完整的（无向）边信息进行插入
            for(int i=0; i<adjList.size();i++){
                if(i%10000==0)
                    cout<<"v"<<i<<" has "<<adjList[i].size()<<" neighbors"<<endl;
                int Ni = i;
                for(int j=0;j<adjList[Ni].size();j++){
                    int Nj = adjList[Ni][j];
                    int dist = adjDist_list[Ni][j];
                    this->addEdge(Ni,Nj,dist);
                    this->totalEdgeWeight += weight;
                    this->numEdges++;
                    ++edgesAdded;
                }

            }


            if (this->numEdges != edgesAdded) {
                std::cerr << "Number of edges read from file does not match number of edges given in header" << std::endl;
                exit(1);
            }

        } else {
            std::cerr << "Cannot open graph file " << grFilePath << std::endl;
            exit(1);
        }


        //this->fixShortestPathOverflow();

        this->computeCoordinateRanges();  //计算整副路网图的MBR

#ifdef LasVegas_Source
        this->checkNonUniqueCoordinates();
        //this->removeNonUniqueNodes();
        this->removeIsolatedRegions();
        std::cout<<"完成 removeIsolatedRegions"<<endl;

        this->removeIsolatedNodes();
        std::cout<<"完成 removeIsolatedNodes"<<endl;

#endif

#ifdef Twitter_Source
        this->checkNonUniqueCoordinates();
        //this->removeNonUniqueNodes();
        this->removeIsolatedRegions();
        std::cout<<"完成 removeIsolatedRegions"<<endl;

        this->removeIsolatedNodes();
        std::cout<<"完成 removeIsolatedNodes"<<endl;

#endif
        this->outputStandardFormatFiles(outputFilePrefix);
    }



        void addNode(NodeID node, Coordinate x, Coordinate y);
        void addEdge(NodeID source, NodeID target, EdgeWeight weight);
        bool isValidEdge(NodeID source, NodeID target, int rawWeight, EdgeWeight weight);
        bool checkValidNode(NodeID node);
        EdgeWeight correctEdgeWeightByEuclidDist(NodeID source, NodeID target, EdgeWeight weight);
        void fixShortestPathOverflow();
        LongPathDistance computePseudoDiameter(int numRounds = 1);
        void computeCoordinateRanges();
        void computeCoordinateTranslations(std::vector<Point>& coordinates, int& minX, int& minY);
        void checkNonUniqueCoordinates();

        void outputStandardFormatFiles_old(std::string outputFilePrefix){  //原来tenda的方式
            std::string newCoFilePath = outputFilePrefix + "/" + this->networkName + "-" + this->edgeType + "." + constants::NODE_EXT;

            std::ofstream outputCoordFile(newCoFilePath, std::ios::out | std::ios::trunc);
            if (outputCoordFile.is_open()) {
                outputCoordFile << this->networkName << " " << this->edgeType << " " << this->numNodes << "\n";
                for (std::size_t i = 0; i < this->nodes.size(); ++i) {
                    outputCoordFile << this->nodes[i] << " " << this->coordinates[i].first << " " << this->coordinates[i].second << "\n";
                }
            } else {
                std::cerr << "Cannot open coordinate output file " << newCoFilePath << std::endl;
            }

            std::string newGrFilePath = outputFilePrefix + "/" + networkName + "-" + edgeType + "." + constants::EDGE_EXT;

            std::ofstream outputGrFile(newGrFilePath, std::ios::out | std::ios::trunc);
            if (outputGrFile.is_open()) {
                outputGrFile << this->networkName << " " << this->edgeType << " " << this->numNodes << " " << this->numEdges << "\n";
                for (std::size_t i = 0; i < this->nodes.size(); ++i) {
                    for (std::size_t j = 0; j < this->neighbours[i].size(); ++j) {
                        int scaled_dist = this->neighbours[i][j].second;
                        //float dist =  scaled_dist *1.0 / scale_edge;
                        outputGrFile << this->nodes[i] << " " << this->neighbours[i][j].first << " " << scaled_dist << "\n";
                    }
                }
            } else {
                std::cerr << "Cannot open graph output file " << newGrFilePath << std::endl;
            }
        }



        void printTransformationStatistics();
        void removeNonUniqueNodes();
        void removeIsolatedNodes();
        void removeIsolatedRegions();
        
        // Graph Information
        unsigned int numNodes = 0, numEdges = 0;
        std::string networkName = "", edgeType = "";
        
        // Error Detected/Corrected Statistics
        int numSelfLoops = 0, numParallelEdges = 0, numNegativeEdgeWeights = 0, numEuclideanCorrections = 0, numCorrectionsOverFactorTen = 0;
        int numNonUniqueNodesRemoved = 0, numDisconnectedNodesRemoved = 0, numNonUniqueNodeEdgesRemoved = 0, numDisconnectedNodeEdgesRemoved = 0;
        int numRegionsFound = 0, numIsolatedRegionNodesRemoved = 0, numIsolatedRegionEdgesRemoved = 0;
        unsigned long long totalEdgeWeight = 0, weightAddedByEuclideanCorrections = 0;
        int xTranslation = 0, yTranslation = 0;
        bool wasDeflated = false;
        int deflationFactor = 0;
        unsigned long long oldTotalEdgeWeight = 0;

        // Graph
        std::vector<NodeID> nodes;
        std::vector<std::vector<NodeEdgeWeightPair>> neighbours;
        std::vector<CoordinatePair> coordinates;
        std::unordered_set<NodeID> nonUniqueNodes;
        
};

#endif // _TRANSFORMINPUTCOMMAND_H
