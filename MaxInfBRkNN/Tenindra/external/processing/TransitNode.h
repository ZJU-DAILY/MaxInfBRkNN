/* Original Work Copyright (C) 2012 Lingkun Wu, Xiaokui Xiao, Dingxiong Deng, Gao Cong, Andy Diwen Zhu, Shuigeng Zhou
 * Modified Work Copyright (C) 2015 Tenindra Abeywickrama
 *
 * This file is part of Road Network kNN Experimental Evaluation.
 * 
 * The following file is a derivative work of the code from 
 * http://sourceforge.net/projects/ntu-sp-exp/ developed for the paper below.
 * The authors have requested users of this code to cite the paper.
 * 
 * Lingkun Wu, Xiaokui Xiao, Dingxiong Deng, Gao Cong, Andy Diwen Zhu, Shuigeng Zhou
 * Shortest Path and Distance Queries on Road Networks: An Experimental Evaluation
 * PVLDB 5(5): 406-417 (2012)
 *
 * That work is in turn a derivative work of the code from
 * Contraction Hierarchies (http://algo2.iti.kit.edu/english/routeplanning.php),
 * licensed under AGPLv3. Thus this work is also licensed under that license.
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

#ifndef TRANSIT_H
#define TRANSIT_H

#include "../datastr/graph/graph.h"

#include "../datastr/graph/UpdateableGraph.h"
#include "DijkstraCH.h"
#include "../many/manytomany.h"

#include <unordered_map>

#include "../../queue/BinaryMaxHeap.h"
#include "../../queue/BinaryMinHeap.h"
#include "../../processing/StaticRtree.h"

//#define  DEBUG

namespace processing{

	typedef datastr::graph::SearchGraph MyGraph;
	typedef DijkstraCHManyToManyFW DijkstraManyToManyFW;
	typedef DijkstraCHManyToManyBW DijkstraManyToManyBW;
	

	template <typename Graph>
	class TransitNode {
	public:
		
		const static NodeID MAX_CROSS_FOR_EDGE = 1000;
		const static NodeID MAX_TRANSITNODE_FOR_CELL = 20000;
		const static NodeID EACH_CELL_TRANSIT_NUM = 30;

		////for statistics
		//int m_iLevelOne;
		//int m_iLocal;

		TransitNode(Graph* _graph, NodeID _ascale, CoordinateType _range): graph(_graph), ascale(_ascale), range(_range)
		{
			cellLength = range / ascale;
			NodeID totalNodes = graph->noOfNodes();
			avgNodesPerCell = totalNodes / (ascale * ascale);
			allCellNum = ascale * ascale;
			// allocate each cell contained node allocate space for each cell Transit node
			//eachCellContainedNode = (NodeID **) malloc(sizeof(NodeID *) * allCellNum);
			eachCellTransitNode = (NodeID **) malloc(sizeof(NodeID *) * allCellNum);
			for (NodeID i = 0; i < allCellNum ; i++){
				/*eachCellContainedNode[i] = (NodeID *) malloc(sizeof(NodeID) * avgNodesPerCell);
				eachCellContainedNode[i][0] = avgNodesPerCell;
				eachCellContainedNode[i][1] = 2;*/
				
				eachCellTransitNode[i] = (NodeID *) malloc(sizeof(NodeID) * EACH_CELL_TRANSIT_NUM);
				eachCellTransitNode[i][0] = EACH_CELL_TRANSIT_NUM;
				eachCellTransitNode[i][1] = 2;
			}
			
			//eachCellTransitNode.resize(ascale * ascale);
			// for the n * n graph, there are (2n^2 + 2*n) edges		
			initHashMap();
		}
			
		
		TransitNode(Graph* _graph, datastr::graph::UpdateableGraph *_upd, NodeID _ascale, CoordinateType _range): graph(_graph), ascale(_ascale), range(_range), dijkstraCHTest(graph), dijkstraBIDTest(_upd)
		{
			upd = _upd;
			cellLength = range / ascale;
			NodeID totalNodes = graph->noOfNodes();
			avgNodesPerCell = totalNodes / (ascale * ascale);
			if (avgNodesPerCell <= 100){
				avgNodesPerCell = 100;
			}
			allCellNum = ascale * ascale;
			// allocate each cell contained node allocate space for each cell Transit node
			//eachCellContainedNode = (NodeID **) malloc(sizeof(NodeID *) * allCellNum);
			eachCellTransitNode = (NodeID **) malloc(sizeof(NodeID *) * allCellNum);
			for (NodeID i = 0; i < allCellNum ; i++){
				/*eachCellContainedNode[i] = (NodeID *) malloc(sizeof(NodeID) * avgNodesPerCell);
				eachCellContainedNode[i][0] = avgNodesPerCell;
				eachCellContainedNode[i][1] = 2;*/
				
				eachCellTransitNode[i] = (NodeID *) malloc(sizeof(NodeID) * EACH_CELL_TRANSIT_NUM);
				eachCellTransitNode[i][0] = EACH_CELL_TRANSIT_NUM;
				eachCellTransitNode[i][1] = 2;
			}
			
			//eachCellTransitNode.resize(ascale * ascale);
			// for the n * n graph, there are (2n^2 + 2*n) edges		
			initHashMap();
		}

	
		//Constructor, read the transit node, each cell transit node, each cell contained node from file
		//_g is search graph file, 
		//_tfile is transit node file 
		TransitNode(Graph *_g, datastr::graph::UpdateableGraph *_upd, ifstream &_tfile): graph(_g), upd(_upd), dijkstraCHTest(graph), dijkstraBIDTest(_upd)
		{
			deserialize(_tfile);
		}
		
        TransitNode() {};
        
        ~TransitNode() {
            for (NodeID i = 0; i < allCellNum ; ++i) {
                if (eachCellTransitNode != NULL) {
                    free(eachCellTransitNode[i]);
                }
                if (eachCellContainedNode != NULL) {
                    free(eachCellContainedNode[i]);
                }
                if (eachCellTransitMapping != NULL) {
                    free(eachCellTransitMapping[i]);
                }
            }
            free(eachCellTransitNode);
            free(eachCellContainedNode);
            free(eachCellTransitMapping);
            delete[] eachCellDistTable;
            // Note: graph and upd maybe used again so they deleted outside of TNR
        }
        
        // Delayed initialization
        void loadIndex(Graph *_g, datastr::graph::UpdateableGraph *_upd, ifstream &_tfile) {
            graph = _g;
            upd = _upd;
            dijkstraCHTest.loadGraph(_g);
            dijkstraBIDTest.loadGraph(_upd);
            deserialize(_tfile);
        }
        
        void getKNNs(StaticRtree& rtree, unsigned int k, NodeID queryNodeID, std::vector<NodeID>& kNNs, 
                     std::vector<EdgeWeight>& kNNDistances, Coordinate queryNodeX, Coordinate queryNodeY, 
                     std::vector<CoordinateType>& xCoords, std::vector<CoordinateType>& yCoords)
        {
            // Retrieve kNN by euclidean distance
            std::vector<NodeID> euclideanKNNs;
            std::vector<EuclideanDistanceSqrd> euclideanKNNDistances;
            BinaryMinHeap<EuclideanDistanceSqrd,RtreeDataTuple> heap = rtree.getKNNs(k,queryNodeX,queryNodeY,euclideanKNNs,euclideanKNNDistances);
            // Note: We keep heap so that we incrementally retrieve further euclidean NNs

            // We compute the network distances to each of these
            BinaryMaxHeap<EdgeWeight,NodeID> knnCandidates;
            EdgeWeight spDist, Dk, euclidDist;
            for (std::size_t i = 0; i < euclideanKNNs.size(); ++i) {
                spDist = this->transitShortestPathQuery(queryNodeID, euclideanKNNs[i], xCoords, yCoords);
                knnCandidates.insert(euclideanKNNs[i],spDist);
            }

            // While the euclidean distance to the next euclidean nearest neigbour
            // is smaller than than network distance to the current kth neighbour
            // we can potentially find a closer nearest neighbour. Keep searching
            // until this lower bound exceed the kth neighbour network distance.
            Dk = knnCandidates.getMaxKey();
            NodeID nextEuclidNN;
            EuclideanDistanceSqrd currEuclidDistSqrd;
            while (true) {
                if (rtree.getNextNearestNeighbour(heap,queryNodeX,queryNodeY,nextEuclidNN,currEuclidDistSqrd)) {
                    // Note: If there were less than k objects to begin with then getNextNearestNeighbour
                    // would return false (as the heap is empty) and we not reach here
                    euclidDist = std::floor(std::sqrt(currEuclidDistSqrd)); // Floor as it's a lowerbound
                    if (euclidDist < Dk) {
                        spDist = this->transitShortestPathQuery(queryNodeID, nextEuclidNN, xCoords, yCoords);
                        if (spDist < Dk) {
                            // Only insert if it is a better kNN candidate
                            knnCandidates.insert(nextEuclidNN,spDist);
                            knnCandidates.extractMaxElement();
                            Dk = knnCandidates.getMaxKey();
                        }
                    } else {
                        break;
                    }
                } else {
                    // This mean no nearest neighbours were found (we have reported
                    // all objects) so we can stop the search
                    break;
                }
            }

            knnCandidates.populateKNNs(kNNs,kNNDistances);
        }
        
        void getKNNsByTravelTimes(StaticRtree& rtree, unsigned int k, NodeID queryNodeID, std::vector<NodeID>& kNNs, 
                     std::vector<EdgeWeight>& kNNDistances, Coordinate queryNodeX, Coordinate queryNodeY, 
                     std::vector<CoordinateType>& xCoords, std::vector<CoordinateType>& yCoords, double maxGraphSpeed)
        {
            // Retrieve kNN by euclidean distance
            std::vector<NodeID> euclideanKNNs;
            std::vector<EuclideanDistanceSqrd> euclideanKNNDistances;
            BinaryMinHeap<EuclideanDistanceSqrd,RtreeDataTuple> heap = rtree.getKNNs(k,queryNodeX,queryNodeY,euclideanKNNs,euclideanKNNDistances);
            // Note: We keep heap so that we incrementally retrieve further euclidean NNs

            // We compute the network distances to each of these
            BinaryMaxHeap<EdgeWeight,NodeID> knnCandidates;
            EdgeWeight spDist, Dk, minTimeByEuclid;
            for (std::size_t i = 0; i < euclideanKNNs.size(); ++i) {
                spDist = this->transitShortestPathQuery(queryNodeID, euclideanKNNs[i], xCoords, yCoords);
                knnCandidates.insert(euclideanKNNs[i],spDist);
            }

            // While the euclidean distance to the next euclidean nearest neigbour
            // is smaller than than network distance to the current kth neighbour
            // we can potentially find a closer nearest neighbour. Keep searching
            // until this lower bound exceed the kth neighbour network distance.
            Dk = knnCandidates.getMaxKey();
            NodeID nextEuclidNN;
            EuclideanDistanceSqrd currEuclidDistSqrd;
            while (true) {
                if (rtree.getNextNearestNeighbour(heap,queryNodeX,queryNodeY,nextEuclidNN,currEuclidDistSqrd)) {
                    // Note: If there were less than k objects to begin with then getNextNearestNeighbour
                    // would return false (as the heap is empty) and we not reach here
                    minTimeByEuclid = std::ceil(std::sqrt(currEuclidDistSqrd)/maxGraphSpeed); 
                    if (minTimeByEuclid < Dk) {
                        spDist = this->transitShortestPathQuery(queryNodeID, nextEuclidNN, xCoords, yCoords);
                        if (spDist < Dk) {
                            // Only insert if it is a better kNN candidate
                            knnCandidates.insert(nextEuclidNN,spDist);
                            knnCandidates.extractMaxElement();
                            Dk = knnCandidates.getMaxKey();
                        }
                    } else {
                        break;
                    }
                } else {
                    // This mean no nearest neighbours were found (we have reported
                    // all objects) so we can stop the search
                    break;
                }
            }

            knnCandidates.populateKNNs(kNNs,kNNDistances);
        }
        
		// init the BorderEdgeHashMap, each edge is associated with an Bucket
		void initHashMap()
		{
			NodeID edgetotal = 2 * ascale * (ascale + 1);
			edgeBucket.resize(edgetotal);
			for (NodeID u = 0; u < edgetotal; u++){
				edgeBucket[u].reserve(MAX_CROSS_FOR_EDGE);
			}
			//every horizontal border edge with a bucket, insert from left to right
			bool firstTime = true;
			SpaceVertex beg, end;
			beg.set(0, 0);
			end.set(ascale, 0);
			for (NodeID i = 0; i < ascale + 1 ; i++){
				insertEdgeBucketByLine(beg, end, firstTime);	
				beg.increseYPosition();
				end.increseYPosition();
			}
			//every vertical border edge with a bucket, insert from bottom to top
			firstTime = false;
			beg.set(0, 0);
			end.set(0, ascale);
			for (NodeID i = 0; i < ascale + 1 ; i++){
				insertEdgeBucketByLine(beg, end, firstTime);
				beg.increseXPosition();
				end.increseXPosition();
			}
		}
		
	private:
		
		struct Cell{
			int xvertex;
			int yvertex;
			Cell()
			{

			}
			Cell(int _x, int _y): xvertex(_x), yvertex(_y)
			{

			}
			/*int getMappedValue()
			{
				return xvertex * ascale + yvertex;
			}*/
		};
		// get the vertex u's cell position
		void getCell(NodeID u, Cell &cell, vector <CoordinateType> &xcord, vector <CoordinateType> &ycord)
		{
			cell.xvertex = (NodeID) (xcord[u] / cellLength);		//是否先加括号? earnestwu
			cell.yvertex = (NodeID) (ycord[u] / cellLength);		//是否先加括号? earnestwu
			if (xcord[u] == range){
				cell.xvertex = ascale - 1;
			}
			if (ycord[u] == range)
				cell.yvertex = ascale - 1;
		}

		// if two cell is four grid away, return true, else return false
		bool isFourGridAway(Cell &a, Cell &b)
		{
			int x1, x2;
			x1 = (a.xvertex - b.xvertex);
			if (x1 < 0) x1 = -x1;
			x2 = (a.yvertex - b.yvertex);
			if (x2 < 0) x2 = -x2;
			if ( x1 > 4 || x2 > 4)
				return true;
			else
				return false;
		}

		bool isEightGridAway(Cell &a, Cell &b)
		{
			int x1, x2;
			x1 = (a.xvertex - b.xvertex);
			if (x1 < 0) x1 = -x1;
			x2 = (a.yvertex - b.yvertex);
			if (x2 < 0) x2 = -x2;
			if ( x1 > 8 || x2 > 8)
				return true;
			else
				return false;
		}

		struct EasyEdge{
			NodeID source;
			NodeID target;
			EasyEdge(NodeID _source, NodeID _target):source(_source), target(_target)
			{

			}
			NodeID minid()
			{
				if (source < target) return source;
				else 
					return target;
				//return ((source < target) ? (source : target));
			}
		};
		struct SpaceVertex{
			friend ostream& operator<<( ostream& os, const SpaceVertex& sv ) {
				os << " ( "<< sv.x << " , " <<  sv.y << " )";
				return os;
			}
			int x;
			int y;
			SpaceVertex()
			{

			};
			SpaceVertex(int _x, int _y): x(_x), y(_y)
			{

			}
			void set(int _x, int _y)
			{
				x = _x;
				y = _y;
			}
			void setXPosition(int _x)
			{
				x = _x;
			}
			void setYPosition(int _y)
			{
				y = _y;
			}
			void increseXPosition()
			{
				x++;
			}
			void increseYPosition()
			{
				y++;
			}
			int getXPosition()
			{
				return x;
			}
			int getYPosition()
			{
				return y;
			}
		};
		class BorderEdge{
		public:
			BorderEdge()
			{

			}
			BorderEdge(SpaceVertex _start, SpaceVertex _end): start(_start), end(_end)
			{

			}
			void setBorderEdge(SpaceVertex t1, SpaceVertex t2)
			{
				start =  t1;
				end = t2;
			}
			
			void intToChar(NodeID x, char* s) const
			{
				int i = 0;
				while (x){
					s[i++]  = x % 10 + '0';
					x /= 10;
				}
				while (i < 3)
					s[i++] = '0';
				char t;
				t = s[0];
				s[0] = s[2];
				s[2] = t;
				s[3] ='\0';
			}
			
			//hash function, mapping a border edge with a unique ID as the key of hash_map
			int getValue() const
			{
				int temp;
				//temp = 0;
				//temp += start.x;
				//temp *= 31;
				//temp += start.y;
				//temp *= 31;
				//temp += end.x;
				//temp *= 31;
				//temp += end.y;
				temp = 0;
				temp += start.x;
				temp = temp << 8;
				temp += start.y;
				temp = temp << 8;
				temp += end.x;
				temp = temp << 8;
				temp += end.y;
				return temp;
			}
		
			// this opearator overloading is not necessary
			bool operator==(const BorderEdge &bo)
			{
				return (bo.start.x == start.x) && (bo.start.y == start.y) && (bo.end.x == end.x) && (bo.end.y == end.y);
			}
		//private:
			SpaceVertex start;
			SpaceVertex end;
		};

		/*struct string_less : public binary_function<const string, const string, bool>
		{
		public:
			result_type operator()(const first_argument_type& _Left, const second_argument_type& _Right) const
			{
				return(_Left.compare(_Right) < 0 ? true : false);
			}
		};
		typedef hash_map<string, vector <EasyEdge> *, hash_compare<string, string_less>> BorderEdgeHashMap;*/

		typedef unordered_map<int, vector <EasyEdge> *> BorderEdgeHashMap;

		struct Rectangle{

			Rectangle(CoordinateType x1, CoordinateType y1, CoordinateType x2, CoordinateType y2){
				xleft = x1;
				yleft = y1;
				xright = x2;
				yright = y2;
			}
			Rectangle()
			{

			}
			CoordinateType xleft;
			CoordinateType yleft;
			CoordinateType xright;
			CoordinateType yright;
		};

		//// the road nodes, and transit nodes of each cell and the distance table
		//struct EachCellEntry{
		//	//map <NodeID, NodeID> srcMap;
		//	//map <NodeID, NodeID> transitMap;
		//	NodeID **eachCellNode;
		//	NodeID **eachCellTransitNode;
		//	Matrix <EdgeWeight> cellDistTable;
		//};

		bool isNodeInRectangle(NodeID u, Rectangle &rec, 
			vector <CoordinateType> &xcord, vector <CoordinateType> &ycord)
		{
			CoordinateType ux = xcord[u];
			CoordinateType uy = ycord[u];
			return ((ux >= rec.xleft && ux <= rec.xright) && (uy >= rec.yleft && uy <= rec.yright));
		}

		
		bool isEdgeInRectangle(EasyEdge e, Rectangle &rec, 
			vector <CoordinateType> &xcord, vector <CoordinateType> &ycord)
		{
			NodeID u = e.source;
			NodeID v = e.target;
			bool sourceIn = isNodeInRectangle(u, rec, xcord, ycord);
			bool targetIn = isNodeInRectangle(v, rec, xcord, ycord);
			return ((sourceIn && !targetIn) || (!sourceIn && targetIn));
		}

		

public:

		/**
		*using transit node for shortest path query
		*only return the distance between src and target
		* this is basic for shortest query with edges
		*/
		EdgeWeight transitShortestPathQuery(NodeID src, NodeID trg, vector <CoordinateType> &xcord, vector <CoordinateType> &ycord)
		{
			Cell sCell, tCell;
			//double t1, t2;
			getCell(src, sCell, xcord, ycord);
			getCell(trg, tCell, xcord, ycord);
			//src and trg is four gird away
			if (isFourGridAway(sCell, tCell)){
				//m_iLevelOne++;
				NodeID mappedSrc, mappedTrg;
				int srcPos, trgPos;
				mappedSrc = sCell.xvertex * ascale + sCell.yvertex;
				mappedTrg = tCell.xvertex * ascale + tCell.yvertex;

				srcPos = binarysearch(eachCellContainedNode[mappedSrc], 2, eachCellContainedNode[mappedSrc][1], src);
				trgPos = binarysearch(eachCellContainedNode[mappedTrg], 2, eachCellContainedNode[mappedTrg][1], trg);

				// the inner position of source and target in the cell distance table
				srcPos -= 2;
				trgPos -= 2;		
				//NodeID srcAccess, trgAccess;
				if (srcPos >= 0 && trgPos >= 0){
					//NodeID srcTransit, trgTransit;
					NodeID srcTransitMap, trgTransitMap;
					EdgeWeight srcDist, trgDist, transitDist;
					EdgeWeight dist = INT_MAX;

					for (NodeID u = 2; u < eachCellTransitNode[mappedSrc][1]; u++){
						// matrix 's index is begin with 0
						srcDist = eachCellDistTable[mappedSrc].value(srcPos, u - 2); 
						//get the cell's transit node, index begin from 2
						//srcTransit = eachCellTransitNode[mappedSrc][u];
						srcTransitMap = eachCellTransitMapping[mappedSrc][u];
						for (NodeID v = 2; v < eachCellTransitNode[mappedTrg][1]; v++){
							trgDist = eachCellDistTable[mappedTrg].value(trgPos, v - 2);
							//trgTransit = eachCellTransitNode[mappedTrg][v];
							trgTransitMap = eachCellTransitMapping[mappedTrg][v];

							transitDist = transitDistTable.value(srcTransitMap, trgTransitMap);
							if (srcDist + trgDist + transitDist < dist){
								dist = srcDist + trgDist + transitDist;
							}
						}
					}
					return dist;
				}
			}
			else{	
				//m_iLocal++;
				//using CH
				EdgeWeight res = dijkstraCHTest.bidirSearch(src, trg);
				dijkstraCHTest.clear();
				return res;

				////using Bidijkstra
				//EdgeWeight res = dijkstraBIDTest.bidirSearch(src, trg);
				//dijkstraBIDTest.clear();
				//return res;	
			}		
		}


		
		/**
		 * using transit node for shortest path query, only return the distance between src and trg
		 * this is basic for shortest path query with edge
		 * this function also returns the source access node and target access node
		 * i.e. the route is  source ------> source access node -------> target access node ------> target
		 * @srcAccess the source access node
		 * @trgAccess the target access node
		**/
		EdgeWeight transitShortestPathQuery(NodeID src, NodeID trg, vector <CoordinateType> &xcord, vector <CoordinateType> &ycord,
			NodeID &srcAccess, NodeID &trgAccess)
		{
			Cell sCell, tCell;
			getCell(src, sCell, xcord, ycord);
			getCell(trg, tCell, xcord, ycord);
			//src and trg is four gird away
			if (isFourGridAway(sCell, tCell)){
			/*	cout << "non local query: " << src << " ( " << sCell.xvertex << ", "
					<< sCell.yvertex << " ) " << " to " << trg << " (" 
					<< tCell.xvertex << "," << tCell.yvertex<< " ) " << endl;*/
				
				NodeID mappedSrc, mappedTrg;
				int srcPos, trgPos;
				mappedSrc = sCell.xvertex * ascale + sCell.yvertex;
				mappedTrg = tCell.xvertex * ascale + tCell.yvertex;
				//srcTransit = eachCellTransitNode();
				srcPos = binarysearch(eachCellContainedNode[mappedSrc], 2, eachCellContainedNode[mappedSrc][1], src);
				trgPos = binarysearch(eachCellContainedNode[mappedTrg], 2, eachCellContainedNode[mappedTrg][1], trg);
				//cout << "source pos " << srcPos << " target pos " << trgPos << endl;
				// the inner position of source and target in the cell distance table
				srcPos -= 2;
				trgPos -= 2;
				//NodeID srcAccess, trgAccess;
				if (srcPos >= 0 && trgPos >= 0){
					NodeID srcTransit, trgTransit;
					NodeID srcTransitMap, trgTransitMap;
					EdgeWeight srcDist, trgDist, transitDist;
					EdgeWeight dist = INT_MAX;
					for (NodeID u = 2; u < eachCellTransitNode[mappedSrc][1]; u++){
						// matrix 's index is begin with 0
						srcDist = eachCellDistTable[mappedSrc].value(srcPos, u - 2); 
						//get the cell's transit node, index begin from 2
						//srcTransit = eachCellTransitNode[mappedSrc][u];
						srcTransitMap = eachCellTransitMapping[mappedSrc][u];
						for (NodeID v = 2; v < eachCellTransitNode[mappedTrg][1]; v++){
							trgDist = eachCellDistTable[mappedTrg].value(trgPos, v - 2);
							trgTransitMap = eachCellTransitMapping[mappedTrg][v];
							//trgTransit = eachCellTransitNode[mappedTrg][v];
							//cout << srcTransit << " " << trgTransit << endl;
							//cout << transitNodeMapping[srcTransit] << " " << transitNodeMapping[trgTransit] << endl;
							//cout <<transitDistTable.value(transitNodeMapping[srcTransit], transitNodeMapping[trgTransit]) << endl;
							transitDist = transitDistTable.value(srcTransitMap, trgTransitMap);
							if (srcDist + trgDist + transitDist < dist){
								dist = srcDist + trgDist + transitDist;
								trgAccess = eachCellTransitNode[mappedTrg][v];
								srcAccess = eachCellTransitNode[mappedSrc][u];
							}
						}
					}
					//cout << "srcAccess: " << srcAccess << " srcTransit: " << trgAccess << " " << endl;
					return dist;
				}
			}
			else{
				/*cout << "local query: " << src << " ( " << sCell.xvertex << ", "
					<< sCell.yvertex << " ) " << " to " << trg << " (" 
					<< tCell.xvertex << "," << tCell.yvertex<< " ) " << endl;*/
				EdgeWeight res = dijkstraCHTest.bidirSearch(src, trg);
				dijkstraCHTest.clear();
				return res;

				//EdgeWeight res = dijkstraBIDTest.bidirSearch(src, trg);
				//dijkstraBIDTest.clear();
				//return res;
			}			
		}
		
		////only used for path outputting
		//EdgeWeight distanceToTransitNode(NodeID src, NodeID trgTransit, vector <CoordinateType> &xcord, vector <CoordinateType> &ycord,
		//	NodeID &srcAccess)
		//{
		//	Cell sCell, tCell;
		//	getCell(src, sCell, xcord, ycord);
		//	getCell(trgTransit, tCell, xcord, ycord);
		//	//src and trg is four gird away
		//	if (isFourGridAway(sCell, tCell)){
		//		/*cout << "non local query: " << src << " ( " << sCell.xvertex << ", "
		//			<< sCell.yvertex << " ) " << " to " << trgTransit << " (" 
		//			<< tCell.xvertex << "," << tCell.yvertex<< " ) " << endl;	
		//		cout << src << " to the trgTransit " << trgTransit << endl;*/
		//		NodeID mappedSrc;
		//		int srcPos;
		//		mappedSrc = sCell.xvertex * ascale + sCell.yvertex;
		//		//mappedTrg = tCell.xvertex * ascale + tCell.yvertex;
		//		srcPos = binarysearch(eachCellContainedNode[mappedSrc], 2, eachCellContainedNode[mappedSrc][1], src);
		//		//cout << "source pos " << srcPos << endl;
		//		// the inner position of source and target in the cell distance table
		//		srcPos -= 2;
		//		NodeID trgPos;
		//		//NodeID srcAccess, trgAccess;
		//		if (srcPos >= 0){
		//			NodeID srcTransit;
		//			NodeID srcTransitMap;
		//			EdgeWeight srcDist, transitDist;
		//			EdgeWeight dist = INT_MAX;
		//			for (NodeID u = 2; u < eachCellTransitNode[mappedSrc][1]; u++){
		//				// matrix 's index is begin with 0
		//				srcDist = eachCellDistTable[mappedSrc].value(srcPos, u - 2); 
		//				//get the cell's transit node, index begin from 2
		//				srcTransitMap = eachCellTransitMapping[mappedSrc][u];
		//				//srcTransit = eachCellTransitNode[mappedSrc][u];
		//				trgPos = transitNodeMapping[trgTransit];
		//				transitDist = transitDistTable.value(srcTransitMap, trgPos);
		//				if (srcDist + transitDist < dist){
		//					dist = srcDist + transitDist;
		//					srcAccess = eachCellTransitNode[mappedSrc][u];
		//				}
		//			}
		//			//cout << "srcAccess: " << srcAccess << endl;
		//			return dist;
		//		}
		//	}
		//	else{
		//		/*cout << "local query: " << src << " ( " << sCell.xvertex << ", "
		//			<< sCell.yvertex << " ) " << " to " << trgTransit << " (" 
		//			<< tCell.xvertex << "," << tCell.yvertex<< " ) " << endl;*/
		//		processing::DijkstraCH<datastr::graph::SearchGraph, NormalPQueue, 2, false> dijkstraCHTest(graph);
		//		EdgeWeight res = dijkstraCHTest.bidirSearch(src, trgTransit);
		//		dijkstraCHTest.clear();
		//		//cout << src << " to the trgTransit " << trgTransit << " : " << res << endl;
		//		return res;
		//	}			
		//}





		/**
		* Get the distance from source to one transit node in level one, need less table looking up compared with normal
		* distance query
		* it is similar with levelOne_transitShortestPathQuery
		* only used for path outputting
		*/

		EdgeWeight distanceToTransitNode(datastr::graph::UpdateableGraph* updGraph, NodeID src, NodeID trgTransit, vector <CoordinateType> &xcord, vector <CoordinateType> &ycord,
			NodeID &srcAccess)
		{
			Cell sCell, tCell;
			getCell(src, sCell,  xcord, ycord);
			getCell(trgTransit, tCell, xcord, ycord);
			//src and trg is four gird away
			if (isFourGridAway(sCell, tCell)){
				NodeID mappedSrc;
				int srcPos;
				mappedSrc = sCell.xvertex * ascale + sCell.yvertex;		
				srcPos = binarysearch(eachCellContainedNode[mappedSrc], 2, eachCellContainedNode[mappedSrc][1], src);
				//the inner position of source and target in the cell distance table
				srcPos -= 2;
				NodeID trgPos;
				if (srcPos >= 0){
					NodeID srcTransitMap;
					EdgeWeight srcDist, transitDist;
					EdgeWeight dist = INT_MAX;
					trgPos = transitNodeMapping[trgTransit];
					for (NodeID u = 2; u < eachCellTransitNode[mappedSrc][1]; u++){
						//matrix 's index begins with 0
						srcDist = eachCellDistTable[mappedSrc].value(srcPos, u - 2); 
						//get the cell's transit node, index begins from 2
						srcTransitMap = eachCellTransitMapping[mappedSrc][u];
						//trgPos = transitNodeMapping[trgTransit];
						transitDist = transitDistTable.value(srcTransitMap, trgPos);
						if (srcDist + transitDist < dist){
							dist = srcDist + transitDist;
							srcAccess = eachCellTransitNode[mappedSrc][u];
						}
					}
					return dist;
				}
			}
			else{
				EdgeWeight res = dijkstraCHTest.bidirSearch(src, trgTransit);
				dijkstraCHTest.clear();
				return res;
			}			
		}


	
		//find the target access node and the distance from source to this node in level one
		//similar with levelOne_transitShortestPathQuery
		EdgeWeight transitShortestPathQuery_GetDist(datastr::graph::UpdateableGraph *updGraph, NodeID src, NodeID trg, vector <CoordinateType> &xcord, vector <CoordinateType> &ycord,
			NodeID &srcAccess, NodeID &trgAccess, EdgeWeight &srcToSAccessDistance, EdgeWeight &srcToTAccessDistance)
		{
			Cell sCell, tCell;
			getCell(src, sCell,  xcord, ycord);
			getCell(trg, tCell,  xcord, ycord);

			//src and trg are four girds away in level one
			if (isFourGridAway(sCell, tCell)){
				NodeID mappedSrc, mappedTrg;
				int srcPos, trgPos;
				mappedSrc = sCell.xvertex * ascale + sCell.yvertex;
				mappedTrg = tCell.xvertex * ascale + tCell.yvertex;

				srcPos = binarysearch(eachCellContainedNode[mappedSrc], 2, eachCellContainedNode[mappedSrc][1], src);
				trgPos = binarysearch(eachCellContainedNode[mappedTrg], 2, eachCellContainedNode[mappedTrg][1], trg);
				// the inner position of source and target in the cell distance table
				srcPos -= 2;
				trgPos -= 2;
				//NodeID srcAccess, trgAccess;
				if (srcPos >= 0 && trgPos >= 0){
					NodeID srcTransitMap, trgTransitMap;
					EdgeWeight srcDist, trgDist, transitDist;
					EdgeWeight dist = INT_MAX;
					for (NodeID u = 2; u < eachCellTransitNode[mappedSrc][1]; u++){
						//matrix 's index begins with 0
						srcDist = eachCellDistTable[mappedSrc].value(srcPos, u - 2); 
						//get the cell's transit node, index begins from 2
						srcTransitMap = eachCellTransitMapping[mappedSrc][u];
						for (NodeID v = 2; v < eachCellTransitNode[mappedTrg][1]; v++){
							trgDist = eachCellDistTable[mappedTrg].value(trgPos, v - 2);
							trgTransitMap = eachCellTransitMapping[mappedTrg][v];
							transitDist = transitDistTable.value(srcTransitMap, trgTransitMap);

							if (srcDist + trgDist + transitDist < dist){
								dist = srcDist + trgDist + transitDist;
								trgAccess = eachCellTransitNode[mappedTrg][v];
								srcAccess = eachCellTransitNode[mappedSrc][u];
								srcToSAccessDistance = srcDist;
								srcToTAccessDistance = srcDist + transitDist;
							}
						}
					}
					return dist;
				}
			}
			else{
				srcAccess = SPECIAL_NODEID;
				trgAccess = SPECIAL_NODEID;									
				EdgeWeight res = dijkstraCHTest.bidirSearch(src, trg);
				dijkstraCHTest.clear();
				return res;
			}			
		}

		
		/**
		  *	using dist table output the shortest path with edges
		  * when src and target is larger than eight grids away, we could only use distance table to output the path,
		  * we start from source to the most farthest node which is 4 grids away to the target, then we change from 
		  * target to meet the node we find before
		  * if the src and target is less than eight grids away, the we also start from source to the most farthest node
		  * which is 4 grids away to the target, then from target to the the most farthest nodes which is 4 girds away 
		  * source, then combine this gap together
		  * using sourceToFourGridAwayTarget() to get the most farthest node which is 4 grids away to the target,
		  * using transitShortestPathQueryWithEdge() to get the gap route  
		*/
		EasyPath transitShortestPathQuery_UsingDist(datastr::graph::UpdateableGraph *updGraph, NodeID src, NodeID trg, vector <CoordinateType> &xcord, vector <CoordinateType> &ycord)
		{
			Cell sCell, tCell;
			EasyPath epathF, epathR, epathM;
			NodeID endingVertex;
			getCell(src, sCell, xcord, ycord);
			getCell(trg, tCell, xcord, ycord);
			endingVertex = SPECIAL_NODEID;
			
			if (isEightGridAway(sCell, tCell)){
				epathF = sourceToTargetAccessNode(updGraph, src, trg, true, endingVertex, xcord, ycord);
				//epathF = sourceToFourGridAwayTarget(updGraph, src, trg, true, endingVertex, xcord, ycord);
				endingVertex = epathF.lastNode();
				epathR = sourceToTargetAccessNode(updGraph, trg, src, false, endingVertex, xcord, ycord);
				//epathR = sourceToFourGridAwayTarget(updGraph, trg, src, false, endingVertex, xcord, ycord);		
				if (epathF.lastNode() != epathR.firstNode()){
					epathM = transitShortestPathQueryWithEdge(updGraph, epathF.lastNode(), epathR.firstNode(), xcord, ycord);
					epathF.add(epathM);
				}
				epathF.add(epathR);
			}
			else{
				//if(isFourGridAway(sCell, tCell)){
				//	epathF = sourceToTargetAccessNode(updGraph, src, trg, true, endingVertex, xcord, ycord);
				//	//epathF = sourceToFourGridAwayTarget(updGraph, src, trg, true, endingVertex, xcord, ycord);	
				//	//reverse the path in the function
				//	epathR = sourceToTargetAccessNode(updGraph, trg, src, false, endingVertex, xcord, ycord);
				//	//epathR = sourceToFourGridAwayTarget(updGraph, trg, src, false, endingVertex, xcord, ycord);
				//	epathM = transitShortestPathQueryWithEdge(updGraph, epathF.lastNode(), epathR.firstNode(), xcord, ycord);
				//	epathF.add(epathM);
				//	epathF.add(epathR);
				//}
				//else{				 
					//CHANGEING
					Path p;
					
					EdgeWeight res = dijkstraCHTest.bidirSearch(src, trg);
					dijkstraCHTest.pathTo(p, trg, -1, true, true);
					dijkstraCHTest.clear();

					////need changing if using dijkstra for local query
					//EdgeWeight res = dijkstraBIDTest.bidirSearch(src, trg);
					//dijkstraBIDTest.pathTo(p, trg, -1, true, true);
					//dijkstraBIDTest.clear();
					
					return EasyPath(p);
				//}
			}
			return epathF;
		}
			
		/**
		* path query function used in level one, find the path from source node to its target access node
		* @param updGraph the original graph
		* @param src source node
		* @param trg target node
		* @param xcord x coordinate
		* @param ycord y coordinate
		* using levelOne_transitShortestPathQuery_GetDist to find the target access node and the distance from source to this node
		*/
		EasyPath sourceToTargetAccessNode(datastr::graph::UpdateableGraph *updGraph, NodeID src, NodeID trg, bool forward,NodeID endingVertex, 
			vector <CoordinateType> &xcord, vector <CoordinateType> &ycord)
		{

			NodeID curNode, curTrg;
			NodeID curAccess, srcAccess, trgAccess;

			EdgeWeight curDist, srcToSAccessDistance, srcToTAccessDistance, transitDistance;
			EdgeWeight eweight;
			EdgeWeight x, y, trgAccessDist;

			Cell sCell, tCell, curCell, tAccessCell, curNodeCell;
			getCell(trg, tCell, xcord, ycord);
			getCell(src, sCell, xcord, ycord);

			EasyPath p(src);
			if (src == trg)
				return p;
			curNode =  src;	

			//get the target access node and the distance from source to target access node
			transitShortestPathQuery_GetDist(updGraph, curNode, trg, xcord, ycord, srcAccess, trgAccess, srcToSAccessDistance, srcToTAccessDistance);
			transitDistance = srcToTAccessDistance - srcToSAccessDistance;

			getCell(trgAccess, tAccessCell, xcord, ycord);
			curDist = 0;
			int first, last, edgeindex;
			Edge edge;
			set <NodeID> addNode;
			addNode.insert(src);
			set <NodeID>::iterator sit;
			while (true){
				first = updGraph->firstLevelEdge(curNode);
				last = updGraph->lastEdge(curNode);
				for (edgeindex = first; edgeindex < last; edgeindex++){
					edge = updGraph->edge(edgeindex);
					curTrg = edge.target();
					sit = addNode.find(curTrg);
					if (sit != addNode.end()) continue;
					getCell(curTrg, curCell, xcord, ycord);
					if (!isFourGridAway(curCell, tAccessCell))
						continue;
					eweight = edge.weight();
					//Get the distance from source to target access node
					trgAccessDist = distanceToTransitNode(updGraph, curTrg, trgAccess, xcord, ycord, curAccess);
					x = srcToTAccessDistance - trgAccessDist;
					if ((curDist + eweight) == x){
						p.add(curTrg);
						addNode.insert(curTrg);		
						if (curTrg == endingVertex || curTrg == trgAccess){
							if (!forward)
								p.reverse();
							return p;
						}
						curNode = curTrg;
						curDist += eweight;
						break;
					}
				}
				getCell(curNode, curNodeCell, xcord, ycord);
				if (edgeindex == last || !isFourGridAway(curNodeCell, tAccessCell)){
					if (!forward)
						p.reverse();				
					return p;
				}
			}
		}


		
		/**
		* this function is just used to connect the gap and output the path
		*/
		EasyPath transitShortestPathQueryWithEdge(datastr::graph::UpdateableGraph *updGraph, NodeID src, NodeID trg, vector <CoordinateType> &xcord, vector <CoordinateType> &ycord)
		{
			Path p;
			EasyPath ep;
			//need changing
			//using ch
			EdgeWeight res = dijkstraCHTest.bidirSearch(src, trg);
			dijkstraCHTest.pathTo(p, trg, -1, true, true);
			dijkstraCHTest.clear();

			////using bijijkstra, need change if bidijkstra for local query
			//	
			//EdgeWeight res = dijkstraBIDTest.bidirSearch(src, trg);
			//dijkstraBIDTest.pathTo(p, trg, -1, true, true);	
			//dijkstraBIDTest.clear();
			
			ep.setPath(p);
			return ep;
		}	
		
		void computeAccessNodeForEdge(datastr::graph::UpdateableGraph *updGraph, vector <CoordinateType> &xcord, vector <CoordinateType> &ycord, fstream &out)
		{
			for (NodeID u = 0; u < updGraph->noOfNodes(); u++){
				OtherEdgeID first = updGraph->firstEdge(u);
				OtherEdgeID last = updGraph->lastEdge(u);
				//cout << first << " " << last << endl;
				for (OtherEdgeID edgeindex = first; edgeindex < last; edgeindex++){
					Edge edge = updGraph->edge(edgeindex);
					OtherEdgeID revid = updGraph->reverseEdge(u, edgeindex);
					// the edge has been tested
					if (revid < edgeindex) continue;
					NodeID sou = u;
					NodeID tar = edge.target();
					insertToBorderEdgeHashMap(sou, tar, xcord, ycord, out);
				}
			}
		}
		
		// compute all cell's transit node by calculate the transit node one cell by one cell
		void computeAllCellTransitNode(vector <CoordinateType> &xcord, vector <CoordinateType> &ycord, fstream &out)
		{
			Cell c;
			for (NodeID i = 0; i < ascale; i++){
				for (NodeID j = 0; j < ascale; j++)
				{
					if ((i * ascale + j ) % 1000 == 0)
						cout << "now computing cell ( " << i << " " << j<< " )'s transit nodes" << endl;
					c.xvertex = i;
					c.yvertex = j;
					computeTransitNodeForCell(c, xcord, ycord, out);
				}
			}
			//out << allCellTransitNode.size() << endl;
			//set <NodeID>::iterator sit = allCellTransitNode.begin();
			//int count = 0;
			//while(sit != allCellTransitNode.end()){
			//	out << *sit << " ";
			//	count++;
			//	sit++;
			//	if (count % 10 == 0)	out << endl;
			//}

		}
		// given a cell, compute this cell's transit node
		void computeTransitNodeForCell(Cell &cell, vector <CoordinateType> &xcord, vector <CoordinateType> &ycord, fstream &out)
		{
			vector <NodeID> startBoundNode;
			vector <NodeID> innerBoundNode;
			vector <NodeID> outerBoundNode;
			Cell innerCell(cell.xvertex - 2, cell.yvertex - 2);
			Cell outerCell(cell.xvertex - 4, cell.yvertex - 4);
			
			startBoundNode.reserve(MAX_TRANSITNODE_FOR_CELL);
			innerBoundNode.reserve(MAX_TRANSITNODE_FOR_CELL);
			outerBoundNode.reserve(MAX_TRANSITNODE_FOR_CELL);
#ifdef DEBUG
			//out << "inner" << innerCell.xvertex << " " << innerCell.yvertex << endl;
			//out << "outer" << outerCell.xvertex << " " << outerCell.yvertex << endl;
#endif
			
			getStartVertexAccessNode(cell, startBoundNode, xcord, ycord, out);
			getInnerAccessNode(innerCell, innerBoundNode, xcord, ycord, out);
			getOuterAccessNode(outerCell, outerBoundNode, xcord, ycord, out);
			removeDupAccessNodes(startBoundNode);
			removeDupAccessNodes(innerBoundNode);
			removeDupAccessNodes(outerBoundNode);
#ifdef DEBUG
			out << endl;
			out << "after remove duplicates " << endl;
			out << "start size: " << startBoundNode.size() << endl;
			for (NodeID u = 0; u < startBoundNode.size(); u++)
				out << startBoundNode[u] << " " ;
			out << endl;

			out << "inner size: " << innerBoundNode.size() << endl;
			for (NodeID u = 0; u < innerBoundNode.size(); u++)
				out << innerBoundNode[u] << " " ;
			out << endl;

			out << "outer size: " << outerBoundNode.size() << endl;
			for (NodeID u = 0; u < outerBoundNode.size(); u++)
				out << outerBoundNode[u] << " " ;
			out << endl;
#endif

			selectTransitNodeForCell(cell, startBoundNode, innerBoundNode, outerBoundNode, xcord, ycord);

			//out << endl;
			//out << " cell : (" << cell.xvertex << " " << cell.yvertex << ") 's transitnode: ";
			//NodeID mappedValue = cell.xvertex * ascale + cell.yvertex;
			//out << eachCellTransitNode[mappedValue][1] - 2 << endl;
			//for (NodeID u = 2; u < eachCellTransitNode[mappedValue][1]; u++){
			//	out << eachCellTransitNode[mappedValue][u] << " ";
			//	if ((u - 1) % 10 == 0) out << endl;
			//}
			/*out << eachCellTransitNode[mappedValue].size() << endl;
			set <NodeID>::iterator sit = eachCellTransitNode[mappedValue].begin();
			int count = 0;
			while (sit != eachCellTransitNode[mappedValue].end()){
				out << *sit << " ";
				sit++;
				count++;
				if (count % 10 == 0) out << endl;
			}*/
			//out << endl;
		}
		
		typedef map <pair<NodeID, NodeID>, pair< NodeID, EdgeWeight> > TransitCandidateType;
        
        // Note: Commented out code below because computeAllCellTransitNode_BySweepLine is not being used 
        // anywhere and it is using a function which is causing an error.
        
// 		void computeAllCellTransitNode_BySweepLine(vector <CoordinateType> &xcord, vector <CoordinateType> &ycord, fstream &out)
// 		{
// 			unsigned int i;
// 			SpaceVertex start, end;
// 			TransitCandidateType transitCandidate; 
// 			bool vertical = true;
// 	
// 
// 			// horizon sweep line
// 		
// 			for (i = 3; i < ascale - 2; i++){
// 				cout << "horizon line " << i << endl;
// 				computeHorizonSweepLineTransitNode(i, transitCandidate, xcord, ycord);
// 				addToTransitNodeSet(transitCandidate, xcord, ycord);
// 				transitCandidate.clear();	
// 			}
// 				// vertical sweep line
// 			for (i = 3; i < ascale - 2; i++){
// 				cout << "vertical line " << i << endl;
// 				computeVerticalSweepLineTransitNode(i, transitCandidate, xcord, ycord);
// 				addToTransitNodeSet(transitCandidate, xcord, ycord);
// 				transitCandidate.clear();
// 				
// 			}
// 			
// 			//addToTransitNodeSet(transitCandidate, xcord, ycord);
// 			out << allCellTransitNode.size() << endl;
// 			set <NodeID>::iterator sit = allCellTransitNode.begin();
// 			int count = 0;
// 			while(sit != allCellTransitNode.end()){
// 				out << *sit << " ";
// 				count++;
// 				sit++;
// 				if (count % 10 == 0)	out << endl;
// 			}
// 			cout << "writing all transit node... done " << endl;
// 		}

// 		void computeVerticalSweepLineTransitNode(NodeID line, TransitCandidateType & transitCandidate,
// 			vector <CoordinateType> &xcord, vector <CoordinateType> &ycord)
// 		{
// 			SpaceVertex start, end;
// 			SpaceVertex left1, left2, left3, left4, right1, right2, right3, right4, t1, t2, t3, t4;
// 			int x1, x2, x3, x4, y1, y2, y3, y4;
// 			int i, j, k, section;
// 			NodeID v;
// 			EdgeWeight newDis;
// 			vector <NodeID> leftBoundary, rightBoundary;
// 			vector <NodeID> boundaryNode;
// 			int leftSize, rightSize;
// 			start.set(line, 0);
// 			end.set(line, 1);
// 			
// 
// 			//compute the minimum vl + vr
// 			for (section = 0; section < ascale; section++){	
// 				//cout << section << endl;
// 				BorderEdge bo(start, end);				
// 				start.increseYPosition();
// 				end.increseYPosition();	
// 				vector <EasyEdge> copy(*beHashMap[bo.getValue()]);
// 				//cout  << "candidate size " << copy.size() << endl;
// 				if (copy.empty()) continue;		
// 				leftBoundary.clear();
// 				rightBoundary.clear();
// 				boundaryNode.clear();
// 				// get the left boundary and right boundary node
// 				x1 = start.getXPosition() - 3;
// 				x2 = start.getXPosition() - 2;
// 				x3 = start.getXPosition() + 2;
// 				x4 = start.getXPosition() + 3;
// 				//cout << x1 << " " << x2 << " " << x3 << " " << x4 << endl;
// 				if (start.getYPosition() - 2 > 0){
// 					y1 = start.getYPosition() - 2;
// 				}
// 				else
// 					y1 = 0;
// 				if (end.getYPosition() + 3 <= ascale){
// 					y2 = start.getYPosition() + 3;
// 				}
// 				else
// 					y2 = ascale;
// 				//y1 = 0;
// 				//y2 = ascale;
// 				//cout << "y1 and y2 " << y1 << " " << y2 << endl;
// 				left1.set(x1, y1);
// 				left2.set(x2, y1);
// 				left3.set(x1, y2);
// 				left4.set(x2, y2);
// 				right1.set(x3, y1);
// 				right2.set(x4, y1);
// 				right3.set(x3, y2);
// 				right4.set(x4, y2);
// 				getAccessNodeByLine(left1, left3, leftBoundary, xcord, ycord);
// 				getAccessNodeByLine(left2, left4, leftBoundary, xcord, ycord);
// 				getAccessNodeByLine(right1, right3, rightBoundary, xcord, ycord);
// 				getAccessNodeByLine(right2, right4, rightBoundary, xcord, ycord);
// 				//cout << "line size " << leftBoundary.size() << " " << rightBoundary.size() << endl;
// 				while (y1 <= y2){
// 					t1.set(x1, y1);
// 					t2.set(x2, y1);
// 					t3.set(x3, y1);
// 					t4.set(x4, y1);
// 					BorderEdge b1(t1, t2);
// 					BorderEdge b2(t3, t4);
// 					vector <EasyEdge> e1(*beHashMap[b1.getValue()]);				
// 					vector <EasyEdge> e2(*beHashMap[b2.getValue()]);
// 					//cout << "e1 and e2 size : "  << e1.size() << " " << e2.size() << endl;
// 					for(NodeID u = 0; u < e1.size(); u++){
// 						leftBoundary.push_back(e1[u].minid());
// 					}
// 					for(NodeID u = 0; u < e2.size(); u++){
// 						rightBoundary.push_back(e2[u].minid());
// 					}
// 					t1.increseYPosition();
// 					t2.increseYPosition();
// 					t3.increseYPosition();
// 					t4.increseYPosition();
// 					y1++;
// 				}
// 				removeDupAccessNodes(leftBoundary);
// 				removeDupAccessNodes(rightBoundary);
// 				if (leftBoundary.size() == 0 || rightBoundary.size() == 0)
// 					continue;
// 				//cout << "remove duplicate: line size " << leftBoundary.size() << " " << rightBoundary.size() << endl;		
// 				leftSize = leftBoundary.size();
// 				rightSize = rightBoundary.size();
// 				for (k = 0; k < leftSize; k++)
// 					boundaryNode.push_back(leftBoundary[k]);
// 				for (k = 0; k < rightSize; k++)
// 					boundaryNode.push_back(rightBoundary[k]);
// 				//removeDupAccessNodes(boundaryNode);
// 				vector <NodeID> edgeNodes;
// 				for (NodeID u = 0; u < copy.size(); u++)
// 					edgeNodes.push_back(copy[u].minid());
// 				int earlyStopLevel = 10;
// 				const bool performBucketScans = true;
// 				ManyToMany<MyGraph, DijkstraManyToManyFW, DijkstraManyToManyBW, performBucketScans> mtmOuter(graph, earlyStopLevel);
// 				Matrix <EdgeWeight> ma(copy.size(), boundaryNode.size());
// 				mtmOuter.computeMatrix(edgeNodes, boundaryNode, ma);
// 				//cout << "totall size " << boundaryNode.size() << endl;
// 				processing::DijkstraCH<datastr::graph::UpdateableGraph, NormalPQueue, 2, false> dijkstraCHTest(upd);
// 				TransitCandidateType::iterator mit;
// 				for(NodeID u = 0; u < copy.size(); u++){
// 					EasyEdge ee = copy[u];					
// 					//dijkstraCHTest.localQueryForTransitNode(ee.minid(), boundaryNode);
// 					//cout << "check vertex " << ee.minid() << endl;
// 					for (i = 0; i < leftSize; i++){
// 						for (j = 0; j < rightSize; j++){
// 							newDis = ma.value(u, i) + ma.value(u, leftSize + j);
// 							//newDis = dijkstraCHTest.distanceTo(leftBoundary[i]) + dijkstraCHTest.distanceTo(rightBoundary[j]);
// 							mit = transitCandidate.find(make_pair(leftBoundary[i], rightBoundary[j]));							
// 							if (mit != transitCandidate.end()){
// 								if (newDis < mit->second.second){
// 
// 									mit->second.first = ee.minid();
// 									mit->second.second = newDis;
// 								}
// 							}
// 							else{
// 								transitCandidate.insert(make_pair(make_pair(leftBoundary[i], rightBoundary[j]), make_pair(ee.minid(), newDis)));
// 							}
// 						}
// 					}
// 					//dijkstraCHTest.clear();
// 				}
// 				mtmOuter.clear();
// 								
// 			}
// 		}
// 
// 
// 
// 		void computeHorizonSweepLineTransitNode(NodeID line, TransitCandidateType & transitCandidate,
// 			vector <CoordinateType> &xcord, vector <CoordinateType> &ycord)
// 		{
// 			SpaceVertex start, end;
// 			SpaceVertex bottom1, bottom2, bottom3, bottom4, top1, top2, top3, top4, t1, t2, t3, t4;
// 			int x1, x2, y1, y2, y3, y4;
// 			int i, j, k, section;
// 			NodeID v;
// 			EdgeWeight newDis;
// 			vector <NodeID> topBoundary, bottomBoundary;
// 			vector <NodeID> boundaryNode;
// 			int topSize, bottomSize;
// 			start.set(0, line);
// 			end.set(1, line);
// 			for (section = 0; section < ascale; section++){	
// 				BorderEdge bo(start, end);				
// 				start.increseXPosition();
// 				end.increseXPosition();	
// 				vector <EasyEdge> copy(*beHashMap[bo.getValue()]);
// 				if (copy.empty()) continue;
// 				bottomBoundary.clear();
// 				topBoundary.clear();
// 				boundaryNode.clear();
// 
// 
// 				y1 = start.getYPosition() - 3;
// 				y2 = start.getYPosition() - 2;
// 				y3 = start.getYPosition() + 2;
// 				y4 = start.getYPosition() + 3;
// 				//cout << y1 << " " << y2 << " " << y3 << " " << y4 << endl;
// 				if (start.getXPosition() - 2 > 0){
// 					x1 = start.getXPosition() - 2;
// 				}
// 				else
// 					x1 = 0;
// 				if (end.getXPosition() + 3 <= ascale){
// 					x2 = start.getXPosition() + 3;
// 				}
// 				else
// 					x2 = ascale;
// 				bottom1.set(x1, y1);
// 				bottom2.set(x2, y1);
// 				bottom3.set(x1, y2);
// 				bottom4.set(x2, y2);
// 				top1.set(x1, y3);
// 				top2.set(x2, y3);
// 				top3.set(x1, y4);
// 				top4.set(x2, y4);
// 				getAccessNodeByLine(bottom1, bottom2, bottomBoundary, xcord, ycord);
// 				getAccessNodeByLine(bottom3, bottom4, bottomBoundary, xcord, ycord);
// 				getAccessNodeByLine(top1, top2, topBoundary, xcord, ycord);
// 				getAccessNodeByLine(top3, top4, topBoundary, xcord, ycord);
// 				while (x1 <= x2){
// 					t1.set(x1, y1);
// 					t2.set(x1, y2);
// 					t3.set(x1, y3);
// 					t4.set(x1, y4);
// 					BorderEdge b1(t1, t2);
// 					BorderEdge b2(t3, t4);
// 					vector <EasyEdge> e1(*beHashMap[b1.getValue()]);				
// 					vector <EasyEdge> e2(*beHashMap[b2.getValue()]);
// 					for(NodeID u = 0; u < e1.size(); u++){
// 						bottomBoundary.push_back(e1[u].minid());
// 					}
// 					for(NodeID u = 0; u < e2.size(); u++){
// 						topBoundary.push_back(e2[u].minid());
// 					}
// 					t1.increseXPosition();
// 					t2.increseXPosition();
// 					t3.increseXPosition();
// 					t4.increseXPosition();
// 					x1++;
// 				}
// 				removeDupAccessNodes(bottomBoundary);
// 				removeDupAccessNodes(topBoundary);
// 
// 				if (bottomBoundary.size() ==  0 || topBoundary.size() == 0)
// 					continue;
// 
// 				//cout << "two size " << bottomBoundary.size() << " " << topBoundary.size() << endl;
// 
// 				bottomSize = bottomBoundary.size();
// 				topSize = topBoundary.size();
// 				for (k = 0; k < bottomSize; k++)
// 					boundaryNode.push_back(bottomBoundary[k]);
// 				for (k = 0; k < topSize; k++)
// 					boundaryNode.push_back(topBoundary[k]);
// 				//removeDupAccessNodes(boundaryNode);
// 				
// 				vector <NodeID> edgeNodes;
// 				for (NodeID u = 0; u < copy.size(); u++)
// 					edgeNodes.push_back(copy[u].minid());
// 				int earlyStopLevel = 10;
// 				const bool performBucketScans = true;
// 				ManyToMany<MyGraph, DijkstraManyToManyFW, DijkstraManyToManyBW, performBucketScans> mtmOuter(graph, earlyStopLevel);
// 				Matrix <EdgeWeight> ma(copy.size(), boundaryNode.size());
// 				mtmOuter.computeMatrix(edgeNodes, boundaryNode, ma);
// 				//cout << bottomSize << " " << topSize << endl;
// 				//cout  << "many to many "<< ma << endl;
// 				//processing::DijkstraCH<datastr::graph::UpdateableGraph, NormalPQueue, 2, false> dijkstraCHTest(upd);
// 				TransitCandidateType::iterator mit;
// 				for(NodeID u = 0; u < copy.size(); u++){
// 					EasyEdge ee = copy[u];					
// 					//dijkstraCHTest.localQueryForTransitNode(ee.minid(), boundaryNode);
// 					for (i = 0; i < bottomSize; i++){
// 						for (j = 0; j < topSize; j++){
// 							newDis = ma.value(u, i) + ma.value(u, bottomSize + j);
// 							//newDis = dijkstraCHTest.distanceTo(bottomBoundary[i]) + dijkstraCHTest.distanceTo(topBoundary[j]);
// 							/*if (newDis != (ma.value(u, i) + ma.value(u, bottomSize + j))){
// 								cout << "( i, j ) "<< i << " " << j << endl;
// 								cout << ma.value(u, i) << " " << ma.value(u, bottomSize + j) << endl;
// 								cout << dijkstraCHTest.distanceTo(bottomBoundary[i]) << " " << dijkstraCHTest.distanceTo(topBoundary[j]) << endl;
// 								cout << newDis << " " << ma.value(u, i) + ma.value(u, bottomSize + j) << endl;
// 								exit(1);
// 							}*/
// 							mit = transitCandidate.find(make_pair(bottomBoundary[i], topBoundary[j]));							
// 							if (mit != transitCandidate.end()){
// 								if (newDis < mit->second.second){
// 									mit->second.first = ee.minid();
// 									mit->second.second = newDis;
// 								}
// 							}
// 							else{
// 								transitCandidate.insert(make_pair(make_pair(bottomBoundary[i], topBoundary[j]), make_pair(ee.minid(), newDis)));
// 							}
// 						}
// 					}
// 					//dijkstraCHTest.clear();
// 				}
// 				mtmOuter.clear();
// 								
// 			}
// 		}

		void addToTransitNodeSet(TransitCandidateType &transitCandidate, vector <CoordinateType> &xcord,
			vector <CoordinateType> &ycord)
		{
			NodeID src, trg;
			int mappedSrc, mappedTrg;
			Cell srcCell, trgCell;
			int low, high;
			int pos;
			TransitCandidateType::iterator mit = transitCandidate.begin();
			//cout << "add to transit node " << transitCandidate.size() << endl;
			while (mit != transitCandidate.end()){
				src = mit->first.first;
				trg = mit->first.second;
				getCell(src, srcCell, xcord, ycord);
				getCell(trg, trgCell, xcord, ycord);
				allCellTransitNode.insert(mit->second.first);
				//cout << mit->second.first << endl;
				mappedSrc = srcCell.xvertex * ascale + srcCell.yvertex;
				mappedTrg = trgCell.xvertex * ascale + trgCell.yvertex;
				low = 2;
				high =  eachCellTransitNode[mappedSrc][1];	
				pos = binarysearch(eachCellTransitNode[mappedSrc], low, high, mit->second.first);
				if (pos < 0){
					if (high >= eachCellTransitNode[mappedSrc][0]){
						eachCellTransitNode[mappedSrc] = (NodeID *)realloc(eachCellTransitNode[mappedSrc], (high + EACH_CELL_TRANSIT_NUM) * sizeof(NodeID));
						eachCellTransitNode[mappedSrc][0] = high + EACH_CELL_TRANSIT_NUM;								
					}
					insert(eachCellTransitNode[mappedSrc], -pos - 1 ,mit->second.first, high);
					eachCellTransitNode[mappedSrc][1]++;
				}

				low = 2;
				high =  eachCellTransitNode[mappedTrg][1];	
				pos = binarysearch(eachCellTransitNode[mappedTrg], low, high, mit->second.first);
				if (pos < 0){
					if (high >= eachCellTransitNode[mappedTrg][0]){
						eachCellTransitNode[mappedTrg] = (NodeID *)realloc(eachCellTransitNode[mappedTrg], (high + EACH_CELL_TRANSIT_NUM) * sizeof(NodeID));
						eachCellTransitNode[mappedTrg][0] = high + EACH_CELL_TRANSIT_NUM;								
					}
					insert(eachCellTransitNode[mappedTrg], -pos - 1 ,mit->second.first, high);
					eachCellTransitNode[mappedTrg][1]++;
				}
				mit++;
			}
			cout << "now transit size " << allCellTransitNode.size() << endl;
		}


		// process everty node to the corresponding cell
		void processEachNodeToCell(vector <CoordinateType> &xcord, vector <CoordinateType> &ycord)
		{
			NodeID  totalNodes = graph->noOfNodes();
			// allocate the space
			eachCellContainedNode = (NodeID **) malloc(sizeof(NodeID *) * allCellNum);
			for (NodeID i = 0; i < allCellNum ; i++){
				eachCellContainedNode[i] = (NodeID *) malloc(sizeof(NodeID) * avgNodesPerCell);
				eachCellContainedNode[i][0] = avgNodesPerCell;
				eachCellContainedNode[i][1] = 2;
			}
			Cell c;
			NodeID pos, mapu, len;
			// find every node to the corresponding cell
			for (NodeID u = 0; u < totalNodes; u++){
				getCell(u, c, xcord, ycord);
				mapu = c.xvertex * ascale + c.yvertex;
				pos = eachCellContainedNode[mapu][1];
				len = eachCellContainedNode[mapu][0];
				if (pos >= len){
					eachCellContainedNode[mapu] = (NodeID *)realloc(eachCellContainedNode[mapu], (len + avgNodesPerCell) * sizeof(NodeID));
					eachCellContainedNode[mapu][0] += avgNodesPerCell;
				}
				eachCellContainedNode[mapu][pos] = u;
				eachCellContainedNode[mapu][1]++;
			}

			/*fstream out("tnoutput/eachVertex_a=256.txt", ios::out);
			for (NodeID cellIndex = 0; cellIndex < allCellNum; cellIndex++){
				Cell c;
				c.xvertex = cellIndex / ascale;
				c.yvertex = cellIndex % ascale;
				out << " cell : (" << c.xvertex << " " << c.yvertex << ") 's containing node: "
					<< eachCellContainedNode[cellIndex][1] - 2 << endl;
				for (NodeID u = 2; u < eachCellContainedNode[cellIndex][1]; u++){
					out << eachCellContainedNode[cellIndex][u] << " ";
					if ((u - 1) % 100 == 0) out << endl;
				}
				out << endl;
				out << "transit node " << eachCellTransitNode[cellIndex][1] - 2 << endl;
				for (NodeID u = 2; u < eachCellTransitNode[cellIndex][1]; u++){
					out << eachCellTransitNode[cellIndex][u] << " ";
					if ((u - 1) % 10 == 0) out << endl;
				}
				out << endl;
			}
			out << endl;
			out.close();*/
			cout << "writing each cell containing node over " << endl; 
		}

		void buildDistanceTable()
		{
			//build the distance table of each transit node
			const bool performBucketScans = true;
			int earlyStopLevel = 10;
			
			set <NodeID>::iterator sit = allCellTransitNode.begin();
			//build the  transit node mapping
			int count = 0;
			while (sit != allCellTransitNode.end()){
				transitNodeMapping.insert(make_pair(*sit, count));
				count++;
				sit++;
			}
			//the internal position of each transit node for each cell
			eachCellTransitMapping = (NodeID **) malloc(sizeof(NodeID *) * allCellNum);
			for (unsigned int u = 0; u < allCellNum; u++){
				int size = eachCellTransitNode[u][0];

				eachCellTransitMapping[u] = (NodeID *) malloc(sizeof(NodeID) * size);
				eachCellTransitMapping[u][0] = eachCellTransitNode[u][0];
				eachCellTransitMapping[u][1] = eachCellTransitNode[u][1];
				for (int i = 2; i < size; i++){
					eachCellTransitMapping[u][i] = transitNodeMapping[eachCellTransitNode[u][i]];
				}
			}
			//build level one distance table
			NodeID allTransitNum = allCellTransitNode.size();
			transitDistTable.setRowAndCol(allTransitNum, allTransitNum);
			transitDistTable.init(0);
			ManyToMany<MyGraph, DijkstraManyToManyFW, DijkstraManyToManyBW, performBucketScans> mtm(graph, earlyStopLevel);
			vector<NodeID> allCellTransitCopy(allCellTransitNode.begin(), allCellTransitNode.end());
			mtm.computeMatrix(allCellTransitCopy, allCellTransitCopy, transitDistTable);		
			//build the each cell's distance table, from the cell's inner node to the access node
			eachCellDistTable = new Matrix <EdgeWeight>[allCellNum];
			for (NodeID cellIndex = 0; cellIndex < allCellNum; cellIndex++){
				ManyToMany<MyGraph, DijkstraManyToManyFW, DijkstraManyToManyBW, performBucketScans> mtminner(graph, earlyStopLevel);
				NodeID srcSize = eachCellContainedNode[cellIndex][1] - 2;
				NodeID transitSize = eachCellTransitNode[cellIndex][1] - 2;
				vector <NodeID> srcNode(srcSize);
				vector <NodeID> transitNode(transitSize);
				for (NodeID u = 0; u < srcSize; u++){
					srcNode[u] = eachCellContainedNode[cellIndex][u + 2];						
				}
				for (NodeID u = 0; u < transitSize; u++){
					transitNode[u] = eachCellTransitNode[cellIndex][u + 2];						
				}
				//initialize this cell's distancetable
				eachCellDistTable[cellIndex].setRowAndCol(srcSize, transitSize);
				eachCellDistTable[cellIndex].init(0);
				// compute the matrix
				mtminner.computeMatrix(srcNode, transitNode, eachCellDistTable[cellIndex]);
			}
			cout << "build each cell distance table finished" << endl;
		}
		
		
		void getThreeGridAwayTransitNodes(int cellX, int cellY, set <NodeID> &localTransit)
		{
			int x1, y1, x2, y2;
			x1 = ( ((cellX - 3) > 0) ? (cellX - 3) : 0 );
			y1 = ( ((cellY - 3) > 0) ? (cellY - 3) : 0 );
			x2 = ( ((cellX + 4) < ascale) ? (cellX + 4) : ascale);
			y2 = ( ((cellY + 4) < ascale) ? (cellY + 4) : ascale);
			int i, j, k;
			int mappedValue;
			for (i = x1; i < x2; i++){
				for (j = y1; j < y2; j++){
					//cout << " tackle cell " << i << " * " << j << " ";
					mappedValue = i * ascale + j;
					//cout << eachCellTransitNode[mappedValue][1] - 2 << endl;
					for (k = 2; k < eachCellTransitNode[mappedValue][1]; k++){
						localTransit.insert(eachCellTransitNode[mappedValue][k]);
					}
				}
			}
			/*set <NodeID>::iterator sit = localTransit.begin();
			while (sit != localTransit.end()){
				cout << *sit << " ";
				sit++;
			}
			cout << endl;*/
		}
		/* using floyd-warshal algorithm for all pairs shortest algorithm */
		void buildDistanceTableByWarshal()
		{
			int i, j, k, mappedValue;
			vector <CompleteEdge> edgeList, combinedEdgeList;
			edgeList.reserve(graph->noOfEdges());
			set <NodeID> localTransit;
			set <NodeID>::iterator sit;
			int cnt;
			vector <NodeID> localSrc, localTrg, localContained;
			int localSize;
			Matrix <EdgeWeight> ma;

			eachCellDistTable = new Matrix <EdgeWeight>[allCellNum];
			//for (NodeID cellIndex = 0; cellIndex < allCellNum; cellIndex++){
			//	ManyToMany<MyGraph, DijkstraManyToManyFW, DijkstraManyToManyBW, performBucketScans> mtminner(graph, earlyStopLevel);
			//	NodeID srcSize = eachCellContainedNode[cellIndex][1] - 2;
			//	NodeID transitSize = eachCellTransitNode[cellIndex][1] - 2;
			//	vector <NodeID> srcNode(srcSize);
			//	vector <NodeID> transitNode(transitSize);
			//	for (NodeID u = 0; u < srcSize; u++){
			//		srcNode[u] = eachCellContainedNode[cellIndex][u + 2];						
			//	}
			//	for (NodeID u = 0; u < transitSize; u++){
			//		transitNode[u] = eachCellTransitNode[cellIndex][u + 2];						
			//	}
			//	//initialize this cell's distancetable
			//	eachCellDistTable[cellIndex].setRowAndCol(srcSize, transitSize);
			//	eachCellDistTable[cellIndex].init(0);
			//	// compute the matrix
			//	mtminner.computeMatrix(srcNode, transitNode, eachCellDistTable[cellIndex]);
			//}
			for (i = 0; i < ascale; i++){
				for (j = 0; j < ascale; j++){
					mappedValue = i * ascale + j;
					if (mappedValue % 100 == 0)
						cout << "transit graph: have computed cell " << mappedValue << endl;
					localTransit.clear();
					getThreeGridAwayTransitNodes(i, j, localTransit);
					localSize = eachCellTransitNode[mappedValue][1] - 2;
					//cout << "local size " << localSize << " " << localTransit.size() << endl;
					if (localSize > 0){
						ma.clear();
						ma.setRowAndCol(localSize, localTransit.size());
						ma.init(0);
						//cout << "init matrix end" << endl;
						localSrc.resize(localSize);
						for (NodeID u = 0; u < localSize; u++){						
							localSrc[u] = eachCellTransitNode[mappedValue][u + 2];
							//cout << localSrc[u] << " ";
						}
						localTrg.resize(localTransit.size());
						cnt = 0;
						sit = localTransit.begin();
						while (sit != localTransit.end()){
							localTrg[cnt] = *sit;
							sit++;
							cnt++;
						}
						
						NodeID srcSize = eachCellContainedNode[mappedValue][1] - 2;
						//cout << " source size " << srcSize << endl;
						localContained.resize(srcSize);
						for (NodeID u = 0; u < srcSize; u++){						
							localContained[u] = eachCellContainedNode[mappedValue][u + 2];
							//cout << localContained[u] << " ";
						}
						//cout << endl;
						eachCellDistTable[mappedValue].setRowAndCol(srcSize, localSize);
						eachCellDistTable[mappedValue].init(0);
						localManyToMany_DistTable(localSrc, localTrg, ma, localContained, eachCellDistTable[mappedValue]);
										
						for (NodeID u = 0; u < localSrc.size(); u++){
							for (NodeID v = 0; v < localTrg.size(); v++){
								CompleteEdge edge1(transitNodeMapping[localSrc[u]], transitNodeMapping[localTrg[v]], ma.value(u, v), false, true, false);
								CompleteEdge edge2(transitNodeMapping[localTrg[v]], transitNodeMapping[localSrc[u]], ma.value(u, v), false, false, true);
								edgeList.push_back(edge1);
								edgeList.push_back(edge2);
							}
						}
					}
				}
			}
			sort(edgeList.begin(), edgeList.end());
			removeDuplicates(edgeList, combinedEdgeList);

			/*fstream combinedFile("tnoutput/combinedEdge.ch", ios::out | ios::binary);
			VectorSerializer<CompleteEdge, unsigned int, ComplexSerializer<CompleteEdge>>::serialize(combinedFile,combinedEdgeList);
			combinedFile.close();*/

			/*fstream combinedFile("tnoutput/combinedEdge.ch", ios::in | ios::binary);
			VectorSerializer<CompleteEdge, unsigned int, ComplexSerializer<CompleteEdge>>::deserialize(combinedFile,combinedEdgeList);
			combinedFile.close();*/

			int nodeNum = allCellTransitNode.size();
			//// standard all pair shortest algorithm
			//Matrix <EdgeWeight> trDistRef(allCellTransitNode.size(), allCellTransitNode.size());
			//trDistRef.init(Weight::MAX_VALUE);
			//fstream edgeFile("tnoutput/cEdge.txt", ios::out);
			//for (i = 0; i < combinedEdgeList.size(); i++){
			//	//edgeFile << combinedEdgeList[i].source() << " " << combinedEdgeList[i].target() << " " << combinedEdgeList[i].weight() << endl;
			//	trDistRef.set(combinedEdgeList[i].source(), combinedEdgeList[i].target(), combinedEdgeList[i].weight());
			//}
			//edgeFile.close();
			//cout << "begin floyed algorith" << endl;
			//for (i = 0; i < nodeNum; i++)
			//	trDistRef.set(i, i, 0);
			//for (k = 0; k < nodeNum; k++){
			//	for (i = 0; i < nodeNum; i++){
			//		for (j = 0; j < nodeNum; j++){
			//			//cout  << i << " " << j << " " << trDistRef.value(i, k)  << " "<< trDistRef.value(k, j) << " " << trDistRef.value(i, j) << endl;
			//			if (trDistRef.value(i, k) != Weight::MAX_VALUE && trDistRef.value(k, j) != Weight::MAX_VALUE){
			//				if (trDistRef.value(i, k) + trDistRef.value(k, j) < trDistRef.value(i, j))
			//					trDistRef.set(i, j, trDistRef.value(i, k) + trDistRef.value(k, j));
			//			}
			//		}
			//	}
			//	if (k % 10 == 0)
			//		cout << k << endl;
			//}

			vector <LevelID> nodeLevels(nodeNum, nodeNum);
			datastr::graph::UpdateableGraph *trUpd = new datastr::graph::UpdateableGraph(combinedEdgeList, nodeLevels);
			cout << "new transit node graph has been built " << endl;
			//Matrix <EdgeWeight> trDistRef(allCellTransitNode.size(), allCellTransitNode.size());
			transitDistTable.setRowAndCol(allCellTransitNode.size(), allCellTransitNode.size());
			transitDistTable.init(0);
			

			processing::DijkstraCH<datastr::graph::UpdateableGraph, NormalPQueue, 2, false> dijkstraCHTest(trUpd);
			for (i = 0; i < trUpd->noOfNodes(); i++){
				dijkstraCHTest.searchWithoutTarget(i);
				for(j = 0; j < trUpd->noOfNodes(); j++){
					if (i == j)
						transitDistTable.set(i, j, 0);
					else
						transitDistTable.set(i, j, dijkstraCHTest.distanceTo(j));
				}
				dijkstraCHTest.clear();
				if ( i % 100 == 0)
					cout << i << endl;
			}
			cout << " distance tabel finished " << endl;
			
		}

		void printEdgeBucket(fstream &out)
		{
			NodeID edgetotal = 2 * ascale * (ascale + 1);
			//output the horizontal line
			bool firstTime = true;
			SpaceVertex beg, end;
			beg.set(0, 0);
			end.set(ascale, 0);
			for (NodeID i = 0; i < ascale + 1; i++){
				printEdgeBucketByLine(out, beg, end, firstTime);	
				beg.increseYPosition();
				end.increseYPosition();
			}

			firstTime = false;
			beg.set(0, 0);
			end.set(0, ascale);
			for (NodeID i = 0; i < ascale + 1; i++){
				printEdgeBucketByLine(out, beg, end, firstTime);
				beg.increseXPosition();
				end.increseXPosition();
			}
		}
		void printEdgeBucketByLine(fstream &out,SpaceVertex startV, SpaceVertex endV, bool firstTime)
		{
			SpaceVertex currentNode, nextNode;
			NodeID startPosition, endPosition;
			BorderEdge be;
			/*NodeID buckctStart;
			if (firstTime)
				buckctStart = 0;
			else
				buckctStart = ascale * (ascale + 1);*/
			bool horizon = false;
			if (startV.getYPosition() == endV.getYPosition()){
				horizon = true;
			}
			if (horizon){
				startPosition = startV.getXPosition();
				endPosition = endV.getXPosition();
			}
			else{
				startPosition = startV.getYPosition();
				endPosition = endV.getYPosition();
			}

			currentNode.set(startV.getXPosition(), startV.getYPosition());			
			//NodeID mappedBucket;
			for (NodeID i = startPosition; i < endPosition; i++){
				if (horizon){
					nextNode.set(currentNode.getXPosition() + 1, currentNode.getYPosition());
					//mappedBucket = buckctStart + currentNode.getYPosition() * ascale + currentNode.getXPosition();
				}		
				else{
					nextNode.set(currentNode.getXPosition(), currentNode.getYPosition() + 1);
					//mappedBucket = buckctStart + currentNode.getXPosition() * ascale + currentNode.getYPosition();
				}
				be.setBorderEdge(currentNode, nextNode);
				vector <EasyEdge> b(*beHashMap[be.getValue()]);

				out << currentNode << " to " << nextNode << endl;
				out << "bucket size : "<< b.size() << endl;
				for (NodeID s = 0; s < b.size(); s++)
					out  << s << " (" << b[s].source << " ---> " << b[s].target << " ) ";
				out << endl;
				
				/*pair <NodeID, vector <EasyEdge> *> bucketPair(be.getValue(), &edgeBucket[mappedBucket]);
				beHashMap.insert(bucketPair);*/

				if (horizon)
					currentNode.increseXPosition();
				else
					currentNode.increseYPosition();
			}
		}

		void printOneCell(fstream &out)
		{
			SpaceVertex currentNode, nextNode;
			currentNode.set(102,127);
			nextNode.set(102,128);
			BorderEdge be;
			be.setBorderEdge(currentNode, nextNode);
			vector <EasyEdge> b(*beHashMap[be.getValue()]);

			out << currentNode << " to " << nextNode << endl;
			out << "bucket size : "<< b.size() << endl;
			for (NodeID s = 0; s < b.size(); s++){			
				out  << s << " (" << b[s].source << " ---> " << b[s].target << " ) ";
				if (s % 10 == 0) out << endl;
			}
			out << endl;

		}
private:
		//  asssociate all edges of the line from startV to endV  with an EdgeBucket
		void insertEdgeBucketByLine(SpaceVertex startV, SpaceVertex endV, bool firstTime)
			{
				SpaceVertex currentNode, nextNode;
				NodeID startPosition, endPosition;
				BorderEdge be;
				NodeID buckctStart;
				if (firstTime)
					buckctStart = 0;
				else
					buckctStart = ascale * (ascale + 1);
				bool horizon = false;
				if (startV.getYPosition() == endV.getYPosition()){
					horizon = true;
				}
				if (horizon){
					startPosition = startV.getXPosition();
					endPosition = endV.getXPosition();
				}
				else{
					startPosition = startV.getYPosition();
					endPosition = endV.getYPosition();
				}

				currentNode.set(startV.getXPosition(), startV.getYPosition());			
				NodeID mappedBucket;
				for (NodeID i = startPosition; i < endPosition; i++){
					if (horizon){
						nextNode.set(currentNode.getXPosition() + 1, currentNode.getYPosition());
						mappedBucket = buckctStart + currentNode.getYPosition() * ascale + currentNode.getXPosition();
					}		
					else{
						nextNode.set(currentNode.getXPosition(), currentNode.getYPosition() + 1);
						mappedBucket = buckctStart + currentNode.getXPosition() * ascale + currentNode.getYPosition();
					}
					be.setBorderEdge(currentNode, nextNode);
					vector <EasyEdge> b;
#ifdef DEBUG
					//cout << be.getValue() << endl;
#endif
					pair <int, vector <EasyEdge> *> bucketPair(be.getValue(), &edgeBucket[mappedBucket]);
					beHashMap.insert(bucketPair);

					if (horizon)
						currentNode.increseXPosition();
					else
						currentNode.increseYPosition();
				}
			}


//		// insert an edge to the  border edge Hash Map, when the edge intersects with one border edge, it is 
//		// associated with this border edge's hash map
//		void insertToBorderEdgeHashMap(NodeID u, NodeID v, vector <CoordinateType> &xcord, vector <CoordinateType> &ycord, fstream &out)
//		{
//			
//			Cell startCell, endCell, temp;
//			getCell(u, startCell,xcord, ycord);
//			getCell(v, endCell, xcord, ycord);
//
//
//			//always assure startcell is in the left of endcell
//			NodeID t;
//			if (xcord[u] > xcord[v]){
//				// swap startcell, endcell
//				temp = startCell;
//				startCell = endCell;
//				endCell = temp;
//				// swap u, v
//				t = u;
//				u = v;
//				v = t;
//			}
//
////#ifdef DEBUG
////			out << "edge: "<< u << " " << v << endl;
////			out << "cell position " << endl;
////			out << startCell.xvertex << " " << startCell.yvertex << " ";
////			out << endCell.xvertex << " " << endCell.yvertex << endl;
////#endif
//
//			NodeID nextXBorder, nextYBorder;
//			CoordinateType nextXBorderCor, nextYBorderCor;
//			CoordinateType nextX, nextY;
//			SpaceVertex startPo, endPo;
//			// calulate the line 's ration and b--(y = k * x + b)
//			// we should consider the line is vertical with the coordinate line
//			if((xcord[v] == xcord[u]) || (ycord[v] == ycord[u])){
//				// this is the vertical line
//				if (xcord[v] == xcord[u]){
////#ifdef DEBUG
////					out << "this is vertical line " << endl;
////#endif
//					if (ycord[u] < ycord[v]){
//						nextXBorder = startCell.xvertex + 1;
//						nextYBorder = startCell.yvertex + 1;
//					}
//					else{
//						nextXBorder = startCell.xvertex + 1;
//						nextYBorder = startCell.yvertex;
//					}
//					while (true){
//						nextYBorderCor = nextYBorder * cellLength;
//						if(nextYBorder < 0 || nextYBorder > ascale) break;
//						if (ycord[u] < ycord[v]){
//							if (ycord[v] < nextYBorderCor) break;
//						}
//						else{
//							if (ycord[v] > nextYBorderCor) break;
//						}
//						endPo.set(nextXBorder, nextYBorder);
//						startPo.set(nextXBorder - 1, nextYBorder);
//						if (ycord[u] < ycord[v]){
//							nextYBorder++;			
//						}
//						else{
//							nextYBorder--;			
//						}
//						if (startPo.x > endPo.x || (startPo.x == endPo.x && startPo.y > endPo.y)){
//							SpaceVertex sv = startPo;
//							startPo = endPo;
//							endPo = sv;
//						}
//						
//						BorderEdge be(startPo, endPo);
//						(*beHashMap[be.getValue()]).push_back(EasyEdge(u, v));
//		//#ifdef DEBUG
//		//				/*cout << be.start.getXPosition() << " " << be.start.getYPosition() << endl;
//		//				cout << be.end.getXPosition() << " " << be.end.getYPosition() << endl;*/
//		//				out << be.getValue() << endl;
//		//				out << nextXBorder << " " << nextYBorder << endl;
//		//#endif
//					}
//				}
//				//horizonal line
//				if (ycord[u] == ycord[v]){
////#ifdef DEBUG
////					out << "this is horizon line" << endl;
////#endif
//					nextXBorder = startCell.xvertex + 1;
//					nextYBorder = startCell.yvertex + 1;
//					while (true){
//						nextXBorderCor = nextXBorder * cellLength;
//						if(nextXBorder < 0 || nextXBorder > ascale) break;
//						if (xcord[v] < nextXBorderCor) break;
//						endPo.set(nextXBorder, nextYBorder);
//						startPo.set(nextXBorder, nextYBorder - 1);
//						nextXBorder++;
//						if (startPo.x > endPo.x || (startPo.x == endPo.x && startPo.y > endPo.y)){
//							SpaceVertex sv = startPo;
//							startPo = endPo;
//							endPo = sv;
//						}
//						
//						BorderEdge be(startPo, endPo);
//						(*beHashMap[be.getValue()]).push_back(EasyEdge(u, v));
//		//#ifdef DEBUG
//		//				/*cout << be.start.getXPosition() << " " << be.start.getYPosition() << endl;
//		//				cout << be.end.getXPosition() << " " << be.end.getYPosition() << endl;*/
//		//				out << be.getValue() << endl;
//		//				out << nextXBorder << " " << nextYBorder << endl;
//		//#endif
//						
//					}
//					
//				}
//			}
//			else{
//				CoordinateType ratio = (ycord[v] - ycord[u]) / (xcord[v] - xcord[u]);
//				CoordinateType b = ycord[v] - ratio * xcord[v];
//
//	/*#ifdef DEBUG
//				out << " ratio = "<< ratio << " " << b << endl;
//	#endif*/
//	
//				if (ycord[u] < ycord[v]){
//					nextXBorder = startCell.xvertex + 1;
//					nextYBorder = startCell.yvertex + 1;
//				}
//				else{
//					nextXBorder = startCell.xvertex + 1;
//					nextYBorder = startCell.yvertex;
//				}
//	/*#ifdef DEBUG
//				out << nextXBorder << " " << nextYBorder << endl;
//	#endif
//		*/	
//				while (true){
//					nextXBorderCor = nextXBorder * cellLength;
//					nextYBorderCor = nextYBorder * cellLength;
//					//cout << nextXBorderCor << " " << nextYBorderCor << endl;
//					endPo.set(nextXBorder, nextYBorder);
//					if (nextXBorder < 0 || nextXBorder > ascale || nextYBorder < 0 || nextYBorder > ascale)
//						break;
//					if (ycord[u] < ycord[v]){
//						if (xcord[v] < nextXBorderCor && ycord[v] < nextYBorderCor) break;
//					}
//					else{
//						if (xcord[v] < nextXBorderCor && ycord[v] > nextYBorderCor) {					
//							break;
//						}
//					}
//					nextY = ratio * nextXBorder * cellLength + b;
//					nextX = (nextYBorder * cellLength - b) / ratio;
//					// check the top case
//					if (ycord[u] < ycord[v]){
//						// when the vertex is in the inter
//						if (nextY >= nextYBorderCor){
//							startPo.set(nextXBorder - 1, nextYBorder);
//							nextYBorder++;
//						}
//						else{
//							startPo.set(nextXBorder, nextYBorder - 1);
//							nextXBorder++;
//						}
//
//					}
//					// check the bottom case
//					else{
//						if (nextY >= nextYBorderCor){
//							startPo.set(nextXBorder, nextYBorder + 1);
//							nextXBorder++;
//						}
//						else{
//							startPo.set(nextXBorder - 1, nextYBorder);
//							nextYBorder--;
//						}
//					}
//
//					// assure startPo is in the left or in the bottom of the endPo
//					if (startPo.x > endPo.x || (startPo.x == endPo.x && startPo.y > endPo.y)){
//						SpaceVertex sv = startPo;
//						startPo = endPo;
//						endPo = sv;
//					}
//					
//					BorderEdge be(startPo, endPo);
//					(*beHashMap[be.getValue()]).push_back(EasyEdge(u, v));
//	//#ifdef DEBUG
//	//				/*cout << be.start.getXPosition() << " " << be.start.getYPosition() << endl;
//	//				cout << be.end.getXPosition() << " " << be.end.getYPosition() << endl;*/
//	//				out << be.getValue() << endl;
//	//				out << nextXBorder << " " << nextYBorder << endl;
//	//#endif
//					}
//				}
////#ifdef DEBUG
////			out << endl;
////#endif
//
//			}

		void insertToBorderEdgeHashMap(NodeID u, NodeID v, vector <CoordinateType> &xcord, vector <CoordinateType> &ycord, fstream &out)
		{
			
			Cell startCell, endCell, temp;
			getCell(u, startCell,xcord, ycord);
			getCell(v, endCell, xcord, ycord);


			//always assure startcell is in the left of endcell
			NodeID t;
			if (xcord[u] > xcord[v]){
				// swap startcell, endcell
				temp = startCell;
				startCell = endCell;
				endCell = temp;
				// swap u, v
				t = u;
				u = v;
				v = t;
			}

			NodeID nextXBorder, nextYBorder;
			CoordinateType nextXBorderCor, nextYBorderCor;
			CoordinateType nextX, nextY;
			SpaceVertex startPo, endPo;
			// calulate the line 's ration and b--(y = k * x + b)
			// we should consider the line is vertical with the coordinate line
			if((xcord[v] == xcord[u]) || (ycord[v] == ycord[u])){
				// this is the vertical line
				if (xcord[v] == xcord[u]){
					if (ycord[u] < ycord[v]){
						nextXBorder = startCell.xvertex + 1;
						nextYBorder = startCell.yvertex + 1;
					}
					else{
						nextXBorder = startCell.xvertex + 1;
						nextYBorder = startCell.yvertex;
					}
					while (true){
						if(nextYBorder < 0 || nextYBorder > ascale) break;
						nextYBorderCor = nextYBorder * cellLength;
						if (ycord[u] < ycord[v]){
							if (ycord[v] < nextYBorderCor) break;
						}
						else{
							if (ycord[v] > nextYBorderCor) break;
						}
						endPo.set(nextXBorder, nextYBorder);
						startPo.set(nextXBorder - 1, nextYBorder);
						if (ycord[u] < ycord[v]){
							nextYBorder++;			
						}
						else{
							nextYBorder--;			
						}
						//if (startPo.x > endPo.x || (startPo.x == endPo.x && startPo.y > endPo.y)){
						//	SpaceVertex sv = startPo;
						//	startPo = endPo;
						//	endPo = sv;
						//}
						
						BorderEdge be(startPo, endPo);
						(*beHashMap[be.getValue()]).push_back(EasyEdge(u, v));
					}
				}
				//horizonal line
				if (ycord[u] == ycord[v]){
					nextXBorder = startCell.xvertex + 1;
					nextYBorder = startCell.yvertex + 1;
					while (true){
						if(nextXBorder < 0 || nextXBorder > ascale) break;
						nextXBorderCor = nextXBorder * cellLength;
						if (xcord[v] < nextXBorderCor) break;
						endPo.set(nextXBorder, nextYBorder);
						startPo.set(nextXBorder, nextYBorder - 1);
						nextXBorder++;
						//if (startPo.x > endPo.x || (startPo.x == endPo.x && startPo.y > endPo.y)){
						//	SpaceVertex sv = startPo;
						//	startPo = endPo;
						//	endPo = sv;
						//}
						
						BorderEdge be(startPo, endPo);
						(*beHashMap[be.getValue()]).push_back(EasyEdge(u, v));
						
					}
					
				}
			}
			else{
				CoordinateType ratio = (ycord[v] - ycord[u]) / (xcord[v] - xcord[u]);
				CoordinateType b = ycord[v] - ratio * xcord[v];
	
				if (ycord[u] < ycord[v]){
					nextXBorder = startCell.xvertex + 1;
					nextYBorder = startCell.yvertex + 1;
				}
				else{
					nextXBorder = startCell.xvertex + 1;
					nextYBorder = startCell.yvertex;
				}
				while (true){
					if (nextXBorder < 0 || nextXBorder > ascale || nextYBorder < 0 || nextYBorder > ascale)
						break;
					nextXBorderCor = nextXBorder * cellLength;
					nextYBorderCor = nextYBorder * cellLength;
					//cout << nextXBorderCor << " " << nextYBorderCor << endl;
					//endPo.set(nextXBorder, nextYBorder);
					if (ycord[u] < ycord[v]){
						if (xcord[v] < nextXBorderCor && ycord[v] < nextYBorderCor) break;
					}
					else{
						if (xcord[v] < nextXBorderCor && ycord[v] > nextYBorderCor) break;
					}
					nextY = ratio * nextXBorder * cellLength + b;
					nextX = (nextYBorder * cellLength - b) / ratio;
					// check the top case
					if (ycord[u] < ycord[v]){
						endPo.set(nextXBorder, nextYBorder);
						// when the vertex is in the inter
						if (nextY >= nextYBorderCor){
							startPo.set(nextXBorder - 1, nextYBorder);
							nextYBorder++;
						}
						else{
							startPo.set(nextXBorder, nextYBorder - 1);
							nextXBorder++;
						}

					}
					// check the bottom case
					else{
						if (nextY >= nextYBorderCor){
							//startPo.set(nextXBorder, nextYBorder + 1);
							startPo.set(nextXBorder, nextYBorder);
							endPo.set(nextXBorder, nextYBorder + 1);
							nextXBorder++;
						}
						else{
							startPo.set(nextXBorder - 1, nextYBorder);
							endPo.set(nextXBorder, nextYBorder);
							nextYBorder--;
						}
					}

					//// assure startPo is in the left or in the bottom of the endPo
					//if (startPo.x > endPo.x || (startPo.x == endPo.x && startPo.y > endPo.y)){
					//	SpaceVertex sv = startPo;
					//	startPo = endPo;
					//	endPo = sv;
					//}
					
					BorderEdge be(startPo, endPo);
					(*beHashMap[be.getValue()]).push_back(EasyEdge(u, v));
					}
				}
			}
	

		// get the start vertex's own cell's bounding box and get the access node
		void getStartVertexAccessNode(Cell &cell, vector <NodeID> &startAccessNode, 
			vector <CoordinateType> &xcord, vector <CoordinateType> &ycord, fstream &out)
			{
//#ifdef DEBUG
//				out << "start access node " << endl;
//#endif
				// the cell's bounding box's length is one
				getAccessNode(cell, startAccessNode, 1, xcord, ycord, out);
//#ifdef DEBUG
//				out << "start size: " << startAccessNode.size() << endl;
//				for (NodeID u = 0; u < startAccessNode.size(); u++)
//					out << startAccessNode[u] << " " ;
//				out << endl;
//#endif
			}

		// get the cell's inner bounding box's access node, if we want to adjust a, b
		// we can set the inner access box's cell length
		void getInnerAccessNode(Cell &cell, vector <NodeID> &innerAccessNode,
			vector <CoordinateType> &xcord, vector <CoordinateType> &ycord, fstream &out)
			{
//#ifdef DEBUG
//				out << "inner node" << endl;
//#endif
				getAccessNode(cell, innerAccessNode, 5, xcord, ycord, out);
//#ifdef DEBUG
//				out << "inner size: " << innerAccessNode.size() << endl;
//				for (NodeID u = 0; u < innerAccessNode.size(); u++)
//					out << innerAccessNode[u] << " " ;
//				out << endl;
//#endif
			}

		// get the cell's outer bounding box's access node, if we want to adjust a, b
		// we need set the outer access bounding box's cell length
		void getOuterAccessNode(Cell &cell, vector <NodeID> &outerAccessNode,
			vector <CoordinateType> &xcord, vector <CoordinateType> &ycord, fstream &out)
			{
//#ifdef DEBUG
//				out << "outer node" << endl;
//#endif
				getAccessNode(cell, outerAccessNode, 9, xcord, ycord, out);
//#ifdef DEBUG
//				out << "outer size: " << outerAccessNode.size() << endl;
//				for (NodeID u = 0; u < outerAccessNode.size(); u++)
//					out << outerAccessNode[u] << " " ;
//				out << endl;
//#endif
			}

		// given a cell, and the cell's length, get the access node of this cell
		void getAccessNode(Cell &cell, vector <NodeID> &accessNode, NodeID cellnum,
			vector <CoordinateType> &xcord, vector <CoordinateType> &ycord, fstream &out)
			{
				SpaceVertex lb, rb, rt, lt;
				SpaceVertex tStart, tEnd;
				SpaceVertex currentNode, nextNode;
				BorderEdge be;
				Rectangle rec;
				lb.set(cell.xvertex, cell.yvertex);
				rb.set(cell.xvertex + cellnum, cell.yvertex);
				rt.set(cell.xvertex + cellnum, cell.yvertex + cellnum);
				lt.set(cell.xvertex, cell.yvertex + cellnum);
				//cout << lb << " " << rb << " " << rt << " " << lt << " " << endl;

				if (lb.getXPosition() < 0){
					rec.xleft = 0;
				}
				else{
					rec.xleft = lb.getXPosition() * cellLength;
				}
				if (lb.getYPosition() < 0)
					rec.yleft = 0;
				else{
					rec.yleft = lb.getYPosition() * cellLength;
				}

				if (rt.getXPosition() < (int)ascale)
					rec.xright = rt.getXPosition() * cellLength;
				else{
					rec.xright = ascale * cellLength;
				}
				if (rt.getYPosition() < (int)ascale)
					rec.yright = rt.getYPosition() * cellLength;
				else{
					rec.yright = ascale * cellLength;
				}

				// get the bottom edge
				NodeID ybottom = lb.getYPosition();
				if (ybottom >=0 && ybottom <= ascale){
					tStart = lb;
					tEnd = rb;
					if (lb.getXPosition() < 0) tStart.setXPosition(0);
					if (rb.getXPosition() > (int)ascale) tEnd.setXPosition(ascale);
					getAccessNodeByLine(tStart, tEnd, rec, accessNode, xcord, ycord, out);
				}
				// get the top edge
				NodeID ytop = lt.getYPosition();
				if (ytop >=0 && ytop <= ascale){
					tStart = lt;
					tEnd = rt;
					//cout << lt.getXPosition() << endl;
					//cout << tStart << endl;
					if (lt.getXPosition() < 0) tStart.setXPosition(0);
					if (rt.getXPosition() > (int)ascale) tEnd.setXPosition(ascale);
					getAccessNodeByLine(tStart, tEnd, rec, accessNode, xcord, ycord, out);
				}

				// get the left edge
				NodeID xleft = lb.getXPosition();
				if (xleft >=0 && xleft <= ascale){
					tStart = lb;
					tEnd = lt;
					if (lb.getYPosition() < 0) tStart.setYPosition(0);
					if (lt.getYPosition() > (int)ascale) tEnd.setYPosition(ascale);
					getAccessNodeByLine(tStart, tEnd, rec, accessNode, xcord, ycord, out);
				}

				NodeID xright = rb.getXPosition();
				if (xright >=0 && xright <= ascale){
					tStart = rb;
					tEnd = rt;
					if (rb.getYPosition() < 0) tStart.setYPosition(0);
					if (rt.getYPosition() > (int)ascale) tEnd.setYPosition(ascale);
					//cout << rb.getYPosition() << " "<< ascale << endl;
					//cout << "last tEnd : "<< tEnd << endl;
					getAccessNodeByLine(tStart, tEnd, rec, accessNode, xcord, ycord,out);
				}

			}

		// given the start vertex and the end vertex, get the access node of this line, used by 
		// the function getAccessNode()
		void getAccessNodeByLine(SpaceVertex startV, SpaceVertex endV, 
			Rectangle &rec, vector <NodeID> &accessNode,
			vector <CoordinateType> &xcord, vector <CoordinateType> &ycord, fstream &out)
			{
//#ifdef DEBUG
//				out << startV << " " << endV << endl;
//#endif
				SpaceVertex currentNode, nextNode;
				NodeID startPosition, endPosition;
				BorderEdge be;
				bool horizon = false;
				if (startV.getYPosition() == endV.getYPosition()){
					horizon = true;
				}
				if (horizon){
					startPosition = startV.getXPosition();
					endPosition = endV.getXPosition();
				}
				else{
					startPosition = startV.getYPosition();
					endPosition = endV.getYPosition();
				}

				currentNode.set(startV.getXPosition(), startV.getYPosition());			
				for (NodeID i = startPosition; i < endPosition; i++){
					if (horizon){
						nextNode.set(currentNode.getXPosition() + 1, currentNode.getYPosition());
					}		
					else{
						nextNode.set(currentNode.getXPosition(), currentNode.getYPosition() + 1);
					}
					be.setBorderEdge(currentNode, nextNode);
					vector <EasyEdge> copy(*beHashMap[be.getValue()]);
					for(NodeID u = 0; u < copy.size(); u++){
						EasyEdge ee = copy[u];
//#ifdef DEBUG
//						out << ee.source << " " << ee.target << endl;
//#endif
						if (isEdgeInRectangle(ee, rec, xcord, ycord))
							accessNode.push_back(ee.minid());
					}

					if (horizon)
						currentNode.increseXPosition();
					else
						currentNode.increseYPosition();
				}
			}
		
// 		void getAccessNodeByLine(SpaceVertex startV, SpaceVertex endV, 
// 			vector <NodeID> &accessNode,
// 			vector <CoordinateType> &xcord, vector <CoordinateType> &ycord)
// 		{
// 
// 			SpaceVertex currentNode, nextNode;
// 			NodeID startPosition, endPosition;
// 			BorderEdge be;
// 			bool horizon = false;
// 			if (startV.getYPosition() == endV.getYPosition()){
// 				horizon = true;
// 			}
// 			if (horizon){
// 				startPosition = startV.getXPosition();
// 				endPosition = endV.getXPosition();
// 			}
// 			else{
// 				startPosition = startV.getYPosition();
// 				endPosition = endV.getYPosition();
// 			}
// 
// 			currentNode.set(startV.getXPosition(), startV.getYPosition());			
// 			for (NodeID i = startPosition; i < endPosition; i++){
// 				if (horizon){
// 					nextNode.set(currentNode.getXPosition() + 1, currentNode.getYPosition());
// 				}		
// 				else{
// 					nextNode.set(currentNode.getXPosition(), currentNode.getYPosition() + 1);
// 				}
// 				be.setBorderEdge(currentNode, nextNode);
// 				vector <EasyEdge> copy(*beHashMap[be.getValue()]);
// 				for(NodeID u = 0; u < copy.size(); u++){
// 					EasyEdge ee = copy[u];
// 
// 					if (isEdgeInRectangle(ee, rec, xcord, ycord))
// 					accessNode.push_back(ee.minid());
// 				}
// 
// 				if (horizon)
// 					currentNode.increseXPosition();
// 				else
// 					currentNode.increseYPosition();
// 			}
// 		}
		void removeDupAccessNodes(vector <NodeID> &node)
		{
			sort(node.begin(), node.end());
			node.erase( unique( node.begin(), node.end() ), node.end() );	
		}


		int binarysearch(NodeID *&arrays, int l, int h, NodeID v){
			int cur=l;
			int low=l;
			int high=h;
			while(low<high){
				cur=low+((high-low)/2);
				if(arrays[cur]==v)
					return cur;
				else{
					if(arrays[cur]>v)
						high=cur;
					else
						low=cur+1;
				}
			}
			return -low-1;
		}
		//insert v into arrays, make sure the index of arrays is not out of range
		void insert(NodeID *&arrays, int pos, int v, int length){
			int k=0;
			for(k=length;k>pos;k--)
				arrays[k]=arrays[k-1];
			arrays[pos]=v;
		}
		void localManyToMany(vector <NodeID> &src, vector <NodeID> &trg, Matrix <EdgeWeight> & ma)
		{
			processing::DijkstraCH<datastr::graph::UpdateableGraph, NormalPQueue, 2, false> dijkstraCHTest(upd);
			for (NodeID u = 0; u < src.size(); u++){
				//cout  << " source " << src[u] << endl;
				dijkstraCHTest.localQueryForTransitNode(src[u], trg);
				for (NodeID j = 0; j < trg.size(); j++){
					//cout << trg[j] << " " << dijkstraCHTest.distanceTo(trg[j]) << endl;
					ma.set(u, j, dijkstraCHTest.distanceTo(trg[j]));
				}
				dijkstraCHTest.clear();
			}
		}
		void localManyToMany_DistTable(vector <NodeID> &src, vector <NodeID> &trg, Matrix <EdgeWeight> & ma, vector <NodeID> &localContained, Matrix <EdgeWeight> & distTable)
		{
			processing::DijkstraCH<datastr::graph::UpdateableGraph, NormalPQueue, 2, false> dijkstraCHTest(upd);
			for (NodeID u = 0; u < src.size(); u++){
				//cout  << " source " << src[u] << endl;
				dijkstraCHTest.localQueryForTransitNode(src[u], trg);
				for (NodeID j = 0; j < trg.size(); j++){
					//cout << trg[j] << " " << dijkstraCHTest.distanceTo(trg[j]) << endl;
					ma.set(u, j, dijkstraCHTest.distanceTo(trg[j]));					
				}
				for (NodeID k = 0; k < localContained.size(); k++){
					distTable.set(k, u, dijkstraCHTest.distanceTo(localContained[k]));
				}
				dijkstraCHTest.clear();
			}
		}

        // Note: Commented out code below because downstream function is causing error
        // and localManyToManyUsingTree is never used anywhere.
    
//      void localManyToManyUsingTree(Rectangle &rec, vector <CoordinateType> &xcord, vector <CoordinateType> &ycord, vector <NodeID> &src, 
//          vector <NodeID> &trg, Matrix <EdgeWeight> & ma)
//      {
//          processing::DijkstraCH<datastr::graph::UpdateableGraph, NormalPQueue, 2, false> dijkstraCHTest(upd);
//          for (NodeID u = 0; u < src.size(); u++){
//              //cout  << " source "<< src[u] << endl;
//              dijkstraCHTest.localQueryForTransitNode_USINGTREE(src[u], trg, rec.xleft, rec.yleft, rec.xright, rec.yright, xcord, ycord);
//              for (NodeID j = 0; j < trg.size(); j++){
//                  //cout << trg[j] << " " << dijkstraCHTest.distanceTo(trg[j]) << endl;
//                  ma.set(u, j, dijkstraCHTest.distanceTo(trg[j]));
//              }
//              dijkstraCHTest.clear();
//          }
//      }

		//given the bound node, select the cell's transit node
		void selectTransitNodeForCell(Cell &cell, vector <NodeID> &startBoundNode, vector <NodeID> &innerBoundNode,
			vector <NodeID> &outerBoundNode, vector <CoordinateType> &xcord, vector <CoordinateType> &ycord)
		{
			NodeID startSize = startBoundNode.size();
			NodeID outerSize = outerBoundNode.size();
			NodeID innerSize = innerBoundNode.size();
			if (outerSize == 0 || startSize == 0 || innerSize == 0){
				return;
			}

			Matrix <EdgeWeight> startToOuter(startSize, outerSize);
			startToOuter.init(Weight::MAX_VALUE);

			Matrix <EdgeWeight> innerToOuter(outerSize, innerSize);
			Matrix <EdgeWeight> innerToStart(startSize, innerSize);
			const bool performBucketScans = true;
			int earlyStopLevel = 10;
			ManyToMany<MyGraph, DijkstraManyToManyFW, DijkstraManyToManyBW, performBucketScans> mtmOuter(graph, earlyStopLevel);
			ManyToMany<MyGraph, DijkstraManyToManyFW, DijkstraManyToManyBW, performBucketScans> mtmStart(graph, earlyStopLevel);

			mtmOuter.computeMatrix(outerBoundNode, innerBoundNode, innerToOuter);
			mtmStart.computeMatrix(startBoundNode, innerBoundNode, innerToStart);
			int mappedValue = cell.xvertex * ascale + cell.yvertex;
			EdgeWeight d;
			EdgeWeight f_min;
			EdgeWeight temp_d;
			int temp_k;
			for (NodeID u = 0; u < startSize; u++){
				for (NodeID v = 0; v < outerSize; v++){			
					d = startToOuter.value(u, v);
					for (NodeID k = 0; k < innerSize; k++){							
						temp_d =  innerToOuter.value(v, k) + innerToStart.value(u, k);
						if (temp_d < d)
						{
							d = temp_d;
							temp_k = k;
							f_min = innerToStart.value(u, k);
						}
						//a shortest path may pass many transit nodes, for each shortest path
						//we chose the first vertex of the inner boundary nodes for which the 
						//shortest path passes
						if (temp_d == d){
							if (innerToStart.value(u, k) < f_min)
								temp_k = k;
							f_min = innerToStart.value(u, k);		
						}
					}
					allCellTransitNode.insert(innerBoundNode[temp_k]);
					// insert to the cell, use binary search
					int low = 2;
					int high =  eachCellTransitNode[mappedValue][1];
					int pos = binarysearch(eachCellTransitNode[mappedValue], low, high, innerBoundNode[temp_k]);
					if (pos < 0){
						if (high >= (int)eachCellTransitNode[mappedValue][0]){
							eachCellTransitNode[mappedValue] = (NodeID *)realloc(eachCellTransitNode[mappedValue], (high + EACH_CELL_TRANSIT_NUM) * sizeof(NodeID));
							//eachCellTransitNode[mappedValue][0] += avgNodesPerCell;								
							eachCellTransitNode[mappedValue][0] = high + EACH_CELL_TRANSIT_NUM;								
						}
						insert(eachCellTransitNode[mappedValue], -pos - 1 ,innerBoundNode[temp_k], high);
						eachCellTransitNode[mappedValue][1]++;
					}
				}
			}

			//Matrix <EdgeWeight> innerToOuter(outerSize, innerSize);
			//innerToOuter.init(Weight::MAX_VALUE);
			//Matrix <EdgeWeight> innerToStart(startSize, innerSize);
			//innerToStart.init(Weight::MAX_VALUE);
			//
			///*SpaceVertex lbouter, rtouter, lbinner, rtinner;
			//lbouter.set(cell.xvertex - 4, cell.yvertex - 4);
			//rtouter.set(cell.xvertex + 5, cell.yvertex + 5);
			//lbinner.set(cell.xvertex - 2, cell.yvertex - 2);
			//rtinner.set(cell.xvertex + 3, cell.yvertex + 3);
			//Rectangle outerRec(lbouter.getXPosition() * cellLength, lbouter.getYPosition() * cellLength, rtouter.getXPosition() * cellLength, rtouter.getYPosition() * cellLength);
			//Rectangle innerRec(lbinner.getXPosition() * cellLength, lbinner.getYPosition() * cellLength, rtinner.getXPosition() * cellLength, rtinner.getYPosition() * cellLength);
			//*/
			////Matrix <EdgeWeight> innerToOuterRef(innerSize, outerSize);
			//////localManyToManyUsingTree(outerRec, xcord, ycord, innerBoundNode, outerBoundNode, innerToOuterRef);
			////localManyToMany(innerBoundNode, outerBoundNode, innerToOuterRef);
			////for (NodeID i = 0; i < outerSize; i++){
			////	for (NodeID j = 0; j < innerSize; j++){
			////		//cout << innerToOuterRef.value(j, i) << endl;
			////		innerToOuter.set(i, j, innerToOuterRef.value(j, i));
			////	}
			////}
			//localManyToMany(startBoundNode, innerBoundNode, innerToStart);
			//localManyToMany(outerBoundNode, innerBoundNode, innerToOuter);
			////localManyToManyUsingTree(innerRec, xcord, ycord, startBoundNode, innerBoundNode, innerToStart);
			////cout << innerToStart << endl;
	

			////const bool performBucketScans = true;
			////int earlyStopLevel = 10;
			////ManyToMany<MyGraph, DijkstraManyToManyFW, DijkstraManyToManyBW, performBucketScans> mtmOuter(graph, earlyStopLevel);
			////ManyToMany<MyGraph, DijkstraManyToManyFW, DijkstraManyToManyBW, performBucketScans> mtmStart(graph, earlyStopLevel);
			////mtmOuter.computeMatrix(outerBoundNode, innerBoundNode, innerToOuter);
			////mtmStart.computeMatrix(startBoundNode, innerBoundNode, innerToStart);
			//		
			//NodeID mappedValue = cell.xvertex * ascale + cell.yvertex;
			//EdgeWeight d;
			////EdgeWeight d1;
			////EdgeWeight d2;
			//EdgeWeight temp_d;
			//NodeID temp_k;
			//EdgeWeight f_min;
			//for (NodeID u = 0; u < startSize; u++){
			//	for (NodeID v = 0; v < outerSize; v++){
			//		d = startToOuter.value(u, v);
			//		for (NodeID k = 0; k < innerSize; k++){							
			//			//if ((innerToOuter.value(v, k) == Weight::MAX_VALUE) || (innerToStart.value(u, k) == Weight::MAX_VALUE))
			//			//	continue;
			//			//cout << u << " " << v << " " << k << endl;
			//			temp_d =  innerToOuter.value(v, k) + innerToStart.value(u, k);
			//			//cout << innerToOuter.value(v, k)<< " " << innerToStart.value(u, k) << " "<< temp_d << endl;
			//			//cout << (temp_d < d) << endl;
			//			if (temp_d < d)
			//			{
			//				d = temp_d;
			//				temp_k = k;
			//				f_min = innerToStart.value(u, k);
			//			}
			//			if (temp_d == d){
			//				if (innerToStart.value(u, k) < f_min)
			//					temp_k = k;
			//				f_min = innerToStart.value(u, k);		
			//			}
			//		}
			//		if (d  < Weight::MAX_VALUE){
			//			allCellTransitNode.insert(innerBoundNode[temp_k]);
			//			//eachCellTransitNode[mappedValue].insert(innerBoundNode[temp_k]);
			//			// insert to the cell, use binary search
			//			int low = 2;
			//			int high =  eachCellTransitNode[mappedValue][1];
			//			int pos = binarysearch(eachCellTransitNode[mappedValue], low, high, innerBoundNode[temp_k]);
			//			//cout << pos << endl;
			//			if (pos < 0){
			//				if (high >= (int)eachCellTransitNode[mappedValue][0]){
			//					eachCellTransitNode[mappedValue] = (NodeID *)realloc(eachCellTransitNode[mappedValue], (high + EACH_CELL_TRANSIT_NUM) * sizeof(NodeID));
			//					eachCellTransitNode[mappedValue][0] += EACH_CELL_TRANSIT_NUM;								
			//				}
			//				insert(eachCellTransitNode[mappedValue], -pos - 1 ,innerBoundNode[temp_k], high);
			//				eachCellTransitNode[mappedValue][1]++;
			//			}	
			//		}
			//	}
			//}
		}
		
		//void selectTransitNodeForCell(Cell &cell, vector <NodeID> &startBoundNode, vector <NodeID> &innerBoundNode,
		//	vector <NodeID> &outerBoundNode, vector <CoordinateType> &xcord, vector <CoordinateType> &ycord)
		//{
		//	NodeID startSize = startBoundNode.size();
		//	NodeID outerSize = outerBoundNode.size();
		//	NodeID innerSize = innerBoundNode.size();
		//	if (outerSize == 0 || startSize == 0 || innerSize == 0){
		//		return;
		//	}
		//	//cout << startSize << " "  << outerSize << " "  << innerSize << endl; 
		//	Matrix <EdgeWeight> startToOuter(startSize, outerSize);
		//	startToOuter.init(Weight::MAX_VALUE);
		//	
		//	Matrix <EdgeWeight> innerToOuter(outerSize, innerSize);
		//	Matrix <EdgeWeight> innerToStart(startSize, innerSize);

		//	const bool performBucketScans = true;
		//	int earlyStopLevel = 10;
		//	ManyToMany<MyGraph, DijkstraManyToManyFW, DijkstraManyToManyBW, performBucketScans> mtmOuter(graph, earlyStopLevel);
		//	ManyToMany<MyGraph, DijkstraManyToManyFW, DijkstraManyToManyBW, performBucketScans> mtmStart(graph, earlyStopLevel);
		//	mtmOuter.computeMatrix(outerBoundNode, innerBoundNode, innerToOuter);
		//	mtmStart.computeMatrix(startBoundNode, innerBoundNode, innerToStart);
		//			
		//	NodeID mappedValue = cell.xvertex * ascale + cell.yvertex;
		//	EdgeWeight d;
		//	//EdgeWeight d1;
		//	//EdgeWeight d2;
		//	EdgeWeight temp_d;
		//	NodeID temp_k;
		//	for (NodeID u = 0; u < startSize; u++){
		//		for (NodeID v = 0; v < outerSize; v++){
		//			d = startToOuter.value(u, v);
		//			for (NodeID k = 0; k < innerSize; k++){							
		//				temp_d =  innerToOuter.value(v, k) + innerToStart.value(u, k);
		//				if (temp_d < d)
		//				{
		//					d = temp_d;
		//					temp_k = k;
		//				}
		//			}
		//			allCellTransitNode.insert(innerBoundNode[temp_k]);
		//			//eachCellTransitNode[mappedValue].insert(innerBoundNode[temp_k]);
		//			// insert to the cell, use binary search
		//			int low = 2;
		//			int high =  eachCellTransitNode[mappedValue][1];
		//			int pos = binarysearch(eachCellTransitNode[mappedValue], low, high, innerBoundNode[temp_k]);
		//			//cout << pos << endl;
		//			if (pos < 0){
		//				if (high >= eachCellTransitNode[mappedValue][0]){
		//					eachCellTransitNode[mappedValue] = (NodeID *)realloc(eachCellTransitNode[mappedValue], (high + EACH_CELL_TRANSIT_NUM) * sizeof(NodeID));
		//					eachCellTransitNode[mappedValue][0] += EACH_CELL_TRANSIT_NUM;								
		//				}
		//				insert(eachCellTransitNode[mappedValue], -pos - 1 ,innerBoundNode[temp_k], high);
		//				eachCellTransitNode[mappedValue][1]++;
		//			}						
		//		}
		//	}
		//}
		//
		public:
			/** Serializes the transit node to the given stream. */
			/**************************************************************************
			 * write ascale, range to the stream
			 * write each cell transit nodes
			 * write each cell contained nodes
			 * write transit distance table
			 * write each cell distance table
			*************************************************************************/
			void serialize(ostream& out) {	
				VERBOSE( cout << "serialize transit node begin... ascale: "  << ascale 
					<< " Range: " << range << " " << endl );
				// first write ascale, range
				PrimitiveSerializer<NodeID>::serialize(out, ascale);
				PrimitiveSerializer<CoordinateType>::serialize(out, range);

				//write all cell transit node
				SetSerializer<NodeID, NodeID>::serialize(out, allCellTransitNode);
				cout << "there are " << allCellTransitNode.size() << " transit nodes in total. " << endl;

				//write each cell transit nodes
				for (NodeID u = 0; u < allCellNum; u++){
					NodeID size = eachCellTransitNode[u][0];
					PrimitiveSerializer<NodeID>::serialize(out, size);
					for (NodeID i = 0; i < size; i++){
						PrimitiveSerializer<NodeID>::serialize(out, eachCellTransitNode[u][i]);
					}
				}

				//write each cell contained nodes
				for (NodeID u = 0; u < allCellNum; u++){
					NodeID size = eachCellContainedNode[u][0];
					PrimitiveSerializer<NodeID>::serialize(out, size);
					for (NodeID i = 0; i < size; i++){
						PrimitiveSerializer<NodeID>::serialize(out, eachCellContainedNode[u][i]);
					}
				}
				//write transit node matrix
				transitDistTable.serialize(out);
				
				//write each cell matrix
				for (NodeID u = 0; u < allCellNum; u++){
					eachCellDistTable[u].serialize(out) ;
				}
				VERBOSE( cout << "done." << endl );
				
			}

			/** Deserializes the graph from the given stream. */
			void deserialize(istream& in) {
// 				VERBOSE( cout << "transit node::deserialize " << flush );
				//fstream out("tnoutput/a.txt", ios::out);
				

				// first read ascale, range
				PrimitiveSerializer<NodeID>::deserialize(in, ascale);
				PrimitiveSerializer<CoordinateType>::deserialize(in, range);
// 				cout << ascale << " " << range << endl;
				
				cellLength = range / ascale;
				NodeID totalNodes = graph->noOfNodes();
				avgNodesPerCell = totalNodes / (ascale * ascale);
				if (avgNodesPerCell == 0){
					//avgNodesPerCell == 100;	//应该为=？earnestwu
					avgNodesPerCell = 100;
				}
				allCellNum = ascale * ascale;

				//read all cell transit node
				SetSerializer<NodeID, NodeID>::deserialize(in, allCellTransitNode);
// 				cout << "read all cell transit nodes " << allCellTransitNode.size() << " in all " << endl;
				

				//read each cell transit nodes
				NodeID size;		
				eachCellTransitNode = (NodeID **) malloc(sizeof(NodeID *) * allCellNum);
				for (NodeID u = 0; u < allCellNum; u++){
					PrimitiveSerializer<NodeID>::deserialize(in, size);
					eachCellTransitNode[u] = (NodeID *) malloc(sizeof(NodeID) * size);
					for (NodeID i = 0; i < size; i++){
						PrimitiveSerializer<NodeID>::deserialize(in, eachCellTransitNode[u][i]);
					}
				}		
// 				cout << "each cell transit node read finished " << endl;

				//read each cell contained nodes
				eachCellContainedNode = (NodeID **) malloc (sizeof(NodeID *) * allCellNum);
				for (NodeID u = 0; u < allCellNum; u++){
					PrimitiveSerializer<NodeID>::deserialize(in, size);
					eachCellContainedNode[u] = (NodeID *) malloc(sizeof(NodeID *) * size);
					for (NodeID i = 0; i < size; i++){
						PrimitiveSerializer<NodeID>::deserialize(in, eachCellContainedNode[u][i]);
					}
				}
// 				cout << "each cell contained node read finished " << endl;

				
				// build the  transit node mapping
				transitNodeMapping.clear();
				set <NodeID>::iterator sit = allCellTransitNode.begin();
				int count = 0;
				while (sit != allCellTransitNode.end()){
					transitNodeMapping.insert(make_pair(*sit, count));
					count++;
					sit++;
				}
				
				//read transit node matrix
				transitDistTable.deserialize(in);
				
// 				cout << "read all pair transit node distance table, finished " << endl;
				

				eachCellDistTable = new Matrix <EdgeWeight> [allCellNum];
				//read each cell matrix
				for (NodeID u = 0; u < allCellNum; u++){
					eachCellDistTable[u].deserialize(in) ;
				}

				eachCellTransitMapping = (NodeID **) malloc(sizeof(NodeID *) * allCellNum);
				for (NodeID u = 0; u < allCellNum; u++){
					NodeID size = eachCellTransitNode[u][0];
					eachCellTransitMapping[u] = (NodeID *) malloc(sizeof(NodeID) * size);
					eachCellTransitMapping[u][0] = eachCellTransitNode[u][0];
					eachCellTransitMapping[u][1] = eachCellTransitNode[u][1];
					for (NodeID i = 2; i < size; i++){
						eachCellTransitMapping[u][i] = transitNodeMapping[eachCellTransitNode[u][i]];
					}
				}
				
// 				cout << "reading each matrix over " << endl;
// 				VERBOSE( cout << " ascale: "  << ascale << " Range: " << range << " " << endl );
// 				VERBOSE( cout << "done..." << endl );
			} 
			
			void deserialize_withMatrix(istream& in) {
// 				VERBOSE( cout << "transit node::deserialize " << flush );
				//fstream out("tnoutput/a.txt", ios::out);


				// first read ascale, range
				PrimitiveSerializer<NodeID>::deserialize(in, ascale);
				PrimitiveSerializer<CoordinateType>::deserialize(in, range);
				//out << ascale << " " << range << endl;

				cellLength = range / ascale;
				NodeID totalNodes = graph->noOfNodes();
				avgNodesPerCell = totalNodes / (ascale * ascale);
				allCellNum = ascale * ascale;

				//read all cell transit node
				SetSerializer<NodeID, NodeID>::deserialize(in, allCellTransitNode);
				//read each cell transit nodes
				NodeID size;

				eachCellTransitNode = (NodeID **) malloc(sizeof(NodeID *) * allCellNum);
				for (NodeID u = 0; u < allCellNum; u++){
					PrimitiveSerializer<NodeID>::deserialize(in, size);
					eachCellTransitNode[u] = (NodeID *) malloc(sizeof(NodeID) * size);
					for (NodeID i = 0; i < size; i++){
						PrimitiveSerializer<NodeID>::deserialize(in, eachCellTransitNode[u][i]);
					}
				}


				//read each cell contained nodes
				eachCellContainedNode = (NodeID **) malloc (sizeof(NodeID *) * allCellNum);
				for (NodeID u = 0; u < allCellNum; u++){
					PrimitiveSerializer<NodeID>::deserialize(in, size);
					eachCellContainedNode[u] = (NodeID *) malloc(sizeof(NodeID *) * size);
					for (NodeID i = 0; i < size; i++){
						PrimitiveSerializer<NodeID>::deserialize(in, eachCellContainedNode[u][i]);
					}
				}


				// build the  transit node mapping
				//fstream out("tnoutput/alltn.txt", ios::out);
				transitNodeMapping.clear();
				set <NodeID>::iterator sit = allCellTransitNode.begin();
				int count = 0;
				while (sit != allCellTransitNode.end()){
					transitNodeMapping.insert(make_pair(*sit, count));
					count++;
					sit++;
				}

				eachCellDistTable = new Matrix <EdgeWeight>[allCellNum];
				//read each cell matrix
				for (NodeID u = 0; u < allCellNum; u++){
					eachCellDistTable[u].deserialize(in) ;
				}

// 				cout << "Reading each matrix over " << endl;

				const bool performBucketScans = true;
				int earlyStopLevel = 10;
				NodeID allTransitNum = allCellTransitNode.size();
				//read transit node matrix
				//transitDistTable.deserialize(in);
				transitDistTable.setRowAndCol(allTransitNum, allTransitNum);
				transitDistTable.init(0);
				ManyToMany<MyGraph, DijkstraManyToManyFW, DijkstraManyToManyBW, performBucketScans> mtm(graph, earlyStopLevel);
				vector<NodeID> allCellTransitCopy(allCellTransitNode.begin(), allCellTransitNode.end());
				mtm.computeMatrix(allCellTransitCopy, allCellTransitCopy, transitDistTable);
// 				VERBOSE( cout << " ascale: "  << ascale << " Range: " << range << " " << endl );
// 				VERBOSE( cout << "done..." << endl );
			} 


	
private:
		Graph* graph = NULL;
		datastr::graph::UpdateableGraph *upd = NULL;

		NodeID ascale;
		CoordinateType range;

        //one case: use ch for local query
        processing::DijkstraCH<datastr::graph::SearchGraph, NormalPQueue, 2, true> dijkstraCHTest;

        //the other case: use bidijkstra for local query
        processing::DijkstraCH<datastr::graph::UpdateableGraph, NormalPQueue, 2, true> dijkstraBIDTest;

		//can be computed from graph, ascale, range
		NodeID avgNodesPerCell;
		NodeID allCellNum;
		CoordinateType cellLength;

		//each border edge is associated with an edgeBucket
		vector < vector <EasyEdge> > edgeBucket;
		// border edge hash map stores the mapping from  every border edge to the corresponding edge bucket
		BorderEdgeHashMap beHashMap;
		//all cell's transit node
		set <NodeID> allCellTransitNode;
		map <NodeID, NodeID> transitNodeMapping;
		
		//vector <EachCellEntry> cellEntry;

		// each cell's transit node
		NodeID **eachCellTransitNode = NULL;
		//vector < set <NodeID> > eachCellTransitNode;		
		//each cell's original node
		NodeID **eachCellContainedNode = NULL;
		
		Matrix <EdgeWeight> transitDistTable;
		Matrix<EdgeWeight> *eachCellDistTable = NULL;

		NodeID **eachCellTransitMapping = NULL;

};
}

typedef processing::TransitNode <datastr::graph::SearchGraph> TransitNodeTest;

#endif