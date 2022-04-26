/* Copyright (C) 2012 Lingkun Wu, Xiaokui Xiao, Dingxiong Deng, Gao Cong, Andy Diwen Zhu, Shuigeng Zhou
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

#include "CoordinateLoader.h"

#include <limits>
#include <iostream>
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

namespace CoordinateLoader
{
    void serializeCoordinate(std::string outputFilePath, std::vector<CoordinateType>& _xcord, std::vector<CoordinateType>& _ycord) {
        std::ofstream outputFS(outputFilePath, std::ios::out | std::ios::binary);
        if (outputFS.is_open()) {
            boost::archive::binary_oarchive oa(outputFS);
            oa << _xcord << _ycord;
        } else {
            std::cerr << "Cannot open binary output file " << outputFilePath << std::endl;
            exit(1);
        }
    }
    
    void deserializeCoordinate(std::string inputFilePath, std::vector<CoordinateType>& _xcord, std::vector<CoordinateType>& _ycord) {
        std::ifstream inputFS(inputFilePath, std::ios::in | std::ios::binary);
        if (inputFS.is_open()) {
            boost::archive::binary_iarchive ia(inputFS);
            ia >> _xcord >> _ycord;
        } else {
            std::cerr << "Cannot open binary input file " << inputFilePath << std::endl;
            exit(1);
        }
    }
    
    CoordinateType getCoordinateRange(std::vector<CoordinateType>& _xcord, std::vector<CoordinateType>& _ycord) {
        CoordinateType g_dMinX = std::numeric_limits<int>::max();
        CoordinateType g_dMaxX = std::numeric_limits<int>::min();
        CoordinateType g_dMinY = std::numeric_limits<int>::max();
        CoordinateType g_dMaxY = std::numeric_limits<int>::min();

        CoordinateType x, y;
        for (NodeID u = 0; u < _xcord.size(); u++){
            //std::cout << v << x << y << std::endl;
            x = _xcord[u];
            y = _ycord[u];

            if (x < g_dMinX)
                g_dMinX = x;
            if (x > g_dMaxX)
                g_dMaxX = x;
            if (y < g_dMinY)
                g_dMinY = y;
            if (y > g_dMaxY)
                g_dMaxY = y;
        }
        CoordinateType g_xRange = g_dMaxX - g_dMinX;
        CoordinateType g_yRange = g_dMaxY - g_dMinY;
        CoordinateType g_Range = (g_xRange < g_yRange) ? g_yRange : g_xRange;
        CoordinateType g_minRange = (g_xRange < g_yRange) ? g_xRange : g_yRange;
        
//         std::cout << "Coordinate Range: " << std::endl;
//         std::cout << "x range: " << g_dMaxX - g_dMinX << " y range: " << g_dMaxY - g_dMinY << std::endl;
//         std::cout << "range : " << g_Range << " minRange: " << g_minRange << std::endl;
//         std::cout << "Min (X,Y) = " << "(" << static_cast<int>(g_dMinX) << "," << static_cast<int>(g_dMinY) << ")" << std::endl;
//         std::cout << "Max (X,Y) = " << "(" << static_cast<int>(g_dMaxX) << "," << static_cast<int>(g_dMaxY) << ")" << std::endl;
        return g_Range;
    }
    
    // Also returns the maximum coordinate range
    void readCoordinate(std::string coFilePath, std::vector<CoordinateType>& _xcord, std::vector<CoordinateType>& _ycord) {
        std::ifstream coFS(coFilePath, std::ios::in);
        if (coFS.is_open()) {
            NodeID numNodes;
            coFS >> numNodes;
            CoordinateType xmin, xmax, ymin, ymax;
            xmin = std::numeric_limits<int>::max();
            xmax = std::numeric_limits<int>::min();;
            ymin = std::numeric_limits<int>::max();
            ymax = std::numeric_limits<int>::min();;

            _xcord.resize(numNodes);
            _ycord.resize(numNodes);
            NodeID v;
            CoordinateType x, y;
            for (NodeID u = 0; u < numNodes; u++){
                coFS >> v >> x >> y;
                //std::cout << v << x << y << endl;
                _xcord[v] = x;
                _ycord[v] = y;

                if (x < xmin)
                    xmin = x;
                if (x > xmax)
                    xmax = x;
                if (y < ymin)
                    ymin = y;
                if (y > ymax)
                    ymax = y;
            }
            CoordinateType g_xRange = xmax - xmin;
            CoordinateType g_yRange = ymax - ymin;
            CoordinateType g_Range = (g_xRange < g_yRange) ? g_yRange : g_xRange;
            CoordinateType g_minRange = (g_xRange < g_yRange) ? g_xRange : g_yRange;
            
//             std::cout << "Coordinate Range: " << std::endl;
//             std::cout << "x range: " << xmax - xmin << " y range: " << ymax - ymin << std::endl;
//             std::cout << "range : " << g_Range << " minRange: " << g_minRange << std::endl;
//             std::cout << "Min (X,Y) = " << "(" << static_cast<int>(xmin) << "," << static_cast<int>(ymin) << ")" << std::endl;
//             std::cout << "Max (X,Y) = " << "(" << static_cast<int>(xmax) << "," << static_cast<int>(ymax) << ")" << std::endl;
        } else {
            std::cerr << "Cannot open text coordinate file " << coFilePath << std::endl;
            exit(1);
        }
        coFS.close();
    }
    
};
