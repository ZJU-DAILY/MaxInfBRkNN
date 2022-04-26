/* Original Work Copyright (C) 2005, 2006, 2007, 2008 Robert Geisberger, Dominik Schultes, Peter Sanders, Universitaet Karlsruhe (TH)
 * Modified Work Copyright (C) 2012 Lingkun Wu, Xiaokui Xiao, Dingxiong Deng, Gao Cong, Andy Diwen Zhu, Shuigeng Zhou
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

#ifndef EDGEWEIGHT_H
#define EDGEWEIGHT_H

#include <limits>

typedef unsigned int EdgeWeight;
typedef int SignedEdgeWeight;


/**
 * Encapsulates several special values of edge weights which are used as flags.
 * The concrete implementation depends on the chosen type of the edge weights.
 * Here: edge weight type = unsigned integer
 */
class Weight
{
public:
    static const SignedEdgeWeight SIGNED_MIN_VALUE = -__INT_MAX__ - 1;
    static const EdgeWeight SIGNED_MAX_VALUE = (EdgeWeight)__INT_MAX__;
    
    /**
     * max value of an edge weight (due to the chosen type)
     * (used to represent 'infinity')
     */
    static const EdgeWeight MAX_VALUE = __INT_MAX__ * 2U + 1; 
    static const EdgeWeight MIN_VALUE = 0;
    
    /**
     * max integer value of an edge weight which is actually used;
     * the highest values can be used as special values (flags);
     * all edge weights between 0 and MAX_INTEGER have no special meaning;
     * MAX_INTEGER has to be less than MAX_VALUE
     */
    static const EdgeWeight MAX_INTEGER = MAX_VALUE - 1;
};

const EdgeWeight Weight::MAX_VALUE;

typedef unsigned long long Checksum;

//#include "datastr/graph/weightCompress.h"
//typedef WeightCompressorIdentity WeightCompressor;

#endif // EDGEWEIGHT_H
