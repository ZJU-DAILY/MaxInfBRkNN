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

#ifndef ELIMINATIONWEIGHT_H
#define ELIMINATIONWEIGHT_H

#include <limits>

/**
 * Encapsulates several special values of the elimination weights or 
 * priority terms for the node ordering. This class is required by the
 * BinaryHeap priority queue.
 * Here: elimination weight type = double
 */
class EliminationWeight
{
public:
    typedef double Type;
    
    /**
     * max value of an edge weight (due to the chosen type)
     * (used to represent 'infinity')
     */
    /*static const Type MAX_VALUE = __DBL_MAX__;
    static const Type MIN_VALUE = -__DBL_MAX__;*/
	static const Type MAX_VALUE;
	static const Type MIN_VALUE;
    
};

const double EliminationWeight::MAX_VALUE = __DBL_MAX__;
const double EliminationWeight::MIN_VALUE = -__DBL_MAX__;


#endif // ELIMINATIONWEIGHT_H
