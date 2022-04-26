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

#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <ctime>
//#include <sys/time.h>


/** Returns a timestamp ('now') in seconds (incl. a fractional part). */
inline clock_t timestamp() {
   /* struct timeval tp;
    gettimeofday(&tp, NULL);
    return double(tp.tv_sec) + tp.tv_usec / 1000000.;*/
	return clock()/CLOCKS_PER_SEC*1000; // Return current time in ms
}


/**
 * Provides methods that can be used to display the
 * progress of a procedure.
 */
class Percent
{
 private:
    typedef int percentInt;
    typedef unsigned int valueInt;
    
 public:
    /**
     * Constructor.
     * @param maxValue the value that corresponds to 100%
     *                 (or 0% if reverse counting is activated)
     * @param reverse instead of raising from 0 to maxValue,
     *                the value falls from maxValue to 0
     * @param step the progress is shown in steps of 'step' percent
     */
    Percent(valueInt maxValue, bool reverse = false, percentInt step = 2) {
        reinit(maxValue, reverse, step);
    }

    /** Reinitializes this object. */
    void reinit(valueInt maxValue, bool reverse = false, percentInt step = 2) {
        _maxValue = maxValue;
        _intervalPercent = _maxValue / 100;
        _nextThreshold = _intervalPercent;
        _reverse = reverse;
        _lastPercent = 0;
        _step = step;
        
        if (reverse) _nextThreshold = _maxValue - _intervalPercent;
    }

    /** If there has been significant progress, display it. */
    void printStatus(valueInt currentValue) {
        if ( ! _reverse) {
            if (currentValue >= _nextThreshold) {
                _nextThreshold += _intervalPercent;
                printPercent( currentValue / (double)_maxValue * 100 );
            }
            if (currentValue + 1 == _maxValue) finish();            
        }
        else {
            if (currentValue <= _nextThreshold) {
                _nextThreshold -= _intervalPercent;
                printPercent( 100 - (currentValue / (double)_maxValue * 100) );
            }
        }
    }
        
 private:
    valueInt _maxValue;
    valueInt _intervalPercent;
    valueInt _nextThreshold;
    bool _reverse;
    percentInt _lastPercent;
    percentInt _step;    

    /** Displays the new progress. */
    void printPercent(double percent) {
        while (percent >= _lastPercent+_step) {
            _lastPercent+=_step;
            if (_lastPercent % 10 == 0) {
                std::cout << " " << _lastPercent << "% ";
            }
            else {
                std::cout << ".";
            }
            std::cout.flush();
        }
    }

    void finish() {
        std::cout << " 100%" << std::endl;
    }
};


#endif // UTILS_H
