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

#ifndef _STOPWATCH_H
#define _STOPWATCH_H

#include <string>
#include <chrono>

class StopWatch {
    
public:
    StopWatch();
    void start();
    void stop();
    void pause();
    void resume();
    void reset();
    double getTimeNs();
    double getTimeUs();
    double getTimeMs();
    
private:
    bool paused;
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::high_resolution_clock::time_point p_t1;
    std::chrono::high_resolution_clock::time_point p_t2;
    std::chrono::nanoseconds durationNs;
    std::chrono::nanoseconds pausedDurationNs;
};

#endif // _STOPWATCH_H
