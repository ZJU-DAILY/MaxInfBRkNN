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

#include "StopWatch.h"

#include <sstream>

StopWatch::StopWatch()
{
    this->durationNs = this->durationNs.zero();
    this->pausedDurationNs = this->pausedDurationNs.zero();
}

double StopWatch::getTimeNs()
{
    // Assumes resume() has been called before stop() already
    return (this->durationNs-this->pausedDurationNs).count();
}

double StopWatch::getTimeUs()
{
    return this->getTimeNs()/1000;
}

double StopWatch::getTimeMs()
{
    return this->getTimeNs()/1000000;
}

void StopWatch::start()
{
    this->t1 = std::chrono::high_resolution_clock::now();
}

void StopWatch::stop()
{
    this->t2 = std::chrono::high_resolution_clock::now();
    this->durationNs = std::chrono::duration_cast<std::chrono::nanoseconds>(this->t2 - this->t1);
}

void StopWatch::pause()
{
    this->p_t1 = std::chrono::high_resolution_clock::now();
}

void StopWatch::resume()
{
    this->p_t2 = std::chrono::high_resolution_clock::now();
    // Accumulate amount of time that been paused 
    this->pausedDurationNs += std::chrono::duration_cast<std::chrono::nanoseconds>(this->p_t2 - this->p_t1);

}

void StopWatch::reset()
{
    this->durationNs = this->durationNs.zero();
    this->pausedDurationNs = this->pausedDurationNs.zero();
}
