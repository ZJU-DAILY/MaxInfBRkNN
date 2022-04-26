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

#ifndef _GEOMETRY_H
#define _GEOMETRY_H

#include "../common.h"

#include <vector>

namespace geometry
{
    bool intervalsIntersect(DistanceBound lowerBound1, DistanceBound upperBound1, DistanceBound lowerBound2, DistanceBound upperBound2);
    void unionIntervals(DistanceBound& lowerBound, DistanceBound& upperBound, std::vector<DistanceBound>& lowerBounds, std::vector<DistanceBound>& upperBounds);
    void getRightTopCorner(Coordinate left, Coordinate bottom, int width, Coordinate& right, Coordinate& top);
    bool isOverlapping(Coordinate left1, Coordinate bottom1, Coordinate right1, Coordinate top1,
                       Coordinate left2, Coordinate bottom2, Coordinate right2, Coordinate top2);
    bool getIntersectingRectangle(Coordinate left1, Coordinate bottom1, Coordinate right1, Coordinate top1,
                                  Coordinate left2, Coordinate bottom2, Coordinate right2, Coordinate top2, 
                                  Coordinate& left, Coordinate& bottom, Coordinate& right, Coordinate& top);
    EuclideanDistance getMinDistToRectangle(Coordinate x, Coordinate y, Coordinate left, Coordinate bottom, Coordinate right, Coordinate top);
    EuclideanDistance getMaxDistToRectangle(Coordinate x, Coordinate y, Coordinate left, Coordinate bottom, Coordinate right, Coordinate top);
    void translateCoordinates(Coordinate& x, Coordinate& y, int xTranslation, int yTranslation);
    EuclideanDistance getEuclideanDist(Coordinate x1, Coordinate y1, Coordinate x2, Coordinate y2);
    double getEuclideanDist_jins(Coordinate x1, Coordinate y1, Coordinate x2, Coordinate y2);
    double getEuclideanDist_BRI(float x1, float y1, float x2, float y2);

    unsigned int maxPower2(unsigned int x);
    bool contains(Coordinate left, Coordinate bottom, Coordinate right, Coordinate top, Coordinate x, Coordinate y);
}

#endif // _GEOMETRY_H