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

#include "geometry.h"

#include <algorithm>
#include <cmath>
#include <assert.h>

namespace geometry
{
    
    bool intervalsIntersect(DistanceBound lowerBound1, DistanceBound upperBound1, DistanceBound lowerBound2, DistanceBound upperBound2) {
        return lowerBound1 <= upperBound2 && lowerBound2 <= upperBound1; 
    }

    void unionIntervals(DistanceBound &lowerBound, DistanceBound &upperBound, std::vector<DistanceBound>& lowerBounds, std::vector<DistanceBound>& upperBounds) {
        // The paper defines union as the tightest interval that contains all other intervals
        assert (lowerBounds.size() > 0 && "No intervals provided");
        assert (lowerBounds.size() == upperBounds.size() && "Unequal number of upper and lower bounds to union");
        if (lowerBounds.size() > 0) {
            lowerBound = lowerBounds[0];
            upperBound = upperBounds[0];
            for (std::size_t i = 1; i < lowerBounds.size(); ++i) {
                if (lowerBound > lowerBounds[i]) {
                    lowerBound = lowerBounds[i];
                }
                if (upperBound < upperBounds[i]) {
                    upperBound = upperBounds[i];
                }
            }
        }
    }
    
    void getRightTopCorner(Coordinate left, Coordinate bottom, int width, Coordinate& right, Coordinate& top)
    {
        // Width and height are the same since we are working with square region
        // And we use fact that left (resp. bottom) is smaller than right (resp. top)
        right  = left + width;
        top = bottom + width;
    }
    
    bool isOverlapping(Coordinate left1, Coordinate bottom1, Coordinate right1, Coordinate top1, 
                       Coordinate left2, Coordinate bottom2, Coordinate right2, Coordinate top2)
    {
        if (right1 <= left2 || right2 <= left1 || top2 <= bottom1 || top1 <= bottom2) {
            return false;
        } else {
            return true;
        }
    }

    bool getIntersectingRectangle(Coordinate left1, Coordinate bottom1, Coordinate right1, Coordinate top1, 
                                  Coordinate left2, Coordinate bottom2, Coordinate right2, Coordinate top2, 
                                  Coordinate& left, Coordinate& bottom, Coordinate& right, Coordinate& top)
    {
        // Return false and does not set valid values if input rectangles do not overlap
        left = left1 > left2 ? left1 : left2; // max(left1,left)
        right = right1 < right2 ? right1 : right2; // min(right1,right2)
        if (left > right) {
            return false;
        }
        bottom = bottom1 > bottom2 ? bottom1 : bottom2; // max(bottom1,bottom2)
        top = top1 < top2 ? top1 : top2; // min(top1,top2)
        if (bottom > top) {
            return false;
        }
        // Note: If they overlap at a single point true will be returned
        return true;
    }

    EuclideanDistance getMinDistToRectangle(Coordinate x, Coordinate y, Coordinate left, Coordinate bottom, Coordinate right, Coordinate top)
    {
        EuclideanDistance dist;

        // First get nearest point on the rectangle (not necessarily a corner)
        Coordinate nearestX, nearestY;
        bool xWithinPerimeter = false;
        bool yWithinPerimeter = false;
        // We use the fact that left is always smaller than right and 
        // similarly for bottom and top
        if (x >= left && x <= right) {
            // This means the x value of the nearest point is the same as
            // query point, so the min dist line will be perpendicular
            xWithinPerimeter = true;
            nearestX = x;
        } else if (x > right) {
            // We already know left < right by defintion
            nearestX = right;
        } else /*if (x < left)*/ {
            // Since left < right we cannot have a situation where
            // x > right and x < left at the same time. Similarly for 
            nearestX = left;
        }
        
        if (y >= bottom && y <= top) {
            // This means the y value of the nearest point is the same as
            // query point, so the min dist line will be perpendicular
            yWithinPerimeter = true;
            nearestY = y;
        } else if (y > top) {
            // We already know bottom < top by defintion
            nearestY = top;
        } else /*if (y < bottom)*/ {
            nearestY = bottom;
        }
        
        //std::cout << "Nearest Point = (" << nearestX << "," << nearestY << ")" << std::endl;
        
        // If the min dist line is perpendicular we can avoid an
        // expensive Euclidean distance calculation
        if (xWithinPerimeter) {
            dist = std::abs(nearestY-y);
        } else if (yWithinPerimeter) {
            dist = std::abs(nearestX-x);
        }else {
            dist = std::sqrt(std::pow(nearestX-x,2) + std::pow(nearestY-y,2));
        }

        return dist;
    }

    EuclideanDistance getMaxDistToRectangle(Coordinate x, Coordinate y, Coordinate left, Coordinate bottom, Coordinate right, Coordinate top)
    {
        EuclideanDistance dist;
        
        // First find the furthest point on the rectangle
        // Unlike minimum distance case max dist will always be on a corner

        Coordinate furthestX, furthestY;

        // Again we use fact that left (resp. bottom) is smaller than right (resp. top)
        EuclideanDistance midX = left + (right-left)/2, midY = bottom + (top-bottom)/2;

        if (x < midX) {
            furthestX = right;
        } else {
            furthestX = left;
        }
        
        if (y < midY) {
            furthestY = top;
        } else {
            furthestY = bottom;
        }
        
        //std::cout << "Furthest Point = (" << furthestX << "," << furthestY << ")" << std::endl;

        // And the max dist will never be a perpendicular distance
        dist = std::sqrt(std::pow(furthestX-x,2) + std::pow(furthestY-y,2));
        
        return dist;
        
    }
    
    void translateCoordinates(Coordinate& x, Coordinate& y, int xTranslation, int yTranslation)
    {
        x = x - xTranslation;
        y = y - yTranslation;
    }
    
    EuclideanDistance getEuclideanDist(Coordinate x1, Coordinate y1, Coordinate x2, Coordinate y2)
    {
        return std::sqrt(std::pow(x2-x1,2) + std::pow(y2-y1,2));;
    }

    double getEuclideanDist_jins(Coordinate x1, Coordinate y1, Coordinate x2, Coordinate y2){
        double edge_dist = std::sqrt(std::pow(x2-x1,2) + std::pow(y2-y1,2));
        return edge_dist;
    }

    double getEuclideanDist_BRI(float x1, float y1, float x2, float y2){
        float edge_dist = std::sqrt(std::pow(x2-x1,2) + std::pow(y2-y1,2));
        return edge_dist;
    }
    
    unsigned int maxPower2(unsigned int x) {
        assert(x > 0);
        x = x | (x >> 1);
        x = x | (x >> 2);
        x = x | (x >> 4);
        x = x | (x >> 8);
        x = x | (x >>16);
        return x - (x >> 1);
    }
    
    bool contains(Coordinate left, Coordinate bottom, Coordinate right, Coordinate top, Coordinate x, Coordinate y)
    {
        if (x >= left && x <= right && y >= bottom && y <= top) {
            return true;
        } else {
            return false;
        }
    }

}