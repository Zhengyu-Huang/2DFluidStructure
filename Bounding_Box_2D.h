//
// Created by icme-huang on 8/29/16.
//

#ifndef INC_2DFLUIDSTRUCTURE_BOUNDING_BOX_2D_H
#define INC_2DFLUIDSTRUCTURE_BOUNDING_BOX_2D_H

#include "Vector2D.h"
#include "Array.h"
class Bounding_Box_2D {
public:
    Vector2D<double> minCorner, maxCorner;
public:
    Bounding_Box_2D(const Array<Vector2D<double>> &points){
        minCorner = points[0];
        maxCorner = points[0];
        for (int i = 1; i < points.len; i++) {
            for(int d = 0; d < 2 ; d++)
            if(points[i][d] < minCorner[d])
                minCorner[d] = points[i][d];
            else if(points[i][d] > maxCorner[d])
                maxCorner[d] = points[i][d];
        }

    }
    int Max_Length_Axis() const{return ((maxCorner[0] - minCorner[0]) < (maxCorner[1] - minCorner[1]) ? 1: 0); }

    ~Bounding_Box_2D(){}

};


#endif //INC_2DFLUIDSTRUCTURE_BOUNDING_BOX_2D_H
