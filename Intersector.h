//
// Created by icme-huang on 8/28/16.
//

#include "Vector2D.h"


#ifndef INC_2DFLUIDSTRUCTURE_INTERSECTOR_H
#define INC_2DFLUIDSTRUCTURE_INTERSECTOR_H

class Segment_Hierachy;

class Intersector {
public:
    /*Fluid information*/
    Vector2D<double> *X;

    /*Structure information*/
    Vector2D<double> *Xs;         //current structure boundary position
    Vector2D<double> *Xs_n;       //structure boundary position of last step
    Segment_Hierachy * segmentHierachy;


};


#endif //INC_2DFLUIDSTRUCTURE_INTERSECTOR_H
