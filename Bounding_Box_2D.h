//
// Created by icme-huang on 8/29/16.
//

#ifndef INC_2DFLUIDSTRUCTURE_BOUNDING_BOX_2D_H
#define INC_2DFLUIDSTRUCTURE_BOUNDING_BOX_2D_H

#include "Vector2D.h"
#include "Array.h"

class Bounding_Box_2D {
public:
    /*minCorner: vector (x_min,y_min)
     *maxCorner: vector (x_max,y_max)
     */
    Vector2D<double> minCorner, maxCorner;
public:
    Bounding_Box_2D():minCorner(0.0,0.0),maxCorner(0.0,0.0){}
    /* Initialize Bounding box with an array of points*/
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
    Bounding_Box_2D(const Vector2D<double> &minCornerInput, const Vector2D<double> &maxCornerInput)
            :minCorner(minCornerInput),maxCorner(maxCornerInput){}

    friend std::ostream& operator<<(std::ostream& cout,const Bounding_Box_2D& box){
        cout<<"[" << box.minCorner[0] <<"," << box.maxCorner[0] <<"] " <<"[" << box.minCorner[1] <<"," << box.maxCorner[1] <<"]"<<std::endl;
        return cout;
    };

    /* return the direction of large edge of the bounding box , 0 means x axis , 1 means y asix */
    int Max_Length_Axis() const{return ((maxCorner[0] - minCorner[0]) < (maxCorner[1] - minCorner[1]) ? 1: 0); }

    /* enlarge bounding box*/
    void Enlarge(double thickness) {maxCorner += thickness;minCorner -=thickness;}

    /* This is a static function, return a bounding box containing two different bounding boxes */
    static Bounding_Box_2D Combine(const Bounding_Box_2D &box1, const Bounding_Box_2D &box2){
        return Bounding_Box_2D(Vector2D<double>::Componentwise_Min(box1.minCorner,box2.minCorner),Vector2D<double>::Componentwise_Max(box1.maxCorner,box2.maxCorner));
    }
    ~Bounding_Box_2D(){}

};



#endif //INC_2DFLUIDSTRUCTURE_BOUNDING_BOX_2D_H
