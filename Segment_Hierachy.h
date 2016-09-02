//
// Created by icme-huang on 8/28/16.
//

#ifndef INC_2DFLUIDSTRUCTURE_SEGMENT_HIERACHY_H
#define INC_2DFLUIDSTRUCTURE_SEGMENT_HIERACHY_H


#include "Vector2D.h"
#include "Array.h"
#include "KD_Tree.h"

class Segment_Hierachy {
public:
    Array<Vector2D<double>>  points;
    Array<Vector2D<int>>     edgeSet;

 //Segment_Hierachy(TRIANGLE_MESH& triangle_mesh_input,GEOMETRY_PARTICLES<VECTOR<T,3> >& particles_input,ARRAY<TRIANGLE_3D<T> >& triangle_list_input,const bool update_boxes=true,const int triangles_per_group=0);

    virtual ~Segment_Hierachy();
private:

    void Initialization_Hierarchy_Using_KDTree();

};


#endif //INC_2DFLUIDSTRUCTURE_TRIANGLE_HIERACHY_H
