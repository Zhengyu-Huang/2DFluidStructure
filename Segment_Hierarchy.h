//
// Created by icme-huang on 8/28/16.
//

#ifndef INC_2DFLUIDSTRUCTURE_SEGMENT_HIERACHY_H
#define INC_2DFLUIDSTRUCTURE_SEGMENT_HIERACHY_H


#include "Vector2D.h"
#include "Array.h"
#include "KD_Tree.h"
#include <stack>

class Segment_Hierarchy {
public:
    Array<Vector2D<double>>  points; // coordinate of all points
    Array<Vector2D<int>>     edgeSet;// edge set, in terms of point id

    int numLeaves,root;              // numLeaves: number of leaves of the bounding box hierarchy tree
                                     // root: id of root box, root+1 is the number of boxes in box hierarchy
                                     // root = 2*numLeaves - 2
    Array<Bounding_Box_2D> boxHierarchy; // boxHierarchy has tree structure, which have numLeaves leaf boxes for each edge,
                                         // their id are 0,1,2...numLeaves-1, and parent bounding boxes with id numLeaves,
                                         // numLeaves+1 ... root, the parent bounding boxe is simply formulated by combination
                                         // of its two children box

    Array<int> parents;              // parents[i] is the id of node i's parent
    Array<Vector2D<int>> children;   // children[i][0] children[i][1] are the children id of node (i + numLeaves)

    /* constructor of Segment Hierarchy with points and edgeSet */
    Segment_Hierarchy(const Array<Vector2D<double>> & points,const Array<Vector2D<int>> & edgeSet,const bool update_boxes=true);

    /* Print information of the segment hierarchy */
    void Print_Info();

    mutable std::stack<int> traversalStack; /* used for traverse the segment hierarchy tree */

    /* output is candidates
     * output all these leaf boxes, which might intersect box.
     * For the box-box intersection, which is lazy intersection dectection
     * the testBox is enlarged by thickness
     */
    void Intersection_List(Bounding_Box_2D &testBox, Array<int> &candidates,double thickness);
    void Intersection_List_Helper(int root, Bounding_Box_2D box, Array<int> &candidates,double thicknessParameter);

    /* if box is leaf box return true, else return false */
    bool Leaf(const int box) const {return box < numLeaves;}

    ~Segment_Hierarchy(){};
private:
    /* Initialize numLeaves, root, parents, children
     * 1. construct Kd tree by using the center points of each edge
     * 2. initialize box hierarchy tree with the same structure of the kd tree
     */
    void Initialize_Hierarchy_Using_KD_Tree();
    int  Initialize_Hierarchy_Using_KD_Tree_Helper(KD_Tree_Node*);

    /* Initialize boxHierarchy
     * using the tree structure in parents, children
     * initialize boxHierarchy, boxHierarchy[i] is the bounding box of subtree
     * with root i
     * For some special case, it is necessary to enlarge each leaf bounding
     * box by enlargeThickness
     */
    void Update_Boxes(double enlargeThickness = 0.0);
    void Update_Leaf_Boxes(double enlargeThickness = 0.0);
    void Update_Nonleaf_Boxes();
};

//##################################################################################



#endif //INC_2DFLUIDSTRUCTURE_TRIANGLE_HIERACHY_H
