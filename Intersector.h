//
// Created by icme-huang on 8/28/16.
//




#ifndef INC_2DFLUIDSTRUCTURE_INTERSECTOR_H
#define INC_2DFLUIDSTRUCTURE_INTERSECTOR_H
#include "Vector2D.h"
#include "Array.h"
#include "Segment_Hierarchy.h"
#include "map"


/* for each fluid edge, save these structure edges, closest one to the start node and the end node */
struct IntersectionResult {
    int intersectTimes;        // intersection time with structure boundary edges
    int structureSegmentID[2]; //  intersected structure edge IDs
    double structureSegmentCoord[2]; // structure edge intersecting position, Intersect at alpha*startNode + (1-alpha)*endNode
    double edgeCoord[2]; //  fluid edge intersecting poition, Intersect at alpha*startNode + (1-alpha)*endNode
    IntersectionResult():structureSegmentID{-1,-1},intersectTimes(0){}
    ~IntersectionResult(){}
};

class Segment_Hierachy;
//should be in the problem descriptor
class Intersector {
public:


    int xFluidNum;  //number of fluid nodes in each row
    int yFluidNum; //number of fluid nodes in each column
    int numXf;     //number of fluid nodes(numXf = xFluidNum*yFluidNum)
    Vector2D<double> *Xf; //pointer to fluid node coordinate, the fluid node is labeled from left to right, from bottom to top,
    // node at i row, j column is labeled as i*xFluidNum + j
    Vector2D<double> *boxMax, *boxMin; //Bounding box of each node and its neighbors, bounding box is boxMin[i], boxMax[i]


    int numXs;                          //number of structure boundary nodes
    Array<Vector2D<double>> Xs;         //current structure boundary position
    Array<Vector2D<double>> Xs_n;       //structure boundary position of last step
    int numEdge;                        //number of structure boundary edges
    Array<Vector2D<int>> edgeSet;       //node pairs of the structure boundary edge, if the structure is 2D, the edge should be counterclockwise.

    /* Intersector information */
    Segment_Hierarchy * segmentHierarchy; //Segment hierarchy of structure boundary edges
    int * status;                         //  current fluid status, 1,2,3.. means active with different fluid id, 0 means inactive
    int * status_n;                       //  fluid state of last step
    bool * sweptNode;                      //  whether fluid is swept by the structure during last time step
    bool * occlude;                       //  true means the node is in the structure boundary with some tolerance, otherwise false
    bool * occlude_n;
    int * forwardMapping;                   //  if the bounding box of fluid node i, lazy intersects the bounding boxes of structure boundary edges,
    int * backwardMapping;                 //  node i -> forwardMapping[i], and backwardMapping[forwardMapping[i]] = i
    Array<Array<int>> candidates;         // and these structure boundary edge labels are in candidates[forwardMapping[i]] array

    bool* xEdgeIsIntersect;               // each horizontal edge has the same label as its left node,if it intersects structure boundary, it is true
    bool* yEdgeIsIntersect;                // each vertical edge has the same label as its lower node,if it intersects structure boundary, it is true
    std::map<int,IntersectionResult>  xEdgeResult; // for each horizontal edge, if it intersects structure boundary, it is maps to IntersectionResult,

    std::map<int,IntersectionResult>  yEdgeResult; // for each vertical edge, if it intersects structure boundary, it is maps to IntersectionResult,
    // each vertical edge has the same label as its lower node


    /* numXs : number of structure boundary nodes
     * Xs : coordinate of structure boudary nodes Points(i) = (Xs[2*i],Xs[2*i+1])
     * numEdge : number of structure boudary edge
     * ij : endpoint pairs Edge(i) = (Points(2*i), Points(2*i+1))
     * thickness: intersection detection tolerance
     * Initialize Structure info, Intersector info, Fluid info
     */
    Intersector(int xFluidNum, int yFluidNum, Vector2D<double> *Xf, int numXs, double *Xs, int numEdge, int* ij, double thickness);
    ~Intersector();


    /* Xs : coordinate of structure boudary nodes Points(i) = (Xs[2*i],Xs[2*i+1])
    * thickness: intersection detection tolerance in real length
    * reset structure boundary position, thickness and compute intersector*/
    void Recompute(double *Xs, double thickness);

    /*Print intersector infomation*/
    void Print_Intersector_Info();


private:
    /* update array boxMax and boxMin
     * for bounding box of each fluid node*/
    void Find_Node_Bounding_Boxes();

    /* update Array<Array<int>> candidates, and int * forwardMapping, and int* intersectOrNot
     * for these lazy-intersected fluid nodes, whose bounding box intersects some structure boundary edge bounding box),
     * intersectOrNot[i] = true, otherwise false
     * return the number of these lazy-intersected nodes */
    int  Compute_Close_Segments(bool * intersectOrNot, double thickness);

    /* update std::map<int,IntersectionResult> xEdgeResult and std::map<int,IntersectionResult> yEdgeResult
     * and bool* occlude
     * input intersectOrNot is the result in int  Compute_Close_Segments(bool * intersectOrNot, double thickness);
     * return the number of intersected fluid edges */
    int Find_Intersect_Result(const bool * intersectOrNot,double thickness);

    /* update edgesRes of fluid edge (startNode, endNode) and bool* occlude*/
    bool Intersect_Result_Helper(int startNode, int endNode, IntersectionResult& edgesRes,double thickness);

    /* check the distance between p and segment v1v2
     * if distance < thckness return true
     * the closest point in the segment is v1 + *beta*(v2 -v1)
     * */
    bool Point_Inside_Segment(const Vector2D<double> &p, const Vector2D<double> &v1, const Vector2D<double> &v2,double thickness,double *beta) const;

    /* compute the intersection point of segment v1v2 and u1u2
     * update alpha and beta
     * the intersection is v1+alpha(v2-v1),  and u1 + beta(u2-u1)
     * return true if two segments intersect, the distance between two segments are less than thickness
     * */
    bool Ray_Interaction(double * alpha, double * beta, const Vector2D<double> &v1, const Vector2D<double> &v2,
                          const Vector2D<double> &u1, const Vector2D<double> &u2,double thickness) const;

    /* update bool* sweptNode
     * if node i is swept by Xs_n and Xs, namely node i is is quad Xs[j] Xs[j+1] and Xs_n[j+1] Xs_n[j]
     * swept[i] = ture
     * */
    void Compute_Swetp_Nodes(const int intersectNum, double thickness);
    /* return whether node p is in quad v1 v2 v2_n v1_n with tolerance thickness
     * if the distance of the node and the quad is < thickness
     * regard the node as in the quad
     * */
    bool Inside_Quad(const Vector2D<double> &p, const Vector2D<double> &v1, const Vector2D<double> &v2,
                    const Vector2D<double> &v1_n, const Vector2D<double> &v2_n,double thickness);



    /* update *int status
     * for these nodes that are not swept, status[i] = status_n[i]
     * update the status of ther other nodes by flood fill method
     * */
    void Find_Status();
    /* update *int status
     * Assume node(i,0) is in fluid(active), and
     * update the status of ther other nodes by flood fill method
     * */
    void Find_Status_Using_FloodFill();



};


#endif //INC_2DFLUIDSTRUCTURE_INTERSECTOR_H
