//
// Created by icme-huang on 8/29/16.
//


#ifndef INC_2DFLUIDSTRUCTURE_KD_TREE_H
#define INC_2DFLUIDSTRUCTURE_KD_TREE_H

#include "KD_Tree_Node.h"
#include "Bounding_Box_2D.h"
#include <vector>

class KD_Tree
{
public:
    KD_Tree_Node* rootNode;  // A pointer to KD tree root node

    Pointer_Pool pool;       // Store all these nodes of KD Tree sequentially, in preorder(DLR)



    KD_Tree() :rootNode(0)
    {}

    ~KD_Tree() {}


public:

    /* Create left balanced KD tree by an array of points
     * left balanced KD tree has full left subtree, 2^k leaves
     */
    void Create_KD_Tree(const Array<Vector2D<double>> &pointsToBalance);

    /* Print information of each node in preorder */
    void Print_Tree_Info();
private:

    /* Create Sub KD tree of a subset of the points array, point[permutationArray[firstIndex]], point[permutationArray[firstIndex+1]]
     * ......point[permutationArray[lastIndex]]
     * box: bounding box of the targe points(this sub KD tree)
     * points: points array
     * cell: the pointer to the root of the KD tree
     */
    void Create_KD_Tree_Helper(KD_Tree_Node* cell,const int firstIndex,const int lastIndex,
                                        const Array<Vector2D<double>> &points,int* permutationArray,const Bounding_Box_2D &box);

    /* The left balances subKdTree contains node Id from firstIndex to the lastIndex
     * return the partition index of the first nodes in the right subtree of the subKDTree
     */
    int Choose_Partition_Index(const int firstIndex,const int lastIndex);

    /* Stably sort the the subarray of permutationArray from firstIndex to the lastIndex
     * make sure any firstIndex <= i < j <= lastIndes
     * points[permutationArray[i]][axis] <= points[permutationArray[j]][axis]
     */
    void Median_Split(const int partitionIndex,const int firstIndex,const int lastIndex,const Array<Vector2D<double>> &points,
                               int* permutationArray,const int axis);

    /* class used for stably sort in Median_Split function
     * return points[i][axis] <= points[j][axis]
     */
    class Sort_Helper_Less_Equal {
    public:
        Array<Vector2D<double>> const * points;
        int axis;
        double splitValue;
        Sort_Helper_Less_Equal(const Array<Vector2D<double>>* p,int a){
            points=p;
            axis=a;
        }
        bool operator()(int i,int j){return (*points)[i][axis]< (*points)[j][axis];}
        ~Sort_Helper_Less_Equal(){}
    };

    /* Print the information of sub KD tree, with root pointer node*/
    void Print_Tree_Info_Helper(KD_Tree_Node* node);

};

/*
void KD_Tree_Test() {
    std::cout << "Start KD Tree Test" << std::endl;

    int segmentNumber = 7;
    Array<Vector2D<double>> centroids(segmentNumber);

    centroids[0] = Vector2D<double>(2.0, 3.0);
    centroids[1] = Vector2D<double>(4.0, 7.0);
    centroids[2] = Vector2D<double>(5.0, 4.0);
    centroids[3] = Vector2D<double>(7.0, 2.0);
    centroids[4] = Vector2D<double>(8.0, 1.0);
    centroids[5] = Vector2D<double>(9.0, 3.0);
    centroids[6] = Vector2D<double>(9.0, 6.0);
    KD_Tree kdTree;

    kdTree.Create_KD_Tree(centroids);
    kdTree.Print_Tree_Info();
}
 */
#endif //INC_2DFLUIDSTRUCTURE_KD_TREE_H
