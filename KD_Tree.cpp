//
// Created by icme-huang on 8/29/16.
//
#include <iostream>
#include "KD_Tree.h"
#include "Integer_Log.h"

/* Create left balanced(left subtree is always full) KD Tree without using internal nodes
 */

void KD_Tree::Create_KD_Tree(const Array<Vector2D<double>> &pointsToBalance)
{

    assert(pointsToBalance.Length() >0);

    int len = pointsToBalance.Length();

    pool.Allocate_Empty_Pool(2*len-1);

    rootNode=pool.New();

    int* permutationArray = new int[len];


    for(int i = 0; i < pointsToBalance.Length(); i++) permutationArray[i] = i;


    Bounding_Box_2D box = Bounding_Box_2D(pointsToBalance);

    Create_KD_Tree_Helper(rootNode,0,len-1,pointsToBalance,permutationArray,box);

}

void KD_Tree::Create_KD_Tree_Helper(KD_Tree_Node* cell,const int firstIndex,const int lastIndex,
                                    const Array<Vector2D<double>> &points,int* permutationArray,const Bounding_Box_2D &box) {
    if (lastIndex == firstIndex) {
        cell->splitAxis = -1; //for leaf cell, the split Axis is -1
        cell->nodeIndex = permutationArray[firstIndex];
        return;
    }
    int partitionIndex = Choose_Partition_Index(firstIndex, lastIndex);

    cell->splitAxis = box.Max_Length_Axis();

    Median_Split(partitionIndex, firstIndex, lastIndex, points, permutationArray, cell->splitAxis);
    cell->splitValue = points[permutationArray[partitionIndex]][cell->splitAxis];
    cell->nodeIndex = -1;

    assert(partitionIndex > firstIndex);
    cell->left = pool.New();
    Bounding_Box_2D leftBox(box);
    leftBox.maxCorner[cell->splitAxis] = cell->splitValue;
    Create_KD_Tree_Helper(cell->left, firstIndex, partitionIndex - 1, points, permutationArray,
                                         leftBox);
    if (partitionIndex <= lastIndex) {
        cell->right = pool.New();
        Bounding_Box_2D rightBox(box);
        rightBox.minCorner[cell->splitAxis] = cell->splitValue;
        Create_KD_Tree_Helper(cell->right, partitionIndex, lastIndex, points, permutationArray,
                                             rightBox);
    }
}

int KD_Tree::Choose_Partition_Index(const int firstIndex,const int lastIndex){
    int elements=lastIndex-firstIndex+1,filledLeftSubtreeElements=(1<<Integer_Log(elements));
    return firstIndex+ std::min(elements-filledLeftSubtreeElements/2,filledLeftSubtreeElements); // this index goes to the right subtree

}

void KD_Tree::Median_Split(const int partitionIndex,const int firstIndex,const int lastIndex,const Array<Vector2D<double>> &points,
                           int* permutationArray,const int axis){

    Sort_Helper_Less_Equal sortHelperLessEqual(&points,axis);

    std::stable_sort(permutationArray + firstIndex, permutationArray + lastIndex + 1,sortHelperLessEqual);

}
/*preorder traversal
 * Data LeftChild RightChild
 */
void KD_Tree::Print_Tree_Info(){

    KD_Tree_Node * Node;
    for(int i =0; i < pool.next; i++){
        Node = pool.pools + i;
        std::cout << "splitAxis is " << Node->splitAxis << " splitValue is " << Node->splitValue
        << " nodeIndes is " << Node->nodeIndex <<std::endl;

    }

}

