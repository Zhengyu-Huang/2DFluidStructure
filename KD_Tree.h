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
    KD_Tree_Node* rootNode;

    Pointer_Pool pool;// Store all these pointer of KD Tree sequentially



    KD_Tree() :rootNode(0)
    {}

    ~KD_Tree() {}


public:

//#####################################################################

    void Create_KD_Tree(const Array<Vector2D<double>> &pointsToBalance);
    void Print_Tree_Info();
private:
    void Create_KD_Tree_Helper(KD_Tree_Node* cell,const int firstIndex,const int lastIndex,
                                        const Array<Vector2D<double>> &points,int* permutationArray,const Bounding_Box_2D &box);

    int Choose_Partition_Index(const int firstIndex,const int lastIndex);

    void Median_Split(const int partitionIndex,const int firstIndex,const int lastIndex,const Array<Vector2D<double>> &points,
                               int* permutationArray,const int axis);


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

    void Print_Tree_Info_Helper(KD_Tree_Node* Node);

//#####################################################################
};


#endif //INC_2DFLUIDSTRUCTURE_KD_TREE_H
