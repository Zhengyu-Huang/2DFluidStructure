//
// Created by icme-huang on 8/29/16.
//

#ifndef INC_2DFLUIDSTRUCTURE_KD_TREE_NODE_H
#define INC_2DFLUIDSTRUCTURE_KD_TREE_NODE_H

#include "NonCopyable.h"
#include "Vector2D.h"

class KD_Tree_Node :public NonCopyable{
public:
    int splitAxis;     // -1 means leaf, 0 meas x axis, 1 means y axis
    double splitValue;
    int nodeIndex;
    KD_Tree_Node* left;//pointer to left child
    KD_Tree_Node* right;//pointer to right child

    KD_Tree_Node()
            :splitAxis(0),splitValue(0.),left(0),right(0)
    {}

    ~KD_Tree_Node()
    {}

};



/*  Store all these nodes of KD Tree sequentially, in preorder(DLR) */
class Pointer_Pool :public NonCopyable{

public:
    KD_Tree_Node* pools; // pointer to all KD tree nodes array
    int blockSize;       // the size of pools
    int next;            // the number of active nodes(pools+next is the address of next available position)


    Pointer_Pool()
            :pools(0),blockSize(0.),next(0)
    {}
    /*allocate new pools with size m */
    void Allocate_Empty_Pool(int m){
        blockSize = m;
        pools = new KD_Tree_Node[m];
        next = 0;
    }
    /* return the next available position in pools to store next node */
    KD_Tree_Node* New(){
        assert(next < blockSize);
        return pools + (next++);
    }


    ~Pointer_Pool() {if(blockSize) delete[] pools; }

};

#endif //INC_2DFLUIDSTRUCTURE_KD_TREE_NODE_H
