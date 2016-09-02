//
// Created by icme-huang on 8/29/16.
//

#ifndef INC_2DFLUIDSTRUCTURE_KD_TREE_NODE_H
#define INC_2DFLUIDSTRUCTURE_KD_TREE_NODE_H

#include "NonCopyable.h"
#include "Vector2D.h"

class KD_Tree_Node :public NonCopyable{
public:
    int splitAxis; // 0 means leaf
    double splitValue;
    int nodeIndex;
    KD_Tree_Node* left;
    KD_Tree_Node* right;

    KD_Tree_Node()
            :splitAxis(0),splitValue(0.),left(0),right(0)
    {}

    ~KD_Tree_Node()
    {}

};

class Pointer_Pool :public NonCopyable{

public:
    KD_Tree_Node* pools;
    int blockSize;
    int next;


    Pointer_Pool()
            :pools(0),blockSize(0.),next(0)
    {}

    void Allocate_Empty_Pool(int m){
        blockSize = m;
        pools = new KD_Tree_Node[m];
        next = 0;
    }

    KD_Tree_Node* New(){
        assert(next < blockSize);
        return pools + (next++);
    }


    ~Pointer_Pool() {if(blockSize) delete[] pools; }

};

#endif //INC_2DFLUIDSTRUCTURE_KD_TREE_NODE_H
