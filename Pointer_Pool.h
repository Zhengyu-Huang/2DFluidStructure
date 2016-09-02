//
// Created by icme-huang on 8/29/16.
//

#ifndef INC_2DFLUIDSTRUCTURE_POINTER_POOL_H
#define INC_2DFLUIDSTRUCTURE_POINTER_POOL_H
#include "NonCopyable.h"
#include "Array.h"
#include "KD_Tree_Node.h"


class Pointer_Pool:public NonCopyable{
    private:
        Array<KD_Tree_Node>* pools;
        int blockSize;
        int current;

    public:

        Pointer_Pool(int n):blockSize(n),current(0){
        }

        ~Pointer_Pool() {
            Delete_All();
        }

        KD_Tree_Node* New() // memory will be freed by Delete_All (do not use delete!)
        {
            assert(current < blockSize - 1);
            current++;
            return pools[pools.len]+(current++);
        }

        void Delete_All()
        {
            for(int i=1;i<=pools.len;i++)
                delete[] pools[i];
            pools.Remove_All();
            current=blockSize;
        }

//#####################################################################
    };




#endif //INC_2DFLUIDSTRUCTURE_POINTER_POOL_H
