//
// Created by icme-huang on 8/29/16.
//

#ifndef INC_2DFLUIDSTRUCTURE_ARRAY_H
#define INC_2DFLUIDSTRUCTURE_ARRAY_H


#include <cassert>
#include <boost/type_traits/is_class.hpp>

template<class T>
class Array {
private:
    T* basePointer;


public:
    int len;
    int bufferSize;
    bool usingExternallyAllocatedPointer;

    Array():usingExternallyAllocatedPointer(false),basePointer(0),len(0), bufferSize(0)
    {}
    Array(int m):usingExternallyAllocatedPointer(false),len(m), bufferSize(m) {basePointer = new T[m];}




    int Length() const {return len;}

    ~Array()    {if(!usingExternallyAllocatedPointer)  delete[] basePointer; }

    const T& operator[](const int i) const {

        assert(i>=0 && i<len);
        return basePointer[i];
    }
    T& operator[](const int i) {

        assert(i>=0 && i<len);
        return basePointer[i];
    }

    void Ensure_Enough_Space(const int newLen,const bool copyExistingElements=true) {
        if(bufferSize < newLen)
            Resize_Helper(4*newLen/3+2,false,copyExistingElements);
    }



    int Append(const T& element) {
        Ensure_Enough_Space(len+1);
        len++;
        (*this)[len-1]= element;//because the new number is in id len-1
        return len;
    }
    void Remove_All() // if elements are non-primitive this may waste memory
    {len = 0;}

private:
    void Resize_Helper(const int newSize,const bool initializeNewElements=true,const bool copyExistingElements=true){

        T* p = new T[newSize];
        len = (len < newSize ? len: newSize);
        if(copyExistingElements) for(int i=0;i < len; i++) p[i]=basePointer[i];
        if(!boost::is_class<T>::value && initializeNewElements) for(int i = len;i < newSize;i++) p[i]=T();
        delete[] basePointer;
        basePointer=p;
        bufferSize = newSize;
    }

};


#endif //INC_2DFLUIDSTRUCTURE_ARRAY_H
