//
// Created by icme-huang on 8/28/16.
//

#ifndef INC_2DFLUIDSTRUCTURE_VECTOR2D_H
#define INC_2DFLUIDSTRUCTURE_VECTOR2D_H
#include <cassert>
template<class T>
class Vector2D {
public:
    T v[2];
    Vector2D() {v[0] = v[1] = 0; }
    Vector2D(T x, T y) {v[0] = x; v[1] = y; }
    Vector2D(const Vector2D<T> & vec) {v[0] = vec.v[0]; v[1] = vec.v[1]; }
    T &operator[](int i) { return v[i]; }
    T operator[](int i) const { return v[i]; }

    Vector2D<T> operator+(const Vector2D<T> & v2) const
    {
        Vector2D<T> sum;
        sum.v[0] = v[0] + v2[0];
        sum.v[1] = v[1] + v2[1];
        return sum;
    }


    Vector2D<T> operator/(const T cst)
    {
        Vector2D<T> vNew;
        vNew[0] = v[0] / cst;
        vNew[1] = v[1] / cst;
        return vNew;
    }

    Vector2D<T>& operator=(const Vector2D<T> vec)
    {
        for(int i=0;i<2;i++) v[i]=vec.v[i];
        return *this;
    }




};


#endif //INC_2DFLUIDSTRUCTURE_VECTOR2D_H
