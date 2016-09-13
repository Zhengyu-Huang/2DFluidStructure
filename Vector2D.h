//
// Created by icme-huang on 8/28/16.
//

#ifndef INC_2DFLUIDSTRUCTURE_VECTOR2D_H
#define INC_2DFLUIDSTRUCTURE_VECTOR2D_H
#include <cassert>
#include<iostream>
template <class T>
inline T min(const T a, const T b)
{return (a<b)?a:b;}
template <class T>
inline T max(const T a, const T b)
{return (a>b)?a:b;}

template<class T>
class Vector2D {
public:
    T v[2];
    Vector2D() {v[0] = v[1] = 0; }
    Vector2D(T x, T y) {v[0] = x; v[1] = y; }
    Vector2D(const Vector2D<T> & vec) {v[0] = vec.v[0]; v[1] = vec.v[1]; }
    T &operator[](int i) { return v[i]; }
    T operator[](int i) const { return v[i]; }

    static Vector2D<T> Componentwise_Min(const Vector2D<T>& v1,const Vector2D<T>& v2)
    {return Vector2D<T>(min(v1.v[0],v2.v[0]),min(v1.v[1],v2.v[1]));}

    static Vector2D<T> Componentwise_Max(const Vector2D<T>& v1,const Vector2D<T>& v2)
    {return Vector2D<T>(max(v1.v[0],v2.v[0]),max(v1.v[1],v2.v[1]));}

    bool All_Less_Equal(const Vector2D<T>& v1) const
    {return v[0]<=v1[0] && v[1]<=v1[1];}

    bool All_Greater_Equal(const Vector2D<T>& v1) const
    {return v[0]>=v1[0] && v[1]>=v1[1];}

    void Get(T& element1,T& element2) const
    {element1=v[0];element2=v[1];}


    Vector2D<T> operator+(const Vector2D<T> & v2) const
    {
        Vector2D<T> sum;
        sum.v[0] = v[0] + v2[0];
        sum.v[1] = v[1] + v2[1];
        return sum;
    }

    Vector2D<T> operator+(const T&a) const
    {
        Vector2D<T> sum;
        sum.v[0] = v[0] + a;
        sum.v[1] = v[1] + a;
        return sum;
    }

    Vector2D<T> operator-(const Vector2D<T> & v2) const
    {
        Vector2D<T> sum;
        sum.v[0] = v[0] - v2[0];
        sum.v[1] = v[1] - v2[1];
        return sum;
    }

    Vector2D<T> operator-(const T&a) const
    {
        Vector2D<T> sum;
        sum.v[0] = v[0] - a;
        sum.v[1] = v[1] - a;
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

    Vector2D<T> &operator+=(const T&a)
    {
        v[0] = v[0] + a;
        v[1] = v[1] + a;
        return *this;
    }
    Vector2D<T> &operator-=(const T&a)
    {

        v[0] = v[0] - a;
        v[1] = v[1] - a;
        return *this;
    }
    friend std::ostream& operator<<(std::ostream& cout,const Vector2D<T>& vec){
        cout << vec.v[0] <<" " << vec.v[1] << std::endl;
        return cout;
    }




};


#endif //INC_2DFLUIDSTRUCTURE_VECTOR2D_H
