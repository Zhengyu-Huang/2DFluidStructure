//
// Created by Zhengyu on 9/18/16.
//

#ifndef INC_2DFLUIDSTRUCTURE_MATRIXVEC_H
#define INC_2DFLUIDSTRUCTURE_MATRIXVEC_H
#include "assert.h"


template<int dim> class MatrixVec;
/*************************************************************************************************************************/
template<class T>
class Expr {
public:
    int len;
    T x;
    Expr(T v) : x(v) { len = x.size(); }
    Expr(T v, int l) : x(v) { len = l; }
    double operator[] (int i) const { return x[i]; }
    int size() const { return len; }

};
template<class T1, class T2>
class Diff {
    T1 a;
    T2 b;
    int len;
public:
    Diff(T1 aa, T2 bb, int l) : a(aa), b(bb) { len = l; }
    double operator[](int i) const { return a[i]-b[i]; }
    int size() const { return len; }
};

template<class T1, class T2>
class Sum {
    T1 a;
    T2 b;
    int len;
public:
    Sum(T1 aa, T2 bb, int l) : a(aa), b(bb) { len = l; }
    double operator[](int i) const { return a[i]+b[i]; }
    int size() const { return len; }

};
template<class T>
class OuterProd {
    double y;
    T a;
    int len;
public:
    OuterProd(double yy, T aa, int l) : y(yy), a(aa) { len = l; }
    double operator[](int i) const { return y*a[i]; }
    int size() const { return len; }

};

/* Operator overloading for - */
template<int dim>
inline Expr<Diff<double *, double *>> operator-(const MatrixVec<dim> &v1, const MatrixVec<dim> &v2)
{
    return Expr<Diff<double *, double *>> ( Diff<double *, double *>(reinterpret_cast<double *>(v1.v), reinterpret_cast<double *>(v2.v), v1.size()*dim) );
}

template<class T, int dim>
inline Expr<Diff<T, double*>> operator-(const Expr<T> &x, const MatrixVec<dim>  &v)
{
    return Expr<Diff<T, double *>>( Diff<T, double *>(x.x, reinterpret_cast<double *>(v.v), x.size()) );

}

template<class T, int dim>
inline Expr<Diff<double *, T>> operator-(const MatrixVec<dim>  &v, const Expr<T> &x)
{
    return Expr<Diff<double *, T>>( Diff<double *, T>(x.x, reinterpret_cast<double *>(v.v), x.size()) );

}

template<class T1, class T2>
inline Expr<Diff<T1, T2>> operator-(const Expr<T1> &x,  const Expr<T2> &y)
{
    return Expr<Diff<T1, T2>>( Diff<T1, T2>(x.x, y.x, x.size()) );

}


/* Operator overloading for + */
template<int dim>
inline Expr<Sum<double *, double *>> operator+(const MatrixVec<dim> &v1, const MatrixVec<dim> &v2)
{
    return Expr<Sum<double *, double *>> ( Sum<double *, double *>(reinterpret_cast<double *>(v1.v), reinterpret_cast<double *>(v2.v), v1.size()*dim) );
}

template<class T, int dim>
inline Expr<Sum<T, double*>> operator+(const Expr<T> &x, const MatrixVec<dim>  &v)
{
    return Expr<Sum<T, double *>>( Sum<T, double *>(x.x, reinterpret_cast<double *>(v.v), x.size()) );
}

template<class T, int dim>
inline Expr<Sum<double *, T>> operator+(const MatrixVec<dim>  &v, const Expr<T> &x)
{
    return Expr<Sum<double *, T>>( Sum<double *, T>(x.x, reinterpret_cast<double *>(v.v), x.size()) );
}

template<class T1, class T2>
inline Expr<Sum<T1, T2>> operator+(const Expr<T1> &x,  const Expr<T2> &y)
{
    return Expr<Sum<T1, T2>>( Sum<T1, T2>(x.x, y.x, x.size()) );

}

/*operator overloading for * */
template<int dim>
inline Expr<OuterProd<double *>> operator*(double y, const MatrixVec<dim> &v)
{
    return Expr<OuterProd<double *>> ( OuterProd<double *>(y, reinterpret_cast<double *>(v.v), v.size()*dim) );
}

template<int dim>
inline Expr<OuterProd<double *>> operator*( const MatrixVec<dim> &v,double y)
{
    return Expr<OuterProd<double *>> ( OuterProd<double *>(y, reinterpret_cast<double *>(v.v), v.size()*dim) );
}

template<class T>
inline Expr<OuterProd<T>> operator*(double y, const Expr<T> &x)
{

    return Expr<OuterProd<T>>( OuterProd<T>(y, x.x, x.size()) );

}

template<class T>
inline Expr<OuterProd<T>> operator*(const Expr<T> &x,double y)
{

    return Expr<OuterProd<T>>( OuterProd<T>(y, x.x, x.size()) );

}


/******************************************************************************/
template<int dim>
class MatrixVec {
public:
    int xNum;  //number of fluid nodes in each row
    int yNum; //number of fluid nodes in each column
    int numXf;     //number of fluid nodes(numXf = xFluidNum*yFluidNum)
    double (*v)[dim];
    bool usingExternallyAllocatedPointer;

    MatrixVec(int xFluidNum, int yFluidNum, double (*)[dim] = 0);

    ~MatrixVec();

    int numCol() const { return xNum; }
    int numRow() const { return yNum; }
    int size() const { return numXf; }

    double* operator()(int i, int j) {
        assert(i >= 0 && i < xNum && j >= 0 && j < yNum);
        return v[i * xNum + j];
    }

    MatrixVec<dim> &operator=(const MatrixVec<dim> &);
    template<class T>
    MatrixVec<dim> &operator=(const Expr<T> &expr);

};

template<int dim>
MatrixVec<dim>::MatrixVec(int xFluidNum, int yFluidNum, double (*vv)[dim]) :xNum(xFluidNum),yNum(yFluidNum),numXf(xFluidNum*yFluidNum) {
    if (vv) {
        usingExternallyAllocatedPointer = true;
        v = vv;
    } else {
        usingExternallyAllocatedPointer = false;
        v = new double[numXf][dim];

    }
}

template<int dim>
MatrixVec<dim>::~MatrixVec() {
    if (!usingExternallyAllocatedPointer && v) delete [] v;
}



template<int dim>
MatrixVec<dim> & MatrixVec<dim>::operator=(const MatrixVec<dim> & y) {
    const double *yy = reinterpret_cast<double *>(y.v);
    double *vv = reinterpret_cast<double *>(this->v);
//#pragma omp parallel for
    for (int i = 0; i < numXf * dim; ++i) vv[i] = yy[i];
    return *this;
}

template<int dim>
template<class T>
inline
MatrixVec<dim> & MatrixVec<dim>::operator=(const Expr<T> &expr)
{

    const T &x = expr.x;
    double *vv = reinterpret_cast<double *>(this->v);

//#pragma omp parallel for
    for (int i = 0; i < numXf * dim; ++i) vv[i] = x[i];
    return *this;

}



//------------------------------------------------------------------------------



#endif //INC_2DFLUIDSTRUCTURE_MATRIXVEC_H
