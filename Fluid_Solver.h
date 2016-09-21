//
// Created by Zhengyu on 9/17/16.
//

#ifndef INC_2DFLUIDSTRUCTURE_FLUID_SOLVER_H
#define INC_2DFLUIDSTRUCTURE_FLUID_SOLVER_H

#include "intersector.h"
#include "MatrixVec.h"
template<int dim>
class Fluid_Solver {
private:
    int xFluidNum;  //number of fluid nodes in each row
    int yFluidNum; //number of fluid nodes in each column
    int numXf;     //number of fluid nodes(numXf = xFluidNum*yFluidNum)
    Vector2D<double> *Xf;
    double *controlVolume;
    Intersector intersector;


    MatrixVec<dim> U0;
    MatrixVec<dim> k1;
    MatrixVec<dim> k2;


    bool RK2;
    bool FE;

public:

    Fluid_Solver();
    ~Fluid_Solver();

    int Solve_NonLinear_System(MatrixVec<dim> &U, int);

private:


    void Update_Intersector(MatrixVec<dim> &U); // Common part to the two following functions.

    void Solve_Using_RK2(MatrixVec<dim> &U, double t0, MatrixVec<dim> &Ubc);


    void computeRKUpdate(MatrixVec<dim>& Ulocal, MatrixVec<dim>& dU, int it);

    //void computeRKUpdateHH(DistSVec<double,dim>& Ulocal,DistVec<double>& dHH) ;
};




#endif //INC_2DFLUIDSTRUCTURE_FLUID_SOLVER_H
