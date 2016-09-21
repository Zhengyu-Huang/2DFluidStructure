//
// Created by Zhengyu on 9/17/16.
//

#include "Fluid_Solver.h"
/*
template <int dim>

int Fluid_Solver::Solve_NonLinear_System(MatrixVec<dim> &U, int) {
    DistSVec<double,dim> Ubc(this->getVecInfo());

    commonPart(U);
    if(RK4)     solveNLAllRK4(U,t0,Ubc);
    else if(FE) solveNLAllFE(U,t0,Ubc);
    else        solveNLAllRK2(U,t0,Ubc);

    this->updateBoundaryExternalState();
}

void solveNLSystemOneBlock(DistSVec<double,dim> &U);
//  void solveNLSystemTwoBlocks(DistSVec<double,dim> &U);


void commonPart(DistSVec<double,dim> &U); // Common part to the two following functions.

tmeplate<int dim>
void Fluid_Solver::Solve_Using_RK2(MatrixVec<dim> &U, double t0, MatrixVec<dim> &Ubc) {
    computeRKUpdate(U, k1, 1);
    this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);
    U0 = U - k1;


    this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
    this->checkSolution(U0);




    computeRKUpdate(U0, k2, 1);

    this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);

    U = U - 1.0 / 2.0 * (k1 + k2);

    this->spaceOp->applyExtrapolationToSolutionVector(U, Ubc);

    this->spaceOp->applyBCsToSolutionVector(U, this->distLSS);

    this->checkSolution(U);

}

void computeRKUpdate(MatrixVec<dim>& Ulocal, MatrixVec<dim>& dU, int it);

void computeRKUpdateHH(DistSVec<double,dim>& Ulocal,
                       DistVec<double>& dHH) ;
};
*/