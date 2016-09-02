#include <iostream>
//#include <Eigen/Dense>
//using Eigen::MatrixXd;
#include"KD_Tree.h"
using namespace std;
#include <vector>
#include "Vector2D.h"
#include "Array.h"

void KD_Tree_Test(){
    std::cout <<"Start KD Tree Test" << std::endl;

    int segmentNumber = 7;
    Array<Vector2D<double>> centroids(segmentNumber);

    centroids[0] = Vector2D<double>(2.0,3.0);
    centroids[1] = Vector2D<double>(4.0,7.0);
    centroids[2] = Vector2D<double>(5.0,4.0);
    centroids[3] = Vector2D<double>(7.0,2.0);
    centroids[4] = Vector2D<double>(8.0,1.0);
    centroids[5] = Vector2D<double>(9.0,3.0);
    centroids[6] = Vector2D<double>(9.0,6.0);
    KD_Tree kdTree;

    kdTree.Create_KD_Tree(centroids);
    kdTree.Print_Tree_Info();


}
int main() {
    KD_Tree_Test();




    return 0;
}
