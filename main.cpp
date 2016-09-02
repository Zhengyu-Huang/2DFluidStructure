#include <iostream>
//#include <Eigen/Dense>
//using Eigen::MatrixXd;
#include"Segment_Hierarchy.h"
using namespace std;
#include <vector>
#include "Vector2D.h"
#include "Array.h"
#include<math.h>
void Segment_Hierarchy_Test() {
    const double PI = 3.1415926535898;
    std::cout << "Segment_Hierarchy_Test" << std::endl;

    int segmentNumber = 4;
    Array<Vector2D<double>> points(segmentNumber);
    Array<Vector2D<int>> edgeSet(segmentNumber);
    for(int i = 0; i < segmentNumber; i++){
        points[i] = Vector2D<double>(cos(i*2.0*PI/segmentNumber),sin(i*2.0*PI/segmentNumber));
        edgeSet[i] = Vector2D<int>(i, (i+1)%segmentNumber);

    }

    Segment_Hierarchy segmentHierarchy(points,edgeSet);
    segmentHierarchy.Print_Info();






}

int main() {
    //KD_Tree_Test();
    Segment_Hierarchy_Test();




    return 0;
}
