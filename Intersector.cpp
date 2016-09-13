//
// Created by icme-huang on 8/28/16.
//

#include "Intersector.h"

Intersector::Intersector(int inputXFluidNum, int inputYFluidNum, Vector2D<double> *inputXf, int inputNumXs,
                         double * inputXs, int inputNumEdge, int * inputEdge, double thickness)
        :xFluidNum(inputXFluidNum), yFluidNum(inputYFluidNum), Xf(inputXf), numXs(inputNumXs), numEdge(inputNumEdge),
         Xs(inputNumXs),Xs_n(inputNumXs),edgeSet(inputNumEdge) {
    //Initialize structure boundary info

    for(int i = 0; i<numEdge; i++) edgeSet[i] = Vector2D<int>(inputEdge[2*i],inputEdge[2*i+1]);
    //Initialize fluid info
    numXf = xFluidNum*yFluidNum;
    //Initialize intersector information

    status = new int[numXf];
    status_n = new int[numXf];
    sweptNode = new int[numXf];
    forwardMapping = new int[numXf];
    occlude = new bool[numXf];
    boxMax = new Vector2D<double>[numXf];
    boxMin = new Vector2D<double>[numXf];
    Recompute(inputXs, thickness);


}

Intersector::~Intersector() {
    delete[] status;
    delete[] status_n;
    delete[] sweptNode;
    delete[] forwardMapping;
    delete[] boxMax;
    delete[] boxMin;


}

void Intersector::Find_Node_Bounding_Boxes(){
    for(int i = 0; i < yFluidNum; i++)
        for(int j = 0; j < xFluidNum; j++){
            if(i != 0)/*bottom*/
                boxMin[i*xFluidNum + j][1] = Xf[(i-1)*xFluidNum + j][1];
            else
                boxMin[i*xFluidNum + j][1] = Xf[i*xFluidNum + j][1];
            if(i != yFluidNum - 1)/*top*/
                boxMax[i*xFluidNum + j][1] = Xf[(i+1)*xFluidNum + j][1];
            else
                boxMax[i*xFluidNum + j][1] = Xf[i*xFluidNum + j][1];
            if(j != 0)/*left*/
                boxMin[i*xFluidNum + j][0] = Xf[i*xFluidNum + j - 1][0];
            else
                boxMin[i*xFluidNum + j][0] = Xf[i*xFluidNum + j][0];
            if(j != xFluidNum - 1)/*right*/
                boxMax[i*xFluidNum + j][0] = Xf[i*xFluidNum + j + 1][0];
            else
                boxMax[i*xFluidNum + j][0] = Xf[i*xFluidNum + j][0];

        }

}
void Intersector::Recompute(double *inputXs, double thickness) {
    // for hasCloseTriangle
    bool* intersectOrNot = new bool[numXf];
    if(segmentHierarchy) delete segmentHierarchy;

    for(int i = 0; i<numXs; i++) Xs[i] = Vector2D<double>(inputXs[2*i],inputXs[2*i+1]);
    segmentHierarchy = new Segment_Hierarchy(Xs,edgeSet);


    Find_Node_Bounding_Boxes();
    Compute_Close_Segments(intersectOrNot, thickness);
    Find_Intersect_Result(intersectOrNot,thickness);


}



int Intersector::Compute_Close_Segments(bool * intersectOrNot,double thickness)
{
    int numCloseNodes = 0;
    int shrunkIndex = 0;
    for(int i=0;i<numXf;++i){
        Array<int> cand;
        Bounding_Box_2D box(boxMin[i],boxMax[i]);
        segmentHierarchy->Intersection_List(box,cand,thickness);

        intersectOrNot[i]= cand.Length() > 0 ? true:false;
        if(intersectOrNot[i]) {
            forwardMapping[i] = shrunkIndex;
            ++numCloseNodes;
            candidates.Append(cand);
            shrunkIndex++;

        }
    }
    return shrunkIndex;
}





int Intersector::Find_Intersect_Result(const bool * intersectOrNot,double thickness)
{

    int intersectedEdgeCount=0;
    /* consider horizontal edges */
    for(int i = 0; i < yFluidNum; i++)
        for(int j = 0; j < xFluidNum - 1; j++){
            // node pair (i*xFluidNum + j, i*xFluidNum + j + 1)
            if(intersectOrNot[i*xFluidNum+j] && intersectOrNot[i*xFluidNum + j + 1]) {
                IntersectionResult edgeRes;
                if(Intersect_Result_Helper(i*xFluidNum + j, i*xFluidNum + j + 1, edgeRes, thickness)){
                    xEdgeResult[i*xFluidNum  + j] = edgeRes;
                    intersectedEdgeCount++;

                }

            }



    }

    /* consider vertical edges */
    for(int i = 0; i < yFluidNum - 1; i++)
        for(int j = 0; j < xFluidNum; j++){
            // node pair (i*xFluidNum + j, (i+1)*xFluidNum + j)
            if(intersectOrNot[i*xFluidNum+j] && intersectOrNot[i*xFluidNum + j + 1]) {
                IntersectionResult edgeRes;
                if(Intersect_Result_Helper(i*xFluidNum + j, (i+1)*xFluidNum + j, edgeRes, thickness)) {
                    yEdgeResult[i * xFluidNum + j] = edgeRes;
                    intersectedEdgeCount++;
                }
            }



        }


    return intersectedEdgeCount;
}




bool Intersector::Intersect_Result_Helper(int startNode, int endNode, IntersectionResult& edgesRes,double thickness){

    Vector2D<double> startNodeCoord = Xf[startNode], endNodeCoord = Xf[endNode];
    double alpha = 0.0, beta = 0.0;
    int candidateStartNode = forwardMapping[startNode], candidateEndNode = forwardMapping[endNode];
    for (int j = 0; j < candidates[candidateStartNode].Length(); ++j)
        for (int k = 0; k < candidates[candidateEndNode].Length(); ++k)
            if (candidates[candidateStartNode][j] == candidates[candidateEndNode][k]) {
                int segmentID = candidates[candidateStartNode][j];
                Vector2D<double> startCandidateNodeCoord = Xs[edgeSet[segmentID][0]], endCandidateNodeCoord = Xs[edgeSet[segmentID][1]];
                if(Line_Interaction(&alpha, &beta, startNodeCoord,endNodeCoord,startCandidateNodeCoord,endCandidateNodeCoord,thickness)) {
                    if (alpha < thickness)
                        occlude[startNode] = true;
                    if (alpha > 1-thickness)
                        occlude[endNode] = true;
                    if(edgesRes.structureSegmentID[0] < 0) {
                        edgesRes.structureSegmentCoord[0] = beta;
                        edgesRes.structureSegmentID[0] = segmentID;
                        edgesRes.edgeCoord[0] = alpha;
                    }else if(alpha < edgesRes.edgeCoord[0]){
                        edgesRes.structureSegmentCoord[0] = beta;
                        edgesRes.structureSegmentID[0] = segmentID;
                        edgesRes.edgeCoord[0] = alpha;

                    }

                    if(edgesRes.structureSegmentID[1] < 0) {
                        edgesRes.structureSegmentCoord[1] = beta;
                        edgesRes.structureSegmentID[1] = segmentID;
                        edgesRes.edgeCoord[1] = alpha;
                    }else if(alpha > edgesRes.edgeCoord[1]){
                        edgesRes.structureSegmentCoord[1] = beta;
                        edgesRes.structureSegmentID[1] = segmentID;
                        edgesRes.edgeCoord[1] = alpha;

                    }
                    edgesRes.intersectTimes++;


                } else break;


            }
    return edgesRes.intersectTimes;
}



bool Intersector::Line_Interaction(double * alpha, double * beta, const Vector2D<double> &v1, const Vector2D<double> &v2,
                      const Vector2D<double> &u1, const Vector2D<double> &u2,double thickness) const {
   double v21x = v2[0] - v1[0],v21y = v2[1] - v1[1], u21x = u2[0] - u1[0], u21y = u2[1] - u1[1];
   double delta = v21x*u21y - u21x*v21y;
   *alpha = (u21y * (u1[0] - v1[0]) - u21x*(u1[1] - v1[1]))/delta;
   *beta = -(-v21y * (u1[0] - v1[0]) + v21x*(u1[1] - v1[1]))/delta;
    if(*alpha > -thickness && *alpha < 1.0+thickness && *beta > -thickness && *beta < 1.0 + thickness)
        return true;
    else
        return false;

}


void Intersector::Print_Intersector_Info(){
    for (int i = 0; i < numXf; i++)
        if(occlude[i])
            std::cout <<"node " << i << " is occluded" << std::endl;
    ;
    std::cout << "HORIZONTAL EDGE" << std::endl;
    for(std::map<int,IntersectionResult>::iterator it = xEdgeResult.begin(); it!=xEdgeResult.end(); ++it){
        std::cout<<"Edge: "<< it->first << " intersects: "<<it->second.structureSegmentID[0] <<" at its " <<it->second.edgeCoord[0]
        << " at " << it->second.structureSegmentCoord[0] <<std::endl;
        if(it->second.intersectTimes > 1)
        std::cout<<"Edge: "<< it->first << " intersects: "<<it->second.structureSegmentID[1] <<" at its " <<it->second.edgeCoord[1]
        << " at " << it->second.structureSegmentCoord[1] <<std::endl;
    }

    std::cout << "VERTICAL EDGE" << std::endl;
    for(std::map<int,IntersectionResult>::iterator it = yEdgeResult.begin(); it!=yEdgeResult.end(); ++it){
        std::cout<<"Edge: "<< it->first << " intersects: "<<it->second.structureSegmentID[0] <<" at its " <<it->second.edgeCoord[0]
        << " at " << it->second.structureSegmentCoord[0] <<std::endl;
        if(it->second.intersectTimes > 1)
            std::cout<<"Edge: "<< it->first << " intersects: "<<it->second.structureSegmentID[1] <<" at its " <<it->second.edgeCoord[1]
            << " at " << it->second.structureSegmentCoord[1] <<std::endl;
    }

}