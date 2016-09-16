//
// Created by icme-huang on 8/28/16.
//

#include "Intersector.h"
#include "queue"
#include "math.h"

Intersector::Intersector(int inputXFluidNum, int inputYFluidNum, Vector2D<double> *inputXf, int inputNumXs,
                         double * inputXs, int inputNumEdge, int * inputEdge, double thickness)
        :xFluidNum(inputXFluidNum), yFluidNum(inputYFluidNum), Xf(inputXf), numXs(inputNumXs), numEdge(inputNumEdge),
         Xs(inputNumXs),Xs_n(inputNumXs),edgeSet(inputNumEdge) {
    //Initialize structure boundary info
    for(int i = 0; i<numXs; i++) Xs[i] = Vector2D<double>(inputXs[2*i],inputXs[2*i+1]);
    for(int i = 0; i<numEdge; i++) edgeSet[i] = Vector2D<int>(inputEdge[2*i],inputEdge[2*i+1]);
    //Initialize fluid info
    numXf = xFluidNum*yFluidNum;
    //Initialize intersector information
    segmentHierarchy = new Segment_Hierarchy(Xs,edgeSet);

    status = new int[numXf];
    status_n = new int[numXf];
    sweptNode = new bool[numXf];
    forwardMapping = new int[numXf];
    backwardMapping = new int[numXf];
    xEdgeIsIntersect = new bool[numXf];
    yEdgeIsIntersect = new bool[numXf];
    occlude = new bool[numXf];
    occlude_n = new bool[numXf];
    boxMax = new Vector2D<double>[numXf];
    boxMin = new Vector2D<double>[numXf];


    bool* intersectOrNot = new bool[numXf];
    Find_Node_Bounding_Boxes();
    Compute_Close_Segments(intersectOrNot, thickness);
    Find_Intersect_Result(intersectOrNot,thickness);
    Find_Status_Using_FloodFill();



}

Intersector::~Intersector() {
    delete[] status;
    delete[] status_n;
    delete[] occlude;
    delete[] occlude_n;
    delete[] sweptNode;
    delete[] forwardMapping;
    delete[] backwardMapping;
    delete[] boxMax;
    delete[] boxMin;
    delete[] xEdgeIsIntersect;
    delete[] yEdgeIsIntersect;


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

    if(segmentHierarchy) delete segmentHierarchy;

    for(int i = 0; i<numXs; i++) Xs[i] = Vector2D<double>(inputXs[2*i],inputXs[2*i+1]);
    segmentHierarchy = new Segment_Hierarchy(Xs,edgeSet);


    bool* intersectOrNot = new bool[numXf];
    Find_Node_Bounding_Boxes();
    int intersectNum = Compute_Close_Segments(intersectOrNot, thickness);

    for(int i = 0; i < numXf; i++) occlude_n[i] = occlude[i];
    Find_Intersect_Result(intersectOrNot,thickness);
    Compute_Swetp_Nodes(intersectNum, thickness);
    for(int i = 0; i < numXf; i++) status_n[i] = status[i];
    Find_Status();


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
            backwardMapping[shrunkIndex] = i;
            ++numCloseNodes;
            candidates.Append(cand);
            shrunkIndex++;
        }

    }
    return shrunkIndex;
}





int Intersector::Find_Intersect_Result(const bool * intersectOrNot,double thickness)
{
    for(int i = 0; i < numXf; i++){
        occlude[i] = false;
        xEdgeIsIntersect[i] = false;
        yEdgeIsIntersect[i] = false;
    }
    xEdgeResult.clear();
    yEdgeResult.clear();


    int intersectedEdgeCount=0;
    /* consider horizontal edges */
    for(int i = 0; i < yFluidNum; i++)
        for(int j = 0; j < xFluidNum - 1; j++){
            // node pair (i*xFluidNum + j, i*xFluidNum + j + 1)
            if(intersectOrNot[i*xFluidNum+j] && intersectOrNot[i*xFluidNum + j + 1]) {
                IntersectionResult edgeRes;
                if(Intersect_Result_Helper(i*xFluidNum + j, i*xFluidNum + j + 1, edgeRes, thickness)){
                    xEdgeResult[i*xFluidNum  + j] = edgeRes;
                    xEdgeIsIntersect[i*xFluidNum  + j] = true;
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
                    yEdgeIsIntersect[i*xFluidNum  + j] = true;
                    intersectedEdgeCount++;
                }
            }



        }


    return intersectedEdgeCount;
}




bool Intersector::Intersect_Result_Helper(int startNode, int endNode, IntersectionResult& edgesRes,double thickness){

    Vector2D<double> startNodeCoord = Xf[startNode], endNodeCoord = Xf[endNode];

    int candidateStartNode = forwardMapping[startNode], candidateEndNode = forwardMapping[endNode];
    Array<int> cand;
    for (int j = 0; j < candidates[candidateStartNode].Length(); ++j)
        for (int k = 0; k < candidates[candidateEndNode].Length(); ++k)
            if (candidates[candidateStartNode][j] == candidates[candidateEndNode][k]) {
                cand.Append(candidates[candidateStartNode][j]);
                break;
            }


    //check whether startNode and endNode are occluded
    for(int i = 0; i < cand.Length(); i++) {
        int segmentID = cand[i];
        double beta;
        Vector2D<double> startCandidateNodeCoord = Xs[edgeSet[segmentID][0]], endCandidateNodeCoord = Xs[edgeSet[segmentID][1]];
        if (Point_Inside_Segment(startNodeCoord, startCandidateNodeCoord, endCandidateNodeCoord, thickness,&beta)) {
            occlude[startNode] = true;
            edgesRes.structureSegmentCoord[0] = beta;
            edgesRes.structureSegmentID[0] = segmentID;
            edgesRes.edgeCoord[0] = 0.0;
        }
        if (Point_Inside_Segment(endNodeCoord, startCandidateNodeCoord, endCandidateNodeCoord, thickness, &beta)) {
            occlude[endNode] = true;
            edgesRes.structureSegmentCoord[1] = beta;
            edgesRes.structureSegmentID[1] = segmentID;
            edgesRes.edgeCoord[1] = 1.0;
        }
    }

    double alpha = 0.0, beta = 0.0;
    for(int i = 0; i < cand.Length(); i++) {
        int segmentID = cand[i];
        Vector2D<double> startCandidateNodeCoord = Xs[edgeSet[segmentID][0]], endCandidateNodeCoord = Xs[edgeSet[segmentID][1]];
        if(!occlude[startNode] && Ray_Interaction(&alpha, &beta, startNodeCoord,endNodeCoord,startCandidateNodeCoord,endCandidateNodeCoord,thickness))
        if(edgesRes.structureSegmentID[0] < 0 || alpha < edgesRes.edgeCoord[0]) {
            edgesRes.structureSegmentCoord[0] = beta;
            edgesRes.structureSegmentID[0] = segmentID;
            edgesRes.edgeCoord[0] = alpha;
            edgesRes.intersectTimes++;
        }
        if(!occlude[endNode] &&Ray_Interaction(&alpha, &beta, endNodeCoord,startNodeCoord,startCandidateNodeCoord,endCandidateNodeCoord,thickness))

        if(edgesRes.structureSegmentID[1] < 0 || 1 - alpha > edgesRes.edgeCoord[1]) {
            edgesRes.structureSegmentCoord[1] = beta;
            edgesRes.structureSegmentID[1] = segmentID;
            edgesRes.edgeCoord[1] = 1 - alpha;
            edgesRes.intersectTimes++;
        }



    }

    return occlude[endNode] || occlude[startNode] ||edgesRes.intersectTimes;
}

bool Intersector::Point_Inside_Segment(const Vector2D<double> &p, const Vector2D<double> &v1, const Vector2D<double> &v2, double thickness, double *beta) const{
    double dist = 0.0;
    double cross = (v2[0] - v1[0]) * (p[0] - v1[0]) + (v2[1] - v1[1]) * (p[1] - v1[1]);
    if (cross <= 0) { // PA1A2 is obtuse angle
        dist = sqrt((p[0] - v1[0]) * (p[0] - v1[0]) + (p[1] - v1[1]) * (p[1] - v1[1]));
        *beta = 0.0;
        return dist > thickness ? false:true;
    }

    double d2 = (v2[0] - v1[0]) * (v2[0] - v1[0]) + (v2[1] - v1[1]) * (v2[1] - v1[1]);
    if (cross >= d2) { // PA2A1 is obtuse angle
        dist =sqrt((p[0] - v2[0]) * (p[0] - v2[0]) + (p[1] - v2[1]) * (p[1] - v2[1]));
        *beta = 1.0;
        return dist > thickness ? false:true;
    }

    double r = cross / d2;
    *beta = r;
    double px = v1[0] + (v2[0] - v1[0]) * r;
    double py = v1[1] + (v2[1] - v1[1]) * r;
    dist = sqrt((p[0] - px) * (p[0] - px) + (p[1] - py) * (p[1] - py));
    return dist > thickness ? false:true;

}

bool Intersector::Ray_Interaction(double * alpha, double * beta, const Vector2D<double> &v1, const Vector2D<double> &v2,
                                  const Vector2D<double> &u1, const Vector2D<double> &u2,double thickness) const {

    double v21x = v2[0] - v1[0],v21y = v2[1] - v1[1], u21x = u2[0] - u1[0], u21y = u2[1] - u1[1];
    double vlength = sqrt(v21x*v21x + v21y*v21y),ulength = sqrt(u21x*u21x + u21y*u21y);
    double delta = v21x*u21y - u21x*v21y;
    //consider the case:parallel lines
    if(delta == 0.0) {
        std::cout <<"Parallel lines in Ray Interaction" << std::endl;
        double alpha1,alpha2;
        bool result1, result2;
        result1 = Point_Inside_Segment(u1, v1, v2, thickness,&alpha1);
        result2 = Point_Inside_Segment(u1, v1, v2, thickness,&alpha2);

        if(result1 && result2)
            alpha1 < alpha2? *alpha = alpha1, *beta = 0.0 :*alpha = alpha2, *beta = 1.0 ;
        else {
            if (result1){
                *alpha = alpha1;
                *beta = 0.0;
            }
            if (result2){
                *alpha = alpha2;
                *beta = 1.0;
            }
        }




        return result1 || result2;
    }


    *alpha = (u21y * (u1[0] - v1[0]) - u21x*(u1[1] - v1[1]))/delta;
    *beta = -(-v21y * (u1[0] - v1[0]) + v21x*(u1[1] - v1[1]))/delta;

    if(*alpha > -thickness/vlength && *alpha < 1.0+thickness/vlength && *beta > -thickness/ulength && *beta < 1.0 + thickness/ulength)
        return true;
    else
        return false;

}

void Intersector::Compute_Swetp_Nodes(const int intersectNum, double thickness){


    for(int i = 0; i < numXf; i++) sweptNode[i] = false;
    for(int i = 0 ; i < intersectNum; i++){
        int fluidNode = backwardMapping[i];
        const Vector2D<double>& fluidNodePosition = Xf[fluidNode];
        for(int j=0; j<=candidates[i].Length() && !sweptNode[fluidNode];++j){
            int segmentID = candidates[i][j];
            Vector2D<double> v1 = Xs[edgeSet[segmentID][0]], v2 = Xs[edgeSet[segmentID][1]];
            Vector2D<double> v1_n = Xs_n[edgeSet[segmentID][0]], v2_n = Xs_n[edgeSet[segmentID][1]];
            sweptNode[fluidNode] = Inside_Quad(fluidNodePosition, v1,v2,v1_n,v2_n, thickness);
        }
    }
}
/* consider a ray with start point at p, toward x = +oo
 * if the number of intersections between the ray and edges is odd p is inside the quad
 * otherwise outside
 */
bool Intersector::Inside_Quad(const Vector2D<double> &p, const Vector2D<double> &v1, const Vector2D<double> &v2,
                             const Vector2D<double> &v1_n, const Vector2D<double> &v2_n,double thickness){

    {
        int nCount = 4;
        Vector2D<double> Points[4] = { v1, v2, v2_n, v1_n };
        int nCross = 0;
        for (int i = 0; i < nCount; i++)
        {
            Vector2D<double> pStart = Points[i];
            Vector2D<double> pEnd   = Points[(i + 1) % nCount];
            double len = sqrt((pEnd[0] - pStart[0])*(pEnd[0] - pStart[0]) + (pEnd[1] - pStart[1])*(pEnd[1] - pStart[1]));

            double delta = pEnd[1] - pStart[1];
            //consider the case:parallel lines
            if(delta == 0.0) {
                std::cout <<"Parallel lines in InsideQuad" << std::endl;
                if (pStart[1] == p[1] && (p[0] - pStart[0])*(p[0] - pEnd[0]) < 0.0)
                    nCross++;
                continue;
            }

            double alpha = (p[1] - pStart[1])/delta;
            if(alpha  < -thickness/len || alpha > 1 + thickness/len) continue;
            double x = pStart[0] + alpha*(pEnd[0] - pStart[0]);
            if(x > p[0] - thickness) nCross++;

        }

        // odd number of intersections -> outside, even number of intersections -> inside
        return (nCross % 2 == 1);
    }
}

void Intersector::Find_Status(){

    // Easy stuff first
//#pragma omp parallel for
    for(int i = 0;numXf;++i)
        if(occlude[i]) status[i] = 0;
        else if(!sweptNode[i]&&!occlude_n[i]) status[i]=status_n[i];



    // Next handle swept nodes
    int iteration_count=0;
    int flags[2] = {1,1}; // flags = {needs_iteration, detected_change}
    while(flags[0]>0 && flags[1]>0){
        flags[0]=0;
        flags[1]=0;
        ++iteration_count;
//#pragma omp parallel for

        for(int i = 0; i < yFluidNum;++i)
            for(int j = 0; j < xFluidNum; j++){
                int myId = i*xFluidNum + j;
                if(status[myId] == -1){
                    int stat=-1;
                    int neighbor = (i-1)*xFluidNum + j;
                    if(i!=0 && !yEdgeIsIntersect[neighbor] && status[neighbor] > 0){
                        stat=status[neighbor];
                        flags[1]=1;
                        break;
                    }

                    neighbor = (i+1)*xFluidNum + j;
                    if(i!=yFluidNum-1 && !yEdgeIsIntersect[myId] && status[neighbor] > 0){
                        stat=status[neighbor];
                        flags[1]=1;
                        break;
                    }
                    neighbor = i*xFluidNum + j - 1;
                    if(j!=0 && !yEdgeIsIntersect[neighbor] && status[neighbor] > 0){
                        stat=status[neighbor];
                        flags[1]=1;
                        break;
                    }
                    neighbor = i*xFluidNum + j + 1;
                    if(j!=xFluidNum - 1 && !yEdgeIsIntersect[myId] && status[neighbor] > 0){
                        stat=status[neighbor];
                        flags[1]=1;
                        break;
                    }
                    if(stat == -1) flags[0]=1;
                    else status[i*xFluidNum + j]=stat;
                }
            }
    }
    // Finish any remaining untouched nodes
#pragma omp parallel for
    for(int i=0;i<numXf;++i)
        if(status[i]==-1)
            status[i]=0;



}


void Intersector::Find_Status_Using_FloodFill() {

    for(int i = 0; i < numXf; i++) status[i] = 0; //inactive nodes
    // Perform local floodFill
    std::queue<Vector2D<int>> floodFill;
    for(int i = 0 ; i < yFluidNum; i++) {
        floodFill.push(Vector2D<int>(i, 0));
        status[i*xFluidNum] = 1;
    }

    while(!floodFill.empty()){
        Vector2D<int> currentNode = floodFill.front();
        floodFill.pop();
        int i = currentNode[0], j = currentNode[1], myId = i*xFluidNum + j;
        int neighbor = (i-1)*xFluidNum + j;
        if(i!= 0 && status[neighbor] <= 0 && !occlude[neighbor] && !yEdgeIsIntersect[neighbor]){
            status[neighbor] = status[myId];
            floodFill.push(Vector2D<int>(i-1, j));

        }

        neighbor = (i+1)*xFluidNum + j;
        if(i!=yFluidNum-1 && status[neighbor] <= 0 && !occlude[neighbor] &&  !yEdgeIsIntersect[myId] ){
            status[neighbor] = status[myId];
            floodFill.push(Vector2D<int>(i+1, j));
        }
        neighbor = i*xFluidNum + j - 1;
        if(j!=0 && status[neighbor] <= 0 && !occlude[neighbor] &&  !xEdgeIsIntersect[neighbor]){
            status[neighbor] = status[myId];
            floodFill.push(Vector2D<int>(i, j-1));
        }
        neighbor = i*xFluidNum + j + 1;
        if(j!=xFluidNum - 1 && status[neighbor] <= 0 && !occlude[neighbor] &&  !xEdgeIsIntersect[myId]){
            status[neighbor] = status[myId];
            floodFill.push(Vector2D<int>(i, j+1));
        }



    }

}
void Intersector::Print_Intersector_Info(){
    for (int i = 0; i < numXf; i++)
        if(occlude[i])
            std::cout <<"node " << i << " is occluded" << std::endl;

    for (int i = 0; i < numXf; i++)
        std::cout <<"status of node " << i << " is " << status[i] << std::endl;


    std::cout << "HORIZONTAL EDGE" << std::endl;
    for(std::map<int,IntersectionResult>::iterator it = xEdgeResult.begin(); it!=xEdgeResult.end(); ++it){
        std::cout<<"1Edge: "<< it->first << " intersects: "<<it->second.structureSegmentID[0] <<" at its " <<it->second.edgeCoord[0]
        << " at " << it->second.structureSegmentCoord[0] <<std::endl;
        if(it->second.structureSegmentID[0] != -1)
            std::cout<<"2Edge: "<< it->first << " intersects: "<<it->second.structureSegmentID[1] <<" at its " <<it->second.edgeCoord[1]
            << " at " << it->second.structureSegmentCoord[1] <<std::endl;
    }

    std::cout << "VERTICAL EDGE" << std::endl;
    for(std::map<int,IntersectionResult>::iterator it = yEdgeResult.begin(); it!=yEdgeResult.end(); ++it){
        std::cout<<"1Edge: "<< it->first << " intersects: "<<it->second.structureSegmentID[0] <<" at its " <<it->second.edgeCoord[0]
        << " at " << it->second.structureSegmentCoord[0] <<std::endl;
        if(it->second.structureSegmentID[0] != -1)
            std::cout<<"2Edge: "<< it->first << " intersects: "<<it->second.structureSegmentID[1] <<" at its " <<it->second.edgeCoord[1]
            << " at " << it->second.structureSegmentCoord[1] <<std::endl;
    }

}

