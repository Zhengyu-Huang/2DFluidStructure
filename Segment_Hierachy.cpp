//
// Created by icme-huang on 8/28/16.
//

#include "Segment_Hierachy.h"

void Segment_Hierachy::KDTree_Initialization() {
    int segmentNumber = edgeSet.Length();
    Array<Vector2D<double>> centroids(segmentNumber);
    for(int i = 0; i < segmentNumber; i++) centroids[i] = (points[edgeSet[i][0]] + points[edgeSet[i][1]]) /2.0;




}