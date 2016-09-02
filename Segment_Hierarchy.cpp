//
// Created by icme-huang on 8/28/16.
//

#include "Segment_Hierarchy.h"
Segment_Hierarchy::Segment_Hierarchy(const Array<Vector2D<double>> &pointsInput,const Array<Vector2D<int>> &edgeSetInput,
                                   const bool update_boxes):points(pointsInput),edgeSet(edgeSetInput){

    if(pointsInput.len) {
        Initialize_Hierarchy_Using_KD_Tree();
        if (update_boxes) Update_Boxes(); else {
            numLeaves = 0;
            root = 0;
        }
    }

}


void Segment_Hierarchy::Initialize_Hierarchy_Using_KD_Tree() {
    int segmentNumber = edgeSet.Length();
    Array<Vector2D<double>> centroids(segmentNumber);
    for(int i = 0; i < segmentNumber; i++) centroids[i] = (points[edgeSet[i][0]] + points[edgeSet[i][1]]) /2.0;
    KD_Tree kdTree;

    kdTree.Create_KD_Tree(centroids);
    numLeaves=edgeSet.len;
    parents.Resize(numLeaves);
    children.Remove_All();
    root=this->Initialize_Hierarchy_Using_KD_Tree_Helper(kdTree.rootNode);
    assert(root==2*numLeaves-2);
    boxHierarchy.Resize(root+1);
    kdTree.Print_Tree_Info();
    //boxRadius.Resize(root);

}

int Segment_Hierarchy::Initialize_Hierarchy_Using_KD_Tree_Helper(KD_Tree_Node * node) {
    if(!node->left&&!node->right)return node->nodeIndex;
    int leftChild=Initialize_Hierarchy_Using_KD_Tree_Helper(node->left);
    if(!node->right)return leftChild;
    int rightChild=Initialize_Hierarchy_Using_KD_Tree_Helper(node->right);
    children.Append(Vector2D<int>(leftChild,rightChild));
    parents.Append(0);
    return parents[leftChild]=parents[rightChild]=children.len+numLeaves-1;

}
void Segment_Hierarchy::Update_Boxes(double enlargeThickness){
    Update_Leaf_Boxes( );
    Update_Nonleaf_Boxes();
}

void Segment_Hierarchy::Update_Leaf_Boxes(double enlargeThickness){
    for(int i=0;i<numLeaves;i++)

        boxHierarchy[i]=Bounding_Box_2D(Vector2D<double>::Componentwise_Min(points[edgeSet[i][0]],points[edgeSet[i][1]]),
                                        Vector2D<double>::Componentwise_Max(points[edgeSet[i][0]],points[edgeSet[i][1]]));

    if(enlargeThickness)
        for(int i=0;i<numLeaves;i++)
            boxHierarchy[i].Enlarge(enlargeThickness);
}

void Segment_Hierarchy::Update_Nonleaf_Boxes(){
    for(int i=numLeaves;i<boxHierarchy.len;i++)
        boxHierarchy[i]=Bounding_Box_2D::Combine(boxHierarchy[children[i-numLeaves][0]],boxHierarchy[children[i-numLeaves][1]]);
}

void Segment_Hierarchy::Print_Info(){
    std::cout<<"Parents"<<std::endl;
    for(int i = 0; i < parents.len; i++)
        std::cout << parents[i] << " ";
    std::cout<<std::endl;

    std::cout<<"Children"<<std::endl;
    for(int i = 0; i < children.len; i++)
        std::cout << children[i];
    std::cout<<std::endl;

    std::cout<<"Box Hierarchy"<<std::endl;
    for(int i = 0; i < boxHierarchy.len; i++)
        std::cout << boxHierarchy[i];
    std::cout<<std::endl;
}