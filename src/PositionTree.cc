#include "interface/PositionTree.h"

PositionTree::PositionTree(uint64* idx, TTree* tree, int nPlanes)
{
    tree_ = tree ? tree : new TTree();

    n_planes = nPlanes;
    
    index=idx;
}

void PositionTree::Init()
{
    //---global branches
    tree_->Branch("index", index, "index/l");
    n_hitsX = new int[n_planes];
    n_hitsY = new int[n_planes];
    X = new float[n_planes];
    Y = new float[n_planes];
    //---position branches
    tree_->Branch("n_planes", &n_planes, "n_planes/I");
    tree_->Branch("n_hitsX", n_hitsX, "n_hitsX[n_planes]/I");
    tree_->Branch("n_hitsY", n_hitsY, "n_hitsY[n_planes]/I");
    tree_->Branch("X", X, "X[n_planes]/F");
    tree_->Branch("Y", Y, "Y[n_planes]/F");
}
