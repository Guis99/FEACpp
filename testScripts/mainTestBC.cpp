#include "..\Solvers\Solvers.hpp"
#include "..\Meshing\Meshing.hpp"

#include <iostream>
#include <vector>

std::vector<double> bc1(std::array<double,2>* startpoint, std::array<double,2>* endpoint, int allocSize) {
    std::vector<double> out;
    out.reserve(allocSize);
    std::cout<<"--------"<<std::endl;
    for (int i=0; i<allocSize; i++) {
        auto pushvar = *(startpoint+i);
        double x = pushvar[0]; double y = pushvar[1];
        std::cout<<x<<", "<<y<<std::endl;
        out.push_back(-x*(x-4));
    }

    return out;
}

int main() {
    std::vector<double> xdivs;
    int nxElem = 2;
    xdivs.reserve(nxElem);
    for (int i = 0; i < nxElem; i++) {
        xdivs.push_back(2);
    }

    std::vector<double> ydivs;
    int nyElem = 2;
    ydivs.reserve(nyElem);

    for (int i = 0; i < nyElem; i++) {
        ydivs.push_back(2);
    }

    int xdeg = 2; int ydeg = 2;

    int widthX = nxElem*xdeg+1;
    int widthY = nyElem*ydeg+1;

    Meshing::BasicMesh::BasicMesh2D mesh(xdeg, ydeg, xdivs, ydivs, 0, 0);  
    std::cout<<mesh.nNodes()<<std::endl; 

    std::vector<int> boundaryNodes = mesh.getBoundaryNodes();
    std::vector<int> freeNodes = mesh.getFreeNodes();
    std::vector<std::array<double,2>> boundaryNodePos = mesh.posOfNodes(boundaryNodes);
    int numBoundaryNodes = boundaryNodes.size();

    for (int i=0; i<numBoundaryNodes; i++) {
        std::cout<<boundaryNodes[i]<<", "<<boundaryNodePos[i][0]<<", "<<boundaryNodePos[i][1]<<std::endl;
    }

    std::vector<bcFunc> DirichletBcs;
    DirichletBcs.resize(4);

    DirichletBcs[0] = bc1;
    DirichletBcs[1] = bc1;
    DirichletBcs[2] = bc1;
    DirichletBcs[3] = bc1;

    auto Bcs = Solvers::MatrixAssembly::EvalBoundaryCond(mesh, boundaryNodes, DirichletBcs);
    std::cout<<"lol"<<std::endl;
    std::cout<<Bcs<<std::endl;
}