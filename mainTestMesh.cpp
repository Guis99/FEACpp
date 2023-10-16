#include "Meshing/include/BasicMesh2D.hpp"

#include <iostream>
#include <vector>

int main() {
    std::vector<double> xdivs;
    int nxElem = 1;
    xdivs.reserve(nxElem);
    for (int i = 0; i < nxElem; i++) {
        xdivs.push_back(2);
    }

    std::vector<double> ydivs;
    int nyElem = 1;
    ydivs.reserve(nyElem);

    for (int i = 0; i < nyElem; i++) {
        ydivs.push_back(2);
    }

    int xdeg = 10; int ydeg = 10;

    Meshing::BasicMesh::BasicMesh2D mesh(xdeg, ydeg, xdivs, ydivs, -1, -1);

    std::cout<<mesh.nNodes()<<std::endl<<mesh.nElements()<<std::endl;

    int nodeIdx; std::cin>>nodeIdx;
    Meshing::BasicMesh::Element elem = mesh.Elements[nodeIdx];

    
    std::vector<std::array<double, 2>> nodePos = mesh.posInElem(nodeIdx);

    for (int y = 0; y < ydeg+1; y++) {
        for (int x = 0; x < ydeg+1; x++) {
            int ind = (xdeg+1)*y + x;
            std::cout<<nodePos[ind][0]<<", "<<nodePos[ind][1]<<" | ";
            // std::cout<<elem.Nodes[ind]<<" | ";
        }
        std::cout<<std::endl;
    }
}