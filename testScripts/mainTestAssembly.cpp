#include "..\Solvers\Solvers.hpp"
#include "..\Meshing\Meshing.hpp"


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

    int xdeg = 20; int ydeg = 3;

    int widthX = nxElem*xdeg+1;
    int widthY = nyElem*ydeg+1;

    Meshing::BasicMesh::BasicMesh2D mesh(xdeg, ydeg, xdivs, ydivs, 0, 0);   

    Eigen::SparseMatrix K = Solvers::MatrixAssembly::StiffnessMatrix(mesh, 2);
    Eigen::SparseMatrix M = Solvers::MatrixAssembly::MassMatrix(mesh,2);
    Eigen::SparseMatrix F = Solvers::MatrixAssembly::AssembleFVec(mesh, 2);

    std::cout<<K<<std::endl<<M<<std::endl<<F<<std::endl;
}