#include "..\Solvers\Solvers.hpp"
#include "..\Meshing\Meshing.hpp"

#include <iostream>
#include <vector>

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

    Eigen::SparseMatrix K = Solvers::MatrixAssembly::StiffnessMatrix(mesh, 2);
    Eigen::SparseMatrix M = Solvers::MatrixAssembly::MassMatrix(mesh, 2);
    Eigen::SparseMatrix F = Solvers::MatrixAssembly::AssembleFVec(mesh, 2);

    // std::cout<<K<<std::endl<<M<<std::endl<<F<<std::endl;
    std::vector<int> boundaryNodes = mesh.getBoundaryNodes();
    std::vector<int> freeNodes = mesh.getFreeNodes();

    Eigen::SparseMatrix ns = Solvers::MatrixAssembly::GetNullSpace(mesh, boundaryNodes, 1);
    Eigen::SparseMatrix cs = Solvers::MatrixAssembly::GetNullSpace(mesh, freeNodes, 0);

    std::cout<<ns<<std::endl<<cs<<std::endl;
}