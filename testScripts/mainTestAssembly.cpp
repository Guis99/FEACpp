#include "..\Solvers\Solvers.hpp"
#include "..\Meshing\Meshing.hpp"

#include <iostream>
#include <vector>

int main() {
    int nxElem;
    int nyElem;
    int xdeg; int ydeg;

    std::cin>>nxElem;
    std::cin>>nyElem;
    std::cin>>xdeg;
    std::cin>>ydeg;

    std::vector<double> xdivs;
    xdivs.reserve(nxElem);
    for (int i = 0; i < nxElem; i++) {
        xdivs.push_back(2);
    }

    std::vector<double> ydivs;
    ydivs.reserve(nyElem);

    for (int i = 0; i < nyElem; i++) {
        ydivs.push_back(2);
    }

    int widthX = nxElem*xdeg+1;
    int widthY = nyElem*ydeg+1;

    Meshing::BasicMesh::BasicMesh2D mesh(xdeg, ydeg, xdivs, ydivs, 0, 0);   
    std::cout<<mesh.nNodes()<<std::endl;

    Eigen::SparseMatrix K = Solvers::MatrixAssembly::StiffnessMatrix(mesh, 2);
    std::cout<<"one done"<<std::endl;
    Eigen::SparseMatrix M = Solvers::MatrixAssembly::MassMatrix(mesh, 2);
    std::cout<<"two done"<<std::endl;
    Eigen::SparseMatrix F = Solvers::MatrixAssembly::AssembleFVec(mesh, 2);
    std::cout<<"three done"<<std::endl;

    int nNodes = mesh.nNodes();
    Eigen::SparseMatrix<double> Kr(nNodes,nNodes);
    Eigen::SparseMatrix<double> Mr(nNodes,nNodes);
    Eigen::SparseMatrix<double> Fr(nNodes,1);

    // Solvers::MatrixAssembly::AssembleMatrices(mesh,Mr,Kr,Fr,2,2,2);
    std::cout<<nNodes<<"x"<<nNodes<<" matrix generated"<<std::endl;

    // std::cout<<K-Kr<<std::endl<<M-Mr<<std::endl<<F-Fr<<std::endl;


    // std::vector<int> boundaryNodes = mesh.getBoundaryNodes();
    // std::vector<int> freeNodes = mesh.getFreeNodes();

    // Eigen::SparseMatrix ns = Solvers::MatrixAssembly::GetNullSpace(mesh, boundaryNodes, 1);
    // Eigen::SparseMatrix cs = Solvers::MatrixAssembly::GetNullSpace(mesh, freeNodes, 0)

    auto boundaryNodes = mesh.getBoundaryNodes();
    auto freeNodes = mesh.getFreeNodes();

    SpD ns(nNodes,boundaryNodes.size());
    SpD cs(nNodes,freeNodes.size());

    auto boundaryNodesInterm = boundaryNodes;
    std::sort(boundaryNodesInterm.begin(), boundaryNodesInterm.end());

    std::vector<int> boundaryNodesSorted;
    boundaryNodesSorted.resize(boundaryNodes.size());

    std::vector<int> moveIndices;
    moveIndices.reserve(boundaryNodes.size());

    for (int i=0; i<boundaryNodes.size(); i++) {
        auto elmToMove = boundaryNodes[i];
        int lb = 0; int ub = boundaryNodes.size()-1;
        int mid = (lb+ub)/2;
        while (elmToMove != boundaryNodesInterm[mid]) {
            mid = (lb+ub)/2;
            if (boundaryNodesInterm[mid] > elmToMove) {
                ub = mid;
            }
            else {
                lb = mid+1;
            }
        }
        boundaryNodesSorted[mid] = elmToMove;
        moveIndices.push_back(mid);
    }

    for (int i=0; i<boundaryNodes.size(); i++) {
        std::cout<<boundaryNodes[i]<<"-> "<<boundaryNodesSorted[i]<<std::endl;
    }

    // Solvers::MatrixAssembly::GetExtensionMatrices(mesh,boundaryNodesSorted,freeNodes,ns,cs);

    // std::cout<<ns<<std::endl<<cs<<std::endl;
}