#include <iostream>
#include <algorithm>
#include "..\include\MatrixAssembly.hpp"
#include "..\..\Meshing\Meshing.hpp"

using namespace Solvers;

Eigen::SparseMatrix<double> MatrixAssembly::MassMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double c) {
    
}

Eigen::SparseMatrix<double> MatrixAssembly::StiffnessMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double k) {
    int xdeg = inputMesh.xdeg; int ydeg = inputMesh.ydeg;
    int nNodes = inputMesh.nNodes(); int nElements = inputMesh.nElements();
    std::vector<double> gaussPointsX = Utils::genGaussPoints(xdeg);
    std::vector<double> gaussPointsY = Utils::genGaussPoints(ydeg);
 
    Eigen::Matrix(double, xdeg, xdeg) Bx; Eigen::Matrix(double, ydeg, ydeg) By;

    double AxInitializer[(xdeg+1)*(xdeg+1)]; double AyInitializer[(ydeg+1)*(ydeg+1)];
    int AxInitIdx = 0; int AyInitIdx = 0;

    for (int k=0; k<xdeg+1; k++) {
        std::vector<double> xPartials = Utils::numDeriv(.00001, k, gaussPointsX, gaussPointsX);
        for (int kk=0; kk<xdeg+1;kk++) {
            AxInitializer[AxInitIdx++] = xPartials[kk];
        }
    }
    
    for (int k=0; k<ydeg+1; k++) {
        std::vector<double> yPartials = Utils::numDeriv(.00001, k, gaussPointsY, gaussPointsY);
        for (int kk=0; kk<ydeg+1;kk++) {
            AyInitializer[AyInitIdx++] = yPartials[kk];
        }
    }

    Eigen::Matrix(double, xdeg+1, xdeg+1) Ax(AxInitializer); 
    Eigen::Matrix(double, ydeg+1, ydeg+1) Ay(AyInitializer);

    Ax.transpose();
    Ay.transpose();

    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nElements*(xdeg+1)*(ydeg+1)*(xdeg+1)*(ydeg+1));

    for (int i=0; i<nElements; i++) {

    }

    Eigen::SparseMatrix<double> mat(nNodes,nNodes);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
    
}

Eigen::SparseMatrix<double> MatrixAssembly::ConvectionMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double b) {
    
}

Eigen::SparseMatrix<double> MatrixAssembly::AssembleFVec(Meshing::BasicMesh::BasicMesh2D &inputMesh, double f) {
    int xdeg = inputMesh.xdeg; int ydeg = inputMesh.ydeg;
    int nNodes = inputMesh.nNodes(); int nElements = inputMesh.nElements();
    std::vector<double> gaussPointsX = Utils::genGaussPoints(xdeg);
    std::vector<double> gaussPointsY = Utils::genGaussPoints(ydeg);

    for (int k=0; k<xdeg+1; k++) {
        std::vector<double> xPartials = Utils::numDeriv(.00001, k, gaussPointsX, gaussPointsX);
    }
    
    for (int k=0; k<ydeg+1; k++) {
        std::vector<double> yPartials = Utils::numDeriv(.00001, k, gaussPointsY, gaussPointsY);
    }

    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nElements*(xdeg+1)*(ydeg+1));
}