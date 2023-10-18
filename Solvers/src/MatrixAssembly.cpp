#include <iostream>
#include <algorithm>
#include "..\include\MatrixAssembly.hpp"

using namespace Solvers;

Eigen::SparseMatrix<double> MatrixAssembly::MassMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double c) {
    int xdeg = inputMesh.xdeg; int ydeg = inputMesh.ydeg;
    int numXNodes = xdeg + 1; int numYNodes = ydeg + 1;
    int numElemNodes = numXNodes * numYNodes;
    int nNodes = inputMesh.nNodes(); int nElements = inputMesh.nElements();
    std::vector<double> gaussPointsX = Utils::genGaussPoints(xdeg);
    std::vector<double> gaussPointsY = Utils::genGaussPoints(ydeg);

    // Generate quadrature weight matrices
    // Get integrals of basis functions
    std::vector<double> integX = Utils::integrateLagrange(gaussPointsX);
    std::vector<double> integY = Utils::integrateLagrange(gaussPointsY);

    // Copy data into matrix
    Eigen::Map<Eigen::VectorXd> integXMat(integX.data(),numXNodes);
    Eigen::Map<Eigen::VectorXd> integYMat(integY.data(),numYNodes); 

    // get elem. matrix
    Eigen::MatrixXd weightMat(numElemNodes, numElemNodes);
    weightMat << Eigen::kroneckerProduct((Eigen::MatrixXd)integXMat.asDiagonal(),(Eigen::MatrixXd)integYMat.asDiagonal());
    weightMat *= c;

    // Turn to vector because matrix is diagonal
    Eigen::VectorXd weightVector = weightMat.diagonal();

    // Initialize i,j,v triplet list for sparse matrix
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nElements * numElemNodes * numElemNodes);

    // Integrate over all elements
    Eigen::VectorXd localElemMat(numElemNodes);
    for (auto &elm : inputMesh.Elements) {
        double Lx = elm.getWidth(); double Ly = elm.getHeight(); // Jacobian factors
        // calculate local matrix
        localElemMat = weightVector*Lx*Ly/4;
        // Get nodes in element
        std::vector<int> nodesInElm = elm.Nodes;
        // Generate i,j,v triplets
        for (int i=0; i<numElemNodes; i++) {
            tripletList.emplace_back(nodesInElm[i], nodesInElm[i], localElemMat(i));
        }     
    }

    // Declare and construct sparse matrix from triplets
    Eigen::SparseMatrix<double> mat(nNodes,nNodes);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

Eigen::SparseMatrix<double> MatrixAssembly::StiffnessMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double k) {
    int xdeg = inputMesh.xdeg; int ydeg = inputMesh.ydeg;
    int numXNodes = xdeg + 1; int numYNodes = ydeg + 1;
    int numElemNodes = numXNodes * numYNodes;
    int nNodes = inputMesh.nNodes(); int nElements = inputMesh.nElements();
    std::vector<double> gaussPointsX = Utils::genGaussPoints(xdeg);
    std::vector<double> gaussPointsY = Utils::genGaussPoints(ydeg);


    // Generate derivative matrices
    std::vector<double> AxInitializer; std::vector<double> AyInitializer;
    AxInitializer.reserve(numXNodes*numXNodes); AyInitializer.reserve(numYNodes*numYNodes);

    double *AxInitIdx = AxInitializer.data(); 
    double *AyInitIdx = AyInitializer.data(); 

    // Generate derivatives for each basis function, copy to full array
    for (int k=0; k<numXNodes; k++) { 
        std::vector<double> xPartials = Utils::numDeriv(.00001, k, gaussPointsX, gaussPointsX);
        std::copy(xPartials.begin(), xPartials.end(), AxInitIdx);
        AxInitIdx += numXNodes;
    }
    
    for (int k=0; k<numYNodes; k++) {
        std::vector<double> yPartials = Utils::numDeriv(.00001, k, gaussPointsY, gaussPointsY);
        std::copy(yPartials.begin(), yPartials.end(), AyInitIdx);
        AyInitIdx += numYNodes;
    }

    // map derivative values to matrix
    Eigen::Map<Eigen::MatrixXd> Ax(AxInitializer.data(), numXNodes, numXNodes); 
    Eigen::Map<Eigen::MatrixXd> Ay(AyInitializer.data(), numYNodes, numYNodes);

    std::cout<<Ax<<std::endl;

    // Generate quadrature weight matrices
    // Get integrals of basis functions
    std::vector<double> integX = Utils::integrateLagrange(gaussPointsX);
    std::vector<double> integY = Utils::integrateLagrange(gaussPointsY);

    Eigen::Map<Eigen::VectorXd> integXMat(integX.data(),numXNodes);
    Eigen::Map<Eigen::VectorXd> integYMat(integY.data(),numYNodes);

    Eigen::MatrixXd weightMat(numElemNodes, numElemNodes);
    weightMat << Eigen::kroneckerProduct((Eigen::MatrixXd)(integXMat.asDiagonal()),(Eigen::MatrixXd)(integYMat.asDiagonal()));

    std::cout<<weightMat<<std::endl;
    // Generate mass matrices
    Eigen::MatrixXd Bx; Bx.setIdentity(numXNodes, numXNodes);
    Eigen::MatrixXd By; By.setIdentity(numYNodes, numYNodes);

    Eigen::MatrixXd coeffMat(numElemNodes, numElemNodes);
    coeffMat.setIdentity(); coeffMat *= k;
    
    // Get element-wise matrix intermediates
    Eigen::MatrixXd combinedX(numElemNodes, numElemNodes);
    combinedX << Eigen::kroneckerProduct(By, Ax);
    Eigen::MatrixXd combinedY(numElemNodes, numElemNodes);
    combinedY << Eigen::kroneckerProduct(Ay, Bx);

    // Initialize i,j,v triplet list for sparse matrix
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nElements * numElemNodes * numElemNodes);

    // Integrate over all elements
    Eigen::MatrixXd localElemMat(numElemNodes, numElemNodes);
    for (auto &elm : inputMesh.Elements) {
        double Lx = elm.getWidth(); double Ly = elm.getHeight(); // Jacobian factors
        // calculate local matrix
        localElemMat = combinedX*coeffMat*weightMat*combinedX.transpose()*Ly/Lx +
                        combinedY*coeffMat*weightMat*combinedY.transpose()*Lx/Ly;
        
        // Get nodes in element
        std::vector<int> nodesInElm = elm.Nodes;
        // Generate i,j,v triplets
        for (int j=0; j<numElemNodes; j++) {
            for (int i=0; i<numElemNodes; i++) {
                tripletList.emplace_back(nodesInElm[i],nodesInElm[j],localElemMat(i,j));
            }
        }
    }

    // Declare and construct sparse matrix from triplets
    Eigen::SparseMatrix<double> mat(nNodes,nNodes);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

Eigen::SparseMatrix<double> MatrixAssembly::AssembleFVec(Meshing::BasicMesh::BasicMesh2D &inputMesh, double f) {
        int xdeg = inputMesh.xdeg; int ydeg = inputMesh.ydeg;
    int numXNodes = xdeg + 1; int numYNodes = ydeg + 1;
    int numElemNodes = numXNodes * numYNodes;
    int nNodes = inputMesh.nNodes(); int nElements = inputMesh.nElements();
    std::vector<double> gaussPointsX = Utils::genGaussPoints(xdeg);
    std::vector<double> gaussPointsY = Utils::genGaussPoints(ydeg);

    // Generate quadrature weight matrices
    // Get integrals of basis functions
    std::vector<double> integX = Utils::integrateLagrange(gaussPointsX);
    std::vector<double> integY = Utils::integrateLagrange(gaussPointsY);

    Eigen::Map<Eigen::VectorXd> integXMat(integX.data(),numXNodes);
    Eigen::Map<Eigen::VectorXd> integYMat(integY.data(),numYNodes);

    Eigen::MatrixXd weightMat(numElemNodes, numElemNodes);
    weightMat << Eigen::kroneckerProduct((Eigen::MatrixXd)integXMat.asDiagonal(),(Eigen::MatrixXd)integYMat.asDiagonal());

    Eigen::VectorXd sourceVec = f * weightMat.diagonal();

    // Initialize i,j,v triplet list for sparse matrix
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nElements * numElemNodes * numElemNodes);

    // Integrate over all elements
    Eigen::VectorXd localElemMat(numElemNodes);
    for (auto &elm : inputMesh.Elements) {
        double Lx = elm.getWidth(); double Ly = elm.getHeight(); // Jacobian factors
        // calculate local matrix
        localElemMat = sourceVec*Lx*Ly/4;
        // Get nodes in element
        std::vector<int> nodesInElm = elm.Nodes;
        // Generate i,j,v triplets
        for (int i=0; i<numElemNodes; i++) {
            tripletList.emplace_back(nodesInElm[i], 0, localElemMat(i));
        }  
    }

    // Declare and construct sparse matrix from triplets
    Eigen::SparseMatrix<double> mat(nNodes,1);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}