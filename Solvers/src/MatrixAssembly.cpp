#include <iostream>
#include <algorithm>
#include "..\include\MatrixAssembly.hpp"

using namespace Solvers;

// typedef Eigen::SparseMatrix<double> SpD;
// typedef Eigen::MatrixXd DD;
// typedef Eigen::VectorXd DvD;

DD MatrixAssembly::GenerateQuadWeights(std::vector<double> &gpX, std::vector<double> &gpY, int numXNodes, int numYNodes, int numElemNodes) {
    // Generate quadrature weight matrices
    // Get integrals of basis functions
    std::vector<double> integX = Utils::integrateLagrange(gpX);
    std::vector<double> integY = Utils::integrateLagrange(gpY);

    Eigen::Map<DvD> integXMat(integX.data(),numXNodes);
    Eigen::Map<DvD> integYMat(integY.data(),numYNodes);

    DD weightMat(numElemNodes, numElemNodes);
    weightMat << Eigen::kroneckerProduct((DD)integXMat.asDiagonal(),(DD)integYMat.asDiagonal());

    return weightMat;
}

SpD MatrixAssembly::MassMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double c) {
    int xdeg = inputMesh.xdeg; int ydeg = inputMesh.ydeg;
    int numXNodes = xdeg + 1; int numYNodes = ydeg + 1;
    int numElemNodes = numXNodes * numYNodes;
    int nNodes = inputMesh.nNodes(); int nElements = inputMesh.nElements();
    std::vector<double> gaussPointsX = Utils::genGaussPoints(xdeg);
    std::vector<double> gaussPointsY = Utils::genGaussPoints(ydeg);

    // get base elem matrix
    DD weightMat = MatrixAssembly::GenerateQuadWeights(gaussPointsX, gaussPointsY, numXNodes, numYNodes, numElemNodes);
    weightMat *= c;

    // Turn to vector because matrix is diagonal
    DvD weightVector = weightMat.diagonal();

    // Initialize i,j,v triplet list for sparse matrix
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nElements * numElemNodes * numElemNodes);

    // Integrate over all elements
    DvD localElemMat(numElemNodes);
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
    SpD mat(nNodes,nNodes);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

SpD MatrixAssembly::StiffnessMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double k) {
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
    Eigen::Map<DD> Ax(AxInitializer.data(), numXNodes, numXNodes); 
    Eigen::Map<DD> Ay(AyInitializer.data(), numYNodes, numYNodes);

    // Generate quadrature weight matrices
    DD weightMat = MatrixAssembly::GenerateQuadWeights(gaussPointsX, gaussPointsY, numXNodes, numYNodes, numElemNodes);

    // Generate mass matrices
    DD Bx; Bx.setIdentity(numXNodes, numXNodes);
    DD By; By.setIdentity(numYNodes, numYNodes);

    DD coeffMat(numElemNodes, numElemNodes);
    coeffMat.setIdentity(); coeffMat *= k;
    
    // Get element-wise matrix intermediates
    DD combinedX(numElemNodes, numElemNodes);
    combinedX << Eigen::kroneckerProduct(By, Ax);
    DD combinedY(numElemNodes, numElemNodes);
    combinedY << Eigen::kroneckerProduct(Ay, Bx);

    // Initialize i,j,v triplet list for sparse matrix
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nElements * numElemNodes * numElemNodes);

    // Integrate over all elements
    DD localElemMat(numElemNodes, numElemNodes);
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
    SpD mat(nNodes,nNodes);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

SpD MatrixAssembly::AssembleFVec(Meshing::BasicMesh::BasicMesh2D &inputMesh, double f) {
    int xdeg = inputMesh.xdeg; int ydeg = inputMesh.ydeg;
    int numXNodes = xdeg + 1; int numYNodes = ydeg + 1;
    int numElemNodes = numXNodes * numYNodes;
    int nNodes = inputMesh.nNodes(); int nElements = inputMesh.nElements();
    std::vector<double> gaussPointsX = Utils::genGaussPoints(xdeg);
    std::vector<double> gaussPointsY = Utils::genGaussPoints(ydeg);

    // Generate quadrature weight matrices
    DD weightMat = MatrixAssembly::GenerateQuadWeights(gaussPointsX, gaussPointsY, numXNodes, numYNodes, numElemNodes);

    // Turn weight mat int vector and mult. by source since diagonal
    DvD sourceVec = f * weightMat.diagonal();

    // Initialize i,j,v triplet list for sparse matrix
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nElements * numElemNodes * numElemNodes);

    // Integrate over all elements
    DvD localElemMat(numElemNodes);
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
    SpD mat(nNodes,1);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

SpD MatrixAssembly::GetNullSpace(Meshing::BasicMesh::BasicMesh2D &inputMesh, std::vector<int> &boundaryNodes, bool nullSpace) {
    int nNodes = inputMesh.nNodes();
    int nNonZeroes;

    if (nullSpace) {
        nNonZeroes = boundaryNodes.size();
        }
    else {
        nNonZeroes = nNodes - boundaryNodes.size();
        }
    
    SpD nullSpaceMat(nNodes, nNonZeroes);

    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nNonZeroes);

    for (int i=0; i<nNonZeroes; i++) {
        tripletList.emplace_back(inputMesh.Nodes[i], i, 1.0);
    }

    return nullSpaceMat;
}

DvD MatrixAssembly::EvalBoundaryCond(Meshing::BasicMesh::BasicMesh2D &inputMesh, std::vector<int> &boundaryNodes, DvD (*DirichletBcs[4]) (std::vector<std::array<double, 2>>&)) {
    int xdeg = inputMesh.xdeg; int ydeg = inputMesh.ydeg;
    int numXElems = inputMesh.xOffsets.size(); int numYElems = inputMesh.yOffsets.size();
    int xWidth = xdeg*numXElems; int yWidth = ydeg*numYElems;
    
    int numBoundaryNodes = boundaryNodes.size();

    std::vector<std::array<double,2>> boundaryNodePos = inputMesh.posOfNodes(boundaryNodes);

    std::vector<double> boundaryNodeValues;
    boundaryNodeValues.reserve(numBoundaryNodes);

    for (int i=0; i<4; i++) {
        auto currFunc = *(DirichletBcs[i]);
    }

    Eigen::Map<DvD> boundaryNodeValuesVec(boundaryNodeValues.data(), numBoundaryNodes, 1);
    return (DvD)boundaryNodeValuesVec;
}

DvD MatrixAssembly::ApplyBoundaryCond(SpD &StiffnessMatrix, SpD &fVec, SpD &columnSpace, SpD &nullSpace, DvD &boundaryVals) {
    // TODO: Make sure boundary nodes are properly rearranged again
    SpD A11 = columnSpace.transpose() * StiffnessMatrix * columnSpace;
    SpD A12 = columnSpace.transpose() * StiffnessMatrix * nullSpace;
    SpD F11 = columnSpace.transpose() * fVec;
    SpD F21 = nullSpace.transpose() * fVec;

    // TODO: find proper sparse solving workflow
    DvD x;
    x = columnSpace * x;
    x += boundaryVals;
    return x;
}