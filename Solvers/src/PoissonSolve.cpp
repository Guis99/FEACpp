#include "..\include\MatrixAssembly.hpp"
#include <iostream>

using namespace Solvers;

DD MatrixAssembly::PoissonSolve(Meshing::BasicMesh::BasicMesh2D &inputMesh,
                double c,
                double k,
                double f) {
    
    std::vector<int> boundaryNodes = inputMesh.getBoundaryNodes();
    std::vector<int> freeNodes = inputMesh.getFreeNodes();
    int nNodes = inputMesh.nNodes();

    std::cout<<"Assembling LHS matrix"<<std::endl;
    SpD KMatrix = MatrixAssembly::StiffnessMatrix(inputMesh, k);
    std::cout<<"Assembling RHS vector"<<std::endl;
    SpD FMatrix = MatrixAssembly::AssembleFVec(inputMesh, f);
    std::cout<<"Assembling boundary condition vector"<<std::endl;
    DvD boundaryVals = MatrixAssembly::EvalDirichletBoundaryCond(inputMesh, boundaryNodes);

    SpD nullSpace(nNodes, boundaryNodes.size());
    SpD columnSpace(nNodes, freeNodes.size());
    std::cout<<"Assembling null-orth matrix"<<std::endl;
    MatrixAssembly::GetExtensionMatrices(inputMesh, boundaryNodes, freeNodes, nullSpace, columnSpace);
    std::cout<<"Solving system with "<<freeNodes.size()<<" nodes"<<std::endl;
    DD x = MatrixAssembly::ComputeSolutionStationary(KMatrix, FMatrix, columnSpace, nullSpace, boundaryVals);
    return x;
}