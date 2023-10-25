#include "..\include\MatrixAssembly.hpp"
#include <iostream>

using namespace Solvers;

DD MatrixAssembly::PoissonSolve(Meshing::BasicMesh::BasicMesh2D &inputMesh,
                std::vector<bcFunc> DirichletBcs,
                double c,
                double k,
                double f) {
    
    std::vector<int> boundaryNodes = inputMesh.getBoundaryNodes();
    std::vector<int> freeNodes = inputMesh.getFreeNodes();
    int nNodes = inputMesh.nNodes();

    SpD KMatrix = MatrixAssembly::StiffnessMatrix(inputMesh, k);
    SpD FMatrix = MatrixAssembly::AssembleFVec(inputMesh, f);
    DvD boundaryVals = MatrixAssembly::EvalBoundaryCond(inputMesh, boundaryNodes, DirichletBcs);

    SpD nullSpace(nNodes, boundaryNodes.size());
    SpD columnSpace(nNodes, freeNodes.size());
    MatrixAssembly::GetExtensionMatrices(inputMesh, boundaryNodes, freeNodes, nullSpace, columnSpace);

    DD x = MatrixAssembly::ComputeSolutionStationary(KMatrix, FMatrix, columnSpace, nullSpace, boundaryVals);
    return x;
}