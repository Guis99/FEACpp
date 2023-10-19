#include "..\..\..\3PDep\eigen-3.4.0\Eigen\Core"
#include "..\..\..\3PDep\eigen-3.4.0\Eigen\Sparse"
#include "..\..\..\3PDep\eigen-3.4.0\unsupported\Eigen\KroneckerProduct"
#include "..\..\Meshing\Meshing.hpp"

#ifndef SOLVERS
#define SOLVERS

typedef Eigen::SparseMatrix<double> SpD;
typedef Eigen::MatrixXd DD;
typedef Eigen::VectorXd DvD;

namespace Solvers {
    namespace MatrixAssembly {
        DD GenerateQuadWeights(std::vector<double> &gpX, std::vector<double> &gpY, int numXNodes, int numYNodes, int numElemNodes);
        SpD MassMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double c);
        SpD StiffnessMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double k);
        SpD AssembleFVec(Meshing::BasicMesh::BasicMesh2D &inputMesh, double f);
        SpD GetNullSpace(Meshing::BasicMesh::BasicMesh2D &inputMesh, std::vector<int> &boundaryNodes, bool nullSpace);
        DvD EvalBoundaryCond(Meshing::BasicMesh::BasicMesh2D &inputMesh, std::vector<int> &boundaryNodes, DvD (*DirichletBcs[4]) (std::vector<std::array<double, 2>>&));
        DvD ApplyBoundaryCond(SpD &StiffnessMatrix, SpD &fVec, SpD &columnSpace, SpD &nullSpace, DvD &boundaryVals);
    }
}

#endif