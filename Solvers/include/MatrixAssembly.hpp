#include "..\..\..\3PDep\eigen-3.4.0\Eigen\Core"
#include "..\..\..\3PDep\eigen-3.4.0\Eigen\Sparse"
#include "..\..\..\3PDep\eigen-3.4.0\unsupported\Eigen\KroneckerProduct"
#include "..\..\Meshing\Meshing.hpp"


namespace Solvers {
    namespace MatrixAssembly {
        Eigen::MatrixXd GenerateQuadWeights(std::vector<double> &gpX, std::vector<double> &gpY, int numXNodes, int numYNodes, int numElemNodes);
        Eigen::SparseMatrix<double> MassMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double c);
        Eigen::SparseMatrix<double> StiffnessMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double k);
        Eigen::SparseMatrix<double> AssembleFVec(Meshing::BasicMesh::BasicMesh2D &inputMesh, double f);
    }
}