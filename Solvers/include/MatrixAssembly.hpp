#include "..\..\..\3PDep\eigen-3.4.0\Eigen\Core"
#include "..\..\..\3PDep\eigen-3.4.0\Eigen\Sparse"
#include "..\..\..\3PDep\eigen-3.4.0\unsupported\Eigen\KroneckerProduct"

#include "..\Utils\Utils.hpp"
#include "..\Meshing\Meshing.hpp"

namespace Solvers {
    namespace MatrixAssembly {
        Eigen::SparseMatrix<double> MassMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double c);
        Eigen::SparseMatrix<double> StiffnessMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double k);
        Eigen::SparseMatrix<double> ConvectionMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double b);
        Eigen::SparseMatrix<double> AssembleFVec(Meshing::BasicMesh::BasicMesh2D &inputMesh, double f);
    }
}