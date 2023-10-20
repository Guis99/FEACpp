#include "..\..\..\3PDep\eigen-3.4.0\Eigen\Core"
#include "..\..\..\3PDep\eigen-3.4.0\Eigen\Sparse"
#include "..\..\..\3PDep\eigen-3.4.0\unsupported\Eigen\KroneckerProduct"
#include "..\..\Meshing\Meshing.hpp"

#ifndef SOLVERS
#define SOLVERS

typedef Eigen::SparseMatrix<double> SpD;
// Dynamically-sized matrix of doubles
typedef Eigen::MatrixXd DD;
// Dynamically-sized vector of doubles
typedef Eigen::VectorXd DvD;
// pointer to function returning std::vector<float>
typedef std::vector<double> (*bcFunc) (std::array<double,2>*,std::array<double,2>*,int); 

namespace Solvers {
    namespace MatrixAssembly {
        DD GenerateQuadWeights(std::vector<double> &gpX, std::vector<double> &gpY, int numXNodes, int numYNodes, int numElemNodes);
        void AssembleMatrices(Meshing::BasicMesh::BasicMesh2D &inputMesh,
                    SpD &MassMatrix,
                    SpD &StiffnessMatrix,
                    SpD &SourceVector,
                    double c,
                    double k,
                    double f);
        SpD MassMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double c);
        SpD StiffnessMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double k);
        SpD AssembleFVec(Meshing::BasicMesh::BasicMesh2D &inputMesh, double f);
        SpD GetNullSpace(Meshing::BasicMesh::BasicMesh2D &inputMesh, std::vector<int> &nodesToAdd, bool nullSpace);
        DvD EvalBoundaryCond(Meshing::BasicMesh::BasicMesh2D &inputMesh, std::vector<int> &boundaryNodes, std::vector<bcFunc>);
        DvD ComputeSolution(SpD &StiffnessMatrix, SpD &fVec, SpD &columnSpace, SpD &nullSpace, DvD &boundaryVals);
    }
}

#endif