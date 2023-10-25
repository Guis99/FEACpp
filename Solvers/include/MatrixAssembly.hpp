#include "..\..\Dependencies\Eigen\Core"
#include "..\..\Dependencies\Eigen\Sparse"
#include "..\..\Dependencies\unsupported\Eigen\KroneckerProduct"
#include "..\..\Meshing\Meshing.hpp"

#ifndef SOLVERS
#define SOLVERS

typedef Eigen::SparseMatrix<double> SpD;
// Dynamically-sized matrix of doubles
typedef Eigen::MatrixXd DD;
// Dynamically-sized vector of doubles
typedef Eigen::VectorXd DvD;
// pointer to function returning std::vector<float>
typedef std::vector<double> (*bcFunc) (std::array<double,2>*,int); 

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
        // SpD GetNullSpace(Meshing::BasicMesh::BasicMesh2D &inputMesh, std::vector<int> &nodesToAdd, bool nullSpace);
        void GetExtensionMatrices(Meshing::BasicMesh::BasicMesh2D &inputMesh,
                                        std::vector<int> &boundaryNodes, 
                                        std::vector<int> &freeNodes,
                                        SpD &nullSpace,
                                        SpD &columnSpace);
        DvD EvalBoundaryCond(Meshing::BasicMesh::BasicMesh2D &inputMesh, std::vector<int> &boundaryNodes, std::vector<bcFunc>);
        DvD ComputeSolutionStationary(SpD &StiffnessMatrix, SpD &fVec, SpD &columnSpace, SpD &nullSpace, DvD &boundaryVals);
        DD PoissonSolve(Meshing::BasicMesh::BasicMesh2D &inputMesh,
                std::vector<bcFunc> DirichletBcs,
                double c,
                double k,
                double f);
    }
}

#endif