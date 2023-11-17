#include "..\Solvers\Solvers.hpp"
#include "..\Meshing\Meshing.hpp"

#include <iostream>
#include <fstream>
#include <vector>

using namespace Solvers;

int main() {
    int nxElem;
    int nyElem;
    int xdeg; int ydeg;

    std::cin>>nxElem;
    std::cin>>nyElem;
    std::cin>>xdeg;
    std::cin>>ydeg;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::vector<double> xdivs;
    xdivs.reserve(nxElem);
    for (int i = 0; i < nxElem; i++) {
        xdivs.push_back(4.0/nxElem);
    }

    std::vector<double> ydivs;
    ydivs.reserve(nyElem);

    for (int i = 0; i < nyElem; i++) {
        ydivs.push_back(4.0/nyElem);
    }

    int widthX = nxElem*xdeg+1;
    int widthY = nyElem*ydeg+1;

    Meshing::BasicMesh::BasicMesh2D mesh(xdeg, ydeg, xdivs, ydivs, 0, 0);   

    double c = 1;
    double k = 1;
    double f = 0;

    double timeStep;
    int numTimeSteps;

    auto allNodePos = mesh.allNodePos();
    auto startpoint = allNodePos.data();
    int allocSize = allNodePos.size();
    auto ICVec = Utils::EvalSymbolicBC(startpoint, allocSize, "Specify initial condition");
    Eigen::Map<DvD> InitialCondition(ICVec.data(), allocSize, 1);
    DvD initialCondition = (DvD)InitialCondition;

    std::cout<<"Specify time-step: "<<std::endl;
    std::cin>>timeStep;
    std::cout<<"Specify number of time-steps: "<<std::endl;
    std::cin>>numTimeSteps;

    std::vector<int> boundaryNodes = mesh.getBoundaryNodes();
    std::vector<int> freeNodes = mesh.getFreeNodes();
    int nNodes = mesh.nNodes();

    std::cout<<"Assembling LHS matrix"<<std::endl;
    SpD KMatrix = MatrixAssembly::StiffnessMatrix(mesh, k);
    std::cout<<"Assembling Mass matrix"<<std::endl;
    SpD MMatrix = MatrixAssembly::MassMatrix(mesh, c);
    std::cout<<"Assembling RHS vector"<<std::endl;
    SpD FMatrix = MatrixAssembly::AssembleFVec(mesh, f);
    std::cout<<"Assembling boundary condition vector"<<std::endl;
    DvD boundaryVals = MatrixAssembly::EvalDirichletBoundaryCond(mesh, boundaryNodes);

    SpD nullSpace(nNodes, boundaryNodes.size());
    SpD columnSpace(nNodes, freeNodes.size());
    std::cout<<"Assembling null-orth matrix"<<std::endl;
    MatrixAssembly::GetExtensionMatrices(mesh, boundaryNodes, freeNodes, nullSpace, columnSpace);
    std::cout<<"Solving system with "<<freeNodes.size()<<" nodes"<<std::endl;
    auto solns = MatrixAssembly::ComputeSolutionTimeDependent1stOrder(KMatrix, MMatrix, FMatrix, columnSpace, nullSpace, boundaryVals,
                                                                initialCondition, timeStep, numTimeSteps);

    std::vector<double> xGrid;
    std::vector<double> yGrid;
    xGrid.reserve(widthX);
    yGrid.reserve(widthY);

    for (int i=0; i<widthX; i++) {
        auto xcoords = mesh.posOfNodes(std::vector<int>{i});
        xGrid.push_back(xcoords[0][0]);
    } 

    for (int i=0; i<widthY; i++) {
        int yIdx = i*widthX;
        auto ycoords = mesh.posOfNodes(std::vector<int>{yIdx});
        yGrid.push_back(ycoords[0][1]);
    } 

    Eigen::Map<DD> xOffsets(xGrid.data(), widthX, 1);
    Eigen::Map<DD> yOffsets(yGrid.data(), widthY, 1);

    std::ofstream fileX("x.txt");
    std::ofstream fileY("y.txt");
    std::ofstream fileZ("z.txt");

    if (fileZ.is_open())
    {
        for (auto x : solns) {
            fileZ << x;
        }
    }
    if (fileX.is_open())
    {
        fileX << xOffsets;
    }
    if (fileY.is_open())
    {
        fileY << yOffsets;
    }
}