#include "..\Solvers\Solvers.hpp"
#include "..\Meshing\Meshing.hpp"

#include <iostream>
#include <fstream>
#include <vector>

std::vector<double> bc1(std::array<double,2>* startpoint, int allocSize) {
    std::vector<double> out;
    out.reserve(allocSize);
    std::cout<<"--------"<<std::endl;
    for (int i=0; i<allocSize; i++) {
        auto pushvar = *(startpoint+i);
        double x = pushvar[0]; double y = pushvar[1];
        std::cout<<x<<", "<<y<<std::endl;
        out.push_back(-x*(x-4) + y*(y-4));
    }

    return out;
}

int main() {
    int nxElem;
    int nyElem;
    int xdeg; int ydeg;

    std::cin>>nxElem;
    std::cin>>nyElem;
    std::cin>>xdeg;
    std::cin>>ydeg;

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
    std::cout<<mesh.nNodes()<<std::endl;

    std::vector<bcFunc> DirichletBcs;
    DirichletBcs.resize(4);

    DirichletBcs[0] = bc1;
    DirichletBcs[1] = bc1;
    DirichletBcs[2] = bc1;
    DirichletBcs[3] = bc1;

    double c = 1;
    double k = 1;
    double f;

    std::cin>>f;

    DD x = Solvers::MatrixAssembly::PoissonSolve(mesh, DirichletBcs, c, k, f);

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
        fileZ << x;
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