#include "..\Solvers\Solvers.hpp"
#include "..\Meshing\Meshing.hpp"
#include "..\Dependencies\MathParser\MathParser.hpp"

#include <iostream>
#include <fstream>
#include <vector>

std::vector<double> bc1(std::array<double,2>* startpoint, int allocSize) {
    MathParser::QuickArray<double> x;
    MathParser::QuickArray<double> y;

    x.reserve(allocSize); y.reserve(allocSize);

    for (int i=0; i<allocSize; i++) {
        auto coord = *(startpoint+i);
        x.push_back(coord[0]);
        y.push_back(coord[1]);
    }

    MathParser::InitMaps();
    MathParser::SetVariable("x", x);
    MathParser::SetVariable("y", y);
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::string str;
    std::cout<<"here1"<<std::endl;
    std::cout<<"Specify Boundary Condition: ";
    std::cout<<"here2"<<std::endl;
    std::getline(std::cin,str);
    std::cout<<"here3"<<std::endl;
    std::cout<<std::endl;

    std::cout<<"here4"<<std::endl;
    std::cout<<str<<std::endl;
    auto result = MathParser::ParseText(str);
    std::cout<<"here5"<<std::endl;

    auto out = result.release();

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
    double f = 0;

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

    x.resize(widthX,widthY);
    std::cout<<x<<std::endl;
}