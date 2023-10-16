#include "..\include\BasicMesh2D.hpp"
#include "..\..\Utils\include\LagrangeInterpolant.hpp"

#include <iostream>

Meshing::BasicMesh::Element::Element(int EID, std::vector<int> Nodes, double Boundaries[4]) {
    this->EID = EID;
    this->Nodes = Nodes;
    this->Boundaries;
}

Meshing::BasicMesh::Node::Node(int NID, std::array<double, 2> Position, int nClass) {
    this->NID = NID;
    this->Position = Position;
    this->nClass = nClass;
}

Meshing::BasicMesh::BasicMesh2D::BasicMesh2D(int xdeg, int ydeg, 
                    std::vector<double> xdiv, std::vector<double> ydiv, 
                    double xstart, double ystart) {
    this->xdeg = xdeg;
    this->ydeg = ydeg;
    this->ydiv = ydiv;
    this->xdiv = xdiv;

    int nElemX = xdiv.size();
    int nElemY = ydiv.size();

    xOffsets.reserve(nElemX+1);
    yOffsets.reserve(nElemY+1);

    xOffsets[0] = xstart;
    yOffsets[0] = ystart;

    for (int i = 1; i < nElemX+1; i++) {
        xOffsets[i] = xOffsets[i-1]+xdiv[i-1];
    }

    for (int i = 1; i < nElemY+1; i++) {
        yOffsets[i] = yOffsets[i-1]+ydiv[i-1];
    }

    std::vector<double> xSpacing = Utils::genGaussPoints(xdeg);
    std::vector<double> ySpacing = Utils::genGaussPoints(ydeg);

    int widthX = nElemX*xdeg+1;
    int widthY = nElemY*ydeg+1;

    Elements.reserve(nElemX*nElemY);
    Nodes.reserve(widthX*widthY);

    int xIndLocal;
    int xElemInd;
    int yIndLocal = -2; 
    int yElemInd = 0;

    int NID = -1; int nClass;
    for (int y = 0; y < widthY; y++) {
        yIndLocal++;
        xIndLocal = -2;
        xElemInd = 0;
        double ay = yOffsets[yElemInd]; double by = yOffsets[yElemInd+1];
        for (int x = 0; x < widthX; x++) {
            xIndLocal++; NID++;
            double ax = xOffsets[xElemInd]; double bx = xOffsets[xElemInd+1];

            if (y==0 || y==widthY-1 || x==0 || x==widthX-1) {
                nClass = 1;
            }
            else {
                nClass = 0;
            }

            std::array<double,2> pos;
            pos[0] = transformPoint(xSpacing[xIndLocal+1],ax,bx);
            pos[1] = transformPoint(ySpacing[yIndLocal+1],ay,by);
            Nodes.emplace_back(NID, pos, nClass);
            if (xIndLocal == xdeg-1) {
                xElemInd++;
                xIndLocal = -1;
            }
        }
        if (yIndLocal == ydeg-1) {
            yElemInd++;
            yIndLocal = -1;
        }
    }
    int EID = -1;
    for (int j = 0; j < nElemY; j++) {
        for (int i = 0; i < nElemX; i++) {
            EID++;
            double bounds[4] = {xOffsets[i],xOffsets[i+1],yOffsets[j],yOffsets[j+1]};
            std::vector<int> currNodes;
            currNodes.reserve((xdeg+1)*(ydeg+1));

            for (int yy = (i)*xdeg+1; yy < (i+1)*xdeg+2; yy++) {
                for (int xx = (j)*ydeg+1; xx < (j+1)*ydeg+2; xx++) {
                    currNodes.push_back((yy-1)*widthX+xx-1);
                }
            }

            Elements.emplace_back(EID, currNodes, bounds);
        }
    }
}


std::vector<std::array<double, 2>> Meshing::BasicMesh::BasicMesh2D::allNodePos() {
    std::vector<std::array<double, 2>> out;
    out.reserve(Nodes.size());

    for (const auto &node : Nodes) {
        out.push_back(node.Position);
    }

    return out;
}

std::vector<std::array<double, 2>> Meshing::BasicMesh::BasicMesh2D::posInElem(int ElemID) {
    std::vector<std::array<double, 2>> out;
    out.reserve((xdeg+1)*(ydeg+1));

    for (const auto &node : Elements[ElemID].Nodes) {
        std::cout<<Nodes[node].Position[0]<<", "<<Nodes[node].Position[1]<<std::endl;
        out.push_back(Nodes[node].Position);
    }

    return out;
}

int Meshing::BasicMesh::BasicMesh2D::nNodes() {
    return Nodes.size();
}

int Meshing::BasicMesh::BasicMesh2D::nElements() {
    return Elements.size();
}

double Meshing::BasicMesh::BasicMesh2D::transformPoint(double x, double a, double b) {
    double result = a + ((b - a) / 2.0) * (x + 1.0);
    return result;
}


