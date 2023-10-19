#include "..\include\BasicMesh2D.hpp"

#include <iostream>

Meshing::BasicMesh::Element::Element(int EID, std::vector<int> Nodes, std::vector<double> Boundaries) {
    this->EID = EID;
    this->Nodes = Nodes;
    this->Boundaries = Boundaries;
}

double Meshing::BasicMesh::Element::getWidth() {
    return Boundaries[1] - Boundaries[0];
}

double Meshing::BasicMesh::Element::getHeight() {
    return Boundaries[3] - Boundaries[2];
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
            std::vector<int> currNodes;
            currNodes.reserve((xdeg+1)*(ydeg+1));

            for (int yy = j*ydeg; yy < (j+1)*ydeg+1; yy++) {
                for (int xx = i*xdeg; xx < (i+1)*xdeg+1; xx++) {
                    currNodes.push_back(yy*widthX+xx);
                }
            }
            std::vector<double> bounds = {xOffsets[i],xOffsets[i+1],yOffsets[j],yOffsets[j+1]};
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
        out.push_back(Nodes[node].Position);
    }

    return out;
}

std::vector<std::array<double, 2>> Meshing::BasicMesh::BasicMesh2D::posOfNodes(std::vector<int> NodeIds) {
    std::vector<std::array<double, 2>> out;
    out.reserve(NodeIds.size());

    for (const auto nodeID : NodeIds) {
        out.push_back(Nodes[nodeID].Position);
    }

    return out;
}

std::vector<int> Meshing::BasicMesh::BasicMesh2D::getBoundaryNodes() {
    int numXElems = xOffsets.size();
    int numYElems = yOffsets.size();
    int xWidth = xdeg*numXElems; int yWidth = ydeg*numYElems;
    int numBoundaryNodes = 2*xWidth + 2*yWidth;
    std::vector<int> boundaryNodes; boundaryNodes.reserve(numBoundaryNodes);
    
    for (int i=0; i<xWidth; i++) {
        boundaryNodes[i] = i;
        boundaryNodes[i+xWidth+yWidth] = this->nNodes() - i - 1;
    }

    for (int i=0; i<yWidth; i++) {
        boundaryNodes[i+xWidth] = i*(xWidth+1)+xWidth;
        boundaryNodes[i+2*xWidth+yWidth] = (yWidth-i)*(xWidth+1);
    }
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


