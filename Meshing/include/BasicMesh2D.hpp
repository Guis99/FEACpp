#include <vector>
#include <array>

namespace Meshing {
    namespace BasicMesh {
        class Element {
            public:
                int EID;
                std::vector<Node> Nodes;
                double Boundaries[4];

                Element(int EID, std::vector<Node> Nodes, double Boundaries[4]);
        };

        class Node {
            public:
                int NID;
                std::array<double, 2> Position;
                int nClass;

                Node(int NID, std::array<double, 2> Position, int nClass);
        };

        class BasicMesh2D {
            public:
                int xdeg;
                int ydeg;
                std::vector<double> xdiv;
                std::vector<double> ydiv;
                std::vector<Node> Nodes;
                std::vector<Element> Elements;
                std::vector<double> xOffsets;
                std::vector<double> yOffsets;

                BasicMesh2D(int xdeg, int ydeg, std::vector<double> xdiv, std::vector<double> ydiv, double xstart, double ystart);

                int numNodes() { return this->Nodes.size(); }
                int numElems() { return this->Elements.size(); }

                std::vector<std::array<double, 2>> allNodePos();
                std::vector<std::array<double, 2>> posInElem(int ElemID);
        };
    }
}