#include "..\..\..\3PDep\eigen-3.4.0\Eigen\Eigen"

#include "..\Utils\Utils.hpp"
#include "..\Meshing\Meshing.hpp"

namespace Solvers {
    namespace MatrixAssembly {
        void MassMatrix();
        void StiffnessMatrix();
        void ConvectionMatrix();
    }
}