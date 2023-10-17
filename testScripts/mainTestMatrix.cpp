#include "..\..\3PDep\eigen-3.4.0\Eigen\Eigen"
// #include "..\..\3PDep\eigen-3.4.0\Eigen\LU"
// #include "..\..\3PDep\eigen-3.4.0\Eigen\Sparse"

#include <iostream>
#include <vector>

int main() {
    // Basic linear solve
    Eigen::Matrix<double, 3, 3> A;
    A << .4, .2, 0.0,
        .43, 0.0, .1,
        0.0, -.25, .7;

    Eigen::Matrix<double, 3, 1> b;
    b << 1.6,
        1.82,
        .7;

    Eigen::Matrix<double, 3, 1> x = A.lu().solve(b);

    for (int i = 0; i < 3; i++) {
        std::cout<<x[i]<<std::endl;
    }

    std::cout<<A<<std::endl;
    std::cout<<x<<std::endl;

    // Sparse matrix construction - the spicy stuff
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(20);

    tripletList.emplace_back(0,0,1);
    tripletList.emplace_back(1,1,1);
    tripletList.emplace_back(2,2,1);
    tripletList.emplace_back(3,3,1);
    tripletList.emplace_back(4,4,1);

    tripletList.emplace_back(0,1,1);
    tripletList.emplace_back(1,0,1);
    tripletList.emplace_back(1,2,1);
    tripletList.emplace_back(2,1,1);
    tripletList.emplace_back(2,3,1);
    tripletList.emplace_back(3,2,1);
    tripletList.emplace_back(3,4,1);
    tripletList.emplace_back(4,3,1);

    tripletList.emplace_back(2,2,1);
    tripletList.emplace_back(3,3,1);
    tripletList.emplace_back(2,3,1);
    tripletList.emplace_back(3,2,1);
    tripletList.emplace_back(1,0,1);
    tripletList.emplace_back(1,2,1);
    tripletList.emplace_back(4,3,1);

    Eigen::SparseMatrix<double> mat(5,5);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    std::cout<<mat<<std::endl;

    // Test matrix data import

    float data[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16};
    Eigen::Matrix<float,4,4> bigMat(data);

    for (int i=0;i<4;i++) {
        for (int j=0;j<4;j++) {
            std::cout<<bigMat(i,j)<<" ";
        }
        std::cout<<std::endl;
    }

    // std::cout<<data<<std::endl;
}