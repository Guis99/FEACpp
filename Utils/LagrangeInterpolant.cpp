#include <iostream>
#include <algorithm>
#include "Utils.hpp"

std::vector<double> Utils::genGaussPoints(int degree) {
    std::vector<double> gaussPoints;
    gaussPoints.reserve(degree+1);
    for (int i=0; i < degree+1; i++) {
        gaussPoints[i] = -cos(i*PI/degree);
    }

    return gaussPoints;
}

std::vector<double> Utils::evalLagrangeInterp(int k, std::vector<double> evalPoints, std::vector<double> &gaussPoints) {
    int nNodes = gaussPoints.size();
    int nPoints = evalPoints.size();

    std::vector<double> out;
    out.reserve(nPoints);  

    for (int i = 0; i < nPoints; i++) {
        double curr = 1;
        double point = evalPoints[i];

        for (int j = 0; j < nNodes; j++) {
            if (j != k) {
            curr *= (point-gaussPoints[j])/(gaussPoints[k]-gaussPoints[j]);
            }
        }
        out[i] = curr;
    }

    return out;

}

std::vector<double> Utils::numDeriv(double h, std::vector<double> evalPoints) {
    std::vector<double> out;

    return out;
}