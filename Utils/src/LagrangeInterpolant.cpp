#include <iostream>
#include <algorithm>
#include "..\include\LagrangeInterpolant.hpp"

std::vector<double> Utils::genGaussPoints(int degree) {
    std::vector<double> gaussPoints;
    gaussPoints.reserve(degree+1);

    for (int i=0; i < degree+1; i++) {
        gaussPoints[i] = -cos(i*PI/degree);
    }

    return gaussPoints;
}

std::vector<double> Utils::evalLagrangeInterp(int k, std::vector<double> evalPoints, std::vector<double> &gaussPoints) {
    // Evaluate the kth Lagrange interpolant at the given locations
    int nNodes = gaussPoints.size();
    int nPoints = evalPoints.size();

    std::vector<double> out;
    out.reserve(nPoints);  

    double curr;
    for (int i = 0; i < nPoints; i++) {
        curr = 1.0;
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

std::vector<double> Utils::numDeriv(double h, int k, std::vector<double> &evalPoints, std::vector<double> &gaussPoints) {
    // Fourth-order four-point derivative approximation given by
    // (-f(x + 2h) + 8f(x + h) − 8f(x -h) + f(x 2h))/12h

    std::vector<double> out;
    double denom = 1/(12*h);

    for (int i = 0; i < evalPoints.size(); i++) {
        double currPoint = evalPoints[i];
        std::vector<double> centeredPoints = {currPoint+2*h,currPoint+h,currPoint-h,currPoint-2*h};
        std::vector<double> funcVals = Utils::evalLagrangeInterp(k, centeredPoints, gaussPoints);

        out[i] = (-funcVals[0]+8*funcVals[1]-8*funcVals[2]+funcVals[3])/denom;
    }

    return out;
}

std::vector<double> Utils::integrateLagrange(std::vector<double> &gaussPoints) {
    int nInt = gaussPoints.size();

    std::vector<double> evalPoints;
    std::vector<double> out;
    std::vector<double> h;
    evalPoints.reserve(2*nInt-1);
    out.reserve(nInt);
    h.reserve(nInt-1);

    for (int i = 0; i < nInt-1; i++) {
        evalPoints[2*i] = gaussPoints[i];
        h[i] = (gaussPoints[i+1]-gaussPoints[i])/2;
        evalPoints[2*i+1] = gaussPoints[i]+h[i];
    }
    evalPoints.push_back(gaussPoints[nInt-1]);

    double currInt;
    for (int k = 0; k < nInt; k++) {
        currInt = 0.0;
        std::vector<double> vals = Utils::evalLagrangeInterp(k, evalPoints, gaussPoints);
        for (int i = 0; i < nInt-1; i++) {
            currInt += h[i]*(vals[2*i]+4*vals[2*i+1]+vals[2*i+2]);
        }
        out[k] = currInt/3;
    }

    return out;
}