#include <iostream>
#include <algorithm>
#include "..\include\LagrangeInterpolant.hpp"
#include <cmath>

std::vector<double> Utils::genGaussPoints(int degree) {
    std::vector<double> gaussPoints;
    gaussPoints.reserve(degree+1);

    for (int i=0; i < degree+1; i++) {
        gaussPoints.push_back(-cos(i*PI/degree));
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

    // use step size h <= .001 for best results

    std::vector<double> out;
    out.reserve(evalPoints.size());
    double denom = 12*h;

    for (int i = 0; i < evalPoints.size(); i++) {
        double currPoint = evalPoints[i];
        std::vector<double> centeredPoints = {currPoint+2*h,currPoint+h,currPoint-h,currPoint-2*h};
        std::vector<double> funcVals = Utils::evalLagrangeInterp(k, centeredPoints, gaussPoints);

        out.push_back((-funcVals[0]+8*funcVals[1]-8*funcVals[2]+funcVals[3])/denom);
    }

    return out;
}

std::vector<double> Utils::integrateLagrange(std::vector<double> &gaussPoints) {
    // Integrates given basis polynomials with Simpson's rule
    int nInt = gaussPoints.size();

    std::vector<double> evalPoints;
    std::vector<double> out;
    int numDiv = 1000;
    double h = 2/double(numDiv);

    evalPoints.resize(numDiv);
    out.reserve(nInt);

    double currLoc = -1;
    for (int i = 0; i < numDiv; i++) {
        evalPoints[i] = currLoc;
        currLoc += h;
    }
    evalPoints.push_back(1);

    double currInt;
    for (int k = 0; k < nInt; k++) {
        currInt = 0.0;
        std::vector<double> vals = Utils::evalLagrangeInterp(k, evalPoints, gaussPoints);
        for (int i = 0; i < evalPoints.size()-2; i+=2) {
            currInt += h*(vals[i]+4*vals[i+1]+vals[i+2]);
        }
        out.push_back(currInt/3);
    }

    return out;
}

std::vector<double> Utils::ReshuffleNodeVals(std::vector<int> RmOrder, std::vector<int> CwOrder, std::vector<double> shuffleArray) {
    int numObjs = CwOrder.size();
    std::vector<int> shuffledIdxs;
    shuffledIdxs.reserve(numObjs);
    for (int i=0; i<numObjs; i++) {
        int elmToMove = CwOrder[i];
        int lb = 0; int ub = numObjs - 1;
        int mid = (lb+ub)/2;
        while (elmToMove != RmOrder[mid]) {
            mid = (lb+ub)/2;
            if (RmOrder[mid] > elmToMove) {
                ub = mid;
            }
            else {
                lb = mid+1;
            }
        }
        shuffledIdxs.push_back(mid);
    }

    std::vector<double> out;
    out.resize(numObjs);
    for (int i=0; i<numObjs; i++) {
        out[shuffledIdxs[i]] = shuffleArray[i];
    }

    return out;
}