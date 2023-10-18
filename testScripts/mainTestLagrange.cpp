#include "..\Utils\include\LagrangeInterpolant.hpp"
#include <iostream>
#include <vector>

int main() {
    std::cout<<"ayyyyy lmao"<<std::endl;
    int deg; double h;
    std::cout<<"order: ";
    std::cin>>deg;
    std::cout<<std::endl<<"step for derivative: ";
    std::cin>>h;

    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);
    std::vector<double> integ = Utils::integrateLagrange(gaussPoints);

    for (double i : gaussPoints) {
        std::cout<<i<<", ";
    }
    std::cout<<std::endl;

    for (double i : integ) {
        std::cout<<i<<", ";
    }
    std::cout<<std::endl;

    for (int k=0;k<deg+1;k++) {
        std::vector<double> derivs = Utils::numDeriv(h,k,gaussPoints,gaussPoints);
        for (double i : derivs) {
            std::cout<<i<<"|";
        }
        std::cout<<std::endl<<"-------------------------------------------------------------------------------------------"<<std::endl;
    }
    
    return 0;
}