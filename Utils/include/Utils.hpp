#define PI 3.1415926535897932384

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>

enum TokenType {
    NUMBER,
    OPERATOR,
    OPEN_PAREN,
    CLOSE_PAREN,
    INDEP_VARIABLE,
    MATH_FUNCTION
    // Add more categories as needed
};

struct Token {
    TokenType type;
    std::string value;
};
 
namespace Utils {
    std::vector<double> evalLagrangeInterp(int k, std::vector<double> evalPoints, std::vector<double> &gaussPoints);
    std::vector<double> genGaussPoints(int degree);
    std::vector<double> numDeriv(double h, int k, std::vector<double> &evalPoints, std::vector<double> &gaussPoints);
    std::vector<double> integrateLagrange(std::vector<double> &gaussPoints);
    std::vector<double> ReshuffleNodeVals(std::vector<int> RmOrder, std::vector<int> CwOrder, std::vector<double> shuffleArray);
    std::vector<double> ParseExpression(std::string inputString);
    std::vector<Token> TokenizeString(std::string &inputString);
    std::vector<Token> ShuntingYard(std::vector<Token> &tokens);
}