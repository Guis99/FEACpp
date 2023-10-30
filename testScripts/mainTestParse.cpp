#include "../Utils/include/Utils.hpp"

int main() {
    std::string str;
    std::getline(std::cin,str);
    str.erase(remove_if(str.begin(), str.end(), ::isspace), str.end());
    std::cout<<str<<std::endl;

    auto out = Utils::ParseExpression(str);
}