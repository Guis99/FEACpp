#include "..\include\Utils.hpp"

#include <regex>
#include <stack>

std::vector<double> Utils::ParseExpression(std::string inputString) {
    
    std::vector<double> out;
    auto tokens = Utils::TokenizeString(inputString);

    for (auto str : tokens) {
        std::cout<<str.value<<std::endl;
    }

    std::cout<<tokens.size()<<std::endl;
    return out;
}

std::vector<Token> Utils::TokenizeString(std::string &inputString) {
    std::vector<Token> tokens;
    std::string variablePattern = "[xy]";
    std::string mathPattern = "sin|cos|tan";
    std::string regexString = R"((\d+)|([-+*/^])|([\(\{\[])|([\)\}\]])|()" + 
                                variablePattern + R"()|()" + 
                                mathPattern + R"())";
    
    std::cout<<regexString<<std::endl;
    std::regex regexPattern(regexString);
    std::smatch match;

    auto pos = inputString.cbegin();
    while (std::regex_search(pos, inputString.cend(), match, regexPattern)) {
        for (size_t i = 1; i < match.size(); ++i) {
            if (!match[i].str().empty()) {
                Token token;
                token.value = match[i].str();
                switch (i) {
                    case 1: token.type = NUMBER; break;
                    case 2: token.type = OPERATOR; break;
                    case 3: token.type = OPEN_PAREN; break;
                    case 4: token.type = CLOSE_PAREN; break;
                    case 5: token.type = INDEP_VARIABLE; break;
                    case 6: token.type = MATH_FUNCTION; break;
                }
                tokens.push_back(token);
            }
        }
        pos = match.suffix().first;
    }
    return tokens;
}

std::vector<Token> Utils::ShuntingYard(std::vector<Token> &tokens) {
    std::stack<double> valStack;
    std::stack<Token> opStack;

    for (auto token : tokens) {
        switch (token.type) {
            case NUMBER: break;
            case OPERATOR: break;
            case OPEN_PAREN: break;
            case CLOSE_PAREN: break;
            case INDEP_VARIABLE: break;
            case MATH_FUNCTION: break;
        }
    }



    return RpnVec;
}

