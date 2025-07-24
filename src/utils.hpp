#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <sbml/SBMLDocument.h>

#include <deque>
#include <vector>
#include <memory>

#define eprintf(...) fprintf(stdout, __VA_ARGS__)
#define TODO(msg) assert(false && msg)

#define control(bool_exp) do {\
    if(!(bool_exp)) { \
        eprintf("[FATAL ERROR] %s:%d: assertion failed: "#bool_exp"\n", __FILE__, __LINE__); \
        exit(1);\
    } \
} while(0)

class MathMLIterator {
    std::deque<libsbml::ASTNode*> frontier;

public:
    MathMLIterator(libsbml::ASTNode *head) {
        this->frontier.push_back(head);
    }

    libsbml::ASTNode *next() {
        if(frontier.empty()) return NULL;
        libsbml::ASTNode *result = frontier[0];
        frontier.pop_front();
        for(u_int i=0; i < result->getNumChildren(); i++) {
            frontier.push_back(result->getChild(i));
        }
        return result;
    }
};


class DFSExplorer {
    std::vector<const libsbml::XMLNode*> frontier;

public:
    DFSExplorer(const libsbml::XMLNode *head) {
        this->frontier.push_back(head);
    }

    const libsbml::XMLNode *next() {
        if(frontier.empty()) return NULL;
        const libsbml::XMLNode *result = frontier.back();
        frontier.pop_back();
        for(u_int i = 0; i < result->getNumChildren(); ++i) {
            this->frontier.push_back(&result->getChild(i));
        }
        return result;
    }
};


#endif // UTILS_HPP_