#ifndef TISSUE_HPP_
#define TISSUE_HPP_

#include <unordered_map>
#include <string>
#include <optional>

class Tissue {
    const char *name;

    std::unordered_map<std::string,double> protein_to_concentration;

public:

    Tissue(const char *name): name(name) {}

    std::optional<double> get_protein_concentraion(std::string name) const {
        return std::nullopt;
    }
    
};

#endif // TISSUE_HPP_ 