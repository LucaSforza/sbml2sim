#ifndef CORE_CONVERTOR_HPP_
#define CORE_CONVERTOR_HPP_

#include <stdio.h>
#include <assert.h>

#include <sbml/SBMLDocument.h>

#include <string>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <optional>

#include "utils.hpp"

#define eprintf(...) fprintf(stdout, __VA_ARGS__)
#define TODO(msg) assert(false && msg)

#define SBO_INHIBITOR 20
#define SBO_ACTIVATOR 21
#define SBO_ENZYME 13

#define EPSILON (1e-6)

using SpeciesToId = std::unordered_map<std::string, std::string>; 
using IdToSpecies = std::unordered_map<std::string, std::vector<std::string>>; 

// species id -> UniProdId
using ProteinToId  = SpeciesToId;
using IdToProtein  = IdToSpecies;
using CompoundToId = SpeciesToId;
using IdToCompound = IdToSpecies;

class SpeciesInformation {
    ProteinToId p_to_id;
    IdToProtein id_to_p;
    CompoundToId c_to_id;
    IdToCompound id_to_c;
public:

    bool is_protein(const std::string& species_id) const {
        return (p_to_id.find(species_id) != p_to_id.end());
    }

    bool is_compound(const std::string& species_id) const {
        return (c_to_id.find(species_id) != c_to_id.end());
    }

    bool is_registered(const std::string& species_id) const {
        return is_protein(species_id) || is_compound(species_id);
    }

    std::optional<std::string> get_compound(const std::string& chebi_id) {
        auto it = id_to_c.find(chebi_id);
        if (it != id_to_c.end() && !it->second.empty()) {
            return it->second[0]; // TODO: generilize and checks for compartement
        }
        return std::nullopt;
    } 

    std::optional<std::string> get_protein(const std::string& uniprod_id) {
        auto it = id_to_p.find(uniprod_id);
        if (it != id_to_p.end() && !it->second.empty()) {
            return it->second[0]; // TODO: generilize and checks for compartement
        }
        return std::nullopt;
    }

    void register_protein(const std::string& species_id, const std::string& uniprot_id) {
        if (is_registered(species_id)) {
            throw std::runtime_error("species is registered, cannot register another time");
        }
        p_to_id[species_id] = uniprot_id;
        id_to_p[uniprot_id].push_back(species_id);
    }

    void register_compound(const std::string& species_id, const std::string& chebi_id) {
        if (is_registered(species_id)) {
            throw std::runtime_error("species is registered, cannot register another time");
        }
        c_to_id[species_id] = chebi_id;
        id_to_c[chebi_id].push_back(species_id);
    }

    // Methods taking libsbml::Species*
    bool is_protein(const libsbml::Species *species) const {
        return is_protein(species->getId());
    }

    bool is_compound(const libsbml::Species *species) const {
        return is_compound(species->getId());
    }

    bool is_registered(const libsbml::Species *species) const {
        return is_registered(species->getId());
    }

    void register_protein(const libsbml::Species *species, const std::string &uniprot_id) {
        register_protein(species->getId(), uniprot_id);
    }

    void register_compound(const libsbml::Species *species, const std::string &chebi_id) {
        register_compound(species->getId(), chebi_id);
    }

    const ProteinToId& get_protein_to_id() const { return p_to_id; }
    const IdToProtein& get_id_to_protein() const { return id_to_p; }
    const CompoundToId& get_compound_to_id() const { return c_to_id; }
    const IdToCompound& get_id_to_compound() const { return id_to_c; }

    SpeciesInformation clone_per_tissue(const char **tissues, size_t n_tissue) const {
        SpeciesInformation cloned;
        for (size_t t = 0; t < n_tissue; ++t) {
            std::string tissue = tissues[t];
            // Clone proteins
            for (const auto& [species_id, uniprot_id] : p_to_id) {
                std::string new_species_id = tissue+"_"+species_id;
                cloned.p_to_id[new_species_id] = uniprot_id;
                cloned.id_to_p[uniprot_id].push_back(new_species_id);
            }
            // Clone compounds
            for (const auto& [species_id, chebi_id] : c_to_id) {
                std::string new_species_id = tissue+"_"+species_id ;
                cloned.c_to_id[new_species_id] = chebi_id;
                cloned.id_to_c[chebi_id].push_back(new_species_id);
            }
        }
        return cloned;
    }
};

enum SpeciesTipology {
    UNKNOWN,
    DRUG,
    PROTEIN,
    COMPOUND,
    COMPLEX,
    SET,
};

bool is_drug(const libsbml::Species *species) {
    if(!species->isSetNotes()) return false;
    std::string notes = species->getNotesString();
    return(notes.find("ProteinDrug") != std::string::npos);
}

bool is_protein(const libsbml::Species *species) {
    if(!species->isSetNotes()) return false;
    std::string notes = species->getNotesString();
    return(notes.find("This is a protein") != std::string::npos);
}

bool is_compound(const libsbml::Species *species) {
    if(!species->isSetNotes()) return false;
    std::string notes = species->getNotesString();
    return(notes.find("This is a small compound") != std::string::npos);
}

bool is_complex(const libsbml::Species *species) {
    if(!species->isSetNotes()) return false;
    std::string notes = species->getNotesString();
    return(notes.find("Reactome Complex") != std::string::npos);
}

bool is_set(const libsbml::Species *species) {
    if(!species->isSetNotes()) return false;
    std::string notes = species->getNotesString();
    return(notes.find("Reactome DefinedSet") != std::string::npos);
}

SpeciesTipology get_tipology(const libsbml::Species *species) {
    if(is_drug(species)) return SpeciesTipology::DRUG;
    if(is_protein(species)) return SpeciesTipology::PROTEIN;
    if(is_compound(species)) return SpeciesTipology::COMPOUND;
    if(is_complex(species)) return SpeciesTipology::COMPLEX;
    if(is_set(species)) return SpeciesTipology::SET;
    return SpeciesTipology::UNKNOWN;
}

const std::string UNIPROT_PREFIX = "https://identifiers.org/uniprot:";
const std::string CHEBI_PREFIX   = "https://identifiers.org/CHEBI:";
const std::string IUPHAR_PREFIX  = "https://identifiers.org/iuphar.ligand:";

bool is_drug(const std::string &url) {
    return url.find(IUPHAR_PREFIX) != std::string::npos;    
}

bool is_protein(const std::string &url) {
    return url.find(UNIPROT_PREFIX) != std::string::npos;
}

bool is_compound(const std::string &url) {
    return url.find(CHEBI_PREFIX) != std::string::npos;
}

SpeciesTipology get_tipology(const std::string &url) {
    if(is_protein(url)) return SpeciesTipology::PROTEIN;
    if(is_compound(url)) return SpeciesTipology::COMPOUND;
    if(is_drug(url)) return SpeciesTipology::DRUG;
    return SpeciesTipology::UNKNOWN;
}

std::string get_drug_id(const std::string &url) {
    return url.substr(IUPHAR_PREFIX.length());
}

std::string get_protein_id(const std::string &url) {
    return url.substr(UNIPROT_PREFIX.length());
}

std::string get_compound_id(const std::string &url) {
    return url.substr(CHEBI_PREFIX.length());
}

std::vector<std::string> extract_information_from_species(const libsbml::Species *species, const std::string& tipology_information) {
    std::vector<std::string> result;
    const libsbml::XMLNode* annotation = species->getAnnotation();
    if(annotation == NULL) return std::vector<std::string>();

    DFSExplorer explorer(annotation);

    const libsbml::XMLNode *next = NULL;

    while((next = explorer.next()) != NULL) {
        if (next->getPrefix() == "bqbiol" && 
            next->getName() == tipology_information) {
                next = explorer.next();
                if(next->getName() == "Bag") {
                    for (u_int p = 0; p < next->getNumChildren(); ++p) {
                        const libsbml::XMLNode& liNode = next->getChild(p);
                        if (liNode.getName() != "li") continue;
                        result.push_back(liNode.getAttributes().getValue(0));
                    }
                    break;
                } else {
                    eprintf("[WARNING] not found the bag for species %s\n", species->getId().c_str());
                }
            }
    }

    return result;
}

std::vector<std::string> extract_is_information_from_species(const libsbml::Species *species) {
    return extract_information_from_species(species, "is");
}

std::vector<std::string> extract_has_part_information_from_species(const libsbml::Species *species) {
    return extract_information_from_species(species, "hasPart");
}

std::string extract_protein_id(const libsbml::Species *species) {
    return extract_is_information_from_species(species)[0].substr(UNIPROT_PREFIX.length());
}

std::string extract_compounds_id(const libsbml::Species *species) {
    return extract_is_information_from_species(species)[0].substr(CHEBI_PREFIX.length());
}

/*
    @return true iff the document has a fatal error
*/
bool check_error(libsbml::SBMLDocument *document) {
    u_int errors;
    bool seriousErrors = false;
    if((errors = document->getNumErrors()) > 0) {
        u_int numReadErrors   = 0;
        u_int numReadWarnings = 0;
        std::string  errMsgRead      = "";
 
        if (errors > 0) {
            for (u_int i = 0; i < errors; i++) {
                if (document->getError(i)->isFatal() || document->getError(i)->isError()) {
                    seriousErrors = true;
                    ++numReadErrors;
                    break;
                } else ++numReadWarnings;
            }
            if(seriousErrors) {
                eprintf("====== PARSING ERROR ======\n\n");
            } else {
                eprintf("====== WARNINGS WHILE PARSING ======\n\n");
            }
            std::ostringstream oss;
            document->printErrors(oss);
            errMsgRead = oss.str();
            std::cout << errMsgRead;
            if(seriousErrors) {
                eprintf("====== END ERRORS ======\n\n");
            } else {
                eprintf("====== END WARNINGS ======\n\n");
            }
            eprintf("errors: %d\nwarnings: %d\n\n", numReadErrors, numReadWarnings);
        }
    }
    return false;
}

double random_kinetic_constant(void) {
    // Limiti fisici tipici per costanti cinetiche: 1e-6 a 1e6 (adimensionali o in unità appropriate)
    double min_exp = -6.0;
    double max_exp = 6.0;
    double scale = static_cast<double>(rand()) / RAND_MAX;
    double x = min_exp + scale * (max_exp - min_exp);
    return pow(10.0, x);
}

// Removes all reactions from the model where the given species appears as a reactant, product, or modifier.
void remove_reactions_with_species(libsbml::Model *model, const std::string& species_id) {
    // Collect indices of reactions to remove
    std::vector<u_int> reactions_to_remove;
    for (u_int i = 0; i < model->getNumReactions(); ++i) {
        libsbml::Reaction* reaction = model->getReaction(i);
        bool found = false;
        // Check reactants
        for (u_int j = 0; j < reaction->getNumReactants(); ++j) {
            if (reaction->getReactant(j)->getSpecies() == species_id) {
                found = true;
                break;
            }
        }
        // Check products
        if (!found) {
            for (u_int j = 0; j < reaction->getNumProducts(); ++j) {
                if (reaction->getProduct(j)->getSpecies() == species_id) {
                    found = true;
                    break;
                }
            }
        }
        // Check modifiers
        if (!found) {
            for (u_int j = 0; j < reaction->getNumModifiers(); ++j) {
                if (reaction->getModifier(j)->getSpecies() == species_id) {
                    // Remove only the modifier, not the whole reaction
                    reaction->removeModifier(j);
                    // Adjust index after removal
                    --j;
                }
            }
        }
        if (found) {
            reactions_to_remove.push_back(i);
        }
    }
    // Remove reactions in reverse order to avoid index shifting
    for (auto it = reactions_to_remove.rbegin(); it != reactions_to_remove.rend(); ++it) {
        model->removeReaction(*it);
    }
}

void eliminate_drugs(libsbml::Model *model) {
    // Iterate over a copy of species IDs to avoid iterator invalidation
    std::vector<std::string> species_ids_to_remove;

    for (u_int i = 0; i < model->getNumSpecies(); ++i) {
        libsbml::Species *s = model->getSpecies(i);
        if (!s) continue;

        // Check for "ProteinDrug" in notes or annotation
        bool is_drug = false;

        // Check notes for "ProteinDrug"
        if (s->isSetNotes()) {
            std::string notes = s->getNotesString();
            if (notes.find("ProteinDrug") != std::string::npos) {
                is_drug = true;
            }
        }

        // Check annotation for IUPHAR ligand (drug) reference
        // TODO: maybe in a DefinedSet not all the reference are drugs
        if (!is_drug && s->isSetAnnotation()) {
            for(std::string& resource : extract_has_part_information_from_species(s)) {
                if (resource.find("https://identifiers.org/iuphar.ligand") != std::string::npos) {
                    is_drug = true;
                    break;
                }
            }
        }

        if (is_drug) {
            species_ids_to_remove.push_back(s->getId());
        }
    }

    // Remove all reactions and the species itself for each drug
    for (const auto& id : species_ids_to_remove) {
        // Remove all reactions where this species appears
        remove_reactions_with_species(model, id);
        // Remove the species itself
        model->removeSpecies(id);
    }
}

bool has_version(const libsbml::Species *s) {
    if(!s->isSetAnnotation()) return false;
    DFSExplorer explorer(s->getAnnotation());
    const libsbml::XMLNode *node = NULL;
    while((node = explorer.next()) != NULL) {
        if(node->getName() == "hasVersion") {
            return true;
        }
    }
    return false;
}

bool is_mutant(const libsbml::Species *s) {
    if(!s->isSetNotes()) return false;
    const std::string& notes = s->getNotesString();
    return notes.find("mutant") != std::string::npos;
}

void register_atomic_species(libsbml::Model *model, SpeciesInformation& info) {
    for (u_int i = 0; i < model->getNumSpecies(); ++i) {
        libsbml::Species *s = model->getSpecies(i);
        if(has_version(s) || is_mutant(s)) continue; // specie fosforata o mutata, non ci interessa
        if(is_protein(s)) {
            std::string id = extract_protein_id(s);
            info.register_protein(s, id);
        } else if(is_compound(s)) {
            std::string id = extract_compounds_id(s);
            info.register_compound(s, id);
        }
    }
}

void delete_species(libsbml::Species *species, libsbml::Model *model) {
    // Remove all reactions where this species appears
    remove_reactions_with_species(model, species->getId());

    // Remove the species itself from the model
    model->removeSpecies(species->getId());
}

void delete_complex(libsbml::Species *complex, libsbml::Model *model) {
    model->removeSpecies(complex->getId());
}

void delete_set(libsbml::Species *set, libsbml::Model *model) {
    model->removeSpecies(set->getId());
}

struct ComplexInformation {
    std::string id;
    SpeciesTipology tipology;
    u_int stoichiometry;
};

std::vector<ComplexInformation> get_complex_information(libsbml::Species *complex) {

    std::vector<std::string> data =  extract_has_part_information_from_species(complex);

    std::vector<ComplexInformation> result;
    if (!complex->isSetNotes()) return result;

    std::string notes = complex->getNotesString();
    // Cerca la parte tra parentesi tonde che contiene la lista degli elementi
    size_t start = notes.find('(');
    size_t end = notes.find(')', start);
    if (start == std::string::npos || end == std::string::npos || end <= start) {
        eprintf("[FATAL ERROR] Complex don't have informations about the complex");
        exit(1);
    }
    std::string content = notes.substr(start + 1, end - start - 1);

    // Split per virgola
    std::stringstream ss(content);
    std::string item;
    while (std::getline(ss, item, ',')) {
        // Rimuovi spazi iniziali/finali
        size_t first = item.find_first_not_of(" \t");
        size_t last = item.find_last_not_of(" \t");
        if (first == std::string::npos) continue;
        item = item.substr(first, last - first + 1);

        // Cerca la forma "2xO14944" oppure solo "O14944"
        size_t x_pos = item.find('x');
        u_int stoichiometry = 1;
        std::string id;
        if (x_pos != std::string::npos && x_pos > 0) {
            // Prima della x c'è la stoichiometria
            std::string stoich_str = item.substr(0, x_pos);
            try {
                stoichiometry = static_cast<u_int>(std::stoul(stoich_str));
            } catch (...) {
                stoichiometry = 1;
            }
            id = item.substr(x_pos + 1);
        } else {
            id = item;
        }
        // Rimuovi eventuali spazi residui
        size_t id_first = id.find_first_not_of(" \t");
        size_t id_last = id.find_last_not_of(" \t");
        if (id_first != std::string::npos)
            id = id.substr(id_first, id_last - id_first + 1);

        if (!id.empty()) {
            // Determina la tipologia usando l'URL corrispondente in data
            if (result.size() >= data.size()) {
                eprintf("[FATAL ERROR] Mismatch between parsed elements and annotation data in complex %s\n", complex->getId().c_str());
                exit(1);
            }
            auto it = std::find_if(result.begin(), result.end(), [&](const ComplexInformation& ci) {
                return ci.id == id;
            });
            if (it != result.end()) {
                it->stoichiometry += stoichiometry;
                continue;
            }
            const std::string& url = data[result.size()];
            SpeciesTipology tipology;
            if (is_protein(url)) {
                tipology = SpeciesTipology::PROTEIN;
                // id = codice dopo l'ultimo ':'
                std::string extracted_id = get_protein_id(url);
                if (id != extracted_id) {
                    eprintf("[FATAL ERROR] Mismatch between id '%s' and url '%s' in complex %s\n", id.c_str(), url.c_str(), complex->getId().c_str());
                    exit(1);
                }
            } else if (is_compound(url)) {
                tipology = SpeciesTipology::COMPOUND;
                std::string extracted_id = get_compound_id(url);
                if (id != extracted_id) {
                    eprintf("[FATAL ERROR] Mismatch between id '%s' and url '%s' in complex %s\n", id.c_str(), extracted_id.c_str(), complex->getId().c_str());
                    exit(1);
                }
            } else {
                eprintf("[FATAL ERROR] Unknown tipology for complex element url: %s\n", url.c_str());
                exit(1);
            }
            result.push_back({id, tipology, stoichiometry});
        } else {
            eprintf("[FATAL ERROR] the string is empty");
            exit(1);
        }
    }

    return result;
}

std::vector<libsbml::Reaction*> get_reaction_where_appear_species(libsbml::Species *species, libsbml::Model *model) {
    std::vector<libsbml::Reaction*> reactions_to_modify;

    // Collect all reactions where the complex appears as reactant, product, or modifier
    for (u_int i = 0; i < model->getNumReactions(); ++i) {
        libsbml::Reaction* reaction = model->getReaction(i);
        bool found = false;
        // Check reactants
        for (u_int j = 0; j < reaction->getNumReactants(); ++j) {
            if (reaction->getReactant(j)->getSpecies() == species->getId()) {
                found = true;
                break;
            }
        }
        // Check products
        if (!found) {
            for (u_int j = 0; j < reaction->getNumProducts(); ++j) {
                if (reaction->getProduct(j)->getSpecies() == species->getId()) {
                    found = true;
                    break;
                }
            }
        }
        // Check modifiers
        if (!found) {
            for (u_int j = 0; j < reaction->getNumModifiers(); ++j) {
                if (reaction->getModifier(j)->getSpecies() == species->getId()) {
                    found = true;
                    break;
                }
            }
        }
        if (found) {
            reactions_to_modify.push_back(reaction);
        }
    }
    return reactions_to_modify;
}


libsbml::Species *create_protein(libsbml::Model *model,const std::string &id, SpeciesInformation &info, const std::string& compartement) {
    libsbml::Species *protein = model->createSpecies();
    protein->setNotes("Genereted by sbml2sim. This is a protein.", true);
    protein->setAnnotation(
    "<annotation><rdf:RDF><rdf:Description><bqbiol:is><rdf:Bag><rdf:li rdf:resource=\""+UNIPROT_PREFIX+id+"\" /></rdf:Bag></bqbiol:is></rdf:Description></rdf:RDF></annotation>");
    protein->setSBOTerm(297);
    protein->setConstant(false);
    protein->setBoundaryCondition(false);
    protein->setHasOnlySubstanceUnits(false);
    protein->setId("species_protein_"+id);
    protein->setCompartment(compartement);
    protein->setName(id);
    info.register_protein(protein, id);
    return protein;
}

libsbml::Species *create_compound(libsbml::Model *model,const std::string &id, SpeciesInformation &info,const std::string& compartement) {
    libsbml::Species *compound = model->createSpecies();
    compound->setNotes("Genereted by sbml2sim. This is a small compound.", true);
    compound->setAnnotation(
    "<rdf:RDF><rdf:Description><bqbiol:is><rdf:Bag><rdf:li rdf:resource=\""+CHEBI_PREFIX+id+"\" /></rdf:Bag></bqbiol:is></rdf:Description></rdf:RDF>");
    compound->setSBOTerm(247);
    compound->setConstant(false);
    compound->setBoundaryCondition(false);
    compound->setHasOnlySubstanceUnits(false);
    compound->setId("species_compound_"+id);
    compound->setCompartment(compartement);
    compound->setName(id);
    info.register_compound(compound, id);
    return compound;
}

void expand_set(libsbml::Species *set, libsbml::Model *model, SpeciesInformation &result) {
    std::vector<std::string> data = extract_has_part_information_from_species(set);
    std::vector<libsbml::Reaction*> reactions_to_modify = get_reaction_where_appear_species(set, model);
    if(reactions_to_modify.size() == 0) {
        // nothing to do
        return;
    }
    assert(data.size() > 0);
    std::vector<std::string> species;
    for(const std::string &url: data) {
        if(is_protein(url)) {
            std::string id = get_protein_id(url);
            std::optional<std::string> p = result.get_protein(id);
            if(!p.has_value()) {
                libsbml::Species *protein = create_protein(model, id, result, set->getCompartment());
                species.push_back(protein->getId());
            } else {
                // TODO: se la proteina è la stessa, ma i compartimenti sono diversi allora bisogna comunque creare una nuova specie
                species.push_back(p.value());
            }
        } else if(is_compound(url)) {
            std::string id = get_compound_id(url);
            std::optional<std::string> c = result.get_compound(id);
            if(!c.has_value()) {
                libsbml::Species *compound = create_compound(model, id, result, set->getCompartment());
                species.push_back(compound->getId());
            } else {
                species.push_back(c.value());
            }
        } else {
            eprintf("[FATAL ERROR] url %s is neither a protein nor a compound\n", url.c_str());
            exit(1);
        }
    }
    assert(species.size() > 0);

    
    for(libsbml::Reaction *reaction : reactions_to_modify) {
        // Conta quante volte il set appare come reactant, product, o modifier
        int reactant_count = 0, product_count = 0, modifier_count = 0;
        for (u_int j = 0; j < reaction->getNumReactants(); ++j)
            if (reaction->getReactant(j)->getSpecies() == set->getId()) reactant_count++;
        for (u_int j = 0; j < reaction->getNumProducts(); ++j)
            if (reaction->getProduct(j)->getSpecies() == set->getId()) product_count++;
        for (u_int j = 0; j < reaction->getNumModifiers(); ++j)
            if (reaction->getModifier(j)->getSpecies() == set->getId()) modifier_count++;

        // Verifica che il set appaia al massimo una volta per reactant/product/modifier
        if (reactant_count + product_count + modifier_count > 1 ) {
            eprintf("[FATAL ERROR] Set %s appears more than once in reaction %s\n", set->getId().c_str(), reaction->getId().c_str());
            exit(1);
        }

        if(modifier_count > 0) {
            // Crea una nuova reazione in cui, per ogni modificatore che è il set,
            // viene sostituito con tutti i membri del set come modificatori.
            libsbml::Reaction* new_reaction = model->createReaction();
            new_reaction->setId(reaction->getId() + "_set_expanded_mod_" + set->getId());

            // Copia reactants
            for (u_int j = 0; j < reaction->getNumReactants(); ++j) {
                const libsbml::SpeciesReference* sr = reaction->getReactant(j);
                libsbml::SpeciesReference* new_sr = new_reaction->createReactant();
                new_sr->setSpecies(sr->getSpecies());
                new_sr->setSBOTerm(sr->getSBOTerm());
                new_sr->setId(sr->getId());
                new_sr->setStoichiometry(sr->getStoichiometry());
            }
            // Copia products
            for (u_int j = 0; j < reaction->getNumProducts(); ++j) {
                const libsbml::SpeciesReference* sr = reaction->getProduct(j);
                libsbml::SpeciesReference* new_sr = new_reaction->createProduct();
                new_sr->setSpecies(sr->getSpecies());
                new_sr->setSBOTerm(sr->getSBOTerm());
                new_sr->setId(sr->getId());
                new_sr->setStoichiometry(sr->getStoichiometry());
            }
            // Copia modifiers, sostituendo il set con tutte le specie del vettore
            for (u_int j = 0; j < reaction->getNumModifiers(); ++j) {
                const libsbml::ModifierSpeciesReference* mr = reaction->getModifier(j);
                if (mr->getSpecies() == set->getId()) {
                    for (const std::string& member : species) {
                        libsbml::ModifierSpeciesReference* new_mr = new_reaction->createModifier();
                        new_mr->setSpecies(member);
                        new_mr->setId(new_reaction->getId() + "_m" + member);
                        new_mr->setSBOTerm(mr->getSBOTerm());
                    }
                } else {
                    libsbml::ModifierSpeciesReference* new_mr = new_reaction->createModifier();
                    new_mr->setSpecies(mr->getSpecies());
                    new_mr->setId(mr->getId());
                    new_mr->setSBOTerm(mr->getSBOTerm());
                }
            }
            new_reaction->setReversible(reaction->getReversible());
        } else {

        // Per ogni specie del set, crea una nuova reazione in cui il set è sostituito dalla specie
        int id = 1;
        for(const std::string& s : species) {
            std::string id_string = std::to_string(id);
            libsbml::Reaction* new_reaction = model->createReaction();
            new_reaction->setId(reaction->getId() + "_set_expanded_" + set->getId()+"_id_"+id_string);

            // Copy reactants
            for (u_int j = 0; j < reaction->getNumReactants(); ++j) {
                const libsbml::SpeciesReference* sr = reaction->getReactant(j);
                if (sr->getSpecies() == set->getId()) {
                    libsbml::SpeciesReference* new_sr = new_reaction->createReactant();
                    new_sr->setSpecies(s);
                    new_sr->setId(new_reaction->getId() + "_" + s);
                    new_sr->setSBOTerm(sr->getSBOTerm());
                    new_sr->setStoichiometry(sr->getStoichiometry());
                } else {
                    libsbml::SpeciesReference* new_sr = new_reaction->createReactant();
                    new_sr->setSpecies(sr->getSpecies());
                    new_sr->setSBOTerm(sr->getSBOTerm());
                    new_sr->setId(sr->getId()+"_r"+id_string);
                    new_sr->setStoichiometry(sr->getStoichiometry());
                }
            }
            // Copy products
            for (u_int j = 0; j < reaction->getNumProducts(); ++j) {
                const libsbml::SpeciesReference* sr = reaction->getProduct(j);
                if (sr->getSpecies() == set->getId()) {
                    libsbml::SpeciesReference* new_sr = new_reaction->createProduct();
                    new_sr->setSpecies(s);
                    new_sr->setId(new_reaction->getId()+"_"+s);
                    new_sr->setSBOTerm(sr->getSBOTerm());
                    new_sr->setStoichiometry(sr->getStoichiometry());
                } else {
                    libsbml::SpeciesReference* new_sr = new_reaction->createProduct();
                    new_sr->setSpecies(sr->getSpecies());
                    new_sr->setId(sr->getId()+"_p"+id_string);
                    new_sr->setSBOTerm(sr->getSBOTerm());
                    new_sr->setStoichiometry(sr->getStoichiometry());
                }
            }
            // Copy reversibility
            new_reaction->setReversible(reaction->getReversible());
            id++;
        }
    }
}

    // Remove the original reactions
    for (libsbml::Reaction* reaction : reactions_to_modify) {
        model->removeReaction(reaction->getId());
    }
}

void eliminate_abstractions(libsbml::Model *model, SpeciesInformation &result) {
    u_int num_species = model->getNumSpecies();
    for (u_int i = 0; i < num_species; ++i) {
        libsbml::Species *s = model->getSpecies(i);
        if(is_set(s)) {
            expand_set(s,model, result);
            delete_set(s,model);
            i--;
            num_species--;
        }
    }
}

void eliminate_useless_species(libsbml::Model *model) {
    // Collect species IDs to remove
    std::vector<std::string> species_to_remove;
    for (u_int i = 0; i < model->getNumSpecies(); ++i) {
        libsbml::Species *species = model->getSpecies(i);
        bool found = false;
        // Check all reactions for this species
        for (u_int j = 0; j < model->getNumReactions(); ++j) {
            libsbml::Reaction *reaction = model->getReaction(j);
            // Reactants
            for (u_int k = 0; k < reaction->getNumReactants(); ++k) {
                if (reaction->getReactant(k)->getSpecies() == species->getId()) {
                    found = true;
                    break;
                }
            }
            // Products
            if (!found) {
                for (u_int k = 0; k < reaction->getNumProducts(); ++k) {
                    if (reaction->getProduct(k)->getSpecies() == species->getId()) {
                        found = true;
                        break;
                    }
                }
            }
            // Modifiers
            if (!found) {
                for (u_int k = 0; k < reaction->getNumModifiers(); ++k) {
                    if (reaction->getModifier(k)->getSpecies() == species->getId()) {
                        found = true;
                        break;
                    }
                }
            }
            if (found) break;
        }
        if (!found) {
            species_to_remove.push_back(species->getId());
        }
    }
    // Remove species
    for (const std::string &id : species_to_remove) {
        model->removeSpecies(id);
    }
}

/**
 * Registra le specie di questo modello SBML, elimina tutte le Drugs da questo modello 
 * Alla prima iterazione in cui si leggono le specie:
 * Se una specie è una normale proteina, viene aggiunta alle proteine
 * 
 * Se una specie è un chemical compound allora viene aggiunto ai compound
 * 
 * Se una specie è un farmaco, allora rimuovilo ed elimina tutte le reazioni in cui appare
 * 
 * Alla seconda iterazione:
 * Se una specie è un Reactome Complex allora prima verifica che non ci siano farmaci in questo complesso
 * Se sono farmaci allora elimina questa specie e tutte le reazioni che compaiono in questa specie
 * Altrimenti per ogni elemento nel complesso:
 *     Verifica se esiste già come specie, altrimenti creala.
 * Poi per ogni reazione in cui appare il complesso sostituiscilo con le singole specie all'interno del complesso
 * 
 * Se una specie è un DefinedSet allora per ogni elemento del set esegui queste operazioni:
 * Se un elemento è un farmaco ignoralo
 * Se un elemento non esiste già allora crealo
 * Per ogni reazione in cui appare questo DefinedSet creane una nuova in cui appare il singolo elemento
 * Finite queste operazioni elemina la specie e elimina le reazioni in cui appare
 */
SpeciesInformation register_all_species(libsbml::Model *model) {
    SpeciesInformation result;

    eliminate_drugs(model);
    eliminate_useless_species(model);
    register_atomic_species(model,result);

    return result;
}

using Inputs = std::unordered_set<std::string>;

bool is_input(const libsbml::Species *species, const libsbml::Model *model) {
    for (u_int i = 0; i < model->getNumReactions(); ++i) {
        const libsbml::Reaction* reaction = model->getReaction(i);
        for (u_int j = 0; j < reaction->getNumProducts(); ++j) {
            if (reaction->getProduct(j)->getSpecies() == species->getId()) {
                return false;
            }
        }
    }
    return true;
}

Inputs collect_all_inputs(const libsbml::Model *model) {
    Inputs result;
    for(u_int i = 0; i < model->getNumSpecies(); ++i) {
        const libsbml::Species *species = model->getSpecies(i);
        if(is_input(species, model)) {
            result.insert(species->getId());
        }
    }
    return result;
}

using Outputs = std::unordered_set<std::string>;

bool is_output(const libsbml::Species *species, const libsbml::Model *model) {
    for (u_int i = 0; i < model->getNumReactions(); ++i) {
        const libsbml::Reaction* reaction = model->getReaction(i);
        for (u_int j = 0; j < reaction->getNumReactants(); ++j) {
            if (reaction->getReactant(j)->getSpecies() == species->getId()) {
                return false;
            }
        }
        // Check if species is a modifier in this reaction
        for (u_int j = 0; j < reaction->getNumModifiers(); ++j) {
            if (reaction->getModifier(j)->getSpecies() == species->getId()) {
                return false;
            }
        }
    }
    return true;
}

Outputs collect_all_outputs(libsbml::Model *model) {
    Outputs result;
    for(u_int i = 0; i < model->getNumSpecies(); ++i) {
        libsbml::Species *species = model->getSpecies(i);
        if(is_output(species, model)) {
            result.insert(species->getId());
            species->setInitialConcentration(0.0);
            species->setBoundaryCondition(true);
        }
    }
    return result;
}

void make_all_input_costant_species(libsbml::Model *model,const Inputs& inputs) {
    for(const std::string& species_id : inputs) {
        libsbml::Species *s = model->getSpecies(species_id);
        assert(s != NULL);
        s->setBoundaryCondition(true);
    }
}

void create_a_fake_reaction_for_all_outputs(libsbml::Model *model, const Outputs& outputs) {
    for(const std::string &species_id : outputs) {
        libsbml::Species *s = model->getSpecies(species_id);
        assert(s != NULL);
        libsbml::Reaction *r = model->createReaction();
        r->setId("degredete_output_"+species_id);
        r->setName("Degredate the output species of the pathway: "+species_id);
        r->setReversible(false);
        r->setCompartment(s->getCompartment());
        libsbml::SpeciesReference* sr = r->createReactant();
        sr->setSpecies(s->getId());
        sr->setId("sr_output_" + species_id);
        sr->setStoichiometry(1.0);
        sr->setSBOTerm(s->getSBOTerm());
        libsbml::KineticLaw* kl = r->createKineticLaw();
        libsbml::Parameter *p = model->createParameter();
        p->setId("output_"+species_id);
        p->setConstant(true);
        p->setValue(1.0); // default
        std::string formula = p->getId() + "*" + sr->getId();
        kl->setFormula(formula);
    }
}

#endif // CORE_CONVERTOR_HPP_