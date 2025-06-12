#ifndef SBMLDOC_HPP_
#define SBMLDOC_HPP_

#include "core_convertor.hpp"

#include <deque>

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

/**
 * @class SBMLDoc
 * @brief A utility class for loading, manipulating, and simulating SBML (Systems Biology Markup Language) documents.
 *
 * This class provides a high-level interface for reading SBML files, modifying kinetic laws and parameters,
 * generating random initial concentrations for species, simulating the model, and saving the modified SBML document.
 * It also supports extracting and displaying gene associations for species.
 *
 * Key features:
 * - Loads an SBML document from a file and parses its model.
 * - Adds or modifies kinetic laws and calculates the total number of kinetic constants.
 * - Allows setting the value of kinetic constants by index.
 * - Generates random initial concentrations for floating species.
 * - Simulates the SBML model and saves the results to a file.
 * - Extracts and prints the mapping between species and their associated genes.
 * - Saves the modified SBML document to a specified file.
 *
 * Dependencies:
 * - libSBML for SBML parsing and manipulation.
 * - RoadRunner for model simulation.
*/
class SBMLDoc {

    libsbml::SBMLDocument *doc;
    libsbml::Model *model;
    int total_kinetic_constant;
    Genes genes;

    rr::RoadRunner rr;

public:

    static SBMLDoc replicate_model_per_tissue(const char *file_path, const char **tissues, size_t n_tissues) {
        return SBMLDoc::replicate_model_per_tissue(libsbml::readSBML(file_path), tissues, n_tissues);
    }

    static SBMLDoc replicate_model_per_tissue(libsbml::SBMLDocument *doc, const char **tissues, size_t n_tissues) {

        if(check_error(doc)) {
            throw std::runtime_error("Error parsing SBML document");
        }

        libsbml::SBMLDocument *doc_result = new libsbml::SBMLDocument();
        libsbml::Model *new_model = doc_result->createModel();
        assert(new_model != NULL);
        libsbml::Model *model = doc->getModel();

        // crea un compartimento per ogni tessuto
        if(model->getNumCompartments() == 0) TODO("compartments are 0");
        for(u_int i=0; i < model->getNumCompartments(); i++) {
            libsbml::Compartment *comp = model->getCompartment(i);
            for(size_t j=0; j < n_tissues; ++j) {
                // TODO: controlla se i singoli parametri sono impostati
                std::string tissue = std::string(tissues[j]);
                libsbml::Compartment *new_comp = new_model->createCompartment();
                new_comp->setConstant(comp->getConstant());
                new_comp->setName(tissue+" "+comp->getName());
                new_comp->setId(tissue+"_"+comp->getId());
                new_comp->setMetaId(tissue+"_"+comp->getMetaId());
                new_comp->setSBOTerm(comp->getSBOTerm());
                new_comp->setAnnotation(comp->getAnnotationString());
            }
        }

        // crea una specie per ogni tessuto del corpo assegnandogli il giusto compartimento
        for(u_int i=0; i < model->getNumSpecies(); ++i) {
            libsbml::Species *s = model->getSpecies(i);
            for(size_t j=0; j < n_tissues; ++j) {
                std::string tissue = std::string(tissues[j]);
                libsbml::Species *new_s = new_model->createSpecies();
                new_s->setBoundaryCondition(s->getBoundaryCondition());
                new_s->setCompartment(tissue+"_"+s->getCompartment());
                new_s->setConstant(s->getConstant());
                new_s->setHasOnlySubstanceUnits(s->getHasOnlySubstanceUnits());
                new_s->setId(tissue+"_"+s->getId());
                new_s->setMetaId(tissue+"_"+s->getMetaId());
                new_s->setName(tissue+" "+s->getName());
                if(s->isSetNotes()) {
                    new_s->setNotes(s->getNotes());
                }
                new_s->setAnnotation(s->getAnnotationString());
            }
        }

        // crea le reazioni facendo attenzione ad assegnare la stessa costante cinetica per le reazioni identiche 
        // ma che si trovano in tessuti diversi

        int kinetic_constants = 0;

        for(u_int i=0; i < model->getNumReactions(); ++i) {
            libsbml::Reaction* r = model->getReaction(i);
            libsbml::KineticLaw* kl = r->getKineticLaw();
            if(kl == NULL) {
                kl = add_kinetic_law(model, r, false, &kinetic_constants);
            } else {
                TODO("aggiungere le costanti cinetiche rispetto alla legge cinetica già esistente");
            }
            // eprintf("[INFO] reaction: %s\n", r->getId().c_str());
            // fflush(stderr);
            for(size_t j=0; j < n_tissues; ++j) {        
                std::string tissue = std::string(tissues[j]);
                libsbml::Reaction *new_r = new_model->createReaction();
                new_r->setCompartment(tissue+"_"+r->getCompartment());
                if(r->isSetNotes()) {
                    new_r->setNotes(r->getNotes());
                }
                // if(r->isSetAnnotation()) {
                //     new_r->setAnnotation(r->getAnnotation());
                // }
                if(r->isSetFast()) {
                    new_r->setFast(r->getFast());
                }
                new_r->setId(tissue+"_"+r->getId());
                new_r->setMetaId(tissue+"_"+r->getMetaId());
                new_r->setName(tissue+" "+r->getName());
                new_r->setReversible(r->getReversible());
                // eprintf("    reactants\n");
                // fflush(stderr);
                for(u_int reactant=0; reactant < r->getNumReactants(); reactant++) {
                    libsbml::SpeciesReference *sr = r->getReactant(reactant);
                    libsbml::SpeciesReference *new_sr = new_r->createReactant();
                    if(sr->isSetConstant()) {
                        new_sr->setConstant(sr->getConstant());
                    }
                    new_sr->setId(tissue+"_"+sr->getId());
                    new_sr->setSBOTerm(sr->getSBOTerm());
                    new_sr->setSpecies(tissue+"_"+sr->getSpecies());
                    new_sr->setStoichiometry(sr->getStoichiometry());
                }
                // eprintf("    product\n");
                // fflush(stderr);
                for(u_int product=0; product < r->getNumProducts(); product++) {
                    libsbml::SpeciesReference *sr = r->getProduct(product);
                    libsbml::SpeciesReference *new_sr = new_r->createProduct();
                    new_sr->setConstant(sr->getConstant());
                    new_sr->setId(tissue+"_"+sr->getId());
                    new_sr->setSBOTerm(sr->getSBOTerm());
                    new_sr->setSpecies(tissue+"_"+sr->getSpecies());
                    new_sr->setStoichiometry(sr->getStoichiometry());
                }
                // eprintf("    modifier\n");
                // fflush(stderr);
                for(u_int modifier=0;modifier < r->getNumModifiers(); modifier++) {
                    libsbml::ModifierSpeciesReference *sr = r->getModifier(modifier);
                    libsbml::ModifierSpeciesReference *new_sr = new_r->createModifier();
                    new_sr->setId(tissue+"_"+sr->getId());
                    new_sr->setSBOTerm(sr->getSBOTerm());
                    new_sr->setSpecies(tissue+"_"+sr->getSpecies());
                }
                // printf("    done\n");
                // fflush(stderr);

                // aggiungi la legge cinetica
                libsbml::KineticLaw *new_kl = new_r->createKineticLaw(); 
                libsbml::ASTNode *new_head = kl->getMath()->deepCopy();
                MathMLIterator iterator(new_head);
                libsbml::ASTNode *next;
                while((next = iterator.next()) != NULL) {
                    assert(!next->isUnknown());
                    if(next->getType() == libsbml::AST_NAME) {
                        const char *name = next->getName();
                        if(name[0] == 'k') {
                            // printf("[INFO] costante cinetica: %s\n", name);
                            // costante cinetica
                            // crea la costante come parametro del new_model se non esiste
                            libsbml::Parameter *constant = NULL;
                            if((constant = model->getListOfParameters()->get(name)) != NULL) {
                                if(new_model->getListOfParameters()->get(name) == NULL) {
                                    libsbml::Parameter *parameter = new_model->createParameter();
                                    parameter->setId(constant->getId());
                                    parameter->setValue(constant->getValue());
                                    parameter->setConstant(constant->getConstant());
                                }
                            } else {
                                eprintf("[FATAL ERROR] kinetic constant doen't exists: %s", name);
                                exit(1);
                            }
                        } else {
                            // prodotto, substrato oppure modificatore
                            // printf("[INFO] altro: %s\n", name);
                            // modifica il nome appendendo davanti il nome del tessuto
                            next->setName((tissue+"_"+std::string(name)).c_str());
                        }
                    }
                }
                new_kl->setMath(new_head);
            }
        }
        // eprintf("[INFO] generation complete\n");
        // fflush(stderr);

        SBMLDoc result = SBMLDoc();
        result.doc = doc_result;
        result.model = new_model;
        result.total_kinetic_constant = kinetic_constants;
        result.genes = extract_species_genes(new_model);
        add_time(new_model);
        add_avg_calculations_only_for_proteins(new_model, result.genes);


        // eprintf("[INFO] returning\n");
        // fflush(stderr);
        return result;
    }

    enum Flags {
        all_convience_rate_law = 1 << 0,
        avg_for_only_proteins = 1 << 1
    };

    SBMLDoc() { }

    /**
     * @param file_path path to the .sbml file
     * @param all_convience_rate_law true if you want to generate only approximate kinetic laws
    */
    SBMLDoc(const char *file_path, int flags) {
        doc = libsbml::readSBML(file_path);
        if(check_error(doc)) {
            throw std::runtime_error("Error parsing SBML document");
        }
        model = doc->getModel();
        total_kinetic_constant = add_kinetic_laws(model, flags & all_convience_rate_law);
        genes = extract_species_genes(model);
        add_time(model);
        if(flags & avg_for_only_proteins) {
            add_avg_calculations_only_for_proteins(model, genes);
        } else {
            add_avg_calculations(model);
        }

    }
    
    ~SBMLDoc() {
        delete this->doc;
    }
    
    /**
    * @returns number of kinetic constants in the document
    */
    int number_of_kinetic_constants() const {
        return this->total_kinetic_constant;
    }

    /**
     * Sets the value of a kinetic constant parameter in the SBML model.
     *
     * @param id The index of the kinetic constant to set (indices range from 0 to n - 1, where n is obtained using number_of_kinetic_constants()).
     * @param value The new value to assign to the specified kinetic constant.
     *
     * This method updates the value of the kinetic constant parameter identified by the given index in the SBML model.
    */
    void set_kinetic_constants(int id, double value) {
        this->model->getParameter(id)->setValue(value);
    }
    
    /**
     * @brief Saves the converted SBML document to the specified output file.
     * @param output_path The path to the file where the SBML document will be saved.
     * @return true if the file was successfully written; false otherwise.
    */
    bool save_converted_file(const char *output_path) const {
        return libsbml::writeSBML(this->doc, output_path);
    }

    /**
     * @brief for all the species generate a random start concentration
    */
    void random_start_concentration() {
        for (u_int i = 0; i < model->getNumSpecies(); ++i) {
            libsbml::Species *species = model->getSpecies(i);
            // Only set for floating species (not boundary or constant)
            if (!species->getBoundaryCondition() && !species->getConstant()) {
                double min_exp = 0;
                double max_exp = 3.0;
                double scale = static_cast<double>(rand()) / RAND_MAX;
                double x = min_exp + scale * (max_exp - min_exp);
                species->setInitialConcentration(x);
            }
        }
    }

    void proteins_start_random_concentration() {
        double min_exp = 0;
        double max_exp = 3.0;
        for (u_int i = 0; i < model->getNumSpecies(); ++i) {
            libsbml::Species *species = model->getSpecies(i);
            if(!species->getBoundaryCondition() && !species->getConstant() && this->is_protein(species->getId().c_str())) {
                // Only set for floating species (not boundary or constant)
                double scale = static_cast<double>(rand()) / RAND_MAX;
                double x = min_exp + scale * (max_exp - min_exp);
                species->setInitialConcentration(x);
            }
        }
    }
    
    /**
     * @brief simulate the model and save results
     * @param output_file path for savings the results
     * @param horizon end time of the simulation
     */
    void simulate(const char *output_file, double horizon) {
        std::string sbmlStr = libsbml::writeSBMLToStdString(this->doc);
        rr.load(sbmlStr);
        rr::SimulateOptions options = rr::SimulateOptions();
        options.output_file = output_file;
        options.duration = horizon;
        rr.simulate(&options);
    }

    /**
     * @brief Prints the mapping between species and their associated genes to the standard output.
     * @note This function does not take any parameters.
    */
    void dump_genes_data(void) const {
        for (const auto& pair : genes) {
            const std::string& species_id = pair.first;
            std::vector<std::string> gene_ids;
            gene_ids.push_back(pair.second);
            std::cout << "Species: " << species_id << " Genes: ";
            for (size_t i = 0; i < gene_ids.size(); ++i) {
            std::cout << gene_ids[i];
            if (i != gene_ids.size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }
    }

    /**
     * @brief Retrieves the gene data associated with this object.
     *
     * @return A constant reference to the collection of genes.
    */
    const Genes& get_genes_data(void) const {
        return genes;
    }

    /**
     * @brief Checks if the given species name corresponds to a protein.
     *
     * Determines whether the specified species name exists in the set of genes.
     *
     * @param specie_name The name of the species to check.
     * @return true if the species name is found in the set of genes (i.e., is a protein), false otherwise.
    */
    bool is_protein(const char *specie_name) const {
        return genes.find(specie_name) != genes.end();
    }
};

#endif // SBMLDOC_HPP_