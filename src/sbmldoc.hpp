#ifndef SBMLDOC_HPP_
#define SBMLDOC_HPP_

#include "core_convertor.hpp"

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

    /**
     * @param file_path path to the .sbml file
     * @param all_convience_rate_law true if you want to generate only approximate kinetic laws
    */
    SBMLDoc(const char *file_path, bool all_convience_rate_law = false) {
        doc = libsbml::readSBML(file_path);
        if(check_error(doc)) {
            throw std::runtime_error("Error parsing SBML document");
        }
        model = doc->getModel();
        total_kinetic_constant = add_kinetic_laws(model, all_convience_rate_law);
        add_avg_calculations(model, false); // TODO: fai in modo che calcoli l'avg solo per le proteine
        genes = extract_species_genes(model);

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
            const std::vector<std::string>& gene_ids = pair.second;
            std::cout << "Species: " << species_id << " Genes: ";
            for (size_t i = 0; i < gene_ids.size(); ++i) {
            std::cout << gene_ids[i];
            if (i != gene_ids.size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }
    }
};

#endif // SBMLDOC_HPP_