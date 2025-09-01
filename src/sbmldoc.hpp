#ifndef SBMLDOC_HPP_
#define SBMLDOC_HPP_

#include "core_convertor.hpp"
#include "avg.hpp"
#include "kinetic_laws.hpp"
#include "utils.hpp"



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
    u_int total_kinetic_constant;
    SpeciesInformation infos;
    Inputs inputs;
    Outputs outputs;

public:

    static SBMLDoc *replicate_model_per_tissue(const char *file_path, const char **tissues, size_t n_tissues) {
        return SBMLDoc::replicate_model_per_tissue(libsbml::readSBML(file_path), tissues, n_tissues);
    }

    // crea le leggi cinetiche se non esistono
    static SBMLDoc *replicate_model_per_tissue(libsbml::SBMLDocument *doc, const char **tissues, size_t n_tissues) {
        SBMLDoc document = SBMLDoc(doc);
        return document.replicate_model_per_tissue(tissues, n_tissues);
    }

    static SBMLDoc *normal_sbml(const char *file_path) {
        SBMLDoc *out = new SBMLDoc(file_path);
        make_all_input_costant_species(out->model, out->inputs);
        return out;
    }

    static SBMLDoc *input_output_normal_species(const char *file_path) {
        SBMLDoc *out = new SBMLDoc(file_path);
        libsbml::Parameter *s = out->model->createParameter();
        s->setId("scale_parameter");
        s->setValue(1.0);
        s->setConstant(true);
        create_a_fake_reaction_for_all_outputs(out->model, out->outputs);
        create_a_fake_reaction_for_all_inputs(out->model, out->inputs);
        return out;
    }

    static SBMLDoc *input_normal_species(const char *file_path) {
        SBMLDoc *out = new SBMLDoc(file_path);
        libsbml::Parameter *s = out->model->createParameter();
        s->setId("scale_parameter");
        s->setValue(1.0);
        s->setConstant(true);
        create_a_fake_reaction_for_all_inputs(out->model, out->inputs);
        return out;
    }

    SBMLDoc() { }

    /**
     * @param file_path path to the .sbml file
    */
    SBMLDoc(const char *file_path): SBMLDoc(libsbml::readSBML(file_path)) {}

    SBMLDoc(libsbml::SBMLDocument *doc) {
        if(check_error(doc)) {
            throw std::runtime_error("Error parsing SBML document");
        }
        this->doc = doc;
        model = this->doc->getModel();
        this->total_kinetic_constant = 0;
        this->infos = register_all_species(model);
        this->inputs = collect_all_inputs(model);
        this->outputs = collect_all_outputs(model);
        // create_a_fake_reaction_for_all_outputs(model, outputs);
        // create_a_fake_reaction_for_all_inputs(model, inputs);
        for (const auto &input : this->inputs) {
            std::cout << "[INPUT] " << input << std::endl;
        }
        for (const auto &output : this->outputs) {
            std::cout << "[OUTPUT] " << output << std::endl;
        }
    }
    
    ~SBMLDoc() {
        delete this->doc;
    }

    void add_kinetic_laws_if_not_exists(bool all_convience_rate_law = false) {
        total_kinetic_constant = add_kinetic_laws(this->model, all_convience_rate_law);
    }

    void add_time_to_model() {
        add_time(this->model);
    }

    void add_avg_calculations_for_all_species() {
        add_avg_calculations(this->model, this->outputs);
    }

    void add_avg_calculation_for_all_proteins() {
        add_avg_calculations_only_for_proteins(this->model, this->infos, this->outputs);
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
        libsbml::Parameter* p = model->getParameter(id);
        if (strncmp( p->getId().c_str(), "output_", 7) != 0) {
            p->setValue(value);
        } else {
            printf("[FATAL ERROR] This is a constant output %s\n", p->getId().c_str());
            exit(1);
        }
    }

    // This is a parameter, so the id is refered to the id of the parameter
    void set_constant_output(int id, double value) {
        libsbml::Parameter* p = model->getParameter(id);
        if (strncmp( p->getId().c_str(), "output_", 7) == 0) {
            p->setValue(value);
        } else {
            printf("[FATAL ERROR] This is a kinetic constant %s\n", p->getId().c_str());
            exit(1);
        }
    }
    

    void input_start_random_concentration() {
        for(const std::string &species_id : this->inputs) {
            libsbml::Species* s = model->getSpecies(species_id);
            assert(s->getSBMLDocument() == this->doc);
            double min_exp = -10;
            double max_exp = -6;
            double scale = static_cast<double>(rand()) / RAND_MAX; // TODO: non è statisticamente affidabile UNIX random engine
            double x = min_exp + scale * (max_exp - min_exp);
            s->setInitialConcentration(pow(10, x));
        }
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
            if(this->outputs.find(species->getId()) != this->outputs.end()) continue;
            double min_exp = -10;
            double max_exp = -6;
            double scale = static_cast<double>(rand()) / RAND_MAX;
            double x = min_exp + scale * (max_exp - min_exp);
            species->setInitialConcentration(pow(10, x));
        }
    }

    void proteins_start_random_concentration() {
        double min_exp = -10;
        double max_exp = -6;
        for (u_int i = 0; i < model->getNumSpecies(); ++i) {
            libsbml::Species *species = model->getSpecies(i);
            if(this->is_protein(species->getId().c_str())) {
                double scale = static_cast<double>(rand()) / RAND_MAX;
                double x = min_exp + scale * (max_exp - min_exp);
                species->setInitialConcentration(pow(10,x));
            }
        }
    }

    void small_compound_start_random_concentration() {
        double min_exp = -10;
        double max_exp = -6;
        for (u_int i = 0; i < model->getNumSpecies(); ++i) {
            libsbml::Species *species = model->getSpecies(i);
            if(this->infos.is_compound(species->getId())) {
                double scale = static_cast<double>(rand()) / RAND_MAX;
                double x = min_exp + scale * (max_exp - min_exp);
                species->setInitialConcentration(pow(10,x));
            }
        }
    }

    void set_zero_output_costant() {
        for (u_int i = 0; i < this->model->getNumParameters(); ++i) {
            libsbml::Parameter *p = this->model->getParameter(i);
            // prendi tutti i parametri output e che non siano non costanti
            if (!p->getConstant() || strncmp(p->getId().c_str(), "output_", 7) != 0)
                continue;

            this->set_constant_output(i, 0);
        }
    }

    void random_kinetic_costant_value() {
        printf("[INFO] total kinetic constant: %d\n", this->total_kinetic_constant);
        int count = 0;
        for (u_int i = 0; i < this->model->getNumParameters(); ++i) {
            libsbml::Parameter *p = this->model->getParameter(i);
            // prendi tutti i parametri non output e che non siano non costanti
            if (!p->getConstant() || strncmp(p->getId().c_str(), "output_", 7) == 0)
                continue;
            count++;
            double min_exp = -6;
            double max_exp = 6;
            double scale = static_cast<double>(rand()) / RAND_MAX;
            double x = min_exp + scale * (max_exp - min_exp);
            this->set_kinetic_constants(i, pow(10, x));
        }
        assert(count == this->total_kinetic_constant);
    }
    
    /**
     * @brief simulate the model and save results
     * @param output_file path for savings the results
     * @param horizon end time of the simulation
     */
    void simulate(const char *output_file, double horizon) {
       eprintf("DEPRECATED simualate SBMLDoc\n");
       exit(1); 
    }

    /**
     * @brief Retrieves the gene data associated with this object.
     *
     * @return A constant reference to the collection of genes.
    */
    const ProteinToId &get_proteins_data(void) const {
        return this->infos.get_protein_to_id();
    }

    const CompoundToId &get_compound_data(void) const {
        return this->infos.get_compound_to_id();
    }

    /**
     * @brief Checks if the given species name corresponds to a protein.
     *
     * Determines whether the specified species name exists in the set of proteins.
     *
     * @param specie_name The name of the species to check.
     * @return true if the species name is found in the set of proteins (i.e., is a protein), false otherwise.
    */
    bool is_protein(const char *specie_name) const {
        return this->infos.is_protein(specie_name);
    }

    u_int get_num_compartements() const {
        return model->getNumCompartments();
    }

    void set_initial_concentration(const char *species_id, double value) {
        libsbml::Species *s = model->getSpecies(species_id);
        if(s == NULL) { 
            eprintf("[FATAL ERROR] species %s doens't exists\n", species_id);
            exit(1);
        }

        assert(s->setInitialConcentration(value) == libsbml::LIBSBML_OPERATION_SUCCESS);
    }

    void set_volume_compartement(u_int id_compartement, double volume) {
        libsbml::Compartment *comp = model->getCompartment(id_compartement);
        if(comp == NULL) throw std::runtime_error("Compartement with id: "+std::to_string(id_compartement)+ "doens't exists");
        if(!comp->isSetSpatialDimensions()) {
            comp->setSpatialDimensions(3.0);
        }
        comp->setVolume(volume);
    }

    const char* get_name_compartement(u_int id_compartement) const {
        libsbml::Compartment *comp = model->getCompartment(id_compartement);
        if(comp == NULL) throw std::runtime_error("Compartement with id: "+std::to_string(id_compartement)+ "doens't exists");
        return comp->getName().c_str();
    }

    void set_parameter(const char *id_parameter, double value) {
        libsbml::Parameter *param = model->getParameter(id_parameter);
        if (param == nullptr) {
            fprintf(stderr, "[FATAL ERROR] Parameter with id '%s' does not exist\n", id_parameter);
            exit(1);
        }
        param->setValue(value);
    }

    void remove_all_assigment_rules() {
        for (unsigned int i = 0; i < model->getNumRules(); ) {
            libsbml::Rule* rule = model->getRule(i);
            if (rule->isAssignment()) {
                model->removeRule(i);
                delete rule;
                // Do not increment i, as the next rule shifts into this index
            } else {
                ++i;
            }
        }
    }

    void assigment_rule_for_inputs() {

        // Create 3 constant parameters: input_constant_f, input_constant_k_1, input_constant_k_2
        std::vector<std::string> param_names = {"f", "k_1", "k_2"};
        std::vector<std::string> param_ids;
        for (const auto &name : param_names) {
            std::string param_id = "input_constant_" + name;
            param_ids.push_back(param_id);
            if (!model->getParameter(param_id)) {
                libsbml::Parameter *param = model->createParameter();
                param->setId(param_id);
                param->setConstant(true);
                param->setValue(0.0); // Default value, can be changed later
            }
        }

        for (const std::string &input_id : this->inputs) {
            libsbml::Species *species = model->getSpecies(input_id);
            if (!species) {
                std::cerr << "[FATAL ERROR] Species not found: " << input_id << std::endl;
                exit(1);
            }
            // Create the assignment rule
            libsbml::AssignmentRule *rule = model->createAssignmentRule();
            rule->setVariable(input_id);
            if(species->getInitialConcentration() == 0.0) {
                eprintf("[FATAL ERROR] initial concentration for input species %s is zero", input_id.c_str());
                exit(1);
            } else if(std::isnan(species->getInitialConcentration())) {
                eprintf("[FATAL ERROR] initial concentration for input species %s is NaN\n", input_id.c_str());
                exit(1);
            }
            std::ostringstream oss;
            oss.precision(15);
            oss << std::scientific << species->getInitialConcentration();
            std::string inital_concentration = oss.str();
            std::string omega = "(2*pi*input_constant_f)";
            std::string phi = "((input_constant_k_1*pi)/4)";
            std::string variance_const = "((input_constant_k_2*"+inital_concentration+")/2)";
            assert(rule->setFormula(inital_concentration+"+"+variance_const+"*sin("+omega+"*get_time"+"+"+phi+")") == libsbml::LIBSBML_OPERATION_SUCCESS);
        }
    }

    SBMLDoc *replicate_model_per_tissue(const char **tissues, size_t n_tissues) const {
        libsbml::SBMLDocument *doc_result = new libsbml::SBMLDocument();
        libsbml::Model *new_model = doc_result->createModel();
        assert(new_model != NULL);

        // crea un compartimento per ogni tessuto
        if(model->getNumCompartments() == 0) TODO("compartments are 0");
        printf("[INFO] clone compartements\n");
        fflush(stdout);
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
                if(comp->isSetSize()) {
                    new_comp->setSize(comp->getSize());
                }
                if(comp->isSetSpatialDimensions()) {
                    new_comp->setSpatialDimensions(comp->getSpatialDimensions());
                }
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
        printf("[INFO] end clone compartements\n");
        fflush(stdout);
        // crea le reazioni facendo attenzione ad assegnare la stessa costante cinetica per le reazioni identiche 
        // ma che si trovano in tessuti diversi

        printf("[INFO] cloning reactions\n");
        fflush(stdout);
        for(u_int i=0; i < model->getNumReactions(); ++i) {
            libsbml::Reaction* r = model->getReaction(i);
            libsbml::KineticLaw* kl = r->getKineticLaw();
            if(kl == NULL) {
                eprintf("[FATAL ERROR] replicate_model_per_tissue: SBML document doen't have kinetic law for reaction %s\n", r->getId().c_str());
                exit(1);
            }
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
                        if(name[0] == 'k' || strncmp(name, "output_", 7) == 0) {
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
                                eprintf("[FATAL ERROR] kinetic constant doen't exists: %s\n", name);
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
        printf("[INFO] end cloning reactions\n");
        fflush(stdout);
        // eprintf("[INFO] generation complete\n");
        // fflush(stderr);

        SBMLDoc *result = new SBMLDoc();
        assert(result != NULL);
        result->doc = doc_result;
        
        result->model = new_model;
        
        result->total_kinetic_constant = this->number_of_kinetic_constants();
        
        result->infos = this->infos.clone_per_tissue(tissues, n_tissues);

        for(const std::string&species_id : this->inputs) {
            for(size_t i=0; i < n_tissues; ++i) {
                std::string tissue = std::string(tissues[i]);
                result->inputs.insert(tissue+"_"+species_id);
            }
        }
        for(const std::string&species_id : this->outputs) {
            for(size_t i=0; i < n_tissues; ++i) {
                std::string tissue = std::string(tissues[i]);
                result->outputs.insert(tissue+"_"+species_id);
            }
        } 
        
        // eprintf("[INFO] returning\n");
        // fflush(stderr);
        return result;
    }

    double get_volume_compartement(const char *id) const {
        return this->model->getCompartment(id)->getSize();
    }

    const char *get_compartement(const char *species_id) const {
        libsbml::Species *s = this->model->getSpecies(species_id);
        if(s == NULL) return NULL;
        return s->getCompartment().c_str();
    }

    std::string convert_to_sbml_string() const {
        return libsbml::writeSBMLToString(this->doc);
    }

    std::vector<std::string> get_kinetic_constants() const {
        std::vector<std::string> kinetic_constants;
        for (u_int i = 0; i < model->getNumParameters(); ++i) {
            libsbml::Parameter* param = model->getParameter(i);
            const std::string& id = param->getId();
            if (!id.empty() && id[0] == 'k') {
                kinetic_constants.push_back(id);
            }
        }
        return kinetic_constants;
    }

    std::vector<std::string> get_output_constants() const {
        std::vector<std::string> output_constants;
        for (u_int i = 0; i < model->getNumParameters(); ++i) {
            libsbml::Parameter* param = model->getParameter(i);
            const std::string& id = param->getId();
            if (!id.empty() && strncmp(id.c_str(), "output_", 7) == 0) {
                output_constants.push_back(id);
            }
        }
        return output_constants;
    }

    std::vector<std::string> get_input_constants() const {
        std::vector<std::string> input_constants;
        for (u_int i = 0; i < model->getNumParameters(); ++i) {
            libsbml::Parameter* param = model->getParameter(i);
            const std::string& id = param->getId();
            if (!id.empty() && strncmp(id.c_str(), "input_", 6) == 0) {
                input_constants.push_back(id);
            }
        }
        return input_constants;
    }

    std::vector<std::string> get_input_species() const {
        std::vector<std::string> input_species;
        input_species.reserve(this->inputs.size());
        for (u_int i = 0; i < model->getNumSpecies(); ++i) {
            libsbml::Species *s = model->getSpecies(i);
            if (s == nullptr) continue;
            if (this->inputs.find(s->getId()) != this->inputs.end()) {
                input_species.push_back(s->getId());
            }
        }
        assert(input_species.size() == this->inputs.size());
        return input_species;
    }

    bool is_output(const char *species_id) const {
        return this->outputs.find(species_id) != this->outputs.end();
    }

    bool is_input(const char *species_id) const {
        return this->inputs.find(species_id) != this->inputs.end();
    }

    void set_outputs_constants(void) {
        for(const std::string& sid : this->outputs) {
            libsbml::Species *s = model->getSpecies(sid);
            s->setBoundaryCondition(true);
        }
    }

    void set_outputs_variable(void) {
        for(const std::string& sid : this->outputs) {
            libsbml::Species *s = model->getSpecies(sid);
            s->setBoundaryCondition(false);
        }
        // TODO: create a method for creating a fake reaction for all outputs
        // create_a_fake_reaction_for_all_outputs(model, outputs);
        add_avg_for_outputs(model, outputs);
    }
};

#endif // SBMLDOC_HPP_