#include "sbmldoc.hpp"
#include "parallel_simulation.hpp"
#include "rr_simulator.hpp"

extern "C" {

    SBMLDoc *SBMLDoc_new(const char *file_path) {
        return new SBMLDoc(file_path);
    }

    int SBMLDoc_number_of_kinetic_costant(const SBMLDoc *_this) {
        return _this->number_of_kinetic_constants();
    }

    bool SBMLDoc_save_converted_file(const SBMLDoc *_this, const char *output_path) {
        return _this->save_converted_file(output_path);
    }

    void SBMLDoc_set_kinetic_constants(SBMLDoc *_this, int id, double value) {
        _this->set_kinetic_constants(id, value);
    }

    void SBMLDoc_simulate(SBMLDoc *_this, const char *output_file, double duration) {
        _this->simulate(output_file, duration);
    }

    void SBMLDoc_random_start_concentration(SBMLDoc *_this) {
        _this->random_start_concentration();
    }
    
    void SBMLDoc_delete(SBMLDoc *_this) {
        delete _this;
    }

    const ProteinToId* SBMLDoc_get_proteins_data(const SBMLDoc *_this) {
        return &(_this->get_proteins_data());
    }

    const CompoundToId* SBMLDoc_get_compounds_data(const SBMLDoc *_this) {
        return &(_this->get_compound_data());
    }

    SpeciesToId::const_iterator *SpeciesToId_iterator(const SpeciesToId *_this) {
        return new SpeciesToId::const_iterator(_this->cbegin());
    }

    void SpeciesToId_delete_iterator(SpeciesToId::const_iterator *it) {
        delete it;
    }

    std::pair<std::string, std::string> *SpeciesToId_iterator_next(SpeciesToId::const_iterator *it) {
        if(it == nullptr) {
            return nullptr;
        }
        try {
            auto &pair = *(*it);
            auto *result = new std::pair<std::string, std::string>(pair.first, pair.second);
            ++(*it);
            return result;
        } catch (...) {
            return nullptr;
        }
    }

    bool SpeciesToId_iterator_end(const SpeciesToId *_this,SpeciesToId::const_iterator *it) {
        return(_this->cend() == *it);
    }

    void Pair_delete(std::pair<std::string, std::string> *p) {
        delete p;
    }

    const char *Pair_first_c_str(const std::pair<std::string, std::string> *p) {
        return p->first.c_str();
    }

    const char *Pair_second_c_str(const std::pair<std::string, std::string> *p) {
        return p->second.c_str();
    }
    
    bool SBMLDoc_is_protein(const SBMLDoc *_this, const char *specie) {
        return _this->is_protein(specie);
    }

    void SBMLDoc_random_protein_concentrations(SBMLDoc *_this) {
        _this->proteins_start_random_concentration();
    }

    SBMLDoc *replicate_model_per_tissue(const char *file_path, const char **tissues, size_t n_tissues) {
        return SBMLDoc::replicate_model_per_tissue(file_path, tissues, n_tissues);
    }

    SBMLDoc *SBMLDoc_replicate_model_per_tissue(SBMLDoc *_this, const char **tissues, size_t n_tissue) {
        return _this->replicate_model_per_tissue(tissues, n_tissue);
    }

    void SBMLDoc_add_kinetic_laws_if_not_exists(SBMLDoc *_this) {
        _this->add_kinetic_laws_if_not_exists();
    }

    void SBMLDoc_add_time_to_model(SBMLDoc *_this) {
        _this->add_time_to_model();
    }

    void SBMLDoc_add_avg_calculations_for_all_species(SBMLDoc *_this) {
        _this->add_avg_calculations_for_all_species();
    }

    void SBMLDoc_add_avg_calculation_for_all_proteins(SBMLDoc *_this) {
        _this->add_avg_calculation_for_all_proteins();
    }

    void SBMLDoc_random_kinetic_costant_value(SBMLDoc *_this) {
        _this->random_kinetic_costant_value();
    }

    void set_seed(unsigned int seed) {
        srand(seed);
    }

    u_int SBMLDoc_get_num_compartements(const SBMLDoc *_this) {
        return _this->get_num_compartements();
    }

    void SBMLDoc_set_volume_compartement(SBMLDoc *_this, u_int id_compartement, double volume) {
        _this->set_volume_compartement(id_compartement, volume);
    }

    const char *SBMLDoc_get_name_compartement(const SBMLDoc *_this, u_int id_compartement) {
        return _this->get_name_compartement(id_compartement);
    }

    void SBMDoc_small_compound_start_random_concentration(SBMLDoc *_this) {
        _this->small_compound_start_random_concentration();
    }

    double SBMLDoc_get_volume_compartement(const SBMLDoc *_this, const char *id) {
        return _this->get_volume_compartement(id);
    }

    const char *SBMLDoc_get_compartement(const SBMLDoc *_this, const char *species_id) {
        return _this->get_compartement(species_id);
    }

    void SBMLDoc_set_initial_concentration(SBMLDoc *_this, const char *species_id, double value) {
        _this->set_initial_concentration(species_id, value);
    }

    void SBMLDoc_input_start_random_concentration(SBMLDoc *_this) {
        _this->input_start_random_concentration();
    }

    void SBMLDoc_set_zero_output_costant(SBMLDoc *_this) {
        _this->set_zero_output_costant();
    }

    void SBMLDoc_set_parameter(SBMLDoc *_this, const char *id_parameter, double value) {
        _this->set_parameter(id_parameter, value);
    }

    void SBMLDoc_assigment_rule_for_inputs(SBMLDoc *_this) {
        _this->assigment_rule_for_inputs();
    }

    void SBMLDoc_remove_all_assigment_rules(SBMLDoc *_this) {
        _this->remove_all_assigment_rules();
    }

    // DEALLOCATE THE STRING
    char *SBMLDoc_convert_to_sbml_string(const SBMLDoc *_this ) {
        return strdup(_this->convert_to_sbml_string().c_str());
    }

    void deallocate_string(char *string) {
        free(string);
    }

    const std::vector<std::string>* SBMLDoc_get_kinetic_constants(const SBMLDoc* _this) {
        return new std::vector<std::string>(_this->get_kinetic_constants());
    }

    const std::vector<std::string>* SBMLDoc_get_output_constants(const SBMLDoc* _this) {
        return new std::vector<std::string>(_this->get_output_constants());
    }

    void SBMLDoc_delete_string_vector(const std::vector<std::string>* vec) {
        delete vec;
    }

    size_t SBMLDoc_string_vector_size(const std::vector<std::string>* vec) {
        return vec->size();
    }

    const char* SBMLDoc_string_vector_get(const std::vector<std::string>* vec, size_t idx) {
        if (!vec || idx >= vec->size()) return nullptr;
        return (*vec)[idx].c_str();
    }

    bool SBMLDoc_is_output(const SBMLDoc *_this, const char *species_id) {
        return _this->is_output(species_id);
    }

    bool SBMLDoc_is_input(const SBMLDoc *_this, const char *species_id) {
        return _this->is_input(species_id);
    }

    void SBMLDoc_set_outputs_constants(SBMLDoc *_this) {
        _this->set_outputs_constants();
    }

    void SBMLDoc_set_outputs_variable(SBMLDoc *_this) {
        _this->set_outputs_variable();
    }

    Simulator *rr_Simulator_create(const SBMLDoc *doc) {
        return new rr_Simulator(doc);
    }

    void Simulator_delete(Simulator *_this) {
        delete _this;
    }

    void Simulator_set_parameter(Simulator *_this, const char *id, double value) {
        _this->set_parameter(id, value);
    }

    ParallelSimulator *ParallelSimulator_create(int workers) {
        return new ParallelSimulator(workers);
    }

    void ParallelSimulator_delete(ParallelSimulator *_this) {
        delete _this;
    }

    void ParallelSimulator_add_worker(ParallelSimulator *_this, Simulator *worker) {
        _this->add_worker(worker);
    }

    void ErrorHandler_add_real_concentration(ErrorHandler *_this, const char *id, double value) {
        _this->add_real_concentration(id, value);
    }

    void ErrorHandler_order_real_concentration(ErrorHandler *_this) {
        _this->order_real_concentration();
    }

    void ErrorHandler_delete(ErrorHandler *_this) {
        delete _this;
    }

    OrderingError *OrderingError_create() {
        return new OrderingError();
    }

    Fitness *ParallelSimulator_simulate(ParallelSimulator *_this, ErrorHandler *handler) {
        return _this->simulate(handler);
    }

    bool Fitness_is_error(Fitness *_this, int index) {
        return _this[index].error;
    }

    double Fitness_fitness(Fitness *_this, int index) {
        return _this[index].fitness;
    }

    void Fitness_free(Fitness *_this) {
        free(_this);
    }

    // TODO: get_parameter

    // TODO: binding
    /* void unchoosed_ids(const char **ids, size_t num_ids) {
        std::vector<ParameterId> result;
        for (size_t i = 0; i < num_ids; ++i) {
            result.emplace_back(ids[i]);
        }
        unchoosed_ids(result);
    } */

    void Simulator_random_start_concentrations(Simulator *_this) {
        _this->random_start_concentrations();
    }

    void rr_Simulator_set_unknown_id(rr_Simulator *_this, const char *id) {
        _this->set_unknown_id(id);
    }
}