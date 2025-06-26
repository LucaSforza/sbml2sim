#include "sbmldoc.hpp"

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
}