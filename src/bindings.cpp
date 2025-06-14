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

    void SBMLDoc_dump_proteins_data(const SBMLDoc *_this) {
        _this->dump_proteins_data();
    }
    
    void SBMLDoc_delete(SBMLDoc *_this) {
        delete _this;
    }

    const Proteins* SBMLDoc_get_proteins_data(const SBMLDoc *_this) {
        return &(_this->get_proteins_data());
    }

    Proteins::const_iterator *Proteins_iterator(const Proteins *_this) {
        return new Proteins::const_iterator(_this->cbegin());
    }

    void Proteins_delete_iterator(Proteins::const_iterator *it) {
        delete it;
    }

    std::pair<std::string, std::string> *Proteins_iterator_next(Proteins::const_iterator *it) {
        try {
            auto &pair = *(*it);
            auto *result = new std::pair<std::string, std::string>(pair.first, pair.second);
            ++(*it);
            return result;
        } catch (...) {
            return nullptr;
        }
    }

    void Pair_delete(std::pair<std::string, std::vector<std::string>> *p) {
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
        return SBMLDoc::replicate_model_per_tissue(_this, tissues, n_tissue);
    }

}