#include "sbmldoc.hpp"

extern "C" {

    SBMLDoc *SBMLDoc_new(const char *file_path, int flags) {
        return new SBMLDoc(file_path, flags);
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

    void SBMLDoc_dump_genes_data(const SBMLDoc *_this) {
        _this->dump_genes_data();
    }
    
    void SBMLDoc_delete(SBMLDoc *_this) {
        delete _this;
    }

    const Genes* SBMLDoc_get_genes_data(const SBMLDoc *_this) {
        return &(_this->get_genes_data());
    }

    Genes::const_iterator *Genes_proteins_iterator(const Genes *_this) {
        return new Genes::const_iterator(_this->cbegin());
    }

    void Genes_delete_proteins_iterator(Genes::const_iterator *it) {
        delete it;
    }

    std::pair<std::string, std::vector<std::string>> *Genes_proteins_iterator_next(Genes::const_iterator *it) {
        eprintf("[FATAL ERROR] function Genes_proteins_iterator_next is deprecated");
        exit(1);
        // try {
        //     auto &pair = *(*it);
        //     auto *result = new std::pair<std::string, std:std::string>(pair.first, pair.second);
        //     ++(*it);
        //     return result;
        // } catch (...) {
        //     return nullptr;
        // }
    }

    void Pair_delete(std::pair<std::string, std::vector<std::string>> *p) {
        delete p;
    }

    const char* Pair_first_c_str(const std::pair<std::string, std::vector<std::string>> *p) {
        return p->first.c_str();
    }

    const char** Pair_second_as_cstr_array(const std::pair<std::string, std::vector<std::string>> *p, size_t *size) {
        if (size) {
            *size = p->second.size();
        }
        const char **arr = new const char*[p->second.size()];
        for (size_t i = 0; i < p->second.size(); ++i) {
            arr[i] = p->second[i].c_str();
        }
        return arr;
    }

    void Pair_delete_cstr_array(const char **arr) {
        delete[] arr;
    }
    
    bool SBMLDoc_is_protein(const SBMLDoc *_this, const char *specie) {
        return _this->is_protein(specie);
    }

    void SBMLDoc_random_protein_concentrations(SBMLDoc *_this) {
        _this->proteins_start_random_concentration();
    }

    SBMLDoc SBMLDoc_replicate_model_per_tissue(const char *file_path, const char **tissues, size_t n_tissues) {
        return SBMLDoc::replicate_model_per_tissue(file_path, tissues, n_tissues);
    }

}