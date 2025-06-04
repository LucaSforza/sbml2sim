#include "sbmldoc.hpp"

extern "C" {

    SBMLDoc *SBMLDoc_new(const char *file_path, bool all_convienece_law) {
        return new SBMLDoc(file_path, all_convienece_law);
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

}