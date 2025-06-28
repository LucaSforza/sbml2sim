#ifndef AVG_HPP_
#define AVG_HPP_

#include "core_convertor.hpp"

void add_time(libsbml::Model *model) {
    libsbml::Parameter *time = model->createParameter();
    time->setId("get_time");
    time->setName("Auxiliary variable used in place of time due to naming restrictions");
    time->setConstant(false);
    time->setValue(0.0);
    libsbml::RateRule *time_rule = model->createRateRule();
    time_rule->setVariable("get_time");
    time_rule->setFormula("1");
}

void add_avg_calculations_only_for_proteins(libsbml::Model *model, const SpeciesInformation& infos, Outputs &outputs ) {
    u_int num_species = model->getNumSpecies();

    for(u_int i = 0; i < num_species; ++i) {
        libsbml::Species *s = model->getSpecies(i);
        
        if(!infos.is_protein(s) || outputs.find(s->getId()) != outputs.end() || s->getBoundaryCondition() || s->getConstant()) continue;
        
        std::string avg_param_id = "avg_" + s->getId();
        libsbml::Parameter* avgSpecies = model->createParameter();
        avgSpecies->setId(avg_param_id);
        avgSpecies->setName("Average of " + s->getId());
        if(!std::isnan(s->getInitialConcentration())) {
            avgSpecies->setValue(s->getInitialConcentration());
        } else {
            avgSpecies->setValue(0.0);
        }
        avgSpecies->setConstant(false);
        libsbml::RateRule* avg_rate_rule = model->createRateRule();
        avg_rate_rule->setVariable(avg_param_id);

        // (x - avg_x)/(time + EPSILON)
        avg_rate_rule->setFormula("(" + s->getId() + " - " + avg_param_id+ ")/(get_time + " + std::to_string(EPSILON) + ")");
    }
}

void add_avg_calculations(libsbml::Model *model, Outputs &outputs) {

    u_int num_species = model->getNumSpecies();

    for(u_int i = 0; i < num_species; ++i) {
        libsbml::Species *s = model->getSpecies(i);
        if(s->getBoundaryCondition() || s->getConstant()) continue;
        std::string avg_param_id = "avg_" + s->getId();
        libsbml::Parameter* avgSpecies = model->createParameter();
        avgSpecies->setId(avg_param_id);
        avgSpecies->setName("Average of " + s->getId());
        avgSpecies->setValue(0.0);
        avgSpecies->setConstant(false);
        libsbml::RateRule* avg_rate_rule = model->createRateRule();
        avg_rate_rule->setVariable(avg_param_id);

        // (x - avg_x)/(time + EPSILON)
        avg_rate_rule->setFormula("(" + s->getId() + " - " + avg_param_id+ ")/(get_time + " + std::to_string(EPSILON) + ")");
    }
}

void add_avg_for_outputs(libsbml::Model *model, Outputs &outputs) {
    for (const std::string& output : outputs) {
        const std::string& species_id = output;
        libsbml::Species* s = model->getSpecies(species_id);
        if (!s || s->getBoundaryCondition() || s->getConstant()) {
            eprintf("[FATAL ERROR] the output %s is constant", s->getId().c_str());
            exit(1);
        }

        std::string avg_param_id = "avg_" + species_id;
        libsbml::Parameter* avgSpecies = model->createParameter();
        avgSpecies->setId(avg_param_id);
        avgSpecies->setName("Average of " + species_id);
        avgSpecies->setValue(0.0);
        avgSpecies->setConstant(false);

        libsbml::RateRule* avg_rate_rule = model->createRateRule();
        avg_rate_rule->setVariable(avg_param_id);
        avg_rate_rule->setFormula("(" + s->getId() + " - " + avg_param_id+ ")/(get_time + " + std::to_string(EPSILON) + ")");
    }
}

#endif // AVG_HPP_