#ifndef CORE_CONVERTOR_HPP_
#define CORE_CONVERTOR_HPP_

#include <stdio.h>
#include <assert.h>

#include <sbml/SBMLDocument.h>
#include <rr/rrRoadRunner.h>

#include <string>
#include <chrono>
#include <map>
#include <vector>

#define eprintf(...) fprintf(stdout, __VA_ARGS__)
#define TODO(msg) assert(false && msg)

#define SBO_INHIBITOR 20
#define SBO_ACTIVATOR 21
#define SBO_ENZYME 13

#define EPSILON (1e-6)


// species id -> UniProdId
using Proteins = std::map<std::string, std::string>;

Proteins extract_proteins_ids(const libsbml::Model* model) {
    Proteins results;
    const std::string UNIPROT_PREFIX = "https://identifiers.org/uniprot:";

    for (unsigned int i = 0; i < model->getNumSpecies(); ++i) {
        const libsbml::Species* species = model->getSpecies(i);
        const std::string species_id = species->getId();
        const libsbml::XMLNode* annotation = species->getAnnotation();

        bool found = false;
        if (!annotation) continue;

        std::string uniprot_id;
        
        // Cerca direttamente nei nodi XML
        for (unsigned int j = 0; j < annotation->getNumChildren(); ++j) {
            const libsbml::XMLNode& rdfNode = annotation->getChild(j);
            if (rdfNode.getName() != "RDF") continue;

            for (unsigned int k = 0; k < rdfNode.getNumChildren(); ++k) {
                const libsbml::XMLNode& descNode = rdfNode.getChild(k);
                if (descNode.getName() != "Description") continue;
                
                for (unsigned int m = 0; m < descNode.getNumChildren(); ++m) {
                    const libsbml::XMLNode& hasPartNode = descNode.getChild(m);
                    if (hasPartNode.getPrefix() != "bqbiol" || 
                        hasPartNode.getName() != "is") continue;

                    for (unsigned int n = 0; n < hasPartNode.getNumChildren(); ++n) {
                        const libsbml::XMLNode& bagNode = hasPartNode.getChild(n);
                        if (bagNode.getName() != "Bag") continue;

                        for (unsigned int p = 0; p < bagNode.getNumChildren(); ++p) {
                            const libsbml::XMLNode& liNode = bagNode.getChild(p);
                            if (liNode.getName() != "li") continue;

                            const std::string resource = liNode.getAttributes().getValue(0);
                            if (resource.find(UNIPROT_PREFIX) != std::string::npos) {
                                uniprot_id = resource.substr(UNIPROT_PREFIX.length());
                                found = true;
                                break;
                            }
                        }
                        if(found) break;
                    }
                    if(found) break;
                }
                if(found) break;
            }
            if(found) break;
        }

        if (found) {
            results[species_id] = uniprot_id;
        }
    }
    return results;
}

/*
    @return true iff the document has a fatal error
*/
bool check_error(libsbml::SBMLDocument *document) {
    u_int errors;
    bool seriousErrors = false;
    if((errors = document->getNumErrors()) > 0) {
        unsigned int numReadErrors   = 0;
        unsigned int numReadWarnings = 0;
        std::string  errMsgRead      = "";
 
        if (errors > 0) {
            for (unsigned int i = 0; i < errors; i++) {
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
    return seriousErrors;
}

double random_kinetic_constant(void) {
    // Limiti fisici tipici per costanti cinetiche: 1e-6 a 1e6 (adimensionali o in unità appropriate)
    double min_exp = -6.0;
    double max_exp = 6.0;
    double scale = static_cast<double>(rand()) / RAND_MAX;
    double x = min_exp + scale * (max_exp - min_exp);
    return pow(10.0, x);
}

std::string create_hill_pos_function(libsbml::Model *model, libsbml::ModifierSpeciesReference *modifier, u_int h, int *kinetic_constant_added) {
    assert(modifier->getSBOTerm() == SBO_ACTIVATOR);
    assert(h > 0);
    std::string param_k_regulator = "k_activator_"+modifier->getId();
    libsbml::Parameter *p = model->createParameter();
    p->setId(param_k_regulator);
    p->setValue(1.0);
    p->setConstant(true);
    *kinetic_constant_added += 1;
    if(h == 1) {
        return "(("+ modifier->getSpecies() +")/("+param_k_regulator+"+"+modifier->getSpecies()+"))";
    } else {
        std::string h_str = std::to_string(h);
        return "(("+ modifier->getSpecies() +"^"+h_str+")/(("+param_k_regulator+"^"+h_str+")+("+modifier->getSpecies()+"^"+h_str+")))";
    }

}

std::string create_hill_neg_function(libsbml::Model *model, libsbml::ModifierSpeciesReference *modifier, u_int h, int *kinetic_constant_added) {
    assert(modifier->getSBOTerm() == SBO_INHIBITOR);
    assert(h > 0);
    std::string param_k_regulator = "k_inhibitor_"+modifier->getId();
    libsbml::Parameter *p = model->createParameter();
    p->setId(param_k_regulator);
    p->setConstant(true);
    p->setValue(1.0);
    *kinetic_constant_added += 1;
    if(h == 1) {
        return "(("+ param_k_regulator +")/("+param_k_regulator+"+"+modifier->getSpecies()+"))";
    } else {
        std::string h_str = std::to_string(h);
        return "(("+ param_k_regulator +"^"+h_str+")/(("+param_k_regulator+"^"+h_str+")+("+modifier->getSpecies()+"^"+h_str+")))";
    }
}

std::string create_hill_function(libsbml::Model *model, libsbml::ModifierSpeciesReference *modifier, u_int h, int *kinetic_constant_added) {

    int sbo = modifier->getSBOTerm();

    if(sbo == SBO_ACTIVATOR) {
        return create_hill_pos_function(model, modifier, h, kinetic_constant_added);
    } else if(sbo == SBO_INHIBITOR) {
        return create_hill_neg_function(model, modifier, h, kinetic_constant_added);
    } else {
        eprintf("[FATAL ERROR] the modifier is not an activator on an inhibitor, SBO: %d\n", sbo);
        exit(5); 
    }

}

std::string create_modifier_law(libsbml::Model *model, libsbml::Reaction *r, int *kinetic_constant_added) {
    
    if(r->getNumModifiers() <= 0) {
        return "1";
    }

    libsbml::Parameter *p;
    std::string enzyme;
    bool found_enzyme = false;
    std::string regulation;
    bool found_regulator = false;

    for(int i = 0; i < r->getNumModifiers(); i++) {
        libsbml::ModifierSpeciesReference *m = r->getModifier(i);
        int sbo = m->getSBOTerm();
        if(sbo == -1) {
            eprintf("[FATAL ERROR] SBO not present");
            exit(3);
        }
        if(sbo == SBO_ENZYME) {
            // enzyme
            if(found_enzyme) {
                enzyme += "+";
            }
            found_enzyme = true;
            enzyme += m->getSpecies();
        } else if(sbo == SBO_ACTIVATOR || sbo == SBO_INHIBITOR) {
            // regulator
            if(found_regulator) {
                regulation += "*";
            }
            found_regulator = true;
            regulation += create_hill_function(model, m, 1, kinetic_constant_added);
        } else {
            eprintf("[FATAL ERROR] SBO not recognised: %d", sbo);
            exit(4);
        }

    }

    if(!found_enzyme) {
        eprintf("[WARNING] for reaction %s was not found the enzyme\n", r->getId().c_str());
        enzyme = "1";
    }

    if(!found_regulator) {
        regulation = "1";
    }

    return "("+enzyme+")*("+regulation+")";
}

std::string add_convinience_kinetic_law_reversible(libsbml::Model *model, libsbml::Reaction *r, int *kinetic_constant_added) {

    assert(r->getNumProducts() > 0);

    libsbml::Parameter *p;
    std::string law;
    std::string reaction_part;
    std::string product_part;
    std::string den_substrate_part;
    std::string den_product_part;
    
    std::string modifiers = create_modifier_law(model, r, kinetic_constant_added);

    // create the kinetic constant forcat
    std::string paramNameforcat = "k_forcat_"+r->getId();
    p = model->createParameter();
    p->setId(paramNameforcat);
    p->setConstant(true);
    p->setValue(1.0); // default
    reaction_part += paramNameforcat;
    *kinetic_constant_added += 1;
    // create the kinetic constant backcat
    std::string paramNamebackcat = "k_backcat_"+r->getId();
    p = model->createParameter();
    p->setId(paramNamebackcat);
    p->setConstant(true);
    p->setValue(1.0); // default
    product_part += paramNamebackcat;
    *kinetic_constant_added += 1;

    // Create kinetic constants for each substrate
    for (u_int j = 0; j < r->getNumReactants(); ++j) {
        const libsbml::SpeciesReference* sr = r->getReactant(j);
        std::string paramName = "k_" + r->getId() + "_"+ sr->getSpecies();
        p = model->createParameter();
        p->setId(paramName);
        p->setConstant(true);
        p->setValue(1.0);
        *kinetic_constant_added += 1;
        int stoichiometry = sr->getStoichiometry();
        reaction_part += "*(("+sr->getId()+"/"+paramName+")^"+std::to_string(stoichiometry)+")";
        if(j != 0) den_substrate_part += "*";
        den_substrate_part += "( 1 +("+sr->getId()+"/"+paramName+")";
        for(int s = 2; s <= stoichiometry; s++) {
            den_substrate_part += "+(("+sr->getId()+"/"+paramName+")^"+std::to_string(s)+")";
        }
        den_substrate_part += ")";

    }
    // Create kinetic constants for each product
    for (u_int j = 0; j < r->getNumProducts(); ++j) {
        const libsbml::SpeciesReference* sr = r->getProduct(j);
        std::string paramName = "k_" + r->getId() + "_" + sr->getSpecies();
        p = model->createParameter();
        p->setId(paramName);
        p->setConstant(true);
        p->setValue(1.0);
        *kinetic_constant_added += 1;
        int stoichiometry = sr->getStoichiometry();
        product_part += "*(("+sr->getId()+"/"+paramName+")^"+std::to_string(stoichiometry)+")";
        if(j != 0) den_product_part += "*";
        den_product_part += "( 1 +("+sr->getId()+"/"+paramName+")";
        for(int s = 2; s <= stoichiometry; s++) {
            den_product_part += "+(("+sr->getId()+"/"+paramName+")^"+std::to_string(s)+")";
        }
        den_product_part += ")";
    }

    // Create kinetic law

    law = modifiers+"*((("+reaction_part+")" +"-"+ "("+product_part+"))/(("+den_substrate_part+")"+"+ ("+den_product_part+") - 1))";

    return law;
}

std::string add_convinience_kinetic_law_irreversible(libsbml::Model *model, libsbml::Reaction *r, int *kinetic_constant_added) {

    libsbml::Parameter *p;
    std::string law;
    std::string reaction_part;
    std::string den_substrate_part;
    
    std::string modifiers = create_modifier_law(model, r, kinetic_constant_added);

    // create the kinetic constant forcat
    std::string paramNameforcat = "k_forcat_"+r->getId();
    p = model->createParameter();
    p->setId(paramNameforcat);
    p->setValue(1.0); // default
    p->setConstant(true);
    reaction_part += paramNameforcat;
    *kinetic_constant_added += 1;

    // Create kinetic constants for each substrate
    for (u_int j = 0; j < r->getNumReactants(); ++j) {
        const libsbml::SpeciesReference* sr = r->getReactant(j);
        std::string paramName = "k_" + r->getId() + "_"+ sr->getSpecies();
        p = model->createParameter();
        p->setId(paramName);
        p->setValue(1.0); // default value
        p->setConstant(true);
        *kinetic_constant_added += 1;
        int stoichiometry = sr->getStoichiometry();
        reaction_part += "*(("+sr->getId()+"/"+paramName+")^"+std::to_string(stoichiometry)+")";
        if(j != 0) den_substrate_part += "*";
        den_substrate_part += "( 1 +("+sr->getId()+"/"+paramName+")";
        for(int s = 2; s <= stoichiometry; s++) {
            den_substrate_part += "+(("+sr->getId()+"/"+paramName+")^"+std::to_string(s)+")";
        }
        den_substrate_part += ")";

    }

    // Create kinetic law

    law = modifiers+"*(("+reaction_part+")/("+den_substrate_part+"))";

    return law;
}

std::string add_convinience_kinetic_law(libsbml::Model *model, libsbml::Reaction *r, int *kinetic_constant_added) {
    if(r->getReversible()) {
        return add_convinience_kinetic_law_reversible(model, r, kinetic_constant_added);
    } else {
        return add_convinience_kinetic_law_irreversible(model, r, kinetic_constant_added);
    }
}

// this function ignores modifiers
std::string add_michelis_menten_kinetic_law_irreversible(libsbml::Model *model, libsbml::Reaction *r, int *kinetic_constant_added) {
    std::string substrate_part;
    std::string law;
    libsbml::Parameter *p;
    
    std::string paramNameforcat = "k_forcat_"+r->getId();
    p = model->createParameter();
    p->setId(paramNameforcat);
    p->setConstant(true);
    p->setValue(1.0); // default
    substrate_part += paramNameforcat;
    *kinetic_constant_added += 1;

    for(u_int i = 0; i < r->getNumReactants(); i++) {
        const libsbml::SpeciesReference* sr = r->getReactant(i);
        int stoichiometry = sr->getStoichiometry();
        assert(stoichiometry != 0);
        if(stoichiometry == 1) {
            substrate_part += "* "+sr->getSpecies();
        } else {
            substrate_part += "* ("+sr->getSpecies()+"^"+std::to_string(stoichiometry)+")";
        }
    } 

    law = "("+substrate_part + ")";
    return law;
}

// this function ignores modifiers
std::string add_michelis_menten_kinetic_law_reversible(libsbml::Model *model, libsbml::Reaction *r, int *kinetic_constant_added) {
    std::string substrate_part = add_michelis_menten_kinetic_law_irreversible(model, r, kinetic_constant_added);
    std::string product_part;
    std::string law;
    libsbml::Parameter *p;

    std::string paramNamebackcat = "k_backcat_"+r->getId();
    p = model->createParameter();
    p->setId(paramNamebackcat);
    p->setConstant(true);
    p->setValue(1.0); // default
    product_part += paramNamebackcat;
    *kinetic_constant_added += 1;

    for(u_int i = 0; i < r->getNumProducts(); i++) {
        const libsbml::SpeciesReference* sr = r->getProduct(i);
        int stoichiometry = sr->getStoichiometry();
        assert(stoichiometry != 0);
        if(stoichiometry == 1) {
            product_part += "* "+sr->getSpecies();
        } else {
            product_part += "* ("+sr->getSpecies()+"^"+std::to_string(stoichiometry)+")";
        }
    }

    law = "(("+substrate_part + ")-(" + product_part+"))";
    return law;
}

std::string add_michelis_menten_kinetic_law(libsbml::Model *model, libsbml::Reaction *r, int *kinetic_constant_added) {

    std::string law;
    
    if(r->getReversible()) {
        law = add_michelis_menten_kinetic_law_reversible(model, r, kinetic_constant_added);
    } else {
        law = add_michelis_menten_kinetic_law_irreversible(model, r, kinetic_constant_added);
    }

    if(r->getNumModifiers() > 0) {
        assert(r->getNumModifiers() == 1);
        libsbml::ModifierSpeciesReference *modifier = r->getModifier(0);
        return create_hill_function(model, modifier, 10, kinetic_constant_added) + "*" + law;
    } else {
        return law;
    }
    
}

// @return number of kinetic constants added
libsbml::KineticLaw *add_kinetic_law(libsbml::Model *model, libsbml::Reaction *r, bool all_convience_rate_law, int *kinetic_constant_added) {
    std::string law;
    if(r->getNumModifiers() == 0 && !all_convience_rate_law) {
        law = add_michelis_menten_kinetic_law(model, r, kinetic_constant_added);
    } else if (r->getNumModifiers() == 1 && !all_convience_rate_law){
        libsbml::ModifierSpeciesReference *modifier = r->getModifier(0);
        if(modifier->getSBOTerm() == SBO_ACTIVATOR || modifier->getSBOTerm() == SBO_INHIBITOR) {
            law = add_michelis_menten_kinetic_law(model,r, kinetic_constant_added);
        } else {
            law = add_convinience_kinetic_law(model, r, kinetic_constant_added);
        }
    } else {
        law = add_convinience_kinetic_law(model, r, kinetic_constant_added);
    }

    libsbml::KineticLaw *kl = r->createKineticLaw();

    kl->setFormula(law);

    return kl;
}

// @return number of kinetic costant added
int add_kinetic_laws(libsbml::Model *model, bool all_convience_rate_law) {
    int total_kinetic_constant_added = 0;

    for(u_int i = 0; i < model->getNumReactions(); i++) {
        libsbml::Reaction *r = model->getReaction(i);
        libsbml::KineticLaw *kl = r->getKineticLaw();
        if(kl == NULL) {
            add_kinetic_law(model, r, all_convience_rate_law, &total_kinetic_constant_added);
        }
    }

    if(total_kinetic_constant_added != model->getNumParameters()) {
        eprintf("added: %d\n", total_kinetic_constant_added);
        eprintf("real: %d\n", model->getNumParameters());
        fflush(stderr);
        throw std::runtime_error("SBML has more parameters than those added");
    }

    return total_kinetic_constant_added;
}

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

void add_avg_calculations_only_for_proteins(libsbml::Model *model, Proteins &proteins) {
    u_int num_species = model->getNumSpecies();

    for(u_int i = 0; i < num_species; ++i) {
        libsbml::Species *s = model->getSpecies(i);
        
        if(proteins.find(s->getId()) == proteins.end()) continue;
        
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

void add_avg_calculations(libsbml::Model *model) {

    u_int num_species = model->getNumSpecies();

    for(u_int i = 0; i < num_species; ++i) {
        libsbml::Species *s = model->getSpecies(i);
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


#endif // CORE_CONVERTOR_HPP_