#ifndef KINETIC_LAWS_HPP_
#define KINETIC_LAWS_HPP_

#include "core_convertor.hpp"

std::string create_hill_pos_function(libsbml::Model *model, libsbml::ModifierSpeciesReference *modifier, u_int h, int *kinetic_constant_added) {
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

#define SBO_STIMULATOR 459

std::string create_hill_function(libsbml::Model *model, libsbml::ModifierSpeciesReference *modifier, u_int h, int *kinetic_constant_added) {

    int sbo = modifier->getSBOTerm();

    if(sbo == SBO_ACTIVATOR || sbo == SBO_ENZYME || sbo == SBO_STIMULATOR) {
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
        // assert(r->getNumModifiers() == 1);
        std::string hill_component = "";
        for(u_int i =0; i < r->getNumModifiers(); ++i ) {
            if(i != 0) {
                hill_component += "*";
            }
            hill_component += create_hill_function(model, r->getModifier(i), 10, kinetic_constant_added);
        }

        return hill_component + "*" + law;
    } else {
        return law;
    }
    
}

// @return number of kinetic constants added
libsbml::KineticLaw *add_kinetic_law(libsbml::Model *model, libsbml::Reaction *r, bool all_convience_rate_law, int *kinetic_constant_added) {
    std::string law;
    // if(r->getNumModifiers() == 0 && !all_convience_rate_law) {
    //     law = add_michelis_menten_kinetic_law(model, r, kinetic_constant_added);
    // } else if (r->getNumModifiers() == 1 && !all_convience_rate_law){
    //     libsbml::ModifierSpeciesReference *modifier = r->getModifier(0);
    //     if(modifier->getSBOTerm() == SBO_ACTIVATOR || modifier->getSBOTerm() == SBO_INHIBITOR) {
    //         law = add_michelis_menten_kinetic_law(model,r, kinetic_constant_added);
    //     } else {
    //         law = add_convinience_kinetic_law(model, r, kinetic_constant_added);
    //     }
    // } else {
    //     law = add_convinience_kinetic_law(model, r, kinetic_constant_added);
    // }
    law = add_michelis_menten_kinetic_law(model, r, kinetic_constant_added);

    libsbml::KineticLaw *kl = r->createKineticLaw();

    assert(kl->setFormula(law) == libsbml::LIBSBML_OPERATION_SUCCESS);

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
    // printf("num parameters: %d, added: %d\n", model->getNumParameters(), total_kinetic_constant_added);
    // fflush(stdout);
    // if(total_kinetic_constant_added != model->getNumParameters()) {
    //     eprintf("added: %d\n", total_kinetic_constant_added);
    //     eprintf("real: %d\n", model->getNumParameters());
    //     fflush(stderr);
    //     throw std::runtime_error("SBML has more parameters than those added");
    // }

    return total_kinetic_constant_added;
}

#endif // KINETIC_LAWS_HPP_