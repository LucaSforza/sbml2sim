#include "sbmldoc.hpp"

#define shift(argc,argv) (assert(argc > 0), argc--,*(argv++))

int main(int argc, const char **argv) {

    const char *program_name = shift(argc, argv);
    const char *sbml_path = shift(argc, argv);

    printf("[INFO] %s: %s\n",program_name, sbml_path);

    SBMLDoc doc = SBMLDoc(sbml_path);
    doc.save_converted_file("prova.sbml");
    srand(time(NULL));
    doc.random_start_concentration();
    doc.add_kinetic_laws_if_not_exists();
    doc.random_kinetic_costant_value();
    doc.save_converted_file("prova.sbml");
    doc.simulate("out.csv", 1000.0);
}

