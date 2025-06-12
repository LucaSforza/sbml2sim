#include "sbmldoc.hpp"

#define shift(argc,argv) (assert(argc > 0), argc--,*(argv++))

int main(int argc, const char **argv) {

    const char *program_name = shift(argc, argv);
    const char *sbml_path = shift(argc, argv);

    printf("[INFO] %s: %s\n",program_name, sbml_path);

    SBMLDoc doc = SBMLDoc::replicate_model_per_tissue(sbml_path, argv, argc);
    doc.save_converted_file("prova.sbml");
    srand(time(NULL));
    doc.random_start_concentration();
    doc.simulate("out.csv", 100.0);
}

