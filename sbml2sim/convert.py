import bioutils
import s2s

import sys

PROGRAM_NAME = sys.argv[0]

def parse_args() -> tuple[str, str]:
    if len(sys.argv) < 2:
        print(f"Usage: {PROGRAM_NAME} <input_file_path> [output_file_path]")
        sys.exit(1)
    input_file_path = sys.argv[1]
    output_file_path = None
    if len(sys.argv) >= 3:
        output_file_path = sys.argv[2]
    else:
        output_file_path = "output.sbml"
    return input_file_path, output_file_path

# TODO: aggiungi la possibilità di mettere proteomiche e costanti cinetiche

def main():
    input_file, output_file = parse_args()
    sbml = s2s.SBMLDoc_input_output(input_file)
    volume_cell_breast_cancer_cell = 1.76 * 10**12  # nanometers
    bioutils.set_compartement_size(sbml, volume_cell_breast_cancer_cell)
    sbml.add_kinetic_laws_if_not_exists()
    sbml.add_time_to_model()
    sbml.add_avg_calculations_for_all_species()
    sbml.set_outputs_variable()
    sbml.save_converted_file(output_file)

if __name__ == "__main__":
    main()