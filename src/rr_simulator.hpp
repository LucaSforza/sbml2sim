#ifndef RR_SIMULATOR_
#define RR_SIMULATOR_

#include <rr/rrRoadRunner.h>

#include "parallel_simulation.hpp"

class rr_Simulator : public Simulator {

    rr::RoadRunner simulator;
    ErrorRule *error_rule = NULL;

public:
    virtual ~rr_Simulator() override = default;

    rr_Simulator(const SBMLDoc *doc): simulator(doc->convert_to_sbml_string()) {
        rr::SimulateOptions options;
        options.start = 0;
        options.duration = 100;
        simulator.setSimulateOptions(options);
        // Seleziona solo le colonne che iniziano per 'avg_'
        // TODO:
        /* std::vector<std::string> allColumns = simulator.getSelectionList();
        std::vector<std::string> selectedColumns;
        for (const auto& col : allColumns) {
            if (col.rfind("avg_", 0) == 0) { // inizia con 'avg_'
                selectedColumns.push_back(col);
            }
        }
        simulator.setSelections(selectedColumns); */
    }

    void set_parameter(const char *id, double value) override {
        this->simulator.setGlobalParameterByName(id, value);
    }

    std::optional<SimulationResult> simulate(const std::unordered_set<SpeciesId>& ids) override {

        const ls::DoubleMatrix *result = NULL;

        const std::string prefix = "avg_";

        try {
            result = simulator.simulate();
        } catch (const std::exception& e) {
            return std::nullopt;
        }

        assert(result != NULL);

        SimulationResult simResult;
        simResult.reserve(ids.size());


        // Ottieni l'indice dell'ultima riga (fine simulazione)
        int lastRow = result->numRows() - 1;

        const std::vector<std::string>& cols = result->getColNames();

        for (int i = 0; i < result->numCols(); ++i) {
            const std::string& col = cols[i];
            if(col.compare(0, prefix.size(), prefix) == 0) {
                const std::string& clean_name = col.substr(prefix.size());
                if(ids.find(clean_name) != ids.end()) {
                    double value = (*result)(lastRow, i);
                    simResult.emplace_back(clean_name, value);
                }
            }
        }

        return simResult;
    } // simulale()


};

class Search_Solutions_Velocity: public rr_Simulator {
    std::unordered_map<ParameterId, bool> choosed_ids;
    std::vector<ParameterId> unchoosed_ids;
public:

    void set_unchoosed_ids(std::vector<ParameterId> ids) {
        this->unchoosed_ids = ids;
    }

    void choose_paramater_velocity(const ParameterId& id, bool velocity) {
        choosed_ids[id] = velocity;
    }
    
    double simulate_error(const std::unordered_set<SpeciesId>& ids, ErrorRule *rule) override {
        // TODO: simula tutte le combinazioni apparte quelle prescelte
        TODO("simulate_error");
    }

};

#endif // RR_SIMULATOR_