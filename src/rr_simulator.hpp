#ifndef RR_SIMULATOR_
#define RR_SIMULATOR_

#include <rr/rrRoadRunner.h>

#include "parallel_simulation.hpp"

class rr_Simulator : public Simulator {

    rr::RoadRunner simulator;

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

protected:
    const rr::RoadRunner& get_simulator() {
        return this->simulator;
    }
};

struct ParameterResult {
    const char *paramater_id;
    bool velocity;
};

class Search_Solutions_Velocity: public rr_Simulator {
    std::unordered_map<ParameterId, bool> choosed_ids;
    std::vector<ParameterId> unchoosed_ids;

    std::unordered_map<ParameterId, bool> best_parameters;
    std::vector<std::vector<bool>> best_assigns;
public:

    void set_unchoosed_ids(std::vector<ParameterId> ids) {
        this->unchoosed_ids = ids;
    }

    void choose_paramater_velocity(const ParameterId& id, bool velocity) {
        choosed_ids[id] = velocity;
    }

    ParameterResult **get_best_assigns() const {
        ParameterResult **results = (ParameterResult**)malloc(sizeof(ParameterResult)*best_assigns.size());
        for(size_t i = 0; i < best_assigns.size(); i++) {
            ParameterResult *one_result = (ParameterResult*)malloc(sizeof(ParameterResult)*unchoosed_ids.size());
            results[i] = one_result;
        }

        for (size_t j = 0; j < best_assigns.size(); j++) {
            for (size_t i = 0; i < unchoosed_ids.size(); ++i) {
                results[j][i].paramater_id = strdup(unchoosed_ids[i].c_str()); // TODO: not duplicate??
                results[j][i].velocity = best_assigns[j][i];
            }
        }
        return results;
    }
    
    double simulate_error(const std::unordered_set<SpeciesId>& ids, ErrorRule *rule) override {
        // Imposta anche i parametri già scelti
        for (const auto& [id, val] : choosed_ids) {
            this->set_parameter(id.c_str(), val ? 1e2 : 1e-2);
        }

        std::vector<bool> curr_assign(unchoosed_ids.size(), false);
        
        // Numero di combinazioni: 2^(unchoosed_ids.size())
        size_t num_combinations = 1ULL << unchoosed_ids.size();

        for (size_t mask = 0; mask < num_combinations; ++mask) {
            // Imposta i parametri booleani secondo la combinazione corrente
            for (size_t i = 0; i < unchoosed_ids.size(); ++i) {
                bool value = (mask >> i) & 1;
                curr_assign[i] = value;
                this->set_parameter(unchoosed_ids[i].c_str(), value ? 1e2 : 1e-2);
            }

            std::optional<SimulationResult> result = simulate(ids);
            if (result.has_value()) {
                double error = 1e6*rule->error(*result);
                if(error < 1e-12) { // if(error == 0.0)
                    best_assigns.push_back(curr_assign);
                }
            }
        }
        return 0.0;
    }

};

#endif // RR_SIMULATOR_