#ifndef RR_SIMULATOR_
#define RR_SIMULATOR_


#include <rr/rrRoadRunner.h>
#include <rr/rrLogger.h>

#include "parallel_simulation.hpp"

class rr_Simulator : public Simulator {

    rr::RoadRunner simulator;

    std::unordered_set<SpeciesId> known_ids;

public:
    virtual ~rr_Simulator() override = default;

    rr_Simulator(const SBMLDoc *doc): simulator(doc->convert_to_sbml_string()) {
        rr::SimulateOptions options;
        options.start = 0;
        options.steps = 100;
        options.duration = 100;
        simulator.setSimulateOptions(options);
        rr::Logger::disableConsoleLogging();
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

    ParameterResult *get_all_parameters() override {
        // Ottieni tutti i parametri globali
        std::vector<std::string> allParams = simulator.getGlobalParameterIds();
        std::vector<ParameterResult> selectedParams;

        const std::string prefix = "k_";
        for (const auto& param : allParams) {
            if (param.compare(0, prefix.size(), prefix) == 0) {
                ParameterResult pr;
                pr.paramater_id = strdup(param.c_str());
                pr.constant = simulator.getGlobalParameterByName(param);
                selectedParams.push_back(pr);
            }
        }

        // Alloca array C-style
        ParameterResult *result = nullptr;
        if (!selectedParams.empty()) {
            result = (ParameterResult*)malloc(sizeof(ParameterResult) * selectedParams.size());
            for (size_t i = 0; i < selectedParams.size(); ++i) {
                result[i] = selectedParams[i];
            }
        }
        return result;
    }

    std::expected<SimulationResult,double> simulate(const std::unordered_set<SpeciesId>& ids) override {

        const ls::DoubleMatrix *result = NULL;

        const std::string prefix = "avg_";

        try {
            result = simulator.simulate();
        } catch (const rr::IntegratorException& ie) {
            double crash_time = simulator.getCurrentTime();
            return std::unexpected(crash_time);
        } catch (const std::exception& e) {
            double crash_time = simulator.getCurrentTime();
            eprintf("[SIMULATION ERROR] %s crashed at: %lf\n", e.what(), crash_time);
            exit(1);
        } catch (...) {
            double crash_time = simulator.getCurrentTime();
            eprintf("[SIMULATION ERROR] unknown exception crashed at: %lf\n", crash_time);
            exit(1);
        }


        assert(result != NULL);

        std::unordered_map<SpeciesId, double> constrained;
        constrained.reserve(ids.size());
        std::unordered_map<SpeciesId, double> not_constrained;
        std::unordered_map<SpeciesId, double> old_values;
        // Ottieni l'indice dell'ultima riga (fine simulazione)
        int lastRow = result->numRows() - 1;

        // Prendi l'id quando il tempo è all'80%
        const std::vector<std::string>& cols = result->getColNames();
        int eightyPercentRow = static_cast<int>(result->numRows() * 0.4); // TODO: param
        if (eightyPercentRow >= result->numRows()) eightyPercentRow = result->numRows() - 1;
        // Salva in old_values tutte le specie che iniziano per 'avg_' all'80% del tempo
        for (int i = 0; i < result->numCols(); ++i) {
            const std::string& col = cols[i];
            if (col.compare(0, prefix.size(), prefix) == 0) {
                const std::string& clean_name = col.substr(prefix.size());
                double value = (*result)(eightyPercentRow, i);
                old_values[clean_name] = value;
            }
        }

        for (int i = 0; i < result->numCols(); ++i) {
            const std::string& col = cols[i];
            if(col.compare(0, prefix.size(), prefix) == 0) {
                const std::string& clean_name = col.substr(prefix.size());
                if(ids.find(clean_name) != ids.end()) {
                    double value = (*result)(lastRow, i);
                    constrained[clean_name] = value;
                } else {
                    // TODO: controllare che non sia di output
                    double value = (*result)(lastRow, i);
                    not_constrained[clean_name] = value;
                }
            }
        }

        // Calcola minimi e massimi per ogni "avg_" durante tutta la simulazione
        std::vector<std::pair<SpeciesId, double>> min_values;
        std::vector<std::pair<SpeciesId, double>> max_values;

        int nRows = result->numRows();
        int nCols = result->numCols();
        min_values.reserve(nCols);
        max_values.reserve(nCols);

        for (int c = 0; c < nCols; ++c) {
            const std::string& col = cols[c];
            if (col.compare(0, prefix.size(), prefix) != 0) continue;

            const SpeciesId clean_name = col.substr(prefix.size());
            // Inizializza con il valore nella prima riga
            double v0 = (*result)(0, c);
            double minv = v0;
            double maxv = v0;
            // Scorri tutte le righe per trovare min e max
            for (int r = 1; r < nRows; ++r) {
                double v = (*result)(r, c);
                if (v < minv) minv = v;
                if (v > maxv) maxv = v;
            }
            min_values.emplace_back(clean_name, minv);
            max_values.emplace_back(clean_name, maxv);
        }

        SimulationResult simResult;

        simResult.constrained = constrained;
        simResult.not_constrained = not_constrained;
        simResult.old_values  = old_values;
        simResult.mins = min_values;
        simResult.maxs = max_values; 

        return simResult;
    } // simulale()

    void random_start_concentrations() override {
        const std::vector<std::string>& speciesIds = simulator.getFloatingSpeciesIds();
        for(const std::string& id : speciesIds) {
            if (known_ids.find(id) == known_ids.end()) {
                double min_exp = -10;
                double max_exp = -6;
                double scale = static_cast<double>(rand()) / RAND_MAX;
                double x = min_exp + scale * (max_exp - min_exp);
                simulator.setValue("init(" + id + ")", pow(10, x));
            }
        }
    }
    void set_initial_concentration(const char *species_id, double value) override {
        if (!species_id) return;
        std::string id(species_id);
        // set initial concentration via RoadRunner "init(name)"
        simulator.setValue("init(" + id + ")", value);
        known_ids.insert(id);
    }

    void set_known_id(SpeciesId id) {
        this->known_ids.insert(id);
    }

protected:
    const rr::RoadRunner& get_simulator() {
        return this->simulator;
    }
};


class Search_Solutions_Velocity: public rr_Simulator {
    // not a boolean velocity
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
                results[j][i].constant = best_assigns[j][i] ? 1e3 : 1e-3; // TODO: choose
            }
        }
        return results;
    }
    
    double simulate_error(const std::unordered_set<SpeciesId>& ids,const ErrorHandler *handler, int *errors) override {
        // Imposta anche i parametri già scelti
        for (const auto& [id, val] : choosed_ids) {
            this->set_parameter(id.c_str(), val ? 1e3 : 1e-3);
        }

        std::vector<bool> curr_assign(unchoosed_ids.size(), false);
        
        // Numero di combinazioni: 2^(unchoosed_ids.size())
        size_t num_combinations = 1ULL << unchoosed_ids.size();

        for (size_t mask = 0; mask < num_combinations; ++mask) {
            // Imposta i parametri booleani secondo la combinazione corrente
            for (size_t i = 0; i < unchoosed_ids.size(); ++i) {
                bool value = (mask >> i) & 1;
                curr_assign[i] = value;
                this->set_parameter(unchoosed_ids[i].c_str(), value ? 1e3 : 1e-3);
            }

            std::expected<SimulationResult,double> result = simulate(ids);
            if (result) {
                double error = 1e6*handler->error(*result);
                if(error < 1e-12) { // if(error == 0.0)
                    best_assigns.push_back(curr_assign);
                }
            }
        }
        return 0.0;
    } 

};

#endif // RR_SIMULATOR_