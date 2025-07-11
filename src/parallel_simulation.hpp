#ifndef PARALLEL_SIMULATION_
#define PARALLEL_SIMULATION_

#include <vector>
#include <unordered_set>
#include <optional>

#include <assert.h>
#include <omp.h>
#include <stdio.h>

#include "core_convertor.hpp"
#include "utils.hpp"

using SimulationResult = std::vector<std::pair<SpeciesId, double>>;

struct Fitness {
    // TODO: padding
    int simulation_id;
    double fitness;
    bool error;

    Fitness(int id, double fitness): simulation_id(id), fitness(fitness), error(false) { }

    Fitness(int id, bool error): simulation_id(id), error(error) {}
};

bool comp(const std::pair<SpeciesId, double>& a, const std::pair<SpeciesId, double>& b) {
    return a.second < b.second;
}

class ErrorHandler {
    std::unordered_set<SpeciesId> constrained_species;
    std::vector<std::pair<SpeciesId, double>> real_conc;
protected:

    const std::vector<std::pair<SpeciesId, double>>& get_real_concs() const {
        return this->real_conc;
    }

public:
    virtual ~ErrorHandler() = default;
    virtual double error(SimulationResult& sim_result) const = 0;

    void add_real_concentration(const char *id, double value) {
        control(constrained_species.insert(id).second);
        this->real_conc.push_back(std::pair(id, value));
    }

    void order_real_concentration() {
        std::sort(this->real_conc.begin(), this->real_conc.end(), comp);
    }

    const std::unordered_set<SpeciesId>& get_constrained_species() const {
        return this->constrained_species;
    }
};

class OrderingError: public ErrorHandler {
public:

    ~OrderingError() override = default;

    OrderingError() = default;

    double error(SimulationResult& result) const override {
        std::sort(result.begin(), result.end(), comp);
        auto real_conc = this->get_real_concs();
        control(result.size() == real_conc.size());
        double err = 0.0;
        for (size_t i = 0; i < result.size(); ++i) {
            if(result[i].first != real_conc[i].first) {
                if(result[i].second < -1e-6) {
                    err += 100.0;
                }
                double a = std::log10(std::abs((result[i].second + LITTLE_EPSILON) / (real_conc[i].second + LITTLE_EPSILON)));
                err += a * a;
            }
        }
        return err;
    }

    
};

struct ParameterResult {
    const char *paramater_id;
    double constant;
};


class Simulator {
public:
    // TODO: random start concentration for input that are not known
    virtual ~Simulator() {};
    virtual void set_parameter(const char *id, double value) = 0;
    virtual ParameterResult *get_all_parameters() = 0;
    virtual std::optional<SimulationResult> simulate(const std::unordered_set<SpeciesId>& ids) = 0;

    virtual double simulate_error(const std::unordered_set<SpeciesId>& ids, ErrorHandler *rule) {
        std::optional<SimulationResult> result = this->simulate(ids);
        if(result.has_value()) {
            return rule->error(result.value());
        } else {
            return NAN;
        }
    }
};

class ParallelSimulator {
    // TODO: align, false sharing
    std::vector<Simulator*> sims;

    int workers;

    double ordering_error(SimulationResult& result) const {
        TODO("ordering_error");
    }
public:
    // ATTENTION: we are not freeing the elements
    ~ParallelSimulator() = default;

    ParallelSimulator(int workers): workers(workers) {}

    /***
     * @note we will not free this istance  
     */
    void add_worker(Simulator *sim) {
        this->sims.push_back(sim);
    }

    /**
     * @returns error
     */
    Fitness *simulate(const ErrorHandler *handler) const {
        if(this->workers <= 0) {
            eprintf("[FATAL ERROR] %s:%d: assertion failed this->workers <= 0");
            exit(1);
        }

        omp_set_num_threads(this->workers);

        Fitness *r = (Fitness*)malloc(sizeof(Fitness)*this->sims.size());
        
        #pragma omp parallel for schedule(dynamic)
        for(int i = 0; i < this->sims.size(); ++i) {
            // printf("Thread %d is running\n", omp_get_thread_num());
            std::optional<SimulationResult> result = this->sims[i]->simulate(handler->get_constrained_species());
            // TODO: choose the error
            if(result.has_value()) {
                double error = handler->error(result.value());
                // TODO: false sharing
                r[i] = Fitness(i, error);
            } else {
                r[i] = Fitness(i, true);
            }
        }
        
        return r;
    }
};

#endif // PARALLEL_SIMULATION_