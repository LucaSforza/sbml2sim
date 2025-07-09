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

class Simulator {
public:
    // TODO: random start concentration for input that are not known
    virtual ~Simulator() {};
    virtual void set_parameter(const char *id, double value) = 0;
    virtual std::optional<SimulationResult> simulate(const std::unordered_set<SpeciesId>& ids) = 0;
};

class ParallelSimulator {
    // TODO: align, false sharing
    Simulator **sims;
    std::unordered_set<SpeciesId> species;
    std::vector<std::pair<SpeciesId, double>> real_conc;
    int lenght;
    int workers;

    double ordering_error(SimulationResult& result) const {
        std::sort(result.begin(), result.end(), comp);
        assert(result.size() == real_conc.size());
        double err = 0.0;
        for (size_t i = 0; i < result.size(); ++i) {
            if(result[i].first != real_conc[i].first) {
                double a = std::log10(result[i].second + LITTLE_EPSILON);
                double b = std::log10(real_conc[i].second + LITTLE_EPSILON);
                double diff = a - b;
                err += diff * diff;
            }
        }
        return err;
    }

    static bool comp(const std::pair<SpeciesId, double>& a, const std::pair<SpeciesId, double>& b) {
        return a.second < b.second;
    }
public:
    // ATTENTION: we are not freeing the elements
    ~ParallelSimulator() {
        free(sims);
    }

    ParallelSimulator(int workers): workers(workers), lenght(0) {
        this->sims = (Simulator**)malloc(sizeof(Simulator*)*workers);
    }

    /***
     * @note we will not free this istance  
     */
    void add_worker(Simulator *sim) {
        sims[this->lenght++] = sim;
    }

    void add_real_concentration(const char *id, double value) {
        assert(this->species.insert(id).second);
        this->real_conc.push_back(std::pair(id, value));
    }

    void order_real_concentration() {
        std::sort(this->real_conc.begin(), this->real_conc.end(), comp);
    }

    /**
     * @returns error
     */
    double simulate(int *crashed) const {
        assert(this->lenght == this->workers);
        assert(this->workers > 0);
        int total_crashed = 0;
        double error = 0.0;

        omp_set_num_threads(workers);
        
        #pragma omp parallel for schedule(dynamic) reduction(+: error) reduction(+: total_crashed)
        for(int i = 0; i < this->workers; ++i) {
            std::optional<SimulationResult> result = this->sims[i]->simulate(this->species);
            // TODO: choose the error
            if(result.has_value()) {
                error += this->ordering_error(result.value());
            } else {
                total_crashed++;
            }
        }
        
        *crashed = total_crashed;
        return error;
    }
};

#endif // PARALLEL_SIMULATION_