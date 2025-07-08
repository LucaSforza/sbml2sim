#ifndef PARALLEL_SIMULATION_
#define PARALLEL_SIMULATION_

#include <vector>
#include <assert.h>
#include <omp.h>

#include "core_convertor.hpp"
#include "utils.hpp"

using Parameters = std::map<std::string, bool>;
using SimulationResult = std::vector<std::pair<SpeciesId, float>>;

class Simulator {
public:
    virtual ~Simulator();
    virtual SimulationResult simulate();
};

class ParallelSimulator {
    // TODO: align, false sharing
    Simulator **sims;
    int lenght;
    int workers;
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

    SimulationResult *simulate() {
        assert(this->lenght == this->workers);
        assert(this->workers > 0);
        SimulationResult *results = (SimulationResult*)malloc(sizeof(SimulationResult)*this->workers);
        #pragma omp parallel for schedule(dynamic)
        for(int i = 0; i < this->workers; ++i) {
            SimulationResult result = this->sims[i]->simulate();
            // false sharing?
            results[i] = result;
        }
        return results;
    }
};

#endif // PARALLEL_SIMULATION_