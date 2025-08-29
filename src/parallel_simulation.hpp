#ifndef PARALLEL_SIMULATION_
#define PARALLEL_SIMULATION_

#include <vector>
#include <unordered_set>
#include <optional>
#include <expected>

#include <assert.h>
#include <omp.h>
#include <stdio.h>

#include "core_convertor.hpp"
#include "utils.hpp"

#include "parallel_simulation.hpp"
#include "loss_function.hpp"


struct ParameterResult {
    const char *paramater_id;
    double constant;
};


class Simulator {
public:
    virtual ~Simulator() {};
    virtual void set_parameter(const char *id, double value) = 0;
    virtual void set_initial_concentration(const char *species_id, double value) = 0;
    virtual ParameterResult *get_all_parameters() = 0;
    // ritorna un simulation result oppure il tempo in cui il simulatore va in crash
    virtual std::expected<SimulationResult,double> simulate(const std::unordered_set<SpeciesId>& ids) = 0;
    virtual double simulate_error(const std::unordered_set<SpeciesId>& ids,const ErrorHandler *handler, int *errors) {
        std::expected<SimulationResult,double> result = this->simulate(ids);
        if(result) {
            *errors = 0;
        } else {
            *errors = 1;
        }
        return handler->error(result);
    }

    virtual void random_start_concentrations() = 0;
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
            int errors = 0;
            double fitness_sum = 0.0;
            int error_count = 0;
            auto result = this->sims[i]->simulate(handler->get_constrained_species());
            if(!result) error_count += 1;
            
            fitness_sum += handler->error(result);
    
            r[i] = Fitness(i, error_count > 0, fitness_sum);
        }
        
        return r;
    }

    void simulate_only(const ErrorHandler *handler) const {
        if(this->workers <= 0) {
            eprintf("[FATAL ERROR] %s:%d: assertion failed this->workers <= 0");
            exit(1);
        }

        omp_set_num_threads(this->workers);

        #pragma omp parallel for schedule(dynamic)
        for(int i = 0; i < this->sims.size(); ++i) {
            this->sims[i]->simulate_error(handler->get_constrained_species(), handler, NULL);
        }
    }
};

#endif // PARALLEL_SIMULATION_