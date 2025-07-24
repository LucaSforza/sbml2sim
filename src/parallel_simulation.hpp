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

struct SimulationResult {
    std::unordered_map<SpeciesId, double> old_values;
    std::vector<std::pair<SpeciesId, double>> constrained;
    std::vector<std::pair<SpeciesId, double>> not_constrained;
};

struct Fitness {
    // TODO: padding
    int simulation_id;
    double fitness;
    bool error;

    Fitness(int id, double fitness): simulation_id(id), fitness(fitness), error(false) { }

    Fitness(int id, bool error): simulation_id(id), error(error) {}

    Fitness(int id, bool error, double fitness): simulation_id(id), fitness(fitness), error(error) {}
};

bool comp(const std::pair<SpeciesId, double>& a, const std::pair<SpeciesId, double>& b) {
    return a.second < b.second;
}

class ErrorHandler {
    std::unordered_set<SpeciesId> constrained_species;
    std::vector<std::pair<SpeciesId, double>> real_conc;
    std::unordered_set<SpeciesId> outputs;
protected:

    const std::vector<std::pair<SpeciesId, double>>& get_real_concs() const {
        return this->real_conc;
    }

    const std::unordered_set<SpeciesId>& get_outputs() const {
        return this->outputs;
    }

public:
    virtual ~ErrorHandler() = default;
    virtual double error(std::expected<SimulationResult,double> sim_result) const = 0;
    virtual void add_output(SpeciesId id) {
        this->outputs.insert(id);
    }

    double get_real_conc(SpeciesId id) const {
        for (const auto& pair : this->real_conc) {
            if (pair.first == id) {
                return pair.second;
            }
        }
        throw std::runtime_error("SpeciesId not found in real_conc");
    }

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

#define WATER_CONC 0.5

class OrderingError: public ErrorHandler {
public:

    ~OrderingError() override = default;

    OrderingError() = default;

    double error(std::expected<SimulationResult,double> sim_result) const override {
        double err = 0.0;
        if(sim_result) {
            SimulationResult simResult = *sim_result;
            std::vector<std::pair<SpeciesId,double>>& result = simResult.constrained;
            std::sort(result.begin(), result.end(), comp);
            auto real_conc = this->get_real_concs();
            control(result.size() == real_conc.size());
            for (size_t i = 0; i < result.size(); ++i) {
                bool not_good = false;    
                if(result[i].second > WATER_CONC || result[i].second < -1e-6) {
                    // unstable
                    not_good = true;
                }
                if(result[i].first != real_conc[i].first || not_good) {
                    double a = std::log10(std::abs((result[i].second + LITTLE_EPSILON) / (this->get_real_conc(result[i].first)+ LITTLE_EPSILON)));
                    if(not_good) {
                        err += 1e4*a*a;
                    } else {
                        err += a * a;
                    }
                }
            }

            const auto& outputs = this->get_outputs();
            
            for (const auto& nc : simResult.not_constrained) {
                if (outputs.find(nc.first) == outputs.end() && (nc.second > WATER_CONC || nc.second < -1e-6)) {
                    // not bounded
                    double a = std::log10(std::abs((nc.second + LITTLE_EPSILON)));
                    err += a;
                }
            } 


        } else {
            // TODO: scegli l'orizzonte di simulazione
            double error = sim_result.error();
            if(error < 1e-12) {
                error = 1e-12;
            }
            err = 1e6*(1/((error)/100.0));
        }
        return err;
    }

    
};

class ClassicalError: public ErrorHandler {
public:
    ~ClassicalError() override = default;

    ClassicalError() = default;

    double error(std::expected<SimulationResult,double> sim_result) const override {
        double err = 0.0;
        if(sim_result) {
            SimulationResult simResult = *sim_result;
            std::vector<std::pair<SpeciesId,double>>& result = simResult.constrained;
            std::sort(result.begin(), result.end(), comp);
            auto real_conc = this->get_real_concs();
            
            control(result.size() == real_conc.size());
            for (size_t i = 0; i < result.size(); ++i) {
                bool not_good = false;    
                if(result[i].second > WATER_CONC || result[i].second < -1e-6) {
                    // unstable
                    not_good = true;
                }
                double a = std::log10(std::abs((result[i].second) / (this->get_real_conc(result[i].first))));
                if(not_good) {
                    err += 1e4*a*a;
                } else {
                    err += a * a;
                }

                double b = (std::log10(simResult.old_values[result[i].first]) - std::log10(result[i].second));
                err += b*b;
            }

            const auto& outputs = this->get_outputs();
            
            for (const auto& nc : simResult.not_constrained) {
                if (outputs.find(nc.first) == outputs.end() && (nc.second > WATER_CONC || nc.second < -1e-6)) {
                    // not bounded
                    double a = std::log10(std::abs((nc.second)));
                    err += a;
                    double b = (std::log10(simResult.old_values[nc.first]) - std::log10(nc.second));
                    err += b*b;
                }
            } 


        } else {
            // TODO: scegli l'orizzonte di simulazione
            double error = sim_result.error();
            if(error < 1e-12) {
                error = 1e-12;
            }
            err = 1e6*(1/((error)/100.0));
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
            int restarts = 1; // TODO: cambiare il numero di restarts

            for (int restart = 0; restart < restarts; ++restart) {
                // this->sims[i]->random_start_concentrations();
                int errors = 0;
                double fitness = this->sims[i]->simulate_error(handler->get_constrained_species(), handler, &errors);
                fitness_sum += fitness;
                error_count += errors;
            }

            double avg_fitness = fitness_sum / (double)restarts;
            r[i] = Fitness(i, error_count > 0, avg_fitness);
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