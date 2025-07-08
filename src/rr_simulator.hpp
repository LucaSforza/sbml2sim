#ifndef RR_SIMULATOR_
#define RR_SIMULATOR_

#include <rrRoadRunner.h>

#include "parallel_simulation.hpp"

class rr_Simulator : public Simulator {

    rr::RoadRunner simulator;

public:
    virtual ~rr_Simulator() override = default;

    rr_Simulator(SBMLDoc *doc) {
        simulator(doc->convert_to_sbml_string());
    }

    void set_parameter(const char *id, double value) override {
        this->simulator.setGlobalParameterByName(id, value);
    }

    SimulationResult simulate(const std::unordered_map<SpeciesId>& ids) override {

    }
};

#endif // RR_SIMULATOR_