#ifndef RR_SIMULATOR_
#define RR_SIMULATOR_

#include <rrRoadRunner.h>

#include "parallel_simulation.hpp"

class rr_Simulator : public Simulator {

    rr::RoadRunner simulator; 

public:
    virtual ~rr_Simulator() override = default;
    SimulationResult simulate(Parameters& param) {

    }
};

#endif // RR_SIMULATOR_