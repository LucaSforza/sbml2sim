from abc import ABC, abstractmethod
import nevergrad as ng

import time

from typing import Any

from bioutils import ParameterId

import s2s

class Optimizer(ABC):

    def __init__(self, model: s2s.Simulator, simulator: s2s.ParallelSimulator, loss_function: s2s.ErrorHandler):
        super().__init__()
        self.model = model
        self.simulator = simulator
        # TODO: generalizza per computazione parallela
        self.simulator.add_worker(model)
        self.loss_function = loss_function
        
    @abstractmethod
    def minimize(self, args: Any) -> dict[str, dict[ParameterId, float] | float]:
        ...

attempt = 0

class NevergradOpt(Optimizer):

    def __init__(self, model: s2s.Simulator, simulator: s2s.ParallelSimulator, loss_function: s2s.ErrorHandler):
        super().__init__(model, simulator, loss_function)

    def minimize(self, args: dict[str,ng.p.Parameter | str | int]) -> Any:
        global attempt
        # args: {
        #   "params": ng.p !!mandatory!!
        #   "optimizer": str (optional) (name of nevergrad optimizer class, e.g. "OnePlusOne")
        #   "budget": int (optional)
        # }
        cfg = args or {}
        parametrization: ng.p.Parameter = cfg["params"]
        budget: int = cfg.get("budget", 3000)
        optimizer_name: str = cfg.get("optimizer", "OnePlusOne")

        optimizer_cls = getattr(ng.optimizers, optimizer_name, ng.optimizers.OnePlusOne)
        opt = optimizer_cls(parametrization=parametrization, budget=budget)

        def objective(params):
            attempt += 1
            self.model.set_parameters(params.value)
            start_time = time.time()
            # TODO: generalizza per esecuzione parallela
            (loss, err) = self.simulator.simulate(self.loss_function)[0]
            print(f"[INFO] attempt: {attempt}")
            print(f"[INFO] loss: {loss}")
            str_err = "WITH" if err else "WITHOUT"
            print(f"[INFO] termineted {str_err} error in {time.time() - start_time}")
            return loss

        recommendation = opt.minimize(objective)
        attempt = 0
        return recommendation.value