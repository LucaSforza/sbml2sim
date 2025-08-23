from abc import ABC, abstractmethod
import math
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
    def minimize(self, args: Any) -> tuple[dict[str, dict[ParameterId, float]], float]:
        ...

attempt = 0
best_value = math.inf

class NevergradOpt(Optimizer):

    def __init__(self, model: s2s.Simulator, simulator: s2s.ParallelSimulator, loss_function: s2s.ErrorHandler):
        super().__init__(model, simulator, loss_function)

    def minimize(self, args: dict[str,ng.p.Parameter | str | int]) -> tuple[dict[str, dict[ParameterId, float]], float]:
        
        # args: {
        #   "params": ng.p !!mandatory!!
        #   "optimizer": str (optional) (name of nevergrad optimizer class, e.g. "OnePlusOne")
        #   "budget": int (optional)
        # }
        cfg = args or {}
        parametrization: ng.p.Parameter = cfg["params"]
        budget: int = cfg.get("budget", 3000)
        optimizer_name: str = cfg.get("optimizer", "NGOpt")

        def objective(params):
            global attempt, best_value
            attempt += 1
            # print(f"[INFO] parameters: {params}")
            self.model.set_parameters(params)
            start_time = time.time()
            # TODO: generalizza per esecuzione parallela
            (loss, err) = self.simulator.simulate(self.loss_function)[0]
            print(f"[INFO] attempt: {attempt}")
            print(f"[INFO] loss: {loss}")
            str_err = "WITH" if err else "WITHOUT"
            if loss < best_value:
                best_value = loss
            print(f"[INFO] termineted {str_err} error in {time.time() - start_time}")
            if err:
                return math.inf
            else:
                return loss

        optimizer_cls = getattr(ng.optimizers, optimizer_name, ng.optimizers.OnePlusOne)
        if optimizer_name == "BayesianOptimization":
            opt = ng.optimizers._BO(parametrization=parametrization, budget=budget,init_budget=20, utility_kind="ei")
        else:
            opt = optimizer_cls(parametrization=parametrization, budget=budget, num_workers=1)

        recommendation = opt.minimize(objective)
        return recommendation.value, best_value