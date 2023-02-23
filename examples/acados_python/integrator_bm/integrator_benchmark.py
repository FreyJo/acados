# -*- coding: future_fstrings -*-
#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#

from acados_template import AcadosSim, AcadosSimSolver, make_object_json_dumpable
from pendulum_model import export_pendulum_ode_model
from utils import plot_pendulum, latexify_plot
import numpy as np
import json

JAC_REUSE = 0
STAGE_RANGE = range(1, 9)

def run_experiment(num_stages: int = 3, integrator_type='IRK') -> np.ndarray:
    sim = AcadosSim()

    # export model
    model = export_pendulum_ode_model()

    # set model_name
    sim.model = model

    Tf = 0.1
    nx = model.x.size()[0]
    nu = model.u.size()[0]
    N = 20000

    # set simulation time
    sim.solver_options.T = Tf
    # set options
    sim.solver_options.integrator_type = integrator_type
    sim.solver_options.num_stages = num_stages
    sim.solver_options.num_steps = 3
    sim.solver_options.collocation_type = "GAUSS_RADAU_IIA"
    sim.solver_options.sim_method_jac_reuse = JAC_REUSE
    sim.solver_options.newton_iter = 10 # for implicit integrator
    sim.solver_options.newton_tol = 1e-10

    # create
    acados_integrator = AcadosSimSolver(sim)

    simX = np.zeros((N+1, nx))
    x0 = np.array([5.863919e-02, 3.176782e+00, -2.557677e-01, 3.518985e+00])
    u0 = np.array([0.0])
    acados_integrator.set("u", u0)

    # simX[0,:] = x0
    cpu = np.zeros((N,))

    for i in range(N):
        # set initial state
        acados_integrator.set("x", x0)
        # initialize IRK
        if sim.solver_options.integrator_type == 'IRK':
            acados_integrator.set("xdot", np.zeros((nx,)))

        # solve
        status = acados_integrator.solve()
        # get solution
        # simX[i+1,:] = acados_integrator.get("x")
        cpu[i] = acados_integrator.get("time_tot")

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    S_forw = acados_integrator.get("S_forw")
    print("S_forw, sensitivities of simulation result wrt x,u:\n", S_forw)
    del acados_integrator

    return cpu

    # # plot results
    # plot_pendulum(np.linspace(0, N*Tf, N+1), 10, np.repeat(u0, N), simX)

def get_run_id(num_stages: int) -> str:
    return f'num_stages_{num_stages}'

def get_results_filename(integrator_type: str, jac_reuse) -> str:
    return f'integrator_benchmark_{integrator_type.lower()}_jac_reuse_{jac_reuse}.json'


def main_benchmark(integrator_type):
    results = {}
    for ns in STAGE_RANGE:
        cpu = run_experiment(num_stages=ns, integrator_type=integrator_type)
        run_id = get_run_id(ns)
        results[run_id] = cpu
    filename = get_results_filename(integrator_type, jac_reuse=JAC_REUSE)
    with open(filename, 'w') as f:
        json.dump(results, f, default=make_object_json_dumpable)

def main_eval(integrator_type):
    filename = get_results_filename(integrator_type, jac_reuse=JAC_REUSE)
    with open(filename, 'r') as f:
        results = json.load(f)
    for ns in STAGE_RANGE:
        run_id = get_run_id(ns)
        print(f'num_stages: {ns}, mean: {np.mean(results[run_id])*1e3:.3f} ms, min: {np.min(results[run_id])*1e3:.3f} ms, max: {np.max(results[run_id])*1e3:.3f} ms')

import matplotlib.pyplot as plt
INTEGRATOR_TYPES = ['IRK', 'GNSF']
def main_eval_plot():
    results = {}
    latexify_plot()
    for integrator_type in INTEGRATOR_TYPES:
        filename = get_results_filename(integrator_type, jac_reuse=JAC_REUSE)
        with open(filename, 'r') as f:
            results[integrator_type] = json.load(f)
        means = []
        for ns in STAGE_RANGE:
            run_id = get_run_id(ns)
            means = np.append(means, np.mean(results[integrator_type][run_id])*1e3)
        plt.plot(STAGE_RANGE, means, label=integrator_type)
    plt.xlabel('num stages')
    plt.ylabel('mean runtime [ms]')
    plt.title(f'with jac_reuse ={JAC_REUSE}')
    plt.grid()
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main_benchmark(integrator_type='IRK')
    main_benchmark(integrator_type='GNSF')
    main_eval_plot()
