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

# TODO:
# - test on single problem
#   - without slacks to use more QP solvers
#   - vary horizon length!

# - shift in closed loop
# - vary tolerances?

import numpy as np

from nominal_control import run_nominal_control_closed_loop
from utils import get_chain_params, load_results_from_json

import matplotlib.pyplot as plt

from plot_utils import latexify

chain_params = get_chain_params()
chain_params["save_results"] = True
chain_params["show_plots"] = False

N_MASSES = list(range(3,9))
# N_MASSES = [8]
# QP_SOLVERS = ["FULL_CONDENSING_DAQP", "FULL_CONDENSING_HPIPM"]
QP_SOLVERS = ["FULL_CONDENSING_HPIPM", "FULL_CONDENSING_DAQP", "FULL_CONDENSING_QPOASES", "PARTIAL_CONDENSING_HPIPM"]
# QP_SOLVERS = ["PARTIAL_CONDENSING_HPIPM"]

def run_benchmark():
    for qp_solver in QP_SOLVERS:
        for n_mass in N_MASSES:
            # adjust parameters wrt experiment
            chain_params["n_mass"] = n_mass
            chain_params["qp_solver"] = qp_solver
            print(chain_params)

            run_nominal_control_closed_loop(chain_params)



def eval_benchmark():
    plt.figure()
    latexify()
    id = 'closed_loop'
    for qp_solver in QP_SOLVERS:
        timings_qp = np.zeros(len(N_MASSES))
        sqp_iter = np.zeros(len(N_MASSES))
        for i, n_mass in enumerate(N_MASSES):
            # adjust parameters wrt experiment
            chain_params["n_mass"] = n_mass
            chain_params["qp_solver"] = qp_solver

            results = load_results_from_json(chain_params, id=id)
            timings_qp[i] = np.mean(np.array(results['timings_qp']) / np.array(results['sqp_iter']))
            # sqp_iter[i] = results['sqp_iter']
            # print(f"sqp_iter {results['sqp_iter']}")

        plt.plot(N_MASSES, timings_qp, label=qp_solver)
    plt.legend()
    plt.xlabel('number of masses')
    plt.ylabel('mean CPU time / QP')
    plt.yscale('log')
    plt.xticks(N_MASSES)
    plt.grid()
    plt.show()

if __name__ == "__main__":
    # run_benchmark()
    eval_benchmark()