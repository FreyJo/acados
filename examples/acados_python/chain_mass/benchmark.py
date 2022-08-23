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

import numpy as np

from nominal_control import run_nominal_control_closed_loop, run_nominal_control_open_loop
from utils import get_chain_params, load_results_from_json

import matplotlib.pyplot as plt

from plot_utils import latexify


N_MASSES = list(range(3,9))
# N_MASSES = [8]
# QP_SOLVERS = ["FULL_CONDENSING_DAQP", "FULL_CONDENSING_HPIPM"]
QP_SOLVERS = ["FULL_CONDENSING_HPIPM", "FULL_CONDENSING_DAQP", "FULL_CONDENSING_QPOASES", "PARTIAL_CONDENSING_HPIPM"]
# QP_SOLVERS = ["PARTIAL_CONDENSING_HPIPM"]
ID = 'closed_loop'
ID = 'open_loop'

# N_ocp_list = [8, 9, 11, 12]
N_ocp_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100]
# N_ocp_list = [10, 15, 20, 25, 30]
# N_ocp_list = [40, 50, 60, 70, 80, 90, 100]

def run_benchmark():
    chain_params = get_chain_params()
    for qp_solver in QP_SOLVERS:
        for n_mass in N_MASSES:
            # adjust parameters wrt experiment
            chain_params["n_mass"] = n_mass
            chain_params["qp_solver"] = qp_solver
            print(chain_params)
            if ID == 'closed_loop':
                run_nominal_control_closed_loop(chain_params)
            elif ID == 'open_loop':
                run_nominal_control_open_loop(chain_params)

def eval_benchmark():
    plt.figure()
    latexify()
    chain_params = get_chain_params()

    for qp_solver in QP_SOLVERS:
        timings_qp = np.zeros(len(N_MASSES))
        sqp_iter = np.zeros(len(N_MASSES))
        for i, n_mass in enumerate(N_MASSES):
            # adjust parameters wrt experiment
            chain_params["n_mass"] = n_mass
            chain_params["qp_solver"] = qp_solver

            results = load_results_from_json(chain_params, id='closed_loop')
            timings_qp[i] = np.mean(np.array(results['timings_qp']))
            # timings_qp[i] = np.mean(np.array(results['timings_qp']) / np.array(results['sqp_iter']))
            # sqp_iter[i] = results['sqp_iter']
            # print(f"sqp_iter {results['sqp_iter']}")
        label = qp_solver.replace('_', ' ')
        # label.replace('PARTIAL_CONDENSING', 'PARTIAL_CONDENSING N' + qp):

        plt.plot(N_MASSES, timings_qp, label=label)
    plt.legend()
    plt.xlabel('number of masses')
    plt.ylabel('mean CPU time / QP')
    plt.yscale('log')
    plt.xticks(N_MASSES)
    plt.grid()
    plt.show()

def run_benchmark_N():
    chain_params = get_chain_params()
    for qp_solver in QP_SOLVERS:
        for N_ocp in N_ocp_list:
            chain_params['N'] = N_ocp
            chain_params["qp_solver"] = qp_solver
            run_nominal_control_open_loop(chain_params)

def eval_benchmark_N():
    plt.figure()
    latexify()
    chain_params = get_chain_params()

    for qp_solver in QP_SOLVERS:
        timings_qp = np.zeros(len(N_ocp_list))
        sqp_iter = np.zeros(len(N_ocp_list))
        for i, N_ocp in enumerate(N_ocp_list):
            # adjust parameters wrt experiment
            chain_params["N"] = N_ocp
            chain_params["qp_solver"] = qp_solver

            results = load_results_from_json(chain_params, id='open_loop')
            timings_qp[i] = np.mean(np.array(results['timings_qp']) / np.array(results['sqp_iter']))
            # sqp_iter[i] = results['sqp_iter']
            # print(f"sqp_iter {results['sqp_iter']}")
        label = qp_solver.replace('_', ' ')
        plt.plot(N_ocp_list, timings_qp, label=label)
    plt.legend()
    plt.xlabel('N')
    plt.ylabel('mean CPU time / QP')
    plt.yscale('log')
    plt.xticks(N_ocp_list)
    plt.grid()
    plt.show()

if __name__ == "__main__":
    # run_benchmark()
    eval_benchmark()
    # run_benchmark_N()
    eval_benchmark_N()