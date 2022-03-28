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

import sys
import os
sys.path.insert(0, os.path.join('..','getting_started','common'))

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver
from pendulum_model import export_pendulum_ode_model
from utils import plot_pendulum
import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt

# create ocp object to formulate the OCP
ocp = AcadosOcp()

# set model
model = export_pendulum_ode_model()
ocp.model = model
print(f"model {model.__dict__}")

Tf = 1.0
nx = model.x.size()[0]
nu = model.u.size()[0]
ny = nx + nu
ny_e = nx
N = 20

# set dimensions
ocp.dims.N = N
# NOTE: all dimensions but N are now detected automatically in the Python
#  interface, all other dimensions will be overwritten by the detection.

# set cost module
ocp.cost.cost_type = 'LINEAR_LS'
ocp.cost.cost_type_e = 'LINEAR_LS'

Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
R = 2*np.diag([1e-0])

ocp.cost.W = scipy.linalg.block_diag(Q, R)

ocp.cost.W_e = Q

ocp.cost.Vx = np.zeros((ny, nx))
ocp.cost.Vx[:nx,:nx] = np.eye(nx)

Vu = np.zeros((ny, nu))
Vu[4,0] = 1.0
ocp.cost.Vu = Vu

ocp.cost.Vx_e = np.eye(nx)

ocp.cost.yref  = np.zeros((ny, ))
ocp.cost.yref_e = np.zeros((ny_e, ))

# set constraints
Fmax = 10
x0 = np.array([0.0, np.pi/8, 0.0, 0.0])
ocp.constraints.constr_type = 'BGH'
ocp.constraints.x0 = x0

ocp.constraints.lbu = np.array([-Fmax])
ocp.constraints.ubu = np.array([+Fmax])
ocp.constraints.idxbu = np.array([0])

ocp.solver_options.levenberg_marquardt = 1e0

xmax = 1.1
ocp.constraints.lbx = np.array([-xmax])
ocp.constraints.ubx = np.array([+xmax])
ocp.constraints.idxbx = np.array([0])

ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
ocp.solver_options.integrator_type = 'ERK'
ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI

ocp.solver_options.qp_solver_cond_N = N

# set prediction horizon
ocp.solver_options.tf = Tf

acados_ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_ocp_' + model.name + '.json')
acados_integrator = AcadosSimSolver(ocp, json_file = 'acados_ocp_' + model.name + '.json')

Nsim = 50

relaxed_iterate = 'relaxed_sol.json'
time_tot_relaxed = np.zeros((Nsim,))
time_tot_exact = np.zeros((Nsim,))

tol_relax = 1e-3
tol_exact = 1e-6
tau_tol_factor = 1e-2

Nplot = 2
for irun, RELAXED in enumerate([True, False]):
    qp_iter_relaxed = np.zeros((Nsim,))
    qp_iter_exact = np.zeros((Nsim,))

    acados_ocp_solver.reset()
    simX = np.ndarray((Nsim+1, nx))
    simU = np.ndarray((Nsim, nu))

    xcurrent = x0
    simX[0,:] = xcurrent

    # closed loop
    for i in range(Nsim):

        # set new ocp constraints
        acados_ocp_solver.set(0, "lbx", xcurrent)
        acados_ocp_solver.set(0, "ubx", xcurrent)

        if RELAXED:
            compute_new_relaxed = True
            if i > 0:
                # check if relaxed solution has to be recomputed
                acados_ocp_solver.load_iterate(filename=relaxed_iterate)
                residuals = acados_ocp_solver.get_residuals(recompute=True)
                if residuals[0] < tol_relax and all(residuals[2:4] < tol_relax):
                    compute_new_relaxed = False

            if compute_new_relaxed:
                # solve relaxed
                acados_ocp_solver.options_set('tol_comp', tol_relax)
                acados_ocp_solver.options_set('tol_stat', tol_relax)
                acados_ocp_solver.options_set('tol_ineq', tol_relax)
                acados_ocp_solver.options_set('qp_tol_ineq', tol_relax)
                # acados_ocp_solver.options_set('tol_eq', tol_relax)
                # acados_ocp_solver.options_set('qp_warm_start', 2)
                acados_ocp_solver.options_set('qp_tau_min', tol_relax*tau_tol_factor) # minimum value of barrier parameter in HPIPM

                # acados_ocp_solver.options_set('qp_mu0', 1e-2) # initial value for complementarity slackness

                if i > 0:
                    acados_ocp_solver.load_iterate(filename=relaxed_iterate)
                status = acados_ocp_solver.solve()

                print("stats relaxed OCP")
                acados_ocp_solver.print_statistics()
                if status != 0:
                    raise Exception('acados acados_ocp_solver returned status {}. Exiting.'.format(status))
                acados_ocp_solver.store_iterate(filename=relaxed_iterate, overwrite=True)
                qp_iter_relaxed[i] = sum(acados_ocp_solver.get_stats('qp_iter'))
                time_tot_relaxed[i] = sum(acados_ocp_solver.get_stats('time_tot'))

        # solve exact
        acados_ocp_solver.options_set('tol_comp', tol_exact)
        acados_ocp_solver.options_set('tol_stat', tol_exact)
        acados_ocp_solver.options_set('tol_ineq', tol_exact)
        acados_ocp_solver.options_set('qp_warm_start', 0)
        acados_ocp_solver.options_set('qp_tol_ineq', tol_exact)
        acados_ocp_solver.options_set('qp_tau_min', tol_exact*tau_tol_factor)

        # acados_ocp_solver.options_set('tol_eq', 1e-6)
        status = acados_ocp_solver.solve()
        acados_ocp_solver.print_statistics()
        if status != 0:
            raise Exception('acados acados_ocp_solver returned status {}. Exiting.'.format(status))
        qp_iter_exact[i] = sum(acados_ocp_solver.get_stats('qp_iter'))
        time_tot_exact[i] = sum(acados_ocp_solver.get_stats('time_tot'))

        simU[i,:] = acados_ocp_solver.get(0, "u")
        print("stats exact OCP")

        # simulate system
        acados_integrator.set("x", xcurrent)
        acados_integrator.set("u", simU[i,:])

        status = acados_integrator.solve()
        if status != 0:
            raise Exception('acados integrator returned status {}. Exiting.'.format(status))

        # update state
        xcurrent = acados_integrator.get("x")
        simX[i+1,:] = xcurrent



    print(f'qp_iter relaxed solver: {sum(qp_iter_relaxed)}, exact: {sum(qp_iter_exact)}, total: {sum(qp_iter_relaxed) + sum(qp_iter_exact)}')
    print(f'time_tot relaxed solver: {sum(time_tot_relaxed):.3f}, exact: {sum(time_tot_exact):.3f}')

    print(f'sum qp_iter: {qp_iter_relaxed + qp_iter_exact}')
    print(f'qp_iter relaxed: {qp_iter_relaxed}')
    print(f'qp_iter exact: {qp_iter_exact}')


    plt.subplot(Nplot, 1, irun+1)
    plt.bar(range(Nsim), qp_iter_relaxed)
    plt.bar(range(Nsim), qp_iter_exact, bottom=qp_iter_relaxed)
    plt.ylabel('qp iter')
    plt.title(f"simulation with relaxed solution {RELAXED}")
plt.show()

import pdb; pdb.set_trace()

# plot results
plot_pendulum(np.linspace(0, Tf/N*Nsim, Nsim+1), Fmax, simU, simX)