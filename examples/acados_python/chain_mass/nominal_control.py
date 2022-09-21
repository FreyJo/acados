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

# authors: Florian Messerer, Jonathan Frey
#
# This implementation of the nonlinear chain example follows the problem formulation:
# https://github.com/dkouzoup/hanging-chain-acado
# and the publication
# Recent Advances in Quadratic Programming Algorithmsfor Nonlinear Model Predictive Control
# https://cdn.syscop.de/publications/Kouzoupis2018.pdf
#
# The problem therein is extended by applying disturbances at each sampling time.
# These disturbances are the parameters of the model described in
# export_disturbed_chain_mass_model.py



import os
import numpy as np
import scipy.linalg

from acados_template import AcadosOcp, AcadosOcpSolver

from export_chain_mass_model import export_chain_mass_model
from export_chain_mass_integrator import export_chain_mass_integrator

from plot_utils import *
from utils import compute_steady_state, sampleFromEllipsoid, save_results_as_json
import matplotlib.pyplot as plt


def export_chain_mass_ocp_solver(chain_params):
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # chain parameters
    n_mass = chain_params["n_mass"]
    M = chain_params["n_mass"] - 2 # number of intermediate masses
    Ts = chain_params["Ts"]
    N = chain_params["N"]
    with_wall = chain_params["with_wall"]
    yPosWall = chain_params["yPosWall"]
    m = chain_params["m"]
    D = chain_params["D"]
    L = chain_params["L"]
    perturb_scale = chain_params["perturb_scale"]

    nlp_iter = chain_params["nlp_iter"]
    nlp_tol = chain_params["nlp_tol"]
    seed = chain_params["seed"]
    qp_solver = chain_params["qp_solver"]

    np.random.seed(seed)

    nparam = 3*M
    W = perturb_scale * np.eye(nparam)

    # export model
    model = export_chain_mass_model(n_mass, m, D, L)

    # set model
    ocp.model = model

    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = nx + nu
    Tf = N * Ts

    # initial state
    xPosFirstMass = np.zeros((3,1))
    xEndRef = np.zeros((3,1))
    xEndRef[0] = L * (M+1) * 6

    xrest = compute_steady_state(n_mass, m, D, L, xPosFirstMass, xEndRef)

    x0 = xrest

    # set dimensions
    ocp.dims.N = N

    # set cost module
    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    Q = 2*np.diagflat( np.ones((nx, 1)) )
    q_diag = np.ones((nx,1))
    strong_penalty = M+1
    q_diag[3*M] = strong_penalty
    q_diag[3*M+1] = strong_penalty
    q_diag[3*M+2] = strong_penalty
    Q = 2*np.diagflat( q_diag )

    R = 2*np.diagflat( 1e-2 * np.ones((nu, 1)) )

    ocp.cost.W = scipy.linalg.block_diag(Q, R)
    ocp.cost.W_e = Q

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[nx:nx+nu, :] = np.eye(nu)
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)

    yref = np.vstack((xrest, np.zeros((nu,1)))).flatten()
    ocp.cost.yref = yref
    ocp.cost.yref_e = xrest.flatten()

    # set constraints
    umax = 1*np.ones((nu,))

    ocp.constraints.constr_type = 'BGH'
    ocp.constraints.lbu = -umax
    ocp.constraints.ubu = umax
    ocp.constraints.x0 = x0.reshape((nx,))
    ocp.constraints.idxbu = np.array(range(nu))

    # wall constraint
    if with_wall:
        nbx = M + 1
        Jbx = np.zeros((nbx,nx))
        for i in range(nbx):
            Jbx[i, 3*i+1] = 1.0

        ocp.constraints.Jbx = Jbx
        ocp.constraints.lbx = yPosWall * np.ones((nbx,))
        ocp.constraints.ubx = 1e9 * np.ones((nbx,))

        # slacks
        ocp.constraints.Jsbx = np.eye(nbx)
        L2_pen = 1e3
        L1_pen = 1
        ocp.cost.Zl = L2_pen * np.ones((nbx,))
        ocp.cost.Zu = L2_pen * np.ones((nbx,))
        ocp.cost.zl = L1_pen * np.ones((nbx,))
        ocp.cost.zu = L1_pen * np.ones((nbx,))


    # solver options
    ocp.solver_options.qp_solver = qp_solver
    ocp.solver_options.qp_solver_iter_max = 1000
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'IRK'
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI
    ocp.solver_options.nlp_solver_max_iter = nlp_iter

    ocp.solver_options.sim_method_num_stages = 2
    ocp.solver_options.sim_method_num_steps = 2
    ocp.solver_options.qp_solver_cond_N = N//4
    ocp.solver_options.qp_tol = nlp_tol
    ocp.solver_options.tol = nlp_tol
    # ocp.solver_options.nlp_solver_tol_eq = 1e-9

    # set prediction horizon
    ocp.solver_options.tf = Tf

    acados_ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_ocp_' + model.name + '.json')

    return acados_ocp_solver


def run_nominal_control_closed_loop(chain_params):

    # chain parameters
    n_mass = chain_params["n_mass"]
    M = chain_params["n_mass"] - 2 # number of intermediate masses
    Ts = chain_params["Ts"]
    Tsim = chain_params["Tsim"]
    u_init = chain_params["u_init"]
    yPosWall = chain_params["yPosWall"]
    m = chain_params["m"]
    D = chain_params["D"]
    L = chain_params["L"]
    perturb_scale = chain_params["perturb_scale"]

    save_results = chain_params["save_results"]
    show_plots = chain_params["show_plots"]
    seed = chain_params["seed"]

    np.random.seed(seed)

    nparam = 3*M
    W = perturb_scale * np.eye(nparam)

    acados_ocp_solver = export_chain_mass_ocp_solver(chain_params)
    nx = acados_ocp_solver.acados_ocp.dims.nx
    nu = acados_ocp_solver.acados_ocp.dims.nu

    acados_integrator = export_chain_mass_integrator(n_mass, m, D, L)

    #%% get initial state by disturbing the rest position a bit
    x0 = acados_ocp_solver.acados_ocp.constraints.lbx_0
    xcurrent = x0.reshape((nx,))
    for i in range(5):
        acados_integrator.set("x", xcurrent)
        acados_integrator.set("u", u_init)

        status = acados_integrator.solve()
        if status != 0:
            raise Exception('acados integrator returned status {}. Exiting.'.format(status))

        # update state
        xcurrent = acados_integrator.get("x")

    #%% actual simulation
    N_sim = int(np.floor(Tsim/Ts))
    simX = np.ndarray((N_sim+1, nx))
    simU = np.ndarray((N_sim, nu))
    wall_dist = np.zeros((N_sim,))

    timings = np.zeros((N_sim,))
    timings_qp = np.zeros((N_sim,))
    timings_lin = np.zeros((N_sim,))
    sqp_iter = np.zeros((N_sim,))
    simX[0,:] = xcurrent

    # closed loop
    for i in range(N_sim):

        # solve ocp
        acados_ocp_solver.set(0, "lbx", xcurrent)
        acados_ocp_solver.set(0, "ubx", xcurrent)

        status = acados_ocp_solver.solve()
        timings[i] = acados_ocp_solver.get_stats("time_tot")[0]
        timings_qp[i] = acados_ocp_solver.get_stats("time_qp")[0]
        timings_lin[i] = acados_ocp_solver.get_stats("time_lin")[0]
        sqp_iter[i] = acados_ocp_solver.get_stats("sqp_iter")[0]

        if status != 0:
            acados_ocp_solver.print_statistics()
            raise Exception('acados acados_ocp_solver returned status {} in time step {}. Exiting.'.format(status, i))

        simU[i,:] = acados_ocp_solver.get(0, "u")
        print("control at time", i, ":", simU[i,:])

        # simulate system
        acados_integrator.set("x", xcurrent)
        acados_integrator.set("u", simU[i,:])

        pertubation = sampleFromEllipsoid(np.zeros((nparam,)), W)
        acados_integrator.set("p", pertubation)

        status = acados_integrator.solve()
        if status != 0:
            raise Exception('acados integrator returned status {}. Exiting.'.format(status))

        # update state
        xcurrent = acados_integrator.get("x")
        simX[i+1,:] = xcurrent

        # xOcpPredict = acados_ocp_solver.get(1, "x")
        # print("model mismatch = ", str(np.max(xOcpPredict - xcurrent)))
        yPos = xcurrent[range(1,3*M+1,3)]
        wall_dist[i] = np.min(yPos - yPosWall)
        print("time i = ", str(i), " dist2wall ", str(wall_dist[i]))

    print("dist2wall (minimum over simulation) ", str(np.min(wall_dist)))

    results = {"simX": simX, "simU": simU, "timings": timings, "timings_qp": timings_qp, "timings_lin": timings_lin, "sqp_iter": sqp_iter}

    #%% plot results
    if os.environ.get('ACADOS_ON_CI') is None and show_plots:
        plot_chain_control_traj(simU)
        plot_chain_position_traj(simX, yPosWall=yPosWall)
        plot_chain_velocity_traj(simX)

        # animate_chain_position(simX, xPosFirstMass, yPosWall=yPosWall)
        # animate_chain_position_3D(simX, xPosFirstMass)

        plt.show()

    if save_results:
        save_results_as_json(results, chain_params, id='closed_loop')

    return


def run_nominal_control_open_loop(chain_params):

    # chain parameters
    n_mass = chain_params["n_mass"]
    M = chain_params["n_mass"] - 2 # number of intermediate masses
    u_init = chain_params["u_init"]
    yPosWall = chain_params["yPosWall"]
    m = chain_params["m"]
    D = chain_params["D"]
    L = chain_params["L"]
    perturb_scale = chain_params["perturb_scale"]

    N = chain_params["N"]
    N_run = chain_params["N_run"]
    save_results = chain_params["save_results"]
    show_plots = chain_params["show_plots"]
    seed = chain_params["seed"]

    np.random.seed(seed)

    nparam = 3*M
    W = perturb_scale * np.eye(nparam)

    acados_ocp_solver = export_chain_mass_ocp_solver(chain_params)
    nx = acados_ocp_solver.acados_ocp.dims.nx
    nu = acados_ocp_solver.acados_ocp.dims.nu

    acados_integrator = export_chain_mass_integrator(n_mass, m, D, L)

    #%% get initial state by disturbing the rest position a bit
    xrest = acados_ocp_solver.acados_ocp.constraints.lbx_0
    xcurrent = xrest.reshape((nx,))
    for i in range(5):
        acados_integrator.set("x", xcurrent)
        acados_integrator.set("u", u_init)

        status = acados_integrator.solve()
        if status != 0:
            raise Exception('acados integrator returned status {}. Exiting.'.format(status))

        # update state
        xcurrent = acados_integrator.get("x")

    #%% solve actual problem
    timings = np.zeros((N_run,))
    timings_qp = np.zeros((N_run,))
    timings_lin = np.zeros((N_run,))
    sqp_iter = np.zeros((N_run,))

    simX = np.ndarray((N+1, nx))
    simU = np.ndarray((N, nu))

    # initialize
    acados_ocp_solver.set(0, "x", xcurrent)
    acados_ocp_solver.store_iterate('default_init.json', overwrite=True)

    acados_ocp_solver.set(0, "lbx", xcurrent)
    acados_ocp_solver.set(0, "ubx", xcurrent)
    # experiment loop
    for i in range(N_run):

        acados_ocp_solver.load_iterate('default_init.json')

        # solve ocp

        status = acados_ocp_solver.solve()
        timings[i] = acados_ocp_solver.get_stats("time_tot")[0]
        timings_qp[i] = acados_ocp_solver.get_stats("time_qp")[0]
        timings_lin[i] = acados_ocp_solver.get_stats("time_lin")[0]
        sqp_iter[i] = acados_ocp_solver.get_stats("sqp_iter")[0]

        if status != 0:
            acados_ocp_solver.print_statistics()
            raise Exception('acados acados_ocp_solver returned status {} in time step {}. Exiting.'.format(status, i))

        if i == 0:
            for i in range(N+1):
                simX[i,:] = acados_ocp_solver.get(i, "x")
            for i in range(N):
                simU[i,:] = acados_ocp_solver.get(i, "u")

    results = {"simX": simX, "simU": simU, "timings": timings, "timings_qp": timings_qp, "timings_lin": timings_lin, "sqp_iter": sqp_iter}

    #%% plot results
    if os.environ.get('ACADOS_ON_TRAVIS') is None and show_plots:
        plot_chain_control_traj(simU)
        plot_chain_position_traj(simX, yPosWall=yPosWall)
        plot_chain_velocity_traj(simX)

        # animate_chain_position(simX, xPosFirstMass, yPosWall=yPosWall)
        # animate_chain_position_3D(simX, xPosFirstMass)

        plt.show()

    if save_results:
        save_results_as_json(results, chain_params, id='open_loop')

    return