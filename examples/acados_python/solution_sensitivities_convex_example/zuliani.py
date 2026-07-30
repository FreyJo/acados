#
# Copyright (c) The acados authors.
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

import numpy as np
import casadi as ca
from acados_template import AcadosOcp, AcadosMultiphaseOcp, AcadosModel, AcadosOcpSolver, latexify_plot
import matplotlib.pyplot as plt
latexify_plot()


def export_parametric_nlp() -> AcadosMultiphaseOcp:

    x1 = ca.SX.sym("x1")
    x2 = ca.SX.sym("x2")
    x3 = ca.SX.sym("x3")
    theta = ca.SX.sym("p_global", 1)

    ocp1 = AcadosOcp()

    ocp1.model.x = ca.vertcat(x1, x2)
    ocp1.model.cost_expr_ext_cost_0 = theta / 2 * (x1**2) + x2
    ocp1.model.disc_dyn_expr = - (x1 + x2)
    ocp1.cost.cost_type_0 = "EXTERNAL"
    ocp1.model.name = "zuliani_phase_1"

    ocp2 = AcadosOcp()
    ocp2.model.name = "zuliani_phase_2"
    ocp2.model.x = x3
    ocp2.model.cost_expr_ext_cost_e = x3
    ocp2.cost.cost_type_e = "EXTERNAL"
    ocp2.model.disc_dyn_expr = x3 # actually not used


    ocp1.model.p_global = theta
    ocp2.model.p_global = theta

    mocp = AcadosMultiphaseOcp(N_list=[1, 0])

    mocp.set_phase(ocp1, 0)
    mocp.set_phase(ocp2, 1)

    mocp.solver_options.qp_solver = "FULL_CONDENSING_HPIPM"
    mocp.solver_options.hessian_approx = "EXACT"
    mocp.solver_options.N_horizon = 1
    mocp.solver_options.tf = 1.0
    mocp.mocp_opts.integrator_type = ["DISCRETE", "DISCRETE"]

    mocp.p_global_values = np.zeros((1,))
    # mocp.solver_options.with_solution_sens_wrt_params = True
    # mocp.solver_options.with_value_sens_wrt_params = True
    mocp.solver_options.nlp_solver_ext_qp_res = 1

    return mocp

def solve_and_compute_sens(p_test, tau):
    np_test = p_test.shape[0]

    ocp = export_parametric_nlp()
    ocp.solver_options.tau_min = tau
    ocp.solver_options.qp_solver_t0_init = 0
    ocp.solver_options.nlp_solver_ext_qp_res = 1
    ocp.solver_options.nlp_solver_max_iter = 2 # QP should converge in one iteration

    ocp_solver = AcadosOcpSolver(ocp, json_file="parameter_augmented_acados_ocp.json", verbose=False)

    sens_x = np.zeros(np_test)
    nx = 3
    solution = np.zeros((nx, np_test))

    ocp_solver.set(0, 'x', np.array([2.0, 1.0]))
    ocp_solver.set(1, 'x', np.array([2.0]))

    for i, p in enumerate(p_test):
        p_val = np.array([p])

        ocp_solver.set_p_global_and_precompute_dependencies(p_val)
        status = ocp_solver.solve()
        sol = np.hstack((ocp_solver.get(0, "x"), ocp_solver.get(1, "x")))
        print(f"p={p}, solution={sol}")
        solution[:, i] = sol

        # ocp_solver.print_statistics()
        # qp = ocp_solver.get_last_qp()
        # print(f"qp = {qp}")

        if status != 0:
            ocp_solver.print_statistics()
            raise Exception(f"OCP solver returned status {status} at {i}th p value {p}, {tau=}.")
            # print(f"OCP solver returned status {status} at {i}th p value {p}, {tau=}.")
            # breakpoint()
        # compare with analytic solution
        if np.abs(sol[0] - (1 / p)) > 1e-6:
            raise Exception(f"solution does not match analytic solution at {i}th p value {p}, {tau=}.")
        if np.abs((sol[1] + sol[2]) + (1 / p)) > 1e-6:
            raise Exception(f"solution does not match analytic solution at {i}th p value {p}, {tau=}.")

        sens_x[i] = 0.0  # placeholder for sensitivity value
        # status = ocp_solver.setup_qp_matrices_and_factorize()
        # if status != 0:
        #     ocp_solver.print_statistics()
        #     raise Exception(f"OCP solver returned status {status} in setup_qp_matrices_and_factorize at {i}th p value {p}, {tau=}.")

        # # Calculate the policy gradient
        # out_dict = ocp_solver.eval_solution_sensitivity(0, "p_global", return_sens_x=True, return_sens_u=False)
        # sens_x[i] = out_dict['sens_x'].item()

    return solution, sens_x

def main():
    delta_p = 0.002
    p_test = np.arange(0.9, 1.1, delta_p)
    sens_list = []
    labels_list = []
    sol_list = []
    tau = 1e-6
    solution, sens_x = solve_and_compute_sens(p_test, tau)

    # Compare to numerical gradients
    sens_x_fd = np.gradient(solution, delta_p)
    test_tol = 1e-2
    median_diff = np.median(np.abs(sens_x - sens_x_fd))
    # print(f"Median difference between policy gradient obtained by acados and via FD is {median_diff} should be < {test_tol}.")
    # # test: check median since derivative cannot be compared at active set changes
    # assert median_diff <= test_tol

    sens_list.append(sens_x)
    labels_list.append(r"$\tau = 10^{-6}$")
    sol_list.append(solution)

    tau_vals = [1e-4, 1e-3, 1e-2]
    tau_vals = [1e-4]
    for tau in tau_vals:
        sol_tau, sens_x_tau = solve_and_compute_sens(p_test, tau)
        sens_list.append(sens_x_tau)
        labels_list.append(r"$\tau = 10^{" + f"{int(np.log10(tau))}" + r"}$")
        # labels_list.append(r"$\tau =" + f"{tau}" + r"$")
        sol_list.append(sol_tau)

    plot_solution_sensitivities_results(p_test, sol_list, sens_list, labels_list,
                 title=None, parameter_name=r"$\theta$")
    plot_solution_sensitivities_results(p_test, sol_list, sens_list, labels_list,
                 title=None, parameter_name=r"$\theta$", horizontal_plot=True)

def plot_solution_sensitivities_results(p_test, sol_list, sens_list, labels_list, title=None, parameter_name="", fig_filename=None, horizontal_plot=False):
    p_min = p_test[0]
    p_max = p_test[-1]
    linestyles = ["--", "-.", "--", ":", "-.", ":"]

    nsub = 2
    if horizontal_plot:
        _, ax = plt.subplots(nrows=1, ncols=nsub, sharex=False, figsize=(12, 3.0))
    else:
        _, ax = plt.subplots(nrows=nsub, ncols=1, sharex=True, figsize=(6.5,5))

    isub = 0
    # plot analytic solution
    for i, sol in enumerate(sol_list):
        ax[isub].plot(p_test, sol, label=labels_list[i], linestyle=linestyles[i])
    ax[isub].set_xlim([p_test[0], p_test[-1]])
    ax[isub].set_ylabel(r"solution $x^{\star}$")
    if title is not None:
        ax[isub].set_title(title)
    ax[isub].legend()

    isub += 1

    # plot numerical sensitivities
    for i, sens_x_tau in enumerate(sens_list):
        ax[isub].plot(p_test, sens_x_tau, label=labels_list[i], color=f"C{i}", linestyle=linestyles[i])
    ax[isub].set_xlim([p_test[0], p_test[-1]])
    ax[isub].set_ylabel(r"derivative $\partial_\theta x^{\star}$")
    # ax[isub].legend(ncol=2)

    for i in range(nsub):
        ax[i].grid(True)
        if horizontal_plot:
            ax[i].set_xlabel(f"{parameter_name}")
    ax[-1].set_xlabel(f"{parameter_name}")

    plt.tight_layout()

    if fig_filename is not None:
        plt.savefig(fig_filename)
        print(f"stored figure as {fig_filename}")
    plt.show()


if __name__ == "__main__":
    main()

    # to plot only analytic solution
    # plot_solution_sensitivities_results([-2, 2], [], [], [], parameter_name=r"$\theta$", fig_filename="solution_sens_non_ocp_analytic.pdf")
