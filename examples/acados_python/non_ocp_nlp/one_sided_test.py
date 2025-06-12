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

from acados_template import AcadosOcp, AcadosOcpSolver, ACADOS_INFTY, AcadosOcpFlattenedIterate
import numpy as np
import casadi as ca


def create_solver(solver_name: str, nlp_solver_type: str = 'SQP_WITH_FEASIBLE_QP',
                  variant: int = 0):

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    nx = 2
    # set model
    ocp.model.name = f"qp_{solver_name}"
    ocp.model.x = ca.SX.sym('x', nx, 1)

    ny = nx

    # discretization
    N = 0

    ocp.cost.W_e = 2*np.diag([1e1, 1e1])

    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    ocp.cost.Vx_e = np.eye((nx))
    ocp.cost.yref_e = np.ones((ny, ))

    # set constraints
    xmax = 2.0
    ocp.constraints.lbx_e = -xmax * np.ones((nx,))
    ocp.constraints.ubx_e = +xmax * np.ones((nx,))
    ocp.constraints.idxbx_e = np.arange(nx)

    # define soft nonlinear constraint
    scale_h = 1.0
    radius = 1.0
    ocp.model.con_h_expr_e = scale_h * (ocp.model.x[0]**2 + ocp.model.x[1]**2)
    ocp.constraints.lh_e = -1000 * np.ones((1,))
    ocp.constraints.lh_e = -ACADOS_INFTY * np.ones((1,))
    ocp.constraints.uh_e = scale_h * radius**2 * np.ones((1,))

    # soften
    ocp.constraints.idxsh_e = np.array([0])
    ocp.cost.zl_e = np.array([0.0])
    ocp.cost.zu_e = np.array([1e4])
    ocp.cost.Zl_e = np.array([0.0])
    ocp.cost.Zu_e = np.array([1e2])
    if variant == 0:
        # works as expected
        ocp.constraints.lsh_e = -ACADOS_INFTY * np.ones((1,))
    elif variant == 1:
        # does not work as expected -> HPIPM issue?
        # slack is not determined as cost is 0
        # res comp does not converge
        ocp.constraints.lsh_e = -0 * np.ones((1,))
    elif variant == 2:
        # works as expected
        ocp.constraints.lsh_e = -0 * np.ones((1,))
        ocp.cost.Zl_e = 1e2 * np.ones((1,))
    elif variant == 3:
        # doesnt work: res_stat in NLP solver does not converge
        # -> ignore contributions of masked slacked constraints?
        ocp.constraints.lsh_e = -ACADOS_INFTY * np.ones((1,))
        ocp.cost.Zl_e = 1e2 * np.ones((1,))
    elif variant == 4:
        # doesnt work: res_stat in NLP solver does not converge
        # -> ignore contributions of masked slacked constraints?
        ocp.cost.Zl_e = 1e2 * np.ones((1,))
        ocp.constraints.lsh_e = -ACADOS_INFTY * np.ones((1,))

    # set options
    solver_options = ocp.solver_options
    solver_options.N_horizon = N

    solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    qp_tol = 1e-7
    solver_options.qp_tol = qp_tol
    solver_options.qp_solver_ric_alg = 1
    solver_options.qp_solver_mu0 = 1e4
    solver_options.qp_solver_warm_start = 1
    solver_options.qp_solver_iter_max = 400
    solver_options.hessian_approx = 'GAUSS_NEWTON'
    solver_options.nlp_solver_type = nlp_solver_type
    solver_options.globalization = 'FUNNEL_L1PEN_LINESEARCH'
    solver_options.nlp_solver_ext_qp_res = 1

    # create ocp solver
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)

    return ocp, ocp_solver


def call_solver(ocp_solver: AcadosOcpSolver) -> AcadosOcpFlattenedIterate:
    # solve
    status = ocp_solver.solve()
    ocp_solver.print_statistics()

    sqp_iter = ocp_solver.get_stats('sqp_iter')
    if status != 0:
        # raise RuntimeError(f"acados returned status {status} after {sqp_iter} SQP iterations.")
        print(f'acados returned status {status}.')

    print(f"cost function value = {ocp_solver.get_cost()} after {sqp_iter} SQP iterations")
    sol = ocp_solver.store_iterate_to_flat_obj()
    return sol


def main():
    print("Reference ...")
    ocp, ocp_solver = create_solver("2", nlp_solver_type="SQP", variant=3)
    sol = call_solver(ocp_solver)
    print(f"Reference solution: {sol}")



if __name__ == '__main__':
    main()

