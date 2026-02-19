import json
import sys
import os
sys.path.insert(0, '../getting_started')
import numpy as np
import casadi as ca
from typing import Union

from acados_template import AcadosOcp, AcadosOcpSolver
from pendulum_model import export_pendulum_ode_model

def formulate_ocp(Tf: float = 1.0, N: int = 20)-> AcadosOcp:
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    # add dummy control
    model.u = ca.vertcat(model.u, ca.SX.sym('dummy_u'))
    ocp.model = model

    # set h constraints
    ocp.model.con_h_expr_0 = ca.norm_2(model.x)
    ocp.constraints.lh_0 = np.array([0])
    ocp.constraints.uh_0 = np.array([3.16])

    nx = model.x.rows()
    nu = model.u.rows()

    # set prediction horizon
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = Tf

    # cost matrices
    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2, 1e-2])

    # path cost
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.model.cost_y_expr = ca.vertcat(model.x, model.u)
    ocp.cost.yref = np.zeros((nx+nu,))
    ocp.cost.W = ca.diagcat(Q_mat, R_mat).full()

    # terminal cost
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.cost.yref_e = np.zeros((nx,))
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.W_e = Q_mat

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0, np.pi, 0, 0])  # initial state
    ocp.constraints.idxbx_0 = np.array([0, 1, 2, 3])

    # set partial bounds for state
    ocp.constraints.lbx = np.array([-10,-10])
    ocp.constraints.ubx = np.array([10,10])
    ocp.constraints.idxbx = np.array([0,3])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'EXACT' # 'GAUSS_NEWTON', 'EXACT'
    ocp.solver_options.regularize_method = 'CONVEXIFY' # 'NO_REGULARIZATION', 'PROJECT', 'CONVEXIFY'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.globalization = 'MERIT_BACKTRACKING' # turns on globalization

    return ocp

def main():
    N_horizon = 20
    Tf = 1.0
    ocp = formulate_ocp(Tf, N_horizon)

    initial_iterate = ocp.create_default_initial_iterate()

    ## solve using acados
    # create acados solver
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    ocp_solver.set_iterate(initial_iterate)
    # solve with acados
    status = ocp_solver.solve()

    # test dump with default qp_type - compare Python and C backends
    ocp_solver.dump_last_qp_to_json(filename='last_qp_pendulum_python.json', overwrite=True, qp_type='default', backend='Python')
    ocp_solver.dump_last_qp_to_json(filename='last_qp_pendulum_C.json', overwrite=True, qp_type='default', backend='C')
    #compare the two json files
    with open('last_qp_pendulum_python.json', 'r') as f:
        python_json = json.load(f)
    with open('last_qp_pendulum_C.json', 'r') as f:
        C_json = json.load(f)

    assert python_json.keys() == C_json.keys(), f"Key mismatch: {python_json.keys() ^ C_json.keys()}"
    for k in python_json:
        np.testing.assert_allclose(python_json[k], C_json[k], atol=1e-6, equal_nan=True, err_msg=f"Error in {k}")
    os.remove('last_qp_pendulum_python.json')
    os.remove('last_qp_pendulum_C.json')
    print("dump default QP: Python and C backends match.")

    # test dump with scaled qp_type
    ocp_scaled = formulate_ocp(Tf, N_horizon)
    ocp_scaled.solver_options.qpscaling_scale_objective = 'OBJECTIVE_GERSHGORIN'
    ocp_scaled.solver_options.qpscaling_scale_constraints = 'INF_NORM'
    ocp_solver_scaled = AcadosOcpSolver(ocp_scaled, verbose=False)
    ocp_solver_scaled.solve()

    ocp_solver_scaled.dump_last_qp_to_json(filename='last_qp_pendulum_scaled_C.json', overwrite=True, qp_type='scaled', backend='C')
    with open('last_qp_pendulum_scaled_C.json', 'r') as f:
        scaled_C_json = json.load(f)
    assert len(scaled_C_json) > 0, "Scaled QP JSON should not be empty."

    scaled_qp = ocp_solver_scaled.get_last_scaled_qp()
    assert len(scaled_qp) > 0, "get_last_scaled_qp should return non-empty dict."
    os.remove('last_qp_pendulum_scaled_C.json')
    print("dump scaled QP: C backend works and get_last_scaled_qp works.")


if __name__ == "__main__":
    main()
