import json
import sys
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

def test_dump_default_qp():
    """Test dumping the default QP from both Python and C backends and compare results."""
    N_horizon = 20
    Tf = 1.0
    ocp = formulate_ocp(Tf, N_horizon)

    initial_iterate = ocp.create_default_initial_iterate()

    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    ocp_solver.set_iterate(initial_iterate)
    ocp_solver.solve()
    ocp_solver.dump_last_qp_to_json(filename='last_qp_pendulum_python.json', overwrite=True, backend='Python')
    ocp_solver.dump_last_qp_to_json(filename='last_qp_pendulum_C.json', overwrite=True, backend='C')

    with open('last_qp_pendulum_python.json', 'r') as f:
        python_json = json.load(f)
    with open('last_qp_pendulum_C.json', 'r') as f:
        C_json = json.load(f)

    assert python_json.keys() == C_json.keys(), f"Key mismatch: {python_json.keys() ^ C_json.keys()}"
    for k in python_json:
        np.testing.assert_allclose(python_json[k], C_json[k], atol=1e-6, equal_nan=True, err_msg=f"Error in {k}")
    print("test_dump_default_qp passed.")


def test_dump_scaled_qp():
    """Test dumping the scaled QP requires qpscaling to be active."""
    N_horizon = 20
    Tf = 1.0
    ocp = formulate_ocp(Tf, N_horizon)
    ocp.solver_options.qpscaling_scale_objective = 'OBJECTIVE_GERSHGORIN'
    ocp.solver_options.qpscaling_scale_constraints = 'INF_NORM'

    initial_iterate = ocp.create_default_initial_iterate()

    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    ocp_solver.set_iterate(initial_iterate)
    ocp_solver.solve()
    ocp_solver.dump_last_qp_to_json(filename='last_scaled_qp_C.json', overwrite=True, backend='C', qp_type='scaled')
    print("test_dump_scaled_qp passed.")


def test_dump_qp_type_validation():
    """Test that invalid qp_type values raise appropriate errors."""
    N_horizon = 20
    Tf = 1.0
    ocp = formulate_ocp(Tf, N_horizon)

    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    ocp_solver.solve()

    # 'relaxed' requires SQP_WITH_FEASIBLE_QP
    try:
        ocp_solver.dump_last_qp_to_json(filename='should_not_exist.json', overwrite=True, backend='C', qp_type='relaxed')
        raise AssertionError("Expected ValueError for qp_type='relaxed' with SQP solver")
    except ValueError as e:
        print(f"Correctly raised ValueError: {e}")

    # 'scaled' requires qpscaling to be active
    try:
        ocp_solver.dump_last_qp_to_json(filename='should_not_exist.json', overwrite=True, backend='C', qp_type='scaled')
        raise AssertionError("Expected ValueError for qp_type='scaled' without qpscaling")
    except ValueError as e:
        print(f"Correctly raised ValueError: {e}")
    print("test_dump_qp_type_validation passed.")


def main():
    test_dump_default_qp()
    test_dump_scaled_qp()
    test_dump_qp_type_validation()

if __name__ == "__main__":
    main()
