import json
import sys
import os
sys.path.insert(0, '../getting_started')
import numpy as np
import casadi as ca

from acados_template import AcadosOcp, AcadosOcpSolver
from pendulum_model import export_pendulum_ode_model

def formulate_ocp(Tf: float = 1.0, N: int = 20) -> AcadosOcp:
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()

    # set prediction horizon
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = Tf

    # cost matrices
    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2])

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

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP'

    return ocp


def test_dump_qp_default():
    """Test dump_last_qp_to_json with default qp_type, comparing Python and C backends."""
    N_horizon = 20
    Tf = 1.0
    ocp = formulate_ocp(Tf, N_horizon)

    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    ocp_solver.solve()

    python_file = 'test_qp_default_python.json'
    c_file = 'test_qp_default_C.json'

    ocp_solver.dump_last_qp_to_json(filename=python_file, overwrite=True, qp_type='default', backend='Python')
    ocp_solver.dump_last_qp_to_json(filename=c_file, overwrite=True, qp_type='default', backend='C')

    with open(python_file, 'r') as f:
        python_json = json.load(f)
    with open(c_file, 'r') as f:
        c_json = json.load(f)

    assert python_json.keys() == c_json.keys(), f"Key mismatch: {python_json.keys() ^ c_json.keys()}"
    for k in python_json:
        np.testing.assert_allclose(python_json[k], c_json[k], atol=1e-6, equal_nan=True, err_msg=f"Error in field {k}")

    # cleanup
    os.remove(python_file)
    os.remove(c_file)
    print("test_dump_qp_default passed.")


def test_dump_qp_scaled():
    """Test dump_last_qp_to_json with scaled qp_type when QP scaling is enabled."""
    N_horizon = 20
    Tf = 1.0
    ocp = formulate_ocp(Tf, N_horizon)

    # enable QP scaling
    ocp.solver_options.qpscaling_scale_objective = 'OBJECTIVE_GERSHGORIN'
    ocp.solver_options.qpscaling_scale_constraints = 'INF_NORM'

    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    ocp_solver.solve()

    c_file = 'test_qp_scaled_C.json'
    ocp_solver.dump_last_qp_to_json(filename=c_file, overwrite=True, qp_type='scaled', backend='C')

    with open(c_file, 'r') as f:
        c_json = json.load(f)

    assert len(c_json) > 0, "Scaled QP JSON should not be empty."

    # also test get_last_scaled_qp
    scaled_qp = ocp_solver.get_last_scaled_qp()
    assert len(scaled_qp) > 0, "get_last_scaled_qp should return non-empty dict."

    # cleanup
    os.remove(c_file)
    print("test_dump_qp_scaled passed.")


def main():
    test_dump_qp_default()
    test_dump_qp_scaled()
    print("\nAll dump QP tests passed!")


if __name__ == "__main__":
    main()
