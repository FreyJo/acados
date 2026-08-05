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

"""
Tests that verify AcadosOcpBatchSolver and AcadosSimBatchSolver produce the same
results as calling the corresponding non-batch solvers individually.
"""

import sys
sys.path.insert(0, '../pendulum_on_cart/common')

import numpy as np
import scipy.linalg
import casadi as ca

from acados_template import (
    AcadosOcp, AcadosOcpSolver, AcadosOcpBatchSolver,
    AcadosSim, AcadosSimSolver, AcadosSimBatchSolver,
)
from pendulum_model import export_pendulum_ode_model

TOL = 1e-7
N_BATCH = 5


# ---------------------------------------------------------------------------
# OCP helpers
# ---------------------------------------------------------------------------

def setup_ocp(tol: float = TOL) -> AcadosOcp:
    ocp = AcadosOcp()
    ocp.model = export_pendulum_ode_model()

    Tf = 1.0
    N = 20
    ocp.solver_options.N_horizon = N

    Q_mat = 2 * np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2 * np.diag([1e-2])
    cost_W = scipy.linalg.block_diag(Q_mat, R_mat)

    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.cost.W_e = Q_mat
    ocp.cost.W = cost_W
    ocp.model.cost_y_expr = ca.vertcat(ocp.model.x, ocp.model.u)
    ocp.model.cost_y_expr_e = ocp.model.x
    ocp.cost.yref = np.zeros(ocp.model.cost_y_expr.shape).flatten()
    ocp.cost.yref_e = np.zeros(ocp.model.cost_y_expr_e.shape).flatten()

    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])
    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'IRK'
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.nlp_solver_tol_stat = tol
    ocp.solver_options.nlp_solver_tol_eq = tol
    ocp.solver_options.nlp_solver_tol_ineq = tol
    ocp.solver_options.nlp_solver_tol_comp = tol
    ocp.solver_options.tf = Tf
    ocp.solver_options.with_batch_functionality = True

    return ocp


def solve_ocp_individually(x0_batch: np.ndarray, tol: float = TOL):
    """Solve OCPs one-by-one and return collected results."""
    ocp = setup_ocp(tol=tol)
    solver = AcadosOcpSolver(ocp, verbose=False)

    N_batch = x0_batch.shape[0]
    nx = ocp.dims.nx
    N_horizon = ocp.solver_options.N_horizon

    x_flat = []
    u_flat = []
    x_stage1 = []  # x at stage 1 for all batch members
    u_stage0 = []  # u at stage 0 for all batch members

    for i in range(N_batch):
        solver.solve_for_x0(x0_bar=x0_batch[i], fail_on_nonzero_status=False)
        x_flat.append(solver.get_flat('x'))
        u_flat.append(solver.get_flat('u'))
        x_stage1.append(solver.get(1, 'x'))
        u_stage0.append(solver.get(0, 'u'))

    return {
        'x_flat': np.array(x_flat),
        'u_flat': np.array(u_flat),
        'x_stage1': np.array(x_stage1),
        'u_stage0': np.array(u_stage0),
    }


def solve_ocp_batch(x0_batch: np.ndarray, tol: float = TOL):
    """Solve OCPs using the batch solver and return the same results."""
    ocp = setup_ocp(tol=tol)
    batch_solver = AcadosOcpBatchSolver(ocp, N_batch_init=x0_batch.shape[0], verbose=False)

    batch_solver.constraints_set(0, 'lbx', x0_batch)
    batch_solver.constraints_set(0, 'ubx', x0_batch)
    batch_solver.solve()

    x_flat = batch_solver.get_flat('x')
    u_flat = batch_solver.get_flat('u')
    x_stage1 = batch_solver.get(1, 'x')
    u_stage0 = batch_solver.get(0, 'u')

    return {
        'x_flat': x_flat,
        'u_flat': u_flat,
        'x_stage1': x_stage1,
        'u_stage0': u_stage0,
        'batch_solver': batch_solver,
    }


# ---------------------------------------------------------------------------
# Sim helpers
# ---------------------------------------------------------------------------

def setup_sim() -> AcadosSim:
    sim = AcadosSim()
    sim.model = export_pendulum_ode_model()
    sim.solver_options.T = 0.2
    sim.solver_options.integrator_type = 'IRK'
    sim.solver_options.num_stages = 3
    sim.solver_options.num_steps = 3
    sim.solver_options.newton_iter = 10
    return sim


def simulate_individually(x_batch: np.ndarray, u_batch: np.ndarray):
    """Simulate each (x, u) pair with a separate integrator call."""
    sim = setup_sim()
    integrator = AcadosSimSolver(sim, verbose=False)

    N_batch = x_batch.shape[0]
    nx = x_batch.shape[1]
    x_next = np.zeros((N_batch, nx))
    for i in range(N_batch):
        x_next[i] = integrator.simulate(x=x_batch[i], u=u_batch[i])
    return x_next


def simulate_batch(x_batch: np.ndarray, u_batch: np.ndarray):
    """Simulate each (x, u) pair using AcadosSimBatchSolver."""
    N_batch = x_batch.shape[0]
    sim = setup_sim()
    batch_integrator = AcadosSimBatchSolver(sim, N_batch=N_batch, verbose=False)

    for n in range(N_batch):
        batch_integrator.sim_solvers[n].set('x', x_batch[n])
        batch_integrator.sim_solvers[n].set('u', u_batch[n])

    batch_integrator.solve()

    nx = x_batch.shape[1]
    x_next = np.zeros((N_batch, nx))
    for n in range(N_batch):
        x_next[n] = batch_integrator.sim_solvers[n].get('x')
    return x_next


# ---------------------------------------------------------------------------
# Test functions
# ---------------------------------------------------------------------------

def test_ocp_batch_get_matches_individual():
    """AcadosOcpBatchSolver.get() should match individual AcadosOcpSolver.get()."""
    x0_batch = np.tile(np.array([0.0, np.pi, 0.0, 0.0]), (N_BATCH, 1))
    # vary initial states slightly so each problem is different
    rng = np.random.default_rng(0)
    x0_batch += 0.1 * rng.standard_normal(x0_batch.shape)

    ref = solve_ocp_individually(x0_batch)
    batch = solve_ocp_batch(x0_batch)

    err_x = np.max(np.abs(ref['x_stage1'] - batch['x_stage1']))
    err_u = np.max(np.abs(ref['u_stage0'] - batch['u_stage0']))

    assert err_x < TOL * 100, f"get(1, 'x'): max error {err_x} exceeds tolerance"
    assert err_u < TOL * 100, f"get(0, 'u'): max error {err_u} exceeds tolerance"
    print(f"test_ocp_batch_get_matches_individual: PASSED (err_x={err_x:.2e}, err_u={err_u:.2e})")


def test_ocp_batch_get_flat_matches_individual():
    """AcadosOcpBatchSolver.get_flat() should match individual AcadosOcpSolver.get_flat()."""
    x0_batch = np.tile(np.array([0.0, np.pi, 0.0, 0.0]), (N_BATCH, 1))
    rng = np.random.default_rng(1)
    x0_batch += 0.1 * rng.standard_normal(x0_batch.shape)

    ref = solve_ocp_individually(x0_batch)
    batch = solve_ocp_batch(x0_batch)

    err_x = np.max(np.abs(ref['x_flat'] - batch['x_flat']))
    err_u = np.max(np.abs(ref['u_flat'] - batch['u_flat']))

    assert err_x < TOL * 100, f"get_flat('x'): max error {err_x} exceeds tolerance"
    assert err_u < TOL * 100, f"get_flat('u'): max error {err_u} exceeds tolerance"
    print(f"test_ocp_batch_get_flat_matches_individual: PASSED (err_x={err_x:.2e}, err_u={err_u:.2e})")


def test_ocp_batch_eval_solution_sensitivity_wrt_initial_state():
    """
    AcadosOcpBatchSolver.eval_solution_sensitivity(with_respect_to='initial_state')
    should produce the same sensitivities as individual AcadosOcpSolver calls.
    """
    x0_batch = np.tile(np.array([0.0, np.pi, 0.0, 0.0]), (N_BATCH, 1))
    rng = np.random.default_rng(2)
    x0_batch += 0.1 * rng.standard_normal(x0_batch.shape)

    # --- individual solvers ---
    ocp = setup_ocp()
    solver = AcadosOcpSolver(ocp, verbose=False)

    nx = ocp.dims.nx
    nu = ocp.dims.nu
    ref_sens_u = np.zeros((N_BATCH, nu, nx))
    ref_sens_x = np.zeros((N_BATCH, nx, nx))

    for i in range(N_BATCH):
        solver.solve_for_x0(x0_bar=x0_batch[i], fail_on_nonzero_status=False)
        solver.setup_qp_matrices_and_factorize()
        d = solver.eval_solution_sensitivity(
            0,
            with_respect_to='initial_state',
            return_sens_x=True,
            return_sens_u=True,
        )
        ref_sens_u[i] = d['sens_u']
        ref_sens_x[i] = d['sens_x']

    # --- batch solver ---
    ocp_b = setup_ocp()
    batch_solver = AcadosOcpBatchSolver(ocp_b, N_batch_init=N_BATCH, verbose=False)

    batch_solver.constraints_set(0, 'lbx', x0_batch)
    batch_solver.constraints_set(0, 'ubx', x0_batch)
    batch_solver.solve()
    batch_solver.setup_qp_matrices_and_factorize()

    batch_d = batch_solver.eval_solution_sensitivity(
        0,
        with_respect_to='initial_state',
        return_sens_x=True,
        return_sens_u=True,
    )

    err_u = np.max(np.abs(ref_sens_u - batch_d['sens_u']))
    err_x = np.max(np.abs(ref_sens_x - batch_d['sens_x']))

    assert err_u < TOL * 100, f"eval_solution_sensitivity sens_u: max error {err_u} exceeds tolerance"
    assert err_x < TOL * 100, f"eval_solution_sensitivity sens_x: max error {err_x} exceeds tolerance"
    print(f"test_ocp_batch_eval_solution_sensitivity_wrt_initial_state: PASSED (err_u={err_u:.2e}, err_x={err_x:.2e})")


def test_ocp_batch_eval_solution_sensitivity_wrt_p_global():
    """
    AcadosOcpBatchSolver.eval_solution_sensitivity(with_respect_to='p_global') should
    produce the same sensitivities as individual AcadosOcpSolver calls.
    Uses the parametric convex OCP from solution_sensitivities_convex_example.
    """
    sys.path.insert(0, '../solution_sensitivities_convex_example')
    from setup_parametric_ocp import export_parametric_ocp

    nx, nu = 4, 2
    learnable_params = ["A", "Q", "b"]

    x0_batch = 0.1 * (-1) ** np.arange(nx)
    x0_batch = np.tile(x0_batch, (N_BATCH, 1))
    rng = np.random.default_rng(5)
    x0_batch += 0.05 * rng.standard_normal(x0_batch.shape)

    # --- individual solvers ---
    ocp_ref = export_parametric_ocp(nx, nu, learnable_params=learnable_params)
    ocp_ref.code_gen_options.with_solution_sens_wrt_params_forw = True
    ocp_ref.solver_options.qp_solver_ric_alg = 0
    solver_ref = AcadosOcpSolver(ocp_ref, verbose=False)

    n_p_global = ocp_ref.p_global_values.shape[0]
    ref_sens_u = np.zeros((N_BATCH, nu, n_p_global))
    ref_sens_x = np.zeros((N_BATCH, nx, n_p_global))

    for i in range(N_BATCH):
        solver_ref.reset()
        solver_ref.solve_for_x0(x0_bar=x0_batch[i], fail_on_nonzero_status=False)
        solver_ref.setup_qp_matrices_and_factorize()
        d = solver_ref.eval_solution_sensitivity(
            0,
            with_respect_to='p_global',
            return_sens_x=True,
            return_sens_u=True,
        )
        ref_sens_u[i] = d['sens_u']
        ref_sens_x[i] = d['sens_x']

    # --- batch solver ---
    ocp_b = export_parametric_ocp(nx, nu, learnable_params=learnable_params)
    ocp_b.code_gen_options.with_solution_sens_wrt_params_forw = True
    ocp_b.solver_options.qp_solver_ric_alg = 0
    ocp_b.solver_options.with_batch_functionality = True

    batch_solver = AcadosOcpBatchSolver(ocp_b, N_batch_init=N_BATCH, verbose=False)
    batch_solver.constraints_set(0, 'lbx', x0_batch)
    batch_solver.constraints_set(0, 'ubx', x0_batch)
    for n in range(N_BATCH):
        batch_solver.ocp_solvers[n].reset()
    batch_solver.solve()
    batch_solver.setup_qp_matrices_and_factorize()

    batch_d = batch_solver.eval_solution_sensitivity(
        0,
        with_respect_to='p_global',
        return_sens_x=True,
        return_sens_u=True,
    )

    err_u = np.max(np.abs(ref_sens_u - batch_d['sens_u']))
    err_x = np.max(np.abs(ref_sens_x - batch_d['sens_x']))

    assert err_u < TOL * 100, f"eval_solution_sensitivity (p_global) sens_u: max error {err_u} exceeds tolerance"
    assert err_x < TOL * 100, f"eval_solution_sensitivity (p_global) sens_x: max error {err_x} exceeds tolerance"
    print(f"test_ocp_batch_eval_solution_sensitivity_wrt_p_global: PASSED (err_u={err_u:.2e}, err_x={err_x:.2e})")


def test_ocp_batch_eval_adjoint_solution_sensitivity():
    """
    AcadosOcpBatchSolver.eval_adjoint_solution_sensitivity() should produce the same
    adjoint sensitivities as individual AcadosOcpSolver calls.
    Uses the parametric convex OCP from solution_sensitivities_convex_example.
    """
    sys.path.insert(0, '../solution_sensitivities_convex_example')
    from setup_parametric_ocp import export_parametric_ocp

    nx, nu = 4, 2
    learnable_params = ["A", "Q", "b"]

    x0_batch = 0.1 * (-1) ** np.arange(nx)
    x0_batch = np.tile(x0_batch, (N_BATCH, 1))
    rng = np.random.default_rng(6)
    x0_batch += 0.05 * rng.standard_normal(x0_batch.shape)

    seed_x_val = np.ones((nx, 1))
    seed_u_val = np.ones((nu, 1))

    # --- individual solvers ---
    ocp_ref = export_parametric_ocp(nx, nu, learnable_params=learnable_params)
    ocp_ref.code_gen_options.with_solution_sens_wrt_params_adj = True
    ocp_ref.solver_options.qp_solver_ric_alg = 0
    solver_ref = AcadosOcpSolver(ocp_ref, verbose=False)

    n_p_global = ocp_ref.p_global_values.shape[0]
    ref_adj = np.zeros((N_BATCH, 1, n_p_global))

    for i in range(N_BATCH):
        solver_ref.reset()
        solver_ref.solve_for_x0(x0_bar=x0_batch[i], fail_on_nonzero_status=False)
        solver_ref.setup_qp_matrices_and_factorize()
        adj = solver_ref.eval_adjoint_solution_sensitivity(
            seed_x=[(1, seed_x_val)],
            seed_u=[(0, seed_u_val)],
        )
        ref_adj[i] = adj

    # --- batch solver ---
    ocp_b = export_parametric_ocp(nx, nu, learnable_params=learnable_params)
    ocp_b.code_gen_options.with_solution_sens_wrt_params_adj = True
    ocp_b.solver_options.qp_solver_ric_alg = 0
    ocp_b.solver_options.with_batch_functionality = True

    batch_solver = AcadosOcpBatchSolver(ocp_b, N_batch_init=N_BATCH, verbose=False)
    batch_solver.constraints_set(0, 'lbx', x0_batch)
    batch_solver.constraints_set(0, 'ubx', x0_batch)
    for n in range(N_BATCH):
        batch_solver.ocp_solvers[n].reset()
    batch_solver.solve()
    batch_solver.setup_qp_matrices_and_factorize()

    # Batch seeds have shape (n_batch, dim, n_seeds)
    batch_seed_x = np.tile(seed_x_val[np.newaxis, :, :], (N_BATCH, 1, 1))
    batch_seed_u = np.tile(seed_u_val[np.newaxis, :, :], (N_BATCH, 1, 1))

    batch_adj = batch_solver.eval_adjoint_solution_sensitivity(
        seed_x=[(1, batch_seed_x)],
        seed_u=[(0, batch_seed_u)],
    )

    err = np.max(np.abs(ref_adj - batch_adj))
    assert err < TOL * 100, f"eval_adjoint_solution_sensitivity: max error {err} exceeds tolerance"
    print(f"test_ocp_batch_eval_adjoint_solution_sensitivity: PASSED (err={err:.2e})")


def test_sim_batch_matches_individual():
    """AcadosSimBatchSolver results should match individual AcadosSimSolver results."""
    rng = np.random.default_rng(4)
    x_batch = rng.standard_normal((N_BATCH, 4))
    u_batch = rng.standard_normal((N_BATCH, 1))

    x_next_ref = simulate_individually(x_batch, u_batch)
    x_next_batch = simulate_batch(x_batch, u_batch)

    err = np.max(np.abs(x_next_ref - x_next_batch))
    assert err < 1e-10, f"SimBatchSolver: max error {err} exceeds tolerance"
    print(f"test_sim_batch_matches_individual: PASSED (err={err:.2e})")


def test_sim_batch_num_threads_property():
    """Verify num_threads_in_batch_solve property getter/setter on AcadosSimBatchSolver."""
    sim = setup_sim()
    batch_integrator = AcadosSimBatchSolver(sim, N_batch=2, num_threads_in_batch_solve=1, verbose=False)

    assert batch_integrator.num_threads_in_batch_solve == 1
    batch_integrator.num_threads_in_batch_solve = 2
    assert batch_integrator.num_threads_in_batch_solve == 2
    print("test_sim_batch_num_threads_property: PASSED")


def test_ocp_batch_num_threads_property():
    """Verify num_threads_in_batch_solve property getter/setter on AcadosOcpBatchSolver."""
    ocp = setup_ocp()
    batch_solver = AcadosOcpBatchSolver(ocp, N_batch_init=2, num_threads_in_batch_solve=1, verbose=False)

    assert batch_solver.num_threads_in_batch_solve == 1
    batch_solver.num_threads_in_batch_solve = 2
    assert batch_solver.num_threads_in_batch_solve == 2
    print("test_ocp_batch_num_threads_property: PASSED")


def test_ocp_batch_exceeding_n_batch_max_raises():
    """Requesting more samples than N_batch_max should raise ValueError."""
    x0 = np.array([[0.0, np.pi, 0.0, 0.0]])
    ocp = setup_ocp()
    N_batch_max = 2
    batch_solver = AcadosOcpBatchSolver(ocp, N_batch_init=N_batch_max, verbose=False)

    batch_solver.constraints_set(0, 'lbx', np.tile(x0, (N_batch_max, 1)))
    batch_solver.constraints_set(0, 'ubx', np.tile(x0, (N_batch_max, 1)))
    batch_solver.solve()

    try:
        batch_solver.get_flat('x', N_batch_max + 1)
        raise AssertionError("Expected ValueError was not raised")
    except ValueError:
        pass
    print("test_ocp_batch_exceeding_n_batch_max_raises: PASSED")


if __name__ == '__main__':
    test_ocp_batch_get_matches_individual()
    test_ocp_batch_get_flat_matches_individual()
    test_ocp_batch_eval_solution_sensitivity_wrt_initial_state()
    test_ocp_batch_eval_solution_sensitivity_wrt_p_global()
    test_ocp_batch_eval_adjoint_solution_sensitivity()
    test_sim_batch_matches_individual()
    test_sim_batch_num_threads_property()
    test_ocp_batch_num_threads_property()
    test_ocp_batch_exceeding_n_batch_max_raises()
    print("\nAll batch solver tests passed.")
