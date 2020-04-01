#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schoels, Jonas Schlagenhauf, Moritz Diehl
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

import sys, os, json
import numpy as np

from ctypes import *

from copy import deepcopy


class AcadosOcpSolver:
    """
    class to interact with the acados ocp solver C object
    """
    def __init__(self, acados_ocp, json_file='acados_ocp_nlp.json', use_compiled_solver=False):

        self.solver_created = False
        model = acados_ocp["model"]

        self.shared_lib_name = 'c_generated_code/libacados_ocp_solver_' + model["name"] + '.so'

        # get
        self.shared_lib = CDLL(self.shared_lib_name)
        self.shared_lib.acados_create()
        self.solver_created = True

        self.shared_lib.acados_get_nlp_opts.restype = c_void_p
        self.nlp_opts = self.shared_lib.acados_get_nlp_opts()

        self.shared_lib.acados_get_nlp_dims.restype = c_void_p
        self.nlp_dims = self.shared_lib.acados_get_nlp_dims()

        self.shared_lib.acados_get_nlp_config.restype = c_void_p
        self.nlp_config = self.shared_lib.acados_get_nlp_config()

        self.shared_lib.acados_get_nlp_out.restype = c_void_p
        self.nlp_out = self.shared_lib.acados_get_nlp_out()

        self.shared_lib.acados_get_nlp_in.restype = c_void_p
        self.nlp_in = self.shared_lib.acados_get_nlp_in()

        self.shared_lib.acados_get_nlp_solver.restype = c_void_p
        self.nlp_solver = self.shared_lib.acados_get_nlp_solver()

        self.acados_ocp = acados_ocp


    def solve(self, rti_phase=0):
        """
        solve the ocp with current input
        :param rti_phase: 0 = preparation + feedback, 1 = preparation only,
         2 = feedback only (if SQP_RTI is used, otherwise only 0 (default) is allowed)
        """
        if isinstance(rti_phase, int) == False or rti_phase < 0 or rti_phase > 2: 
            raise Exception('AcadosOcpSolver.solve(): argument \'rti_phase\' can ' 
                'take only values 0, 1, 2 for SQP-RTI-type solvers')
        self.shared_lib.acados_solve.argtypes = [c_int]

        status = self.shared_lib.acados_solve(rti_phase)
        return status


    def get(self, stage_, field_):
        """
        get the last solution of the solver:
            :param stage: integer corresponding to shooting node
            :param field_: string in ['x', 'u', 'z', 'pi']
        """

        out_fields = ['x', 'u', 'z', 'pi']
        field = field_
        field = field.encode('utf-8')

        if (field_ not in out_fields):
            raise Exception('AcadosOcpSolver.get(): {} is an invalid argument.\
                    \n Possible values are {}. Exiting.'.format(out_fields))

        self.shared_lib.ocp_nlp_dims_get_from_attr.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p]
        self.shared_lib.ocp_nlp_dims_get_from_attr.restype = c_int

        dims = self.shared_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage_, field)

        out = np.ascontiguousarray(np.zeros((dims,)), dtype=np.float64)
        out_data = cast(out.ctypes.data, POINTER(c_double))

        self.shared_lib.ocp_nlp_out_get.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
        self.shared_lib.ocp_nlp_out_get(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage_, field, out_data)

        return out


    def print_statistics(self):
        stat = self.get_stats("statistics")

        if self.acados_ocp["solver_options"]["nlp_solver_type"] == 'SQP':
            print('\niter\tres_stat\tres_eq\t\tres_ineq\tres_comp\tqp_stat\tqp_iter')
            if stat.shape[0]>7:
                print('\tqp_res_stat\tqp_res_eq\tqp_res_ineq\tqp_res_comp')
            for jj in range(stat.shape[1]):
                print('{:d}\t{:e}\t{:e}\t{:e}\t{:e}\t{:d}\t{:d}'.format( \
                     int(stat[0][jj]), stat[1][jj], stat[2][jj], \
                     stat[3][jj], stat[4][jj], int(stat[5][jj]), int(stat[6][jj])))
                if stat.shape[0]>7:
                    print('\t{:e}\t{:e}\t{:e}\t{:e}'.format( \
                        stat[7][jj], stat[8][jj], stat[9][jj], stat[10][jj]))
            print('\n')
        elif self.acados_ocp["solver_options"]["nlp_solver_type"] == 'SQP_RTI':
            print('\niter\tqp_stat\tqp_iter')
            if stat.shape[0]>3:
                print('\tqp_res_stat\tqp_res_eq\tqp_res_ineq\tqp_res_comp')
            for jj in range(stat.shape[1]):
                print('{:d}\t{:d}\t{:d}'.format( int(stat[0][jj]), int(stat[1][jj]), int(stat[2][jj])))
                if stat.shape[0]>3:
                    print('\t{:e}\t{:e}\t{:e}\t{:e}'.format( \
                         stat[3][jj], stat[4][jj], stat[5][jj], stat[6][jj]))
            print('\n')

        return


    def get_stats(self, field_):
        """
        get the information of the last solver call:
            :param field_: string in ['time_tot', 'time_lin', 'time_qp', 'time_reg', 'sqp_iter']
        """

        fields = ['time_tot',  # total cpu time previous call
                  'time_lin',  # cpu time for linearization
                  'time_qp',   # cpu time qp solution
                  'time_qp_solver_call',  # cpu time inside qp solver (without converting the QP)
                  'time_reg',  # cpu time regularization
                  'sqp_iter',  # number of SQP iterations
                  'statistics',  # table with info about last iteration
                  'stat_m',
                  'stat_n',
                ]

        field = field_
        field = field.encode('utf-8')
        if (field_ not in fields):
            raise Exception('AcadosOcpSolver.get_stats(): {} is not a valid argument.\
                    \n Possible values are {}. Exiting.'.format(fields, fields))

        if field_ in ['sqp_iter', 'stat_m', 'stat_n']:
            out = np.ascontiguousarray(np.zeros((1,)), dtype=np.int64)
            out_data = cast(out.ctypes.data, POINTER(c_int64))

        elif field_ == 'statistics':
            sqp_iter = self.get_stats("sqp_iter")
            stat_m = self.get_stats("stat_m")
            stat_n = self.get_stats("stat_n")

            min_size = min([stat_m, sqp_iter+1])

            out = np.ascontiguousarray(
                        np.zeros( (stat_n[0]+1, min_size[0]) ), dtype=np.float64)
            out_data = cast(out.ctypes.data, POINTER(c_double))

        else:
            out = np.ascontiguousarray(np.zeros((1,)), dtype=np.float64)
            out_data = cast(out.ctypes.data, POINTER(c_double))

        self.shared_lib.ocp_nlp_get.argtypes = [c_void_p, c_void_p, c_char_p, c_void_p]
        self.shared_lib.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, out_data)

        return out

    # Note: this function should not be used anymore, better use cost_set, constraints_set
    def set(self, stage_, field_, value_):

        cost_fields = ['y_ref', 'yref']
        constraints_fields = ['lbx', 'ubx', 'lbu', 'ubu']
        out_fields = ['x', 'u', 'pi']

        # cast value_ to avoid conversion issues
        value_ = value_.astype(float)

        field = field_
        field = field.encode('utf-8')

        stage = c_int(stage_)

        # treat parameters separately
        if field_ is 'p':
            self.shared_lib.acados_update_params.argtypes = [c_int, POINTER(c_double)]
            self.shared_lib.acados_update_params.restype = c_int
            value_data = cast(value_.ctypes.data, POINTER(c_double))
            self.shared_lib.acados_update_params(stage, value_data, value_.shape[0])
        else:
            if (field_ not in constraints_fields) and \
                    (field_ not in cost_fields) and (field_ not in out_fields):
                raise Exception("AcadosOcpSolver.set(): {} is not a valid argument.\
                    \nPossible values are {} and {}. Exiting.".format(field, \
                    cost_fields, constraints_fields, out_fields))

            self.shared_lib.ocp_nlp_dims_get_from_attr.argtypes = \
                [c_void_p, c_void_p, c_void_p, c_int, c_char_p]
            self.shared_lib.ocp_nlp_dims_get_from_attr.restype = c_int

            dims = self.shared_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, \
                self.nlp_dims, self.nlp_out, stage_, field)

            if value_.shape[0] != dims:
                msg = 'AcadosOcpSolver.set(): mismatching dimension for field "{}" '.format(field_)
                msg += 'with dimension {} (you have {})'.format(dims, value_.shape[0])
                raise Exception(msg)

            value_data = cast(value_.ctypes.data, POINTER(c_double))
            value_data_p = cast((value_data), c_void_p)

            if field_ in constraints_fields:
                self.shared_lib.ocp_nlp_constraints_model_set.argtypes = \
                    [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
                self.shared_lib.ocp_nlp_constraints_model_set(self.nlp_config, \
                    self.nlp_dims, self.nlp_in, stage, field, value_data_p)
            elif field_ in cost_fields:
                self.shared_lib.ocp_nlp_cost_model_set.argtypes = \
                    [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
                self.shared_lib.ocp_nlp_cost_model_set(self.nlp_config, \
                    self.nlp_dims, self.nlp_in, stage, field, value_data_p)
            elif field_ in out_fields:
                self.shared_lib.ocp_nlp_out_set.argtypes = \
                    [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
                self.shared_lib.ocp_nlp_out_set(self.nlp_config, \
                    self.nlp_dims, self.nlp_out, stage, field, value_data_p)

        return


    def cost_set(self, stage_, field_, value_):
        """
        set numerical data in the cost module of the solver:
            :param stage_: integer corresponding to shooting node
            :param field_: string, e.g. 'yref', 'W'
            :param value_: of appropriate size
        """
        # cast value_ to avoid conversion issues
        value_ = value_.astype(float)

        field = field_
        field = field.encode('utf-8')

        stage = c_int(stage_)
        self.shared_lib.ocp_nlp_cost_dims_get_from_attr.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p, POINTER(c_int)]
        self.shared_lib.ocp_nlp_cost_dims_get_from_attr.restype = c_int

        dims = np.ascontiguousarray(np.zeros((2,)), dtype=np.intc)
        dims_data = cast(dims.ctypes.data, POINTER(c_int))

        self.shared_lib.ocp_nlp_cost_dims_get_from_attr(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage_, field, dims_data)

        value_shape = value_.shape
        if len(value_shape) == 1:
            value_shape = (value_shape[0], 0)

        if value_shape != tuple(dims):
            raise Exception('acados_solver.set(): mismatching dimension', \
                ' for field "{}" with dimension {} (you have {})'.format( \
                field_, tuple(dims), value_shape))

        value_data = cast(value_.ctypes.data, POINTER(c_double))
        value_data_p = cast((value_data), c_void_p)

        self.shared_lib.ocp_nlp_cost_model_set.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
        self.shared_lib.ocp_nlp_cost_model_set(self.nlp_config, \
            self.nlp_dims, self.nlp_in, stage, field, value_data_p)

        return


    def constraints_set(self, stage_, field_, value_):
        """
        set numerical data in the constraint module of the solver:
        Parameters:
            :param stage_: integer corresponding to shooting node
            :param field_: string, e.g. 'lbx'
            :param value_: of appropriate size
        """
        # cast value_ to avoid conversion issues
        value_ = value_.astype(float)

        field = field_
        field = field.encode('utf-8')

        stage = c_int(stage_)
        self.shared_lib.ocp_nlp_constraint_dims_get_from_attr.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p, POINTER(c_int)]
        self.shared_lib.ocp_nlp_constraint_dims_get_from_attr.restype = c_int

        dims = np.ascontiguousarray(np.zeros((2,)), dtype=np.intc)
        dims_data = cast(dims.ctypes.data, POINTER(c_int))

        self.shared_lib.ocp_nlp_constraint_dims_get_from_attr(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage_, field, dims_data)

        value_shape = value_.shape
        if len(value_shape) == 1:
            value_shape = (value_shape[0], 0)

        if value_shape != tuple(dims):
            raise Exception('acados_solver.set(): mismatching dimension' \
                ' for field "{}" with dimension {} (you have {})'.format(field_, tuple(dims), value_shape))

        value_data = cast(value_.ctypes.data, POINTER(c_double))
        value_data_p = cast((value_data), c_void_p)

        self.shared_lib.ocp_nlp_constraints_model_set.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
        self.shared_lib.ocp_nlp_constraints_model_set(self.nlp_config, \
            self.nlp_dims, self.nlp_in, stage, field, value_data_p)

        return


    def options_set(self, field_, value_):
        """
        set options of the solver:
        Parameters:
            :param field_: string, e.g. 'print_level'
            :param value_: of type int
        """
        # cast value_ to avoid conversion issues
        if isinstance(value_, int) is not True:
            raise Exception('solver options must be of type int. You have {}'.format())

        field = field_
        field = field.encode('utf-8')

        value_ctypes = c_int(value_)

        self.shared_lib.ocp_nlp_solver_opts_set.argtypes = \
            [c_void_p, c_void_p, c_char_p, c_void_p]
        self.shared_lib.ocp_nlp_solver_opts_set(self.nlp_config, \
            self.nlp_opts, field, byref(value_ctypes))

        return


    def __del__(self):
        if self.solver_created:
            self.shared_lib.acados_free()
            del self.shared_lib

        # NOTE: DLL cannot be easily unloaded!!!
        # see https://stackoverflow.com/questions/359498/how-can-i-unload-a-dll-using-ctypes-in-python
        # while isLoaded(self.shared_lib_name):
        #     dlclose(handle)
