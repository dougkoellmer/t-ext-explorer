
import sympy
import math

import tmc_view
import tmc_utils

import traceback

from tmc_consts import *
from tmc_types import *
from tmc_utils import *
from tmc_view import *


from ai_code.compatible_basis import *
from ai_code.symbolic import *
from ai_code.reduction import *
from ai_code.fields import *

import ai_code.reduce_skew as reduce_skew

import tmc_phi_psi

import pprint

from ai_code import *


def reduce_to_finite_field(expression, p, n):

    if hasattr(expression, 'as_ordered_terms'):
        expression = expression.subs(Q_SYMBOL, p**n)
        expression = expression.subs(A_SYMBOL, reduce_skew.alpha)
        expression = expression.subs(THETA_SYMBOL, reduce_skew.theta)
        expression = expression.subs(TAU_SYMBOL, reduce_skew.tau)
        expression = expression.subs(TAU_SYMBOL_NON_COMM, reduce_skew.tau)

    return reduce_skew.reduce_skew(expression, p, n)

def convert_parenthesized_int_to_q(expr):

    new_expr = expr.replace(
        ParenthesizedInteger, 
        # lambda arg: sympy.UnevaluatedExpr(arg)
        lambda arg: arg
    )

    return new_expr

def convert_twist_symbol_to_q(expr):
    clean_expr = expr.replace(
        lambda x: isinstance(x, TwistSymbol),
        lambda x: x.base ** x.index
    )

    return clean_expr

def ore_matrix_to_R_matrix(ore_matrix: OreMatrix, phi_psi_options: PhiPsiOptions):
    raw: sympy.Matrix = None

    if isinstance(ore_matrix, list) and len(ore_matrix) > 0:
        raw = sp.Matrix(ore_matrix)
    else:
        raw = ore_matrix.display()

    raw = raw.replace(TAU_SYMBOL, TAU_SYMBOL_NON_COMM)
    # raw = toggle_commutativity_for_all_symbols_in_matrix(raw, change_to_commutative=True, keep_tau_non_comm=True)
    raw = raw.applyfunc(tmc_utils.move_tau_to_the_right_without_twist)

    leaf_count = 0
    byte_count = 0

    if phi_psi_options.calc_leaf_and_byte_counts == True:
        leaf_count = tmc_utils.mathematica_leaf_count(raw)
        byte_count = tmc_utils.mathematica_matrix_byte_count(raw)

    total_number_of_reductions_made = -1
    elapsed_time = -1.0

    return R_Matrix(raw, [], [], [], [], [], [], [], leaf_count, byte_count, total_number_of_reductions_made, elapsed_time)


def compare_matrices(name, matrix_a: sympy.Matrix, matrix_b: sympy.Matrix, phi_psi_options: PhiPsiOptions) -> bool:
    
    if matrix_a.shape != matrix_b.shape:
        print(f"Size mismatch: Matrix A has '{matrix_a.shape}', Matrix B has '{matrix_b.shape}'")
        return False

    all_matrix_entries_match = True
    matrix_a = toggle_commutativity_for_all_symbols_in_matrix(matrix_a, change_to_commutative=True, keep_tau_non_comm=False)

    if phi_psi_options.execute_symbolically == False:
        matrix_b = toggle_commutativity_for_all_symbols_in_matrix(matrix_b, change_to_commutative=True, keep_tau_non_comm=False)
        matrix_a = matrix_a.subs(A_SYMBOL, ALPHA_SYMBOL)


    matrix_a = matrix_a.applyfunc(convert_parenthesized_int_to_q)
    matrix_b = matrix_b.applyfunc(convert_twist_symbol_to_q)

    for row in range(matrix_a.rows):
        for col in range(matrix_a.cols):
            val_a = matrix_a[row, col]
            val_b = matrix_b[row, col]

            if not val_a.equals(val_b):
                all_matrix_entries_match = False
                position = (row+1, col+1)
                print(f"{name} mismatch, position ({position}): Matrix A has '{val_a}', Matrix B has '{val_b}'")

    return all_matrix_entries_match


def convert_2d_list_to_matrix(list_2d):
    
    aliases = {"theta": THETA_SYMBOL, "tau": TAU_SYMBOL, "θ": THETA_SYMBOL, "τ": TAU_SYMBOL, "FJTwist": FJTwist}

    return sp.Matrix([[sp.sympify(entry, locals=aliases) for entry in row] for row in list_2d])

def run_reduction(phi_psi_input: PhiPsiInput, phi_psi_options: PhiPsiOptions):

    phi_module: T_Module = phi_psi_input.Phi_module
    psi_module: T_Module = phi_psi_input.Psi_module

    Phi_input = phi_module.as_ai_input
    Psi_input = psi_module.as_ai_input
    
    response_dict = tmc_view.phi_psi_input_to_dict(phi_psi_input)

    # print("AI REDUCTION")
    # print(f"{Phi_input}")
    # print(f"{Psi_input}")

    parameters = None

    if phi_psi_options.execute_symbolically == True:
        parameters = FJSymbolicParameters(theta=THETA_SYMBOL, tau=TAU_SYMBOL)
    else:
        p = phi_psi_options.prime_characteristic
        n = phi_psi_options.extension_degree

        Phi_input = Phi_input.subs(A_SYMBOL, ALPHA_SYMBOL)
        Psi_input = Psi_input.subs(A_SYMBOL, ALPHA_SYMBOL)

        parameters = FJParametersPN(p, n, theta=THETA_SYMBOL, tau=TAU_SYMBOL, alpha=ALPHA_SYMBOL)

    data = ReductionData(Phi_input, Psi_input, parameters)

    report = FJReport(data)
    large_sequence = FJExactSequence(data)

    # print(f"DATA")
    # pprint.pp(vars(data))
    # print(f"REPORT")
    # print(f"{report}")
    

    Ext_OreMatrix: OreMatrix = report["Pi"] if phi_psi_options.execute_symbolically == True else large_sequence["Pi"]
    Ext0_OreMatrix: OreMatrix = large_sequence["Pi0"]
    K_OreMatrix: OreMatrix = large_sequence["K"]
    R_OreMatrix: OreMatrix = large_sequence["R"]
    P_OreMatrix: OreMatrix = large_sequence["P"]
    S_OreMatrix: OreMatrix = large_sequence["S"]

    Ext_Matrix = ore_matrix_to_R_matrix(Ext_OreMatrix, phi_psi_options)

    Ext0_Matrix = None
    K_Matrix = None
    R_MatrixInstance = None
    P_Matrix = None
    S_Matrix = None

    if phi_psi_options.only_generate_ext == False:
        Ext0_Matrix = ore_matrix_to_R_matrix(Ext0_OreMatrix, phi_psi_options)
        K_Matrix = ore_matrix_to_R_matrix(K_OreMatrix, phi_psi_options)
        R_MatrixInstance = ore_matrix_to_R_matrix(R_OreMatrix, phi_psi_options)
        P_Matrix = ore_matrix_to_R_matrix(P_OreMatrix, phi_psi_options)
        S_Matrix = ore_matrix_to_R_matrix(S_OreMatrix, phi_psi_options)

    response_dict["Ext_matrix"] = tmc_view.r_matrix_to_dict(Ext_Matrix)
    response_dict["Ext0_matrix"] = tmc_view.r_matrix_to_dict(Ext0_Matrix)
    response_dict["K_matrix"] = tmc_view.r_matrix_to_dict(K_Matrix)
    response_dict["R_matrix"] = tmc_view.r_matrix_to_dict(R_MatrixInstance)
    response_dict["P_matrix"] = tmc_view.r_matrix_to_dict(P_Matrix)
    response_dict["S_matrix"] = tmc_view.r_matrix_to_dict(S_Matrix)

    while True:
        tmc_utils.check_connection_still_alive()


    rightmost_basis: RightmostBasis = tmc_phi_psi.get_rightmost_basis(phi_module, psi_module)
    initial_analysis: InitialAnalysis = InitialAnalysis(rightmost_basis, None, -1)

    response_dict = response_dict | tmc_view.initial_analysis_to_dict(initial_analysis, phi_module, psi_module)

    total_reduction_count = 0
    total_leaf_count = 0
    total_byte_count = 0
    if phi_psi_options.only_generate_ext == False:
        total_reduction_count = Ext_Matrix.total_reduction_count + Ext0_Matrix.total_reduction_count + K_Matrix.total_reduction_count
        total_leaf_count = Ext_Matrix.leaf_count + Ext0_Matrix.leaf_count + K_Matrix.leaf_count
        total_byte_count = Ext_Matrix.byte_count + Ext0_Matrix.byte_count + K_Matrix.byte_count
    else:
        total_reduction_count = Ext_Matrix.total_reduction_count
        total_leaf_count = Ext_Matrix.leaf_count
        total_byte_count = Ext_Matrix.byte_count

    response_dict["total_reduction_count"] = total_reduction_count
    response_dict["total_leaf_count"] = total_leaf_count
    response_dict["total_byte_count"] = total_byte_count

    return PhiPsiReductionOutputFromAiCode(Ext_Matrix, Ext0_Matrix, K_Matrix, R_Matrix, P_Matrix, S_Matrix, response_dict)

def rightmost_compatible_basis(t_module_matrices):
    combined_matrix = sympy.Matrix.hstack(*t_module_matrices)

    combined_matrix_as_list = combined_matrix.tolist()

    try:
        start_time = int(time.time_ns())
        
        parameters = FJSymbolicParameters(theta=THETA_SYMBOL, tau=TAU_SYMBOL)

        result = compatible_basis_report(combined_matrix_as_list, parameters.field)

        end_time = int(time.time_ns())

        if result['Found'] == True:
            basis_matrix_as_list = result['BasisMatrix']

            for i in range(len(basis_matrix_as_list)):
                for k in range(len(basis_matrix_as_list[i])):
                    basis_matrix_as_list[i][k] = basis_matrix_as_list[i][k].as_expr()

            # print(f"{type(basis_matrix_as_list)} {basis_matrix_as_list}")

            basis_matrix = sympy.Matrix(basis_matrix_as_list)
            
            basis_matrix_inverse = basis_matrix.inv()
            basis_matrix_inverse_non_comm = basis_matrix_inverse.applyfunc(lambda matrix_elem: toggle_commutativity_for_all_symbols(matrix_elem, change_to_commutative=False))


            t_module_matrix_subset = []
            column_lists = []
            sorted_column_lists = []
            flag = []
            basis_subset = None
            last_m = -1
            total_seconds_to_generate = (end_time-start_time) / 1e9

            to_return = RightmostBasis(t_module_matrix_subset, column_lists, sorted_column_lists, flag, basis_subset, basis_matrix, basis_matrix_inverse, basis_matrix_inverse_non_comm, last_m, total_seconds_to_generate)

            return to_return
        else:
            return None
    except Exception as e:
        # traceback.print_exc()
        print(f"Could not find compatible basis due to exception: {e}")
        return None

