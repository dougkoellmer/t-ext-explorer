

import sympy

from tmc_consts import *
from tmc_utils import *
import tmc_utils
import copy
import tmc_phi_psi
import tmc_ai_bridge

import itertools

import random

def roll_dice(sides, start = 1):
    return random.randint(start, sides)

def validate_simplification(expression, expected_result_expression):
    simplified_expression = tmc_utils.simplify_each_term_in_expression(expression)

    latex_to_return = ""

    if simplified_expression != expected_result_expression:
        print(f"TEST_ERROR: original={expression}   simplified={simplified_expression}   expected={expected_result_expression}")
        latex_to_return = "\\textcolor{red}{" + sympy.latex(expression) + " \\neq " + sympy.latex(expected_result_expression) + "}";
    else:
        print(f"TEST_SUCCESS: original={expression}   simplified={simplified_expression}.  expected={expected_result_expression}")
        latex_to_return = sympy.latex(expression) + " = " + sympy.latex(simplified_expression)

    return latex_to_return

def validate_tau_swap(expression, expected_result_expression):
    swapped_expression = move_tau_to_the_right_of_other_symbols(expression)
    swapped_expression = toggle_commutativity_for_all_symbols(swapped_expression, change_to_commutative=True)

    latex_to_return = ""

    if swapped_expression != expected_result_expression:
        print(f"TEST_ERROR: original={expression}   swapped={swapped_expression}   expected={expected_result_expression}")
        latex_to_return = "\\textcolor{red}{" + sympy.latex(expression) + " \\neq " + sympy.latex(expected_result_expression) + "}";
    else:
        print(f"TEST_SUCCESS: original={expression}   swapped={swapped_expression}.  expected={expected_result_expression}")
        latex_to_return = sympy.latex(expression) + " = " + sympy.latex(swapped_expression)

    return latex_to_return

def binary_combinations_matrices(rows, cols) -> list[sympy.Matrix]:
    total_elements = rows * cols
    
    # flat_combinations = itertools.product([0, A_SYMBOL], repeat=total_elements)
    flat_combinations = itertools.product([0, 1], repeat=total_elements)
    
    matrices = []
    for comb in flat_combinations:
        M = sympy.Matrix(rows, cols, comb)
        matrices.append(M)
        
    return matrices

def split_matrix_into_blocks(M, block_rows, block_cols):
    total_rows, total_cols = M.shape
    
    # Check if the block dimensions perfectly divide the matrix dimensions
    if total_rows % block_rows != 0 or total_cols % block_cols != 0:
        raise ValueError("Matrix dimensions are not perfectly divisible by the block dimensions.")
        
    num_row_blocks = total_rows // block_rows
    num_col_blocks = total_cols // block_cols
    
    blocks = []
    for r in range(num_row_blocks):
        for c in range(num_col_blocks):
            row_start = r * block_rows
            row_end = (r + 1) * block_rows
            col_start = c * block_cols
            col_end = (c + 1) * block_cols
            
            submatrix = M[row_start:row_end, col_start:col_end]
            blocks.append(submatrix)
            
    return blocks

def generate_binary_matrices(n: int, d: int) -> list[list[sympy.Matrix]]:
    block_matrices = binary_combinations_matrices(d, n*d)
    to_return = []
    for ith_block_matrix in block_matrices:
        if ith_block_matrix.is_zero_matrix:
            pass
        else:
            ith_block_matrix_split = split_matrix_into_blocks(ith_block_matrix, d, d)
            if len(ith_block_matrix_split) > 1:
                all_the_rest_are_zeros = True
                for i in range(1, len(ith_block_matrix_split)):
                    ith_sub_block_matrix = ith_block_matrix_split[i]
                    if ith_sub_block_matrix.is_zero_matrix:
                        pass
                    else:
                        all_the_rest_are_zeros = False
                        break
                
                if all_the_rest_are_zeros == False:
                    to_return.append(ith_block_matrix_split)
            else:
                to_return.append(ith_block_matrix_split)
    
    return to_return


def multiply_matrices_by_tau(matrices: list[sympy.Matrix]):
    to_return = []
    for i in range(len(matrices)):
        ith_matrix = matrices[i]
        if i == 0:
            to_return.append(ith_matrix)
        else:
            tau_to_the_pow = TAU_SYMBOL_NON_COMM**i
            ith_matrix_times_tau = ith_matrix * tau_to_the_pow
            to_return.append(ith_matrix_times_tau)
        
    return to_return


def add_thetas_to_diagonal(matrix: sympy.Matrix, force = False):
    for m in range(matrix.rows):
        if matrix[m, m] == 0 or force == True:
            matrix[m, m] = THETA_SYMBOL
        else:
            matrix[m, m] = matrix[m, m] * THETA_SYMBOL

def fuzz_test():
    sympy.init_printing(use_unicode=True)
    sympy.init_printing(order='none')

    show_reduction_steps = False
    replace_c_with_tau = True
    replace_q_with_integer = True
    group_coefficients = True
    strip_integer_exponents_for_symbolic_execution = True
    execute_symbolically = True
    prime_characteristic = 3
    extension_degree = 2
    reduce_log_output = True
    only_generate_ext = True
    reduce_latex_output = True
    calc_leaf_and_byte_counts = False

    phi_psi_options = PhiPsiOptions(
        show_reduction_steps,
        replace_c_with_tau,
        replace_q_with_integer,
        group_coefficients,
        strip_integer_exponents_for_symbolic_execution,
        execute_symbolically,
        prime_characteristic,
        extension_degree,
        reduce_log_output,
        only_generate_ext,
        reduce_latex_output,
        calc_leaf_and_byte_counts
    )

    
    fuzz_test_random(phi_psi_options)
    # fuzz_test_dimensions(phi_psi_options)

def print_fuzz_test_error(phi_expr, psi_expr, error_message):
    print(f"FUZZ TEST ERROR: {error_message}")
    sympy.pprint(phi_expr)
    sympy.pprint(psi_expr)

def set_matrix_to_random_values(matrix):
    for i in range(matrix.rows):
        for j in range(matrix.cols):
            rand_int = random.randint(0, 3)
            if rand_int == 0:
                matrix[i, j] = 0
            elif rand_int == 1:
                matrix[i, j] = 1
            elif rand_int == 2:
                matrix[i, j] = THETA_SYMBOL
            elif rand_int == 3:
                matrix[i, j] = A_SYMBOL
            

def generate_random_t_module(matrix_count, d, phi_psi_options: PhiPsiOptions):

    matrices = []
    for i in range(matrix_count):
        ith_matrix = sympy.zeros(d, d)
        
        set_matrix_to_random_values(ith_matrix)

        if i == 0:
            add_thetas_to_diagonal(ith_matrix, True)
        
        matrices.append(ith_matrix)

    matrices = multiply_matrices_by_tau(matrices)
    phi_expr = sympy.MatAdd(*matrices, evaluate=False)

    return phi_expr


def generate_random_phi_psi_pair(phi_psi_options: PhiPsiOptions) -> PhiPsiInput:
    psi_matrix_count = roll_dice(3, 2)
    phi_matrix_count = roll_dice(4, psi_matrix_count+1)

    d_phi = roll_dice(2)
    d_psi = roll_dice(2)

    phi_expr = generate_random_t_module(phi_matrix_count, d_phi, phi_psi_options)
    psi_expr = generate_random_t_module(psi_matrix_count, d_psi, phi_psi_options)

    phi_psi_input = tmc_phi_psi.validate_input(phi_expr, psi_expr, phi_psi_options)

    return phi_psi_input


def fuzz_test_random(phi_psi_options: PhiPsiOptions):
    limit = 10000

    attempt_count = 0
    invalid_input_count = 0
    basis_count = 0
    error_count = 0
    match_count = 0
    mismatch_count = 0

    def fuzz_test_random_error(msg, phi_psi_input: PhiPsiInput):
        print("FUZZ ERROR: " + msg)
        tmc_view.log_input_in_various_ways(phi_psi_input.Phi_module.simplified, phi_psi_input.Psi_module.simplified)
    
    for i in range(limit):
        try:
            attempt_count += 1
            # print_fuzz_test_error(phi_expr, psi_expr, f"Exception: {e}")
            phi_psi_input = generate_random_phi_psi_pair(phi_psi_options)
            phi_module = phi_psi_input.Phi_module
            psi_module = phi_psi_input.Psi_module

            phi_module_for_log = sympy.MatAdd(*phi_module.expanded_matrices, evaluate=False)
            psi_module_for_log = sympy.MatAdd(*psi_module.expanded_matrices, evaluate=False)

            print(f"attempt_count={attempt_count-1}, invalid_input_count={invalid_input_count}, basis_count={basis_count}, error_count={error_count}, match_count={match_count}, mismatch_count={mismatch_count}")
            sympy.pprint(phi_module_for_log)
            sympy.pprint(psi_module_for_log)
        
            if len(phi_psi_input.input_errors) == 0:
                only_do_ai_reduction = True
                rightmost_basis = tmc_phi_psi.get_rightmost_basis(phi_module, psi_module)

                rightmost_basis_expr_for_log = sympy.MatAdd(*phi_module.expanded_matrices, evaluate=False)

                if rightmost_basis is not None:
                    basis_count += 1

                    tmc_view.log_input_in_various_ways(phi_module.simplified, psi_module.simplified)
                
                    if only_do_ai_reduction == True:

                        ai_reduction_output: PhiPsiReductionOutputFromAiCode = tmc_ai_bridge.run_reduction(phi_module, psi_module, phi_psi_options)

                        if ai_reduction_output.Ext_Matrix is not None and ai_reduction_output.Ext_Matrix.raw.rows > 0:
                            print("REDUCED")
                            pass
                        else:
                            unexpected_error_count = unexpected_error_count + 1
                            fuzz_test_random_error("ai_reduction_output.Ext_matrix came back None/empty", phi_psi_input)
                    else:
                        phi_psi_reduction_output: PhiPsiReductionOutput = tmc_phi_psi.run_analysis(phi_psi_input, phi_psi_options, rightmost_basis)

                        if phi_psi_reduction_output is not None:
                            if phi_psi_reduction_output.found_rightmost_basis == True:

                                if phi_psi_reduction_output.Ext_matrix is not None:
                                    # phi_psi_reduction_output.Ext_matrix.raw[0, 0] = -2
                                    ai_reduction_output: PhiPsiReductionOutputFromAiCode = tmc_ai_bridge.run_reduction(phi_psi_input.Phi_module, phi_psi_input.Psi_module, phi_psi_options)
                                    match = tmc_ai_bridge.compare_matrices("W_Matrix", phi_psi_reduction_output.Ext_matrix.raw, ai_reduction_output.Ext_Matrix.raw, phi_psi_options)

                                    if match == True:
                                        match_count += 1
                                        pass
                                    else:
                                        mismatch_count += 1
                                        fuzz_test_random_error("Original code and AI code don't agree.", phi_psi_input)
                                else:
                                    error_count += 1
                                    fuzz_test_random_error("phi_psi_reduction_output.Ext_matrix came back None", phi_psi_input)
                            else:
                                # No rightmost basis, but we already checked above anyway so should never get here, treat as error.
                                error_count += 1
                                fuzz_test_random_error("Should have already determined that there isn't a rightmost basis.", phi_psi_input)
                        else:
                            error_count += 1
                            fuzz_test_random_error("phi_psi_reduction_output came back None", phi_psi_input)
                else:
                    # print("NO RIGHTMOST BASIS")
                    # sympy.pprint(rightmost_basis_expr_for_log)
                    pass # Expected, no rightmost basis
            else:
                if phi_module.has_critical_error == False and psi_module.has_critical_error == False:
                    invalid_input_count += 1
                    pass
                else:
                    error_count += 1
                    fuzz_test_random_error(f"Strange input: f{phi_psi_input.input_errors}", phi_psi_input)
        except Exception as e:
            error_count += 1
            fuzz_test_random_error(f"Exception {e}", phi_psi_input)

def fuzz_test_dimensions(phi_psi_options: PhiPsiOptions):

    MIN_MATRICES_PHI = 3
    MIN_MATRICES_PSI = 2
    MAX_MATRICES_PHI = 5
    MAX_MATRICES_PSI = 4

    MIN_MATRICES_PHI = 3
    MIN_MATRICES_PSI = 2
    MAX_MATRICES_PHI = 3
    MAX_MATRICES_PSI = 2

    total_count = 0
    total_that_had_rightmost_basis = 0
    total_mismatches_found = 0
    unexpected_error_count = 0
    skip_count = 0

    start_time = int(time.time_ns())

    MIN_DIMENSION = 1
    MAX_DIMENSION = 2

    for phi_dimension in range(MIN_DIMENSION, MAX_DIMENSION+1):
        for psi_dimension in range(MIN_DIMENSION, MAX_DIMENSION+1):
            for i in range(MIN_MATRICES_PHI, MAX_MATRICES_PHI+1):
                number_of_phi_terms = i
                print(f"FUZZ_LOG Generating Phi combinations for phi_max_degree={i}, phi_d={phi_dimension}, psi_d={psi_dimension}")
                ith_phi_equations = generate_binary_matrices(number_of_phi_terms, phi_dimension)
                # print(ith_phi_equations)

                all_psi_modules = []

                for z in range(len(ith_phi_equations)):
                    ith_phi_equation = ith_phi_equations[z]
                    add_thetas_to_diagonal(ith_phi_equation[0])
                    ith_phi_equation = multiply_matrices_by_tau(ith_phi_equation)
                    current_time = int(time.time_ns())
                    total_seconds_so_far = (current_time-start_time) / 1e9
                    phi_expr = sympy.MatAdd(*ith_phi_equation, evaluate=False)
                    phi_module: T_Module = tmc_view.parse_t_module_latex(phi_expr, "Phi", phi_psi_options)

                    print(f"FUZZ_LOG Checking for Phi combination {z+1} of {len(ith_phi_equations)}, seconds_running={total_seconds_so_far}, total_count={total_count}, total_errors={unexpected_error_count}, skip_count={skip_count}, total_with_basis={total_that_had_rightmost_basis}")
                    sympy.pprint(phi_expr)

                    position_in_psi_modules_array = 0
                    
                    for k in range(MIN_MATRICES_PSI, MAX_MATRICES_PSI+1):
                        number_of_psi_terms = min(k, number_of_phi_terms-1)

                        ith_psi_equations = None
                        have_to_generate_modules = False
                        generated_psi_modules = []

                        if position_in_psi_modules_array >= len(all_psi_modules):
                            ith_psi_equations = generate_binary_matrices(number_of_psi_terms, psi_dimension)
                            have_to_generate_modules = True
                        else:
                            ith_psi_equations = all_psi_modules[position_in_psi_modules_array]

                        print(f"ith_psi_equations={len(ith_psi_equations)}")
                        
                        position_in_psi_modules_array = position_in_psi_modules_array + 1

                        for y in range(len(ith_psi_equations)):

                            dice_roll = 1
                            if psi_dimension == 3:
                                dice_roll = 1000
                            elif psi_dimension == 2:
                                dice_roll = 20

                            if roll_dice(dice_roll) != 1:
                                skip_count = skip_count + 1
                                continue
                            
                            ith_psi_equation = ith_psi_equations[y]

                            psi_module: T_Module = None
                            if have_to_generate_modules == True:
                                add_thetas_to_diagonal(ith_psi_equation[0])
                                ith_psi_equation = multiply_matrices_by_tau(ith_psi_equation)
                                psi_expr = sympy.MatAdd(*ith_psi_equation, evaluate=True)
                                psi_module = tmc_view.parse_t_module_latex(psi_expr, "Psi", phi_psi_options)
                                generated_psi_modules.append(psi_module)
                                # print("generating_psi")
                            else:
                                psi_module = ith_psi_equation
                                # print("using_cached_psi")

                            psi_expr = psi_module.simplified
                            
                            # print("COMPARING:")
                            # sympy.pprint(phi_expr)
                            # sympy.pprint(psi_expr)

                            try:

                                phi_psi_input: PhiPsiInput = tmc_phi_psi.validate_input(phi_module, psi_module, phi_psi_options)

                                if len(phi_psi_input.input_errors) > 0:
                                    if phi_module.has_critical_error == False and psi_module.has_critical_error == False:
                                        # print(f"INPUT ERROR: f{phi_psi_input.input_errors}")
                                        # sympy.pprint(phi_expr)
                                        # sympy.pprint(psi_expr)

                                        # Expected once in a while from the way we're generating phi/psi combinations.
                                        # E.g. First matrix may not be nilpotent, or all zeros in one matrix may make its effective degree not high enough, etc.
                                        pass
                                    else:
                                        print("OH NO")
                                        sympy.pprint(phi_expr)
                                        sympy.pprint(psi_expr)
                                        unexpected_error_count = unexpected_error_count + 1
                                        print_fuzz_test_error(phi_expr, psi_expr, f"Input error(s): {phi_psi_input.input_errors}")
                                else:
                                    total_count = total_count + 1

                                    only_do_ai_reduction = False
                                    rightmost_basis = tmc_phi_psi.get_rightmost_basis(phi_module, psi_module)

                                    rightmost_basis_expr_for_log = sympy.MatAdd(*phi_module.expanded_matrices, evaluate=False)

                                    if rightmost_basis is not None:
                                        total_that_had_rightmost_basis = total_that_had_rightmost_basis + 1
                                    
                                        if only_do_ai_reduction == True:

                                            ai_reduction_output: PhiPsiReductionOutputFromAiCode = tmc_ai_bridge.run_reduction(phi_module, psi_module, phi_psi_options)

                                            if ai_reduction_output.Ext_Matrix is not None and ai_reduction_output.Ext_Matrix.raw.rows > 0:
                                                pass
                                            else:
                                                unexpected_error_count = unexpected_error_count + 1
                                                print_fuzz_test_error(phi_expr, psi_expr, "ai_reduction_output.Ext_matrix came back None/empty")
                                        else:
                                            phi_psi_reduction_output: PhiPsiReductionOutput = tmc_phi_psi.run_analysis(phi_psi_input, phi_psi_options, rightmost_basis)

                                            if phi_psi_reduction_output is not None:
                                                if phi_psi_reduction_output.found_rightmost_basis == True:

                                                    if phi_psi_reduction_output.Ext_matrix is not None:
                                                        # phi_psi_reduction_output.Ext_matrix.raw[0, 0] = -2
                                                        ai_reduction_output: PhiPsiReductionOutputFromAiCode = tmc_ai_bridge.run_reduction(phi_psi_input.Phi_module, phi_psi_input.Psi_module, phi_psi_options)
                                                        match = tmc_ai_bridge.compare_matrices("W_Matrix", phi_psi_reduction_output.Ext_matrix.raw, ai_reduction_output.Ext_Matrix.raw, phi_psi_options)

                                                        if match == True:
                                                            # print("MATCH!")
                                                            pass
                                                        else:
                                                            total_mismatches_found = total_mismatches_found + 1
                                                            print_fuzz_test_error(phi_expr, psi_expr, "Original code and AI code don't agree.")
                                                    else:
                                                        unexpected_error_count = unexpected_error_count + 1
                                                        print_fuzz_test_error(phi_expr, psi_expr, "phi_psi_reduction_output.Ext_matrix came back None")
                                                else:
                                                    # No rightmost basis, but we already checked above anyway so should never get here, treat as error.
                                                    unexpected_error_count = unexpected_error_count + 1
                                                    print_fuzz_test_error(phi_expr, psi_expr, "Should have already determined that there isn't a rightmost basis.")
                                            else:
                                                unexpected_error_count = unexpected_error_count + 1
                                                print_fuzz_test_error(phi_expr, psi_expr, "phi_psi_reduction_output came back None")
                                    else:
                                        # print("NO RIGHTMOST BASIS")
                                        # sympy.pprint(rightmost_basis_expr_for_log)
                                        pass # Expected, no rightmost basis
                            except Exception as e:
                                unexpected_error_count = unexpected_error_count + 1
                                print_fuzz_test_error(phi_expr, psi_expr, f"Exception: {e}")

                        if have_to_generate_modules == True:
                            all_psi_modules.append(generated_psi_modules)

    print(f"TOTAL VALID INPUTS: {total_count}")
    print(f"INPUTS WITH BASIS: {total_that_had_rightmost_basis}")
    print(f"MISMATCHING OUTPUT COUNT: {total_mismatches_found}")
    print(f"UNEXPECTED ERROR COUNT: {unexpected_error_count}")
    


def compare_output(Phi_latex, Psi_latex):
    phi_psi_input = tmc_phi_psi.validate_input(Phi_latex, Psi_latex, phi_psi_options)

    phi_psi_reduction_output: PhiPsiReductionOutput = tmc_phi_psi.run_analysis(phi_psi_input, phi_psi_options)

    
    print(f"STARTING MATRIX COMPARISON")
    ai_reduction: PhiPsiReductionOutputFromAiCode = tmc_ai_bridge.run_reduction(phi_psi_input.Phi_module, phi_psi_input.Psi_module, phi_psi_options)
    tmc_ai_bridge.compare_matrices("Ext_Matrix", phi_psi_reduction_output.Ext_matrix.raw, ai_reduction.Ext_Matrix.raw, phi_psi_options)
    # tmc_ai_bridge.compare_matrices("V_Matrix", phi_psi_reduction_output.V_matrix.raw, ai_reduction.V_matrix)
    # tmc_ai_bridge.compare_matrices("W_V_Matrix", phi_psi_reduction_output.W_V_matrix.raw, ai_reduction.W_V_matrix)
    print(f"ENDING MATRIX COMPARISON")

def run_tests():
    latex_array = []

    if FEATURE_FLAG__RUN_UNIT_TESTS == False:
        return latex_array

    fuzz_test()

    return

def run_twist_tests():

    p_n_test_expr = C_SYMBOL**Q_SYMBOL**4
    p_n_test_expr_after = tmc_utils.simplify_twisted_modular(p_n_test_expr, 3, 2)
    print(f"simplify_twisted_modular: {p_n_test_expr}, {p_n_test_expr_after}")

    p_n_test_expr = C_SYMBOL**Q_SYMBOL**2
    p_n_test_expr_after = tmc_utils.simplify_twisted_modular(p_n_test_expr, 3, 2)
    print(f"simplify_twisted_modular: {p_n_test_expr}, {p_n_test_expr_after}")

    p_n_test_expr = C_SYMBOL**Q_SYMBOL**3
    p_n_test_expr_after = tmc_utils.simplify_twisted_modular(p_n_test_expr, 3, 2)
    print(f"simplify_twisted_modular: {p_n_test_expr}, {p_n_test_expr_after}")

    p_n_test_expr = C_SYMBOL**(Q_SYMBOL+2)
    p_n_test_expr_after = tmc_utils.simplify_twisted_modular(p_n_test_expr, 3, 2)
    print(f"simplify_twisted_modular: {p_n_test_expr}, {p_n_test_expr_after}")

    p_n_test_expr = -C_SYMBOL*TAU_SYMBOL**2
    p_n_test_expr_after = tmc_utils.simplify_twisted_modular(p_n_test_expr, 3, 2)
    print(f"simplify_twisted_modular: {p_n_test_expr}, {p_n_test_expr_after}")

    c_replacement_test = (C_SYMBOL**Q_SYMBOL**2)**Q_SYMBOL**2
    c_replacement_test_simple = sympy.simplify(c_replacement_test)
    print(f"c_replacement_test={c_replacement_test}, c_replacement_test_simple={c_replacement_test_simple}")

    if TAU_SYMBOL_NON_COMM == sympy.symbols("tau", commutative=False):
        print("TAU EQUAL FOR COMM=False (Expected)")
    else:
        print("TAU NOT EQUAL FOR COMM=False (Not Expected)")

    if TAU_SYMBOL_NON_COMM == sympy.symbols("tau", commutative=True):
        print("TAU EQUAL FOR COMM=True COMPARISON (Not Expected)")
    else:
        print("TAU NOT EQUAL FOR COMM=True COMPARISON (Expected)")

    # Some local shorthand to make tests easier to read. Introduced after experience with tests below which are hard to read.
    b = B_SYMBOL_FOR_TESTS
    c = C_SYMBOL
    q = Q_SYMBOL
    tau = TAU_SYMBOL_NON_COMM


    print(f"{9 == sympy.sympify('9')}")

    latex_array.append(validate_simplification(B_SYMBOL_FOR_TESTS**(-Q_SYMBOL-1)*B_SYMBOL_FOR_TESTS**(Q_SYMBOL+1), 1))
    latex_array.append(validate_simplification(B_SYMBOL_FOR_TESTS**(-Q_SYMBOL-1)*B_SYMBOL_FOR_TESTS**(Q_SYMBOL+1)*TAU_SYMBOL_NON_COMM**2, TAU_SYMBOL_NON_COMM**2))
    # latex_array.append(validate_simplification(b*b**(-q - 1)*c**(q**2)*tau, TAU_SYMBOL_NON_COMM))

    latex_array.append(validate_tau_swap( TAU_SYMBOL_NON_COMM*(A_SYMBOL_NON_COMM**-1)*C_SYMBOL_NON_COMM, A_SYMBOL**(-1*Q_SYMBOL)*C_SYMBOL**Q_SYMBOL*TAU_SYMBOL_NON_COMM ))

    latex_array.append(validate_tau_swap( TAU_SYMBOL_NON_COMM * non_comm_number(9), 9**Q_SYMBOL*TAU_SYMBOL_NON_COMM ))
    latex_array.append(validate_tau_swap( C_SYMBOL_NON_COMM*TAU_SYMBOL_NON_COMM, C_SYMBOL*TAU_SYMBOL_NON_COMM ))
    latex_array.append(validate_tau_swap( TAU_SYMBOL_NON_COMM*C_SYMBOL_NON_COMM, C_SYMBOL**Q_SYMBOL*TAU_SYMBOL_NON_COMM ))
    latex_array.append(validate_tau_swap( TAU_SYMBOL_NON_COMM**2*C_SYMBOL_NON_COMM, C_SYMBOL**Q_SYMBOL**2*TAU_SYMBOL_NON_COMM**2 ))
    latex_array.append(validate_tau_swap( TAU_SYMBOL_NON_COMM**2*C_SYMBOL_NON_COMM*TAU_SYMBOL_NON_COMM, C_SYMBOL**Q_SYMBOL**2*TAU_SYMBOL_NON_COMM**3 ))
    latex_array.append(validate_tau_swap( TAU_SYMBOL_NON_COMM**2*C_SYMBOL_NON_COMM**Q_SYMBOL**1, C_SYMBOL**Q_SYMBOL**3*TAU_SYMBOL_NON_COMM**2 ))
    latex_array.append(validate_tau_swap( TAU_SYMBOL_NON_COMM**2*C_SYMBOL_NON_COMM**Q_SYMBOL**1, C_SYMBOL**Q_SYMBOL**3*TAU_SYMBOL_NON_COMM**2 ))
    latex_array.append(validate_tau_swap( TAU_SYMBOL_NON_COMM**2*C_SYMBOL_NON_COMM**Q_SYMBOL**2, C_SYMBOL**Q_SYMBOL**4*TAU_SYMBOL_NON_COMM**2 ))
    
    latex_array.append(validate_tau_swap( TAU_SYMBOL_NON_COMM*C_SYMBOL_NON_COMM*TAU_SYMBOL_NON_COMM, C_SYMBOL**Q_SYMBOL*TAU_SYMBOL_NON_COMM**2 ))
    latex_array.append(validate_tau_swap( TAU_SYMBOL_NON_COMM*C_SYMBOL_NON_COMM*TAU_SYMBOL_NON_COMM**2, C_SYMBOL**Q_SYMBOL*TAU_SYMBOL_NON_COMM**3 ))

    latex_array.append(validate_tau_swap( (TAU_SYMBOL_NON_COMM*C_SYMBOL_NON_COMM*TAU_SYMBOL_NON_COMM**2)/2, (C_SYMBOL**Q_SYMBOL*TAU_SYMBOL_NON_COMM**3)/2 ))


    latex_array.append(validate_tau_swap( TAU_SYMBOL_NON_COMM**2*C_SYMBOL_NON_COMM*A_SYMBOL_NON_COMM, (C_SYMBOL**Q_SYMBOL**2)*(A_SYMBOL**Q_SYMBOL**2)*TAU_SYMBOL_NON_COMM**2 ))
    latex_array.append(validate_tau_swap( TAU_SYMBOL_NON_COMM**2*C_SYMBOL_NON_COMM*(A_SYMBOL_NON_COMM**Q_SYMBOL**2), (C_SYMBOL**Q_SYMBOL**2)*(A_SYMBOL**Q_SYMBOL**4)*TAU_SYMBOL_NON_COMM**2 ))


    # latex_array.append(validate_tau_swap( TAU_SYMBOL_NON_COMM*5*c_symbol, 5**q_symbol*c_symbol**q_symbol*TAU_SYMBOL_NON_COMM ))

    return latex_array