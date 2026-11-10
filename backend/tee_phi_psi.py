
import sympy
import math

import tee_view
import tee_utils
import tee_ai_bridge

import traceback

from tee_consts import *
from tee_types import *
from tee_utils import *
from tee_view import *

from copy import deepcopy

from sympy.polys.galoistools import gf_irreducible
from sympy.polys.domains import ZZ

import time

def validate_input(Phi_latex_or_expr_or_module: str, Psi_latex_or_expr_or_module: str, phi_psi_options: PhiPsiOptions) -> tee_types.PhiPsiInput:
    input_errors = []
    module_name_phi = "<math-span>\\Phi</math-span>"
    module_name_psi = "<math-span>\\Psi</math-span>"

    # print(f"validating1={Psi_latex_or_expr_or_module}")


    Phi_module: T_Module = tee_view.parse_t_module_latex(Phi_latex_or_expr_or_module, module_name_phi, phi_psi_options)
    Psi_module: T_Module = tee_view.parse_t_module_latex(Psi_latex_or_expr_or_module, module_name_psi, phi_psi_options)

    
    if Phi_module.max_deg > -1 and Phi_module.max_deg <= 1:
        if isinstance(Phi_latex_or_expr_or_module, tee_types.T_Module) == False:
            Phi_module.parse_errors.append(f"{module_name_phi} must have degree of at least 2")
        else:
            input_errors.append(f"{module_name_phi} must have degree of at least 2")
        
    
    if Phi_module.max_deg > -1 and Psi_module.max_deg == 0:
        if isinstance(Psi_latex_or_expr_or_module, tee_types.T_Module) == False:
            Psi_module.parse_errors.append(f"{module_name_psi} must have degree of at least 1")
        else:
            input_errors.append(f"{module_name_phi} must have degree of at least 2")

    input_errors += Phi_module.parse_errors
    input_errors += Psi_module.parse_errors

    if len(input_errors) > 0:
        is_max_degree_error = True
        return tee_types.PhiPsiInput(Phi_module, Psi_module, input_errors, is_max_degree_error)
    else:
        if Psi_module.max_deg >= Phi_module.max_deg:
            last_error = f"{module_name_phi} must be of greater degree than {module_name_psi}"
            if isinstance(Phi_latex_or_expr_or_module, tee_types.T_Module) == False:
                Phi_module.parse_errors.append(last_error)
            input_errors.append(last_error)

        return tee_types.PhiPsiInput(Phi_module, Psi_module, input_errors)


def run_analysis(phi_psi_input: PhiPsiInput, phi_psi_options: PhiPsiOptions, known_rightmost_basis: RightmostBasis = None):

    Phi_module: T_Module = phi_psi_input.Phi_module
    Psi_module: T_Module = phi_psi_input.Psi_module

    if phi_psi_options.reduce_latex_output == False:
        response_dict = tee_view.phi_psi_input_to_dict(phi_psi_input)
    else:
        response_dict = {}

    Ext_matrix: R_Matrix = None
    V_matrix: R_Matrix = None
    W_V_matrix: R_Matrix = None

    if len(phi_psi_input.input_errors) > 0:
        to_return = PhiPsiReductionOutput(Ext_matrix, V_matrix, W_V_matrix, response_dict, False)
        return to_return

    if phi_psi_options.reduce_log_output == False:
        tee_view.log_input_in_various_ways(Phi_module.simplified, Psi_module.simplified)

    try:
        initial_analysis: InitialAnalysis = run_initial_analysis(Phi_module, Psi_module, phi_psi_options, known_rightmost_basis)

        if initial_analysis is None:
            response_dict["error"] = "Could not find rightmost basis."
            to_return = PhiPsiReductionOutput(Ext_matrix, V_matrix, W_V_matrix, response_dict, False)
            return to_return

        Ext_matrix: R_Matrix = get_r_matrix(Phi_module, Psi_module, initial_analysis, PhiPsiSubvariant.ORIGINAL_LOGIC, phi_psi_options)

        total_reduction_count = 0
        total_leaf_count = 0
        total_byte_count = 0
        if phi_psi_options.only_generate_ext == False:
            V_matrix: R_Matrix = get_r_matrix(Phi_module, Psi_module, initial_analysis, PhiPsiSubvariant.CONSTANT_TERM_LOGIC, phi_psi_options)
            W_V_matrix: R_Matrix = get_r_matrix(Phi_module, Psi_module, initial_analysis, PhiPsiSubvariant.EMBEDDING_LOGIC, phi_psi_options, Ext_matrix.Eij_tau_k_matrices)
            total_reduction_count = Ext_matrix.total_reduction_count + V_matrix.total_reduction_count + W_V_matrix.total_reduction_count
            total_leaf_count = Ext_matrix.leaf_count + V_matrix.leaf_count + W_V_matrix.leaf_count
            total_byte_count = Ext_matrix.byte_count + V_matrix.byte_count + W_V_matrix.byte_count
        else:
            total_reduction_count = Ext_matrix.total_reduction_count
            total_leaf_count = Ext_matrix.leaf_count
            total_byte_count = Ext_matrix.byte_count

        if phi_psi_options.reduce_latex_output == False:
            response_dict = response_dict | tee_view.initial_analysis_to_dict(initial_analysis, Phi_module, Psi_module)

        if phi_psi_options.reduce_latex_output == False:
            Ext_matrix_dict = tee_view.r_matrix_to_dict(Ext_matrix)

            if phi_psi_options.only_generate_ext == False:
                V_matrix_dict = tee_view.r_matrix_to_dict(V_matrix)
                W_V_matrix_dict = tee_view.r_matrix_to_dict(W_V_matrix)
        else:
            Ext_matrix_dict = {}
            V_matrix_dict = {}
            W_V_matrix_dict = {}

        response_dict["Ext_matrix"] = Ext_matrix_dict

        if phi_psi_options.only_generate_ext == False:
            response_dict["V_matrix"] = V_matrix_dict
            response_dict["W_V_matrix"] = W_V_matrix_dict

        response_dict["total_reduction_count"] = total_reduction_count
        response_dict["total_leaf_count"] = total_leaf_count
        response_dict["total_byte_count"] = total_byte_count

        # print(f"response_dict={response_dict}")

    except Exception as e:
        traceback.print_exc()
        catch_all_exception_message = f"Error: {e}"

        response_dict["error"] = catch_all_exception_message

    to_return = PhiPsiReductionOutput(Ext_matrix, V_matrix, W_V_matrix, response_dict, True)

    return to_return

def get_rightmost_basis(Phi_module: T_Module, Psi_module: T_Module) -> RightmostBasis:
    number_of_Phi_matrices = len(Phi_module.expanded_matrices)
    number_of_Psi_matrices = len(Psi_module.expanded_matrices)
    matrices_to_use_for_determining_rightmost_basis = []

    for i in range(number_of_Phi_matrices):
        if i >= number_of_Psi_matrices:
            matrices_to_use_for_determining_rightmost_basis.append(Phi_module.expanded_matrices[i])

    rightmost_basis: RightmostBasis = tee_ai_bridge.rightmost_compatible_basis(matrices_to_use_for_determining_rightmost_basis)

    return rightmost_basis

def run_initial_analysis(Phi_module: T_Module, Psi_module: T_Module, phi_psi_options: PhiPsiOptions, known_rightmost_basis: RightmostBasis = None) -> InitialAnalysis:

    start_time = int(time.time_ns())

    rightmost_basis: RightmostBasis = known_rightmost_basis if known_rightmost_basis is not None else get_rightmost_basis(Phi_module, Psi_module)

    if rightmost_basis is not None:
        Eij_and_Aij_matrices: Eij_And_Aij_Matrices = get_Eij_and_Aij_matrices(Phi_module.simplified, Psi_module.simplified, rightmost_basis, phi_psi_options)

        end_time = int(time.time_ns())
        elapsed_time = (end_time-start_time) / 1e9

        initial_analysis: InitialAnalysis = InitialAnalysis(rightmost_basis, Eij_and_Aij_matrices, elapsed_time)

        return initial_analysis
    # else:
    #     raise ValueError("Could not find rightmost basis.")

    return None

def get_r_matrix(Phi_module: T_Module, Psi_module: T_Module, initial_analysis: InitialAnalysis, subvariant: PhiPsiSubvariant, phi_psi_options: PhiPsiOptions, basis_matrices_from_previous = None) -> R_Matrix :
    rightmost_basis = initial_analysis.rightmost_basis
    Eij_and_Aij_matrices = initial_analysis.Eij_and_Aij_matrices

    reduction_limit = FEATURE_FLAG__REDUCTION_LIMIT_EXT if subvariant == PhiPsiSubvariant.ORIGINAL_LOGIC else FEATURE_FLAG__REDUCTION_LIMIT_OTHER

    stripping_integers = False
    replacing_q_with_p_n = False
    further_process_q_k = False
    if phi_psi_options.strip_integer_exponents_for_symbolic_execution and phi_psi_options.execute_symbolically:
        further_process_q_k = True
        stripping_integers = True
    elif phi_psi_options.execute_symbolically == False and phi_psi_options.prime_characteristic >= 0 and phi_psi_options.extension_degree >= 0:
        further_process_q_k = True
        replacing_q_with_p_n = True

    p_val = phi_psi_options.prime_characteristic
    n_val = phi_psi_options.extension_degree
    
    # Used repeatedly below for debug output.
    U_inverse_matrix_latex = ""
    if phi_psi_options.reduce_latex_output == False:
        U_inverse_matrix_latex = matrix_to_latex(rightmost_basis.basis_matrix_inverse)
    
    degree_k_start_values = []
    degree_k_finish_values = []
    Eij_tau_k_matrices_to_return = []
    Psi_Eij_c_tau_k_matrices_to_return = []

    if basis_matrices_from_previous is not None:
        Eij_tau_k_matrices_to_return = basis_matrices_from_previous

    decompositions = []
    
    reduction_matrices = []
    final_reduction_matrices = []
    reduction_counts = []

    r_matrix_col_count = 0

    start_time = int(time.time_ns())

    total_number_of_reductions_made = 0

    i_times_k = 0

    for i in range(len(Eij_and_Aij_matrices.Eij_matrices)):

        check_connection_still_alive()

        ith_Aij_Phi_minus_Psi_Aij_matrix = Eij_and_Aij_matrices.Aij_Phi_minus_Psi_Aij_matrices[i]
        ith_Eij_matrix = Eij_and_Aij_matrices.Eij_matrices[i]
        ith_Nij = Eij_and_Aij_matrices.Nij_list[i]
        contains_constant_term = Eij_and_Aij_matrices.contains_constant_term_list[i]

        degree_k_start = 0 if subvariant == PhiPsiSubvariant.ORIGINAL_LOGIC else 1
        degree_k_finish = ith_Nij #EXCLUSIVE: so the final degree is ith_Nij-1 (at least as first defined, may be changed right below)

        using_constant_term_logic = (subvariant == PhiPsiSubvariant.CONSTANT_TERM_LOGIC or subvariant == PhiPsiSubvariant.EMBEDDING_LOGIC) and contains_constant_term

        # if subvariant == PhiPsiSubvariant.CONSTANT_TERM_LOGIC:
        #     print(f"i={i} contains_constant_term={contains_constant_term}")

        if using_constant_term_logic:
            degree_k_finish = ith_Nij+1 #AGAIN, EXCLUSIVE: so the final degree is ith_Nij here.
        
        degree_k_start_values.append(degree_k_start)
        degree_k_finish_values.append(degree_k_finish)

        r_matrix_col_count += degree_k_finish - degree_k_start

        kth_Eij_tau_k_matrices = Eij_and_Aij_matrices.Eij_tau_k_matrices_wrapper[i]
        kth_Eij_c_tau_k_matrices = Eij_and_Aij_matrices.Eij_c_tau_k_matrices_wrapper[i]
        Psi_Eij_c_tau_k_matrices = Eij_and_Aij_matrices.Psi_Eij_c_tau_k_matrices_wrapper[i];

        for degree_k in range(degree_k_start, degree_k_finish):

            check_connection_still_alive()

            kth_Eij_tau_k_with_effective_elem = kth_Eij_c_tau_k_matrices[degree_k] if subvariant == PhiPsiSubvariant.EMBEDDING_LOGIC else kth_Eij_tau_k_matrices[degree_k]
            kth_Psi_Eij_c_tau_k_with_effective_elem = Psi_Eij_c_tau_k_matrices[degree_k]

            kth_Eij_tau_k = kth_Eij_tau_k_with_effective_elem.raw
            ongoing_reduction_matrix_with_effective_elem = kth_Psi_Eij_c_tau_k_with_effective_elem if subvariant != PhiPsiSubvariant.EMBEDDING_LOGIC else kth_Eij_tau_k_with_effective_elem
            ongoing_reduction_matrix = ongoing_reduction_matrix_with_effective_elem.raw
            

            if basis_matrices_from_previous is None:
                Eij_tau_k_matrices_to_return.append(kth_Eij_tau_k_with_effective_elem)
            
            if subvariant == PhiPsiSubvariant.EMBEDDING_LOGIC:
                Psi_Eij_c_tau_k_matrices_to_return.append(kth_Eij_tau_k_with_effective_elem)
            else:
                Psi_Eij_c_tau_k_matrices_to_return.append(kth_Psi_Eij_c_tau_k_with_effective_elem)

            ith_decompositions = []
            ith_reduction_matrices = []

            need_to_reduce = True

            iteration_count_circuit_breaker = 100; # Set to some high number to start, will be estimated better below.
            iteration_count = 0
            reduction_count = 0

            reduction_latex_array = []

            ongoing_reduction_matrix__max_deg_tuple = None

            while need_to_reduce == True:

                check_connection_still_alive()

                iteration_count = iteration_count + 1

                if iteration_count >= iteration_count_circuit_breaker:
                    error_message = "Reached maximum iteration circuit breaker. Probably hit an infinite loop case."
                    print(error_message)
                    raise ValueError(f"May have hit infinite loop circuit_breaker={iteration_count_circuit_breaker}")
                    break
                
                min_estimated_reduction_count = 0

                
                # COMMENTING OUT if/else block. This was original logic, so we don't have to search the matrix for the highest degree on the first iteration, since we "know" it.
                # However a case came up with Ext0 where the element/term to be reduced first was not in the current Eij position. It may be that for Ext (aka W) it can
                # indeed be assumed that Eij position is the first element/term to be reduced (if it needs to be). Possible minor optimization.
                # if iteration_count == 1:
                #     # First reduction, so we can use pre-processed value.
                #     ongoing_reduction_matrix__max_deg_tuple = ongoing_reduction_matrix_with_effective_elem.max_deg_tuple
                # else:
                #     # Already done one or more reductions, which updates the value for `kth_Psi_Eij_c_tau_k`,
                #     # so need to search afresh for the max degree. 
                #     preferred_element_ordinal = ongoing_reduction_matrix__max_deg_tuple.zero_based_element_ordinal
                #     ongoing_reduction_matrix__max_deg_tuple = tee_utils.get_max_degree_of_tau_in_matrix(ongoing_reduction_matrix, subvariant, Eij_and_Aij_matrices, preferred_element_ordinal)

                # NEW LOGIC: Always search the matrix fresh for the term to be reduced.
                require_tau_to_be_rightmost_factor = True
                ongoing_reduction_matrix__max_deg_tuple = tee_utils.get_max_degree_of_tau_in_matrix(ongoing_reduction_matrix, require_tau_to_be_rightmost_factor, subvariant, Eij_and_Aij_matrices)

                ongoing_reduction_matrix__element_ordinal = ongoing_reduction_matrix__max_deg_tuple.zero_based_element_ordinal
                ongoing_reduction__Nij = Eij_and_Aij_matrices.Nij_list[ongoing_reduction_matrix__element_ordinal]
                ongoing_reduction__contains_constant_term = Eij_and_Aij_matrices.contains_constant_term_list[ongoing_reduction_matrix__element_ordinal]

                max_deg_found = ongoing_reduction__Nij

                if subvariant == PhiPsiSubvariant.CONSTANT_TERM_LOGIC:
                    if ongoing_reduction__contains_constant_term:
                        need_to_reduce = ongoing_reduction_matrix__max_deg_tuple.max_degree > ongoing_reduction__Nij
                        if need_to_reduce:
                            min_estimated_reduction_count = ongoing_reduction_matrix__max_deg_tuple.max_degree - ongoing_reduction__Nij
                    else:
                        need_to_reduce = ongoing_reduction_matrix__max_deg_tuple.max_degree >= ongoing_reduction__Nij
                        if need_to_reduce:
                            min_estimated_reduction_count = 1+(ongoing_reduction_matrix__max_deg_tuple.max_degree - ongoing_reduction__Nij)
                elif subvariant == PhiPsiSubvariant.ORIGINAL_LOGIC:
                    need_to_reduce = ongoing_reduction_matrix__max_deg_tuple.max_degree >= ongoing_reduction__Nij
                    if need_to_reduce:
                        min_estimated_reduction_count = 1+(ongoing_reduction_matrix__max_deg_tuple.max_degree - ongoing_reduction__Nij)
                elif subvariant == PhiPsiSubvariant.EMBEDDING_LOGIC:
                    if ongoing_reduction_matrix__max_deg_tuple.max_degree >= ongoing_reduction__Nij:
                        ongoing_reduction_matrix = kth_Eij_tau_k
                        need_to_reduce = True
                        max_deg_found = ongoing_reduction__Nij
                        min_estimated_reduction_count = 1 # Should only require one reduction
                    else:
                        need_to_reduce = False
                else:
                    raise ValueError(f"Unhandled subvariant={subvariant}")
                
                if iteration_count == 1:
                    kth_Psi_Eij_c_tau_k__row = kth_Psi_Eij_c_tau_k_with_effective_elem.effective_zero_based_row
                    kth_Psi_Eij_c_tau_k__col = kth_Psi_Eij_c_tau_k_with_effective_elem.effective_zero_based_col
                    kth_Psi_Eij_c_tau_k__expression = ongoing_reduction_matrix[kth_Psi_Eij_c_tau_k__row, kth_Psi_Eij_c_tau_k__col]
                    kth_Psi_Eij_c_tau_k__term_count = len(kth_Psi_Eij_c_tau_k__expression.as_ordered_terms())
                    
                    iteration_count_circuit_breaker = min_estimated_reduction_count * kth_Psi_Eij_c_tau_k__term_count
                    iteration_count_circuit_breaker = iteration_count_circuit_breaker + 1

                    # Give some more breathing room just in case.
                    #  TODO: See if worst-case loop iteration count can be proven exactly and use that as circuit breaker.
                    iteration_count_circuit_breaker = iteration_count_circuit_breaker * 10

                if need_to_reduce == True:

                    total_number_of_reductions_made = total_number_of_reductions_made + 1
                    reduction_count = reduction_count + 1

                    tau_degree_surplus = ongoing_reduction_matrix__max_deg_tuple.max_degree - (max_deg_found)

                    Eij_matrix = Eij_and_Aij_matrices.Eij_matrices[ongoing_reduction_matrix__max_deg_tuple.zero_based_element_ordinal]

                    if FEATURE_FLAG__INCLUDE_INTERMEDIATE_REDUCTION_MATRICES:
                        if len(ith_reduction_matrices) > 0 :
                            # Only append the matrix to be decomposed if it's not the first reduction pass.
                            ith_reduction_matrices.append(ongoing_reduction_matrix)

                    new_Aij_matrix = MU_SYMBOL_NON_COMM * TAU_SYMBOL_NON_COMM**tau_degree_surplus * Eij_matrix.raw * rightmost_basis.basis_matrix_inverse_non_comm
                    new_Aij_matrix = new_Aij_matrix.applyfunc(tee_utils.move_tau_to_the_right_of_other_symbols)
                    new_Aij_matrix = new_Aij_matrix.applyfunc(lambda matrix_elem: toggle_commutativity_for_all_symbols(matrix_elem, change_to_commutative=True))

                    if phi_psi_options.reduce_latex_output == False:
                        Eij_matrix_latex = matrix_to_latex(Eij_matrix.raw)

                    element_position_text = ""

                    if phi_psi_options.show_reduction_steps == True:
                        if len(reduction_latex_array) > 0:
                            reduction_latex_array.append("<br><br>")

                        element_position_text = f"{(ongoing_reduction_matrix__max_deg_tuple.zero_based_row+1), (ongoing_reduction_matrix__max_deg_tuple.zero_based_col+1), }"
                        reduction_latex_array.append("<div style=\"text-align:left; text-decoration:underline; font-weight:bold;\">Reduction " + str(iteration_count) + "</div>")
                        context_text = ""
                        if subvariant == PhiPsiSubvariant.CONSTANT_TERM_LOGIC:
                            context_text = "contains constant term" if ongoing_reduction__contains_constant_term == True else "does not contain constant term"
                        
                        if context_text:
                            context_text = " ("+context_text+")"
                        # context_text += f", must reduce to {max_deg_found}"
                        
                        reduction_latex_array.append("\\mathbf{$N_{"+element_position_text+"}$ =\\text{ }"+str(ongoing_reduction__Nij)+" \\text{"+context_text+"}")
                        reduction_latex_array.append("\\mathbf{$A_{"+element_position_text+"}$ = \\mu\\cdot\\tau^k\cdot$E_{"+element_position_text+"}$\cdot$U^{-1}$}")
                        reduction_latex_array.append(f"= \\mu\\cdot\\tau^{tau_degree_surplus}\\cdot"+Eij_matrix_latex+"\\cdot"+U_inverse_matrix_latex)
                        reduction_latex_array.append(f"= " + matrix_to_latex(new_Aij_matrix))

                    new_Aij_Phi = tee_utils.multiply_matrices_non_comm(new_Aij_matrix, Phi_module.simplified, finish=False)
                    new_Psi_Aij = tee_utils.multiply_matrices_non_comm(Psi_module.simplified, new_Aij_matrix, finish=False)

                    if phi_psi_options.show_reduction_steps == True:
                        reduction_latex_array.append("\\mathbf{\Phi_t = " + matrix_to_latex(Phi_module.simplified))
                        reduction_latex_array.append("\\mathbf{\Psi_t = " + matrix_to_latex(Psi_module.simplified))

                    if phi_psi_options.show_reduction_steps == True:
                        reduction_latex_array.append("\\mathbf{$A_{"+element_position_text+"}$\cdot\Phi_t - \Psi_t\cdot$A_{"+element_position_text+"}$} = " + matrix_to_latex(new_Aij_Phi) + " - " + matrix_to_latex(new_Psi_Aij))

                    new_Aij_Phi = new_Aij_Phi.applyfunc(tee_utils.move_tau_to_the_right_of_other_symbols)
                    new_Psi_Aij = new_Psi_Aij.applyfunc(tee_utils.move_tau_to_the_right_of_other_symbols)

                    if phi_psi_options.show_reduction_steps == True:
                        reduction_latex_array.append("= " + matrix_to_latex(new_Aij_Phi) + " - " + matrix_to_latex(new_Psi_Aij) + "\\textbf{(apply twist)}")

                    move_tau = False # We already applied the twist to new_Psi_Aij above.
                    new_Aij_Phi = tee_utils.finish_multiply_matrices_non_comm(new_Aij_Phi, move_tau)
                    new_Psi_Aij = tee_utils.finish_multiply_matrices_non_comm(new_Psi_Aij, move_tau)
                    
                    if phi_psi_options.show_reduction_steps == True:
                        reduction_latex_array.append("= " + matrix_to_latex(new_Aij_Phi) + " - " + matrix_to_latex(new_Psi_Aij) +  "\\textbf{(simplify terms)}")

                    new_Aij_Phi__minus__Psi_Aij = new_Aij_Phi - new_Psi_Aij
                    new_Aij_Phi__minus__Psi_Aij = new_Aij_Phi__minus__Psi_Aij.expand()

                    if phi_psi_options.show_reduction_steps == True:
                        reduction_latex_array.append("= " + matrix_to_latex(new_Aij_Phi__minus__Psi_Aij) +  "\\textbf{(after subtraction)}")

                    # reduction_latex_array.append("= " + matrix_to_latex(new_Aij_Phi__minus__Psi_Aij))

                    if FEATURE_FLAG__INCLUDE_INTERMEDIATE_REDUCTION_MATRICES:
                        ith_reduction_matrices.append(new_Aij_Phi__minus__Psi_Aij.copy())

                    # reduction_latex_array.append("= " + matrix_to_latex(new_Aij_Phi__minus__Psi_Aij))

                    commutative_true = True
                    new_Aij_Phi__minus__Psi_Aij = new_Aij_Phi__minus__Psi_Aij.applyfunc(lambda matrix_elem: tee_utils.toggle_commutativity_for_all_symbols(matrix_elem, commutative_true))
                    new_Aij_Phi__minus__Psi_Aij = tee_utils.convert_matrix_to_simple_terms(new_Aij_Phi__minus__Psi_Aij)

                    if phi_psi_options.group_coefficients == True:
                        new_Aij_Phi__minus__Psi_Aij = new_Aij_Phi__minus__Psi_Aij.applyfunc(tee_utils.group_coefficients_by_tau)

                        if phi_psi_options.show_reduction_steps == True:
                            reduction_latex_array.append(" = " + matrix_to_latex(new_Aij_Phi__minus__Psi_Aij) + "\\textbf{(group coefficients by} \\tau \\textbf{)}")

                    
                    non_comm_coefficient = tee_utils.toggle_commutativity_for_all_symbols(ongoing_reduction_matrix__max_deg_tuple.coefficient, True)
                    new_Aij_Phi__minus__Psi_Aij__effective_expression = new_Aij_Phi__minus__Psi_Aij[ongoing_reduction_matrix__max_deg_tuple.zero_based_row, ongoing_reduction_matrix__max_deg_tuple.zero_based_col]
                    new_Aij_Phi__minus__Psi_Aij__max_deg_tuple = get_max_degree_of_tau_in_expression(new_Aij_Phi__minus__Psi_Aij__effective_expression)
                    # print("new_Aij_Phi__minus__Psi_Aij__max_deg_tuple.coefficient=")
                    # sympy.pprint(new_Aij_Phi__minus__Psi_Aij__max_deg_tuple.coefficient)
                    solve_for_mu_equation = sympy.Eq(new_Aij_Phi__minus__Psi_Aij__max_deg_tuple.coefficient, non_comm_coefficient)
                    mu_value_list = sympy.solve(solve_for_mu_equation, MU_SYMBOL)

                    mu_value = None
                    if len(mu_value_list) > 0:
                        if len(mu_value_list) > 1:
                            print(f"WARNING: Found more than one solution for mu in {solve_for_mu_equation}")
                        mu_value = mu_value_list[0]
                    else:
                        raise ValueError(f"Could not solve for mu in {solve_for_mu_equation}")
                    
                    new_Aij_Phi__minus__Psi_Aij = new_Aij_Phi__minus__Psi_Aij.applyfunc(lambda matrix_elem: matrix_elem.subs(MU_SYMBOL, mu_value))

                    if phi_psi_options.show_reduction_steps == True:
                        # reduction_latex_array.append("= " + matrix_to_latex(new_Aij_Phi) + " - " + matrix_to_latex(new_Psi_Aij) + "\\textbf{(substitute mu with} "+sympy.latex(mu_value)+" \\textbf{)}}")
                        reduction_latex_array.append(" = " + matrix_to_latex(new_Aij_Phi__minus__Psi_Aij) + "\\textbf{(substitute} \\mu \\textbf{ with } "+sympy.latex(mu_value)+" \\textbf{)}}")

                    strip_integer_exponents_for_symbolic_execution = False
                    new_Aij_Phi__minus__Psi_Aij = tee_utils.convert_matrix_to_simple_terms(new_Aij_Phi__minus__Psi_Aij, strip_integer_exponents_for_symbolic_execution)
                    if phi_psi_options.show_reduction_steps == True:
                        reduction_latex_array.append(" = " + matrix_to_latex(new_Aij_Phi__minus__Psi_Aij) + "\\textbf{(simplify terms after substitution)}")

                    if further_process_q_k == True:
                        new_Aij_Phi__minus__Psi_Aij = tee_utils.convert_matrix_to_simple_terms(new_Aij_Phi__minus__Psi_Aij, strip_integer_exponents_for_symbolic_execution=stripping_integers, p=p_val, n=n_val)

                        if phi_psi_options.show_reduction_steps == True:
                            if stripping_integers == True:
                                reduction_latex_array.append(" = " + matrix_to_latex(new_Aij_Phi__minus__Psi_Aij) + "\\textbf{(strip exponents from integers)}")
                            # elif replacing_q_with_p_n == True:
                            #     reduction_latex_array.append(" = " + matrix_to_latex(new_Aij_Phi__minus__Psi_Aij) + "\\textbf{(replace q with p^n ("+str(p_val)+"^"+str(n_val)+") )}")


                    if phi_psi_options.group_coefficients == True:
                        new_Aij_Phi__minus__Psi_Aij = new_Aij_Phi__minus__Psi_Aij.applyfunc(tee_utils.group_coefficients_by_tau)

                        if phi_psi_options.show_reduction_steps == True:
                            reduction_latex_array.append(" = " + matrix_to_latex(new_Aij_Phi__minus__Psi_Aij) + "\\textbf{(AGAIN group coefficients by} \\tau \\textbf{)}")

                    reduced_matrix = ongoing_reduction_matrix - new_Aij_Phi__minus__Psi_Aij
                    reduced_matrix = reduced_matrix.expand()
                    reduced_matrix = tee_utils.convert_matrix_to_simple_terms(reduced_matrix, strip_integer_exponents_for_symbolic_execution)

                    if phi_psi_options.show_reduction_steps == True:
                        reduction_latex_array.append("\\mathbf{"+matrix_to_latex(ongoing_reduction_matrix) + " - " + matrix_to_latex(new_Aij_Phi__minus__Psi_Aij) + "}")
                        reduction_latex_array.append(" = " + matrix_to_latex(reduced_matrix) + "\\textbf{(after subtraction)}")

                    if further_process_q_k == True and replacing_q_with_p_n == True:
                        reduced_matrix = tee_utils.convert_matrix_to_simple_terms(reduced_matrix, strip_integer_exponents_for_symbolic_execution=False, p=p_val, n=n_val)

                        if phi_psi_options.show_reduction_steps == True:
                            reduction_latex_array.append(" = " + matrix_to_latex(reduced_matrix) + "\\textbf{(again replace q with p^n ("+str(p_val)+"^"+str(n_val)+") )}")

                    if phi_psi_options.group_coefficients == True:
                        reduced_matrix = reduced_matrix.applyfunc(tee_utils.group_coefficients_by_tau)

                        if phi_psi_options.show_reduction_steps == True:
                            reduction_latex_array.append(" = " + matrix_to_latex(reduced_matrix) + "\\textbf{(group coefficients by} \\tau \\textbf{)}")


                    if FEATURE_FLAG__INCLUDE_INTERMEDIATE_REDUCTION_MATRICES:
                        ith_reduction_matrices.append(new_Aij_Phi__minus__Psi_Aij)

                    ongoing_reduction_matrix = reduced_matrix
                    
                else: # if need_to_reduce == False:
                    
                    frobenius_expanded_final_reduction = sympy.zeros(ongoing_reduction_matrix.rows, ongoing_reduction_matrix.cols)

                    for row_index in range(ongoing_reduction_matrix.rows):
                        for col_index in range(ongoing_reduction_matrix.cols):
                            element = ongoing_reduction_matrix[row_index, col_index]

                            check_connection_still_alive()

                            if element == 0:
                                continue

                            all_degrees_of_tau = tee_utils.get_all_degrees_of_tau_in_expression(element)

                            coefficients_by_degrees = {}

                            for degree_tuple in all_degrees_of_tau:
                                existing_list = coefficients_by_degrees.get(degree_tuple.max_degree)

                                if existing_list:
                                    pass
                                else:
                                    existing_list = []
                                    coefficients_by_degrees[degree_tuple.max_degree] = existing_list
                                
                                existing_list.append(degree_tuple.coefficient)

                            # Start off with "zero" expression.
                            element__frobenius_expanded_expression = sympy.Integer(0)

                            coefficients_by_degrees_iterator = coefficients_by_degrees.items();

                            if phi_psi_options.show_reduction_steps == True:
                                # Use `sorted()` here just so it shows tau-pows in expected order.
                                coefficients_by_degrees_iterator = sorted(coefficients_by_degrees_iterator)
                            
                            for degree_of_tau, list_of_coefficients_for_this_degree in coefficients_by_degrees_iterator:
                                sum_of_coefficients = sum(list_of_coefficients_for_this_degree)

                                if phi_psi_options.group_coefficients == True:
                                    # "Undo" expansion for better readability
                                    if hasattr(sum_of_coefficients, "expand") and callable(getattr(sum_of_coefficients, "expand")):
                                        sum_of_coefficients = tee_utils.expand_frobenius_2(sum_of_coefficients, Q_SYMBOL)
                                        sum_of_coefficients = sum_of_coefficients.expand()

                                if phi_psi_options.show_reduction_steps == True:
                                    element__frobenius_expanded_expression += (sum_of_coefficients) * TAU_SYMBOL_NON_COMM**degree_of_tau

                                if phi_psi_options.replace_c_with_tau == True:
                                    sum_of_coefficients = tee_utils.replace_c_with_tau(sum_of_coefficients)

                                if replacing_q_with_p_n:
                                    # print(f"reduce_to_finite_field start {sum_of_coefficients}, p={phi_psi_options.prime_characteristic}, n={phi_psi_options.extension_degree}")
                                    sum_of_coefficients = tee_ai_bridge.reduce_to_finite_field(sum_of_coefficients, phi_psi_options.prime_characteristic, phi_psi_options.extension_degree)
                                    # print(f"reduce_to_finite_field end {sum_of_coefficients}")

                                if phi_psi_options.replace_q_with_integer == True:
                                    sum_of_coefficients = tee_utils.replace_q_power_with_integer(sum_of_coefficients)

                                ith_decomposition = Psi_Eij_c_tau_k_decomposition(row_index, col_index, degree_of_tau, sum_of_coefficients)

                                ith_decompositions.append(ith_decomposition)

                            frobenius_expanded_final_reduction[row_index, col_index] = element__frobenius_expanded_expression

                    # If not showing reduction steps in UI, then we save a few cycles by not evaluating/rendering
                    # the frobenius-expanded version of the indvidual reduction. The final Ext/Ext0 matrix will be
                    # the same regardless.
                    if phi_psi_options.show_reduction_steps == True:
                        ith_reduction_matrices.append(frobenius_expanded_final_reduction)
                        final_reduction_matrices.append(frobenius_expanded_final_reduction)

                        reduction_latex_array.append("<br><br>")
                        reduction_latex_array.append("<div style=\"text-align:left; text-decoration:underline; font-weight:bold;\">Final Frobenius-expanded Reduction</div>")
                        reduction_latex_array.append(" = "+matrix_to_latex(frobenius_expanded_final_reduction)+ "")
                    else:
                        ith_reduction_matrices.append(ongoing_reduction_matrix)
                        final_reduction_matrices.append(ongoing_reduction_matrix)

                if reduction_limit > 0 and iteration_count >= reduction_limit:
                    if need_to_reduce == True:
                        final_reduction_matrices.append(ongoing_reduction_matrix)
                    need_to_reduce = False

            decompositions.append(ith_decompositions)

            reduction_counts.append(reduction_count)
            if reduction_latex_array:
                reduction_matrices.append(reduction_latex_array)
            else:
                reduction_matrices.append([])
            
            i_times_k += 1
    
    r_matrix_row_count = len(Eij_tau_k_matrices_to_return) if subvariant == PhiPsiSubvariant.EMBEDDING_LOGIC else r_matrix_col_count
    r_matrix = sympy.zeros(r_matrix_row_count, r_matrix_col_count)

    curr_col_index_in_r_matrix = 0
    # It's possible to tighten up this nest of loops and even merge it with the above nest. Certainly it could be "less code".
    # But it's arguable whether other lower-code methods would be faster overall.
    # `n` for every nested loop here is generally going to be low in the grand scheme of things (few hundred max)
    # and there is nothing heavy-weight going on. So for now it's easier to set the R matrix elements as its own separate thing.
    # It makes it a bit more readable, robust, and maintainable this way.
    for ith_list_of_decompositions in decompositions:

        for ith_decomposition in ith_list_of_decompositions:

            needle = 1

            if ith_decomposition.degree_of_tau > TAU_DEGREE_OF_ZERO:
                needle = TAU_SYMBOL_NON_COMM**ith_decomposition.degree_of_tau

            for index_in_Eij_tau_k_matrices in range(len(Eij_tau_k_matrices_to_return)):
                ith_Eij_tau_k_matrix = Eij_tau_k_matrices_to_return[index_in_Eij_tau_k_matrices]

                element_in_ith_Eij_matrix = ith_Eij_tau_k_matrix.raw[ith_decomposition.row_index, ith_decomposition.col_index]
                if element_in_ith_Eij_matrix == needle:
                    col_index_in_r = curr_col_index_in_r_matrix
                    row_index_in_r = index_in_Eij_tau_k_matrices

                    final_entry = ith_decomposition.sum_of_coefficients

                    if hasattr(final_entry, 'as_ordered_terms'):
                        # Have to do this finaly simplify() for some reason, to get rid of random (-1) factors at the front of some terms.
                        simplified_terms = [(term.doit() if hasattr(term, 'doit') else term) for term in final_entry.as_ordered_terms()]
                        final_entry = sympy.Add(*simplified_terms)
                    elif hasattr(final_entry, 'doit'):
                        final_entry = final_entry.doit()

                    r_matrix[row_index_in_r, col_index_in_r] = final_entry

        curr_col_index_in_r_matrix += 1

    end_time = int(time.time_ns())
    elapsed_time = (end_time-start_time) / 1e9

    if phi_psi_options.reduce_log_output == False:
        print(f"reduction count for {subvariant}: {total_number_of_reductions_made}")
        print(f"execution time for {subvariant}: {elapsed_time} seconds")

    # leaf_count = sum(1 for _ in sympy.preorder_traversal(r_matrix))
    leaf_count = 0
    byte_count = 0
    if phi_psi_options.calc_leaf_and_byte_counts == True:
        leaf_count = tee_utils.mathematica_leaf_count(r_matrix)
        byte_count = tee_utils.mathematica_matrix_byte_count(r_matrix)

    if phi_psi_options.reduce_log_output == False:
        print(f"leaf_count={leaf_count}")
        print(f"byte_count={byte_count}")

    return R_Matrix(r_matrix, Eij_tau_k_matrices_to_return, Psi_Eij_c_tau_k_matrices_to_return, reduction_matrices, final_reduction_matrices, degree_k_start_values, degree_k_finish_values, reduction_counts, leaf_count, byte_count, total_number_of_reductions_made, elapsed_time)


def get_Eij_and_Aij_matrices(phi_module: sympy.Matrix, psi_module:sympy.Matrix, rightmost_basis: RightmostBasis, phi_psi_options: PhiPsiOptions) -> Eij_And_Aij_Matrices:
    initial_Eij_matrices = []
    initial_Aij_matrices = []
    initial_Aij_Phi__minus__Psi_Aij_matrices = []
    Nij_list = []
    contains_constant_term_list = []
    Aij_Phi_matrices: list[sympy.Matrix] = []
    Psi_Aij_matrices: list[sympy.Matrix] = []

    Eij_tau_k_matrices_wrapper: list[list[MatrixWithEffectiveElement]] = []
    Eij_c_tau_k_matrices_wrapper: list[list[MatrixWithEffectiveElement]] = []
    Psi_Eij_c_tau_k_matrices_wrapper: list[list[MatrixWithEffectiveElement]] = []

    Eij_row_count = psi_module.rows
    Eij_col_count = phi_module.rows
    Eij_element_count = Eij_row_count * Eij_col_count

    for i in range(Eij_element_count):

        check_connection_still_alive()

        zero_based_row = math.floor(i / Eij_col_count)
        zero_based_col = i % Eij_col_count

        ith_Eij_matrix_raw = sympy.zeros(Eij_row_count, Eij_col_count)
        ith_Eij_matrix_raw[zero_based_row, zero_based_col] = 1
        initial_Eij_matrices.append(MatrixWithEffectiveElement(ith_Eij_matrix_raw, zero_based_row, zero_based_col, 0, i, None))

        ith_Aij_matrix = MU_SYMBOL * ith_Eij_matrix_raw * rightmost_basis.basis_matrix_inverse
        initial_Aij_matrices.append(ith_Aij_matrix)

        # Only have to do the tau -> c^q^k conversion for Psi*Aij
        ith_Aij_Phi = tee_utils.multiply_matrices_non_comm(ith_Aij_matrix, phi_module, finish=True, move_tau=False)
        ith_Psi_Aij = tee_utils.multiply_matrices_non_comm(psi_module, ith_Aij_matrix, finish=True, move_tau=True)

        Aij_Phi_matrices.append(ith_Aij_Phi)
        Psi_Aij_matrices.append(ith_Psi_Aij)

        ith_Aij_Phi__minus__Psi_Aij_matrix = ith_Aij_Phi - ith_Psi_Aij
        ith_Aij_Phi__minus__Psi_Aij_matrix = ith_Aij_Phi__minus__Psi_Aij_matrix.expand()

        initial_Aij_Phi__minus__Psi_Aij_matrices.append(ith_Aij_Phi__minus__Psi_Aij_matrix)

        ith_Nij = get_max_degree_of_tau_in_matrix(ith_Aij_Phi__minus__Psi_Aij_matrix).max_degree

        Nij_list.append(int(ith_Nij))

        contains_constant_term = tee_utils.does_matrix_contain_constant_term(ith_Aij_Phi__minus__Psi_Aij_matrix)
        contains_constant_term_list.append(contains_constant_term)

        degree_k_start = 0
        degree_k_finish = ith_Nij+1

        Eij_tau_k_matrices = []
        Eij_c_tau_k_matrices = []
        Psi_Eij_c_tau_k_matrices = []

        for degree_k in range(degree_k_start, degree_k_finish):

            check_connection_still_alive()

            kth_Eij_tau_k = ith_Eij_matrix_raw * TAU_SYMBOL_NON_COMM**degree_k
            kth_Eij_c_tau_k = ith_Eij_matrix_raw * (C_SYMBOL*TAU_SYMBOL_NON_COMM**degree_k)
            kth_Psi_Eij_c_tau_k = tee_utils.multiply_matrices_non_comm(psi_module, kth_Eij_c_tau_k)

            if phi_psi_options.group_coefficients == True:
                # I don't think this is strictly necessary, but applying just in case.
                kth_Psi_Eij_c_tau_k = kth_Psi_Eij_c_tau_k.applyfunc(tee_utils.group_coefficients_by_tau)

            kth_Psi_Eij_c_tau_k = tee_utils.convert_matrix_to_simple_terms(kth_Psi_Eij_c_tau_k, strip_integer_exponents_for_symbolic_execution=False, p=phi_psi_options.prime_characteristic, n=phi_psi_options.extension_degree)

            kth_Eij_tau_k__max_deg_tuple = tee_utils.get_max_degree_of_tau_in_matrix(kth_Eij_tau_k)
            kth_Eij_c_tau_k__max_deg_tuple = tee_utils.get_max_degree_of_tau_in_matrix(kth_Eij_c_tau_k)
            kth_Psi_Eij_c_tau_k__max_deg_tuple = tee_utils.get_max_degree_of_tau_in_matrix(kth_Psi_Eij_c_tau_k)

            Eij_tau_k_matrices.append(MatrixWithEffectiveElement(kth_Eij_tau_k, zero_based_row, zero_based_col, degree_k, i, kth_Eij_tau_k__max_deg_tuple))
            Eij_c_tau_k_matrices.append(MatrixWithEffectiveElement(kth_Eij_c_tau_k, zero_based_row, zero_based_col, degree_k, i, kth_Eij_c_tau_k__max_deg_tuple))
            Psi_Eij_c_tau_k_matrices.append(MatrixWithEffectiveElement(kth_Psi_Eij_c_tau_k, zero_based_row, zero_based_col, degree_k, i, kth_Psi_Eij_c_tau_k__max_deg_tuple))
        
        Eij_tau_k_matrices_wrapper.append(Eij_tau_k_matrices)
        Eij_c_tau_k_matrices_wrapper.append(Eij_c_tau_k_matrices)
        Psi_Eij_c_tau_k_matrices_wrapper.append(Psi_Eij_c_tau_k_matrices)
        

    return Eij_And_Aij_Matrices(
        Eij_row_count,
        Eij_col_count,
        initial_Eij_matrices,
        initial_Aij_matrices,
        initial_Aij_Phi__minus__Psi_Aij_matrices,
        Nij_list,
        contains_constant_term_list,
        Aij_Phi_matrices,
        Psi_Aij_matrices,

        Eij_tau_k_matrices_wrapper,
        Eij_c_tau_k_matrices_wrapper,
        Psi_Eij_c_tau_k_matrices_wrapper
    )
