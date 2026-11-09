
import sympy
import math
import re

from sympy.parsing.sympy_parser import parse_expr

from latex2sympy2_extended import latex2sympy

# from sympy.parsing.latex import parse_latex_lark

from sympy.parsing.latex import parse_latex
from sympy import UnevaluatedExpr

from sympy.core.parameters import evaluate

import tmc_types
import tmc_utils
import tmc_ai_bridge

from tmc_consts import *
from tmc_types import *


LATEX__MATRIX_STRETCH = "\\renewcommand{\\arraystretch}{1.5}"
LATEX__MATRIX_STRETCH_BIGGER = "\\renewcommand{\\arraystretch}{2}"
LATEX__H_SPACE = "\\hspace{5px}"

def r_matrix_to_dict(r_matrix:R_Matrix):

    Eij_tau_k_matrices_latex = []
    Psi_Eij_c_tau_k_matrices_latex = []
    reduction_matrices_latex = []
    final_reduction_matrices_latex = []

    to_return = {
        "Eij_tau_k_matrices_latex": list_of_matrices_to_list_of_latex_strings(r_matrix.Eij_tau_k_matrices),
        "Psi_Eij_c_tau_k_matrices_latex": list_of_matrices_to_list_of_latex_strings(r_matrix.Psi_Eij_c_tau_k_matrices),
        # "reduction_matrices_latex": list_of_matrices_to_list_of_latex_strings(r_matrix.reduction_matrices),
        "reduction_matrices_latex": r_matrix.reduction_matrices,
        "final_reduction_matrices_latex": list_of_matrices_to_list_of_latex_strings(r_matrix.final_reduction_matrices),
        "r_matrix": r_matrix_to_2d_array_of_latex(r_matrix.raw),
        "degree_k_start_values": r_matrix.degree_k_start_values,
        "degree_k_finish_values": r_matrix.degree_k_finish_values,
        "reduction_counts": r_matrix.reduction_counts,
        "total_reduction_count": r_matrix.total_reduction_count,
        "total_seconds_to_generate": r_matrix.total_seconds_to_generate,
        "leaf_count": r_matrix.leaf_count,
        "byte_count": r_matrix.byte_count
    }

    return to_return

def phi_psi_input_to_dict(phi_psi_input: PhiPsiInput):
    Phi_module_simplified_comm = phi_psi_input.Phi_module.simplified
    Psi_module_simplified_comm = phi_psi_input.Psi_module.simplified

    # Phi_module_simplified_comm = phi_psi_input.Phi_module.simplified.subs(TAU_SYMBOL_NON_COMM, TAU_SYMBOL)
    # Psi_module_simplified_comm = phi_psi_input.Psi_module.simplified.subs(TAU_SYMBOL_NON_COMM, TAU_SYMBOL)
    Phi_simplified_latex = matrix_to_latex_for_simplify_and_expand_t_module(Phi_module_simplified_comm)
    Psi_simplified_latex = matrix_to_latex_for_simplify_and_expand_t_module(Psi_module_simplified_comm)

    return {
        "was_able_to_condense_phi": True if Phi_module_simplified_comm is not None else False,
        "was_able_to_condense_psi": True if Psi_module_simplified_comm is not None else False,
        "was_able_to_validate_phi": True if len(phi_psi_input.Phi_module.parse_errors) == 0 else False,
        "was_able_to_validate_psi": True if len(phi_psi_input.Psi_module.parse_errors) == 0 else False,
        "phi_module_condensed_latex": Phi_simplified_latex,
        "psi_module_condensed_latex": Psi_simplified_latex,
        "input_errors": phi_psi_input.input_errors
    }

def initial_analysis_to_dict(initial_analysis: InitialAnalysis, Phi_module: T_Module, Psi_module: T_Module):

    U_matrix_original_latex = LATEX__MATRIX_STRETCH + matrix_to_latex(initial_analysis.rightmost_basis.basis_matrix)
    U_matrix_inverse_latex = LATEX__MATRIX_STRETCH + matrix_to_latex(initial_analysis.rightmost_basis.basis_matrix_inverse)

    Eij_matrices_latex = []
    Aij_matrices_latex = []
    Aij_Phi__minus__Psi_Aij_matrices_latex = []

    if initial_analysis.Eij_and_Aij_matrices is not None:
        for i in range(len(initial_analysis.Eij_and_Aij_matrices.Eij_matrices)):

            ith_Eij_matrix = initial_analysis.Eij_and_Aij_matrices.Eij_matrices[i]
            row = ith_Eij_matrix.effective_zero_based_row
            col = ith_Eij_matrix.effective_zero_based_col
            # row = math.floor(i / initial_analysis.Eij_and_Aij_matrices.Eij_col_count)
            # col = i % initial_analysis.Eij_and_Aij_matrices.Eij_col_count

            
            ith_Aij_matrix = initial_analysis.Eij_and_Aij_matrices.Aij_matrices[i]
            ith_Aij_Phi__minus__Psi_Aij_matrix = initial_analysis.Eij_and_Aij_matrices.Aij_Phi_minus_Psi_Aij_matrices[i];

            ith_Aij_Phi_matrix = initial_analysis.Eij_and_Aij_matrices.Aij_Phi_matrices[i];
            ith_Psi_Aij_matrix = initial_analysis.Eij_and_Aij_matrices.Psi_Aij_matrices[i];

            ith_Eij_matrix_latex = LATEX__H_SPACE + "$E_{"+str(row+1)+","+str(col+1)+"}$=" + LATEX__MATRIX_STRETCH + matrix_to_latex(ith_Eij_matrix.raw)
            ith_Aij_matrix_latex = " $A_{"+str(row+1)+","+str(col+1)+"}$=" + LATEX__MATRIX_STRETCH + matrix_to_latex(ith_Aij_matrix)
            ith_Aij_Phi__minus__Psi_Aij_matrix_latex = LATEX__MATRIX_STRETCH_BIGGER + matrix_to_latex(ith_Aij_Phi__minus__Psi_Aij_matrix)

            ith_Aij_Phi_matrix_latex = LATEX__MATRIX_STRETCH_BIGGER + matrix_to_latex(ith_Aij_Phi_matrix)
            ith_Psi_Aij_matrix_latex = LATEX__MATRIX_STRETCH_BIGGER + matrix_to_latex(ith_Psi_Aij_matrix)

            final_Aij_Phi__minus__Psi_Aij_matrix_latex = ith_Aij_Phi_matrix_latex + " - " + ith_Psi_Aij_matrix_latex + " = " + ith_Aij_Phi__minus__Psi_Aij_matrix_latex

            Eij_matrices_latex.append(ith_Eij_matrix_latex)
            Aij_matrices_latex.append(ith_Aij_matrix_latex)
            Aij_Phi__minus__Psi_Aij_matrices_latex.append(final_Aij_Phi__minus__Psi_Aij_matrix_latex)

    return {
        "util_matrix_latex": U_matrix_original_latex,
        "inverse_matrix_latex": U_matrix_inverse_latex,
        "e_matrices_latex": Eij_matrices_latex,
        "a_matrices_latex": Aij_matrices_latex,
        "multiply_by_a_latex": Aij_Phi__minus__Psi_Aij_matrices_latex,
        "max_degrees": initial_analysis.Eij_and_Aij_matrices.Nij_list if initial_analysis.Eij_and_Aij_matrices is not None else [],
        "contains_constant_term_list": initial_analysis.Eij_and_Aij_matrices.contains_constant_term_list if initial_analysis.Eij_and_Aij_matrices is not None else [],
        "seconds_to_generate_preprocessing": initial_analysis.total_seconds_to_generate,
        "seconds_to_generate_rightmost_basis": initial_analysis.rightmost_basis.total_seconds_to_generate
    }

def dimensional_analysis_matrices_to_2d_arrays_of_latex(matrices:list[sympy.Matrix], matrix_subset: MatrixSubset):
    to_return = []

    for i in range(len(matrices)):
        ith_matrix = matrices[i]
        to_return.append(dimensional_analysis_matrix_to_2d_array_of_latex(ith_matrix, i, matrix_subset))

    return to_return


def dimensional_analysis_matrix_to_2d_array_of_latex(matrix:sympy.Matrix, matrix_index: int, matrix_subset: MatrixSubset):
    to_return = []

    for zero_based_row in range(matrix.rows):
        new_row = []

        for zero_based_col in range(matrix.cols):
            element = matrix[zero_based_row, zero_based_col]
            element_latex = sympy.latex(element, order='none')

            is_highlighted_col = False

            for ith_col_index in matrix_subset.col_indeces:
                if ith_col_index.matrix_index == matrix_index and ith_col_index.col_index == zero_based_col:
                    is_highlighted_col = True

            if is_highlighted_col:
                element_latex = "\\textcolor{red}{\\bm{" + element_latex + "}}"

            new_row.append(element_latex)

        to_return.append(new_row)
    
    return to_return

def r_matrix_to_2d_array_of_latex(r_matrix: sympy.Matrix):
    to_return = []

    for zero_based_row in range(r_matrix.rows):
        new_row = []

        for zero_based_col in range(r_matrix.cols):
            element = r_matrix[zero_based_row, zero_based_col]
            if element == 0:
                new_row.append("0")
            else:
                # element_latex = "\\textcolor{orange}{\\mathbf{" + sympy.latex(element) + "}}"
                # element_latex = "\\textcolor{orange}{" + sympy.latex(element) + "}"
                element_latex = sympy.latex(element)
                new_row.append(element_latex)

        to_return.append(new_row)
    
    return to_return

def matrix_to_latex_for_simplify_and_expand_t_module(matrix: sympy.Matrix):
    if matrix is not None:
        return sympy.latex(matrix, mat_str='pmatrix', mat_delim="", order='none')
    else:
        return "\\hspace{5px}\\text{(Could not simplify)}";
    # return sympy.latex(matrix, mat_str='pmatrix', mat_delim="", order='rev-lex')
    # return sympy.latex(matrix, mat_str='pmatrix', mat_delim="")

def matrix_to_latex(matrix: sympy.Matrix):
    return sympy.latex(matrix, mat_delim="(")

# def matrix_to_latex(matrix: sympy.Matrix):
#     return sympy.latex(matrix, mat_str='pmatrix')

def matrix_to_latex_square_brackets(matrix):
    return sympy.latex(matrix, mat_delim="[")

def list_of_matrices_to_list_of_latex_strings(list_of_matrices):
    to_return = []
    for ith in list_of_matrices:
        to_return.append(matrix_to_latex(ith))

    return to_return

def list_of_matrices_to_list_of_latex_strings(list_of_matrices, include_matrix_stretch = True) -> list[str]:
    if len(list_of_matrices) <= 0:
        return []

    to_return = []

    latex_stretch = LATEX__MATRIX_STRETCH if include_matrix_stretch else ""

    for ith in list_of_matrices:

        ith_json = "";

        if isinstance(ith, MatrixWithEffectiveElement):
            ith_json = {
                "raw": latex_stretch + matrix_to_latex(ith.raw),
                "effective_zero_based_row": ith.effective_zero_based_row,
                "effective_zero_based_col": ith.effective_zero_based_col,
                "degree_of_tau": ith.degree_of_tau,
            }
        elif isinstance(ith, list):
            ith_json = list_of_matrices_to_delimited_latex_string(ith)
        else:
            if ith != None:
                ith_json = latex_stretch + matrix_to_latex(ith)
            else:
                ith_json = ""

        to_return.append(ith_json)
    
    return to_return

def list_of_matrices_to_delimited_latex_string(list_of_matrices, delimiter = ",", space_before_delimiter = False, latex_h_space=LATEX__H_SPACE, include_matrix_stretch=True):
    if len(list_of_matrices) <= 0:
        return ""

    to_return = LATEX__MATRIX_STRETCH if include_matrix_stretch else ""
    add_delimiter = False

    index = 0

    for ith in list_of_matrices:
        if add_delimiter:
            # to_return += "\\Rightarrow" + latex_h_space
            if space_before_delimiter:
                to_return += latex_h_space
            to_return += delimiter + latex_h_space
        else:
            add_delimiter = True

        to_return += matrix_to_latex(ith)

        index = index + 1
    
    return to_return

def safe_parse_with_symbolic_d(latex_str):
    # Match any LaTeX command (\delta, \frac, etc.) OR the single character 'd'
    pattern = r'(\\[a-zA-Z]+)|(d)'
    
    def replace(match):
        # If group 1 matched, it's a safe LaTeX command. Return it unchanged.
        if match.group(1):
            return match.group(1)
        # Otherwise, it's a standalone variable 'd'. Wrap it.
        return '{d}'
        
    sanitized = re.sub(pattern, replace, latex_str)
    to_return = ""

    # print(f"sanitized={sanitized}")


    # to_return = parse_latex(sanitized)

    with sympy.evaluate(False):
        to_return = latex2sympy(sanitized)

    # sanitized = sanitized.replace("\\", "\\\\")
    # to_return = parse_expr(sanitized, evaluate=False)

    # to_return = parse_latex_lark(sanitized)

    return to_return

def convert_t_module_latex_matrix_to_sympy_matrix(ith_matrix, variant: Variant) -> [sympy.Matrix]:
    d = len(ith_matrix)
    ith_matrix_symbolic = sympy.zeros(d, d)

    for row_index, row in enumerate(ith_matrix):
        for col_index, element_latex in enumerate(row):
            element_latex = element_latex.replace("^{}", "")
            # element_latex = element_latex.replace("d", "{d}")
            # print(f"element_latex={element_latex}")
            # element = latex2sympy(element_latex)
            element = safe_parse_with_symbolic_d(element_latex)
            # element = parse_latex_lark(element_latex)
            # print(f"element_latex={element_latex}, element={element}")

            if variant == Variant.VARIANT__DIAGONAL:
                element = tmc_utils.toggle_commutativity_for_all_symbols(element, False)
            
            ith_matrix_symbolic[row_index, col_index] = element
    
    return ith_matrix_symbolic

def log_input_in_various_ways(phi_simplified: sympy.Matrix, psi_simplified: sympy.Matrix):
    Phi_as_sage_string = f"{phi_simplified.tolist()}"
    Psi_as_sage_string = f"{psi_simplified.tolist()}"

    Phi_as_mathematica_string = Phi_as_sage_string.replace("[", "{")
    Phi_as_mathematica_string = Phi_as_mathematica_string.replace("]", "}")
    Phi_as_mathematica_string = Phi_as_mathematica_string.replace("theta", "Θ")
    Phi_as_mathematica_string = Phi_as_mathematica_string.replace("tau", "τ")

    Psi_as_mathematica_string = Psi_as_sage_string.replace("[", "{")
    Psi_as_mathematica_string = Psi_as_mathematica_string.replace("]", "}")
    Psi_as_mathematica_string = Psi_as_mathematica_string.replace("theta", "Θ")
    Psi_as_mathematica_string = Psi_as_mathematica_string.replace("tau", "τ")
    
    print("SAGE:")
    print(f"Phi={Phi_as_sage_string}")
    print(f"Psi={Psi_as_sage_string}")

    print("MATHEMATICA:")
    print(f"Phi={Phi_as_mathematica_string}")
    print(f"Psi={Psi_as_mathematica_string}")

    print("LATEX:")
    print(f"Phi={matrix_to_latex(phi_simplified)}")
    print(f"Phi={matrix_to_latex(psi_simplified)}")

def convert_t_module_latex_list_to_list_of_matrices(t_module_latex_array, variant: Variant) -> [sympy.Matrix]:
    to_return = []

    for ith_matrix in t_module_latex_array:
        ith_matrix_symbolic = convert_t_module_latex_matrix_to_sympy_matrix(ith_matrix, variant)
        
        to_return.append(ith_matrix_symbolic)
    
    return to_return

def parse_t_module_latex(t_module_latex_or_expr, module_name, phi_psi_options: PhiPsiOptions) -> tmc_types.T_Module:

    if isinstance(t_module_latex_or_expr, tmc_types.T_Module):
        return t_module_latex_or_expr

    has_critical_error = False

    stripping_integers = False
    replacing_q_with_p_n = False
    further_process_q_k = False
    if phi_psi_options.strip_integer_exponents_for_symbolic_execution and phi_psi_options.execute_symbolically:
        further_process_q_k = True
        stripping_integers = True
    elif phi_psi_options.execute_symbolically == False and phi_psi_options.prime_characteristic >= 0 and phi_psi_options.extension_degree >= 0:
        further_process_q_k = True
        replacing_q_with_p_n = True

    ALLOWED_SYMBOLS_IN_INPUT_IF_REPLACING_Q_WITH_P_N = [THETA_SYMBOL, TAU_SYMBOL, A_SYMBOL]


    d = 0
    max_deg = -1
    expanded_matrices = []
    t_module_matrix = None
    t_module_expression = None
    t_module_expression_simplified = None
    parse_errors = []

    if t_module_latex_or_expr:
        try:
            # print(f"module_name={module_name} parsing...")
            t_module_expression = safe_parse_with_symbolic_d(t_module_latex_or_expr) if isinstance(t_module_latex_or_expr, str) else t_module_latex_or_expr
            

            try:
                t_module_expression_simplified = t_module_expression.doit()

                # print(f"module_name={module_name}, t_module_expression_simplified={t_module_expression_simplified}")

                t_module_matrix = None

                if isinstance(t_module_expression_simplified, sympy.Matrix):
                    t_module_matrix = t_module_expression_simplified
                    d = t_module_matrix.rows

                    if t_module_matrix.rows == t_module_matrix.cols:
                        pass
                    else:
                        parse_errors.append(f"{module_name} matrix must be square")
                else:
                    d = 1
                    t_module_matrix = sympy.Matrix([[t_module_expression_simplified]])
                    
                t_module_matrix_expanded = t_module_matrix.expand()

                # print(f"module_name={module_name}, t_module_matrix_expanded={t_module_matrix_expanded}")

                expanded_matrices = tmc_utils.expand_t_module_matrix(t_module_matrix_expanded)

                # print(f"module_name={module_name}, expanded_matrices={expanded_matrices}")

                max_deg = len(expanded_matrices)-1

                k = 0
                expanded_matrices_times_tau = []

                for ith_matrix in expanded_matrices:

                    expanded_matrices_times_tau.append(ith_matrix * TAU_SYMBOL_NON_COMM**k)

                    free_symbols_in_matrix = ith_matrix.free_symbols

                    for ith_free_symbol in free_symbols_in_matrix:
                        if replacing_q_with_p_n == True:
                            if ith_free_symbol not in ALLOWED_SYMBOLS_IN_INPUT_IF_REPLACING_Q_WITH_P_N:
                                has_critical_error = True
                                parse_errors.append(f"{module_name} can only contain the following symbols, because p and n are defined: {ALLOWED_SYMBOLS_IN_INPUT_IF_REPLACING_Q_WITH_P_N}")
                        elif ith_free_symbol in DISALLOWED_SYMBOLS_IN_INPUT:
                            has_critical_error = True
                            parse_errors.append(f"{module_name} cannot contain any of the following symbols: {DISALLOWED_SYMBOLS_IN_INPUT}")
                            break

                    if k == 0:
                        if tmc_utils.is_valid_first_t_module_matrix(ith_matrix):
                            pass
                        else:
                            parse_errors.append(module_name + " must be nilpotent, i.e. <math-span>(M_0 - \\theta \\cdot I)^d = 0</math-span>")

                    k += 1

                t_module_matrix = sympy.MatAdd(*expanded_matrices_times_tau).doit()

                # DO NOT DO THIS. The reason we expand the matrices and resimplify them, is to make sure tau
                # is on the right, without twists or anything.
                # t_module_matrix = t_module_matrix_expanded

            except Exception as e:
                has_critical_error = True
                parse_errors.append(f"Could not evaluate expression for {module_name}: {e}")
                print(f"{e}")
        except Exception as e:
            has_critical_error = True
            parse_errors.append(f"Could not parse latex for {module_name}: {e}")
    else:
        has_critical_error = True
        parse_errors.append(f"{module_name} cannot be empty")

    t_module_matrix_as_list = None
    t_module_as_ai_input = None
    if t_module_matrix is not None:
        t_module_matrix_as_2d_list = t_module_matrix.tolist()

        t_module_as_ai_input = tmc_ai_bridge.convert_2d_list_to_matrix(t_module_matrix_as_2d_list)
        t_module_as_ai_input = tmc_utils.toggle_commutativity_for_all_symbols_in_matrix(t_module_as_ai_input, change_to_commutative=True, keep_tau_non_comm=False)

    return tmc_types.T_Module(d, max_deg, expanded_matrices, t_module_matrix, t_module_matrix_as_2d_list, t_module_as_ai_input, t_module_expression, t_module_expression_simplified, parse_errors, has_critical_error)
