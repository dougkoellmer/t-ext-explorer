
import sympy

import sys

from tee_consts import *
from tee_types import *
from itertools import combinations
import time
import tee_consts
import tee_view

import threading
import json

thread_local = threading.local()

import hashlib

import socket

def non_comm_number(number):
     return sympy.Symbol(str(number), commutative=False)

# This probably won't match exactly what Mathematica returns, but is at least
# the same general idea.
def mathematica_matrix_byte_count(matrix):
    # Start with the memory of the matrix container shell itself
    total_bytes = sys.getsizeof(matrix)
    
    # Add the memory of the flat internal list container holding the data
    # Use M.flat() instead of the deprecated M._mat
    total_bytes += sys.getsizeof(matrix.flat())
    
    # Optionally track of object IDs we've already counted to avoid double-counting 
    # shared references (like duplicate 0s, 1s, and symbols).
    # For now, keep this disabled to better-match Mathematica's approach.
    account_for_shared_memory = False
    seen_ids = set() if account_for_shared_memory == True else None
    
    # Recursively traverse every unique sub-expression in the matrix
    for element in matrix.flat():
        for node in sympy.preorder_traversal(element):
            if seen_ids is not None and id(node) not in seen_ids:
                seen_ids.add(id(node))
                total_bytes += sys.getsizeof(node)
                
    return total_bytes

# This should match Mathematica output for LeafCount, but don't rely on it absolutely.
# It was tested on a few examples and things matched, but again no absolute guarantees here.
def mathematica_leaf_count(expr):
    # If it's a matrix, convert it to a nested Python list first
    if hasattr(expr, 'tolist'):
        expr = expr.tolist()
        
    if isinstance(expr, list):
        # 1 leaf for the 'head' (the list itself) + the leaves of its items
        return 1 + sum(mathematica_leaf_count(item) for item in expr)
    
    # For standalone SymPy expressions, preorder_traversal behaves exactly like LeafCount
    return sum(1 for _ in sympy.preorder_traversal(expr))

def is_valid_first_t_module_matrix(matrix):
    """
    Checks if a matrix M can be the first matrix of a t-module.
    Condition: M - theta * I must be nilpotent.
    """
    if not matrix.is_square:
        return False
        
    d = matrix.shape[0]
    I = sympy.eye(d)
    
    nilpotent_test = (matrix - THETA_SYMBOL * I)**d
    
    return nilpotent_test.expand().is_zero_matrix

def toggle_commutativity_for_all_symbols_in_matrix(matrix, change_to_commutative, keep_tau_non_comm=True):
    matrix = matrix.applyfunc(lambda matrix_elem: toggle_commutativity_for_all_symbols(matrix_elem, change_to_commutative=change_to_commutative, is_inside_power=False, keep_tau_non_comm=keep_tau_non_comm))

    return matrix

def toggle_commutativity_for_all_symbols(expression, change_to_commutative, is_inside_power=False, keep_tau_non_comm=True):
    if is_inside_power == False and isinstance(expression, sympy.Number) and change_to_commutative == False:
        if FEATURE_FLAG__APPLY_TWIST_TO_NUMBERS == True:
            return non_comm_number(expression)
        else:
            return expression
    elif isinstance(expression, sympy.Symbol):
        if change_to_commutative == True:
            if keep_tau_non_comm == True and expression == TAU_SYMBOL_NON_COMM:
                return expression # An exception is made for TAU_SYMBOL_NON_COMM. It cannot be made commutative.
            else:
                if FEATURE_FLAG__APPLY_TWIST_TO_NUMBERS == True:
                    try:
                        # Attempt to parse the symbol's name back into a real SymPy number, which is always commutative
                        return sympy.sympify(expression.name)
                    except (sympy.SympifyError, ValueError):
                        return sympy.symbols(expression.name, commutative=change_to_commutative)
                else:
                    return sympy.symbols(expression.name, commutative=change_to_commutative)
        else:
            return sympy.symbols(expression.name, commutative=change_to_commutative)
    elif isinstance(expression, sympy.Pow):
        base, exp = expression.as_base_exp()
        return toggle_commutativity_for_all_symbols(base, change_to_commutative, is_inside_power, keep_tau_non_comm)**exp
    elif isinstance(expression, sympy.Mul):
        expression_factors = expression.as_ordered_factors()
        new_factors = []

        for ith_factor in expression_factors:
            ith_new_factor = toggle_commutativity_for_all_symbols(ith_factor, change_to_commutative, is_inside_power, keep_tau_non_comm)
            
            new_factors.append(ith_new_factor)
        
        return sympy.Mul(*new_factors, evaluate=True)

    elif isinstance(expression, sympy.Add):
        expression_terms = expression.as_ordered_terms()
        new_terms = []
        for ith_term in expression_terms:
            ith_new_term = toggle_commutativity_for_all_symbols(ith_term, change_to_commutative, is_inside_power, keep_tau_non_comm)
            
            new_terms.append(ith_new_term)
        
        return sum(new_terms)
    else:
        # For all other types of expressions, return the original expression.
        return expression

def check_connection_still_alive():
    if FEATURE_FLAG__RUN_UNIT_TESTS:
        return

    current_epoch_seconds = int(time.time())

    request_timeout = tee_consts.get_request_timeout_in_secs(thread_local.response.path, thread_local.is_localhost)

    if request_timeout > 0:
        seconds_since_processing_started = current_epoch_seconds - thread_local.processing_start_time
        if seconds_since_processing_started > request_timeout:
            raise ValueError(f"Exceeded maximum processing time of {request_timeout} seconds. You may consider running the server locally on your machine, where there is no restriction.")

    if FEATURE_FLAG__USE_JSON_CHUNKING:
        last_check_time = getattr(thread_local, "last_check_time", 0)
        seconds_since_last_check = current_epoch_seconds - last_check_time

        if seconds_since_last_check > MINIMUM_SECONDS_BETWEEN_CHUNKED_CONNECTION_LIVENESS_CHECKS:
            print("Checking connection still alive.")
            response = thread_local.response
            response.wfile.write(b"1\r\n\n\r\n")
            # Triggers BrokenPipeError if user closes page or navigates away.
            # NOTE: This can take maybe a minute to actually throw an exception after page is closed.
            # Better late than never!
            response.wfile.flush()
            print("Connection still alive.")

            thread_local.last_check_time = current_epoch_seconds

def get_final_json_chunk_bytes(json_dict):
    json_bytes = json.dumps(json_dict).encode('utf-8')
        
    return json_bytes

def expand_t_module_matrix(t_module: sympy.Matrix):
    t_module = t_module.applyfunc(lambda matrix_elem: toggle_commutativity_for_all_symbols(matrix_elem, change_to_commutative=False))

    max_deg_tuple = get_max_degree_of_tau_in_matrix(t_module, require_tau_to_be_rightmost_factor=False)
    max_deg = max_deg_tuple.max_degree

    to_return = []

    t_module = t_module.applyfunc(lambda matrix_elem: toggle_commutativity_for_all_symbols(matrix_elem, change_to_commutative=True, keep_tau_non_comm=False))

    for deg in range(max_deg + 1):
        # Extract the coefficient matrix for tau**deg
        coeff_mat = t_module.applyfunc(lambda x: sympy.Poly(x, TAU_SYMBOL).coeff_monomial(TAU_SYMBOL**deg))
        
        to_return.append(coeff_mat)
        
    return to_return

def write_final_json_chunk(response, json_bytes, use_json_chunking=FEATURE_FLAG__USE_JSON_CHUNKING):
    if use_json_chunking:

        hex_len = f"{len(json_bytes):X}\r\n".encode('utf-8')
        response.wfile.write(hex_len + json_bytes + b"\r\n")
        
        response.wfile.write(b"0\r\n\r\n")
        response.wfile.flush()
    else:
        response.wfile.write(json_bytes)

def multiply_matrices_non_comm(matrix1, matrix2, finish=True, move_tau=True):
    matrix1 = matrix1.applyfunc(lambda matrix_elem: toggle_commutativity_for_all_symbols(matrix_elem, change_to_commutative=False))
    matrix2 = matrix2.applyfunc(lambda matrix_elem: toggle_commutativity_for_all_symbols(matrix_elem, change_to_commutative=False))

    resulting_matrix = matrix1 * matrix2;

    resulting_matrix = resulting_matrix.expand()

    if finish == True:
        resulting_matrix = finish_multiply_matrices_non_comm(resulting_matrix, move_tau)

    return resulting_matrix

def finish_multiply_matrices_non_comm(resulting_matrix, move_tau = True):
    if move_tau == True:
        resulting_matrix = resulting_matrix.applyfunc(move_tau_to_the_right_of_other_symbols)

    resulting_matrix = resulting_matrix.applyfunc(lambda matrix_elem: toggle_commutativity_for_all_symbols(matrix_elem, change_to_commutative=True))
    resulting_matrix = convert_matrix_to_simple_terms(resulting_matrix)

    return resulting_matrix


def get_degree_of_tau(tau_symbol_or_tau_to_some_degree):
    if tau_symbol_or_tau_to_some_degree == TAU_SYMBOL_NON_COMM:
        return 1
    elif isinstance(tau_symbol_or_tau_to_some_degree, sympy.Pow):
        base, exp = tau_symbol_or_tau_to_some_degree.as_base_exp()

        if base == TAU_SYMBOL_NON_COMM:
            return exp
    
    return -1


def move_tau_to_the_right_without_twist(expr, var = TAU_SYMBOL_NON_COMM):
    if expr.is_Mul:
        fixed_args = [move_tau_to_the_right_without_twist(arg, var) for arg in expr.args]
        
        others = [f for f in fixed_args if not f.has(var)]
        x_terms = [f for f in fixed_args if f.has(var)]
        
        return sympy.Mul(*(others + x_terms), evaluate=False)
        
    elif expr.is_Add:
        return sympy.Add(*[move_tau_to_the_right_without_twist(arg, var) for arg in expr.args], evaluate=False)
        
    elif expr.is_Pow:
        return sympy.Pow(move_tau_to_the_right_without_twist(expr.base, var), expr.exp, evaluate=False)
        
    return expr

def move_tau_to_the_right_of_other_symbols_in_multiply_expression(expression):
    if isinstance(expression, sympy.Mul) == False:
        return expression

    found_something_to_move = True

    while found_something_to_move:

        found_something_to_move = False

        expression_factors = expression.as_ordered_factors()

        for i, ith_arg in enumerate(expression_factors):
            
            degree_of_tau = get_degree_of_tau(ith_arg)

            if degree_of_tau > TAU_DEGREE_OF_ZERO:
                if i < len(expression_factors)-1 :
                    next_arg = expression_factors[i+1]
                    degree_of_tau_for_next_arg = get_degree_of_tau(next_arg)

                    if degree_of_tau_for_next_arg <= TAU_DEGREE_OF_ZERO:

                        next_non_tau_arg = next_arg
                        next_non_tau_with_twist = next_non_tau_arg

                        if isinstance(next_non_tau_arg, sympy.Pow):
                            first_base, first_exp = next_non_tau_arg.as_base_exp()

                            next_non_tau_with_twist = first_base**(first_exp * Q_SYMBOL**degree_of_tau)
                        else:
                            next_non_tau_with_twist = next_non_tau_arg**Q_SYMBOL**degree_of_tau
                            
                        
                        expression_factors[i] = next_non_tau_with_twist
                        expression_factors[i+1] = ith_arg
                    else: # means next_arg is also TAU_SYMBOL_NON_COMM
                        expression_factors[i] = TAU_SYMBOL_NON_COMM**(degree_of_tau + degree_of_tau_for_next_arg)
                        del expression_factors[i+1]

                    expression = sympy.Mul(*expression_factors, evaluate=True)
                    found_something_to_move = True
                    break

    return expression


def expand_frobenius_2(expr, q_symbol):
    """
    Distributes any exponent that is a power of q (e.g., q, q^2, q^3) 
    over terms inside an addition: (A + B)^(q^k) -> A^(q^k) + B^(q^k)
    """
    def is_q_power(exponent):
        # Matches q itself, or q**k
        if exponent == q_symbol:
            return True
        if exponent.is_Pow and exponent.base == q_symbol:
            return True
        return False

    query = lambda x: x.is_Pow and x.base.is_Add and is_q_power(x.exp)
    
    transform = lambda x: sympy.Add(*[term**x.exp for term in x.base.args])
    
    return expr.replace(query, transform)

# Manually distributes an outer exponent over a sum, 
# and applies standard power rules to simplify the nested powers.
def expand_frobenius(expression, power):

    # NOTE: Initial implementation here using expression.replace proved the concept but can't handle things like nested parentheses.
    # return expression.replace(
    #     lambda expr: expr.is_Pow and expr.exp == power and expr.base.is_Add,
    #     lambda expr: sympy.Add(*[term**power for term in expr.base.args])
    # )

    def worker(e):
        if isinstance(e, sympy.Pow) and e.exp == power and isinstance(e.base, sympy.Add):
            # Distribute the power over each term in the sum
            new_terms = []
            for term in e.base.args:
                # Apply the power to each factor inside the term, then simplify rules like (x^a)^b = x^(a*b)
                distributed_term = sympy.powdenest(term**power, force=True)
                new_terms.append(distributed_term)
            return sympy.Add(*new_terms)
        return e

    return expression.replace(lambda e: isinstance(e, sympy.Pow) and e.exp == power, worker)


def move_tau_to_the_right_of_other_symbols(expression):
    if isinstance(expression, sympy.Symbol) or isinstance(expression, sympy.Pow) or isinstance(expression, sympy.Number):
        return expression
    elif isinstance(expression, sympy.Mul):
        return move_tau_to_the_right_of_other_symbols_in_multiply_expression(expression)
    elif isinstance(expression, sympy.Add):
        new_terms = []
        expression_terms = expression.as_ordered_terms()
        for ith_term in expression_terms:
            ith_new_term = move_tau_to_the_right_of_other_symbols(ith_term)
            
            new_terms.append(ith_new_term)
        
        return sum(new_terms)
    else:
        # For all other types of expressions, return the original expression.
        return expression

def is_positive_int(value):
    is_pos_sympy_int = value.is_integer and value.is_positive

    if is_pos_sympy_int:
        return True
    elif isinstance(value, int) and type(value) is not bool and value > 0:
        return True
    else:
        return False


def replace_q_power_with_integer(expression, is_inside_power = False, simplifying_non_symbolic_bases = False, base_is_symbolic_if_inside_power=True):
    
    if hasattr(expression, "replace"):
        # print(f"ERERE1 non_symbolic={simplifying_non_symbolic_bases}: {expression}")

        new_expr = expression

        def transform_exponent_tree(exponent_expr):
            step1 = exponent_expr.replace(
                lambda e: e.is_Pow and e.base == Q_SYMBOL,
                lambda e: ParenthesizedInteger(e.exp)
            )
            
            step2 = step1.replace(
                Q_SYMBOL,
                ParenthesizedInteger(1)
            )
            
            return step2

        def transform_exponent_tree_strip_ints(exponent_expr):
            step1 = exponent_expr.replace(
                lambda e: e.is_Pow and e.base == Q_SYMBOL,
                lambda e: 1
            )
            
            if hasattr(step1, "replace"):
                step2 = step1.replace(
                    Q_SYMBOL,
                    1
                )
            else:
                step2 = step1
            
            return step2

        def outer_transformer(x):
            if not x.base.free_symbols and simplifying_non_symbolic_bases == True: # ( (x.is_Pow and x.base == Q_SYMBOL) or x == Q_SYMBOL):
                new_exponent = transform_exponent_tree_strip_ints(x.exp)
                return x.base ** new_exponent
            else:
                if simplifying_non_symbolic_bases == True:
                    return x
                else:
                    # Send the exponent tree to our sequential transformer
                    new_exponent = transform_exponent_tree(x.exp)
                    return x.base ** new_exponent
            
        new_expr = new_expr.replace(
            lambda x: x.is_Pow and x.base != C_SYMBOL,
            outer_transformer
        )

        # print(f"ERERE2 non_symbolic={simplifying_non_symbolic_bases}: {new_expr}")

        return new_expr
    else:
        return expression

    # if isinstance(expression, sympy.Symbol):
    #     if is_inside_power and expression == Q_SYMBOL:
    #         if simplifying_non_symbolic_bases == True:
    #             if base_is_symbolic_if_inside_power == False:
    #                 return 1
    #             else:
    #                 return expression
    #         else:
    #             return ParenthesizedInteger(1)
    #     else:
    #         return expression
    # elif isinstance(expression, sympy.Pow):
    #     base_1, exp_1 = expression.as_base_exp()

    #     if is_inside_power:
    #         if base_1 == Q_SYMBOL:
    #             if simplifying_non_symbolic_bases == True:
    #                 if base_is_symbolic_if_inside_power == False:
    #                     return 1
    #                 else:
    #                     return expression
    #             else:
    #                 if is_positive_int(exp_1):
    #                     return ParenthesizedInteger(exp_1)
    #                 else:
    #                     # Probably shouldn't be able to get here.
    #                     return expression
    #         else:
    #             return expression
    #     else:
    #         if base_1 == C_SYMBOL:
    #             return expression
    #         else:
    #             if exp_1 == Q_SYMBOL:
    #                 if simplifying_non_symbolic_bases == True:
    #                     if not base_1.free_symbols:
    #                         if exp_1.is_negative:
    #                             return 1 / base_1
    #                         else:
    #                             return base_1
    #                     else:
    #                         return expression
    #                 else:
    #                     return base_1**ParenthesizedInteger(1)
    #             elif isinstance(exp_1, sympy.Pow):
    #                 base_2, exp_2 = exp_1.as_base_exp()

    #                 if base_2 == Q_SYMBOL:
    #                     if simplifying_non_symbolic_bases == True:
    #                         if not base_1.free_symbols:
    #                             if base_2.is_negative:
    #                                 return 1 / base_1
    #                             else:
    #                                 return base_1
    #                         else:
    #                             return expression
    #                     else:
    #                         if is_positive_int(exp_2):
    #                             return base_1**ParenthesizedInteger(exp_2)
    #                         else:
    #                             # Probably shouldn't be able to get here.
    #                             return expression
    #                 else:
    #                     # Probably shouldn't be able to get here.
    #                     return expression
    #             else:
    #                 if not base_1.free_symbols:
    #                     base_is_symbolic_if_inside_power = False
    #                 return base_1**replace_q_power_with_integer(exp_1, is_inside_power=True, simplifying_non_symbolic_bases=simplifying_non_symbolic_bases, base_is_symbolic_if_inside_power=base_is_symbolic_if_inside_power)
    # elif isinstance(expression, sympy.Mul):
    #     new_factors = []
    #     expression_factors = expression.as_ordered_factors()
    #     for ith_factor in expression_factors:
    #         ith_new_factor = replace_q_power_with_integer(ith_factor, is_inside_power, simplifying_non_symbolic_bases=simplifying_non_symbolic_bases, base_is_symbolic_if_inside_power=base_is_symbolic_if_inside_power)
            
    #         new_factors.append(ith_new_factor)
        
    #     return sympy.Mul(*new_factors, evaluate=False)
    # elif isinstance(expression, sympy.Add):
    #     new_terms = []
    #     expression_terms = expression.as_ordered_terms()
    #     for ith_term in expression_terms:
    #         ith_new_term = replace_q_power_with_integer(ith_term, is_inside_power, simplifying_non_symbolic_bases=simplifying_non_symbolic_bases, base_is_symbolic_if_inside_power=base_is_symbolic_if_inside_power)
            
    #         new_terms.append(ith_new_term)
        
    #     return sum(new_terms)
    # else:
    #     # For all other types of expressions, return the original expression.
    #     return expression

def replace_c_with_tau(expression):
    # print(f"replace_c_with_tau={expression}")
    if isinstance(expression, sympy.Symbol):
        # print(f"replace_c_with_tau__symbol={expression}")
        if expression == C_SYMBOL:
            # print(f"replace_c_with_tau__symbol__1")
            # c^q^0 = tau^0 = 1
            return 1
        else:
            return expression
    elif isinstance(expression, sympy.Pow):
        base_1, exp_1 = expression.as_base_exp()

        if base_1 == C_SYMBOL:
            if exp_1 == Q_SYMBOL:
                # c^q^1 = tau^1 = tau
                return TAU_SYMBOL_NON_COMM
            elif isinstance(exp_1, sympy.Pow):
                base_2, exp_2 = exp_1.as_base_exp()

                if base_2 == Q_SYMBOL:
                    # print(f"replace_c_with_tau__c^q^k={expression}")
                    if is_positive_int(exp_2):
                        # print(f"replace_c_with_tau__c^q^k__doing_replacement")
                        return TAU_SYMBOL_NON_COMM**exp_2
                    else:
                        # Probably shouldn't be able to get here.
                        return expression
                else:
                    # Probably shouldn't be able to get here.
                    return expression
            else:
                # Probably shouldn't be able to get here.
                return expression
        else:
            return expression
    elif isinstance(expression, sympy.Mul):
        new_factors = []
        expression_factors = expression.as_ordered_factors()
        for ith_factor in expression_factors:
            ith_new_factor = replace_c_with_tau(ith_factor)
            
            new_factors.append(ith_new_factor)
        
        return sympy.Mul(*new_factors, evaluate=True)
    elif isinstance(expression, sympy.Add):
        new_terms = []
        expression_terms = expression.as_ordered_terms()
        for ith_term in expression_terms:
            ith_new_term = replace_c_with_tau(ith_term)
            
            new_terms.append(ith_new_term)
        
        return sum(new_terms)
    else:
        # For all other types of expressions, return the original expression.
        return expression


def get_coefficient_of_factor_with_tau_on_end(expression_factor, tau_degree):
    tau_to_the_pow = TAU_SYMBOL_NON_COMM**tau_degree
    return expression_factor.subs(tau_to_the_pow, 1)

def get_all_degrees_of_tau_in_expression(expression, require_tau_to_be_rightmost_factor = True):
    to_return = []

    if expression == 0:
        pass
    elif isinstance(expression, sympy.Add):
        expression_terms = expression.as_ordered_terms()

        for ith_term in expression_terms:
            to_return = to_return + get_all_degrees_of_tau_in_expression(ith_term, require_tau_to_be_rightmost_factor)
    elif isinstance(expression, sympy.Mul):
        # For this case in general, previous implementation simply looked at the last factor in the expression to get the degree of tau.
        # Looking at it later, more robust is to look through each factor, for tau. However this function as a whole pretty much assumes
        # that tau is already on the end, in order to get the coefficient. So an error is thrown if that's not the case.
        expression_factors = expression.as_ordered_factors()

        highest_tau_degree_found_in_expression = TAU_DEGREE_OF_ZERO

        last_factor = expression_factors[-1] # get last factor
        
        for ith_factor in expression_factors:
            ith_factor_degree = get_degree_of_tau(ith_factor)

            if ith_factor_degree > highest_tau_degree_found_in_expression:
                highest_tau_degree_found_in_expression = ith_factor_degree

                if ith_factor == last_factor:
                    pass # all good
                else:
                    if require_tau_to_be_rightmost_factor == True:
                        raise ValueError(f"Expected tau factor {ith_factor} to be the last factor in the expression {expression}.");

        if highest_tau_degree_found_in_expression > TAU_DEGREE_OF_ZERO:
            coefficient = get_coefficient_of_factor_with_tau_on_end(expression, highest_tau_degree_found_in_expression)
            to_return.append(MaxDegreeOfTauInExpression(highest_tau_degree_found_in_expression, coefficient))
        else:
            to_return.append(MaxDegreeOfTauInExpression(TAU_DEGREE_OF_ZERO, expression))
    else:
        degree_of_tau = get_degree_of_tau(expression)
        if degree_of_tau > TAU_DEGREE_OF_ZERO:
            to_return.append(MaxDegreeOfTauInExpression(degree_of_tau, COEFFICIENT_OF_ONE))
        else:
            to_return.append(MaxDegreeOfTauInExpression(TAU_DEGREE_OF_ZERO, expression))
            
    return to_return

def does_matrix_contain_any_zero_columns(matrix:sympy.Matrix) -> bool:
    return any(matrix.col(i).is_zero_matrix for i in range(matrix.shape[1]))

def list_of_matrix_subset_col_indeces_to_ints(indeces:list[MatrixSubsetCol]) -> list[int]:
    raw = []
    for ith in indeces:
        raw.append(ith.col_index)

    return raw

def get_matrix_from_matrix_subset(t_module_matrices: list[sympy.Matrix], matrix_subset:MatrixSubset) -> sympy.Matrix:
    d = t_module_matrices[0].cols
    to_return = sympy.zeros(d, d)

    for ith in matrix_subset.col_indeces:
        matrix_from_t_module = t_module_matrices[ith.matrix_index]

        to_return[:, ith.col_index] = matrix_from_t_module[:, ith.col_index]

    return to_return

def get_max_degree_of_tau_in_expression(expression, require_tau_to_be_rightmost_factor = True):

    all_degrees = get_all_degrees_of_tau_in_expression(expression, require_tau_to_be_rightmost_factor)

    if len(all_degrees) > 0:
        highest:MaxDegreeOfTauInExpression = None
        for ith in all_degrees:
            # print(f"ith_degree.max_degree={ith.max_degree}, ith_degree.coefficient={ith.coefficient}")

            if highest == None or ith.max_degree > highest.max_degree:
                highest = ith
        
        return highest
    else:
        return MaxDegreeOfTauInExpression(-1, None)

def group_coefficients_by_tau(pre_expanded_expression):

    # NOTE: Cannot use this simple approach right below for non-comm symbols, so have to make it more complex further below.
    # Obviously tau_powers array would be generated dynamically given a max degree. This was just a test.
    # tau_powers = [TAU_SYMBOL_NON_COMM, TAU_SYMBOL_NON_COMM**2, TAU_SYMBOL_NON_COMM**3, TAU_SYMBOL_NON_COMM**4]
    # grouped_expression = sympy.collect(expression, tau_powers)

    grouped_dict = {}

    # WARNING: callers are responsible for making sure the expression is already exanded.
    for term in sympy.Add.make_args(pre_expanded_expression):
        # args_cnc() splits terms into [commutative] and [non-commutative]
        c_part, nc_part = term.args_cnc()
        
        nc_key = sympy.Mul(*nc_part) if nc_part else 1
        c_coeff = sympy.Mul(*c_part) if c_part else 1
        
        grouped_dict[nc_key] = grouped_dict.get(nc_key, 0) + c_coeff

    # Reconstruct using evaluate=False to force coeff on the left and tau on the right
    collected_terms = []
    for power, coeff in grouped_dict.items():
        # coeff = sympy.cancel(coeff)

        if power == 1:
            # This is a pure constant term (e.g. 5)
            collected_terms.append(coeff)
        else:
            # evaluate=False guarantees 'coeff' is left, 'tau to the pow' is right
            collected_terms.append(sympy.Mul(coeff, power, evaluate=False))

    collected_expr = sympy.Add(*collected_terms)

    return collected_expr

def get_max_degree_of_tau_in_matrix(matrix, require_tau_to_be_rightmost_factor = True, subvariant = PhiPsiSubvariant.CONSTANT_TERM_LOGIC, Eij_and_Aij_matrices: Eij_And_Aij_Matrices = None, preferred_element_ordinal: int = -1):
    Nij_list = []
    contains_constant_term_list = []

    if Eij_and_Aij_matrices is not None:
        Nij_list = Eij_and_Aij_matrices.Nij_list
        contains_constant_term_list = Eij_and_Aij_matrices.contains_constant_term_list

    max_deg = TAU_DEGREE_OF_ZERO
    max_deg_row = -1
    max_deg_col = -1
    coefficient = -1
    element_ordinal = -1

    max_deg_delta = -1

    current_element_ordinal = 0
    for zero_based_row in range(matrix.rows):
        for zero_based_col in range(matrix.cols):
            element = matrix[zero_based_row, zero_based_col]
            max_deg_in_expression = get_max_degree_of_tau_in_expression(element, require_tau_to_be_rightmost_factor)

            if len(Nij_list) > 0:
                ith_Nij = Nij_list[current_element_ordinal]
                ith_contains_constant_term = contains_constant_term_list[current_element_ordinal]
                Nij_base = ith_Nij

                if ith_contains_constant_term and subvariant == PhiPsiSubvariant.CONSTANT_TERM_LOGIC:
                    Nij_base += 1

                delta = max_deg_in_expression.max_degree - Nij_base

                if delta >= 0 and max_deg_in_expression.max_degree > max_deg:
                # if delta > max_deg_delta:
                # if delta > max_deg_delta or (current_element_ordinal == preferred_element_ordinal and delta >= 0):
                    max_deg_delta = delta
                    max_deg = max_deg_in_expression.max_degree
                    coefficient = max_deg_in_expression.coefficient
                    max_deg_row = zero_based_row
                    max_deg_col = zero_based_col
                    element_ordinal = current_element_ordinal

                    # if current_element_ordinal == preferred_element_ordinal:
                    #     break
            else:
                if max_deg_in_expression.max_degree > max_deg:
                    max_deg = max_deg_in_expression.max_degree
                    coefficient = max_deg_in_expression.coefficient
                    max_deg_row = zero_based_row
                    max_deg_col = zero_based_col
                    element_ordinal = current_element_ordinal
            
            current_element_ordinal += 1

    return MaxDegreeOfTauInMatrix(max_deg, max_deg_row, max_deg_col, element_ordinal, coefficient)

def does_matrix_contain_constant_term(matrix):
    
    for zero_based_row in range(matrix.rows):
        for zero_based_col in range(matrix.cols):
            element = matrix[zero_based_row, zero_based_col]
            all_degrees_of_tau = get_all_degrees_of_tau_in_expression(element)

            for ith_degree_of_tau in all_degrees_of_tau:
                if ith_degree_of_tau.max_degree == 0:
                    return True
            
    return False

def sub_if_int_exponent(old_sym, new_val):
    def condition_and_replace(node):
        if isinstance(node, sympy.Pow):
            base, exp = node.args
            if isinstance(base, sympy.Integer) and exp.has(old_sym):
                new_exp = exp.subs(old_sym, new_val)
                return sympy.Pow(base, new_exp)
        return node
    return condition_and_replace

def simplify_each_term_in_expression(expression: sympy.Expr, strip_integer_exponents_for_symbolic_execution: bool = False, p = -1, n = -1) -> sympy.Expr:
    expression_terms = expression.as_ordered_terms()
    new_terms = []

    process_further = False
    
    if strip_integer_exponents_for_symbolic_execution == True:
        process_further = True
    elif p >= 0 or n >= 0:
        process_further = True

    for ith_term in expression_terms:
        try:
            tau: MaxDegreeOfTauInExpression = get_max_degree_of_tau_in_expression(ith_term)
            tau_to_the_pow = TAU_SYMBOL_NON_COMM**tau.max_degree
            ith_new_term = ith_term.subs(tau_to_the_pow, 1)

            if process_further == True:
                # print(f"before={ith_new_term}")

                def expand_term(term):
                    # print(f"before_sub={term}")
                    term = expand_frobenius_2(term, Q_SYMBOL)
                    term = sympy.expand_power_base(term)
                    # term = sympy.powsimp(term)
                    # term = term.expand(force=True, power_exp=False)
                    # print(f"mid_sub={term}")
                    # term = term.replace(lambda node: isinstance(node, sympy.Pow), lambda node: sub_if_int_exponent(Q_SYMBOL, 1)(node))
                    if strip_integer_exponents_for_symbolic_execution == True:
                        term = replace_q_power_with_integer(term, is_inside_power=False, simplifying_non_symbolic_bases=True)
                    # elif p >= 0 or n >= 0:
                    #     # print(f"here: {p}, {n}")
                    #     term = simplify_twisted_modular(term, p, n)
                    # print(f"after_sub={term}")

                    return term

                if isinstance(ith_new_term, sympy.Mul):
                    ith_new_term_factors = ith_new_term.as_ordered_factors()
                    ith_new_term_factors_new = []
                    for ith_new_term_factor in ith_new_term_factors:
                        ith_new_term_factor = expand_term(ith_new_term_factor)
                        ith_new_term_factors_new.append(ith_new_term_factor)
                    
                    ith_new_term = sympy.Mul(*ith_new_term_factors_new, evaluate=True)
                else:
                    ith_new_term = expand_term(ith_new_term)

                # print(f"after={ith_new_term}")
                

            ith_new_term = sympy.simplify(ith_new_term)
            # ith_new_term = sympy.factor(ith_term)
            ith_new_term = ith_new_term * tau_to_the_pow
        except Exception as e:
            print(f"WARNING: Could not simplify individual expression term in matrix: {ith_term}")
            print(f"e={e}")
            ith_new_term = ith_term
        
        new_terms.append(ith_new_term)
    
    new_expression = sum(new_terms)

    return new_expression

# def simplify_twisted_modular(expr, p_val, n_val):
#     sp = sympy
#     """
#     Simplifies twisted polynomial expressions over GF(p^n),
#     correctly expanding Frobenius powers (base^q) over additions and multiplications.
#     """
#     p_sym = sp.Symbol('p')
#     n_sym = sp.Symbol('n')
#     q_sym = sp.Symbol('q')
    
#     exponent_modulus = p_val**n_val - 1

#     def expand_frobenius(node):
#         """Recursively distributes exponents that match q (or p^k) over Add/Mul nodes."""
#         if not node.is_Pow:
#             return node
            
#         base, exponent = node.base, node.exp
        
#         # Evaluate if the exponent represents a Frobenius power (q = p^n)
#         exponent_eval = exponent.subs(q_sym, p_sym**n_sym).subs({p_sym: p_val, n_sym: n_val})
        
#         # If the exponent is a multiple of p or evaluates to q, apply Freshman's Dream
#         if exponent_eval % p_val == 0 or exponent_eval == p_val**n_val:
#             # Case 1: (A + B)^q = A^q + B^q (mod p)
#             if base.is_Add:
#                 return base.func(*[expand_frobenius(sp.Pow(arg, exponent)) for arg in base.args])
            
#             # Case 2: (A * B)^q = A^q * B^q
#             if base.is_Mul:
#                 return base.func(*[expand_frobenius(sp.Pow(arg, exponent)) for arg in base.args])
                
#         return node

#     def reduce_node(node):
#         """Performs standard modular exponent reduction and coefficient reduction."""
#         if node.is_Pow:
#             base = node.base
#             exponent = node.exp
#             exponent_subbed = exponent.subs(q_sym, p_sym**n_sym).subs({p_sym: p_val, n_sym: n_val})
            
#             # Reduce exponents modulo p^n - 1
#             mod_exponent = sp.Mod(exponent_subbed, exponent_modulus)
#             mod_exponent = sp.simplify(mod_exponent)
#             return base ** mod_exponent
            
#         elif node.is_Integer:
#             return sp.Mod(node, p_val)
            
#         return node

#     # Step 1: Push Frobenius exponents down to the individual variables/constants
#     frobenius_expanded = expr.bottom_up(expand_frobenius)
    
#     # Step 2: Apply modular reduction to the exponents and integer constants
#     processed_expr = frobenius_expanded.bottom_up(reduce_node)
    
#     # Step 3: Final cleanup for any residual structural 'q' variables
#     final_expr = processed_expr.subs(q_sym, p_val**n_val)
    
#     return sp.simplify(final_expr)
    
def simplify_twisted_modular(expr, p_val, n_val):
    """
    Simplifies twisted polynomial expressions over the finite field GF(p^n).
    Reduces tower exponents modulo (p^n - 1) and base coefficients modulo p.
    """
    p_sym = sympy.Symbol('p')
    n_sym = sympy.Symbol('n')
    q_sym = Q_SYMBOL
    
    exponent_modulus = p_val**n_val - 1

    def transform(node):
        if node.is_Pow:
            base = node.base
            exponent = node.exp
            
            # Substitute q = p^n, then substitute specific p and n values into the exponent
            exponent_subbed = exponent.subs(q_sym, p_sym**n_sym).subs({p_sym: p_val, n_sym: n_val})
            
            # Reduce exponent modulo (p^n - 1)
            mod_exponent = sympy.Mod(exponent_subbed, exponent_modulus)
            mod_exponent = sympy.simplify(mod_exponent)
            
            return base ** mod_exponent
        
        elif node.is_Integer:
            return sympy.Mod(node, p_val)
        
        return node

    if hasattr(expr, "replace"):
        processed_expr = expr.replace(lambda x: True, transform)
    else:
        processed_expr = expr
    
    if hasattr(processed_expr, "subs"):
        final_expr = processed_expr.subs(q_sym, p_val**n_val)
    else:
        final_expr = processed_expr
    
    return sympy.simplify(final_expr)

# WARNING: Should only use this method after ensuring that tau twist has been applied.
# E.g. move_tau_to_the_right_of_other_symbols
def convert_matrix_to_simple_terms(matrix:sympy.Matrix, strip_integer_exponents_for_symbolic_execution: bool = False, p = -1, n = -1):
    
    hit_error = False

    try:
        # matrix = sympy.expand(matrix)
        # matrix = sympy.simplify(matrix)

        matrix = matrix.applyfunc(lambda matrix_elem: simplify_each_term_in_expression(matrix_elem, strip_integer_exponents_for_symbolic_execution, p, n) )
    except Exception as e:
        print(f"WARNING: Could not simplify matrix, trying per-element simplification instead.")
        print(f"{tee_view.matrix_to_latex(matrix)}")
        hit_error = True

    # if hit_error == True:
    #     print(f"post-factored Matrix:")
    #     print(f"{tee_view.matrix_to_latex(matrix)}")
    matrix = sympy.expand(matrix)

    return matrix

def unique_hash(input_string):
    encoded_string = input_string.encode('utf-8')
    hash_object = hashlib.sha256(encoded_string)
    
    return hash_object.hexdigest()

def is_running_on_localhost(request):
    bound_address = request.server.server_address[0]
        
    is_bound_locally = bound_address in ['localhost', '127.0.0.1', '::1', '0.0.0.0']

    return is_bound_locally