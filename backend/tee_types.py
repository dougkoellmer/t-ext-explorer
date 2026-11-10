
import sympy
from typing import NamedTuple, Any, Optional
from enum import Enum
import threading
from collections import OrderedDict

class PhiPsiOptions(NamedTuple):
    show_reduction_steps: bool
    replace_c_with_tau: bool
    replace_q_with_integer: bool
    group_coefficients: bool
    strip_integer_exponents_for_symbolic_execution: bool
    execute_symbolically: bool
    prime_characteristic: int # the p in p^n
    extension_degree: int # n in p^n
    reduce_log_output: bool = False
    only_generate_ext: bool = False
    reduce_latex_output: bool = False
    calc_leaf_and_byte_counts: bool = True

class MaxDegreeOfTauInMatrix(NamedTuple):
    max_degree: int
    zero_based_row: int
    zero_based_col: int
    zero_based_element_ordinal: int
    coefficient: sympy.Expr

class MatrixWithEffectiveElement(NamedTuple):
    raw: sympy.Matrix
    effective_zero_based_row: int
    effective_zero_based_col: int
    degree_of_tau: int
    zero_based_element_ordinal: int
    max_deg_tuple: MaxDegreeOfTauInMatrix

class PhiPsiSubvariant(Enum):
    ORIGINAL_LOGIC = 1
    CONSTANT_TERM_LOGIC = 2
    EMBEDDING_LOGIC = 3

class Variant(Enum):
    VARIANT__PHI_PSI_EMBEDDING = 1
    VARIANT__DIAGONAL = 2
    VARIANT__RIGHTMOST_BASIS = 3

class T_Module(NamedTuple):
    d: int
    max_deg: int
    expanded_matrices: list[sympy.Matrix]
    simplified: sympy.Matrix
    simplified_as_2d_list: list
    as_ai_input: sympy.Matrix
    expression: sympy.Expr
    expression_simplified: sympy.Expr
    parse_errors: list[str]
    has_critical_error: bool

class PhiPsiInput(NamedTuple):
    Phi_module: T_Module
    Psi_module: T_Module
    input_errors: list[str]
    is_max_degree_error: bool = False

class MaxDegreeOfTauInExpression(NamedTuple):
    max_degree: int
    coefficient: sympy.Expr

class Eij_And_Aij_Matrices(NamedTuple):
    Eij_row_count: int
    Eij_col_count: int
    Eij_matrices: list[MatrixWithEffectiveElement]
    Aij_matrices: list[sympy.Matrix]
    Aij_Phi_minus_Psi_Aij_matrices: list[sympy.Matrix]
    Nij_list: list[int]
    contains_constant_term_list: list[bool]

    # For debugging...
    Aij_Phi_matrices: list[sympy.Matrix]
    Psi_Aij_matrices: list[sympy.Matrix]

    Eij_tau_k_matrices_wrapper: list[list[MatrixWithEffectiveElement]]
    Eij_c_tau_k_matrices_wrapper: list[list[MatrixWithEffectiveElement]]
    Psi_Eij_c_tau_k_matrices_wrapper: list[list[MatrixWithEffectiveElement]]

class Psi_Eij_c_tau_k_decomposition(NamedTuple):
    row_index: int
    col_index: int
    degree_of_tau: int
    sum_of_coefficients: sympy.Expr

class MatrixSubsetCol(NamedTuple):
    matrix_index: int
    col_index: int
    vector: sympy.Matrix

class MatrixSubset(NamedTuple):
    raw: sympy.Matrix
    col_indeces: list[MatrixSubsetCol]
    ordinal: int

    def col_indeces_raw(self) -> list[int]:
        raw = []
        for ith in self.col_indeces:
            raw.append(ith.col_index)

        return raw

class MatrixFlagElem(NamedTuple):
    raw: sympy.Matrix
    matrices: list[sympy.Matrix]
    rank: int

class RightmostBasis(NamedTuple):
    t_module_matrix_subset: list[sympy.Matrix]
    column_lists: list[MatrixSubset]
    sorted_column_lists: list[MatrixSubset]
    flag: list[MatrixFlagElem]
    basis_subset: MatrixSubset
    basis_matrix: sympy.Matrix
    basis_matrix_inverse: sympy.Matrix
    basis_matrix_inverse_non_comm: sympy.Matrix
    last_m: int
    total_seconds_to_generate: float

class InitialAnalysis(NamedTuple):
    rightmost_basis: RightmostBasis
    Eij_and_Aij_matrices: Eij_And_Aij_Matrices
    total_seconds_to_generate: float

class R_Matrix(NamedTuple):
    raw: sympy.Matrix
    Eij_tau_k_matrices: list[MatrixWithEffectiveElement]
    Psi_Eij_c_tau_k_matrices: list[MatrixWithEffectiveElement]
    reduction_matrices: list[str]
    final_reduction_matrices: list[sympy.Matrix]
    degree_k_start_values: list[int]
    degree_k_finish_values: list[int]
    reduction_counts: list[int]
    leaf_count: int
    byte_count: int
    total_reduction_count: int
    total_seconds_to_generate: float

class PhiPsiReductionOutput(NamedTuple):
    Ext_matrix: R_Matrix
    V_matrix: R_Matrix
    W_V_matrix: R_Matrix
    response_dict: dict
    found_rightmost_basis: bool

class PhiPsiReductionOutputFromAiCode(NamedTuple):
    Ext_Matrix: R_Matrix
    Ext0_Matrix: R_Matrix
    K_Matrix: R_Matrix
    R_Matrix: R_Matrix
    P_Matrix: R_Matrix
    S_Matrix: R_Matrix
    response_dict: dict

class ThreadSafeLRUCache:
    def __init__(self, capacity: int):
        if capacity <= 0:
            raise ValueError("Capacity must be positive")
        self.capacity = capacity
        self.cache = OrderedDict()
        # A reentrant lock allows the thread holding the lock to re-acquire it
        self.lock = threading.RLock()

    def get(self, key: Any) -> Optional[Any]:
        """Retrieve an item and move it to the end to track recency."""
        with self.lock:
            if key not in self.cache:
                return None
            # Move key to the end (marked as Most Recently Used)
            self.cache.move_to_end(key)
            return self.cache[key]

    def put(self, key: Any, value: Any) -> None:
        """Insert or update an item, evicting the oldest if over capacity."""
        with self.lock:
            if key in self.cache:
                # Update existing key and mark as MRU
                self.cache.move_to_end(key)
            self.cache[key] = value
            
            # Evict the oldest item if capacity is exceeded
            if len(self.cache) > self.capacity:
                # last=False pops the first item (Least Recently Used)
                self.cache.popitem(last=False)

    def clear(self) -> None:
        """Thread-safe clearance of all cache data."""
        with self.lock:
            self.cache.clear()

# class ParenthesizedInteger(sympy.Function):
#     def _latex(self, printer):
#         return r"\left({}\right)".format(printer._print(self.args[0]))


class ParenthesizedInteger(sympy.Function):
    def _latex(self, printer, **kwargs):
        # 1. Print the inner argument wrapped in parentheses
        base_latex = r"\left({}\right)".format(printer._print(self.args[0]))
        
        # 2. Check if SymPy passed an exponent to append
        if 'exp' in kwargs and kwargs['exp'] is not None:
            expr_latex = printer._print(kwargs['exp'])
            return f"{base_latex}^{{{expr_latex}}}"
            
        return base_latex