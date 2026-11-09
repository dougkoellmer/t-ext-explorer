"""FJReduction-SymPy 0.4.4: exact and formal Algorithms 1--13."""
import sympy as _sp
from .fields import (Field, Element, Parameters, FJError, prime_field,
                     finite_field, function_field, rational_field)
from .symbolic import (TwistSymbol, FJTwist, SymbolicField, SymbolicElement,
                       SymbolicParameters, FJSymbolicParameters,
                       SpecializeFrobenius, FJSymbolicConditions,
                       FJSpecializationReport)
from .compatible_basis import (preprocess, build_column_lists, order_column_lists,
                               compatible_ordering, dfs_selection,
                               compatible_basis, compatible_basis_report)
from .ore import OreMatrix, skew_multiply, skew_equal
from .reduction import ReductionData, construct_fj_basis
from .linear import normalize_matrix as _normalize_matrix, rank as _rank

__version__ = "0.4.4"
FJReductionVersion = __version__


def field_matrix(matrix, field):
    """Readable SymPy Matrix from a matrix of exact field elements."""
    if not matrix:
        return _sp.zeros(0, 0)
    return _sp.Matrix(len(matrix), len(matrix[0]),
                      lambda i, j: field.expr(matrix[i][j]))


def IsIndependent(vectors, column, field=None):
    field = field or rational_field()
    values = _normalize_matrix(list(vectors)+[list(column)], field)
    return _rank(values) == len(values)


def _data(*args):
    if len(args) == 1 and isinstance(args[0], ReductionData):
        return args[0]
    if len(args) == 3:
        return ReductionData(*args)
    raise FJError("InvalidArguments", "Use data, or phi, psi, parameters.")


# Mathematica-style names preserve recognizable entry points and argument order.
# Unlike Wolfram Failure objects, input errors raise FJError with a stable code.
CBPrimeField = prime_field
CBFiniteExtensionField = finite_field
CBRationalFunctionField = function_field
Preprocess = preprocess
BuildColumnLists = build_column_lists
OrderColumnLists = order_column_lists
CompatibleOrdering = compatible_ordering
DFSSelection = dfs_selection
CompatibleBasis = compatible_basis
CompatibleBasisReport = compatible_basis_report
FJParameters = Parameters
FJReductionData = ReductionData
ConstructFJBasis = construct_fj_basis
TMult = skew_multiply
SkewEqualQ = skew_equal


def CBFiniteField(p, generator, modulus):
    """Mathematica-compatible ordering of the finite-field constructor."""
    return finite_field(p, modulus, generator)


def CBNormalize(field, value):
    return field(value)


def CBInverse(field, value):
    return field(value).inverse()


def SkewNormalize(matrix, parameters):
    return OreMatrix.from_matrix(matrix, parameters)


def FJFrobenius(value, k, parameters):
    return parameters.frobenius(value, k)


def FJInnerBiderivation(u, data):
    return data.inner(u)


def FJReduce(v, data):
    return data.reduce(v)


def FJReduceReport(v, data):
    return data.reduce_report(v)


def FJExtension(*args):
    return _data(*args).extension()


# ===== BEGIN: computation of the true Ext_0 and exact sequence =====
# Public functions; implementation is in fj_sympy/sequence.py.

def FJReduceExtTau0(v, data):
    return data.reduce_zero(v)


def FJReduceExtTau0Report(v, data):
    return data.reduce_zero_report(v)


def FJExtTau0(*args):
    return _data(*args).ext_tau0()


def FJTransition(*args):
    return _data(*args).transition()


def GaPart(data):
    return data.ga_part()


def GaPartReport(data):
    return data.constant_data()


# ===== BEGIN: assembly and verification of the exact sequence (Algorithm 13) =====

def FJExactSequence(*args):
    data = _data(*args)
    report = data.exact_sequence()
    if isinstance(data.parameters, SymbolicParameters):
        report = dict(report, FrobeniusMode="Symbolic",
                      SpecializationConditions=FJSymbolicConditions(data.parameters))
    return report


def FJReport(*args):
    return FJExactSequence(*args)


# ===== END: assembly and verification of the exact sequence (Algorithm 13) =====
# ===== END: computation of the true Ext_0 and exact sequence =====


def FJGeneralReport(*args):
    data = _data(*args)
    report = data.general_report()
    if isinstance(data.parameters, SymbolicParameters):
        report = dict(report, FrobeniusMode="Symbolic",
                      SpecializationConditions=FJSymbolicConditions(data.parameters))
    return report


def FJGeneralReduceExtTau0Report(v, data):
    return data.general_reduce_zero_report(v)


def FJGeneralReduceExtTau0(v, data):
    return data.general_reduce_zero_report(v)["ReducedMatrix"]


def FJEvaluate(a, x, parameters=None):
    if parameters is not None:
        a = OreMatrix.from_matrix(a, parameters)
    return a.evaluate(x)


# ===== BEGIN: true Ext_0 construction (Algorithms 10-12) =====
# Public entry points for Algorithms 10-12.

def FJConstantData(data):
    return data.constant_data()


def FJSequenceCoordinates(data):
    return data.sequence_coordinates()


def FJSequenceActions(pi, data):
    return data.sequence_actions(pi)


# ===== END: true Ext_0 construction (Algorithms 10-12) =====

# Partial basis entry points below.
def FJPartialExactSequence(*args):
    return _data(*args).partial_exact_sequence()


def FJPartialExtTau0(*args):
    return _data(*args).partial_ext_tau0()


def FJPartialTransition(*args):
    return _data(*args).partial_transition()


def FJPartialReduceExtTau0Report(v, data):
    return data.partial_reduce_zero_report(v)


def FJPartialReduceExtTau0(v, data):
    return data.partial_reduce_zero(v)


def FJPartialGaPart(data):
    return data.partial_ga_part()


__all__ = [
    "TwistSymbol", "FJTwist", "SymbolicField", "SymbolicElement",
    "SymbolicParameters", "FJSymbolicParameters", "SpecializeFrobenius",
    "FJSymbolicConditions", "FJSpecializationReport",
    "Field", "Element", "Parameters", "FJError", "OreMatrix", "ReductionData",
    "prime_field", "finite_field", "function_field", "rational_field", "field_matrix",
    "preprocess", "build_column_lists", "order_column_lists", "compatible_ordering",
    "dfs_selection", "compatible_basis", "compatible_basis_report",
    "skew_multiply", "skew_equal", "construct_fj_basis", "FJReductionVersion",
    "CBPrimeField", "CBFiniteField", "CBFiniteExtensionField", "CBRationalFunctionField",
    "CBNormalize", "CBInverse",
    "Preprocess", "BuildColumnLists", "OrderColumnLists", "CompatibleOrdering",
    "DFSSelection", "CompatibleBasis", "CompatibleBasisReport", "IsIndependent",
    "FJParameters", "FJReductionData", "ConstructFJBasis", "TMult", "SkewEqualQ",
    "SkewNormalize", "FJFrobenius", "FJInnerBiderivation", "FJReduce",
    "FJReduceReport", "FJExtension", "FJReduceExtTau0", "FJReduceExtTau0Report",
    "FJExtTau0", "FJTransition", "GaPart", "GaPartReport", "FJExactSequence",
    "FJReport", "FJGeneralReport", "FJGeneralReduceExtTau0Report",
    "FJGeneralReduceExtTau0", "FJEvaluate", "FJConstantData", "FJSequenceCoordinates",
    "FJSequenceActions", "FJPartialExactSequence", "FJPartialExtTau0",
    "FJPartialTransition", "FJPartialReduceExtTau0Report", "FJPartialReduceExtTau0", "FJPartialGaPart"]

# from .example_reports import (FJExampleCompute,FJExampleReport,FJPartialBasisComparison,
#                               FJPartialBasisReport,FJJ2Example,ReportView)
# __all__ += ['FJExampleCompute','FJExampleReport','FJPartialBasisComparison',
#             'FJPartialBasisReport','FJJ2Example','ReportView']
