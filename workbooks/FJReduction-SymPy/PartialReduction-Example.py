# # Partial basis: coordinate comparison and J2
# Edit the input pair below for your own examples. For the supplied 14-dimensional pair, both coordinate systems give the same Ext0. The separate J2 pair demonstrates failure of the individual count.

from pathlib import Path
import sys
project_dir = Path.cwd()
if not (project_dir / "fj_sympy").is_dir():
    raise RuntimeError("Start Jupyter in the extracted FJReduction-SymPy folder.")
sys.path.insert(0, str(project_dir))
import sympy as sp
import fj_sympy as fj
from IPython.display import display
from ExampleReports import *
print("FJReduction version:", fj.FJReductionVersion)
assert fj.FJReductionVersion == "0.4.5"
theta, tau, a, b, c = sp.symbols("theta tau a b c")
def parse_module(rows):
    aliases = {"theta": theta, "tau": tau, "θ": theta, "τ": tau, "FJTwist": fj.FJTwist}
    return sp.Matrix([[sp.sympify(entry, locals=aliases) for entry in row] for row in rows])


# EDIT ONLY THESE TWO MATRICES for your own example.
Phi = [[theta + tau**4, tau], [0, theta + tau**3]]
Psi = [[theta + a*tau, 0], [1 + tau, theta + b*tau]]


# Clear old results before computing the current inputs.
example_result = data = sequence = None
calculation_error = None
try:
    Phi_input, Psi_input = parse_module(Phi), parse_module(Psi)
    parameters = fj.FJSymbolicParameters(theta=theta, tau=tau)
    # ===== BEGIN: computation of the true Ext_0 and exact sequence =====
    example_result = FJExampleCompute(Phi_input, Psi_input, parameters)
    data, sequence = example_result["Data"], example_result["Sequence"]
    # ===== END: computation of the true Ext_0 and exact sequence =====
except (fj.FJError, ValueError, TypeError, ZeroDivisionError) as error:
    calculation_error = error
    print("The selected hypotheses do not permit this input:", error)
else:
    print("Dimensions (Ext, Ext0, quotient):", tuple(sequence[k] for k in ("FullDimension", "ExtTau0Dimension", "H")))
    assert sequence["ExactSequenceVerified"]


if example_result is not None:
    display(FJExampleReport(example_result))

partial_comparison = None
if example_result is None:
    print("First compute a valid input in section 3.")
else:
    partial_comparison = FJPartialBasisComparison(data, sequence)
    display(FJPartialBasisReport(partial_comparison))


# ===== BEGIN: computation of the true Ext_0 and exact sequence =====
j2_theta,j2_tau = sp.symbols("theta tau")
j2Counterexample = FJJ2Example(j2_theta, j2_tau)
# ===== END: computation of the true Ext_0 and exact sequence =====
display(j2Counterexample["Output"])
summary = j2Counterexample["Summary"]
assert tuple(summary[k] for k in ("FullDimension", "TrueExt0Dimension", "H", "IndividualZeroCount", "InvalidCountPrediction")) == (8, 6, 2, 1, 7)
assert summary["PartialBasisStatus"] == "Unavailable"
assert j2Counterexample["Result"]["Sequence"]["ExactSequenceVerified"]

