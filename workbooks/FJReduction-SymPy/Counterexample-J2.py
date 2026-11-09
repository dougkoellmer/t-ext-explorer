# # J2 counterexample
# Computed example with Psi upper-right 1+tau. The subsection is in PartialReduction-Example.tex.

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


# ===== BEGIN: computation of the true Ext_0 and exact sequence =====
j2_theta,j2_tau = sp.symbols("theta tau")
j2Counterexample = FJJ2Example(j2_theta, j2_tau)
# ===== END: computation of the true Ext_0 and exact sequence =====
display(j2Counterexample["Output"])
summary = j2Counterexample["Summary"]
assert tuple(summary[k] for k in ("FullDimension", "TrueExt0Dimension", "H", "IndividualZeroCount", "InvalidCountPrediction")) == (8, 6, 2, 1, 7)
assert summary["PartialBasisStatus"] == "Unavailable"
assert j2Counterexample["Result"]["Sequence"]["ExactSequenceVerified"]

