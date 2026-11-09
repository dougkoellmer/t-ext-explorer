# # My examples: automatic finite fields, Ext and Ext0
# 
# Version 0.4.5. Edit p and n in section 1a and Phi/Psi in section 2. Run All after changing the field. The symbolic filename is retained; choose mode="symbolic" for formal Frobenius.

# ## 1. Load the package
# Keep the notebook and source files in the extracted package folder. Restart the kernel after installing a new version.

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
theta, tau, alpha, a, b, c = sp.symbols("theta tau alpha a b c")
def parse_module(rows):
    aliases = {"theta": theta, "tau": tau, "θ": theta, "τ": tau, "FJTwist": fj.FJTwist, "alpha": alpha, "α": alpha}
    return sp.Matrix([[sp.sympify(entry, locals=aliases) for entry in row] for row in rows])


# ## 1a. Choose your field automatically
# Set any prime p and positive integer n. Concrete mode constructs F_(p**n)(theta) automatically; no modulus input is needed. For n>1 use alpha for the finite-field generator. Other coefficient symbols such as a and b are transcendental variables, so the field is enlarged to F_(p**n)(theta,a,b,...). Use numerical or alpha-based coefficients to stay in F_(p**n)(theta). After changing the settings, Run All. Finite-field construction is performed inside section 3 before every calculation.

# FIELD SETTINGS: edit only p and n to choose F_(p**n).
p = 3
n = 1
mode = "concrete"  # "concrete" uses q=p**n; "symbolic" keeps formal Frobenius.


# ## 2. Your input matrices
# Use coefficients on the left of tau, `**` for powers and `*` for multiplication. The source must have an FJ basis, nilpotent constant parts and deg(Psi) smaller than its minimum jump degree. General coefficient symbols are not Frobenius-fixed. Strings such as `"theta + u*tau**3"` introduce new symbols automatically.

# EDIT ONLY THESE TWO MATRICES for your own example.
Phi = [[theta + tau**4, tau], [0, theta + tau**3]]
Psi = [[theta + a*tau, 0], [1 + tau, theta + b*tau]]


# ## 3. Compute Ext, the true Ext0 and exact sequence
# BEGIN/END comments mark the calculation. Each input change creates fresh parameters. Invalid input clears prior results and reports the actual condition.

# Clear old results before computing the current inputs.
example_result = data = sequence = None
calculation_error = None
try:
    Phi_input, Psi_input = parse_module(Phi), parse_module(Psi)
    if mode == "concrete":
        extra_symbols = sorted(
            (Phi_input.free_symbols | Psi_input.free_symbols) - {theta, tau, alpha},
            key=str)
        parameters = fj.FJParametersPN(p, n, theta=theta, tau=tau,
                                      alpha=alpha, coefficient_symbols=extra_symbols)
        print(f"Field settings: p={p}, n={n}, q={parameters.q}")
        print("Coefficient field:", parameters.field)
        if parameters.field.extension_modulus is not None:
            print("Automatic modulus:", parameters.field.extension_modulus.as_expr())
        if extra_symbols:
            print("Additional coefficient variables:", extra_symbols,
                  "(transcendental; Frobenius acts by q-th powers)")
    elif mode == "symbolic":
        parameters = fj.FJSymbolicParameters(theta=theta, tau=tau)
        print("Formal Frobenius mode: p and n are not specialized.")
    else:
        raise ValueError('mode must be "concrete" or "symbolic"')
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


# ## 4. Descriptive output
# Ext is displayed first, followed by the true Ext0 and exact sequence. The report numbers basis representatives and explains matrix shapes, roles and every verification identity. Raw objects remain in `example_result["Sequence"]`.

if example_result is not None:
    display(FJExampleReport(example_result))
else:
    print("No valid result:", calculation_error)

# ## 5a. Partial basis for your input
# When applicable, W and its additive-polynomial inverse verify the same subgroup in different coordinates. Otherwise, the precise failed conditions are explained; the general Ext0 result remains available.

partial_comparison = None
if example_result is None:
    print("First compute a valid input in section 3.")
else:
    partial_comparison = FJPartialBasisComparison(data, sequence)
    display(FJPartialBasisReport(partial_comparison))


# ## 5b. J2 counterexample with upper-right 1+tau in Psi
# Here Phi=theta I2+J2+tau^2 I2 and Psi=theta I2+J2+J2 tau. The direction count is 1, but h=2 and dim Ext0=6. The candidate count 7 is invalid. The quotient action is nonscalar. For every c, the candidate parameters (0,c,-c^q,0,0,0,c) give the zero class.

# ===== BEGIN: computation of the true Ext_0 and exact sequence =====
j2_theta,j2_tau = sp.symbols("theta tau")
j2_parameters = (fj.FJParametersPN(p, n, theta=j2_theta, tau=j2_tau)
                 if mode == "concrete" else None)
j2Counterexample = FJJ2Example(j2_theta, j2_tau, parameters=j2_parameters)
# ===== END: computation of the true Ext_0 and exact sequence =====
display(j2Counterexample["Output"])
summary = j2Counterexample["Summary"]
assert tuple(summary[k] for k in ("FullDimension", "TrueExt0Dimension", "H", "IndividualZeroCount", "InvalidCountPrediction")) == (8, 6, 2, 1, 7)
assert summary["PartialBasisStatus"] == "Unavailable"
assert j2Counterexample["Result"]["Sequence"]["ExactSequenceVerified"]

