# # FJReduction with Python and SymPy
# **Version 0.4.5 — Algorithms 1–13, with explanations and step-by-step comments.**
# 
# This notebook computes compatible jump bases, normal representatives,
# extension actions and the article's exact sequence. It uses exact arithmetic
# and runs in an ordinary Python kernel, without SageMath or Mathematica.
# 
# Extract the complete ZIP. In its folder, run
# **python -m pip install -r requirements.txt**, then
# **python -m notebook GettingStarted.ipynb**. Select a Python 3 kernel and run
# the cells in order.
# 
# The input presentation of $\Phi$ stays fixed throughout. The permutation
# $\sigma$ gives an inspection/elimination order, without changing $\Phi$.
# 
# **For a fixed characteristic and numerical $q$, start with this notebook.**
# For generic formulas with formal Frobenius, use `MyExamples-Symbolic.ipynb`.
# 
# **For a fixed characteristic and numerical $q$, start with this notebook.**
# For generic formulas with formal Frobenius, use `MyExamples-Symbolic.ipynb`.

# Step 1: locate the extracted package beside this notebook.
from pathlib import Path
import sys

project_dir = Path.cwd()
if not (project_dir / "fj_sympy").is_dir():
    candidate = project_dir / "FJReduction-SymPy"
    if candidate.is_dir():
        project_dir = candidate
    else:
        raise RuntimeError("Start Jupyter in the extracted FJReduction-SymPy folder.")
sys.path.insert(0, str(project_dir))

# Step 2: load SymPy, the package, and ordinary Jupyter display functions.
import sympy as sp
from IPython.display import display, Math
from fj_sympy import *

sp.init_printing(use_latex="mathjax")
print("FJReduction-SymPy:", FJReductionVersion)
print("SymPy:", sp.__version__)
print("Project folder:", project_dir.name)

# ## 1. Choose p and n automatically
# 
# Edit only p and n in the next cell. The package constructs F_(p**n)(theta), chooses an irreducible polynomial automatically and sets q=p**n. For n>1, alpha denotes the finite-field generator. Run All after changing either setting.
# 
# Sections 9 and 10 are independent demonstrations with their own explicitly fixed fields. Section 11 and the final J2 example use your selected p and n.

# FIELD SETTINGS: edit these two numbers.
p = 3
n = 1

theta, tau, alpha = sp.symbols("θ τ α")
parameters = FJParametersPN(p, n, theta=theta, tau=tau, alpha=alpha)
L = parameters.field
q = parameters.q
print(f"Field settings: p={p}, n={n}, q={q}")
print("Coefficient field:", L)
if L.extension_modulus is not None:
    print("Automatic modulus:", L.extension_modulus.as_expr())
display(Math(fr"\tau c=c^{{{q}}}\tau"))

# ## 2. The running example
# 
# The highest coefficient matrix of $\Phi_t$ is singular. The second
# coefficient block supplies the missing pivot. The target has dimension one
# and degree one.

# Enter the ORIGINAL presentations. No conjugation or isomorphic replacement
# will be performed by the basis-search or reduction routines.
phi = sp.Matrix([[theta + tau**2, 0],
                 [tau**3, theta + tau**2 + tau**3]])
psi = sp.Matrix([[theta + tau]])
display(Math(r"\Phi_t="), phi)
display(Math(r"\Psi_t="), psi)

# Read the positive coefficient blocks A_1, A_2, A_3.
phi_ore = OreMatrix.from_matrix(phi, parameters)
blocks = [phi_ore.coefficient(j) for j in range(1, phi_ore.degree + 1)]
for j, block in enumerate(blocks, start=1):
    display(Math(fr"A_{j}="), field_matrix(block, L))

# ## 3. Algorithms 1–6: find a compatible jump basis
# 
# Algorithm 1 finds the maximal suffix of full rank and its rank-increment
# quotas. Algorithms 3–4 form and order candidate lists. Algorithm 5 searches
# with quota and independence pruning; Algorithm 2 checks each complete
# selection using a precedence graph. Algorithm 6 combines these steps.
# 
# All reported column and block tags are one-based, as in the article.

# Step 1 (Algorithm 1): obtain the retained suffix, flag and block quotas.
found, trim, flag, required = Preprocess(blocks, L)
print("Full rank:", found)
print("First retained block:", trim["FirstBlock"])
print("Flag:", flag)
print("Required:", required)

# Step 2 (Algorithms 3--4): list nonzero columns and sort by list length.
candidate_lists = BuildColumnLists(trim)
ordered_lists = OrderColumnLists(candidate_lists)
print("Candidate lists:", candidate_lists)
print("Ordered lists:", ordered_lists)

# Step 3 (Algorithms 5 and 2): search and certify a compatible ordering.
found, basis, sigma = DFSSelection(trim, ordered_lists, required)
print("Compatible basis:", [(i, j, tuple(L.expr(x) for x in c))
                            for i, j, c in basis])
print("Permutation sigma:", sigma)

# Step 4 (Algorithm 6): the master routine gives the same result.
basis_report = CompatibleBasisReport(blocks, L)
assert basis_report["Found"] and basis_report["Sigma"] == sigma
display(Math(r"C="), field_matrix(basis_report["BasisMatrix"], L))

# ## 4. Algorithm 7: prepare the reduction data
# 
# The selected degrees are $(j_1,j_2)=(2,3)$, so
# $N=\min(j_1,j_2)=2$ and $J=j_1+j_2=5$. The required gap is
# $\deg_\tau\Psi_t=1<N$.
# 
# Preparation checks the nilpotent constant parts, verifies actual successive
# rank jumps and constructs the basic inner directions $E_{i,s}C^{-1}$.
# This uses $C^{-1}$ for elimination without changing $\Phi$.

# Step 1: construct and certify the reduction data for the original inputs.
data = ReductionData(phi, psi, parameters)

# Step 2: inspect the basis and bounds.
print("Compatible order:", [s + 1 for s in data.sigma])
print("Jump degrees:", data.jumps)
print("N = minimum jump degree:", data.minimum_jump)
print("J = sum of jump degrees:", sum(data.jumps))
print("Full coordinate family E:", data.basis())
assert data.phi == OreMatrix.from_matrix(phi, parameters)
display(Math(r"C^{-1}="), field_matrix(data.cb_inverse, L))

# ## 5. Algorithm 8: reduction with an inner witness
# 
# We reduce a representative $V$. Each step chooses the largest offending
# degree, then follows $\sigma$ and the target-row order. It subtracts an
# inner biderivation and records its coefficient and shift.
# 
# The report verifies
# $V-V_{\mathrm{red}}=U\Phi_t-\Psi_t U$ exactly.

# Step 1: choose a representative exceeding both column bounds.
V = sp.Matrix([[tau**4 + theta*tau, tau**4 + 1]])
reduction = FJReduceReport(V, data)

# Step 2: display the resulting normal form and accumulated inner witness.
display(Math(r"V="), V)
display(Math(r"V_{\mathrm{red}}="), reduction["ReducedMatrix"])
display(Math(r"U="), reduction["InnerWitness"])

# Step 3: explain each cancellation in the trace.
for number, step in enumerate(reduction["Steps"], start=1):
    print(f"Step {number}: degree {step['Degree']}, target row {step['Row']}, "
          f"original column {step['Column']}, shift {step['Shift']}, "
          f"coefficient {L.expr(step['Coefficient'])}")

# Step 4: independently recheck the reported identity and idempotence.
V_ore = OreMatrix.from_matrix(V, parameters)
assert V_ore - reduction["ReducedMatrix"] == data.inner(reduction["InnerWitness"])
assert data.reduce(reduction["ReducedMatrix"]) == reduction["ReducedMatrix"]
print("Inner witness verified:", reduction["VerifiedDifference"])

# ## 6. Algorithm 9: the full extension action
# 
# For every coordinate generator, the algorithm applies $\Psi_t$ to a
# **formal** coefficient $c$ and reduces the result. A term $b c^{q^h}$
# becomes the entry $b\tau^h$ in $\Pi_t$. Using only $c=1$ would lose these
# Frobenius indices.

# Compute the action in the full normal coordinate order E.
Pi = FJExtension(data)
display(Math(r"\Pi_t="), Pi)
print("Full dimension:", Pi.rows)

# Check a numeric coordinate vector by direct representative reduction.
x = [L(theta), L(1), L(2), L(theta + 1), L(0)]
direct_coordinates = data.coordinates(data.reduce(data.psi @ data.representative(x)))
assert direct_coordinates == Pi.evaluate(x)
print("Direct reduction agrees with the additive action:", True)

# ## 7. Correct Ext0 and the exact sequence: new Algorithms 10–13
# 
# Algorithm 10 forms the column-vectorized map $D$ of $H\mapsto HN_\Phi-N_\Psi H$ and computes $h=ed-\operatorname{rank}D$. Independent image columns form $B$, and $[B\ Z]^{-1}$ has rows $U,P_0$.
# 
# Algorithm 11 inserts image constants and positive normal monomials into the existing full coordinates with $K$. It also returns $R,P$. Algorithm 12 computes $\Pi_{0,t}=R\Pi_tK$ with Ore products and the general quotient action $Q_t$. Algorithm 13 checks ranks, exactness, morphisms and Anderson nilpotence. The quotient need not be scalar.

# ===== BEGIN: computation of the true Ext_0 and exact sequence =====
# ===== BEGIN: true Ext_0 construction (Algorithms 10-12) =====
# Algorithms 10-11: constant data and coordinate maps.
constant_info = FJConstantData(data)
coordinates0 = FJSequenceCoordinates(data)
print("h = ed - rank D:", constant_info["H"])
display(Math("D="), field_matrix(constant_info["D"], L))
print("Correct Ext0 basis:", coordinates0["ExtTau0Basis"])
# Full reduction tests P x=0 and returns R x.
zero_reduction = FJReduceExtTau0Report([[tau**4, tau**4]], data)
display(Math(r"V_{0,\mathrm{normal}}="), zero_reduction["ReducedMatrix"])
display(Math(r"V_{0,\mathrm{zero}}="), zero_reduction["ZeroConstantRepresentative"])
assert zero_reduction["VerifiedDifference"] and zero_reduction["ZeroConstantWitnessVerified"]
# Algorithm 12: restricted action and constant inclusion.
Pi0 = FJExtTau0(data)
K = FJTransition(data)
display(Math(r"\Pi_{0,t}="), Pi0)
display(Math("K="), K)

# ===== END: true Ext_0 construction (Algorithms 10-12) =====
# Blok główny kończy się po FJExactSequence w następnej komórce.


# Step 4: assemble and verify the full article exact sequence.
# ===== BEGIN: exact sequence assembly and checks (Algorithm 13) =====
sequence = FJExactSequence(data)
# ===== END: exact sequence assembly and checks (Algorithm 13) =====
# ===== END: computation of the true Ext_0 and exact sequence =====
P = sequence["P"]
display(Math(r"P="), P)
print("Dimensions (M, M0, h):",
      (sequence["FullDimension"], sequence["ExtTau0Dimension"], sequence["H"]))

# Step 5: inspect all six morphism/additive-inverse checks.
for key in ("InclusionVerified", "QuotientVerified", "CompositionZeroVerified",
            "InjectivityVerified", "SurjectivityVerified", "KernelImageVerified",
            "AndersonExt0Verified", "AndersonQuotientVerified"):
    print(key + ":", sequence[key])
assert sequence["ExactSequenceVerified"]
display(Math(r"0\longrightarrow \mathbb G_a^3"
             r"\xrightarrow{K}\mathbb G_a^5"
             r"\xrightarrow{P}\mathbb G_a^2\longrightarrow0"))
display(Math("R="), sequence["R"])
display(Math("Q_t="), sequence["QuotientMatrix"])


# ## 8. The fourteen-dimensional example
# 
# The dimensions remain $(M,M_0,h)=(14,12,2)$. The new inclusion $K$ is constant. Its image is the same Ext0 subgroup as the partial basis shifted-coordinate construction, but its basis and action coordinates differ. The old matrices are preserved only in the final “Partial reduction” example.

# Step 1: enter the original source and target.
phi_large = sp.Matrix([[theta + tau**4, tau],
                       [0, theta + tau**3]])
psi_large = sp.Matrix([[theta + theta**2*tau, 0],
                       [1 + tau, theta + theta*tau]])
display(Math(r"\Phi_t="), phi_large)
display(Math(r"\Psi_t="), psi_large)

# Step 2: prepare data and compute all action and exact-sequence matrices.
large_data = ReductionData(phi_large, psi_large, parameters)
# ===== BEGIN: computation of the true Ext_0 and exact sequence =====
# Algorithms 10-13 for this example.
large_sequence = FJExactSequence(large_data)
# ===== END: computation of the true Ext_0 and exact sequence =====
print("Jump degrees:", large_sequence["JumpDegrees"])
print("Compatible order:", large_sequence["Sigma"])
print("Dimensions:", (large_sequence["FullDimension"],
                      large_sequence["ExtTau0Dimension"], large_sequence["H"]))
print("epsilon:", large_sequence["Epsilon"])
assert large_sequence["ExactSequenceVerified"]

# Ext, then correct Ext0 and its exact-sequence maps.
for label, key in ((r"\Pi_t", "Pi"), (r"\Pi_{0,t}", "Pi0"), ("K", "K"), ("R", "R"), ("P", "P"), ("Q_t", "QuotientMatrix")):
    display(Math(label + "="), large_sequence[key])
assert large_sequence["K"].degree == 0


# ## 9. Finite extensions and $\mathbb F_9(\theta)$
# 
# We construct $\mathbb F_9=\mathbb F_3[\alpha]/(\alpha^2+1)$ from its modulus.
# The ordinary SymPy expression GF(9) is not used: it represents the integers
# modulo nine, which are not a field.
# 
# With $q=9$, $\alpha^9=\alpha$ and $\tau\theta=\theta^9\tau$.

# Step 1: construct the finite field by its irreducible polynomial.
F9 = finite_field(3, alpha**2 + 1)
L9 = function_field(F9, theta)
parameters9 = Parameters(L9, q=9, theta=theta, tau=tau)
assert not F9(alpha**2 + 1)
assert parameters9.frobenius(alpha) == L9(alpha)

# Step 2: verify the q=9 Ore multiplication rule.
twisted = TMult([[tau]], [[theta]], parameters9)
display(Math(r"\tau\theta="), twisted)
assert SkewEqualQ(twisted, [[theta**9*tau]], parameters9)

# Step 3: compute a scalar extension with an algebraic coefficient.
data9 = ReductionData([[theta + alpha*tau + tau**2]], [[theta + tau]], parameters9)
sequence9 = FJExactSequence(data9)
display(Math(r"\Pi_t="), sequence9["Pi"])
print("Exact sequence over F9(theta):", sequence9["ExactSequenceVerified"])

# ## 10. Matching nonzero nilpotents
# 
# For $N_\Phi=N_\Psi=J_2$, individual zero directions give the former count 1, whereas the correct rank calculation gives $h=2$. The quotient has a nonzero nilpotent action. The current algorithm handles this directly.

parameters3 = Parameters(prime_field(3), q=3, theta=1, tau=tau)
general_data = ReductionData([[1+tau**2, 1], [0, 1+tau**2]], [[1, 1], [0, 1]], parameters3)
# ===== BEGIN: computation of the true Ext_0 and exact sequence =====
# Algorithms 10-13 for this example.
general = FJExactSequence(general_data)
# ===== END: computation of the true Ext_0 and exact sequence =====
print("Former individual count (diagnostic only):", general["PseudocodeGaCount"])
print("Correct h:", general["H"])
print("Dimensions:", (general["FullDimension"], general["ExtTau0Dimension"], general["H"]))
display(Math("Q_t="), general["QuotientMatrix"])
assert general["H"] == 2 and general["ExactSequenceVerified"]


# ## 11. Enter your own example
# 
# Edit Phi_t and Psi_t below. This cell rebuilds the field using p and n from section 1. Use alpha-based coefficients for elements of the selected finite field when n>1. Additional coefficient symbols are treated as transcendental variables. The source must admit a compatible jump basis and the target degree must be less than every selected jump degree. Keep the original presentations, with constant parts theta I plus nilpotents.

# Step 2: replace these matrices with your own original Phi_t and Psi_t.
my_phi = [[theta + tau**3, tau, 0],
          [0, theta + tau**2, 0],
          [0, 0, theta + tau**2]]
my_psi = [[theta + tau]]

# Rebuild the field using the selected p and n and all input coefficient variables.
my_extra_symbols = sorted(
    (sp.Matrix(my_phi).free_symbols | sp.Matrix(my_psi).free_symbols)
    - {theta, tau, alpha}, key=str)
my_parameters = FJParametersPN(p, n, theta=theta, tau=tau, alpha=alpha,
                               coefficient_symbols=my_extra_symbols)
my_field = my_parameters.field
print(f"Own example: p={p}, n={n}, q={my_parameters.q}")
print("Coefficient field:", my_field)

# Step 3: compute and inspect the report.
my_data = ReductionData(my_phi, my_psi, my_parameters)
# ===== BEGIN: computation of the true Ext_0 and exact sequence =====
# Algorithms 10-13 for this example.
my_report = FJReport(my_data)
# ===== END: computation of the true Ext_0 and exact sequence =====
print("Jump degrees:", my_report["JumpDegrees"])
print("Dimensions:", (my_report["FullDimension"], my_report["ExtTau0Dimension"],
                      my_report["H"]))
display(Math(r"\Pi_t="), my_report["Pi"])
display(Math(r"\Pi_{0,t}="), my_report["Pi0"])
assert my_report["ExactSequenceVerified"]

# ## 12. Verify the package
# 
# This runs the regression suite, including independent basis-search oracles
# and entry-by-entry article matrix comparisons. It writes a JSON report.
# 
# Use @ on OreMatrix objects for multiplication. A.display() returns a display
# view; ordinary multiplication of those displayed SymPy matrices does not
# apply the Frobenius rule.

# Run the tests with THIS kernel's Python interpreter.
import subprocess
completed = subprocess.run([sys.executable, str(project_dir / "run_tests.py")],
                           cwd=project_dir, text=True, capture_output=True)
print(completed.stdout)
print(completed.stderr)
assert completed.returncode == 0, "See the test output above."

# Read the recorded counts.
import json
test_summary = json.loads((project_dir / "Validation" / "test-results.json").read_text())
print("Tests:", test_summary["tests"])
print("Failures:", test_summary["failures"])
print("Errors:", test_summary["errors"])
assert test_summary["all_passed"]

# ## Full descriptive report
# Numbered representatives, matrix sizes and roles, and exact checked identities.

from IPython.display import display
from ExampleReports import *
# ===== BEGIN: computation of the true Ext_0 and exact sequence =====
walkthrough_result = FJExampleCompute(data.phi, data.psi, data.parameters)
# ===== END: computation of the true Ext_0 and exact sequence =====
display(FJExampleReport(walkthrough_result))

# ## Partial basis
# An explicit Ore coordinate isomorphism is checked when the candidate is valid. The general rank construction applies also when the candidate fails.

walkthrough_comparison = FJPartialBasisComparison(walkthrough_result["Data"],walkthrough_result["Sequence"])
display(FJPartialBasisReport(walkthrough_comparison))

# ## J2 counterexample
# The independent pair uses Psi upper-right 1+tau. Computed dimensions are 8/6/2; the individual count is 1. This final example does not overwrite the tutorial input.

# ===== BEGIN: computation of the true Ext_0 and exact sequence =====
j2_theta,j2_tau = sp.symbols("theta tau")
j2_parameters = FJParametersPN(p, n, theta=j2_theta, tau=j2_tau)
j2Counterexample = FJJ2Example(j2_theta, j2_tau, parameters=j2_parameters)
# ===== END: computation of the true Ext_0 and exact sequence =====
display(j2Counterexample["Output"])
summary = j2Counterexample["Summary"]
assert tuple(summary[k] for k in ("FullDimension", "TrueExt0Dimension", "H", "IndividualZeroCount", "InvalidCountPrediction")) == (8, 6, 2, 1, 7)
assert summary["PartialBasisStatus"] == "Unavailable"
assert j2Counterexample["Result"]["Sequence"]["ExactSequenceVerified"]

