# FJReduction-SymPy 0.4.5

Automatic finite fields: choose a prime p and a positive integer n; q=p**n.

Open **MyExamples-Symbolic.ipynb** (section 1a) or **GettingStarted.ipynb** / **GettingStarted-Concrete.ipynb** (section 1). Change p and n, then **Run All**. MyExamples defaults to concrete mode; choose mode="symbolic" to retain formal Frobenius.

The package chooses a reproducible irreducible polynomial automatically. For n>1, alpha is the finite-field generator and alpha**q=alpha under field arithmetic. For n=1, use coefficients in the prime field; there is no algebraic alpha generator. Additional symbols in your matrices are transcendental coefficient variables, with q-power Frobenius. Thus symbolic a and b give F_q(theta,a,b), while coefficients expressed in alpha give F_q(theta).

Programmatic use:

```python
import sympy as sp
import fj_sympy as fj
theta, tau, alpha = sp.symbols("theta tau alpha")
p, n = 3, 2
parameters = fj.FJParametersPN(p, n, theta=theta, tau=tau, alpha=alpha)
print(parameters.field)
print(parameters.field.extension_modulus)
```

Use TMult for Ore multiplication and field coercion or FJFrobenius to reduce finite-field expressions. Plain SymPy powers do not apply the modulus automatically. Restart the kernel after replacing an older package. Very large p or n can make field construction and exact Frobenius calculations expensive.

GettingStarted sections 9 and 10 deliberately demonstrate fixed fields independently. Your own example and the final J2 example use your selected field. MyExamples section 5b also uses the selected field in concrete mode.

The full Ext, corrected Ext0, exact-sequence and partial-basis routines retain their existing behavior. BEGIN/END comments identify Ext0 calculations. MyExamples clears old results when field settings or hypotheses are invalid and gives an error. Existing validation files are historical evidence from earlier releases; current automatic-field checks are in Validation/automatic-fields-0.4.5.json and Validation/test-results.json.

## True Ext0 and exact sequence

D represents H -> H NPhi - NPsi H using column vectorization.
The quotient dimension is h=ed-rank D and dim Ext0=M-h. Its normal-coordinate
basis consists of image-D constants followed by the original positive monomials.
K is the constant inclusion, R recovers coordinates on its image, and P is the
quotient projection. The quotient action Q need not be scalar.

The report lists input matrices, numbered representatives, matrix sizes and
roles, the general quotient action and the identity behind every True.
Raw data stay in result["Sequence"]. R and AdditiveSection are additive coordinate
maps; they do not assert a splitting as t-modules. Use TMult for Ore products.

**BEGIN/END comments explicitly mark the true Ext_0 and exact-sequence
calculations in all example notebooks, their Python companions, and fj_sympy/sequence.py (Algorithms 10-13), fj_sympy/reduction.py and fj_sympy/__init__.py.**
The report adapter marks the actual FJExactSequence call as well.

## Partial basis and the J2 counterexample

Partial basis uses shifted monomial ranges determined by epsilon. When its
conditions hold, W=R T_partial and WInverse=L_partial K are explicitly computed.
The report checks both inverse identities, T_partial=K W and the action
intertwining identity using Ore multiplication. This proves the same subgroup
in different coordinates; it is not a different Ext0.

The revised counterexample is

```
Phi_t = [[theta+tau**2, 1], [0, theta+tau**2]]
Psi_t = [[theta, 1+tau], [0, theta]]
```

Both constant nilpotents are J2. The individual direction count is 1, but h=2;
Ext dimension is 8 and true Ext0 dimension is 6. The candidate count 7 is
invalid. Q=[[theta,0],[1,theta]] is nonscalar. For every c,

[c E11 tau^2] + [c E22 tau^2] - [c^q E12 tau] = 0.

Partial basis returns an explained applicability rejection for this pair.
The general rank-based calculation succeeds. The J2 cell leaves your input
matrices unchanged. The article subsection gives the full/restricted actions,
constant maps and exact sequence, including characteristic two.

## Reusable example/report API

```python
import fj_sympy as fj
import sympy as sp
theta,tau = sp.symbols("theta tau")
from ExampleReports import (FJExampleCompute, FJExampleReport,
    FJPartialBasisComparison, FJPartialBasisReport, FJJ2Example)
from IPython.display import display
p = fj.FJSymbolicParameters(theta=theta,tau=tau)
result = FJExampleCompute([[theta+tau**2]], [[theta+tau]], p)
display(FJExampleReport(result))
comparison = FJPartialBasisComparison(result["Data"], result["Sequence"])
display(FJPartialBasisReport(comparison))
j2 = FJJ2Example(theta,tau)
display(j2["Output"])
```

Use concrete FJParameters instead when fixing the field and q.
FJExampleCompute retains mathematical exceptions with their native error codes.
Generic symbolic specialization must preserve denominators, pivot conditions
and the compatible basis; reported conditions have their stated scope.

## Verification

Run `python run_tests.py` from the extracted folder. The release contains **120 tests**,
including seven new report/example regressions. For notebook checks, run
`python Validation/validate_notebooks_044.py`.
It executes every code cell in fresh processes, captures rich output, and tests
four changed-input scenarios. Validation reports record the actual engine,
source hashes and scope. Jupyter browser UI interactions are not claimed tested.

The SymPy checks execute this Python implementation, including the verified full and restricted J2 action matrices, formal Frobenius dependence and zero-constant witnesses.

## Source and article

ALGORITHMS.md maps the mathematical entry points. The original concrete examples
and partial basis APIs remain available. PartialReduction-Example.tex is the
English J2 subsection, identical to the revised Mathematica article fragment.
See existing copyright/attribution files for author credits and terms.
