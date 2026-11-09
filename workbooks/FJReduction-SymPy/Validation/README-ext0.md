# Validation for 0.4.5

The full Python/SymPy implementation was executed. All 109 tests passed,
including nine new tests in tests/test_exact_sequence.py and the existing
100 regressions. No tests were skipped.

The new tests cover:
- Independent Kronecker formula for D, h=2 versus former count 1 for J2.
- Nonscalar quotients and exact sequences over prime, extension and function
  fields, including positive-degree targets and nonidentity jump bases.
- Empty positive blocks and zero-dimensional Ext0.
- Removable nonzero constants, both inner witnesses, and quotient obstructions.
- Formal coefficient families and required Frobenius twists on K.
- Symbolic/concrete nonscalar comparisons in characteristics 2,3,5.
- Preservation of the old partial matrices and an invertible Ore coordinate
  change to the new Ext0 model.
- New Algorithm 10--12 entry points and invalid action shapes.

Validation/test-results.json records actual counts, timing and source hashes.
Validation/test-output.txt contains the complete successful test output.
The partial basis matrix-entry regressions remain attached to partial_* methods;
the same field cases independently verify the current exact sequence.

verify_notebooks.py executed 50 code cells across five fresh-IPython scenarios:
concrete walkthrough, default symbolic walkthrough, separate partial example,
new coefficient names introduced by changing only Phi/Psi, and a user pair
with matching J2 nilpotents and a nonscalar quotient. All passed. Executed
outputs are stored in the delivered notebooks. GettingStarted.ipynb is the
compatibility copy of the verified concrete walkthrough. The Python companions
contain the same inputs. Notebook JSON validation and Python syntax checks passed.

This is real Python code-cell execution with captured rich display, not a
browser-based Jupyter UI inspection. No Mathematica or SageMath execution is
claimed. full-ext-preservation.json confirms that the full reduction, normal
basis and Ext source block is byte-identical to the supplied archive.
Partial basis reference-provenance.json refers to the old partial-coordinate
formulas and should be read together with the current source-hash reports.

Earlier comment-only update: begin/end markers were added to the notebook code,
Python companions and source modules. comment-markers.json certifies that all
Python and affected notebook code-cell ASTs are unchanged, and records current
source hashes. Executed notebook outputs are preserved. The preceding 109-test
report records the validated code before these comments; no test rerun was
necessary for this update.
