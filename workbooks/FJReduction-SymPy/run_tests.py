"""Run all tests and write a compact reproducible validation report."""
from pathlib import Path
import hashlib
import json
import platform
import sys
import time
import unittest
import sympy

ROOT = Path(__file__).resolve().parent


def main():
    sys.path.insert(0, str(ROOT))
    suite = unittest.defaultTestLoader.discover(str(ROOT/"tests"), pattern="test_*.py")
    start = time.perf_counter()
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted((ROOT/"fj_sympy").glob("*.py"))}
    report = {
        "implementation": "FJReduction-SymPy 0.4.5",
        "python": platform.python_version(), "sympy": sympy.__version__,
        "tests": result.testsRun, "failures": len(result.failures), "errors": len(result.errors),
        "skipped": len(result.skipped), "all_passed": result.wasSuccessful(),
        "elapsed_seconds": round(time.perf_counter()-start, 3),
        "independent_checks": [
            "512 binary matrix comparisons over Q and F2 against a rank-jump oracle",
            "36 random dimension-three basis existence comparisons",
            "Every partial basis partial-coordinate matrix entry across 11 field/q configurations; current exact sequence ranks, Ore identities and Anderson conditions",
            "Numeric coefficient-mode reduction versus formal additive action columns",
            "Six Ore polynomial morphism/inverse identities for each article exact sequence",
            "18 random reduction witness, idempotence and additive action cases",
            "Formal Frobenius laws, rational normalization and Ore associativity",
            "Independent symbolic/concrete comparisons for q=2,3,5 and q=9",
            "Recorded denominator/pivot conditions and rejected specializations"],
        "source_sha256": hashes,
        "scope": "Executes this SymPy implementation. Reference outputs come from the supplied Mathematica 0.2.0 tests and article formulas; Mathematica/Sage are not rerun here."}
    (ROOT/"Validation").mkdir(exist_ok=True)
    (ROOT/"Validation/test-results.json").write_text(json.dumps(report, indent=2)+"\n")
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
