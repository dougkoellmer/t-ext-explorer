"""Regression tests and independent mathematical checks.

The existence oracle enumerates permutations and actual rank jumps, using
SymPy DomainMatrix/Matrix for rank rather than the production graph or DFS.
The large example is compared entry by entry with the article's LaTeX.
"""
import copy
from itertools import permutations, product
from pathlib import Path
import random
import re
import unittest
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from fj_sympy import *
from fj_sympy.linear import identity, zero, rank, inverse, multiply

alpha, theta, tau = sp.symbols("α θ τ")
ROOT = Path(__file__).resolve().parents[1]


def oracle_rank(columns, d, p):
    if not columns:
        return 0
    matrix = [[v[i] for v in columns] for i in range(d)]
    if p == 0:
        return sp.Matrix(matrix).rank()
    return DomainMatrix.from_list_sympy(d, len(columns), matrix).convert_to(sp.GF(p)).rank()


def oracle_jump_basis(blocks, p):
    d = len(blocks[0])
    for sigma in permutations(range(d)):
        columns, tags = [], []
        for j in range(len(blocks)-1, -1, -1):
            for s in sigma:
                column = [row[s] for row in blocks[j]]
                if oracle_rank(columns+[column], d, p) > len(columns):
                    columns.append(column)
                    tags.append((s, j+1))
        if len(columns) == d and len({s for s, _ in tags}) == d:
            return True
    return False


def expected_article_matrices(p):
    """Parse the four displayed matrices, independently of reduction."""
    source = (ROOT/"tests/fixtures/article-large-example.tex").read_text()
    th = p.field.expr(p.theta)
    a, b, delta = th**2, th, th**p.q-th
    values = {"0": 0, "1": 1, "-1": -1, "theta": th,
              "tau": tau, "-tau": -tau, "atau": a*tau, "-atau": -a*tau,
              "btau": b*tau, "-btau": -b*tau, "btau^2": b*tau**2,
              "alpha": a*a**p.q*tau**2, "beta": b*b**p.q*tau**2,
              "gamma": (b+a**p.q)*tau**2, "xi": -a*delta*tau,
              "omega": (a+b-delta)*tau}
    patterns = {
        "Pi": r'\\Pi_t=\s*\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}',
        "Pi0": r'\\Pi_\{0,t\}=\s*\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}',
        "T": r' T=\s*\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}',
        "P": r' P=\s*\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}'}
    result = {}
    for name, pattern in patterns.items():
        raw = re.search(pattern, source, re.S).group(1)
        rows = [[re.sub(r'\s|\\', '', cell) for cell in row.split('&')]
                for row in raw.split(r'\\') if row.strip()]
        result[name] = OreMatrix.from_matrix([[values[cell] for cell in row] for row in rows], p)
    return result


class FieldsTests(unittest.TestCase):
    def test_prime_fields_exact(self):
        for p in (2, 3, 5, 7):
            f = prime_field(p)
            self.assertFalse(f(p))
            for i in range(1, p):
                self.assertEqual(f(i)*f(i).inverse(), f.one)

    def test_F4_and_F9_all_units(self):
        for p, polynomial in ((2, alpha**2+alpha+1), (3, alpha**2+1)):
            f = finite_field(p, polynomial)
            self.assertFalse(f(polynomial))
            for coefficients in product(range(p), repeat=2):
                value = f(coefficients[0]+coefficients[1]*alpha)
                if value:
                    self.assertEqual(value*value.inverse(), f.one)
                    self.assertEqual(value**(p**2), value)

    def test_fraction_field_extension_inverses(self):
        for p, polynomial in ((2, alpha**2+alpha+1), (3, alpha**2+1)):
            f = function_field(finite_field(p, polynomial), theta)
            for expression in (theta+alpha, (theta+alpha)/(alpha+1),
                               (theta**2+alpha*theta+1)/(theta+1)):
                a = f(expression)
                self.assertEqual(a/a, f.one)
                self.assertEqual(a*a.inverse(), f.one)
                self.assertEqual(f(a.as_expr()), a)

    def test_fraction_unit_normalization(self):
        f = function_field(prime_field(3), theta)
        self.assertEqual(f(2)/f(2), f.one)
        self.assertEqual(f(2*theta)/f(2), f(theta))
        self.assertEqual(f((theta**2-1)/(theta-1)), f(theta+1))

    def test_multivariable_extension(self):
        x = sp.Symbol("x")
        f = function_field(finite_field(3, alpha**2+1), theta, x)
        a = f((alpha*theta+x)/(theta+x+1))
        self.assertEqual(a/a, f.one)

    def test_invalid_fields(self):
        for p in (1, 4, 9):
            with self.assertRaises(FJError):
                prime_field(p)
        with self.assertRaises(FJError):
            finite_field(3, alpha**2-1)

    def test_float_and_zero_denominator(self):
        f = prime_field(3)
        for value in (0.5, sp.Rational(1, 3)):
            with self.assertRaises(FJError):
                f(value)

    def test_parameters(self):
        for q in (1, 4, 9):
            with self.assertRaises(FJError):
                Parameters(prime_field(3), q, 1, tau)
        p = Parameters(finite_field(3, alpha**2+1), 9, alpha, tau)
        self.assertEqual(p.frobenius(alpha), p.field(alpha))
        with self.assertRaises(FJError):
            p.frobenius(1, -1)


class CompatibleBasisTests(unittest.TestCase):
    def test_reported_example(self):
        a = [[1, 0, 0, 0], [0, 1, 1, 1]]
        report = compatible_basis_report(a, prime_field(3))
        self.assertTrue(report["Found"])
        self.assertEqual(report["Sigma"], [2, 1])
        self.assertEqual(report["Basis"], [(1, 1, (1, 0)), (2, 2, (0, 1))])

    def test_theta_column_preserved(self):
        f = function_field(prime_field(3), theta)
        report = compatible_basis_report([[theta, 0, 0, 0], [0, 1, 1, 1]], f)
        self.assertEqual(report["Basis"][0][2][0], f(theta))
        self.assertEqual(report["Sigma"], [2, 1])

    def test_maximal_trim_original_indices(self):
        report = compatible_basis_report([[9, 4, 1, 0, 0, 0], [7, 2, 0, 1, 0, 0]])
        self.assertEqual(report["BlockChoices"], [2, 2])
        self.assertEqual(report["Required"], [0, 2, 0])
        self.assertEqual(report["Flag"], [(2, 2), (3, 0)])

    def test_block_and_horizontal_inputs(self):
        f = prime_field(3)
        blocks = [[[1, 0], [0, 1]], [[0, 0], [1, 1]]]
        self.assertEqual(compatible_basis(blocks, f),
                         compatible_basis([[1, 0, 0, 0], [0, 1, 1, 1]], f))

    def test_stable_sort_and_column_lists(self):
        _, trim, _, _ = preprocess([[1, 0, 0, 0], [0, 1, 1, 1]])
        self.assertEqual(build_column_lists(trim),
                         [[(1, 2, (0, 1)), (1, 1, (1, 0))],
                          [(2, 2, (0, 1)), (2, 1, (0, 1))]])
        self.assertEqual(order_column_lists([[1, 2], [3], [4, 5], [6]]),
                         [[3], [6], [1, 2], [4, 5]])

    def test_cycle_and_rank_deficiency(self):
        for a in ([[1, 1, 0, 0], [0, 1, 1, 1]],
                  [[1, 2, 3, 4], [0, 0, 0, 0]],
                  [[1, 0, 0, 0], [0, 0, 1, 0]]):
            self.assertFalse(compatible_basis(a)[0])

    def test_no_independent_transversal(self):
        a = [[1, 1, 0, 1, 1, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1]]
        self.assertFalse(compatible_basis(a)[0])

    def test_backtracking_leaf_regression(self):
        a = [[0, 0, 0, 0, 0, -1, 1, -1, 0],
             [1, 0, -1, 0, 0, 0, 1, 0, -1],
             [-1, -1, 0, 1, 0, 1, 0, 0, 0]]
        self.assertFalse(compatible_basis(a, prime_field(3))[0])

    def test_exhaustive_512_field_comparisons(self):
        # 256 binary 2x4 inputs, independently checked over Q and F_2.
        for entries in product((0, 1), repeat=8):
            a = [list(entries[:4]), list(entries[4:])]
            blocks = [[[row[j] for j in range(k, k+2)] for row in a] for k in (0, 2)]
            for p in (0, 2):
                f = rational_field() if not p else prime_field(p)
                self.assertEqual(compatible_basis(a, f)[0], oracle_jump_basis(blocks, p),
                                 msg=f"{a}, characteristic {p}")

    def test_random_three_dimensional_oracle(self):
        rng = random.Random(1463)
        for p in (2, 3, 5):
            for _ in range(12):
                blocks = [[[rng.randrange(p) for _ in range(3)] for _ in range(3)] for _ in range(3)]
                self.assertEqual(compatible_basis(blocks, prime_field(p))[0],
                                 oracle_jump_basis(blocks, p))

    def test_inputs_are_not_mutated(self):
        a = [[1, 0, 0, 0], [0, 1, 1, 1]]
        original = copy.deepcopy(a)
        compatible_basis(a)
        self.assertEqual(a, original)

    def test_input_errors(self):
        for a in ([], [[1, 0, 1], [0, 1, 0]], [[1], [0, 1]], [[0.5]]):
            with self.assertRaises(FJError):
                compatible_basis(a)

    def test_invalid_quota(self):
        _, trim, _, _ = preprocess([[1, 0], [0, 1]])
        with self.assertRaises(FJError):
            dfs_selection(trim, build_column_lists(trim), [1])


class OreTests(unittest.TestCase):
    def setUp(self):
        self.f = function_field(prime_field(3), theta)
        self.p = Parameters(self.f, 3, theta, tau)

    def test_noncommutativity(self):
        left = TMult([[tau]], [[theta]], self.p)
        right = TMult([[theta]], [[tau]], self.p)
        self.assertEqual(left.display(), sp.Matrix([[theta**3*tau]]))
        self.assertNotEqual(left, right)

    def test_associativity(self):
        a, b, c = (OreMatrix.from_matrix(m, self.p) for m in
                   ([[tau, 1], [2, tau**2]], [[theta+tau, 2], [tau, 1]],
                    [[tau, 0], [theta, 2+tau]]))
        self.assertEqual((a@b)@c, a@(b@c))

    def test_skew_normalization(self):
        self.assertTrue(SkewEqualQ([[3*tau**2+4*tau-1]], [[tau+2]], self.p))

    def test_formal_frobenius(self):
        f = self.f
        formal = OreMatrix.monomial([[(f.one,)]], 0, self.p, True)
        twist = OreMatrix.from_matrix([[tau]], self.p) @ formal
        self.assertEqual(twist.coefficient(1)[0][0], (f.zero, f.one))
        self.assertEqual(twist.evaluate_formal(theta),
                         OreMatrix.from_matrix([[theta**3*tau]], self.p))

    def test_additive_evaluation(self):
        m = OreMatrix.from_matrix([[tau, theta], [0, tau**2]], self.p)
        self.assertEqual(m.evaluate([theta, 1]), [self.f(theta**3+theta), self.f.one])

    def test_non_polynomial_rejected(self):
        with self.assertRaises(FJError):
            OreMatrix.from_matrix([[1/tau]], self.p)

    def test_coefficient_copy_and_immutability(self):
        m = OreMatrix.from_matrix([[theta+tau]], self.p)
        coeff = m.coefficient(0)
        coeff[0][0] = self.f.zero
        self.assertEqual(m.coefficient(0)[0][0], self.f(theta))


class ReductionTests(unittest.TestCase):
    def setUp(self):
        self.f = prime_field(3)
        self.p = Parameters(self.f, 3, 1, tau)
        self.data = ReductionData([[1+tau**2]], [[1+tau]], self.p)

    def test_scalar_reduction_and_trace(self):
        r = self.data.reduce_report([[tau**7+2*tau**3+2]])
        self.assertTrue(SkewEqualQ(r["ReducedMatrix"], [[2]], self.p))
        self.assertTrue(r["VerifiedDifference"])
        degrees = [s["Degree"] for s in r["Steps"]]
        self.assertEqual(degrees, sorted(degrees, reverse=True))

    def test_witness_idempotence(self):
        v = OreMatrix.from_matrix([[2*tau**6+tau**3]], self.p)
        r = self.data.reduce_report(v)
        self.assertEqual(v-r["ReducedMatrix"], self.data.inner(r["InnerWitness"]))
        self.assertEqual(self.data.reduce(r["ReducedMatrix"]), r["ReducedMatrix"])

    def test_scalar_matrices_and_inverse_maps(self):
        r = self.data.exact_sequence()
        self.assertEqual(r["Pi"], OreMatrix.from_matrix([[1, 0], [tau, 1+tau**2]], self.p))
        self.assertEqual(r["Pi0"], OreMatrix.from_matrix([[1+tau**2]], self.p))
        self.assertEqual(r["T"].display(), sp.Matrix([[0], [1]]))
        self.assertEqual(r["InclusionLeftInverse"].display(), sp.Matrix([[0, 1]]))
        self.assertEqual(r["AdditiveSection"].display(), sp.Matrix([[1], [0]]))

    def test_zero_reduction(self):
        r = self.data.reduce_zero_report([[tau**2]])
        self.assertEqual(r["ReducedMatrix"].display(), sp.Matrix([[tau]]))
        self.assertTrue(r["VerifiedDifference"])
        with self.assertRaises(FJError):
            self.data.reduce_zero([[1+tau**2]])

    def test_original_phi_and_input_preserved(self):
        phi, psi = [[1+tau**2, 0], [tau**3, 1+tau**2+tau**3]], [[1+tau]]
        originals = copy.deepcopy((phi, psi))
        data = ReductionData(phi, psi, self.p)
        data.exact_sequence()
        self.assertEqual((phi, psi), originals)
        self.assertEqual(data.phi, OreMatrix.from_matrix(originals[0], self.p))
        self.assertEqual(data.psi, OreMatrix.from_matrix(originals[1], self.p))

    def test_dimension_zero(self):
        r = ReductionData([[1+tau]], [[1]], self.p).exact_sequence()
        self.assertEqual((r["Pi0"].shape, r["T"].shape, r["P"].shape),
                         ((0, 0), (1, 0), (1, 1)))
        self.assertTrue(r["ExactSequenceVerified"])

    def test_degree_nilpotent_shape_errors(self):
        cases = [([[1+tau**2]], [[1+tau**2]]),
                 ([[2+tau**2]], [[1+tau]]),
                 ([[1+tau, 0, 0]], [[1]]),
                 ([[1+tau, 0], [tau**2, 1]], [[1]])]
        for phi, psi in cases:
            with self.assertRaises(FJError):
                ReductionData(phi, psi, self.p)
        with self.assertRaises(FJError):
            self.data.reduce([[tau**2, 0]])

    def test_dependent_constant_directions_general_branch(self):
        phi = [[1+tau**2, 1], [0, 1+tau**2]]
        psi = [[1, 1], [0, 1]]
        d = ReductionData(phi, psi, self.p)
        a = d.constant_analysis()
        self.assertEqual((a["H"], a["KernelDimension"], a["ConstantMapRank"]), (1, 2, 2))
        with self.assertRaises(FJError) as caught:
            d.partial_exact_sequence()
        self.assertEqual(caught.exception.code, "ArticleCoordinateConditionFailed")
        r = d.general_report()
        self.assertEqual((r["FullDimension"], r["ExtTau0Dimension"], r["QuotientDimension"]), (8, 6, 2))
        self.assertIsNone(r["GaDimension"])
        self.assertTrue(d.general_reduce_zero_report([[tau**2, 0], [0, tau**2]])["ReducedMatrix"].is_zero)

    def test_general_zero_positive_target(self):
        d = ReductionData([[1+tau**2, 1], [0, 1+tau**2]],
                          [[1+tau, 1], [0, 1+tau]], self.p)
        v = [[tau**4, tau**3], [2*tau**2, tau**5]]
        r = d.general_reduce_zero_report(v)
        self.assertTrue(r["VerifiedDifference"])
        self.assertEqual(d.reduce(r["ReducedMatrix"]), r["FullReducedRepresentative"])
        self.assertEqual(d.general_reduce_zero_report(r["ReducedMatrix"])["ReducedMatrix"],
                         r["ReducedMatrix"])


def _field_cases():
    f4, f9 = finite_field(2, alpha**2+alpha+1), finite_field(3, alpha**2+1)
    return [
        ("F2", prime_field(2), 2, 1), ("F3", prime_field(3), 3, 1),
        ("F5", prime_field(5), 5, 1), ("F4_q2", f4, 2, alpha),
        ("F4_q4", f4, 4, alpha), ("F9_q3", f9, 3, alpha),
        ("F9_q9", f9, 9, alpha),
        ("F3_theta", function_field(prime_field(3), theta), 3, theta),
        ("F4_theta", function_field(f4, theta), 4, theta),
        ("F9_theta_q3", function_field(f9, theta), 3, theta),
        ("F9_theta_q9", function_field(f9, theta), 9, theta)]


class CrossFieldTests(unittest.TestCase):
    pass


def _make_field_test(f, q, th, case):
    def check(self):
        p = Parameters(f, q, th, tau)
        if case == "running":
            phi, psi = [[th+tau**2, 0], [tau**3, th+tau**2+tau**3]], [[th+tau]]
        elif case == "large":
            phi = [[th+tau**4, tau], [0, th+tau**3]]
            psi = [[th+th**2*tau, 0], [1+tau, th+th*tau]]
        else:
            phi = [[th+tau**2, 1], [0, th+tau**2]]
            psi = [[th+tau]]
        d = ReductionData(phi, psi, p)
        r = d.exact_sequence()
        self.assertTrue(r["ExactSequenceVerified"])
        if case == "running":
            self.assertEqual((r["FullDimension"], r["ExtTau0Dimension"], r["H"]), (5, 3, 2))
            self.assertEqual(r["Sigma"], [2, 1])
        elif case == "large":
            self.assertEqual((r["FullDimension"], r["ExtTau0Dimension"], r["H"]), (14, 12, 2))
            # Every entry of Pi, Pi0, T and P, including nonconstant twists.
            for key, expected in expected_article_matrices(p).items():
                self.assertEqual(d.partial_exact_sequence()[key], expected, msg=key)
            self.assertTrue(r["K"].degree == 0)
        else:
            self.assertEqual((r["FullDimension"], r["ExtTau0Dimension"], r["H"]), (4, 3, 1))
        # Independently check each action column at a nonconstant/nontrivial value.
        x = f(theta+alpha+1) if f.variables and f.extension_modulus else (
            f(theta+1) if f.variables else f(alpha+1) if f.extension_modulus
            else f(1 if f.characteristic == 2 else 2))
        for col, triple in enumerate(d.basis(one_based=False)):
            actual = d.coordinates(d.reduce(d.psi @ d.generator(triple, value=x)))
            vector = [f.zero]*r["FullDimension"]
            vector[col] = x
            self.assertEqual(actual, r["Pi"].evaluate(vector))
        v = OreMatrix.from_matrix([[tau**5+tau, 2*tau**4] for _ in range(d.e)], p)
        reduced = d.reduce_report(v)
        self.assertTrue(reduced["VerifiedDifference"])
    return check


for name, f, q, th in _field_cases():
    for case in ("running", "large", "nilpotent_phi"):
        setattr(CrossFieldTests, f"test_{name}_{case}", _make_field_test(f, q, th, case))


class ArticleAdditionalTests(unittest.TestCase):
    def test_three_dimensional_example_6_1(self):
        p = Parameters(function_field(prime_field(3), theta), 3, theta, tau)
        d = ReductionData([[theta+tau**3, tau, 0], [0, theta+tau**2, 0],
                          [0, 0, theta+tau**2]], [[theta+tau]], p)
        expected = [[theta, 0, 0, 0, 0, 0, 0],
                    [tau, theta, tau**2, 0, 0, 0, 0],
                    [0, tau, theta, 0, 0, 0, 0],
                    [0, 0, 0, theta, 0, 0, 0],
                    [0, 0, -tau, tau, theta+tau**2, 0, 0],
                    [0, 0, 0, 0, 0, theta, 0],
                    [0, 0, 0, 0, 0, tau, theta+tau**2]]
        self.assertTrue(SkewEqualQ(d.extension(), expected, p))
        expected0 = [[theta, tau**2, 0, 0], [tau, theta, 0, 0],
                     [0, -tau, theta+tau**2, 0], [0, 0, 0, theta+tau**2]]
        self.assertTrue(SkewEqualQ(d.ext_tau0(), expected0, p))
        self.assertTrue(d.exact_sequence()["ExactSequenceVerified"])

    def test_random_reduction_witnesses_and_actions(self):
        rng = random.Random(20261009)
        for case in range(18):
            f, pchar = prime_field(3), 3
            p = Parameters(f, pchar, 1, tau)
            n = 2+case % 2
            phi = [[(1 if i == j else 0)+sum(
                    rng.randrange(3)*tau**k for k in range(1, n))
                    +(tau**n if i == j else 0) for j in range(2)] for i in range(2)]
            psi = [[(1 if i == j else 0)+rng.randrange(3)*tau for j in range(2)] for i in range(2)]
            d = ReductionData(phi, psi, p)
            pi = d.extension()
            v = [[sum(rng.randrange(3)*tau**k for k in range(n+3)) for _ in range(2)] for _ in range(2)]
            r = d.reduce_report(v)
            self.assertTrue(r["VerifiedDifference"])
            self.assertEqual(d.reduce(r["ReducedMatrix"]), r["ReducedMatrix"])
            x = [f(rng.randrange(3)) for _ in d.basis()]
            actual = d.coordinates(d.reduce(d.psi @ d.representative(x)))
            self.assertEqual(actual, pi.evaluate(x))


if __name__ == "__main__":
    unittest.main(verbosity=2)
