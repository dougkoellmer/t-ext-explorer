"""Formal Frobenius laws and independent concrete specialization checks."""
import pickle
import unittest
import sympy as sp
from fj_sympy import *

th, tau, a, b, c, x = sp.symbols('θ τ a b c x')


class FormalFrobeniusTests(unittest.TestCase):
    def setUp(self):
        self.p = FJSymbolicParameters(th, tau)
        self.f = self.p.field

    def test_twist_composition(self):
        self.assertEqual(FJFrobenius(FJFrobenius(a, 2, self.p), 3, self.p), FJTwist(a, 5))
        self.assertEqual(FJTwist(FJTwist(a, 2), 3), FJTwist(a, 5))

    def test_identity(self):
        self.assertEqual(FJTwist(a, 0), a)
        self.assertEqual(FJFrobenius((a+b)/(1-a), 0, self.p), (a+b)/(1-a))

    def test_rational_homomorphism(self):
        self.assertEqual(FJFrobenius((th-a)/(1+b), 2, self.p),
                         (FJTwist(th, 2)-FJTwist(a, 2))/(1+FJTwist(b, 2)))

    def test_product_inverse_and_power(self):
        self.assertEqual(FJFrobenius(a**2/(b-th), 3, self.p),
                         FJTwist(a, 3)**2/(FJTwist(b, 3)-FJTwist(th, 3)))

    def test_theta_twists(self):
        self.assertEqual(FJFrobenius(th, 1, self.p), FJTwist(th, 1))
        self.assertNotEqual(FJFrobenius(th, 1, self.p), th)

    def test_integer_constants(self):
        self.assertEqual(FJFrobenius(sp.Rational(2, 3), 5, self.p), sp.Rational(2, 3))

    def test_characteristic_arithmetic(self):
        p = FJSymbolicParameters(th, tau, characteristic=3)
        self.assertEqual(p.field(3*a), 0)
        self.assertEqual(p.field(sp.Rational(1, 2)), -1)
        self.assertEqual(p.field(2*a+2*a), a)

    def test_forbidden_denominator(self):
        p = FJSymbolicParameters(th, tau, characteristic=2)
        with self.assertRaises(FJError):
            p.field(sp.Rational(1, 2))

    def test_fixed_coefficients(self):
        p = FJSymbolicParameters(th, tau, fixed_symbols=(b,))
        self.assertEqual(FJFrobenius(a+b, 2, p), FJTwist(a, 2)+b)
        self.assertEqual(p.field(FJTwist(b, 3)), b)

    def test_rational_normalization(self):
        self.assertEqual(self.f((a*a-b*b)/(a-b)), a+b)
        self.assertEqual(self.f(2*a)/self.f(2*b), a/b)
        u = FJTwist(a, 2)
        self.assertEqual(self.f((u*u-b*b)/(u-b)), u+b)

    def test_exact_zero_and_nonzero(self):
        self.assertFalse(self.f((a*a-b*b)/(a-b)-a-b))
        self.assertTrue(self.f(FJTwist(a, 1)-a))

    def test_latex_plain_and_pickle(self):
        atom = FJTwist(th, 2)
        self.assertEqual(str(atom), 'θ^(2)')
        self.assertEqual(sp.latex(atom), r'θ^{(2)}')
        self.assertEqual(pickle.loads(pickle.dumps(atom)), atom)
        # Symbols sharing a printed name but having different assumptions stay distinct.
        self.assertNotEqual(FJTwist(sp.Symbol('a', positive=True), 1), FJTwist(a, 1))
        self.assertNotEqual(FJTwist(a, 1), sp.Symbol('a__FJtwist_1'))

    def test_invalid_inputs(self):
        for r in [-1, True, sp.Rational(1, 2), sp.Symbol('r')]:
            with self.assertRaises(FJError):
                FJFrobenius(a, r, self.p)
        for value in [0.5, sp.sin(a), sp.sqrt(a), tau, FJTwist(tau, 1)]:
            with self.assertRaises(FJError):
                self.f(value)
        with self.assertRaises(FJError):
            FJSymbolicParameters(characteristic=9)
        with self.assertRaises(FJError):
            FJTwist(a+b, 2)

    def test_skew_rule(self):
        product = TMult([[tau]], [[a*tau]], self.p)
        self.assertEqual(product.display(), sp.Matrix([[FJTwist(a, 1)*tau**2]]))

    def test_skew_associativity(self):
        A, B, C = [OreMatrix.from_matrix(v, self.p) for v in
                   ([[a+tau, b*tau]], [[th*tau], [1+a*tau**2]], [[b+tau]])]
        self.assertEqual((A @ B) @ C, A @ (B @ C))

    def test_formal_parameter_display(self):
        A = OreMatrix([[[(self.f(a), self.f(b))]]], 1, 1, self.p, formal=True)
        self.assertEqual(A.display(), sp.Matrix([[a*c+b*FJTwist(c, 1)]]))
        self.assertEqual(A.evaluate_formal(x).display(),
                         sp.Matrix([[a*x+b*FJTwist(x, 1)]]))

    def test_conditions_excluded_primes(self):
        self.f(2).inverse()
        self.f(3*a).inverse()
        conditions = FJSymbolicConditions(self.p)
        self.assertEqual(conditions['ExcludedPrimes'], [2, 3])
        self.assertIn(3*a, conditions['NonzeroConditions'])

    def test_specialization_denominator_vanishes(self):
        self.f(FJTwist(a, 1)-a).inverse()
        concrete = Parameters(function_field(prime_field(3), th), 3, th, tau)
        bad = FJSpecializationReport(self.p, concrete, {a: 1})
        good = FJSpecializationReport(self.p, concrete, {a: th})
        self.assertFalse(bad['RecordedConditionsSatisfied'])
        self.assertTrue(good['RecordedConditionsSatisfied'])

    def test_fixed_symbol_specialization(self):
        p = FJSymbolicParameters(th, tau, fixed_symbols=(b,))
        concrete = Parameters(function_field(prime_field(3), th), 3, th, tau)
        self.assertFalse(FJSpecializationReport(p, concrete, {b: th})['RecordedConditionsSatisfied'])
        self.assertTrue(FJSpecializationReport(p, concrete, {b: 1})['RecordedConditionsSatisfied'])

    def test_characteristic_mismatch(self):
        p = FJSymbolicParameters(th, tau, characteristic=3)
        concrete = Parameters(function_field(prime_field(5), th), 5, th, tau)
        self.assertFalse(FJSpecializationReport(p, concrete)['RecordedConditionsSatisfied'])

    def test_specialization_commutes_with_frobenius(self):
        value = (th-a)/(1+b)
        for prime, q, base in [(2, 2, prime_field(2)), (3, 3, prime_field(3)),
                               (3, 9, finite_field(3, x**2+1, x)), (5, 5, prime_field(5))]:
            with self.subTest(q=q):
                concrete = Parameters(function_field(base, th), q, th, tau)
                assignments = {a: th+1, b: th**2}
                specialized = SpecializeFrobenius(value, concrete, assignments)
                self.assertEqual(SpecializeFrobenius(FJFrobenius(value, 2, self.p),
                                                    concrete, assignments),
                                 concrete.frobenius(specialized, 2))

    def test_specialization_commutes_with_ore(self):
        A = OreMatrix.from_matrix([[a+th*tau, b*tau]], self.p)
        B = OreMatrix.from_matrix([[tau], [a+th*tau]], self.p)
        concrete = Parameters(function_field(prime_field(3), th), 3, th, tau)
        assignments = {a: th+1, b: th**2}
        self.assertEqual(SpecializeFrobenius(A @ B, concrete, assignments),
                         SpecializeFrobenius(A, concrete, assignments) @
                         SpecializeFrobenius(B, concrete, assignments))

    def test_twisted_assignments_rejected(self):
        concrete = Parameters(function_field(prime_field(3), th), 3, th, tau)
        with self.assertRaises(FJError):
            SpecializeFrobenius(a, concrete, {FJTwist(a, 1): 1})


class SymbolicIntegrationTests(unittest.TestCase):
    def test_running_exact_sequence(self):
        p = FJSymbolicParameters(th, tau)
        phi = [[th+tau**2, 0], [tau**3, th+tau**2+tau**3]]
        psi = [[th+tau]]
        d = ReductionData(phi, psi, p)
        report = FJExactSequence(d)
        self.assertEqual((report['FullDimension'], report['ExtTau0Dimension'], report['H']), (5, 3, 2))
        self.assertTrue(report['ExactSequenceVerified'])
        red = FJReduceReport([[a*tau**5+b, th*tau**4]], d)
        self.assertTrue(red['VerifiedDifference'])
        self.assertEqual(d.reduce(red['ReducedMatrix']), red['ReducedMatrix'])

    def test_symbolic_pivots_and_reduction(self):
        p = FJSymbolicParameters(th, tau)
        d = ReductionData([[th+a*tau**3]], [[th+b*tau]], p)
        red = FJReduceReport([[c*tau**5]], d)
        self.assertTrue(red['VerifiedDifference'])
        self.assertEqual(d.reduce(red['ReducedMatrix']), red['ReducedMatrix'])
        self.assertTrue(FJExactSequence(d)['ExactSequenceVerified'])
        concrete = Parameters(function_field(prime_field(3), th), 3, th, tau)
        assignments = {a: th+1, b: th, c: th**2}
        dd = ReductionData(SpecializeFrobenius(d.phi, concrete, assignments),
                           SpecializeFrobenius(d.psi, concrete, assignments), concrete)
        direct = dd.reduce([[th**2*tau**5]])
        self.assertEqual(SpecializeFrobenius(red['ReducedMatrix'], concrete, assignments), direct)

    def test_fourteen_dimensional_specialization(self):
        p = FJSymbolicParameters(th, tau)
        phi = [[th+tau**4, tau], [0, th+tau**3]]
        psi = [[th+th**2*tau, 0], [1+tau, th+th*tau]]
        d = ReductionData(phi, psi, p)
        report = FJExactSequence(d)
        self.assertEqual((report['FullDimension'], report['ExtTau0Dimension'], report['H']), (14, 12, 2))
        self.assertTrue(report['ExactSequenceVerified'])
        self.assertTrue(report['Pi'].display().has(FJTwist(th, 1)))
        for prime in [2, 3, 5]:
            with self.subTest(prime=prime):
                concrete = Parameters(function_field(prime_field(prime), th), prime, th, tau)
                direct = FJExactSequence(phi, psi, concrete)
                self.assertTrue(FJSpecializationReport(p, concrete)['RecordedConditionsSatisfied'])
                for key in ['Pi', 'Pi0', 'T', 'P']:
                    self.assertEqual(SpecializeFrobenius(report[key], concrete), direct[key])

    def test_fixed_characteristic_exact_sequence(self):
        p = FJSymbolicParameters(th, tau, characteristic=3)
        d = ReductionData([[th+tau**2, 0], [tau**3, th+tau**2+tau**3]], [[th+tau]], p)
        self.assertTrue(FJExactSequence(d)['ExactSequenceVerified'])

    def test_general_nilpotent_branch(self):
        p = FJSymbolicParameters(th, tau)
        d = ReductionData([[th+tau**2, 1], [0, th+tau**2]], [[th, 1], [0, th]], p)
        report = FJGeneralReport(d)
        self.assertTrue(report['InclusionVerified'])
        self.assertTrue(report['QuotientVerified'])
