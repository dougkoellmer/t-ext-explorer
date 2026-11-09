"""Regressions for the article's new Algorithms 10--13."""
import unittest
import sympy as sp
from fj_sympy import *
from fj_sympy.linear import identity, multiply

th, tau, a, b = sp.symbols('theta tau a b')


class ExactSequenceTests(unittest.TestCase):
    def nilpotent_pair(self, parameters, degree=2, positive_target=False, leading=None):
        theta = parameters.field.expr(parameters.theta)
        lead = sp.eye(2) if leading is None else sp.Matrix(leading)
        phi = sp.Matrix([[theta, 1], [0, theta]]) + lead*tau**degree
        psi = sp.Matrix([[theta, 1], [0, theta]]) + (sp.eye(2)*tau if positive_target else sp.zeros(2))
        return ReductionData(phi, psi, parameters)

    def test_column_vectorization_and_quotient(self):
        p = Parameters(prime_field(3), 3, 1, tau)
        d = self.nilpotent_pair(p)
        info = FJConstantData(d)
        J = sp.Matrix([[0, 1], [0, 0]])
        expected = sp.kronecker_product(J.T, sp.eye(2))-sp.kronecker_product(sp.eye(2), J)
        self.assertEqual(info['D'], [[p.field(v) for v in row] for row in expected.tolist()])
        self.assertEqual((info['H'], info['ConstantMapRank'], info['PseudocodeZeroCount']), (2, 2, 1))
        self.assertEqual(field_matrix(info['QuotientMatrix'], p.field), sp.Matrix([[1, 0], [1, 1]]))
        self.assertEqual(GaPart(d)[0], 2)
        self.assertEqual(FJPartialGaPart(d)[0], 1)
        r = FJExactSequence(d)
        self.assertEqual((r['FullDimension'], r['ExtTau0Dimension'], r['H']), (8, 6, 2))
        self.assertTrue(r['ExactSequenceVerified'])
        self.assertFalse(r['ScalarAdditiveQuotient'])

    def test_nonscalar_quotient_across_fields(self):
        alpha = sp.Symbol('alpha')
        fields = [(prime_field(2), 2, 1), (prime_field(3), 3, 1),
                  (prime_field(5), 5, 1), (function_field(finite_field(3, alpha**2+1, alpha), th), 9, th)]
        for f, q, theta in fields:
            with self.subTest(q=q):
                p = Parameters(f, q, theta, tau)
                d = self.nilpotent_pair(p, positive_target=True, leading=[[1, 1], [0, 1]])
                r = d.exact_sequence()
                self.assertEqual((r['FullDimension'], r['ExtTau0Dimension'], r['H']), (8, 6, 2))
                self.assertTrue(r['AndersonExt0Verified'] and r['AndersonQuotientVerified'])
                self.assertEqual(r['Pi'] @ r['K'], r['K'] @ r['Pi0'])
                self.assertEqual(r['P'] @ r['Pi'], r['QuotientMatrix'] @ r['P'])

    def test_empty_blocks(self):
        p = Parameters(prime_field(3), 3, 1, tau)
        r = self.nilpotent_pair(p, degree=1).exact_sequence()
        self.assertEqual((r['FullDimension'], r['ExtTau0Dimension'], r['H']), (4, 2, 2))
        r = FJExactSequence([[1+tau]], [[1]], p)
        self.assertEqual((r['Pi0'].shape, r['K'].shape, r['R'].shape), ((0, 0), (1, 0), (0, 1)))
        self.assertTrue(r['ExactSequenceVerified'])

    def test_removable_constant_and_two_witnesses(self):
        p = Parameters(prime_field(3), 3, 1, tau)
        d = self.nilpotent_pair(p)
        v = OreMatrix.from_matrix([[0, 1], [0, 0]], p)
        r = FJReduceExtTau0Report(v, d)
        self.assertEqual(r['ReducedMatrix'], v)
        self.assertEqual(r['ZeroConstantRepresentative'], OreMatrix.from_matrix([[2*tau**2, 0], [0, 0]], p))
        self.assertEqual(v-r['ReducedMatrix'], d.inner(r['InnerWitness']))
        self.assertEqual(v-r['ZeroConstantRepresentative'], d.inner(r['ZeroConstantInnerWitness']))
        self.assertTrue(r['ZeroConstantWitnessVerified'])
        with self.assertRaises(FJError) as error:
            d.reduce_zero([[0, 0], [1, 0]])
        self.assertEqual(error.exception.code, 'ClassNotInExtZero')

    def test_formal_parameter_representative(self):
        p = FJSymbolicParameters(th, tau)
        d = ReductionData([[th+tau**2, a], [0, th+tau**2]], [[th+b*tau]], p)
        v = d.generator((0, 1, 0), formal=True)
        r = d.reduce_zero_report(v)
        self.assertTrue(r['VerifiedDifference'] and r['ZeroConstantWitnessVerified'])
        self.assertEqual(v-r['ZeroConstantRepresentative'], d.inner(r['ZeroConstantInnerWitness']))

    def test_nonfixed_inclusion_requires_twists(self):
        p = FJSymbolicParameters(th, tau)
        d = ReductionData([[th+tau**2, a], [0, th+tau**2]], [[th+b*tau]], p)
        r = d.exact_sequence()
        self.assertEqual(r['Pi0'], r['R'] @ r['Pi'] @ r['K'])
        R, K, f = r['R'].coefficient(0), r['K'].coefficient(0), p.field
        wrong = OreMatrix([multiply(multiply(R, C, f), K, f) for C in r['Pi'].coefficients],
                          r['Pi0'].rows, r['Pi0'].cols, p)
        self.assertNotEqual(r['Pi0'], wrong)
        self.assertTrue(r['Pi0'].display().has(FJTwist(a, 1)))

    def test_symbolic_nonscalar_specialization(self):
        p = FJSymbolicParameters(th, tau)
        phi, psi = [[th+tau**2, a], [0, th+tau**2]], [[th+tau, b], [0, th+tau]]
        r = FJExactSequence(phi, psi, p)
        for prime in (2, 3, 5):
            with self.subTest(prime=prime):
                cp = Parameters(function_field(prime_field(prime), th), prime, th, tau)
                assignments = {a: th, b: th+1}
                direct = FJExactSequence(SpecializeFrobenius(r["Data"].phi, cp, assignments),
                                         SpecializeFrobenius(r["Data"].psi, cp, assignments), cp)
                for key in ['Pi', 'Pi0', 'K', 'R', 'P', 'QuotientMatrix']:
                    self.assertEqual(SpecializeFrobenius(r[key], cp, assignments), direct[key])

    def test_partial_example_is_preserved_and_equivalent(self):
        p = FJSymbolicParameters(th, tau)
        d = ReductionData([[th+tau**4, tau], [0, th+tau**3]],
                          [[th+th**2*tau, 0], [1+tau, th+th*tau]], p)
        r, old = d.exact_sequence(), d.partial_exact_sequence()
        self.assertIs(r['Pi'], old['Pi'])
        self.assertEqual(r['K'].degree, 0)
        self.assertEqual(old['T'].degree, 1)
        W = r['R'] @ old['T']
        V = old['InclusionLeftInverse'] @ r['K']
        self.assertEqual(r['K'] @ W, old['T'])
        self.assertEqual(r['Pi0'] @ W, W @ old['Pi0'])
        ident = OreMatrix.monomial(identity(p.field, 12), 0, p)
        self.assertEqual(W @ V, ident)
        self.assertEqual(V @ W, ident)

    def test_algorithm_entry_points_and_bad_action_shape(self):
        p = FJSymbolicParameters(th, tau)
        d = ReductionData([[th+tau**2]], [[th+tau]], p)
        c = FJSequenceCoordinates(d)
        acts = FJSequenceActions(FJExtension(d), d)
        self.assertEqual(acts['Pi0'], FJExtTau0(d))
        self.assertEqual(c['K'], FJTransition(d))
        with self.assertRaises(FJError):
            FJSequenceActions([[th]], d)
