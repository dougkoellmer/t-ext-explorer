"""Automatic p,n fields, Frobenius/Ore laws, and reuse after a field change."""
import unittest
import sympy as sp
import fj_sympy as fj
from fj_sympy.example_reports import FJExampleCompute, FJJ2Example

class AutomaticFieldTests(unittest.TestCase):
    def setUp(self):
        self.theta, self.tau, self.alpha, self.u = sp.symbols('theta tau alpha u')

    def test_fields_and_generator_relations(self):
        for p,n in [(2,1),(3,1),(5,1),(2,2),(2,3),(3,2),(5,2),(3,3)]:
            with self.subTest(p=p,n=n):
                par=fj.FJParametersPN(p,n,self.theta,self.tau,self.alpha)
                self.assertEqual(par.q,p**n)
                self.assertEqual(par.field.characteristic,p)
                self.assertFalse(par.field(p))
                self.assertEqual(par.field.constant_degree,n)
                if n>1:
                    poly=par.field.extension_modulus
                    self.assertEqual(poly.degree(),n)
                    self.assertTrue(poly.is_irreducible)
                    self.assertFalse(par.field(poly.as_expr()))
                    self.assertEqual(par.frobenius(self.alpha),par.field(self.alpha))
                    self.assertEqual(fj.finite_field_pn(p,n,self.alpha).extension_modulus,poly)

    def test_rational_frobenius_and_ore_multiplication(self):
        for p,n in [(2,2),(3,2),(5,2)]:
            with self.subTest(p=p,n=n):
                par=fj.FJParametersPN(p,n,self.theta,self.tau,self.alpha,
                                     coefficient_symbols=(self.u,))
                c=(self.alpha+self.theta)/(self.theta+self.u)
                product=fj.TMult([[self.tau]],[[c]],par)
                expected=par.field(c)**par.q
                self.assertEqual(product.coefficient(1)[0][0],expected)
                self.assertEqual(par.frobenius(self.u),par.field(self.u**par.q))

    def test_change_field_recomputes_example_and_j2(self):
        for p,n in [(3,1),(2,2),(3,2),(5,1)]:
            with self.subTest(p=p,n=n):
                par=fj.FJParametersPN(p,n,self.theta,self.tau,self.alpha)
                result=FJExampleCompute([[self.theta+self.tau**2]],[[self.theta+self.tau]],par)
                self.assertTrue(result['Sequence']['ExactSequenceVerified'])
                self.assertEqual(result['Data'].parameters.q,p**n)
                j2=FJJ2Example(self.theta,self.tau,parameters=par)
                self.assertEqual(j2['Summary']['TrueExt0Dimension'],6)
                self.assertTrue(j2['Result']['Sequence']['ExactSequenceVerified'])

    def test_invalid_settings_are_rejected(self):
        for p,n in [(4,1),(9,2),(True,2),(3,0),(3,-1),(3,1.5),(3,True)]:
            with self.subTest(p=p,n=n):
                with self.assertRaises(fj.FJError):
                    fj.FJParametersPN(p,n,self.theta,self.tau,self.alpha)

if __name__=='__main__':
    unittest.main()
