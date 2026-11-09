"""Verified Mathematica J2 example and descriptive report regressions."""
import unittest
import sympy as sp
import fj_sympy as fj
from fj_sympy import FJExampleCompute,FJExampleReport,FJPartialBasisComparison,FJPartialBasisReport,FJJ2Example
_InputError=fj.FJError

def _symbols():return sp.symbols('theta tau c')
def _parameters(t,u):return fj.FJSymbolicParameters(theta=t,tau=u)
def _plain(x):return str(x.display())

class ExampleReportTests(unittest.TestCase):
    def setUp(self):
        self.theta,self.tau,self.c=_symbols()
        self.p=_parameters(self.theta,self.tau)

    def test_j2_actions_match_verified_mathematica(self):
        t,u=self.theta,self.tau;r=FJJ2Example(t,u);s=r['Result']['Sequence'];q=r['Summary']
        self.assertEqual(tuple(q[k] for k in ['FullDimension','TrueExt0Dimension','H','IndividualZeroCount','InvalidCountPrediction']),(8,6,2,1,7))
        self.assertEqual(q['PartialBasisStatus'],'Unavailable')
        pi=[[t,0,0,0,1,0,0,0],[0,t,0,0,u,1,0,0],[0,0,t,0,0,-u,1,0],[0,0,0,t,0,0,u,1],[0,0,0,0,t,0,0,0],[0,0,0,0,0,t,0,0],[0,0,0,0,0,0,t,0],[0,0,0,0,0,0,0,t]]
        pi0=[[t,1,0,0,-u,0],[0,t,0,0,0,0],[0,0,t,0,1,0],[0,u,0,t,0,1],[0,0,0,0,t,0],[0,0,0,0,0,t]]
        p=s['Parameters']
        self.assertTrue(fj.SkewEqualQ(s['Pi'],pi,p));self.assertTrue(fj.SkewEqualQ(s['Pi0'],pi0,p))
        self.assertTrue(fj.SkewEqualQ(s['QuotientMatrix'],[[t,0],[1,t]],p))
        self.assertTrue(s['ExactSequenceVerified'])

    def test_formal_frobenius_dependence_and_zero_constant_witness(self):
        t,u,c=self.theta,self.tau,self.c
        data=fj.FJReductionData([[t+u**2,1],[0,t+u**2]],[[t,1+u],[0,t]],self.p)
        twist=fj.FJTwist(c,1)
        v=[[c*u**2,-twist*u],[0,c*u**2]]
        self.assertTrue(fj.SkewEqualQ(fj.FJReduce(v,data),[[0,0],[0,0]],self.p))
        w=fj.FJReduceExtTau0Report([[-c,0],[0,c]],data)
        self.assertTrue(w['ZeroConstantWitnessVerified'])
        self.assertTrue(fj.SkewEqualQ(w['ZeroConstantRepresentative'],[[twist*u,0],[-c*u**2,0]],self.p))

    def test_partial_basis_polynomial_isomorphism_for_large_pair(self):
        t,u=self.theta,self.tau
        result=FJExampleCompute([[t+u**4,u],[0,t+u**3]],[[t+t**2*u,0],[1+u,t+t*u]],self.p)
        comparison=FJPartialBasisComparison(result['Data'],result['Sequence'])
        self.assertEqual(comparison['Status'],'Available');self.assertTrue(comparison['SameSubgroupVerified'])
        self.assertIn('tau',_plain(comparison['W']))
        view=FJPartialBasisReport(comparison)
        self.assertIn('W^(-1)',view.text)

    def test_zero_dimensional_ext0_preserves_empty_matrix_shapes(self):
        t,u=self.theta,self.tau;r=FJExampleCompute([[t+u]],[[t]],self.p)
        self.assertEqual(r['Sequence']['ExtTau0Dimension'],0)
        view=FJExampleReport(r)
        self.assertIn('Matrix size: 0 x 1',view.text)
        self.assertIn('Empty matrix (0 x 1)',view._repr_html_())

    def test_invalid_input_retains_mathematical_exception(self):
        t,u=self.theta,self.tau
        with self.assertRaises(_InputError):FJExampleCompute([[t]],[[t+u]],self.p)

    def test_own_scalar_example_and_raw_inputs(self):
        t,u=self.theta,self.tau;phi=[[t+u**2]];psi=[[t+u]]
        r=FJExampleCompute(phi,psi,self.p)
        self.assertIs(r['InputPhi'],phi);self.assertIs(r['InputPsi'],psi)
        self.assertEqual(tuple(r['Sequence'][k] for k in ['FullDimension','ExtTau0Dimension','H']),(2,1,1))
        self.assertTrue(r['Sequence']['ExactSequenceVerified'])

    def test_rich_report_explains_expected_partial_rejection(self):
        r=FJJ2Example(self.theta,self.tau);view=r['Output'];html=view._repr_html_()
        self.assertIn('Frobenius',view.text);self.assertIn('rank D',view.text)
        self.assertIn('linearly dependent',view.text);self.assertIn('τ',html)
        self.assertNotIn('history',view.text.lower());self.assertNotIn('Î',html)
        self.assertIn('Why partial basis is unavailable',html)


if __name__=='__main__':unittest.main()
