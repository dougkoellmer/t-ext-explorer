import sympy as sp
import fj_sympy as fj
_InputError = fj.FJError

def _parts(data):
    return data.parameters, data.e, data.d

def _nil(data, which):
    return data.nphi if which == 'Phi' else data.npsi

def _symbolic(theta,tau):
    return fj.FJSymbolicParameters(theta=theta,tau=tau)

def _diagnostics(data):
    if isinstance(data.parameters,fj.SymbolicParameters):
        return {'Scope':'Recorded denominators and inverted pivots; conservative conditions, not a complete classification of all q.',**fj.FJSymbolicConditions(data.parameters)}
    return {'Scope':'The specified concrete field and q.','Characteristic':data.field.characteristic,'q':data.parameters.q}

def _product(a,b,p):
    return fj.TMult(a,b,p)

def _equal(a,b,p):
    return fj.SkewEqualQ(a,b,p)

def _identity(n,p):
    return fj.OreMatrix.from_matrix(sp.eye(n),p)

def _unit(data,i,s,k):
    out=sp.zeros(data.e,data.d);out[i-1,s-1]=data.parameters.tau**k
    return out

def _matrix_data(value,data,shape=None):
    if isinstance(value,fj.OreMatrix):value=value.display()
    if isinstance(value,sp.MatrixBase):return int(value.rows),int(value.cols),value.tolist()
    rows=list(value);n,m=shape or (len(rows),len(rows[0]) if rows else 0)
    def expr(x):
        return data.field.expr(x) if isinstance(x,(fj.Element,fj.SymbolicElement)) else x
    return int(n),int(m),[[expr(x) for x in row] for row in rows]

"""English example reports. Arithmetic is supplied by the package adapter.

Raw matrices and exceptions remain available in the returned dictionaries.
The renderer uses HTML and plain text, without changing mathematical objects.
"""
from html import escape
import re

__all__ = ['FJExampleCompute', 'FJExampleReport', 'FJPartialBasisComparison',
           'FJPartialBasisReport', 'FJJ2Example', 'ReportView']


def _notation(value):
    text = escape(str(value))
    text = re.sub(r'\b(theta|tau)\b', lambda m: {'theta':'θ','tau':'τ'}[m[0]], text)
    text = re.sub(r'FJTwist\((\w+),\s*(\d+)\)', lambda m: m[1]+'<sup>('+m[2]+')</sup>', text)
    text = re.sub(r'\^\((\d+)\)', lambda m: '<sup>('+m[1]+')</sup>', text)
    text = re.sub(r'(?:\*\*|\^)(-?\d+)', lambda m: '<sup>'+m[1]+'</sup>', text)
    return text.replace('*', '·')


class ReportView:
    """Rich Jupyter output and readable plain text for the same report."""
    def __init__(self):
        self._html, self._text = [], []

    def heading(self, text):
        self._html.append('<h3>'+escape(text)+'</h3>')
        self._text.append('\n'+text)
        return self

    def paragraph(self, text):
        self._html.append('<p>'+escape(str(text))+'</p>')
        self._text.append(str(text))
        return self

    def table(self, rows, header=None):
        rows = ([header] if header else [])+list(rows)
        body = []
        for i, row in enumerate(rows):
            tag = 'th' if header and i == 0 else 'td'
            body.append('<tr>'+''.join('<'+tag+' style="padding:5px 10px;border:1px solid #aaa;text-align:left">'
                                      +escape(str(v))+'</'+tag+'>' for v in row)+'</tr>')
            self._text.append(' | '.join(str(v) for v in row))
        self._html.append('<div style="overflow-x:auto"><table style="border-collapse:collapse">'+''.join(body)+'</table></div>')
        return self

    def matrix(self, label, description, value, data, shape=None):
        rows, cols, entries = _matrix_data(value, data, shape)
        self.heading(label).paragraph(description).paragraph(f'Matrix size: {rows} x {cols}')
        self._html.append(_matrix_html(rows, cols, entries))
        self._text += ['['+', '.join(str(v) for v in row)+']' for row in entries]
        return self

    def basis(self, basis, data, zero=False):
        rows = []
        for number, term in enumerate(basis, 1):
            if zero and term['Kind'] == 'ConstantImage':
                label, value = 'Constant in im D', term['Matrix']
            else:
                i, s, k = term['Position'] if zero else term
                label, value = f'E({i},{s}) tau^{k}', _unit(data, i, s, k)
            n, m, entries = _matrix_data(value, data)
            rows.append('<tr><td>'+str(number)+'</td><td>'+escape(label)+'</td><td>'+_matrix_html(n,m,entries)+'</td></tr>')
            self._text.append(f'{number}: {label}: '+str(entries))
        self._html.append('<div style="overflow-x:auto"><table class="fj-basis"><tr><th>Coordinate</th><th>Kind / index</th><th>Representative matrix</th></tr>'+''.join(rows)+'</table></div>')
        return self

    @property
    def text(self):
        return '\n'.join(self._text)

    def __repr__(self):
        return self.text

    def _repr_html_(self):
        return '<div class="fj-example-report">'+''.join(self._html)+'</div>'


def _matrix_html(rows, cols, entries):
    if not rows or not cols:
        return '<p>Empty matrix ('+str(rows)+' x '+str(cols)+')</p>'
    body = ''.join('<tr>'+''.join('<td style="padding:4px 10px;white-space:nowrap;text-align:center">'+_notation(v)+'</td>' for v in row)+'</tr>' for row in entries)
    return '<div style="overflow-x:auto"><table style="border-left:2px solid;border-right:2px solid;margin:6px 0">'+body+'</table></div>'


def FJExampleCompute(phi, psi, parameters):
    """Compute your input pair; mathematical errors retain their native type."""
    data = fj.FJReductionData(phi, psi, parameters)
    # ===== BEGIN: computation of the true Ext_0 and exact sequence =====
    sequence = fj.FJExactSequence(data)
    # ===== END: computation of the true Ext_0 and exact sequence =====
    return {'InputPhi': phi, 'InputPsi': psi, 'Data': data,
            'Sequence': sequence, 'Characteristics': _diagnostics(data)}


def FJExampleReport(result):
    """Inputs, indexed bases, matrix roles/shapes and exact checked identities."""
    if not isinstance(result, dict) or 'Sequence' not in result:
        raise ValueError('First compute a result with FJExampleCompute.')
    s, data = result['Sequence'], result['Data']
    p, e, d = _parts(data)
    a, m, m0, h = s['ConstantAnalysis'], s['FullDimension'], s['ExtTau0Dimension'], s['H']
    v = ReportView().heading('Input and conventions')
    v.matrix('Phi_t', 'Input: the action of t on Phi.', result['InputPhi'], data)
    v.matrix('Psi_t', 'Input: the action of t on Psi.', result['InputPsi'], data)
    v.paragraph(f'd = dim Phi = {d}; e = dim Psi = {e}. Coordinates are column vectors.')
    v.paragraph('An entry c tau^r is the additive map x |-> c x^(q^r). FJTwist[c,r] denotes the formal coefficient c^(q^r). Use TMult for Ore products, including Frobenius twists.')
    v.matrix('C: compatible basis', 'Columns are selected compatible basis vectors. Jump degrees give the normal degree bounds.', s['BasisMatrix'], data)
    v.table([('sigma',s['Sigma']),('Jump degrees j_s',s['JumpDegrees']),('Selected vectors / blocks',s['CompatibleBasis'])])
    v.heading('Ext: full normal coordinates').paragraph(f'M = e * sum j_s = {m}. E(i,s) has a single 1 at (i,s); retain 0 <= r < j_s. The numbering is the row/column order of Pi.')
    v.basis(s['FullBasis'], data)
    v.matrix('Pi: action of t on Ext', 'Column j is the image of coordinate j. Powers of tau encode Frobenius powers of that coordinate.', s['Pi'], data)
    v.heading('True Ext0: kernel of the constant obstruction').paragraph('Ext0 consists of classes admitting a representative with zero constant term. A normal representative may have a nonzero constant in im D, removable by an inner biderivation.')
    v.matrix('N_Phi','Constant coefficient of Phi_t minus theta I_d.',_nil(data,'Phi'),data)
    v.matrix('N_Psi','Constant coefficient of Psi_t minus theta I_e.',_nil(data,'Psi'),data)
    v.paragraph('T(H) = H N_Phi - N_Psi H. Column vectorization stacks H11,...,He1,H12,...,Hed. D = Transpose[N_Phi] tensor I_e - I_d tensor N_Psi.')
    v.matrix('D: removable constants','D vec(H) = vec(T(H)); constants modulo im D give the obstruction to membership in Ext0.',s['D'],data)
    v.table([('ed',e*d),('rank D',a['ConstantMapRank']),('h = ed-rank D',h),('dim Ext0 = M-h',m0),('Quotient dimension',h)])
    for key, desc, shape in [
        ('B','Independent columns of D, forming a basis of im D.',(e*d,a['ConstantMapRank'])),
        ('Z','Representatives of the constant quotient; [B | Z] is invertible.',(e*d,h)),
        ('U','First rows of [B | Z]^(-1); U B=I and U Z=0.',(a['ConstantMapRank'],e*d)),
        ('P0','Last rows of [B | Z]^(-1); P0 B=0 and P0 Z=I.',(h,e*d))]:
        v.matrix(key,desc,s[key],data,shape)
    v.heading('Ordered Ext0 basis').paragraph('First image-D constants, then all positive normal monomials in their original order. This determines columns of K and rows of R and Pi0.')
    v.basis(s['ExtTau0Basis'], data, zero=True)
    for key, desc, shape in [
        ('Pi0','Pi0 = R o Pi o K using Ore multiplication, in the basis above.',(m0,m0)),
        ('K','x=K y converts Ext0 coordinates to full normal Ext coordinates; rank K=M-h.',(m,m0)),
        ('R','y=R x for x in ker P; R K=I. No t-module morphism property is asserted for R on all of Ext.',(m0,m)),
        ('P','P x=P0 vec(F(0)) for the normal representative F. A class is in Ext0 if and only if P x=0.',(h,m)),
        ('AdditiveSection','S is an additive-group section: P S=I and K R+S P=I. No splitting as t-modules is asserted.',(m,h))]:
        v.matrix('S' if key=='AdditiveSection' else key,desc,s[key],data,shape)
    v.heading('Quotient basis')
    for number, value in enumerate(s['QuotientBasis'],1):
        v.matrix(f'Quotient representative {number}','Constant representative in Mat(e,d)/im D.',value,data)
    v.matrix('Q_t: quotient action','Q_t=theta I_h+P0 (I_d tensor N_Psi) Z. Scalarity is checked rather than assumed.',s['QuotientMatrix'],data,(h,h))
    v.table([('Q_t is theta I_h',s['ScalarAdditiveQuotient'])])
    v.heading('Exact sequence').paragraph(f'0 -> (Ga^{m0},Pi0) --K--> (Ga^{m},Pi) --P--> (Ga^{h},Q_t) -> 0')
    identities=[('K is a morphism','Pi o K = K o Pi0 (Ore)','InclusionVerified'),
        ('P is a morphism','P o Pi = Q_t o P (Ore)','QuotientVerified'),
        ('Composition','P K=0','CompositionZeroVerified'),
        ('Injectivity','R K=I_(M-h)','InjectivityVerified'),
        ('Surjectivity','P S=I_h','SurjectivityVerified'),
        ('Kernel = image','rank K=M-h, rank P=h, K R+S P=I_M','KernelImageVerified'),
        ('Anderson for Ext0',f'(Pi0(0)-theta I)^{e}=0','AndersonExt0Verified'),
        ('Anderson for Q',f'(Q_t-theta I)^{e}=0','AndersonQuotientVerified'),
        ('Overall','All checks above','ExactSequenceVerified')]
    v.heading('Identities and conditions checked').table([(title,identity,s[key]) for title,identity,key in identities],['Check','Identity / condition','Result'])
    v.heading('Characteristic and specialization diagnostics')
    for key, value in result['Characteristics'].items():
        v.paragraph(f'{key}: {value}')
    v.paragraph('Specializing generic coefficients requires preserving the reported ranks, compatible basis and nonvanishing conditions. The diagnostic scope is stated above.')
    return v


def FJPartialBasisComparison(data, sequence=None):
    """Compare partial basis with the true Ext0; retain an explained rejection."""
    s = fj.FJExactSequence(data) if sequence is None else sequence
    a, p = s['ConstantAnalysis'], _parts(data)[0]
    c = {'Data':data,'CorrectSequence':s,'H':a['H'],'PartialZeroCount':a['PseudocodeZeroCount'],
         'Epsilon':a['Epsilon'],'CoordinateKernelCondition':a['CoordinateKernelCondition'],
         'ScalarAdditiveQuotient':a['ScalarAdditiveQuotient']}
    try:
        partial = fj.FJPartialExactSequence(data)
    except _InputError as error:
        c.update(Status='Unavailable', PartialResult=error, Reason=str(error))
        return c
    w = _product(s['R'],partial['T'],p)
    inverse = _product(partial['InclusionLeftInverse'],s['K'],p)
    ident = _identity(s['ExtTau0Dimension'],p)
    ok = all([_equal(s['Pi'],partial['Pi'],p),
        _equal(_product(s['K'],w,p),partial['T'],p),
        _equal(_product(w,inverse,p),ident,p),_equal(_product(inverse,w,p),ident,p),
        _equal(_product(s['Pi0'],w,p),_product(w,partial['Pi0'],p),p)])
    if not ok:
        raise ValueError('The explicit partial basis coordinate isomorphism failed.')
    c.update(Status='Available',PartialResult=partial,W=w,WInverse=inverse,SameSubgroupVerified=ok)
    return c


def FJPartialBasisReport(c):
    v=ReportView().heading('Partial basis: comparison and applicability')
    v.paragraph('For U=E(i,s) C^(-1), epsilon(i,s)=0 if U N_Phi-N_Psi U=0 and 1 otherwise. The candidate partial basis retains 1 <= r < j_s+epsilon(i,s). Full reduction gives the additive-polynomial inclusion T_partial. Counting zero entries misses kernel directions that are linear combinations.')
    s,data=c['CorrectSequence'],c['Data']
    v.matrix('epsilon','A zero records an individual coordinate intertwiner, not all linear combinations.',c['Epsilon'],data)
    v.table([('Number of zero entries',c['PartialZeroCount']),('Correct h',c['H']),('True Ext0 dimension',s['ExtTau0Dimension']),('Status',c['Status'])])
    if c['Status']=='Unavailable':
        v.table([('Kernel spanned by epsilon=0 directions',c['CoordinateKernelCondition']),('Q is scalar',c['ScalarAdditiveQuotient'])])
        tag=getattr(c['PartialResult'],'code',getattr(c['PartialResult'],'tag',''))
        reason=('The nonzero images of individual directions are linearly dependent. The direction count does not determine the full rank or kernel dimension.' if tag=='ArticleCoordinateConditionFailed' else 'The quotient action is nonscalar; this partial basis construction requires a scalar quotient.' if tag=='ArticleScalarQuotientConditionFailed' else c['Reason'])
        v.heading('Why partial basis is unavailable').paragraph(reason)
        v.paragraph('The general Ext0 result remains valid. Raw applicability details are in comparison["PartialResult"]. No action or exact sequence is assigned to this invalid candidate basis.')
        return v
    partial=c['PartialResult']
    v.heading('Ordered partial basis representatives').basis(partial['ExtTau0Basis'],data)
    for label,key,desc in [('Pi0_partial','Pi0','Action of t in partial basis coordinates.'),('T_partial','T','Full normal coordinates of each reduced candidate with a formal scalar; tau acts on that scalar by Frobenius.'),('P_partial','P','Projection returned with the partial basis construction.')]:
        v.matrix(label,desc,partial[key],data)
    v.matrix('W = R o T_partial','Change from partial basis to the true Ext0 coordinates; T_partial=K o W.',c['W'],data)
    v.matrix('W^(-1)','Inverse change L_partial o K, using the additive left inverse of T_partial.',c['WInverse'],data)
    v.table([('Same subgroup',c['SameSubgroupVerified']),('Exact sequence in partial basis coordinates',partial['ExactSequenceVerified'])])
    v.paragraph('Checked by Ore multiplication: common Pi; T_partial=K o W; W o W^(-1)=W^(-1) o W=I; Pi0 o W=W o Pi0_partial. These are different coordinates for the same Ext0.')
    return v


def FJJ2Example(theta, tau, parameters=None):
    """Computed J2 example with upper-right 1+tau in Psi_t; user inputs unchanged."""
    p=parameters if parameters is not None else _symbolic(theta,tau)
    result=FJExampleCompute([[theta+tau**2,1],[0,theta+tau**2]],[[theta,1+tau],[0,theta]],p)
    s=result['Sequence'];comparison=FJPartialBasisComparison(result['Data'],s)
    summary={'FullDimension':s['FullDimension'],'TrueExt0Dimension':s['ExtTau0Dimension'],
        'H':s['H'],'IndividualZeroCount':comparison['PartialZeroCount'],
        'InvalidCountPrediction':s['FullDimension']-comparison['PartialZeroCount'],
        'PartialBasisStatus':comparison['Status'],'QuotientMatrix':s['QuotientMatrix']}
    v=ReportView().heading('J2 counterexample: the individual direction count fails')
    v.paragraph('Phi_t=theta I2+J2+tau^2 I2; Psi_t=theta I2+J2+J2 tau. This separate pair has Ext dimension 8.')
    v.paragraph('For H={{x,y},{z,w}}, HJ2-J2H={{-z,x-w},{0,z}}. The kernel is span{I2,J2}. Only E12 is an individual kernel direction; E11+E22=I2 supplies the missing direction.')
    v.table([(k,val) for k,val in summary.items() if k!='QuotientMatrix'],['Computed quantity','Result'])
    v.paragraph('For every c, [c E11 tau^2]+[c E22 tau^2]-[c^q E12 tau]=0. The candidate parameters (0,c,-c^q,0,0,0,c) represent an inner biderivation and the zero class. The candidate count 7 is not a valid Ext0 dimension.')
    for view in (FJExampleReport(result),FJPartialBasisReport(comparison)):
        v._html.extend(view._html);v._text.extend(view._text)
    return {'Result':result,'Comparison':comparison,'Summary':summary,'Output':v}
