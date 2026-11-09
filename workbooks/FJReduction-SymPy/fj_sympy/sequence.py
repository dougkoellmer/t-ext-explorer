"""New Algorithms 10--13 from main.tex (2026-10-12).
Column vectorization; unchanged full normal-coordinate order.
"""
from .fields import FJError
from .linear import zero, identity, unit, transpose, flatten, add, scale, multiply, rref, rank, inverse
from .ore import OreMatrix, linear_add, linear_scale


# ===== BEGIN: computation of the true Ext_0 and exact sequence =====
# Article pseudocode implementation; full Ext comes from data.extension().

def _vec(matrix):
    return flatten(transpose(matrix))


def _unvec(vector, e, d):
    return [[vector[s*e+i] for s in range(d)] for i in range(e)]


# ===== BEGIN: true Ext_0 construction (Algorithms 10-12) =====

def constant_data(data):
    """Algorithm 10: D, h, B, Z, U, P0, over the chosen exact field."""
    cached = getattr(data, '_constant_data', None)
    if cached is not None:
        return cached
    f, e, d = data.field, data.e, data.d
    n = e*d
    units = [unit(f, e, d, i, s) for s in range(d) for i in range(e)]
    columns = [_vec(add(multiply(u, data.nphi, f), scale(multiply(data.npsi, u, f), -f.one))) for u in units]
    D = transpose(columns)
    pivots = rref(D)[1]
    r, h = len(pivots), n-len(pivots)
    B = [[D[i][j] for j in pivots] for i in range(n)]
    pivot_rows = rref(transpose(B))[1] if r else []
    complement = [j for j in range(n) if j not in pivot_rows]
    Z = [[f.one if i == j else f.zero for j in complement] for i in range(n)]
    G = inverse([b+z for b, z in zip(B, Z)], f)
    U, P0 = G[:r], G[r:]
    images = [_vec(multiply(data.npsi, _unvec([Z[i][j] for i in range(n)], e, d), f)) for j in range(h)]
    nil = multiply(P0, transpose(images), f, cols=h)
    Q = add(scale(identity(f, h), data.parameters.theta), nil)
    # Former coordinate count is retained only as a diagnostic.
    oldcols = [flatten(add(multiply(u, data.nphi, f), scale(multiply(data.npsi, u, f), -f.one))) for u in data.units]
    epsilon = [[int(any(oldcols[i*d+s])) for s in range(d)] for i in range(e)]
    count = sum(not any(c) for c in oldcols)
    result = {'D': D, 'ConstantTermMap': D, 'Vectorization': 'ColumnMajor',
              'ConstantMapRank': r, 'H': h, 'KernelDimension': h, 'QuotientDimension': h,
              'B': B, 'Z': Z, 'U': U, 'P0': P0, 'ImageBasis': B,
              'ImageCoordinates': U, 'QuotientProjection': P0, 'QuotientSection': Z,
              'ImageIndices': pivots, 'ImageUnits': [units[j] for j in pivots],
              'QuotientNilpotentPart': nil, 'QuotientMatrix': Q,
              'Epsilon': epsilon, 'S': epsilon, 'PseudocodeZeroCount': count,
              'CoordinateKernelCondition': count == h,
              'ScalarAdditiveQuotient': not any(flatten(nil))}
    data._constant_data = result
    return result


def sequence_coordinates(data):
    """Algorithm 11: bases and constant K, R, P in original full coordinates."""
    cached = getattr(data, '_sequence_coordinates', None)
    if cached is not None:
        return cached
    a, f, p = constant_data(data), data.field, data.parameters
    full = data.basis(one_based=False)
    e, d, m = data.e, data.d, len(full)
    r, h = a['ConstantMapRank'], a['H']
    m0 = m-h
    constants = {(i, s): j for j, (i, s, k) in enumerate(full) if k == 0}
    Ic = [constants[i, s] for s in range(d) for i in range(e)]
    positive = [j for j, t in enumerate(full) if t[2] > 0]
    K, R, P, section = zero(f, m, m0), zero(f, m0, m), zero(f, h, m), zero(f, m, h)
    for j, ix in enumerate(Ic):
        for k in range(r):
            K[ix][k], R[k][ix] = a['B'][j][k], a['U'][k][j]
        for k in range(h):
            P[k][ix], section[ix][k] = a['P0'][k][j], a['Z'][j][k]
    for j, ix in enumerate(positive):
        K[ix][r+j] = R[r+j][ix] = f.one
    basis0 = [{'Kind': 'ConstantImage', 'Matrix': _unvec([a['B'][i][j] for i in range(e*d)], e, d),
               'InnerUnit': a['ImageUnits'][j]} for j in range(r)]
    basis0 += [{'Kind': 'Monomial', 'Position': (full[j][0]+1, full[j][1]+1, full[j][2])} for j in positive]
    basisq = [_unvec([a['Z'][i][j] for i in range(e*d)], e, d) for j in range(h)]
    const = lambda mat, rows, cols: OreMatrix.monomial(mat, 0, p, shape=(rows, cols))
    result = {'FullBasis': data.basis(), 'ExtTau0Basis': basis0, 'QuotientBasis': basisq,
              'K': const(K, m, m0), 'R': const(R, m0, m), 'P': const(P, h, m),
              'AdditiveSection': const(section, m, h), 'FullDimension': m,
              'ExtTau0Dimension': m0, 'ConstantIndices': [i+1 for i in Ic], 'PositiveIndices': [i+1 for i in positive]}
    data._sequence_coordinates = result
    return result


def sequence_actions(pi, data):
    """Algorithm 12: R C_a K^(a), with successive one-step twists of K."""
    a, c, f, p = constant_data(data), sequence_coordinates(data), data.field, data.parameters
    pi = OreMatrix.from_matrix(pi, p)
    m, m0, h = c['FullDimension'], c['ExtTau0Dimension'], a['H']
    if pi.shape != (m, m) or pi.formal:
        raise FJError('InvalidActionShape', 'Pi must be an ordinary Ore action on the full normal basis.')
    twk, R = c['K'].coefficient(0), c['R'].coefficient(0)
    coefficients, images = [], []
    degree = max(0, pi.degree)
    for k in range(degree+1):
        image = multiply(pi.coefficient(k), twk, f, cols=m0)
        images.append(image)
        coefficients.append(multiply(R, image, f, cols=m0))
        if k < degree:
            twk = [[p.frobenius(x, 1) for x in row] for row in twk]
    return {'Pi0': OreMatrix(coefficients, m0, m0, p),
            'QuotientMatrix': OreMatrix.monomial(a['QuotientMatrix'], 0, p, shape=(h, h)),
            'InclusionAction': OreMatrix(images, m, m0, p)}


# ===== END: true Ext_0 construction (Algorithms 10-12) =====


def _power(matrix, n, f):
    result, base = identity(f, len(matrix)), matrix
    while n:
        if n & 1:
            result = multiply(result, base, f, cols=len(matrix))
        n //= 2
        if n:
            base = multiply(base, base, f, cols=len(matrix))
    return result


# ===== BEGIN: assembly and verification of the exact sequence (Algorithm 13) =====

def exact_sequence(data):
    """Algorithm 13: assembly, ranks, Ore identities and Anderson checks."""
    a, c, f, p = constant_data(data), sequence_coordinates(data), data.field, data.parameters
    full = data.extension()
    actions = sequence_actions(full, data)
    pi0, Q = actions['Pi0'], actions['QuotientMatrix']
    K, R, P, section = c['K'], c['R'], c['P'], c['AdditiveSection']
    m, m0, h = c['FullDimension'], c['ExtTau0Dimension'], a['H']
    ident = lambda n: OreMatrix.monomial(identity(f, n), 0, p, shape=(n, n))
    rankK, rankP = rank(K.coefficient(0)), rank(P.coefficient(0))
    nil0 = add(pi0.coefficient(0), scale(identity(f, m0), -p.theta))
    checks = {'InclusionVerified': actions['InclusionAction'] == K @ pi0,
              'QuotientVerified': P @ full == Q @ P,
              'CompositionZeroVerified': (P @ K).is_zero,
              'InjectivityVerified': R @ K == ident(m0) and rankK == m0,
              'SurjectivityVerified': P @ section == ident(h) and rankP == h,
              'KernelImageVerified': (P @ K).is_zero and rankK == m0 and rankP == h
                   and m0 == m-h and K @ R + section @ P == ident(m),
              'AndersonExt0Verified': not any(flatten(_power(nil0, data.e, f))),
              'AndersonQuotientVerified': not any(flatten(_power(a['QuotientNilpotentPart'], data.e, f)))}
    if not all(checks.values()):
        raise FJError('ExactSequenceVerificationFailed', 'An Algorithm 13 rank, Ore, exactness or Anderson check failed.')
    return {**c, **checks, 'Version': '0.4.5', 'Data': data, 'Parameters': p,
            'BasisMatrix': data.cb, 'BasisInverse': data.cb_inverse, 'CompatibleBasis': data.basis_report['Basis'],
            'Sigma': [s+1 for s in data.sigma], 'JumpDegrees': data.jumps[:],
            'MinimumJumpDegree': data.minimum_jump, 'N': data.thresholds,
            'Pi': full, 'Pi0': pi0, 'QuotientMatrix': Q, 'T': K, 'Transition': K,
            'QuotientProjection': P, 'InclusionLeftInverse': R, 'H': h, 'QuotientDimension': h,
            'GaDimension': h if a['ScalarAdditiveQuotient'] else None,
            'Epsilon': a['Epsilon'], 'S': a['Epsilon'], 'PseudocodeGaCount': a['PseudocodeZeroCount'],
            'ConstantAnalysis': a, 'ScalarAdditiveQuotient': a['ScalarAdditiveQuotient'],
            'CoordinateKernelCondition': a['CoordinateKernelCondition'], 'InclusionRank': rankK, 'ProjectionRank': rankP,
            'D': a['D'], 'B': a['B'], 'Z': a['Z'], 'U': a['U'], 'P0': a['P0'],
            'ExtTau0Mode': 'ImageDPlusPositiveNormalMonomials', 'CanonicalFlagVerified': True,
            'ExactSequenceVerified': True}


# ===== END: assembly and verification of the exact sequence (Algorithm 13) =====


def _apply(matrix, vector, data, formal=False):
    if formal:
        out = []
        for row in matrix:
            value = ()
            for coefficient, x in zip(row, vector):
                value = linear_add(value, linear_scale(x, coefficient), data.field)
            out.append(value)
        return out
    return [row[0] for row in multiply(matrix, [[x] for x in vector], data.field, cols=1)]


# ===== BEGIN: true Ext_0 coordinates and zero-constant witness =====

def reduce_zero_report(v, data):
    """Full reduction, P x = 0, Ext0 coordinates R x, and a constant lift."""
    a, c, f, p = constant_data(data), sequence_coordinates(data), data.field, data.parameters
    v = OreMatrix.from_matrix(v, p)
    result = data.reduce_report(v)
    reduced, x = result['ReducedMatrix'], data.coordinates(result['ReducedMatrix'])
    quotient = _apply(c['P'].coefficient(0), x, data, v.formal)
    if any(quotient):
        raise FJError('ClassNotInExtZero', 'The full reduced class has nonzero quotient coordinates P x.')
    coords = _apply(c['R'].coefficient(0), x, data, v.formal)
    if _apply(c['K'].coefficient(0), coords, data, v.formal) != x:
        raise FJError('ExtZeroCoordinatesFailed', 'K (R x) did not recover the full class coordinates.')
    imagecoords = _apply(a['U'], _vec(reduced.coefficient(0)), data, v.formal)
    H = [[() if v.formal else f.zero for _ in range(data.d)] for _ in range(data.e)]
    for u, value in zip(a['ImageUnits'], imagecoords):
        for i in range(data.e):
            for s in range(data.d):
                H[i][s] = linear_add(H[i][s], linear_scale(value, u[i][s]), f) if v.formal else H[i][s]+value*u[i][s]
    ho = OreMatrix.monomial(H, 0, p, v.formal, shape=(data.e, data.d))
    zero_rep = reduced-data.inner(ho)
    witness = result['InnerWitness']+ho
    if any(flatten(zero_rep.coefficient(0))) or v-zero_rep != data.inner(witness):
        raise FJError('ExtZeroWitnessFailed', 'The zero-constant witness identity failed.')
    return {**result, 'FullReducedRepresentative': reduced, 'Coordinates': coords,
            'QuotientCoordinates': quotient, 'ZeroConstantRepresentative': zero_rep,
            'ConstantLift': ho, 'ZeroConstantInnerWitness': witness,
            'ZeroConstantWitnessVerified': True, 'ReductionMode': 'FullReductionThenImageDCoordinates'}

# ===== END: true Ext_0 coordinates and zero-constant witness =====
# ===== END: computation of the true Ext_0 and exact sequence =====
