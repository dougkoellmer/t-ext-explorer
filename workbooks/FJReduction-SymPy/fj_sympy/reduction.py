"""Algorithms 7--13 and the checked extension exact sequence.

The original Phi and Psi are immutable OreMatrix objects. The jump-basis
inverse is used only to construct inner biderivations, never to replace
Phi with a conjugate or an isomorphic t-module.
"""
from .fields import FJError
from .linear import (zero, identity, unit, transpose, flatten, add, scale,
                     multiply, rank, rref, inverse)
from .compatible_basis import compatible_basis_report
from .ore import OreMatrix, linear_scale, linear_add, linear_columns_to_ore


def _nilpotent(matrix, f):
    power = identity(f, len(matrix))
    for _ in range(len(matrix)):
        power = multiply(power, matrix, f)
    return not any(flatten(power))


def _linear_map(matrix, vector, f):
    """Apply an L-linear map to Frobenius-index coefficient tuples."""
    out = []
    for row in matrix:
        entry = ()
        for coefficient, value in zip(row, vector):
            entry = linear_add(entry, linear_scale(value, coefficient), f)
        out.append(entry)
    return out


class ReductionData:
    """Prepare and certify the fixed input used by Algorithms 7--13."""
    def __init__(self, phi, psi, parameters):
        self.parameters, self.field = parameters, parameters.field
        p, f = parameters, self.field
        self.phi = OreMatrix.from_matrix(phi, p)
        self.psi = OreMatrix.from_matrix(psi, p)
        self.d, self.e = self.phi.rows, self.psi.rows
        if (not self.d or not self.e or self.phi.cols != self.d
                or self.psi.cols != self.e or self.phi.degree < 1
                or self.phi.formal or self.psi.formal):
            raise FJError("InvalidModuleShape", "Use square ordinary Phi and Psi, with deg(Phi)>0.")
        # Step 1: check the nilpotent parts of the given constant matrices.
        self.nphi = add(self.phi.coefficient(0), scale(identity(f, self.d), -p.theta))
        self.npsi = add(self.psi.coefficient(0), scale(identity(f, self.e), -p.theta))
        if not _nilpotent(self.nphi, f) or not _nilpotent(self.npsi, f):
            raise FJError("InvalidNilpotentPart", "The constant part must be theta I + nilpotent.")
        # Step 2 (Algorithm 7): find the tagged jump basis without changing Phi.
        self.basis_report = compatible_basis_report(
            [self.phi.coefficient(k) for k in range(1, self.phi.degree+1)], f)
        if not self.basis_report["Found"]:
            raise FJError("NoCompatibleBasis", "The positive blocks have no compatible jump basis.")
        self.cb = self.basis_report["BasisMatrix"]
        self.cb_inverse = inverse(self.cb, f)
        self.sigma = [s-1 for s in self.basis_report["Sigma"]]
        self.jumps = self.basis_report["BlockChoices"][:]
        self.minimum_jump = min(self.jumps)
        if self.psi.degree >= self.minimum_jump:
            raise FJError("DegreeConditionFailed", "Require deg_tau(Psi) < min(j_s).")
        self.thresholds = [self.jumps[:] for _ in range(self.e)]
        # Step 3: independently verify actual successive rank jumps.
        columns, jumps = [], []
        for k in range(self.phi.degree, 0, -1):
            matrix = self.phi.coefficient(k)
            for s in self.sigma:
                vector = tuple(row[s] for row in matrix)
                if rank(transpose(columns+[vector])) > len(columns):
                    columns.append(vector)
                    jumps.append((s+1, k, vector))
        if sorted(jumps) != self.basis_report["Basis"]:
            raise FJError("FlagVerificationFailed", "The selected columns fail the rank-jump definition.")
        # Step 4: prepare U_(i,s)=E_(i,s) C^-1 and verify each leading pivot.
        self.units = [multiply(unit(f, self.e, self.d, i, s), self.cb_inverse, f)
                      for i in range(self.e) for s in range(self.d)]
        for i in range(self.e):
            for s in range(self.d):
                u = OreMatrix.monomial(self.units[i*self.d+s], 0, p)
                g = self.inner(u)
                if g.degree != self.jumps[s] or g.coefficient(g.degree)[i][s] != f.one:
                    raise FJError("PivotVerificationFailed", "A basic inner pivot is not monic.")
        self._analysis, self._full, self._zero = None, None, None

    def inner(self, u):
        """Compute U Phi - Psi U with the Ore multiplication rule."""
        u = OreMatrix.from_matrix(u, self.parameters)
        return u @ self.phi - self.psi @ u

    def basis(self, zero_part=False, one_based=True):
        """Order coordinates by target row, ORIGINAL column, then tau degree."""
        shifts = self.constant_analysis()["Epsilon"] if zero_part else None
        triples = [(i, s, r) for i in range(self.e) for s in range(self.d)
                   for r in range(1 if zero_part else 0,
                                  self.jumps[s]+(shifts[i][s] if zero_part else 0))]
        return [(i+1, s+1, r) for i, s, r in triples] if one_based else triples

    def generator(self, triple, formal=False, value=None):
        """Internal zero-based triple -> c E_(i,s) tau^r (or a numeric input)."""
        i, s, r = triple
        f = self.field
        z = () if formal else f.zero
        matrix = [[z for _ in range(self.d)] for _ in range(self.e)]
        matrix[i][s] = (f.one,) if formal else f(1 if value is None else value)
        return OreMatrix.monomial(matrix, r, self.parameters, formal)

    def coordinates(self, a, basis=None):
        basis = self.basis(one_based=False) if basis is None else basis
        return [a.coefficient(r)[i][s] for i, s, r in basis]

    def representative(self, coordinates):
        """Build a full normal representative from an ordinary coordinate vector."""
        basis = self.basis(one_based=False)
        if len(coordinates) != len(basis):
            raise FJError("VectorShape", "The full coordinate vector has the wrong length.")
        out = OreMatrix.zero(self.e, self.d, self.parameters)
        for t, x in zip(basis, coordinates):
            out = out + self.generator(t, value=x)
        return out

    def _reduce(self, v, thresholds, record=False):
        """Algorithm 8; Algorithm 10 uses the same loop with shifted bounds."""
        if v.shape != (self.e, self.d):
            raise FJError("RepresentativeShape", "Use dim(Psi) rows and dim(Phi) columns.")
        p, f, a = self.parameters, self.field, v
        witness = OreMatrix.zero(self.e, self.d, p, v.formal)
        steps, count = [], 0
        bound = self.e * self.d * max(0, v.degree+1)
        while True:
            # Step 1: locate the largest degree violating a column bound.
            # Ties follow sigma, then increasing target row.
            positions = [(a.entry_degree(i, s), order, i, s)
                         for order, s in enumerate(self.sigma) for i in range(self.e)
                         if a.entry_degree(i, s) >= thresholds[i][s]]
            if not positions:
                break
            r, order, i, s = min(positions, key=lambda t: (-t[0], t[1], t[2]))
            shift = r - self.jumps[s]
            mu = a.coefficient(r)[i][s]
            # Step 2: the twist of C^-1 is necessary because tau^shift moves
            # past the source coefficient matrices.
            u = self.units[i*self.d+s]
            matrix = [[linear_scale(mu, p.frobenius(x, shift)) if v.formal
                       else mu * p.frobenius(x, shift) for x in row] for row in u]
            h = OreMatrix.monomial(matrix, shift, p, v.formal)
            g = self.inner(h)
            # Step 3: certify cancellation and no reintroduced earlier pivot.
            leading = g.coefficient(r)
            if (g.degree > r or leading[i][s] != mu or any(
                    leading[ii][ss] for ii in range(self.e) for ss in self.sigma[:order])):
                raise FJError("TriangularPivotFailed", "The cancellation is not triangular.")
            # Step 4: subtract the inner biderivation and accumulate its witness.
            a, witness = a-g, witness+h
            count += 1
            if count > bound:
                raise FJError("ReductionProgressFailed", "The certified reduction bound was exceeded.")
            if record:
                steps.append({"Degree": r, "Row": i+1, "Column": s+1,
                              "Shift": shift, "Coefficient": mu})
        return a, witness, steps

    def reduce_report(self, v):
        """Algorithm 8 with a step-by-step trace and exact inner witness."""
        v = OreMatrix.from_matrix(v, self.parameters)
        reduced, witness, steps = self._reduce(v, self.thresholds, record=True)
        verified = v-reduced == self.inner(witness)
        if not verified:
            raise FJError("WitnessVerificationFailed", "V-V_red is not inner(U).")
        return {"ReducedMatrix": reduced, "InnerWitness": witness,
                "Steps": steps, "VerifiedDifference": verified}

    def reduce(self, v):
        v = OreMatrix.from_matrix(v, self.parameters)
        return self._reduce(v, self.thresholds)[0]

    def _extension(self, basis, thresholds):
        """Build action columns with a FORMAL c, not by setting c=1."""
        columns = []
        for triple in basis:
            # Step 1: apply Psi_t to c times one coordinate generator.
            v = self.psi @ self.generator(triple, formal=True)
            # Step 2: reduce and read sum b_h c^(q^h) in every coordinate.
            reduced = self._reduce(v, thresholds)[0]
            columns.append(self.coordinates(reduced, basis))
        # Step 3: turn each c^(q^h) into tau^h in the action matrix.
        return linear_columns_to_ore(columns, len(basis), len(basis), self.parameters)

    def extension(self):
        """Algorithm 9: the action Pi_t in the full normal family E."""
        if self._full is None:
            self._full = self._extension(self.basis(one_based=False), self.thresholds)
        return self._full

    def constant_analysis(self):
        """Constant map, epsilon, actual kernel, quotient coordinates/action."""
        if self._analysis is not None:
            return self._analysis
        f, size = self.field, self.e*self.d
        # Step 1: evaluate U NPhi - NPsi U on all basic constant directions.
        columns = [flatten(add(multiply(u, self.nphi, f),
                               scale(multiply(self.npsi, u, f), -f.one)))
                   for u in self.units]
        constant_map = transpose(columns)
        reduced, pivots = rref(constant_map)
        image_rank, h = len(pivots), size-len(pivots)
        zero_columns = [j for j, col in enumerate(columns) if not any(col)]
        epsilon = [[int(i*self.d+s not in zero_columns) for s in range(self.d)]
                   for i in range(self.e)]
        # Step 2: complete an image basis to a basis of the constant space.
        complete = [columns[j] for j in pivots]
        complement = []
        for j in range(size):
            vector = [f.one if l == j else f.zero for l in range(size)]
            if rank(transpose(complete+[vector])) > len(complete):
                complete.append(vector)
                complement.append(vector)
        binv = inverse(transpose(complete), f)
        projection, image_coordinates = binv[image_rank:], binv[:image_rank]
        section = transpose(complement) if h else zero(f, size, 0)
        # Step 3: compute the quotient action, instead of assuming it is scalar.
        images = [flatten(multiply(self.npsi,
                  [v[i*self.d:(i+1)*self.d] for i in range(self.e)], f))
                  for v in complement]
        nilpotent = multiply(projection, transpose(images), f, cols=h) if h else []
        quotient = add(scale(identity(f, h), self.parameters.theta), nilpotent)
        # Step 4: record the full kernel, including dependent nonzero directions.
        kernel_vectors = []
        for j in range(size):
            if j in pivots:
                continue
            vector = [f.zero] * size
            vector[j] = f.one
            for row, pivot in enumerate(pivots):
                vector[pivot] = -reduced[row][j]
            kernel_vectors.append(vector)
        kernel_matrices = []
        for vector in kernel_vectors:
            matrix = zero(f, self.e, self.d)
            for u, c in zip(self.units, vector):
                matrix = add(matrix, scale(u, c))
            kernel_matrices.append(matrix)
        self._analysis = {
            "ConstantTermMap": constant_map, "ConstantMapRank": image_rank,
            "ImageIndices": pivots, "ImageUnits": [self.units[j] for j in pivots],
            "ImageBasis": transpose([columns[j] for j in pivots]) if pivots else zero(f, size, 0),
            "ImageCoordinates": image_coordinates, "QuotientProjection": projection,
            "QuotientSection": section, "QuotientMatrix": quotient,
            "QuotientNilpotentPart": nilpotent, "KernelDimension": h,
            "KernelCoefficientVectors": kernel_vectors, "KernelMatrices": kernel_matrices,
            "Epsilon": epsilon, "S": epsilon, "H": len(zero_columns),
            "PseudocodeZeroCount": len(zero_columns),
            "CoordinateKernelCondition": len(zero_columns) == h,
            "ScalarAdditiveQuotient": not any(flatten(nilpotent))}
        return self._analysis

    def partial_ga_part(self):
        """Algorithm 13: literal zero-column count h and epsilon_(i,s)."""
        a = self.constant_analysis()
        return a["H"], [row[:] for row in a["Epsilon"]]

    def _require_article(self, scalar=False):
        a = self.constant_analysis()
        if not a["CoordinateKernelCondition"]:
            raise FJError("ArticleCoordinateConditionFailed",
                          "The nonzero constant directions are dependent. "
                          "Use general_report/general_reduce_zero_report explicitly.")
        if scalar and not a["ScalarAdditiveQuotient"]:
            raise FJError("ArticleScalarQuotientConditionFailed",
                          "The quotient action is not theta times the identity.")
        return a

    def _zero_thresholds(self):
        epsilon = self.constant_analysis()["Epsilon"]
        return [[self.jumps[s]+epsilon[i][s] for s in range(self.d)] for i in range(self.e)]

    def partial_transition(self):
        """Algorithm 12: reduce each E0 generator into the full family E."""
        self._require_article()
        basis0, columns = self.basis(True, False), []
        for triple in basis0:
            g = self.generator(triple, formal=True)
            columns.append(self.coordinates(self.reduce(g)))
        return linear_columns_to_ore(columns, len(self.basis()), len(basis0), self.parameters)

    def partial_ext_tau0(self):
        """Algorithm 11: Pi0_t in 1 <= r < j_s+epsilon_(i,s)."""
        self._require_article()
        if self._zero is None:
            self._zero = self._extension(self.basis(True, False), self._zero_thresholds())
        return self._zero

    def partial_reduce_zero_report(self, v):
        """Algorithm 10: zero-constant normal form, trace and inner witness."""
        self._require_article()
        v = OreMatrix.from_matrix(v, self.parameters)
        if any(flatten(v.coefficient(0))):
            raise FJError("NonzeroConstantRepresentative", "Ext_tau,0 input must have zero constant term.")
        reduced, witness, steps = self._reduce(v, self._zero_thresholds(), record=True)
        verified = not any(flatten(reduced.coefficient(0))) and v-reduced == self.inner(witness)
        if not verified:
            raise FJError("ExtZeroWitnessFailed", "The zero-constant witness check failed.")
        return {"ReducedMatrix": reduced, "FullReducedRepresentative": self.reduce(reduced),
                "Coordinates": self.coordinates(reduced, self.basis(True, False)),
                "InnerWitness": witness, "Steps": steps, "VerifiedDifference": True,
                "ReductionMode": "ManuscriptCoordinateKernel"}

    def partial_reduce_zero(self, v):
        return self.partial_reduce_zero_report(v)["ReducedMatrix"]

    def _projection(self):
        a, basis = self.constant_analysis(), self.basis(one_based=False)
        matrix = [[a["QuotientProjection"][j][i*self.d+s] if r == 0 else self.field.zero
                   for i, s, r in basis] for j in range(a["KernelDimension"])]
        return OreMatrix.monomial(matrix, 0, self.parameters,
                                  shape=(a["KernelDimension"], len(basis)))

    def _article_inverse_maps(self):
        """Construct additive L and S; S is not asserted to be t-linear."""
        a, f, p = self.constant_analysis(), self.field, self.parameters
        full_basis, zero_basis = self.basis(one_based=False), self.basis(True, False)
        columns = []
        for triple in full_basis:
            g = self.generator(triple, formal=True)
            if triple[2] == 0:
                # Step 1: split the constant term into image and complement.
                x = flatten(g.coefficient(0))
                y = _linear_map(a["ImageCoordinates"], x, f)
                z = _linear_map(a["QuotientProjection"], x, f)
                h = [[() for _ in range(self.d)] for _ in range(self.e)]
                for u, coefficient in zip(a["ImageUnits"], y):
                    for i in range(self.e):
                        for s in range(self.d):
                            h[i][s] = linear_add(h[i][s], linear_scale(coefficient, u[i][s]), f)
                c = _linear_map(a["QuotientSection"], z, f)
                c = [c[i*self.d:(i+1)*self.d] for i in range(self.e)]
                # Step 2: subtract the image's inner lift and the complement.
                v = g-self.inner(OreMatrix.monomial(h, 0, p, True))-OreMatrix.monomial(c, 0, p, True)
                if any(flatten(v.coefficient(0))):
                    raise FJError("ConstantLiftVerificationFailed", "The lift did not remove the constant term.")
                g = self._reduce(v, self._zero_thresholds())[0]
            columns.append(self.coordinates(g, zero_basis))
        left = linear_columns_to_ore(columns, len(zero_basis), len(full_basis), p)
        section_matrix = [[a["QuotientSection"][i*self.d+s][j] if r == 0 else f.zero
                           for j in range(a["KernelDimension"])] for i, s, r in full_basis]
        section = OreMatrix.monomial(section_matrix, 0, p,
                                    shape=(len(full_basis), a["KernelDimension"]))
        return left, section

    def partial_exact_sequence(self):
        """Article coordinates with six independently checked polynomial identities."""
        a = self._require_article(scalar=True)
        full, pi0, transition, projection = self.extension(), self.partial_ext_tau0(), self.partial_transition(), self._projection()
        left, section = self._article_inverse_maps()
        p, f, h = self.parameters, self.field, a["H"]
        quotient = OreMatrix.monomial(scale(identity(f, h), p.theta), 0, p, shape=(h, h))
        id_full = OreMatrix.monomial(identity(f, full.rows), 0, p, shape=full.shape)
        id_zero = OreMatrix.monomial(identity(f, pi0.rows), 0, p, shape=pi0.shape)
        id_quotient = OreMatrix.monomial(identity(f, h), 0, p, shape=(h, h))
        checks = {
            "InclusionVerified": full @ transition == transition @ pi0,
            "QuotientVerified": projection @ full == quotient @ projection,
            "CompositionZeroVerified": (projection @ transition).is_zero,
            "InjectivityVerified": left @ transition == id_zero,
            "SurjectivityVerified": projection @ section == id_quotient,
            "KernelImageVerified": transition @ left + section @ projection == id_full}
        if not all(checks.values()) or pi0.rows+h != full.rows:
            raise FJError("ExactSequenceVerificationFailed", "A morphism or additive inverse identity failed.")
        return {**self._report_base(), **checks, "Pi": full, "Pi0": pi0,
                "T": transition, "P": projection, "Transition": transition,
                "QuotientProjection": projection, "QuotientMatrix": quotient,
                "InclusionLeftInverse": left, "AdditiveSection": section,
                "FullDimension": full.rows, "ExtTau0Dimension": pi0.rows,
                "QuotientDimension": h, "GaDimension": h, "H": h,
                "ExtTau0Basis": self.basis(True), "ExtTau0Mode": "ManuscriptCoordinateKernel",
                "ExactSequenceVerified": True}

    def _report_base(self):
        a = self.constant_analysis()
        return {"Version": "0.2.0", "Data": self, "Parameters": self.parameters,
                "BasisMatrix": self.cb, "BasisInverse": self.cb_inverse,
                "CompatibleBasis": self.basis_report["Basis"],
                "Sigma": [s+1 for s in self.sigma], "JumpDegrees": self.jumps[:],
                "MinimumJumpDegree": self.minimum_jump, "N": self.thresholds,
                "Epsilon": a["Epsilon"], "S": a["Epsilon"], "FullBasis": self.basis(),
                "PseudocodeGaCount": a["H"], "ConstantAnalysis": a,
                "ScalarAdditiveQuotient": a["ScalarAdditiveQuotient"],
                "CoordinateKernelCondition": a["CoordinateKernelCondition"],
                "CanonicalFlagVerified": True}

    def _general_zero_structure(self):
        """Explicit fallback: constant image directions plus positive monomials."""
        a, f, p = self.constant_analysis(), self.field, self.parameters
        basis = self.basis(one_based=False)
        columns, descriptions = [], []
        for j, index in enumerate(a["ImageIndices"]):
            columns.append([a["ImageBasis"][i*self.d+s][j] if r == 0 else f.zero
                            for i, s, r in basis])
            descriptions.append({"Kind": "ConstantImage",
                                 "InnerPosition": (index//self.d+1, index%self.d+1)})
        for j, triple in enumerate(basis):
            if triple[2] > 0:
                columns.append([f.one if l == j else f.zero for l in range(len(basis))])
                descriptions.append({"Kind": "Monomial",
                                     "Position": (triple[0]+1, triple[1]+1, triple[2])})
        dim = len(columns)
        k = transpose(columns) if dim else zero(f, len(basis), 0)
        transition = OreMatrix.monomial(k, 0, p, shape=(len(basis), dim))
        if not dim:
            return transition, OreMatrix.zero(0, 0, p), descriptions, [], []
        rows = rref(transpose(k))[1]
        inv = inverse([k[i] for i in rows], f)
        action = self.extension() @ transition
        pi0 = OreMatrix([multiply(inv, [m[i] for i in rows], f) for m in action.coefficients],
                        dim, dim, p)
        if action != transition @ pi0:
            raise FJError("ExtZeroStabilityFailed", "The general kernel is not invariant.")
        return transition, pi0, descriptions, rows, inv

    def general_report(self):
        """General constant-term construction, including non-scalar quotients."""
        a, p = self.constant_analysis(), self.parameters
        if a["CoordinateKernelCondition"]:
            transition, pi0, descriptions = self.partial_transition(), self.partial_ext_tau0(), self.basis(True)
            mode = "ManuscriptCoordinateKernel"
        else:
            transition, pi0, descriptions, _, _ = self._general_zero_structure()
            mode = "GeneralConstantTermKernel"
        full, projection, h = self.extension(), self._projection(), a["KernelDimension"]
        quotient = OreMatrix.monomial(a["QuotientMatrix"], 0, p, shape=(h, h))
        inclusion_ok = full @ transition == transition @ pi0
        quotient_ok = projection @ full == quotient @ projection
        composition_ok = (projection @ transition).is_zero
        if not (inclusion_ok and quotient_ok and composition_ok):
            raise FJError("ExactSequenceVerificationFailed", "A general kernel/quotient identity failed.")
        return {**self._report_base(), "Pi": full, "Pi0": pi0, "T": transition,
                "P": projection, "Transition": transition, "QuotientProjection": projection,
                "QuotientMatrix": quotient, "FullDimension": full.rows,
                "ExtTau0Dimension": pi0.rows, "QuotientDimension": h,
                "GaDimension": h if a["ScalarAdditiveQuotient"] else None,
                "ExtTau0Basis": descriptions, "ExtTau0Mode": mode,
                "InclusionVerified": inclusion_ok, "QuotientVerified": quotient_ok,
                "CompositionZeroVerified": composition_ok}

    def general_reduce_zero_report(self, v):
        """General fallback with an actual zero-constant representative/witness."""
        a = self.constant_analysis()
        if a["CoordinateKernelCondition"]:
            return self.partial_reduce_zero_report(v)
        p, f = self.parameters, self.field
        v = OreMatrix.from_matrix(v, p)
        if any(flatten(v.coefficient(0))):
            raise FJError("NonzeroConstantRepresentative", "Ext_tau,0 input must have zero constant term.")
        full_report = self.reduce_report(v)
        reduced = full_report["ReducedMatrix"]
        fullcoords = self.coordinates(reduced)
        transition, _, descriptions, rows, inv = self._general_zero_structure()
        coordinates = [row[0] for row in multiply(inv, [[fullcoords[i]] for i in rows], f, cols=1)]
        reconstructed = [row[0] for row in multiply(transition.coefficient(0),
                           [[x] for x in coordinates], f, cols=1)]
        if reconstructed != fullcoords:
            raise FJError("ExtZeroCoordinatesFailed", "The class is outside the constant-term kernel.")
        # Lift the image coordinates back to a constant inner biderivation.
        h = zero(f, self.e, self.d)
        for u, c in zip(a["ImageUnits"], coordinates[:a["ConstantMapRank"]]):
            h = add(h, scale(u, c))
        h = OreMatrix.monomial(h, 0, p)
        zero_representative = reduced-self.inner(h)
        witness = full_report["InnerWitness"]+h
        if any(flatten(zero_representative.coefficient(0))) or v-zero_representative != self.inner(witness):
            raise FJError("ExtZeroWitnessFailed", "The general zero-constant witness failed.")
        return {"ReducedMatrix": zero_representative, "FullReducedRepresentative": reduced,
                "Coordinates": coordinates, "InnerWitness": witness,
                "VerifiedDifference": True, "Steps": full_report["Steps"],
                "ReductionMode": "GeneralConstantTermKernel"}

    # ===== BEGIN: computation of the true Ext_0 and exact sequence =====
    # Delegacja do fj_sympy/sequence.py. Historyczne metody partial_* są powyżej.
    def constant_data(self):
        from .sequence import constant_data
        return constant_data(self)

    def sequence_coordinates(self):
        from .sequence import sequence_coordinates
        return sequence_coordinates(self)

    def sequence_actions(self, pi):
        from .sequence import sequence_actions
        return sequence_actions(pi, self)

    def exact_sequence(self):
        from .sequence import exact_sequence
        return exact_sequence(self)

    def ga_part(self):
        a = self.constant_data()
        return a["H"], [row[:] for row in a["Epsilon"]]

    def transition(self):
        return self.sequence_coordinates()["K"]

    def ext_tau0(self):
        return self.sequence_actions(self.extension())["Pi0"]

    def reduce_zero_report(self, v):
        from .sequence import reduce_zero_report
        return reduce_zero_report(v, self)

    def reduce_zero(self, v):
        return self.reduce_zero_report(v)["ReducedMatrix"]

    # ===== END: computation of the true Ext_0 and exact sequence =====


def construct_fj_basis(phi, parameters):
    """Algorithm 7 without needing a target."""
    ph = OreMatrix.from_matrix(phi, parameters)
    if ph.rows < 1 or ph.cols != ph.rows or ph.degree < 1:
        raise FJError("InvalidPhi", "Use a nonempty square Phi of positive degree.")
    report = compatible_basis_report([ph.coefficient(k) for k in range(1, ph.degree+1)],
                                    parameters.field)
    if not report["Found"]:
        raise FJError("NoCompatibleBasis", "No compatible jump basis exists.")
    return report["BasisMatrix"], report["Sigma"]
