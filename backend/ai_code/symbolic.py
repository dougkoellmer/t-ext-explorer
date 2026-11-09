"""Formal q-Frobenius on a rational difference field.

Unspecified characteristic means generic rational formulas with a formal
endomorphism, NOT the power map on a characteristic-zero field. Specialize
the twists to q-power maps in positive characteristic before using a result
over a concrete field. Recorded denominator/pivot conditions must survive.
"""
from dataclasses import dataclass
from functools import lru_cache
import sympy as sp
from .fields import Element, FJError


def _index(r):
    if isinstance(r, bool) or not isinstance(r, (int, sp.Integer)) or r < 0:
        raise FJError("InvalidTwist", "Use a nonnegative integer twist index.")
    return int(r)


class TwistSymbol(sp.Symbol):
    """An atomic coefficient a^(r), distinct from an ordinary power."""
    def __new__(cls, base, r):
        r = _index(r)
        base = sp.sympify(base)
        if isinstance(base, TwistSymbol):
            r, base = r + base.index, base.base
        if not isinstance(base, sp.Symbol) or base.is_commutative is not True:
            raise FJError("InvalidTwistBase", "Twist atoms require a commutative Symbol.")
        if not r:
            return base
        obj = sp.Symbol.__xnew__(cls, f"{base.name}__FJtwist_{r}", commutative=True)
        obj.base, obj.index = base, r
        return obj

    def _hashable_content(self):
        return super()._hashable_content() + (self.base, self.index)

    def _sympystr(self, printer):
        return f"{printer.doprint(self.base)}^({self.index})"

    def _latex(self, printer):
        return printer.doprint(self.base) + r"^{(" + str(self.index) + ")}"

    def __getnewargs_ex__(self):
        return (self.base, self.index), {}


def FJTwist(base, r=1):
    """Return a twist atom; use FJFrobenius for a rational expression."""
    return TwistSymbol(base, r)


@lru_cache(maxsize=32768)
def _normalize(expression, characteristic):
    atoms = tuple(sorted(expression.atoms(sp.Symbol), key=sp.default_sort_key))
    domain = sp.QQ if characteristic is None else sp.GF(characteristic)
    try:
        numerator, denominator = expression.as_numer_denom()
        if atoms:
            # Domain rational functions give exact equality after cancellation.
            field = domain.frac_field(*atoms)
            a, b = field.from_sympy(sp.expand(numerator)), field.from_sympy(sp.expand(denominator))
            if not b:
                raise ZeroDivisionError("Zero denominator")
            value = a / b
            if value:
                lc = value.denom.LC
                value = field.field.new(value.numer.quo_ground(lc),
                                        value.denom.quo_ground(lc))
            return field.to_sympy(value)
        a, b = domain.from_sympy(numerator), domain.from_sympy(denominator)
        if not b:
            raise ZeroDivisionError("Zero denominator")
        return domain.to_sympy(a / b)
    except Exception as exc:
        raise FJError("InvalidSymbolicCoefficient",
                      "Use rational functions of symbols and FJTwist atoms; "
                      "the denominator must be nonzero in the chosen characteristic.") from exc


@dataclass(frozen=True, eq=False)
class SymbolicElement:
    field: object
    value: sp.Expr

    def __bool__(self):
        return self.value != 0

    def __eq__(self, other):
        try:
            return self.field(self.value - self.field(other).value).value == 0
        except (FJError, TypeError):
            return False

    def __add__(self, other):
        return self.field(self.value + self.field(other).value)

    __radd__ = __add__

    def __neg__(self):
        return self.field(-self.value)

    def __sub__(self, other):
        return self + (-self.field(other))

    def __rsub__(self, other):
        return self.field(other) - self

    def __mul__(self, other):
        return self.field(self.value * self.field(other).value)

    __rmul__ = __mul__

    def inverse(self):
        if not self:
            raise ZeroDivisionError("Cannot invert zero in the formal field.")
        self.field.require_nonzero(self.value)
        return self.field(1 / self.value)

    def __truediv__(self, other):
        return self * self.field(other).inverse()

    def __rtruediv__(self, other):
        return self.field(other) / self

    def __pow__(self, exponent):
        if isinstance(exponent, bool) or not isinstance(exponent, (int, sp.Integer)):
            raise FJError("InvalidExponent", "Use an integer exponent.")
        if exponent < 0:
            return self.inverse() ** (-exponent)
        return self.field(self.value ** int(exponent))

    def as_expr(self):
        return self.value

    def _repr_latex_(self):
        return "$" + sp.latex(self.value) + "$"

    def __repr__(self):
        return str(self.value)


class SymbolicField:
    """Rational coefficients with automatically discovered twist generators.

    characteristic=None computes generic formulas over Q. A prime fixes
    coefficient arithmetic modulo that prime while retaining formal twists.
    """
    def __init__(self, characteristic=None, fixed_symbols=(), forbidden_symbols=()):
        if characteristic is not None and (
                isinstance(characteristic, bool) or not isinstance(characteristic, int)
                or not sp.isprime(characteristic)):
            raise FJError("InvalidPrime", "Use None or a prime characteristic.")
        if any(not isinstance(x, sp.Symbol) or isinstance(x, TwistSymbol)
               for x in fixed_symbols):
            raise FJError("InvalidFixedSymbols", "Use ordinary coefficient symbols.")
        self.characteristic = characteristic
        self.fixed_symbols = frozenset(fixed_symbols)
        self.forbidden_symbols = frozenset(forbidden_symbols)
        self._conditions = set()

    def __call__(self, value):
        if isinstance(value, SymbolicElement):
            if value.field is self:
                return value
            value = value.as_expr()
        elif isinstance(value, Element):
            # Do not silently lift a finite-field residue to a generic formula.
            raise FJError("ParameterMismatch", "Specialize explicitly between coefficient modes.")
        if isinstance(value, bool):
            value = int(value)
        expression = sp.sympify(value)
        if expression.has(sp.Float):
            raise FJError("ApproximateInput", "Use exact integers or Rational.")
        if any(x in self.forbidden_symbols or
               (isinstance(x, TwistSymbol) and x.base in self.forbidden_symbols)
               for x in expression.atoms(sp.Symbol)):
            raise FJError("InvalidCoefficient", "tau and q are not coefficient generators.")
        functions = expression.atoms(sp.Function)
        if functions or expression.has(sp.I) or any(
                not power.exp.is_Integer for power in expression.atoms(sp.Pow)):
            raise FJError("InvalidSymbolicCoefficient", "Only rational functions are supported.")
        # Canonicalize user-supplied twists of Frobenius-fixed coefficients.
        expression = expression.xreplace({x: x.base for x in expression.atoms(TwistSymbol)
                                          if x.base in self.fixed_symbols})
        denominator = expression.as_numer_denom()[1]
        normalized = _normalize(expression, self.characteristic)
        self.require_nonzero(denominator)
        return SymbolicElement(self, normalized)

    def require_nonzero(self, expression):
        """Record denominators and inverted pivots before cancellation."""
        numerator, denominator = sp.sympify(expression).as_numer_denom()
        for part in (numerator, denominator):
            if part not in (sp.S.One, -sp.S.One):
                self._conditions.add(part)

    @property
    def zero(self):
        return SymbolicElement(self, sp.S.Zero)

    @property
    def one(self):
        return SymbolicElement(self, sp.S.One)

    def expr(self, value):
        return self(value).as_expr()

    def frobenius(self, value, r=1):
        r = _index(r)
        expression = self(value).as_expr()
        if not r:
            return self(expression)
        # Homomorphism on sums, products, inverses, and integer powers.
        substitutions = {}
        for atom in expression.atoms(sp.Symbol):
            base = atom.base if isinstance(atom, TwistSymbol) else atom
            if base not in self.fixed_symbols:
                substitutions[atom] = FJTwist(atom, r)
        return self(expression.xreplace(substitutions))

    def conditions(self):
        return tuple(sorted(self._conditions, key=sp.default_sort_key))

    def __repr__(self):
        base = "generic Q" if self.characteristic is None else f"F_{self.characteristic}"
        return f"Formal rational difference field ({base})"


class SymbolicParameters:
    """Parameters accepted by every FJ/Ore routine without a numerical q."""
    symbolic_frobenius = True

    def __init__(self, theta=None, tau=None, characteristic=None, fixed_symbols=()):
        theta = sp.Symbol("θ") if theta is None else sp.sympify(theta)
        tau = sp.Symbol("τ") if tau is None else tau
        if not isinstance(tau, sp.Symbol) or isinstance(tau, TwistSymbol):
            raise FJError("InvalidTau", "Use an ordinary unassigned Symbol for tau.")
        self.q = sp.Symbol("q", integer=True, positive=True)
        self.tau = tau
        self.field = SymbolicField(characteristic, fixed_symbols, (tau, self.q))
        self.theta = self.field(theta)

    def frobenius(self, value, k=1):
        return self.field.frobenius(value, k)

    def parameter_monomial(self, c, h):
        # The separate additive coordinate is never treated as Frobenius-fixed.
        return FJTwist(c, h)


FJSymbolicParameters = SymbolicParameters


def SpecializeFrobenius(value, parameters, substitutions=None):
    """Evaluate twist atoms as a^(q^r) using concrete Parameters.

    Substitutions specify only UNTWiSTED generators. Twisted values are
    determined by q; arbitrary independent assignments are not allowed.
    """
    if getattr(parameters, "symbolic_frobenius", False):
        raise FJError("ConcreteParametersRequired", "Use Parameters with numerical q=p^f.")
    if isinstance(value, OreMatrixType()):
        if value.formal:
            raise FJError("FormalEvaluation", "Specialize action matrices, not c-linear families.")
        from .ore import OreMatrix
        return OreMatrix([[[SpecializeFrobenius(x, parameters, substitutions)
                            for x in row] for row in m] for m in value.coefficients],
                         value.rows, value.cols, parameters)
    if isinstance(value, sp.MatrixBase):
        return value.applyfunc(lambda x: parameters.field.expr(
            SpecializeFrobenius(x, parameters, substitutions)))
    expression = value.as_expr() if hasattr(value, "as_expr") else sp.sympify(value)
    substitutions = substitutions or {}
    if any(not isinstance(x, sp.Symbol) or isinstance(x, TwistSymbol)
           for x in substitutions):
        raise FJError("InvalidSubstitution", "Assign untwisted coefficient symbols only.")
    replacements = {}
    for atom in expression.atoms(sp.Symbol):
        if isinstance(atom, TwistSymbol):
            base = substitutions.get(atom.base, atom.base)
            replacements[atom] = parameters.field.expr(parameters.frobenius(base, atom.index))
        elif atom in substitutions:
            replacements[atom] = parameters.field.expr(substitutions[atom])
    return parameters.field(expression.xreplace(replacements))


def OreMatrixType():
    # Deferred import avoids a module cycle during package initialization.
    from .ore import OreMatrix
    return OreMatrix


def FJSymbolicConditions(parameters):
    """Conservative generic specialization conditions; not a proof for all q."""
    field = parameters.field
    if not isinstance(field, SymbolicField):
        raise FJError("SymbolicParametersRequired", "Use FJSymbolicParameters.")
    primes = set()
    for condition in field.conditions():
        atoms = tuple(sorted(condition.atoms(sp.Symbol), key=sp.default_sort_key))
        coefficients = (sp.Poly(condition, *atoms, domain=sp.QQ).coeffs()
                        if atoms else [sp.Rational(condition)])
        integers = []
        for coefficient in coefficients:
            primes.update(sp.factorint(int(coefficient.q)))
            integers.append(abs(int(coefficient.p)))
        content = sp.igcd(*integers) if len(integers) > 1 else (integers[0] if integers else 0)
        if content:
            primes.update(sp.factorint(int(content)))
    return {"Characteristic": field.characteristic,
            "GenericMode": field.characteristic is None,
            "NonzeroConditions": list(field.conditions()),
            "ExcludedPrimes": sorted(primes),
            "ExcludedQRule": "Exclude all powers of an excluded prime, and any q/parameter "
                             "specialization at which a recorded condition vanishes. "
                             "Recompute the compatible basis and report after specialization."}


def FJSpecializationReport(symbolic_parameters, concrete_parameters, substitutions=None):
    """Check the recorded conditions for ONE concrete field/q/parameter choice."""
    conditions = FJSymbolicConditions(symbolic_parameters)
    failures = []
    symbolic_field = symbolic_parameters.field
    if (symbolic_field.characteristic is not None and
            symbolic_field.characteristic != concrete_parameters.field.characteristic):
        failures.append({"Condition": "characteristic", "Reason": "Characteristic mismatch"})
    for symbol in symbolic_field.fixed_symbols:
        try:
            value = (substitutions or {}).get(symbol, symbol)
            if concrete_parameters.frobenius(value) != concrete_parameters.field(value):
                failures.append({"Condition": symbol, "Reason": "Not Frobenius-fixed"})
        except (FJError, ZeroDivisionError) as exc:
            failures.append({"Condition": symbol, "Reason": str(exc)})
    for condition in conditions["NonzeroConditions"]:
        try:
            if not SpecializeFrobenius(condition, concrete_parameters, substitutions):
                failures.append({"Condition": condition, "Reason": "Specializes to zero"})
        except (FJError, ZeroDivisionError) as exc:
            failures.append({"Condition": condition, "Reason": str(exc)})
    return {"Q": concrete_parameters.q,
            "Characteristic": concrete_parameters.field.characteristic,
            "RecordedConditionsSatisfied": not failures,
            "FailedConditions": failures,
            "RecomputeReportAfterSpecialization": True}
