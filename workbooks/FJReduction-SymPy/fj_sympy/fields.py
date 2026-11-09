"""Exact coefficient fields backed by SymPy domains.

All arithmetic stays in the chosen characteristic. SymPy expressions are
used only at the input/display boundary, never to infer a matrix's field.
"""
from dataclasses import dataclass
from functools import lru_cache
import random
import sympy as sp
from sympy.polys.agca.extensions import FiniteExtension
from sympy.polys.polyclasses import DMP


def _normalize(domain, value):
    """Canonicalize units in rational denominators, also inside extensions."""
    if isinstance(domain, FiniteExtension):
        coefficients = [_normalize(domain.domain, c) for c in value.rep.to_list()]
        return domain.new(DMP(coefficients, domain.domain))
    if getattr(domain, "is_FractionField", False) and value:
        # SymPy 1.14 can leave 2/2 distinct from 1 over GF(3)(theta).
        # A monic denominator removes this unit ambiguity exactly.
        numerator, denominator = value.numer, value.denom
        leading = denominator.LC
        numerator = numerator.quo_ground(leading)
        denominator = denominator.quo_ground(leading)
        return domain.field.new(numerator, denominator)
    return value


def _poly_trim(a):
    a = list(a)
    while a and not a[-1]:
        a.pop()
    return a


def _poly_sub(a, b, ground):
    return _poly_trim([(a[i] if i < len(a) else ground.zero)
                       - (b[i] if i < len(b) else ground.zero)
                       for i in range(max(len(a), len(b)))])


def _poly_mul(a, b, ground):
    out = [ground.zero for _ in range(max(0, len(a)+len(b)-1))]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] = out[i+j] + x*y
    return _poly_trim(out)


def _poly_divmod(a, b, ground):
    a, b = _poly_trim(a), _poly_trim(b)
    if not b:
        raise ZeroDivisionError("Polynomial division by zero.")
    out = [ground.zero for _ in range(max(0, len(a)-len(b)+1))]
    while a and len(a) >= len(b):
        shift, coefficient = len(a)-len(b), a[-1]/b[-1]
        out[shift] = out[shift] + coefficient
        a = _poly_sub(a, [ground.zero]*shift+[coefficient*x for x in b], ground)
    return _poly_trim(out), a


@dataclass(frozen=True, eq=False)
class Element:
    """Exact domain element with canonical normalization after each operation."""
    field: object
    value: object

    def __post_init__(self):
        object.__setattr__(self, "value", _normalize(self.field.domain, self.value))

    def __bool__(self):
        return bool(self.value)

    def __eq__(self, other):
        try:
            other = self.field(other)
            return self.value == other.value
        except (FJError, TypeError):
            return False

    def __add__(self, other):
        return Element(self.field, self.value + self.field(other).value)

    __radd__ = __add__

    def __neg__(self):
        return Element(self.field, -self.value)

    def __sub__(self, other):
        return self + (-self.field(other))

    def __rsub__(self, other):
        return self.field(other) - self

    def __mul__(self, other):
        return Element(self.field, self.value * self.field(other).value)

    __rmul__ = __mul__

    def inverse(self):
        if not self:
            raise ZeroDivisionError("Cannot invert zero in a field.")
        domain = self.field.domain
        if not isinstance(domain, FiniteExtension) or not self.field.variables:
            return Element(self.field, domain.one / self.value)
        # Extended Euclid over normalized ground-field coefficients. This
        # bypasses SymPy's finite-extension inversion over rational functions.
        ground = Field(domain.domain, self.field.characteristic)
        modulus = [Element(ground, v) for v in reversed(domain.mod.to_list())]
        value = [Element(ground, v) for v in reversed(self.value.rep.to_list())]
        r0, r1, s0, s1 = modulus, value, [], [ground.one]
        while r1:
            quotient, remainder = _poly_divmod(r0, r1, ground)
            r0, r1 = r1, remainder
            s0, s1 = s1, _poly_sub(s0, _poly_mul(quotient, s1, ground), ground)
        if len(r0) != 1:
            raise FJError("NonInvertible", "The extension modulus does not define a field.")
        _, result = _poly_divmod([v/r0[0] for v in s0], modulus, ground)
        return Element(self.field, domain.new(DMP(
            [v.value for v in reversed(result)], domain.domain)))

    def __truediv__(self, other):
        return self * self.field(other).inverse()

    def __rtruediv__(self, other):
        return self.field(other) / self

    def __pow__(self, exponent):
        if not isinstance(exponent, int):
            raise FJError("InvalidExponent", "Use an integer field exponent.")
        if exponent < 0:
            return self.inverse() ** (-exponent)
        # Binary powering normalizes after every multiplication.
        result, base = self.field.one, self
        while exponent:
            if exponent & 1:
                result = result * base
            exponent //= 2
            if exponent:
                base = base * base
        return result

    def as_expr(self):
        return self.field.domain.to_sympy(self.value)

    def __repr__(self):
        return str(self.as_expr())


class FJError(ValueError):
    """A mathematical/input condition failed; code identifies it."""
    def __init__(self, code, message):
        self.code = code
        super().__init__(f"{code}: {message}")


@dataclass(frozen=True)
class Field:
    domain: object
    characteristic: int
    constant_degree: int = 1
    extension_modulus: object = None
    variables: tuple = ()

    def __call__(self, value):
        """Coerce exactly; reject floats and denominators zero in the field."""
        if isinstance(value, Element):
            if value.field == self:
                return value
            value = value.as_expr()
        if self.domain.of_type(value):
            try:
                return Element(self, self.domain.convert(value))
            except Exception as exc:
                raise FJError("InvalidCoefficient", str(exc)) from exc
        if isinstance(value, bool):
            value = int(value)
        if hasattr(value, "as_expr"):
            value = value.as_expr()
        try:
            expression = sp.sympify(value)
            if expression.has(sp.Float):
                raise FJError("ApproximateInput", "Use exact integers or Rational.")
            numerator, denominator = expression.as_numer_denom()
            a = Element(self, self.domain.from_sympy(sp.expand(numerator)))
            b = Element(self, self.domain.from_sympy(sp.expand(denominator)))
            if not b:
                raise FJError("ZeroDenominator", "The denominator is zero in L.")
            return a / b
        except FJError:
            raise
        except Exception as exc:
            raise FJError("InvalidCoefficient",
                          f"{value!s} cannot be interpreted in {self.domain}.") from exc

    @property
    def zero(self):
        return Element(self, self.domain.zero)

    @property
    def one(self):
        return Element(self, self.domain.one)

    def expr(self, value):
        """Return a readable SymPy expression representing an exact element."""
        return self(value).as_expr()

    def __repr__(self):
        return str(self.domain)


def rational_field():
    """Q is useful for Algorithms 1--6 only."""
    return Field(sp.QQ, 0)


def prime_field(p):
    if not isinstance(p, int) or isinstance(p, bool) or not sp.isprime(p):
        raise FJError("InvalidPrime", "p must be prime; GF(9) is not F_9.")
    return Field(sp.GF(p), p)


def finite_field(p, modulus, generator=None):
    """Construct F_(p^r) = F_p[a]/(modulus), checking irreducibility."""
    base = prime_field(p)
    if isinstance(modulus, sp.Poly):
        generator = generator or modulus.gens[0]
        expression = modulus.as_expr()
    else:
        expression = sp.sympify(modulus)
        symbols = tuple(expression.free_symbols)
        if generator is None:
            if len(symbols) != 1:
                raise FJError("InvalidModulus", "Specify the extension generator.")
            generator = symbols[0]
    poly = sp.Poly(expression, generator, domain=base.domain).monic()
    if poly.degree() < 1 or not poly.is_irreducible:
        raise FJError("ReducibleModulus", "Use a nonconstant irreducible polynomial.")
    return Field(FiniteExtension(poly), p, poly.degree(), poly)


@lru_cache(maxsize=128, typed=True)
def finite_field_pn(p, n, generator=None):
    """Build F_(p**n), choosing a reproducible irreducible modulus automatically.

    p is a prime and n is a positive integer. The same p, n and generator
    use the same polynomial across calls and runs; no global RNG is changed.
    For n=1 the result is the prime field, with no algebraic generator.
    """
    base = prime_field(p)
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise FJError("InvalidExtensionDegree", "n must be a positive integer.")
    if generator is not None and not isinstance(generator, sp.Symbol):
        raise FJError("InvalidGenerator", "Use a SymPy Symbol for alpha.")
    if n == 1:
        return base
    generator = generator if generator is not None else sp.Symbol("alpha")
    rng = random.Random(f"FJReduction-finite-field:{p}:{n}")
    while True:
        coefficients = [1] + [rng.randrange(p) for _ in range(n-1)] + [rng.randrange(1, p)]
        modulus = sp.Poly.from_list(coefficients, generator, domain=base.domain)
        if modulus.is_irreducible:
            return finite_field(p, modulus, generator)


def FJParametersPN(p, n, theta=None, tau=None, alpha=None, coefficient_symbols=()):
    """Automatically construct F_(p**n)(theta, ...) and its q-power Frobenius."""
    theta = theta if theta is not None else sp.Symbol("theta")
    tau = tau if tau is not None else sp.Symbol("tau")
    alpha = alpha if alpha is not None else sp.Symbol("alpha")
    constants = finite_field_pn(p, n, alpha)
    field = function_field(constants, theta, *coefficient_symbols)
    return Parameters(field, q=p**n, theta=theta, tau=tau)


def function_field(base, *variables):
    """Construct rational functions, including F_(p^r)(theta).

    For an extension base, build F_p(theta)[a]/(modulus). This avoids
    SymPy's incomplete coercions for FiniteExtension.frac_field(theta).
    """
    if not variables or any(not isinstance(v, sp.Symbol) for v in variables):
        raise FJError("InvalidVariables", "Supply distinct SymPy symbols.")
    all_variables = base.variables + tuple(variables)
    if len(set(all_variables)) != len(all_variables):
        raise FJError("InvalidVariables", "Function-field variables must be distinct.")
    if base.extension_modulus is None:
        root = sp.QQ if base.characteristic == 0 else sp.GF(base.characteristic)
        domain = root.frac_field(*all_variables)
    else:
        poly = base.extension_modulus
        if poly.gens[0] in all_variables:
            raise FJError("InvalidVariables", "Separate algebraic and transcendental generators.")
        root = sp.GF(base.characteristic).frac_field(*all_variables)
        domain = FiniteExtension(sp.Poly(poly.as_expr(), poly.gens[0], domain=root))
    return Field(domain, base.characteristic, base.constant_degree,
                 base.extension_modulus, all_variables)


@dataclass(frozen=True)
class Parameters:
    field: Field
    q: int
    theta: object
    tau: sp.Symbol = sp.Symbol("τ")

    def __post_init__(self):
        # Step 1: check that q is a positive power of char(L).
        p, remainder, exponent = self.field.characteristic, self.q, 0
        if p < 2 or not isinstance(remainder, int) or isinstance(remainder, bool) or remainder < 2:
            raise FJError("InvalidFrobenius", "Use positive characteristic and q=p^f.")
        while remainder % p == 0:
            exponent += 1
            remainder //= p
        if remainder != 1:
            raise FJError("InvalidFrobenius", "q must be a power of char(L).")
        # Step 2: F_(p^f) must embed into the specified constant field.
        if self.field.constant_degree % exponent:
            raise FJError("MissingConstantField", "The coefficient field must contain F_q.")
        if not isinstance(self.tau, sp.Symbol):
            raise FJError("InvalidTau", "tau must be an unassigned SymPy Symbol.")
        if self.tau in self.field.variables or (
                self.field.extension_modulus is not None
                and self.tau == self.field.extension_modulus.gens[0]):
            raise FJError("InvalidTau", "tau cannot also be a coefficient-field generator.")
        object.__setattr__(self, "theta", self.field(self.theta))

    def frobenius(self, value, k=1):
        if not isinstance(k, int) or k < 0:
            raise FJError("NegativeTwist", "Frobenius indices must be nonnegative integers.")
        return self.field(value) ** (self.q ** k)
