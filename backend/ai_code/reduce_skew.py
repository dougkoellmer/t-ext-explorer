import sympy as sp
from functools import lru_cache
from random import Random

# Input symbols preserve multiplication order.
theta, tau, alpha = sp.symbols(
    "theta tau a", commutative=False
)

# Internal, commuting coefficient-field symbols.
_T, _A = sp.symbols("theta alpha")

@lru_cache(maxsize=32)
def _ore_field_data(p, n):
    if not isinstance(p, (int, sp.Integer)) or not sp.isprime(p):
        raise ValueError("p must be prime.")
    if not isinstance(n, (int, sp.Integer)) or n < 1:
        raise ValueError("n must be a positive integer.")

    p, n = int(p), int(n)
    K = sp.GF(p).frac_field(_T)
    f = sp.Poly(_A, _A, domain=K)

    if n > 1:
        rng = Random(0)
        while True:
            coefficients = (
                [1]
                + [rng.randrange(p) for _ in range(n - 1)]
                + [rng.randrange(1, p)]
            )
            g = sp.Poly.from_list(
                coefficients, _A, modulus=p
            )
            if g.is_irreducible:
                f = sp.Poly(g.as_expr(), _A, domain=K)
                break

    return K, f


def field_modulus(p, n):
    """Return the defining polynomial for alpha."""
    f = _ore_field_data(p, n)[1]
    return f.as_expr().xreplace({_A: alpha})


def reduce_skew(expr, p, n):
    """Apply Frobenius, then reduce all left coefficients."""
    K, f = _ore_field_data(p, n)
    q = int(p)**int(n)
    expr = sp.sympify(expr)

    if expr.atoms(sp.Float):
        raise ValueError("Use exact numbers, not floats.")
    if expr.free_symbols - {theta, tau, alpha}:
        raise ValueError(
            "Recreate expr using the new theta, tau, alpha."
        )

    def C(a):
        return sp.Poly(a, _A, domain=K).rem(f).as_expr()

    def constant(a):
        a = C(a)
        return {0: a} if a != 0 else {}

    def add(P, Q):
        out = dict(P)
        for k, b in Q.items():
            out[k] = C(out.get(k, 0) + b)
        return {k: a for k, a in out.items() if a != 0}

    def multiply(P, Q):
        out = {}
        for i, a in P.items():
            for j, b in Q.items():
                # FIRST: tau^i b = b^(q^i) tau^i.
                # The q-Frobenius fixes alpha.
                twisted_b = b.subs(_T, _T**(q**i))

                # THEN: reduce the left coefficient.
                out[i + j] = C(
                    out.get(i + j, 0) + a*twisted_b
                )

        return {k: a for k, a in out.items() if a != 0}

    def parse(e):
        if e == theta:
            return constant(_T)
        if e == alpha:
            return constant(_A)
        if e == tau:
            return {1: sp.Integer(1)}

        if e.is_Rational:
            if int(e.q) % int(p) == 0:
                raise ZeroDivisionError(
                    "Denominator is divisible by p."
                )
            return constant(e)

        if e.is_Add:
            out = {}
            for a in e.args:
                out = add(out, parse(a))
            return out

        if e.is_Mul:
            out = {0: sp.Integer(1)}
            for a in e.args:
                out = multiply(out, parse(a))
            return out

        if e.is_Pow:
            base, exponent = e.args
            if exponent.is_Integer is not True:
                raise ValueError("Only integer powers are allowed.")

            base = parse(base)
            exponent = int(exponent)

            if exponent < 0:
                if not base:
                    raise ZeroDivisionError(
                        "Denominator vanishes in the field."
                    )
                if set(base) != {0}:
                    raise ValueError(
                        "Only coefficient-field elements can be inverted."
                    )

                inverse = sp.invert(
                    sp.Poly(base[0], _A, domain=K), f
                )
                base = constant(inverse.as_expr())
                exponent = -exponent

            out = {0: sp.Integer(1)}
            while exponent:
                if exponent & 1:
                    out = multiply(out, base)
                base = multiply(base, base)
                exponent //= 2

            return out

        raise ValueError("Unsupported expression: " + str(e))

    coefficients = parse(expr)
    to_input = {_T: theta, _A: alpha}

    return sp.Add(*(
        a.xreplace(to_input)*tau**k
        for k, a in sorted(coefficients.items(), reverse=True)
    ))