"""Matrices over L{tau}, and sparse Frobenius-linear formal parameters.

An ordinary entry is sum a_k tau^k. A formal coefficient is a tuple
(b_0,...,b_h) meaning sum b_h c^(q^h). We never expand c^(q^h)
as a dense polynomial. This is essential when q or h is large.
"""
from dataclasses import dataclass
import sympy as sp
from .fields import FJError


def linear_trim(values):
    values = tuple(values)
    while values and not values[-1]:
        values = values[:-1]
    return values


def linear_add(a, b, f):
    return linear_trim([(a[h] if h < len(a) else f.zero)
                        + (b[h] if h < len(b) else f.zero)
                        for h in range(max(len(a), len(b)))])


def linear_scale(a, scalar):
    return linear_trim([v * scalar for v in a])


def linear_frobenius(a, k, p):
    return () if not a else (p.field.zero,) * k + tuple(
        p.frobenius(v, k) for v in a)


@dataclass(frozen=True, init=False)
class OreMatrix:
    """Immutable coefficient data. Use @ for skew multiplication."""
    coefficients: tuple
    rows: int
    cols: int
    parameters: object
    formal: bool

    def __init__(self, coefficients, rows, cols, parameters, formal=False):
        if rows < 0 or cols < 0:
            raise FJError("MatrixShape", "Dimensions cannot be negative.")
        f = parameters.field
        values = []
        for matrix in coefficients:
            if len(matrix) != rows or any(len(row) != cols for row in matrix):
                raise FJError("MatrixShape", "Coefficient matrix has the wrong shape.")
            values.append(tuple(tuple(linear_trim(f(v) for v in x) if formal else f(x)
                                      for x in row) for row in matrix))
        while values and not any(x for row in values[-1] for x in row):
            values.pop()
        object.__setattr__(self, "coefficients", tuple(values))
        object.__setattr__(self, "rows", rows)
        object.__setattr__(self, "cols", cols)
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "formal", formal)

    @classmethod
    def from_matrix(cls, matrix, parameters):
        """Parse displayed polynomials; tau here is a display placeholder."""
        if isinstance(matrix, cls):
            if matrix.parameters != parameters:
                raise FJError("ParameterMismatch", "Use the same coefficient field and q.")
            return matrix
        if isinstance(matrix, sp.MatrixBase):
            rows, cols = matrix.shape
            matrix = matrix.tolist()
        else:
            rows, cols = len(matrix), (len(matrix[0]) if matrix else 0)
            if any(len(row) != cols for row in matrix):
                raise FJError("MatrixShape", "The displayed matrix must be rectangular.")
        polys = []
        try:
            for row in matrix:
                polys.append([sp.Poly(sp.expand(sp.sympify(v)), parameters.tau,
                                     domain="EX") for v in row])
        except (sp.PolynomialError, TypeError) as exc:
            raise FJError("InvalidSkewEntry", "Every entry must be polynomial in tau.") from exc
        degree = max([int(v.degree()) for row in polys for v in row if not v.is_zero] + [-1])
        coefficients = [[[parameters.field(polys[i][j].nth(k)) for j in range(cols)]
                         for i in range(rows)] for k in range(degree+1)]
        return cls(coefficients, rows, cols, parameters)

    @classmethod
    def zero(cls, rows, cols, parameters, formal=False):
        return cls([], rows, cols, parameters, formal)

    @classmethod
    def monomial(cls, matrix, degree, parameters, formal=False, shape=None):
        if not isinstance(degree, int) or degree < 0:
            raise FJError("InvalidDegree", "Use a nonnegative integer tau degree.")
        rows, cols = shape or (len(matrix), len(matrix[0]) if matrix else 0)
        z = () if formal else parameters.field.zero
        zero = [[z for _ in range(cols)] for _ in range(rows)]
        return cls([zero for _ in range(degree)] + [matrix], rows, cols, parameters, formal)

    @property
    def shape(self):
        return self.rows, self.cols

    @property
    def degree(self):
        return len(self.coefficients) - 1

    @property
    def is_zero(self):
        return not self.coefficients

    def coefficient(self, k):
        """Return a fresh matrix, so callers cannot alter stored coefficients."""
        if 0 <= k < len(self.coefficients):
            return [list(row) for row in self.coefficients[k]]
        z = () if self.formal else self.parameters.field.zero
        return [[z for _ in range(self.cols)] for _ in range(self.rows)]

    def entry_degree(self, i, s):
        return next((k for k in range(self.degree, -1, -1)
                     if self.coefficients[k][i][s]), -1)

    def _compatible(self, b):
        if self.parameters != b.parameters:
            raise FJError("ParameterMismatch", "The skew matrices use different parameters.")

    def __add__(self, b):
        self._compatible(b)
        if self.shape != b.shape or self.formal != b.formal:
            raise FJError("MatrixShape", "Addition requires equal shapes and coefficient modes.")
        f, values = self.parameters.field, []
        for k in range(max(self.degree, b.degree)+1):
            a_k, b_k = self.coefficient(k), b.coefficient(k)
            values.append([[linear_add(a_k[i][j], b_k[i][j], f) if self.formal
                            else a_k[i][j] + b_k[i][j]
                            for j in range(self.cols)] for i in range(self.rows)])
        return OreMatrix(values, self.rows, self.cols, self.parameters, self.formal)

    def scale(self, scalar):
        scalar = self.parameters.field(scalar)
        return OreMatrix([[[linear_scale(v, scalar) if self.formal else v * scalar
                            for v in row] for row in matrix]
                          for matrix in self.coefficients],
                         self.rows, self.cols, self.parameters, self.formal)

    def __neg__(self):
        return self.scale(-1)

    def __sub__(self, b):
        return self + (-b)

    def __matmul__(self, b):
        """Skew product: (a tau^k)(b tau^h)=a b^(q^k) tau^(k+h)."""
        self._compatible(b)
        if self.cols != b.rows:
            raise FJError("MatrixShape", "Inner dimensions of a skew product do not agree.")
        if self.formal and b.formal:
            raise FJError("NonlinearParameterProduct", "Two formal parameters cannot be multiplied.")
        p, f, formal = self.parameters, self.parameters.field, self.formal or b.formal
        if self.is_zero or b.is_zero:
            return OreMatrix.zero(self.rows, b.cols, p, formal)
        z = () if formal else f.zero
        values = [[[z for _ in range(b.cols)] for _ in range(self.rows)]
                  for _ in range(self.degree + b.degree + 1)]
        for k, a_k in enumerate(self.coefficients):
            for h, b_h in enumerate(b.coefficients):
                for i in range(self.rows):
                    for l in range(self.cols):
                        x = a_k[i][l]
                        if not x:
                            continue
                        for j in range(b.cols):
                            y = b_h[l][j]
                            if not y:
                                continue
                            # Step 1: twist the RIGHT coefficient, including c.
                            y = linear_frobenius(y, k, p) if b.formal else p.frobenius(y, k)
                            # Step 2: multiply by the left coefficient.
                            if self.formal:
                                product = linear_scale(x, y)
                            elif b.formal:
                                product = linear_scale(y, x)
                            else:
                                product = x * y
                            # Step 3: accumulate at tau degree k+h.
                            old = values[k+h][i][j]
                            values[k+h][i][j] = linear_add(old, product, f) if formal else old + product
        return OreMatrix(values, self.rows, b.cols, p, formal)

    def display(self, parameter=None):
        """SymPy display matrix; multiply the OreMatrix objects, not this view."""
        f, tau = self.parameters.field, self.parameters.tau
        c = parameter if parameter is not None else sp.Symbol("c")
        def expression(value):
            if not self.formal:
                return f.expr(value)
            monomial = getattr(self.parameters, "parameter_monomial",
                               lambda c, h: c ** (self.parameters.q ** h))
            return sum((f.expr(v) * monomial(c, h)
                        for h, v in enumerate(value)), sp.S.Zero)
        return sp.Matrix(self.rows, self.cols, lambda i, j: sum(
            (expression(m[i][j]) * tau**k for k, m in enumerate(self.coefficients)),
            sp.S.Zero))

    def _repr_latex_(self):
        return self.display()._repr_latex_()

    def __repr__(self):
        return f"OreMatrix{self.shape}, degree={self.degree}\n{self.display()}"

    def evaluate(self, vector):
        """Evaluate an additive-polynomial matrix at a vector over L."""
        if self.formal:
            raise FJError("FormalEvaluation", "Use evaluate_formal for a formal parameter.")
        if len(vector) != self.cols:
            raise FJError("VectorShape", "The coordinate vector has the wrong length.")
        p, f = self.parameters, self.parameters.field
        x = [f(v) for v in vector]
        return [sum((m[i][j] * p.frobenius(x[j], k)
                     for k, m in enumerate(self.coefficients)
                     for j in range(self.cols)), f.zero) for i in range(self.rows)]

    def evaluate_formal(self, value):
        if not self.formal:
            raise FJError("FormalEvaluation", "This matrix has no formal parameter.")
        p, f = self.parameters, self.parameters.field
        x = f(value)
        return OreMatrix([[[sum((v * p.frobenius(x, h) for h, v in enumerate(entry)),
                               f.zero) for entry in row] for row in m]
                          for m in self.coefficients], self.rows, self.cols, p)


def linear_columns_to_ore(columns, rows, cols, parameters):
    """Interpret c^(q^h) coefficients as tau^h in an additive map."""
    degree = max([len(entry)-1 for column in columns for entry in column] + [-1])
    f = parameters.field
    coefficients = [[[columns[j][i][h] if h < len(columns[j][i]) else f.zero
                     for j in range(cols)] for i in range(rows)]
                    for h in range(degree+1)]
    return OreMatrix(coefficients, rows, cols, parameters)


def skew_multiply(a, b, parameters):
    return OreMatrix.from_matrix(a, parameters) @ OreMatrix.from_matrix(b, parameters)


def skew_equal(a, b, parameters):
    return OreMatrix.from_matrix(a, parameters) == OreMatrix.from_matrix(b, parameters)
