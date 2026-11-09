"""Dense linear algebra over exact SymPy domain elements.

Elimination keeps characteristic and extension relations fixed.
No generic Matrix.rank() is allowed to guess the field.
"""
from .fields import FJError


def zero(f, rows, cols):
    return [[f.zero for _ in range(cols)] for _ in range(rows)]


def identity(f, n):
    return [[f.one if i == j else f.zero for j in range(n)] for i in range(n)]


def transpose(a):
    return [list(row) for row in zip(*a)] if a else []


def flatten(a):
    return [x for row in a for x in row]


def add(a, b):
    return [[x + y for x, y in zip(u, v)] for u, v in zip(a, b)]


def scale(a, x):
    return [[v * x for v in row] for row in a]


def multiply(a, b, f, cols=None):
    cols = len(b[0]) if cols is None and b else (cols or 0)
    return [[sum((a[i][k] * b[k][j] for k in range(len(b))), f.zero)
             for j in range(cols)] for i in range(len(a))]


def rref(a):
    """Return a new reduced matrix and zero-based pivot-column indices."""
    m = [list(row) for row in a]
    pivots, row = [], 0
    for col in range(len(m[0]) if m else 0):
        # Step 1: locate a nonzero pivot using exact field equality.
        pivot = next((i for i in range(row, len(m)) if m[i][col]), None)
        if pivot is None:
            continue
        # Step 2: normalize using inversion in the same field.
        m[row], m[pivot] = m[pivot], m[row]
        divisor = m[row][col]
        m[row] = [v / divisor for v in m[row]]
        # Step 3: eliminate this column in every other row.
        for i in range(len(m)):
            if i != row and m[i][col]:
                coefficient = m[i][col]
                m[i] = [x - coefficient * y for x, y in zip(m[i], m[row])]
        pivots.append(col)
        row += 1
        if row == len(m):
            break
    return m, pivots


def rank(a):
    return len(rref(a)[1])


def inverse(a, f):
    n = len(a)
    if any(len(row) != n for row in a):
        raise FJError("MatrixShape", "Inversion requires a square matrix.")
    if not n:
        return []
    reduced, pivots = rref([u + v for u, v in zip(a, identity(f, n))])
    if pivots[:n] != list(range(n)):
        raise FJError("SingularMatrix", "The matrix is singular over the chosen field.")
    return [row[n:] for row in reduced]


def unit(f, e, d, i, s):
    return [[f.one if (a, b) == (i, s) else f.zero
             for b in range(d)] for a in range(e)]


def normalize_matrix(a, f):
    if hasattr(a, "tolist"):
        a = a.tolist()
    if not isinstance(a, (list, tuple)) or not a:
        raise FJError("MatrixShape", "Use a nonempty rectangular matrix.")
    try:
        cols = len(a[0])
        if cols < 1 or any(len(row) != cols for row in a):
            raise FJError("MatrixShape", "Use a nonempty rectangular matrix.")
        return [[f(v) for v in row] for row in a]
    except TypeError as exc:
        raise FJError("MatrixShape", "Use rows of coefficients.") from exc
