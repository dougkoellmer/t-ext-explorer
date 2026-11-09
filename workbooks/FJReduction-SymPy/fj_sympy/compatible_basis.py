"""Algorithms 1--6, preserving the Mathematica search and tie-breaking.

Public source-column and block tags are ONE-based, as in the manuscript.
Internal Python offsets are zero-based. Input matrices are never mutated.
"""
from .fields import FJError, rational_field
from .linear import normalize_matrix, transpose, rank, inverse, multiply


def preprocess(a, field=None):
    """Algorithm 1: full rank, maximal suffix, flag and block quotas."""
    f = field or rational_field()
    # Step 1: accept either [A1,...,An] or [A1 | ... | An].
    if isinstance(a, (list, tuple)) and a and isinstance(a[0], (list, tuple)) and a[0] and isinstance(a[0][0], (list, tuple)):
        blocks = [normalize_matrix(b, f) for b in a]
        d = len(blocks[0])
        if any(len(b) != d or any(len(row) != d for row in b) for b in blocks):
            raise FJError("BlockShape", "Every coefficient block must be d by d.")
        matrix = [sum((b[i] for b in blocks), []) for i in range(d)]
    else:
        matrix = normalize_matrix(a, f)
        d = len(matrix)
        if len(matrix[0]) % d:
            raise FJError("BlockShape", "The number of columns must be a multiple of d.")
        blocks = [[row[j:j+d] for row in matrix] for j in range(0, len(matrix[0]), d)]
    n = len(blocks)
    if rank(matrix) < d:
        return False, None, [], []
    # Step 2: read suffix ranks from right to left, stopping at full rank.
    required, ranks, previous = [0] * n, [None] * n, 0
    for j in range(n - 1, -1, -1):
        suffix = [sum((b[i] for b in blocks[j:]), []) for i in range(d)]
        current = rank(suffix)
        ranks[j] = current
        required[j] = current - previous
        if current == d:
            # Step 3: retain ORIGINAL block numbers; removed quotas stay zero.
            trim = {"Blocks": blocks[j:], "Dimension": d, "FirstBlock": j+1,
                    "OriginalBlockCount": n, "BlockIndices": list(range(j+1, n+1)),
                    "Matrix": suffix, "SuffixRanks": ranks, "Field": f}
            return True, trim, [(r+1, ranks[r]) for r in range(j, n)], required
        previous = current
    raise FJError("InternalPreprocessError", "No full-rank suffix was found.")


def build_column_lists(trim):
    """Algorithm 3: nonzero candidates (original position, block, vector)."""
    d, k, n = trim["Dimension"], trim["FirstBlock"], trim["OriginalBlockCount"]
    lists = [[] for _ in range(d)]
    for i in range(d):
        for j in range(n, k-1, -1):
            column = tuple(row[i] for row in trim["Blocks"][j-k])
            if any(column):
                lists[i].append((i+1, j, column))
    return lists


def order_column_lists(lists):
    """Algorithm 4: stable increasing-cardinality sort."""
    return sorted(lists, key=len)


def _valid_candidate(candidate, trim):
    try:
        i, j, vector = candidate
        d, k, n = trim["Dimension"], trim["FirstBlock"], trim["OriginalBlockCount"]
        return (isinstance(i, int) and isinstance(j, int) and 1 <= i <= d
                and k <= j <= n and tuple(vector) == tuple(
                    row[i-1] for row in trim["Blocks"][j-k]))
    except (TypeError, ValueError):
        return False


def compatible_ordering(trim, basis):
    """Algorithm 2: precedence graph, then lexicographic topological order."""
    d, f = trim["Dimension"], trim["Field"]
    basis = sorted(basis, key=lambda t: t[0])
    if (len(basis) != d or [t[0] for t in basis] != list(range(1, d+1))
            or not all(_valid_candidate(t, trim) for t in basis)):
        raise FJError("InvalidBasis", "Choose one tagged source column per original position.")
    cb = transpose([t[2] for t in basis])
    if rank(cb) != d:
        raise FJError("DependentBasis", "The selected vectors must be independent.")
    inv, choices, edges = inverse(cb, f), [t[1] for t in basis], set()
    # Step 1: express each input block in the selected basis.
    # C^-1 is a computational coordinate tool: Phi is NOT conjugated.
    for j, block in zip(trim["BlockIndices"], trim["Blocks"]):
        selected = {s for s in range(d) if choices[s] == j}
        coordinates = multiply(inv, block, f)
        # Step 2: a nonzero coefficient in selected row l forces l before i.
        for i in set(range(d)) - selected:
            for l in selected:
                if coordinates[l][i]:
                    edges.add((l, i))
    # Step 3: deterministic topological sort; a cycle means incompatibility.
    remaining, sigma = set(range(d)), []
    while remaining:
        available = [v for v in remaining if not any(
            end == v and start in remaining for start, end in edges)]
        if not available:
            return False, []
        v = min(available)
        sigma.append(v+1)
        remaining.remove(v)
    return True, sigma


def dfs_selection(trim, ordered_lists, required):
    """Algorithm 5: quota and independence pruning, with full backtracking."""
    d, n = trim["Dimension"], trim["OriginalBlockCount"]
    if (len(required) != n or any(not isinstance(v, int) or v < 0 for v in required)
            or sum(required) != d or any(required[:trim["FirstBlock"]-1])):
        raise FJError("InvalidRequired", "Quotas must sum to d and vanish on removed blocks.")
    if len(ordered_lists) != d:
        raise FJError("InvalidColumnLists", "Supply exactly d column lists.")
    if any(not candidates for candidates in ordered_lists):
        return False, [], []
    if (not all(_valid_candidate(t, trim) for ls in ordered_lists for t in ls)
            or any(len({t[0] for t in ls}) != 1 for ls in ordered_lists)
            or sorted(ls[0][0] for ls in ordered_lists) != list(range(1, d+1))):
        raise FJError("InvalidColumnLists", "Lists must cover each original position once.")
    used, partial = [0] * n, []

    def search(level):
        if level == d:
            # Sort a COPY; partial remains indexed by search depth.
            basis = sorted(partial, key=lambda t: t[0])
            found, sigma = compatible_ordering(trim, basis)
            return (basis, sigma) if found else None
        for candidate in ordered_lists[level]:
            _, block, vector = candidate
            if used[block-1] >= required[block-1]:
                continue
            if rank([t[2] for t in partial] + [vector]) != level + 1:
                continue
            partial.append(candidate)
            used[block-1] += 1
            result = search(level+1)
            if result is not None:
                return result
            # Undo BOTH state updates before exploring the next candidate.
            used[block-1] -= 1
            partial.pop()
        return None

    answer = search(0)
    return (True, answer[0], answer[1]) if answer is not None else (False, [], [])


def compatible_basis_report(a, field=None):
    """Algorithm 6: orchestrate Algorithms 1--5 and expose their outputs."""
    found, trim, flag, required = preprocess(a, field)
    if not found:
        return {"Found": False, "Basis": [], "Sigma": [], "BasisMatrix": [],
                "BlockChoices": [], "Required": [], "Flag": [], "TrimmedData": None,
                "Status": "RankDeficient"}
    lists = build_column_lists(trim)
    found, basis, sigma = dfs_selection(trim, order_column_lists(lists), required)
    return {"Found": found, "Basis": basis, "Sigma": sigma,
            "BasisMatrix": transpose([t[2] for t in basis]) if found else [],
            "BlockChoices": [t[1] for t in basis] if found else [],
            "Required": required, "Flag": flag, "TrimmedData": trim,
            "Status": "Found" if found else "NoCompatibleSelection"}


def compatible_basis(a, field=None):
    report = compatible_basis_report(a, field)
    return report["Found"], report["Basis"], report["Sigma"]
