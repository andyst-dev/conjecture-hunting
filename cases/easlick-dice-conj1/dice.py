"""Graph G_n of Easlick (arXiv:2610.40321): vertices = n-sided dice with faces in [n]
(multisets of size n from {1..n}, standard die included); arc X -> Y iff
#{(i,j): x_i > y_j} == floor(n^2/2) + 1  (uniform, minimal winning advantage, Eq. (1)).

All arithmetic is integer (exact). Counts use
    beats(x, y) = sum_v m_x(v) * #{j : y_j < v}
where m_x(v) is the multiplicity of face v in x.
"""
import itertools
import numpy as np
import scipy.sparse as sp


def dice(n):
    """All non-decreasing n-tuples over [n], in lexicographic order."""
    return [t for t in itertools.combinations_with_replacement(range(1, n + 1), n)]


def target(n):
    return n * n // 2 + 1


def mult_matrix(D, n):
    M = np.zeros((len(D), n), dtype=np.int64)
    for k, t in enumerate(D):
        for v in t:
            M[k, v - 1] += 1
    return M


def below_matrix(M):
    """C[k, v-1] = #{faces of die k strictly below v}."""
    C = np.zeros_like(M)
    C[:, 1:] = np.cumsum(M, axis=1)[:, :-1]
    return C


def beats_count(x, y):
    """Direct definition, used for cross-checks."""
    return sum(1 for a in x for b in y if a > b)


def build_graph(n, block=2000, rows=None, cols=None):
    """CSR adjacency of G_n (optionally the subgraph induced on index lists)."""
    D = dice(n)
    M = mult_matrix(D, n)
    C = below_matrix(M)
    T = target(n)
    R = np.arange(len(D)) if rows is None else np.asarray(rows)
    K = R if cols is None else np.asarray(cols)
    Ct = C[K].T.copy()
    src_parts, dst_parts = [], []
    for s in range(0, len(R), block):
        blk = M[R[s:s + block]] @ Ct
        i, j = np.nonzero(blk == T)
        src_parts.append((i + s).astype(np.int32))
        dst_parts.append(j.astype(np.int32))
    src = np.concatenate(src_parts)
    dst = np.concatenate(dst_parts)
    A = sp.csr_matrix((np.ones(len(src), dtype=np.int8), (src, dst)), shape=(len(R), len(K)))
    return D, A
