"""
This file contains some common helper functions needed to create trinomial
graphs, such as interfacing with flint, computing resultants, checking if
something is a power of 2, and so on.
"""

from flint import fmpq_poly


def dres(n, i, j):
    """
    check if the edge {i, j} exists in the nth trinomial graph, meaning that a
    certain resultant is a signed power of 2.
    """
    p = trinom(n, i)
    q = trinom(n, j)
    r = abs(int(p.resultant(q)))
    return is_dyadic(r)


def dres_poly(f, g):
    return is_dyadic(abs(f.resultant(g)))


def are_scalable(n, i, j):
    p = trinom(n, i)
    q = trinom(n, j)
    return are_scalable_poly(p, q)


def are_scalable_poly(f, g):
    """
    check if the denominators of the bezout cofactors of f and g are powers of
    2.

    this should be equivalent to dyadicallyResolve(f, g), but this is here as a
    sanity check.
    """
    g, r, s = f.xgcd(g)

    if g != 1:
        return False

    return is_dyadic(r.denom()) and is_dyadic(s.denom())


def trinom(n, k):
    x = fmpq_poly([0, 1])
    return x**n - x**k + 1


def is_dyadic(x):
    if x == 0:
        return False

    while x & 1 == 0:
        x >>= 1

    return x == 1
