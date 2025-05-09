#!/usr/bin/env python3

"""
graph.py - construct trinomial graphs

This file contains the main graph creation functions, including the parallel
processing boilerplate.
"""

import sys
import itertools
import multiprocessing as mp
from math import comb as binomial
from functools import cache
from collections import deque

import networkx as nx
from utils import trinom, dres, are_scalable
from tqdm import tqdm
import pynauty


def make_graph(n, cond="congruence"):
    """
    make the trinomial graph on {1, 2, ..., n - 1}, which has the edge {i, j}
    if and only if the resultant of x^n - x^i + 1 and x^n - x^j + 1 is +-a
    power of 2.
    """
    G = nx.Graph()

    match cond:
        case "congruence":
            resultant_count = 0
            with mp.Pool() as p:
                points = ((n, i) for i in range(1, n))
                rets = p.imap_unordered(mp_congruence_visit, points)
                for edges, count in tqdm(rets, total=n - 1):
                    resultant_count += count
                    G.add_edges_from(edges)

                return G
        case "scalable":
            check = mp_scalable
        case "resultant":
            check = mp_resultant
        case "relprime":
            check = mp_relprime
        case _:
            raise ValueError(
                "cond must be one of: 'congruence', 'scalable', "
                "'resultant', 'relprime'"
            )

    pairs = itertools.combinations(range(1, n), 2)
    points = ((n, i, j) for i, j in pairs)
    with mp.Pool() as p:
        total = binomial(n - 1, 2)
        res = p.imap_unordered(check, points, chunksize=max(1, total // 1000))
        res = tqdm(res, total=binomial(n - 1, 2))
        edges = ((i, j) for i, j, scalable in res if scalable)
        G.add_edges_from(edges)

    return G


def pow2_indepdence_check(n, i, j):
    v = ord2(i - j)
    if v > 0:
        power = pow(2, v)
        if i % power not in {n % power, 0}:
            return False

    return True


def count_resultants(n):
    """
    count the number of resultants computed by the naive "check every pair"
    approach after excluding our power of 2 independent sets.
    """
    check = pow2_indepdence_check
    pairs = itertools.combinations(range(1, n), 2)
    return sum(check(n, i, j) for i, j in pairs)


def mp_scalable(args):
    n, i, j = args

    if not pow2_indepdence_check(n, i, j):
        return i, j, False

    return i, j, are_scalable(*args)


def mp_resultant(args):
    n, i, j = args

    if not pow2_indepdence_check(n, i, j):
        return i, j, False

    return i, j, dres(*args)


def mp_relprime(args):
    n, i, j = args
    p = trinom(n, i)
    q = trinom(n, j)
    return i, j, p.gcd(q) == 1


def mp_congruence_visit(args):
    n, i = args
    visit = set(range(1, n - i))

    ret_edges = deque()
    resultant_count = 0
    while visit:
        d = min(visit)
        resultant_count += 1
        if not dres(n, i, i + d):
            visit -= set(range(d, n - i, d))
        else:
            edges = ((k, k + d) for k in range(i, n - d, d))
            ret_edges.extend(edges)
            visit -= {d}

    return ret_edges, resultant_count


def make_graph_nauty(n):
    G = make_graph(n)
    N = pynauty.Graph(n - 1)

    for k in range(1, n):
        N.connect_vertex(k - 1, [v - 1 for v in G[k].keys()])

    return N


def ord2(n):
    k = 0
    while n & 1 == 0:
        k += 1
        n >>= 1

    return k


def largestClique(n):
    """
    find the cliques of largest size in the graph S(n).
    if minimize is true, then choose the clique with the fewest average number
    of terms in the bezout cofactors.
    """
    G = make_graph(n)
    cliques = nx.clique.find_cliques(G)

    return max(cliques, key=len)


def largestCliques(n):
    G = make_graph(n)
    cs = nx.clique.find_cliques(G)
    size = largestCliqueSize(n)
    return [sorted(c) for c in cs if len(c) == size]


def largestMidCliques(n):
    """
    find the cliques of maximal size which contain the vertex n / 2. n must be
    even.
    """
    if n % 2 != 0:
        raise ValueError("need even n for a mid clique")

    G = make_graph(n)
    cs = list(nx.clique.find_cliques(G))
    max_size = max(map(len, cs))
    return [c for c in cs if len(c) == max_size and n // 2 in c]


@cache
def largestCliqueSize(n):
    return len(largestClique(n))


def graphColors(n):
    G = make_graph(n)
    colors = nx.coloring.greedy_color(G, "independent_set")
    # this is horredenously inefficient
    inv_colors = {
        color: sorted([v for v, c in colors.items() if c == color])
        for color in set(colors.values())
    }

    return inv_colors


@cache
def nColors(n):
    return len(graphColors(n))


def main(argv):
    try:
        n = int(sys.argv[1])
    except:
        print("trouble parsing n\narguments:", sys.argv[1:], file=sys.stderr)
        sys.exit(1)

    if n < 3:
        print("n must be >= 3")
        sys.exit(1)

    G = make_graph(n)
    cs = list(nx.clique.find_cliques(G))
    cs = [sorted(c) for c in cs]

    def get_largest(cond):
        restricted_cs = [c for c in cs if cond(c)]
        if not restricted_cs:
            return [], 0
        max_size = max(map(len, restricted_cs))
        max_restricted = [c for c in restricted_cs if len(c) == max_size]
        return max_restricted, max_size

    max_cs, max_size = get_largest(lambda L: True)

    def is_centered(L):
        return set(L) == set(n - x for x in L)

    max_centered, max_centered_size = get_largest(is_centered)
    max_half, max_half_size = get_largest(lambda L: all(x <= n // 2 for x in L))

    print("# clique number, maximal symmetric clique, maximal clique <= n/2")
    print(f"{max_size}, {max_centered_size}, {max_half_size}")
    print()

    print("maximum cliques")
    print("-----------")
    for c in max_cs:
        print(c)

    print()

    print("maximal symmetric cliques")
    print("-----------")
    for c in max_centered:
        print(c)

    print()

    print("maximial skew cliques")
    print("-----------")
    for c in max_half:
        print(c)


if __name__ == "__main__":
    main(sys.argv)
