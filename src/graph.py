#!/usr/bin/env python3

import multiprocessing as mp
import networkx as nx
import itertools
from search import pointToPoly, areScalable, areScalableRes, areRelPrime, areBinaryScalable, dyadicCofactors
from math import factorial, log2, floor
from tqdm import tqdm
from functools import cache
import pynauty
import sys

def mp_scalable(args):
    n, logfloor, i, j = args

    v = min(ord2(i - j), logfloor)
    if v > 0:
        power = pow(2, v)
        if i % power not in {n % power, 0}:
            return i, j, False

    p = pointToPoly(n, 1, i)
    q = pointToPoly(n, 1, j)
    return i, j, areScalable(p, q)

def mp_resultant(args):
    n, logfloor, i, j = args

    v = min(ord2(i - j), logfloor)
    if v > 0:
        power = pow(2, v)
        if i % power not in {n % power, 0}:
            return i, j, False

    p = pointToPoly(n, 1, i)
    q = pointToPoly(n, 1, j)
    return i, j, areScalableRes(p, q)

def mp_binscal(args):
    n, i, j = args
    p = pointToPoly(n, 1, i)
    q = pointToPoly(n, 1, j)
    return i, j, areBinaryScalable(p, q)

def mp_relprime(args):
    n, i, j = args
    p = pointToPoly(n, 1, i)
    q = pointToPoly(n, 1, j)
    return i, j, areRelPrime(p, q)

def mp_dyadic(args):
    n, i, j = args
    p = pointToPoly(n, 1, i)
    q = pointToPoly(n, 1, j)
    return i, j, dyadicCofactors(p, q)

def binomial(n, k):
    return (factorial(n) // factorial(k)) // factorial(n - k)

def ord2(n):
    k = 0
    while n & 1 == 0:
        k += 1
        n >>= 1

    return k

def makeGraph(n, cond="resultant"):
    G = nx.Graph()
    G.add_nodes_from(range(1, n))
    pairs = itertools.combinations(range(1, n), 2)
    points = ((n, i, j) for i, j in pairs)

    match cond:
        case "scalable":
            check = mp_scalable
            # the scalable check now has a shortcut if we know floor(log2(n)),
            # so pass that to its function.
            points = ((n, floor(log2(n)), i, j) for n, i, j in points)
        case "resultant":
            check = mp_resultant
            points = ((n, floor(log2(n)), i, j) for n, i, j in points)
        case "relprime":
            check = mp_relprime
        case "dyadic":
            check = mp_dyadic
        case "binary scalable":
            check = mp_binscal
        case _:
            raise ValueError("cond must be one of: 'scalable', 'resultant', 'relprime', 'dyadic', 'binary scalable'")

    with mp.Pool() as p:
        total = binomial(n-1, 2)
        res = p.imap_unordered(check, points, chunksize=max(1, total//1000))
        res = tqdm(res, total=binomial(n-1, 2))
        edges = ((i, j) for i, j, scalable in res if scalable)
        G.add_edges_from(edges)

    return G

def makeGraphNauty(n):
    G = makeGraph(n)
    N = pynauty.Graph(n-1)

    for k in range(1, n):
        N.connect_vertex(k-1, [v-1 for v in G[k].keys()])

    return N

def largestClique(n, minimize=False):
    """
    find the cliques of largest size in the graph S(n).
    if minimize is true, then choose the clique with the fewest average number
    of terms in the bezout cofactors.
    """
    G = makeGraph(n)
    cliques = nx.clique.find_cliques(G)

    if minimize:
        def weight(clique):
            avg = 0
            for i, j in itertools.combinations(clique, 2):
                p = pointToPoly(n, 1, i)
                q = pointToPoly(n, 1, j)
                g, r, s = p.xgcd(q)
                avg -= sum(r.coeffs()) / len(r.coeffs())

            return avg / binomial(len(clique), 2)

        key = weight
    else:
        key = len

    return max(cliques, key=key)

def largestCliques(n):
    G = makeGraph(n)
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

    G = makeGraph(n)
    cs = list(nx.clique.find_cliques(G))
    max_size = max(map(len, cs))
    return [c for c in cs if len(c) == max_size and n // 2 in c]

@cache
def largestCliqueSize(n):
    return len(largestClique(n))

def drawColors(n):
    G = makeGraph(n)
    colors = nx.coloring.greedy_color(G, "independent_set")
    unique_colors = set(colors.values())
    graph_color_to_mpl_color = dict(zip(unique_colors, mpc.TABLEAU_COLORS))
    node_colors = [graph_color_to_mpl_color[colors[n]] for n in G.nodes()]

    pos = nx.shell_layout(G)
    nx.draw(
        G,
        pos,
        with_labels=True,
        node_color=node_colors,
        font_color="#333333",
    )

def graphColors(n):
    G = makeGraph(n)
    colors = nx.coloring.greedy_color(G, "independent_set")
    # this is horredenously inefficient
    inv_colors = {color: sorted([v for v, c in colors.items() if c == color]) for color in set(colors.values())}

    return inv_colors

@cache
def nColors(n):
    return len(graphColors(n))

if __name__ == "__main__":
    try:
        n = int(sys.argv[1])
    except:
        print("trouble parsing n\narguments:", sys.argv[1:], file=sys.stderr)
        sys.exit(1)

    if n < 3:
        print("n must be >= 3")
        sys.exit(1)

    G = makeGraph(n, cond="resultant")
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
    max_half, max_half_size = get_largest(lambda L: all(x <= n//2 for x in L))

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
