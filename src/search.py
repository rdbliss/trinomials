import itertools
import multiprocessing as mp
from tqdm import tqdm
from flint import fmpq_poly, fmpz_poly
from collections import Counter
from math import factorial

def search(maxcoeff, setsize, n, extraCheck=None):
    vals = range(1, maxcoeff + 1, 2)
    coeffs = itertools.product(vals, repeat=setsize)
    ks = itertools.combinations_with_replacement(range(n), setsize)
    points = itertools.product(coeffs, ks)
    points = (p for p in points if max(Counter(zip(*p)).values()) == 1)

    if extraCheck:
        points = filter(extraCheck, points)

    points = ((n, cs, ks) for cs, ks in points)
    work = int(((maxcoeff + 1) // 2)**setsize * (n**setsize / factorial(setsize)))

    goods = []
    with mp.Pool() as p:
        # break the work into 1000 ~equally sized chunks
        results = tqdm(p.imap_unordered(check, points, chunksize=max(1, work//1000)), total=work)
        for res in results:
            if res:
                goods.append(res)
                break

    real_goods = []
    for point in goods:
        n, cs, ks = point
        polys = tuple(pointToPoly(n, c, k) for c, k in zip(cs, ks))
        real_goods.append(polys)
        print(point)
        print(polys)

    return real_goods

def isDyadic(x):
    while x & 1 == 0:
        x >>= 1

    return x == 1

def polyIsDyadic(p):
    return isDyadic(p.denom())

def areScalable(f, g):
    g, r, s = f.xgcd(g)

    if g != 1:
        return False

    return isDyadic(r.denom()) and isDyadic(s.denom())

def areScalableRes(f, g):
    return isDyadic(abs(int(f.resultant(g))))

def areRelPrime(f, g):
    g, r, s = f.xgcd(g)

    return g == 1

def dyadicCofactors(f, g):
    g, r, s = f.xgcd(g)

    return isDyadic(r.denom()) and isDyadic(s.denom())

def areBinaryScalable(f, g):
    g, r, s = f.xgcd(g)

    return g == 1 and r.denom() == 1 and s.denom() == 1

def pointToPoly(n, c, k):
    if k == 0:
        cs = [1 - c] + [0] * (n - 1) + [1]
    else:
        cs = [1] + [0] * (k - 1) + [-c] + [0] * (n - k - 1) + [1]

    return fmpq_poly(cs)

def check(argTuple):
    n, coeffs, ks = argTuple
    fs = [pointToPoly(n, c, k) for c, k in zip(coeffs, ks)]

    pairs = itertools.combinations(fs, 2)
    if all(areScalable(f, g) for f, g in pairs):
        return argTuple

    return None

def binomial(n, k):
    return (factorial(n) // factorial(k)) // factorial(n - k)

def mp_scalablePoint(args):
    n, i, j = args
    p = pointToPoly(n, 1, i)
    q = pointToPoly(n, 1, j)
    return n, areScalable(p, q)

def countScalablePairs(N):
    def pointGenerator():
        for k in range(1, N+1):
            pairs = itertools.combinations(range(1, k), 2)
            points = ((k, i, j) for i, j in pairs)
            yield from points

    total = (N * (N**2 - 3 * N + 2)) // 6
    points = tqdm(pointGenerator(), total=total)

    results = [0] * N
    with mp.Pool() as p:
        res = p.imap_unordered(mp_scalablePoint, points)
        for n, scalable in res:
            results[n-1] += scalable

    return results
