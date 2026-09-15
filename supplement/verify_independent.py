#!/usr/bin/env python3
"""Second implementation of the scalar certificate verifier.

No import from certify.py or verify_certificate.py. Box coordinates, weights,
coefficients and comparisons use fractions.Fraction exactly. Only powers and
square roots are enclosed on a 160-bit dyadic grid, with directed integer rounding.
This cross-checks the same analytic reduction; it is not a proof-assistant check.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json
import sys
import time

BITS = 160
T = 1 << BITS
COORD = 1 << 100
GAMMA = Q(5, 3)


def need(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def ceildiv(a: int, b: int) -> int:
    need(b > 0, 'nonpositive divisor')
    return (a + b - 1) // b


def halfpower(x: Q, twice_exponent: int, upper: bool) -> Q:
    """Enclose x**(twice_exponent/2) for x >= 0 by monotone arithmetic."""
    need(x >= 0 and twice_exponent >= 0, 'invalid power')
    rnd = ceildiv if upper else lambda a, b: a // b
    base = rnd(x.numerator * T, x.denominator)
    result, z, exponent = T, base, twice_exponent // 2
    while exponent:
        if exponent % 2:
            result = rnd(result * z, T)
        exponent //= 2
        if exponent:
            z = rnd(z * z, T)
    if twice_exponent % 2:
        root = isqrt(base * T)
        if upper and root * root < base * T:
            root += 1
        result = rnd(result * root, T)
    return Q(result, T)


@lru_cache(maxsize=500)
def mesh(m_upper: int, subdivision: int):
    need(subdivision > 0, 'empty subdivision')
    kmax = m_upper.bit_length() + 3
    coarse = [0] + [COORD >> k for k in range(kmax, 0, -1)] + [COORD]
    vertices = {0, COORD}
    for lo, hi in zip(coarse, coarse[1:]):
        for j in range(subdivision + 1):
            z = lo + (hi - lo) * j // subdivision
            vertices.update((z, COORD-z))
    points = tuple(Q(z, COORD) for z in sorted(vertices))
    weights = []
    # Integrate the two linear interpolation basis functions exactly.
    for l, r in zip(points, points[1:]):
        z = r-l
        pairs = []
        for power in (1, 2):
            total = (r**(power+1)-l**(power+1))/(power+1)
            moment = (r**(power+2)-l**(power+2))/(power+2)
            pairs.append(((r*total-moment)/z, (moment-l*total)/z))
        weights.append(tuple(pairs))
    return points, tuple(weights)


def certify_h(ml: int, mh: int, al: Q, ah: Q, h: Q) -> str:
    dl, du = 1-ah*ah, 1-al*al
    if h <= 1/(1+Q(mh, 2)*du):
        return 'elementary'
    need(dl > 0, 'non-elementary polar bound at a=1')
    upper = (halfpower(1+du*h, mh+2, True)
             -halfpower(al, 2*(mh+2), False)
             -Q(ml+2, 2)*dl*(1+h))
    need(upper < 0, 'uncertified polar lower bound')
    return 'chord'


def certify_box(ml: int, mh: int, al: Q, ah: Q, h: Q, sub: int) -> str:
    if h > 2*ah-ah*ah:
        return 'polar'
    d = (ah*ah+h)/2
    need(d < 1, 'invalid denominator bound')
    points, weights = mesh(mh, sub)
    vals1, vals2 = [], []
    for t in points:
        envelope = max(Q(0), 1-h*t-al*al*t*(1-t))
        value = halfpower(envelope, ml-1, True)
        vals1.append(value)
        vals2.append(value/(1-d*t))
    k1 = sum((w[0][0]*vals1[j]+w[0][1]*vals1[j+1]
              for j, w in enumerate(weights)), Q(0))
    k2 = sum((w[1][0]*vals2[j]+w[1][1]*vals2[j+1]
              for j, w in enumerate(weights)), Q(0))
    f = halfpower(1-h, ml, True)
    rr = GAMMA*ah**2*mh*k1
    if f+rr+Q(mh, mh+1)*ah < 1:
        return 'direct'
    ee = halfpower(1-h, ml+1, True)
    cc = GAMMA*ah**3*mh*(mh+1)*k2/2
    alpha = Q(mh, 2)*(1-al*al)
    if 2*cc <= al and 0 < alpha <= Q(ml-1, 2):
        ratio = Q(mh, 2)/alpha*halfpower(alpha/(1+alpha), ml+1, True)
        if ratio < Q(1, 4):
            return 'boundary'
    if 2*cc <= ah and ah+ee < 1:
        return 'center-monotone'
    if 2*cc > ah and ee < 1 and 4*(ah+cc)**3 < 27*cc*(1-ee)**2:
        return 'center-cubic'
    raise ValueError('no exclusion condition is satisfied')


def verify(path: Path) -> dict:
    raw = path.read_bytes()
    doc = json.loads(raw)
    need(doc['precision_bits'] == 100, 'coordinate scale mismatch')
    need(doc['m_start'] == 5 and doc['m_stop_exclusive'] == 1000000, 'wrong domain')
    cover = defaultdict(list)
    counts, polar = Counter(), Counter()
    started = time.monotonic()
    for number, row in enumerate(doc['records']):
        need(len(row) == 7, 'invalid record')
        ml, mh, als, ahs, hs, tag, subdivision = row
        need(isinstance(ml, int) and isinstance(mh, int) and 5 <= ml <= mh < 1000000, 'invalid degree')
        al, ah, h = (Q(int(s), COORD) for s in (als, ahs, hs))
        need(0 <= al < ah <= 1 and 0 <= h < 1, 'invalid box')
        try:
            polar[certify_h(ml, mh, al, ah, h)] += 1
            actual = certify_box(ml, mh, al, ah, h, subdivision)
        except Exception as exc:
            raise ValueError(f'record {number}: {exc}') from exc
        counts[actual] += 1
        cover[ml, mh].append((al, ah))
        if number % 1000 == 0:
            print(f'verified {number+1} records', file=sys.stderr, flush=True)
    previous = 4
    for (ml, mh), ranges in sorted(cover.items()):
        need(ml == previous+1, 'degree gap or overlap')
        endpoint = Q(0)
        for lo, hi in sorted(ranges):
            need(lo == endpoint, 'parameter gap or overlap')
            endpoint = hi
        need(endpoint == 1, 'missing boundary')
        previous = mh
    need(previous == 999999, 'missing final degree')
    return {'status': 'PASS', 'implementation': 'separate exact-Fraction checker; 160-bit directed power bounds',
            'imports_original_checker': False, 'm_range': [5,999999], 'degree_blocks': len(cover),
            'parameter_boxes': len(doc['records']), 'criteria': dict(counts), 'polar_bounds': dict(polar),
            'certificate_sha256': sha256(raw).hexdigest(), 'seconds': round(time.monotonic()-started, 3)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', nargs='?', type=Path, default=Path('certificate.json'))
    args = parser.parse_args()
    print(json.dumps(verify(args.certificate), indent=2))
