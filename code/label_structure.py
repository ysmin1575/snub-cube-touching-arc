"""Independent snub-cube reconstruction and structural audit.

Run: python label_structure.py
Dependencies: sympy, python-flint (mpmath is a SymPy dependency).
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import permutations, product, combinations
from pathlib import Path
import json
import sympy as sp
from flint import arb, fmpq, ctx

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "out"
ctx.prec = 256


@dataclass(frozen=True)
class Q:
    """Exact element a+b*t+c*t^2, with t^3=t^2+t+1."""
    a: F = F(0)
    b: F = F(0)
    c: F = F(0)

    @staticmethod
    def of(x):
        return x if isinstance(x, Q) else Q(F(x))

    def __add__(self, other):
        o = Q.of(other)
        return Q(self.a + o.a, self.b + o.b, self.c + o.c)

    __radd__ = __add__

    def __neg__(self):
        return Q(-self.a, -self.b, -self.c)

    def __sub__(self, other):
        return self + -Q.of(other)

    def __rsub__(self, other):
        return Q.of(other) - self

    def __mul__(self, other):
        o = Q.of(other)
        z = [F(0)] * 5
        for i, x in enumerate((self.a, self.b, self.c)):
            for j, y in enumerate((o.a, o.b, o.c)):
                z[i + j] += x * y
        # t^4=2*t^2+2*t+1; t^3=t^2+t+1.
        return Q(z[0] + z[3] + z[4], z[1] + z[3] + 2*z[4], z[2] + z[3] + 2*z[4])

    __rmul__ = __mul__

    def inverse(self):
        assert self != ZERO
        t = sp.Symbol("t")
        expr = self.a + self.b*t + self.c*t*t
        z = sp.Poly(sp.invert(expr, t**3-t*t-t-1, t), t)
        return Q(*(F(z.nth(i)) for i in range(3)))

    def __truediv__(self, other):
        return self * Q.of(other).inverse()

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        z = ONE
        for _ in range(n):
            z = z*self
        return z

    def ball(self, t_ball):
        return aq(self.a) + aq(self.b)*t_ball + aq(self.c)*t_ball*t_ball

    def coeff(self):
        return [str(x) for x in (self.a, self.b, self.c)]


ZERO, ONE, T = Q(), Q(F(1)), Q(F(0), F(1))
INV_T = T*T-T-1


def aq(x):
    x = F(x)
    return arb(fmpq(x.numerator, x.denominator))


def root_ball():
    def p(x):
        return x*x*x-x*x-x-1
    lo, hi = F(1), F(2)
    for _ in range(240):
        mid = (lo+hi)/2
        if p(mid) < 0:
            lo = mid
        else:
            hi = mid
    assert p(lo) < 0 < p(hi)
    # p'=(3t+1)(t-1)>0 on (1,2), while both stationary values
    # p(-1/3)=-22/27 and p(1)=-2 are negative: unique real root.
    assert p(F(-1, 3)) == F(-22, 27) and p(F(1)) == -2
    return arb(fmpq(((lo+hi)/2).numerator, ((lo+hi)/2).denominator),
               fmpq(((hi-lo)/2).numerator, ((hi-lo)/2).denominator)), lo, hi


TB, TLO, THI = root_ball()


def sign(x):
    if x == ZERO:
        return 0
    v = x.ball(TB)
    if v > 0:
        return 1
    if v < 0:
        return -1
    raise ArithmeticError(f"Unresolved sign: {x}")


def dot(v, w):
    return sum((a*b for a, b in zip(v, w)), ZERO)


def add(v, w):
    return tuple(a+b for a, b in zip(v, w))


def sub(v, w):
    return tuple(a-b for a, b in zip(v, w))


def scale(a, v):
    return tuple(a*x for x in v)


def cross(v, w):
    return (v[1]*w[2]-v[2]*w[1], v[2]*w[0]-v[0]*w[2], v[0]*w[1]-v[1]*w[0])


def vertices():
    base = (ONE, INV_T, T)
    result = []
    for perm in permutations(range(3)):
        parity = sum(perm[i] > perm[j] for i in range(3) for j in range(i+1, 3)) % 2
        for signs in product((-1, 1), repeat=3):
            if sum(s > 0 for s in signs) % 2 == parity:
                result.append(tuple(signs[i]*base[perm[i]] for i in range(3)))
    assert len(result) == len(set(result)) == 24
    return result


RAW = vertices()
C = (T*T-T-2, T*T-T, ZERO)
A4 = (ZERO, ZERO, ONE)
B = cross(C, A4)
C2 = dot(C, C)


def reflect(v, u=C):
    return sub(v, scale(2*dot(v, u)/dot(u, u), u))


PINNED = {0: 1, 6: 8, 2: 20, 11: 23, 8: 16, 3: 0, 20: 21,
          21: 3, 23: 19, 17: 22}
REMAINING = iter(i for i in range(24) if i not in PINNED.values())
LABEL_TO_RAW = [PINNED[i] if i in PINNED else next(REMAINING) for i in range(24)]
V = [RAW[i] for i in LABEL_TO_RAW]


def area_coeff(a, d, p):
    """Signed projected area in the plane coordinates (B+s*C, a4).

    A common positive horizontal normalization is omitted.
    Returns its exact constant and linear coefficient in s.
    """
    e, q = sub(d, a), sub(p, a)
    return tuple(dot(e, w)*q[2]-e[2]*dot(q, w) for w in (B, C))


def area(a, d, p, s):
    f0, f1 = area_coeff(a, d, p)
    return f0+s*f1


def labeling_audit():
    """Task 1: exact index constraints and certified geometric pinning."""
    assert len(set(LABEL_TO_RAW)) == 24
    assert sub(V[0], V[6]) == C
    assert reflect(V[2]) == V[2] and reflect(V[11]) == V[11]
    f0, f1 = area_coeff(V[8], V[3], reflect(V[20]))
    g0, g1 = area_coeff(V[21], V[23], reflect(V[17]))
    assert (f0, f1) == (g0, g1)
    assert sign(f0) == 1 and sign(f1) == -1
    sstar = -f0/f1
    assert area(V[8], V[3], reflect(V[20]), sstar) == ZERO
    for a, d in ((8, 3), (21, 23)):
        assert all(sign(area(V[a], V[d], V[j], sstar)) == 1
                   for j in range(24) if j not in (a, d))
    # The exact root lies between the two edge endpoints, with positive
    # coordinates at vertices 8 and 3 in the triangle (8,3,0).
    denom = area(V[8], V[3], V[0], sstar)
    bary8 = area(V[3], V[0], reflect(V[20]), sstar)/denom
    bary3 = area(V[0], V[8], reflect(V[20]), sstar)/denom
    bary0 = area(V[8], V[3], reflect(V[20]), sstar)/denom
    assert sign(denom) == sign(bary8) == sign(bary3) == 1
    assert bary0 == ZERO and bary8+bary3 == ONE
    # The fixed-vertex edges named in the corner discussion are actual
    # oriented supporting edges at an interior point of the claimed arc.
    for a, d in ((2, 8), (23, 11)):
        assert all(sign(area(V[a], V[d], V[j], Q(F(1, 10)))) >= 0
                   for j in range(24))
    s0 = F("0.16071324478583")
    assert sstar.ball(TB) > aq(s0)
    # atan is an Arb interval operation; compare against the exact paper
    # bracket in theta, rather than equating rounded decimal constants.
    theta = sstar.ball(TB).atan()
    theta0 = aq(s0).atan()
    assert theta < theta0+aq(F(1, 10**14))
    h = lambda v: (-v[0], -v[1], v[2])
    assert h(V[2]) == V[11] and h(V[20]) == V[17]
    assert h(V[8]) == V[23] and h(V[3]) == V[21]
    data = {
        "status": "VERIFIED: a labeling consistent with all named PDF index constraints",
        "original_SPEC_labeling_recovered": False,
        "nonuniqueness": "The 14 unnamed labels are assigned deterministically; v0/v6 also admit another ordered chord pair.",
        "enumeration": "lexicographic permutation outer loop; signs (-1,+1) inner loop; #positive signs parity equals permutation parity",
        "field": "Q[t]/(t^3-t^2-t-1); coefficients are [constant,t,t^2] rational strings",
        "t_interval_exact_rational": [str(TLO), str(THI)],
        "t_arb": str(TB),
        "label_to_raw": LABEL_TO_RAW,
        "vertices": [[x.coeff() for x in v] for v in V],
        "c": [x.coeff() for x in C],
        "b": [x.coeff() for x in B],
        "threshold_s_exact": sstar.coeff(),
        "threshold_s_arb": str(sstar.ball(TB)),
        "threshold_theta_arb": str(theta),
        "edge_area_coefficients": [f0.coeff(), f1.coeff()],
        "threshold_barycentric_exact": [bary8.coeff(), bary3.coeff(), bary0.coeff()],
        "pinned_indices": sorted(PINNED),
        "proof": [
            "v0-v6=c holds as an exact Q(t) equality.",
            "K0 fixes v2 and v11 exactly.",
            "The two named oriented edge tests have identical exact affine polynomials in s=tan(theta), with positive constant and negative slope certified by Arb.",
            "At their exact common zero, both named edges support the shadow strictly away from their endpoints.",
            "At that zero the image of K0*v20 has positive barycentric coordinates at v8 and v3 in triangle (8,3,0), and zero at v0.",
            "The named fixed-vertex edges support the shadow at s=1/10.",
            "Half-turn maps (v2,v20,v8,v3) to (v11,v17,v23,v21), reversing oriented endpoints as required.",
            "Arb verifies theta0 < theta_star < theta0+10^-14.",
        ],
    }
    OUT.mkdir(exist_ok=True)
    (OUT/"labeling.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    return data


def mat_apply(m, v):
    return tuple(dot(row, v) for row in m)


def structure_audit():
    """Task 2: exact structure, and an exhaustive upper bound on G."""
    x = sp.Symbol("t")
    assert sp.Poly(x**3-x*x-x-1, x).is_irreducible
    assert dot(B, B) == C2 and dot(B, C) == dot(C, A4) == dot(B, A4) == ZERO
    assert sign(C2) == 1
    parallel = [(i, j) for i, j in combinations(range(24), 2)
                if cross(sub(V[i], V[j]), C) == (ZERO, ZERO, ZERO)]
    assert len(parallel) == 5 and (0, 6) in parallel
    lookup = {v:i for i, v in enumerate(V)}
    kappa = {i: lookup[reflect(v)] for i, v in enumerate(V) if reflect(v) in lookup}
    fixed = [i for i, j in kappa.items() if i == j]
    swapped = [(i, j) for i, j in kappa.items() if i < j]
    assert len(kappa) == 12 and len(swapped) == 5 and fixed == [2, 11]
    assert all(kappa[kappa[i]] == i for i in kappa)
    assert all(tuple(-x for x in v) not in lookup for v in V)

    # An arbitrary orthogonal symmetry is determined by the images of
    # three independent vertices. Enumerate EVERY possible ordered image
    # triple in V^3 satisfying their exact Gram matrix. This gives the
    # required upper bound, not merely a subgroup of signed permutations.
    base_indices = next(tri for tri in combinations(range(24), 3)
                        if dot(V[tri[0]], cross(V[tri[1]], V[tri[2]])) != ZERO)
    p, q, r = (V[i] for i in base_indices)
    detu = dot(p, cross(q, r))
    idetu = detu.inverse()
    inverse_rows = [scale(idetu, z) for z in (cross(q, r), cross(r, p), cross(p, q))]
    gram = [[dot(a, b) for b in (p, q, r)] for a in (p, q, r)]
    dots = [[dot(a, b) for b in V] for a in V]
    gram_candidates, proper, improper, matrices = 0, [], [], []
    for i in range(24):
        for j in range(24):
            if dots[i][j] != gram[0][1]:
                continue
            for k in range(24):
                if dots[i][k] != gram[0][2] or dots[j][k] != gram[1][2]:
                    continue
                gram_candidates += 1
                images = (V[i], V[j], V[k])
                m = tuple(tuple(sum((images[a][b]*inverse_rows[a][c] for a in range(3)), ZERO)
                                for c in range(3)) for b in range(3))
                mm = [[sum((m[a][i]*m[a][j] for a in range(3)), ZERO)
                       for j in range(3)] for i in range(3)]
                assert mm == [[ONE if i == j else ZERO for j in range(3)] for i in range(3)]
                detm = dot(m[0], cross(m[1], m[2]))
                assert detm in (ONE, -ONE)
                mapped = [mat_apply(m, v) for v in V]
                if not all(v in lookup for v in mapped):
                    continue
                perm = [lookup[v] for v in mapped]
                matrices.append(m)
                (proper if detm == ONE else improper).append(perm)
    assert len(set(matrices)) == len(matrices)
    assert len(proper) == 24 and len(improper) == 0
    # All 24 matrices also visibly are signed coordinate permutations,
    # but this observation is not used to bound the full symmetry group.
    assert all(all(sum(z != ZERO for z in row) == 1 and all(z in (ZERO, ONE, -ONE) for z in row)
                   for row in m) for m in matrices)
    assert all(all(sum(m[i][j] != ZERO for i in range(3)) == 1 for j in range(3)) for m in matrices)

    # Symbolic rational matrix identity for every real eps,delta. This
    # proves Ku*K0=N(rK0)/D(rK0), so also verifies the cross-product sign.
    eps, delta = sp.symbols("eps delta", real=True)
    sy = lambda z: sp.Rational(z.a.numerator,z.a.denominator)+sp.Rational(z.b.numerator,z.b.denominator)*x+sp.Rational(z.c.numerator,z.c.denominator)*x*x
    modpoly = x**3-x*x-x-1

    def reduce_expr(z):
        num, den = sp.fraction(sp.cancel(z))
        pn = sp.Poly(num, x, domain=sp.QQ.frac_field(eps, delta)).rem(sp.Poly(modpoly, x, domain=sp.QQ.frac_field(eps, delta)))
        pd = sp.Poly(den, x, domain=sp.QQ.frac_field(eps, delta)).rem(sp.Poly(modpoly, x, domain=sp.QQ.frac_field(eps, delta)))
        assert not pd.is_zero
        return pn.is_zero or sp.cancel(pn.as_expr()) == 0

    cs, bs = sp.Matrix([sy(z) for z in C]), sp.Matrix([sy(z) for z in B])
    a4 = sp.Matrix([0,0,1])
    u = cs+eps*bs+delta*a4
    c2 = (cs.T*cs)[0]
    r0 = -eps*a4+delta*bs/c2
    ident = sp.eye(3)
    kk0 = ident-2*u*u.T/(u.T*u)[0]
    k0 = ident-2*cs*cs.T/c2
    skew = sp.Matrix([[0,-r0[2],r0[1]],[r0[2],0,-r0[0]],[-r0[1],r0[0],0]])
    d0 = 1+(r0.T*r0)[0]
    n0 = (1-(r0.T*r0)[0])*ident+2*r0*r0.T+2*skew
    residual = kk0*k0-n0/d0
    assert all(reduce_expr(z) for z in residual)
    gibbs_cross = cs.cross(u)/(cs.T*u)[0]-r0
    assert all(reduce_expr(z) for z in gibbs_cross)

    result = {
        "status": "VERIFIED",
        "field_irreducible": True,
        "b_norm_equals_c_norm_exact": True,
        "parallel_unordered_chords": parallel,
        "reflection_vertex_map": {str(i):j for i,j in kappa.items()},
        "reflection_fixed": fixed,
        "reflection_swapped_pairs": swapped,
        "rotation_group_order": len(proper),
        "improper_symmetry_count": len(improper),
        "minus_I_maps_any_vertex_to_vertex": False,
        "group_upper_bound_method": "exhaustive exact Gram-preserving ordered image triples of three independent vertices, followed by exact vertex-set preservation",
        "independent_vertex_triple": base_indices,
        "gram_candidate_count": gram_candidates,
        "group_vertex_permutations": proper,
        "gibbs_identity_exact": "Ku*K0=N(-eps*a4+delta*b/|c|^2)/D(-eps*a4+delta*b/|c|^2)",
        "gibbs_axis_orientation": "c cross u; right-handed rotation angle 2 angle(c,u), with identity at u=c",
        "gibbs_denominators": "<c,u>=|c|^2>0; |u|^2=|c|^2*(1+eps^2)+delta^2>0; D=1+|r|^2>0",
    }
    (OUT/"structure.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


if __name__ == "__main__":
    result = labeling_audit()
    print("Task 1 VERIFIED: consistent labeling saved to", "out/labeling.json")
    print("theta*:", result["threshold_theta_arb"])
    print("label-to-raw:", LABEL_TO_RAW)
    result = structure_audit()
    print("Task 2 VERIFIED:", {k:result[k] for k in ("parallel_unordered_chords", "reflection_fixed", "reflection_swapped_pairs", "rotation_group_order", "gram_candidate_count", "gibbs_identity_exact")})
