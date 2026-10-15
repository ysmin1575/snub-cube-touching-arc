"""Single-process exact-field/Arb check of Tom 7's snub-cube tie.

Run from the repository root: python -B code/tom_arc_certify.py
Only out/tom_arc_certificate.json is written. No search or Lean build.
The field primitives are reused from the independently audited arc_certify.
"""
from pathlib import Path
from itertools import permutations, product
import json
import time
import arc_certify as ac

OUT = Path(__file__).resolve().parents[1] / "out"
Z, O, T = ac.ZERO, ac.ONE, ac.T


def scale(a, v):
    return tuple(ac.mul(a, b) for b in v)


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def reflect(v, n):
    return sub(v, scale(ac.div(2 * ac.dot(v, n), ac.dot(n, n)), n))


def S(v):
    return (-v[0], v[2], v[1])


def tom_vertices():
    seed = (T, O, -ac.mul(T, T))
    result = []
    for perm in permutations(range(3)):
        parity = sum(perm[i] > perm[j] for i in range(3) for j in range(i + 1, 3)) % 2
        for signs in product((-1, 1), repeat=3):
            if (-1) ** parity * signs[0] * signs[1] * signs[2] == 1:
                result.append(tuple(signs[i] * seed[perm[i]] for i in range(3)))
    assert len(result) == len({tuple(map(str, v)) for v in result}) == 24
    return result


def main():
    start = time.monotonic()
    ac.ctx.prec = 384
    V = tom_vertices()
    t2 = ac.mul(T, T)
    x, e, y = (O, Z, t2), (-t2, Z, O), (Z, O, Z)
    k = O - ac.div(O, T)
    assert ac.sign(k) > 0
    assert ac.dot(x, e) == Z and ac.dot(x, x) == ac.dot(e, e)
    assert S(x) == tuple(ac.div(a, k) for a in ac.C)
    assert S(e) == tuple(ac.div(a, k) for a in ac.B)
    assert S(y) == (Z, Z, O)
    # Regenerate Tom's orbit separately, then match exact coordinates.
    mapped = [scale(T, S(v)) for v in ac.V]
    assert {tuple(map(str, v)) for v in V} == {tuple(map(str, v)) for v in mapped}
    V = mapped  # coordinate matching, not Tom's source slot numbering
    K = [reflect(v, x) for v in V]
    assert K == [scale(T, S(reflect(v, ac.C))) for v in ac.V]

    # Cayley w=mu*y has N=(1-mu^2)I+2mu^2 yy^T+2mu[y]_cross.
    # For EVERY vertex the horizontal numerator is
    # (1+mu^2)*(e-mu*x).v; the vertical coordinate is v_y.
    # These are coefficient equalities over Q(t), not sampled rotations.
    for v in V:
        ycross = (v[2], Z, -v[0])
        n2 = sub(scale(2 * ac.dot(y, v), y), v)
        lhs = [ac.dot(e, v), 2 * ac.dot(e, ycross) + ac.dot(x, v),
               ac.dot(e, n2) + 2 * ac.dot(x, ycross), ac.dot(x, n2)]
        rhs = [ac.dot(e, v), -ac.dot(x, v), ac.dot(e, v), -ac.dot(x, v)]
        assert lhs == rhs
        assert [v[1], 2 * ycross[1], n2[1]] == [v[1], Z, v[1]]

    def area(i, j, p):
        edge, point = sub(V[j], V[i]), sub(p, V[i])
        return tuple(ac.mul(ac.dot(edge, n), point[1]) -
                     ac.mul(edge[1], ac.dot(point, n)) for n in (e, x))

    hull = [2, 0, 18, 17, 23, 22, 4, 14]
    star, end = 2 * O - T, O / 3
    counts = {"hole": 0, "plug": 0}
    multiplier = ac.div(t2, k)
    for plug, points in ((False, V), (True, K)):
        for i, j in zip(hull, hull[1:] + hull[:1]):
            for h, p in enumerate(points):
                a, d = area(i, j, p)
                assert (a, d) == tuple(ac.mul(multiplier, z) for z in ac.area(i, j, h, plug))
                assert ac.sign(a) >= 0 and ac.sign(a + ac.mul(star, d)) >= 0
                counts["plug" if plug else "hole"] += 1
    # Geometric hull checks cover all mu in [0,star], including zero.
    ac.certificate(hull, Z, star, False)
    ac.certificate(hull, Z, star, True)
    a, d = area(2, 0, K[8])
    factor = 2 * ac.mul(ac.mul(T, t2), T - O)
    assert d == -factor and a == ac.mul(factor, star)
    assert ac.sign(factor) > 0
    # The critical outer edge remains supporting ALL THE WAY to mu=1/3.
    for p in V:
        aa, dd = area(2, 0, p)
        assert ac.sign(aa) >= 0 and ac.sign(aa + ac.mul(end, dd)) >= 0
    # Consequently every mu in (star,1/3] fails closed containment.
    theta = ac.evaluate(star).atan()
    checks = {}
    for name, mu in (("mu=0.165", ac.arb("0.165")),
                     ("eps=0.165 rad", ac.arb("0.165").tan()),
                     ("eps=0.1594 rad", ac.arb("0.1594").tan())):
        assert mu > ac.evaluate(star) and mu < ac.evaluate(end)
        det = ac.evaluate(a) + mu * ac.evaluate(d)
        assert det < 0
        # Divide by the positive edge length in orthonormal screen coordinates.
        edge = sub(V[0], V[2])
        horizontal_scale = (ac.evaluate(ac.dot(x, x)) * (1 + mu * mu)).sqrt()
        dx = (ac.evaluate(ac.dot(edge, e)) + mu * ac.evaluate(ac.dot(edge, x))) / horizontal_scale
        length = (dx * dx + ac.evaluate(edge[1]) ** 2).sqrt()
        signed_distance = det / horizontal_scale / length
        checks[name] = {"mu_ball": mu.str(70), "signed_area_ball": det.str(70),
                        "signed_perpendicular_distance_ball": signed_distance.str(70)}
    # The audit's 40-place decimal interval is retained without rounding.
    lo = "0.1593506261845691446953081463393635724078"
    hi = "0.1593506261845691446953081463393635724079"
    assert theta > ac.arb(lo) and theta < ac.arb(hi)
    result = {
        "status": "VERIFIED", "precision_bits": ac.ctx.prec,
        "tom_vertices": "O orbit of (t,1,-t^2); unnormalized",
        "map": "audit v = S(tom v)/t, S(X,Y,Z)=(-X,Z,Y), det(S)=1",
        "view": "u=x-mu*e1, x=(1,0,t^2), e1=(-t^2,0,1)",
        "rotation": "Cayley w=mu*e_y; R_y(2*atan(mu))=H_u*H_x=K_u*K_x",
        "parameter_relation": "mu=tan(eps)=tan(theta); eps=theta; audit chart epsilon=-mu",
        "mu_star_exact": ac.expression(star),
        "mu_star_ball": ac.evaluate(star).str(70),
        "eps_star_enclosure": [lo, hi],
        "degree_ball": (theta * 180 / ac.arb.pi()).str(65),
        "relative_rotation_endpoint_ball": (2 * theta).str(65),
        "domain": "0<=mu<=1/3, equivalently 0<=eps<=atan(1/3)",
        "containment_iff": "mu<=2-t; clearance is zero there (equal vertical supports)",
        "all_affine_endpoint_checks": counts,
        "critical_vertices": {"edge_start": [ac.expression(z) for z in V[2]],
                              "edge_end": [ac.expression(z) for z in V[0]],
                              "plug_before_Kx": [ac.expression(z) for z in V[8]]},
        "critical_determinant": "2*t^3*(t-1)*(2-t-mu)",
        "certified_failures": checks,
        "scope": "This checks the touching family only, not Tom 7's full Lean certificates or a global non-Rupert theorem.",
        "elapsed_seconds": time.monotonic() - start,
    }
    assert result["elapsed_seconds"] < 300
    OUT.mkdir(exist_ok=True)
    (OUT / "tom_arc_certificate.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
