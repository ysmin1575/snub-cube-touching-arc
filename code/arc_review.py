"""Independent review of arc_certify.py, without changing that script."""
from fractions import Fraction as F
from pathlib import Path
import json
import arc_certify as ac
import label_structure as ls

OUT = Path(__file__).resolve().parents[1] / "out"


def from_poly(p):
    z = [F(str(x)) for x in p.coeffs()]
    z += [F(0)]*(3-len(z))
    assert len(z) == 3
    return ls.Q(*z)


def main():
    # Independently regenerated exact vertices; match root enumeration by
    # coordinates, rather than assuming its index order or trusting JSON.
    root_vertices = [tuple(from_poly(x) for x in row) for row in ac.V]
    assert set(root_vertices) == set(ls.V)
    root_to_label = [ls.V.index(row) for row in root_vertices]
    hull = [root_to_label[i] for i in (2,0,18,17,23,22,4,14)]
    assert len(set(hull)) == 8
    star = 2-ls.T
    assert ls.sign(star) > 0
    counts = {"nonzero_height_difference_edges":0,"horizontal_noncollapse_edges":0,
              "strict_turns_at_zero":0,"zero_turns_at_zero":0}
    for n,(a,b) in enumerate(zip(hull,hull[1:]+hull[:1])):
        # Exact equality of every determinant coefficient to the independent
        # script's coefficient, for both the hole and the reflected plug.
        for k in range(24):
            rk = root_to_label.index(k)
            for plug in (False,True):
                ours = ls.area_coeff(ls.V[a],ls.V[b],ls.reflect(ls.V[k]) if plug else ls.V[k])
                theirs = tuple(from_poly(z) for z in ac.area((2,0,18,17,23,22,4,14)[n],
                                                               (2,0,18,17,23,22,4,14)[(n+1)%8],rk,plug))
                assert ours == theirs
                assert ls.sign(ours[0]) >= 0
                assert ls.sign(ours[0]+star*ours[1]) >= 0
        # Explicit full-interval noncollapse: if the edge is horizontal,
        # its affine horizontal separation has the same nonzero sign at
        # the two endpoints. Otherwise its height separation is nonzero.
        e = ls.sub(ls.V[b],ls.V[a])
        if e[2] != ls.ZERO:
            counts["nonzero_height_difference_edges"] += 1
        else:
            x0 = ls.dot(e,ls.B)
            x1 = x0+star*ls.dot(e,ls.C)
            assert ls.sign(x0) == ls.sign(x1) != 0
            counts["horizontal_noncollapse_edges"] += 1
        c = hull[(n+2)%8]
        q0,q1 = ls.area_coeff(ls.V[a],ls.V[b],ls.V[c])
        s0,s1 = ls.sign(q0),ls.sign(q0+star*q1)
        assert s0 >= 0 and s1 > 0
        counts["strict_turns_at_zero" if s0 > 0 else "zero_turns_at_zero"] += 1
    # At s=0, distinct vertices can have identical shadow coordinates;
    # these eight selected vertices must remain distinct there as well.
    for n,a in enumerate(hull):
        for b in hull[:n]:
            dz=ls.V[a][2]-ls.V[b][2]
            dx=ls.dot(ls.sub(ls.V[a],ls.V[b]),ls.B)
            assert dz != ls.ZERO or dx != ls.ZERO
    # Exact critical value and relevant named edges agree with task 1.
    a,d = ls.area_coeff(ls.V[8],ls.V[3],ls.reflect(ls.V[20]))
    assert -a/d == star
    assert tuple(from_poly(z) for z in ac.area(2,0,8,True)) == (a,d)
    assert tuple(root_to_label[i] for i in (2,0,8)) == (8,3,20)
    assert tuple(root_to_label[i] for i in (18,17,5)) == (21,23,17)
    report = {
        "status": "VERIFIED: independent review and exact additional geometric checks",
        "root_to_label": root_to_label,
        "hull_in_labels": hull,
        "extra_checks": counts,
        "all_384_affine_determinants_independently_equal": True,
        "all_edges_nondegenerate_for_entire_positive_interval": True,
        "all_turns_strict_for_positive_s_and_nonnegative_at_zero": True,
        "critical_value_exact": "tan(theta_star)=2-t",
        "rotation_equivalence": "eps=-tan(theta),delta=0 gives Gibbs r=tan(theta)*a4, hence positive right-handed rotation +2theta. Ku leaves every shadow point fixed. No translation is used.",
        "review_finding": "Full-interval noncollapse, distinct hull vertices, strict cyclic turns and all supporting halfplanes certify the actual bounded convex hull; no certificate flaw remains.",
    }
    (OUT/"arc_review.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()
