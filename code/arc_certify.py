"""Independent exact-field and Arb certificate for Theorem 4.2.

No numerical hull is trusted: every listed edge is checked against all points.
fmpq_poly arithmetic is in Q[t]/(t^3-t^2-t-1). Arb certifies every nonzero sign.
"""
from pathlib import Path
import itertools
import json
import time
from flint import arb, ctx, fmpq, fmpq_poly

OUT = Path(__file__).resolve().parents[1] / 'out'
F = fmpq_poly([-1, -1, -1, 1])
ZERO = fmpq_poly([])
ONE = fmpq_poly([1])
T = fmpq_poly([0, 1])
ctx.prec = 384

def mul(a, b):
    return (a*b) % F

def div(a, b):
    g, x, _ = b.xgcd(F)
    assert g == ONE
    return mul(a, x)

def dot(a, b):
    return sum((mul(x, y) for x, y in zip(a, b)), ZERO)

def root_ball():
    lo, hi = fmpq(1), fmpq(2)
    for _ in range(300):
        mid = (lo+hi)/2
        value = mid**3-mid**2-mid-1
        if value < 0:
            lo = mid
        else:
            hi = mid
    assert lo**3-lo**2-lo-1 < 0 < hi**3-hi**2-hi-1
    # p'= (3t+1)(t-1); both stationary values are negative.
    assert fmpq(-1,3)**3-fmpq(-1,3)**2-fmpq(-1,3)-1 == fmpq(-22,27)
    ball = arb(str((lo+hi)/2), str((hi-lo)/2))
    assert 3*ball**2-2*ball-1 > 0
    return ball

TB = root_ball()

def evaluate(p):
    r = arb(0)
    for a in reversed(p.coeffs()):
        r = r*TB + arb(str(a))
    return r

def sign(p):
    if p == ZERO:
        return 0
    value = evaluate(p)
    if value > 0:
        return 1
    if value < 0:
        return -1
    raise ArithmeticError('Unresolved sign: '+str(p)+' '+str(value))

def expression(p):
    return [str(x) for x in p.coeffs()] + ['0']*(3-len(p))

def vertices():
    result=[]
    for signs in itertools.product((-1, 1), repeat=3):
        v=[signs[i]*x for i,x in enumerate((ONE, T*T-T-ONE, T))]
        for perm in itertools.permutations(range(3)):
            parity=sum(perm[i]>perm[j] for i in range(3) for j in range(i+1, 3))%2
            if parity == sum(x==1 for x in signs)%2:
                result.append([v[k] for k in perm])
    assert len(result)==24 and len({tuple(map(str,v)) for v in result})==24
    return result

V=vertices()
C=[T*T-T-2*ONE, T*T-T, ZERO]
B=[C[1], -C[0], ZERO]
CV=[dot(v,C) for v in V]
BV=[dot(v,B) for v in V]
Z=[v[2] for v in V]

def area(i,j,k,plug=False):
    """Oriented determinant A+s D, using (b.s v+s c.v,z)."""
    a=mul(BV[j]-BV[i],Z[k]-Z[i])-mul(Z[j]-Z[i],BV[k]-BV[i])
    ck=-CV[k] if plug else CV[k]
    d=mul(CV[j]-CV[i],Z[k]-Z[i])-mul(Z[j]-Z[i],ck-CV[i])
    return a,d

def affine_at(poly,s):
    return poly[0]+mul(poly[1],s)

def certificate(hull,lo,hi,plug):
    # The candidate boundary is a simple, strictly convex cyclic polygon:
    # all vertices are distinct throughout the interval, every turn is positive,
    # and (below) every supporting halfplane contains every hole vertex.
    for i,j in itertools.combinations(hull,2):
        if Z[i]!=Z[j]:
            assert sign(mul(Z[i]-Z[j],Z[i]-Z[j]))>0
        else:
            xdiff=(BV[i]-BV[j],CV[i]-CV[j])
            assert sign(affine_at(xdiff,lo))*sign(affine_at(xdiff,hi))>0
    for pos,i in enumerate(hull):
        turn=area(hull[pos-1],i,hull[(pos+1)%len(hull)])
        assert sign(affine_at(turn,lo))>0 and sign(affine_at(turn,hi))>0
    counts={'positive_both':0,'zero_lo':0,'zero_hi':0,'identical_zero':0}
    records=[]
    for i,j in zip(hull,hull[1:]+hull[:1]):
        assert sign(mul(BV[j]-BV[i],BV[j]-BV[i])+mul(Z[j]-Z[i],Z[j]-Z[i]))>0
        for k in range(24):
            a,d=area(i,j,k,plug)
            lv,hv=affine_at((a,d),lo),affine_at((a,d),hi)
            ls,hs=sign(lv),sign(hv)
            assert ls>=0 and hs>=0,(i,j,k,plug,str(a),str(d),ls,hs)
            category=('identical_zero' if a==ZERO and d==ZERO else
                      'zero_lo' if ls==0 else 'zero_hi' if hs==0 else 'positive_both')
            counts[category]+=1
            records.append(dict(edge=[i,j],vertex=k,A=expression(a),D=expression(d),
                                signs=[ls,hs]))
    # Nonzero area, together with the verified supporting cyclic boundary,
    # ensures the intersection is exactly the bounded convex hull.
    a,d=area(hull[0],hull[1],hull[2])
    assert sign(affine_at((a,d),lo))>0 and sign(affine_at((a,d),hi))>0
    return counts,records

def enclosure(x, digits=65):
    m,r,e=x.mid_rad_10exp(digits)
    from decimal import Decimal,localcontext
    with localcontext() as dc:
        dc.prec=digits+30
        scale=Decimal(10)**int(e)
        return [format(Decimal(int(m-r))*scale,'f'),format(Decimal(int(m+r))*scale,'f')]

def main():
    start=time.monotonic()
    # Candidate only. Its correctness is established by exhaustive edge checks.
    hull=[2,0,18,17,23,22,4,14]
    a,d=area(2,0,8,True)
    assert sign(a)>0 and sign(d)<0
    factor=2*mul(T-ONE,T-ONE)
    assert d==-factor and a==mul(factor,2*ONE-T)
    star=div(-a,d)
    sb=evaluate(star)
    theta=sb.atan()
    report_interval=['0.1593506261845691446953081463393635724078',
                     '0.1593506261845691446953081463393635724079']
    assert theta-arb(report_interval[0])>0 and arb(report_interval[1])-theta>0
    theta0=arb('0.16071324478583').atan()
    assert theta-theta0 > 0
    assert theta0+arb('1e-14')-theta > 0
    other=area(18,17,5,True)
    assert affine_at(other,star)==ZERO
    hole_counts,hole_records=certificate(hull,ZERO,star,False)
    plug_counts,plug_records=certificate(hull,ZERO,star,True)
    # Extend the supporting edge beyond the first failure with rigorous bounds.
    end=(theta0+arb('1e-6')).tan()
    for k in range(24):
        a0,d0=area(2,0,k)
        v0=affine_at((a0,d0),star)
        vend=evaluate(a0)+end*evaluate(d0)
        assert sign(v0)>=0
        if a0==ZERO and d0==ZERO:
            assert vend.contains(0)
        else:
            assert vend > 0
    assert evaluate(d)*(theta0+arb('1e-14')).tan()+evaluate(a) < 0
    assert evaluate(d)*arb('0.16071324478583')+evaluate(a) > 0
    # Reproduce the PDF's ancillary support, derivative and triangle claims
    # throughout I=[theta0-1e-9,theta0+1e-6], rather than at one point.
    theta_I=theta0+arb('0.0000004995','0.0000005005')
    s_I=theta_I.tan()
    assert theta_I.sin()>0 and theta_I.cos()>0
    normalized_derivative=-(sb*theta_I.sin()+theta_I.cos())
    assert normalized_derivative<0
    for k in range(24):
        if k not in (2,0):
            aa,dd=area(2,0,k)
            assert evaluate(aa)+s_I*evaluate(dd)>0
    v0=V.index([-ONE,T*T-T-ONE,T])
    for test in (area(2,0,v0),area(0,v0,8,True),area(v0,2,8,True)):
        assert evaluate(test[0])+s_I*evaluate(test[1])>0
    # Opposite sense: search for one supporting edge and a strictly failing plug.
    opposite=None
    minus_end=-end
    for i,j in itertools.permutations(range(24),2):
        valid=True
        for k in range(24):
            aa,dd=area(i,j,k)
            if sign(aa)<0:
                valid=False; break
            if aa==ZERO and dd==ZERO:
                continue
            if not evaluate(aa)+minus_end*evaluate(dd)>=0:
                valid=False; break
        if not valid:
            continue
        for k in range(24):
            aa,dd=area(i,j,k,True)
            if aa==ZERO and sign(dd)>0:
                opposite=dict(edge=[i,j],plug=k,A=expression(aa),D=expression(dd),
                              D_ball=evaluate(dd).str(65))
                break
        if opposite:
            break
    assert opposite is not None
    assert opposite['edge']==[14,2] and opposite['plug']==2
    assert area(14,2,2,True)==(ZERO,4*(T*T-T-ONE))
    # Exact half-turn exchange of the two critical obstructions.
    negxy=lambda v:[-v[0],-v[1],v[2]]
    assert negxy(V[2])==V[17] and negxy(V[0])==V[18] and negxy(V[8])==V[5]
    assert area(2,0,8,True)==area(18,17,5,True)
    result=dict(status='VERIFIED',precision_bits=ctx.prec,t_ball=TB.str(80),
                s_star_Qt=expression(star),s_star_ball=sb.str(70),
                theta_star_ball=theta.str(70),theta_star_enclosure=enclosure(theta),
                report_enclosure=report_interval,
                degree_ball=(theta*180/arb.pi()).str(60),
                theta_minus_theta0=(theta-theta0).str(60),hull=hull,
                critical_area_factor_Qt=expression(factor),
                critical_area_formula='2*(t-1)^2*(2-t-s)',
                PDF_near_threshold_support_derivative_triangle='VERIFIED on I',
                hole_counts=hole_counts,plug_counts=plug_counts,
                opposite=opposite,elapsed_seconds=time.monotonic()-start,
                proof='All determinants affine in s; exact field zeros, Arb endpoint signs. Convex interpolation covers every s in [0,s*]. A supporting edge and negative slope exclude all s in (s*,tan(theta0+1e-6)].')
    (OUT/'arc_certificate.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    (OUT/'arc_halfplanes.json').write_text(json.dumps(dict(hole=hole_records,plug=plug_records),indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
