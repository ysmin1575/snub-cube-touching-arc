"""Independent continuum coverage and numerical-constant audit of snubs.pdf.

Run: python code/coverage_audit.py
Writes only out/coverage.json and out/coverage.md.
All numeric constants are exact rational/dyadic inputs, with Arb sign checks.
Symbolic inverse-map identities and a mathematical density proof cover the
continuum; sampled points are not used to establish coverage.
This script does not recreate the unavailable branch-and-bound certificates.
"""

from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import json
import platform
import time

from flint import arb, ctx
import sympy as sp


START = time.perf_counter()
BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "out"
SOURCE_DOI = "10.5281/zenodo.23118432"
ctx.prec = 256
CHECKS = []


def fl(text):
    return Q.from_float(float(text))


def ball(q):
    q = Q(q)
    return arb(str(q.numerator)) / arb(str(q.denominator))


def positive(name, q):
    q = Q(q)
    enclosure = ball(q)
    assert q > 0 and enclosure > 0, (name, q, enclosure)
    CHECKS.append({"name": name, "exact_margin": str(q),
                   "arb_margin": str(enclosure), "result": "positive"})


def equal(name, left, right):
    assert left == right, (name, left, right)
    CHECKS.append({"name": name, "left": str(left), "right": str(right),
                   "result": "exact equality"})


def symbolic(name, left, right):
    residuals = [sp.cancel(a - b) for a, b in zip(left, right)]
    assert all(a == 0 for a in residuals), (name, residuals)
    CHECKS.append({"name": name, "method": "exact symbolic identity",
                   "residuals": [str(a) for a in residuals]})


# Inputs below are stated formulas, not bundled third-party text.

ks, wf, s1, s2 = Q(3, 100), Q(3, 2), Q(37, 100), Q(1, 5)
c03, c06, c012 = fl("0.03"), fl("0.06"), fl("0.012")
f137, f120 = Q(137, 100), Q(6, 5)

# Constants needed for the complete inverse-map proof, including every bound
# where rounding can affect an endpoint. The rational equalities are exact.
positive("A3: fl(0.03) < kappa_s = 3/100", ks - c03)
positive("A3 W0: (137/100)/(3/2) = 137/150 < 1", 1 - f137 / wf)
positive("A3 W0: (6/5)/(3/2) = 4/5 < 1", 1 - f120 / wf)
positive("A3 Shigh: fl(1.4) > 137/100", fl("1.4") - f137)
positive("A3 Shigh: fl(1.25) > 6/5", fl("1.25") - f120)
equal("A3 covered region: fl(0.75)/(3/2) = 1/2", fl("0.75") / wf, Q(1, 2))
equal("A3 covered region: fl(0.5) = 1/2", fl("0.5"), Q(1, 2))
equal("A4 sheared first coordinate: 1 + 37/100 = 137/100", 1 + s1, f137)
equal("A4 sheared second coordinate: 1 + 1/5 = 6/5", 1 + s2, f120)
equal("A4 doubling: fl(0.06) = 2 fl(0.03)", c06, 2 * c03)
positive("A4 bound 1.37 under float convention is large enough", fl("1.37") - f137)
positive("A4 displayed fl(1.2) bound is 0.2 ulp too small", f120 - fl("1.2"))
equal("Tube epsilon lower endpoint uses identical dyadic constant", fl("-0.18"), -fl("0.18"))
equal("Tube epsilon upper endpoint uses identical dyadic constant", fl("-0.02"), -fl("0.02"))
positive("Tube parameter interval has positive width", fl("0.18") - fl("0.02"))
positive("Tube transverse endpoint fl(0.012) is positive", c012)
positive("Corner-to-tube epsilon overlap has positive width", c03 - fl("0.02"))
positive("Corner transverse radius exceeds tube transverse radius", c03 - c012)
positive("Shigh covered |d|<=1/2 fits |d|<=1", 1 - Q(1, 2))
positive("Shigh covered rho<=3/2 fits rho<=2", 2 - Q(3, 2))
positive("Shigh covered |f1'|<=3/4 fits SP first bound", f137 - Q(3, 4))
positive("Shigh covered |f2'|<=3/4 fits SP second bound", f120 - Q(3, 4))
equal("W0/Ws hemisphere cap upper endpoint", fl("1.0"), Q(1))
equal("A upper ratio endpoint", fl("2.0"), Q(2))

# Explicit symbolic inverse maps. Sup-norm selections and inequalities are
# proved in the prose below rather than assumed from samples.
eps, delta, rx, ry, rz, m, R, z, b0, b1, c2 = sp.symbols(
    "eps delta rx ry rz m R z b0 b1 c2", nonzero=True)
target = (eps, delta, rx, ry, rz)
rK = lambda e, d: sp.Matrix((d*b0/c2, d*b1/c2, -e))

for pivot in range(4):
    for sign in (-1, 1):
        # A chosen pivot guarantees max(|delta|, |r-rK|)=z and z>0.
        normed = sp.Matrix((delta, *(sp.Matrix((rx, ry, rz))-rK(eps, delta)))) / z
        reconstructed_delta = z * normed[0]
        reconstructed_r = rK(eps, reconstructed_delta) + z * sp.Matrix(normed[1:])
        symbolic(f"A2 tube inverse face {pivot},{sign}",
                 (eps, reconstructed_delta, *reconstructed_r), target)

for edge in range(2):
    for sign in (-1, 1):
        for face in range(3):
            for face_sign in (-1, 1):
                # eps/m,delta/m is an edge; r/R is a cube face.
                symbolic(f"A4 inverse edge {edge},{sign} face {face},{face_sign}",
                         (m*eps/m, m*delta/m,
                          m*(R/m)*rx/R, m*(R/m)*ry/R, m*(R/m)*rz/R), target)
for face in range(3):
    for sign in (-1, 1):
        symbolic(f"A4 B inverse face {face},{sign}",
                 (R*eps/R, R*delta/R, R*rx/R, R*ry/R, R*rz/R), target)

mu, d, rho, f1, f2, y, sigma = sp.symbols("mu d rho f1 f2 y sigma", nonzero=True)


def S(a, q, p, v, w):
    return (-a, a*q, a*p*(v+q*sp.Rational(37, 100)),
            a*p*(w+q*sp.Rational(1, 5)), a*p)


special = S(mu, d, rho, f1, f2)
symbolic("A3 W0 continuum inverse",
         S(sp.Rational(3, 100)*y*mu/(sp.Rational(3, 100)*y),
           y*d/y, rho, sp.Rational(3, 2)*y*f1/(sp.Rational(3, 2)*y),
           sp.Rational(3, 2)*y*f2/(sp.Rational(3, 2)*y)), special)
symbolic("A3 Ws continuum inverse on lower face sigma=1-rho",
         S(sp.Rational(3, 100)*y*mu/(sp.Rational(3, 100)*y),
           y*(1-rho)*d/(y*(1-rho)), 1-(1-rho),
           sp.Rational(3, 2)*y*(1-rho)*f1/(sp.Rational(3, 2)*y*(1-rho)),
           sp.Rational(3, 2)*y*(1-rho)*f2/(sp.Rational(3, 2)*y*(1-rho))), special)
symbolic("A3 C continuum inverse",
         S(mu, sigma*d/sigma, 1+sigma*(rho-1)/sigma,
           sp.Rational(3, 2)*sigma*f1/(sp.Rational(3, 2)*sigma),
           sp.Rational(3, 2)*sigma*f2/(sp.Rational(3, 2)*sigma)), special)
symbolic("A4 special inverse shear",
         S(m, delta/m, R/m, rx/R-(delta/m)*sp.Rational(37, 100),
           ry/R-(delta/m)*sp.Rational(1, 5)), (-m, delta, rx, ry, R))

# Finite combinatorial counts, not evidence that any leaf was certified.
family_counts = {"tube": 2*4, "A": (2*2*2*3-1)*2, "B": 2*3,
                 "W0": 1+2*3, "Ws": 1+2*3, "C": 2*4-1, "Shigh": 1}
equal("Appendix family total", sum(family_counts.values()), 82)
main_sd = {"outside": 93116, "global": 231708, "local": 9116}
main_sc = {"outside": 117356, "global": 1062181, "local": 105704, "Z": 33919}
equal("Snub dodecahedron leaf total", sum(main_sd.values()), 333940)
equal("Snub cube leaf total", sum(main_sc.values()), 1319160)
equal("Snub cube Z subtotal", 27421+6498, main_sc["Z"])
zoom_leaves = 2170487 + 220
binary_nodes = 2*zoom_leaves - 82
equal("Conditional full-binary forest node total", binary_nodes, 4341332)
positive("4341332 rounds to 4.3 million: lower rounding boundary", binary_nodes-4250000)
positive("4341332 rounds to 4.3 million: upper rounding boundary", 4350000-binary_nodes)
positive("Reported largest family 2.2M can fit overall 4.341332M", binary_nodes-2200000)

PROOF = r"""
1. A.1 개방성과 정리 5.3의 조밀성 단계.
차트의 u는 u·c=|c|²>0이므로 0이 아니고 D=1+|r|²>=1이다.
u×b도 0이 아니므로 e1=(u×b)/|u×b|, e2=(u/|u|)×e1로
연속적인 그림자 평면 좌표계를 정할 수 있다. 유한한 꼭짓점의 투영,
그 볼록껍질 및 지지함수는 매개변수에 관하여 균등하게 연속이다.
실제 엄격 포함 K+t⊂int H가 있으면 단위 방향 w 전체에 대해
hH(w)-hK(w)-t·w의 최소값은 양수이다. 단위 원은 콤팩트하므로
이 양의 여유는 충분히 작은 매개변수 변화 뒤에도 유지된다.
따라서 통과점 집합은 열려 있고 비통과점 집합은 닫혀 있다.
각도 좌표에서는 고정 평면에 투영하므로 같은 증명이 바로 적용된다.

E={r=0}∪{δ=0,r=rK0(ε,0)}는 R^5의 닫힌 집합이며 각 성분은
각각 차원 2, 1인 아핀 부분공간이므로 내부가 없다. TUBE는
(ε,δ,r-rK0)의 가역 아핀 좌표에서 비퇴화 닫힌 상자이고 CORNER도
비퇴화 닫힌 상자이다. 두 영역 모두 내부의 폐포다. 영역의 모든
점은 내부에서 E를 피하는 점열의 극한이다. 따라서 영역\E의
비통과가 증명되었다면 닫힘에 의해 전체 영역도 비통과다.
경계가 있다는 이유로 조밀성 논증을 반박할 근거는 없다.

2. A.2 튜브의 연속체 덮개.
a=(δ,r-rK0(ε,δ)), z=||a||∞로 두면 예외 곡선 밖에서
0<z<=fl(0.012)이다. a/z는 [-1,1]^4의 적어도 한 면에 놓인다.
그 면의 고정 좌표를 제외한 세 좌표를 자유변수로 하고
y0=ε,y1=z로 두면 정확히 원래 점을 얻는다. 면 모서리의 중복은
덮개에 문제가 없고 모든 거리변수 y1은 양수이다.

3. A.3 특수 차트의 연속체 덮개.
0<µ<=fl(0.03), |d|<=1, 0<ρ<=2, |f1'|<=137/100,
|f2'|<=6/5를 가정한다. ρ<=1/2이면
a=max(µ/(3/100),|d|,|f1'|/(3/2),|f2'|/(3/2))이다.
0<a<=1이고 (mhat,hhat)=(µ/((3/100)a),d/a,
f1'/((3/2)a),f2'/((3/2)a))는 [0,1]×[-1,1]^3의 단위
sup-sphere 위에 있다. mhat>0이므로 필요한 거리변수도 모두 양수다.
y0=a,y1=ρ로 W0에 넣으면 원래 S를 정확히 복원한다.

ρ>1/2이고 Shigh의 covered region 밖이면 해당 점 자체를
Shigh 좌표로 사용한다. 정확한 137/100,6/5는 각각 float64
1.4,1.25보다 작으므로 반올림 때문에 생기는 빈틈이 없다.
covered region 안이면 |d|<=1/2, |ρ-1|<=1/2, |f'|<=3/4이다.
z=(d,f1'/(3/2),f2'/(3/2),ρ-1), σ=||z||∞라 두면
제외된 arc point 이외에는 0<σ<=1/2이다.
ρ-1=-σ이면 a=max(µ/(3/100),|d|/σ,|f1'|/((3/2)σ),
|f2'|/((3/2)σ))∈(0,1]이고
(mhat,hhat)=(µ/((3/100)a),d/(aσ),f1'/((3/2)aσ),
f2'/((3/2)aσ))가 단위 sup-sphere에 놓인다.
y0=a,y1=σ인 Ws가 점을 복원하며 mhat,a,σ는 양수다.
그렇지 않으면 z/σ의 sup-norm이 1에 도달하는 면 중 하나는
제외된 g3=-1 면이 아니다. y0=µ,y1=σ로 해당 C 면에 넣으면
원래 점을 복원한다. ties도 Ws를 우선 사용하거나 다른 C 면을
선택하면 처리된다. C의 거리변수 µ,σ도 양수다.

Shigh의 covered region은 SP의 일부이므로 covered leaf를
W0/Ws/C에 위와 같이 넘기는 것에는 순환 의존이 없다.
실제 220개의 covered cells 전부가 그 지역 안인지 여부는
인증서가 없어 확인하지 못한다.

4. A.4 코너의 연속체 덮개.
m=max(|ε|,|δ|), R=||r||∞>0으로 둔다. R<=2m이면 m>0이며
y0=m,y2=R/m∈(0,2], (ε/m,δ/m)는 정사각형의 한 변,
r/R는 [-1,1]^3의 한 면에 놓인다. y2<=1/2 또는 y2>=1/2인
A chart가 적용된다. 일반 edge-face pair는 24개 중 23개다.
특수 pair ε/m=-1,rz/R=1이면 d=δ/m, f=(rx/R,ry/R)이고
x=S(m,d,R/m,f-ds)이다. |f1-ds1|<=1+37/100=137/100,
|f2-ds2|<=1+1/5=6/5이므로 정확한 SP 정의에 들어간다.
여기에는 반올림된 1.2가 아니라 정확한 6/5를 사용해야 한다.
R>2m이면 y0=R<=fl(0.06), |ε/R|,|δ/R|<1/2이므로 B chart다.
m=0이어도 이 분기는 그대로 작동한다. r=0과 특수 pair의
(d,f',ρ)=(0,0,1)만 남으며 앞의 닫힘 논증이 처리한다.

이 덮개 증명은 실제 비통과 인증서가 존재한다는 것을 증명하지
않는다. 정리 5.3에 적용하려면 모든 family cell의 유효성,
비음수 가중치, 엄격 양의 witness, 정확한 나눗셈 및 완전한
tiling이 별도로 증명되어야 한다. PDF의 Lemma 5.2는 그중
비음수 가중치 조건을 논리적으로 보장하지 못한다.
"""

FINDINGS = [
    {"id": "COV-01", "severity": "critical", "location": "Lemma 5.2 / Theorem 5.3, pp.6-7",
     "status": "PDF criterion has a mathematical gap",
     "finding": "가중치가 공통 부호라는 조건에서 λk>=0를 추론할 수 없다. 공통 음수 가중치와 양의 weighted witness는 엄격 포함과도 양립한다. 따라서 덮개 자체가 정확해도 적힌 인증 기준으로 정리 5.3을 얻을 수 없다.",
     "repair": "실제 인증서마다 λk>=0를 직접 검사하거나 공통 부호를 양수로 정규화한 가중치에 조건 (3)을 다시 검사해야 한다. 이 결함은 일반 보조정리의 반례이며 snub 자체가 Rupert라는 반례는 아니다.",
     "evidence": "독립 lemmas_audit.py의 exact rectangular-box counterexample; 본 스크립트는 기존 파일을 읽거나 수정하지 않는다."},
    {"id": "COV-02", "severity": "high", "location": "Sections 3,5.3,5.4,6 / supplement [11]",
     "status": "not checkable from supplied PDF",
     "finding": "SPEC.md, SPEC_SC.md, main/zoom certificates, producer 및 independent verifier가 제공되지 않았다. 실제 82개 family의 모든 leaf validity, weights, exact divisibility, positivity와 tilings; 220개 covered cell의 포함; chart gluing과 주 인증서의 완전한 root-domain coverage를 실행 검증할 수 없다.",
     "repair": "정확한 revision 또는 hashes와 함께 사양, 데이터/재생성 스크립트, 두 검증기를 제공하고 모든 검증 결과 및 정리의 선행조건을 확인해야 한다.",
     "evidence": "PDF explicitly delegates test lists, scales, file formats, sphere covers, trigonometric details and certificates to [11]."},
    {"id": "COV-03", "severity": "medium", "location": "Lemma 3.1 / local leaves, p.4",
     "status": "supplement-dependent domain and coverage",
     "finding": "principal rotation의 범위, r0<π, Fu의 실제 enclosure, S²×[0,tan(r0/2)]의 sphere/T cover 및 local witness nonnegativity를 확인할 자료가 없다. 이 내용이 없으면 local leaf 전체의 통과 배제를 재검증할 수 없다.",
     "repair": "명시된 principal-angle/r0 조건과 local certificate data를 제공한다.",
     "evidence": "Exact local formula is independently audited elsewhere; PDF alone does not instantiate Fu/r0/cover."},
    {"id": "COV-04", "severity": "low", "location": "Lemma A.4 second bullet, p.9",
     "status": "literal float64 inconsistency, no uncovered region after exact repair",
     "finding": "§5.3의 'decimal constants below'를 Appendix까지 적용하면 |f2-ds2|<=1.2라는 문장은 거짓이다. d=-1,f2=1에서 값은 정확한6/5이며 fl(1.2)는 그보다 1/22517998136852480 작다. A3가 SP에 정확한6/5를 쓰고 Shigh가1.25까지 넓어졌으므로 실제 덮개에 구멍이 생기지는 않는다.",
     "repair": "A4의1.37,1.2를 정확한137/100,6/5로 적는다.",
     "evidence": "Exact rational bound and 256-bit Arb sign check."},
    {"id": "COV-05", "severity": "low", "location": "Lemma 2.2, p.3",
     "status": "coordinate convention insufficiently explicit",
     "finding": "Ry의 회전 방향을 명시하지 않았다. 표시된 X(θ,φ)와 맞으려면 일반적인 active Ry와 반대 부호의 convention이 필요하다. 일관된 convention을 선택하면 reduction 자체는 성립한다.",
     "repair": "Ry 행렬을 본문 또는 사양에 명시한다.",
     "evidence": "Independent symbolic audit of A^T e3 in lemmas_audit.py."},
]

SUMMARY = {
    "coverage_continuum": "A1-A4 and density proof correct conditional on valid complete certificates and corrected nonnegative weights",
    "theorem_5_3": "not independently established: flawed signed-weight criterion and unavailable certificates",
    "family_counts": family_counts,
    "snub_dodecahedron_leaves": main_sd,
    "snub_cube_leaves": main_sc,
    "zoom_leaves": zoom_leaves,
    "binary_forest_nodes_if_full_binary": binary_nodes,
    "node_count_scope": "The paper does not specify binary branching here; 4341332 is a consistency calculation under that usual model, not reconstructed actual nodes.",
    "precision": ctx.prec,
    "exact_dyadic_constants": {a: str(fl(a)) for a in
        ("0.03", "0.06", "0.012", "0.18", "0.02", "1.4", "1.25", "1.37", "1.2", "0.75", "0.5")},
    "counterexample_to_A4_literal_float_bound": {
        "mu": str(c03), "d": "-1", "rho": "1", "f": ["0", "1"],
        "f2_minus_d_s2": "6/5", "float_1.2": str(fl("1.2")),
        "margin": str(f120-fl("1.2")),
        "inside_CORNER": "epsilon=-mu, delta=-mu, r=(0,mu,mu); all magnitudes <=fl(0.03)<fl(0.06)",
        "outside_exceptional_arc": True,
        "is_still_covered": True,
    },
    "reviewed_source_text_sha256": 'f6988b1e63ba16a3f7493b2c0fb0a5451da355d454a6e568637b5c25a045c319',
    "reviewed_source_pdf_sha256": 'b79ee5f3477e5e88b4967c50672e06293250f40875c7c633b74b1823f560786a',
    "python": platform.python_version(),
    "runtime_seconds": round(time.perf_counter()-START, 3),
    "checks": CHECKS,
    "findings": FINDINGS,
    "continuum_proof": PROOF.strip(),
}

OUT.mkdir(exist_ok=True)
(OUT / "coverage.json").write_text(json.dumps(SUMMARY, ensure_ascii=False, indent=2), encoding="utf-8")

md = """# 덮개와 전체 증명 구조의 독립 감사

A.1–A.4의 역매개화와 정리 5.3의 폐포·조밀성 단계는 연속체 전체에 대해 성립합니다. 점 샘플링을 근거로 삼지 않았습니다. 다만 정리 5.3과 두 snub의 비-Rupert 결론을 재검증했다는 뜻은 아닙니다. Lemma 5.2의 부호 조건은 잘못 적혀 있고, 실제 인증서와 검증 사양은 제공되지 않았습니다.

| 우선순위 | 위치 | 판정 |
|---|---|---|
| Critical | Lemma 5.2 → Theorem 5.3 | 공통 부호로 비음수 가중치를 추론하는 결함. 실제 모든 witness의 양수 정규화와 조건 (3) 재검사가 필요합니다. |
| High | §§3,5–6 / [11] | 데이터·사양·검증기가 없어 실제 전역 leaf 배제와 tiling을 확인할 수 없습니다. |
| Medium | Lemma 3.1 / local leaves | principal-angle/r0 조건, Fu, sphere/T 덮개는 별도 확인이 필요합니다. |
| Low | A.4 | fl(1.2)는 정확한6/5보다 작습니다. 정확한 분수로 고치면 덮개는 유지됩니다. |
| Low | Lemma 2.2 | Ry 부호 convention을 명시해야 합니다. 일관된 선택이면 reduction이 성립합니다. |

## 수치와 개수

82개 family의 구성은 8+46+6+7+7+7+1=82로 일치합니다. Snub dodecahedron leaf는93,116+231,708+9,116=333,940이며, snub cube leaf는117,356+1,062,181+105,704+33,919=1,319,160입니다. Z leaf는27,421+6,498=33,919입니다. 이는 합계 일치만 확인한 것이며 실제 leaf의 수나 타일링을 확인한 결과는 아닙니다.

Zoom leaf는2,170,487+220=2,170,707입니다. 82개의 완전 이진 tree라는 조건에서는 node가2×2,170,707−82=4,341,332로, 논문의 '4.3 million' 반올림 표현과 일치합니다. PDF는 해당 위치에서 이진 분할을 명시하지 않으므로 이 계산은 조건부 정합성 검사입니다. '가장 큰 W0 face m 2.2 million nodes'도 전체 개수와 모순되지 않지만 실제 node 수는 확인할 수 없습니다.

fl(0.03)<3/100의 정확한 차이는1/900719925474099200입니다. fl(1.4)>137/100, fl(1.25)>6/5를 exact dyadic rational과256-bit Arb로 검증했습니다. Shigh의 넓어진1.25 경계는 SP의 정확한6/5 경계를 포함합니다. fl(0.06)=2fl(0.03)이고, covered region의0.75/(3/2)=1/2도 정확합니다. 전체 margin과 Arb enclosure는 coverage.json에 기록했습니다.

## 연속체 덮개 증명

""" + PROOF.strip() + "\n\n## 검증 범위\n\n" + \
"공개판은 원문의 해당 식을 재구성하여, 정확한 수치 상수의 부등식, 각 역매개화의 symbolic identity, family/leaf 합계를 검증했습니다. 실제82개 family의 certificate 통과, 220개 covered leaf의 위치, sphere covering, main tiling 및Z-leaf의 interval enclosure는 검증하지 못했습니다. Theorem4.2와 일반 보조정리의 세부 독립 검증은 별도 감사 산출물에 있습니다.\n"
(OUT / "coverage.md").write_text(md, encoding="utf-8")
print(json.dumps({"checks": len(CHECKS), "family_count": sum(family_counts.values()),
                  "runtime_seconds": SUMMARY["runtime_seconds"],
                  "outputs": ["out/coverage.json", "out/coverage.md"]}, ensure_ascii=False))
