"""Exact symbolic audit of the separation, reduction and zoom lemmas.

Run: python code/lemmas_audit.py
Only out/lemmas.json and out/lemmas.md are written.
"""

from itertools import product
from pathlib import Path
import json
import platform
import time

import sympy as sp


AUDIT = Path(__file__).resolve().parents[1]
SOURCE = "https://doi.org/10.5281/zenodo.23118432"
OUT = AUDIT / "out"
START = time.perf_counter()
CHECKS = []


def zero(name, expression):
    values = list(expression) if isinstance(expression, sp.MatrixBase) else [expression]
    residuals = [sp.factor(sp.expand(value)) for value in values]
    assert all(value == 0 for value in residuals), (name, residuals)
    CHECKS.append({"name": name, "method": "exact symbolic identity", "residuals": [str(v) for v in residuals]})


def exact(name, condition, values):
    assert condition, (name, values)
    CHECKS.append({"name": name, "method": "exact integer/rational evaluation", "values": values})


def cross2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def skew(a):
    return sp.Matrix([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])


def rz(angle):
    return sp.Matrix([[sp.cos(angle), -sp.sin(angle), 0], [sp.sin(angle), sp.cos(angle), 0], [0, 0, 1]])


def ry(angle):
    # The sign convention required by the paper's stated X(theta,phi).
    return sp.Matrix([[sp.cos(angle), 0, -sp.sin(angle)], [0, 1, 0], [sp.sin(angle), 0, sp.cos(angle)]])


def trig_zero(name, expression):
    values = list(expression) if isinstance(expression, sp.MatrixBase) else [expression]
    residuals = [sp.trigsimp(value) for value in values]
    assert all(value == 0 for value in residuals), (name, residuals)
    CHECKS.append({"name": name, "method": "exact trigonometric identity", "residuals": [str(v) for v in residuals]})


# Source locations are citations; exact checks do not require source text.

# Lemma 2.1: the balancing identity is unrestricted; separation also needs positivity.
planar_vectors = [sp.Matrix(sp.symbols(f"w{i}x w{i}y")) for i in range(1, 4)]
planar_weights = [cross2(planar_vectors[(i + 1) % 3], planar_vectors[(i + 2) % 3]) for i in range(3)]
zero("Lemma 2.1 three-vector balancing", sum((planar_weights[i] * planar_vectors[i] for i in range(3)), sp.zeros(2, 1)))

# Lemma 5.1: every coordinate of the basis matrix remains symbolic.
epsilon, delta = sp.symbols("epsilon delta")
B = sp.Matrix(3, 3, sp.symbols("B11 B12 B13 B21 B22 B23 B31 B32 B33"))
nu_coordinates = [sp.Matrix(sp.symbols(f"alpha{i} xi{i} zeta{i}")) for i in range(1, 4)]
u_coordinates = sp.Matrix([1, epsilon, delta])
u_generic = B * u_coordinates
nu_generic = [B * coordinate for coordinate in nu_coordinates]
deltas = []
for k in range(3):
    i, j = (k + 1) % 3, (k + 2) % 3
    ai, xi, zi = nu_coordinates[i]
    aj, xj, zj = nu_coordinates[j]
    d = xi * zj - zi * xj - epsilon * (ai * zj - zi * aj) + delta * (ai * xj - xi * aj)
    deltas.append(d)
    zero(f"Lemma 5.1 coordinate determinant Delta{k + 1}", sp.Matrix.hstack(u_coordinates, nu_coordinates[i], nu_coordinates[j]).det() - d)
    zero(f"Lemma 5.1 arbitrary-basis determinant Delta{k + 1}", sp.Matrix.hstack(u_generic, nu_generic[i], nu_generic[j]).det() - B.det() * d)

nu_determinant = sp.Matrix.hstack(*nu_coordinates).det()
zero("Lemma 5.1 generic vector Cramer identity", sum((deltas[i] * nu_generic[i] for i in range(3)), sp.zeros(3, 1)) - nu_determinant * u_generic)
zero("Lemma 5.1 generic physical cross-product balancing", sum((deltas[i] * u_generic.cross(nu_generic[i]) for i in range(3)), sp.zeros(3, 1)))
F = sp.symbols("F", nonzero=True)
zero("Lemma 5.1 common factor rational identity", sum((deltas[i] / F * u_generic.cross(nu_generic[i]) for i in range(3)), sp.zeros(3, 1)))

# Lemma 2.2: the four moves, the Euler frame and the final alpha normal form.
I = sp.eye(3)
e3 = sp.Matrix([0, 0, 1])
projection = sp.Matrix([[1, 0, 0], [0, 1, 0]])
N = sp.diag(1, -1, -1)
Z = sp.Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
reflection = sp.diag(1, -1)
zero("Lemma 2.2 N proper orthogonal", N.T * N - I)
zero("Lemma 2.2 N determinant +1", N.det() - 1)
zero("Lemma 2.2 N planar reflection", projection * N - reflection * projection)
zero("Lemma 2.2 N negates view", N.T * e3 + e3)
zero("Lemma 2.2 Z proper orthogonal", Z.T * Z - I)
zero("Lemma 2.2 Z determinant +1", Z.det() - 1)
theta, phi, beta, beta1, beta2, alpha = sp.symbols("theta phi beta beta1 beta2 alpha", real=True)
frame = Z * ry(phi) * rz(-theta)
X = sp.Matrix([sp.cos(theta) * sp.sin(phi), sp.sin(theta) * sp.sin(phi), sp.cos(phi)])
trig_zero("Lemma 2.2 Euler frame view X", frame.T * e3 - X)
trig_zero("Lemma 2.2 Euler frame proper orthogonal", frame.T * frame - I)
trig_zero("Lemma 2.2 Euler frame determinant +1", frame.det() - 1)
trig_zero("Lemma 2.2 common Rz preserves view", rz(beta).T * e3 - e3)
planar_rotation = rz(beta)[:2, :2]
trig_zero("Lemma 2.2 common Rz planar isometry", projection * rz(beta) - planar_rotation * projection)
trig_zero("Lemma 2.2 alpha normal form", rz(beta1 - beta2) * Z * ry(phi) * rz(-theta) - Z * rz(beta1 - beta2) * ry(phi) * rz(-theta))
A = sp.Matrix(3, 3, sp.symbols("A11 A12 A13 A21 A22 A23 A31 A32 A33"))
G = sp.Matrix(3, 3, sp.symbols("G11 G12 G13 G21 G22 G23 G31 G32 G33"))
zero("Lemma 2.2 right-group action on view", (A * G).T * e3 - G.T * (A.T * e3))
zero("Lemma 2.2 common N action on both views", (N * A).T * e3 + A.T * e3)

# Every SO(3) stabilizer of e3 is a planar SO(2) matrix, hence Rz(beta).
v, w = sp.symbols("v w", real=True)
stabilizer = sp.Matrix([[v, -w, 0], [w, v, 0], [0, 0, 1]])
zero("Lemma 2.2 SO(2) stabilizer orthogonality", stabilizer.T * stabilizer - sp.diag(v * v + w * w, v * v + w * w, 1))
zero("Lemma 2.2 SO(2) stabilizer determinant", stabilizer.det() - (v * v + w * w))

# Lemma 3.1: exact Rodrigues identity and the coefficient norm of its error term.
eta = sp.Matrix(sp.symbols("eta1 eta2 eta3", real=True))
a = sp.Matrix(sp.symbols("a1 a2 a3", real=True))
U = sp.Matrix(sp.symbols("U1 U2 U3", real=True))
T = sp.symbols("T", real=True)
sin_s = 2 * T / (1 + T * T)
cos_s = (1 - T * T) / (1 + T * T)
E = cos_s * I + sin_s * skew(eta) + (1 - cos_s) * eta * eta.T
functional = -eta.dot(a.cross(U)) + T * (a.dot(U) - eta.dot(a) * eta.dot(U))
zero("Lemma 3.1 Rodrigues support difference", (1 + T * T) * ((a - E * a).dot(U) - sin_s * functional))
coefficient = eta.cross(U) + T * (U - eta * eta.dot(U))
zero("Lemma 3.1 linear functional coefficient", functional - a.dot(coefficient))
zero("Lemma 3.1 norm identity modulo eta norm", coefficient.dot(coefficient) - (1 + T * T) * (U.dot(U) - eta.dot(U) ** 2) - (eta.dot(eta) - 1) * (U.dot(U) + T * T * eta.dot(U) ** 2))

# Section 4.3: r parallel to a4 fixes a4 exactly, not merely to first order.
rx, ry_symbol, rz_symbol = sp.symbols("rx ry rz", real=True)
r = sp.Matrix([rx, ry_symbol, rz_symbol])
D = 1 + r.dot(r)
numerator = (1 - r.dot(r)) * I + 2 * r * r.T + 2 * skew(r)
zero("Gibbs matrix proper orthogonal", numerator.T * numerator - D * D * I)
zero("Gibbs matrix determinant", numerator.det() - D ** 3)
zero("Section 4.3 no-fit plane fixed a4", numerator.subs({rx: 0, ry_symbol: 0}) * e3 - D.subs({rx: 0, ry_symbol: 0}) * e3)
zero("Section 4.3 horizontal view orthogonal a4", sp.Matrix([sp.symbols("cx"), sp.symbols("cy"), 0]).dot(e3))

# A fully exact counterexample to the common-sign certificate criterion.
# The chosen basis obeys b = c x a4, as in the paper.
vertices = [sp.Matrix([sx, 2 * sy, 3 * sz]) for sx, sy, sz in product([-1, 1], repeat=3)]
Ain = I
Aout = sp.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
chart_rotation = Aout.T
zero("Counterexample Ain in SO(3)", Ain.T * Ain - I)
zero("Counterexample Aout in SO(3)", Aout.T * Aout - I)
zero("Counterexample Aout determinant +1", Aout.det() - 1)
zero("Counterexample chart rotation Gibbs r=(1,1,1)", numerator.subs({rx: 1, ry_symbol: 1, rz_symbol: 1}) - 4 * chart_rotation)
inner_images = [projection * Ain * p for p in vertices]
outer_images = [projection * Aout * p for p in vertices]
inner_extents = [(min(p[i] for p in inner_images), max(p[i] for p in inner_images)) for i in range(2)]
outer_extents = [(min(p[i] for p in outer_images), max(p[i] for p in outer_images)) for i in range(2)]
exact("Counterexample strict containment", all(outer_extents[i][0] < inner_extents[i][0] < inner_extents[i][1] < outer_extents[i][1] for i in range(2)), {"inner": [[str(x) for x in pair] for pair in inner_extents], "outer": [[str(x) for x in pair] for pair in outer_extents], "translation": [0, 0]})

clockwise = [sp.Matrix([1, 0]), sp.Matrix([0, -1]), sp.Matrix([-1, 1])]
negative_weights = [cross2(clockwise[(i + 1) % 3], clockwise[(i + 2) % 3]) for i in range(3)]
support_gaps = [max(p.dot(direction) for p in inner_images) - max(p.dot(direction) for p in outer_images) for direction in clockwise]
zero("Counterexample clockwise balancing", sum((negative_weights[i] * clockwise[i] for i in range(3)), sp.zeros(2, 1)))
exact("Counterexample common negative weights with positive sum", negative_weights == [-1, -1, -1] and support_gaps == [-1, -1, -2] and sum(negative_weights[i] * support_gaps[i] for i in range(3)) == 4, {"weights": [str(x) for x in negative_weights], "support_gaps": [str(x) for x in support_gaps], "weighted_sum": "4"})

c = sp.Matrix([1, 0, 0])
a4 = e3
b = c.cross(a4)
chart_u = c
counter_nus = [-e3, sp.Matrix([0, 1, 0]), sp.Matrix([0, -1, 1])]
counter_basis = sp.Matrix.hstack(c, b, a4)
counter_coordinates = [counter_basis.inv() * nu for nu in counter_nus]
counter_deltas = []
for i in range(3):
    ci, cj = counter_coordinates[(i + 1) % 3], counter_coordinates[(i + 2) % 3]
    counter_deltas.append(ci[1] * cj[2] - ci[2] * cj[1])
counter_directions = [chart_u.cross(nu) for nu in counter_nus]
plug_vertices = [sp.Matrix([1, 2, 3]), sp.Matrix([1, 2, 3]), sp.Matrix([-1, -2, 3])]
hole_vertices = [sp.Matrix([1, 2, 3]), sp.Matrix([1, 2, 3]), sp.Matrix([1, -2, -3])]
counter_values = [4 * (chart_rotation * plug_vertices[i] - hole_vertices[i]).dot(counter_directions[i]) for i in range(3)]
validity_values = [[(hole_vertices[i] - vertex).cross(chart_u).dot(counter_nus[i]) for vertex in vertices] for i in range(3)]
exact("Counterexample all validity inequalities", all(v >= 0 for row in validity_values for v in row), {"values": [[str(v) for v in row] for row in validity_values]})
zero("Counterexample Delta balancing in paper basis", sum((counter_deltas[i] * counter_directions[i] for i in range(3)), sp.zeros(3, 1)))
exact("Counterexample zoom conditions 2 and 3", counter_deltas == [-1, -1, -1] and counter_values == [-4, -4, -8] and sum(counter_deltas[i] * counter_values[i] for i in range(3)) == 16, {"Delta": [str(v) for v in counter_deltas], "g": [str(v) for v in counter_values], "weighted_sum": "16", "mF": "1", "family": "constant Phi=(epsilon,delta,r)=(0,0,1,1,1) on Y=[1,2]^5", "distance_variables": ["y0>0"]})

result = {
    "source": str(SOURCE),
    "source_sha256": 'f6988b1e63ba16a3f7493b2c0fb0a5451da355d454a6e568637b5c25a045c319',
    "python": platform.python_version(),
    "sympy": sp.__version__,
    "arithmetic": "exact symbolic polynomials, exact trigonometric identities, integers and rationals; no numerical sign test",
    "checks": CHECKS,
    "findings": [
        {"lemma": "2.1", "status": "PASS", "scope": "Balanced separation is correct when all weights are nonnegative and at least one positive weight has a nonzero direction. The three-vector identity itself does not enforce positivity."},
        {"lemma": "2.2", "status": "PASS_WITH_ROTATION_CONVENTION", "scope": "All four reduction moves and the Euler form hold for a finite proper rotation group. Ry(phi) uses [[cos,0,-sin],[0,1,0],[sin,0,cos]] to match the paper's X(theta,phi). The opposite active-rotation convention negates the horizontal components of X. Common proper N realizes simultaneous view negation."},
        {"lemma": "3.1", "status": "PASS_WITH_RANGE_REQUIREMENT", "scope": "Rodrigues expression and error coefficient norm verified exactly. The min-over-maximizers step uses sin(s)>0 and finite tan(s/2), thus requires the local range 0<s<pi."},
        {"lemma": "5.1", "status": "PASS", "scope": "Verified with all nine entries of the basis matrix symbolic. Common division preserves balance where the factor is nonzero; exact polynomial quotients extend the identity to the zero set."},
        {"lemma": "5.2", "status": "LOGICAL_GAP", "scope": "Common sign does not imply nonnegative weights. The criterion with common negative weights can accept a strict passage. The rectangular-box example refutes that generic criterion, not the snub-specific no-Rupert conclusion or the actual certificate data."},
        {"lemma": "4.3 W plane", "status": "PASS", "scope": "On delta=0 and r parallel to a4, both shadows have identical supports at plus/minus a4, so balanced weights (1,1) rule out all translated strict containments."},
    ],
    "elapsed_seconds": round(time.perf_counter() - START, 6),
}

report = r"""# 본문 보조정리의 정확한 대수 검증

검증 명령: `python code/lemmas_audit.py`

설치된 SymPy로 정수·유리수와 기호 다항식 항등식을 검증했습니다. 부호를 부동소수점 근삿값으로 판정한 항목은 없습니다. 검증 대상은 본문 보조정리이며, 보충자료의 전체 셀·분할 인증서는 이 스크립트의 대상이 아닙니다.

| 대상 | 판정 | 확인 사항 |
|---|---|---|
| Lemma 2.1 | 타당 | 균형식과 비음수 가중치로 평행이동 항이 소거되고 엄격부등식이 유지됩니다. |
| Lemma 2.2 | 회전 규약을 맞추면 타당 | 네 가지 변환, Euler 틀, 동시 시선 반전과 두 Dirichlet 셀로의 축소가 성립합니다. |
| Lemma 3.1 | 국소 각도 범위에서 타당 | Rodrigues 식과 오차항의 계수 노름이 정확히 일치합니다. `0<s<π` 범위가 필요합니다. |
| Lemma 5.1 | 타당 | 일반적인 3×3 기저행렬의 9개 성분을 모두 기호로 두어 균형식을 검증했습니다. |
| Lemma 5.2 | 증명에 부호 조건 누락 | “공통 부호”만으로 `λk≥0`가 따라오지 않습니다. 공통 음수이면 판정 기준 자체가 일반적으로 실패합니다. |
| §4.3의 W 평면 | 타당한 정확한 장애물 | `±a4` 방향의 지지함수가 각각 같으므로 어떤 평행이동으로도 엄격 포함할 수 없습니다. |

## Lemma 2.1: 균형 분리와 필요한 부호

엄격 포함이면 모든 `w≠0`에서

`hin(w)+t·w<hout(w)`

가 성립하고 `w=0`에서는 양변이 0입니다. 따라서 `λk≥0`, `Σλkwk=0`이고 적어도 하나의 `λk>0, wk≠0`이면 `Σλk(hin−hout)(wk)<0`가 됩니다. 본문의 가정 `Σλk(hin−hout)(wk)≥0`와 모순입니다. 이 증명은 정확합니다.

세 벡터의 `λ1=w2×w3`, `λ2=w3×w1`, `λ3=w1×w2`에 대해서는 두 좌표의 `Σλiwi`가 모두 기호적으로 0입니다. 다만 이 항등식은 가중치의 비음수성을 보장하지 않습니다. 모두 음수이면 모두에 −1을 곱해 비음수로 만들 수 있지만, 동시에 지지함수 합의 부호도 바뀝니다.

## Lemma 2.2: 회전군 축소와 키랄성

`Z`와 `N=diag(1,−1,−1)`은 모두 행렬식 +1의 직교행렬입니다. 네 단계의 근거는 다음과 같습니다.

본문의 `X(θ,φ)`와 일치하도록 `Ry(φ)=[[cosφ,0,−sinφ],[0,1,0],[sinφ,0,cosφ]]` 규약을 사용했습니다. 반대 부호의 능동 회전 규약을 쓰면 `Fᵀe3`의 수평 두 성분이 본문의 X와 반대가 됩니다. 본문은 Ry의 행렬을 제시하지 않으므로 명세·구현이 어느 규약인지 확인해야 합니다. 이는 매개변수 규약 문제이며 축소 원리의 실패를 의미하지 않습니다.

1. `A→Ag`에서 `gP=P`이므로 그림자가 그대로이며 시선은 `gᵀAᵀe3`가 됩니다.
2. 두 행렬에 공통으로 `N`을 왼쪽에서 곱하면 그림자는 같은 평면 반사 `H=diag(1,−1)`를 받고 평행이동은 `Ht`로 바뀝니다. `πN=Hπ`이고 `Nᵀe3=−e3`이므로 두 시선이 함께 반전됩니다. 실제 3차원 운동은 적절한 회전입니다.
3. 공통 왼쪽 인자 `Rz(β)`는 그림자 평면과 평행이동을 같은 각도로 회전시키며 시선을 바꾸지 않습니다.
4. `F=ZRy(φ)Rz(−θ)`는 `Fᵀe3=X(θ,φ)`를 만족합니다. `Aᵀe3=X`이면 `AFᵀ∈SO(3)`이고 `AFᵀe3=e3`입니다. `e3`를 고정하는 `SO(3)` 행렬은 평면의 `SO(2)` 회전 `Rz(β)`이므로 `A=Rz(β)F`입니다. 구면 극점에서도 이 안정자 논증이 성립합니다.

`K=G∪(−G)`는 유한 직교군입니다. 유한 궤도에서 `c*·kᵀuout`를 최대화하는 `k=sg`를 선택하면 `kᵀuout`이 `Cell2`에 들어갑니다. `s=−1`이면 공통 `N`을 적용한 뒤 바깥쪽에 `g`를 오른쪽에서 곱하면 됩니다. 이어 안쪽에만 `G`의 원소를 곱해 `Cell1`에 넣습니다. 마지막 공통 `Rz`로 바깥 회전각을 없애면 안쪽의 상대각 `α`만 남습니다. `Z`는 `Rz`와 가환하므로 본문의 최종 형식과 일치합니다. 반사 대칭이나 `−I∈G`를 요구하지 않으므로 키랄 입체에도 적용됩니다.

## Lemma 5.1: 일반 기저에서의 정확한 항등식

기저행렬을 `B=[c b a4]`, `U=(1,ε,δ)ᵀ`, `Vk=(αk,ξk,ζk)ᵀ`로 두었습니다. 모든 `Bij`를 서로 독립적인 기호로 두고 다음을 직접 전개해 0임을 확인했습니다.

`det(BU,BVi,BVj)−det(B)Δk=0`

`ΣΔk BVk−det(V1,V2,V3)BU=0`

`ΣΔk((BU)×(BVk))=0`

따라서 고정된 수치 기저에서만 성립하는 검사가 아닙니다. 실제 기저에서 `det(B)≠0`이면 본문의 행렬식 몫 표현도 유효합니다. 공통 인자로 나누면 그 인자가 0이 아닌 곳에서 균형식이 유지됩니다. 모든 가중치의 정확한 다항식 나눗셈이 가능하면 정역의 소거법칙으로 몫의 균형식이 0인 점까지 다항식 항등식으로 확장됩니다. 이때 인자의 부호는 균형식에 영향을 주지 않지만 분리부등식의 방향에는 영향을 줍니다.

## Lemma 5.2: 공통 음수 부호에 대한 정확한 반례

조건 (2)는 가중치에 “공통 부호”만 요구합니다. 증명의 “So λk(y)≥0”는 그 조건으로부터 도출되지 않습니다. 이 문제를 보여 주는 일반적인 정확한 예는 직육면체

`P=[−1,1]×[−2,2]×[−3,3]`

입니다. `Ain=I`, `Aout(x,y,z)=(y,z,x)`는 모두 적절한 회전이고 그림자는 각각

`Sin=[−1,1]×[−2,2]`, `Sout=[−2,2]×[−3,3]`

입니다. 평행이동 0으로 엄격 포함합니다. 시계방향 삼중항

`w1=(1,0), w2=(0,−1), w3=(−1,1)`

은 정확히 `λ=(−1,−1,−1)`, `Σλiwi=0`을 만족합니다. 지지함수 차이는 `(−1,−1,−2)`이므로 `Σλi(hin−hout)(wi)=4>0`입니다. 비음수 조건을 공통 부호 조건으로 바꿨다면 실제 통과를 잘못 배제하게 됩니다.

같은 예를 §5의 다항식 형식에도 직접 넣었습니다. `c=e1`, `a4=e3`, `b=c×a4=−e2`, `ε=δ=0`, `r=(1,1,1)`로 두면 `D=4`, `R(r)(x,y,z)=(z,x,y)`입니다. `u=e1` 방향의 구멍 그림자는 `(y,z)` 좌표에서 `[−2,2]×[−3,3]`, 플러그 그림자는 `[−1,1]×[−2,2]`로 엄격 포함합니다.

| k | νk | u×νk | 플러그 p | 구멍 q | Δk | gk |
|---|---|---|---|---|---:|---:|
| 1 | `(0,0,−1)` | `(0,1,0)` | `(1,2,3)` | `(1,2,3)` | −1 | −4 |
| 2 | `(0,1,0)` | `(0,0,1)` | `(1,2,3)` | `(1,2,3)` | −1 | −4 |
| 3 | `(0,−1,1)` | `(0,−1,−1)` | `(−1,−2,3)` | `(1,−2,−3)` | −1 | −8 |

8개 꼭짓점 모두에 대해 조건 (1)의 지지 꼭짓점 유효성을 정확히 검사했습니다. 조건 (2)의 공통 음수 부호가 성립하고 `ΣΔkgk=16>0`입니다. 상수 다항식 가족 `Φ=(0,0,1,1,1)`, `Y=[1,2]^5`, 거리 변수 `y0`, `mF=1`이면 나눗셈 조건도 자동으로 만족하고 모든 거리 변수가 양수입니다. 따라서 명시된 일반 인증 기준은 실제 엄격 통과를 받아들일 수 있습니다.

이 예는 직육면체에 대한 일반 기준의 실패입니다. 스너브 큐브의 통과를 구성한 것이 아니므로 스너브 큐브가 Rupert가 아니라는 결론 자체를 반박하지 않습니다. 실제 보충 인증서에서 모든 가중치를 비음수로 정규화했는지는 별도로 확인해야 합니다. 본문을 고치려면 조건 (2)에 `λk≥0`를 명시하고 조건 (3)을 같은 정규화된 가중치로 평가해야 합니다. 공통 음수 가중치를 −1로 바꾸는 경우 양수 합도 음수 합으로 바뀌므로 조건 (3)을 재검사해야 합니다.

## Lemma 3.1: Rodrigues 식과 오차항

`|η|=1`, `T=tan(s/2)`에 대해 `sin s=2T/(1+T²)`, `cos s=(1−T²)/(1+T²)`를 대입하면

`(a−Ea)·U=sin s[−η·(a×U)+T(a·U−(η·a)(η·U))]`

가 정확히 성립합니다. 대괄호의 `a`에 대한 계수는 `v=η×U+T(U−η(η·U))`이고

`|v|²=(1+T²)(|U|²−(η·U)²)≤(1+T²)|U|²`

입니다. `|Ain(x)p−Ain(c)p|≤ein|p|`와 Cauchy–Schwarz를 합치면 본문의 오차항 `ein|p||U|√(1+T²)`가 나옵니다. 이 계수 노름 항등식은 `|η|²−1`로 정확히 나누어 떨어지는 잔차까지 확인했습니다.

최댓값을 얻는 `q`를 골라 `Fu`의 최솟값으로 내려가는 단계에서는 `sin s>0`가 필요합니다. 또한 `T`가 유한해야 하므로 표준적인 국소 범위 `0<s<π`로 해석해야 합니다. 본문은 `0<s≤r0`만 쓰므로 보충 명세에서 `r0<π`를 확인하거나 그 조건을 본문에 명시하면 됩니다.

## §4.3: W 평면의 정확한 비통과 증명

`δ=0`이면 `u=c+εb`가 수평이므로 `a4⊥u`입니다. `r∥a4`인 Gibbs 회전은 `a4`를 정확히 고정합니다. 따라서 같은 입체의 구멍과 플러그는 그림자 평면의 `+a4`, `−a4`에서 지지함수가 각각 같습니다. 두 방향에 가중치 `(1,1)`을 주면 균형식과 지지함수 차이의 합 0이 성립하고 Lemma 2.1이 모든 평행이동에 대한 엄격 통과를 배제합니다. 수치 관찰 `k≈1.78`이나 2차 근사는 이 증명에 쓰지 않았습니다.
"""

OUT.mkdir(exist_ok=True)
(OUT / "lemmas.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
(OUT / "lemmas.md").write_text(report, encoding="utf-8")
print(json.dumps({"checks": len(CHECKS), "elapsed_seconds": result["elapsed_seconds"], "status": "exact checks passed; Lemma 5.2 common-sign gap documented", "outputs": ["out/lemmas.json", "out/lemmas.md"]}, ensure_ascii=False))
