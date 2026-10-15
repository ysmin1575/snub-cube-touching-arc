# On the touching arc of the snub cube: a closed form and an independent check

**Sumin Yu — Independent researcher, Republic of Korea**  
**Version 1.0 · 2026-10-15**  
Note text and PDF: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Abstract

The snub cube admits an initial arc of poses in which one congruent shadow is contained in another with boundary contact. L, W established this arc and enclosed its endpoint numerically in Theorem 4.2 of their preprint. We give the compact expression \(\theta_*=\arctan(2-t)\), where \(t\) is the tribonacci constant, and an independent certificate using exact arithmetic in \(\mathbb Q(t)\) and Arb interval arithmetic. A change of coordinates identifies this arc with the tie described in Tom 7's repository; its angular parameter has the same endpoint, whereas its Cayley parameter ends at \(2-t\). Thus the repository's approximate value \(0.165\) overstates the initial touching range under either interpretation. Two separate remarks concern a weight-sign inference and a binary64 bound in the L, W text. The scope is this arc and these textual points; no global non-Rupert certificate is reproduced here.

## 1. Setup and statement

Let \(t\in(1,2)\) be the unique real root of \(t^3-t^2-t-1=0\), and put \(q=1/t=t^2-t-1\). Let \(V\) comprise the even permutations of \((\pm1,\pm q,\pm t)\) with an even number of plus signs and the odd permutations with an odd number of plus signs. The 24-point convex hull \(P=\operatorname{conv}V\) is one chirality of the snub cube. A uniform scaling or a simultaneous orthogonal change of coordinates does not affect the containment question.

Write

\[
a_4=e_z,\quad C=(q-1,q+1,0),\quad B=C\times a_4=(q+1,1-q,0),
\]

\[
u_\theta=\widehat C\cos\theta-\widehat B\sin\theta,\qquad
K_n=I-\frac{2nn^{\mathsf T}}{n\cdot n},\qquad K_0=K_C.
\]

Here a hat denotes normalization, and \(\pi_u\) denotes orthogonal projection onto \(u^\perp\). Rotations are active, right-handed rotations. In particular, the \(xz\)-block of \(R_y(\alpha)\) is \(\left(\begin{smallmatrix}\cos\alpha&\sin\alpha\\-\sin\alpha&\cos\alpha\end{smallmatrix}\right)\). We ask for **closed containment at zero translation**, not containment in the interior. For unit \(u\), define the clearance of a fixed pose by \(\min_{d\perp u,\ |d|=1}(h_{\rm hole}(d)-h_{\rm plug}(d))\). Thus closed containment gives nonnegative clearance, and contact gives zero clearance.

**Theorem.** For \(0\leq\theta\leq\arctan(1/3)\),

\[
\pi_{u_\theta}(K_0P)\subseteq\pi_{u_\theta}(P)
\quad\Longleftrightarrow\quad
0\leq\theta\leq\boxed{\arctan(2-t)}.
\]

Equivalently the plug can be \(R_z(2\theta)P\), with no translation. The clearance is exactly zero throughout this initial arc. Every angle immediately beyond its endpoint, up to \(\arctan(1/3)\), fails closed containment for this pose. The angular domain is explicit because the family is periodic; this is not an assertion that containment occurs only once among all real angles.

## 2. Proof and certificate

Set \(s=\tan\theta\). Since \(|B|=|C|\), coordinates for the screen, after a positive horizontal rescaling, are

\[
Q_s(v)=(B\cdot v+sC\cdot v,v_z),\qquad
Q'_s(v)=(B\cdot v-sC\cdot v,v_z).
\]

The latter are the coordinates of \(K_0v\). The rescaling is invertible and hence preserves containment and the signs of oriented areas. In a labeling consistent with every named vertex constraint in [LW], the outer boundary has cyclic order

\[
(8,3,21,23,11,13,6,2).
\]

This order is a candidate used by the certificate, not a numerically trusted hull. Its eight edges are tested against all 24 outer vertices and all 24 reflected plug vertices. Each of these 384 determinants has the affine form \(A+sD\) over \(\mathbb Q(t)\). Nonnegative values at \(s=0\) and \(s=2-t\) establish every inequality throughout the interval by linear interpolation. Exact zero tests are performed in \(\mathbb Q[t]/(t^3-t^2-t-1)\); nonzero signs are enclosed using Arb. Distinctness of the boundary vertices, noncollapsing edges, and strictly positive cyclic turns are also checked throughout the interval, establishing that the supporting polygon is the actual convex hull. A second implementation of the field arithmetic recomputes all 384 coefficients.

The first failing test is the reflected vertex \(p=(-t,q,-1)\), labeled 20, against the oriented edge

\[
a=(-t,-1,-q),\qquad b=(-1,-q,-t),
\]

labeled \(8\to3\). Direct reduction of its determinant gives

\[
\det\bigl(Q_s(b)-Q_s(a),Q'_s(p)-Q_s(a)\bigr)
=2(t-1)^2(2-t-s).
\]

The positive prefactor \(2(t-1)^2\) makes the root exactly \(s_*=2-t\). The half-turn \(\operatorname{diag}(-1,-1,1)\) gives the identical obstruction at reflected vertex 17 against edge \(21\to23\). The first edge remains supporting for **every** \(s\in[0,1/3]\), as its outer-vertex affine inequalities are nonnegative at both endpoints. Its plug determinant is strictly negative for \(s\in(2-t,1/3]\), proving the converse.

Finally \(K_{u_\theta}K_0=R_z(2\theta)\), and \(\pi_{u_\theta}K_{u_\theta}=\pi_{u_\theta}\). This proves the rotation equivalence. The plug rotation preserves \(v_z\), and \(a_4\perp u_\theta\), so both vertical support values are equal; the contained shadow always touches. This observation also rules out strict containment by translating this particular pair of shadows, since translation cannot reduce their equal vertical widths.

At 384-bit working precision the certified enclosure is

\[
\begin{aligned}
\theta_*\in[\,&0.1593506261845691446953081463393635724078,\\
&0.1593506261845691446953081463393635724079\,]\,\mathrm{rad}.
\end{aligned}
\]

It is approximately \(9.13011834314267604063121785735^\circ\). The decimal interval has exact rational endpoints and width \(10^{-40}\); it is not a confidence interval or a sample-based estimate.

## 3. Comparison with the two sources

In [LW, Theorem 4.2], \(\theta_0=\arctan(0.16071324478583)\) and \(\theta_*\in[\theta_0,\theta_0+10^{-14}]\). Our check gives \(\theta_*-\theta_0\approx8.64418\times10^{-15}\), within that bracket. The preprint already establishes the arc and characterizes its endpoint by the first signed-area zero. Our contribution is its simplification, an independent containment certificate, and the explicit comparison below. To our knowledge, the closed form is not stated in the supplied October 3 version of [LW] or the inspected repository version of [Tom7]. This is a claim about those versions, not a literature-wide priority claim.

Tom 7's exact frame is

\[
x=(1,0,t^2),\quad e_1=(-t^2,0,1),\quad Z=e_y.
\]

The wall in question is the \(XZ\)-plane of outer views. On side 1 of the cap chart, the tie has

\[
u_T=x-\mu e_1,\qquad w=\mu e_y,\qquad
\phi=\arctan\mu,\qquad R_T=R_y(2\phi).
\]

Because \(x\perp e_1\) and \(|x|=|e_1|\), \(\phi\) is the angle past \(u^*=x/|x|\), namely the angle denoted \(\varepsilon\) in the markdown explanation. With \(H_n=2nn^{\mathsf T}/|n|^2-I=-K_n\), the half-turn convention used by the repository, \(H_{u_T}H_x=K_{u_T}K_x=R_y(2\phi)\). Individually \(H_u\) negates a screen projection; \(K_u\) fixes it. Their products coincide. The inner viewing direction is \(R_T^{\mathsf T}u_T\), explaining the opposite-angle inner view in the C++ comment.

The proper signed permutation \(S(X,Y,Z)=(-X,Z,Y)\), followed by division by \(t\), maps Tom's solid to \(P\). If \(k=1-q>0\), then

\[
Sx=C/k,\qquad Se_1=B/k,\qquad Se_y=a_4.
\]

Consequently \(\phi=\varepsilon_{\rm Tom}=\theta\), \(\mu=\tan\theta\), and the L, W chart variable is \(\varepsilon_{\rm LW}=-\mu\). These are the **same family**. Directly in Tom's unnormalized screen coordinates \(((e_1+\mu x)\cdot v,v_Y)\), the critical edge vertices and plug vertex are

\[
(t^2,-1,-t),\quad(t,-t^2,-1),\quad(t^2,-t,1),
\]

and their determinant is

\[
2(t-1)-2t(t+1)\mu=2t(t+1)(2-t-\mu).
\]

Thus \(\mu_*=2-t=0.160713244785838867448\ldots\) and \(\varepsilon_{\rm Tom,*}=\theta_*\); the relative rotation ends at \(2\theta_*=0.318701252369138\ldots\) radians. The markdown states about \(0.165\) **radians** (`snub-cube.md:57–60`), while the C++ comment states \(\mu\) up to about \(0.165\) (`search/sc/sccapcert.cc:23–26`). The new certificate proves failure at both \(\mu=0.165\) and \(\varepsilon_{\rm Tom}=0.165\) radians: the critical signed determinants are respectively enclosed near \(-0.0447731222\) and \(-0.0605847034\). Hence changing between angle and Cayley parameter does not resolve either approximate bound. This corrects a descriptive range, without showing any failure of the global cap/tube exclusion certificates. Their \(\mu\leq1/3\) domain need not consist entirely of touching poses.

All 503 tracked files of [Tom7] were enumerated and all 500 text files searched, including Lean, C++, generated data and Markdown. The three binary images were excluded from textual search. The scan did not cover repository history, external `nopert229` prototypes, linked videos/papers, or unprovided [LW] supplements; `novelty.md` records the evidence and limits.

## 4. Two remarks on [LW]

**Weight signs in Lemma 5.2.** A common sign does not imply that the weights are nonnegative. Flipping three negative weights also flips their weighted support gap, so it requires rechecking the strict inequality. An explicit counterexample to the general algebraic criterion is

\[
P_b=[-1,1]\times[-2,2]\times[-3,3],\quad u=e_1,\quad
R(x,y,z)=(z,x,y),\quad r=(1,1,1),\quad D=4.
\]

In the \((y,z)\) screen the hole is \([-2,2]\times[-3,3]\) and the plug is \([-1,1]\times[-2,2]\), strictly contained at translation zero. Choose \(c=e_1,a_4=e_3,b=-e_2\), \(\varepsilon=\delta=0\), and these tests:

| \(k\) | \(\nu_k\) | \(u\times\nu_k\) | plug vertex \(p_k\) | hole support vertex \(q_k\) | \(\Delta_k\) | \(g_k\) |
|---|---|---|---|---|---:|---:|
| 1 | \((0,0,-1)\) | \((0,1,0)\) | \((1,2,3)\) | \((1,2,3)\) | −1 | −4 |
| 2 | \((0,1,0)\) | \((0,0,1)\) | \((1,2,3)\) | \((1,2,3)\) | −1 | −4 |
| 3 | \((0,-1,1)\) | \((0,-1,-1)\) | \((-1,-2,3)\) | \((1,-2,-3)\) | −1 | −8 |

Here \(g_k=\langle N(r)p_k-Dq_k,u\times\nu_k\rangle\), \(\sum\Delta_k(u\times\nu_k)=0\), and \(\sum\Delta_kg_k=16>0\). Each hole support test is valid. The constant chart \(\Phi=(0,0,1,1,1)\) on \([1,2]^5\), with distance variable \(y_0\) and \(m_F=1\), also meets the divisibility conditions. Thus the common-sign inference and its corresponding general criterion fail. This box is not a snub cube and is not a counterexample to either snub solid's non-Rupert conclusion. The needed correction is to require nonnegative weights and evaluate condition (3) using those same weights. Whether the actual [LW] certificates already do so is outside this check.

**Appendix A.4 and binary64 constants.** Section 5.3 interprets decimal constants as binary64 values. However

\[
\operatorname{fl}(1.2)=\frac{5404319552844595}{4503599627370496}<\frac65,
\quad \frac65-\operatorname{fl}(1.2)=\frac1{22517998136852480}.
\]

For the special pair with \(s_2=1/5\), \(f_2=1,d=-1\) gives \(|f_2-ds_2|=6/5\), violating the displayed bound if read as \(\leq\operatorname{fl}(1.2)\). For example \(\mu=1/100\), \(\varepsilon=\delta=-\mu\), \(r=(0,\mu,\mu)\) is a corresponding CORNER point. The exact fraction \(6/5\) repairs this wording. It is already admitted by the relevant exact chart bound, and the enlarged Shigh bound \(5/4\) is larger. This discrepancy is 0.2 ulp at 1.2; it does not introduce a new hole in the continuum cover checked by our audit.

## 5. Reproducibility, scope and author statement

The public repository contains `code/` (six mathematical scripts), `data/` (reference JSON certificates and source provenance), and `note/` (this note and the source-review record). From the repository root run:

```text
python -B self_test.py
python -B run_all.py
```

Both commands run the checks sequentially in one Python process with numerical thread limits set to 1. The self-test checks SHA-256 integrity, reruns the mathematics and compares eight JSON certificates, and checks the five-page PDF. Generated results go to `out/`. No upstream PDF, clone, Lean build, C++ search or network access is needed for these checks. Historical input hashes identify the reviewed sources; source-document inspection is a separate procedure in the README and `novelty.md`.

The recorded environment is Python 3.13.15, NumPy 2.5.2, mpmath 1.3.0, SymPy 1.14.0, python-flint 0.9.0, and pypdf 6.20.1. NumPy is recorded for environment parity; the certificate path does not use it. mpmath is a SymPy dependency; pypdf is used by the release self-test. Script hashes, measured runtime, and exact claim-to-command mapping are in the repository README, `SHA256SUMS`, and `data/reproducibility.json`.

The reviewed preprint is DOI `10.5281/zenodo.23118432`. The noperts snapshot is commit `0824ef598e8637481848092ff283d815cf0c04e8`. Its recorded source search covers that snapshot only, not repository history or a literature-wide priority claim. This verification concerns Theorem 4.2 and the local family through \(\mu=1/3\), plus the two textual remarks. It does not certify all poses/translations, either global non-Rupert result, or a global Nieuwland constant.

Data and code availability. The version 1.0 release is available at [https://github.com/ysmin1575/snub-cube-touching-arc](https://github.com/ysmin1575/snub-cube-touching-arc) and archived on Zenodo; its version-specific code DOI is `CODE_DOI_PENDING` (to be filled in after the release).

Use of AI. The author chose to audit the snub-cube preprint, directed the project and decided at each stage what to check and what to publish. AI systems were used extensively as research assistants: Claude (Anthropic) for planning, literature search, writing and running code, and drafting; GPT-6.1 Sol (OpenAI, via Codex) for independent derivations, the verification scripts and the first draft of this note. Every mathematical claim in this note is checked by the published scripts in exact or outward-rounded interval arithmetic, so the reader does not need to trust any AI system: running self_test.py reproduces all certified statements.

**References**

- **[LW]** L, W (creator name as registered in Zenodo). *The Snub Cube and Snub Dodecahedron Are Not Rupert*. Preprint dated 3 October 2026; the PDF title includes “the” before “Snub Dodecahedron.” [DOI: 10.5281/zenodo.23118432](https://doi.org/10.5281/zenodo.23118432), especially Theorem 4.2, Lemma 5.2 and Appendix A.4.
- **[Tom7]** Tom 7. *noperts*, snub-cube description and sources, inspected at commit `0824ef598e8637481848092ff283d815cf0c04e8`, 14 October 2026; the snub-cube status entry reports verification on 13 October. [Pinned description](https://github.com/tom7/noperts/blob/0824ef598e8637481848092ff283d815cf0c04e8/snub-cube.md#L53-L60), [C++ chart description](https://github.com/tom7/noperts/blob/0824ef598e8637481848092ff283d815cf0c04e8/search/sc/sccapcert.cc#L6-L26), [exact Lean frame](https://github.com/tom7/noperts/blob/0824ef598e8637481848092ff283d815cf0c04e8/proof/SnubCube/Caps.lean#L35-L64).
