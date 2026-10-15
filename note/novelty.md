# Audit of the touching-arc endpoint in the reviewed sources

Audit date: 2026-10-14 (Asia/Seoul). This note checks source statements and the scope of a possible novelty claim. It does not compute or prove the proposed endpoint formula.

## Verdict

The reviewed manuscript already proves the touching arc and characterizes its endpoint as a unique signed-area zero, with a rigorous numerical bracket of angular width `10^-14`. The compact expression `theta_* = arctan(2 - t)` is not explicitly stated in the manuscript text examined here. No explicit occurrence of that expression, its direct `Q(t)` coefficient encoding `(2, -1, 0)`, or its characteristic decimal values was found in the tracked text of the reviewed Tom 7 `noperts` snapshot. Relevant snub-cube Lean definitions and C++ source were also read to distinguish their exact model, approximate discussion of the tie, and certificate-domain constants.

This supports a statement limited to those reviewed versions. It does not establish that the formula has never been stated elsewhere, that it is absent from all of Tom 7's work, or that it is absent from unreviewed supplementary material. Nor does source absence validate the formula; its mathematical derivation must be checked separately.

Suggested wording, once the formula has been proved:

> The touching arc and a rigorous numerical enclosure of its endpoint are established in [L,W, Theorem 4.2]. We have not found the compact formula `theta_* = arctan(2 - t)` stated explicitly in the October 3, 2026 manuscript or in Tom 7's `noperts` repository at commit `0824ef598e8637481848092ff283d815cf0c04e8`; we give an algebraic derivation here.

If the bibliography label `[Tom7]` denotes the 2025 paper or the video rather than this repository, the above audit does not yet justify substituting that label for the repository citation.

## Reviewed versions and bibliographic information

- Manuscript: `upstream/snubs.txt`, the existing extraction from `upstream/snubs.pdf`. The PDF title is **The Snub Cube and the Snub Dodecahedron Are Not Rupert** (text lines 1–2). It is dated **October 3, 2026** (line 41). No author names occur in the extraction. The live Zenodo API record captured by the main audit in `historical Zenodo API snapshot` identifies exactly one creator entry with literal name **`L, W`**, publication date **2026-10-03**, and metadata title **The Snub Cube and Snub Dodecahedron Are Not Rupert**. The omission of the second “the” is a metadata/PDF title difference, not a different paper. The version DOI is **10.5281/zenodo.23118432** ([record API](https://zenodo.org/api/records/23118432), [DOI](https://doi.org/10.5281/zenodo.23118432)). The API's `files[0].checksum` is `md5:5fe9b7c898534f37d6b5d4a19897150a`; a lightweight `Get-FileHash -Algorithm MD5` check of the local PDF returned the same value. Thus the local PDF matches the record's published file. The registered creator string is verified; it does not establish two individual authors or expand their identities. The bibliography occupies text lines 399–416. Reference [11] names `SPEC.md` and `SPEC SC.md` as supplementary certificate specifications (lines 415–416); lines 344–345 say that specifications, producers, certificates or regeneration scripts, and both verifiers form the supplementary material.
- Repository: existing clone `upstream/noperts`, remote `https://github.com/tom7/noperts.git`. Audited HEAD: **`0824ef598e8637481848092ff283d815cf0c04e8`**. Both author and committer timestamps returned by `git log -1` are **`2026-10-14T04:01:45+00:00`**, with author **Tom 7** and subject **“Septopert rebuild: one thread per core by default; disk note for the new data”**. `git status --short` returned no entries during the audit. The repository README identifies itself as “Tom 7, 10 Oct 2026” at line 2; this is distinct from the HEAD timestamp.
- The README links to Tom 7's 2025 SIGBOVIK paper and subsequent video at lines 13–15. Those linked works were not reviewed by this subtask. Repository comments also refer to `nopert229/notes/SNUB-CUBE.md` and `nopert229/sccap.cc`, which are not tracked in this clone (`search/sc/sccapcert.cc:1–4`; `proof/SnubCube/Caps.lean:14–16`). Their contents are outside this audit's scope.

## What the manuscript states

| Location in `upstream/snubs.txt` | Source statement and its implication |
|---|---|
| 14–20 | The abstract states a touching arc of approximately 9.13 degrees and describes an exact algebraic chart over `Q(t)`. This establishes that the arc itself is prior work. |
| 187–194 | The manuscript defines `t^3 = t^2 + t + 1`, gives `1/t = t^2 - t - 1`, fixes the vertex convention, and sets `a4 = e3`, `c = (t^2 - t - 2, t^2 - t, 0)`, and `b = c × a4`. These are model/chart constants, not a stated endpoint formula. |
| 209–213 | The reflection relation explains why the special plug has the same shadow as the mirror image, with contact tests vanishing identically on the special surface. |
| 214–220 | Theorem 4.2 uses `epsilon = -tan(theta)` and `theta0 = arctan(0.16071324478583)`. It places `theta_*` in `[theta0, theta0 + 10^-14]`, reports `theta_* = 0.1593506261845...` radians, approximately `9.130118` degrees, proves containment for `0 <= theta <= theta_*`, and failure just past that endpoint. The exponent is **−14**, not −4. |
| 232–243 | Near the endpoint, the active image is that of vertex `v20` against the oriented edge `(v8, v3)`, with a corresponding test for `v17` against `(v21, v23)`. The proof encloses the derivative/signs and characterizes each threshold as a unique zero of a signed-area function. It takes the minimum of the two zeros and notes their equality by half-turn symmetry. This is an implicit exact characterization with a numerical enclosure; the text does not simplify that zero to `arctan(2 - t)`. |

Thus a new derivation should credit both the existence of the touching family and the rigorous endpoint characterization to the manuscript. The incremental contribution, if correct, is the simplified algebraic expression and its derivation.

## What the repository states and encodes

All repository line references below are relative to `upstream/noperts` and to the audited commit.

| Location | Source statement and interpretation |
|---|---|
| `snub-cube.md:53–60` | Gives the zero-margin view using `p = (-t^2, ±t, 1)` and describes a one-parameter exact tie with relative rotation `R_y(2 epsilon)`. The prose says clearance zero for `epsilon` up to about `0.165 rad`. It does not name an exact terminal angle. |
| `search/sc/sccapcert.cc:6–26` | Defines the exact model as the O-orbit of `(t, 1, -t^2)`, `x = (1, 0, t^2)`, and the charts. At the tie, the relative rotation is about `e_y` by `2 atan(mu)`; the comment says zero clearance for **`mu`** up to about `0.165`. This is a different parameter from the markdown's **angle** `epsilon`. The two approximate statements should not be treated as interchangeable exact endpoint bounds. |
| `proof/SnubCube/Caps.lean:18–21, 35–64` | Explicit exact frame constants include `L2 = (4 - 4t + 12t^2)/7`; the cap uses `mu0 = 1/8`, and the tube uses `tieMu0 = 1/3`. The seed is `(t, 1, -t^2)` in exact `TribQ` coefficients. These constants are model or certificate-domain parameters, rather than the last touching value. |
| `proof/SnubCube/Caps.lean:68–79` | The cap and tube claims exclude **Rupert poses**. They are not claims that every point in the tube is a touching passage. In particular, `tieMu0 = 1/3` should not be reported as the touching endpoint. |
| `proof/Pipeline/KCapTree.lean:635–645` | Defines the tube as a neighborhood of the special family using bounds on the view and Cayley-vector frame coordinates. Its upper `mu` bound is a certificate-domain bound. |
| `search/sc/sccapcert.cc:185, 223–228, 252–286` | The actual setup uses rational cap/tube bounds `1/8` and `1/3`, constructs the orbit of the exact seed, and computes the frame constant `L2` from edge-on supporting triangles. No named endpoint calculation appears in this setup. |
| `search/sc/sc_model.cc:147–184, 232–259, 277–280` | Constructs the normalized/shrunken model from the same exact seed and enumerates the zero-margin directions `e_a × p`. The `ustar` export is a **view direction**, not the angular length of the touching arc. |
| `search/sc/sc_model.txt:2–28, 92` | Stores the normalized vertices and `ustar = (0.2834724367428329595549687093248158376139, 0, 0.9589803843703376638300331996548725128970)`. The numeric coordinate constants should not be mistaken for a terminal angle. |
| `proof/SnubCube/ModelData.lean:13–49` | Stores the rational base point, slot elements, rational vertices, and normalization scale. It does not define an endpoint constant. |
| `proof/Pipeline/TribQ.lean:10–18, 110–146` | Defines the field representation as `a + b t + c t^2`. Thus an explicit reduced representation of `2 - t` would have coefficients `(2, -1, 0)`. The searched direct encodings were not found. This does not exclude expressions that would simplify to the same value only after further algebra. |

The markdown's verification reports at `snub-cube.md:119–130` are source reports of a successful Lean verification and reproduction. This subtask did not run or independently reproduce them.

## Search scope and negative evidence

`rg` was attempted first but is unavailable in the shell. The fallback was `git grep` over the complete tracked snapshot, not a file-extension-limited grep. `git ls-files` returned **503 tracked files**. `git grep -l -I -E '.' -- .` returned **500 nonempty text files**. The only three excluded files were `jcreed-sagemath.png`, `septopert-dry-salt-bed.jpg`, and `septopert.gif`. The tracked ASCII STL was included. Lean proofs and generated data, all `search/sc` C++, generic search support, other solids' code/data, Markdown, scripts, and configuration text were within the textual search scope.

Search families included:

- `arctan`, `atan`, `0.1593`, `0.1607`, `0.165`, `9.13`, `tie`, `wall`, `ustar`, `u*`, `endpoint`, and `closed form`;
- direct `2 - t`/`2 − t`, `2 - tribonacci`, and reversed `tribonacci - 2` spellings;
- explicit `TOf(2) - t`, `TSub(TOf(2), t)`, `TribQ.sub (TribQ.ofRat 2) TribQ.t`, and `ofRat 2 - TribQ.t` forms;
- direct coefficient forms corresponding to `(2, -1, 0)`, including `TQ{Q(2), Q(-1), Q(0)}`, `tq 2 (-1) 0`, Lean tuple notation, and list/array notation;
- the characteristic decimal prefixes `0.160713`, `0.159350`, and `9.130118`.

The searches for direct endpoint expressions/encodings and characteristic decimal prefixes returned no matches (`git grep` exit status 1). The broader searches located the cited tie prose, generic trigonometric uses, and unrelated hits for other solids. Source sections carrying the snub-cube constants and tie definitions were read directly, as listed above.

The full manuscript extraction was searched for trigonometric expressions, the candidate formula, closed-form terminology, and bibliographic metadata; the complete arc proof at lines 209–243 was read. It contains the decimal arctangent enclosure and implicit unique-zero argument described above.

Negative search results show absence of the searched explicit statements in this textual snapshot. They do not constitute a formal proof of absence of every algebraically equivalent expression, an audit of all history/branches, a literature-wide priority search, or an inspection of binary figures, external prototypes, videos, the 2025 Tom 7 paper, or the manuscript's supplementary files.

## Replayable search commands and recorded results

These PowerShell commands are recorded from the lightweight source audit above. Run from `upstream/noperts` with the audited commit checked out. They read files and Git metadata only. The broader searches produce irrelevant hits as well as the cited snub-cube sections; an absence claim uses the targeted searches and the direct source reading together.

```powershell
git rev-parse HEAD
git log -1 --format='%H%n%aI%n%cI%n%an%n%s'
git status --short
$tracked = @(git ls-files)
'TRACKED_COUNT ' + $tracked.Count
$textFiles = @(git grep -l -I -E '.' -- .)
'NONEMPTY_TEXT_COUNT ' + $textFiles.Count
Compare-Object -ReferenceObject $tracked -DifferenceObject $textFiles |
  Where-Object SideIndicator -eq '<=' | ForEach-Object InputObject
```

Recorded result: HEAD `0824ef598e8637481848092ff283d815cf0c04e8`; the timestamps and author listed above; no status entries; `TRACKED_COUNT 503`; `NONEMPTY_TEXT_COUNT 500`; only `jcreed-sagemath.png`, `septopert-dry-salt-bed.jpg`, and `septopert.gif` absent from the text-file list.

```powershell
git grep -n -I -i -E '2[[:space:]]*-[[:space:]]*t|arctan|atan|0\.1593|0\.1607|0\.165|9\.13|u\*|u_star|ustar|closed.form' -- .
git grep -n -I -i -E 'tie|wall|ustar|u\*|0\.165|endpoint' -- .
```

Recorded relevant results include `snub-cube.md:53–60`, `search/sc/sccapcert.cc:23–26`, `search/sc/sc_model.cc:237–280`, `search/sc/sc_model.txt:92`, and the Lean cap/chart sources listed above. Broad `u*` and substring searches also return multiplication expressions, ordinary words, generic trigonometry, and unrelated solids. The output was inspected and filtered for relevant source locations; broad-hit presence is not an endpoint-formula match.

```powershell
git grep -n -I -E '2[[:space:]]*[-−][[:space:]]*(tribonacci|t)([^A-Za-z0-9_]|$)|tribonacci[[:space:]]*[-−][[:space:]]*2|TOf\(2\)[[:space:]]*-[[:space:]]*t|TQ[[:space:]]*\{[[:space:]]*Q\(2\)[[:space:]]*,[[:space:]]*Q\(-1\)[[:space:]]*,[[:space:]]*Q\(0\)|tq[[:space:]]+2[[:space:]]+\(?-1\)?[[:space:]]+0|⟨[[:space:]]*2[[:space:]]*,[[:space:]]*[-−]1[[:space:]]*,[[:space:]]*0[[:space:]]*⟩|\[[[:space:]]*2[[:space:]]*,[[:space:]]*-1[[:space:]]*,[[:space:]]*0[[:space:]]*\]' -- .
'EXACT_PATTERN_EXIT ' + $LASTEXITCODE
git grep -n -I -E '0\.160713|0\.159350|9\.130118|0\.1593|0\.1607' -- .
'NUMERIC_PATTERN_EXIT ' + $LASTEXITCODE
git grep -n -I -E 'TSub\([[:space:]]*TOf\(2\)[[:space:]]*,[[:space:]]*t|TribQ\.sub[[:space:]]*\([[:space:]]*TribQ\.ofRat[[:space:]]+2[[:space:]]*\)[[:space:]]+TribQ\.t|ofRat[[:space:]]+2[[:space:]]*-[[:space:]]*TribQ\.t|TQ[[:space:]]*\{[[:space:]]*Q\(2\)[[:space:]]*,[[:space:]]*Q\(0\)[[:space:]]*-[[:space:]]*Q\(1\)|[⟨\[{][[:space:]]*\(?2[[:space:]]*,[[:space:]]*\(?-1\)?[[:space:]]*,[[:space:]]*0' -- .
'EXPANDED_PATTERN_EXIT ' + $LASTEXITCODE
```

Recorded result for each targeted search: no output and exit status 1 (`EXACT_PATTERN_EXIT 1`, `NUMERIC_PATTERN_EXIT 1`, `EXPANDED_PATTERN_EXIT 1`). Exit 1 here means no matching lines; the calls did not return a Git/regular-expression error. These syntactic patterns are intentionally documented so the limits of the negative evidence are reviewable.
