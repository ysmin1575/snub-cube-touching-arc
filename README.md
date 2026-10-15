# On the touching arc of the snub cube: a closed form and an independent check

**Sumin Yu** — Independent researcher, Republic of Korea  
**Version 1.0 · 2026-10-15**

This repository accompanies *On the touching arc of the snub cube: a closed form and an independent check*. It independently verifies the local closed-containment claims of Theorem 4.2 in L, W's preprint, [DOI 10.5281/zenodo.23118432](https://doi.org/10.5281/zenodo.23118432).

For the tribonacci constant `t`, the unique real root of `t^3-t^2-t-1=0`, the first endpoint is **theta_star = arctan(2-t)**. The certificate checks all 384 affine half-plane determinants in exact `Q(t)` arithmetic and certifies nonzero signs with Arb. An independent field implementation checks those coefficients. The note proves closed containment, with zero clearance, for the specified zero-translation family up to this endpoint and failure after it through `tan(theta)=1/3`.

The same family is identified in [Tom 7's noperts](https://github.com/tom7/noperts/tree/0824ef598e8637481848092ff283d815cf0c04e8). Its angular parameter ends at `theta_star` and its Cayley parameter at `2-t`. Both descriptive values `0.165` exceed the touching interval under their respective interpretations. This is a correction to the descriptive range, not a demonstrated failure of Tom 7's global certificates. Separate remarks exhibit the common-negative-weight issue in Lemma 5.2 and the exact `6/5` versus binary64 `1.2` discrepancy in Appendix A.4 (0.2 ulp).

Scope: a touching family and two textual remarks. This package does not verify every pose or translation, either global non-Rupert theorem, the full Lean/C++ certificates, or a global Nieuwland constant. The closed-form novelty statement concerns the supplied preprint version and pinned repository snapshot only. The historical search report is `note/novelty.md`; absence of a searched expression is not literature-wide priority evidence.

## Requirements and quick start

The versions actually installed and used/recorded for this release are Python **3.13.15**, NumPy **2.5.2**, mpmath **1.3.0**, SymPy **1.14.0**, python-flint **0.9.0**, and pypdf **6.20.1**. Use Python 3.13. NumPy is pinned for environment parity; the released mathematical checks do not import it. mpmath is a SymPy dependency. python-flint supplies the exact field and Arb arithmetic; pypdf checks the release PDF. SciPy, Git, TeX, upstream PDFs, network access and CI are not required for the installed offline self-test.

From the repository root, after selecting Python 3.13:

```text
python -B -m pip install -r requirements.txt
python -B self_test.py
```

The dependency installation requires network access unless wheels are already available. Both commands return exit code **0** on success. The self-test limits numerical threads to 1, executes all six checks sequentially in **one Python process**, verifies the inventory hashes and PDF, and compares eight recomputed JSON certificates with the reference data. It creates only `out/` results. No subprocess, worker, long search or overnight job is launched. Do not use Python `-O` or `PYTHONOPTIMIZE`; assertions form part of the checks and the runner rejects optimized execution.

Expected summary:

```text
HASH CHECK PASS: ... files
PRIVACY SCAN PASS: ... files; no candidates; no third-party PDFs
label_structure.py PASS (... s)
arc_certify.py PASS (... s)
arc_review.py PASS (... s)
lemmas_audit.py PASS (... s)
coverage_audit.py PASS (... s)
tom_arc_certify.py PASS (... s)
SELF-TEST PASS: 6 scripts, 8 certificates, 5 PDF pages; ... s
```

A failed assertion, missing dependency, unreadable input or hash mismatch returns a nonzero exit code and traceback. No failure is interpreted as a successful certificate. The two-minute target concerns the installed self-test, not package installation, TeX compilation or upstream source inspection. The release-preparation mathematical suite took about **8.3 seconds** (about 8.7 seconds including interpreter startup) on the author's Windows machine while another computation was running; the fresh-copy self-test measurement is recorded in the author's release report. Runtime is hardware dependent. `data/reproducibility.json` records the actual baseline run; `out/reproducibility.json` records each new invocation.

## Reproduce each claim

```text
python -B run_all.py
```

This runs the complete mathematical suite and returns 0 only if all six tasks complete. Each task can also be isolated using these exact commands; output setup and thread limits remain the same:

| Claim or check | Command from repository root | Expected output in `out/` |
|---|---|---|
| Exact vertices, consistent labels, symmetry group, reflection and rotation identities | `python -B run_all.py --only label_structure` | `labeling.json`, `structure.json`; verified labeling, 24 proper symmetries, zero improper symmetries |
| Endpoint `2-t`, angular enclosure, preprint bracket, 384 supporting inequalities, local derivative/triangle statements | `python -B run_all.py --only arc_certify` | `arc_certificate.json`, `arc_halfplanes.json`; `status=VERIFIED`, `s_star_Qt=["2","-1","0"]` |
| Independent field implementation, all 384 coefficients, hull geometry and named critical edge | `python -B run_all.py --only arc_review` | `arc_review.json`; all determinant coefficients equal; hull labels `[8,3,21,23,11,13,6,2]` |
| Lemma 5.2 negative weights, strict box containment, valid supports, balancing and witness 16; equal vertical support obstruction | `python -B run_all.py --only lemmas_audit` | `lemmas.json`, `lemmas.md`; exact checks pass and the generic criterion gap is documented |
| Appendix A.4 exact/binary64 bound, repaired chart bounds and conditional inverse-map identities | `python -B run_all.py --only coverage_audit` | `coverage.json`, `coverage.md`; exact margins and identities pass; global certificate remains outside scope |
| Tom-frame map, angular/Cayley relation, endpoint, containment/failure through `mu=1/3`, failures at both interpretations of 0.165 | `python -B run_all.py --only tom_arc_certify` | `tom_arc_certificate.json`; `status=VERIFIED`, 192 hole and 192 plug tests |
| All above plus shipped-file hashes, reference certificate equality and five-page PDF | `python -B self_test.py` | `self_test.json`; `status=PASS`, exit 0, runtime below 120 seconds |
| File inventory and privacy candidate scan | `python -B release_check.py` | `release_scan.json`; no candidates, exit 0 |

The endpoint interval has exact decimal rational bounds:

```text
0.1593506261845691446953081463393635724078
0.1593506261845691446953081463393635724079
```

The note's approximately `9.130118343142676` degrees and the two negative Tom-frame determinants, approximately `-0.0447731222` and `-0.0605847034`, appear in the Arb certificates. The A.4 discrepancy is exactly `1/22517998136852480`, or `1/5` of a binary64 ulp at 1.2. The self-test compares mathematical contents exactly; only elapsed time and interpreter/version observations are excluded from reference comparison. It recomputes assertions rather than trusting `status` labels alone.

The public scripts were copied and adapted from the original audit: output paths are portable, redundant source-document presence checks were removed, and a numerical discovery helper was omitted. Exact mathematical assertions were retained. `data/provenance.json` distinguishes original source/script hashes from release hashes. The unchanged historical source hashes identify the reviewed preprint; they do not claim the offline runner opens a bundled copy. Author-supplied source/document provenance is evidence for what was reviewed, not a substitute for checking the cited sources.

## Inspect source statements and repeat the limited novelty search

Read Theorem 4.2, Lemma 5.2, Section 5.3 and Appendix A.4 via the [preprint DOI](https://doi.org/10.5281/zenodo.23118432). The reviewed PDF SHA-256 is `b79ee5f3477e5e88b4967c50672e06293250f40875c7c633b74b1823f560786a`. No third-party PDF or extracted full text is redistributed. Download it separately if desired and use Python's `hashlib.file_digest` to compare it with `data/provenance.json`. The source version is dated 3 October 2026; do not silently substitute a later version.

After separately saving that PDF as `upstream/snubs.pdf`, this exact command prints its SHA-256 and returns 0 only if it matches the reviewed version:

```text
python -B -c "from pathlib import Path; import hashlib,json; p=Path('upstream/snubs.pdf'); h=hashlib.sha256(p.read_bytes()).hexdigest(); print(h); assert h==json.loads(Path('data/provenance.json').read_text(encoding='utf-8'))['reviewed_preprint']['pdf_sha256']"
```

Tom's description and exact frame are cited at these fixed source locations: [snub-cube.md, lines 53–60](https://github.com/tom7/noperts/blob/0824ef598e8637481848092ff283d815cf0c04e8/snub-cube.md#L53-L60), [sccapcert.cc, lines 6–26](https://github.com/tom7/noperts/blob/0824ef598e8637481848092ff283d815cf0c04e8/search/sc/sccapcert.cc#L6-L26), and [Caps.lean](https://github.com/tom7/noperts/blob/0824ef598e8637481848092ff283d815cf0c04e8/proof/SnubCube/Caps.lean#L35-L79).

Optional source replay requires Git and network access and is not part of the timed self-test. In a fresh directory, obtain the fixed snapshot:

```text
git init upstream/noperts
git -C upstream/noperts remote add origin https://github.com/tom7/noperts.git
git -C upstream/noperts fetch --depth 1 origin 0824ef598e8637481848092ff283d815cf0c04e8
git -C upstream/noperts checkout --detach FETCH_HEAD
```

Then change into `upstream/noperts` and execute the exact PowerShell search commands in `note/novelty.md`. Recorded inventory: 503 tracked files, 500 nonempty text files, three binary images excluded. Targeted endpoint-expression and decimal-prefix searches returned **Git exit 1 (no matching lines)**; an exit greater than 1 is a search error. Broad searches have matches and generally return 0. This replay checks syntactic negative evidence only. Original manual reading of the manuscript and relevant source sections remains a separate scholarly review step; neither the self-test nor grep proves semantic absence of every equivalent expression.

## Build the note PDF

The standalone `note/note.tex` is the authoritative typesetting source; update `note/note.md` consistently when editing. The release PDF was rebuilt with **Tectonic 0.17.0**, then inspected using PyMuPDF **1.26.7**. These are optional document tools, not mathematical dependencies. With Tectonic on PATH:

```text
tectonic --untrusted --outdir note note/note.tex
```

Alternatively, with a suitable XeLaTeX installation:

```text
xelatex -interaction=nonstopmode -halt-on-error -output-directory=note note/note.tex
xelatex -interaction=nonstopmode -halt-on-error -output-directory=note note/note.tex
```

Expected exit code: 0; `note/note.pdf` has five A4 pages. Compilation may populate the TeX cache or download packages depending on your installation. Inspect all pages after author edits. The reference suite checks the PDF's page count and required reproduction commands but does not certify all layout or prose changes automatically.

After intentional source, metadata or PDF edits, refresh the file inventory and rerun the self-test:

```text
python -B release_check.py --refresh
python -B self_test.py
```

If mathematical assertions change, first review them independently and rerun `run_all.py`; review and update the eight reference JSON certificates before accepting a new baseline. Never refresh a changed mathematical baseline solely to silence a failed comparison.

## Citation and archive

Sumin Yu (2026), *On the touching arc of the snub cube: a closed form and an independent check*, version 1.0, released 2026-10-15. Affiliation: Independent researcher, Republic of Korea. Structured citation metadata is in `CITATION.cff`.

Repository: [https://github.com/ysmin1575/snub-cube-touching-arc](https://github.com/ysmin1575/snub-cube-touching-arc). The release is archived on Zenodo; the version-specific code DOI is shown on the [GitHub release page](https://github.com/ysmin1575/snub-cube-touching-arc/releases/tag/v1.0) and Zenodo badge. The GitHub release archives the code through the Zenodo integration. The code DOI and any separate note DOI are not yet known; no DOI is assigned in this package. The cited L, W DOI identifies the reviewed preprint, not this release.

Use of AI. The author chose to audit the snub-cube preprint, directed the project and decided at each stage what to check and what to publish. AI systems were used extensively as research assistants: Claude (Anthropic) for planning, literature search, writing and running code, and drafting; GPT-6.1 Sol (OpenAI, via Codex) for independent derivations, the verification scripts and the first draft of this note. Every mathematical claim in this note is checked by the published scripts in exact or outward-rounded interval arithmetic, so the reader does not need to trust any AI system: running self_test.py reproduces all certified statements.

## Licenses and disclosure

The original code is licensed under the MIT License in `LICENSE`: Copyright (c) 2026 Sumin Yu. The note text and PDF (`note/note.md`, `note/note.tex`, and `note/note.pdf`) are licensed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/). These are separate licenses for the code and note; the note is not dual-licensed under MIT.

No upstream source checkout, third-party preprint, full extracted preprint text, external executable, API credential or local machine path is part of the public inventory. Public authorship and source-attribution names are intentional. The release scanner checks file text and PDF text/metadata for candidate paths, email addresses, common credential formats and private-key headers. It reports candidates without displaying secret values. Pattern scanning cannot guarantee exhaustive detection or establish copyright ownership. Dependencies and cited works retain their own licenses; dependency packages are not vendored here. `out/`, virtual environments, third-party replays and build caches are excluded from the repository inventory.

## File inventory and SHA-256

The table lists every shipped file except `README.md` and `SHA256SUMS`. `SHA256SUMS` additionally includes the complete README hash. A checksum manifest cannot include a stable hash of itself; its integrity is supplied by the immutable GitHub/Zenodo archive or an independently recorded archive hash. `python -B release_check.py --refresh` updates both the table and manifest together. Do not convert line endings on upload.

<!-- FILES START -->

| File | SHA-256 |
|---|---|
| `.gitattributes` | `705fd4d6451a31d36b3df7de96f83f30ac976c9b4a6d1e51671d8e2f33e2d0da` |
| `.gitignore` | `31b5d963ff5d50801b70e9d926f1f3d1f9dd51c8b6fa082ef5253d4a1a33714c` |
| `CITATION.cff` | `b32ee916c13594850ff6ff3c2e2372a1cf77c37ae0671899ce32468f8efa2c1b` |
| `code/arc_certify.py` | `63f40ff09cf6f14b2873c8396fee1de6729ffe8cbd2722dec0707b848f4f6ff4` |
| `code/arc_review.py` | `458113f97bec151fc664aa1c748be338af706f709d9a872ddb076101030bdd23` |
| `code/coverage_audit.py` | `fb35a73c7186c7507cc74f66bb168df4e080aecc7c29ed873619334f27b6325e` |
| `code/label_structure.py` | `e53ac828a8bef5f77e273ef724c5ce1bcf6cc59e8e80e4aefd4c2c4dc5e35adb` |
| `code/lemmas_audit.py` | `83ef9113b00adde8a6ece40ab071aa0a339e78ada8ef7193d256d08cf904937e` |
| `code/tom_arc_certify.py` | `2ab58bf27df7eae9a75665bd2692c1a64f3bd9405d72ed4c4eaa50425b0ce9ef` |
| `data/arc_certificate.json` | `8a26c9f9e203fb4e1721915d6cf81c8e63a8bb6add236682b82f8beb7d1d53fa` |
| `data/arc_halfplanes.json` | `d596d2f5921213f08f54ee9578129b102900bae61949bb7f4432cc647b517eb7` |
| `data/arc_review.json` | `1203c0700b31c642d9a9c5b9257176a40781f964b8a307cef931849524c0b1f6` |
| `data/coverage.json` | `966606443c0fb0fc000bf52367c47a5027a6cf434599a75b1ac3e674630d79b9` |
| `data/labeling.json` | `24a161a8d93bb15b2240dbadf52ed82c99e4989fbfc57a93c5c7c9a68de441d3` |
| `data/lemmas.json` | `a69cb30ece948f395fd9665e0ec688a1f313aef46582b7a4df16ab20fc8a7569` |
| `data/provenance.json` | `4b36145659deee7bfdebe859f4cbd81aae9542a8e2ea27d20cdd5e5fb9c5be51` |
| `data/reproducibility.json` | `6a8fe6a68715b72c564ba036baf5e873be563db60ff5730806ba34a2df1acfab` |
| `data/structure.json` | `278092f54ef3e0b36969f429870b2b109abdc46e2827c1b990122dad5fe32096` |
| `data/tom_arc_certificate.json` | `27f40dace966f2c96a401d8f38b492dbb318b74d8ab69c16dbd56d068ccf7b9a` |
| `LICENSE` | `cee17c4a566cb6671cbed70397c7f596b36db8a38dc9664fa8aa92d64d93ab9f` |
| `note/note.md` | `797c54af0267336f6b9a21eb37c36180f74f2677d7218228ca42d1748258b947` |
| `note/note.pdf` | `db9af9a3fb91b2f5ff63340a26791d2265a27c1d4460fad0dd79bc8d1c98f403` |
| `note/note.tex` | `506768264b06a26e3dcd8a7efbfa87f3f1dcec79d9e13158e260391266a97b2b` |
| `note/novelty.md` | `28a4723bc5a598880e5c976d83624ac1390cd04f07038abde85166c112efb2ee` |
| `release_check.py` | `7e8d9ae3313ff0393f4d1cf208faf1e07ec7ce6db4f4eb57d784dcad2b681c7c` |
| `requirements.txt` | `5ea23b28cd67c6e10d21b06e0695756c2ae58f155a6e09bafa42f033a0e747f1` |
| `run_all.py` | `82b89e0a031ebd07e17b4ae4f6d83ba63f5136c0f34117d86c352f1133ad2d83` |
| `self_test.py` | `6df957f9bfb8af4c5146a85c236f258d3bc8c4b09010032aa79ef9fdf11db963` |

<!-- FILES END -->
