"""Offline self-test: integrity, six mathematical checks, PDF/source consistency."""
from fractions import Fraction
import json
from pathlib import Path
import time
import run_all

ROOT = Path(__file__).resolve().parent


def semantic(value):
    # Runtime and interpreter metadata are observations, not certificate claims.
    ignored = {'elapsed_seconds', 'runtime_seconds', 'python', 'sympy'}
    if isinstance(value, dict):
        return {k: semantic(v) for k, v in value.items() if k not in ignored}
    if isinstance(value, list):
        return [semantic(v) for v in value]
    return value


def main():
    if not __debug__:
        raise RuntimeError('Do not use -O: mathematical assertions must execute.')
    start = time.monotonic()
    from release_check import check_hashes, scan
    check_hashes()
    scan()
    run_all.main()
    names = ('labeling', 'structure', 'arc_certificate', 'arc_halfplanes', 'arc_review', 'lemmas', 'coverage', 'tom_arc_certificate')
    for name in names:
        actual = json.loads((ROOT / 'out' / (name + '.json')).read_text(encoding='utf-8'))
        expected = json.loads((ROOT / 'data' / (name + '.json')).read_text(encoding='utf-8'))
        assert semantic(actual) == semantic(expected), name + ' differs from reference certificate'
    arc = json.loads((ROOT / 'out/arc_certificate.json').read_text())
    tom = json.loads((ROOT / 'out/tom_arc_certificate.json').read_text())
    assert arc['s_star_Qt'] == tom['mu_star_exact'] == ['2', '-1', '0']
    assert tom['all_affine_endpoint_checks'] == {'hole': 192, 'plug': 192}
    assert arc['report_enclosure'] == tom['eps_star_enclosure']
    discrepancy = Fraction(6, 5) - Fraction.from_float(1.2)
    assert discrepancy == Fraction(1, 22517998136852480)
    assert discrepancy / Fraction(1, 2**52) == Fraction(1, 5)
    from pypdf import PdfReader
    pdf = PdfReader(ROOT / 'note/note.pdf')
    text = '\n'.join(p.extract_text() for p in pdf.pages)
    assert len(pdf.pages) == 5
    assert '23118432' in text
    assert 'python -B self_test.py' in text and 'python -B run_all.py' in text
    # Verify the reproduction commands in both published note sources.
    for name in ('note.md', 'note.tex'):
        source = (ROOT / 'note' / name).read_text(encoding='utf-8')
        assert 'python -B self_test.py' in source and 'python -B run_all.py' in source
    elapsed = time.monotonic() - start
    assert elapsed < 120, 'Self-test exceeded the two-minute target'
    result = {'status': 'PASS', 'seconds': elapsed, 'mathematical_scripts': 6,
              'matched_json_certificates': 8, 'pdf_pages': len(pdf.pages), 'exit_code': 0}
    (ROOT / 'out/self_test.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print('SELF-TEST PASS: 6 scripts, 8 certificates, 5 PDF pages; ' + format(elapsed, '.3f') + ' s')


if __name__ == '__main__':
    main()
