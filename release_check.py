"""Refresh/verify the publication inventory; scan for accidental disclosure."""
import argparse
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent


def files():
    return sorted(p for p in ROOT.rglob('*') if p.is_file() and not any(
        part in {'out', '.git', '.venv', '__pycache__', 'upstream'} for part in p.relative_to(ROOT).parts))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def refresh():
    inventory = [p for p in files() if p.name not in {'README.md', 'SHA256SUMS'}]
    rows = '\n'.join('| `' + p.relative_to(ROOT).as_posix() + '` | `' + sha(p) + '` |' for p in inventory)
    readme = ROOT / 'README.md'
    text = readme.read_text(encoding='utf-8')
    before, tail = text.split('<!-- FILES START -->', 1)
    _, after = tail.split('<!-- FILES END -->', 1)
    readme.write_text(before + '<!-- FILES START -->\n\n| File | SHA-256 |\n|---|---|\n' + rows + '\n\n<!-- FILES END -->' + after, encoding='utf-8')
    (ROOT / 'SHA256SUMS').write_text(''.join(sha(p) + '  ' + p.relative_to(ROOT).as_posix() + '\n' for p in files() if p.name != 'SHA256SUMS'), encoding='utf-8')
    print('MANIFEST REFRESHED')


def check_hashes():
    expected = {}
    for line in (ROOT / 'SHA256SUMS').read_text(encoding='utf-8').splitlines():
        digest, name = line.split('  ', 1)
        assert not Path(name).is_absolute() and '..' not in Path(name).parts
        expected[name] = digest
    actual = {p.relative_to(ROOT).as_posix(): sha(p) for p in files() if p.name != 'SHA256SUMS'}
    assert actual == expected, 'Inventory or SHA-256 mismatch; after deliberate edits run --refresh'
    print('HASH CHECK PASS: ' + str(len(actual)) + ' files')


def scan():
    from pypdf import PdfReader
    # Candidate detection complements content review, not a proof of ownership.
    patterns = {
        'windows_absolute_path': re.compile(r'(?<![A-Za-z0-9])[A-Za-z]:[\\/]'),
        'unix_home_path': re.compile(r'/(?:home|Users)/[A-Za-z0-9_.-]+/'),
        'email': re.compile(r'[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}'),
        'github_token': re.compile('gh[pousr]_' + r'[A-Za-z0-9]{20,}'),
        'api_key': re.compile('sk-' + r'[A-Za-z0-9_-]{20,}'),
        'aws_key': re.compile('AKIA' + r'[A-Z0-9]{16}'),
        'private_key': re.compile('-----BEGIN ' + '(?:RSA |EC |OPENSSH )?PRIVATE KEY-----'),
    }
    findings = []
    scanned = files()
    for p in scanned:
        if p.suffix == '.pdf':
            doc = PdfReader(p)
            text = '\n'.join(page.extract_text() for page in doc.pages) + json.dumps(dict(doc.metadata or {}), ensure_ascii=False)
        else:
            text = p.read_text(encoding='utf-8')
        for kind, pattern in patterns.items():
            if pattern.search(text):
                findings.append({'file': p.relative_to(ROOT).as_posix(), 'kind': kind})
    assert not findings, 'Disclosure candidates: ' + json.dumps(findings)
    assert [p.relative_to(ROOT).as_posix() for p in scanned if p.suffix == '.pdf'] == ['note/note.pdf']
    assert not any(p.suffix.lower() in {'.exe', '.dll', '.zip', '.pem', '.key'} for p in scanned)
    result = {'status': 'PASS', 'files_scanned': len(scanned), 'findings': findings,
              'intentional_public_names': ['Sumin Yu', 'L, W and Tom 7 (source citations)'],
              'third_party_pdf_count': 0,
              'limitations': 'Pattern scan plus content review; authorship names are intentional. Not an exhaustive secret-detection or copyright-clearance guarantee.'}
    out = ROOT / 'out'
    out.mkdir(exist_ok=True)
    (out / 'release_scan.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print('PRIVACY SCAN PASS: ' + str(len(scanned)) + ' files; no candidates; no third-party PDFs')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refresh', action='store_true', help='Refresh hashes after intentional author edits')
    args = parser.parse_args()
    if args.refresh:
        refresh()
    check_hashes()
    scan()
