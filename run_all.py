"""Recompute six public mathematical checks sequentially in one process."""
import contextlib
import argparse
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import runpy
import sys
import time

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
SCRIPTS = ('label_structure', 'arc_certify', 'arc_review', 'lemmas_audit', 'coverage_audit', 'tom_arc_certify')


def main(only=None):
    if not __debug__:
        raise RuntimeError('Assertions are required: do not use python -O.')
    out = ROOT / 'out'
    out.mkdir(exist_ok=True)
    records = []
    start = time.monotonic()
    for name in ((only,) if only else SCRIPTS):
        path = ROOT / 'code' / (name + '.py')
        capture = io.StringIO()
        tick = time.monotonic()
        with contextlib.redirect_stdout(capture), contextlib.redirect_stderr(capture):
            runpy.run_path(str(path), run_name='__main__')
        seconds = time.monotonic() - tick
        (out / (name + '.log')).write_text(capture.getvalue(), encoding='utf-8')
        records.append({'script': 'code/' + path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'seconds': seconds, 'exit_code': 0})
        print(name + '.py PASS (' + format(seconds, '.3f') + ' s)')
    result = {'status': 'PASS', 'python': sys.version.split()[0],
              'packages': {n: importlib.metadata.version(n) for n in ('numpy', 'mpmath', 'sympy', 'python-flint', 'pypdf')},
              'execution': 'one process; sequential; numerical thread limits 1; no subprocess or worker',
              'seconds': time.monotonic() - start, 'checks': records}
    (out / 'reproducibility.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--only', choices=SCRIPTS, help='Run one check with the same thread limits and output setup')
    main(parser.parse_args().only)
