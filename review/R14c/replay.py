#!/usr/bin/env python3
"""Run R14c in a temporary directory and compare frozen outputs (stdlib only)."""
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
CHECKS = (
    ('sturm_certificate.py', 'sturm_certificate.log'),
    ('verify.py', 'execution.log'),
    ('conservative_check.py', 'conservative.log'),
    ('pointcount_check.py', 'pointcounts.log'),
    ('prediction_check.py', 'predictions.log'),
)
if not __debug__:
    raise RuntimeError('Assertions are required; do not use python -O')
with TemporaryDirectory(prefix='j023-r14c-') as tmp:
    for script, log in CHECKS:
        output = subprocess.check_output([sys.executable, str(ROOT / script)], cwd=tmp)
        if output != (ROOT / log).read_bytes():
            raise RuntimeError(f'Output mismatch: {script}')
        print(f'PASS {script}')
    for name in ('results.json', 'pointcounts.json', 'predictions.json'):
        if (Path(tmp) / name).read_bytes() != (ROOT / name).read_bytes():
            raise RuntimeError(f'Output mismatch: {name}')
        print(f'PASS {name}')
print('PASS: 5 computations and 3 result files; no historical data modified')
