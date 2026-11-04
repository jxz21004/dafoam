#!/usr/bin/env python3
"""Frozen L3 connectivity/AD checks using the locally installed DAFoam package."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np
import dafoam

p = argparse.ArgumentParser()
p.add_argument('--source-case', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--hisa-repo', type=Path, required=True)
p.add_argument('--hisa-lib', type=Path, required=True)
a = p.parse_args()
repo = Path(__file__).resolve().parents[1]
installed = repo / 'build/hisa-local/python'
assert Path(dafoam.__file__).resolve().is_relative_to(installed), dafoam.__file__
out = a.output.resolve()
out.mkdir(parents=True, exist_ok=False)
case = out / 'case'
case.mkdir()
for name in ('0', 'constant', 'system'):
    shutil.copytree(a.source_case / name, case / name)
shutil.copy2(a.source_case / 'final_state.npz', case / 'final_state.npz')
reference = np.load(case / 'final_state.npz')['residual'].copy()
scripts = a.hisa_repo.resolve() / 'tests/highOrderNewton'
for script in ('audit_connectivity.py', 'check_parity.py', 'check_derivatives.py'):
    command = [sys.executable, str(scripts / script), '--case', str(case), '--library-dir', str(a.hisa_lib.resolve())]
    if script == 'audit_connectivity.py':
        command += ['--dafoam-library-dir', str(repo / 'build/hisa-local/lib')]
    with (out / (script + '.log')).open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
connectivity = json.loads((case / 'connectivity_report.json').read_text())
derivative = json.loads((case / 'derivative_report.json').read_text())
report = {
    'python_package': str(Path(dafoam.__file__).resolve()),
    'loaded_dafoam': connectivity['loaded_dafoam'],
    'max_column_relative_error': max(v['relative_column_error'] for v in connectivity['columns']),
    'missing_sa_entries': sum(v['sa_missing_count'] for v in connectivity['columns']),
    'native_ad_residual_relative_difference': derivative['native_ad_residual_relative_difference'],
    'max_best_block_derivative_error': max(v['relative_error'] for v in derivative['best_by_block']),
    'residual_max_change_from_previous_version': float(abs(np.load(case / 'final_state.npz')['residual'] - reference).max()),
    'scope': 'Five Ux columns and 16 random AD/FD block projections at frozen L3 state; not steady convergence or aerodynamic gradient validation.'
}
report['passed'] = (report['missing_sa_entries'] == 0 and report['max_column_relative_error'] < 1e-4
                    and report['native_ad_residual_relative_difference'] < 1e-12
                    and report['max_best_block_derivative_error'] < 1e-4)
(out / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2), flush=True)
assert report['passed'], 'Validation gates failed'
