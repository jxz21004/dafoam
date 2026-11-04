#!/usr/bin/env python3
"""Audit bounded-update acceptance, CFL feedback, cost and convergence."""
import argparse
import json
import re
from pathlib import Path
import numpy as np

p = argparse.ArgumentParser()
p.add_argument('--root', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
number = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'
results = []
for group in ['dafoam_L3_limited_baseline', 'dafoam_L3_limited_compare']:
    for row in json.loads((a.root / group / 'summary.json').read_text()):
        assert row['returncode'] == 0 and row['initial_state_max_difference'] == 0
        log = (Path(row['case']) / 'primal.log').read_text()
        entries = re.findall(r'LimitedUpdate: alpha (\S+) steadyRatio (\S+) pseudoRatio (\S+) residualEvaluations (\d+) cfl (\S+) nextCfl (\S+)', log)
        history = [dict(zip(['alpha', 'steady_ratio', 'pseudo_ratio', 'evaluations', 'cfl', 'next_cfl'], map(float, e))) for e in entries]
        if row['limited_update']:
            assert len(history) == row['report']['primal_iterations']
            for i, h in enumerate(history):
                assert 0 <= h['alpha'] <= 1
                assert h['evaluations'] <= row['update_max_trials'] + 2
                if h['alpha'] > 0:
                    assert h['pseudo_ratio'] <= 1 - 1e-4 * h['alpha'] + 1e-12
                factor = 0.5 if h['alpha'] < 0.25 else (1.2 if h['alpha'] >= 0.8 and h['steady_ratio'] < 0.99 else 1)
                assert np.isclose(h['next_cfl'], min(100000, max(0.1, h['cfl'] * factor)))
                if i:
                    assert np.isclose(h['cfl'], history[i-1]['next_cfl'])
            drifts = [float(v) for v in re.findall(r'repeatError (\S+)', log)]
            assert len(drifts) == len(history) and max(drifts) < 1e-8
            row['max_zero_step_drift'] = max(drifts)
        outer = re.findall(r'GMRES iteration:\s*0\s+Residual:\s*(.*?)\s+Krylov steps:', log)
        row['last_max_outer_residual'] = max(map(float, re.findall(number, outer[-1])))
        row['flow_threshold_met'] = row['last_max_outer_residual'] < 1e-8
        row['updates_accepted'] = sum(h['alpha'] > 0 for h in history)
        row['updates_rejected'] = sum(h['alpha'] == 0 for h in history)
        row['extra_residual_evaluations'] = sum(int(h['evaluations']) for h in history)
        row['alpha_range'] = [min(h['alpha'] for h in history), max(h['alpha'] for h in history)] if history else None
        row['accepted_steps_with_increased_steady_merit'] = sum(h['alpha'] > 0 and h['steady_ratio'] > 1 for h in history)
        row['max_extra_evaluations_per_step'] = max((int(h['evaluations']) for h in history), default=0)
        results.append(row)
        print(row['limited_update'], row['cfl'], row['report']['raw_residual_l2'], row['report']['primal_seconds'], row['peak_rss_kb'], row['last_max_outer_residual'], row['updates_rejected'])
# The disabled branch must preserve the previously validated trajectory.
old = (a.root / 'dafoam_L3_mf_off_baseline/cfl5_0/primal.log').read_text()
new = (a.root / 'dafoam_L3_limited_baseline/cfl5_0/primal.log').read_text()
pattern = r'GMRES iteration:.*'
old_gmres = re.findall(pattern, old)
new_gmres = re.findall(pattern, new)
assert old_gmres and old_gmres == new_gmres[:len(old_gmres)]
a.output.mkdir(parents=True, exist_ok=True)
(a.output / 'results.json').write_text(json.dumps({'default_off_first_20_iterations_identical': True, 'results': results}, indent=2) + '\n')
