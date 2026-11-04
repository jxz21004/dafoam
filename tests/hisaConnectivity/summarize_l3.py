#!/usr/bin/env python3
"""Summarize native HiSA outer residuals separately from GMRES inner residuals."""
import argparse
import csv
import json
from pathlib import Path
import re
import statistics

p = argparse.ArgumentParser()
p.add_argument('case', type=Path)
a = p.parse_args()
text = (a.case / 'primal.log').read_text()
number = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'
rows = []
iteration = 0
cfl = None
for line in text.splitlines():
    m = re.search(r'pseudoTime: iteration\s+(\d+)', line)
    if m:
        iteration = int(m[1])
    m = re.search(r'Pseudo Courant No:\s*(' + number + ')', line)
    if m:
        cfl = float(m[1])
    m = re.search(r'GMRES iteration:\s*0\s+Residual:\s*(.*?)\s+Krylov steps:', line)
    if m:
        vals = [float(v) for v in re.findall(number, m[1])]
        assert len(vals) == 5
        rows.append(dict(iteration=iteration, cfl=cfl, rho=vals[0], rhoUx=vals[1], rhoUy=vals[2], rhoUz=vals[3], rhoE=vals[4], max_outer_residual=max(vals), krylov_steps=0))
    m = re.search(r'GMRES iteration:\s*[1-9]\d*.*Krylov steps:\s*(\d+)', line)
    if m and rows:
        rows[-1]['krylov_steps'] = int(m[1])  # cumulative across restart cycles
assert rows, 'No outer residuals found'
sa = re.findall(r'nuTilda initRes:\s*(' + number + r')\s+finalRes:\s*(' + number + ')', text)
report = json.loads((a.case / 'primal_report.json').read_text())
last = rows[-1]
summary = dict(threshold=1e-8, residual_definition='HiSA normalized outer initial residual (not the GMRES final inner residual and not a strict initial-relative residual)',
               flow_threshold_met=last['max_outer_residual'] < 1e-8,
               sa_last_initial_residual=float(sa[-1][0]) if sa else None,
               sa_last_final_residual=float(sa[-1][1]) if sa else None,
               iterations=report['primal_iterations'], outer_samples=len(rows), seconds=report['primal_seconds'],
               gmres_subspace=40, gmres_relative_tolerance=0.05, preconditioner='LUSGS; LaxFriedrichs inviscid + laplacian viscous + pseudo-time diagonal',
               cfl_cap=100000, cfl_observed_max=max(v['cfl'] for v in rows), last=last,
               total_krylov_steps=sum(v['krylov_steps'] for v in rows),
               last100_mean_outer=statistics.mean(v['max_outer_residual'] for v in rows[-100:]),
               samples={str(i): next((v for v in rows if v['iteration']==i), None) for i in [1,100,1000,2000,5000,10000]},
               fixed_definition_relative_l2=report['fixed_definition_relative_l2'], raw_residual_l2=report['raw_residual_l2'],
               raw_residual_max=report['raw_residual_max'], temperature_range=[report['temperature_min'],report['temperature_max']], functions=report['functions'],
               solver_converged_message=bool(re.search(r'pseudoTime: converged in',text)), gradient_verified=False)
summary['converged'] = summary['flow_threshold_met'] and bool(sa) and summary['sa_last_initial_residual'] < 1e-8
with (a.case/'convergence.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
(a.case/'convergence_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))

import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/hisa-mpl-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axes=plt.subplots(2,1,figsize=(8,6),sharex=True)
x=[v['iteration'] for v in rows]
axes[0].semilogy(x,[v['max_outer_residual'] for v in rows],label='Max outer flow residual')
axes[0].axhline(1e-8,color='black',linestyle='--',label='Requested threshold')
axes[0].set_ylabel('HiSA normalized residual');axes[0].legend()
axes[1].semilogy(x,[v['cfl'] for v in rows],label='Actual CFL')
axes[1].axhline(100000,color='black',linestyle='--',label='CFL upper limit')
axes[1].set_ylabel('CFL');axes[1].set_xlabel('Outer iteration');axes[1].legend()
for ax in axes:ax.grid(True,alpha=.3)
fig.tight_layout();fig.savefig(a.case/'convergence.png',dpi=160)
