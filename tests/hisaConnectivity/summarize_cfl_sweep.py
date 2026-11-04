#!/usr/bin/env python3
"""Compare fixed-CFL tests using the same starting residual and wall-time basis."""
import argparse,json,os,re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--runs',type=Path,nargs='+',required=True);p.add_argument('--continuation',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
rows=[]
for path in a.runs:
 for row in json.loads((path/'summary.json').read_text()):
  assert row['returncode']==0 and row['initial_state_max_difference']==0
  r=row['report'];rows.append(dict(cfl=row['cfl'],iterations=r['primal_iterations'],residual_l2_initial=r['initial_raw_residual_l2'],residual_l2_final=r['raw_residual_l2'],residual_ratio=r['fixed_definition_relative_l2'],seconds=r['primal_seconds'],krylov_steps=row['total_krylov_steps'],decades_per_second=row['decades_reduction_per_primal_second'],actual_cfl_range=row['actual_cfl_range']))
rows.sort(key=lambda r:r['cfl']);assert max(r['residual_l2_initial'] for r in rows)==min(r['residual_l2_initial'] for r in rows)
cont=json.loads((a.continuation/'summary.json').read_text())[0];r=cont['report'];first=next(v for v in rows if v['cfl']==5)
text=(Path(cont['case'])/'primal.log').read_text();values=re.findall(r'GMRES iteration:\s*0\s+Residual:\s*(.*?)\s+Krylov steps:',text);numbers=r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?';outer=[max(map(float,re.findall(numbers,v))) for v in values]
continuation=dict(cfl=5,additional_iterations=r['primal_iterations'],total_iterations=first['iterations']+r['primal_iterations'],additional_seconds=r['primal_seconds'],total_seconds=first['seconds']+r['primal_seconds'],residual_l2_final=r['raw_residual_l2'],residual_ratio_over_continuation=r['fixed_definition_relative_l2'],residual_ratio_from_common_start=r['raw_residual_l2']/first['residual_l2_initial'],total_krylov_steps=first['krylov_steps']+cont['total_krylov_steps'],last_normalized_outer=outer[-1],last100_mean_outer=sum(outer[-100:])/len(outer[-100:]),flow_threshold_met=outer[-1]<1e-8)
report=dict(scope='Fixed-CFL diagnostic from identical stalled L3 state; not an ANK implementation. Timing is a single serial run per setting, not a statistical benchmark.',sweep=rows,continuation=continuation)
a.output.mkdir(parents=True,exist_ok=True);(a.output/'comparison.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
os.environ.setdefault('MPLCONFIGDIR','/tmp/hisa-mpl-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axes=plt.subplots(1,2,figsize=(10,4))
x=[str(v['cfl']) for v in rows]
axes[0].bar(x,[v['residual_ratio'] for v in rows]);axes[0].axhline(1,color='black',linestyle='--');axes[0].set_ylabel('Final / initial residual L2 (500 steps)')
axes[1].bar(x,[v['seconds'] for v in rows]);axes[1].set_ylabel('Primal wall time (s)')
for ax in axes:ax.set_xlabel('Fixed CFL');ax.grid(True,axis='y',alpha=.3)
fig.tight_layout();fig.savefig(a.output/'comparison.png',dpi=160)
