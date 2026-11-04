#!/usr/bin/env python3
"""Validate the bounded matrix-free experiment and retain unsuccessful solves."""
import argparse,json,re
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
num=r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'
old=a.root/'dafoam_L3_mf_old_baseline/cfl5_0';off=a.root/'dafoam_L3_mf_off_baseline/cfl5_0'
x=np.load(old/'final_state.npz');y=np.load(off/'final_state.npz');default=dict(state_max_difference=float(abs(x['state']-y['state']).max()),residual_max_difference=float(abs(x['residual']-y['residual']).max()))
check=(a.root/'dafoam_L3_mf_check/cfl5_0/primal.log').read_text();checks=[float(v) for v in re.findall(r'MatrixFreeCheck:.*relativeDifference (\S+)',check)];drifts=[float(v) for v in re.findall(r'repeatError (\S+)',check)]
assert default['state_max_difference']==0 and default['residual_max_difference']==0
assert len(checks)==3 and max(checks)<1e-4 and max(drifts)<1e-8
results=[]
for group in ['dafoam_L3_mf_off_baseline','dafoam_L3_mf_compare']:
 for row in json.loads((a.root/group/'summary.json').read_text()):
  assert row['returncode']==0 and row['initial_state_max_difference']==0
  text=(Path(row['case'])/'primal.log').read_text();trials=re.findall(r'MatrixFreeFlow: alpha (\S+) meritRatio (\S+) residualEvaluations (\d+) repeatError (\S+)',text)
  row['flow_updates_accepted']=sum(float(v[0])>0 for v in trials)
  row['flow_updates_rejected']=sum(float(v[0])==0 for v in trials)
  row['residual_evaluations_in_matrix_free_path']=sum(int(v[2]) for v in trials)
  row['max_repeat_error']=max([float(v[3]) for v in trials],default=0)
  row['line_search_history']=[dict(alpha=float(v[0]),merit_ratio=float(v[1]),evaluations=int(v[2])) for v in trials]
  row['cycles_at_limit']=len(re.findall(r'GMRES iteration:\s*3\s+Residual:',text)) if row['matrix_free'] else 0
  linear=[];initial=None
  for iteration,values in re.findall(r'GMRES iteration:\s*(\d+)\s+Residual:\s*(.*?)\s+Krylov steps:',text):
   values=np.array(list(map(float,re.findall(num,values))))
   if int(iteration)==0:initial=values
   elif initial is not None:
    linear.append(bool(np.all((values<1e-12)|(values/np.maximum(initial,1e-300)<0.05))))
  row['linear_solves_meeting_tolerance']=sum(linear)
  row['linear_solves_not_meeting_tolerance']=len(linear)-sum(linear)

  outer=re.findall(r'GMRES iteration:\s*0\s+Residual:\s*(.*?)\s+Krylov steps:',text)
  row['last_max_outer_residual']=max(map(float,re.findall(num,outer[-1])))
  row['flow_threshold_met']=row['last_max_outer_residual']<1e-8
  results.append(row)
summary=dict(default_path_regression=default,directional_checks={'scope':'First Krylov direction in each of three outer iterations, epsilon versus epsilon/2, frozen turbulence','max_relative_difference':max(checks),'max_zero_step_residual_drift':max(drifts),'passed':True},scope='20-step experiments at identical state; MF includes residual-based backtracking and is NOT a fully coupled flow/SA Newton method.',results=results)
a.output.mkdir(parents=True,exist_ok=True);(a.output/'results.json').write_text(json.dumps(summary,indent=2)+'\n')
for r in results:print(r['matrix_free'],r['cfl'],r['report']['fixed_definition_relative_l2'],r['report']['primal_seconds'],r['flow_updates_accepted'],r['cycles_at_limit'],r['last_max_outer_residual'])
