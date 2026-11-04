#!/usr/bin/env python3
"""Assess first-order convergence and second-order continuation separately."""
import argparse,json,re
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
num=r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'
rows=[]
for group in ['dafoam_L3_firstorder_2000','dafoam_L3_first_to_second_2000','dafoam_L3_direct_second_2000']:
    r=json.loads((a.root/group/'summary.json').read_text())[0]
    assert r['returncode']==0 and r['initial_state_max_difference']==0
    case=Path(r['case']);s=(case/'primal.log').read_text()
    outer=[max(map(float,re.findall(num,v))) for v in re.findall(r'GMRES iteration:\s*0\s+Residual:\s*(.*?)\s+Krylov steps:',s)]
    sa=re.findall(r'nuTilda initRes: (\S+) finalRes: (\S+) nIters: (\d+)',s)
    z=np.load(case/'final_state.npz');initial=np.load(case/'initial_residual.npz')['residual'];n=len(initial)//6
    r.update(group=group,outer_solves_logged=len(outer),last_max_outer_residual=outer[-1],min_max_outer_residual=min(outer),last_200_outer_range=[min(outer[-200:]),max(outer[-200:])],flow_normalized_threshold_met=outer[-1]<1e-8,last_sa_linear=dict(initial=float(sa[-1][0]),final=float(sa[-1][1]),iterations=int(sa[-1][2])),initial_sa_l2=float(np.linalg.norm(initial[5*n:])),final_sa_l2=float(np.linalg.norm(z['residual'][5*n:])),final_flow_l2=float(np.linalg.norm(z['residual'][:5*n])))
    rows.append(r)
first=np.load(Path(rows[0]['case'])/'final_state.npz')['state']
restart=np.load(Path(rows[1]['case'])/'initial_state.npz')['state']
assert np.array_equal(first,restart)
a.output.mkdir(parents=True,exist_ok=True)
(a.output/'results.json').write_text(json.dumps(dict(results=rows,first_to_second_state_max_difference=float(abs(first-restart).max()),two_stage_primal_seconds=sum(r['report']['primal_seconds'] for r in rows[:2])),indent=2)+'\n')
for r in rows:print(r['group'],r['report']['primal_iterations'],r['last_max_outer_residual'],r['report']['raw_residual_l2'],r['final_sa_l2'],r['report']['primal_seconds'])
