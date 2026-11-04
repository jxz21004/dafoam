#!/usr/bin/env python3
"""Keep adjoint linear accuracy and finite-iteration FD stability separate."""
import argparse,json,re
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
base=a.root/'dafoam_case/validation/dafoam_L3_direct_second_2000/cfl5_0';root=a.root/'dafoam_case/validation/dafoam_L3_gradient_h1e4'
reference=np.load(base/'final_state.npz')['state'];reports={};history={};outer={}
num=r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'
for name in ['plus','minus','zero']:
    case=root/name;reports[name]=json.loads((case/'primal_report.json').read_text());assert reports[name]['primal_iterations']==1000
    assert np.array_equal(np.load(case/'initial_state.npz')['state'],reference)
    s=(case/'run.log').read_text();history[name]={};t=None
    for line in s.splitlines():
        m=re.match(r'Time = (\S+)',line)
        if m:t=int(float(m[1]))
        m=re.match(r'(cd|cl): (\S+) final:',line)
        if m and t is not None:history[name].setdefault(t,{})[m[1]]=float(m[2])
    v=re.findall(r'GMRES iteration:\s*0\s+Residual:\s*(.*?)\s+Krylov steps:',s)
    outer[name]=max(map(float,re.findall(num,v[-1])))
fd_history={str(t):{k:(history['plus'][t][k]-history['minus'][t][k])/2e-4 for k in ['cd','cl']} for t in sorted(history['plus'].keys()&history['minus'].keys())}
fd={k:(reports['plus']['functions'][k]-reports['minus']['functions'][k])/2e-4 for k in ['cd','cl']}
midpoint={k:0.5*(reports['plus']['functions'][k]+reports['minus']['functions'][k])-reports['zero']['functions'][k] for k in ['cd','cl']}
adj={}
for name in ['adjoint','adjoint_full']:
    path=root/name/'adjoint_report.json'
    if path.exists():
        r=json.loads(path.read_text());r['base_state_max_difference']=float(abs(np.load(root/name/'final_state.npz')['state']-reference).max());assert r['base_state_max_difference']==0
        for k,v in r['functions'].items():
            v['difference_from_fd_1000']=v['d_dalpha']-fd[k]
            v['linear_accuracy_pass']=v['ksp_reason']>0 and v['true_linear_relative_residual']<=1e-8
        adj[name]=r
    run_file=root/(name+'_runs.json')
    if run_file.exists():adj.setdefault(name,{})['run']=json.loads(run_file.read_text())
result=dict(variable='alpha_deg',central_difference_step_deg=1e-4,base_report=json.loads((base/'primal_report.json').read_text()),fd_reports=reports,fd_last_max_outer_residual=outer,fd_1000=fd,fd_by_continuation_iteration=fd_history,midpoint_minus_zero=midpoint,adjoints=adj,gradient_verified=False,limitation='FD perturbed states remain unconverged; derivatives depend on continuation length. A small adjoint linear residual would not alone validate the steady gradient.')
a.output.mkdir(parents=True,exist_ok=True);(a.output/'results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(fd_1000=fd,midpoint_minus_zero=midpoint,outer=outer,adjoints={k:v.get('functions',{}) for k,v in adj.items()}),indent=2))
