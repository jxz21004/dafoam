#!/usr/bin/env python3
"""Check first-order alpha gradients with converged paired FD states."""
import argparse,json,re
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
base=a.root/'dafoam_case/validation/dafoam_L3_firstorder_2000/cfl5_0';root=a.root/'dafoam_case/validation/dafoam_L3_first_gradient_h1e4'
reference=np.load(base/'final_state.npz')['state'];reports={};history={};outer={};sa={}
num=r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'
for name in ['plus','minus','zero']:
    case=root/name;reports[name]=json.loads((case/'primal_report.json').read_text())
    assert np.array_equal(np.load(case/'initial_state.npz')['state'],reference)
    s=(case/'run.log').read_text();history[name]={};t=None
    for line in s.splitlines():
        m=re.match(r'Time = (\S+)',line)
        if m:t=int(float(m[1]))
        m=re.match(r'(cd|cl): (\S+) final:',line)
        if m and t is not None:history[name].setdefault(t,{})[m[1]]=float(m[2])
    v=re.findall(r'GMRES iteration:\s*0\s+Residual:\s*(.*?)\s+Krylov steps:',s)
    outer[name]=max(map(float,re.findall(num,v[-1])));assert outer[name]<1e-10
    vals=re.findall(r'nuTilda initRes: (\S+) finalRes: (\S+) nIters: (\d+)',s)
    sa[name]=dict(initial=float(vals[-1][0]),final=float(vals[-1][1]))
fd_history={str(t):{k:(history['plus'][t][k]-history['minus'][t][k])/2e-4 for k in ['cd','cl']} for t in sorted(history['plus'].keys()&history['minus'].keys())}
fd={k:(reports['plus']['functions'][k]-reports['minus']['functions'][k])/2e-4 for k in ['cd','cl']}
midpoint={k:0.5*(reports['plus']['functions'][k]+reports['minus']['functions'][k])-reports['zero']['functions'][k] for k in ['cd','cl']}
adj=json.loads((root/'adjoint_full/adjoint_report.json').read_text())
assert np.array_equal(np.load(root/'adjoint_full/final_state.npz')['state'],reference)
for name,r in adj['functions'].items():
    r['fd']=fd[name];r['relative_difference_to_fd']=abs(r['d_dalpha']-fd[name])/max(abs(fd[name]),1e-30)
    r['independent_linear_pass_1e8']=r['ksp_reason']>0 and r['true_linear_relative_residual']<=1e-8
result=dict(variable='alpha_deg',step=1e-4,base_report=json.loads((base/'primal_report.json').read_text()),fd_reports=reports,fd_outer_residuals=outer,fd_sa_linear=sa,fd=fd,fd_history=fd_history,midpoint_minus_zero=midpoint,adjoint=adj,scope='First-order L3 alpha gradients at one FD step only; not a second-order or geometry-gradient validation.')
a.output.mkdir(parents=True,exist_ok=True);(a.output/'results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(fd=fd,midpoint=midpoint,adjoint=adj['functions']),indent=2))
