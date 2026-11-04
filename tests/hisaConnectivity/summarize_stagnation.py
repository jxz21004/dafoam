#!/usr/bin/env python3
"""Audit frozen-SA and actual-update-direction experiments."""
import argparse
import json
import re
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
results=[]
for group in ['dafoam_L3_sa_audit','dafoam_L3_sa_control','dafoam_L3_sa_frozen','dafoam_L3_direction_ready','dafoam_L3_direction_control']:
    row=json.loads((a.root/group/'summary.json').read_text())[0]
    assert row['returncode']==0 and row['initial_state_max_difference']==0
    case=Path(row['case']); initial=np.load(case/'initial_state.npz')['state']; final=np.load(case/'final_state.npz'); n=len(initial)//6
    row['sa_state_max_change']=float(abs(final['state'][5*n:]-initial[5*n:]).max())
    row['sa_residual_l2']=float(np.linalg.norm(final['residual'][5*n:]))
    row['flow_residual_l2']=float(np.linalg.norm(final['residual'][:5*n]))
    if row['freeze_sa']:assert row['sa_state_max_change']==0
    log=(case/'primal.log').read_text()
    row['sa_audit']=[dict(zip(['lagged_l2','refreshed_l2','full_l2','refresh_difference_l2','boundary_difference_l2','phi_difference_l2'],map(float,v))) for v in re.findall(r'SAAudit: laggedL2 (\S+) refreshedL2 (\S+) fullL2 (\S+) refreshDifferenceL2 (\S+) boundaryDifferenceL2 (\S+) phiDifferenceL2 (\S+)',log)]
    row['direction_derivatives']=[dict(zip(['h','relative_error','true_derivative_norm'],map(float,v))) for v in re.findall(r'DirectionAuditJ: h (\S+) relativeError (\S+) trueDerivativeNorm (\S+)',log)]
    row['direction_steps']=[dict(zip(['alpha','predicted_ratio','actual_ratio'],map(float,v))) for v in re.findall(r'DirectionAuditStep: alpha (\S+) predictedRatio (\S+) actualRatio (\S+)',log)]
    drifts=list(map(float,re.findall(r'DirectionAuditRepeat: (\S+)',log)))
    if row['direction_audit']:
        assert len(drifts)==3 and max(drifts)<1e-8
        assert len(row['direction_derivatives'])==6 and len(row['direction_steps'])==12
        for j in range(0,6,2):
            x,y=row['direction_derivatives'][j:j+2]
            assert np.isclose(x['relative_error'],y['relative_error'],rtol=1e-6)
        row['max_zero_step_drift']=max(drifts)
    results.append(row)
regressions=[]
for x,y in [('dafoam_L3_direction_ready','dafoam_L3_direction_control'),('dafoam_L3_sa_audit','dafoam_L3_mf_off_baseline')]:
    u=np.load(a.root/x/'cfl5_0/final_state.npz');v=np.load(a.root/y/'cfl5_0/final_state.npz')
    errors={k:float(abs(u[k]-v[k]).max()) for k in ['state','residual']}
    assert all(z==0 for z in errors.values())
    regressions.append(dict(cases=[x,y],max_differences=errors))
a.output.mkdir(parents=True,exist_ok=True)
(a.output/'results.json').write_text(json.dumps(dict(results=results,regressions=regressions),indent=2)+'\n')
print('Passed initial-state, frozen-SA, diagnostic-restoration and step-size consistency checks.')
