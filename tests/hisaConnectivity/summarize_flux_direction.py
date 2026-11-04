#!/usr/bin/env python3
"""Compare flux derivatives at the same state and saved conservative direction."""
import argparse,gzip,hashlib,json,re
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
ref=a.root/'dafoam_L3_fixed_direction_reference/cfl5_0'
def hashes(case):
    result={}
    for name in ['auditDeltaRho','auditDeltaE','auditDeltaMom']:
        f=case/'constant'/name
        data=f.read_bytes() if f.exists() else gzip.decompress(f.with_suffix('.gz').read_bytes())
        result[name]=hashlib.sha256(data).hexdigest()
    return result
reference_hashes=hashes(ref);rows=[]
for label in ['ausm_high','ausm_first','lf_central','lf_geometric']:
    row=json.loads((a.root/('dafoam_L3_fixed_direction_v2_'+label)/'summary.json').read_text())[0]
    assert row['returncode']==0 and row['initial_state_max_difference']==0
    case=Path(row['case']);assert hashes(case)==reference_hashes
    text=(case/'primal.log').read_text()
    pairs=re.findall(r'DirectionAuditJ: h (\S+) relativeError (\S+) trueDerivativeNorm (\S+) densityError (\S+) momentumError (\S+) energyError (\S+)',text)
    checks=[dict(zip(['h','relative_error','true_derivative_norm','density_error','momentum_error','energy_error'],map(float,x))) for x in pairs]
    assert len(checks)==2
    for k in checks[0]:
        if k!='h':assert np.isclose(checks[0][k],checks[1][k],rtol=1e-6)
    drift=float(re.search(r'DirectionAuditRepeat: (\S+)',text)[1]);assert drift<1e-8
    row.update(label=label,derivative_checks=checks,zero_step_drift=drift,direction_file_hashes=hashes(case))
    rows.append(row)
u=np.load(ref/'final_state.npz');v=np.load(Path(rows[0]['case'])/'final_state.npz')
restore={k:float(abs(u[k]-v[k]).max()) for k in ['state','residual']};assert all(x==0 for x in restore.values())
a.output.mkdir(parents=True,exist_ok=True)
(a.output/'results.json').write_text(json.dumps(dict(results=rows,stored_direction_roundtrip=restore),indent=2)+'\n')
for r in rows:print(r['label'],r['derivative_checks'][0])
