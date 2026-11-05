#!/usr/bin/env python3
"""Summarize the current-HiSA standalone RAE2822 smoke test."""
from pathlib import Path
import hashlib,json,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('/mdouser/test_aerostructure');out=Path(__file__).resolve().parent
case=root/'reports/hisa_airfoil_baseline/current_rae2822/simulation'
log=Path('/tmp/rae_current.log').read_text()
rows=np.loadtxt(case/'postProcessing/forceCoeffsCompressible/0/coefficient.dat')
res=[]
for line in log.splitlines():
 if 'GMRES iteration: 0 ' in line:
  part=line.split('Residual:')[1].split('Krylov')[0];vals=[float(x) for x in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?',part)];res.append(max(vals))
lib=root/'dafoam_case/validation/hisa_direction_build/lib/libHisaSolver.so'
r={'case':str(case),'library':str(lib),'library_sha256':hashlib.sha256(lib.read_bytes()).hexdigest(),'cells':23552,'iterations':int(re.findall(r'pseudoTime: iteration (\d+)',log)[-1]),'native_converged':'pseudoTime: converged' in log,'last_normalized_outer_max':res[-1],'last_force_time':float(rows[-1,0]),'CD':float(rows[-1,1]),'CL':float(rows[-1,4]),'last_100_steps_CD_range':[float(rows[rows[:,0]>=rows[-1,0]-100,1].min()),float(rows[rows[:,0]>=rows[-1,0]-100,1].max())],'last_100_steps_CL_range':[float(rows[rows[:,0]>=rows[-1,0]-100,4].min()),float(rows[rows[:,0]>=rows[-1,0]-100,4].max())],'last_timing':re.findall(r'ExecutionTime = .*',log)[-1],'final_log_has_End':log.rstrip().endswith('End'),'external_accuracy_validated':False}
(out/'summary.json').write_text(json.dumps(r,indent=2)+'\n')
np.savetxt(out/'coefficients.csv',rows[:,[0,1,4]],delimiter=',',header='iteration,CD,CL',comments='')
fig,ax=plt.subplots(1,3,figsize=(12,3.5),layout='constrained');ax[0].semilogy(np.arange(1,len(res)+1),res);ax[0].axhline(1e-8,color='k',ls='--');ax[0].set(ylabel='Native normalized outer max',xlabel='Iteration');ax[1].plot(rows[:,0],rows[:,1]);ax[1].set(xlabel='Iteration',ylabel='CD');ax[2].plot(rows[:,0],rows[:,4]);ax[2].set(xlabel='Iteration',ylabel='CL');fig.savefig(out/'convergence.png',dpi=160);print(json.dumps(r,indent=2))
