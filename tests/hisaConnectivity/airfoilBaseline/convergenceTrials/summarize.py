#!/usr/bin/env python3
from pathlib import Path
import re,json,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parent;root=out.parents[4]
results={};fig,ax=plt.subplots(1,3,figsize=(12,3.5),layout='constrained')
for name,logname in [('cfl5','cfl5'),('minmod_cfl5','minmod'),('tight_linear','tight'),('gauss_gradient','gauss'),('gauss_cfl50','gauss50'),('gauss_cfl500','gauss500'),('gauss_cfl5000','gauss5000'),('gauss_tight','gausstight')]:
 case=root/'reports/hisa_airfoil_baseline/convergence'/name
 logpath=case/'log.hisa'
 if not logpath.exists():logpath=Path('/tmp/rae_conv_'+logname+'.log')
 log=logpath.read_text()
 rows=np.loadtxt(case/'postProcessing/forceCoeffsCompressible/0/coefficient.dat');res=[]
 for line in log.splitlines():
  if 'GMRES iteration: 0 ' in line:res.append(max(float(v) for v in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?',line.split('Residual:')[1].split('Krylov')[0])))
 results[name]={'steps':int(re.findall(r'pseudoTime: iteration (\d+)',log)[-1]),'native_converged':'pseudoTime: converged' in log,'final_residual':res[-1],'force_sample_step':int(rows[-1,0]),'CD':float(rows[-1,1]),'CL':float(rows[-1,4]),'last_SA_solve':re.findall(r'smoothSolver:  Solving for nuTilda, .*',log)[-1],'timing':re.findall(r'ExecutionTime = .*',log)[-1]}
 ax[0].semilogy(np.arange(1,len(res)+1),res,label=name);ax[1].plot(rows[:,0],rows[:,1],label=name);ax[2].plot(rows[:,0],rows[:,4],label=name)
 folder=out/name;folder.mkdir(exist_ok=True)
 for n in ['fvSolution','fvSchemes','controlDict']:(folder/n).write_text((case/'system'/n).read_text().rstrip()+'\n')
 np.savetxt(folder/'history.csv',rows[:,[0,1,4]],delimiter=',',header='iteration,CD,CL',comments='')
 (folder/'final_excerpt.log').write_text('\n'.join(x.rstrip() for x in log.splitlines()[-35:]).rstrip()+'\n')
for case_name in ['tight_linear','gauss_gradient','gauss_cfl50','gauss_cfl500','gauss_cfl5000','gauss_tight']:
 raw=[]
 case=root/'reports/hisa_airfoil_baseline/convergence'/case_name
 for folder in sorted((p for p in case.iterdir() if p.name.isdigit()),key=lambda p:int(p.name)):
  row={'step':int(folder.name)}
  for name in ['rhoR','rhoUR','rhoER']:
   if not (folder/name).exists():continue
   ss=(folder/name).read_text();m=re.search(r'internalField\s+nonuniform\s+List<\w+>\s+\d+\s*\((.*?)\n\)',ss,re.S);assert m
   vals=np.fromstring(m[1].replace('(',' ').replace(')',' '),sep=' ')
   row[name+'_L2']=float(np.linalg.norm(vals));row[name+'_max']=float(abs(vals).max())
  if len(row)>1:raw.append(row)
 results[case_name+'_raw_residual_history']=raw
for a,title in zip(ax,['Native normalized outer max','CD','CL']):a.set(xlabel='Continuation iteration',ylabel=title);a.legend(fontsize=8)
fig.savefig(out/'comparison.png',dpi=160)
(out/'results.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,indent=2))
