#!/usr/bin/env python3
"""Overlay saved first-order and wVanLeer test histories."""
from pathlib import Path
import json,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parent
fig,axs=plt.subplots(1,3,figsize=(12,3.5),layout='constrained')
for folder,log,label in [(out.parent/'a279_m0734',Path('/tmp/rae279.log'),'wVanLeer'),(out,Path('/tmp/rae279_firstorder.log'),'firstOrder')]:
 data=np.loadtxt(folder/'coefficients.csv',delimiter=',',skiprows=1)
 residual=[]
 for line in log.read_text().splitlines():
  if 'GMRES iteration: 0 ' in line:
   vals=[float(v) for v in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?',line.split('Residual:')[1].split('Krylov')[0])];residual.append(max(vals))
 axs[0].semilogy(np.arange(1,len(residual)+1),residual,label=label)
 axs[1].plot(data[:,0],data[:,1],label=label);axs[2].plot(data[:,0],data[:,2],label=label)
for ax,title in zip(axs,['Native outer normalized max','CD','CL']):ax.set(xlabel='Iteration',ylabel=title);ax.legend()
axs[0].axhline(1e-8,color='k',ls='--');fig.savefig(out/'order_comparison.png',dpi=160)
a=json.loads((out/'summary.json').read_text());b=json.loads((out.parent/'a279_m0734/summary.json').read_text())
result={k:{'first_order':a[k],'wVanLeer':b[k],'relative_difference':a[k]/b[k]-1} for k in ['CL','CD']}
(out/'comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
