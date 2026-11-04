#!/usr/bin/env python3
"""Bounded fixed-CFL comparisons from identical saved L3 fields."""
import argparse
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('--source',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--source-time',default='10000')
p.add_argument('--hisa-repo',type=Path,required=True)
p.add_argument('--hisa-lib',type=Path,required=True)
p.add_argument('--steps',type=int,default=500)
p.add_argument('--cfl',type=float,nargs='+',default=[58,100,300,1000])
a=p.parse_args();source=a.source.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
runner=a.hisa_repo.resolve()/'tests/wingAero/run_case.py'
reference=np.load(source/'final_state.npz')['state']
results=[]
for cfl in a.cfl:
 case=out/('cfl'+str(cfl).replace('.','_'));case.mkdir()
 for name in ['constant','system']:shutil.copytree(source/name,case/name)
 # Keep only volume fields, not saved pseudo-time/CFL control state.
 (case/'0').mkdir()
 for f in (source/a.source_time).iterdir():
  if f.is_file() and not f.name.startswith('pseudoCo'):shutil.copy2(f,case/'0'/f.name)
 path=case/'system/controlDict';s=path.read_text()
 for key,value in [('endTime',a.steps),('writeInterval',a.steps)]:
  s,n=re.subn(r'(?m)^(\s*'+key+r'\s+)[^;]+;',lambda m:m[1]+str(value)+';',s);assert n==1
 path.write_text(s)
 path=case/'system/fvSolution';s=path.read_text()
 for key,value in [('pseudoCoNum',cfl),('pseudoCoNumMax',100000),('pseudoCoNumMaxIncreaseFactor',1)]:
  s,n=re.subn(r'(?m)^(\s*'+key+r'\s+)[^;]+;',lambda m:m[1]+str(value)+';',s);assert n==1
 for key,value in [('resetPseudo','true'),('pseudoCoNumMinDecreaseFactor','1')]:
  if re.search(r'(?m)^\s*'+key+r'\s+',s):
   s=re.sub(r'(?m)^(\s*'+key+r'\s+)[^;]+;',lambda m:m[1]+value+';',s)
  else:
   s=s.replace('pseudoTime\n{','pseudoTime\n{\n    '+key+' '+value+';')
 path.write_text(s)
 start=time.monotonic()
 with (case/'primal.log').open('w') as log:
  try:
   proc=subprocess.run([sys.executable,str(runner),'--case',str(case),'--library-dir',str(a.hisa_lib.resolve()),'--restart'],stdout=log,stderr=subprocess.STDOUT,timeout=180)
   code=proc.returncode
  except subprocess.TimeoutExpired:code='timeout'
 row={'cfl':cfl,'returncode':code,'wall_seconds_including_startup':time.monotonic()-start,'case':str(case)}
 text=(case/'primal.log').read_text();actual=[float(v) for v in re.findall(r'Pseudo Courant No:\s*(\S+)',text)]
 row['actual_cfl_range']=[min(actual),max(actual)] if actual else None
 if code==0:
  report=json.loads((case/'primal_report.json').read_text());row['report']=report
  # Accumulate the last cumulative Krylov count of each outer solve.
  counts=[]
  for line in text.splitlines():
   m=re.search(r'GMRES iteration:\s*(\d+).*Krylov steps:\s*(\d+)',line)
   if m:
    if int(m[1])==0:counts.append(0)
    else:counts[-1]=int(m[2])
  row['total_krylov_steps']=sum(counts)
  row['decades_reduction_per_primal_second']=-float(np.log10(report['fixed_definition_relative_l2']))/report['primal_seconds']
  row['initial_state_max_difference']=float(np.max(abs(np.load(case/'initial_state.npz')['state']-reference)))
  assert row['initial_state_max_difference'] < 1e-12, 'Restart state changed'
  assert actual and all(abs(v-cfl)<1e-10 for v in actual), 'CFL was not held fixed'
 results.append(row);(out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
 print(json.dumps(row),flush=True)
