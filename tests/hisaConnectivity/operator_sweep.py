#!/usr/bin/env python3
"""Compare native HiSA implicit operators and update bounds at a common state."""
import argparse,json,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--source-time',default='500');p.add_argument('--output',type=Path,required=True);p.add_argument('--hisa-repo',type=Path,required=True);p.add_argument('--hisa-lib',type=Path,required=True);p.add_argument('--steps',type=int,default=300);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
variants=[('baseline5','LaxFriedrichs',None,5,False),('central5','centralLaxFriedrichs',None,5,False),('discrete5','discreteLaxFriedrichs','centralLaxFriedrichs',5,False),('discrete100','discreteLaxFriedrichs','centralLaxFriedrichs',100,False),('discrete100_bounded','discreteLaxFriedrichs','centralLaxFriedrichs',100,True),('baseline100_bounded','LaxFriedrichs',None,100,True)]
results=[]
for name,jac,pc,cfl,bound in variants:
 command=[sys.executable,str(Path(__file__).with_name('cfl_sweep.py')),'--source',str(a.source),'--source-time',a.source_time,'--output',str(a.output/name),'--hisa-repo',str(a.hisa_repo),'--hisa-lib',str(a.hisa_lib),'--steps',str(a.steps),'--cfl',str(cfl),'--jacobian',jac]
 if pc:command+=['--pc-jacobian',pc]
 if bound:command+=['--bounding']
 with (a.output/(name+'.log')).open('w') as log:
  result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
 path=a.output/name/'summary.json';row={'variant':name,'driver_returncode':result.returncode}
 if path.exists():row.update(json.loads(path.read_text())[0])
 results.append(row);(a.output/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
 print(name,row.get('returncode'),row.get('report',{}).get('fixed_definition_relative_l2'),flush=True)
