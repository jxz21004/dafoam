#!/usr/bin/env python3
"""Diagnostic alpha gradients at a nonconverged base; no accuracy gate is bypassed silently."""
import argparse,json,re,shutil,subprocess,sys,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=['adjoint','adjoint_full','fd'],required=True);a=p.parse_args()
root=a.root.resolve(); source=root/'dafoam_case/validation/dafoam_L3_direct_second_2000/cfl5_0'; out=root/'dafoam_case/validation/dafoam_L3_gradient_h1e4';out.mkdir(exist_ok=True)
runner=root/'Hisa4DAFoam/tests/wingAero/run_case.py';lib=root/'dafoam_case/validation/hisa_direction_build/lib'
variants=[(a.mode,2.)] if a.mode.startswith('adjoint') else [('plus',2.0001),('minus',1.9999),('zero',2.)]
results=[]
for name,alpha in variants:
    case=out/name;case.mkdir()
    for folder in ['constant','system']:shutil.copytree(source/folder,case/folder)
    (case/'0').mkdir()
    for f in (source/'2000').iterdir():
        if f.is_file() and not f.name.startswith('pseudoCo'):shutil.copy2(f,case/'0'/f.name)
    shutil.copy2(source/'final_state.npz',case/'final_state.npz')
    pth=case/'system/controlDict';s=pth.read_text()
    for key,value in [('endTime',1000),('writeInterval',1000)]:
        s,n=re.subn(r'(?m)^(\s*'+key+r'\s+)[^;]+;',lambda m:m[1]+str(value)+';',s);assert n==1
    pth.write_text(s)
    cmd=[sys.executable,str(runner),'--case',str(case),'--library-dir',str(lib),'--restart','--alpha',str(alpha)]
    if name.startswith('adjoint'):cmd+=['--probe','--adjoint','--allow-unconverged-adjoint','--adjoint-restart','40']
    if name=='adjoint_full':cmd+=['--adjoint-full-lu']
    start=time.monotonic()
    with (case/'run.log').open('w') as log:
        try:code=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=600 if name.startswith('adjoint') else 180).returncode
        except subprocess.TimeoutExpired:code='timeout'
    r=dict(name=name,alpha_deg=alpha,returncode=code,wall_seconds=time.monotonic()-start,case=str(case));results.append(r)
    (out/(a.mode+'_runs.json')).write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(r),flush=True)
