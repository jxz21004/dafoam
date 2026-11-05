#!/usr/bin/env python3
"""Run the user's two-point static-aeroelastic models in isolated case folders."""
import argparse,hashlib,json,os,re,shutil,sys,time,traceback
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--solver',choices=['adflow','hisa'],required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
root=a.root.resolve();work=a.output.resolve();work.mkdir(parents=True,exist_ok=False)
source=root/('asopt2.py' if a.solver=='adflow' else 'asopt2-hisa.py')
s=source.read_text().split('prob = om.Problem()')[0]
# Both original scripts have identical Aitken settings except the iteration cap.
s=re.sub(r'maxiter = (?:50|100),\n(\s*)iprint = 2,\n(\s*)atol = 1e-4,',r'maxiter = 100,\n\1iprint = 2,\n\2atol = 1e-4,',s)
if a.solver=='hisa':
    s=s.replace('f"mesh_aero.{MPhysVariables.Aerodynamics.Surface.Mesh.COORDINATES}"', '"mesh_aero.x_aero0"')
(work/'model_source.py').write_text(s)
for name in ['tacs_setup.py','dafoam_group.py','optrecord.py']:shutil.copy2(root/name,work/name)
(work/'mphys_input_files').symlink_to(root/'mphys_input_files',target_is_directory=True)
(work/'output').mkdir()
if a.solver=='hisa':
    case=work/'dafoam_case';template=root/'dafoam_case/validation/dafoam_L3_firstorder_2000/cfl5_0';case.mkdir()
    for name in ['constant','system','0']:shutil.copytree(template/name,case/name)
    path=case/'system/controlDict';cfg=path.read_text()
    for key,value in [('endTime',4000),('writeInterval',4000)]:
        cfg,n=re.subn(r'(?m)^(\s*'+key+r'\s+)[^;]+;',lambda m:m[1]+str(value)+';',cfg);assert n==1
    path.write_text(cfg)
# Preserve each script's physics and structural definitions; run analysis only.
os.chdir(work);sys.path.insert(0,str(work));sys.argv=[str(source),'--level','L3','--task','run']
ns={'__name__':'comparison_model','__file__':str(source)};start=time.monotonic()
report={'solver':a.solver,'source_script':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'case':str(work),'mpi_ranks':1,'coupling_maxiter':100,'coupling_atol':1e-4,'coupling_rtol':1e-8}
try:
    exec(compile(s,str(work/'model_source.py'),'exec'),ns)
    import numpy as np
    ns['inputdir']=str(root/'mphys_input_files')+'/'
    ns['outputdir']=str(work/'output')+'/'
    om=ns['om'];prob=om.Problem(model=ns['AerostructuralFlightPoint']())
    prob.setup(mode='rev');prob.final_setup()
    report['setup_seconds']=time.monotonic()-start
    t=time.monotonic();prob.run_model();report['solve_seconds']=time.monotonic()-t;report['completed']=True
    vals={};arrays={}
    for name,meta in prob.model.list_outputs(out_stream=None,prom_name=True):
        if name.startswith(('cruise.','maneuver.')):
            v=np.asarray(meta['val']);prom=meta.get('prom_name',name)
            if v.size==1:vals[prom]=float(v.flat[0])
            if any(token in name for token in ['u_struct','u_aero','f_aero','f_struct','displacement']):
                arrays[name]=v;vals[name]={'shape':list(v.shape),'norm':float(np.linalg.norm(v)),'min':float(v.min()),'max':float(v.max())}
    np.savez(work/'coupled_fields.npz',**arrays);report['outputs']=vals
    report['coupling_iterations']={k:getattr(prob.model,k).coupling.nonlinear_solver._iter_count for k in ['cruise','maneuver']}
    if any(n == 0 for n in report['coupling_iterations'].values()):
        raise RuntimeError('No coupling iteration executed; results are not valid')
    if not all(np.isfinite(v).all() for v in arrays.values()):
        raise RuntimeError('Nonfinite coupled fields; results are not valid')
except BaseException as e:
    report['completed']=False;report['error']=str(e);traceback.print_exc()
finally:
    report['wall_seconds']=time.monotonic()-start
    report['modules']={k:getattr(sys.modules[k],'__file__',None) for k in ['adflow','dafoam','tacs','funtofem','mphys','openmdao'] if k in sys.modules}
    paths=set(line.split()[-1] for line in Path('/proc/self/maps').read_text().splitlines() if any(k in line for k in ['libadflow','libHisaSolver','libDASolver']) and line.split()[-1].startswith('/'))
    if 'adflow' in sys.modules:
        package = Path(sys.modules['adflow'].__file__).parent
        report['installed_adflow_libraries'] = {str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in package.glob('*.so')}
    report['libraries']={x:hashlib.sha256(Path(x).read_bytes()).hexdigest() for x in paths}
    (work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('COMPARISON_REPORT '+json.dumps(report),flush=True)
if not report['completed']:sys.exit(1)
