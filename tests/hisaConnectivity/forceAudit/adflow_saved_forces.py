#!/usr/bin/env python3
"""Evaluate force components from an existing ADflow CGNS solution; no flow solve."""
import argparse
import json
from pathlib import Path
from adflow import ADFLOW
from baseclasses import AeroProblem
from mpi4py import MPI
p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
case=a.case.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
options=json.loads((case/'options.json').read_text());options['restartFile']=str(case/'rigid2_000_vol.cgns');options['outputDirectory']=str(out)
solver=ADFLOW(options=options,comm=MPI.COMM_WORLD)
ap=AeroProblem(name='rigid2',mach=.85,altitude=10000.,alpha=2.,areaRef=45.5,chordRef=3.25,evalFuncs=['cd','cdp','cdv','cdm','cl','clp','clv'])
solver.setAeroProblem(ap)
# Refresh boundary stresses at the loaded state, without a state update.
solver.getResidual(ap)
funcs={};solver.evalFunctions(ap,funcs)
ref=json.loads((case/'report.json').read_text())['functions']
delta={k:funcs[k]-ref[k] for k in ref}
report={'functions':funcs,'difference_from_saved_report':delta,'flow_iterations_executed':0,'source_case':str(case),'reference_area':45.5,'dynamic_pressure':float(ap.q)}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('SAVED_FORCE_REPORT',json.dumps(report))
assert max(abs(v) for v in delta.values())<1e-8, 'Reloaded forces do not reproduce saved state'
