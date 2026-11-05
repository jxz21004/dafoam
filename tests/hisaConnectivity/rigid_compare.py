#!/usr/bin/env python3
"""Run an isolated alpha=2 L3 rigid case with the existing solver settings."""
import argparse
import ast
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--solver', choices=['adflow', 'hisa'], required=True)
parser.add_argument('--hisa-flux', choices=['AUSMPlusUp', 'JST'], default='AUSMPlusUp')
parser.add_argument('--adflow-relative-tolerance', type=float, default=1e-8)
parser.add_argument('--hisa-restart-fields', type=Path)
parser.add_argument('--hisa-adiabatic', action='store_true')
parser.add_argument('--hisa-cfl', type=float, default=5.)
args = parser.parse_args()
if args.hisa_cfl <= 0:
    parser.error('--hisa-cfl must be positive')
root = args.root.resolve()
work = args.output.resolve()
work.mkdir(parents=True, exist_ok=False)

if args.solver == 'hisa':
    template = root / 'dafoam_case/validation/dafoam_L3_firstorder_2000/cfl5_0'
    for folder in ['constant', 'system']:
        shutil.copytree(template / folder, work / folder)
    fields = args.hisa_restart_fields.resolve() if args.hisa_restart_fields else template / '0'
    shutil.copytree(fields, work / '0',
                    ignore=shutil.ignore_patterns('uniform', 'pseudoCo*')
                    if args.hisa_restart_fields else None)
    if args.hisa_adiabatic:
        path = work / '0/T'
        compressed = not path.exists()
        if compressed:
            path = path.with_suffix('.gz')
        opener = gzip.open if compressed else open
        with opener(path, 'rt') as stream:
            text = stream.read()
        text, count = re.subn(r'(?m)(^\s*wing\s*\n\s*\{)[^}]*}',
                              r'\1\n        type zeroGradient;\n    }', text)
        assert count == 1
        with opener(path, 'wt') as stream:
            stream.write(text)
    path = work / 'system/fvSolution'
    text, count = re.subn(r'(?m)^(\s*pseudoCoNum\s+)[^;]+;',
                          lambda m: m[1] + str(args.hisa_cfl) + ';', path.read_text())
    assert count == 1
    path.write_text(text)
    path = work / 'system/controlDict'
    text = path.read_text()
    for key in ['endTime', 'writeInterval']:
        text, count = re.subn(r'(?m)^(\s*' + key + r'\s+)[^;]+;', r'\g<1>4000;', text)
        assert count == 1
    path.write_text(text)
    if args.hisa_flux == 'JST':
        path = work / 'system/fvSchemes'
        text, count = re.subn(r'(?m)^fluxScheme\s+[^;]+;',
                              'fluxScheme JST;\njst_k2 0.5;\njst_k4 0.02;\nsensorName pressure;',
                              path.read_text())
        assert count == 1
        path.write_text(text)
    cmd = [sys.executable, str(root / 'Hisa4DAFoam/tests/wingAero/run_case.py'),
           '--case', str(work), '--library-dir',
           str(root / 'dafoam_case/validation/hisa_direction_build/lib'), '--alpha', '2']
    if args.hisa_restart_fields:
        cmd.append('--restart')
    (work / 'run_config.json').write_text(json.dumps({
        'flux': args.hisa_flux, 'cfl': args.hisa_cfl,
        'adiabatic_wall': args.hisa_adiabatic,
        'restart_fields': str(fields), 'restart': bool(args.hisa_restart_fields),
        'command': cmd}, indent=2) + '\n')
    result = subprocess.run(cmd, check=False)
    sys.exit(result.returncode)

# Extract only the original option dictionary, without initializing TACS/MELD.
source = root / 'asopt2.py'
tree = ast.parse(source.read_text())
assignments = [node for node in ast.walk(tree) if isinstance(node, ast.Assign)
               and any(isinstance(t, ast.Name) and t.id == 'aero_options' for t in node.targets)]
assert len(assignments) == 1
scope = {'inputdir': str(root / 'mphys_input_files') + '/',
         'outputdir': str(work) + '/', 'args': argparse.Namespace(level='L3'),
         'forcesAsTractions': False}
options = eval(compile(ast.Expression(assignments[0].value), str(source), 'eval'), scope)
options['L2ConvergenceRel'] = args.adflow_relative_tolerance
(work / 'options.json').write_text(json.dumps(options, indent=2) + '\n')
os.chdir(work)
from mpi4py import MPI
from adflow import ADFLOW
import adflow
from baseclasses import AeroProblem

assert MPI.COMM_WORLD.size == 1
ap = AeroProblem(name='rigid2', mach=.85, altitude=10000., alpha=2.,
                 areaRef=45.5, chordRef=3.25, evalFuncs=['cl', 'cd'])
start = time.monotonic()
solver = ADFLOW(options=options, comm=MPI.COMM_WORLD)
setup_seconds = time.monotonic() - start
start = time.monotonic()
solver(ap)
solve_seconds = time.monotonic() - start
funcs = {}
solver.evalFunctions(ap, funcs)
norms = [float(v) for v in solver.getResNorms()]
report = {'solver': 'adflow', 'alpha': 2., 'rigid': True, 'mpi_ranks': 1,
          'functions': funcs, 'setup_seconds': setup_seconds, 'solve_seconds': solve_seconds,
          'resnorm_initial_start_final': norms,
          'final_over_initial': norms[2] / norms[0],
          'solve_failed': bool(ap.solveFailed), 'fatal_fail': bool(ap.fatalFail),
          'freestream': {k: float(getattr(ap, k)) for k in ['P', 'T', 'rho', 'V', 'mu', 'nu', 'q']},
          'module': adflow.__file__, 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest()}
(work / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
print('RIGID_REPORT ' + json.dumps(report), flush=True)
if ap.solveFailed or ap.fatalFail:
    sys.exit(1)
