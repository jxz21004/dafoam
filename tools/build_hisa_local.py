#!/usr/bin/env python3
"""Build native/ADR DAFoam and Python bindings without changing system libraries."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess

p = argparse.ArgumentParser()
p.add_argument('--root', type=Path, default=Path('/home/dafoamuser/dafoam'))
p.add_argument('--hisa-lib', type=Path, required=True)
p.add_argument('--jobs', type=int, default=4)
a = p.parse_args()
repo = Path(__file__).resolve().parents[1]
root = a.root.resolve()
out = repo / 'build/hisa-local'
lib = out / 'lib'
lib.mkdir(parents=True, exist_ok=True)
for name in ('adjoint', 'include', 'pyDASolvers'):
    source = repo / 'src' / name
    target = out / 'src' / name
    for f in source.rglob('*'):
        if f.is_file() and not f.is_symlink() and not any(part.startswith('linux64') or part == 'lnInclude' for part in f.parts):
            dest = target / f.relative_to(source)
            dest.parent.mkdir(parents=True, exist_ok=True)
            data = f.read_bytes().replace(b'\r\n', b'\n')
            if not dest.exists() or dest.read_bytes() != data:
                dest.write_bytes(data)
f = out / 'src/adjoint/Make/files'
s = f.read_text()
s = s.replace('$(DAFOAM_ROOT_PATH)/OpenFOAM/sharedLibs/libDASolver', str(lib / 'libDASolver'))
f.write_text(s)
f = out / 'src/pyDASolvers/Make/options'
f.write_text(f.read_text().replace('LIB_LIBS =', 'LIB_LIBS = -L' + str(lib) + ' '))
env = dict(os.environ, DAFOAM_ROOT_PATH=str(root), OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
for mode, tree in (('', 'OpenFOAM-v2506'), ('ADR', 'OpenFOAM-AD')):
    command = '''source "$1" "WM_AD_MODE=$2"
set -e
export LD_LIBRARY_PATH="$3:$4:$LD_LIBRARY_PATH"
wmakeLnInclude "$5/src/adjoint"
wmake -j"$6" libso "$5/src/adjoint"
cd "$5/src/pyDASolvers"
cythonize -3 pyDASolvers.pyx
wmake -j"$6" libso
'''
    with (out / ('build-' + (mode or 'native') + '.log')).open('w') as log:
        subprocess.run(['bash', '-c', command, 'build', str(root / 'OpenFOAM' / tree / 'etc/bashrc'), mode,
                        str(lib), str(a.hisa_lib.resolve()), str(out), str(a.jobs)], env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    import sysconfig
    destination = repo / 'dafoam/libs' / mode
    destination.mkdir(exist_ok=True)
    shutil.copy2(out / 'src/pyDASolvers' / ('pyDASolvers' + mode + '.so'),
                 destination / ('pyDASolvers' + sysconfig.get_config_var('EXT_SUFFIX')))
subprocess.run(['python', '-m', 'pip', 'install', '--no-deps', '--no-build-isolation',
                '--target', str(out / 'python'), '--upgrade', str(repo)], check=True)
import shlex
(out / 'activate.sh').write_text(
    'export PYTHONPATH=' + shlex.quote(str(out / 'python')) + ':${PYTHONPATH:-}\n'
    'export LD_LIBRARY_PATH=' + shlex.quote(str(lib) + ':' + str(a.hisa_lib.resolve())) + ':${LD_LIBRARY_PATH:-}\n')
print('Installed to', out, flush=True)
