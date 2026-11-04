# HiSA SA connectivity compatibility test

`DAStateInfoHisaFoam.C` extends only the HiSA `nuTildaRes` second-ring
connectivity with `U`. High-order face mass flux depends on reconstructed
velocities; SA convection uses that flux. The shared SA model and other
solvers' connectivity remain unchanged. This does not alter the SA equation
or fix primal nonlinear convergence by itself.

Local build (native and reverse AD, including the Python extensions):

```sh
python tools/build_hisa_local.py --hisa-lib /absolute/path/to/patched/hisa/lib
source /home/dafoamuser/dafoam/OpenFOAM/OpenFOAM-v2506/etc/bashrc
source build/hisa-local/activate.sh
```

The builder normalizes Windows line endings in its build copy, writes shared
libraries to `build/hisa-local/lib`, and installs the Python package to
`build/hisa-local/python`. It does not replace the system DAFoam installation.
It uses the existing OpenFOAM, PETSc, and unchanged DAFoam auxiliary turbulence
libraries. The HiSA library must include the previously validated `hisaLimited`
AD parity fix for the supplied test cases. Forward AD is not built by this tool.

Validation (from this repository, with the local environment active):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python tools/validate_hisa_local.py \
  --source-case ../dafoam_case/validation/hisa_limited_fix/connectivity \
  --output ../dafoam_case/validation/dafoam_hisa_scoped \
  --hisa-repo ../Hisa4DAFoam \
  --hisa-lib ../dafoam_case/validation/hisa_limited_fix/build/lib
```

Use a fresh output directory. No old coloring files are copied. The check
verifies the installed Python package path and actual native/ADR DAFoam library
paths, tests five assembled velocity columns against uncolored forward
finite differences, checks native/ADR residual parity, and compares 16 random
AD block projections against central differences over multiple step sizes.
These checks do not establish full-matrix correctness, steady aerodynamic
gradient accuracy, or convergence on L2/L3. Other solver source isolation is
checked separately; a full other-solver runtime regression suite is not implied.


## Validated result (2026-11-04)

Built and installed DAFoam 5.1.1 from clone baseline `3b3795d`, with this HiSA-only
connectivity extension, using OpenFOAM-v2506 (native) and OpenFOAM-AD (ADR).
Both matching Python extensions were rebuilt. The loaded package and library
paths/hashes are recorded in `summary.json`.

- Five checked Ux columns: zero missing SA entries; maximum relative error
  4.3675e-6 against identical uncolored forward differences (step 1e-8).
- Native/ADR residual relative difference: 1.0351e-16.
- Sixteen random block projections: maximum best-step AD/central-FD relative
  error 8.8227e-6 (gate 1e-4).
- Frozen-state native residual unchanged from the previous installation/test.
- L3 (3024 cells), uniform initialization: 20 primal iterations completed,
  primalFail=0, finite residuals, temperature range 197.42–292.65 K. This is an
  execution smoke test, not convergence: raw residual L2 remains 1.2632e6.
- Coloring regenerated: 457 colors; reported matrix assembly CPU time 18.25 s
  (excludes initialization and coloring).

Shared SA, SimpleFoam state information, and RhoSimpleFoam state information
were checked against this clone's HEAD with line endings normalized and are
unchanged. No full runtime regression of other solvers was performed.

## L3 stagnation diagnostics

See [SA refresh, frozen-SA and update-direction checks](stagnationDiagnostics/README.md)
for the follow-up to the unsuccessful bounded-update experiment. The measured
approximate-Jacobian prediction mismatch is a diagnostic finding, not a fix.

The [fixed-direction flux/reconstruction audit](fluxDirection/README.md) further
separates reconstruction changes from the AUSM/LF operator mismatch.

[First-order convergence and second-order continuation](orderStartup/README.md)
records a first-order outer-tolerance success and the subsequent second-order
stagnation, with a direct second-order control.

[Adjoint gradients at the residual plateau](gradientAtPlateau/README.md) compare
alpha sensitivities against central differences with a 1e-4-degree step,
including independent adjoint residual checks and finite-iteration FD drift.
