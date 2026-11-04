# Native HiSA implicit-operator tests, 2026-11-04

All initial 300-step runs start from exactly the same L3 checkpoint, after the
previous CFL=5 500-step run. Initial fixed-definition residual L2: 80068.8245.
The high-order AUSMPlusUp/wVanLeer residual, mesh, turbulence model, GMRES
subspace 40 and relative tolerance 0.05 are unchanged. Only the implicit
operator, its preconditioner matrix, and explicitly named bounds/CFL differ.

| Operator | Global CFL | Bounds | Final L2 / initial L2 | Primal seconds |
|---|---:|---|---:|---:|
| LaxFriedrichs | 5 | off | 0.94866 | 7.02 |
| centralLaxFriedrichs | 5 | off | 0.05573 | 6.78 |
| discreteLaxFriedrichs, fixed | 5 | off | 0.39391 | 9.26 |
| LaxFriedrichs | 100 | on | 5.64415 | see results.json |
| discreteLaxFriedrichs, fixed | 100 | off | abort at iteration 13 | — |
| discreteLaxFriedrichs, fixed | 100 | on | 56.20884 | 11.75 |

All discrete tests use a separate centralLaxFriedrichs/laplacian matrix for
LUSGS, as required by that operator's source documentation. Bounds are the
existing local density/internal-energy bounds, 0.8–1.2 times previous values,
with local pseudo-time damping; this is not residual-based line search.
Global CFL is fixed; when bounds are enabled the local effective CFL varies.

The central operator uses arithmetic face weights (0.5) for this inviscid
linearization, whereas the original LaxFriedrichs uses geometric weights.
The observed benefit does not prove it is the exact high-order Jacobian.
Both the main matrix and LUSGS matrix change in the central test; their
individual contributions have not been separated.

From the central 300-step checkpoint (L2=4462.6159), three further 1000-step
runs gave:

| CFL | Final L2 | Final maximum normalized outer residual | Seconds |
|---:|---:|---:|---:|
| 5 | 3766.8210 | 6.9424e-5 | 23.11 |
| 100 | 24199.3240 | 4.0154e-4 | 23.49 |
| 1000 | 41498.1918 | 7.7220e-4 | 23.98 |

None reached 1e-8. central+CFL5 is the best tested configuration here, with
strong initial improvement but subsequent stagnation. A larger CFL alone is
still not effective. No steady adjoint/gradient acceptance is claimed.

## Confirmed optional-operator compatibility bug

Before repair, discreteLaxFriedrichs aborts in OpenFOAM-v2506 with:
`Operator + : undefined for oriented and unoriented types`.
The geometric and reconstructed interpolation weights were marked oriented.
A face interpolation coefficient should be unoriented. Setting w, wu, we,
and wrho unoriented (under the same OpenFOAM version guard already used in
LaxFriedrichs) eliminates this compatibility failure. Only that branch changes.
Native/ADR builds succeeded; the corrected branch completes 300 steps at CFL5.
Native/ADR frozen-state residual relative difference: 1.2097e-16.
At CFL100 the subsequent abort is a different failure: nonphysical flow states
lead to the thermodynamic temperature inversion failing after large updates.
It is not the original orientation exception. See results.json for both errors.

Patched libraries are isolated under
`dafoam_case/validation/hisa_discrete_fix/build/lib`; default system libraries
and the existing workspace activation script were not replaced.

## Reproduction

With the workspace DAFoam environment active, use operator_sweep.py for the
initial six configurations. Use cfl_sweep.py with --jacobian and
--pc-jacobian for individual cases; --bounding enables the bounded variant.
Use a fresh output directory. For the repaired discrete branch, prepend its
library directory to LD_LIBRARY_PATH and pass that same directory to
--hisa-lib. The runner validates actual HiSA library paths and hashes.

Runtime directories under dafoam_case/validation are named by keys in
results.json. The initial failures and later numerical failures are retained,
not counted as successful converged runs. Timing is a single serial sample
per configuration and excludes initialization; no statistical speedup claim.
