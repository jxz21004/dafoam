# L3 stagnation diagnostics (2026-11-04)

The most direct evidence is a large mismatch between the assembled approximate
flow Jacobian and the true high-order residual derivative along the actual
update. Freezing SA did not remove stagnation. This local experiment identifies
a priority for investigation, not a proven implementation bug or a convergence
fix. All diagnostic options default to false.

## Controlled setup

Serial L3, 3024 cells, Mach 0.85, alpha 2 degrees. Start every case from
`dafoam_L3_central_continue/cfl5_0/1000`, raw residual L2 3766.8210377503738.
High-order AUSMPlusUp/wVanLeer, hisaLimited viscous scheme; centralLaxFriedrichs
approximate Jacobian, laplacian viscous Jacobian, LUSGS, GMRES nKrylov 40,
relative tolerance 0.05, absolute tolerance 1e-12. CFL held at 5, cap 100000.
Matrix-free and limited-update modes are off.

Native and ADR libraries compiled in isolated workspace directories:

- DAFoam: `dafoam/build/hisa-diagnostics/lib`
- HiSA: `dafoam_case/validation/hisa_direction_build/lib`

DAFoam adds only a HiSA-specific `hisaDiagnostics` dictionary in fvSolution:

```
hisaDiagnostics { freezeSA false; auditSA true; }
```

`freezeSA` skips the segregated SA correction; it is a diagnostic, not a
production turbulence closure. `auditSA` samples only the last requested
iteration, after the SA correction, so its refreshing of fluxes cannot change
the preceding trajectory. Shared SA source is unchanged. The native tests
use the existing workspace Python bindings. Library hashes are in results.json.
The SA-only runs loaded the previous unmodified ADR library; both diagnostic
ADR libraries subsequently compiled successfully, and the direction runs used
the completed pair. No new AD derivative or adjoint test was performed.

## SA equation and refresh audit

`DASpalartAllmaras::correct` and `calcResiduals(isPC=0)` share assembly, including
relaxation. The former solves the assembled system and updates nuTilda, while
the latter reassembles the nonlinear equation at the supplied state. The
linear solver's normalized initial/final residual is not the per-volume
nonlinear residual returned by `getResiduals`. A small linear residual alone
cannot prove the nonlinear SA residual has vanished.

In the ordinary HiSA path, mass flux is computed before the conservative flow
update; after updateFields, the SA solve uses this lagged flux. DAFoam's
getResiduals refreshes boundary/intermediate fields and fluxes. At the end of
20 steps, without another SA solve:

| SA residual definition at the same nuTilda | L2 |
|---|---:|
| Reassembled with lagged primal flux | 7.24840559 |
| After refreshing flow flux | 7.41975286 |
| After full DAFoam boundary/intermediate refresh | 7.41975286 |

The flux refresh changes the residual vector by L2 1.18399518; the additional
boundary refresh changes it by only 1.55e-11. Thus the refresh difference is
real, but the lagged-flux nonlinear residual is also nonzero. It does not
explain away the observed residual. Lagged segregated coupling is not by
itself a coding error; at a converged state the lag must disappear.
The audited final state and residual exactly match the prior ordinary 20-step
run, verifying this end-only diagnostic did not alter its result.

## Frozen-SA comparison

| 300 further iterations | Final raw residual L2 | SA block L2 | Solve time |
|---|---:|---:|---:|
| Normal SA updates | 4807.30887 | 7.19585 | 7.38 s |
| Frozen nuTilda | 5451.33082 | 16.48561 | 7.17 s |

The frozen nuTilda array is bitwise unchanged. Both totals exceed the initial
3766.82; freezing SA does not eliminate the flow stagnation on this interval.
The SA residual in a frozen experiment is not a convergence target, since
that equation is intentionally not solved. This result does not establish
that SA coupling is irrelevant in every regime.

## Actual update direction

Enable `directionAudit true` in flowSolver/GMRES. In the ordinary approximate
solve, save the base state W0 and the computed update deltaW. Apply the actual
assembled matrix A to deltaW, then infer Japprox*deltaW = D*deltaW - A*deltaW,
where D is the frozen pseudo-time diagonal and R has native sign -div(flux).
Compare with central differences of R(W0 +/- h*deltaW), with turbulence frozen.
Use density, momentum and energy scaling identical to the limited-update
experiment. This norm differs from the aggregate raw residual above.

For the first of three steps:

| Fraction of actual update | Predicted residual norm / initial | Actual / initial |
|---|---:|---:|
| 1 | 0.371591 | 0.968040 |
| 0.5 | 0.661206 | 0.981427 |
| 0.1 | 0.930776 | 0.995891 |
| 0.01 | 0.993056 | 0.999580 |

The relative directional derivative error
`||Japprox*deltaW - Jtrue*deltaW|| / ||Jtrue*deltaW||` is 4.69282, 4.59156,
and 4.29388 for the three steps. For h=1e-3 and h=5e-4 these reported errors
agree within 1e-6 relative. Zero-step residual drift is below 1e-12. Thus the
mismatch persists close to W0; finite full-step nonlinearity alone cannot
explain it. These checks concern only the actual update directions at three
nearby states, not the entire Jacobian or all operating conditions.

Every diagnostic restores the solved conservative state and original fluxes,
residuals and bounding flag before the ordinary update continues. The final
state and residual exactly match the corresponding three-step run with the
diagnostic disabled. The probe is not a new iteration scheme.

The approximate operator is intentionally not an exact high-order Jacobian.
The measured mismatch warrants isolating convection/reconstruction, viscous
and boundary contributions next. It does not prove an algebraic sign bug,
that LUSGS is unsuitable, or that full matrix-free NK would be cost-effective.

## Reproduce and audit

Source the native OpenFOAM environment and `dafoam/build/hisa-local/activate.sh`.
Prepend the diagnostic DAFoam and desired HiSA library directories to
LD_LIBRARY_PATH. Limit BLAS/OpenMP to one thread. Wait for all library builds
to finish before running (an initial direction trial failed to import an
incomplete library during linking; it produced no numerical result and is
excluded).

Use `cfl_sweep.py` with the common source/time above, --jacobian
centralLaxFriedrichs and --cfl 5, plus:

- `--audit-sa --steps 20`: `dafoam_L3_sa_audit`
- `--audit-sa --steps 300`: `dafoam_L3_sa_control`
- `--audit-sa --freeze-sa --steps 300`: `dafoam_L3_sa_frozen`
- `--direction-audit --steps 3`: `dafoam_L3_direction_ready`
- `--steps 3`: `dafoam_L3_direction_control`

The first three used the existing HiSA limited-update library with its options
disabled; the last two used hisa_direction_build. All outputs are fresh
subdirectories of dafoam_case/validation. Retain full logs there. Audit with:

```
python dafoam/tests/hisaConnectivity/summarize_stagnation.py \
  --root dafoam_case/validation \
  --output dafoam/tests/hisaConnectivity/stagnationDiagnostics
```
