# Frozen-turbulence matrix-free flow experiment, 2026-11-04

A default-off HiSA option replaces GMRES's iterative matrix products with
central differences of the actual high-order native flow residual. Central
Lax-Friedrichs/laplacian matrices are retained for LUSGS preconditioning and
for the legacy residual normalization. The main linear action is
`A v = frozenDdt*v - (R(W+h*v)-R(W-h*v))/(2*h)` in strong form. Native R is the
negative steady flux balance. There is no extra cell-volume factor.

The flow variables are conservative rho, rhoU, rhoE. SA is still advanced
separately by DAFoam; turbulence is held fixed during the flow solve. This is
not a fully coupled six-variable Newton method. Thermodynamics, boundary
conditions and high-order fluxes are reevaluated for each directional sample.
The pseudo-time diagonal remains frozen throughout each linear solve.

After the linear solve, a scaled steady-flow residual norm is checked with
backtracking (at most 12 trials). Failed steps restore the starting flow state;
SA can still change afterwards. This is a preliminary steady-residual decrease
test, not ADflow's unsteady-residual ANK line search. Positivity checks and the
native temperature bounds protect the tested ideal-gas case. Tests cover only
serial, steady, stationary-mesh L3. No general parallel or other-thermodynamics
validation is claimed. No steady adjoint gradients were computed.

## Checks and results

- Both native and ADR libraries compiled in an isolated directory.
- New-library option-off versus old-library 20-step state/residual differences:
  exactly zero.
- Zero-step residual reproducibility: maximum relative drift 9.63e-13 in the
  initial 3-step check.
- First Krylov direction per outer step: epsilon versus epsilon/2 derivative
  difference at most 3.72e-8 (epsilon=1e-6; this is not a full Jacobian audit).
- All comparisons below start from exactly the same saved state, with initial
  fixed-definition residual L2 3766.8210. GMRES subspace=40, relative tolerance
  0.05, absolute tolerance 1e-12. MF tests cap GMRES at 3 restart cycles/120
  Krylov steps per outer iteration. Baseline retains its previous maxIter=60
  but only needs 3 Krylov steps per iteration.

| Variant | CFL | Final L2 / initial L2 | Seconds | Krylov steps | Accepted flow updates |
|---|---:|---:|---:|---:|---:|
| Existing central/LUSGS | 5 | 1.17553 | 0.72 | 60 | no line search |
| MF + backtracking | 5 | 0.87185 | 8.28 | 252 | 2 / 20 |
| MF + backtracking | 100 | 0.64198 | 65.54 | 2400 | 11 / 20 |
| MF + backtracking | 1000 | 0.72353 | 61.37 | 2400 | 12 / 20 |

Both baseline and MF+CFL5 met the inner tolerance in all 20 solves. Neither
MF+CFL100 nor MF+CFL1000 met it in any of the 20 solves within the chosen cap.
Extra residual evaluations within the MF path: 1770, 14593, 14565 respectively
(includes state restoration and trial evaluations; excludes normal outer/SA
work). The finite-difference direction uses two samples and restores the base
state; the existing GMRES early-exit check also evaluates a product after each
Krylov direction, accounting for much of the cost.

No case reached the requested outer 1e-8 threshold. This short, single-run
comparison does not establish a convergence-time advantage. MF changes both
operator and update acceptance, so their individual contributions are not
isolated. Rejecting the flow step does not freeze subsequent separate SA
updates, explaining why total residuals can still vary after rejection.

## Configuration and reproduction

In flowSolver/GMRES:

```
inviscidJacobian centralLaxFriedrichs;
viscousJacobian laplacian;
preconditioner LUSGS;
nKrylov 40;
maxIter 3;
matrixFree true;
matrixFreeEpsilon 1e-6;
matrixFreeCheck false;
matrixFreeMaxTrials 12;
```

Use `cfl_sweep.py --matrix-free --mf-check --mf-epsilon 1e-6` for the 3-step
consistency check. Use --steps 20 --cfl 5 100 1000 without --mf-check for the
bounded comparison. Source checkpoint:
`dafoam_case/validation/dafoam_L3_central_continue/cfl5_0`, --source-time 1000.

The isolated library directory is
`dafoam_case/validation/hisa_matrix_free/build/lib`. Prepend it to
LD_LIBRARY_PATH after sourcing the workspace DAFoam activate.sh, and pass the
same path to --hisa-lib. System/default libraries were not replaced. Exact
loaded-library paths, hashes, timings and line-search histories are in
results.json. Runtime logs are under dafoam_case/validation/dafoam_L3_mf_*.

Run summarize_matrix_free.py --root dafoam_case/validation --output OUTPUT to
recheck default-path equality and finite-difference/zero-step gates and
regenerate the report. The method remains experimental and default-off.
