# Bounded approximate flow updates (2026-11-04)

Outcome: the bounded-update mechanism works, but **none of these L3 runs
converges to 1e-8**. Keep the option experimental and disabled by default.
This test does not justify replacing the existing iteration method.

## Method

The native and reverse-AD HiSA libraries were rebuilt in the workspace under
`dafoam_case/validation/hisa_limited_update/build/lib`. Installed libraries
were not changed. The existing centralLaxFriedrichs/laplacian approximate
operator and LUSGS preconditioner are retained; matrix-free Jv is disabled.
SA remains segregated and advances after the flow solve, including after a
rejected flow step. No shared SA model change is made.

For the native sign convention R = -div(flux) + sources, freeze the pseudo-time
diagonal D at the start of the solve and test

    G(alpha) = D * alpha * deltaW - R(W0 + alpha * deltaW).

Use the same scaled norm for G and R, with density, momentum and total-energy
scales frozen from W0. Accept if ||G(alpha)|| <= (1 - 1e-4 alpha) ||R(W0)||.
This is a per-step pseudo-time defect check, not a steady convergence test.

Start alpha at min(1, updateLimit / largest scaled conservative change).
Density and energy use local magnitudes with floors 1e-3 times their mean
magnitude. Momentum uses its local magnitude with a floor of 0.1 times the
momentum scale. Also require the relative specific internal-energy change
<= updateLimit, positive density/internal energy and no temperature clipping.
At most four trials, each unsuccessful trial halving alpha, then restore
W0 if all fail. This limits updates to 20% using the stated scales, not
20% relative to arbitrarily small local momentum components.

CFL halves after rejection or alpha < 0.25; it grows by 1.2 only when
alpha >= 0.8 and the scaled steady flow residual decreases by at least 1%.
Otherwise it stays fixed. Clamp to [0.1, 100000]. The old SER controller is
bypassed when limitedUpdate is enabled. Use resetPseudo false to preserve
this feedback between iterations.

    flowSolver
    {
        solver GMRES;
        GMRES
        {
            inviscidJacobian centralLaxFriedrichs;
            viscousJacobian laplacian;
            preconditioner LUSGS;
            nKrylov 40;
            maxIter 60;
            solverTol 1e-12 (1e-12 1e-12 1e-12) 1e-12;
            solverTolRel 0.05 (0.05 0.05 0.05) 0.05;
            matrixFree false;
            limitedUpdate true;
            updateLimit 0.2;
            updateMaxTrials 4;
        }
    }

## Identical-start comparisons

Serial L3, 3024 cells, high-order AUSMPlusUp/wVanLeer and hisaLimited viscous
scheme. Every run starts from `dafoam_L3_central_continue/cfl5_0/1000`;
initial state arrays match exactly and initial raw steady residual L2 is
3766.8210377503738. Each run performs 1000 additional iterations.

| Method | Initial CFL | Actual CFL range | Final raw steady L2 | Final max normalized outer flow residual | Solve seconds | Peak process RSS (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| Original, fixed CFL | 5 | 5 | 5707.52 | 9.229e-5 | 24.05 | 594.90 |
| Bounded update | 5 | 3.69–99.05 | 6358.74 | 5.366e-5 | 36.93 | 596.35 |
| Bounded update | 100 | 3.20–100 | 4434.56 | 6.333e-5 | 36.88 | 596.98 |

Raw steady L2 is evaluated with the same DAFoam residual definition at the
final state. The native normalized outer monitor uses its own normalization;
it must not be equated with raw L2 or the inner GMRES residual. Neither
monitor establishes the requested convergence. The raw L2 ratios to the
common starting state are 1.5152, 1.6881 and 1.1773, respectively.

The bounded runs accept 930/1000 and 922/1000 flow steps, use 3660 and 3697
Krylov directions (baseline 3000), and add 2438 and 2452 residual evaluations.
There are at most six additional residual evaluations per step: one
zero-step consistency check, four trials and one restoration. Roughly half
of accepted steps increase the scaled steady merit (493 and 511). Thus a
reduced local pseudo-time defect alone does not enforce sustained steady
residual reduction. The tests do not isolate whether the remaining obstacle
is the approximate operator, flux reconstruction or segregated turbulence.

The cost is about 1.54 times baseline, with approximately 1.5–2.1 MiB extra
peak process RSS in this small case. RSS is Linux ru_maxrss in the actual
solver Python process, including libraries/PETSc/DAFoam, measured once per
case. It is not an algorithm-only allocation measurement or a large-mesh
memory scaling result. No additional full Jacobian or larger Krylov basis
was introduced. Neither gradients nor parallel/dynamic-mesh operation were
validated in this experiment.

## Reproduction and checks

Activate the workspace DAFoam build, put the experimental HiSA library
directory first in LD_LIBRARY_PATH, and set OPENBLAS_NUM_THREADS=1 and
OMP_NUM_THREADS=1. Run from the workspace root, using fresh output paths:

```sh
python dafoam/tests/hisaConnectivity/cfl_sweep.py \
  --source dafoam_case/validation/dafoam_L3_central_continue/cfl5_0 \
  --source-time 1000 --output dafoam_case/validation/dafoam_L3_limited_compare \
  --hisa-repo Hisa4DAFoam \
  --hisa-lib dafoam_case/validation/hisa_limited_update/build/lib \
  --jacobian centralLaxFriedrichs --limited-update --steps 1000 --cfl 5 100
```

Baseline uses the same command without --limited-update, with --cfl 5 and
output `dafoam_L3_limited_baseline`. Each child has a 180-second timeout.

```sh
python dafoam/tests/hisaConnectivity/summarize_limited_update.py \
  --root dafoam_case/validation \
  --output dafoam/tests/hisaConnectivity/limitedUpdate
```

The audit checks identical initial states, every accepted pseudo-residual
inequality, trial evaluation budgets, every CFL update and propagation,
and zero-step residual drift < 1e-8. With limitedUpdate disabled, the first
20 iterations' GMRES history is identical to the previously validated
matrix-free-disabled run. Native and ADR compilation passed; this does not
constitute a new AD derivative validation. Library hashes and numerical
reports are stored in results.json; full logs/fields remain in the workspace.
