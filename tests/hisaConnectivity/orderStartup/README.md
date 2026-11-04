# First-order solve and second-order continuation (2026-11-04)

**First-order AUSM reaches the native normalized outer-flow tolerance 1e-8.
Switching that solution back to second order does not eliminate stagnation
in the next 2000 iterations.** The staged route costs about twice as much as
the direct second-order comparison without demonstrating a convergence benefit.
This is a restart experiment from an existing L3 flow, not a fresh uniform-flow
startup validation or proof that continuation can never help other cases.

## Conditions and controls

Serial L3, 3024 cells, Mach 0.85, alpha 2 degrees. Common starting checkpoint:
`dafoam_L3_central_continue/cfl5_0/1000`. Both routes start from exactly the
same saved state; the first-to-second restart also preserves its state exactly.

- AUSMPlusUp flux, centralLaxFriedrichs approximate inviscid Jacobian,
  laplacian viscous Jacobian and LUSGS preconditioner.
- GMRES subspace 40, relative tolerance 0.05, absolute tolerance 1e-12.
- CFL fixed at 5, configured cap 100000. No CFL-growth change in this test.
- Native pseudoTol 1e-8 for all five flow components. SA is solved normally.
- First order uses firstOrder reconstruction for rho/U/T; second order uses
  wVanLeer. Other discretization settings are retained.
- Limited updates, matrix-free products, frozen-SA and direction diagnostics
  are disabled. No source-code changes were required for these solves.
- Each run has a 2000-iteration budget and a 180-second process timeout.

Native/ADR libraries came from `hisa_direction_build`, with DAFoam from the
workspace `hisa-diagnostics` build. Provenance hashes are in results.json.

## Results

| Stage | Reported final iteration/time index | Final max normalized outer flow residual | Independent raw steady L2 | Solve seconds |
|---|---:|---:|---:|---:|
| First order from common checkpoint | 1911 | 9.95334e-9 | 0.318596 | 42.99 |
| Second order from first-order solution | 2000 | 6.78851e-5 | 5842.73 | 46.52 |
| Direct second order from common checkpoint | 2000 | 6.28542e-5 | 4634.70 | 46.12 |

The first-order run stops before its budget after meeting the native outer
criterion. It logs 1910 flow solves and final time index 1911 (the control
loop performs the convergence check on the following iteration). The last
200 logged outer residual maxima decrease from about 2.68e-8 to 9.95e-9.
The first-order raw steady residual decreases from 171760.78 to 0.318596,
a ratio of 1.85488e-6. Thus it has met the **normalized outer 1e-8 criterion**,
not a raw residual 1e-8 criterion or a fixed-initial L2-relative 1e-8 criterion.
The generic run_case.py report's stricter raw acceptance flag remains false.
No gradients were validated.

Reconstruction changes the discrete residual equation. At the unchanged
first-order final state, switching to second order gives initial raw steady
L2 168547.56; it is expected to differ greatly from the first-order final
0.318596. Second order reduces it to 5842.73 but stalls: the final 200 outer
residual maxima range from 5.04e-5 to 9.42e-5. Direct second order ranges from
4.93e-5 to 9.16e-5 over the same final window. Neither achieves 1e-8.
Do not compare the residual drop across the discretization switch as if it
were a single unchanged equation.

The staged route takes 89.51 seconds in total, versus 46.12 seconds for direct
second order. GMRES direction totals are 6665 for first order, 6066 for the
second-order continuation and 6000 for direct second order. The direct and
staged paths have different total iteration budgets; both second-order stages
have the same 2000-step budget. Final values at an oscillating plateau are
not evidence of a statistically significant advantage for either path.

## SA checks

| Stage | Last SA linear initial residual | Last SA linear final residual | Independent final SA block L2 |
|---|---:|---:|---:|
| First order | 1.39518e-11 | 6.88839e-14 | 0.00263951 |
| First to second order | 9.73233e-9 | 5.85000e-11 | 9.88004 |
| Direct second order | 6.91962e-9 | 5.34644e-11 | 6.88317 |

All last SA linear solves take one iteration. These linear normalized metrics
and independently reassembled nonlinear residuals have different definitions.
The first-order lift/drag are CL=0.290475 and CD=0.163780; they belong to the
first-order discrete solution and must not substitute for a converged
second-order aerodynamic result.

## Reproduction

Activate native OpenFOAM and the workspace DAFoam Python package. Prepend
`dafoam/build/hisa-diagnostics/lib` and
`dafoam_case/validation/hisa_direction_build/lib` to LD_LIBRARY_PATH.
Use OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1.

Run cfl_sweep.py with:

```
--hisa-repo Hisa4DAFoam
--hisa-lib dafoam_case/validation/hisa_direction_build/lib
--jacobian centralLaxFriedrichs --flux-scheme AUSMPlusUp --steps 2000 --cfl 5
```

Case-specific arguments (outputs under dafoam_case/validation):

1. `--source dafoam_case/validation/dafoam_L3_central_continue/cfl5_0
   --source-time 1000 --reconstruction firstOrder`
   output `dafoam_L3_firstorder_2000`.
2. `--source dafoam_case/validation/dafoam_L3_firstorder_2000/cfl5_0
   --source-time 1911 --reconstruction wVanLeer`
   output `dafoam_L3_first_to_second_2000`.
3. The original common source/time, `--reconstruction wVanLeer`,
   output `dafoam_L3_direct_second_2000`.

Use fresh output directories. Recompute the audit with:

```
python dafoam/tests/hisaConnectivity/summarize_order_startup.py \
  --root dafoam_case/validation \
  --output dafoam/tests/hisaConnectivity/orderStartup
```

The script checks successful execution, exact restart states, separates
native outer-flow and SA linear monitors from independent residual blocks,
and retains results including failures to meet the convergence target.
