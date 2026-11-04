# Fixed-CFL benefit comparison, 2026-11-04

All seven 500-step tests started from the identical saved L3 state after the
10000-step run. Each recorded initial state matches that checkpoint exactly;
initial fixed-definition residual L2 is 438957.5037118538 in every case.
High-order AUSMPlusUp/wVanLeer, mesh, SA model, GMRES subspace 40, relative
linear tolerance 0.05, and LUSGS preconditioning were retained.

The CFL cap remains 100000. For this diagnostic only, resetPseudo=true and
both CFL increment/decrement factors equal one hold the actual CFL at the
specified value. Log checks confirm constant CFL. This is not an ANK update
implementation and does not replace the default adaptive CFL policy.

| CFL | Final / initial residual L2 | Primal seconds | Total Krylov steps |
|---:|---:|---:|---:|
| 5 | 0.1824 | 11.58 | 1501 |
| 10 | 0.5428 | 11.62 | 1505 |
| 20 | 0.8125 | 13.10 | 2498 |
| 58 | 1.0041 | 16.19 | 3500 |
| 100 | 1.0260 | 15.08 | 3500 |
| 300 | 1.0552 | 14.66 | 3500 |
| 1000 | 1.0536 | 14.79 | 3500 |

Timing is one serial sample per setting, excludes process startup, and is not
a statistical speed benchmark. Krylov counts and residual reductions provide
a more direct explanation of the lower-CFL benefit here.

CFL 5 was continued for another 1500 steps from its 500-step checkpoint.
Total: 2000 steps, 45.85 s primal time, 6001 Krylov steps. Residual L2 became
74874.8237 (0.17057 of the common initial value). The extra 1500 steps reduced
residual L2 by only 6.49%; the final HiSA normalized maximum outer residual was
7.3909e-4. It did not meet the requested 1e-8 threshold.

Conclusion: increasing CFL is not beneficial from this stalled state. Lowering
it improves the residual level and reduces linear work, but does not establish
steady convergence. There is no evidence here to justify unconditional CFL
growth or claim that a complete ANK strategy would or would not converge.
No adjoint/steady aerodynamic gradient test was run.

Reproduce with the workspace DAFoam environment active:

```sh
python tests/hisaConnectivity/cfl_sweep.py \
 --source ../dafoam_case/validation/dafoam_L3_cfl100000_tol1e8 \
 --output NEW_OUTPUT --hisa-repo ../Hisa4DAFoam \
 --hisa-lib ../dafoam_case/validation/hisa_limited_fix/build/lib \
 --cfl 5 10 20 58 100 300 1000 --steps 500
```

The runner saves logs, states, actual CFL ranges, residual metrics and timings.
Use --source-time 500 and --source pointing to a completed case to continue
its checkpoint. Runtime evidence for this test is under dafoam_case/validation:
`dafoam_L3_fixed_cfl`, `dafoam_L3_low_cfl`, `dafoam_L3_cfl5_continue`.
`summary.json` in each run contains full per-case reports; `comparison.json`
and `comparison.png` here consolidate the results.
