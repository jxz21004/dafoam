# L3 high-order convergence test, 2026-11-04

Result: **not converged** after 10000 iterations (318.79 s primal wall time).
The requested threshold was 1e-8 for each HiSA normalized outer flow residual;
this is the initial residual of each outer linear system, not its final inner
GMRES residual and not a strictly initial-relative residual.

Used the workspace DAFoam install built from commit 825ad51 and the validated
HiSA native/ADR libraries. L3: 3024 cells, Mach 0.85, altitude 10000 m, alpha 2 deg,
uniform initialization, AUSMPlusUp and wVanLeer reconstruction, unchanged mesh.
The exact numerical settings are included alongside this report.

- CFL ceiling: 100000; observed maximum 78.10, final 57.53. The ceiling does not
  force the adaptive CFL to grow when the residual stagnates.
- GMRES subspace: 40; relative tolerance 0.05, absolute tolerance 1e-12;
  LUSGS preconditioner with LaxFriedrichs inviscid/laplacian viscous approximation.
- Total Krylov steps: 69914 (6.9914 per outer iteration).
- Maximum outer residual: 1.8985e-2, final-100 average 1.8955e-2.
- Fixed-definition residual L2 / initial L2: 2.57296e-2.
- SA last initial/final linear residual: 7.6503e-7 / 1.4161e-8.
- `primalFail=0` only reports a successful solver return here; it does not
  override the failed convergence threshold.

No aerodynamic gradient calculation was attempted. No additional flow-physics
or nonlinear algorithm changes were made during this test.

Full runtime files and checkpoints:
`dafoam_case/validation/dafoam_L3_cfl100000_tol1e8` under the common workspace.
The final state remains available there for future diagnostics.

Run the primal using `Hisa4DAFoam/tests/wingAero/run_case.py --case CASE
--library-dir HISA_LIB` after sourcing the workspace install's activate.sh.
From DAFoam, run `python tests/hisaConnectivity/summarize_l3.py CASE` to regenerate
the summary, CSV, and plot from primal.log and primal_report.json.
