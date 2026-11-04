# Adjoint/FD gradients at the L3 residual plateau (2026-11-04)

**Gradient accuracy is not established.** Low-order preconditioning fails to
solve the adjoint adequately. Full-residual LU preconditioning improves it
substantially, but independently recomputed linear residuals remain about
1e-6. The requested central FD at h=1e-4 degrees also depends strongly on
continuation length because its perturbed primal solves remain unconverged.
Do not treat the final FD numbers as an accurate steady-gradient reference.

## Base and variable

Base: `dafoam_L3_direct_second_2000/cfl5_0/2000`, L3 3024 cells, alpha=2 degrees,
Mach 0.85. Native max normalized outer residual is 6.28542e-5; independently
computed raw steady L2 is 4634.7021. High-order AUSMPlusUp/wVanLeer, central LF
approximate primal Jacobian, LUSGS, primal GMRES subspace 40 and fixed CFL 5.

Differentiate CD and CL with respect to alpha in **degrees**, with freestream
speed, thermodynamics and geometry fixed. Force directions include their
explicit alpha dependence through patchV. Plus/minus angles are 2.0001 and
1.9999 degrees. Each perturbed solve, and a zero-perturbation control, starts
from the exact same saved internal state and runs 1000 further iterations.
No turbulence freezing or limited-update/direction diagnostics are enabled.

Adjoint tests do not advance the primal base. They explicitly opt into
`--allow-unconverged-adjoint`; primal acceptance and gradient verification
remain false. A previously unexercised harness omission was corrected:
`initializedRdWTMatrixFree()` is required before creating the matrix-free KSP,
as in DAFoam's normal MPhys workflow. This is a test-driver fix, not a newly
established HiSA or shared-SA solver bug.

## Adjoint solves

Both use the reverse-AD high-order residual transpose operator, GMRES restart
40, max 2000 iterations, relative tolerance 1e-8 and absolute tolerance 1e-13.
The gradient is F_alpha - psi^T R_alpha. Independent verification computes
||J^T psi - F_w|| / ||F_w|| using calcJacTVecProduct with the same state scaling.

| Preconditioner | Function | GMRES iterations | Independent relative linear residual | Gradient per degree |
|---|---|---:|---:|---:|
| Low-order FD matrix, right ASM/ILU(1), RCM, overlap 1 | CD | 2000 | 0.805639 | 0.00730595 |
| Same | CL | 2000 | 0.137013 | 0.10408913 |
| Full-residual FD matrix, right LU | CD | 56 | 1.22509e-6 | 0.01072365 |
| Same | CL | 58 | 9.21142e-7 | 0.09985881 |

The low-order matrix uses calcdRdWT(isPC=1), LF convection with the off-diagonal
tauMC term removed. The full-residual matrix uses calcdRdWT(isPC=0), including
the high-order residual, with the configured state FD increment 1e-7. Both
are colored FD preconditioners; the KSP operator remains reverse AD. The
matrix size is 18144 x 18144, with 2274356 stored nonzeros.

The low-order runs hit the iteration cap (PETSc reason -3). The LU runs report
PETSc reason 2 and monitored relative residuals below 1e-8, but the separate
AD product check gives the values above. This discrepancy is retained, not
hidden by the KSP success flag. Its cause has not been isolated; none passes
the independent 1e-8 gate. The much smaller LU linear errors still do not
establish the accuracy of the steady gradient at a nonzero primal residual.

Low-order solve times are 134.70 s (CD) and 130.50 s (CL), excluding shared
setup. LU takes 182.38 s for CD including first factorization and 8.17 s for
CL with factorization reused; total LU-process wall time is 217.79 s. Full LU
is a small-case diagnostic, not a proposed scalable engineering default.

## Requested finite difference

For each N, use paired histories at the same N:

    derivative(N) = [F(alpha+h, N) - F(alpha-h, N)] / (2h).

| Continuation iterations N | dCD/dalpha | dCL/dalpha |
|---:|---:|---:|
| 20 | 0.00816114 | 0.06599822 |
| 100 | 0.00921596 | 0.10241727 |
| 200 | 0.00990366 | 0.10449337 |
| 500 | 0.00555664 | 0.11753372 |
| 1000 | -0.02270226 | -0.04053179 |

The final positive/negative/zero cases have normalized outer residual maxima
5.7753e-5 / 5.4391e-5 / 5.6602e-5, respectively. They have not converged to
1e-8. Their final raw L2 values are 4848.79 / 4822.93 / 4953.63.

At N=1000, half the sum of the perturbed functions minus the zero-control
function is -7.0753e-6 for CD and +8.6122e-6 for CL. These are comparable to
or larger than the plus-minus differences themselves (-4.5405e-6 and
-8.1064e-6). This is evidence against assuming a clean linear, converged FD
response; it is not a measurement of random roundoff noise alone.

The 200-step numbers happen to be closer to the LU adjoint values. Selecting
that snapshot would be unjustified, because the FD estimate is not stable as
iterations continue. No relative "gradient accuracy percentage" is claimed.
At finite N the FD differentiates the iteration-dependent output, while the
steady adjoint approximates a local implicit derivative at the fixed base;
these are not guaranteed to match away from a steady root.

## Reproduction and artifacts

Use the native OpenFOAM environment, workspace DAFoam Python package, and
prepend `dafoam/build/hisa-diagnostics/lib` and
`dafoam_case/validation/hisa_direction_build/lib` to LD_LIBRARY_PATH.
Set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1. Builds are unchanged from
the prior tests; library provenance is retained in the reports.

From the workspace root, each mode requires fresh corresponding case folders:

```
python dafoam/tests/hisaConnectivity/gradient_plateau.py --root "$PWD" --mode adjoint
python dafoam/tests/hisaConnectivity/gradient_plateau.py --root "$PWD" --mode fd
python dafoam/tests/hisaConnectivity/gradient_plateau.py --root "$PWD" --mode adjoint_full
python dafoam/tests/hisaConnectivity/summarize_gradient_plateau.py \
  --root "$PWD" --output dafoam/tests/hisaConnectivity/gradientAtPlateau
```

Outputs are under `dafoam_case/validation/dafoam_L3_gradient_h1e4`; each adjoint
process has a 600-second cap and each primal FD process a 180-second cap.
Psi, RHS and independent J^T psi arrays are saved there. results.json retains
the full FD history, residual metrics, gradients, run status and provenance.
The summary checks exact internal-state restarts and unchanged adjoint base
states. Original cases are untouched. No AD implementation accuracy or
optimization-ready gradient claim follows from this unsuccessful validation.
