# First-order L3 adjoint gradient comparison (2026-11-04)

**At an alpha step of 1e-4 degrees, first-order CD/CL adjoint gradients agree
with central differences to 5.59e-5 and 9.57e-5 relative, respectively.**
Both independently recomputed adjoint linear residuals are below 1e-8.
This is a positive check of these two first-order alpha derivatives at one
state and step size, not a validation of second-order or geometry gradients.

## Configuration

Base: `dafoam_L3_firstorder_2000/cfl5_0/1911`, the first-order AUSM solution
with native normalized flow residual maximum 9.95334e-9. Its raw steady L2
is 0.318596; normalized and raw residual gates are distinct. L3 has 3024 cells,
Mach 0.85 and alpha=2 degrees. All probes retain firstOrder reconstruction
of rho/U/T. Flow solver: centralLaxFriedrichs approximate Jacobian, laplacian
viscous Jacobian, LUSGS, GMRES subspace 40, CFL fixed at 5 (configured cap
100000). SA is solved normally.

FD uses alpha=2 +/- 1e-4 degrees, plus an alpha=2 zero control, all starting
from the same saved internal state. To reduce finite-difference contamination,
the perturbed/zero solves use a tighter native pseudoTol of 1e-10, with a
3000-step budget. They stop after reaching that criterion rather than after
a fixed number of steps. The adjoint base remains the original 1e-8 state.

Adjoint: reverse-AD first-order residual transpose, GMRES restart 40, max
2000 iterations, relative tolerance 1e-8 and absolute tolerance 1e-13. The
preconditioner is LU of the full colored finite-difference residual Jacobian
transpose (calcdRdWT(isPC=0), state FD increment 1e-7), used only as a small-case
diagnostic. It does not replace the primal LUSGS or install a new default.

## Gradients

All derivatives are per degree, including the explicit change in lift/drag
reference directions with alpha.

| Function | Adjoint | Central FD, h=1e-4 deg | Relative difference | Percent difference |
|---|---:|---:|---:|---:|
| CD | 0.00798135029 | 0.00798090428 | 5.58848e-5 | 0.005588% |
| CL | 0.09970684338 | 0.09969730052 | 9.57184e-5 | 0.009572% |

Relative difference is abs(adjoint-FD)/abs(FD). No universal acceptance
threshold or multi-step FD accuracy proof is inferred from this one comparison.

Both adjoint solves converge in four Krylov iterations (PETSc reason 2).
Independent ||J^T psi - F_w||/||F_w|| is 1.54597e-10 for CD and 5.56300e-11
for CL. Unlike the earlier second-order plateau test, KSP success is confirmed
by the independent product check. CD solve time is 63.15 s including initial
LU factorization, and CL takes 0.716 s with the factor reused. Total adjoint
process wall time including initialization/coloring/matrix assembly is 88.43 s.
This cost is not a production scalability benchmark.

## Finite-difference convergence

| Run | Final time index | Raw steady L2 | CD | CL |
|---|---:|---:|---:|---:|
| alpha + h | 989 | 0.00392687 | 0.163781110569 | 0.290484763269 |
| alpha - h | 996 | 0.00359092 | 0.163779514388 | 0.290464823809 |
| alpha | 983 | 0.00388032 | 0.163780312473 | 0.290474793572 |

All three final native normalized flow residual maxima are below 1e-10.
The midpoint of the perturbed outputs minus the zero-control output is
5.85e-12 for CD and -3.29e-11 for CL, much smaller than the perturbation
signals. At matched intermediate steps, the paired derivatives approach the
final result: at steps 200/500/800/980, dCL/dalpha is
0.09932567/0.09957941/0.09967998/0.09969667. This is qualitatively different
from the drifting/sign-changing second-order FD history.

The remaining discrepancy may include finite primal tolerance, the adjoint
base being less converged than the FD endpoints, and finite-difference error.
This experiment does not separate those contributions. Generic driver raw
acceptance/gradient_verified flags are not forcibly set true; the numerical
comparison and its limited scope are recorded separately here.

## Reproduction

Use the native OpenFOAM environment, workspace DAFoam Python package, and
prepend `dafoam/build/hisa-diagnostics/lib` and
`dafoam_case/validation/hisa_direction_build/lib` to LD_LIBRARY_PATH. Set
OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1. No solver libraries changed in
this experiment. Provenance is recorded in results.json.

Common arguments to gradient_plateau.py from the workspace root:

```
--root "$PWD"
--source dafoam_case/validation/dafoam_L3_firstorder_2000/cfl5_0
--source-time 1911
--output dafoam_case/validation/dafoam_L3_first_gradient_h1e4
```

Run once with `--mode fd --steps 3000 --flow-tolerance 1e-10`, and once with
`--mode adjoint_full`. Use fresh mode subdirectories. Process limits are
180 seconds per FD solve and 600 seconds for the adjoint test.

```
python dafoam/tests/hisaConnectivity/summarize_first_gradient.py \
  --root "$PWD" --output dafoam/tests/hisaConnectivity/firstOrderGradient
```

The audit checks exact initial internal states, unchanged adjoint base state,
and every FD endpoint's normalized flow criterion. It retains independent
linear residual checks, FD history, SA linear residuals and library hashes.
Full logs and adjoint/RHS/product arrays remain in the workspace case folders.
