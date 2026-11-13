# HiSA compatibility and validation

The HiSA-specific `nuTildaRes` second-ring connectivity includes velocity dependencies from reconstructed face mass flux. The shared SA model and other solver connectivity are unchanged by this scoped compatibility fix. It is not itself a primal convergence correction.

Build tooling remains in `tools/build_hisa_local.py`; `tools/validate_hisa_local.py` drives frozen-state connectivity, native/AD parity and derivative checks using the reusable Hisa4DAFoam regression helpers. Supply explicit source case, output and library paths outside this repository. Existing `build/hisa-local` installations remain available and ignored by Git.

Project validation and convergence scans were removed from `tests/hisaConnectivity`. Historical evidence remains in Git history and the cleanup archive. Current RAE2822 baseline and ADflow reference are maintained in `../validation/hisa_adflow_comparison`. New project test scripts and reports must be placed there rather than under DAFoam tests.

## Freestream changes between primal solves

DAHisaFoam supports these `daOptions` defaults, following ADflow:

```python
"infChangeCorrection": True,
"infChangeCorrectionTol": 1e-12,
"infChangeCorrectionType": "offset",  # alternatively "rotate"
```

The current implementation supports characteristic farfield velocity changes
(speed and angle) at fixed freestream pressure/temperature for the ideal-gas,
constant-Cv setup. It does not implement ADflow's entire AeroProblem interface;
thermodynamic changes with correction enabled are rejected explicitly.
Other DAFoam solver classes are unaffected.

Before a native primal solve, compare the updated farfield reference to the
previous reference (initially read from the loaded checkpoint). Offset adds the
freestream velocity increment to cell velocities. Rotate rotates/scales them.
Density is preserved and the freestream total-energy-density increment is
converted consistently into updated cell pressure/temperature. Wall boundary
conditions are reapplied. No random noise is added.

The threshold measures the nondimensional free-stream velocity/energy change,
using velocity scale sqrt(pInf/rhoInf). Unchanged or below-threshold conditions
leave the initial field untouched; `False` disables correction. Native initial
field preparation is outside the residual differentiation path and does not
modify the separate reverse-AD solver state.

After an applied correction, hybrid stopping requires at least
`max(1, primalMinIters)` updates during that call, subject to the existing global
iteration budget. The cumulative iteration count and original residual baseline
are preserved. This startup protection does not replace convergence and gradient
accuracy verification.
