# HiSA compatibility and validation

The HiSA-specific `nuTildaRes` second-ring connectivity includes velocity dependencies from reconstructed face mass flux. The shared SA model and other solver connectivity are unchanged by this scoped compatibility fix. It is not itself a primal convergence correction.

Build tooling remains in `tools/build_hisa_local.py`; `tools/validate_hisa_local.py` drives frozen-state connectivity, native/AD parity and derivative checks using the reusable Hisa4DAFoam regression helpers. Supply explicit source case, output and library paths outside this repository. Existing `build/hisa-local` installations remain available and ignored by Git.

Project validation and convergence scans were removed from `tests/hisaConnectivity`. Historical evidence remains in Git history and the cleanup archive. Current RAE2822 baseline and ADflow reference are maintained in `../validation/hisa_adflow_comparison`. New project test scripts and reports must be placed there rather than under DAFoam tests.
