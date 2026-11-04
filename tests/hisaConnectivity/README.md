# HiSA SA connectivity compatibility test

`DAStateInfoHisaFoam.C` extends only the HiSA `nuTildaRes` second-ring
connectivity with `U`. High-order face mass flux depends on reconstructed
velocities; SA convection uses that flux. The shared SA model and other
solvers' connectivity remain unchanged. This does not alter the SA equation
or fix primal nonlinear convergence by itself.

Local build (native and reverse AD, including the Python extensions):

```sh
python tools/build_hisa_local.py --hisa-lib /absolute/path/to/patched/hisa/lib
source /home/dafoamuser/dafoam/OpenFOAM/OpenFOAM-v2506/etc/bashrc
source build/hisa-local/activate.sh
```

The builder normalizes Windows line endings in its build copy, writes shared
libraries to `build/hisa-local/lib`, and installs the Python package to
`build/hisa-local/python`. It does not replace the system DAFoam installation.
It uses the existing OpenFOAM, PETSc, and unchanged DAFoam auxiliary turbulence
libraries. The HiSA library must include the previously validated `hisaLimited`
AD parity fix for the supplied test cases. Forward AD is not built by this tool.

Validation (from this repository, with the local environment active):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python tools/validate_hisa_local.py \
  --source-case ../dafoam_case/validation/hisa_limited_fix/connectivity \
  --output ../dafoam_case/validation/dafoam_hisa_scoped \
  --hisa-repo ../Hisa4DAFoam \
  --hisa-lib ../dafoam_case/validation/hisa_limited_fix/build/lib
```

Use a fresh output directory. No old coloring files are copied. The check
verifies the installed Python package path and actual native/ADR DAFoam library
paths, tests five assembled velocity columns against uncolored forward
finite differences, checks native/ADR residual parity, and compares 16 random
AD block projections against central differences over multiple step sizes.
These checks do not establish full-matrix correctness, steady aerodynamic
gradient accuracy, or convergence on L2/L3. Other solver source isolation is
checked separately; a full other-solver runtime regression suite is not implied.


## Validated result (2026-11-04)

Built and installed DAFoam 5.1.1 from clone baseline `3b3795d`, with this HiSA-only
connectivity extension, using OpenFOAM-v2506 (native) and OpenFOAM-AD (ADR).
Both matching Python extensions were rebuilt. The loaded package and library
paths/hashes are recorded in `summary.json`.

- Five checked Ux columns: zero missing SA entries; maximum relative error
  4.3675e-6 against identical uncolored forward differences (step 1e-8).
- Native/ADR residual relative difference: 1.0351e-16.
- Sixteen random block projections: maximum best-step AD/central-FD relative
  error 8.8227e-6 (gate 1e-4).
- Frozen-state native residual unchanged from the previous installation/test.
- L3 (3024 cells), uniform initialization: 20 primal iterations completed,
  primalFail=0, finite residuals, temperature range 197.42–292.65 K. This is an
  execution smoke test, not convergence: raw residual L2 remains 1.2632e6.
- Coloring regenerated: 457 colors; reported matrix assembly CPU time 18.25 s
  (excludes initialization and coloring).

Shared SA, SimpleFoam state information, and RhoSimpleFoam state information
were checked against this clone's HEAD with line endings normalized and are
unchanged. No full runtime regression of other solvers was performed.

## L3 stagnation diagnostics

See [SA refresh, frozen-SA and update-direction checks](stagnationDiagnostics/README.md)
for the follow-up to the unsuccessful bounded-update experiment. The measured
approximate-Jacobian prediction mismatch is a diagnostic finding, not a fix.

The [fixed-direction flux/reconstruction audit](fluxDirection/README.md) further
separates reconstruction changes from the AUSM/LF operator mismatch.

[First-order convergence and second-order continuation](orderStartup/README.md)
records a first-order outer-tolerance success and the subsequent second-order
stagnation, with a direct second-order control.

[Adjoint gradients at the residual plateau](gradientAtPlateau/README.md) compare
alpha sensitivities against central differences with a 1e-4-degree step,
including independent adjoint residual checks and finite-iteration FD drift.

[First-order alpha-gradient check](firstOrderGradient/README.md) gives a positive
comparison at h=1e-4 degrees, with independently verified adjoint linear
residuals; this does not validate the stalled second-order gradients.

## 一阶与二阶的配置切换

切换位置是 **OpenFOAM 算例的 `system/fvSchemes`**，使用同一个 Python
求解入口。当前测试中，一阶配置为：

```foam
interpolationSchemes
{
    default          linear;
    reconstruct(rho) firstOrder;
    reconstruct(U)   firstOrder;
    reconstruct(T)   firstOrder;
}
```

二阶将这三行的 `firstOrder` 全部换成 `wVanLeer`。这里的阶数指流场
对流项的重构：一阶使用单元两侧的常值，`wVanLeer` 在光滑区域采用名义
二阶限制重构，在受限区域可能降阶；并不是所有方程一起切换阶数。

也可以通过本目录的 `cfl_sweep.py` 设置：

```sh
--reconstruction firstOrder   # 一阶
--reconstruction wVanLeer     # 二阶
```

这两个参数只是让脚本修改上述三行。切换本身不需要修改 Python 求解逻辑
或 `fvSolution`，也不改变 AUSMPlusUp 通量公式、黏性离散、SA 方程、LUSGS
和 GMRES 设置。伴随仍使用同一套 AD 代码，对当前所选重构的残差求导。

**库依赖：`firstOrder` 是本工作区在 Hisa4DAFoam 中新增的重构选项**
（HiSA 提交 `5e95626`），不能假定原版 OpenFOAM/HiSA 已支持。
当前测试加载的库位于工作区
`dafoam_case/validation/hisa_direction_build/lib`。换环境时须编译并加载包含
该选项的 HiSA 库；运行 DAFoam 伴随测试时，普通版和 ADR 版均需对应更新。
本工作区没有替换系统安装库，因此只修改字典、不加载相应库可能报
`Unknown discretisation scheme firstOrder`。

梯度测试中的容差是另一项独立设置：一阶有限差分测试通过
`gradient_plateau.py --flow-tolerance 1e-10` 收紧了 `fvSolution` 的
`pseudoTime/pseudoTol`；二阶测试的目标为 `1e-8`，但出现停滞。
收紧容差是为了减少有限差分污染，不是切换重构阶数的必要条件。
