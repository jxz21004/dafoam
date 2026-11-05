# L3、迎角 2° 刚性流场对照

2026-11-05：重新从自由来流启动纯气动求解，没有 TACS、MELD 或结构位移。ADflow 使用本地 2.13.1 (`8155e981`)；HiSA 使用前期验证的一阶配置。两者各自达到原生停止条件，但气动力差异显著，说明差异在结构耦合之前就已存在。

## 条件

共同条件：L3（3024 单元）、Mach 0.85、迎角 2°、海拔 10000 m、参考面积 45.5、弦长 3.25、单 MPI rank、BLAS/OMP 单线程。来流 P=26499.8979 Pa、T=223.25209 K、V=254.60282 m/s、rho=0.41350775 kg/m³、q=13402.32335 Pa。

ADflow 直接读 `mphys_input_files/wing_vol_L3.cgns`。HiSA 使用该 CGNS 转换/修复后的 Foam 网格，未施加结构变形；这不意味着两份体网格坐标完全一致。

ADflow 选项从用户 `asopt2.py` 的 `aero_options` 字典提取，保持 RANS/SA、ANK、中心通量加标量耗散等设置，只将 `L2ConvergenceRel` 从 1e-4 收紧为 1e-8。最大 10000 步，绝对阈值仍为 1e-14。`anksubspacesize=-1` 沿用脚本的自动设置，不把它报告为固定 40。完整显式选项见 `adflow_tight_options.json`。

HiSA 为一阶 AUSM（rho/U/T 的 `firstOrder` 重构），GMRES 子空间 40、LUSGS 预处理、centralLaxFriedrichs 对流近似雅可比和 laplacian 黏性近似雅可比。线性相对容差 0.05，CFL 固定 5、配置上限 100000；最大 4000 步、原生 `pseudoTol=1e-8`。`hisa_options.json` 中 `adjEqnOption.gmresRestart=100` 是未执行的伴随配置，不是流场 GMRES 子空间。

## 结果

| 指标 | ADflow | 一阶 HiSA |
|---|---:|---:|
| 原生外层迭代 | 60 | 2171 |
| 求解计时（s） | 9.35 | 39.42 |
| CL | 0.41516875645 | 0.29047468070 |
| CD | 0.03619010882 | 0.16378032800 |
| 各自原生停止条件 | 达到 | 达到 |

HiSA 的 CL 低约 30.03%，CD 为 ADflow 的 4.526 倍。计时来自各自求解调用，不含初始化；此次是单次运行，存在短暂同时执行，不作为严格性能基准。ADflow 日志 `Iter Tot=436` 与外层 `Iter=60` 是不同列，不将 436 作为外层步数。HiSA 报告 `primal_iterations=2172` 来自最终时间索引，实际求解日志明确为 2171 步。

残差定义分别为：

- ADflow `getResNorms()`：初始 378899.8357，最终 5.89609e-5，最终/初始 **1.55611e-10**。设置阈值为 1e-8，最后一步越过该阈值。
- HiSA 原生最后外层归一化分量最大值 **9.93556e-9**，满足原生 1e-8。另行重算的 DAFoam 原始稳态 L2 为 0.406094，初始为 1.70604e7，最终/初始 **2.38033e-8**。
- HiSA 原始分量最大值为 0.0944925，未达到测试脚本额外的“原始残差 max<1e-8”门槛，所以 JSON 的 `accepted=false`。这不等于原生停止失败，也不应将原生 1e-8 宣称为原始绝对残差 1e-8。

原脚本阈值 1e-4 的 ADflow 试跑在相对残差 9.22631e-5 时返回：CL=0.42053975、CD=0.03575996。收紧后仍有可见变化，因此本表使用收紧后的结果。不能拿先前静气弹平衡态系数与本表的刚性系数直接相减归因于求解器。

## 解释与后续定位

本次排除了结构变形作为气动力差异的必要条件，但不是同离散方案的精度验证：HiSA 是一阶，而 ADflow 主流场使用中心通量加耗散；ADflow 的 SA 离散选项为 `first order`，不能因此把整个 ADflow 流场称为一阶。粗网格上一阶耗散是应优先检查的因素，尚不能定量断言它解释全部阻力差异。

壁面及 SA 配置也尚未严格对齐：HiSA 的 nut 使用 `nutUSpaldingWallFunction`，nuTilda 来流取分子运动黏度；ADflow 打印 `eddyVisInfRatio=0.009`、`useft2SA=True`。这些设置名称不能直接互换，后续应逐项核对模型及边界条件、压力/黏性阻力分解和 Cp 分布，再判断离散或实现问题。本次没有做新的伴随梯度测试。

## 复现与证据

运行器为 `../rigid_compare.py`。先按 `../aerostructComparison/README.md` 激活本地 OpenFOAM、DAFoam、ADflow 与 HiSA 库，在工作区根目录执行：

```sh
python dafoam/tests/hisaConnectivity/rigid_compare.py \
  --root /mdouser/test_aerostructure --solver adflow \
  --output /mdouser/test_aerostructure/reports/rigid_comparison/adflow_new
python dafoam/tests/hisaConnectivity/rigid_compare.py \
  --root /mdouser/test_aerostructure --solver hisa \
  --output /mdouser/test_aerostructure/reports/rigid_comparison/hisa_new
```

必须使用新的输出目录；ADflow 默认相对阈值为 1e-8，可用 `--adflow-relative-tolerance` 显式指定。HiSA 使用现有 `Hisa4DAFoam/tests/wingAero/run_case.py`，不加 restart，从自由来流初始化。

结果归档：`results.json`、`adflow_loose_report.json`、两份收敛日志以及显式配置 JSON。完整场文件保存在工作区 `reports/rigid_comparison/{adflow,adflow_tight,hisa}`；完整日志位于 `/tmp/rigid_{adflow,adflow_tight,hisa}.log`。
