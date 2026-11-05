# HiSA JST：2° 刚性 L3 试验

本次只运行 HiSA，沿用已有 ADflow 收敛结果。JST 测试失败：约 1441 步时原生 Foam abort，未获得可比较的收敛 CL/CD。

除通量外，保持此前刚性 HiSA 设置：自由来流初始化、Mach 0.85、迎角 2°、L3 3024 单元、CFL 固定 5、GMRES 子空间 40、LUSGS、centralLaxFriedrichs 近似雅可比、原生 pseudoTol=1e-8、最多 4000 步。新增配置：

```foam
fluxScheme JST;
jst_k2 0.5;
jst_k4 0.02;
sensorName pressure;
```

保留 `hisaForceJSTFlux=0`，原始流场和 AD 路径的日志均明确打印 `Selecting flux scheme JST`。原有 reconstruct(rho/U/T) firstOrder 字典项保留，但 JST 主对流通量直接采用线性中心插值，不经过 AUSM 的这些重构项。本次未修改求解器源码，也未运行伴随。

残差在约 100—1300 步的过程中没有持续下降，典型归一化外层分量仍为 1e-2。随后出现巨大负密度，内层 GMRES 达到 60 次循环/2400 个 Krylov 步仍不收敛；1440 步附近 CL/CD 达到非物理的 1e12 量级。1441 步温度反算超过 100 次迭代后自行 abort（包装进程退出码 250）。没有实施人工终止，没有生成最终 primal_report.json。

不能选取中途曾接近 ADflow 的 CD 当作结果：这些状态未收敛，最终发散。该结果只能说明当前 JST 实现配合这套参数无法完成此次测试，不能证明所有 JST 参数均不可用，也不能将温度反算退出本身认定为根本原因。

下一步可独立测试更小 CFL 或从已收敛的一阶状态切换 JST，以区分启动/时间推进稳定性与离散问题；如果仍然停滞，应检查该实现的耗散尺度、压力传感器及预处理匹配，而非只调大迭代上限。

复现（先使用上级报告中的本地环境）：

```sh
python dafoam/tests/hisaConnectivity/rigid_compare.py \
  --root /mdouser/test_aerostructure --solver hisa --hisa-flux JST \
  --output /mdouser/test_aerostructure/reports/rigid_comparison/hisa_jst_new
```

完整日志 `/tmp/rigid_hisa_jst.log`，运行目录 `reports/rigid_comparison/hisa_jst`。本目录归档配置、失败摘要及关键日志，已有 ADflow 结果保持不变。
