# 一阶 AUSM 壁面通量核对

直接调用当前已编译的 `hisa_direction_build/lib/libHisaSolver.so`，对保存的绝热、关闭壁面函数的一阶状态 `hisa_adiabatic_lowre/4000` 评估通量，不更新流场。该状态尚未达到目标收敛阈值。

## 实测

252个wing面全部满足：

- 壁面质量通量 = 0。
- 壁面对流能量通量 = 0（不代表此程序检查了全部黏性/导热能量项）。
- 动量通量向量与 `p_wall Sf` 的差 = 0。
- 从动量通量反算的压力与p_wall最大误差2.18e-11 Pa。
- p_wall与相邻owner单元压力之差 = 0，符合zeroGradient压力边界。

因此本次没有发现“后处理壁面压力力与原始AUSM壁面动量通量不同”的问题。它也不证明内部面通量、能量方程、网格度量或原生/AD全部残差正确。

源码解释：`reconstructionScheme.C`对非耦合边界的左右状态都取边界值。无滑移壁面U=0，两侧Mach=0，AUSM压力分裂各为1/2；同侧压力相等，质量通量压力修正为零，因此得到零质量通量和p Sf动量通量。`wallFluxAudit.C`调用库而非手工复写该公式，验证了实际计算路径。

## 高压已存在于内部单元

例如最强Cp面（face1，x=0.6901 m、y=1.1994 m）：

- 壁面及owner压力：61414.82 Pa。
- owner温度：253.07 K，密度0.84556 kg/m³。
- owner速度：(31.96,59.85,4.35) m/s。
- owner法向速度约0.0856 m/s；它不是壁面速度，也不表示壁面泄漏。

按该算例Cp=1005 J/(kg K)、R=8314.462618/28.97、来流P/T/V计算，理想气体等熵来流总压约42505.10 Pa。16个壁面面的静压已高于此参考，face1的局部等熵总压估计约63391.47 Pa。这个参考是物理合理性诊断，不是对黏性RANS离散场逐点严格上界的证明；结合未收敛状态，不能仅凭它宣布某条代码有bug。但足以说明异常不只是力积分公式，而应追查内部压力/密度及其离散更新。

下一步应沿前缘异常单元检查相邻内部面质量、动量和能量通量，核对守恒变量与热力学转换，以及总压从远场到前缘的变化。当前未修改求解器，也未进行欧拉滑移壁面对照。

## 复现

激活OpenFOAM v2506后：

```sh
export FOAM_USER_APPBIN=/mdouser/test_aerostructure/dafoam/build/force-audit/bin
export LD_LIBRARY_PATH=/mdouser/test_aerostructure/dafoam_case/validation/hisa_direction_build/lib:$LD_LIBRARY_PATH
wmake dafoam/tests/hisaConnectivity/wallFluxAudit
$FOAM_USER_APPBIN/hisaWallFluxAudit \
 -case /mdouser/test_aerostructure/reports/rigid_comparison/hisa_adiabatic_lowre \
 -time 4000 -lib libcompressibleTurbulenceModels.so
```

程序用保存的p/T/U重建EOS密度和守恒量，对当前一阶通量评估；不等于读取正在运行中的所有内部缓存。编译链接路径显式指向当前验证库。逐面结果见faces.csv，摘要见summary.json，完整日志/tmp/wall_flux_run.log。
