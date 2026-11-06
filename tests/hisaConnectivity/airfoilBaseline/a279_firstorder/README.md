# 同工况一阶重构对照

RAE2822，迎角2.79°、T=255.56 K、Mach0.734、SA，当前HiSA库。只把fvSchemes的`interpolate wVanLeer`改为`interpolate firstOrder`，rho/U/T的reconstruct项沿用该变量；其余网格、来流、边界、GMRES8/LUSGS、CFL初值1上限200、原生1e-8门槛和1250步预算保持一致。两次都从均匀来流启动，不使用高阶末态续算。

| 指标 | wVanLeer | firstOrder |
|---|---:|---:|
| 迭代数 | 1250 | 1250 |
| 最终原生归一化残差最大分量 | 5.39803e-3 | 2.60578e-5 |
| CL | 0.777985707 | 0.707326069 |
| CD | 0.019120276 | 0.046270367 |
| 达到1e-8 | 否 | 否 |

一阶最终残差约低207倍，后期仍下降，没有高阶那样明显的5e-3平台；但本次尚未达到目标收敛条件。相较高阶末步，CL低9.08%、CD高142%（2.42倍）。两者都不是充分收敛结果，不能据此证明谁已具备所需精度。

一阶最后100步CL范围0.70731153—0.70732607，CD范围0.04626925—0.04627037；气动力较稳定仍不代表伴随所需残差达标。此次运行墙钟338秒、CPU329.53秒，非重复性能基准，不据此下阶次耗时结论。

运行目录`reports/hisa_airfoil_baseline/rae2822_a279_m0734_firstorder`，完整日志`/tmp/rae279_firstorder.log`。库SHA、系数历史和配置已归档。`order_comparison.png`为两种重构的直接对照；本次未修改求解器源码。

在本地OpenFOAM环境及当前HiSA库路径下运行独立hisa可复现。全新计算需从0开始，不能直接在已完成的latestTime目录运行并当作初次计算。复现汇总：

```sh
MPLCONFIGDIR=/tmp/mpl-audit python dafoam/tests/hisaConnectivity/airfoilBaseline/summarize.py \
 --case reports/hisa_airfoil_baseline/rae2822_a279_m0734_firstorder \
 --log /tmp/rae279_firstorder.log \
 --output dafoam/tests/hisaConnectivity/airfoilBaseline/a279_firstorder
MPLCONFIGDIR=/tmp/mpl-audit python dafoam/tests/hisaConnectivity/airfoilBaseline/a279_firstorder/compare_orders.py
```
