# 一阶 HiSA 绝热壁面设置及续算

已为独立对照算例设置绝热壁面：`0/T`（或压缩的`T.gz`）中

```foam
wing
{
    type zeroGradient;
}
```

在当前无穿透、无滑移、Fourier导热的静止壁面条件下，温度法向梯度为零对应零法向热通量。不再固定壁温300 K；速度、压力和湍流边界保持原配置。原300 K算例未覆盖。

运行器新增`--hisa-adiabatic`选项，同时支持压缩/未压缩温度文件。它只修改新算例的wing温度边界，不改变来流温度或求解器源码，默认行为仍保留历史配置。

## 试算结果

从已收敛的一阶300 K场启动，L3、迎角2°、AUSMPlusUp一阶、CFL5、GMRES40、LUSGS，新增4000步，71.23秒。

**尚未达到原生1e-8收敛标准**：最终原生归一化外层最大分量约5.077e-6。原始稳态L2为866.04，相对本次初始JST无关的一阶残差16678.54的比值为0.05193。`primal_fail=0`不表示残差合格。

| 指标 | 300 K壁面，已收敛 | 绝热壁面，4000步未收敛 |
|---|---:|---:|
| 压力CD | 0.1547303590 | 0.1547549726 |
| 黏性CD | 0.0090499690 | 0.0089060647 |
| 总CD | 0.1637803280 | 0.1636610373 |
| 总CL | 0.2904746807 | 0.2936600512 |

绝热最后一步的独立壁面力积分再次复现DAFoam输出。目前没有看到巨大阻力被壁温修正消除，但由于新状态未达收敛门槛，不把这张表当作最终收敛差值。壁温不一致是已修正的对照条件问题，尚不能说它就是阻力异常根因。

## 复现

先激活前述本地OpenFOAM/DAFoam/HiSA环境，再运行：

```sh
python dafoam/tests/hisaConnectivity/rigid_compare.py \
  --root /mdouser/test_aerostructure --solver hisa --hisa-adiabatic \
  --hisa-restart-fields /mdouser/test_aerostructure/reports/rigid_comparison/hisa/2172 \
  --output /mdouser/test_aerostructure/reports/rigid_comparison/hisa_adiabatic_new
```

后处理使用`forceAudit`程序，对时间4000运行，再用`summarize.py`积分。结果位于`reports/rigid_comparison/hisa_adiabatic`，完整日志`/tmp/rigid_hisa_adiabatic.log`。本次未重算ADflow，未做伴随验证。
