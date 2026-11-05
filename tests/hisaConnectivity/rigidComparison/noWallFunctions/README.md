# 绝热一阶 HiSA：关闭壁面函数

通过DAFoam实际运行选项`primalBC.useWallFunction=False`关闭，而非仅改文件后被运行器重新开启。`run_case.py`新增`--no-wall-functions`；对照运行器新增`--hisa-no-wall-functions`。历史默认设置保持不变，此新算例明确选择关闭。

运行日志确认nut改为`nutLowReWallFunction`、alphat为`fixedValue`；最终保存的wing/nut为uniform 0。`nutLowReWallFunction`虽含WallFunction名称，实际提供零壁面涡黏度，未使用Spalding壁面律。T仍为zeroGradient绝热，U仍为无滑移。未改动共享SA实现。

从此前绝热算例4000步状态启动，再执行4000步。保持一阶AUSM、CFL5、GMRES40、LUSGS、原生阈值1e-8。最终原生外层最大分量1.04991e-6，未收敛到目标。ADflow没有重算。

| 指标 | 绝热、壁面函数开启（未收敛） | 绝热、关闭（未收敛） |
|---|---:|---:|
| 黏性CD | 0.008906065 | 0.008071909 |
| 压力CD | 0.154754973 | 0.154678917 |
| 总CD | 0.163661037 | 0.162750826 |
| CL | 0.293660051 | 0.295199463 |

黏性阻力下降约9.4%，但仍约为已有ADflow黏性CD=0.003739392的2.16倍。两份HiSA结果都尚未达到收敛阈值，不能视为严格收敛对照或认定壁面函数是全部原因。后处理独立积分见summary.json。

复现（激活前述本地环境）：

```sh
python dafoam/tests/hisaConnectivity/rigid_compare.py \
 --root /mdouser/test_aerostructure --solver hisa \
 --hisa-adiabatic --hisa-no-wall-functions \
 --hisa-restart-fields /mdouser/test_aerostructure/reports/rigid_comparison/hisa_adiabatic/4000 \
 --output /mdouser/test_aerostructure/reports/rigid_comparison/hisa_adiabatic_lowre_new
```

本次算例在reports/rigid_comparison/hisa_adiabatic_lowre；日志/tmp/rigid_hisa_lowre.log。原始壁面函数算例未覆盖，未执行伴随。
