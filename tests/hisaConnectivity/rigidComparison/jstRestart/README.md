# 一阶收敛解启动 JST，CFL=1

本次仅测试 HiSA，ADflow 不重算。L3、迎角 2°，从 `reports/rigid_comparison/hisa/2172` 的一阶 AUSM 收敛场启动，切换 JST（k2=0.5、k4=0.02、pressure sensor），CFL 从 5 降到 1。其余通量外设置保持不变：GMRES 子空间40、LUSGS、centralLF近似雅可比、原生阈值1e-8。

结果：完成新增4000步，没有原生abort，但未收敛。初始状态向量与一阶最终状态逐项核验，最大绝对差为0，因此不是误从自由来流启动。

| 指标 | 结果 |
|---|---:|
| 新增迭代 | 4000 |
| 求解时间 | 64.93 s |
| 初始 JST 原始稳态 L2 | 1.92226349e6 |
| 最终 JST 原始稳态 L2 | 3.78320622e6 |
| 最终/初始 | 1.96810 |
| 最后原生外层归一化分量最大值 | 0.115607 |
| 最后 CL（未收敛） | 0.145273 |
| 最后 CD（未收敛） | 0.056173 |

一阶 AUSM 收敛不意味着切换离散后的 JST 残差也很小。本次 JST 残差没有下降，不能用较低的 CD 宣称准确性改善。报告中 `primal_fail=0` 只表示接口没有标记失败，日志没有原生收敛消息，独立原始残差检查 `accepted=false`；正常退出不等于收敛。

重启目录包含 `uniform` 元数据，导致日志编号沿用旧计数，从2173到6172，确实执行了4000个新增步骤，并非6172步。最终报告的 `primal_iterations=4000` 是本次最终时间值。运行器现已排除重启目录中的 `uniform` 和 `pseudoCo*`，后续独立试验从零计数并使用明确指定的CFL。本次日志确认实际CFL为1。

这次试验同时改变启动场和CFL，不能区分是哪一项避免了此前的爆炸式退出。能确认的是：两项组合尚不足以使这套JST收敛。下一步应审查该JST耗散尺度、边界处理、传感器及近似雅可比的匹配；如果继续降低CFL，应作为受控诊断，不能把残差下降慢预先当作唯一原因。

复现（使用上级报告的本地环境，新输出目录）：

```sh
python dafoam/tests/hisaConnectivity/rigid_compare.py \
  --root /mdouser/test_aerostructure --solver hisa --hisa-flux JST \
  --hisa-cfl 1 \
  --hisa-restart-fields /mdouser/test_aerostructure/reports/rigid_comparison/hisa/2172 \
  --output /mdouser/test_aerostructure/reports/rigid_comparison/hisa_jst_restart_cfl1_new
```

原始算例 `reports/rigid_comparison/hisa_jst_restart_cfl1`，完整日志 `/tmp/rigid_hisa_jst_restart_cfl1.log`。本目录归档配置、结果和日志摘录。本次没有修改求解器源码，也没有执行伴随测试。
