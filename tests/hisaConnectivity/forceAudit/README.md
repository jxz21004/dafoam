# 已收敛的一阶 HiSA 壁面力检查

对象：L3、Mach 0.85、迎角2°，`reports/rigid_comparison/hisa/2172`。本次不重新求流场，不计算ADflow或伴随。

## 核心结果

独立程序读取保存的p、U、nut及热物性，通过OpenFOAM几何和梯度算子重新积分252个wing面：压力力 `p Sf`，黏性力 `Sf & [-(mu+rho nut) dev(2 symm grad U)]`。它不调用DAFoam函数计算接口，但共享OpenFOAM几何/梯度，因而不是对几何或离散正确性的完全独立证明。

| 分量 | CD | CL |
|---|---:|---:|
| 压力 | 0.1547303590 | 0.2907723117 |
| 黏性 | 0.0090499690 | -0.0002976310 |
| 总和 | 0.1637803280 | 0.2904746807 |

与DAFoam报告的差分别仅约2.8e-16、1.6e-15。压力贡献占总阻力94.47%。这表明大阻力主要由该流场的壁面压力产生，当前没有证据显示是CL/CD投影、缩放或求和错误。黏性项相对较小，但此次也不能宣称它已物理验证。

归一化采用q=13402.323345 Pa、S=45.5，阻力方向(cos2°,0,sin2°)，升力方向(-sin2°,0,cos2°)。对比原先options中的定义一致。

绝对压力换成p-p∞后，CL/CD变化小于1e-15。wing面积向量之和约(0,-1.881261,0)：半翼根部未包含在wing积分中，y方向不闭合，因此参考压力会改变侧向力；x/z方向闭合，不能用它解释本次升阻力异常。

## 新发现：壁面热边界不一致

HiSA保存场的 `T/wing` 是 `fixedValue uniform 300`（等温壁面），来流温度223.252 K。读取原CGNS边界后，各翼面类型均为22，即`BCWallViscous`。ADflow源码 `src/partitioning/readCGNSGrid.F90` 将该类型映射为`NSWallAdiabatic`，而非等温壁面。

因此此前HiSA与ADflow的热边界条件并未一致。固定300 K是否足以解释阻力差异，尚未做受控求解，不能直接宣称已找到全部根因。下一步应单独将HiSA改为绝热壁面，其他一阶设置保持不变，重新收敛并重复力分解。

壁面Cp范围为[-1.59943,2.60716]；最高值面的位置已列在results.json中，应作为后续压力分布检查对象。没有在本次直接修改该流场或边界条件。

## 复现

激活OpenFOAM v2506环境，从工作区根目录：

```sh
export FOAM_USER_APPBIN=/mdouser/test_aerostructure/dafoam/build/force-audit/bin
wmake dafoam/tests/hisaConnectivity/forceAudit
$FOAM_USER_APPBIN/hisaForceAudit \
  -case /mdouser/test_aerostructure/reports/rigid_comparison/force_audit \
  -time 2172 -lib libcompressibleTurbulenceModels.so
python dafoam/tests/hisaConnectivity/forceAudit/summarize.py \
  --case reports/rigid_comparison/force_audit \
  --reference reports/rigid_comparison/hisa/primal_report.json
```

后处理算例复制原始constant、system和2172目录，保持原始算例不变。faces.csv保存面中心、面积向量、压力、黏度和力分量；results.json保存积分结果；cgns_boundary_types.json保存通过本地CGNS C接口读取的边界枚举值。程序当前针对串行L3/2°审核，参考压力等常数固定于此算例，不作为通用后处理包发布。

## 与 ADflow 的摩擦阻力对比

从已有 `adflow_tight/rigid2_000_vol.cgns` 重载收敛场，不进行流场迭代。仅在固定状态上调用 `getResidual` 刷新壁面应力，再求 `cdp/cdv`。若省略刷新，黏性力缓存未初始化，首次提取结果无效，不能用于比较。刷新后总CD与原报告相差1.34e-10，总CL相差-3.09e-9，通过复现检查。CGNS读取有54条类型转换警告，因此以复现误差界定后处理精度。

| 状态 | 黏性 CD | 压力 CD | 总 CD |
|---|---:|---:|---:|
| ADflow，绝热、已收敛 | 0.003739392 | 0.032450717 | 0.036190109 |
| HiSA一阶，300 K、已收敛 | 0.009049969 | 0.154730359 | 0.163780328 |
| HiSA一阶，绝热、4000步未收敛 | 0.008906065 | 0.154754973 | 0.163661037 |

黏性CD是有效黏性应力合力沿来流方向的投影，作为摩擦阻力对比指标；不把它与局部皮肤摩擦系数Cf混淆。三者均用S=45.5、q=13402.323345、迎角2°。

HiSA已收敛300 K状态的黏性阻力为ADflow的2.42倍，相差0.00531058，确实也偏大。绝热未收敛状态约为2.38倍，尚不能作为最终收敛比值。从已收敛300 K状态与ADflow的总阻力差分解看，约95.8%的差值来自压力项；摩擦偏大值得检查，但无法解释全部阻力差异。

下一步定位摩擦差异，应核对壁面处理、y+、nut及局部剪切分布；HiSA采用nutUSpaldingWallFunction，尚未完成与ADflow壁面处理的等效性验证。离散、热边界及体网格差异仍存在，不直接将2.42倍判定为黏性应力实现错误。

复现后处理：激活本地ADflow，运行`adflow_saved_forces.py --case reports/rigid_comparison/adflow_tight --output <新的输出目录>`。无`solver(ap)`调用；结果归档为`adflow_components.json`。首次未刷新缓存的失败提取保留在工作区`adflow_force_components`，有效提取为`adflow_force_components_refreshed`。
