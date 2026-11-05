# 表面异常区域定位

比较已有ADflow收敛表面场与HiSA绝热、关闭壁面函数的一阶4000步状态。HiSA原生残差约1.05e-6，尚未达到1e-8，因此定位结果为诊断，不能视为最终收敛差值。没有重新求解流场。

## 几何匹配检查

两边均252个壁面面元。CGNS坐标重建四边形面积加权面心后，全部一一匹配，最大距离1.543e-5 m。不能使用四顶点简单平均：在该粗网格上它会造成最多0.12 m的面心误差及上下表面错配。

进一步将ADflow Cp投影到HiSA面积向量上积分，压力CD=0.03245071656，与ADflow原生压力CD=0.03245071711一致到5.5e-10。这支持当前面对应、Cp字段读取和压力投影正确。

## 压力差集中在哪里

[表面Cp对比](surface_maps.png)使用相同色标，按流体域面法向z分量分为两行，避免上下表面完全重叠。上排Sf_z<0主要为机翼上表面，下排Sf_z>=0还包含其他朝向的面，不能把它严格称为纯下表面。[逐面Cp对照](cp_comparison.png)给出同一面元的数值对应。

- HiSA Cp范围[-1.6206,2.6051]，ADflow为[-0.7829,1.3588]。
- 压力阻力差约0.122228，最大的单面贡献位于前缘及邻近上下表面，不只表现为某一激波位置的轻微平移。
- 例如前缘面(x,y,z)≈(3.2262,5.9639,0.000164)m，HiSA Cp=2.5430，ADflow Cp=1.0934，单面贡献ΔCD≈0.00912。
- 前缘窄条对应NSWallAdiabaticBCZone1，合计ΔCD≈0.03379，约占总压力差27.6%。邻近主表面Zone7、Zone10分别贡献0.05209和0.03436。
- 绝对阻力差贡献最大的10个面，占所有面绝对贡献总和41.5%。这里是绝对贡献比例，不等同于净阻力差比例。

这将优先排查区域缩小到前缘高压区及其上下游邻近面；应继续检查这些面的单元压力、速度、面通量和局部网格度量，尚不能凭表面图确定是数值耗散、壁面通量还是网格问题。

## 摩擦高值区

[HiSA局部Cf图](hisa_friction_map.png)采用壁面应力的切向模长除以q，不是局部阻力投影。最高点集中在翼尖附近：最大Cf≈0.02141，位置(x,y,z)≈(8.6281,14.0072,-0.01880)m；其次约0.01147、0.01006，同样位于y≈14 m。

高Cf不代表对总阻力贡献最大，还需要乘面面积并投影。已有ADflow表面CGNS只包含速度、Cp和Mach，没有局部剪切或Cf，因此本次**没有生成ADflow与HiSA的逐点摩擦差图**。不能把HiSA翼尖高值直接认定为相对ADflow的异常。已有总CDv对比仍见forceAudit记录。

## 数据与复现

`matched_faces.csv`给出每个HiSA面的位置、匹配距离、两边Cp、ΔCp、该面的ΔCD压力贡献及HiSA Cf。`summary.json`含排名和分区贡献。脚本`compare.py`通过本地CGNS C接口读ADF格式数据，再使用面积加权面心匹配；工作区和依赖路径显式固定为本次环境。

```sh
MPLCONFIGDIR=/tmp/mpl-audit python dafoam/tests/hisaConnectivity/surfaceAudit/compare.py
```

输入：`reports/rigid_comparison/adflow_tight/rigid2_000_surf.cgns`与`reports/rigid_comparison/hisa_adiabatic_lowre/force_audit_faces.csv`。
