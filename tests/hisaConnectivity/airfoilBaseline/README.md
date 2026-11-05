# 当前 HiSA：RAE2822 独立翼型简单测试

按用户要求使用当前修改版HiSA，没有切换原始上游版本。采用仓库自带`examples/rae2822`算例，在独立目录复制并运行，不调用DAFoam求解接口，不经过此前CGNS转换/修复程序。

## 实测结论

1250步正常完成，未发生发散退出，CL/CD趋稳，但原生残差停滞于约5.57e-3，**未达到1e-8收敛阈值**。没有与实验Cp或参考解进行误差核对，因此本次只是运行/稳定性检查，不能称为精度验证通过。

| 项目 | 结果 |
|---|---:|
| 网格单元数 | 23552 |
| 最终CL | 0.690557 |
| 最终CD | 0.01383583 |
| 最终原生外层归一化最大分量 | 0.00556942 |
| 迭代步数 | 1250 |
| 程序报告CPU/墙钟时间 | 156.79 / 165 s |
| 最后100步CL范围 | 0.690514—0.6905725 |
| 最后100步CD范围 | 0.01376274—0.01386207 |

曲线见`convergence.png`；原始力系数表中Cd(f)、Cd(r)是前/后分解，不作为摩擦/压力分量使用。这里仅提取总Cd和Cl。

## 配置与隔离

- 运行器：`/home/dafoamuser/dafoam/OpenFOAM/sharedBins/hisa`。
- LD_LIBRARY_PATH优先加载当前`dafoam_case/validation/hisa_direction_build/lib/libHisaSolver.so`；ldd确认绑定路径，库SHA256见summary.json。
- 使用原算例来流U=(233.44,0,9.417) m/s，p=108987.773 Pa，T=255.556 K；约Mach0.73、迎角2.31°。
- AUSMPlusUp，lowMachAusm=false，wVanLeer高阶重构；SA湍流模型。
- GMRES子空间8，LaxFriedrichs近似雅可比，LUSGS，线性相对阈值0.1；CFL自适应，初值1、上限200。
- 原算例的characteristicFarfield、boundaryCorrectedFixedValue和characteristicWallTemperature等边界保持不变。这与此前机翼的简单边界配置不同，但本次未做单变量试验，不能把差异直接归因于某一边界。
- 参考弦长及面积均为0.3048（原算例定义）。

## 网格准备与限制

原算例已提供raetaf.msh，使用gmshToFoam、stitchMesh、changeDictionary、transformPoints按(0.3048,1,0.3048)缩放，再renumberMesh。脚本/字典在副本中仅规范化CRLF。

网格边界改变后生成的旧0/meshPhi会阻止OpenFOAM读取；按原setupMesh中删除该文件的步骤，本次将它移到`meshPhi.beforeBoundaryChange`留存，再运行renumberMesh/checkMesh。没有修改网格点或做质量修复。

checkMesh未全部通过：65个高长宽比单元（最大约6.40e5），319个严重非正交面（最大约89.986°）。具体见mesh_check.log。不能把这次运行描述为优质网格基准，也不能据此排除网格对残差停滞的影响。

## 复现与文件

算例位于`reports/hisa_airfoil_baseline/current_rae2822`。完成网格准备后，simulation/constant/polyMesh复制自mesh，0复制自0.org。激活OpenFOam v2506，设置：

```sh
export LD_LIBRARY_PATH=/mdouser/test_aerostructure/dafoam_case/validation/hisa_direction_build/lib:$LD_LIBRARY_PATH
/home/dafoamuser/dafoam/OpenFOAM/sharedBins/hisa \
 -case /mdouser/test_aerostructure/reports/hisa_airfoil_baseline/current_rae2822/simulation
```

从已有完成目录直接运行会沿latestTime重启；全新复现应重新复制原始算例及0.org。日志为`/tmp/rae_current.log`，结果与后处理在simulation中。归档了fvSchemes/fvSolution/controlDict、系数曲线、网格检查与摘要。`summarize.py`读取上述本地路径生成报告。

这次测试表明当前HiSA在独立翼型算例上可稳定推进，但同样有残差平台。下一步应先核对其Cp/参考结果和残差定义，避免把CL/CD趋稳误判为求解器已验证正确。
