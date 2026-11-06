# RAE2822 高阶收敛试验

使用当前 HiSA 动态库，迎角 2.79°、Mach 0.734、来流温度 255.56 K，其他物理设置、网格、边界条件继承 `../a279_m0734`。本次仅调整算例字典，没有修改或重新编译求解器。未重算 ADflow。

## 本轮结果

高阶配置的原生外层残差由原先约 5.4e-3 改善到 **2.37358e-8**，但所有流场分量低于 1e-8 的严格判据尚未满足，日志没有报告收敛。SA Initial residual 为 **2.51096e-11**，CL=**0.7748341546**，CD=**0.01893123526**。相较原 1250 步高阶结果，CL 变化约 -0.405%，CD 约 -0.989%。这不是实验精度验证。

有效续算路径：原保存场 1200 → CFL=5 运行 1000 步 → CFL=50 运行 1500 步 → CFL=500 运行 1500 步 → CFL=5000 运行 4000 步，总共新增 8000 步。没有把对照分支的步数算入该路径。

最后从 CFL=5000 第 4000 步另做 200 步短对照（gauss_tight），仅将 GMRES 内层相对容差从 0.1 收紧至 1e-3，仍保留子空间 8。外层末值 2.39407e-8，SA 1.85175e-11，CL=0.7748341571，CD=0.01893123529。内层达到约 1e-11，并没有消除外层约 2e-8 平台；该分支 107.18 CPU 秒，原 CFL=5000 阶段 4000 步共 565.61 CPU 秒。这些阶段不等长，时间仅用于记录，不作严格性能基准。

建议保留 Gauss 梯度和原内层容差。若要求进一步收敛或伴随精度，应定位最后平台的方程/单元来源，并用独立残差与梯度检查验证；当前证据不支持继续盲目收紧 GMRES，也不足以断言平台就是双精度舍入误差。全部本轮求解已结束，无后台续算。

## 有效改动

保持 AUSMPlusUp 和 wVanLeer 重构，将 `system/fvSchemes` 中的默认梯度改为：

```foam
gradSchemes
{
    default         Gauss linear;
    grad(nuTilda)   cellLimited Gauss linear 0.5;
}
```

原设置为 `faceLeastSquares linear`。这仍使用高阶重构，没有退回一阶；实际局部精度受限制器和网格影响。默认梯度也可能影响黏性项，不能将试验直接解释为已证实的重构实现 bug。

## 对照和续算

前三项及 `gauss_gradient` 都从原算例保存的第 1200 步独立启动。后续三个 CFL 阶段串行续算，gauss_tight 再从最后阶段续算，各自迭代编号从零开始。图中的横轴是各阶段步数，不是共同起点的总步数。

| 算例 | 改动 | 步数 | 最后外层归一化残差最大值 |
|---|---|---:|---:|
| cfl5 | 原梯度、固定 CFL=5 | 1000 | 1.26348e-3 |
| minmod_cfl5 | CFL=5，Minmod 重构 | 1000 | 1.64739e-3 |
| tight_linear | CFL=5，GMRES 40，内层相对容差 1e-3 | 1000 | 1.48057e-3 |
| gauss_gradient | CFL=5，Gauss 默认梯度 | 1000 | 1.89035e-5 |
| gauss_cfl50 | 接上一行，CFL=50 | 1500 | 4.86684e-6 |
| gauss_cfl500 | 接上一行，CFL=500 | 1500 | 6.53228e-7 |
| gauss_cfl5000 | 接上一行，CFL=5000 | 4000 | 2.37358e-8 |
| gauss_tight | 接上一行，内层容差 1e-3，子空间仍为 8 | 200 | 2.39407e-8 |

除 tight_linear 和 gauss_tight 外，使用 GMRES 子空间 8、maxIter=20、内层相对容差 0.1，LUSGS 预处理，LaxFriedrichs 对流近似雅可比和 laplacian 黏性雅可比。没有通过扩大 GMRES 子空间获得本次改善。

CFL 在各阶段内固定，`resetPseudo true`，增长/下降因子都设为 1。这是人工分阶段续算，不是新实现的自动 CFL 控制算法。

## 残差证据与限制

原 faceLeastSquares 的 tight_linear 对照在第 500 步，未加权密度残差平方和约 97.4% 集中在十个单元，位于前缘近壁 x/c 约 0.0026–0.0043；不在此前 checkMesh 标记的最大长宽比单元组中。小体积单元会放大强形式残差，因此这只是定位证据，不足以单独证明离散错误。

同为 CFL=5，tight_linear 第 1000 步的原始 rhoR L2 约 5.24327e6，Gauss 梯度第 1000 步约 1.15394e4。改善并非只来自归一化分母改变。后续密度残差继续下降，但动量/能量的原始未加权 L2 下降不完全同步，见 results.json。

HiSA 日志的外层残差是每次 GMRES iteration 0 的原生归一化残差，不是相对首次残差，也不是 GMRES 内层最终残差。归一化与伪时间矩阵有关，不能将不同 CFL 阶段直接作为固定归一化的收敛率比较。SA 报告使用 nuTilda 的 Initial residual；Final residual 仅表示本次线性求解精度。

writeObjects 输出 rhoR/rhoUR/rhoER 是当步更新前组装的残差，不是对最终输出场另行重算的残差。即使达到原生 1e-8，也不等价于已经验证伴随梯度精度或实验一致性。

## 文件与复现

完整日志保存在每个算例的 `log.hisa`。完整算例及保存场位于工作区 `reports/hisa_airfoil_baseline/convergence/<name>/`。本目录保存配置、力系数历史、日志末尾、汇总和图，不纳入体网格和大型场文件。

`prepare.py` 创建独立续算副本，不修改输入算例；要求源时刻包含完整流场，例如：

```sh
python dafoam/tests/hisaConnectivity/airfoilBaseline/convergenceTrials/prepare.py \
 --source reports/hisa_airfoil_baseline/rae2822_a279_m0734 --time 1200 \
 --output reports/hisa_airfoil_baseline/convergence/new_trial \
 --cfl 5 --steps 1000 --gauss-gradient
```

载入 OpenFOAM-v2506 环境，并将 `dafoam_case/validation/hisa_direction_build/lib` 置于 LD_LIBRARY_PATH 首位，运行 `hisa -case <新算例>`。本次库 SHA256 为 `964a1065bf44fcd46305726113ff40e82c0a1e66d5036edc03437454e36e15be`。

汇总命令：

```sh
MPLCONFIGDIR=/tmp/mpl-audit python dafoam/tests/hisaConnectivity/airfoilBaseline/convergenceTrials/summarize.py
```

汇总优先读取各算例的 `log.hisa`，不存在时读取 `/tmp/rae_conv_*.log`。
