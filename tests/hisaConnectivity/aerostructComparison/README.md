# L3 静气弹：ADflow 与一阶 HiSA 对照（2026-11-05）

本次安装并实际运行了工作区新克隆的 ADflow。ADflow 的巡航、机动工况均达到现有脚本的静气弹耦合收敛条件；HiSA 首次刚性流场收敛，但第一次结构变形后的流场求解发散。尚未取得 HiSA 静气弹平衡解，也未计算此次静气弹伴随梯度。

## 安装

ADflow 2.13.1，源码提交 `8155e981`。原始克隆含 CRLF，因此复制到 `adflow/build-local/source`，仅在构建副本中规范化源码/构建脚本换行，使用 `config.LINUX_GFORTRAN.mk` 构建。原始源码和系统 ADflow 安装未覆盖。

依赖使用现有 GFortran/MPI、PETSc 3.21.6 (`real-opt`)、CGNS 4.5.0；具体路径和编译库 SHA256 见 `adflow_installation.json`。在构建副本执行 `make`，再执行：

```sh
python -m pip install --no-deps --no-build-isolation \
  --target /mdouser/test_aerostructure/adflow/build-local/python \
  /mdouser/test_aerostructure/adflow/build-local/source
source /mdouser/test_aerostructure/adflow/build-local/activate.sh
```

构建导入检查成功，完整计算实际导入了上述本地 Python 包。ADflow 通过 MExt 动态加载库，本次报告的 `libraries` 字段为空，不能将其解释为已记录运行中动态库哈希；安装产物哈希另记在安装记录中。后续运行的 harness 会另外记录安装目录的 `.so` 哈希。

## 对比条件与边界

- 使用现成的 `asopt2.py` 与 `asopt2-hisa.py` 模型定义；仅运行分析，不运行优化。
- L3 CGNS、Mach 0.85、海拔 10000 m；巡航迎角 2°、机动迎角 5°；参考面积 45.5，参考弦长 3.25。
- 同一 `wingbox.bdf`、`tacs_setup.py`、初始厚度 0.002，以及 MELD 耦合；巡航 1g、机动 2.5g。
- 两者串行、单 MPI rank。耦合均为带 Aitken 的 NLBGS，最大 100 次，绝对容差 1e-4、相对容差 1e-8。
- ADflow 保留原脚本 RANS/ANK 设置，通量离散为 `central plus scalar dissipation`；流场 `L2Convergence=1e-14`、相对阈值 1e-4、最大 10000 步。
- HiSA 使用已验证的一阶 AUSM 设置，GMRES 子空间 40，LUSGS 预处理，centralLF 近似雅可比；流场最大 4000 步，原生归一化外层阈值 1e-8。CFL 实际固定为 5；配置上限仍为 100000，增长因子为 1，因此本次不是大 CFL 试验。
- HiSA 使用转换/修复后的 Foam 体网格；ADflow 使用结构 CGNS。离散、体网格和表面数据布局并不完全相同，不能将本次结果视为同离散精度基准，也不能直接比较两个求解器的残差数值。

harness 修正了生成副本中的一个 MPhys 网格输出连接：本地 DAFoam 输出名为 `x_aero0`，而当前 MPhys 枚举为 `x_aero0_mesh`。原始用户脚本未改动。

## 实测结果

| 求解器/工况 | 耦合结果 | 外层迭代 | CL | CD |
|---|---|---:|---:|---:|
| ADflow 巡航 | 收敛 | 11 | 0.2199381093 | 0.02627823766 |
| ADflow 机动 | 收敛 | 11 | 0.4112184387 | 0.03701378749 |
| HiSA 巡航 | 第二次流场求解发散 | 未收敛 | — | — |
| HiSA 机动 | 前一工况退出，未执行 | — | — | — |

ADflow 两工况合计 `run_model` 约 40.03 s，初始化约 11.70 s。结构质量 380.54588；巡航/机动 KS failure 分别为 1.04347 / 1.84791，意味着这个初始设计不满足对应的 KS≤1 约束，并不表示计算未收敛。完整输出和耦合日志分别见 `adflow_report.json`、`adflow_coupling_excerpt.log`。

HiSA 的首次刚性流场在 2171 步达到原生外层阈值。该状态 CL=0.29047468、CD=0.16378033，**不是静气弹平衡结果，不能与表中的 ADflow 平衡值直接作精度比较**。DAFoam 随后打印的未归一化总残差范数为 0.406094，与原生 1e-8 阈值不是同一指标。

第一次传入结构变形后，网格最小体积仍为正（3.4542e-9），但出现 1 个面方向异常，最大非正交角 87.81°。检查配置允许 1 个异常面，所以仍打印 `Mesh OK`。随后第二次流场求解出现负密度、温度截断，最终在热力学温度反算达到 100 次迭代上限时原生退出。证据见 `hisa_failure_excerpt.log`。

这支持“当前 HiSA 耦合流程的鲁棒性不足”的判断，但还不能证明错误来自 HiSA 核心、网格变形、载荷传递或变形后状态恢复中的哪一环。未得到 HiSA 完整计算时间，因此不报告求解器速度比。下一步应固定首次结构位移，分别检查变形网格与流场重启一致性，再定位发散来源。

额外尝试过 `maxiter=1` 的单次耦合诊断，但实际执行了 0 次耦合迭代，并产生未初始化/非有限输出。这一试验无效，已从可用 harness 中移除；`reports/aerostruct_comparison/hisa_one_pass/report.json` 的 `completed=true` 仅说明函数返回，不代表有效计算，禁止用于比较。

## 复现

在工作区根目录，先加载现有 OpenFOAM/DAFoam 环境：

```sh
export DAFOAM_ROOT_PATH=/home/dafoamuser/dafoam
source /home/dafoamuser/dafoam/OpenFOAM/OpenFOAM-v2506/etc/bashrc
source dafoam/build/hisa-local/activate.sh
source adflow/build-local/activate.sh
export LD_LIBRARY_PATH=/mdouser/test_aerostructure/dafoam/build/hisa-diagnostics/lib:/mdouser/test_aerostructure/dafoam_case/validation/hisa_direction_build/lib:$LD_LIBRARY_PATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENMDAO_REPORTS=0
ulimit -c 0
python dafoam/tests/hisaConnectivity/aerostruct_compare.py \
  --root /mdouser/test_aerostructure --solver adflow \
  --output /mdouser/test_aerostructure/reports/aerostruct_comparison/adflow_new
# 切换 --solver hisa，并使用另一个新的 --output 目录。
```

输出目录必须不存在。运行目录内保存生成的模型脚本、报告和耦合场数组；原生 Foam abort 时 Python 无法写最终报告，必须同时保留标准输出/退出码。当前完整运行目录为 `reports/aerostruct_comparison/adflow`、`hisa_v2`；原始日志为 `/tmp/aerostruct_adflow.log`、`/tmp/aerostruct_hisa_v2.log`。这里归档了关键日志以避免临时文件清理后丢失证据。
