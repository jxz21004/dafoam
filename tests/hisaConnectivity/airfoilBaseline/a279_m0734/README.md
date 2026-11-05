# RAE2822：2.79°、255.56 K、Mach 0.734

使用当前HiSA库，在独立副本从均匀来流启动。保持原翼型网格23552单元、SA模型、AUSMPlusUp/wVanLeer、原生特征边界及1250步预算。未重算ADflow。

按算例Cp=1005、molWeight=28.966计算R=287.0421397、gamma=1.399803603，对应速度235.2062634 m/s，U=(234.9274612,0,11.44876942)。用户未指定压力，沿用原算例108987.773 Pa，rho=1.485727768 kg/m³；弦长0.3048 m，对应Re约6.541e6。力系数参考速度、密度及升阻力方向同时更新，避免沿用旧迎角投影。参考面积仍为0.3048。

1250步正常退出，CPU约154.7秒，墙钟164秒：

- CL=0.777985707。
- CD=0.019120276。
- 原生外层归一化残差最大分量=0.00539803，**未达到1e-8**。
- 最后100步CL范围0.777968—0.778134，CD范围0.0191203—0.0192979。

这些是指定计算预算末尾的结果，不能称为充分收敛解，也尚未与实验数据验证。原网格长宽比/非正交性限制与上一级报告一致。

完整算例：`reports/hisa_airfoil_baseline/rae2822_a279_m0734`；日志`/tmp/rae279.log`。本目录保存来流及力系数配置、结果JSON、系数历史和收敛图。

复现汇总：

```sh
MPLCONFIGDIR=/tmp/mpl-audit python dafoam/tests/hisaConnectivity/airfoilBaseline/summarize.py \
 --case reports/hisa_airfoil_baseline/rae2822_a279_m0734 \
 --log /tmp/rae279.log \
 --output dafoam/tests/hisaConnectivity/airfoilBaseline/a279_m0734
```
