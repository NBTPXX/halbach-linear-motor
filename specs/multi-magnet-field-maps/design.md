# Multi-Magnet Field Maps

Feature Name: multi-magnet-field-maps
Updated: 2026-08-06

## Description

扩展单铁芯开放边界二维磁静态模型。磁铁采用五组连续单元并围绕 DXF 铁芯几何中心对称布置。铁芯使用原始 DXF 轮廓，绕自身中心旋转 180° 后放置于磁铁上方，齿顶与磁铁顶面间距为 1 mm。

## Components and Interfaces

- `cogging_fdm.py`：新增有限磁铁阵列参数、无侧片方案和 14.04 度、26.565 度斜贴几何，导出 PNG 场图。
- `halbach_field_n52.html`：为全部预计算场图提供选择按钮和场景说明。
- `cogging_fdm.py --animate`：以 200 帧计算有限阵列的移动位置能量导数并导出带齿槽力箭头的 GIF 与力-位置曲线；曲线使用 100 个单程位置点。

## Correctness Properties

- 每个场图只有一个由 `CORE_DXF` 表示的铁芯区域。
- 每个场图的磁铁源只落在有限的数组区间内。
- 图上的磁铁轮廓与求解器中的磁铁区域一一对应。
- GIF 中箭头方向与共能量导数的齿槽力符号一致。

## Test Strategy

- 执行 Python 编译检查。
- 生成所有场图并检查输出文件。
- 通过预览页面检查按钮映射。
- 使用图像分析核验单铁芯和多组磁铁轮廓。
