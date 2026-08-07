# 14.04 度磁铁斜贴方案

Feature Name: magnet-skew-14-04
Updated: 2026-08-06

## Description

在既有 10:2 无侧片有限磁铁阵列上，复用斜贴偏移模型生成 14.04 度方案的预计算资源。

## Components and Interfaces

- `cogging_fdm.py --angle 14.04` 生成 14.04 度边缘截面场图、往复动画和齿槽力曲线。
- `halbach_field_n52.html` 为新资源提供方案按钮和资源映射。

## Correctness Properties

- 斜贴偏移由 `16 mm × tan(14.04°)` 计算。
- 动画使用 200 帧 `-24 mm → +24 mm → -24 mm` 往复轨迹，曲线使用 100 个单程位置点。
- 页面按钮名称、静态场图、动画和力曲线指向同一个 14.04 度方案。

## Test Strategy

- 执行 Python 编译检查。
- 检查三个 14.04 度输出文件存在。
- 在预览页面检查按钮资源映射。
