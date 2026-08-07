# FFT 齿槽力分量显示

Feature Name: fft-cogging-components
Updated: 2026-08-07

## Description

对每条 100 点齿槽力曲线执行实数 FFT，保留幅值最大的三个非直流频谱分量，并与总齿槽力曲线叠加显示。

## Components and Interfaces

- `cogging_fdm.py --curve-only --curve <path>`：计算力样本、生成 FFT 分量虚线并输出曲线 PNG。

## Correctness Properties

- FFT 直流分量不参与主要分量排序。
- 分量振幅采用单边谱幅值 `2 × |FFT| / N`。
- 图例显示每个分量的空间周期与振幅。

## Test Strategy

- 执行 Python 编译检查。
- 生成七条曲线并检查标题、虚线和图例。
- 使用图片分析核验曲线可读性。
