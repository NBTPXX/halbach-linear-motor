# Requirements Document

## Introduction

在每种磁铁方案的齿槽力曲线中显示主要 FFT 空间谐波分量。

## Requirements

### Requirement 1

**User Story:** AS 电机设计人员, I want 查看齿槽力的主要频率分量, so that 我可以判断波形由哪些空间谐波构成。

#### Acceptance Criteria

1. THE curve generator SHALL calculate the discrete Fourier transform from the 100 single-pass force samples.
2. THE curve renderer SHALL select the three largest non-DC Fourier components by amplitude.
3. THE curve renderer SHALL display every selected component as a dashed line.
4. THE curve renderer SHALL label every selected component with spatial period in millimetres and amplitude in newtons.
