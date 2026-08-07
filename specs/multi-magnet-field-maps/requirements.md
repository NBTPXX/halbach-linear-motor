# Requirements Document

## Introduction

磁场图应保留用户提供的单个 DXF 铁芯，并在铁芯上方展示多组有限长度的磁铁阵列，以便观察端部磁场和行程范围。

## Glossary

- **铁芯**：由用户提供的 DXF 闭合轮廓定义的单个 2D 软磁区域。
- **磁铁组**：沿 Y 方向连续放置的一组 24 mm Halbach 单元。
- **场图方案**：由磁化配比、侧片状态和斜贴角度共同定义的求解结果。

## Requirements

### Requirement 1

**User Story:** AS 电机设计人员, I want 单铁芯与多组磁铁同图显示, so that 我可以评估有限阵列端部场。

#### Acceptance Criteria

1. THE field-map generator SHALL include one DXF iron-core outline in every field-map solution.
2. THE field-map generator SHALL include five finite magnet groups above the iron core in every field-map solution.
3. THE field-map generator SHALL align the center of the five-group magnet array with the geometric center of the DXF iron core.
4. THE field-map renderer SHALL draw the outline of every magnet block included in the solution.
5. THE field-map generator SHALL rotate the DXF iron core by 180 degrees and place the tooth tips 1 mm above the magnet array.
6. THE field-map generator SHALL export a 200-frame GIF for each field-map scenario with the iron core moving in the sequence -24 mm, +24 mm, -24 mm relative to the magnet array.
7. THE GIF renderer SHALL show a force-direction arrow on the iron core, scale the arrow length by the instantaneous cogging-force magnitude, and display the core position and force value.
8. THE field-map generator SHALL export a cogging-force versus core-position curve for each GIF scenario.
9. THE curve generator SHALL calculate the discrete Fourier transform of each force curve and render the three largest non-DC components as dashed lines with spatial-period and amplitude labels.

### Requirement 2

**User Story:** AS 电机设计人员, I want 查看关键磁铁结构, so that 我可以比较配比、侧片和斜贴的影响。

#### Acceptance Criteria

1. THE field-map generator SHALL render 6:6, 8:4, 10:2 and 11:1 Halbach configurations with side magnets.
2. THE field-map generator SHALL render a 10:2 configuration without side magnets.
3. THE field-map generator SHALL render 10:2 configurations with 14.04 degree and 26.6 degree skew angles.
4. WHEN a viewer selects a field-map scenario, THE web page SHALL display the matching generated image, animation, force curve, and scenario label.
