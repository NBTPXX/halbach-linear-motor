# Requirements Document

## Introduction

将 14.04 度斜贴作为独立的 10:2 无侧片磁铁方案加入 Halbach 直线电机分析页面。

## Glossary

- **14.04 度斜贴方案**：磁铁沿叠厚 16 mm 方向形成 14.04 度相对偏移的 10:2 无侧片磁铁阵列。

## Requirements

### Requirement 1

**User Story:** AS 电机设计人员, I want 查看 14.04 度斜贴磁铁方案, so that 我可以与零斜贴和 26.565 度斜贴方案比较齿槽力。

#### Acceptance Criteria

1. THE field-map generator SHALL render one 10:2 no-side-magnet field map with a 14.04 degree skew angle.
2. THE field-map generator SHALL export one 90-frame animation for the 14.04 degree skew scenario with the core moving in the sequence -24 mm, +24 mm, -24 mm.
3. THE field-map generator SHALL export one cogging-force curve for the 14.04 degree skew scenario.
4. WHEN a viewer selects the 14.04 degree skew scenario, THE web page SHALL display the matching field map, animation, force curve, and scenario label.
