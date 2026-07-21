# Final Analysis Lock Audit Report

This document serves as the formal audit record of the final analysis pipeline executed on 2026-07-21 for the study **"Cross-National Heterogeneity in the Protective Effect of Education on Cognitive Decline."**

## 1. Primary Analysis Conclusions

Based on the full analytic sample of 144,642 unique individuals (542,426 longitudinal observations), the final fixed models conclude:
1. 全样本中，中国Mid–Low教育组的第4年瞬时年斜率差约0.22；
2. 该教育梯度显著大于三个比较队列（美国 HRS、英国 ELSA、欧洲 SHARE）；
3. 童年 SES 调整后估计基本不变；
4. 经济共同样本中，财富与收入联合调整（M4）使点估计约衰减10%，但该结果仅作为机制探索，不是正式的中介效应结论；
5. CHARLS 2011 初始队列独立验证支持了主要的教育梯度差异；
6. IPCW 敏感性分析显示主要结论对当前的缺失权重方案不敏感，但这不能证明完全不存在失访或死亡偏倚。

## 2. Imputation & Data Constraints

- **Europe_Unknown Country Codes**: There are 9,719 individuals from Europe with an "Unknown" specific country code. Among them, 358 individuals required imputation for childhood SES. Due to this minor imputation constraint, the M1 model incorporates these within the pooled SHARE sample. No further re-estimation of M1 is required.
- **Algebraic Consistency**: All cross-national algebraic consistency checks (e.g., `(China - USA) == (China_Slope - USA_Slope)`) passed with a residual tolerance of `1e-8`.

## 3. Sensitivity Analyses

### A. CHARLS Cohort Validation
- **Wave 2011 Cohort**: Supported the main educational gradient.
- **Wave 2013 Refreshment**: 无直接 t=4 观测支持，不作为正式敏感性结论。
- **Combined 2011+2013 Cohort**: Validated the primary gradient findings.

### B. Collinearity & Model Instability (M4)
- **M4 Diagnostics**: 设计矩阵满秩（88/88），但高阶交互模型条件指数较高（Kappa ~ 20293）。
- **Conclusion**: 关键 M4 边际对比在 20 套插补中保持稳定，因此 M4 仅作为机制探索敏感性模型，单项交互系数不作实质解释。

## 4. Output Artifacts

The final run produced the following cryptographic manifest for traceability:
- `final_primary_results.csv`
- `final_mechanism_results.csv`
- `final_sensitivity_results.csv`
- `charls_sensitivity.csv`
- `ipcw_audit.csv`
- `final_manifest_sha256.txt`

The pipeline is now locked. No further models, subgroups, or bootstrap iterations will be added.
