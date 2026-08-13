# TWPS SCORING SPECIFICATION — LOCKED BEFORE VALIDATION

## Locked endpoint

- Primary program: `D7_M3_CORE100`
- Higher score means greater activity of the transient wound program.
- All 100 locked genes have positive program direction. Gene direction must not be reversed in response to validation results.

## Bulk-expression scoring

For each dataset independently:

1. Map gene identifiers to approved gene symbols and document the mapping.
2. Intersect the expression matrix with the locked CORE100; do not impute or replace missing genes.
3. For every available locked gene, calculate its z-score across samples within that dataset using a transformation appropriate to the processed platform scale.
4. For every sample, calculate TWPS as the arithmetic mean of the available locked-gene z-scores.

For available locked gene `g` and scored analysis unit `i`, the notation is `z_gi = (x_gi - xbar_g) / s_g`, where `x_gi` is expression of gene `g` in unit `i`, `xbar_g` is the mean expression of gene `g` across all scored units within that dataset, and `s_g` is the standard deviation of gene `g` across all scored units within that dataset. Because this standardization is performed separately within each dataset, absolute TWPS magnitudes must not be compared across datasets or platforms; only within-dataset contrasts, directions and evidence patterns are interpretable.

Gene-wise standardization must not use disease labels, outcomes, or group assignments to select genes or directions. Dataset-specific normalization is permitted only to accommodate platform scale; it may not alter membership or direction using outcome information. A mathematically equivalent platform-appropriate implementation is allowed only when specified prospectively in the relevant future phase protocol.

## Pseudobulk or scRNA validation

Use the same locked positive gene direction with a compatible module-score implementation. Preserve donor or patient as the biological inference unit; cells are not independent biological replicates.

## Coverage and reporting

- Primary analysis requires at least 80 of 100 locked genes.
- If fewer than 80 are available, set `primary_score_status=LOW_COVERAGE` and obtain protocol review before inference.
- For every validation dataset report `N_locked`, `N_available`, and `coverage_fraction`.
- Missing genes must not be replaced by correlated substitutes.

## Preregistered Gate B primary hypothesis

- Dataset: GSE113619
- Primary endpoint: locked D7_M3_CORE100 TWPS
- Primary comparison: within-person change from baseline to post-wounding
- Primary inferential contrast: `(keloid-prone post - keloid-prone baseline) - (healthy post - healthy baseline)`
- Equivalent implementation: group × time interaction with subject pairing preserved
- Expected persistence-model direction: the post-wounding TWPS change in keloid-prone individuals shows greater persistence, less attenuation, or greater residual activity than in healthy controls.

Project continuation has no required P-value threshold. Interpretation must jointly consider effect size, 95% confidence interval, direction, paired consistency, and statistical evidence. No other module may replace the primary endpoint after Gate B results are viewed.

## Locked Gate B sensitivity and secondary analyses

The following must all be reported as sensitivity checks, not alternative primary endpoints:

- D7_M3_FULL
- D7_M3_CORE50
- D7_M3_CORE25

The secondary discovery modules D7_M1, D7_M2, D1_M2, and D1_M1 may later be tested only as clearly labelled secondary exploratory endpoints with appropriate multiple-testing handling. They may not replace the primary TWPS because they validate better.
