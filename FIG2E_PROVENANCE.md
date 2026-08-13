# Figure 2E provenance

FIG2E_STATUS=SOURCE_DATA_NOT_TRACEABLE

The frozen Phase 4B implementation created the 106-row continuous-time model input only as an in-memory object named `ct`. It wrote the frozen slope, confidence interval, P value, and days-metadata completeness table, but did not write the sample-level CORE100-plus-days input table. Reconstructing that table would require recomputing CORE100 from the raw count matrix, which Phase 10B prohibits.

ORIGINAL_FROZEN_IMPLEMENTATION=04_scripts/phase4b/run_gse178411_spectrum.R
FROZEN_SUMMARY=06_results/spectrum/GSE178411_TWPS_REPORT.md
PLOT_RECORD_N=NA
SAMPLE_CHANGE=NO
MODEL_CHANGE=NO
STATISTICAL_RECOMPUTATION=NO
ACTION=OMIT_FIGURE2E
