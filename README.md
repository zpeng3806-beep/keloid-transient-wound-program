# keloid_TWPS_project analysis code and reproducibility package

This candidate contains the final analysis scripts supporting the frozen manuscript evidence. It excludes raw or controlled data, personal filesystem paths, credentials, private notes, superseded Gate B code and publication-only rendering code.

## Runtime and reuse information

- Project name: `keloid_TWPS_project`.
- Operating system: the analyses and final package were run or assembled on macOS on Apple Silicon (`arm64`); the exact macOS version used for every analytical phase is `NOT_TRACEABLE`. The final package was assembled on macOS 26.5.1 (build 25F80).
- Programming languages: R, Python and zsh shell scripts.
- Recorded versions: R 4.6.1; Python 3.13.15 for the final GSE241124 run; MAGMA v1.10 Mac binary. Exact cross-phase Python and several package versions are `NOT_TRACEABLE`; see `SOFTWARE_REPRODUCIBILITY.md`.
- Recorded Python dependencies: numpy 2.5.1 and pandas 3.0.5 for GSE241124. Other required modules and R packages can be identified from script import/library statements, but versions not listed in `SOFTWARE_REPRODUCIBILITY.md` are `NOT_TRACEABLE`.
- Input accessions: GSE241132, GSE113619, GSE181316, GSE178411, GSE241124, GCST90652488 and GCST90652489. Public locations and expected inputs are listed in `DATASET_ACCESSION_MANIFEST.tsv`.
- Usage: obtain the listed source inputs, preserve the directory layout and run the scripts in the order below. Scripts contain input, dimension, hash, coverage and sample-structure checks and will stop on mismatches.
- License: original code and associated repository materials for which the authors have redistribution rights are released under the MIT License; see `LICENSE`.
- Third-party data: original public GEO and GWAS datasets are not redistributed here and remain subject to the terms of their respective source repositories. Retrieve them using the accession numbers listed in `DATASET_ACCESSION_MANIFEST.tsv`.
- Public home page or archived identifier: not assigned.

## Evidence and signature locks

- Final evidence freeze SHA256: `94d1f6fc991eca017f0c53e63bb3710a1d2b72eb5eb187ae7715e8a8bba4ac7d`
- CORE100 SHA256: `6d5a62c080f2d0d9e4077983294d2a6a443a9e0e3c565691a4d3a863cc2afcc3`
- The signature manifest and a publication copy of the frozen scoring specification are in `locks/`. The public copy changes only an internal review-routing label; the scoring rules are unchanged.
- The public scoring specification also states the final unambiguous gene-standardization notation and the prohibition on comparing absolute TWPS magnitudes across independently standardized datasets; this clarification does not alter the scoring implementation.
- Gate B is represented only by `run_gate_b_technical_correction.R` (final corrected 96/100 implementation).
- The final Gate B copy retains the final ARNTL2/BMAL2 mapping and all frozen inferential procedures but omits the internal comparison table against superseded outputs.

## Contents and execution order

Run from a project directory that preserves the paths shown below after obtaining the public inputs listed in `DATASET_ACCESSION_MANIFEST.tsv`.

1. Normal-wound discovery and lock:
   - `Rscript 04_scripts/gateA1/run_gate_a1.R`
   - `Rscript 04_scripts/gateA1/run_gate_a2.R`
   - `python3 04_scripts/gateA1/run_gate_a2_enrichment.py`
   - `python3 04_scripts/gateA1/lock_gate_a3.py`
2. Corrected longitudinal Gate B:
   - `Rscript 04_scripts/gateB/run_gate_b_technical_correction.R`
3. Established-keloid Gate C:
   - `python3 04_scripts/gateC/run_gse181316_gate_c.py`
4. GSE178411 structure and patient-aware spectrum:
   - `python3 04_scripts/phase4a/build_gse178411_structure_audit.py`
   - `Rscript 04_scripts/phase4b/run_gse178411_spectrum.R`
5. GSE241124 spatial confirmation:
   - `python3 04_scripts/phase7b_gse241124_spatial_twps.py`
6. Genetics:
   - build locked MAGMA sets with `phase8b2_build_locked_sets.py` using its required command-line paths;
   - run `zsh 04_scripts/genetics/phase8b2_run_gene_models.sh <project-directory>`;
   - freeze the primary results with `phase8b2_freeze_primary.py`;
   - create secondary evidence tables with `phase8b2_secondary_evidence.py --project <project-directory>`.

The scripts intentionally stop when expected inputs, dimensions, hashes, coverage or sample structures differ from the frozen analysis. This package was assembled without rerunning the analyses.

`MANUSCRIPT_NUMBER_PROVENANCE.tsv` maps the principal manuscript values to their frozen project sources.

One non-scientific multipart boundary and HTTP user-agent label in the enrichment helper were neutralized for public release. The Gate B public copy was also limited to final corrected outputs. Request content, final calculations and analysis logic were not changed.

## Known reproducibility limits

- The upstream GSE113619 technical-replicate merge algorithm is not traceable; the final official normalized biological-sample columns were analyzed as deposited.
- The version/procedure used upstream to harmonize GWAS coordinates to GRCh37 is not traceable. The final MAGMA implementation anchors retained rsIDs to ancestry-matched GRCh37 reference BIM positions.
- Exact cross-phase Python and several package versions were not recorded.
- The exact 106-row GSE178411 continuous-time model input was reconstructed from the frozen raw-to-score implementation and is provided in `reconstructed_inputs/`. Its SHA256 is `8172cfd169f99bc6c6eefea1399add675b80a9aeebc25ec7c18a90feddd63885`; the reconstructed slope, confidence interval and P value exactly match the frozen result.
- Public deposition location and DOI will be added only after verified publication.

Release classification: `READY_FOR_PUBLIC_DEPOSIT_WITH_DECLARED_LIMITATIONS`.

## Main-figure regeneration

From the project root, run `python3 04_scripts/publication_build_scirep_figures.py --main-only` after placing the frozen result tables at their documented paths. This performs presentation-only rendering and does not rerun biological analyses.
