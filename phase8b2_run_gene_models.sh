#!/bin/zsh
set -euo pipefail

project_dir=${1:?project directory required}
cd "$project_dir"

export LC_ALL=C
magma='data/external/genetics/MAGMA/bin/magma_v1.10_mac/magma'
gene_loc='data/external/genetics/MAGMA/gene_location/NCBI37.3/NCBI37.3.gene.loc'
work='05_logs/genetics/phase8b2_magma_inputs'
results='06_results/genetics/magma'
mkdir -p "$work" "$results"

for ancestry in EUR AFR; do
  if [[ "$ancestry" == EUR ]]; then
    accession='GCST90652488'
    reference='data/external/genetics/MAGMA/reference/g1000_eur/g1000_eur'
    sample_n='1282582'
  else
    accession='GCST90652489'
    reference='data/external/genetics/MAGMA/reference/g1000_afr/g1000_afr'
    sample_n='139538'
  fi

  python3 04_scripts/genetics/phase8b2_prepare_magma_inputs.py \
    --sumstats "data/external/genetics/GWAS/${accession}.h.tsv.gz" \
    --pval-output "$work/${ancestry}_all_valid_unique_rsid_p.tsv" \
    --qc-output "$work/${ancestry}_sumstat_preparation_qc.tsv"

  awk 'BEGIN { OFS="\t" } { print $2, $1, $4 }' "${reference}.bim" > "$work/${ancestry}_reference_GRCh37_snp_locations.tsv"

  "$magma" \
    --annotate window=0,0 \
    --gene-loc "$gene_loc" \
    --snp-loc "$work/${ancestry}_reference_GRCh37_snp_locations.tsv" \
    --out "$work/${ancestry}_0kb_all_reference_annotation"

  "$magma" \
    --bfile "$reference" \
    --pval "$work/${ancestry}_all_valid_unique_rsid_p.tsv" "N=${sample_n}" \
    --gene-annot "$work/${ancestry}_0kb_all_reference_annotation.genes.annot" \
    --out "$results/${ancestry}_gene_results"
done
