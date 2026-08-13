#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260811)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grep("^--file=", args)])
root <- normalizePath(file.path(dirname(script_arg), "..", ".."))
res <- file.path(root, "06_results", "gateB", "corrected")
dir.create(res, recursive = TRUE, showWarnings = FALSE)

locked_path <- function(name) file.path(root, "06_results", "gateA", name)
read_locked <- function(name) read.delim(locked_path(name), check.names = FALSE)$gene
fmt <- function(x, digits = 12) format(as.numeric(x), digits = digits, trim = TRUE)

# Lock and frozen-input checks.
lock_hash <- strsplit(system2("shasum", c("-a", "256", locked_path("TWPS_PRIMARY_D7_M3_CORE100.tsv")), stdout = TRUE), " +")[[1]][1]
stopifnot(lock_hash == "6d5a62c080f2d0d9e4077983294d2a6a443a9e0e3c565691a4d3a863cc2afcc3")

normalized_path <- file.path(root, "data", "raw_small", "GSE113619", "GSE113619_RNA-seq_keloids_normalized.csv.gz")
gene_info_path <- file.path(root, "data", "intermediate", "gateB", "Homo_sapiens.gene_info.gz")
norm <- read.csv(gzfile(normalized_path), row.names = 1, check.names = FALSE)
stopifnot(nrow(norm) == 24943L, ncol(norm) == 26L, !anyNA(norm), all(is.finite(as.matrix(norm))), !anyDuplicated(rownames(norm)))

sample_map_path <- file.path(root, "06_results", "gateB", "GSE113619_FINAL_SAMPLE_MAP.tsv")
sample_hash <- strsplit(system2("shasum", c("-a", "256", sample_map_path), stdout = TRUE), " +")[[1]][1]
stopifnot(sample_hash == "cc75dddda0a75df0aa1608b320e274f8b06e9195892869c3411112d842443347")
sample_map <- read.delim(sample_map_path, check.names = FALSE)
paired_map <- sample_map[sample_map$included == "YES", ]
stopifnot(nrow(paired_map) == 24L, length(unique(paired_map$subject_id)) == 12L)
stopifnot(all(table(paired_map$subject_id) == 2L))
stopifnot(sum(!duplicated(paired_map$subject_id) & paired_map$group == "keloid-prone") == 8L)
stopifnot(sum(!duplicated(paired_map$subject_id) & paired_map$group == "healthy") == 4L)

# Apply the frozen approved-symbol collapse.
gi <- read.delim(gzfile(gene_info_path), quote = "", comment.char = "", check.names = FALSE)
approved <- ifelse(gi$Symbol_from_nomenclature_authority == "-", gi$Symbol, gi$Symbol_from_nomenclature_authority)
entrez_to_symbol <- setNames(approved, as.character(gi$GeneID))
symbols <- unname(entrez_to_symbol[rownames(norm)])
mapped <- !is.na(symbols) & symbols != "-" & symbols != ""
collapse_symbol <- function(mat, sym) {
  x <- as.matrix(mat[mapped, , drop = FALSE])
  sym <- sym[mapped]
  sums <- rowsum(x, group = sym, reorder = FALSE)
  counts <- as.numeric(table(factor(sym, levels = rownames(sums))))
  sums / counts
}
norm_symbol <- collapse_symbol(norm, symbols)

# Sole authorized correction: locked ARNTL2 deterministically resolves to
# Entrez 56938, whose approved symbol is BMAL2 in the frozen annotation snapshot.
alias_rows <- gi[gi$GeneID == 56938L, ]
stopifnot(nrow(alias_rows) == 1L, approved[gi$GeneID == 56938L] == "BMAL2")
stopifnot("ARNTL2" %in% strsplit(alias_rows$Synonyms, "\\|")[[1]])
stopifnot("56938" %in% rownames(norm), !"ARNTL2" %in% rownames(norm_symbol))
corrected_symbol <- rbind(norm_symbol, ARNTL2 = as.numeric(norm["56938", ]))
colnames(corrected_symbol) <- colnames(norm)

primary <- read_locked("TWPS_PRIMARY_D7_M3_CORE100.tsv")
core50 <- read_locked("TWPS_SENSITIVITY_D7_M3_CORE50.tsv")
core25 <- read_locked("TWPS_SENSITIVITY_D7_M3_CORE25.tsv")
full <- read_locked("TWPS_SENSITIVITY_D7_M3_FULL.tsv")
signatures <- list(CORE100 = primary, CORE50 = core50, CORE25 = core25, FULL = full)

coverage <- lapply(signatures, function(g) {
  available <- intersect(g, rownames(corrected_symbol))
  list(n = length(available), fraction = length(available) / length(g), available = available)
})
stopifnot(coverage$CORE100$n == 96L, coverage$CORE50$n == 46L, coverage$CORE25$n == 23L, coverage$FULL$n == 464L)

# Every locked CORE100 member is recorded before outcome calculation.
coverage_rows <- lapply(primary, function(gene) {
  direct <- gene %in% rownames(norm_symbol)
  corrected_alias <- identical(gene, "ARNTL2")
  available <- gene %in% rownames(corrected_symbol)
  gene_id <- if (corrected_alias) "56938" else {
    ids <- names(entrez_to_symbol)[entrez_to_symbol == gene & names(entrez_to_symbol) %in% rownames(norm)]
    if (length(ids)) paste(ids, collapse = ";") else ""
  }
  data.frame(
    locked_gene = gene,
    mapped_entrez_id = gene_id,
    mapping_route = if (direct) "approved_symbol" else if (corrected_alias) "verified_alias_ARNTL2_to_BMAL2" else "unavailable",
    duplicate_mapping = if (length(strsplit(gene_id, ";", fixed = TRUE)[[1]]) > 1L && gene_id != "") "YES" else "NO",
    final_available = if (available) "YES" else "NO",
    reason_missing = if (available) "" else "No deterministic approved-symbol or authorized ARNTL2 alias row in the frozen primary matrix"
  )
})
coverage_audit <- do.call(rbind, coverage_rows)
write.table(coverage_audit, file.path(res, "GATE_B_CORE100_COVERAGE_CORRECTED.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

zscore_rows <- function(mat) {
  means <- rowMeans(mat)
  sds <- apply(mat, 1, sd)
  keep <- is.finite(sds) & sds > 0
  z <- sweep(mat[keep, , drop = FALSE], 1, means[keep], "-")
  sweep(z, 1, sds[keep], "/")
}
corrected_z <- zscore_rows(corrected_symbol)
score_signature <- function(genes) {
  available <- intersect(genes, rownames(corrected_z))
  colMeans(corrected_z[available, , drop = FALSE])
}
signature_scores <- lapply(signatures, score_signature)

scores <- paired_map[, c("subject_id", "group", "timepoint")]
for (nm in names(signature_scores)) scores[[nm]] <- unname(signature_scores[[nm]][paired_map$biological_sample_id])
scores$CORE100_available <- coverage$CORE100$n
scores$coverage <- coverage$CORE100$fraction
write.table(scores, file.path(res, "GATE_B_TWPS_SCORES_CORRECTED.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

wide_changes <- function(score_col, prefix = "") {
  ids <- unique(scores$subject_id)
  out <- do.call(rbind, lapply(ids, function(id) {
    d <- scores[scores$subject_id == id, ]
    b <- d[[score_col]][d$timepoint == "baseline"]
    p <- d[[score_col]][d$timepoint == "post-wounding"]
    data.frame(subject_id = id, group = d$group[1], baseline = b, post = p, delta = p - b)
  }))
  names(out)[3:5] <- paste0(prefix, c("baseline_TWPS", "post_TWPS", "delta_TWPS"))
  out
}
changes <- wide_changes("CORE100")
for (nm in c("CORE50", "CORE25", "FULL")) {
  x <- wide_changes(nm, paste0(nm, "_"))
  changes <- merge(changes, x[, c("subject_id", paste0(nm, c("_baseline_TWPS", "_post_TWPS", "_delta_TWPS")))], by = "subject_id", sort = FALSE)
}
changes <- changes[match(unique(scores$subject_id), changes$subject_id), ]
write.table(changes, file.path(res, "GATE_B_SUBJECT_CHANGES_CORRECTED.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

hedges_g <- function(x, y) {
  nx <- length(x); ny <- length(y)
  sp <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  d <- (mean(x) - mean(y)) / sp
  J <- 1 - 3 / (4 * (nx + ny - 2) - 1)
  J * d
}
exact_perm <- function(x, y) {
  vals <- c(x, y); nx <- length(x); obs <- mean(x) - mean(y)
  comb <- combn(seq_along(vals), nx)
  stats <- apply(comb, 2, function(idx) mean(vals[idx]) - mean(vals[-idx]))
  c(two = mean(abs(stats) >= abs(obs) - 1e-14), one = mean(stats >= obs - 1e-14), permutations = length(stats))
}
bootstrap_ci <- function(x, y, B = 10000L) {
  diffs <- numeric(B); gs <- numeric(B)
  for (i in seq_len(B)) {
    xb <- sample(x, length(x), replace = TRUE); yb <- sample(y, length(y), replace = TRUE)
    diffs[i] <- mean(xb) - mean(yb)
    gs[i] <- suppressWarnings(hedges_g(xb, yb))
  }
  list(diff = unname(quantile(diffs, c(.025, .975), na.rm = TRUE)),
       g = unname(quantile(gs[is.finite(gs)], c(.025, .975), na.rm = TRUE)))
}
iqr_string <- function(x) paste(format(quantile(x, c(.25, .75)), digits = 8), collapse = ",")
analyse_delta <- function(sig, delta_col, genes) {
  kp <- changes[[delta_col]][changes$group == "keloid-prone"]
  hc <- changes[[delta_col]][changes$group == "healthy"]
  perm <- exact_perm(kp, hc)
  welch <- t.test(kp, hc, var.equal = FALSE)
  boot <- bootstrap_ci(kp, hc)
  data.frame(
    signature = sig, gene_n_available = coverage[[sig]]$n, coverage = coverage[[sig]]$fraction,
    kp_baseline_mean = mean(changes[[sub("delta", "baseline", delta_col, fixed = TRUE)]][changes$group == "keloid-prone"]),
    kp_post_mean = mean(changes[[sub("delta", "post", delta_col, fixed = TRUE)]][changes$group == "keloid-prone"]),
    kp_delta_mean = mean(kp), kp_delta_sd = sd(kp), kp_delta_median = median(kp), kp_delta_iqr = iqr_string(kp),
    healthy_baseline_mean = mean(changes[[sub("delta", "baseline", delta_col, fixed = TRUE)]][changes$group == "healthy"]),
    healthy_post_mean = mean(changes[[sub("delta", "post", delta_col, fixed = TRUE)]][changes$group == "healthy"]),
    healthy_delta_mean = mean(hc), healthy_delta_sd = sd(hc), healthy_delta_median = median(hc), healthy_delta_iqr = iqr_string(hc),
    difference = mean(kp) - mean(hc), hedges_g = hedges_g(kp, hc), ci_low = boot$g[1], ci_high = boot$g[2],
    bootstrap_difference_low = boot$diff[1], bootstrap_difference_high = boot$diff[2],
    permutation_p = perm["two"], permutation_p_one_sided = perm["one"], permutations = perm["permutations"],
    welch_difference_low = unname(welch$conf.int[1]), welch_difference_high = unname(welch$conf.int[2]), welch_p = welch$p.value
  )
}

analyses <- list(
  CORE100 = analyse_delta("CORE100", "delta_TWPS", primary),
  CORE50 = analyse_delta("CORE50", "CORE50_delta_TWPS", core50),
  CORE25 = analyse_delta("CORE25", "CORE25_delta_TWPS", core25),
  FULL = analyse_delta("FULL", "FULL_delta_TWPS", full)
)
sens <- do.call(rbind, analyses)
rownames(sens) <- NULL
write.table(sens, file.path(res, "GATE_B_SIGNATURE_SENSITIVITY_CORRECTED.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

primary_a <- analyses$CORE100
sens_dir <- sum(vapply(analyses[c("CORE50", "CORE25", "FULL")], function(x) x$difference > 0, logical(1)))
if (primary_a$difference < 0 && sens_dir <= 1) {
  classification <- "CONTRADICTORY"
} else if (primary_a$difference <= 0 && sens_dir < 2) {
  classification <- "NO_SUPPORT"
} else if (primary_a$difference > 0 && abs(primary_a$hedges_g) >= .8 && primary_a$permutation_p < .05 && sens_dir >= 2) {
  classification <- "STRONG_SUPPORT"
} else if (primary_a$difference > 0 && abs(primary_a$hedges_g) >= .5 && sens_dir >= 2) {
  classification <- "MODERATE_SUPPORT"
} else {
  classification <- "WEAK_SUPPORT"
}

report <- c(
  "# Final Gate B Analysis", "",
  "IMPLEMENTATION=", "final corrected ARNTL2/BMAL2 mapping", "",
  "TWPS_CHANGED=", "NO", "", "SUBJECTS_CHANGED=", "NO", "", "INPUT_DATASET_CHANGED=", "NO", "", "PRIMARY_ENDPOINT_CHANGED=", "NO", "",
  "PRIMARY_SIGNATURE=", "D7_M3_CORE100", "", "LOCK_SHA256_MATCH=", "YES", "", "PAIRED_TOTAL=", "12", "", "KELOID_PRONE_N=", "8", "", "HEALTHY_N=", "4", "",
  "CORE100_COVERAGE=", "96/100", "",
  "## Corrected primary result", "",
  paste0("KP_BASELINE_MEAN=", fmt(primary_a$kp_baseline_mean)), paste0("KP_POST_MEAN=", fmt(primary_a$kp_post_mean)), paste0("KP_DELTA=", fmt(primary_a$kp_delta_mean)), "",
  paste0("HEALTHY_BASELINE_MEAN=", fmt(primary_a$healthy_baseline_mean)), paste0("HEALTHY_POST_MEAN=", fmt(primary_a$healthy_post_mean)), paste0("HEALTHY_DELTA=", fmt(primary_a$healthy_delta_mean)), "",
  paste0("DIFFERENCE_IN_DELTA=", fmt(primary_a$difference)), paste0("HEDGES_G=", fmt(primary_a$hedges_g)), paste0("HEDGES_G_95CI=", fmt(primary_a$ci_low), ",", fmt(primary_a$ci_high)),
  paste0("BOOTSTRAP_DIFFERENCE_95CI=", fmt(primary_a$bootstrap_difference_low), ",", fmt(primary_a$bootstrap_difference_high)),
  paste0("PERMUTATION_P_TWO_SIDED=", fmt(primary_a$permutation_p)), paste0("PERMUTATION_P_ONE_SIDED=", fmt(primary_a$permutation_p_one_sided)), paste0("WELCH_P=", fmt(primary_a$welch_p)), "",
  "## Sensitivity", "",
  paste0("CORE25=difference ", fmt(analyses$CORE25$difference), "; Hedges g ", fmt(analyses$CORE25$hedges_g), "; permutation P(two-sided) ", fmt(analyses$CORE25$permutation_p)),
  paste0("CORE50=difference ", fmt(analyses$CORE50$difference), "; Hedges g ", fmt(analyses$CORE50$hedges_g), "; permutation P(two-sided) ", fmt(analyses$CORE50$permutation_p)),
  paste0("FULL=difference ", fmt(analyses$FULL$difference), "; Hedges g ", fmt(analyses$FULL$hedges_g), "; permutation P(two-sided) ", fmt(analyses$FULL$permutation_p)), "",
  "## Evidence classification", "", classification, "",
  "## Gate C propagation audit", "", "GATE_C_SHARED_SCORING_CODE=", "NO", "", "GATE_C_TECHNICAL_REAUDIT_REQUIRED=", "NO"
)
writeLines(report, file.path(res, "GATE_B_REPORT_CORRECTED.md"))

cat("TWPS_LOCK_MATCH=YES\n")
cat("CORRECTED_COVERAGE=96/100\n")
cat("CORRECTED_KP_DELTA=", fmt(primary_a$kp_delta_mean), "\n", sep = "")
cat("CORRECTED_HEALTHY_DELTA=", fmt(primary_a$healthy_delta_mean), "\n", sep = "")
cat("CORRECTED_DIFFERENCE=", fmt(primary_a$difference), "\n", sep = "")
cat("CORRECTED_HEDGES_G=", fmt(primary_a$hedges_g), "\n", sep = "")
cat("CORRECTED_PERMUTATION_P=", fmt(primary_a$permutation_p), "\n", sep = "")
cat("CORRECTED_CLASS=", classification, "\n", sep = "")
