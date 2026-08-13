suppressPackageStartupMessages(library(Matrix))

script_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])
root <- normalizePath(file.path(dirname(script_file), "..", ".."))
meta_file <- file.path(root, "data/staged/GSE241132/GSE241132_cell_metadata.txt.gz")
matrix_root <- file.path(root, "data/staged/GSE241132/extracted")
result_dir <- file.path(root, "06_results/gateA")
figure_dir <- file.path(root, "07_figures/draft/gateA1")
intermediate_dir <- file.path(root, "data/intermediate")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read.delim(gzfile(meta_file), check.names = FALSE, stringsAsFactors = FALSE)
meta$timepoint <- sub("^.*D", "D", meta$orig.ident)
meta$timepoint[meta$timepoint == "D0"] <- "Skin"
time_levels <- c("Skin", "D1", "D7", "D30")
meta$timepoint <- factor(meta$timepoint, levels = time_levels)
meta$donor <- meta$Patient

label_counts <- as.data.frame(table(meta$newCellTypes, meta$newMainCellTypes), stringsAsFactors = FALSE)
label_counts <- label_counts[label_counts$Freq > 0, ]
names(label_counts) <- c("original_label", "main_label", "cells")
label_counts$included_as_fibroblast <- ifelse(label_counts$main_label == "Fibroblast", "YES", "NO")
label_counts$reason <- ifelse(label_counts$main_label == "Fibroblast", "Included by official newMainCellTypes annotation", "Excluded: official main label is not Fibroblast")
label_counts <- label_counts[order(label_counts$included_as_fibroblast, decreasing = TRUE), ]
write.table(label_counts[, c("original_label", "cells", "included_as_fibroblast", "reason")], file.path(result_dir, "GSE241132_FIBROBLAST_LABEL_AUDIT.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

fib <- meta[meta$newMainCellTypes == "Fibroblast", ]
sample_ids <- sort(unique(meta$orig.ident))
pb_list <- list()
sample_map <- list()
feature_ref <- NULL

for (sample_id in sample_ids) {
  sample_dir <- file.path(matrix_root, sample_id)
  barcodes <- readLines(gzfile(file.path(sample_dir, paste0(sample_id, "_barcodes.tsv.gz"))))
  features <- read.delim(gzfile(file.path(sample_dir, paste0(sample_id, "_features.tsv.gz"))), header = FALSE, stringsAsFactors = FALSE)
  if (is.null(feature_ref)) feature_ref <- features
  if (!identical(features[[1]], feature_ref[[1]])) stop("Feature order mismatch: ", sample_id)
  mat <- readMM(gzfile(file.path(sample_dir, paste0(sample_id, "_matrix.mtx.gz"))))
  keep_barcodes <- sub(paste0("^", sample_id, "_"), "", fib$barcode[fib$orig.ident == sample_id])
  idx <- match(keep_barcodes, barcodes)
  if (anyNA(idx)) stop("Fibroblast barcode mismatch: ", sample_id)
  pb_list[[sample_id]] <- as.numeric(Matrix::rowSums(mat[, idx, drop = FALSE]))
  sm <- fib[fib$orig.ident == sample_id, ][1, ]
  sample_map[[sample_id]] <- data.frame(sample_id = sample_id, donor = sm$donor, timepoint = as.character(sm$timepoint), n_fibroblast_cells = length(idx), input_type = "raw integer UMI counts", stringsAsFactors = FALSE)
  rm(mat); gc(verbose = FALSE)
}

pb <- do.call(cbind, pb_list)
colnames(pb) <- names(pb_list)
gene_id <- feature_ref[[1]]
gene_symbol <- feature_ref[[2]]
gene <- gene_symbol
bad <- is.na(gene) | gene == "" | duplicated(gene) | duplicated(gene, fromLast = TRUE)
gene[bad] <- gene_id[bad]
rownames(pb) <- gene
sample_map <- do.call(rbind, sample_map)
sample_map$low_cell_count_warning <- ifelse(sample_map$n_fibroblast_cells < 50, "YES", "NO")
sample_map <- sample_map[order(sample_map$donor, match(sample_map$timepoint, time_levels)), ]
write.table(sample_map, file.path(result_dir, "PSEUDOBULK_SAMPLE_MAP.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sample_map[, c("donor", "timepoint", "n_fibroblast_cells", "sample_id", "low_cell_count_warning")], file.path(result_dir, "FIBROBLAST_DONOR_TIME_COUNTS.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("donor", "timepoint", "n_fibroblast_cells", "library_or_sample_id", "low_cell_count_warning"))

pb <- pb[, sample_map$sample_id, drop = FALSE]
libsize <- colSums(pb)
cpm <- t(t(pb) / libsize * 1e6)
keep <- rowSums(cpm >= 1) >= 3
pb_f <- pb[keep, , drop = FALSE]
logexpr <- log2(pmax(cpm[keep, , drop = FALSE], 0) + 0.5)
pb_out <- data.frame(gene = rownames(pb_f), gene_id = gene_id[keep], gene_symbol = gene_symbol[keep], pb_f, check.names = FALSE)
write.table(pb_out, gzfile(file.path(intermediate_dir, "GSE241132_fibroblast_pseudobulk.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)

donors <- sort(unique(sample_map$donor))
sample_for <- function(d, t) sample_map$sample_id[sample_map$donor == d & sample_map$timepoint == t]

calc_metrics <- function(use_donors) {
  effects <- lapply(use_donors, function(d) {
    x <- sapply(time_levels, function(t) logexpr[, sample_for(d, t)])
    colnames(x) <- time_levels
    x
  })
  d1 <- do.call(cbind, lapply(effects, function(x) x[, "D1"] - x[, "Skin"]))
  d7 <- do.call(cbind, lapply(effects, function(x) x[, "D7"] - x[, "Skin"]))
  d30s <- do.call(cbind, lapply(effects, function(x) x[, "D30"] - x[, "Skin"]))
  med_d1 <- apply(d1, 1, median); med_d7 <- apply(d7, 1, median)
  peak_time <- ifelse(med_d1 >= med_d7, "D1", "D7")
  act <- d1
  act[peak_time == "D7", ] <- d7[peak_time == "D7", ]
  att <- do.call(cbind, lapply(seq_along(use_donors), function(i) {
    x <- effects[[i]]
    ifelse(peak_time == "D1", x[, "D1"] - x[, "D30"], x[, "D7"] - x[, "D30"])
  }))
  peak_abs <- abs(act)
  d30_abs <- abs(d30s)
  data.frame(
    gene = rownames(logexpr), peak_time = peak_time,
    activation_median = apply(act, 1, median), activation_mean = rowMeans(act), activation_sd = apply(act, 1, sd),
    activation_concordance = rowSums(act > 0), activation_denominator = length(use_donors),
    attenuation_median = apply(att, 1, median), attenuation_mean = rowMeans(att), attenuation_sd = apply(att, 1, sd),
    attenuation_concordance = rowSums(att > 0), attenuation_denominator = length(use_donors),
    D30_vs_skin_median = apply(d30s, 1, median),
    peak_to_skin_abs_distance = apply(peak_abs, 1, median),
    D30_to_skin_abs_distance = apply(d30_abs, 1, median),
    return_toward_skin_ratio = apply(d30_abs, 1, median) / (apply(peak_abs, 1, median) + 1e-8),
    stringsAsFactors = FALSE
  )
}

metrics <- calc_metrics(donors)
metrics$gene_id <- gene_id[keep][match(metrics$gene, rownames(logexpr))]
metrics$gene_symbol <- gene_symbol[keep][match(metrics$gene, rownames(logexpr))]
metrics$eligible_positive_peak <- metrics$activation_median > 0
metrics$rank_score <- metrics$activation_median + metrics$attenuation_median
metrics <- metrics[order(-metrics$rank_score, metrics$gene), ]
write.table(metrics, file.path(result_dir, "GENE_TRANSIENT_METRICS.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

tier_set <- function(m, tier) {
  if (tier == 1) keep_t <- m$activation_median >= 0.5 & m$attenuation_median >= 0.5 & m$activation_concordance == m$activation_denominator & m$attenuation_concordance == m$attenuation_denominator
  if (tier == 2) keep_t <- m$activation_median >= 0.5 & m$attenuation_median >= 0.3 & m$activation_concordance >= 2 & m$attenuation_concordance >= 2
  if (tier == 3) keep_t <- m$activation_median >= 0.3 & m$attenuation_median >= 0.3 & m$activation_concordance >= 2 & m$attenuation_concordance >= 2
  m[keep_t & m$eligible_positive_peak, , drop = FALSE]
}
tiers <- lapply(1:3, function(i) tier_set(metrics, i))
for (i in 1:3) write.table(tiers[[i]], file.path(result_dir, paste0("TRANSIENT_TIER", i, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

candidate <- tiers[[2]]
candidate <- candidate[, c("gene", "gene_id", "gene_symbol", "peak_time", "activation_median", "activation_mean", "activation_concordance", "attenuation_median", "attenuation_mean", "attenuation_concordance", "D30_vs_skin_median", "return_toward_skin_ratio", "rank_score")]
con <- file(file.path(result_dir, "TWPS_CANDIDATE_NOT_LOCKED.tsv"), "w")
writeLines(c("# STATUS=NOT_LOCKED", "# VALIDATION_DATA_NOT_VIEWED=YES"), con)
write.table(candidate, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)

full_sets <- split(tiers[[2]]$gene, tiers[[2]]$peak_time)
lodo_rows <- list()
for (omit in donors) {
  sub <- setdiff(donors, omit)
  lm <- calc_metrics(sub)
  lm$eligible_positive_peak <- lm$activation_median > 0
  lm$rank_score <- lm$activation_median + lm$attenuation_median
  lt <- tier_set(lm, 2)
  lt <- lt[order(-lt$rank_score, lt$gene), ]
  for (program in c("D1", "D7")) {
    full <- full_sets[[program]]; if (is.null(full)) full <- character()
    cur <- lt$gene[lt$peak_time == program]
    ov <- intersect(full, cur); un <- union(full, cur)
    lodo_rows[[length(lodo_rows) + 1]] <- data.frame(
      omitted_donor = omit, retained_donors = paste(sub, collapse = "+"), program = program,
      full_n = length(full), lodo_n = length(cur), overlap_n = length(ov),
      jaccard = ifelse(length(un) == 0, 1, length(ov) / length(un)),
      full_top20_overlap_n = length(intersect(head(full, 20), head(cur, 20))),
      lodo_top_genes = paste(head(cur, 10), collapse = ";"), stringsAsFactors = FALSE)
  }
}
lodo <- do.call(rbind, lodo_rows)
write.table(lodo, file.path(result_dir, "LODO_STABILITY.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

png(file.path(figure_dir, "fibroblast_cell_counts.png"), width = 1200, height = 700)
par(family = "sans")
cols <- c(Skin = "grey50", D1 = "#d73027", D7 = "#fc8d59", D30 = "#4575b4")
barplot(sample_map$n_fibroblast_cells, names.arg = paste(sample_map$donor, sample_map$timepoint, sep = "\n"), col = cols[sample_map$timepoint], las = 2, ylab = "Fibroblast cells")
dev.off()

png(file.path(figure_dir, "pseudobulk_sample_correlation.png"), width = 900, height = 800)
par(family = "sans")
cor_mat <- cor(logexpr)
heatmap(cor_mat, Rowv = NA, Colv = NA, scale = "none", margins = c(10, 10))
dev.off()

plot_top <- function(program, filename) {
  genes <- head(tiers[[2]]$gene[tiers[[2]]$peak_time == program], 20)
  png(file.path(figure_dir, filename), width = 1100, height = 800)
  par(family = "sans")
  if (length(genes) == 0) { plot.new(); text(.5, .5, paste("No Tier2", program, "candidates")) } else {
    vals <- sapply(time_levels, function(t) {
      ids <- sample_map$sample_id[match(paste(donors, t), paste(sample_map$donor, sample_map$timepoint))]
      apply(logexpr[genes, ids, drop = FALSE], 1, median)
    })
    vals <- vals - vals[, "Skin"]
    matplot(seq_along(time_levels), t(vals), type = "l", lty = 1, xlab = "Time", ylab = "Median log2 expression change vs Skin", xaxt = "n")
    axis(1, seq_along(time_levels), time_levels); abline(h = 0, lty = 2, col = "grey60")
  }
  dev.off()
}
plot_top("D1", "tier2_D1_top20_trajectories.png")
plot_top("D7", "tier2_D7_top20_trajectories.png")

counts <- sapply(tiers, function(x) c(D1 = sum(x$peak_time == "D1"), D7 = sum(x$peak_time == "D7")))
png(file.path(figure_dir, "candidate_count_sensitivity.png"), width = 900, height = 700)
par(family = "sans")
barplot(counts, beside = TRUE, names.arg = c("Tier1", "Tier2", "Tier3"), col = c("#d73027", "#4575b4"), legend.text = rownames(counts), ylab = "Candidate genes")
dev.off()

summary_file <- file.path(result_dir, "gate_a1_summary.tsv")
summary_rows <- data.frame(
  metric = c("fibroblast_cells", "filtered_genes", "tier1_total", "tier1_D1", "tier1_D7", "tier2_total", "tier2_D1", "tier2_D7", "tier3_total", "tier3_D1", "tier3_D7", "lodo_D1_median_jaccard", "lodo_D7_median_jaccard"),
  value = c(nrow(fib), nrow(logexpr), length(tiers[[1]]$gene), sum(tiers[[1]]$peak_time == "D1"), sum(tiers[[1]]$peak_time == "D7"), length(tiers[[2]]$gene), sum(tiers[[2]]$peak_time == "D1"), sum(tiers[[2]]$peak_time == "D7"), length(tiers[[3]]$gene), sum(tiers[[3]]$peak_time == "D1"), sum(tiers[[3]]$peak_time == "D7"), median(lodo$jaccard[lodo$program == "D1"]), median(lodo$jaccard[lodo$program == "D7"])), stringsAsFactors = FALSE)
write.table(summary_rows, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
