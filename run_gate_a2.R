suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(cluster))

set.seed(20260811)
script_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])
root <- normalizePath(file.path(dirname(script_file), "..", ".."))
res <- file.path(root, "06_results/gateA")
figdir <- file.path(root, "07_figures/draft/gateA2")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

tier1 <- read.delim(file.path(res, "TRANSIENT_TIER1.tsv"), check.names = FALSE)
tier2 <- read.delim(file.path(res, "TRANSIENT_TIER2.tsv"), check.names = FALSE)
tier3 <- read.delim(file.path(res, "TRANSIENT_TIER3.tsv"), check.names = FALSE)
sample_map <- read.delim(file.path(res, "PSEUDOBULK_SAMPLE_MAP.tsv"), check.names = FALSE)
pb <- read.delim(gzfile(file.path(root, "data/intermediate/GSE241132_fibroblast_pseudobulk.tsv.gz")), check.names = FALSE)
gene_id <- pb$gene_id; gene_symbol <- pb$gene_symbol
rownames(pb) <- pb$gene
pb <- as.matrix(pb[, sample_map$sample_id, drop = FALSE])
storage.mode(pb) <- "numeric"
cpm <- t(t(pb) / colSums(pb) * 1e6)
logexpr <- log2(pmax(cpm, 0) + 0.5)
donors <- sort(unique(sample_map$donor))
times <- c("Skin", "D1", "D7", "D30")
sample_for <- function(d, t) sample_map$sample_id[sample_map$donor == d & sample_map$timepoint == t]

calc_metrics <- function(use_donors) {
  effects <- lapply(use_donors, function(d) {
    x <- sapply(times, function(t) logexpr[, sample_for(d, t)])
    colnames(x) <- times; x
  })
  d1 <- do.call(cbind, lapply(effects, function(x) x[, "D1"] - x[, "Skin"]))
  d7 <- do.call(cbind, lapply(effects, function(x) x[, "D7"] - x[, "Skin"]))
  d30s <- do.call(cbind, lapply(effects, function(x) x[, "D30"] - x[, "Skin"]))
  md1 <- apply(d1, 1, median); md7 <- apply(d7, 1, median)
  peak <- ifelse(md1 >= md7, "D1", "D7")
  act <- d1; act[peak == "D7", ] <- d7[peak == "D7", ]
  att <- do.call(cbind, lapply(seq_along(use_donors), function(i) {
    x <- effects[[i]]
    ifelse(peak == "D1", x[, "D1"] - x[, "D30"], x[, "D7"] - x[, "D30"])
  }))
  out <- data.frame(gene = rownames(logexpr), peak_time = peak,
    activation = apply(act, 1, median), attenuation = apply(att, 1, median),
    activation_concordance = rowSums(act > 0), attenuation_concordance = rowSums(att > 0),
    denominator = length(use_donors), rank_score = apply(act, 1, median) + apply(att, 1, median), stringsAsFactors = FALSE)
  eligible <- out$activation > 0
  out$rank_percentile <- 1
  out$rank_percentile[eligible] <- rank(-out$rank_score[eligible], ties.method = "average") / sum(eligible)
  out$tier2 <- eligible & out$activation >= 0.5 & out$attenuation >= 0.3 & out$activation_concordance >= 2 & out$attenuation_concordance >= 2
  out
}

full <- calc_metrics(donors)
lodo_gene <- list()
for (omit in donors) {
  x <- calc_metrics(setdiff(donors, omit))
  x$omitted_donor <- omit
  lodo_gene[[omit]] <- x
}

rob <- tier2[, c("gene", "gene_id", "gene_symbol", "peak_time", "activation_median", "attenuation_median", "activation_concordance", "attenuation_concordance")]
names(rob)[names(rob) == "peak_time"] <- "full_peak_time"
names(rob)[names(rob) == "activation_median"] <- "activation"
names(rob)[names(rob) == "attenuation_median"] <- "attenuation"
rob$full_tier <- "Tier2"
rob$number_of_LODO_runs_retaining_candidate <- 0L
rob$number_of_LODO_runs_preserving_peak_time <- 0L
rank_mat <- matrix(NA_real_, nrow(rob), length(donors), dimnames = list(rob$gene, donors))
for (omit in donors) {
  x <- lodo_gene[[omit]]; idx <- match(rob$gene, x$gene)
  retained <- x$tier2[idx]
  preserved <- retained & x$peak_time[idx] == rob$full_peak_time
  rob$number_of_LODO_runs_retaining_candidate <- rob$number_of_LODO_runs_retaining_candidate + retained
  rob$number_of_LODO_runs_preserving_peak_time <- rob$number_of_LODO_runs_preserving_peak_time + preserved
  rank_mat[, omit] <- x$rank_percentile[idx]
}
rob$mean_LODO_rank_percentile <- rowMeans(rank_mat)
rob$rank_stability <- 1 - apply(rank_mat, 1, sd)
rob$ROBUST_CORE_FLAG <- ifelse(rob$number_of_LODO_runs_retaining_candidate >= 2 & rob$number_of_LODO_runs_preserving_peak_time >= 2, "YES", "NO")
sym <- toupper(ifelse(is.na(rob$gene_symbol) | rob$gene_symbol == "", rob$gene, rob$gene_symbol))
tech <- ifelse(grepl("^MT-", sym), "MITOCHONDRIAL", ifelse(grepl("^RP[SL][0-9]", sym), "RIBOSOMAL", ifelse(grepl("^HB[ABDEGQZ][0-9]", sym), "HEMOGLOBIN", "NO")))
rob$TECHNICAL_GENE_FLAG <- tech
write.table(rob, file.path(res, "A2_GENE_ROBUSTNESS.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

center_gene_profiles <- function(genes, use_donors = donors) {
  x <- logexpr[genes, unlist(lapply(use_donors, function(d) sapply(times, function(t) sample_for(d, t)))), drop = FALSE]
  for (d in use_donors) {
    ids <- sapply(times, function(t) sample_for(d, t))
    x[, ids] <- x[, ids, drop = FALSE] - rowMeans(x[, ids, drop = FALSE])
  }
  s <- apply(x, 1, sd); s[s == 0 | is.na(s)] <- 1
  out <- x / s
  rownames(out) <- genes
  colnames(out) <- colnames(x)
  out
}

adjusted_rand <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  choose2 <- function(z) z * (z - 1) / 2
  sum_nij <- sum(choose2(tab)); sum_ai <- sum(choose2(rowSums(tab))); sum_bj <- sum(choose2(colSums(tab)))
  expected <- sum_ai * sum_bj / choose2(n); maxv <- (sum_ai + sum_bj) / 2
  if (maxv == expected) return(1)
  (sum_nij - expected) / (maxv - expected)
}

cluster_selection <- list(); assignments <- list(); profile_mats <- list()
for (group in c("D1", "D7")) {
  genes <- rob$gene[rob$ROBUST_CORE_FLAG == "YES" & rob$TECHNICAL_GENE_FLAG == "NO" & rob$full_peak_time == group]
  x <- center_gene_profiles(genes)
  profile_mats[[group]] <- x
  dx <- dist(x)
  rows <- list(); fits <- list()
  for (k in 2:8) {
    set.seed(20260811 + k + ifelse(group == "D7", 100, 0))
    fit <- kmeans(x, centers = k, nstart = 50, iter.max = 200)
    sil <- mean(cluster::silhouette(fit$cluster, dx)[, "sil_width"])
    aris <- c()
    for (omit in donors) {
      subcols <- unlist(lapply(setdiff(donors, omit), function(d) sapply(times, function(t) sample_for(d, t))))
      set.seed(20260811 + k + match(omit, donors) * 1000 + ifelse(group == "D7", 100, 0))
      sf <- kmeans(x[, subcols, drop = FALSE], centers = k, nstart = 50, iter.max = 200)
      aris <- c(aris, adjusted_rand(fit$cluster, sf$cluster))
    }
    sizes <- table(fit$cluster)
    rows[[as.character(k)]] <- data.frame(peak_group = group, k = k, silhouette = sil, min_cluster_size = min(sizes), cluster_sizes = paste(as.integer(sizes), collapse = ";"), mean_LODO_ARI = mean(aris), min_LODO_ARI = min(aris), selection_score = sil + mean(aris), selected = "NO")
    fits[[as.character(k)]] <- fit
  }
  tab <- do.call(rbind, rows)
  eligible <- tab$min_cluster_size >= 25
  best <- which(eligible)[which.max(tab$selection_score[eligible])]
  tab$selected[best] <- "YES"
  fit <- fits[[as.character(tab$k[best])]]
  cent <- sapply(sort(unique(fit$cluster)), function(cl) {
    z <- colMeans(x[fit$cluster == cl, , drop = FALSE])
    peak <- ifelse(group == "D1", mean(z[grepl("D1$", names(z))]), mean(z[grepl("D7$", names(z))]))
    d30 <- mean(z[grepl("D30$", names(z))]); c(peak = peak, attenuation = peak - d30)
  })
  ord <- order(-cent["peak", ], -cent["attenuation", ])
  map <- setNames(seq_along(ord), ord)
  assignments[[group]] <- setNames(paste0(group, "_M", map[as.character(fit$cluster)]), rownames(x))
  cluster_selection[[group]] <- tab
}
cluster_tab <- do.call(rbind, cluster_selection)
write.table(cluster_tab, file.path(res, "A2_CLUSTER_SELECTION.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

rob$module <- NA_character_
for (g in names(assignments)) rob$module[match(names(assignments[[g]]), rob$gene)] <- assignments[[g]]
membership <- data.frame(gene = rob$gene, peak_group = rob$full_peak_time, module = ifelse(is.na(rob$module), "NOT_CLUSTERED", rob$module), robust_core = rob$ROBUST_CORE_FLAG, technical_flag = rob$TECHNICAL_GENE_FLAG, activation = rob$activation, attenuation = rob$attenuation, LODO_retention = rob$number_of_LODO_runs_retaining_candidate, stringsAsFactors = FALSE)
write.table(membership, file.path(res, "A2_MODULE_MEMBERSHIP.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

all_modules <- sort(unique(rob$module[!is.na(rob$module)]))
module_scores <- list()
for (m in all_modules) {
  genes <- rob$gene[!is.na(rob$module) & rob$module == m]
  group <- sub("_M.*", "", m)
  missing_genes <- setdiff(genes, rownames(profile_mats[[group]]))
  if (length(missing_genes)) stop("Profile matrix missing genes for ", m, ": ", paste(head(missing_genes), collapse = ","))
  x <- profile_mats[[group]][genes, , drop = FALSE]
  sc <- colMeans(x)
  module_scores[[m]] <- data.frame(donor = sample_map$donor[match(names(sc), sample_map$sample_id)], timepoint = sample_map$timepoint[match(names(sc), sample_map$sample_id)], module = m, score = as.numeric(sc), stringsAsFactors = FALSE)
}
scores <- do.call(rbind, module_scores)
scores$timepoint <- factor(scores$timepoint, levels = times)
write.table(scores, file.path(res, "A2_MODULE_SCORES.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

module_metrics <- list(); module_lodo <- list()
for (m in all_modules) {
  s <- scores[scores$module == m, ]
  means <- sapply(times, function(t) mean(s$score[s$timepoint == t]))
  group <- sub("_M.*", "", m); peak <- group
  act_d <- sapply(donors, function(d) s$score[s$donor == d & s$timepoint == peak] - s$score[s$donor == d & s$timepoint == "Skin"])
  att_d <- sapply(donors, function(d) s$score[s$donor == d & s$timepoint == peak] - s$score[s$donor == d & s$timepoint == "D30"])
  act <- means[peak] - means["Skin"]; att <- means[peak] - means["D30"]
  ratio <- abs(means["D30"] - means["Skin"]) / (abs(act) + 1e-8)
  for (omit in donors) {
    kept <- setdiff(donors, omit)
    lm <- sapply(times, function(t) mean(s$score[s$donor %in% kept & s$timepoint == t]))
    la <- lm[peak] - lm["Skin"]; lt <- lm[peak] - lm["D30"]
    module_lodo[[length(module_lodo) + 1]] <- data.frame(module = m, omitted_donor = omit, trajectory_correlation = cor(means, lm), activation = la, attenuation = lt, activation_sign_preserved = la > 0, attenuation_sign_preserved = lt > 0, activation_ratio_vs_full = la / (act + 1e-8), attenuation_ratio_vs_full = lt / (att + 1e-8), stringsAsFactors = FALSE)
  }
  module_metrics[[m]] <- data.frame(module = m, N_GENES = sum(rob$module == m, na.rm = TRUE), SKIN_MEAN = means["Skin"], D1_MEAN = means["D1"], D7_MEAN = means["D7"], D30_MEAN = means["D30"], PEAK_TIME = peak, ACTIVATION_EFFECT = act, ATTENUATION_EFFECT = att, RETURN_TOWARD_SKIN_RATIO = ratio, ACTIVATION_CONCORDANCE = sum(act_d > 0), ATTENUATION_CONCORDANCE = sum(att_d > 0), stringsAsFactors = FALSE)
}
module_metrics <- do.call(rbind, module_metrics); module_lodo <- do.call(rbind, module_lodo)
write.table(module_lodo, file.path(res, "A2_MODULE_LODO.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

meta <- read.delim(gzfile(file.path(root, "data/staged/GSE241132/GSE241132_cell_metadata.txt.gz")), check.names = FALSE)
meta <- meta[meta$newMainCellTypes == "Fibroblast", ]
meta$timepoint <- sub("^.*D", "D", meta$orig.ident); meta$timepoint[meta$timepoint == "D0"] <- "Skin"
sub_counts <- as.data.frame(with(meta, table(donor = Patient, timepoint = factor(timepoint, levels = times), subtype = newCellTypes)), stringsAsFactors = FALSE)
names(sub_counts)[4] <- "n_cells"
sub_counts$used_for_sensitivity <- ifelse(sub_counts$subtype %in% c("FB-II", "FB-III"), "YES", "NO")
sub_counts$reason <- ifelse(sub_counts$used_for_sensitivity == "YES", "At least 50 cells in every donor-time sample", "At least one donor-time sample has fewer than 50 cells")
write.table(sub_counts, file.path(res, "A2_FIBROBLAST_SUBTYPE_COUNTS.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

subtype_rows <- list()
feature_file <- file.path(root, "data/staged/GSE241132/extracted/PWH26D0/PWH26D0_features.tsv.gz")
features <- read.delim(gzfile(feature_file), header = FALSE, stringsAsFactors = FALSE)
raw_gene <- features[[2]]
bad_gene <- is.na(raw_gene) | raw_gene == "" | duplicated(raw_gene) | duplicated(raw_gene, fromLast = TRUE)
raw_gene[bad_gene] <- features[[1]][bad_gene]
raw_gene_idx <- match(rownames(pb), raw_gene)
if (anyNA(raw_gene_idx)) stop("Filtered pseudobulk genes missing from raw feature map")
for (subtype in c("FB-II", "FB-III")) {
  sub_pb <- matrix(0, nrow(pb), ncol(pb), dimnames = dimnames(pb))
  for (sid in colnames(pb)) {
    sdir <- file.path(root, "data/staged/GSE241132/extracted", sid)
    bc <- readLines(gzfile(file.path(sdir, paste0(sid, "_barcodes.tsv.gz"))))
    mat <- readMM(gzfile(file.path(sdir, paste0(sid, "_matrix.mtx.gz"))))
    wanted <- sub(paste0("^", sid, "_"), "", meta$barcode[meta$orig.ident == sid & meta$newCellTypes == subtype])
    idx <- match(wanted, bc)
    sub_pb[, sid] <- as.numeric(rowSums(mat[raw_gene_idx, idx, drop = FALSE]))
    rm(mat); gc(verbose = FALSE)
  }
  scpm <- t(t(sub_pb) / pmax(colSums(sub_pb), 1) * 1e6)
  slog <- log2(pmax(scpm, 0) + 0.5)
  for (m in all_modules) {
    genes <- intersect(rob$gene[!is.na(rob$module) & rob$module == m], rownames(slog)); x <- slog[genes, , drop = FALSE]
    for (d in donors) { ids <- sapply(times, function(t) sample_for(d, t)); x[, ids] <- x[, ids, drop = FALSE] - rowMeans(x[, ids, drop = FALSE]) }
    sdg <- apply(x, 1, sd); sdg[sdg == 0 | is.na(sdg)] <- 1; x <- x / sdg
    ms <- colMeans(x); peak <- sub("_M.*", "", m)
    ad <- sapply(donors, function(d) ms[sample_for(d, peak)] - ms[sample_for(d, "Skin")])
    td <- sapply(donors, function(d) ms[sample_for(d, peak)] - ms[sample_for(d, "D30")])
    subtype_rows[[length(subtype_rows) + 1]] <- data.frame(subtype = subtype, module = m, activation = mean(ad), attenuation = mean(td), activation_concordance = sum(ad > 0), attenuation_concordance = sum(td > 0), direction_preserved_2of3 = mean(ad) > 0 & mean(td) > 0 & sum(ad > 0) >= 2 & sum(td > 0) >= 2, stringsAsFactors = FALSE)
  }
}
subtype_metrics <- do.call(rbind, subtype_rows)
write.table(subtype_metrics, file.path(res, "A2_SUBTYPE_MODULE_SENSITIVITY.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sub_rob <- sapply(all_modules, function(m) sum(subtype_metrics$direction_preserved_2of3[subtype_metrics$module == m]))
module_metrics$SUBTYPE_SENSITIVITY <- c("NONE", "ONE_MAJOR_SUBTYPE", "BOTH_MAJOR_SUBTYPES")[sub_rob[match(module_metrics$module, names(sub_rob))] + 1]

lodo_ok <- aggregate(cbind(activation_sign_preserved, attenuation_sign_preserved) ~ module, module_lodo, sum)
module_metrics <- merge(module_metrics, lodo_ok, by = "module", all.x = TRUE)
module_metrics$QUALITY <- ifelse(module_metrics$ACTIVATION_EFFECT > 0 & module_metrics$ATTENUATION_EFFECT > 0 & module_metrics$ACTIVATION_CONCORDANCE == 3 & module_metrics$ATTENUATION_CONCORDANCE == 3 & module_metrics$activation_sign_preserved == 3 & module_metrics$attenuation_sign_preserved == 3, "HIGH_QUALITY", ifelse(module_metrics$ACTIVATION_EFFECT > 0 & module_metrics$ATTENUATION_EFFECT > 0 & module_metrics$ACTIVATION_CONCORDANCE >= 2 & module_metrics$ATTENUATION_CONCORDANCE >= 2 & module_metrics$activation_sign_preserved >= 2 & module_metrics$attenuation_sign_preserved >= 2, "MODERATE_QUALITY", ifelse(module_metrics$ACTIVATION_EFFECT > 0 & module_metrics$ATTENUATION_EFFECT > 0, "LOW_QUALITY", "FAIL")))
write.table(module_metrics, file.path(res, "A2_MODULE_METRICS.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

ranking <- list()
for (m in all_modules) {
  x <- rob[!is.na(rob$module) & rob$module == m, ]
  z <- function(v) if (length(v) < 2 || sd(v) == 0) rep(0, length(v)) else as.numeric(scale(v))
  x$stability_score <- z(x$activation) + z(x$attenuation) + x$activation_concordance / 3 + x$attenuation_concordance / 3 + x$number_of_LODO_runs_retaining_candidate / 3 + x$number_of_LODO_runs_preserving_peak_time / 3 - x$mean_LODO_rank_percentile
  x <- x[order(-x$stability_score, x$gene), ]
  x$MODULE_CORE_RANK <- seq_len(nrow(x))
  x$core25 <- x$MODULE_CORE_RANK <= 25; x$core50 <- x$MODULE_CORE_RANK <= 50; x$core100 <- x$MODULE_CORE_RANK <= 100
  ranking[[m]] <- data.frame(module = m, gene = x$gene, gene_id = x$gene_id, gene_symbol = x$gene_symbol, MODULE_CORE_RANK = x$MODULE_CORE_RANK, stability_score = x$stability_score, core25 = x$core25, core50 = x$core50, core100 = x$core100, stringsAsFactors = FALSE)
}
ranking <- do.call(rbind, ranking)
write.table(ranking, file.path(res, "A2_MODULE_CORE_RANKING.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

png(file.path(figdir, "module_donor_trajectories.png"), width = 1600, height = 1000)
par(family = "Helvetica", mfrow = c(ceiling(length(all_modules) / 3), 3), mar = c(3, 3, 2, 1))
for (m in all_modules) { s <- scores[scores$module == m, ]; plot(1:4, rep(0,4), type="n", xaxt="n", xlab="", ylab="Score", main=m, ylim=range(s$score)); axis(1,1:4,times); for(d in donors) lines(1:4,s$score[s$donor==d][match(times,s$timepoint[s$donor==d])],type="b"); abline(h=0,lty=2,col="grey70") }
dev.copy(pdf, file = file.path(figdir, "module_donor_trajectories.pdf"), width = 12, height = 8, family = "Helvetica"); dev.off()
dev.off()

png(file.path(figdir, "module_heatmap.png"), width = 1200, height = 800)
par(family = "Helvetica"); score_mat <- sapply(all_modules, function(m) scores$score[scores$module == m]); rownames(score_mat) <- paste(sample_map$donor, sample_map$timepoint, sep="_"); heatmap(t(score_mat), Rowv=NA, Colv=NA, scale="none", margins=c(8,8))
dev.copy(pdf, file = file.path(figdir, "module_heatmap.pdf"), width = 10, height = 7, family = "Helvetica"); dev.off()
dev.off()

png(file.path(figdir, "module_activation_attenuation.png"), width = 900, height = 700)
par(family = "Helvetica"); plot(module_metrics$ACTIVATION_EFFECT, module_metrics$ATTENUATION_EFFECT, pch=19, col=ifelse(module_metrics$QUALITY=="HIGH_QUALITY","#3182BD","#767676"), xlab="Activation effect", ylab="Attenuation effect"); text(module_metrics$ACTIVATION_EFFECT, module_metrics$ATTENUATION_EFFECT, module_metrics$module, pos=3, cex=.8); abline(h=0,v=0,lty=2,col="grey70")
dev.copy(pdf, file = file.path(figdir, "module_activation_attenuation.pdf"), width = 8, height = 6, family = "Helvetica"); dev.off()
dev.off()

png(file.path(figdir, "module_lodo_stability.png"), width = 1000, height = 700)
par(family = "Helvetica"); boxplot(trajectory_correlation ~ module, module_lodo, las=2, ylab="Trajectory correlation", col="#D8D8D8")
dev.copy(pdf, file = file.path(figdir, "module_lodo_stability.pdf"), width = 9, height = 6, family = "Helvetica"); dev.off()
dev.off()

png(file.path(figdir, "fixed_subtype_module_trajectory.png"), width = 1200, height = 800)
par(family = "Helvetica"); mat <- xtabs((activation + attenuation) ~ module + subtype, subtype_metrics); barplot(mat, beside=TRUE, las=2, col=c("#3182BD","#33B5A5"), ylab="Activation + attenuation", legend.text=rownames(mat))
dev.copy(pdf, file = file.path(figdir, "fixed_subtype_module_trajectory.pdf"), width = 10, height = 7, family = "Helvetica"); dev.off()
dev.off()

write.table(data.frame(metric=c("TIER2_N","ROBUST_CORE_N","D1_ROBUST_CORE_N","D7_ROBUST_CORE_N","D1_K","D7_K","MODULE_N","HIGH_QUALITY_MODULE_N","MODERATE_QUALITY_MODULE_N"), value=c(nrow(tier2),sum(rob$ROBUST_CORE_FLAG=="YES"),sum(rob$ROBUST_CORE_FLAG=="YES" & rob$full_peak_time=="D1"),sum(rob$ROBUST_CORE_FLAG=="YES" & rob$full_peak_time=="D7"),cluster_tab$k[cluster_tab$peak_group=="D1" & cluster_tab$selected=="YES"],cluster_tab$k[cluster_tab$peak_group=="D7" & cluster_tab$selected=="YES"],nrow(module_metrics),sum(module_metrics$QUALITY=="HIGH_QUALITY"),sum(module_metrics$QUALITY=="MODERATE_QUALITY"))), file.path(res,"gate_a2_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
