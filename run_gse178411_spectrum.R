#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260811)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grep("^--file=", args)])
root <- normalizePath(file.path(dirname(script_arg), "..", ".."))
res <- file.path(root, "06_results", "spectrum")
fig <- file.path(root, "07_figures", "draft", "spectrum_GSE178411")
dir.create(res, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
fmt <- function(x, d = 10) format(as.numeric(x), digits = d, trim = TRUE)

# Frozen integrity and inputs.
primary_path <- file.path(root, "06_results/gateA/TWPS_PRIMARY_D7_M3_CORE100.tsv")
lock_hash <- strsplit(system2("shasum", c("-a", "256", primary_path), stdout = TRUE), " +")[[1]][1]
stopifnot(lock_hash == "6d5a62c080f2d0d9e4077983294d2a6a443a9e0e3c565691a4d3a863cc2afcc3")
counts_path <- file.path(root, "data/raw_small/GSE178411/GSE178411_counts.txt.gz")
counts <- read.delim(gzfile(counts_path), check.names = FALSE)
stopifnot(ncol(counts) == 108L, nrow(counts) == 28395L, !anyNA(counts), all(counts >= 0), !anyDuplicated(rownames(counts)))
record_map <- read.delim(file.path(res, "GSE178411_RECORD_MAP.tsv"), check.names = FALSE)
meta <- record_map[match(colnames(counts), record_map$sample_id), ]
stopifnot(nrow(meta) == 108L, !anyNA(meta$patient_id_resolved), all(meta$sample_id == colnames(counts)))

# Deterministic NCBI GeneID -> approved symbol mapping; raw duplicate IDs are summed.
gi_path <- file.path(root, "data/intermediate/gateB/Homo_sapiens.gene_info.gz")
gi <- read.delim(gzfile(gi_path), quote = "", comment.char = "", check.names = FALSE)
approved <- ifelse(gi$Symbol_from_nomenclature_authority == "-", gi$Symbol, gi$Symbol_from_nomenclature_authority)
entrez_to_symbol <- setNames(approved, as.character(gi$GeneID))
symbols <- unname(entrez_to_symbol[rownames(counts)])
mapped <- !is.na(symbols) & symbols != "" & symbols != "-"
symbol_counts <- rowsum(as.matrix(counts[mapped, , drop = FALSE]), group = symbols[mapped], reorder = FALSE)
libsize <- colSums(counts)
logcpm <- log2(sweep(symbol_counts, 2, libsize, "/") * 1e6 + 1)

synonym_to_approved <- list()
for (i in seq_len(nrow(gi))) {
  if (!(as.character(gi$GeneID[i]) %in% rownames(counts))) next
  syns <- strsplit(gi$Synonyms[i], "\\|", fixed = FALSE)[[1]]
  syns <- syns[syns != "-"]
  for (s in syns) synonym_to_approved[[s]] <- unique(c(synonym_to_approved[[s]], approved[i]))
}
resolve_signature <- function(genes) {
  out <- lapply(genes, function(g) {
    if (g %in% rownames(logcpm)) {
      ids <- rownames(counts)[mapped & symbols == g]
      data.frame(locked_gene = g, matched_NCBI_gene_id = paste(ids, collapse = ";"), resolved_approved_symbol = g,
                 mapping_type = "APPROVED_SYMBOL", available = "YES", duplicate_mapping = ifelse(length(ids) > 1, "YES", "NO"), final_action = "INCLUDE")
    } else {
      candidates <- intersect(synonym_to_approved[[g]], rownames(logcpm))
      if (length(candidates) == 1L) {
        ids <- rownames(counts)[mapped & symbols == candidates]
        data.frame(locked_gene = g, matched_NCBI_gene_id = paste(ids, collapse = ";"), resolved_approved_symbol = candidates,
                   mapping_type = "UNIQUE_ANNOTATION_SYNONYM", available = "YES", duplicate_mapping = ifelse(length(ids) > 1, "YES", "NO"), final_action = "INCLUDE_AS_LOCKED_GENE")
      } else {
        data.frame(locked_gene = g, matched_NCBI_gene_id = "", resolved_approved_symbol = "",
                   mapping_type = ifelse(length(candidates) > 1, "AMBIGUOUS_SYNONYM", "UNAVAILABLE"), available = "NO", duplicate_mapping = "NO", final_action = "EXCLUDE_NO_SUBSTITUTE")
      }
    }
  })
  x <- do.call(rbind, out)
  resolved <- x$resolved_approved_symbol[x$available == "YES"]
  stopifnot(!anyDuplicated(resolved))
  x
}
read_genes <- function(f) read.delim(file.path(root, "06_results/gateA", f), check.names = FALSE)$gene
genes <- list(CORE100 = read_genes("TWPS_PRIMARY_D7_M3_CORE100.tsv"), CORE50 = read_genes("TWPS_SENSITIVITY_D7_M3_CORE50.tsv"),
              CORE25 = read_genes("TWPS_SENSITIVITY_D7_M3_CORE25.tsv"), FULL = read_genes("TWPS_SENSITIVITY_D7_M3_FULL.tsv"))
maps <- lapply(genes, resolve_signature)
write.table(maps$CORE100[, c("locked_gene", "matched_NCBI_gene_id", "mapping_type", "available", "duplicate_mapping", "final_action")],
            file.path(res, "GSE178411_CORE100_MAPPING.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
coverage <- vapply(maps, function(x) sum(x$available == "YES"), integer(1))
stopifnot(coverage["CORE100"] >= 80L)

zrows <- function(mat) {
  mu <- rowMeans(mat); s <- apply(mat, 1, sd); keep <- is.finite(s) & s > 0
  z <- sweep(mat[keep, , drop = FALSE], 1, mu[keep], "-")
  sweep(z, 1, s[keep], "/")
}
z <- zrows(logcpm)
score_one <- function(map) {
  syms <- map$resolved_approved_symbol[map$available == "YES"]
  colMeans(z[syms, , drop = FALSE])
}
sample_scores <- as.data.frame(lapply(maps, score_one), check.names = FALSE)
sample_scores$sample_id <- rownames(sample_scores)
sample_scores$patient_id <- meta$patient_id_resolved
sample_scores$state <- meta$state_standardized

# Fixed arithmetic mean within patient × state.
key <- paste(sample_scores$patient_id, sample_scores$state, sep = "||")
ps <- do.call(rbind, lapply(split(seq_len(nrow(sample_scores)), key), function(idx) {
  data.frame(patient_id = sample_scores$patient_id[idx[1]], state = sample_scores$state[idx[1]], n_samples = length(idx),
             CORE100 = mean(sample_scores$CORE100[idx]), CORE50 = mean(sample_scores$CORE50[idx]),
             CORE25 = mean(sample_scores$CORE25[idx]), FULL = mean(sample_scores$FULL[idx]))
}))
rownames(ps) <- NULL
write.table(ps, file.path(res, "GSE178411_PATIENT_STATE_TWPS.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# OLS marginal contrast with patient-clustered CR1 standard error.
hedges <- function(a, b) {
  na <- length(a); nb <- length(b)
  sp <- sqrt(((na - 1) * var(a) + (nb - 1) * var(b)) / (na + nb - 2))
  (1 - 3 / (4 * (na + nb - 2) - 1)) * (mean(a) - mean(b)) / sp
}
cluster_contrast <- function(signature, state_a, state_b, expected_sign, label) {
  d <- ps[ps$state %in% c(state_a, state_b), c("patient_id", "state", signature)]
  names(d)[3] <- "y"
  d$x <- as.numeric(d$state == state_a)
  X <- cbind(1, d$x); y <- d$y
  bread <- solve(crossprod(X)); beta <- bread %*% crossprod(X, y); u <- y - as.vector(X %*% beta)
  clusters <- unique(d$patient_id); meat <- matrix(0, 2, 2)
  for (cl in clusters) { ix <- which(d$patient_id == cl); xu <- crossprod(X[ix, , drop = FALSE], u[ix]); meat <- meat + xu %*% t(xu) }
  N <- nrow(d); G <- length(clusters); p <- ncol(X)
  vc <- (G / (G - 1)) * ((N - 1) / (N - p)) * bread %*% meat %*% bread
  se <- sqrt(vc[2, 2]); crit <- qt(.975, df = G - 1); est <- as.numeric(beta[2]); pval <- 2 * pt(-abs(est / se), df = G - 1)
  pa <- unique(d$patient_id[d$state == state_a]); pb <- unique(d$patient_id[d$state == state_b])
  g <- hedges(d$y[d$state == state_a], d$y[d$state == state_b])
  data.frame(contrast = label, signature = signature, n_patients_A = length(pa), n_patients_B = length(pb), patient_overlap = length(intersect(pa, pb)),
             difference = est, ci_low = est - crit * se, ci_high = est + crit * se, effect_size = g, pvalue = pval,
             expected_direction = ifelse(expected_sign > 0, "A>B", "A<B"), observed_direction = ifelse(est > 0, "A>B", ifelse(est < 0, "A<B", "EQUAL")))
}
contrast_specs <- list(
  c("EARLY_WOUND", "UNINJURED_SKIN", 1, "EARLY_VS_UNINJURED"),
  c("LATE_WOUND", "EARLY_WOUND", -1, "LATE_VS_EARLY"),
  c("HYPERTROPHIC_SCAR", "LATE_WOUND", 1, "HTS_VS_LATE"),
  c("HYPERTROPHIC_SCAR", "UNINJURED_SKIN", 1, "HTS_VS_UNINJURED"))
all_results <- do.call(rbind, lapply(names(genes), function(sig) do.call(rbind, lapply(contrast_specs, function(sp) cluster_contrast(sig, sp[1], sp[2], as.numeric(sp[3]), sp[4])))))
classify <- function(row, expected_sign, sens_rows) {
  correct <- sign(row$difference) == expected_sign
  sens_correct <- sum(sign(sens_rows$difference) == expected_sign)
  if (!correct && sens_correct <= 1) return("CONTRADICTORY")
  if (correct && abs(row$effect_size) >= .8 && row$pvalue < .05 && sens_correct == 3) return("STRONG")
  if (correct && abs(row$effect_size) >= .5 && sens_correct >= 2) return("MODERATE")
  if (correct && sens_correct >= 2) return("WEAK")
  "NONE"
}
for (label in unique(all_results$contrast)) {
  ix <- all_results$contrast == label
  expected <- ifelse(label == "LATE_VS_EARLY", -1, 1)
  primary_row <- all_results[ix & all_results$signature == "CORE100", ]
  sens_rows <- all_results[ix & all_results$signature != "CORE100", ]
  ev <- classify(primary_row, expected, sens_rows)
  all_results$evidence_class[ix] <- ifelse(all_results$signature[ix] == "CORE100", ev, "SENSITIVITY")
}
write.table(all_results, file.path(res, "GSE178411_SIGNATURE_SENSITIVITY.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(all_results, file.path(res, "GSE178411_TWPS_SUMMARY.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Days-since-injury: parse metadata only, then fit one prespecified log-time model.
soft <- readLines(gzfile(file.path(root, "data/raw_small/GSE178411/GSE178411_family.soft.gz")), warn = FALSE)
current <- NA_character_; days <- list()
for (line in soft) {
  if (startsWith(line, "^SAMPLE = ")) current <- sub("^\\^SAMPLE = ", "", line)
  if (!is.na(current) && startsWith(line, "!Sample_characteristics_ch1 = days since injury:")) {
    val <- trimws(sub("^.*days since injury:", "", line)); days[[current]] <- suppressWarnings(as.numeric(val))
  }
}
sample_scores$GSM <- meta$GSM_if_available
sample_scores$days_since_injury <- vapply(sample_scores$GSM, function(g) if (is.null(days[[g]])) NA_real_ else days[[g]], numeric(1))
ct <- sample_scores[is.finite(sample_scores$days_since_injury), ]
continuous_text <- "DESCRIPTIVE_CONTINUOUS_TIME_SUPPORT; insufficient valid observations"
continuous <- NULL
if (nrow(ct) >= 10L) {
  X <- cbind(1, log1p(ct$days_since_injury)); y <- ct$CORE100; bread <- solve(crossprod(X)); beta <- bread %*% crossprod(X, y); u <- y - as.vector(X %*% beta)
  cls <- unique(ct$patient_id); meat <- matrix(0, 2, 2)
  for (cl in cls) { ix <- which(ct$patient_id == cl); xu <- crossprod(X[ix, , drop = FALSE], u[ix]); meat <- meat + xu %*% t(xu) }
  N <- nrow(ct); G <- length(cls); vc <- (G/(G-1))*((N-1)/(N-2))*bread%*%meat%*%bread; se <- sqrt(vc[2,2]); p <- 2*pt(-abs(beta[2]/se),G-1); cr <- qt(.975,G-1)
  continuous <- c(n=nrow(ct), missing=nrow(sample_scores)-nrow(ct), patients=G, slope=beta[2], low=beta[2]-cr*se, high=beta[2]+cr*se, p=p)
  continuous_text <- paste0("DESCRIPTIVE_CONTINUOUS_TIME_SUPPORT; log1p(days) slope=", fmt(beta[2]), "; 95% CI=", fmt(continuous['low']), ",", fmt(continuous['high']), "; P=", fmt(p))
}
day_tab <- as.data.frame(table(state = sample_scores$state, valid_days = is.finite(sample_scores$days_since_injury)))
write.table(day_tab, file.path(res, "GSE178411_DAYS_METADATA_COMPLETENESS.tsv"), sep="\t",row.names=FALSE,quote=FALSE)

# Descriptive sparse states.
desc <- do.call(rbind, lapply(c("CHRONIC_WOUND", "NORMAL_SCAR"), function(st) {
  x <- ps$CORE100[ps$state == st]
  data.frame(state=st, patient_n=length(x), median=median(x), minimum=min(x), maximum=max(x), individual_values=paste(fmt(x),collapse=";"))
}))
write.table(desc, file.path(res, "GSE178411_SPARSE_STATES_DESCRIPTIVE.tsv"), sep="\t",row.names=FALSE,quote=FALSE)

# Draft figures.
state_order <- c("UNINJURED_SKIN","EARLY_WOUND","LATE_WOUND","CHRONIC_WOUND","NORMAL_SCAR","HYPERTROPHIC_SCAR")
png(file.path(fig,"A_patient_CORE100_spectrum.png"),1200,750,res=140)
boxplot(CORE100 ~ factor(state,levels=state_order),data=ps,las=2,ylab="Patient-state CORE100 TWPS",xlab="",col="#7AA6C255")
stripchart(CORE100 ~ factor(state,levels=state_order),data=ps,vertical=TRUE,method="jitter",pch=19,add=TRUE,col="#244A66")
dev.off()
prim <- all_results[all_results$signature=="CORE100" & all_results$contrast %in% c("EARLY_VS_UNINJURED","LATE_VS_EARLY","HTS_VS_LATE"),]
png(file.path(fig,"B_primary_effect_forest.png"),950,650,res=140); y <- 3:1
plot(prim$difference,y,xlim=range(c(prim$ci_low,prim$ci_high,0)),yaxt="n",ylab="",xlab="Patient-clustered difference",pch=19); segments(prim$ci_low,y,prim$ci_high,y); abline(v=0,lty=2); axis(2,y,prim$contrast,las=1); dev.off()
png(file.path(fig,"C_signature_direction_consistency.png"),1000,700,res=140)
mat <- xtabs(difference ~ signature + contrast, data=all_results[all_results$contrast %in% c("EARLY_VS_UNINJURED","LATE_VS_EARLY","HTS_VS_LATE"),]); barplot(mat,beside=TRUE,legend.text=TRUE,las=2,ylab="Difference"); abline(h=0,lty=2); dev.off()
if (!is.null(continuous)) { png(file.path(fig,"D_CORE100_vs_days.png"),900,700,res=140); plot(log1p(ct$days_since_injury),ct$CORE100,pch=19,xlab="log1p(days since injury)",ylab="CORE100 TWPS"); abline(lm(CORE100~log1p(days_since_injury),data=ct),col="red",lwd=2); dev.off() }

getr <- function(label, sig="CORE100") all_results[all_results$contrast==label & all_results$signature==sig,][1,]
r1<-getr("EARLY_VS_UNINJURED"); r2<-getr("LATE_VS_EARLY"); r3<-getr("HTS_VS_LATE"); r4<-getr("HTS_VS_UNINJURED")
counts_state <- function(st) c(samples=sum(meta$state_standardized==st),patients=length(unique(meta$patient_id_resolved[meta$state_standardized==st])))
e1<-r1$evidence_class;e2<-r2$evidence_class;e3<-r3$evidence_class
state_labels <- c(UNINJURED="UNINJURED_SKIN",EARLY="EARLY_WOUND",LATE="LATE_WOUND",HTS="HYPERTROPHIC_SCAR",CHRONIC="CHRONIC_WOUND",NORMAL_SCAR="NORMAL_SCAR")
report <- c("# GSE178411 Wound/Scar Spectrum Validation","","## Integrity","","TWPS_SHA256_MATCH=YES","GEO_SAMPLE_N=108","MATRIX_SAMPLE_N=108","MISSING_MATRIX_SAMPLE=NONE","REASON=Phase4A parser correction: the header has no GeneID placeholder and 12-W is the first sample column","NORMALIZATION=raw counts; library-size CPM followed by log2(CPM+1); no differential-expression fitting","PRIMARY_INFERENCE=patient-state arithmetic means with patient-clustered CR1 standard errors","",
"## Cohort used","",paste0("ANALYZED_SAMPLE_N=",nrow(sample_scores)),paste0("ANALYZED_PATIENT_N=",length(unique(sample_scores$patient_id))),"",
unlist(lapply(names(state_labels),function(lbl){st<-state_labels[[lbl]];x<-counts_state(st);c(paste0(lbl,":"),paste0("samples=",x['samples']),paste0("patients=",x['patients']),"")})),
"## Signature coverage","",paste0("CORE100_AVAILABLE=",coverage['CORE100']),paste0("CORE100_COVERAGE=",coverage['CORE100'],"/100"),"",
"## Primary 1 — Early vs Uninjured","",paste0("DIFFERENCE=",fmt(r1$difference)),paste0("95CI=",fmt(r1$ci_low),",",fmt(r1$ci_high)),paste0("EFFECT_SIZE=",fmt(r1$effect_size)),paste0("P=",fmt(r1$pvalue)),paste0("PATIENT_OVERLAP=",r1$patient_overlap),paste0("DIRECTION=",r1$observed_direction),"",
"## Primary 2 — Late vs Early","",paste0("DIFFERENCE=",fmt(r2$difference)),paste0("95CI=",fmt(r2$ci_low),",",fmt(r2$ci_high)),paste0("EFFECT_SIZE=",fmt(r2$effect_size)),paste0("P=",fmt(r2$pvalue)),paste0("PATIENT_OVERLAP=",r2$patient_overlap),paste0("DIRECTION=",r2$observed_direction),"",
"## Pathological support — HTS vs Late","",paste0("DIFFERENCE=",fmt(r3$difference)),paste0("95CI=",fmt(r3$ci_low),",",fmt(r3$ci_high)),paste0("EFFECT_SIZE=",fmt(r3$effect_size)),paste0("P=",fmt(r3$pvalue)),paste0("PATIENT_OVERLAP=",r3$patient_overlap),paste0("DIRECTION=",r3$observed_direction),"",
"## Secondary — HTS vs Uninjured","",paste0("DIFFERENCE=",fmt(r4$difference)),paste0("95CI=",fmt(r4$ci_low),",",fmt(r4$ci_high)),paste0("EFFECT_SIZE=",fmt(r4$effect_size)),paste0("P=",fmt(r4$pvalue)),"",
"## Sensitivity signatures","","CORE25/CORE50/FULL results are stored without selection in GSE178411_SIGNATURE_SENSITIVITY.tsv.","","## Continuous time","",paste0("N_WITH_VALID_DAYS=",ifelse(is.null(continuous),nrow(ct),continuous['n'])),paste0("N_MISSING_DAYS=",nrow(sample_scores)-nrow(ct)),paste0("RESULT=",continuous_text),"",
"## Chronic wound","","DESCRIPTIVE ONLY; see GSE178411_SPARSE_STATES_DESCRIPTIVE.tsv.","","## Normal scar","","DESCRIPTIVE ONLY; see GSE178411_SPARSE_STATES_DESCRIPTIVE.tsv.","",
"## Secondary locked modules","","NOT_RUN; optional and not required for the three CORE100 questions.","","## Optional covariate sensitivity","","SKIPPED; mixed-model packages were unavailable and the primary result was prespecified as minimally adjusted. No covariate model was selected from outcomes.","",
"## Limitations","","- GEO narrative counts are internally inconsistent; record-level metadata were used.","- The matrix header lacks an explicit GeneID placeholder; all 108 sample columns were retained after deterministic parser correction.","- Repeated samples within patients required patient-state aggregation and patient-clustered inference.","- Burn/wound clinical context is heterogeneous.","- Normal scar has two patients and chronic wound has three patients; both are descriptive only.","- This dataset concerns hypertrophic scar/burn wound biology, not keloid-specific disease replication.","",
"## Evidence","",paste0("TEMPORAL_ACTIVATION_SUPPORT=",e1),paste0("TEMPORAL_ATTENUATION_SUPPORT=",e2),paste0("PATHOLOGICAL_HTS_PERSISTENCE_SUPPORT=",e3),"","TWPS_CHANGED=NO")
writeLines(report,file.path(res,"GSE178411_TWPS_REPORT.md"))
cat("CORE100_COVERAGE=",coverage['CORE100'],"/100\n",sep="");cat("ANALYZED_PATIENT_N=",length(unique(sample_scores$patient_id)),"\n",sep="")
for(r in list(r1,r2,r3))cat(r$contrast," difference=",fmt(r$difference)," g=",fmt(r$effect_size)," P=",fmt(r$pvalue)," class=",r$evidence_class,"\n",sep="")
cat("CONTINUOUS=",continuous_text,"\n",sep="")
