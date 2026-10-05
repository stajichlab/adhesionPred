#!/usr/bin/env Rscript
# Audit of REPORT_2026-10-01_sowgp_depth_vs_repeats.md: confounder adjustment and relatedness.
#
# Q1 (status vs whole-gene depth ratio). Response log(ratio). Strains: tree tips with genome
#    depth >= 10x (as script 41). Contrasts: no-model - full-length, fragment - full-length,
#    no-model - fragment. Reported as fold difference in the ratio, exp(coef), with 95% CI.
# Q2 (depth vs unit count). Full-length strains, genome >= 10x (as scripts 40 and 42), one
#    model per species. Responses: array/flank (all reads), array/flank (MAPQ >= 20), whole-gene
#    ratio. Slope per unit with 95% CI; tests of slope = 0 and, for array/flank, slope = 0.25.
#
# Models, each fitted to the same rows:
#   ols_raw      response ~ predictor (+ species for pooled Q1)
#   ols_adj      + log(genome depth) + sequencing center + read-length class + log(insert size)
#   pgls_cds     ols_raw terms, GLS with a phylogenetic covariance (CDS tree), Pagel lambda by ML
#   pgls_cds_adj ols_adj terms, same covariance
#   pgls_snp_adj ols_adj terms, covariance from the genome-wide SNP tree (fewer strains)
# Covariance: Brownian motion on each species' own subtree (ape::vcv), scaled to a maximum
# diagonal of 1, block-diagonal across species. The species split is a fixed effect, not part
# of the covariance, because the CDS tree's species stem (~0.05) is about 500x its median
# within-species branch (~1e-4) and would otherwise swamp within-species relatedness.
# Pagel lambda multiplies the off-diagonal elements; lambda = 0 is OLS.
#
# Inputs : CIMG_04613.coverage.tsv, sowgp_tree_units.tsv, sowgp_array_depth.tsv,
#          cram_readlen.tsv (script 64), Genotyping samples.csv (Center), the two trees
# Outputs: audit_adjusted_models.tsv, stdout
# Run    : Rscript 61_depth_adjusted_models.R | tee audit_adjusted_models.log
suppressMessages({
  library(ape)
  library(phytools)
  library(sandwich)
})

P <- "/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci"
CDS_TREE <- file.path(P, "Phylogeny/results/msa_filter_cds_ascomycota-buildtree/Cocci_cds.488taxa_ascomycota.fa.part.aicc.contree")
SNP_TREE <- file.path(P, "Genotyping/all_C_immitis_ref_RS/strain_tree/unpruned/Cocciall_v1.All.SNP.mfa.treefile")
SAMPLES <- file.path(P, "Genotyping/all_C_immitis_ref_RS/samples.csv")

# ------------------------------------------------------------------ data
cov <- read.delim("CIMG_04613.coverage.tsv")
tu <- read.delim("sowgp_tree_units.tsv")
tu$strain <- sub("^Coccidioides_(immitis|posadasii)_", "", tu$tip)
arr <- read.delim("sowgp_array_depth.tsv")
rl <- read.delim("cram_readlen.tsv")
smp <- read.csv(SAMPLES, colClasses = "character")
smp$center <- ifelse(smp$Center == "", "unrecorded", smp$Center)

d <- merge(cov, tu[, c("strain", "tip", "status", "clade_species", "n_units")], by = "strain")
d <- merge(d, rl[, c("strain", "modal_len", "median_insert")], by = "strain")
d <- merge(d, smp[, c("Strain", "center")], by.x = "strain", by.y = "Strain")
d <- merge(d, arr[, c("strain", "array_flank", "array_flank_q20", "flank")], by = "strain")
d$rl_class <- cut(d$modal_len, c(0, 101, 126, 151, 251, 301), labels = c("le101", "126", "150", "251", "301"))
d$log_ins <- log(d$median_insert)
d$log_gen <- log(d$genome_mean)
d$species <- factor(d$clade_species)
d$status <- factor(d$status, levels = c("full-length", "fragment only", "no SOWgp gene model"))
cat("joined strains:", nrow(d), "\n")
cat("covariate table, tree strains with genome >= 10x:\n")
q1 <- subset(d, genome_mean >= 10)
print(table(q1$status, q1$center))
print(table(q1$status, q1$rl_class))

# ------------------------------------------------------------------ trees -> covariance
species_blocks <- function(tree_file, d, tipcol) {
  tr <- read.tree(tree_file)
  tr <- midpoint.root(tr)
  tr$edge.length[tr$edge.length < 0] <- 0
  out <- list()
  for (sp in levels(d$species)) {
    tips <- d[[tipcol]][d$species == sp]
    tips <- tips[tips %in% tr$tip.label]
    st <- keep.tip(tr, tips)
    C <- vcv(st)
    out[[sp]] <- C / max(diag(C))
  }
  out
}
blockC <- function(blocks, labels, species) {
  n <- length(labels)
  C <- matrix(0, n, n, dimnames = list(labels, labels))
  for (sp in names(blocks)) {
    i <- which(species == sp)
    C[i, i] <- blocks[[sp]][labels[i], labels[i]]
  }
  C
}
cds_blocks <- species_blocks(CDS_TREE, d, "tip")
snp_blocks <- species_blocks(SNP_TREE, d, "strain")

# ------------------------------------------------------------------ GLS with Pagel lambda
# grid stops at 0.999: tips with zero-length terminal branches make C singular at lambda = 1
gls_lambda <- function(y, X, C, grid = c(seq(0, 0.99, by = 0.01), 0.999)) {
  n <- length(y); p <- ncol(X)
  fit1 <- function(lam) {
    V <- C * lam; diag(V) <- diag(C) * (1 + 1e-9)
    L <- t(chol(V))
    yw <- forwardsolve(L, y); Xw <- forwardsolve(L, X)
    qrw <- qr(Xw)
    b <- qr.coef(qrw, yw); r <- yw - Xw %*% b
    s2ml <- sum(r^2) / n
    ll <- -n / 2 * log(2 * pi * s2ml) - sum(log(diag(L))) - n / 2
    list(lam = lam, b = b, Xw = Xw, r = r, ll = ll)
  }
  fits <- lapply(grid, fit1)
  ll <- sapply(fits, `[[`, "ll")
  best <- fits[[which.max(ll)]]
  s2 <- sum(best$r^2) / (n - p)
  V <- s2 * solve(crossprod(best$Xw))
  list(coef = setNames(drop(best$b), colnames(X)), vcov = V, df = n - p, lambda = best$lam,
       lrt_p = pchisq(2 * (max(ll) - fits[[1]]$ll), 1, lower.tail = FALSE), n = n)
}
ols <- function(y, X) {
  f <- lm.fit(X, y); n <- length(y); p <- ncol(X)
  s2 <- sum(f$residuals^2) / (n - p)
  V <- s2 * chol2inv(qr.R(f$qr))[order(f$qr$pivot), order(f$qr$pivot)]
  dimnames(V) <- list(colnames(X), colnames(X))
  list(coef = setNames(f$coefficients, colnames(X)), vcov = V, df = n - p, lambda = NA, lrt_p = NA, n = n)
}
contrast <- function(fit, w, null = 0) {
  w <- w[names(fit$coef)]; w[is.na(w)] <- 0
  est <- sum(w * fit$coef); se <- sqrt(drop(t(w) %*% fit$vcov %*% w))
  tq <- qt(0.975, fit$df)
  c(est = est, lo = est - tq * se, hi = est + tq * se, p = 2 * pt(-abs((est - null) / se), fit$df))
}
design <- function(form, data) {
  X <- model.matrix(form, data)
  X[, colSums(abs(X)) > 0 & !duplicated(t(X)), drop = FALSE]
}

RES <- list()
add <- function(...) RES[[length(RES) + 1]] <<- data.frame(...)

# ------------------------------------------------------------------ Q1
cat("\n======== Q1 status vs log(ratio) ========\n")
f_raw <- ~ status + species
f_adj <- ~ status + species + log_gen + center + rl_class + log_ins
contr <- list(
  "no-model - full-length" = c("statusno SOWgp gene model" = 1),
  "fragment - full-length" = c("statusfragment only" = 1),
  "no-model - fragment" = c("statusno SOWgp gene model" = 1, "statusfragment only" = -1)
)
run_q1 <- function(label, data, form, C = NULL) {
  X <- design(form, data); y <- log(data$ratio)
  fit <- if (is.null(C)) ols(y, X) else gls_lambda(y, X, C)
  for (k in names(contr)) {
    cc <- contrast(fit, contr[[k]])
    add(question = "Q1", response = "log(ratio)", species = "both", model = label, term = k,
        estimate = exp(cc["est"]), ci_low = exp(cc["lo"]), ci_high = exp(cc["hi"]), p = cc["p"],
        p_slope_0.25 = NA, lambda = fit$lambda, lambda_lrt_p = fit$lrt_p, n = fit$n,
        scale = "fold difference in ratio")
    cat(sprintf("%-14s %-24s fold %.3f (95%% CI %.3f-%.3f) p=%.3g  lambda=%s n=%d\n", label, k,
                exp(cc["est"]), exp(cc["lo"]), exp(cc["hi"]), cc["p"],
                ifelse(is.na(fit$lambda), "-", sprintf("%.2f (LRT p=%.2g)", fit$lambda, fit$lrt_p)), fit$n))
  }
  if ("log_gen" %in% names(fit$coef)) {
    cc <- contrast(fit, c(log_gen = 1))
    cat(sprintf("%-14s %-24s coef %.3f (95%% CI %.3f-%.3f) p=%.3g\n", label, "log genome depth", cc["est"], cc["lo"], cc["hi"], cc["p"]))
  }
  invisible(fit)
}
q1 <- q1[order(q1$species, q1$tip), ]
run_q1("ols_status_only", q1, ~ status)
run_q1("ols_raw", q1, f_raw)
run_q1("ols_depth", q1, ~ status + species + log_gen)
run_q1("ols_center", q1, ~ status + species + center)
run_q1("ols_readlen", q1, ~ status + species + rl_class)
run_q1("ols_insert", q1, ~ status + species + log_ins)
run_q1("ols_adj", q1, f_adj)
cat("median log(ratio) by read-length class and status:\n")
print(round(tapply(q1$ratio, list(q1$rl_class, q1$status), median), 3))
cat("full-length strains only: ratio by read-length class:\n")
fls <- subset(q1, status == "full-length")
print(round(tapply(fls$ratio, fls$rl_class, median), 3))
print(kruskal.test(ratio ~ rl_class, data = fls))
# HC3 robust SE for the adjusted OLS
m <- lm(update(f_adj, log(ratio) ~ .), data = q1)
V <- vcovHC(m, type = "HC3")
for (k in names(contr)) {
  w <- setNames(rep(0, length(coef(m))), names(coef(m))); w[names(contr[[k]])] <- contr[[k]]
  est <- sum(w * coef(m)); se <- sqrt(drop(t(w) %*% V %*% w)); tq <- qt(0.975, m$df.residual)
  p <- 2 * pt(-abs(est / se), m$df.residual)
  add(question = "Q1", response = "log(ratio)", species = "both", model = "ols_adj_HC3", term = k,
      estimate = exp(est), ci_low = exp(est - tq * se), ci_high = exp(est + tq * se), p = p,
      p_slope_0.25 = NA, lambda = NA, lambda_lrt_p = NA, n = nobs(m), scale = "fold difference in ratio")
  cat(sprintf("%-14s %-24s fold %.3f (95%% CI %.3f-%.3f) p=%.3g\n", "ols_adj_HC3", k, exp(est),
              exp(est - tq * se), exp(est + tq * se), p))
}
Ccds <- blockC(cds_blocks, q1$tip, q1$species)
run_q1("pgls_cds", q1, f_raw, Ccds)
run_q1("pgls_cds_adj", q1, f_adj, Ccds)
q1s <- q1[q1$strain %in% unlist(lapply(snp_blocks, rownames)), ]
Csnp <- blockC(snp_blocks, q1s$strain, q1s$species)
run_q1("ols_adj_snpset", q1s, f_adj)
run_q1("pgls_snp_adj", q1s, f_adj, Csnp)
# within-species versions (no pooling)
for (sp in levels(q1$species)) {
  s <- droplevels(q1[q1$species == sp, ])
  run_q1(paste0("ols_raw_", sp), s, ~ status)
  run_q1(paste0("pgls_cds_adj_", sp), s, ~ status + log_gen + center + rl_class + log_ins,
         blockC(cds_blocks[sp], s$tip, s$species))
}

# ------------------------------------------------------------------ Q2
cat("\n======== Q2 depth vs unit count (full-length, genome >= 10x) ========\n")
q2 <- subset(d, status == "full-length" & genome_mean >= 10 & flank > 0)
q2 <- q2[order(q2$species, q2$tip), ]
q2$u <- q2$n_units
for (sp in levels(q2$species)) {
  s <- droplevels(q2[q2$species == sp, ])
  ss <- s[s$strain %in% rownames(snp_blocks[[sp]]), ]
  for (resp in c("array_flank", "array_flank_q20", "ratio")) {
    y <- s[[resp]]
    for (mdl in c("ols_raw", "ols_adj", "pgls_cds", "pgls_cds_adj", "pgls_snp_adj")) {
      form <- if (grepl("adj", mdl)) ~ u + log_gen + center + rl_class + log_ins else ~ u
      dat <- if (mdl == "pgls_snp_adj") ss else s
      X <- design(form, dat); yy <- dat[[resp]]
      fit <- switch(mdl,
        ols_raw = ols(yy, X), ols_adj = ols(yy, X),
        pgls_cds = gls_lambda(yy, X, blockC(cds_blocks[sp], dat$tip, dat$species)),
        pgls_cds_adj = gls_lambda(yy, X, blockC(cds_blocks[sp], dat$tip, dat$species)),
        pgls_snp_adj = gls_lambda(yy, X, blockC(snp_blocks[sp], dat$strain, dat$species)))
      cc <- contrast(fit, c(u = 1))
      p25 <- if (resp == "array_flank") contrast(fit, c(u = 1), null = 0.25)["p"] else NA
      add(question = "Q2", response = resp, species = sp, model = mdl, term = "slope per unit",
          estimate = cc["est"], ci_low = cc["lo"], ci_high = cc["hi"], p = cc["p"], p_slope_0.25 = p25,
          lambda = fit$lambda, lambda_lrt_p = fit$lrt_p, n = fit$n, scale = "response units per repeat unit")
      cat(sprintf("%-9s %-15s %-13s slope %+.3f (95%% CI %+.3f to %+.3f) p0=%.3g p0.25=%s lambda=%s n=%d\n",
                  sp, resp, mdl, cc["est"], cc["lo"], cc["hi"], cc["p"],
                  ifelse(is.na(p25), "-", sprintf("%.3g", p25)),
                  ifelse(is.na(fit$lambda), "-", sprintf("%.2f(LRT p=%.2g)", fit$lambda, fit$lrt_p)), fit$n))
    }
  }
}

R <- do.call(rbind, RES)
rownames(R) <- NULL
# BH within each model, over the tests that match the original families
R$family <- with(R, ifelse(question == "Q1", "F1",
                    ifelse(response == "array_flank" & !is.na(p_slope_0.25), "F2", "F2")))
R$q_bh_model_family <- NA
for (key in unique(paste(R$question, R$model))) {
  i <- which(paste(R$question, R$model) == key)
  R$q_bh_model_family[i] <- p.adjust(R$p[i], "BH")
}
i <- which(!is.na(R$p_slope_0.25))
R$q_bh_slope_0.25 <- NA
for (mdl in unique(R$model[i])) {
  j <- i[R$model[i] == mdl]
  R$q_bh_slope_0.25[j] <- p.adjust(R$p_slope_0.25[j], "BH")
}
write.table(R, "audit_adjusted_models.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nwrote audit_adjusted_models.tsv:", nrow(R), "rows\n")
