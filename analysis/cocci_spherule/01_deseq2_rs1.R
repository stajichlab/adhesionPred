#!/usr/bin/env Rscript
# DESeq2 on the existing C. immitis RS kallisto counts (RS1 = FungiDB-46 / NCBI RefSeq
# GCF_000149335.2 annotation; CIMG_ gene IDs). Two replicates per condition, so the
# test has little power: read padj as a ranking aid, not as proof.
# Usage: Rscript 01_deseq2_rs1.R <counts.csv> <out_dir>
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
gene <- sub("-t26_1$", "", rownames(cts))
# kallisto est_counts are fractional. Sum transcripts per gene, then round for DESeq2.
gene_cts <- rowsum(as.matrix(cts), group = gene)
gene_cts <- round(gene_cts)
cond <- factor(sub("\\.r[0-9]+$", "", colnames(gene_cts)),
               levels = c("Mycelia", "Spherule48H", "Spherule8D"))
coldata <- data.frame(row.names = colnames(gene_cts), condition = cond)
dds <- DESeqDataSetFromMatrix(gene_cts, coldata, design = ~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds, quiet = TRUE)

res_list <- list()
for (s in c("Spherule48H", "Spherule8D")) {
  r <- results(dds, contrast = c("condition", s, "Mycelia"))
  d <- as.data.frame(r)
  d$gene <- rownames(d)
  names(d)[names(d) != "gene"] <- paste0(names(d)[names(d) != "gene"], ".", s)
  res_list[[s]] <- d
}
out <- merge(res_list[[1]], res_list[[2]], by = "gene", all = TRUE)
mean_norm <- as.data.frame(sapply(levels(cond), function(l)
  rowMeans(counts(dds, normalized = TRUE)[, cond == l, drop = FALSE])))
names(mean_norm) <- paste0("norm_count.", names(mean_norm))
mean_norm$gene <- rownames(mean_norm)
out <- merge(out, mean_norm, by = "gene", all.x = TRUE)
out <- out[order(-out$log2FoldChange.Spherule48H), ]
write.table(out, file.path(out_dir, "deseq2_rs1_gene.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("genes tested:", nrow(dds), "\n")
cat("DESeq2 version:", as.character(packageVersion("DESeq2")), "\n")
