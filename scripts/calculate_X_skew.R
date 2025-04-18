library(argparser)
library(Matrix)
library(tidyverse)
library(magrittr)
library(data.table)
library(plyranges)
library(GenomicRanges)
library(tibble)
library(bedr)

parser <- arg_parser("Annotate outlier regions with haplotype methylation deltas")
parser <- add_argument(parser, "--sample", help="sample name")
parser <- add_argument(parser, "--X_region_bed", help="regions on X chromosome to test skew")
parser <- add_argument(parser, "--combined", help="Metafora formatted combined methylation bed")
parser <- add_argument(parser, "--hap1", help="Metafora formatted Haplotype 1 bed")
parser <- add_argument(parser, "--hap2", help="Metafora formatted Haplotype 2 bed")
parser <- add_argument(parser, "--out_tsv", help="where to write the haplotype annotated bed file")
argv <- parse_args(parser)

xbed <- fread(argv$X_region_bed)
colnames(xbed) <- c("chrom", "start", "end", "gene")
xbed$sample <- argv$sample
combined <- fread(argv$combined)
hap1 <- fread(argv$hap1)
hap2 <- fread(argv$hap2)

xbed.gr <- makeGRangesFromDataFrame(xbed, keep.extra.columns =T)
hap1.gr <- makeGRangesFromDataFrame(hap1, keep.extra.columns =T)
hap2.gr <- makeGRangesFromDataFrame(hap2, keep.extra.columns=T)
combined.gr <- makeGRangesFromDataFrame(combined, keep.extra.columns=T)

ol <- findOverlaps(xbed.gr, combined.gr)
CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(xbed.gr),nrow(combined)), x=1)
combined$depth[is.na(combined$depth)] <- 0
xbed$combined_depth <- as.integer((CpG_Identity %*% combined$depth) / rowSums(CpG_Identity))

ol <- findOverlaps(xbed.gr, hap1.gr)
# create Segment x CpG identity matrix to indicate which cpgs belong to which segment
CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(xbed.gr),nrow(hap1)), x=1)
hap1$depth[is.na(hap1$depth)] <- 0
hap1$beta[is.na(hap1$beta)] <- 0
hap1$beta <- ((hap1$beta*hap1$depth) + 1) / (hap1$depth + 2)
xbed$hap1_depth <- as.integer((CpG_Identity %*% hap1$depth) / rowSums(CpG_Identity))
hap1$depth <- hap1$depth + 2
# aggregate betas across each segment
xbed$hap1_beta <- as.matrix(((CpG_Identity %*% (hap1$beta * hap1$depth))) / ((CpG_Identity %*% hap1$depth)))
xbed$hap1_deviance_score <-  qnorm(pbeta(q=xbed$hap1_beta, 0.5*xbed$hap1_depth, 0.5*xbed$hap1_depth))

ol <- findOverlaps(xbed.gr, hap2.gr)
CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(xbed.gr),nrow(hap2)), x=1)
hap2$depth[is.na(hap2$depth)] <- 0
hap2$beta[is.na(hap2$beta)] <- 0
hap2$beta <- ((hap2$beta*hap2$depth) + 1) / (hap2$depth + 2)
xbed$hap2_depth <- as.integer((CpG_Identity %*% hap2$depth) / rowSums(CpG_Identity))
hap2$depth <- hap2$depth + 2
xbed$hap2_beta <- as.matrix(((CpG_Identity %*% (hap2$beta * hap2$depth))) / ((CpG_Identity %*% hap2$depth)))
xbed$hap2_deviance_score <-  qnorm(pbeta(q=xbed$hap2_beta, 0.5*xbed$hap2_depth, 0.5*xbed$hap2_depth))

xbed$phasing_percent <- pmin(xbed$combined_depth / (xbed$hap1_depth + xbed$hap2_depth), 1)
xbed$haplotype_coverage_bias <- log(xbed$hap1_depth/xbed$hap2_depth)
xbed$hap_delta <- xbed$hap1_beta - xbed$hap2_beta

fwrite(xbed, argv$out_tsv, sep="\t")
