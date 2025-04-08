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
parser <- add_argument(parser, "--outlier_bed", help="outlier bed output from combine outlier rule")
parser <- add_argument(parser, "--combined", help="Metafora formatted combined methylation bed")
parser <- add_argument(parser, "--hap1", help="Metafora formatted Haplotype 1 bed")
parser <- add_argument(parser, "--hap2", help="Metafora formatted Haplotype 2 bed")
parser <- add_argument(parser, "--annotated_out", help="where to write the haplotype annotated bed file")
argv <- parse_args(parser)

if(file.info(argv$outlier_bed)$size == 0) {
    fwrite(NULL, argv$annotated_out)
    quit(save = "no")
}

outliers <- fread(argv$outlier_bed)
combined <- fread(argv$combined)
hap1 <- fread(argv$hap1)
hap2 <- fread(argv$hap2)

outliers.gr <- makeGRangesFromDataFrame(outliers, keep.extra.columns =T)
hap1.gr <- makeGRangesFromDataFrame(hap1, keep.extra.columns =T)
hap2.gr <- makeGRangesFromDataFrame(hap2, keep.extra.columns=T)
combined.gr <- makeGRangesFromDataFrame(combined, keep.extra.columns=T)

ol <- findOverlaps(outliers.gr, combined.gr)
CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(outliers.gr),nrow(combined)), x=1)
combined$depth[is.na(combined$depth)] <- 0
outliers$combined_depth <- as.integer((CpG_Identity %*% combined$depth) / rowSums(CpG_Identity))

ol <- findOverlaps(outliers.gr, hap1.gr)
# create Segment x CpG identity matrix to indicate which cpgs belong to which segment
CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(outliers.gr),nrow(hap1)), x=1)
hap1$depth[is.na(hap1$depth)] <- 0
hap1$beta[is.na(hap1$beta)] <- 0
hap1$beta <- ((hap1$beta*hap1$depth) + 1) / (hap1$depth + 2)
outliers$hap1_depth <- as.integer((CpG_Identity %*% hap1$depth) / rowSums(CpG_Identity))
hap1$depth <- hap1$depth + 2
# aggregate betas across each segment
outliers$hap1_beta <- as.matrix(((CpG_Identity %*% (hap1$beta * hap1$depth))) / ((CpG_Identity %*% hap1$depth)))

ol <- findOverlaps(outliers.gr, hap2.gr)
CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(outliers.gr),nrow(hap2)), x=1)
hap2$depth[is.na(hap2$depth)] <- 0
hap2$beta[is.na(hap2$beta)] <- 0
hap2$beta <- ((hap2$beta*hap2$depth) + 1) / (hap2$depth + 2)
outliers$hap2_depth <- as.integer((CpG_Identity %*% hap2$depth) / rowSums(CpG_Identity))
hap2$depth <- hap2$depth + 2
outliers$hap2_beta <- as.matrix(((CpG_Identity %*% (hap2$beta * hap2$depth))) / ((CpG_Identity %*% hap2$depth)))

outliers$phasing_percent <- pmin(outliers$combined_depth / (outliers$hap1_depth + outliers$hap2_depth), 1)
outliers$haplotype_coverage_bias <- log(outliers$hap1_depth/outliers$hap2_depth)
outliers$hap_delta <- outliers$hap1_beta - outliers$hap2_beta

fwrite(outliers, argv$annotated_out, sep="\t")
