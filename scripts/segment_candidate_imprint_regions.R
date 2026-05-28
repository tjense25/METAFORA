library(data.table)
library(plyranges)
library(tidyverse)
library(Matrix)
library(matrixStats)
library(fastseg)
library(argparser)

parser <- arg_parser("Script to calculate haplotype effect at CpGs and segment to find potential imprinting regions")
parser <- add_argument(parser, "--hp1_beta", help = "CpG-level matrix of methylation betas for this block from HP1")
parser <- add_argument(parser, "--hp1_cov", help = "CpG-level depth matrix for this block from HP1")
parser <- add_argument(parser, "--hp2_beta", help = "CpG-level matrix of methylation betas for this block from HP2")
parser <- add_argument(parser, "--hp2_cov", help = "CpG-level depth matrix for this block from HP2")
parser <- add_argument(parser, "--chrom", help = "chrom block to do segmentation over, should match files in other arguments")
parser <- add_argument(parser,"--imprint_out", help = "where to write candidate imprinting regions with sample-level haplotype deltas")

argv <- parse_args(parser)
#argv <- NULL
#argv$hp1_beta <- "../METAFORA_output/HP_1.tissue_PBMC/Population_methylation.hp_1.tissue_PBMC.chrom_chr20.1.betas.mat.gz"
#argv$hp1_cov <- "../METAFORA_output/HP_1.tissue_PBMC/Population_methylation.hp_1.tissue_PBMC.chrom_chr20.1.coverage.mat.gz"
#argv$hp2_beta <- "../METAFORA_output/HP_2.tissue_PBMC/Population_methylation.hp_2.tissue_PBMC.chrom_chr20.1.betas.mat.gz"
#argv$hp2_cov <- "../METAFORA_output/HP_2.tissue_PBMC/Population_methylation.hp_2.tissue_PBMC.chrom_chr20.1.coverage.mat.gz"
#argv$chrom <- "chr20.1"
#argv$imprint_out <- "../METAFORA_output/imprinting_regions.tissue_PBMC/Candidate_imprinting_loci.tissue_PBMC.chrom_chr20.1.abs_haplotype_delta.mat"

main <- function(argv) {
    chrom_block <- argv$chrom
    hp1_betas <- fread(argv$hp1_beta)
    hp2_betas <- fread(argv$hp2_beta)

    hp1_cov <- fread(argv$hp1_cov)
    hp2_cov <- fread(argv$hp2_cov)

    hp1.beta.mat <- as.matrix(hp1_betas[,4:ncol(hp1_betas)])
    hp2.beta.mat <- as.matrix(hp2_betas[,4:ncol(hp2_betas)])
    hp1.depth.mat <- as.matrix(hp1_cov[,4:ncol(hp1_cov)])
    hp2.depth.mat <- as.matrix(hp2_cov[,4:ncol(hp2_cov)])

    hp1.beta.mat[hp1.depth.mat<5] <- NA
    hp2.beta.mat[hp2.depth.mat<5] <- NA

    #laplacian smoothing of betas
    hp1.beta.mat <- (hp1.beta.mat*hp1.depth.mat+1)/(hp1.depth.mat+2)
    hp2.beta.mat <- (hp2.beta.mat*hp2.depth.mat+1)/(hp2.depth.mat+2)

    abs_hap_delt <- abs(hp1.beta.mat - hp2.beta.mat)
    pop_hap_delta <- data.table(chrom=hp1_betas$chromosome,start=hp1_betas$start,end=hp1_betas$end)
    pop_hap_delta$median_hap_delta <- rowMedians(abs_hap_delt,na.rm=T)

    pop_hap_delta <- pop_hap_delta[!is.na(pop_hap_delta$median_hap_delta),]

    segment_hap_delta <- function(pop_hap_delta, chrom_block, alpha = 0.01, minSeg = 10) {
      index <- c(1,which(diff(pop_hap_delta$start) > 1000))
      last_start <- index[length(index)]
      last_end <- nrow(pop_hap_delta)
      block_start <- index[-length(index)]
      block_end <- index[-1] - 1
      block_start <- c(block_start, last_start)
      block_end <- c(block_end, last_end)
      blocks = data.frame(start=block_start, end=block_end)
      
      segments <- Reduce(rbind,lapply(1:nrow(blocks), function(i) { 
        block_beta <- pop_hap_delta[blocks$start[i]:blocks$end[i],]
        segs <- as.data.frame(fastseg(block_beta$median_hap_delta, alpha=alpha, minSeg=minSeg, segMedianT=c(.5,0)))
        segs$start <- block_beta$start[segs$start]
        segs$end <- block_beta$end[segs$end]
        segs$width <- segs$end - segs$start
        segs$seqnames <- gsub(chrom_block, pattern="(\\w+)\\.\\d+",replacement="\\1")
        segs
        }))
      makeGRangesFromDataFrame(segments,keep.extra.columns = T)
    }

    segs <- segment_hap_delta(pop_hap_delta,chrom_block) 
    cand_imprint_segs <- segs %>% filter(seg.mean > .5)

    if(length(cand_imprint_segs)<1) {
        file.remove(argv$imprint_out)
        fwrite(data.table(NULL),file=argv$imprint_out)
        return()
    }

    cand_imprint_segs$seg_id <- paste0(seqnames(cand_imprint_segs),":",start(cand_imprint_segs),"-",end(cand_imprint_segs))
    cpg.gr <- hp1_betas %>% select(chromosome,start,end) %>% makeGRangesFromDataFrame
    ol <- findOverlaps(cand_imprint_segs, cpg.gr)
    CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(cand_imprint_segs),length(cpg.gr)), x=1)

    hp1.beta.mat <- as.matrix(hp1_betas[,4:ncol(hp1_betas)])
    hp2.beta.mat <- as.matrix(hp2_betas[,4:ncol(hp2_betas)])

    hp1.beta.mat[is.na(hp1.beta.mat)] <- 0
    hp2.beta.mat[is.na(hp2.beta.mat)] <- 0
    hp1.depth.mat[is.na(hp1.depth.mat)] <- 0
    hp2.depth.mat[is.na(hp2.depth.mat)] <- 0

    hp1_imprint_beta <- matrix(((CpG_Identity%*%(hp1.beta.mat*hp1.depth.mat))+1)/ ((CpG_Identity%*%hp1.depth.mat)+2), nrow=length(cand_imprint_segs), ncol=ncol(hp1.beta.mat))
    hp1_imprint_depth <- matrix((CpG_Identity%*%hp1.depth.mat)/rowSums(CpG_Identity), nrow=length(cand_imprint_segs), ncol=ncol(hp1.depth.mat))
    hp2_imprint_beta <- matrix(((CpG_Identity%*%(hp2.beta.mat*hp2.depth.mat))+1)/ ((CpG_Identity%*%hp2.depth.mat)+2), nrow=length(cand_imprint_segs), ncol=ncol(hp2.beta.mat))
    hp1_imprint_depth <- matrix((CpG_Identity%*%hp2.depth.mat)/rowSums(CpG_Identity), nrow=length(cand_imprint_segs), ncol=ncol(hp2.depth.mat))

    hp1_imprint_beta[hp1_imprint_depth<5] <- NA
    hp2_imprint_beta[hp1_imprint_depth<5] <- NA

    imprint_hap_delta <- abs(hp1_imprint_beta-hp2_imprint_beta)
    colnames(imprint_hap_delta) <- colnames(hp1.beta.mat)

    cand_imprint.hap_deltas <- data.table(cbind(cand_imprint_segs %>% as.data.frame %>% select(seqnames,start,end,seg_id), imprint_hap_delta))
    fwrite(cand_imprint.hap_deltas, file=argv$imprint_out, sep="\t", scipen=999)
}

main(argv)
