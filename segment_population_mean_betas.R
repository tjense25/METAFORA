library(data.table)
library(tibble)
library(fastseg)
library(GenomicRanges)
library(tidyverse)
library(argparser)
library(Matrix)
library(magrittr)
library(purrr)

parser <- arg_parser("Segment population beta profiles and calculate block beta values and global methylation PCs")
parser <- add_argument(parser, "--chrom", help="chromosome to segment")
parser <- add_argument(parser, "--population_beta", help="population mean beta file")
parser <- add_argument(parser, "--beta_mat", help="input sample beta matrix")
parser <- add_argument(parser, "--depth_mat", help="input sample depth matrix")
parser <- add_argument(parser, "--segment_bed_out", help="output bed file of segmented cpgs")

argv <- parse_args(parser)

breakup_large_segments <- function(segs, threshold=200, smaller_seg_size=100) {
    smaller_segments <- Reduce(rbind, lapply(1:nrow(segs), function(i) { 
            if (segs$num.mark[i] < threshold) {return(segs[i,])}
            N = round(segs$num.mark[i] / smaller_seg_size)
            stepsize = round(segs$num.mark[i] / N)
            starts <- seq(segs$start[i],segs$end[i], stepsize)[1:N]
            ends <- c(starts[2:length(starts)] - 1, segs$end[i])
            tmp <- map_dfr(seq_len(N), ~segs[i,])
            tmp$start <- starts
            tmp$end <- ends
            tmp$num.mark <- tmp$end - tmp$start + 1
            return(tmp)
    } ))
    return(smaller_segments)
}

segment_blocks <- function(pop_mean, chrom, alpha = 0.01, minSeg = 10) {
  index <- c(1,which(diff(pop_mean$start) > 10000))
  last_start <- index[length(index)]
  last_end <- nrow(pop_mean)
  block_start <- index[-length(index)]
  block_end <- index[-1] - 1
  block_start <- c(block_start, last_start)
  block_end <- c(block_end, last_end)
  
  segments <-Reduce(rbind,lapply(1:length(block_start), function(i) { 
    block_beta <- pop_mean[block_start[i]:block_end[i],]
    segs <- as.data.frame(fastseg(block_beta$mean_beta, alpha=alpha, minSeg=minSeg, segMedianT=c(.6,.4)))
    segs <- breakup_large_segments(segs)
    segs$start <- block_beta$start[segs$start]
    segs$end <- block_beta$end[segs$end]
    segs$width <- segs$end - segs$start
    segs$seqnames <- chrom
    segs$seg_id <- paste0(chrom,"_",i,"_",1:nrow(segs)) 
    segs
    }))
  makeGRangesFromDataFrame(segments,keep.extra.columns = T)
}

segment_popbeta <- function(betas, depth, pop_mean, this_chrom) { 
  betas.gr <- makeGRangesFromDataFrame(betas)
  betas.mat <-  as.matrix(betas[,4:ncol(betas)])
  depth.mat <- as.matrix(depth[,4:ncol(depth)])
  # cap depth at 30 so all samples have similar influence on mean (higher coverage samples will not influence more)
  depth.mat[depth.mat > 30] <- 30  
  
  #pop_mean = rowSums(betas.mat*depth.mat, na.rm=T) / rowSums(depth.mat, na.rm=T)
  meth_segs <- segment_blocks(pop_mean, this_chrom)
  
  ol <- findOverlaps(meth_segs, betas.gr)
  # create Segment x CpG identity matrix to indicate which cpgs belong to which segment
  CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(meth_segs),nrow(betas)), x=1)
  
  betas.mat[is.na(betas.mat)] <- 0
  depth.mat[is.na(depth.mat)] <- 0
  
  # aggregate betas across each segment
  seg_beta <- as.matrix(((CpG_Identity %*% (betas.mat*depth.mat))+1) / ((CpG_Identity %*% depth.mat)+2))
  
  data.frame(chrom=this_chrom, start=start(meth_segs), end=end(meth_segs), seg_id=meth_segs$seg_id, num_cpg = meth_segs$num.mark, width=width(meth_segs), 
             cpg_density = meth_segs$num.mark/width(meth_segs), as.data.frame(seg_beta))
}

cat("Reading in data . . . \n")
this_chrom <- argv$chrom
pop_mean <- fread(argv$population_beta) 
betas <- fread(argv$beta_mat)
depths <- fread(argv$depth_mat)
cat("Segmenting population mean betas . . . \n")
segment_betas <- segment_popbeta(betas, depths, pop_mean, this_chrom)
fwrite(segment_betas, file=argv$segment_bed_out, quote=F,row.names=F,col.names=T,sep="\t")

cat("Done! have a lovely day :)\n")
