library(data.table)
library(limma)
library(tibble)
library(fastseg)
library(GenomicRanges)
library(tidyverse)
library(argparser)
library(PCAtools)
library(pbmcapply)
library(Matrix)
library(magrittr)
library(matrixStats)
library(plyranges)
options(dplyr.summarise.inform = FALSE, dplyr.join.inform = FALSE)

parser <- arg_parser("Calculate and segment per sample zscore")
parser <- add_argument(parser, "--tissue", help="tissue this sample was sequenced from")
parser <- add_argument(parser, "--chrom", help="chromosome to call outliers on")
parser <- add_argument(parser, "--beta_mat", help="input sample beta matrix")
parser <- add_argument(parser, "--depth_mat", help="input sample depth matrix")
parser <- add_argument(parser, "--global_meth_pcs", help="output global methylation PCs")
parser <- add_argument(parser, "--outlier_bed", help="output global methylation region bed with zscores to this file")
parser <- add_argument(parser, "--outlier_z_mat", help="output outlier_region x samples matrix of zscores to this file")
parser <- add_argument(parser, "--min_abs_zscore", help="minimum absolute zscore threshold for calling methylation region an outlier", type="integer", default=3)
parser <- add_argument(parser, "--min_seg_size", help="minimum number of cpgs in a region to call a candidate outlier during segmentation", type="integer", default=20)
parser <- add_argument(parser, "--min_abs_delta", help="minimum effect size delta (sample_average_methylation - population_median_methylation) for calling methylation region an outlier", type="numeric", default=0.25)
parser <- add_argument(parser, "--max_depth", help="maxmimum depth of read coverage to consider. (regions higher than this depth will be effectively downsampled)", type="integer", default=30)
parser <- add_argument(parser, "--threads", help="number of threads to use for paralellized chrom block segmentation", default=1)
 
argv <- parse_args(parser)
ncores=as.integer(argv$threads)

MIN_Z_THRESH <- function(D) {1.172304 + .0355*D} #optimized parameters accounting for depth from simulation experiment to acheive 90% power for absolute deltas of 0.25
MAX_Z_THRESH <- function(D) {5.5 + 0.16654*D} #optimized  parameters accounting for depth for zscores observed for 0.9 deltas 

segment_candidate_outliers <- function(pop_mean, betas, depth, this_sample, this_chrom, segment_alpha=.01, min_seg_size = 10, MAX_DEPTH=100) {
  betas.sample <- betas %>% select("chromosome","start",all_of(this_sample))
  colnames(betas.sample) <- c("chromosome","start","sample_beta")
  
  depth.sample <- depth %>% select("chromosome","start",all_of(this_sample))
  colnames(depth.sample) <- c("chromosome", "start", "sample_depth")
  
  meth.sample <- pop_mean %>% left_join(betas.sample,join_by(chromosome,start)) %>% left_join(depth.sample,join_by(chromosome,start)) %>% mutate(cpg_num = 1:n())
  meth.sample$sample_depth %<>% pmin(MAX_DEPTH) #cap depth at specified value to not inflate p-values for high-coverage samples

  D=median(depth.sample$sample_depth,na.rm=T)
  MIN_Z=max(1,MIN_Z_THRESH(D))
  MAX_Z=MAX_Z_THRESH(D)

  #beta correction on sample
  meth.sample$sample_beta <- (meth.sample$sample_beta*meth.sample$sample_depth+1)/(meth.sample$sample_depth+2)
  
  # calculate zscore
  meth.sample %<>% mutate(deviance_score = qnorm(pbeta(q=sample_beta, mean_beta*sample_depth, (1-mean_beta)*sample_depth)))
  meth.sample$deviance_score %<>% pmax(-MAX_Z)
  meth.sample$deviance_score %<>% pmin(MAX_Z)
  meth.sample <- meth.sample[!is.na(meth.sample$sample_beta),]

  # segment zscore profile
  index <- c(1,which(diff(meth.sample$start) > 1000)) #break into blocks of contiguous cpgs in at least 1000bp window
  last_start <- index[length(index)]
  last_end <- nrow(meth.sample) #number of cpgs on chrom
  block_start <- index[-length(index)]
  block_end <- index[-1] - 1
  block_start <- c(block_start, last_start)
  block_end <- c(block_end, last_end)

  blocks <- data.frame(start = block_start, end = block_end)
  blocks <- blocks[blocks$end - blocks$start >= min_seg_size,] #do not segment blocks with less than minimum seg_size contiguous cpgs

  cand.segs <- pbmclapply(1:nrow(blocks), function(i) {  
            tmp.meth <- meth.sample[blocks$start[i]:blocks$end[i],]
            as.data.frame(fastseg(tmp.meth$deviance_score, alpha=segment_alpha, minSeg=10, segMedianT = c(MIN_Z,-MIN_Z))) %>%
                  mutate(start = tmp.meth$start[start], end=tmp.meth$start[end],width=end-start, seqnames=this_chrom,ID=this_sample)
              #seg_id=paste0(this_chrom,"_",this_sample,"_",1:nrow(segs)))
            }, mc.cores=ncores) %>% bind_rows()
  if(is.null(cand.segs)||nrow(cand.segs)==0) { 
        return(list("meth.sample"=meth.sample,"cand.segs"=NULL))
  }
  cand.segs$seg_id=paste0(this_chrom,"_",this_sample,"_",1:nrow(cand.segs))

  cand.segs <- cand.segs[!is.na(cand.segs$`seg.mean`),] #remove cand segments that are NA over all cpgs
  cand.segs <- cand.segs[abs(cand.segs$seg.mean) > MIN_Z,]
  cand.segs <- cand.segs[cand.segs$`num.mark` >= min_seg_size,]
  return(list("meth.sample"=meth.sample,"cand.segs"=cand.segs))
}

call_joint_outliers <- function(outliers.gr, cpgs.gr, beta.mat, depth.mat, covariates=NULL) {
  ol <- findOverlaps(outliers.gr, cpgs.gr)
  # create Segment x CpG identity matrix to indicate which cpgs belong to which segment
  CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(outliers.gr),nrow(beta.mat)), x=1)
  
  # aggregate betas across each segment
  region_beta <- as.matrix(((CpG_Identity %*% (beta.mat*depth.mat))+1) / ((CpG_Identity %*% depth.mat)+2))

  M <- log(region_beta/(1-region_beta))
  Batch=NULL
  if(!is.null(covariates) && "Batch" %in% colnames(covariates)) {
      Batch=covariates$Batch
      covariates <- covariates %>% select(-Batch)

      #if Batch only has 1 group or any Batch group has less than 4 individuals 
      #(not big enough sample size to reliable correct batch effects) set batch to NULL
      if((length(table(Batch)) < 2) || any(table(Batch) <= 5)) {
          Batch=NULL
      }
  }
  M.corrected <- removeBatchEffect(M, batch=Batch, covariates=covariates)
  zscores <- t(scale(t(M.corrected)))
  return(zscores)
}

call_outliers <-function(cand.segs, cpgs.gr, beta.mat, depth.mat, sample_id, MIN_ABS_ZSCORE = 3, covariates=NULL) {
  cands.gr <- makeGRangesFromDataFrame(cand.segs, keep.extra.columns = T)
  
  ol <- findOverlaps(cands.gr, cpgs.gr)
  # create Segment x CpG identity matrix to indicate which cpgs belong to which segment
  CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(cands.gr),nrow(beta.mat)), x=1)
  
  # aggregate betas across each segment
  region_beta <- as.matrix(((CpG_Identity %*% (beta.mat*depth.mat))+1) / ((CpG_Identity %*% depth.mat)+2))
  cand.segs$pop_median <- rowMedians(region_beta,na.rm=T)

  colnames(region_beta) <- colnames(beta.mat)
  if (!is.null(covariates) && "Batch" %in% colnames(covariates)) {
      this_batch <- covariates$Batch[rownames(covariates) == sample_id]

      #create batch specific profile only if batch group has more than 10 samples
      if (sum(covariates$Batch == this_batch) > 10){
          cand.segs$pop_median <- rowMedians(matrix(region_beta[,covariates$Batch==this_batch],nrow=dim(region_beta)[1]),na.rm=T)
      }
  }
  cand.segs$delta <- region_beta[,sample_id] - cand.segs$pop_median

  M <- log(region_beta/(1-region_beta))
  #M.scaled <- t(scale(t(M),center=T,scale=T))
  Batch=NULL
  if(!is.null(covariates) && "Batch" %in% colnames(covariates)) {
      Batch=covariates$Batch
      covariates <- covariates %>% select(-Batch)

      #if Batch only has 1 group or any Batch group has less than 4 individuals 
      #(not big enough sample size to reliable correct batch effects) set batch to NULL
      if((length(table(Batch)) < 2) || any(table(Batch) <= 5)) {
          Batch=NULL
      }
  }
  M.corrected <- removeBatchEffect(M, batch=Batch, covariates=covariates)
  zscores <- t(scale(t(M.corrected)))
  cand.segs$zscore <- zscores[,sample_id]

  zscores <- zscores[abs(cand.segs$zscore) > MIN_ABS_ZSCORE,]
  cand.segs <- cand.segs[abs(cand.segs$zscore) > MIN_ABS_ZSCORE,]
  if (nrow(cand.segs)==0) { return(list("outlier.segs"=NULL, "z.mat"=NULL)) }

  zscores <- matrix(zscores, nrow=nrow(cand.segs), ncol=ncol(beta.mat))
  colnames(zscores) <- colnames(beta.mat)
  rownames(zscores) <- cand.segs$seg_id
  return(list("outlier.segs"=cand.segs,"z.mat"=zscores))
}

outlier_pipeline <- function(pop_mean, betas, depths, cpgs.gr, beta.mat, depth.mat, covariates, this_sample,this_chrom, min_seg_size, sample_id, MAX_DEPTH=100, MIN_ABS_ZSCORE=3, MIN_ABS_DELTA=0.25) {
        cat(paste0("Segmenting and calling outliers for sample: ",this_sample, " on chrom ", this_chrom, " . . . \n"))
        cand.outliers <- segment_candidate_outliers(pop_mean, betas, depths, this_sample=this_sample, this_chrom=this_chrom, min_seg_size=min_seg_size, MAX_DEPTH=MAX_DEPTH)
        meth.sample = cand.outliers[["meth.sample"]]
        cand.segs <- cand.outliers[["cand.segs"]]
        if(nrow(cand.segs)==0||is.null(nrow(cand.segs))) { return(NULL) }
        # calculate region aggregated M values and call zscores across samples
        outliers <- call_outliers(cand.segs, cpgs.gr, beta.mat, depth.mat, sample_id=this_sample, MIN_ABS_ZSCORE=MIN_ABS_ZSCORE, covariates=covariates)
        outlier.segs <- outliers[["outlier.segs"]]
        outlier_z_matrix <- outliers[["z.mat"]]
        outlier.segs
        if(nrow(outlier.segs)==0||is.null(nrow(outlier.segs))) { return(NULL) }
        keep <- (abs(outlier.segs$delta) > MIN_ABS_DELTA)
        outlier_z_matrix <- outlier_z_matrix[keep,]
        outlier.segs <- outlier.segs[keep,]
        if(nrow(outlier.segs)==0||is.null(nrow(outlier.segs))) { return(NULL) }
        return(outlier.segs)
}

main <- function(argv) {
    # Read in data
    cat("Reading in data . . . \n")
    this_chrom <- argv$chrom
    betas <- fread(argv$beta_mat)
    depths <- fread(argv$depth_mat)

    cpgs.gr <- makeGRangesFromDataFrame(betas[,1:3])

    MIN_ABS_ZSCORE <- as.numeric(argv$min_abs_zscore)
    MIN_SEG_SIZE <- as.integer(argv$min_seg_size)
    MIN_ABS_DELTA <- as.numeric(argv$min_abs_delta)
    MAX_DEPTH <- as.numeric(argv$max_depth)
    covariates <- NULL
    if(!is.null(argv$global_meth_pcs)) {
        covariates <- read.table(argv$global_meth_pcs, row.names=1, header=T) 
        covariates
        colnames(betas)
        # REMOVE GLOBAL OUTLIERS + RECOMPUTE POP MEAN
        cols=colnames(betas)[colnames(betas) %in% c("chromosome", "start", "end", rownames(covariates))]
        betas <- betas[,..cols]
        depths <- depths[,..cols]
    }


 
    ## Correct pop mean for specific batch
    beta.mat <- as.matrix(betas[,4:ncol(betas)])
    depth.mat <- as.matrix(depths[,4:ncol(depths)])
    beta.mat[is.na(beta.mat)] <- 0
    depth.mat[is.na(depth.mat)] <- 0
    depth.mat[depth.mat > MAX_DEPTH] <- MAX_DEPTH
    unique_batches <- unique(covariates$Batch)
    if(is.null(unique_batches)) {
        unique_batches <- list(NULL)
    }

    combined_outliers <- NULL
    for (this_batch in unique_batches) {
        cat(paste0("Calling outliers for samples in Batch ", this_batch, " . . . \n"))
        cat("Correcting population mean for sample batch . . .\n")
        tmp.bmat <- beta.mat
        tmp.dmat <- depth.mat
        #make population mean specific to Batch to make less susceptible to batch effects:
        if (sum(covariates$Batch==this_batch) > 20) {
            tmp.bmat <- tmp.bmat[,covariates$Batch == this_batch]
            tmp.dmat <- tmp.dmat[,covariates$Batch == this_batch]
        }
        pop_mean <- (rowSums(tmp.bmat*tmp.dmat) + rowSums(tmp.dmat > 0)) / (rowSums(tmp.dmat) + 2*rowSums(tmp.dmat > 0))
        total_depth <- rowSums(tmp.dmat, na.rm=T)
        pop_mean <- data.table(betas[,1:3], total_depth, mean_beta=pop_mean)
        batch_samples <- colnames(tmp.bmat)
        if (!is.null(this_batch)) {
           batch_samples <- rownames(covariates[covariates$Batch == this_batch,])
        }
        rm(tmp.bmat)
        rm(tmp.dmat)
        gc()

        outliers <- lapply(batch_samples, function(x) {
                   outlier_pipeline(pop_mean,betas,depths,cpgs.gr,beta.mat,depth.mat,covariates,
                                    this_sample=x,this_chrom=this_chrom,min_seg_size=MIN_SEG_SIZE,
                                    MAX_DEPTH=MAX_DEPTH,MIN_ABS_ZSCORE=MIN_ABS_ZSCORE,MIN_ABS_DELTA=MIN_ABS_DELTA) 
                }) %>% bind_rows
        cat("\n")
        combined_outliers <- rbind(combined_outliers, outliers)
    }

    cat("Merging sample-level outliers and computing joint cohort z-scores . . . \n")
    combined_outliers$Tissue <- argv$tissue
    combined_outliers$CHROM_TYPE <- "AUTOSOME"

    outliers.gr <- makeGRangesFromDataFrame(combined_outliers)
    outliers.merged <- outliers.gr %>% GenomicRanges::reduce()
    outliers.merged$ID <- paste0("METAFORA_",1:length(outliers.merged),"_",this_chrom)
    outliers.gr %<>% join_overlap_left(outliers.merged) 
    combined_outliers$MERGE_ID <- outliers.gr$ID

    combined_mat <- call_joint_outliers(outliers.merged, cpgs.gr, beta.mat, depth.mat, covariates=covariates)
    outlier_z_mat <- data.frame(chromosome=seqnames(outliers.merged), start=start(outliers.merged), end=end(outliers.merged), MERGE_ID=outliers.merged$ID, combined_mat)
    ## Save Data
    write.table(combined_outliers, file=argv$outlier_bed, row.names=F, col.names=T, quote=F)
    write.table(outlier_z_mat, file=argv$outlier_z_mat,row.names=F,col.names=T, quote=F)
}
 
main(argv)
