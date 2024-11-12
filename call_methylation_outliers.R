library(data.table)
library(limma)
library(ggrepel)
library(tibble)
library(fastseg)
library(GenomicRanges)
library(tidyverse)
library(scran)
library(argparser)
library(sigmoid)
library(PCAtools)
library(pbmcapply)
library(Matrix)
library(magrittr)
library(cowplot)
library(matrixStats)

parser <- arg_parser("Calculate and segment per sample zscore")
parser <- add_argument(parser, "--sample", help="sample to parse and detect outliers in")
parser <- add_argument(parser, "--chrom", help="chromosome to call outliers on")
parser <- add_argument(parser, "--beta_mat", help="input sample beta matrix")
parser <- add_argument(parser, "--depth_mat", help="input sample depth matrix")
parser <- add_argument(parser, "--global_meth_pcs", help="output global methylation PCs")
parser <- add_argument(parser, "--outlier_bed", help="output global methylation region bed with zscores to this file")
parser <- add_argument(parser, "--outlier_z_mat", help="output outlier_region x samples matrix of zscores to this file")
parser <- add_argument(parser, "--plot_dir", help="where to plot outlier regions", default=NULL)
parser <- add_argument(parser, "--min_abs_zscore", help="minimum absolute zscore threshold for calling methylation region an outlier", type="integer", default=3)
parser <- add_argument(parser, "--min_seg_size", help="minimum number of cpgs in a region to call a candidate outlier during segmentation", type="integer", default=20)
parser <- add_argument(parser, "--min_abs_delta", help="minimum effect size delta (sample_average_methylation - population_median_methylation) for calling methylation region an outlier", type="numeric", default=0.25)
parser <- add_argument(parser, "--max_depth", help="maxmimum depth of read coverage to consider. (regions higher than this depth will be effectively downsampled)", type="integer", default=30)
 
argv <- parse_args(parser)

segment_chrom <- function(meth.sample, this_chrom, zscore_threshold=2, segment_alpha=.01, min_seg_size = 20, median_seg_z = 1) {
  # segment zscore profile
  meth.sample <- meth.sample[meth.sample$chrom == this_chrom,]
  
  index <- c(1,which(diff(meth.sample$start) > 1000)) #break into blocks of contiguous cpgs in at least 1000bp window
  last_start <- index[length(index)]
  last_end <- nrow(meth.sample) #number of cpgs on chrom
  block_start <- index[-length(index)]
  block_end <- index[-1] - 1
  block_start <- c(block_start, last_start)
  block_end <- c(block_end, last_end)

  blocks <- data.frame(start = block_start, end = block_end)
  blocks <- blocks[blocks$end - blocks$start >= min_seg_size,] #do not segment blocks with less than minimum seg_size contiguous cpgs

  cand.segs <-Reduce(rbind,pbmclapply(1:nrow(blocks), function(i) {  
    tmp.meth <- meth.sample[blocks$start[i]:blocks$end[i],]
    as.data.frame(fastseg(tmp.meth$zscore, alpha = segment_alpha, 
            minSeg = min_seg_size, segMedianT = c(median_seg_z,-median_seg_z))) %>%
      mutate(start = tmp.meth$start[start], end=tmp.meth$start[end],
             seqnames=this_chrom)}))
  return(cand.segs)
}

segment_candidate_outliers <- function(pop_mean, betas, depth, this_sample, this_chrom, zscore_threshold=2, segment_alpha=.01, min_seg_size = 20, median_seg_z = 1, MAX_DEPTH=30) {
  betas.sample <- betas %>% select("chromosome","start",all_of(this_sample))
  colnames(betas.sample) <- c("chromosome","start","sample_beta")
  
  depth.sample <- depth %>% select("chromosome","start",all_of(this_sample))
  colnames(depth.sample) <- c("chromosome", "start", "sample_depth")
  
  meth.sample <- pop_mean %>% left_join(betas.sample) %>% left_join(depth.sample) %>% mutate(cpg_num = 1:n())
  meth.sample$sample_depth %<>% pmin(MAX_DEPTH) #cap depth at specified value to not inflate p-values for high-coverage samples

  #beta correction on sample
  meth.sample$sample_beta <- (meth.sample$sample_beta*meth.sample$sample_depth+1)/(meth.sample$sample_depth+2)
  
  # calculate zscore
  meth.sample %<>% mutate(zscore = qnorm(pbeta(q=sample_beta, mean_beta*sample_depth, (1-mean_beta)*sample_depth)))
  meth.sample$zscore %<>% pmax(-6)
  meth.sample$zscore %<>% pmin(6)
  meth.sample <- meth.sample[!is.na(meth.sample$sample_beta),]

  # segment zscore profile
  cand.segs <- Reduce(rbind,lapply(c(this_chrom), function(x) {  
        segs <- segment_chrom(meth.sample, x, zscore_threshold, segment_alpha, min_seg_size, median_seg_z)
        if (nrow(segs) == 0) {return(NULL)}
        return(segs %>% mutate(ID=this_sample, seg_id=paste0(x,"_",this_sample,"_",1:nrow(segs))))
      }))
  cand.segs <- cand.segs[abs(cand.segs$seg.mean) > zscore_threshold,]
  return(list("meth.sample"=meth.sample,"cand.segs"=cand.segs))
}

call_outliers <-function(cand.segs, betas, depth, sample_id, MIN_ABS_ZSCORE = 3, covariates=NULL) {
  betas.gr <- makeGRangesFromDataFrame(betas)
  cands.gr <- makeGRangesFromDataFrame(cand.segs, keep.extra.columns = T)
  
  betas.mat <-  as.matrix(betas[,4:ncol(betas)])
  depth.mat <- as.matrix(depth[,4:ncol(depth)])
  
  ol <- findOverlaps(cands.gr, betas.gr)
  # create Segment x CpG identity matrix to indicate which cpgs belong to which segment
  CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(cands.gr),nrow(betas)), x=1)
  
  betas.mat[is.na(betas.mat)] <- 0
  depth.mat[is.na(depth.mat)] <- 0
  # aggregate betas across each segment
  region_beta <- as.matrix(((CpG_Identity %*% (betas.mat*depth.mat))+1) / ((CpG_Identity %*% depth.mat)+2))
  colnames(region_beta) <- colnames(betas.mat)
  if (!is.null(covariates) && "Batch" %in% colnames(covariates)) {
      this_batch <- covariates$Batch[rownames(covariates) == sample_id]
      cand.segs$pop_median <- rowMedians(matrix(region_beta[,covariates$Batch==this_batch],nrow=dim(region_beta)[1]),na.rm=T)
  }
  else {
      cand.segs$pop_median <- rowMedians(region_beta,na.rm=T)
  }
  cand.segs$delta <- region_beta[,sample_id] - cand.segs$pop_median

  B <- region_beta
  M <- log(B/(1-B))
  #M.scaled <- t(scale(t(M),center=T,scale=T))
  Batch=NULL
  if(!is.null(covariates) && "Batch" %in% colnames(covariates)) {
      Batch=covariates$Batch
      covariates <- covariates %>% select(-Batch)
  }
  M.corrected <- removeBatchEffect(M, batch=Batch, covariates=covariates)
  zscores <- t(scale(t(M.corrected)))
  cand.segs$zscore <- zscores[,sample_id]

  zscores <- zscores[abs(cand.segs$zscore) > MIN_ABS_ZSCORE,]
  cand.segs <- cand.segs[abs(cand.segs$zscore) > MIN_ABS_ZSCORE,]
  if (nrow(cand.segs)==0) { return(list("outlier.segs"=NULL, "z.mat"=NULL)) }

  zscores <- matrix(zscores, nrow=nrow(cand.segs), ncol=ncol(betas.mat))
  colnames(zscores) <- colnames(betas.mat)
  rownames(zscores) <- cand.segs$seg_id
  return(list("outlier.segs"=cand.segs,"z.mat"=zscores))
}

plot_outliers <- function(samp, cand.segs, meth.sample, z.mat, plot_dir) {
  expanded.segs <- cand.segs
  size=pmax(1000,(cand.segs$end - cand.segs$start))
  expanded.segs$start=pmax(1,expanded.segs$start - size)
  expanded.segs$end = expanded.segs$end + size

  betas.gr <- makeGRangesFromDataFrame(meth.sample)
  segs.gr <- makeGRangesFromDataFrame(expanded.segs, keep.extra.columns = T)

  ol <- findOverlaps(betas.gr,segs.gr)
  meth.sample <- meth.sample[queryHits(ol)]
  meth.sample$seg_id <- segs.gr$seg_id[subjectHits(ol)]

  i=1
  for (i in 1:nrow(cand.segs)) {
     seg = cand.segs[i,] 
     seg_id = paste0(seg$seqnames,"_",seg$start,"_",seg$end)
     this_seg = seg$seg_id

     #subset to current segement
     tmp.meth.sample <- meth.sample %>% filter(seg_id == this_seg)
     #smooth betas and zscore for plotting
     tmp.meth.sample$sample_beta.smoothed <- sapply(1:nrow(tmp.meth.sample), function(x) mean(tmp.meth.sample$sample_beta[max(1,x-5):min(x+5,nrow(tmp.meth.sample))], na.rm=T))
     tmp.meth.sample$mean_beta.smoothed <- sapply(1:nrow(tmp.meth.sample), function(x) mean(tmp.meth.sample$mean_beta[max(1,x-5):min(x+5,nrow(tmp.meth.sample))], na.rm=T))
     tmp.meth.sample$zscore.smoothed <- sapply(1:nrow(tmp.meth.sample), function(x) mean(tmp.meth.sample$zscore[max(1,x-5):min(x+5,nrow(tmp.meth.sample))],na.rm=T))
     tmp.meth.sample$sd.smoothed <- sapply(1:nrow(tmp.meth.sample), function(x) mean(tmp.meth.sample$sd_beta[max(1,x-5):min(x+5,nrow(tmp.meth.sample))], na.rm=T))

     betas.plot <- tmp.meth.sample %>% 
       pivot_longer(cols=c("mean_beta.smoothed", "sample_beta.smoothed"), names_to="source", values_to="beta") %>% 
       mutate(beta_lower=beta-ifelse(source=="mean_beta.smoothed",yes=sd.smoothed, 0), 
              beta_higher=beta+ifelse(source=="mean_beta.smoothed",yes=sd.smoothed,0))  %>% 
       ggplot(aes(start, beta, ymin=beta_lower, ymax=beta_higher, color=source)) + geom_point() + geom_line(alpha=.8) + geom_ribbon(alpha=.1) + 
           theme_classic() +
           scale_color_manual(values=c("black","red"))  + 
           ylab("Methylation beta") +
           ggtitle(samp,subtitle=seg_id) +
           theme(legend.position = "none")

     zscore.plot <- tmp.meth.sample %>% filter(seg_id == this_seg) %>%
       ggplot(aes(start, zscore.smoothed)) + geom_point() + geom_line() + theme_classic() +
       geom_segment(data=as.data.frame(cand.segs[i,]),mapping=aes(x=start,xend=end,y=seg.mean,yend=seg.mean),color="red", linewidth=2) +
       ylab("population mean difference zscore") 

     z.df <- data.frame(zscore=as.numeric(z.mat[this_seg,]),
                       sample_id=colnames(z.mat),
                       outlier=(colnames(z.mat)==samp))
     hist <- ggplot(z.df, aes(zscore, fill=outlier)) + geom_histogram(bins=50) +  theme_classic() +
         geom_vline(xintercept=seg$zscore, color="red", linetype="dashed") +
         scale_fill_manual(values=c("grey40", "firebrick")) +
         theme(legend.position="none")

     plot_grid(betas.plot,zscore.plot,hist,ncol=1)
     ggsave(file=paste0(plot_dir, "/",seg_id,'.methylation_outlier_plots.pdf'), width = 4, height=8)
    }
}

main <- function(argv) {
    # Read in data
    cat("Reading in data . . . \n")
    this_sample <- argv$sample
    this_chrom <- argv$chrom
    betas <- fread(argv$beta_mat)
    depths <- fread(argv$depth_mat)

    MIN_ABS_ZSCORE <- argv$min_abs_zscore
    MIN_SEG_SIZE <- argv$min_seg_size
    MIN_ABS_DELTA <- argv$min_abs_delta
    MAX_DEPTH <- argv$max_depth
    covariates <- NULL
    if(!is.null(argv$global_meth_pcs)) {
        covariates <- read.table(argv$global_meth_pcs, row.names=1, header=T) 
        # REMOVE GLOBAL OUTLIERS + RECOMPUTE POP MEAN
        cols=colnames(betas)[colnames(betas) %in% c("chromosome", "start", "end", rownames(covariates))]
        betas <- betas[,..cols]
        depths <- depths[,..cols]
    }

    if(!(this_sample %in% cols)) {
        write.table(NULL, file=argv$outlier_bed, row.names=F, col.names=T, quote=F)
        write.table(NULL, file=argv$outlier_z_mat,row.names=T,col.names=T, quote=F)
        cat(paste0(this_sample, " is a global outlier. WILL NOT DETECT OUTLIERS.\n"))
        return()
    }

    ## Correct pop mean for specific batch
    cat("Correcting population mean for sample batch . . .\n")
    beta.mat <- as.matrix(betas[,4:ncol(betas)])
    depth.mat <- as.matrix(depths[,4:ncol(depths)])
    tmp.bmat <- beta.mat
    tmp.dmat <- depth.mat
    tmp.bmat[is.na(tmp.bmat)] <- 0
    tmp.dmat[is.na(tmp.dmat)] <- 0
    tmp.dmat[tmp.dmat > MAX_DEPTH] <- MAX_DEPTH

    #make population mean specific to Batch to make less susceptible to batch effects:
    if(!is.null(argv$global_meth_pcs) && "Batch" %in% colnames(covariates))  {
        this_batch <- covariates$Batch[rownames(covariates) == this_sample]
        this_batch
        tmp.bmat <- tmp.bmat[,covariates$Batch == this_batch]
        tmp.dmat <- tmp.dmat[,covariates$Batch == this_batch]
    }
    pop_mean <- (rowSums(tmp.bmat*tmp.dmat) + rowSums(tmp.dmat > 0)) / (rowSums(tmp.dmat) + 2*rowSums(tmp.dmat > 0))
    pop_sd <- rowSds(tmp.bmat, na.rm=T)
    total_depth <- rowSums(tmp.dmat, na.rm=T)

    pop_mean <- data.table(betas[,1:3], total_depth, mean_beta=pop_mean, sd_beta=pop_sd)
    cat("Calculating population difference zscore and segment to find candidate outlier region . . .\n")
    # calculate pop difference zscore and segment to find candidate outlier regions

    cand.outliers <- segment_candidate_outliers(pop_mean, betas, depths,this_sample=this_sample, this_chrom=this_chrom, min_seg_size=MIN_SEG_SIZE)
    meth.sample = cand.outliers[["meth.sample"]]
    cand.segs <- cand.outliers[["cand.segs"]]
    cand.segs
    if(nrow(cand.segs)==0||is.null(nrow(cand.segs))) {
        write.table(NULL, file=argv$outlier_bed, row.names=F, col.names=T, quote=F)
        write.table(NULL, file=argv$outlier_z_mat,row.names=T,col.names=T, quote=F)
        return()
    }
    # calculate region aggregated M values and call zscores across samples
    cat("Calculating region aggregated M values and zscores . . . \n")
    outliers <- call_outliers(cand.segs, betas, depths, sample_id=this_sample, MIN_ABS_ZSCORE=MIN_ABS_ZSCORE, covariates=covariates)
    outlier.segs <- outliers[["outlier.segs"]]
    outlier_z_matrix <- outliers[["z.mat"]]
    outlier.segs
    if(nrow(outlier.segs)==0||is.null(nrow(outlier.segs))) {
        write.table(NULL, file=argv$outlier_bed, row.names=F, col.names=T, quote=F)
        write.table(NULL, file=argv$outlier_z_mat,row.names=T,col.names=T, quote=F)
        return()
    }
    keep <- (abs(outlier.segs$zscore) > MIN_ABS_ZSCORE & abs(outlier.segs$delta) > MIN_ABS_DELTA)
    outlier_z_matrix <- outlier_z_matrix[keep,]
    outlier.segs <- outlier.segs[keep,]
    outlier.segs
    if(nrow(outlier.segs)==0||is.null(nrow(outlier.segs))) {
        write.table(NULL, file=argv$outlier_bed, row.names=F, col.names=T, quote=F)
        write.table(NULL, file=argv$outlier_z_mat,row.names=T,col.names=T, quote=F)
        return()
    }

    ## Save Data
    outlier_z_matrix = matrix(outlier_z_matrix, nrow=nrow(outlier.segs))
    rownames(outlier_z_matrix)  <- outlier.segs$seg_id
    colnames(outlier_z_matrix) <- colnames(betas)[4:ncol(betas)]
    write.table(outlier.segs, file=argv$outlier_bed, row.names=F, col.names=T, quote=F)
    write.table(outlier_z_matrix, file=argv$outlier_z_mat,row.names=T,col.names=T, quote=F)

    ## Plot Outlier regions
    cat("Plotting outliers . . . \n")
    if (!is.null(argv$plot_dir)) {
      plot_outliers(this_sample, outlier.segs, meth.sample, outlier_z_matrix, plot_dir=argv$plot_dir)
    }
}
 
main(argv)
