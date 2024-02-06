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
parser <- add_argument(parser, "--chrom", help="chromosome to call outliers")
parser <- add_argument(parser, "--sample", help="sample to parse and detect outliers in")
parser <- add_argument(parser, "--population_beta", help="population mean beta file")
parser <- add_argument(parser, "--beta_mat", help="input sample beta matrix")
parser <- add_argument(parser, "--depth_mat", help="input sample depth matrix")
parser <- add_argument(parser, "--global_meth_pcs", help="output global methylation PCs")
parser <- add_argument(parser, "--outlier_bed", help="output global methylation region bed with zscores to this file")
parser <- add_argument(parser, "--outlier_z_mat", help="output outlier_region x samples matrix of zscores to this file")
parser <- add_argument(parser, "--plot_dir", help="where to plot outlier regions", default="NULL")
 
argv <- parse_args(parser)

segment_chrom <- function(meth.sample, this_chrom, zscore_threshold=2, segment_alpha=.01, min_seg_size = 20, median_seg_z = 1) {
  # segment zscore profile
  index <- c(1,which(diff(meth.sample$start) > 10000))
  last_start <- index[length(index)]
  last_end <- nrow(meth.sample)
  block_start <- index[-length(index)]
  block_end <- index[-1] - 1
  block_start <- c(block_start, last_start)
  block_end <- c(block_end, last_end)
  
  cand.segs <-Reduce(rbind,lapply(1:length(block_start), function(i) {  
    tmp.meth <- meth.sample[block_start[i]:block_end[i],]
    as.data.frame(fastseg(tmp.meth$zscore, alpha = segment_alpha, 
            minSeg = min_seg_size, segMedianT = c(median_seg_z,-median_seg_z))) %>%
      mutate(start = tmp.meth$start[start], end=tmp.meth$start[end],
             seqnames=this_chrom)}))
  cand.segs <- cand.segs[abs(cand.segs$seg.mean) > zscore_threshold,]
  return(cand.segs)
}

segment_candidate_outliers <- function(pop_mean, betas, depth, this_sample, this_chrom, zscore_threshold=2, segment_alpha=.01, min_seg_size = 20, median_seg_z = 1) {
  betas.sample <- betas %>% select("chromosome","start",all_of(this_sample))
  colnames(betas.sample) <- c("chrom","start","sample_beta")
  
  depth.sample <- depth %>% select("chromosome","start",all_of(this_sample))
  colnames(depth.sample) <- c("chrom", "start", "sample_depth")
  
  meth.sample <- pop_mean %>% left_join(betas.sample) %>% left_join(depth.sample) %>% mutate(cpg_num = 1:n())
  meth.sample$sample_depth %<>% pmin(20) #cap depth at 20x coverage

  #beta correction on sample
  meth.sample$sample_beta <- (meth.sample$sample_beta*meth.sample$sample_depth+1)/(meth.sample$sample_depth+2)
  #beta correction on population mean
  N <- ncol(betas) - 3
  meth.sample$mean_beta <- (meth.sample$mean_beta*meth.sample$total_depth + N) / (meth.sample$total_depth + 2*N)
  
  # calculate zscore
  meth.sample %<>% mutate(zscore = qnorm(pbeta(q=sample_beta, mean_beta*sample_depth, (1-mean_beta)*sample_depth)))
  meth.sample$zscore %<>% pmax(-5)
  meth.sample$zscore %<>% pmin(5)

  meth.sample <- meth.sample[!is.na(meth.sample$sample_beta),]
    
  # segment zscore profile
  cand.segs <- segment_chrom(meth.sample, this_chrom, zscore_threshold, segment_alpha, min_seg_size, median_seg_z)
  if (nrow(cand.segs)==0) { 
    cand.segs <- NULL
  } else {
    cand.segs %<>% mutate(ID=this_sample, seg_id=paste0(this_chrom,"_",this_sample,"_",1:nrow(cand.segs)))
  }
  return(list("meth.sample"=meth.sample,"cand.segs"=cand.segs))
}

call_outliers <-function(cand.segs, betas, depth, sample_id, covariates=NULL) {
  betas.gr <- makeGRangesFromDataFrame(betas)
  cands.gr <- makeGRangesFromDataFrame(cand.segs, keep.extra.columns = T)
  
  betas.mat <-  as.matrix(betas[,4:ncol(betas)])
  depth.mat <- as.matrix(depth[,4:ncol(depth)])
  depth.mat[depth.mat > 20] <- 20
  
  ol <- findOverlaps(cands.gr, betas.gr)
  # create Segment x CpG identity matrix to indicate which cpgs belong to which segment
  CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(cands.gr),nrow(betas)), x=1)
  
  betas.mat[is.na(betas.mat)] <- 0
  depth.mat[is.na(depth.mat)] <- 0
  # aggregate betas across each segment
  region_beta <- matrix(((CpG_Identity %*% (betas.mat*depth.mat))+1) / ((CpG_Identity %*% depth.mat)+2), nrow=nrow(cand.segs), ncol=ncol(betas.mat))
  colnames(region_beta) <- colnames(betas.mat)
  cand.segs$delta <- region_beta[,sample_id] - rowMedians(region_beta, na.rm=T)
  M <- logit(region_beta)
  M.scaled <- t(scale(t(M),center=T,scale=T))
  #M.corrected <- removeBatchEffect(M.scaled, covariates=covariates)
  #zscores <- t(scale(t(M.corrected)))
  cand.segs$zscore <- M.scaled[,sample_id]

  M.scaled <- M.scaled[abs(cand.segs$zscore) > 2.5,]
  cand.segs <- cand.segs[abs(cand.segs$zscore) > 2.5,]
  if (nrow(cand.segs)==0) { return(list("outlier.segs"=NULL, "z.mat"=NULL)) }

  M.scaled <- matrix(M.scaled, nrow=nrow(cand.segs), ncol=ncol(betas.mat))
  colnames(M.scaled) <- colnames(betas.mat)
  rownames(M.scaled) <- cand.segs$seg_id
  return(list("outlier.segs"=cand.segs,"z.mat"=M.scaled))
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

     betas.plot <- tmp.meth.sample %>% 
       pivot_longer(cols=c("mean_beta.smoothed", "sample_beta.smoothed"), names_to="source", values_to="beta") %>% 
       ggplot(aes(start, beta, color=source)) + geom_point() + geom_line(alpha=.8) + theme_classic()  +
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
    this_chrom <- argv$chrom
    this_sample <- argv$sample
    pop_mean <- fread(argv$population_beta) 
    betas <- fread(argv$beta_mat)
    depths <- fread(argv$depth_mat)
    global_pcs <- read.table(argv$global_meth_pcs)
    # REMOVE GLOBAL OUTLIERS + RECOMPUTE POP MEAN
    #cols=colnames(betas)[colnames(betas) %in% c("chromosome", "start", "end", rownames(global_pcs))]
    #betas <- betas[,..cols]
    #depths <- depths[,..cols]

    #if(!(this_sample %in% cols)) {
    #    write.table(NULL, file=argv$outlier_bed, row.names=F, col.names=T, quote=F)
    #    write.table(NULL, file=argv$outlier_z_mat,row.names=T,col.names=T, quote=F)
    #    cat(paste0(this_sample, " is a global outlier. WILL NOT DETECT OUTLIERS.\n"))
    #    return()
    #}

    cat("Calculating population difference zscore and segment to find candidate outlier region . . .\n")
    # calculate pop difference zscore and segment to find candidate outlier regions
    cand.outliers <- segment_candidate_outliers(pop_mean, betas, depths,this_sample=this_sample, this_chrom=this_chrom)
    meth.sample <- cand.outliers[["meth.sample"]]
    cand.segs <- cand.outliers[["cand.segs"]]
    if(is.null(nrow(cand.segs))) {
        write.table(NULL, file=argv$outlier_bed, row.names=F, col.names=T, quote=F)
        write.table(NULL, file=argv$outlier_z_mat,row.names=T,col.names=T, quote=F)
        return()
    }

    # calculate region aggregated M values and call zscores across samples
    cat("Calculating region aggregated M values and zscores . . . \n")
    outliers <- call_outliers(cand.segs, betas, depths,sample_id=this_sample)
    outlier.segs <- outliers[["outlier.segs"]]
    outlier_z_matrix <- outliers[["z.mat"]]
    if(is.null(nrow(outlier.segs))) {
        write.table(NULL, file=argv$outlier_bed, row.names=F, col.names=T, quote=F)
        write.table(NULL, file=argv$outlier_z_mat,row.names=T,col.names=T, quote=F)
        return()
    }

    ## Plot Outlier regions
    cat("Plotting outliers . . . \n")
    if (!is.null(argv$plot_dir)) {
      plot_outliers(this_sample, outlier.segs, meth.sample, outlier_z_matrix, plot_dir=argv$plot_dir)
    }

    ## Save Data
    write.table(outlier.segs, file=argv$outlier_bed, row.names=F, col.names=T, quote=F)
    write.table(outlier_z_matrix, file=argv$outlier_z_mat,row.names=T,col.names=T, quote=F)
}
 
main(argv)
