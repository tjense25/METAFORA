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

parser <- arg_parser("Segment population beta profiles and calculate block beta values and global methylation PCs")
parser <- add_argument(parser, "--population_beta", help="population mean beta file")
parser <- add_argument(parser, "--beta_mat", help="input sample beta matrix")
parser <- add_argument(parser, "--depth_mat", help="input sample depth matrix")
parser <- add_argument(parser, "--segment_bed_out", help="output bed file of segmented cpgs")
parser <- add_argument(parser, "--varsegment_bed_out", help="output bed file of segmented cpgs")
parser <- add_argument(parser, "--varseg_plot_file", help="file to write variable segment CV2 plot")
parser <- add_argument(parser, "--pca_plot_dir", help="output global methylation PCs")
parser <- add_argument(parser, "--global_meth_pcs_out", help="output global methylation PCs")

argv <- parse_args(parser)

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
    segs <- as.data.frame(fastseg(block_beta$mean_beta, alpha=alpha, minSeg=minSeg, segMedianT=c(.65,.35)))
    segs$start <- block_beta$start[segs$start]
    segs$end <- block_beta$end[segs$end]
    segs$width <- segs$end - segs$start
    segs$seqnames <- chrom
    segs$seg_id <- paste0(chrom,"_",i,"_",1:nrow(segs)) 
    segs[segs$num.mark >= minSeg,]
    }))
  makeGRangesFromDataFrame(segments,keep.extra.columns = T)
}

segment_popbeta_chrom <- function(betas, depth, pop_mean, this_chrom) { 
  # filter to chromosome
  betas <- betas[betas$chrom == this_chrom,]
  depth <- depth[depth$chrom == this_chrom,]
  pop_mean <- pop_mean[pop_mean$chrom == this_chrom,]
  
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

findVariableSegments <- function(seg_betas, plot_file=NULL) {
   seg_betas <- seg_betas[seg_betas$cpg_density > .005,]
   beta.mat <- as.matrix(seg_betas[,8:ncol(seg_betas)])
   seg_cv2 <- data.frame(mean=rowMeans(beta.mat), var=apply(beta.mat,1,var))
   seg_cv2$CV2=seg_cv2$var/seg_cv2$mean
   cv2.fit <- fitTrendCV2(seg_cv2$mean, seg_cv2$CV2, min.mean = 0.2, ncells=10)
   var_segs <- seg_cv2$CV2 > cv2.fit$trend(seg_cv2$mean) + 3*cv2.fit$trend(seg_cv2$mean)*cv2.fit$std.dev/sqrt(ncol(beta.mat))
   if(!is.null(plot_file)) {
    gp <- ggplot(seg_cv2, aes(abs(mean), CV2,color=var_segs)) + geom_point() + scale_color_manual(values=c("black", "red"))  +
      geom_line(data=data.frame(x=seq(0.01,1,.01),y=cv2.fit$trend(seq(0.01,1,.01))), mapping=aes(x,y,color="blue"),color="blue",linewidth=2) + 
        scale_y_log10()
    ggsave(file=plot_file,plot=gp)
   }
   seg_betas[var_segs,]
}
 

findGlobalMethylationPCs <- function(var_seg_betas, plot_dir=NULL) {
     M <- logit(var_seg_betas[,8:ncol(var_seg_betas)]) #convert beta into M-values
     outlier_samples=NULL
     i <- 1
     while(T) {
      M.scaled <- t(scale(t(M), center=T, scale = T))
      pcs <- pca(M.scaled)
      ntest <- PCAtools::findElbowPoint(pcs$variance)
      outliers <- apply(pcs$rotated[,1:2], 2, function(x) (x < quantile(x,.25)-3*IQR(x)) | (x > quantile(x,.75)+3*IQR(x)))
      current_outliers <- unique(rownames(which(outliers,arr.ind=T)))
      outlier_samples <-  c(outlier_samples,current_outliers)

      if (!is.null(plot_dir) & plot_dir != "") {
          x_lim = max(abs(pcs$rotated[,1]))*1.35
          y_lim = max(abs(pcs$rotated[,2]))*1.35
          plot_title = paste(c("samples",current_outliers,"removed"),collapse=" ")
          if (length(current_outliers) == 0) {plot_title="No Global Outliers"}
         pca_plot <- PCAtools::biplot(pcs) + ggtitle(paste0("outlier detection: iteration ",i)) +
            stat_ellipse(level = .95, geom="polygon", fill="blue", alpha=.2) + 
            geom_text_repel(mapping=aes(label=rownames(pcs$rotated))) +
            ylim(-y_lim,y_lim) + xlim(-x_lim, x_lim) +
            labs(title=plot_title) +
            theme(plot.title=element_text(color="red"))
        ggsave(file=paste0(plot_dir, '/PCA_plot_', i, '.pdf'), plot=pca_plot)
      }
      i <- i+1
      if (any(outliers) == F) { break }
      M <- M[,-which(colnames(M) %in% outlier_samples)]
     }
     cat(paste0("Following samples are global methylation outliers and will be excluded: ",outlier_samples,"\n"))
     return(pcs$rotated[,1:ntest])
}

cat("Reading in data . . . \n")
pop_mean <- fread(argv$population_beta) 
betas <- fread(argv$beta_mat)
depths <- fread(argv$depth_mat)

cat("Segmenting population mean betas . . . \n")
valid_chrom = paste0("chr", seq(1,22))
segment_betas <- Reduce(rbind,pbmclapply(valid_chrom, function(x)  segment_popbeta_chrom(betas, depths, pop_mean, x)))
write.table(segment_betas, file=argv$segment_bed_out, quote=F,row.names=F,col.names=T,sep="\t")

cat("Calculating variable methylation segments . . . \n")
var_seg_betas <- findVariableSegments(segment_betas, plot_file=argv$varseg_plot_file)
write.table(segment_betas, file=argv$varsegment_bed_out, quote=F,row.names=F,col.names=T,sep="\t")
cat("Calculating global methylation PCs . . . \n")
global_meth_PCs <- findGlobalMethylationPCs(var_seg_betas, plot_dir=argv$pca_plot_dir)
write.table(global_meth_PCs, file=argv$global_meth_pcs_out, quote=F, row.names=T, col.names=T, sep="\t")
cat("Done! have a lovely day :)\n")
