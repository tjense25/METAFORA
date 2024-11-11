library(data.table)
library(limma)
library(ggrepel)
library(tibble)
library(tidyverse)
library(scran)
library(argparser)
library(sigmoid)
library(PCAtools)
library(Matrix)
library(magrittr)
library(RNOmni)

parser <- arg_parser("Segment population beta profiles and calculate block beta values and global methylation PCs")
parser <- add_argument(parser, "--segment_betas", help="mean betas of all chromosome segments found in segmentation step")
parser <- add_argument(parser, "--varsegment_bed_out", help="output bed file of variable segmented cpgs")
parser <- add_argument(parser, "--varseg_plot_file", help="file to write variable segment variance plot")
parser <- add_argument(parser, "--pca_plot_dir", help="output global methylation PCs")
parser <- add_argument(parser, "--global_meth_pcs_out", help="where to write global PCs")

argv <- parse_args(parser)

findVariableSegments <- function(seg_betas, plot_file=NULL) {
   seg_betas <- segment_betas
   seg_betas <- seg_betas[seg_betas$cpg_density > .005,]
   beta.mat <- as.matrix(seg_betas[,8:ncol(seg_betas)])

   transformed <- logit(.5 + abs(.5-beta.mat))
   seg_cv2 <- data.frame(mean=rowMeans(transformed), var=apply(transformed,1,var))

   seg_cv2$CV2=seg_cv2$var/seg_cv2$mean

   cv2.fit <- fitTrendCV2(seg_cv2$mean, seg_cv2$CV2,min.mean=.1,ncells=10000)
   var_segs <- seg_cv2$CV2 > cv2.fit$trend(seg_cv2$mean) + 3*cv2.fit$trend(seg_cv2$mean)*cv2.fit$std.dev/sqrt(ncol(beta.mat))

   if(!is.null(plot_file)) {
    gp <- ggplot(seg_cv2, aes(abs(mean), CV2,color=var_segs)) + geom_point() + scale_color_manual(values=c("black", "red"))  +
      geom_line(data=data.frame(x=seq(0.01,max(seg_cv2$mean),.01),y=cv2.fit$trend(seq(0.01,max(seg_cv2$mean),.01))), mapping=aes(x,y,color="blue"),color="blue",linewidth=2)
    ggsave(file=plot_file,plot=gp)
   }
   seg_betas[var_segs,]
}
 

findGlobalMethylationPCs <- function(var_seg_betas, plot_dir=NULL) {
     M <- logit(as.matrix(var_seg_betas[,8:ncol(var_seg_betas)])) #convert beta into M-values
     outlier_samples=NULL
     i <- 1
     while(T) {
      M.scaled <- t(scale(t(M), center=T, scale = T))
      pcs <- pca(M.scaled)
      ntest <- PCAtools::findElbowPoint(pcs$variance)
      outliers <- apply(pcs$rotated[,1:2], 2, function(x) (x < quantile(x,.25)-3*IQR(x)) | (x > quantile(x,.75)+3*IQR(x)))
      current_outliers <- unique(rownames(which(outliers,arr.ind=T)))
      outlier_samples <-  c(outlier_samples,current_outliers)
      outlier_samples

      if (!is.null(plot_dir) & plot_dir != "") {
          x_lim = max(abs(pcs$rotated[,1]))*1.35
          y_lim = max(abs(pcs$rotated[,2]))*1.35
          plot_title = paste(c("samples",current_outliers,"removed"),collapse=" ")
          if (length(current_outliers) == 0) {plot_title="No Global Outliers"}
         pca_plot <- PCAtools::biplot(pcs,max.overlaps=0) + ggtitle(paste0("outlier detection: iteration ",i)) +
            stat_ellipse(level = .95, geom="polygon", fill="blue", alpha=.2) + 
            geom_text_repel(mapping=aes(label=rownames(pcs$rotated))) +
            labs(title=plot_title) +
            theme(plot.title=element_text(color="red"))
        ggsave(file=paste0(plot_dir, '/PCA_plot_', i, '.pdf'), plot=pca_plot, width=15,height=15)
        ggsave(file=paste0('PCA_plot_ranknorm.pdf'), plot=pca_plot, width=15,height=15)
      }
      i <- i+1
      if (any(outliers) == F) { break }
      M <- M[,-which(colnames(M) %in% outlier_samples)]
     }
     cat(paste0("Following samples are global methylation outliers and will be excluded: ",outlier_samples,"\n"))
     return(pcs$rotated[,1:20])
}

cat("reading in data . . .")
segment_betas <- fread(argv$segment_betas)
segment_betas <- fread("/oak/stanford/groups/smontgom/tannerj/ADRC_LRS/output/methylation_results/Meth_segments.GRCh38.tissue_PBMC.segment_betas.ALL_AUTOSOMES.bed")

cat("Calculating variable methylation segments . . . \n")
var_seg_betas <- findVariableSegments(segment_betas, plot_file=argv$varseg_plot_file)
fwrite(var_seg_betas, file=argv$varsegment_bed_out, quote=F,row.names=F,col.names=T,sep="\t")

cat("Calculating global methylation PCs . . . \n")
global_meth_PCs <- findGlobalMethylationPCs(var_seg_betas, plot_dir=argv$pca_plot_dir)
fwrite(global_meth_PCs, file=argv$global_meth_pcs_out, quote=F, row.names=T, col.names=T, sep="\t")
