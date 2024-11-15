library(data.table)
library(cowplot)
library(ggrepel)
library(limma)
library(edgeR)
library(tibble)
library(GenomicRanges)
library(tidyverse)
library(argparser)
library(PCAtools)
library(pbmcapply)
library(Matrix)
library(magrittr)
library(mclust)
library(mvtnorm)
library(matrixStats)

parser <- arg_parser("Script to compute global variation PCs in methylation profiles and auto-detect outliers")
parser <- add_argument(parser, "--seg_beta", help = "comma-separated list of summarized segment betas from segmentation for each autosome")
parser <- add_argument(parser, "--seg_depth", help = "comma-separated list of summarized segment depths from segmentation for each autosome")
parser <- add_argument(parser, "--sex_seg_beta", help = "comma-separated list of summarized segment betas from segmentation for each sex chromosome")
parser <- add_argument(parser, "--sex_seg_depth", help = "comma-separated list of summarized segment depths from segmentation for each sex chromosome")
parser <- add_argument(parser, "--covariates", help="a dataframe of additional covariates to regress for during outlier detection", default=NULL)
parser <- add_argument(parser, "--correlation_plot_out", help="where to write correlation plots for determining correlation outliers")
parser <- add_argument(parser, "--pc_biplot_out", help="where to write PCA biplot")
parser <- add_argument(parser, "--sex_plot_out", help="where to write sex chromosome copy number and PCA biplot")
parser <- add_argument(parser, "--sex_chrom_summary_out", help="where to write data frame of sex chromosome copy number estimates")
parser <- add_argument(parser, "--global_meth_pcs_out", help="where to write global PC covariates file ")
parser <- add_argument(parser, "--chrX_seqname", help="seqname of the X chromosome in the reference")
parser <- add_argument(parser, "--chrY_seqname", help="seqname of the Y chromosome in the reference")

args <- parse_args(parser)

seg_betas <- unlist(strsplit(args$seg_beta,","))
seg_depths <- unlist(strsplit(args$seg_depth,","))
segment_betas <- do.call(rbind, lapply(seg_betas, fread))
segment_depths <- do.call(rbind, lapply(seg_depths, fread))
beta.mat <- as.matrix(segment_betas[,8:ncol(segment_betas)]) 
segment_depths <- as.matrix(segment_depths)

median_depth <- colMedians(segment_depths)

cor.mat <- cor(beta.mat)
mean_cors <- rowMeans(cor.mat)
tukey_outlier_limit <- quantile(mean_cors,.25) - 3*IQR(mean_cors) 

correlation_outliers <- names(which(mean_cors < tukey_outlier_limit))

cor_data <- melt(cor.mat)
custom_colors <- c("darkred", "red", "yellow", "white", "cyan", "blue", "darkblue")
cor_plot <- ggplot(cor_data %>% mutate(outlier=Var1%in%correlation_outliers), aes(Var1, Var2,  fill=value)) + geom_tile(col="white") + 
  scale_fill_gradientn(colors=custom_colors, limits = c(0, 1), name = "Correlation") +
  geom_rect(aes(xmin = as.numeric(factor(Var1)) - 0.5, 
                xmax = as.numeric(factor(Var1)) + 0.5,
                ymin = 0.5, ymax = length(unique(Var2)) + 0.5,
                color=outlier, linewidth=outlier),
            fill = NA) +
  scale_color_manual(values=c("white", "red")) + 
  scale_linewidth_manual(values=c(0,1.5)) +
  theme_minimal() + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom") 
mean_cor_hist <- data.frame(sample=names(mean_cors), mean_cor=mean_cors, outlier=names(mean_cors) %in% correlation_outliers) %>% 
  ggplot(aes(mean_cor, fill=outlier, label=ifelse(outlier,yes=sample,no=""))) + geom_histogram() + theme_minimal() + scale_fill_manual(values=c("grey80", "red")) + geom_text_repel(aes(y=1), color="red") +
  theme(legend.position="none") + ylab("Count") + xlab("mean pairwise correlation")
plot_grid(mean_cor_hist, cor_plot, ncol=1, rel_heights=c(2,8))
ggsave(args$correlation_plot_out, width=10, height=13)

#var_seg_betas <- findVariableSegments(segment_betas, plot_file=argv$varseg_plot_file)
#fwrite(var_seg_betas, file=argv$varsegment_bed_out, quote=F,row.names=F,col.names=T,sep="\t")

#cat("Calculating global methylation PCs . . . \n")
#global_meth_PCs <- findGlobalMethylationPCs(var_seg_betas, plot_dir=argv$pca_plot_dir)
#fwrite(global_meth_PCs, file=argv$global_meth_pcs_out, quote=F, row.names=T, col.names=T, sep="\t")

beta.mat.outliers_removed <- beta.mat[,!colnames(beta.mat) %in% correlation_outliers]
B <- beta.mat.outliers_removed
M <- log(B/(1-B))
Mpcs <- pca(M)
Mpcs$metadata <- data.frame(sample=colnames(beta.mat.outliers_removed), median_depth=median_depth[!colnames(beta.mat) %in% correlation_outliers])
biplot(Mpcs, lab=rownames(Mpcs$rotated))
ggsave(args$pc_biplot_out)

npcs <- findElbowPoint(Mpcs$variance)
hidden_factor_df <- data.frame(Sample_name=rownames(Mpcs$rotated), depth=Mpcs$metadata$median_depth, Mpcs$rotated[,1:npcs])

# add additional covariates if they are provided
if (!is.null(args$covariates)) {
    covariates <- fread(args$covariates)
    hidden_factor_df <- left_join(hidden_factor_df, covariates)
}

#sex  chromosome estimation analysis
sex_chrom_copynumber <- NULL
if(args$sex_seg_beta != "SKIP" && args$sex_seg_depth != "SKIP") {

  sex_seg_betas <- unlist(strsplit(args$sex_seg_beta,","))
  sex_seg_depths <- unlist(strsplit(args$sex_seg_depth,","))
  sex_segment_betas <- do.call(rbind, lapply(sex_seg_betas, fread))
  sex_segment_depths <- do.call(rbind, lapply(sex_seg_depths, fread))
  sex_beta.mat <- as.matrix(sex_segment_betas[,8:ncol(sex_segment_betas)]) 
  sex_segment_depths <- as.matrix(sex_segment_depths)


  chrX_median_depth <- colMedians(sex_segment_depths[sex_segment_betas$chrom == args$chrX_seqname,])
  chrY_median_depth <- colMedians(sex_segment_depths[sex_segment_betas$chrom == args$chrY_seqname,])

  sex_chrom_copynumber <- data.frame(Sample_name = colnames(sex_segment_depths),
                                     chrX = chrX_median_depth / median_depth * 2,
                                     chrY = chrY_median_depth / median_depth * 2)

  # remove correlation outliers
  sex_chrom_copynumber <- sex_chrom_copynumber[!colnames(sex_segment_depths) %in% correlation_outliers,]
  sex_beta.mat <- sex_beta.mat[,!colnames(sex_beta.mat) %in% correlation_outliers]
  sex_segment_depths <- sex_segment_depths[,!colnames(sex_segment_depths) %in% correlation_outliers]

  sex_chrom_clustering <- Mclust(sex_chrom_copynumber[,-1], G=2)
  sex_chrom_centers <- t(cbind(sex_chrom_clustering$parameters$mean, "XX"=c(2,0), "XY"=c(1,1)))
  clustering_distances <- as.matrix(dist(sex_chrom_centers))
  XX_index <- which.min(clustering_distances["XX",1:2])
  XY_index <- which.min(clustering_distances["XY", 1:2])
  sex_chrom_map = c(0,0)
  sex_chrom_map[XX_index] <- "XX"
  sex_chrom_map[XY_index] <- "XY"
  XX_index != XY_index

  sex_chrom_copynumber$cluster_probability <- sapply(1:nrow(sex_chrom_copynumber), function(x) {cluster=sex_chrom_clustering$classification[x]; 
                                                     pmvnorm(lower=as.numeric(sex_chrom_copynumber[x,2:3])-.25,
                                                              upper=as.numeric(sex_chrom_copynumber[x,2:3])+.25, 
                                                              mean=sex_chrom_clustering$parameters$mean[,cluster], 
                                                              sigma=sex_chrom_clustering$parameters$variance$sigma[,,cluster])})

  sex_chrom_copynumber$sex <- sapply(sex_chrom_clustering$classification, function(x) sex_chrom_map[x])
  sex_chrom_copynumber$sex[sex_chrom_copynumber$cluster_probability < .5] <- "ambiguous"

  B <- sex_beta.mat
  M <- log(B/(1-B))
  Mpcs <- pca(M)
  Mpcs$metadata <- data.frame(Sample_name=colnames(sex_beta.mat), sex_chrom_copynumber)

  depth_plot <- ggplot(sex_chrom_copynumber, aes(chrX, chrY, color=sex)) + geom_point(size=3) + ggtitle("Sex Chromosome Copy Number") + theme_minimal()
  gbiplot <- biplot(Mpcs, lab=rownames(Mpcs$rotated), colby="sex") + ggtitle("Sex chromosome methylation profiles")
  plot_grid(depth_plot, gbiplot)
  ggsave(args$sex_plot_out, width=10)

  #filter out individuals with ambiguous sex chromosomes (not XX or XY)
  ambiguous_sex_chrom_outlier <- filter(sex_chrom_copynumber, sex=="ambiguous") %>% pull(Sample_name)

  hidden_factor_df$sex <-  as.integer(sex_chrom_copynumber$sex == "XY")
  hidden_factor_df <- hidden_factor_df[!hidden_factor_df$Sample_name %in% ambiguous_sex_chrom_outlier,]

}
fwrite(sex_chrom_copynumber, args$sex_chrom_summary_out, row.names=F, col.names=T, sep="\t")
fwrite(hidden_factor_df, args$global_meth_pcs_out, row.names=F, col.names=T, sep="\t")
