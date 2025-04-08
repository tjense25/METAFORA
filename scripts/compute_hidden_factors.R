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
parser <- add_argument(parser, "--correlation_summary_out",help="where to write correlation analysis summary files")
parser <- add_argument(parser, "--sex_chrom_summary_out", help="where to write data frame of sex chromosome copy number estimates")
parser <- add_argument(parser, "--global_meth_pcs_out", help="where to write global PC covariates file ")
parser <- add_argument(parser, "--plot_out_dir", help="where to write correlation, PC, and sex plots")
parser <- add_argument(parser, "--chrX_seqname", help="seqname of the X chromosome in the reference")
parser <- add_argument(parser, "--chrY_seqname", help="seqname of the Y chromosome in the reference")

args <- NULL
args$seg_beta <- "../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr1.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr2.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr3.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr4.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr5.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr6.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr7.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr8.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr9.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr10.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr11.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr12.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr13.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr14.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr15.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr16.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr17.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr18.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr19.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr20.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr21.bed,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr22.bed"
args$seg_depth <- "../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr1.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr2.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr3.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr4.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr5.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr6.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr7.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr8.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr9.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr10.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr11.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr12.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr13.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr14.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr15.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr16.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr17.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr18.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr19.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr20.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr21.mat,../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.chrom_chr22.mat"
args$sex_seg_beta <- "SKIP"
args$sex_seg_depth <- "SKIP"         
args$correlation_summary_out <- "../METAFORA_output/Global_Methylation_PCA_tissue_Blood/Mean_pairwise_correlation_summary.txt"         
args$sex_chrom_summary_out <- "../METAFORA_output/Global_Methylation_PCA_tissue_Blood/Sex_chromosome_estimates_summary.txt" 
args$global_meth_pcs_out <- "../METAFORA_output/Global_Methylation_PCA_tissue_Blood/PCA_covariates.txt" 
args$plot_out_dir <- "../METAFORA_output/Global_Methylation_PCA_tissue_Blood/"
args$chrX_seqname <- "chrX"
args$chrY_seqname <- "chrY"
args$covariates <- "../MOTRPAC_pilot.covariates.txt"


args <- parse_args(parser)
plot_out <- args$plot_out_dir
seg_betas <- unlist(strsplit(args$seg_beta,","))
seg_depths <- unlist(strsplit(args$seg_depth,","))
segment_betas <- do.call(rbind, lapply(seg_betas, fread))
segment_depths <- do.call(rbind, lapply(seg_depths, fread))
beta.mat <- as.matrix(segment_betas[,8:ncol(segment_betas)]) 
segment_depths <- as.matrix(segment_depths)

median_depth <- colMedians(segment_depths)

cor.mat <- cor(beta.mat)
mean_cors <- rowMeans(cor.mat)

#find tukey outlier limit. (if its too high set to 0.95 to avoid filtering too many samples)
tukey_outlier_limit <- min(0.95, quantile(mean_cors,.25) - 3*IQR(mean_cors) )

correlation_outliers <- names(which(mean_cors < tukey_outlier_limit))

cor_data <- reshape2::melt(cor.mat) %>% mutate(outlier=(Var1%in%correlation_outliers)|(Var2%in%correlation_outliers))
custom_colors <- c("darkred", "red", "yellow", "white", "cyan", "blue", "darkblue")
cor_plot <- ggplot(cor_data, aes(Var1, Var2,  fill=value,color=outlier)) + geom_tile() + 
  scale_fill_gradientn(colors=custom_colors, limits = c(0, 1), name = "Correlation") +
  scale_color_manual(values=c("white", "red")) + 
  theme_minimal() + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom") 
mean_cor_hist <- data.frame(sample=names(mean_cors), mean_cor=mean_cors, outlier=names(mean_cors) %in% correlation_outliers) %>% 
  ggplot(aes(mean_cor, fill=outlier, label=ifelse(outlier,yes=sample,no=""))) + geom_histogram() + theme_minimal() + scale_fill_manual(values=c("grey80", "red")) + geom_text_repel(aes(y=1), color="red") +
  theme(legend.position="none") + ylab("Count") + xlab("mean pairwise correlation")
plot_grid(mean_cor_hist, cor_plot, ncol=1, rel_heights=c(2,8))
ggsave(paste0(plot_out,"/Correlation_matrix.mean_cor_distribution.correlation_outliers.pdf"), width=10, height=13)


mean_correlation_df <- data.frame(Sample_name=colnames(segment_depths), mean_pairwise_correlation=mean_cors, median_depth)
if (!args$covariates=="SKIP") {
    covariates <- fread(args$covariates)
    mean_correlation_df %<>% left_join(covariates)
    if ("Batch" %in% colnames(mean_correlation_df)) {
        #Look at Batch Plots
        unique_batches <- unique(mean_correlation_df$Batch)
        median_batch_mat <- do.call(cbind, lapply(unique_batches, function(this_batch) rowMedians(beta.mat[,mean_correlation_df$Batch==this_batch])))
        colnames(median_batch_mat) <- unique_batches
        batch_cor_mat <- cor(median_batch_mat) %>% reshape2::melt()
        ggplot(batch_cor_mat, aes(Var1, Var2,  fill=value)) + geom_tile(color="white") + 
          geom_text(aes(label = sprintf("%.4f", value)), size = 3, color="white") +  # Add correlation values
          scale_fill_gradientn(colors=custom_colors, limits = c(0, 1), name = "Correlation") +
          theme_minimal() + xlab("") + ylab("") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom") 
        ggsave(paste0(plot_out, "/Median_batch_profile.batch_correlation_matrix.pdf"))
    }
}
if(!"Batch" %in% colnames(mean_correlation_df)) {
    mean_correlation_df$Batch <- "Full Cohort"
}
mean_correlation_df$outlier <- mean_correlation_df$Sample_name %in% correlation_outliers
ggplot(mean_correlation_df, aes(median_depth, mean_pairwise_correlation, label=ifelse(outlier,Sample_name,""), color=Batch)) + geom_point() + theme_minimal() + geom_text_repel(color="Red") + 
    geom_hline(yintercept=tukey_outlier_limit, linetype="dashed", color="red")
ggsave(paste0(plot_out,"/mean_pairwise_correlation.median_depth.scatter_plot.pdf"))
cor_data$Batch1 <- mean_correlation_df$Batch[match(cor_data$Var1, mean_correlation_df$Sample_name)]
cor_data$Batch2 <- mean_correlation_df$Batch[match(cor_data$Var2, mean_correlation_df$Sample_name)]
ggplot(cor_data, aes(Batch2, value, fill=Batch2)) + geom_violin(alpha=.7) + facet_wrap(~Batch1) + 
    theme_minimal() +
    geom_hline(yintercept=tukey_outlier_limit, linetype="dashed", color="red", linewidth=2) +
    ylab("pairwise sample correlation (r)") + 
    xlab("Batch") + 
    theme(panel.border=element_rect(color="black", fill=NA, linewidth=2), axis.text.x=element_text(hjust=1, angle=90), text=element_text(size=18))
ggsave(paste0(plot_out, "/Pairwise_correlations.across_batches.violin_plots.pdf"), width=15, height=15)
fwrite(mean_correlation_df, args$correlation_summary_out, row.names=F, col.names=T, sep="\t")

beta.mat.outliers_removed <- beta.mat[,!colnames(beta.mat) %in% correlation_outliers]
B <- beta.mat.outliers_removed
M <- log(B/(1-B))
Mpcs <- pca(M)
Mpcs$metadata <- data.frame(Sample_name=colnames(beta.mat.outliers_removed), median_depth=median_depth[!colnames(beta.mat) %in% correlation_outliers])
biplot(Mpcs, lab=rownames(Mpcs$rotated))
ggsave(paste0(plot_out, "/Global_Hidden_factors.PCA_biplot.pdf"))

npcs <- findElbowPoint(Mpcs$variance)
hidden_factor_df <- data.frame(Sample_name=rownames(Mpcs$rotated), depth=Mpcs$metadata$median_depth, Mpcs$rotated[,1:npcs])

# add additional covariates if they are provided
if (!args$covariates=="SKIP") {
    covariates <- fread(args$covariates)
    hidden_factor_df <- left_join(hidden_factor_df, covariates)
    if ("Batch" %in% colnames(covariates)) {
        Mpcs$metadata$Batch <- covariates$Batch[match(Mpcs$metadata$Sample_name, covariates$Sample_name)]
        biplot(Mpcs, lab=rownames(Mpcs$rotated), colby="Batch") + theme(legend.position="right")
        ggsave(paste0(plot_out, "/Global_Hidden_factors.PCA_biplot.pdf"), width=8)
    }
}


#sex  chromosome estimation analysis
sex_chrom_copynumber <- data.frame()
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
  XX_index
  XY_index

  sex_chrom_copynumber$cluster_probability <- sapply(1:nrow(sex_chrom_copynumber), function(x) {cluster=sex_chrom_clustering$classification[x]; 
                                                     pmvnorm(lower=as.numeric(sex_chrom_copynumber[x,2:3])-.5,
                                                              upper=as.numeric(sex_chrom_copynumber[x,2:3])+.5, 
                                                              mean=sex_chrom_clustering$parameters$mean[,cluster], 
                                                              sigma=sex_chrom_clustering$parameters$variance$sigma[,,cluster])})

  sex_chrom_copynumber$sex <- sapply(sex_chrom_clustering$classification, function(x) sex_chrom_map[x])
  sex_chrom_copynumber$sex[sex_chrom_copynumber$cluster_probability < .8] <- "ambiguous"

  B <- sex_beta.mat
  M <- log(B/(1-B))
  Mpcs <- pca(M)
  Mpcs$metadata <- data.frame(Sample_name=colnames(sex_beta.mat), sex_chrom_copynumber)

  depth_plot <- ggplot(sex_chrom_copynumber, aes(chrX, chrY, color=sex)) + geom_point(size=3) + ggtitle("Sex Chromosome Copy Number") + theme_minimal() +
    xlab("X chromosome coverage ratio") + ylab("Y chromosome coverage ratio") +
      theme(panel.border = element_rect(color="black", fill=NA, size=1),
             text=element_text(size=16),
             axis.title = element_text(size=18, face="bold"),
             plot.title=element_text(size=20, face="bold"))
  gbiplot <- biplot(Mpcs, lab=rownames(Mpcs$rotated), colby="sex") + ggtitle("Sex chromosome methylation profiles")
  plot_grid(depth_plot, gbiplot)
  ggsave(paste0(plot_out, "/Sex_chromosome_estimate_plots.pdf"), width=10) 

  #filter out individuals with ambiguous sex chromosomes (not XX or XY)
  ambiguous_sex_chrom_outlier <- filter(sex_chrom_copynumber, sex=="ambiguous") %>% pull(Sample_name)

  hidden_factor_df$sex <-  as.integer(sex_chrom_copynumber$sex == "XY")
  hidden_factor_df <- hidden_factor_df[!hidden_factor_df$Sample_name %in% ambiguous_sex_chrom_outlier,]

}
fwrite(sex_chrom_copynumber, args$sex_chrom_summary_out, row.names=F, col.names=T, sep="\t")
fwrite(hidden_factor_df, args$global_meth_pcs_out, row.names=F, col.names=T, sep="\t")
