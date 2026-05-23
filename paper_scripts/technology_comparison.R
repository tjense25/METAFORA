library(data.table)
library(cowplot)
library(corrplot)
library(ggrepel)
library(limma)
library(edgeR)
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

valid_chrom <- paste0("chr", seq(1,22))
segment_betas <- NULL
segment_depths <- NULL
for (chr in valid_chrom) {
    tmp.seg <- fread(paste0("./output/methylation_results/Population_methylation.GRCh38.tissue_Blood/Meth_segments.GRCh38.tissue_Blood.segment_betas.chrom_",chr,".bed"))
    segment_betas <-  rbind(segment_betas, tmp.seg)
    tmp.depths <- fread(paste0("./output/methylation_results/Population_methylation.GRCh38.tissue_Blood/Meth_segments.GRCh38.tissue_Blood.segment_coverage.chrom_",chr,".mat"))
    segment_depths <- rbind(segment_depths,tmp.depths)
}
beta.mat <- as.matrix(segment_betas[,8:ncol(segment_betas)])
segment_depths <- as.matrix(segment_depths)
# create Segment x CpG identity matrix to indicate which cpgs belong to which segment
tech_meta <- fread("./U08.sample_table.metadat.txt")
seq_tech <- tech_meta$Technology[match(colnames(beta.mat), tech_meta$Sample_name)]
site <- tech_meta$Site[match(colnames(beta.mat), tech_meta$Sample_name)]
median_depth <- colMedians(segment_depths)

pb_mean <- rowMedians(beta.mat[,seq_tech=="PacBio"], na.rm=T)
ONT_mean <- rowMedians(beta.mat[,seq_tech=="ONT_R9"], na.rm=T)
R10_mean <- rowMedians(beta.mat[,seq_tech=="ONT_R10"], na.rm=T)
pb_mean_depth <- rowMeans(segment_depths[,seq_tech=="PacBio"])
ONT_mean_depth <- rowMeans(segment_depths[,seq_tech=="ONT_R9"])
R10_mean_depth <- rowMeans(segment_depths[,seq_tech=="ONT_R10"])
#join_mean <- rowMeans(beta.mat)

#permutations <- t(replicate(n=10000, sample(seq_tech)))
#cors <- do.call(c,pbmclapply(1:nrow(permutations), function(x) {cor(rowSums(beta.mat[,permutations[x,]=="PacBio"]), rowSums(beta.mat[,permutations[x,]=="ONT"]))}))
#print(paste0("pb_vs_r9: ", cor(pb_mean,ONT_mean), " pb_vs_r10: ", cor(pb_mean, R10_mean), " r9_vs_r10: ", cor(ONT_mean,R10_mean)))

cor.mat <- cor(beta.mat)
mean_cors <- rowMeans(cor.mat)
ggplot(NULL, aes(median_depth, mean_cors, label=names(mean_cors), shape=site, color=seq_tech)) + geom_point() + geom_text_repel(color="black") + theme(legend.position="none") + scale_color_manual(values=c("turquoise", "blue", "hotpink")) + scale_y_log10()
ggsave("./meth_plots/median_depth_by_mean_correlation.pdf")

pdf("./meth_plots/corrplot.pdf", width=12, height=12)
corrplot(corr=cor.mat, is.corr=T, order='AOE')
dev.off()

transformed.mat <- -log10(1-cor.mat)
diag(transformed.mat) <- 2
pdf("./meth_plots/corrplot.transformed.pdf", width=12, height=12)
corrplot(corr=transformed.mat, order='AOE', is.corr=F)
dev.off()

correlation_outliers <- names(which(mean_cors < .95))
correlation_outliers


cor.mat_cleaned <- cor(beta.mat[,!colnames(beta.mat) %in% correlation_outliers])
cor.melted <- melt(cor.mat_cleaned)
colnames(cor.melted) <- c("Sample1", "Sample2", "cor")
cor.melted %<>% filter(Sample1 != Sample2)
cor.melted$Tech1 <- tech_meta$Technology[match(cor.melted$Sample1, tech_meta$Sample_name)]
cor.melted$Tech2 <- tech_meta$Technology[match(cor.melted$Sample2, tech_meta$Sample_name)]
cor.melted$Site1 <- tech_meta$Site[match(cor.melted$Sample1, tech_meta$Sample_name)]
cor.melted$Site2 <- tech_meta$Site[match(cor.melted$Sample2, tech_meta$Sample_name)]
cor.melted$SameSite <- cor.melted$Site1==cor.melted$Site2
cor.melted$SameTech <- cor.melted$Tech1 == cor.melted$Tech2

cor.melted
median(cor.melted$cor)
result=t.test(data=cor.melted %>% filter(Tech1=="ONT_R10",Tech2=="ONT_R10"), cor~SameSite, method="wilcox")
mean(result$conf.int)
result=t.test(data=cor.melted %>% filter(Tech1=="PacBio",Tech2=="PacBio"), cor~SameSite, method="wilcox")
mean(result$conf.int)
cor.melted %>% filter(Tech1=="ONT_R10", SameTech) %>% group_by(SameSite) %>% summarize(cor=mean(cor))
cor.melted %>% filter(Tech1=="PacBio", SameTech) %>% group_by(SameSite) %>% summarize(cor=mean(cor))


cor.melted %>% filter(Tech1=="ONT_R10", Tech2=="ONT_R9")
t.test(cor.melted %>% filter(Tech1=="ONT_R10", Tech2%in%c("ONT_R10","ONT_R9")), cor~SameTech, method="wilcox")
cor.melted %>% group_by(Tech1,Tech2) %>% summarize(cor=mean(cor))
cor.melted %>% filter(Tech1=="ONT_R10") %>% group_by(Tech2) %>% summarize(cor=mean(cor))


site_colors=list("GSS"="#d4e09b", "Broad"="#f6f4d2", "PMGRC"="#cbdfbd", "UW"="#f19c79")
table(seq_tech)
tech_colors=list("ONT_R10"="turquoise", "ONT_R9"="blue", "PacBio"="hotpink")
ggplot(cor.melted, aes(Tech1, y=cor, fill=Tech1)) + geom_violin() + facet_wrap(~Tech2) + theme_classic() + scale_fill_manual(values=tech_colors) + ylab("Correlation (pearson r)")  + theme(panel.grid.major.y=element_line())
ggsave("./meth_plots/correlation_violines.pdf")

cor.melted %>% filter(Tech1=="ONT_R10", Tech2=="ONT_R10") %>%
    ggplot(aes(SameSite, y=cor, fill=Tech1, alpha=SameSite)) + geom_violin() + facet_grid(.~Site1) + scale_fill_manual(values=tech_colors) + scale_alpha_manual(values=c(.4,1)) +
    theme_minimal() + theme(legend.position="none")
ggsave("./meth_plots/ont_site_correlation.pdf", width=5)
cor.melted %>% filter(Tech1=="PacBio", Tech2=="PacBio") %>%
    ggplot(aes(SameSite, y=cor, fill=Tech1, alpha=SameSite)) + geom_violin() + facet_grid(.~Site1) + scale_fill_manual(values=tech_colors) + scale_alpha_manual(values=c(.4,1)) +
    theme_minimal() + theme(legend.position="none")
ggsave("./meth_plots/pacbio_site_correlation.pdf", width=2.5)

cor.melted %>% filter(Tech1=="ONT_R10", Tech2=="ONT_R10") %>%
    ggplot(aes(SameSite, y=cor, fill=Tech1, alpha=SameSite)) + geom_violin() + facet_grid(.~Site1) + scale_fill_manual(values=tech_colors) + scale_alpha_manual(values=c(.4,1)) +
    theme_minimal()
ggplot(cor.melted %>% filter(Tech1 == "PacBio", Tech2=="PacBio"), aes(Site1, y=cor, fill=Site1)) + geom_violin() + facet_grid(.~Site2) + scale_fill_manual(values=site_colors)
ggsave("./meth_plots/pacbio_site_correlation.pdf", width=8)
ggplot(cor.melted %>% filter(Tech1 == "ONT_R10", Tech2=="ONT_R10"), aes(Site1, y=cor, fill=Site1)) + geom_violin() + facet_grid(.~Site2) + scale_fill_manual(values=site_colors)
ggsave("./meth_plots/ont_site_correlation.pdf", width=12)


ggplot(NULL, aes(x=rowMeans(cor.mat))) + geom_histogram()
ggsave("./meth_plots/mean_cross_sample_correlation.pdf")

tech_mean_segment_beta <- data.frame(Tech=c(rep("PacBio", length(pb_mean)), rep("ONT_R9",length(ONT_mean)), rep("ONT_R10",length(R10_mean))), Mean_beta=c(pb_mean, ONT_mean, R10_mean), mean_depth=c(pb_mean_depth,ONT_mean_depth, R10_mean_depth))
ggplot(tech_mean_segment_beta, aes(x=Mean_beta, fill=Tech)) + geom_density(alpha=.4) + 
    scale_fill_manual(values=c("turquoise", "blue", "hotpink")) + 
    theme_classic() + theme(panel.grid.major.x=element_line(linewidth=1), panel.grid.minor.x=element_line(linewidth=.5))
ggsave("./meth_plots/tech_mean_segment_beta.density.pdf")
ggplot(tech_mean_segment_beta, aes(x=mean_depth, fill=Tech)) + geom_density(alpha=.4) +
    scale_fill_manual(values=c("turquoise", "blue", "hotpink")) + xlim(0,40) +
    theme_classic() + theme(panel.grid.major.x=element_line(linewidth=1), panel.grid.minor.x=element_line(linewidth=.5))
ggsave("./meth_plots/tech_mean_segment_depth.density.pdf")


cpg_density_vs_depth <- data.frame(cpg_density = segment_betas$num_cpg/(segment_betas$end - segment_betas$start + 1), ONT_mean_depth, pb_mean_depth, R10_mean_depth)
cpg_density_vs_depth$density_quantile <- factor(round(percent_rank(cpg_density_vs_depth$cpg_density),1))
ggplot(cpg_density_vs_depth, aes(density_quantile,ONT_mean_depth)) + geom_violin(fill="blue", alpha=.5) + theme_classic() + ylim(0,20) +
    theme(panel.grid.major.y=element_line(linewidth=1)) + ylab("R9 Segment Mean Depth")
ggsave('ONT_cpg_density_vs_depth.violins.pdf')
ggplot(cpg_density_vs_depth, aes(density_quantile,pb_mean_depth)) + geom_violin(fill="hotpink", alpha=.5)  + ylim(0,35) +  theme_classic() +
    theme(panel.grid.major.y=element_line(linewidth=1)) + ylab("PacBio Segment Mean Depth")
ggsave('PB_cpg_density_vs_depth.violins.pdf')

ggplot(cpg_density_vs_depth, aes(density_quantile,pb_mean_depth)) + geom_violin(fill="turquoise", alpha=.5)  + ylim(0,40) +  theme_classic() + 
    theme(panel.grid.major.y=element_line(linewidth=1)) + ylab("R10 Segment Mean Depth")
ggsave('R10_cpg_density_vs_depth.violins.pdf')


B <- beta.mat
M <- log(B/(1-B))
Mpcs <- pca(M)
Mpcs$metadata <- data.frame(sample=colnames(beta.mat), seq_tech, site, median_depth)
biplot(Mpcs, lab=rownames(Mpcs$rotated), colby="seq_tech", shape="site") + scale_color_manual(values=c("turquoise", "blue", "hotpink"))
ggsave('meth_plots/r10_biplot.pdf')

beta.mat.outliers_removed <- beta.mat[,!colnames(beta.mat) %in% correlation_outliers]
B <- beta.mat.outliers_removed
M <- log(B/(1-B))
Mpcs <- pca(M)
Mpcs$metadata <- data.frame(sample=colnames(beta.mat.outliers_removed), seq_tech=seq_tech[!colnames(beta.mat) %in% correlation_outliers], site=site[!colnames(beta.mat) %in% correlation_outliers], median_depth=median_depth[!colnames(beta.mat) %in% correlation_outliers])
Mpcs$metadata
biplot(Mpcs, lab=rownames(Mpcs$rotated), colby="seq_tech", shape="site") + scale_color_manual(values=c("turquoise","blue","hotpink"))
ggsave('meth_plots/r10_biplot.global_removed.pdf')


median_depth.capped <- pmin(median_depth, 35)
data.frame(Mpcs$metadata, Mpcs$rotated) %>% 
    mutate(median_depth_capped = pmin(median_depth, 35)) %>%
    ggplot(aes(PC1,PC2, color=median_depth_capped, shape=seq_tech)) + geom_point(size=3) + theme_classic() +
    theme(panel.grid.major=element_line())
ggsave('meth_plots/r10_biplot.by_depth.pdf')
pc1_depth <- data.frame(Mpcs$metadata, Mpcs$rotated) %>% 
    ggplot(aes(PC1,median_depth)) + geom_point(size=3, aes(shape=seq_tech)) + theme_classic() + 
    geom_smooth(method="lm", color="red") + 
    theme(panel.grid.major=element_line()) + ylab("Median Depth") + theme(text=element_text(size=20), legend.position="none")
pc2_depth <- data.frame(Mpcs$metadata, Mpcs$rotated) %>% 
    ggplot(aes(PC2,median_depth)) + geom_point(size=3, aes(shape=seq_tech)) + theme_classic() +
    geom_smooth(method="lm", color="red") + 
    theme(panel.grid.major=element_line())+ ylab("Median Depth") + theme(text= element_text(size=20), legend.position="none")
plot_grid(pc1_depth, pc2_depth, ncol=1) 
ggsave('meth_plots/PC_depth.pdf', width=5, height=10)

cor(Mpcs$rotated$PC1, Mpcs$metadata$median_depth)
cor(Mpcs$rotated$PC2, Mpcs$metadata$median_depth)

covars <- data.frame(Technology=Mpcs$metadata$seq_tech, median_depth=Mpcs$metadata$median_depth, Mpcs$rotated[,1:7])
#covars <- data.frame(median_depth=Mpcs$metadata$median_depth)

global_meth_outliers <- correlation_outliers
M.corrected <- removeBatchEffect(M, batch=covars$Technology, covariates = cbind(covars$median_depth))
M.depth_corrected <- removeBatchEffect(M, batch=NULL, covariates = cbind(covars$median_depth))
Mpcs.corrected <- pca(M.corrected)
Mpcs.depth_corrected <- pca(M.depth_corrected)
Mpcs.corrected$metadata <- data.frame(sample=colnames(beta.mat.outliers_removed), seq_tech = seq_tech[!colnames(beta.mat) %in% global_meth_outliers], median_depth=median_depth[!colnames(beta.mat) %in% global_meth_outliers], site=site[!colnames(beta.mat) %in% global_meth_outliers])
Mpcs.depth_corrected$metadata <- data.frame(sample=colnames(beta.mat.outliers_removed), seq_tech = seq_tech[!colnames(beta.mat) %in% global_meth_outliers], median_depth=median_depth[!colnames(beta.mat) %in% global_meth_outliers], site=site[!colnames(beta.mat) %in% global_meth_outliers])
#biplot(Mpcs.corrected, lab="", colby="seq_tech", shape="site") + scale_color_manual(values=c("turquoise", "blue", "hotpink"))
biplot(Mpcs.depth_corrected, lab="", colby="seq_tech", shape="site") + scale_color_manual(values=c("dark blue", "turquoise", "hotpink"))
ggsave("./meth_plots/r10_biplot.global_removed.corrected.pdf")

findElbowPoint(variance=Mpcs.corrected$sdev^2)
colnames(segment_betas)

covars$Batch <- covars$Technology
covars$Technology <- NULL
covars$Site <- NULL
covars
fwrite(segment_betas, "./output/methylation_results/Meth_segments.GRCh38.tissue_Blood.segment_betas.bed")
fwrite(covars, "./output/methylation_results/Gloabal_Methylation_PCA_GRCh38_tissue_Blood/PCA_covariates.txt",row.names=T)

covars_to_correct <- Mpcs$metadata 
lr.meta <- data.frame(sample=colnames(beta.mat), Site=site)
covars_to_correct %<>% left_join(lr.meta)
table(covars_to_correct$seq_tech, covars_to_correct$Site)

beta.mat_corrected <- removeBatchEffect(beta.mat.outliers_removed, batch= covars_to_correct$seq_tech, covariates = covars_to_correct$median_depth)
corrected_PCs <- PCAtools::pca(beta.mat_corrected)
corrected_PCs$metadata <- covars_to_correct
biplot(corrected_PCs, lab=rownames(corrected_PCs$rotated), colby="seq_tech", shape="Site") + scale_color_manual(values=c("dark blue", "turquoise", "hotpink"))
ggsave("meth_plots/corrected_pcs.pdf")

beta.mat_corrected <- as.matrix(beta.mat_corrected)
pb_mean <- rowMedians(beta.mat_corrected[,covars_to_correct$seq_tech=="PacBio"], na.rm=T)
ONT_mean <- rowMedians(beta.mat_corrected[,covars_to_correct$seq_tech=="ONT_R9"], na.rm=T)
R10_mean <- rowMedians(beta.mat_corrected[,covars_to_correct$seq_tech=="ONT_R10"], na.rm=T)

tech_mean_segment_beta <- data.frame(Tech=c(rep("PacBio", length(pb_mean)), rep("ONT_R9",length(ONT_mean)), rep("ONT_R10",length(R10_mean))), Mean_beta=c(pb_mean, ONT_mean, R10_mean), mean_depth=c(pb_mean_depth,ONT_mean_depth, R10_mean_depth))
ggplot(tech_mean_segment_beta, aes(x=Mean_beta, fill=Tech)) + geom_density(alpha=.2) + 
    scale_fill_manual(values=c("turquoise", "blue", "hotpink")) + 
    theme_classic() + theme(panel.grid.major.x=element_line(linewidth=1), panel.grid.minor.x=element_line(linewidth=.5))
ggsave("./meth_plots/corrected.tech_mean_segment_beta.pdf")

covars_to_correct <- data.frame(sample=colnames(beta.mat_corrected), seq_tech, median_depth)
beta.mat_corrected <- removeBatchEffect(beta.mat, batch= covars_to_correct$seq_tech, covariates = covars_to_correct$median_depth)
corrected_PCs <- PCAtools::pca(beta.mat_corrected)

#### copy code for Fibroblast
segment_betas <- NULL
segment_depths <- NULL
for (chr in valid_chrom) {
    tmp.seg <- fread(paste0("./output/methylation_results/Population_methylation.GRCh38.tissue_Fibro/Meth_segments.GRCh38.tissue_Fibro.segment_betas.chrom_",chr,".bed"))
    segment_betas <-  rbind(segment_betas, tmp.seg)
    tmp.depths <- fread(paste0("./output/methylation_results/Population_methylation.GRCh38.tissue_Fibro/Meth_segments.GRCh38.tissue_Fibro.segment_coverage.chrom_",chr,".mat"))
    segment_depths <- rbind(segment_depths,tmp.depths)
}
nrow(segment_betas)

beta.mat <- as.matrix(segment_betas[,8:ncol(segment_betas)])
segment_depths <- as.matrix(segment_depths)
median_depth <- colMedians(segment_depths)

tech_meta <- fread("./Methylation.sample_table.metadat.txt")
seq_tech <- tech_meta$Technology[match(colnames(beta.mat), tech_meta$Sample_name)]
ONT_mean <- rowMedians(beta.mat[,seq_tech=="ONT_R9"])
R10_mean <- rowMedians(beta.mat[,seq_tech=="ONT_R10"])
ONT_mean_depth <- rowMeans(segment_depths[,seq_tech=="ONT_R9"])
R10_mean_depth <- rowMeans(segment_depths[,seq_tech=="ONT_R10"])

tech_compare = data.frame(ONT_R9=ONT_mean, ONT_R10=R10_mean)
ggplot(tech_compare, aes(ONT_R9, ONT_R10)) + geom_point(alpha=.2, color="grey80") + geom_density2d()
ggsave('./meth_plots/Fibro.ONT.tech_compare.geom_density2d.pdf', width=10)

tech_cor <- cor(ONT_mean, R10_mean)
tech_cor

cor.mat <- cor(beta.mat)
pdf("./meth_plots/Fibro.corrplot.pdf", width=12, height=12)
corrplot(corr=cor.mat, is.corr=T)
dev.off()
transformed.mat <- -log10(1-cor.mat)
diag(transformed.mat) <- 2
pdf("./meth_plots/Fibro.corrplot.transformed.pdf", width=12, height=12)
corrplot(corr=transformed.mat, order='AOE', is.corr=F)
dev.off()
#beta.mat_seqcorrected <- removeBatchEffect(beta.mat, batch = seq_tech)
#beta.mat_seqcorrected[beta.mat_seqcorrected<0] <- 0.01
#segment_betas <- data.table(segment_betas[,1:7], beta.mat_seqcorrected)

mean_cor = rowMeans(cor.mat)
ggplot(NULL,(aes(median_depth, mean_cor, color=seq_tech, label=colnames(beta.mat)))) + geom_point() + geom_label_repel()
ggsave('./meth_plots/Fibro.mean_correlation.label_plot.pdf')

tech_mean_segment_beta <- data.frame(Tech=c(rep("ONT_R9",length(ONT_mean)), rep("ONT_R10",length(R10_mean))), Mean_beta=c(ONT_mean, R10_mean), mean_depth=c(ONT_mean_depth, R10_mean_depth))
ggplot(tech_mean_segment_beta, aes(x=Mean_beta, fill=Tech)) + geom_density(alpha=.4) +
    scale_fill_manual(values=c("turquoise", "dark blue"))
ggsave("meth_plots/Fibro.tech_mean_segment_beta.density.pdf")
ggplot(tech_mean_segment_beta, aes(x=mean_depth, fill=Tech)) + geom_density(alpha=.4) +
    scale_fill_manual(values=c("turquoise", "dark blue")) + xlim(0,50)
ggsave("meth_plots/Fibro.tech_mean_segment_depth.density.pdf")

B <- beta.mat
M <- log(B/(1-B))
Mpcs <- pca(M)
Mpcs$metadata <- data.frame(sample=colnames(beta.mat), seq_tech, median_depth)
biplot(Mpcs, lab=rownames(Mpcs$rotated), colby="seq_tech") + scale_color_manual(values=c("dark blue", "turquoise"))
ggsave('meth_plots/Fibro.r10_biplot.pdf')

global_meth_outliers <- c("UDN336336_Fibro", "UDN969133_Fibro", "UDN720761_Fibro", "UDN249098_Fibro", "UDN890454_Fibro", "UDN154285_Fibro", "UDN248063_Fibro")
beta.mat.outliers_removed <- beta.mat[,!colnames(beta.mat) %in% global_meth_outliers]
B <- beta.mat.outliers_removed
M <- log(B/(1-B))
Mpcs <- pca(M)
Mpcs$metadata <- data.frame(sample=colnames(beta.mat.outliers_removed), seq_tech = seq_tech[!colnames(beta.mat) %in% global_meth_outliers], median_depth=median_depth[!colnames(beta.mat) %in% global_meth_outliers])
biplot(Mpcs, lab=rownames(Mpcs$rotated), colby="seq_tech") + scale_color_manual(values=c("dark blue", "turquoise", "hotpink"))
ggsave('meth_plots/Fibro.r10_biplot.global_removed.pdf')

elbow <- findElbowPoint(variance=Mpcs$sdev^2)
covars <- data.frame(Batch=Mpcs$metadata$seq_tech, median_depth=Mpcs$metadata$median_depth, Mpcs$rotated[,1:elbow])
covars_to_correct <- Mpcs$metadata 

M_corrected <- removeBatchEffect(M, batch=covars_to_correct$seq_tech, covariates = covars_to_correct$median_depth)
corrected_PCs <- PCAtools::pca(M_corrected)
corrected_PCs$metadata <- covars_to_correct
biplot(corrected_PCs, lab=rownames(corrected_PCs$rotated), colby="seq_tech") + scale_color_manual(values=c("dark blue", "turquoise"))
ggsave("meth_plots/Fibro.corrected_pcs.pdf")


covars$Batch <- NULL
fwrite(segment_betas, "./output/methylation_results/Meth_segments.GRCh38.tissue_Fibro.segment_betas.bed")
fwrite(covars, "./output/methylation_results/Gloabal_Methylation_PCA_GRCh38_tissue_Fibro/PCA_covariates.txt",row.names=T)

seq_stats <- fread("./output/GREGoR_U07.nanoplot_summarized.txt")
seq_stats %>% left_join(tech_meta) %>% group_by(Technology) %>% summarize(median_n50=median(read_length_n50), median_coverage=median(total_bases/3))
seq_stats %>% left_join(tech_meta) %>% filter(read_length_n50>10) %>% group_by(Technology) %>% summarize(median_n50=median(read_length_n50), median_coverage=median(total_bases/3), minrange=min(read_length_n50),maxrange=max(read_length_n50))
