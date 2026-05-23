library(data.table)
library(tidyverse)
library(matrixStats)
library(ggridges)
library(magrittr)
library(plyranges)

chrM.betas <- fread("./METAFORA_output/Population_methylation.tissue_Blood/Population_methylation.tissue_Blood.chrom_chrM.betas.mat.gz")
chrM.depths <- fread("./METAFORA_output/Population_methylation.tissue_Blood/Population_methylation.tissue_Blood.chrom_chrM.coverage.mat.gz")
covariates <- fread("./U08_GREGoR.covariates.txt")
covariates

beta.mat <- as.matrix(chrM.betas[,4:ncol(chrM.betas)])
depth.mat <- as.matrix(chrM.depths[,4:ncol(chrM.depths)])

estimated_error <- colMeans(beta.mat*depth.mat/(depth.mat+2), na.rm=T)
median_depth = colMedians(depth.mat, na.rm=T)

seqtech_colors=c("PacBio"="hotpink", "ONT_R9"="blue", "ONT_R10"="turquoise")
estimated_error <- data.frame(Sample_name=colnames(beta.mat), estimated_error, median_depth) %>% left_join(covariates)
estimated_error$Batch %<>% factor(levels=c("PacBio", "ONT_R9", "ONT_R10"))
#ggplot(estimated_error, aes(estimated_error,fill=Batch)) + geom_histogram(color="black", alpha=.7, position="identity") + facet_wrap(~Batch,ncol=1, scale="free") + scale_fill_manual(values=seqtech_colors) + xlim(0,.2) + theme_minimal()
ggplot(estimated_error, aes(estimated_error,Batch,fill=Batch)) + geom_density_ridges2(color="black", alpha=.3, scale=5) + scale_fill_manual(values=seqtech_colors) + xlim(0,.2) + theme_minimal() + 
    geom_vline(xintercept=.006,color="lightblue", linewidth=1, linetype="dashed") +
    geom_vline(xintercept=.022, color="hotpink", linewidth=1, linetype="dashed") + 
    geom_vline(xintercept=.138, color="navyblue", linewidth=1, linetype="dashed")
ggsave("estimated_error_profiles.pdf")
estimated_error %>% filter(estimated_error==0.02270735)
estimated_error %>% filter(estimated_error==0.1380456)
estimated_error %>% filter(estimated_error==0.1380456)

estimated_error
#ggplot(estimated_error, aes(median_depth, estimated_error, color=Batch)) + geom_point(size=2) + scale_color_manual(values=seqtech_colors) + ylim(0,.2) + theme_minimal()
#ggsave("estimated_error_profiles.pdf")

outliers_combined <- fread("./METAFORA_output/METAFORA.tissue_Blood.methylation_outliers.combined.tsv")
outlier_count <- outliers_combined %>% group_by(ID) %>% summarize(N_outliers=n()) %>% dplyr::rename(Sample_name=ID)

by(estimated_error$estimated_error, estimated_error$Batch, median)
estimated_error %>% left_join(outlier_count) %>% 
    ggplot(aes(estimated_error, N_outliers, color=Batch)) + geom_point(size=2) + scale_color_manual(values=seqtech_colors) + theme_minimal()
ggsave("estimated_error_profiles.pdf")

autosome_betas <- fread("./METAFORA_output/Population_methylation.tissue_Blood/Population_methylation.tissue_Blood.chrom_chr21.betas.mat.gz")
autosome_depths <- fread("./METAFORA_output/Population_methylation.tissue_Blood/Population_methylation.tissue_Blood.chrom_chr21.coverage.mat.gz")
auto_segs <- fread("./METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.chrom_chr21.bed")

auto_segs.gr <- makeGRangesFromDataFrame(auto_segs)
auto_segs.gr$seg_id <- 1:nrow(auto_segs)

betas.gr <- makeGRangesFromDataFrame(autosome_betas, keep.extra.columns =T)
head(auto_segs)

betas.gr %<>% join_overlap_left(auto_segs.gr)

head(betas.gr) %>% as.data.frame
betas.long <- betas.gr %>% as.data.frame() %>% select(-width,-strand) %>% pivot_longer(cols=-c(seqnames,start,end,seg_id), values_to="beta", names_to="Sample_name")
segment_summary.df <- betas.long %>% group_by(Sample_name, seg_id) %>% summarize(beta_mean=mean(beta,na.rm=T), beta_sd=sd(beta,na.rm=T))
segment_summary.df %>% filter(beta_mean<.05) %>% nrow

segment_median_sds <- segment_summary.df %>% ungroup %>% filter(beta_mean<.05) %>%  group_by(Sample_name) %>% summarize(median_sd = median(beta_sd, na.rm=T))
segment_median_sds %>% left_join(estimated_error) %>% 
    ggplot(aes(estimated_error, median_sd, color=Batch)) + geom_point(size=2) + scale_color_manual(values=seqtech_colors) + theme_minimal()
ggsave("./estimated_error_profiles.pdf")


