library(tidyverse)
library(plyranges)
library(ggbreak)
library(rtracklayer)
library(ggrepel)
library(data.table)
library(magrittr)
library(ggbeeswarm)

outliers <- fread("./METAFORA_output/METAFORA_methylation_outlier_regions.tissue_Blood.ALL_CHROM_COMBINED.haplotype_annotated.gene_track_annotated.bed")
select_cols=c("seg_id","width","num.mark","delta","zscore","Tissue","CHROM_TYPE","ID","PRIORITIZATION_SCORE","combined_depth","haplotype_coverage_bias","hap_delta","ImprintDisordDMR","CpgIslands","ATAC_peaks","ABC_Enhancers", "Promoter")
blood.outs <- fread("./METAFORA_output/METAFORA_methylation_outlier_regions.tissue_Blood.ALL_CHROM_COMBINED.haplotype_annotated.gene_track_annotated.bed") %>% select(all_of(select_cols)) %>% unique
fibro.outs <-fread("./METAFORA_output/METAFORA_methylation_outlier_regions.tissue_Fibroblast.ALL_CHROM_COMBINED.haplotype_annotated.gene_track_annotated.bed") %>% select(all_of(select_cols)) %>% unique
outs <- rbind(blood.outs, fibro.outs)
outs %>% filter(num.mark>30)
prop.table(table(outs$delta > 0))
prop.table(table(outs$zscore > 0))
outs %>% filter(ImprintDisordDMR) %>% pull(unique(ID)) %>% length
outs %>% filter(ImprintDisordDMR)  %>% pull(ID) %>% unique %>% length


meta.dat <- fread("./U08_GREGoR.covariates.txt") %>% dplyr::rename(Technology=Batch)
site_map=list("GS"="GSS", "M1"="Broad", "PM"="PMGRC", "UW"="UW", "UD"="GSS")
meta.dat %<>% mutate(Site=as.character(sapply(substr(Sample_name,1,2),function(x) site_map[x])))
meta.dat$ID <- meta.dat$Sample_name

meta.dat$Sample_name <- NULL

depths.blood <- fread("./METAFORA_output/Global_Methylation_PCA_tissue_Blood/PCA_covariates.txt")  %>% dplyr::rename(ID=Sample_name) %>% select(ID,depth)
depths.fibro <- fread("./METAFORA_output/Global_Methylation_PCA_tissue_Fibroblast/PCA_covariates.txt")  %>% dplyr::rename(ID=Sample_name) %>% select(ID,depth)
covars.blood <- fread("./METAFORA_output/Global_Methylation_PCA_tissue_Blood/PCA_covariates.txt")  %>% dplyr::rename(ID=Sample_name)
covars.fibro <- fread("./METAFORA_output/Global_Methylation_PCA_tissue_Fibroblast/PCA_covariates.txt")  %>% dplyr::rename(ID=Sample_name)
depths <- rbind(depths.blood, depths.fibro)
meta.dat %<>% left_join(depths)


metadat <- fread("./U08_GREGoR.sequencing_stats.metadata.txt")
table(metadat$Technology)
nrow(metadat)
ggplot(metadat,aes(coverage, n50, shape=Site, fill=Technology)) + geom_point(size=2,alpha=.6,color="black") +
    scale_fill_manual(values=c("turquoise","blue","hotpink")) +
    scale_shape_manual(values=c(21,22,23,24)) +
    theme_minimal()
ggsave("U08.seq_stats.summary.pdf",width=5,height=4)
ggplot(covars.blood%>%mutate(depth=pmax(pmin(depth,40),10)), aes(PC1,-PC2,shape=Batch,fill=depth)) + 
    geom_point(size=2.8,color="black",alpha=.8) + theme_minimal() + scale_shape_manual(values=c(ONT_R10=21,ONT_R9=24,PacBio=22)) +
    scale_fill_viridis_c(option = "inferno",end=.9) +
    labs(fill = "Median sequencing depth (×)", shape="Sequencing Technology")
ggsave("Methylation_biplot.technology.pdf")

summary(lm(PC2~depth,data=covars.blood))$coefficients
cor.test(-covars.blood$PC2,covars.blood$depth, method="pearson")
ggplot(covars.blood, aes(depth,-PC2)) + geom_point(aes(shape=Batch,fill=Batch),size=2,color="black",alpha=.8) +
    geom_smooth(method="lm",color="red") +
    theme_classic() +
    scale_shape_manual(values=c(ONT_R10=21,ONT_R9=24,PacBio=22)) +
    scale_fill_manual(values=c("turquoise","blue","hotpink"))
ggsave("Methylation.PC_depth_correlation.pdf", width=4,height=3)

outs %>% left_join(meta.dat) %>% filter(Tissue=="Blood")  %>% filter(CHROM_TYPE=="AUTOSOME",num.mark>30) %>% 
    group_by(ID,Site,Technology) %>% summarize(n=dplyr::n()) %>% 
    ggplot(aes(Site,n,fill=Technology)) + geom_boxplot() + scale_y_log10() + facet_grid(.~Technology) +
    theme_minimal() + theme(panel.border=element_rect(fill=NA),legend.position="none", axis.text.x=element_text(angle=30,hjust=1,vjust=1)) + 
    scale_fill_manual(values=c("turquoise", "blue", "hotpink")) +  
    ylab("Number of Outliers per Sample")
ggsave("./U08.outliers_by_tissue_technology.autosomes.pdf", height=5, width=3)


outliers <- outs %>% left_join(meta.dat) %>% filter(num.mark>30)
ggplot(outliers, aes(delta, abs(zscore), shape=Technology, fill=Tissue)) + geom_point(alpha=.4,size=2,color="black") +
    theme_minimal() + theme(panel.border=element_rect(fill=NA)) + 
    geom_vline(xintercept=-.3, linetype="dashed", color="grey") +
    geom_vline(xintercept=.3, linetype="dashed", color="black") +
    scale_fill_manual(values=c("brown", "navajowhite3")) + 
    facet_wrap(.~Tissue,ncol=1,scales="free_y") +
    scale_shape_manual(values=c(21,24,22)) +
    theme(legend.position="none")
ggsave("outliers.volcano.png",width=4,height=5)

outliers %>% group_by(Tissue, delta>0) %>% summarize(dplyr::n()/sum(outliers$Tissue==Tissue[1]))

num_outliers <- outliers %>% filter(CHROM_TYPE=="AUTOSOME") %>% group_by(ID,Technology,Site,depth,Tissue) %>% summarize(n=dplyr::n(), hyper=sum(zscore>0), hypo=sum(zscore<0), median_length=median(width), median_cpgs=median(num.mark))
outliers
num_outliers %>% pull(n) %>% median
outliers

outliers %>% filter(num.mark>30) %>% group_by(ID,Technology) %>% summarize(n=dplyr::n()) %>% ungroup %>% group_by(Technology) %>% summarize(median(n))
outliers %>% filter(Technology=="PacBio", abs(hap_delta)>.25) %>% group_by(ID) %>% summarize(n=dplyr::n()) %>% pull(n) %>% median

summary(lm(abs(haplotype_coverage_bias) ~ Technology,data=outliers%>%filter(is.finite(haplotype_coverage_bias), Technology != "ONT_R9")))
summary(lm(num.mark ~ Technology,data=outliers%>%filter(is.finite(haplotype_coverage_bias), Technology != "ONT_R9")))
summary(lm(width ~ Technology,data=outliers%>%filter(is.finite(haplotype_coverage_bias), Technology != "ONT_R9")))
summary(lm(abs(delta) ~ Technology,data=outliers%>%filter(is.finite(haplotype_coverage_bias), Technology != "ONT_R9")))
summary(lm(CpgIslands ~ Technology,data=outliers%>%filter(is.finite(haplotype_coverage_bias), Technology != "ONT_R9")))
summary(lm(ABC_Enhancers ~ Technology,data=outliers%>%filter(is.finite(haplotype_coverage_bias), Technology != "ONT_R9")))
summary(lm(ATAC_peaks ~ Technology,data=outliers%>%filter(is.finite(haplotype_coverage_bias), Technology != "ONT_R9")))
summary(lm(PRIORITIZATION_SCORE ~ Technology,data=outliers%>%filter(is.finite(haplotype_coverage_bias), Technology != "ONT_R9")))

#Prioritization Plot
global_outs <- priot_count %>% filter(n>100) %>% pull(ID) %>% unique
priot_count <- rbind( 
    outliers %>% group_by(ID,Technology) %>% summarize(n=dplyr::n()) %>% mutate(filter="All Outliers"),
    outliers %>% filter(abs(delta)>.4) %>% group_by(ID,Technology) %>% summarize(n=dplyr::n()) %>% mutate(filter="| delta | > 0.4"),
    outliers %>% filter(abs(delta)>.5) %>% group_by(ID,Technology) %>% summarize(n=dplyr::n()) %>% mutate(filter="| delta | > 0.5"),
    outliers %>% filter(abs(hap_delta)>.25) %>% group_by(ID,Technology) %>% summarize(n=dplyr::n()) %>% mutate(filter = "| haplotype delta | > 0.25"),
    outliers %>% filter(abs(hap_delta)>.5) %>% group_by(ID,Technology) %>% summarize(n=dplyr::n()) %>% mutate(filter = "| haplotype delta | > 0.5"),
    outliers %>% filter(combined_depth>10, abs(haplotype_coverage_bias) < 1) %>% group_by(ID,Technology) %>% summarize(n=dplyr::n()) %>% mutate(filter= "High-confidence (Depth > 10, balanced hapltoype coverage)"),
    outliers %>% filter(CpgIslands) %>% group_by(ID,Technology) %>% summarize(n=dplyr::n()) %>% mutate(filter="Overlap CpG Island"),
    outliers %>% filter(ATAC_peaks) %>% group_by(ID,Technology) %>% summarize(n=dplyr::n())  %>% mutate(filter="Overlap PBMC ATAC Peak"),
    outliers %>% filter(ABC_Enhancers) %>% group_by(ID,Technology) %>% summarize(n=dplyr::n()) %>% mutate(filter="Overlap ABC Enhancer"),
    outliers %>% filter(Promoter) %>% group_by(ID,Technology) %>% summarize(n=dplyr::n()) %>% mutate(filter="Overlap Promoter"),
    outliers %>% filter(ImprintDisordDMR) %>% group_by(ID,Technology) %>% summarize(n=dplyr::n()) %>% mutate(filter="Overlap Imprinting Disorder DMR"),
    outliers %>% filter(PRIORITIZATION_SCORE>3) %>% group_by(ID,Technology) %>% summarize(n=dplyr::n()) %>% mutate(filter="Prioritization Score >= 3"),
    outliers %>% filter(PRIORITIZATION_SCORE>5) %>% group_by(ID,Technology) %>% summarize(n=dplyr::n()) %>% mutate(filter="Prioritization Score >= 5")
)
priot_count$filter %<>% factor(levels=rev(unique(priot_count$filter)))
priot_count.bar <- priot_count %>% group_by(filter) %>% summarize(median_count=median(n), mean_count=mean(n), n_samples=dplyr::n())
priot_count.bar

ggplot(priot_count %>% filter(!ID %in% global_outs), aes(n,filter)) +
    geom_jitter(aes(color=Technology),height=.3,alpha=.5) +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "crossbar", width = 0.5, color="black") +
    scale_color_manual(values=c("turquoise", "blue", "hotpink")) +  
    theme_minimal()  +
    xlab("Number of outliers per sample meeting filter") +
    xlab("Prioritization Filter")
ggsave("Prioritization_beeswarm_plot.pdf",width=10)
