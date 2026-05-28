library(data.table)
library(cowplot)
library(GenomicRanges)
library(magrittr)
library(plyranges)
library(tidyverse)
library(pbmcapply)
library(Matrix)
library(matrixStats)
library(ggrepel)
library(bedr)

parser <- arg_parser("Script to compute global variation PCs in methylation profiles and auto-detect outliers")
parser <- add_argument(parser, "--seg_beta", help = "comma-separated list of summarized segment betas from segmentation for each autosome")
parser <- add_argument(parser, "--seg_depth", help = "comma-separated list of summarized segment depths from segmentation for each autosome")
parser <- add_argument(parser, "--merged_outliers", help="where to write correlation, PC, and sex plots")

args <- parse_args(parser)

methouts <- fread("../METAFORA_output/METAFORA_methylation_outlier_regions.tissue_PBMC.ALL_CHROM_COMBINED.merged_joint_called_zscore.mat.gz") %>% select(chromosome,start,end,MERGE_ID)
methouts$block <- gsub(methouts$MERGE_ID, pattern="METAFORA_\\d+_(\\w+)", replacement="\\1")
methouts$coord <- paste0(methouts$chromosome,":",methouts$start,"-",methouts$end)

this_block <- "chr1.2"
this_tissue <- "PBMC"
output_dir <- "../METAFORA_output"
calculate_hap_beta <- function(this_block, this_hp, this_tissue) {
    print(this_block)
    methout.gr <- methouts %>% dplyr::filter(block == this_block) %>% makeGRangesFromDataFrame(keep.extra.columns=T)
    betas <- fread(paste0(output_dir, "/HP_",this_hp,".tissue_",this_tissue,"/Population_methylation.hp_",this_hp,".tissue_",this_tissue,".chrom_",this_block,".betas.mat.gz"))
    betas.mat <- as.matrix(betas[,4:ncol(betas)])
    depths <- fread(paste0(output_dir, "/HP_",this_hp,".tissue_",this_tissue,"/Population_methylation.hp_",this_hp,".tissue_",this_tissue,".chrom_",this_block,".coverage.mat.gz"))
    depths.mat <- as.matrix(depths[,4:ncol(depths)])

    betas.gr <- makeGRangesFromDataFrame(betas,keep.extra.columns=F)

    ol <- findOverlaps(methout.gr, betas.gr)
    # create Segment x CpG identity matrix to indicate which cpgs belong to which segment
    CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(methout.gr),nrow(betas)), x=1)
      
    betas.mat[is.na(betas.mat)] <- 0
    depths.mat[is.na(depths.mat)] <- 0

    # aggregate betas across each segment
    methout_beta <- as.matrix(((CpG_Identity %*% (betas.mat*depths.mat))+1) / (CpG_Identity %*% depths.mat+2))
    methout_depth <- as.matrix((CpG_Identity %*% depths.mat+2) / rowSums(CpG_Identity))
    methout_beta[methout_depth<5] <- NA #low depth segments set to NA
    data.frame(methout_beta)
}

hap1_betas <- pbmclapply(unique(methouts$block), function(b) calculate_hap_beta(b,1,this_tissue), mc.cores=8) %>% bind_rows 
hap2_betas <- pbmclapply(unique(methouts$block), function(b) calculate_hap_beta(b,2,this_tissue), mc.cores=8) %>% bind_rows
hap1_betas %<>% as.matrix
hap2_betas %<>% as.matrix
hap_delta <- as.matrix(hap1_betas - hap2_betas)
median_hap_delta <- rowMedians(abs(hap_delta),na.rm=T)
rownames(hap_delta) <- methouts$MERGE_ID

write.table(hap_delta, file="../METAFORA_output/METAFORA_methylation_outlier_regions.tissue_PBMC.ALL_CHROM_COMBINED.merged_joint_called_zscore.haplotype_delta.mat", row.names=T, col.names=T, quote=F)

hap_diff.df <- data.frame(METH_ID=rownames(hap_delta), median_hap_delta)
hap_diff_z <- (abs(hap_delta) - rowMeans(abs(hap_delta),na.rm=T))/(rowSds(abs(hap_delta),na.rm=T))
hap_diff.df <- hap_diff_z %>% melt %>% set_colnames(c("METH_ID","sample","hap_diff_z")) %>% left_join(hap_diff.df)

hap_diff_outlier_count <- hap_diff.df %>% group_by(METH_ID) %>% summarize(HAPDELTA_OUTLIER_COUNT=sum(hap_diff_z < -3,na.rm=T))
hap_diff.df %<>% left_join(hap_diff_outlier_count)
hap_diff.df %>% arrange(hap_diff_z) %>% head
hap_diff.df %>% arrange(hap_diff_z) %>% filter(HAPDELTA_OUTLIER_COUNT<3) %>% head

rownames(hap1_betas) <- methouts$MERGE_ID
rownames(hap2_betas) <- methouts$MERGE_ID
hap_diff.df %<>% left_join(hap1_betas%>%melt%>%set_colnames(c("METH_ID","sample","hap1_beta"))) %>% left_join(hap2_betas%>%melt%>%set_colnames(c("METH_ID","sample","hap2_beta")))

hap_diff.df %>% arrange(hap_diff_z) %>% head
methouts %>% filter(MERGE_ID=="METAFORA_105_chr11.1")
methouts %>% filter(MERGE_ID=="METAFORA_64_chr15.1")

missing_hap_seg_burden <- hap_diff.df %>% group_by(sample) %>% summarize(n=sum(hap_diff_z < -5, na.rm=T)) %>% arrange(desc(n))
mars.meta <- fread("../../MARS_cohort.nanostat.sequencing_summary_stats_table.txt")
mars.meta <- fread("../../MARS_cohort.combined_seqstats_demo_covariates.sample_table.txt") %>% rename(sample=Sample)
colnames(mars.meta)
ggplot(missing_hap_seg_burden %>% mutate(high=n>5) %>% left_join(mars.meta), aes(high, Age)) +geom_jitter(width=.3)
ggplot(missing_hap_seg_burden %>% mutate(high=n>5) %>% left_join(mars.meta), aes(high, Mono)) +geom_jitter(width=.3)
ggplot(missing_hap_seg_burden %>% mutate(high=n>5) %>% left_join(mars.meta), aes(high, CD4T)) +geom_jitter(width=.3)
ggplot(missing_hap_seg_burden %>% mutate(high=n>5) %>% left_join(mars.meta), aes(high, CD8T)) +geom_jitter(width=.3)
ggplot(missing_hap_seg_burden %>% mutate(high=n>5) %>% left_join(mars.meta), aes(high, Bcell)) +geom_jitter(width=.3)
ggplot(missing_hap_seg_burden %>% mutate(high=n>5) %>% left_join(mars.meta), aes(high, AFR)) +geom_jitter(width=.3)
ggplot(missing_hap_seg_burden %>% mutate(high=n>5) %>% left_join(mars.meta), aes(high, AFR)) +geom_jitter(width=.3)
missing_hap_seg_burden %>% filter(n>5) %>% left_join(mars.meta) 
ggsave("./tmp.pdf")

high_samples <- missing_hap_seg_burden %>% filter(n>5) %>% pull(sample)

hap_diff.df %>% filter(hap_diff_z < -5) %>% filter(sample %in% high_samples) %>% pull(median_hap_delta) %>% unique
    
methouts %>% filter(MERGE_ID=="METAFORA_281_chr19.2")
methouts %>% filter(MERGE_ID=="METAFORA_282_chr19.2")
 
ggplot(hap_diff.df)
summary(mars.meta$Age)

hp1_segs <- lapply(list.files("../METAFORA_output/HP_1.tissue_PBMC/", pattern="Meth_segments.hp_1.tissue_PBMC.segment_betas.chrom_.*.bed", full.names=T), fread) %>% bind_rows
hp2_segs <- lapply(list.files("../METAFORA_output/HP_2.tissue_PBMC/", pattern="Meth_segments.hp_2.tissue_PBMC.segment_betas.chrom_.*.bed", full.names=T), fread) %>% bind_rows
sex_chrom_segs <- hp1_segs$chrom %in% c("chrX","chrY")
head(hp1_segs[,1:5])
hp1_seg_betas <- as.matrix(hp1_segs[,5:ncol(hp1_segs)])
rownames(hp1_seg_betas) <- hp1_segs$seg_id
hp2_seg_betas <- as.matrix(hp2_segs[,5:ncol(hp2_segs)])
rownames(hp2_seg_betas) <- hp2_segs$seg_id
segs.meta <- hp1_segs %>% select(chrom,start,end,seg_id)
segs.meta %<>% mutate(block=gsub(seg_id,pattern="(\\w+\\.\\d+)_\\d+_\\d+",replacement="\\1"),coord=paste0(chrom,":",start,"-",end))

hp1_seg_cov <- lapply(list.files("../METAFORA_output/HP_1.tissue_PBMC/", pattern="Meth_segments.hp_1.tissue_PBMC.segment_coverage.chrom_.*.mat", full.names=T), fread) %>% bind_rows %>% as.matrix
hp2_seg_cov <- lapply(list.files("../METAFORA_output/HP_2.tissue_PBMC/", pattern="Meth_segments.hp_2.tissue_PBMC.segment_coverage.chrom_.*.mat", full.names=T), fread) %>% bind_rows %>% as.matrix

hp1_seg_betas[hp1_seg_cov<5] <- NA
hp2_seg_betas[hp2_seg_cov<5] <- NA


percent_missing <- rowSums(is.na(seg_hap_diff))/ncol(seg_hap_diff)
summary(percent_missing)

seg_hap_diff <- abs(hp1_seg_betas - hp2_seg_betas)
seg_hap_diff <- seg_hap_diff[!sex_chrom_segs,] #subset to autosomes
seg_hap_diff <- seg_hap_diff[percent_missing<.2,]
hap_diff_medians <- rowMedians(seg_hap_diff,na.rm=T)
hap_diff_q1 <- apply(seg_hap_diff,1, function(x) quantile(x,.01,na.rm=T)) #get 1 percentile of haplotype difference (99% of samples have hap delta > this)
sum(hap_diff_q1>.5,na.rm=T) #highly "constrained" imprinting loci

ggplot(NULL, aes(hap_diff_q1, fill=hap_diff_q1>.5)) + geom_histogram(alpha=.8) + theme_minimal() + scale_fill_manual(values=c("grey","red")) + scale_y_log10() + 
    geom_vline(xintercept=.5,color="red",linetype="dashed")
ggsave("seg.median_hap_diff.distribution.hist.pdf",height=3)

imprint_cand_regions <- rownames(seg_hap_diff)[hap_diff_q1>.5]
imprint_cand_regions <- imprint_cand_regions[!is.na(imprint_cand_regions)]
names(hap_diff_q1) <- rownames(seg_hap_diff)
imprint_cand_regions <- imprint_cand_regions[order(-rowMeans(seg_hap_diff[imprint_cand_regions,],na.rm=T))]

cand_imprint_delta <- seg_hap_diff[imprint_cand_regions,] %>% melt %>% set_colnames(c("meth_id","sample","hap_delta"))
cand_imprint_delta$meth_id %<>% factor(levels=imprint_cand_regions)

cand_imprint_delta %<>% group_by(meth_id) %>% mutate(hap_diff_z=scale(hap_delta)) %>% ungroup 
outburden <- cand_imprint_delta %>% group_by(sample) %>% summarize(n=sum(hap_diff_z < -3,na.rm=T)) %>% arrange(desc(n))
imprinting_outlier_samps <- outburden %>% filter(n>2) %>% pull(sample)
cand_imprint_delta %<>% mutate(global_imprinting_outlier = (sample %in% imprinting_outlier_samps))
outburden
ggplot(cand_imprint_delta, aes(hap_delta, meth_id, color=sample=="RAD38504640",size=sample=="RAD38504640")) + geom_jitter(height=.3,alpha=.8) + theme_minimal() + scale_color_manual(values=c("black","red")) + 
    geom_vline(xintercept=.5,linetype="dashed",color="red") + scale_size_manual(values=c(1,2))
ggsave("tmp.pdf",height=15)

outburden$rank <- nrow(outburden):1
ggplot(outburden, aes(rank,n,color=n>10,label=sample)) + geom_point() + geom_text_repel() + theme_classic() + scale_color_manual(values=c("black","red"))
ggsave('imprinitng_missegregation_outliers.pdf')

m<-"chr20.1_966_61"
samp <- "RAD38504640"
plot_hap_region <- function(m,samp) {
    this_seg <- segs.meta %>% filter(seg_id==m)
    this_block <- this_seg$block
    this_region  <- this_seg$coord
    this_chrom <- gsub(this_block,pattern="(\\w+)\\.\\d+", replacement="\\1")
    hp1 <- data.table(tabix(this_region, paste0('../METAFORA_output/HP_1.tissue_PBMC/Population_methylation.hp_1.tissue_PBMC.chrom_',this_block,'.betas.mat.gz')))
    hp2 <- data.table(tabix(this_region, paste0('../METAFORA_output/HP_2.tissue_PBMC/Population_methylation.hp_2.tissue_PBMC.chrom_',this_block,'.betas.mat.gz')))
    hp1.mat <- hp1[,4:ncol(hp1)] %>% set_colnames(colnames(hp1_seg_betas)) %>% as.matrix
    hp2.mat <- hp2[,4:ncol(hp2)] %>% set_colnames(colnames(hp2_seg_betas)) %>% as.matrix
    rownames(hp1.mat) <- hp1$V2
    rownames(hp2.mat) <- hp2$V2
    hp_meth <- rbind(
        hp1.mat %>% melt %>% set_colnames(c("cpg_pos","sample","beta")) %>% mutate(hp="hp1"),
        hp2.mat %>% melt %>% set_colnames(c("cpg_pos","sample","beta")) %>% mutate(hp="hp2")
    ) %>% as.data.frame %>% mutate(cpg_pos=as.integer(cpg_pos),beta=as.numeric(beta))
    hp_meth$sample %<>% as.character
    hp_mean <- hp_meth %>% group_by(sample,hp) %>% summarize(missing=sum(is.na(beta))/n(), beta=mean(beta,na.rm=T)) %>% ungroup
    samples_keep <- hp_mean %>% filter(missing<.1,!is.na(beta)) %>% group_by(sample) %>% summarize(n=n()) %>% filter(n==2) %>% pull(sample)
    hp_transform <- hp_mean %>% filter(sample %in% samples_keep) %>% group_by(sample) %>% summarize(max_hap=hp[which.max(beta)], min_hap=hp[which.min(beta)]) %>% pivot_longer(-sample,names_to="hap_cat", values_to="hp")
    hp_meth %<>% left_join(hp_transform) %>% filter(!is.na(hap_cat))
    hp_meth %<>% group_by(sample,hp,hap_cat) %>% mutate(beta.smoothed = sapply(1:n(), function(x) mean(beta[max(1,x-3):min(x+3,n())], na.rm=T)))

    set.seed(123)
    subset_samples <- c(sample(unique(as.character(hp_meth$sample)),50),samp)
    lineplt.by_hap <- ggplot(hp_meth%>%filter(sample %in% subset_samples), aes(cpg_pos,beta.smoothed,shape=hap_cat,group=interaction(sample,hap_cat), color=sample==samp,alpha=sample==samp)) + 
        geom_line() + facet_wrap(~hap_cat,ncol=1) + theme_minimal() + scale_color_manual(values=c("black","red")) + scale_alpha_manual(values=c(.1,.9)) +
        ylab("Methylation beta") + xlab(paste0(this_chrom, " position")) + ylim(0,1) +
        theme(legend.position="none")
    sample.plt <- ggplot(hp_meth%>%filter(sample %in% samp), aes(cpg_pos,beta.smoothed,color=hap_cat)) + 
        geom_line() + theme_minimal() + scale_color_manual(values=c("firebrick3","steelblue3")) + scale_alpha_manual(values=c(.1,.9)) + ylim(0,1) +
        ylab("Methylation beta") + xlab(paste0(this_chrom, " position")) + theme(legend.position="top")

    class(hp1.mat) <- "numeric"
    class(hp2.mat) <- "numeric"
    hap_diff.mat <- abs(hp1.mat-hp2.mat)
    hap_diff.df <- hap_diff.mat %>% melt %>% set_colnames(c("cpg_pos","sample","hap_delta"))
    hap_diff.df %<>% group_by(sample) %>% mutate(hap_delta.smoothed = sapply(1:n(), function(x) mean(hap_delta[max(1,x-3):min(x+3,n())],na.rm=T)))
    hap_diff.lineplt <- ggplot(hap_diff.df %>% filter(sample %in% subset_samples), aes(cpg_pos,hap_delta.smoothed,group=sample,color=sample==samp,alpha=sample==samp)) + 
        geom_line() + theme_minimal() + scale_color_manual(values=c("black","red")) + scale_alpha_manual(values=c(.1,.9)) + 
        ylab("| Haplotype delta |") + xlab(paste0(this_chrom, " position"))  + theme(legend.position="none") + ylim(0,1) 
    hap_diff.hist <- cand_imprint_delta %>% filter(meth_id==m) %>%
        ggplot(aes(hap_delta, fill=sample==samp)) + geom_histogram(bins=50) + theme_minimal() + scale_fill_manual(values=c("grey","red")) + 
        geom_vline(xintercept=cand_imprint_delta%>%filter(meth_id==m,sample==samp)%>%pull(hap_delta), color="red", linetype="dashed") +
        theme(legend.position="none") + xlim(0,1)
    mean.meth.plt <- hp_meth %>% group_by(sample,cpg_pos) %>% summarize(mean_meth=mean(beta.smoothed)) %>% filter(sample %in% subset_samples) %>%
        ggplot(aes(cpg_pos,mean_meth,group=sample,color=sample==samp,alpha=sample==samp)) +
        geom_line() + theme_minimal() + scale_color_manual(values=c("black","red")) + scale_alpha_manual(values=c(.1,.9)) +
        ylab("Methylation beta") + xlab(paste0(this_chrom, " position")) + ylim(0,1) + labs(title=paste0("sample: ",samp), subtitle=paste0("region: ", this_region)) +
        theme(legend.position="none",plot.title = element_text(color = "red"))
    topleft <- plot_grid(mean.meth.plt,hap_diff.lineplt,ncol=1)
    top <- plot_grid(topleft,lineplt.by_hap, nrow=1)
    bottom <- plot_grid(sample.plt,hap_diff.hist,nrow=1)
    plot_grid(top,bottom,rel_heights=c(2,1),ncol=1)
    ggsave(paste0("../constrained_imprinting_segregation_loss/",samp,".",this_region,".hap_delta_plots.pdf"),height=13, width=9)
}


cands.region <- cand_imprint_delta %>% filter(hap_diff_z < -3)
head(cands.region)
lapply(1:nrow(cands.region), function(i) plot_hap_region(m=cands.region$meth_id[i],samp=cands.region$sample[i]))
