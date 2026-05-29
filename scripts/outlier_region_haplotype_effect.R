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

#calculate_hap_beta <- function(this_block, this_hp, this_tissue) {
#    print(this_block)
#    methout.gr <- methouts %>% dplyr::filter(block == this_block) %>% makeGRangesFromDataFrame(keep.extra.columns=T)
#    betas <- fread(paste0(output_dir, "/HP_",this_hp,".tissue_",this_tissue,"/Population_methylation.hp_",this_hp,".tissue_",this_tissue,".chrom_",this_block,".betas.mat.gz"))
#    betas.mat <- as.matrix(betas[,4:ncol(betas)])
#    depths <- fread(paste0(output_dir, "/HP_",this_hp,".tissue_",this_tissue,"/Population_methylation.hp_",this_hp,".tissue_",this_tissue,".chrom_",this_block,".coverage.mat.gz"))
#    depths.mat <- as.matrix(depths[,4:ncol(depths)])
#
#    betas.gr <- makeGRangesFromDataFrame(betas,keep.extra.columns=F)
#
#    ol <- findOverlaps(methout.gr, betas.gr)
#    # create Segment x CpG identity matrix to indicate which cpgs belong to which segment
#    CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(methout.gr),nrow(betas)), x=1)
#      
#    betas.mat[is.na(betas.mat)] <- 0
#    depths.mat[is.na(depths.mat)] <- 0
#
#    # aggregate betas across each segment
#    methout_beta <- as.matrix(((CpG_Identity %*% (betas.mat*depths.mat))+1) / (CpG_Identity %*% depths.mat+2))
#    methout_depth <- as.matrix((CpG_Identity %*% depths.mat+2) / rowSums(CpG_Identity))
#    methout_beta[methout_depth<5] <- NA #low depth segments set to NA
#    data.frame(methout_beta)
#}
#
#hap1_betas <- pbmclapply(unique(methouts$block), function(b) calculate_hap_beta(b,1,this_tissue), mc.cores=8) %>% bind_rows 
#hap2_betas <- pbmclapply(unique(methouts$block), function(b) calculate_hap_beta(b,2,this_tissue), mc.cores=8) %>% bind_rows
#hap1_betas %<>% as.matrix
#hap2_betas %<>% as.matrix
#hap_delta <- as.matrix(hap1_betas - hap2_betas)
#median_hap_delta <- rowMedians(abs(hap_delta),na.rm=T)
#rownames(hap_delta) <- methouts$MERGE_ID
#write.table(hap_delta, file="../METAFORA_output/METAFORA_methylation_outlier_regions.tissue_PBMC.ALL_CHROM_COMBINED.merged_joint_called_zscore.haplotype_delta.mat", row.names=T, col.names=T, quote=F)



cand_imprint_segs <- lapply(list.files("../METAFORA_output/imprinting_regions.tissue_PBMC",pattern="*.mat",full.names=T),fread) %>% bind_rows
blocks<-NULL
for (file in list.files("../METAFORA_output/imprinting_regions.tissue_PBMC",pattern="*.mat",full.names=T)) {
    block <- gsub(file,pattern=".*\\.chrom_(\\w+\\.\\d+)\\..*",replacement="\\1")
    tmp <- fread(file)
    blocks <- c(blocks,rep(block,nrow(tmp)))
}

segs.meta <- cand_imprint_segs[,1:4]
segs.meta$block <- blocks

hap_diff.mat <- as.matrix(cand_imprint_segs[,5:ncol(cand_imprint_segs)])
rownames(hap_diff.mat) <- cand_imprint_segs$seg_id
percent_missing <- rowSums(is.na(hap_diff.mat))/ncol(hap_diff.mat)
hap_diff.mat <- hap_diff.mat[percent_missing<.2,]

hap_diff_q1 <- apply(hap_diff.mat,1, function(x) quantile(x,.02,na.rm=T)) #get 1 percentile of haplotype difference (99% of samples have hap delta > this)
sum(hap_diff_q1>.5,na.rm=T) #highly "constrained" imprinting loci

ggplot(NULL, aes(hap_diff_q1, fill=hap_diff_q1>.25)) + geom_histogram(alpha=.8) + theme_minimal() + scale_fill_manual(values=c("grey","red")) + scale_y_log10() + 
    geom_vline(xintercept=.5,color="red",linetype="dashed")
ggsave("seg.median_hap_diff.distribution.hist.pdf",height=3)


hap_diff.mat <- hap_diff.mat[hap_diff_q1>.5,]
cand_imprint_delta <- hap_diff.mat %>% melt %>% set_colnames(c("meth_id","sample","hap_delta"))
imprint_cand_regions <- rownames(hap_diff.mat)[order(rowMeans(hap_diff.mat,na.rm=T))]
cand_imprint_delta$meth_id %<>% factor(levels=imprint_cand_regions)
cand_imprint_delta$sample %<>% as.character

cand_imprint_delta %<>% group_by(meth_id) %>% mutate(hap_diff_z=scale(hap_delta)) %>% ungroup 
outburden <- cand_imprint_delta %>% group_by(sample) %>% summarize(n=sum(hap_diff_z < -5,na.rm=T)) %>% arrange(desc(n))
imprinting_outlier_samps <- outburden %>% filter(n>5) %>% pull(sample)
cand_imprint_delta %<>% mutate(global_outlier=sample%in% imprinting_outlier_samps, outlier_sample=ifelse(global_outlier, yes=sample, no="not global outlier"))

outburden
ggplot(cand_imprint_delta, aes(hap_delta, meth_id, fill=outlier_sample,size=global_outlier,shape=global_outlier)) + geom_jitter(height=.3,alpha=.8,color="black") + theme_minimal() + scale_fill_manual(values=c("black","firebrick2","steelblue2","forestgreen","goldenrod")) + geom_hline(yintercept=.5+1:nrow(hap_diff.mat)) + theme(panel.border=element_rect(fill=NA)) +
    scale_size_manual(values=c(1,2)) + scale_shape_manual(values=c(21,23))
ggsave("tmp.pdf",height=12)

outburden %<>% mutate(global_outlier=sample%in% imprinting_outlier_samps, outlier_sample=ifelse(global_outlier, yes=sample, no="not global outlier"))
outburden$rank <- nrow(outburden):1
set.seed(123)
ggplot(outburden, aes(rank,n,fill=outlier_sample,label=ifelse(global_outlier,sample,""),shape=global_outlier)) + geom_point(size=2,alpha=.8,color="black") + geom_text_repel() + theme_classic() + 
    scale_fill_manual(values=c("black","firebrick2","steelblue2","forestgreen","goldenrod")) + scale_shape_manual(values=c(21,23)) +
    geom_hline(yintercept=5,color="red",linetype="dashed") +
    ylab("number of imprinting | haplotype delta |  outliers per genome") +
    xlab("rank-order samples")
ggsave('imprinitng_missegregation_outliers.pdf')

plot_hap_region <- function(m,samp) {
    this_seg <- segs.meta %>% filter(seg_id==m)
    this_block <- this_seg$block
    this_region  <- this_seg$seg_id
    this_chrom <- gsub(this_block,pattern="(\\w+)\\.\\d+", replacement="\\1")
    hp1 <- data.table(tabix(this_region, paste0('../METAFORA_output/HP_1.tissue_PBMC/Population_methylation.hp_1.tissue_PBMC.chrom_',this_block,'.betas.mat.gz')))
    hp2 <- data.table(tabix(this_region, paste0('../METAFORA_output/HP_2.tissue_PBMC/Population_methylation.hp_2.tissue_PBMC.chrom_',this_block,'.betas.mat.gz')))
    hp1.mat <- hp1[,4:ncol(hp1)] %>% set_colnames(colnames(hap_diff.mat)) %>% as.matrix
    hp2.mat <- hp2[,4:ncol(hp2)] %>% set_colnames(colnames(hap_diff.mat)) %>% as.matrix
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


cands.region <- cand_imprint_delta %>% filter(hap_diff_z < -5)
lapply(1:nrow(cands.region), function(i) plot_hap_region(m=cands.region$meth_id[i],samp=cands.region$sample[i]))
